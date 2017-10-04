import time
from datetime import datetime

import mmap
import os
from Bio.Seq import Seq
from Bio.Alphabet import generic_dna, generic_protein
from mutalyzer.GenRecord import PList, Locus, Gene, Record
from mutalyzer.dbgb.models import Transcript, Reference

from mutalyzer.config import settings


def get_sequence_mmap(file_path, start, end):
    """
    Sequence retrieval.
    :param file_path: Path towards the sequence file.
    :param start: Start position.
    :param end: End position
    :return: The sequence.
    """
    with open(file_path, "r+b") as f:
        # memory-map the file, size 0 means whole file
        mm = mmap.mmap(f.fileno(), 0)
        return mm[start - 1:end]


def get_reference(accession, version):
    """
    Retrieves the database reference entry for the user provided accession and
    version.
    :param accession: The accession from the user provided description.
    :param version: The version number from the user provided description.
    :return: Database reference entry.
    """
    reference = Reference.query.filter_by(accession=str(accession), version=str(version)) \
        .order_by(Reference.id.asc()) \
        .first()
    return reference


def get_mutalyzer_record(reference, db_transcripts):
    """
    Creates a Mutalyzer specific record from the transcript entries retrieved
    from the gbparser database.
    :param reference: A gbparser database reference entry.
    :param db_transcripts:A gbparser database list of transcript.
    :return: The Mutalyzer record.
    """
    record = Record()
    # Populating the record with the generic information.
    record.source_id = reference.accession
    record.id = reference.accession
    record.source_accession = reference.accession
    record.source_version = reference.version
    record.organism = 'Homo sapiens'
    record.molType = 'g'

    # Extracting the transcripts from the DB entries.
    transcripts = []
    for transcript in db_transcripts:
        my_transcript = {
            'gene': transcript.gene,
            'strand': transcript.strand,
            'transcript_start': transcript.transcript_start,
            'transcript_stop': transcript.transcript_stop,
            'cds_start': transcript.cds_start,
            'cds_stop': transcript.cds_stop,
            'exons': [],
            'transcriptID': transcript.transcript_accession + '.' + transcript.transcript_version,
            'proteinID': transcript.protein_accession + '.' + transcript.protein_version,
            'linkMethod': 'ncbi'
        }
        if transcript.exons and isinstance(transcript.exons, list):
            for exon in transcript.exons:
                exon = {'start': exon.start,
                        'stop': exon.stop}
                my_transcript['exons'].append(exon)
        transcripts.append(my_transcript)

    print("%s: End transcripts formation" % datetime.now())

    # Generating the actual record entries in the Mutalyzer format.
    gene_dict = {}
    for transcript in transcripts:
        if transcript['gene'] in gene_dict:
            gene = gene_dict[transcript['gene']]
        else:
            gene = Gene(transcript['gene'])

        if transcript['strand'] == '+':
            gene.orientation = 1
        if transcript['strand'] == '-':
            gene.orientation = -1

        version = gene.newLocusTag()
        my_transcript = Locus(version)

        my_transcript.mRNA = PList()
        my_transcript.mRNA.location = [transcript['transcript_start'],
                                       transcript['transcript_stop']]
        my_transcript.CDS = PList()
        my_transcript.CDS.location = [transcript['cds_start'],
                                      transcript['cds_stop']]
        if transcript.get('exons') and isinstance(transcript.get('exons'), list):
            for exon in transcript['exons']:
                my_transcript.CDS.positionList.append(exon['start'])
                my_transcript.CDS.positionList.append(exon['stop'])
        my_transcript.transcriptID = transcript['transcriptID']
        my_transcript.proteinID = transcript['proteinID']
        my_transcript.linkMethod = 'ncbi'
        my_transcript.transcribe = True

        gene.transcriptList.append(my_transcript)
        gene_dict[gene.name] = gene

    record.geneList = list(gene_dict.values())

    # Get the sequence.
    seq_path = settings.SEQ_PATH + reference.checksum_sequence + '.sequence'
    try:
        seq = Seq(get_sequence_mmap(seq_path, 1, reference.length + 1),
                  generic_dna)
    except IOError:
        return None
    else:
        record.seq = seq

    return record


def get_description_boundary_positions(description):
    """
    Determines the minimum and maximum positions that appear in the operations
    of an HGVS description.

    For 'NC_000001.11:g.[100T>C;750A>G;2000C>T]' it should return 100 and 2000.

    :param description: Parsed HGVS description.
    :return: Minimum and maximum positions that appear in the description.
    """
    if description.SingleAlleleVarSet:
        variants = [v.RawVar for v in description.SingleAlleleVarSet]
    else:
        variants = [description.RawVar]

    positions = set()
    for variant in variants:
        first_location = last_location = variant.StartLoc.PtLoc
        if variant.EndLoc:
            last_location = variant.EndLoc.PtLoc
        positions.add(int(first_location.Main))
        positions.add(int(last_location.Main))

    return min(positions), max(positions)


def configuration_check(output):
    # Check if the required variables were set into settings.
    try:
        settings.SEQ_PATH
        settings.DATABASE_GB_URI
    except (NameError, AttributeError):
        output.addMessage(
            __file__, 2, 'NCSETTINGS', 'Chromosomal database settings not set.')
        return False

    # Check if the sequence folder exists.
    if not os.path.isdir(settings.SEQ_PATH):
        output.addMessage(
            __file__, 2, 'NCSEQDIR', 'Sequence directory does not exist.')
        return False

    return True


def get_nc_record(record_id, description, parsed_description, output):
    """
    Get an NC record from the gbparser database and transform it into the
    Mutalyzer record format.

    :param record_id: HGVS format record description.
    :param parsed_description: HGVS description object.
    :return: The record in mutalyzer format or None if not found.
    """
    if not configuration_check(output):
        return None

    if parsed_description.RefType in ['p', 'm', 'n']:
        output.addMessage(
            __file__, 4, 'ECHROMCOORD', 'Could not retrieve information for '
                                        'the provided {}. coordinate system.'
                                        .format(parsed_description.RefType))
        return None

    # Get the accession
    accession = record_id.split('.')[0]
    version = record_id.split('.')[1]

    reference = get_reference(accession, version)

    # If no reference present in the database we just return.
    if reference is None:
        return None

    if parsed_description.RefType == 'g':
        # Find the start and end positions.
        position_start, position_end = get_description_boundary_positions(parsed_description)
        # Get the DB transcript entries.
        db_transcripts = get_transcripts(reference, position_start, position_end)
    elif parsed_description.RefType == 'c':
        db_transcripts = Transcript.query.filter_by(reference_id=reference.id).all()
    else:
        return None

    print("%s: End nc_db get_transcripts" % datetime.now())

    return get_mutalyzer_record(reference, db_transcripts)


def get_db_boundaries_positions(reference, position_start, position_end):
    """
    Providing two positions on a chromosome reference this method will query
    the gbparser database and will return the extremity positions of the
    transcript entries found that contain the input positions.

    Example 1:
                      p_s                        p_e
                       |                          |
            |------------------------------------------------|
                                                |---------------------|
           t_1                                                          t_2

    Example 2:
                      p_s                        p_e
                       |                          |
            |------------------|
                                             |---------------------|
           t_1                                                    t_2

    Example 2:
                      p_s                        p_e
                       |                          |
                           |------------------|
                                          |---------------------|
                      t_1                                      t_2

    return: t_1, t_2

    :param reference: Database reference entry.
    :param position_start:
    :param position_end:
    :return:
    """
    p_s = position_start
    p_e = position_end

    positions = {p_s, p_e}

    transcripts = Transcript.query.filter_by(reference_id=reference.id). \
        filter(((Transcript.transcript_start <= p_s) & (Transcript.transcript_stop >= p_e)) |
               ((Transcript.transcript_start >= p_s) & (Transcript.transcript_stop <= p_e)) |
               ((Transcript.transcript_start <= p_s) & (Transcript.transcript_stop >= p_s)) |
               ((Transcript.transcript_start <= p_e) & (Transcript.transcript_stop >= p_e))).all()

    for transcript in transcripts:
        positions.add(transcript.transcript_start)
        positions.add(transcript.transcript_stop)

    p_s = min(positions)
    p_e = max(positions)

    return p_s, p_e


def get_transcripts(reference, position_start, position_end):
    """
    Retrieves the transcripts information from the database for the provided
    reference that are between the provided start and end positions to which
    5000 is subtracted and added, respectively.
    :param reference:
    :param position_start:
    :param position_end:
    :return:
    """
    if position_start > 5000:
        p_s = position_start - 5000
    else:
        p_s = 1

    if position_end < reference.length - 5000:
        p_e = position_end + 5000
    else:
        p_e = reference.length

    p_s, p_e = get_db_boundaries_positions(reference, p_s, p_e)

    transcripts = Transcript.query.filter_by(reference_id=reference.id). \
        filter(((Transcript.transcript_start <= p_s) & (Transcript.transcript_stop >= p_e)) |
               ((Transcript.transcript_start >= p_s) & (Transcript.transcript_stop <= p_e)) |
               ((Transcript.transcript_start <= p_s) & (Transcript.transcript_stop >= p_s)) |
               ((Transcript.transcript_start <= p_e) & (Transcript.transcript_stop >= p_e))).all()

    return transcripts
