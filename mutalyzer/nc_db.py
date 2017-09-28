import time
from datetime import datetime

import mmap
from Bio.Seq import Seq
from Bio.Alphabet import generic_dna, generic_protein
from mutalyzer.GenRecord import PList, Locus, Gene, Record
from mutalyzer.dbgb.models import Transcript, Reference


def get_sequence_mmap(file_path, start, end):
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
    :return: Database entry.
    """
    reference = Reference.query.filter_by(accession=str(accession), version=str(version)) \
        .order_by(Reference.id.asc()) \
        .first()
    return reference


def complete_muta_record(record, db_transcripts):
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


def get_nc_record(record_id, description, output):
    """
    Get an NC record from the gbparser database and transform it into the
    Mutalyzer record format.

    :param record_id: HGVS format record description.
    :param description: HGVS description object.
    :return: The record in mutalyzer format or None if not found.
    """
    if description.RefType in ['p', 'm', 'n']:
        output.addMessage(
            __file__, 4, 'ECHROMCOORD', "Could not retrieve information for "
                                        "the provided '{}.' coordinate system.".format(description.RefType))
        return None

    record = Record()
    accession = record_id.split('.')[0]
    version = record_id.split('.')[1]
    reference = get_reference(accession, version)
    if reference is None:
        return None
    record.source_id = reference.accession
    record.id = reference.accession
    record.source_accession = reference.accession
    record.source_version = reference.version
    record.organism = 'Homo sapiens'
    record.molType = 'g'

    if description.RefType == 'g':
        variant = description.RawVar
        print(variant)
        position_start = position_end = variant.StartLoc.PtLoc
        if variant.EndLoc:
            position_end = variant.EndLoc.PtLoc

        print('first location {}'.format(position_start[0]))
        print('last location {}'.format(position_end[0]))
        # Get the DB transcript entries.
        db_transcripts = get_transcripts(accession, version, position_start[0], position_end[0])
    elif description.RefType == 'c':
        db_transcripts = Transcript.query.filter_by(reference_id=reference.id).all()
    else:
        return None

    complete_muta_record(record, db_transcripts)

    # Get the sequence.
    seq_path = "/home/mlefter/projects/mutalyzer/data/gbparserout/" \
               + reference.checksum_sequence + '.sequence'
    record.seq = Seq(get_sequence_mmap(seq_path, 1, reference.length + 1), generic_dna)

    return record


def get_transcripts(accession, version, position_start, position_end):
    reference = get_reference(accession, version)

    p_s = int(position_start) - 100000
    p_e = int(position_end) + 100000

    transcripts = Transcript.query.filter_by(reference_id=reference.id). \
        filter(((Transcript.transcript_start <= p_s) & (p_s <= Transcript.transcript_stop)) |
               ((Transcript.transcript_start <= p_e) & (p_e <= Transcript.transcript_stop)) |
               ((Transcript.transcript_start >= p_s) & (p_e >= Transcript.transcript_stop))
               ) \
        .all()
    return transcripts
