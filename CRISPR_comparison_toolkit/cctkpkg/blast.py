#!/usr/bin/env python3


import sys
import argparse
import subprocess
import multiprocessing
import re
from collections import defaultdict

from . import sequence_operations, file_handling



class BlastResult():
    """
    A class to store column contents from a blast result in an easily 
    retrieved form.

    Attributes:
        qseqid (str): query (e.g., unknown gene) sequence id
        sseqid (str): subject (e.g., reference genome) sequence id
        pident (float): percentage of identical matches
        length (int): alignment length (sequence overlap)
        mismatch (int): number of mismatches
        gapopen (int): number of gap openings
        qstart (int): start of alignment in query
        qend (int): end of alignment in query
        sstart (int): start of alignment in subject
        send (int): end of alignment in subject
        evalue (str): expect value
        bitscore (float): bit score
        qlen (int): length of query sequence
        slen (int): length of subject sequence
        strand (str): Whether the blast hit was on the top (plus) or
          bottom (minus) strand of the DNA
    """
    def __init__(self, blast_line):
        bits = blast_line.split('\t')
        self.qseqid = bits[0]
        self.sseqid = bits[1]
        self.pident = float(bits[2])
        self.length = int(bits[3])
        self.mismatch = int(bits[4])
        self.gapopen = bits[5]
        self.qstart = int(bits[6])
        self.qend = int(bits[7])
        self.sstart = int(bits[8])
        self.send = int(bits[9])
        self.evalue = bits[10]
        self.bitscore = bits[11]
        self.qlen = int(bits[12])
        self.slen = int(bits[13])
        self.sseq = bits[14]
        self.strand = 'plus' if self.sstart < self.send else 'minus'


def run_blastcmd(db, fstring, batch_locations):
    """
    function to call blastdbcmd in a shell and process the output. Uses
    a batch query of the format provided in the blastdbcmd docs. e.g.
    printf "%s %s %s %s\\n%s %s %s\\n" 13626247 40-80 plus 30 14772189 \
    1-10 minus | blastdbcmd -db GPIPE/9606/current/all_contig \
    -entry_batch -
    
    Args:
        db (str): path to the blast db you want to query.
        fstring (str):  The string to give to printf
        batch_locations (str):  the seqid, locations, and strand of all
          the spacers to be retrieved
    
    Returns:
        (list) List of the sequences of regions returned by blastdbcmd
          in response to the query locations submitted
    
    Raises:
        ERROR running blastdbcmd: Raises an exception when blastdbcmd
          returns something to stderr, prints the error as well as the
          information provided to blastdbcmd and aborts the process.
    """
    x = subprocess.run(
        'printf "{}" {}| blastdbcmd -db {} -entry_batch -\
            '.format(fstring, batch_locations, db), 
        shell=True,
        universal_newlines=True,
        capture_output=True
        ) 
    if x.stderr:
        print("ERROR running blastdbcmd on {} :\n{}".format(
            db, batch_locations, x.stderr))
        sys.exit()
    else:
        return [i for i in x.stdout.split('\n') if '>' not in i and len(i) > 0]


def pool_MP_spacer_finder(array_entries, args, threads, chunksize):
    """
    Manages the multiprocess worker pool and returns to arrays found by
    build_arrays_MP to whatever called it.
    Args:
        array_entries (list): A list of lists of blast_entry class
          instances.
        args (argparse class): All of the argparse options given by the
          user.
        threads (int): Number of threads to use for multiprocessing.
        Chunksize (int):    Approximate size of batches into which pool
          should split array_entries 
    
    Returns:
        (list) List of FoundArray instances representing distinct arrays
    """

    pool = multiprocessing.Pool(processes=threads)

    options = [(array_entry, args) for array_entry in array_entries]
    output = pool.starmap(build_arrays_MP, options, chunksize)
    output = [i for i in output]
    pool.close()
    pool.join()

    return output


def build_arrays_MP(array_entry, args):
    """
    Builds the fstring and batch_locations objects to give to
    run_blastdbcmd and then calls that and returns the output.
    Args:
        array_entry (list): list of blast_entry class instances.
        args (argparse class): All of the argparse options given by the
          user.
    
    Returns:
        (FoundArray) FoundArray instance representing the contiguous
          array identified.
    
    Raises:
        Exception: Raises an exception and prints the exception if one
          occurs.
    """
    try:
        # retrieve spacer seqs in batch for quicker blastdbcmd usage.
        sp_fstring = '' 
        sp_batch_locations = ''
        # same for repeats
        rep_fstring = '' 
        rep_batch_locations = ''
       
        array = file_handling.FoundArray()
        if args.regex_pattern:
            array.genome = re.match(
                args.regex_pattern, array_entry[0].sseqid)[0]
        elif args.assemblies:
            array.genome = which_substring_in_string(
                args.assemblies,
                array_entry[0].sseqid)
            print(array.genome)
        else:
            array.genome = sseqid
        n_reps = len(array_entry)
        array.repeats = [a.sseq for a in array_entry] 

        ### Retrieve repeat sequence using blastdbcmd and fstring

        array.repeat_id = array_entry[0].qseqid
        array.contig = array_entry[0].sseqid
        if array_entry[0].strand == 'plus':
            array.start = array_entry[-1].send
            array.stop = array_entry[0].sstart
            for i, entry in enumerate(array_entry):
                if i+1 != n_reps:
                    sp_fstring += '%s %s %s\n'
                    sp_batch_locations += '{} {}-{} {} '.format(
                        entry.sseqid,
                        entry.send+1,
                        array_entry[i+1].sstart-1,
                        entry.strand
                        )
                rep_fstring += '%s %s %s\n'
                rep_batch_locations += '{} {}-{} {} '.format(
                    entry.sseqid,
                    entry.sstart,
                    entry.send,
                    entry.strand
                    )
        else:
            array.start = array_entry[0].send
            array.stop = array_entry[-1].sstart
            for i, entry in enumerate(reversed(array_entry)):
                if i+1 != n_reps:
                    sp_fstring += '%s %s %s\n'
                    sp_batch_locations += '{} {}-{} {} '.format(
                        entry.sseqid,
                        array_entry[-(i+2)].send+1,
                        entry.sstart-1,
                        entry.strand
                        )
                rep_fstring += '%s %s %s\n'
                rep_batch_locations += '{} {}-{} {} '.format(
                    entry.sseqid,
                    entry.sstart,
                    entry.send,
                    entry.strand
                    )
        array.spacers = run_blastcmd(
            args.blast_db_path,
            sp_fstring,
            sp_batch_locations
            )
        array.repeats = run_blastcmd(
            args.blast_db_path,
            rep_fstring,
            rep_batch_locations
            )

        return array


    except Exception as e:
        print(e)
        sys.exit()


def blastn_to_arrays(args):
    """
    Runs blastn of provided repeats against provided blastdb. Processes
    the blast output to find hits that are good enough to keep.
    Passes those good hits to identify_same_array_hits and stores the
    list of arrays returned to it. Passes that list of arrays to
    pool_MP_spacer_finder and stores the arrays returned to it that now
    have spacer information. If arrays were found, these are returned.
    If none are found an error is raised.
    Args:
        args (argparse class): All of the argparse options given by the user.

    
    Returns:
        (list): List of lists of FoundArray instances representing
          distinct arrays.
    
    Raises:
        No arrays found: If no arrays were found in the blast database
          then this raises an error stating that and exits.
    """ 
    blastn_command = "blastn -query {} -db {} -task blastn-short -outfmt \
    '6 std qlen slen sseq' -num_threads {} -max_target_seqs {} -evalue {} {}\
    ".format(
        args.repeats_file,
        args.blast_db_path,
        args.num_threads,
        args.max_target_seqs,
        args.evalue,
        args.other_blast_options
        )
    blast_run = subprocess.run(
        blastn_command,
        shell=True,
        universal_newlines=True,
        capture_output=True
        )

    if blast_run.stderr:
        print("ERROR running blast on {}:\n{}".format(
            args.blast_db_path, blast_run.stderr))
        sys.exit()
    blast_lines = [
        BlastResult(i) for i in blast_run.stdout.split('\n') if len(i) > 0]
    

    good_hits = []

    for line in blast_lines:
        # Discard any blast hits where match is < 90% of query
        if line.length < 0.9 * line.qlen:
            continue
        # Or <80% ID
        if line.pident < 80:
            continue
        if line.length != line.qlen:
            if line.sstart < line.send:
                # if the blast hit is on the positive strand
                if line.qstart != 1:
                    # if the 5' end of repeat is missing
                    # adjust sstart by the number of bases missing from
                    # repeat 5' end
                    line.sstart = line.sstart - (line.qstart - 1)
                if line.qlen != line.qend:
                    # if the 3' end of the repeat is missing
                    # adjust send by the number of bases missing from
                    # repeat 3' end
                    line.send = line.send + (line.qlen - line.length)
            elif line.sstart > line.send:
                # if the blast hit is on the negative strand
                if line.qstart != 1:
                    # if the 5' end of repeat is missing
                    # adjust sstart by the number of bases missing from
                    # repeat 5' end
                    line.sstart = line.sstart + (line.qstart - 1)
                if line.qlen != line.qend:
                    # if the 3' end of the repeat is missing
                    # adjust send by the number of bases missing from
                    # repeat 3' end
                    line.send = line.send - (line.qlen - line.length)
        
        # # Adjust start and end for reverse orientation matches
        if line.sstart > line.send:
            line.sstart, line.send = line.send, line.sstart
        good_hits.append(line)
    # check_blast_overlap(good_hits)
    good_hits_sorted = sorted(
        good_hits, key = lambda x: (x.sseqid, x.strand, x.sstart))
    array_entries = identify_same_array_hits(good_hits_sorted, args)
    

    if len(array_entries) > 0:
        chunksize = max(int(len(array_entries)//args.num_threads), 1)
        arrays = pool_MP_spacer_finder(
            array_entries,
            args,
            args.num_threads,
            chunksize)

        return arrays

    else:
        print("No arrays identified from the query file {} in the database {}\
            ".format(args.repeats_file, args.blast_db_path))
        sys.exit()


def identify_same_array_hits(blast_entries, args):
    """
    Given a list of blast results, identify which of them are from a contiguous spacer array by checking:
    Are they close enough?
    Are they on the same strand?
    Are they in the same genome, contig?
    If it finds 3 or more repeats close together on the same strand of the same contig, groups them into an array.
    Args:
        blast_entries (list):   List of blast entry classes sorted by contig, then start, then location.
        args (argparse class):  The arguments provided to argparse.
    
    Returns:
        (list) List of lists of blast_entry classes.
    """
    arrays = []
    n_spacers = 0 # count spacers in the current array
    last_contig = False
    last_loc = False # store start coord of last repeat
    last_strand = False # store orientation of last repeat
    for entry in blast_entries:
        if last_loc:
            if last_loc > entry.sstart - args.rep_interval and last_strand == entry.strand and last_contig == entry.sseqid:
                this_array.append(entry)
                last_loc = entry.sstart
                if entry == blast_entries[-1]:
                    if len(this_array) > 2:
                        arrays.append(this_array)
                        last_loc = entry.sstart
                        last_strand = entry.strand
                        last_contig = entry.sseqid
                        this_array = [entry]

                    else:
                        last_loc = entry.sstart
                        last_strand = entry.strand
                        last_contig = entry.sseqid
                        this_array = [entry]
                
            else:
                if len(this_array) > 2:
                    arrays.append(this_array)
                    last_loc = entry.sstart
                    last_strand = entry.strand
                    last_contig = entry.sseqid
                    this_array = [entry]

                else:
                    last_loc = entry.sstart
                    last_strand = entry.strand
                    last_contig = entry.sseqid
                    this_array = [entry]
        else:
            last_loc = entry.sstart
            last_strand = entry.strand
            last_contig = entry.sseqid
            this_array = [entry]

    return arrays


def which_substring_in_string(substrings, string):
    present = []
    for s in substrings:
        if s in string:
            present.append(s)

    # If more than one present, longest will be best
    if len(present) > 1:
        return max(present, key=len)
    else:
        return present[0]


def check_repeat_similarity(repeats_file):
    blastn_command = "blastn -query {} -subject {} -task blastn-short -outfmt \
        '6 std qlen slen sseq'\
        ".format(
            repeats_file,
            repeats_file,
            )
    blast_run = subprocess.run(
        blastn_command,
        shell=True,
        universal_newlines=True,
        capture_output=True
        )

    blast_lines = [
    BlastResult(i) for i in blast_run.stdout.split('\n') if len(i) > 0]
    similar_pairs = []
    for line in blast_lines:
        if (
            (line.length > 0.9*line.qlen or line.length > 0.9*line.slen)
            and (line.pident > 90)
            and (line.qseqid != line.sseqid)):
            
            if not (((line.qseqid, line.sseqid) in similar_pairs)
                or ((line.sseqid, line.qseqid) in similar_pairs)):
                similar_pairs.append((line.qseqid, line.sseqid))
    if len(similar_pairs) > 0:
        sys.stderr.write(
            "The following pairs of repeats in the provided file are "
            "over 90% identical across over 90% of their length. This "
            "program can only identify arrays associated with distinct "
            "repeats.\nSimilar pairs of arrays:\n{}\n".format(
                "\n".join(["\t".join(i) for i in similar_pairs])))
        sys.exit()


def build_parser(parser):

    req_options = parser.add_argument_group("Required arguments")
    req_options.add_argument(
        "-r", "--repeats",
        metavar=" ",
        dest="repeats_file",
        required=True,
        help="fasta format file of your CRISPR repeats"
        )
    req_options.add_argument(
        "-d", "--blastdb",
        metavar=" ",
        dest="blast_db_path",
        required=True,
        help="path to blast db files (not including the file extensions)"
        )
    req_options.add_argument(
        "-o", "--outdir",
        metavar=" ",
        required=True,
        help="path to directory you want output files saved in"
        )

    other_options = parser.add_argument_group("Other inputs")
    other_options.add_argument(
        "-a", "--assemblies",
        metavar=" ",
        required=False,
        help="file containing a list of your assembly names"
        )
    other_options.add_argument(
        "-p", "--regex-pattern",
        metavar=" ",
        dest="regex_pattern",
        required=False,
        help="pattern describing your assembly names"
        )
    other_options.add_argument(
        "-q", "--regex-type",
        metavar=" ",
        required=False,
        default='P',
        help="DEFAULT: P. {E, P} Regex type to be used in step taht calls grep \
        See CCTK documentation for details."
        )
    other_options.add_argument(
        "-t", "--threads",
        metavar=" ",
        dest="num_threads",
        required=False,
        default=1,
        type=int,
        help="DEFAULT: 1. Number of threads to use"
        )
    other_options.add_argument(
        "-i", "--repeat-interval",
        metavar=" ",
        dest="rep_interval",
        required=False,
        default=80,
        type=int,
        help="DEFAULT: 80. Set the expected interval between the start \
        position of your repeats."
        )

    blast_options = parser.add_argument_group("BLASTn settings")
    blast_options.add_argument(
        "-e", "--evalue",
        metavar=" ",
        required=False,
        default='10',
        help="DEFAULT: 10. set the evalue for blastn"
        )
    blast_options.add_argument(
        "-m", "--max-target-seqs",
        metavar=" ",
        required=False,
        default='10000',
        help="DEFAULT: 10000. Set the max_target_seqs option for blastn \
        (see CCTK documentation for details)"
        )
    blast_options.add_argument(
        "-x", "--blast-options",
        metavar=" ",
        dest="other_blast_options",
        required=False,
        default='',
        help="DEFAULT: none. Input additional blastn options. \
        Forbidden options: blastn -query -db -task -outfmt -num_threads \
        -max_target_seqs -evalue"
        )
    
    return parser


def main(args):

    outdir = args.outdir + '/' if args.outdir[-1] != '/' else args.outdir

    # Make sure repeats aren't similar enough to return overlapping hits
    check_repeat_similarity(args.repeats_file)
        
    if args.regex_pattern:
        assembly_names = [a for a in subprocess.run(
            "blastdbcmd -db {} -entry all -outfmt '\%a' | grep -o{} '{}' | \
            sort | uniq".format(
                args.blast_db_path,
                args.regex_type,
                args.regex_pattern
                ),
            shell=True,
            universal_newlines=True,
            capture_output=True
            ).stdout.split('\n') if len(a) > 0]

    elif args.assemblies:
        # store assembly names in args.assemblies for easy delivery
        assembly_names = file_handling.read_assembly_list_file(args.assemblies)
        args.assemblies = assembly_names

    else:
        sys.stderr.write(
            "WARNING!!!\n\n"
            "In order to associate arrays found in different contigs from the "
            "same assembly, you must provide either a list of assembly names "
            "or a regex to identify your assemblies. Each sequence in your "
            "blastdb will be treated as a separate assembly.\n")
        assembly_names = [a for a in subprocess.run(
            "blastdbcmd -db {} -entry all -outfmt '\%a'".format(
                args.blast_db_path
                ),
            shell=True,
            universal_newlines=True,
            capture_output=True
            ).stdout.split('\n') if len(a) > 0]
        args.regex_pattern = ".+"

    # Find the names of all the genomes being searched by looking for 
    # the user-provided regex in the blast-db sequence IDs.
    # for all genomes, start their entry in the AssemblyCRISPRs with
    # the placeholder False. This indicated no CRISPRs found.
    # If a CRISPR array is found later this entry will be overwritten
    # with infor about the arrays.
    all_assemblies_dict = {
        k: file_handling.AssemblyCRISPRs() for k in assembly_names}
    for k,v in all_assemblies_dict.items():
        v.accession = k

    all_arrays = blastn_to_arrays(args)

    for array in all_arrays:
        assembly = all_assemblies_dict[array.genome]
        assembly.array_count+=1
        assembly.has_crispr = True
        assembly.arrays[assembly.array_count] = array

    # unpack dict into list for processing
    all_assemblies = [i for i in all_assemblies_dict.values()]

    (
        non_red_spacer_dict,
        non_red_spacer_id_dict,
        non_red_array_dict,
        non_red_array_id_dict
        ) = sequence_operations.non_redundant_CR(all_assemblies)

    # Output summary info
    sys.stderr.write("Total unique spacers: %i\n" %len(non_red_spacer_id_dict))
    sys.stderr.write("Total unique arrays: %i\n" %len(non_red_array_id_dict))

    # Fill in spacer and array ID info
    sequence_operations.add_ids(
        all_assemblies,
        non_red_spacer_id_dict,
        non_red_array_id_dict)

    # Save array files
    file_handling.write_CRISPR_files(
        all_assemblies,
        non_red_spacer_id_dict,
        non_red_array_id_dict,
        args.outdir)


if __name__ == '__main__':

    parser = argparse.ArgumentParser(
        description="Given fasta CRISPR repeats and a blast db of genomes, \
        pulls out spacer arrays of >= 2 spacers")
    parser = build_parser(parser)

    if len(sys.argv) == 1:
        parser.parse_args(['--help'])
    else:
        args = parser.parse_args()

        main(args) 
