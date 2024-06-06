#!/usr/bin/env python3

import argparse

# Default values for functions
POLYA_DEFAULT = 16
MISMATCH_DEFAULT = 0.1
ADAPTOR_LENGTH_DEFAULT = 0
KEEP_DEFAULT = True

stats_dict = {
    "matched": 0,
    "not_matched": 0,
    "5_prime_polyA": 0,
    "3_prime_polyA": 0,
    "5_prime_polyT": 0,
    "3_prime_polyT": 0,
}


def check_read(
    lines,
    polya_length=POLYA_DEFAULT,
    mismatch=MISMATCH_DEFAULT,
    adaptor_length=ADAPTOR_LENGTH_DEFAULT,
    keep=KEEP_DEFAULT,
):
    """Check whether fastq read has matching polyA/T

    :param lines: 4-line fastq read
    :type lines: str
    :param polya_length: length of polyA/T to match, defaults to POLYA_DEFAULT
    :type polya_length: int, optional
    :param mismatch: mismatch rate in polyA/T, defaults to MISMATCH_DEFAULT
    :type mismatch: float, optional
    :param adaptor_length: length of adapter to pad by, defaults to ADAPTOR_LENGTH_DEFAULT
    :type adaptor_length: int, optional
    :param keep: whether to keep matching read, defaults to KEEP_DEFAULT
    :type keep: bool, optional
    :return: list of lines from the matching read
    :rtype: list
    """
    seq = lines[1]
    match = False
    num_mismatches = int(polya_length * mismatch)
    for i in range(0, adaptor_length + polya_length):
        start_window = seq[i : i + polya_length]
        end_window = seq[-(i + polya_length) : -i]
        if start_window.count("A") >= (polya_length - num_mismatches):
            match = True
            stats_dict["5_prime_polyA"] += 1
            break
        elif start_window.count("T") >= (polya_length - num_mismatches):
            match = True
            stats_dict["5_prime_polyT"] += 1
            break
        elif end_window.count("A") >= (polya_length - num_mismatches):
            match = True
            stats_dict["3_prime_polyA"] += 1
            break
        elif end_window.count("T") >= (polya_length - num_mismatches):
            match = True
            stats_dict["3_prime_polyT"] += 1
            break
    if match:
        stats_dict["matched"] += 1
    else:
        stats_dict["not_matched"] += 1
    if keep and match:
        return lines
    elif not keep and not match:
        return lines
    else:
        return []


def filter_reads(
    input_file,
    output_file,
    polya_length=POLYA_DEFAULT,
    mismatch=MISMATCH_DEFAULT,
    adaptor_length=ADAPTOR_LENGTH_DEFAULT,
    keep=KEEP_DEFAULT,
):
    """Filter an input fastq file and output to new fastq file

    :param input_file: the input fastq filename
    :type input_file: str
    :param output_file: the output fastq filename, if False write to stdout
    :type output_file: str
    :param polya_length: length of polyA/T to match, defaults to POLYA_DEFAULT
    :type polya_length: int, optional
    :param mismatch: mismatch rate in polyA/T, defaults to MISMATCH_DEFAULT
    :type mismatch: float, optional
    :param adaptor_length: length of adapter to pad by, defaults to ADAPTOR_LENGTH_DEFAULT
    :type adaptor_length: int, optional
    :param keep: whether to keep matching read, defaults to KEEP_DEFAULT
    :type keep: bool, optional
    """
    with open(input_file, "r") as file:
        if output_file:
            out = open(output_file, "w")
        while True:
            lines = [file.readline() for _ in range(4)]
            if not lines[0]:
                break
            processed = check_read(
                lines,
                polya_length=polya_length,
                mismatch=mismatch,
                adaptor_length=adaptor_length,
                keep=keep,
            )
            if processed and output_file:
                out.write("".join(processed))
            elif processed and not output_file:
                print("".join(processed))


def format_stats(stats_dict, input_file):
    matched = stats_dict["matched"]
    not_matched = stats_dict["not_matched"]
    total = matched + not_matched
    pc_matched = (matched / total) * 100
    five_prime_polyA = stats_dict["5_prime_polyA"]
    five_prime_polyT = stats_dict["5_prime_polyT"]
    three_prime_polyA = stats_dict["3_prime_polyA"]
    three_prime_polyT = stats_dict["3_prime_polyT"]
    return (
        f"# polya-filter stats on {input_file}\n"
        f"\n"
        f"reads with polyA/T match: {matched} ({pc_matched:.2f}% of total)\n"
        f"reads without polyA/T match: {not_matched} ({100-pc_matched:.2f}% of total)\n"
        f"total reads: {total}\n"
        f"\n"
        f"reads with 5' polyA match: {five_prime_polyA} ({(five_prime_polyA/matched)*100:.2f}% of matched)\n"
        f"reads with 5' polyT match: {five_prime_polyT} ({(five_prime_polyT/matched)*100:.2f}% of matched)\n"
        f"reads with 3' polyA match: {three_prime_polyA} ({(three_prime_polyA/matched)*100:.2f}% of matched)\n"
        f"reads with 3' polyT match: {three_prime_polyT} ({(three_prime_polyT/matched)*100:.2f}% of matched)\n"
    )


if __name__ == "__main__":
    # Set up command line arguments:
    parser = argparse.ArgumentParser()

    parser.add_argument(
        "-i", "--input", type=str, help="input fastq filename", required=True
    )
    parser.add_argument(
        "-o",
        "--output",
        type=str,
        help="output fastq filename, omit to write to stdout",
        default=False,
    )
    parser.add_argument(
        "-p",
        "--polya_length",
        type=int,
        help="minimum lenth of polyA/T to filter by",
        default=POLYA_DEFAULT,
    )
    parser.add_argument(
        "-m",
        "--mismatch",
        type=float,
        help="proportion of base mismatches within polyA/T (0.0-1.0)",
        default=MISMATCH_DEFAULT,
    )
    parser.add_argument(
        "-a",
        "--adaptor_length",
        type=int,
        help="length of adaptor to pad with before polyA/T",
        default=ADAPTOR_LENGTH_DEFAULT,
    )
    parser.add_argument(
        "--keep", action="store_true", help="Keep reads with polyA/T match (default)"
    )
    parser.add_argument(
        "--discard",
        dest="keep",
        action="store_false",
        help="Discard reads with polyA/T match",
    )
    parser.set_defaults(keep=KEEP_DEFAULT)
    parser.add_argument(
        "--stats",
        dest="stats",
        action="store_true",
        help="Write stats to text file",
    )
    parser.set_defaults(stats=False)

    args = parser.parse_args()

    input_file = args.input
    output_file = args.output
    polya_length = args.polya_length
    mismatch = args.mismatch
    adaptor_length = args.adaptor_length
    keep = args.keep
    stats = args.stats

    # Run function on file
    filter_reads(
        input_file,
        output_file,
        polya_length=polya_length,
        mismatch=mismatch,
        adaptor_length=adaptor_length,
        keep=keep,
    )

    # Write stats
    if stats:
        if output_file:
            stats_out_file = output_file.replace(".fastq", ".stats")
        else:
            stats_out_file = input_file.replace(".fastq", ".stats")
        formatted_stats = format_stats(stats_dict, input_file)
        with open(stats_out_file, "w") as out:
            out.write(formatted_stats)
