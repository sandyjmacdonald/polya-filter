#!/usr/bin/env python

import argparse

POLYA_DEFAULT = 16
MISMATCH_DEFAULT = 0.1
ADAPTOR_LENGTH_DEFAULT = 30
KEEP_DEFAULT = True


def check_read(
    lines,
    polya_length=POLYA_DEFAULT,
    mismatch=MISMATCH_DEFAULT,
    adaptor_length=ADAPTOR_LENGTH_DEFAULT,
    keep=KEEP_DEFAULT,
):
    seq = lines[1]
    match = False
    num_mismatches = int(polya_length * mismatch)
    for i in range(0, adaptor_length + polya_length):
        start_window = seq[i : i + polya_length]
        end_window = seq[-(i + polya_length) : -i]
        if start_window.count("A") >= (polya_length - num_mismatches):
            match = True
            break
        elif start_window.count("T") >= (polya_length - num_mismatches):
            match = True
            break
        elif end_window.count("A") >= (polya_length - num_mismatches):
            match = True
            break
        elif end_window.count("T") >= (polya_length - num_mismatches):
            match = True
            break
    if keep and match:
        return lines
    elif not keep and not match:
        return lines
    else:
        return None


def filter_reads(
    input_file,
    output_file,
    polya_length=POLYA_DEFAULT,
    mismatch=MISMATCH_DEFAULT,
    adaptor_length=ADAPTOR_LENGTH_DEFAULT,
    keep=KEEP_DEFAULT,
):
    with open(input_file, "r") as file:
        with open(output_file, "w") as out:
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
                if processed is not None:
                    out.write("".join(processed))


if __name__ == "__main__":
    # Set up command line arguments:
    parser = argparse.ArgumentParser()

    parser.add_argument(
        "-i", "--input", type=str, help="input fastq filename", required=True
    )
    parser.add_argument(
        "-o", "--output", type=str, help="output fastq filename", required=True
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
        "--keep", action="store_true", help="Keep reads with polyA/T match"
    )
    parser.add_argument(
        "--discard",
        dest="keep",
        action="store_false",
        help="Discard reads with polyA/T match",
    )
    parser.set_defaults(keep=KEEP_DEFAULT)

    args = parser.parse_args()

    input_file = args.input
    output_file = args.output
    polya_length = args.polya_length
    mismatch = args.mismatch
    adaptor_length = args.adaptor_length
    keep = args.keep

    filter_reads(
        input_file,
        output_file,
        polya_length=polya_length,
        mismatch=mismatch,
        adaptor_length=adaptor_length,
        keep=keep,
    )
