#!/bin/env python

import argparse
import hashlib
import logging
import os
from datetime import datetime
from pathlib import Path

import pandas as pd


def get_checksum(file_name):
    """
    Get MD5 checksum of a file.
    """
    hash_md5 = hashlib.md5()
    with open(file_name, "rb") as f:
        for chunk in iter(lambda: f.read(4096), b""):
            hash_md5.update(chunk)
    return hash_md5.hexdigest()


def die(message):
    raise RuntimeError(message)


def get_old_sample_from_input(input_path, archive_base, cohort):
    """
    Infer old sample name from input path.

    Expected structure:
      /net/eichler/vol28/projects/long_read_archive/nobackups/{cohort}/{old_sample}/...
    """
    input_path = Path(input_path).resolve()
    cohort_base = archive_base / cohort

    try:
        rel = input_path.relative_to(cohort_base)
    except ValueError:
        die(
            f"Input path is not under expected cohort base:\n"
            f"  input_path: {input_path}\n"
            f"  cohort_base: {cohort_base}"
        )

    if len(rel.parts) < 1:
        die(f"Could not infer old sample from input path: {input_path}")

    return rel.parts[0]


def make_relative_under_old_sample(dest_path, old_base):
    """
    Convert absolute DEST_PATH to path relative to old sample directory.

    This prevents creating paths like:
      NEW_SAMPLE/net/eichler/vol28/...
    """
    dest_path = Path(dest_path).resolve()

    try:
        return dest_path.relative_to(old_base)
    except ValueError:
        die(
            f"DEST_PATH is not under old sample base:\n"
            f"  DEST_PATH: {dest_path}\n"
            f"  old_base:  {old_base}"
        )


def validate_no_nested_net_path(path_series, column_name):
    """
    Guard against accidental nested /net/eichler path creation.
    """
    bad = path_series.astype(str).str.contains(r"(^|/)net/eichler/vol28(/|$)", regex=True)

    if bad.any():
        examples = "\n".join(path_series.loc[bad].astype(str).head(10))
        die(
            f"{column_name} contains nested net/eichler paths. "
            f"This likely indicates path construction bug:\n{examples}"
        )


def safe_hardlink(src, dst, overwrite=False):
    """
    Create hard link safely.
    """
    src = Path(src)
    dst = Path(dst)

    if not src.exists():
        return False

    if dst.exists():
        if overwrite:
            dst.unlink()
        else:
            logging.warning(f"Skipping existing destination: {dst}")
            return False

    dst.parent.mkdir(parents=True, exist_ok=True)
    os.link(src, dst)
    return True


def main():
    logging.basicConfig(
        format="%(levelname)s (%(asctime)s): %(message)s "
               "(Line: %(lineno)d [%(filename)s])",
        datefmt="%m/%d/%Y %I:%M:%S %p",
        level=logging.INFO,
    )

    parser = argparse.ArgumentParser(
        formatter_class=argparse.ArgumentDefaultsHelpFormatter
    )

    parser.add_argument(
        "--input",
        "-i",
        type=str,
        required=True,
        help="Input copy record",
    )
    parser.add_argument(
        "--sample_new",
        "-s",
        type=str,
        required=True,
        help="Correct sample",
    )
    parser.add_argument(
        "--cohort",
        "-c",
        type=str,
        required=True,
        help="Cohort of correct sample",
    )
    parser.add_argument(
        "--output",
        "-o",
        type=str,
        required=False,
        help="Output file with files to remove",
    )
    parser.add_argument(
        "--archive_base",
        type=str,
        default="/net/eichler/vol28/projects/long_read_archive/nobackups",
        help="Base directory of long_read_archive nobackups",
    )
    parser.add_argument(
        "--dry-run",
        action="store_true",
        help="Print planned links without creating files",
    )
    parser.add_argument(
        "--overwrite",
        action="store_true",
        help="Overwrite existing destination hardlinks/files",
    )

    args = parser.parse_args()

    input_path = Path(args.input).resolve()
    archive_base = Path(args.archive_base).resolve()

    if not input_path.exists():
        die(f"Input file does not exist: {input_path}")

    old_sample = get_old_sample_from_input(input_path, archive_base, args.cohort)

    old_base = archive_base / args.cohort / old_sample

    if args.cohort:
        new_prefix = Path("..") / args.cohort / args.sample_new
    else:
        new_prefix = Path(args.sample_new)

    if not args.output:
        output_file = "-".join(
            [args.sample_new, "REMOVE", input_path.name]
        )
    else:
        output_file = args.output

    output_file = output_file.replace(".gz", "")

    logging.info(f"Input copy record: {input_path}")
    logging.info(f"Old sample: {old_sample}")
    logging.info(f"New sample: {args.sample_new}")
    logging.info(f"Old base: {old_base}")
    logging.info(f"New prefix: {new_prefix}")

    df = pd.read_csv(input_path, sep="\t")

    required_cols = {
        "DEST_PATH",
        "SEQ_TYPE",
        "RUN_ID",
        "SOURCE_PATH",
        "SIZE",
        "MOD_TIME",
        "STATUS",
        "MD5",
    }

    missing_cols = required_cols - set(df.columns)
    if missing_cols:
        die(f"Input table is missing required columns: {sorted(missing_cols)}")

    new_files = []

    logging.info("Looking for fast5 and pod5, this might take a while")

    for _, row in df.iterrows():
        dest_path = Path(row["DEST_PATH"])

        if not str(dest_path).endswith("fast5"):
            continue

        search_base = dest_path.parent.parent
        search_folder = dest_path.parent.name
        search_folder_pod = search_folder.replace("fast5", "pod5")

        search_file = dest_path.name
        search_file_pod = search_file.replace("fast5", "pod5")

        for folder_check in [search_folder, search_folder_pod]:
            for file_check in [search_file, search_file_pod]:
                new_file = search_base / folder_check / file_check

                if not new_file.is_file():
                    continue

                new_file_str = str(new_file)

                if new_file_str in df["DEST_PATH"].values:
                    continue

                logging.debug(f"Getting info for NEW_FILE: {new_file}")

                file_stat = new_file.stat()
                md5 = get_checksum(new_file)

                new_files.append(
                    {
                        "SAMPLE": "NA",
                        "SEQ_TYPE": row["SEQ_TYPE"],
                        "RUN_ID": row["RUN_ID"],
                        "SOURCE_PATH": "NA",
                        "DEST_PATH": new_file_str,
                        "SIZE": file_stat.st_size,
                        "MOD_TIME": datetime.fromtimestamp(
                            file_stat.st_mtime
                        ).strftime("%Y-%m-%d %H:%M:%S.%f"),
                        "STATUS": "Copied",
                        "MD5": md5,
                    }
                )

    if new_files:
        new_files_df = pd.DataFrame(new_files)
        df = pd.concat([df, new_files_df], ignore_index=True)

    df = df.reset_index(drop=True).fillna("NA")

    fastq_files = []

    # Preserve original logic:
    # old code used "/".join(args.input.split("/")[0:4]) / fastq / basename
    # Safer version: use archive root's first four components equivalent only if needed.
    #
    # For absolute input:
    #   /net/eichler/vol28/...
    # old split[0:4] produced:
    #   /net/eichler/vol28
    fastq_dir = Path("/net/eichler/vol28") / "fastq" / input_path.name.replace(".tab.gz", "")

    if fastq_dir.exists():
        for path in fastq_dir.rglob("*"):
            if path.is_file():
                fastq_files.append({"DEST_PATH": str(path.resolve())})

    fastq_files_df = pd.DataFrame(fastq_files)

    def make_copy_path_from_dest(dest):
        rel = make_relative_under_old_sample(dest, old_base)
        return str(new_prefix / rel)

    def make_new_path_from_dest(dest):
        rel = make_relative_under_old_sample(dest, old_base)
        return str(Path(args.sample_new) / rel)

    df["COPY_PATH"] = df["DEST_PATH"].apply(make_copy_path_from_dest)
    df["NEW_PATH"] = df["DEST_PATH"].apply(make_new_path_from_dest)

    if len(fastq_files_df) > 0:
        def make_fastq_copy_path(dest):
            dest = Path(dest).resolve()
            parts = [
                args.sample_new if part == old_sample else part
                for part in dest.parts
            ]

            # Avoid absolute path becoming absolute destination.
            # Keep relative version under new_prefix.
            rel_parts = [
                args.sample_new if part == old_sample else part
                for part in dest.parts[1:]
            ]

            return str(new_prefix / Path(*rel_parts))

        def make_fastq_new_path(dest):
            dest = Path(dest).resolve()
            rel_parts = [
                args.sample_new if part == old_sample else part
                for part in dest.parts[1:]
            ]
            return str(Path(args.sample_new) / Path(*rel_parts))

        fastq_files_df["COPY_PATH"] = fastq_files_df["DEST_PATH"].apply(
            make_fastq_copy_path
        )
        fastq_files_df["NEW_PATH"] = fastq_files_df["DEST_PATH"].apply(
            make_fastq_new_path
        )

    validate_no_nested_net_path(df["COPY_PATH"], "COPY_PATH")
    validate_no_nested_net_path(df["NEW_PATH"], "NEW_PATH")

    if len(fastq_files_df) > 0:
        validate_no_nested_net_path(fastq_files_df["COPY_PATH"], "fastq COPY_PATH")
        validate_no_nested_net_path(fastq_files_df["NEW_PATH"], "fastq NEW_PATH")

    logging.info("Planned path examples:")
    logging.info(
        "\n" + df[["DEST_PATH", "COPY_PATH", "NEW_PATH"]]
        .head(10)
        .to_string(index=False)
    )

    migrated_files = [str(input_path)]

    logging.info("Linking files")

    for _, row in df.iterrows():
        src = Path(row["DEST_PATH"])
        dst = Path(row["COPY_PATH"])

        if not src.exists():
            logging.warning(f"Source does not exist, skipping: {src}")
            continue

        if args.dry_run:
            logging.info(f"[dry-run] link {src} -> {dst}")
        else:
            linked = safe_hardlink(src, dst, overwrite=args.overwrite)
            if linked:
                migrated_files.append(str(src))

    if len(fastq_files_df) > 0:
        for _, row in fastq_files_df.iterrows():
            src = Path(row["DEST_PATH"])
            dst = Path(row["COPY_PATH"])

            if not src.exists():
                logging.warning(f"Source does not exist, skipping: {src}")
                continue

            if args.dry_run:
                logging.info(f"[dry-run] link {src} -> {dst}")
            else:
                linked = safe_hardlink(src, dst, overwrite=args.overwrite)
                if linked:
                    migrated_files.append(str(src))

    df["SAMPLE"] = args.sample_new
    df["DEST_PATH"] = df["NEW_PATH"]

    out_cols = [
        "SAMPLE",
        "SEQ_TYPE",
        "RUN_ID",
        "SOURCE_PATH",
        "DEST_PATH",
        "SIZE",
        "MOD_TIME",
        "STATUS",
        "MD5",
    ]

    run_id = input_path.name.replace(".tab.gz", "")

    logging.info("Rewriting fast5 and basecall tsv files")

    if not args.dry_run:
        fast5_table = Path("ont_fast5_table.tsv")
        basecall_table = Path("ont_basecall.tsv")

        if fast5_table.exists():
            fast5_df = pd.read_csv(fast5_table, sep="\t")
            fast5_df = fast5_df.loc[
                ~(
                    (fast5_df["SAMPLE"] == old_sample)
                    & (fast5_df["RUN_ID"] == run_id)
                )
            ]
            fast5_df.to_csv(fast5_table, sep="\t", index=False)
        else:
            logging.warning(f"Missing table, skipping: {fast5_table}")

        if basecall_table.exists():
            bc_df = pd.read_csv(basecall_table, sep="\t")
            bc_df = bc_df.loc[
                ~(
                    (bc_df["SAMPLE"] == old_sample)
                    & (bc_df["RUN_ID"] == run_id)
                )
            ]
            bc_df.to_csv(basecall_table, sep="\t", index=False)
        else:
            logging.warning(f"Missing table, skipping: {basecall_table}")

    copy_rename = new_prefix / input_path.relative_to(old_base)

    validate_no_nested_net_path(
        pd.Series([str(copy_rename)]),
        "copy_rename",
    )

    if args.dry_run:
        logging.info(f"[dry-run] write copy record: {copy_rename}")
        logging.info(f"[dry-run] write remove list: {output_file}")
    else:
        copy_rename.parent.mkdir(parents=True, exist_ok=True)
        df[out_cols].to_csv(copy_rename, sep="\t", index=False)

        with open(output_file, "w") as outfile:
            outfile.write("\n".join(migrated_files) + "\n")

    logging.info(f"Output written to {output_file}")


if __name__ == "__main__":
    main()