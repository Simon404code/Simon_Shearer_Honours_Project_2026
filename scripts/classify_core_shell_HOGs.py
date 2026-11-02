import csv
import argparse


def classify_hogs(input_file, output_file):

    with open(input_file, "r", encoding="utf-8-sig", newline="") as infile:

        reader = csv.reader(infile)

        header = next(reader)

        hog_id_column = header[0]
        haplotype_count = len(header) - 1

        with open(
            output_file,
            "w",
            encoding="utf-8",
            newline=""
        ) as outfile:

            writer = csv.writer(outfile)

            writer.writerow([
                "HOG_ID",
                "Category",
                "Haplotypes_present"
            ])

            for row in reader:

                hog_id = row[0]

                # Count haplotypes where the HOG is present
                present_count = sum(
                    1
                    for value in row[1:]
                    if value.strip() and float(value) > 0
                )

                if present_count == haplotype_count:
                    category = "Core"

                elif present_count == haplotype_count - 1:
                    category = "Soft-core"

                else:
                    category = "Shell"

                writer.writerow([
                    hog_id,
                    category,
                    present_count
                ])

    print(f"Classified HOGs written to: {output_file}")


if __name__ == "__main__":

    parser = argparse.ArgumentParser(
        description="Classify HOGs into Core, Soft-core, and Shell."
    )

    parser.add_argument(
        "--input",
        required=True,
        help="Input HOG count matrix CSV file"
    )

    parser.add_argument(
        "--output",
        required=True,
        help="Output classified HOG CSV file"
    )

    args = parser.parse_args()

    classify_hogs(
        input_file=args.input,
        output_file=args.output
    )