import csv




def flag_conflict(pangolin_output,lineage_of_interest,scorpio_call):
    with open(pangolin_output, "r") as f:
        reader= csv.DictReader(f)
        conflicting = {}
        for row in reader:
            if row["lineage"] == lineage_of_interest:
                if row["scorpio_call"] != scorpio_call:
                    print(row["taxon"],row["lineage"],row["scorpio_call"])

            elif row["scorpio_call"] ==scorpio_call:
                pass
                # print(row["taxon"],row["lineage"],row["scorpio_call"])
flag_conflict("/Users/s1680070/repositories/quokka/lineage_report.csv","B.1.351","cB.1.351")