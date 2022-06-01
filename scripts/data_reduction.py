# read and reduce data size MIMIC
import csv
import sys
import os
import re



eventid = sys.argv[1]
new_file = sys.argv[2]
directory = re.sub("/[a-zA-Z0-9._]*$", "", new_file)

if not os.path.isdir(directory):
    sys.exit("The target directory does not exist")

with open("data/MIMIC/mimic-iv-1.0/icu/chartevents.csv", "r", encoding = "utf-8") as f:
    header = f.readline().split(",")
    counter = 0
    reduced_data = list()
    while f:
        counter += 1
        line_raw = f.readline()
        if line_raw == "":
            break
        line = line_raw.split(",")
        if len(line) < 6:
            break
        if line[5] == str(eventid):
            reduced_data.append([line[0]] + [line[3]] + line[5:7])
            sys.stdout.write("\rLine number: " + str(counter))
        if counter == 329499788:
            break

sys.stdout.write("Done, now writing output\n")
with open(new_file, 'w') as f_new:
      
    # using csv.writer method from CSV package
    write = csv.writer(f_new)
      
    write.writerow([header[0]] + [header[3]] + header[5:7])
    write.writerows(reduced_data)