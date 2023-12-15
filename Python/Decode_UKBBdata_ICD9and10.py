
#Firoj Mahmud, Post-doc karolinska Institute, firoj.mahmud@ki.se
##########This is a pyhton code that uses UKB subset data regarding icd9 and icd10 and decode the icd9 and 10 based on the supplied UKB data code 3 and 87
#for further mail to Fiorj mahmud firoj.mahmud@ki.se

import csv
import os
import time

def decode_codes(encoded_file, decoding_files, decoded_file):
    decoding_info = {}
    column_prefixes = {
        'coding87.tsv': ['diagnosis_main', 'diagnoses_secondary'],
        'coding3.tsv': ['type_of_cancer']
    }

    # Read decoding files and populate decoding_info dictionary
    for decoding_file, prefixes in column_prefixes.items():
        with open(decoding_file, 'r') as file:
            reader = csv.reader(file, delimiter='\t')
            for row in reader:
                for prefix in prefixes:
                    code = prefix + row[0]
                    decode_value = row[1]
                    decoding_info[code] = decode_value

    # Decode the encoded file
    with open(encoded_file, 'r') as input_file, open(decoded_file, 'w', newline='') as output_file:
        reader = csv.reader(input_file)
        writer = csv.writer(output_file)

        header_row = next(reader)
        writer.writerow(header_row)

        total_rows = sum(1 for _ in reader)
        input_file.seek(0)

        start_time = time.time()
        processed_rows = 0

        for row in reader:
            eid = row[0]
            decoded_values = []

            for i, value in enumerate(row[1:]):
                column_header = header_row[i+1]

                for decoding_file, prefixes in column_prefixes.items():
                    for prefix in prefixes:
                        if column_header.startswith(prefix):
                            code = column_header.replace(prefix, '')
                            code_with_prefix = prefix + code

                            if code_with_prefix in decoding_info:
                                decoded_values.append(decoding_info[code_with_prefix])
                            else:
                                decoded_values.append(value)
                            break
                    else:
                        continue
                    break
                else:
                    decoded_values.append(value)

            writer.writerow([eid] + decoded_values)

            processed_rows += 1
            if processed_rows % 1000 == 0:
                elapsed_time = time.time() - start_time
                print(f"Processed {processed_rows}/{total_rows} rows. Elapsed time: {elapsed_time:.2f} seconds.")

        elapsed_time = time.time() - start_time
        print(f"Decoding completed. Total rows processed: {processed_rows}. Elapsed time: {elapsed_time:.2f} seconds.")

# Update the file names accordingly:
encoded_file = 'ukb_Matrix_Firoj_ControlHealthyICD9and10codes_only.csv'  # File to be decoded
decoding_files = ['coding87.tsv', 'coding3.tsv']  # Coding files here
decoded_file = 'Decoded_ukb_Matrix_Firoj_ControlHealthyICD9and10codes_only.csv'  # Decoded file name

decode_codes(encoded_file, decoding_files, decoded_file)
