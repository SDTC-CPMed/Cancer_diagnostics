#!/bin/bash

# Function to print progress messages
#The script written by firoj mahmud (firoj.mahmud@ki.se)
# sh <script> <md5sum> <encrypted_file> <key_file>
print_progress() {
    echo "[$(date +'%Y-%m-%d %H:%M:%S')] $1"
}

# Function to check MD5 sum
check_md5sum() {
    expected_md5="$1"
    calculated_md5=$(md5sum "$2" | awk '{print $1}')
    if [ "$expected_md5" != "$calculated_md5" ]; then
        echo "MD5 sum mismatch! Expected: $expected_md5, Calculated: $calculated_md5"
        exit 1
    else
        echo "MD5 sum check passed."
    fi
}

# Display help message
display_help() {
    echo "Usage: $0 <md5sum> <encrypted_file> <key_file>"
    echo "Check MD5 sum, decrypt, and perform conversions."
    echo "Arguments:"
    echo "  <md5sum>          Expected MD5 sum for the encrypted file."
    echo "  <encrypted_file>  Encrypted file name."
    echo "  <key_file>        Key file name."
    exit 1
}

# Check if required files exist
check_files_exist() {
    if ! [ -f "encoding.dat" ] || ! [ -f "ukbconv" ] || ! [ -f "ukbunpack" ]; then
        echo "Required files (encoding.dat, ukbconv, ukbunpack) not found!"
        exit 1
    fi
}

# Start time
start_time=$(date +%s)

# Check if input arguments are provided
if [ $# -eq 0 ]; then
    display_help
fi

# Check if the user wants to display help
if [ "$1" == "-h" ] || [ "$1" == "--help" ]; then
    display_help
fi

# Check if required files exist
check_files_exist

# Assign input arguments
expected_md5="$1"
enc_file="$2"
key_file="$3"

# Check MD5 sum
print_progress "Checking MD5 sum..."
check_md5sum "$expected_md5" "$enc_file"

# Decrypt the files
print_progress "Decrypting the files..."
./ukbunpack "$enc_file" "$key_file"

# Create CSV files
print_progress "Creating CSV files..."
./ukbconv "${enc_file}_ukb" csv -eencoding.dat

# Create R files
print_progress "Creating R files..."
./ukbconv "${enc_file}_ukb" r -eencoding.dat

# Create docs and HTML files
print_progress "Creating docs and HTML files..."
./ukbconv "${enc_file}_ukb" docs -eencoding.dat

# Calculate and print the running time
end_time=$(date +%s)
running_time=$((end_time - start_time))
print_progress "Job done. Running time: $running_time seconds."

