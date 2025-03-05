# Check if correct number of arguments is provided
if [ "$#" -ne 3 ]; then
    echo "Usage: $0 input.gtf output.gtf prefix"
    exit 1
fi

input_gtf="$1"
output_gtf="$2"
prefix="$3"

awk -v prefix="$prefix" '{
    if ($0 ~ /exon_id "[^"]+"/) {
        match($0, /exon_id "([^"]+)"/, arr);
        exon = arr[1];

        if (exon ~ /^ENSE/) {
            # Skip known exons that start with "ENSE"
            print;
            next;
        }

        if (exon ~ /^[^0-9]/) {
            # Already has a prefix, replace it with given prefix
            gsub(/^.*\./, prefix ".", exon);
        } else {
            # No prefix, add given prefix
            exon = prefix "." exon;
        }

        gsub(/exon_id "[^"]+"/, "exon_id \"" exon "\"");
    }
    print;
}' "$input_gtf" > "$output_gtf"

echo "Modified GTF saved to $output_gtf"
