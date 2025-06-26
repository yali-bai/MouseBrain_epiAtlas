#########    All "Joint_Cabernet" in the following code refers to Joint Cabernet.
#!/usr/bin/bash
#SBATCH -J addCompartment_subclass
#SBATCH --mem=60G
#SBATCH --cpus-per-task=2
#SBATCH --partition=compute_fat
#SBATCH -o ./logs/03-addCompartment_subclass_1.log
#SBATCH -e ./logs/03-addCompartment_subclass_1.log
#SBATCH --time=150:00:00

if [ "$#" -ne 2 ]; then
    echo "Usage: $0 <input_csv_filename> <output_csv_filename>"
    exit 1
fi
# "indir" is a custom input path, and "outdir" is a custom output path.
# indir=""
# outdir=""

input_base_path=${indir}
output_base_path=${outdir}

input_csv="${input_base_path}$1"
output_csv="${output_base_path}$2"

compartment_dir= "${indir}/h5/"
subclass_file="../../../04.data/02.metainfo/01.RNA/01.Young_Mouse/subclass_corresponding_name.csv"

cp "$input_csv" "$output_csv"

while IFS=',' read -r subclass_Joint_Cabernet subclass_Nature; do

    if [[ "$subclass_Joint_Cabernet" == "subclass_Joint_Cabernet" ]]; then
        continue
    fi

    compartment_file="${compartment_dir}${subclass_Nature}_100K_compartments.tsv"

    if [[ -f "$compartment_file" ]]; then
        echo "Processing $compartment_file for subclass $subclass_Joint_Cabernet..."

        awk 'BEGIN{OFS=","} NR>1 {chrom=$1; start=$2; region=chrom"_"int(start/100000); compartment=$5; print region, compartment}' "$compartment_file" > temp_compartment_regions.tsv

        echo "Generated temp_compartment_regions.tsv:"
        head temp_compartment_regions.tsv

        awk -v subclass="$subclass_Joint_Cabernet" -F',' 'NR==FNR {compartment[$1]=$2; next} 
            FNR==1 {print $0 "," subclass"_compartment"; next} 
            FNR>1 {region=$1; if (region in compartment) {print $0 "," compartment[region]} else {print $0 ","}}' \
            temp_compartment_regions.tsv "$output_csv" > temp_output.csv

        mv temp_output.csv "$output_csv"
    else
        echo "Warning: $compartment_file not found!"
    fi
done < "$subclass_file"

rm temp_compartment_regions.tsv

