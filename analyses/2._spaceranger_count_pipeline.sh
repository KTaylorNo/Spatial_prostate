# To run the spaceranger count pipeline for each tissue section.

module load Miniconda3/22.11.1-1

spaceranger count --id="Pipeline_output_directory" \
  --transcriptome=</Path/To/Reference/refdata-gex-GRCh38-2020-A> \
  --fastqs=</Path/To/fastqs> \
  --image=</Path/To/hires_tissue_image> \
  --slide=V10A20-121 \
  --area=A1 \
  --localcores=16 \
  --localmem=128
