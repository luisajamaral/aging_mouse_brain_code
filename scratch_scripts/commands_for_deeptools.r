computeMatrix scale-regions -S ../female_split_bams/bigwigs/Female_L2_3_IT_PPP_Glut_2mo_merge.bw  ../female_split_bams/bigwigs/Female_L2_3_IT_PPP_Glut_9mo_merge.bw ../female_split_bams/bigwigs/Female_L2_3_IT_PPP_Glut_18mo_merge.bw  -R  L2_3_IT_PPP_Glut-H3K9me3-Female_filtered.bed --beforeRegionStartLength 1000 --regionBodyLength 500 --afterRegionStartLength 1000 --skipZeros -o L2_3_IT_PPP_Glut-H3K9me3-Female_filtered.matrix.mat.gz -p 14

plotHeatmap -m L2_3_IT_PPP_Glut-H3K9me3-Female.matrix.mat.gz --referencePoint center --kmeans 20 --boxAroundHeatmaps no --heatmapWidth 25  -out ExampleHeatmap1.png


plotHeatmap -m DG_Glut-H3K9me3-Female.matrix.mat.gz --referencePoint center --kmeans 20 --boxAroundHeatmaps no --heatmapWidth 25  -out SG.png


plotHeatmap -m DG_Glut-H3K9me3-Female.matrix.mat.gz \
  -out heatmap_DG.png \
  --heatmapHeight 10 \
  --zMin 0 --zMax 10 \
  --plotTitle "ATAC Signal Heatmap" \
  --samplesLabel "2mo" "9mo" "18mo"

#system("plotHeatmap -m DG_Glut-H3K9me3-Female_filtered4.matrix.mat.gz   -out heatmap_DG_filtered4.svg --heatmapHeight 10  --zMin 0 --zMax 10 --plotTitle 'ATAC Signal at H3K9me3'   --samplesLabel '2mo' '9mo' '18mo'")
