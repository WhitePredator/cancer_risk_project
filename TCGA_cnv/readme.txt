seg_clinic为tcga的cnv及临床数据

cnv.R为运行code



cnv_result.Rdata为cnv（拷贝数变化）结果：
cancer_count为每个癌种所有样本中基因出现拷贝数变化的频率统计。以blca为例，blca中基因名结尾加_0表示为该基因拷贝数降低的事件，结尾加_1表示该基因拷贝数增加的事件。blca_bind结果为不考虑该cnv变化是增加还是减少，均记为1变化。
cancer_pvalue_list为对cancer_count进行的卡方检验p值及fdr调整过的p值。_0及_1同上述表示。
threshold_up=0.848      ##segment_mean的阈值，超过上界认为拷贝数增加，低于下节认为拷贝数下降。
threshold_low=-0.737   ##阈值参考自doi:  10.1038/nature06358

blca.xlsx为blca的结果，存成xlsx方便查看。（示例作用，更多结果在Rdata文件中）