#qc
adata = sc.read('raw_adata.h5ad')
adata=ov.pp.qc(adata,
                   tresh={'mito_perc': 0.1, 'nUMIs': 500, 'detected_genes': 250},
                   doublets_method='scrublet',
                   batch_key='batch')
adata
adata.var['mt'] = adata.var_names.str.startswith('MT-')
adata.var['ribo'] = adata.var_names.str.startswith(('RPS', 'RPL'))
adata.var['hb'] = adata.var_names.str.contains('^HB[^(P)]'
sc.pp.calculate_qc_metrics(
    adata, qc_vars=['mt','ribo','hb'], inplace=True, log1p=True
)
sc.pl.violin(
    adata,
    ['n_genes_by_counts','total_counts','pct_counts_mt'],
    jitter=0.4,
    multi_panel=True,
)

sc.pl.scatter(adata,'total_counts','n_genes_by_counts', color='pct_counts_mt')

adata.layers['counts'] = adata.X.copy()

sc.pp.normalize_total(adata)
sc.pp.log1p(adata)

sc.pp.highly_variable_genes(adata, n_top_genes=2000, batch_key='batch')
sc.pl.highly_variable_genes(adata)
#dimensionality reduction
sc.tl.pca(adata)

sc.pl.pca_variance_ratio(adata, n_pcs=50, log=True)

sc.pl.pca(
    adata,
    color=['sample_id','sample_name','time_from_onset','infection'],
    dimensions=[(0, 1), (2, 3), (0, 1), (2, 3)],
    ncols=2,
    size=2,
)

#BBKNN
sce.pp.bbknn(adata, batch_key='batch')
sc.tl.umap(adata)
#leiden
sc.tl.leiden(adata, flavor='leidenalg', n_iterations=2)
sc.pl.umap(adata, color=['leiden'])
#save_adata
adata.write('umap_adata.h5ad')
#rank_genes_group
sc.tl.rank_genes_groups(adata, 'leiden', method='wilcoxon')
variable_genes = pd.DataFrame(adata.uns['rank_genes_groups']['names']).head(100)
#dotplot
antivirus_list = ['ISG15','DHX58','IFIT3','OASL','IFITM3']
x = sc.pl.dotplot(adata, antivirus_list, groupby='celltype', layer='scaled',
                  var_group_rotation=0, dendrogram=False, dot_max=1,
                  dot_min=0.01, smallest_dot=10, standard_scale='var', cmap='RdBu_r'
                 )
#score
sc.tl.score_genes(adata, antivirus_score, ctrl_size=50, gene_pool=None, n_bins=25, score_name='score_01_antivirus',
                  random_state=0, copy=False, use_raw=None)