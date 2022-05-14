import scanpy as sc
import sys
import sklearn as sk

#Produces a table that allows one to convert from equiv class ID to gene list.
#Params: data - anndata object representing the equiv class matrix
#filename: string path to gene_eqclass.txt or file that contains info on equiv classes

def id_to_gene_list(data,filename):
    eq_classes = [None]*data.X.toarray().shape[1]
    fo = open(filename, "r")
    fo.readline()
    fo.readline()
    for line in fo:
        toks = line.split('\t')
        eq_classes[int(toks[-1])] = toks[0:-1]
    return eq_classes                

#usage: python make_plots.py [ALEVIN FOLDER PATH] [OUTPUT_FOLDER] [PLOT_FILENAME_BASE]
#produces plots and information about the associated count and equivalence class matrices
def main():
    count_mtx_path = sys.argv[1] + "/quants_mat.mtx"
    geq_mtx_path = sys.argv[1] + "/geqc_counts.mtx"
    output_destination = sys.argv[2]
    plot_root = sys.argv[3]
    txt_file = open(output_destination + "/scores.txt", "w")
    count_mtx = sc.read_mtx(count_mtx_path)
    geq_mtx = sc.read_mtx(geq_mtx_path)
    sc.pl.highest_expr_genes(count_mtx,save="_"+plot_root+"_gene_expr.png")
    sc.pl.highest_expr_genes(geq_mtx,save="_"+plot_root+"_eq_expr.png")
    sc.pp.normalize_total(count_mtx)
    sc.pp.normalize_total(geq_mtx)
    sc.pp.log1p(count_mtx)
    sc.pp.log1p(geq_mtx)
    sc.pp.highly_variable_genes(count_mtx)
    sc.pp.highly_variable_genes(geq_mtx)
    count_mtx = count_mtx[:,count_mtx.var.highly_variable]
    geq_mtx = geq_mtx[:,geq_mtx.var.highly_variable]
    for dist in ['euclidean','l1','l2']:
        sc.pp.neighbors(count_mtx,metric=dist)
        sc.pp.neighbors(geq_mtx,metric=dist)
        sc.tl.umap(count_mtx)
        sc.tl.umap(geq_mtx)
        sc.tl.leiden(count_mtx)
        sc.tl.leiden(geq_mtx)
        sc.tl.louvain(count_mtx)
        sc.tl.louvain(geq_mtx)
        txt_file.write(dist + " Count matrix Leiden silhouette score: " + str(sk.metrics.silhouette_score(count_mtx.X,count_mtx.obs['leiden'].array)) + "\n")
        txt_file.write(dist + " Count matrix Leiden bouldin score: " + str(sk.metrics.davies_bouldin_score(count_mtx.X.toarray(),count_mtx.obs['leiden'].array)) + "\n")
        txt_file.write(dist + " Equiv matrix Leiden silhouette score: " + str(sk.metrics.silhouette_score(geq_mtx.X,geq_mtx.obs['leiden'].array)) + "\n")
        txt_file.write(dist + " Equiv matrix Leiden bouldin score: " + str(sk.metrics.davies_bouldin_score(geq_mtx.X.toarray(),count_mtx.obs['leiden'].array))+"\n")
        txt_file.write(dist + " Count matrix louvain silhouette score: " + str(sk.metrics.silhouette_score(count_mtx.X,count_mtx.obs['louvain'].array)) +"\n")
        txt_file.write(dist + " Count matrix louvain bouldin score: " + str(sk.metrics.davies_bouldin_score(count_mtx.X.toarray(),count_mtx.obs['louvain'].array)) + "\n")
        txt_file.write(dist + " Equiv matrix louvain silhouette score: " + str(sk.metrics.silhouette_score(geq_mtx.X,geq_mtx.obs['louvain'].array)) + "\n")
        txt_file.write(dist + " Equiv matrix louvain bouldin score: " + str(sk.metrics.davies_bouldin_score(geq_mtx.X.toarray(),count_mtx.obs['louvain'].array)) + "\n")
        
        sc.pl.umap(count_mtx,color='leiden',title='Count Matrix Leiden',return_fig=True).savefig(output_destination + "/plots/" + plot_root+'_'+dist+'_gene_leiden_umap.png')
        sc.pl.umap(count_mtx,color='louvain',title='Count Matrix Louvain',return_fig=True).savefig(output_destination + "/plots/" + plot_root+'_'+dist+'_gene_louvain_umap.png')
        sc.pl.umap(geq_mtx,color='leiden',title='Equivalence Class Leiden',return_fig=True).savefig(output_destination + "/plots/" + plot_root+'_'+dist+'_eqc_leiden_umap.png')
        sc.pl.umap(geq_mtx,color='louvain',title = 'Equivalence Class Louvain',return_fig=True).savefig(output_destination + "/plots/" + plot_root+'_'+dist+'_eqc_louvain_umap.png')
if __name__ == '__main__':
    main()