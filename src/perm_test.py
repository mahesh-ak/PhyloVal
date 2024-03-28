import numpy as np
from tqdm import tqdm
from io import StringIO
from Bio import Phylo

def tree_str(tree_node, paren_dist):
    length = max(paren_dist - tree_node['dist'], 0.0001)
    ## leaf
    if len(tree_node['nodes']) == 1 and type(tree_node['nodes'][0]) == str:
        name = tree_node['nodes'][0]#.replace('(','-').replace(')','').replace("'",'`')
        return f"{name}:{length}"
    child_str = []
        
    for child_d in tree_node['nodes']:
        txt = tree_str(child_d, tree_node['dist'])
        child_str.append(txt)
    
    return f"({','.join(child_str)}){round(100*tree_node['support'],1)}:{length}"
    
def perm_cluster(word_mat, langs, num_concepts, compare_func, num_runs=100):
    eps = 1e-12
    num_langs = len(langs)
    num_clusters = num_langs
    score_mean = 0.0
    score_max = 0.0
    p_val_max = 0.0
    ret_dict = {'bi_scores': None, 'bi_pval': None, 'last_score': None, 'last_pval': None}
    dist_mat = np.zeros((num_langs, num_langs), dtype=np.float64)
    rnd_dist_mat = np.zeros((num_runs, num_langs, num_langs, num_concepts), dtype=np.float64)
    clusters = [[i] for i in range(num_langs)]
    tree = {'nodes': [{'nodes':[l], 'support': None, 'dist': 0} for l in langs],'support': None, 'dist': 1}
    
    ## Compute distance matrix according to similarity function (num_langs x num_langs)
    for i in range(num_langs):
        for j in range(0,i):
            dists = [compare_func(w1,w2) for w1,w2 in zip(word_mat[i,:],word_mat[j,:])]
            dists = np.array(dists)
            dists = np.ma.masked_where(dists < 0, dists)
            dist_mat[i,j] = dists.mean()
            dist_mat[j,i] = dist_mat[i,j]
    
    ## Compute distance matrices for random permutations (num_runs x num_langs x num_langs)
    print("Computing random scores:")
    for r in tqdm(range(num_runs)):
        word_mat_rnd = word_mat.copy()
        for l_i in range(num_langs):
            word_mat_rnd[l_i,:] = word_mat_rnd[l_i,:][np.random.permutation(num_concepts)]
        for i in range(num_langs):
            for j in range(0,i):
                for c in range(num_concepts):
                    rnd_dist_mat[r,i,j,c] = compare_func(word_mat_rnd[i,c], word_mat_rnd[j,c])
                    rnd_dist_mat[r,j,i,c] = rnd_dist_mat[r,i,j,c]
    rnd_dist_mat = np.ma.masked_where(rnd_dist_mat < 0, rnd_dist_mat)
    rnd_dist_mat = rnd_dist_mat.mean(axis=-1)
    
    ## Compute initial similarity matrix along with support probabilities
    p_val_mat = np.ones((num_clusters,num_clusters),dtype=np.float64)
    score_mat = -1*np.ones((num_clusters,num_clusters),dtype=np.float64)
    for i in range(num_clusters):
        for j in range(0,i):
            cluster_i = clusters[i]
            cluster_j = clusters[j]
            score_mat[i,j] = dist_mat[cluster_i,:][:,cluster_j].mean()

            scores = rnd_dist_mat[:,cluster_i,:][:,:,cluster_j].mean(axis=(1,2))
            p_val_mat[i,j] = (scores < score_mat[i,j]).astype(np.int64).mean()
            score_mean = scores.mean()
            score_mat[i,j] = (score_mean - score_mat[i,j]) / (score_mean + eps)
        
    # Bilateral similarities
    ret_dict['bi_scores'] = score_mat.copy()
    ret_dict['bi_pval'] = p_val_mat.copy()
    
    while (num_clusters > 1):
            
        i_max, j_max = np.unravel_index(score_mat.argmax(), score_mat.shape)
        c_i, c_j = clusters[i_max], clusters[j_max]
        t_i, t_j = tree['nodes'][i_max], tree['nodes'][j_max]
        clusters = [c for i,c in enumerate(clusters) if i != i_max and i != j_max]
        tree['nodes'] = [c for i,c in enumerate(tree['nodes']) if i != i_max and i != j_max]
        cluster_inds = [i for i in range(num_clusters) if i != i_max and i != j_max]
        
        c = c_i + c_j
        clusters.append(c)
        num_clusters -= 1
        score_max = score_mat[i_max,j_max]
        p_val_max = p_val_mat[i_max,j_max]
        tree['nodes'].append({'nodes': [t_i,t_j],'support': 1-p_val_max,'dist': 1-score_max})
        
        if num_clusters == 1:
            break
            
        new_p_val_mat = np.ones((num_clusters,num_clusters),dtype=np.float64)
        new_score_mat = -1*np.ones((num_clusters,num_clusters),dtype=np.float64)
        
        ## Update new similarity matrix
        new_p_val_mat[:num_clusters-1,:num_clusters-1] = p_val_mat[cluster_inds,:][:,cluster_inds]
        new_score_mat[:num_clusters-1,:num_clusters-1] = score_mat[cluster_inds,:][:,cluster_inds]
        
        for j in range(0,num_clusters-1):
            cluster_j = clusters[j]
            new_score_mat[-1,j] = dist_mat[c,:][:,cluster_j].mean()

            scores = rnd_dist_mat[:,c,:][:,:,cluster_j].mean(axis=(1,2))
            new_p_val_mat[-1,j] = (scores < new_score_mat[-1,j]).astype(np.int64).mean()
            score_mean = scores.mean()
            new_score_mat[-1,j] = (score_mean - new_score_mat[-1,j]) / (score_mean + eps)
        
    
        del p_val_mat, score_mat
        p_val_mat = new_p_val_mat
        score_mat = new_score_mat

    ret_dict['tree'] = tree_str(tree['nodes'][0], tree['dist'])
    #print(ret_dict['tree'])
    newick = Phylo.read(StringIO(ret_dict['tree']), "newick")
    Phylo.draw_ascii(newick)

    ret_dict['last_pval'] = p_val_max
    ret_dict['last_score'] = max(score_max, 0.0)
    print('Final p_val:',p_val_max)
    return ret_dict