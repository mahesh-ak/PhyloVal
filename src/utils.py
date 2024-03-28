from lingpy import *

def create_arr(nos_words):
    p1_dolgo = [''] + [tokens2class(ipa2tokens(word),'dolgo')[0] for word in nos_words[1:]]
    turchin_arr = [''] + [''.join(tokens2class(ipa2tokens(word),'dolgo')[:3]) for word in nos_words[1:]]
    turchin_arr = [''] + [t.replace('V','') if t[0] != 'V' else t[:2] for t in turchin_arr[1:]]
    
    return p1_dolgo, turchin_arr

def identity(id1, id2, arr):
    if id1 == 0 or id2 == 0:
        return -1
    
    return 0 if arr[id1] == arr[id2] else 1

def lexstat_score(id1,id2,wl,method='lexstat'):
    if id1 == 0 or id2 == 0:
        return -1
    score = wl.align_pairs(id1,id2, method=method, pprint=False)[2]
    return score 

def get_cognates_wl(lng1,lng2,word_mat,langs,wl,method='sca'):
    l1 = langs.index(lng1)
    l2 = langs.index(lng2)
    cogs = []
    if method == 'lexstat':
        th = 0.7
    else:
        th = 0.55
    
    for id1, id2 in zip(word_mat[l1,:], word_mat[l2,:]):
        if id1 == 0 or id2 == 0:
            continue
        w1, w2, dist = wl.align_pairs(id1,id2, method=method, pprint=False)
        if dist < th:
            cogs.append((w1,w2,round(dist,3)))
    return sorted(cogs,key=lambda x: x[2])

def get_cognates(lng1,lng2,word_mat,langs,nos_words,arr):
    l1 = langs.index(lng1)
    l2 = langs.index(lng2)
    cogs = []
    for id1, id2 in zip(word_mat[l1,:], word_mat[l2,:]):
        if id1 == 0 or id2 == 0:
            continue
        w1, w2, dist = nos_words[id1], nos_words[id2], identity(id1, id2, arr=arr)
        if dist < 0.5:
            cogs.append((w1,w2,dist))
    return sorted(cogs,key=lambda x: x[2])

def sort_by_similarity(lang, mat, langs, p_mat):
    ind = langs.index(lang)
    dists = mat[ind,:]
    inds = dists.argsort()
    return [(langs[i], round(1-dists[i],3), p_mat[ind,i] < 0.05, p_mat[ind,i]) for i in inds]
