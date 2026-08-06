from itertools import combinations
from pathlib import Path
import numpy as np
from numba import njit

P=11
INV=np.array([0]+[pow(x,P-2,P) for x in range(1,P)],dtype=np.int64)
THETA=np.load(Path(__file__).with_name('theta_ame1811.npy')).astype(np.int64)%P
COMB=np.array(list(combinations(range(17),8)),dtype=np.int64)

@njit(cache=False)
def rank_mod(B):
    M=B.copy()
    rows,cols=M.shape
    r=0
    for c in range(cols):
        piv=-1
        for i in range(r,rows):
            if M[i,c]%P:
                piv=i; break
        if piv<0: continue
        if piv!=r:
            for j in range(cols):
                tmp=M[r,j]; M[r,j]=M[piv,j]; M[piv,j]=tmp
        iv=INV[M[r,c]%P]
        for j in range(c,cols): M[r,j]=(M[r,j]*iv)%P
        for i in range(r+1,rows):
            f=M[i,c]%P
            if f:
                for j in range(c,cols): M[i,j]=(M[i,j]-f*M[r,j])%P
        r+=1
        if r==rows: break
    return r

@njit(cache=False)
def certify_one(T, combs):
    bad=0
    n=17
    mask=np.zeros(n,dtype=np.uint8)
    B=np.empty((8,9),dtype=np.int64)
    for z in range(combs.shape[0]):
        for i in range(n): mask[i]=0
        for i in range(8): mask[combs[z,i]]=1
        cc=0
        for j in range(n):
            if mask[j]==0:
                for i in range(8): B[i,cc]=T[combs[z,i],j]
                cc+=1
        if rank_mod(B)<8: bad+=1
    return bad

# trigger compile on one matrix
T0=np.delete(np.delete(THETA,17,axis=0),17,axis=1)
print('primary deficient:',certify_one(T0,COMB),'/',len(COMB))
stored=np.load(Path(__file__).with_name('theta_ame1711.npy')).astype(np.int64)%P
print('stored theta_ame1711 matches deletion of last vertex:',np.array_equal(stored,T0))
ok=0
for v in range(18):
    Tv=np.delete(np.delete(THETA,v,axis=0),v,axis=1)
    bad=certify_one(Tv,COMB)
    print(f'deletion {v+1:2d}: deficient={bad}')
    ok += bad==0
print(f'all-deletion census: {ok}/18')
assert ok==18
print('ALL CHECKS PASSED')
