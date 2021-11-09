import sys
import numpy as np
import paraparser.paraparser as pp
# from mymath.matrix import visualize
from fileparser.parse47 import parse_47
from filegenerator.genfchk import quicksave

def COT(d, forder=False):
    amoao, bmoao = d['trafomo']
    dim = amoao.shape[0]
    assert bmoao.shape[0] == dim
    admao, bdmao = d['density']
    Na = amoao.dot(admao).dot(amoao.T).trace().round().astype(int)
    Nb = bmoao.dot(bdmao).dot(bmoao.T).trace().round().astype(int)
    assert Na >= Nb
    # print(Na, Nb)
    acoamo = np.eye(dim)
    bcobmo = np.eye(dim)
    D2 = np.zeros(dim)
    sao = d['overlap']
    sab = amoao.dot(sao).dot(bmoao.T)
    # visualize(sab)

    U1, D21, V1 = np.linalg.svd(sab[:Na,:Nb], full_matrices=True)
    # print(D21)
    acoamo[:Na] = U1.T.dot(acoamo[:Na])
    bcobmo[:Nb] = V1.dot(bcobmo[:Nb])
    D2[:Nb] = D21
    sab[:Na] = U1.T.dot(sab[:Na])
    sab[:,:Nb] = sab[:,:Nb].dot(V1.T)
    assert np.allclose(sab, acoamo.dot(amoao).dot(sao).dot(bmoao.T).dot(bcobmo.T))
    # visualize(sab)

    if Na > Nb:
        U2, D22, V2 = np.linalg.svd(sab[Nb:Na,Nb:], full_matrices=True)
        # print(D22)
        acoamo[Nb:Na] = U2.T.dot(acoamo[Nb:Na])
        bcobmo[Nb:] = V2.dot(bcobmo[Nb:])
        D2[Nb:Na] = D22
        sab[Nb:Na] = U2.T.dot(sab[Nb:Na])
        sab[:,Nb:] = sab[:,Nb:].dot(V2.T)
        assert np.allclose(sab, acoamo.dot(amoao).dot(sao).dot(bmoao.T).dot(bcobmo.T))
        # visualize(sab)

    U3, D23, V3 = np.linalg.svd(sab[Na:,Na:], full_matrices=True)
    # print(D23[::-1])
    acoamo[Na:] = U3.T[::-1].dot(acoamo[Na:])
    bcobmo[Na:] = V3[::-1].dot(bcobmo[Na:])
    D2[Na:] = D23[::-1]
    sab[Na:] = U3.T[::-1].dot(sab[Na:])
    sab[:,Na:] = sab[:,Na:].dot(V3[::-1].T)
    assert np.allclose(sab, acoamo.dot(amoao).dot(sao).dot(bmoao.T).dot(bcobmo.T))
    # visualize(sab)

    acoao = acoamo.dot(amoao)
    bcoao = bcobmo.dot(bmoao)

    if forder:
        acoao[:Na] = np.concatenate((acoao[Nb:Na], acoao[:Nb][::-1]))
        bcoao[:Na] = np.concatenate((bcoao[Nb:Na], bcoao[:Nb][::-1]))
        D2[:Na] = np.concatenate((D2[Nb:Na], D2[:Nb][::-1]))

    return acoao, bcoao, D2, D2[max(0,Nb-10):min(dim,Na+10)]

def SNO(d):
    amoao, bmoao = d['trafomo']
    admao, bdmao = d['density']
    sdmao = admao - bdmao
    sdmamo = amoao.dot(sdmao).dot(amoao.T)
    assert np.allclose(sdmamo, sdmamo.T)
    vals, vecs = np.linalg.eigh(sdmamo)
    vecs = vecs.T
    order = abs(vals).argsort()[::-1]
    vals, sedoamo = vals[order], vecs[order]
    sedoao = sedoamo.dot(amoao)
    return sedoao, vals

def genSNO(fn47, mfn=None, mode='SNO', silent=False, forder=False):
    print(fn47)
    mode = mode.upper()
    assert mode in ('SNO', 'COT')
    d = parse_47(fn47, silent=True)
    assert d['dim'][0] == 2

    if mode == 'SNO':
        sedoao, vals = SNO(d)
        if not silent:
           print(vals[:20])
        if mfn is not None:
            quicksave(mfn, sedoao, vals, suffix='_'+mode)

    elif mode == 'COT':
        acoao, bcoao, D2, D2p = COT(d, forder=forder)
        if not silent:
           print(D2p)
        if mfn is not None:
            quicksave(mfn, np.array((acoao, bcoao)), np.array((D2, D2)), suffix='_'+mode)

    else:
        raise UserWarning('Unrecognized keyword: mode=%s' % mode)

def main():
    paras = pp.main()
    fn47 = paras.fns[1]
    fnfchk = paras.fns[2] if paras.fns[2:] else None
    mode = 'SNO'
    if 'COT' in [_.upper() for _ in paras.rest] or 'c' in paras.args:
        mode = 'COT'
    if 'SNO' in [_.upper() for _ in paras.rest] or 's' in paras.args:
        mode = 'SNO'
    paras.kwargs.setdefault('mode', mode)
    if 'f' in paras.args:
        paras.kwargs.setdefault('forder', True)
    genSNO(fn47, fnfchk, **paras.kwargs)

if __name__ == '__main__':
    main()


