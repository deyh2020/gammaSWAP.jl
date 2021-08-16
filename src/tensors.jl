function getγtensor(L)
    h = zeros(2^L,2^L)
    basis = [Spinor([(reverse(digits(i, base=2, pad=L)))],[1]) for i in 0:2^L-1]
    for (i,spinorI) in enumerate(basis)
        for (j,spinorJ) in enumerate(basis)
            h[i, j] = sum( [k == l ? 0 : (-1/2*(spinorI*dipoledipole(k,l,spinorJ))) for k in 1:L, l in 1:L] )
        end
    end
    return h
end

function getjtensor(L,β; betaArray=[])
    betas = betaArray == [] ? [β^i for i in 0:L-2] : betaArray
    length(betas) != L-1 ? error("betas is wrong length.") : nothing

    jTensor = zeros(2^L, 2^L, Int(L*(L-1)/2))
    basis = [Spinor([(reverse(digits(i, base=2, pad=L)))],[1]) for i in 0:2^L-1]
    for (i,spinorI) in enumerate(basis)
        for (j,spinorJ) in enumerate(basis)
            baseIndex = 0
            for k in 1:L-1
                for n in 1:L-k
                    jTensor[i, j, baseIndex + n] = betas[k]*spinorI*σiσj(n,n+k,spinorJ)         
                end
                baseIndex += L-k
            end
        end
    end
    return jTensor
end