function getγtensor(L)
    h = zeros(2^L,2^L)
    gamma = GAMMA
    basis = [Spinor([(reverse(digits(i, base=2, pad=L)))],[1]) for i in 0:2^L-1]
    for (i,spinorI) in enumerate(basis)
        for (j,spinorJ) in enumerate(basis)
            h[i, j] = sum( [k == l ? 0 : (-1/2*(spinorI*dipoledipole(k,l,spinorJ))) for k in 1:L, l in 1:L] )
        end
    end
    return h
end

function getjtensor(L,β)
    jTensor = zeros(2^L, 2^L, Int(L*(L-1)/2))
    gamma = GAMMA
    basis = [Spinor([(reverse(digits(i, base=2, pad=L)))],[1]) for i in 0:2^L-1]
    for (i,spinorI) in enumerate(basis)
        for (j,spinorJ) in enumerate(basis)
            baseIndex = 0
            for k in 1:L-1
                for n in 1:L-k
                    jTensor[i, j, baseIndex + n] = β^(k-1)*spinorI*σiσj(n,n+k,spinorJ)         
                end
                baseIndex += L-k
            end
        end
    end
    return jTensor
end