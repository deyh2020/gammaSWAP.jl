const GAMMA = 6.796381511118120480641144794583e-5;

function Ham!(ham,L,jt,γm,jcouplings,γ)
    @views begin
        for i in 1:Int(L*(L-1)/2)
            ham .= ham .+ jcouplings[i].*jt[:,:,i]
        end
        ham .= ham .+ γ.*γm
        nothing
    end
end

function diagexp!(array, result)
    for i in 1:length(@view result[:,1])
        result[i,i] = exp(array[i,i])
    end
    nothing
end

function zeros!(A)
    for i in eachindex(A)
        A[i] = 0.0
    end
    nothing
end