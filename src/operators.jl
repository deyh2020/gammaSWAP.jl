function σx(n::Int,spinor::Spinor)
    newspinor = deepcopy(spinor)
    for state in newspinor.spins
        state[n] = (state[n] + 1) % 2
    end
    return newspinor
end

function σy(n::Int,spinor::Spinor)
    newspinor = deepcopy(spinor)
    for (i,state) in enumerate(newspinor.spins)
        state[n] = (state[n] + 1) % 2
        newspinor.coefficients[i] *= spinor.spins[i][n] == 1 ? im : -im
    end
    return newspinor
end

function σz(n::Int,spinor::Spinor)
    newspinor = deepcopy(spinor)
    for i in 1:length(newspinor.coefficients)
        newspinor.coefficients[i] *= ((spinor.spins[i][n] == 1) ? 1 : -1)
    end
    return newspinor
end

σiσj(i,j,spinor::Spinor) = σx(i,σx(j, spinor)) + σy(i,σy(j, spinor)) + σz(i,σz(j, spinor))

dipoledipole(i,j,spinor::Spinor) = 1/abs(i-j)^3 * ((-1)*σx(i,σx(j, spinor)) + (-1)*σy(i,σy(j, spinor)) + 2*σz(i,σz(j, spinor)))
