using SpinQubits
using Test

@testset "Fake Tests" begin
	for i in 1:10
		@test i == i + 1 - 1
	end
end

# Testing spinors

# Testing operators

# Testing utils

# Testing tensors

# Testing integration
#=@testset "Noiseless fidelities" begin
	#@test calculateFidelities(4,0.01,0.0,0.0,[0.0, 0.0, 0.0],10, 0.01)[2][end] ≈ 0.0018048111455146731
	#@test calculateFidelities(4,0.01,0.0,0.0,[0.0, 0.0, 0.0],10, 0.01)[2][end] ≈ 1.4824046652317513e-5
#end=#