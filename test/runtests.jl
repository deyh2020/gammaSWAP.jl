using SpinQubits, LinearAlgebra, Test

@testset verbose=true "All tests" begin
	@testset verbose=true "Comparison w/ Mathematica tests" begin
		@testset "Noiseless fidelities" begin
			@test calculateFidelities(4, 0.00, 0.0, 0.0, zeros(3), 10, 3.0)[2][end] ≈ 0.0000148240466523175 rtol=1e-8
			@test calculateFidelities(4, 0.01, 0.0, 0.0, zeros(3), 10, 3.0)[2][end] ≈ 0.0018048111455146731 rtol=1e-8
			@test calculateFidelities(4, 0.05, 0.0, 0.0, zeros(3), 10, 3.0)[2][end] ≈ 0.0396766971366405 rtol=1e-8
		end
		@testset "Noisy fidelities" begin
			@test calculateFidelities(4, 0.01, 0.0, 0.0, [0.01, 0.00, 0.00], 5000, 3.0)[2][end] ≈ 0.004784666103416213 rtol=5e-2
			@test calculateFidelities(4, 0.02, 0.0, 0.0, [0.01, 0.00, 0.00], 5000, 3.0)[2][end] ≈ 0.00951159418717562 rtol=5e-2
		end
	end
	@testset "utils.jl" begin
		@testset "diagexp!" begin
			ringo = Diagonal(ones(3))
			SpinQubits.diagexp!(ringo)
			@test ringo == Diagonal(fill(exp(1.0),3))
			paul = Diagonal(zeros(3))
			SpinQubits.diagexp!(paul)
			@test paul == Diagonal(ones(3))
			for i in 1:10
				startDiag = randn(3)
				resultDiag = exp.(startDiag)
				testMat = Diagonal(startDiag)
				SpinQubits.diagexp!(testMat)
				@test testMat == Diagonal(resultDiag)
			end
		end
		@testset "zeros!" begin
			testMat = zeros(5,5)
			SpinQubits.zeros!(testMat)
			@test testMat == zeros(5,5)
			for i in 1:10
				testMat = randn(5,5)
				SpinQubits.zeros!(testMat)
				@test testMat == zeros(5,5)
			end
		end
	end

	@testset "spinors.jl" begin
		spinor1 = SpinQubits.Spinor([[1,0,0,1],[0,0,1,1]],[2.0,1.0])
		spinor2 = SpinQubits.Spinor([[1,0,0,1],[0,1,0,1]],[im,1.0])
		@testset "isequal" begin
			@test spinor1 != spinor2
			@test SpinQubits.isequal(spinor1,spinor1*1.0)
			@test spinor1 == 1.0*spinor1*1.0
			@test SpinQubits.Spinor([[1,0]],[1]) == SpinQubits.Spinor([[1,0]], [1])
		end
		@testset "+" begin
			@test spinor1 + spinor1 == SpinQubits.Spinor([[1,0,0,1],[0,0,1,1]],[4.0,2.0])
			@test spinor2 + spinor2 == SpinQubits.Spinor([[1,0,0,1],[0,1,0,1]],[2*im,2.0])
			@test spinor1 + spinor2 == SpinQubits.Spinor([[1,0,0,1],[0,0,1,1],[0,1,0,1]],[2.0 + 1.0im,1.0,1.0])
			@test spinor1 + spinor2 == spinor2 + spinor1 
			@test_throws MethodError 3.0 + spinor1
		end
		@testset "*" begin
			@test spinor1 * spinor2 ≈ 2*im
			@test spinor2 * spinor1 ≈ -2*im
			@test spinor1 * spinor1 ≈ 5.0
			@test spinor2 * spinor2 ≈ 2.0
			@test 3.0 * spinor1 == SpinQubits.Spinor([[1,0,0,1],[0,0,1,1]],[6.0,3.0])
			@test im * spinor1 == SpinQubits.Spinor([[1,0,0,1],[0,0,1,1]],[2.0im,im])
			@test (spinor2 * -im) * spinor1 == spinor2 * (im * spinor1)
		end
		@testset "normalize" begin
			@test SpinQubits.normalize(spinor1) == SpinQubits.Spinor([[1,0,0,1],[0,0,1,1]],[2.0/sqrt(5),1.0/sqrt(5)])
			@test SpinQubits.normalize(spinor2) == SpinQubits.Spinor([[1,0,0,1],[0,1,0,1]],[im/sqrt(2),1.0/sqrt(2)])
			@test SpinQubits.normalize(spinor1 + spinor2) == SpinQubits.Spinor([[1,0,0,1],[0,0,1,1],[0,1,0,1]],[(2.0 + im)/sqrt(7), 1.0/sqrt(7), 1.0/sqrt(7)])
		end
		@testset "getInitKet" begin # Sanity check, the order of "basis" is ddd ddu dud duu udd udu uud uuu
			@test SpinQubits.getInitKet(2, true) ≈ [0.0, -1/sqrt(2), 1/sqrt(2), 0.0]
			@test SpinQubits.getInitKet(2, false) ≈ [0.0, 1.0, 0.0, 0.0]
			@test SpinQubits.getInitKet(3, true)  ≈ [0.0, 0.0, 0.0, -1/sqrt(2), 0.0, 1/sqrt(2), 0.0, 0.0]
			@test SpinQubits.getInitKet(3, false) ≈ [0.0, 0.0, 0.0,       1.0,  0.0,       0.0, 0.0, 0.0]
			@test SpinQubits.getInitKet(4, true) ≈ [0.0, 0.0, 0.0,        0.0, 0.0, 0.0, 0.0, -1/sqrt(2), 
										  0.0, 0.0, 0.0, 1/sqrt(2), 0.0, 0.0, 0.0,       0.0]
			@test SpinQubits.getInitKet(4, false) ≈ [0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 1.0, 
										  0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0]
		end
	end

	@testset "tensors.jl" begin
		@testset "J tensors" begin
			@test sum(SpinQubits.getjtensor(2, 0.0),dims=3) ≈ [1.0  0.0  0.0 0.0;
																	0.0 -1.0  2.0 0.0;
																	0.0  2.0 -1.0 0.0;
																	0.0  0.0  0.0 1.0]
			@test sum(SpinQubits.getjtensor(2, 0.1),dims=3) ≈ [1.0  0.0  0.0 0.0;
																0.0 -1.0  2.0 0.0;
																0.0  2.0 -1.0 0.0;
																0.0  0.0  0.0 1.0]
			@test sum(SpinQubits.getjtensor(3, 0.01),dims=3) ≈ [2.01 0.0	0.0	0.0	0.0	0.0	0.0	0.0;
																0.0 -0.01 2.0	0.0	0.02 0.0 0.0 0.0;
																0.0 2 -1.99 0.0 2.0	0.0	0.0	0.0;
																0.0 0.0 0.0 -0.01	0.0 2 0.02	0.0;
																0.0 0.02	2.0 0.0 -0.01 0.0 0.0 0.0;
																0.0 0.0 0.0 2.0 0.0 -1.99 2.0 0.0;
																0.0 0.0 0.0 0.02 0.0 2.0 -0.01	0.0;
																0.0 0.0 0.0 0.0 0.0 0.0 0.0 2.01]
		end
	end
end