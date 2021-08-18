#!/usr/bin/env julia

using SpinQubits, ArgParse	

function parse_commandline()
    s = ArgParseSettings()

    @add_arg_table! s begin
        "--length", "-L"
            arg_type = Int
	    default = 4
        "--singlet", "-s"
            arg_type = Bool
            default = false
        "--sigmas", "-σ"
            arg_type = Vector
	    default = zeros(3)
        "--dirspec"
	    arg_type = String
	    default = ""
	"--nReals"
	    arg_type = Int
	    default = 10000
	"--spacing"
	    arg_type = Float64
	    default = 0.03
    end

    return parse_args(s)
end

args = parse_commandline()

calculateFidelities(4, 0.01, 0.0, 0.0, zeros(3), 10, 0.01) # Compilation run

#saveFidelities( 4, 0.01, 0.1, zeros(3), 10, 0.01; singlet=false)

saveFidelities(args["length"], 0.01, 0.1, args["sigmas"], args["nReals"], args["spacing"];index=args["dirspec"], singlet=args["singlet"])
