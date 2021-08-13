function readmathematica(filename)
    datastr = readline(filename)
    datastr = replace(datastr, "},{"=>"; ")
    datastr = replace(datastr, "}}"=>"]")
    datastr = replace(datastr, "{{"=>"[")
    datastr = replace(datastr, "*^"=>"e")
    datastr = replace(datastr, ","=>" ")
    dataMat = eval(Meta.parse(datastr))
    return dataMat[:,1], dataMat[:,2]
end

function saveFidelities(L,BETA,DISGAM,sigmas,nREALS,SPACING;singlet=false,format="mathematica")

    gamString = rpad(DISGAM,4,"0")
    nRealsPrime = maximum(sigmas) == 0 ? 0 : nREALS
    initString = singlet ? "s" : "up"
    sigString = string("σJ",rpad(sigmas[1],4,"0"),"_","σγ",rpad(sigmas[2],4,"0"),"_","στ",rpad(sigmas[3],4,"0"),"_",lpad(nRealsPrime,5,"0"))

    filename = joinpath(pwd(),"jdata",string(format,"_",L,"_",initString,"_β",rpad(BETA,4,"0"),"_γ",gamString,"_",sigString,"_",SPACING))

    if isfile(filename) && parse(Int,split(strip(read(`wc -c $filename`, String))," ")[1]) > 0
        println("File already found. Skipping.")
    else
        println("No file found. Calculating...")
        f = open(filename,"w")
        theseExponents = collect(range(0.0,3.0,step=SPACING))
        data = calculateFidelities(L,BETA,0.0,DISGAM,sigmas,nREALS,SPACING;singlet)
        println("Done.")
        if format=="mathematica"
            # Mathematica plotting format
            print(f, "{")
            for i in 1:length(theseExponents)
                print(f,"{",data[1][i],",",replace(string(data[2][i]),"e"=>"*^"),"}")
                i == length(theseExponents) ? print(f,"}") : print(f,",")
            end
        elseif format == "julia"
            # Julia plotting format
            write(f, data[1], data[2])
        end
        close(f)
    end
    nothing
end

function plotter!(L,BETA,DISGAM,sigmas,nREALS,SPACING;format="mathematica")

    n = length(range(0.0,3.0,step=SPACING))
    emptyarray = zeros(Float64,n,2)

    gamString = rpad(DISGAM,4,"0")
    nRealsPrime = maximum(sigmas) == 0 ? 0 : nREALS
    initString = singlet ? "s" : "up"
    sigString = string("σJ",rpad(sigmas[1],4,"0"),"_","σγ",rpad(sigmas[2],4,"0"),"_","στ",rpad(sigmas[3],4,"0"),"_",lpad(nRealsPrime,5,"0"))

    filename = joinpath(pwd(),"jdata",string(format,"_",L,"_",initString,"_β",rpad(BETA,4,"0"),"_γ",gamString,"_",sigString,"_",SPACING))
    #println(filename)
    if !isfile(filename)
        println("File is missing.")
        return nothing
    elseif format=="mathematica"
        plotdata = readmathematica(filename)
    else
        f = open(filename,"r")
        data = read!(f,emptyarray)
        close(f)
        plotdata = (data[:,1],data[:,2])
    end


    plot!(plotdata, scale=:log10, size=(500,500))
    return nothing
end