# julia
using ArgParse
using FASTX
using BioSequences
# input parameters function
function parse_commandline()

    s = ArgParseSettings(description = "retrieves the CpG islands using the Gardiner-Garden and Frommer (1987) method and exports the id, start, end, strand, gc content and obs/exp ratio")
    @add_arg_table s begin
        "--in"
            help = "input single or multi-fasta file"
            required = true
        "--gc"
            help = "min gc% content"
            arg_type = Int
            default = 50
            required = false
        "--win"
            help = "window to slice a sequence"
            arg_type = Int
            required = true
        "--step"
            help = "step to slice a sequence"
            arg_type = Int
            required = true
        "--ratio"
            help = "min obs/exp ratio"
            arg_type = Float64
            default = 0.6
            required = false
        "--out"
            help = "output txt file"
            required = true
    end
    return parse_args(s)
end

function main()
    parsed_args = parse_commandline()
    println(parsed_args)
# create function for gc content calculation
    function gc_content(seq)
        gc = count("C",seq) + count("G",seq) + count("S",seq)
        return round(gc * 100 / (count("C",seq) + count("G",seq) + count("S",seq) + count("W",seq) + count("T",seq) + count("A",seq)), digits = 2)
    end
    # calculate ratio
    function ratio(seq)
        # calculate obs value
        obs = count("CG",seq)
        exp = count("C",seq) * count("G",seq) / parsed_args["win"]
        return round(obs/exp,  digits = 2)
    end
    # main
    open(parsed_args["out"],"a") do io
        println(io,"id","\t","start","\t","end","\t","strand","\t","%GC","\t","obs","\t","obs/exp")
    end
    reader = open(FASTA.Reader, parsed_args["in"])
    open(parsed_args["out"],"a") do io
        for record in reader
            for i in range(1, step= parsed_args["step"], length(FASTA.sequence(record)) - parsed_args["win"] + 1)
                seq = SubString(String(FASTA.sequence(record)), i,Int(i + parsed_args["win"]))
                if  gc_content(seq) >= parsed_args["gc"] && ratio(seq) >= parsed_args["ratio"]  
                    println(io,FASTA.identifier(record),"\t",i,"\t",Int(i + parsed_args["win"]),"\t","+","\t",gc_content(seq),"\t",count("CG",seq),"\t",ratio(seq))
                end
            end
        end
    end 
    close(reader)
    reader = open(FASTA.Reader, parsed_args["in"])
    open(parsed_args["out"],"a") do io
        for record in reader
            for i in range(1, step= parsed_args["step"], length(FASTA.sequence(record)) - parsed_args["win"] + 1)
                seq = SubString(String(reverse_complement!(FASTA.sequence(record))), i,Int(i + parsed_args["win"]))
                if  gc_content(seq) >= parsed_args["gc"] && ratio(seq) >= parsed_args["ratio"]  
                    println(io,FASTA.identifier(record),"\t",Int(-1 * (i - length(FASTA.sequence(record)) + parsed_args["win"])),"\t",Int(-1 * (i - length(FASTA.sequence(record)))),"\t","-","\t",gc_content(seq),"\t",count("CG",seq),"\t",ratio(seq))
                end
            end
        end
    end               
    close(reader)
end

main()
