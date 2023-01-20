using ArgParse

include("bc.jl")

function run_gap(_ARGS)
   appfolder = dirname(@__FILE__)

   function parse_commandline()
      s = ArgParseSettings(usage="##### BC-CVRP arguments #####")
      @add_arg_table! s begin
         "--instance","-i"
            help = "Instance file path"
            arg_type = String
            default = "GeneralizedAssignment/data/gapC-5-100.txt"
         "--upcutoff","-u"
            help = "upper cutoff value for the solution cost"
            arg_type = Float64
            default = 1e6
      end
      return parse_args(_ARGS, s)
   end
   app = parse_commandline()

   println("Application parameters:")
   for (arg,val) in app
      println(val === nothing ? "  $arg  =>  nothing" : "  $arg  =>  $val")
   end

   # solve the problem an get the solution
   sol, cost, total_cuts = solve_gap(app["instance"], app["upcutoff"])

   println("########################################################")
   if cost != 0 # Is there a solution?
      println("Solution $sol")
      println("Cost $cost")
   else
      println("Solution not found")
   end
   println("########################################################")

   return sol, cost, total_cuts
end
