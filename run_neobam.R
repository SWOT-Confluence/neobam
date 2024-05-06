#' Run neoBAM
#'
#' Execute neoBAM on input files and write to container output directory.
#'
#' Commandline arguments (optional):
#' reach_files
#'    - name of JSON file that contains associated reach file data


# Functions
source("/app/neobam/input.R")
source("/app/neobam/neobam_functions.R")
source("/app/neobam/output.R")
source("/app/neobam/process.R")

# source("neobam/input.R")
# source("neobam/neobam_functions.R")
# source("neobam/output.R")
# source("neobam/process.R")

# Constants
# IN_DIR = file.path("/nas/cee-water/cjgleason/SWOT_Q_UMASS/mnt",  "input")
# OUT_DIR = file.path("/nas/cee-water/cjgleason/SWOT_Q_UMASS/mnt", "output")
IN_DIR = file.path("/mnt/data/input")
OUT_DIR = file.path("/mnt/data/output")
STAN_FILE = file.path("/app", "neobam", "neobam_stan_engine.stan")
# STAN_FILE = file.path( "neobam", "neobam_stan_engine.stan")


#' Create output data structure for invalid observations
#'
#' @param nt number of time steps
#'
#' @return named list of discharge and posteriors
create_invalid_out = function(nt) {
  nt_vector = rep(NA_real_, nt)
  base_discharge = list(nt_vector, nt_vector, nt_vector)
  base_posteriors = list(
    r = list(mean=NA_real_, sd=NA_real_),
    logn = list(mean=NA_real_, sd=NA_real_),
    logWb = list(mean=NA_real_, sd=NA_real_),
    logDb = list(mean=NA_real_, sd=NA_real_)
  )
  return(list(discharge=base_discharge, posteriors=list(base_posteriors, base_posteriors, base_posteriors)))
}

#' Execute neoBAM
main = function() {


    # Identify reach files to process
  start = Sys.time()
  args = commandArgs(trailingOnly=TRUE)
  set_json = ifelse(identical(args, character(0)), "/mnt/data/input/metrosets.json", args[1])
  set_index = strtoi(Sys.getenv("AWS_BATCH_JOB_ARRAY_INDEX")) + 1


  # Get Input
  in_data = get_input(sos_file,set_json,set_index)
    
  
   
    if(typeof(in_data) != "list"){
       return('dummy')}else{
    
  is_valid = check_valid(in_data)
        
    
   
  # Process
  if (is_valid == TRUE) {
  
    neobam_output = process_data(in_data, STAN_FILE)
    out_data = neobam_output
 
                
  } else {
   
    
   return('dummy')
  }
}
    
   

  # Write output
  # write_output(out_data, neobam_output$posteriors, neobam_output$posterior_Q, OUT_DIR)
  end = Sys.time()
  print(paste("Total execution time for set", set_index, ":", (end - start), "seconds."))


    return(list('neobam_output'=out_data,'time'=in_data$time,'reach_id'=in_data$reach_id))
}

neobam_output=main()
