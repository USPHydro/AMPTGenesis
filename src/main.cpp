#include "./AMPTGenesis.cpp"


#include <string>
#include <iostream>
#include <fstream>

#include "boost/program_options.hpp" //Needed for parsing input file

namespace po = boost::program_options;

AMPTGenesis *genesis_ptr; 

po::variables_map get_input_parameters(int ac, char* av[]){

    // Declare a group of options that will be
    // allowed only on command line
    std::string conf_file;
    po::options_description generic("Generic options");
    generic.add_options()
    //  ("version,v", "print version string")
        ("help,h", "produce help message")
        ("genesis-config,p", po::value<std::string>(&conf_file)->default_value("config_grid.cfg"),
                "name of a file with configuration.")
        ;

    // Declare the supported options.
    po::options_description npoints("Number of points in the grid");
    npoints.add_options()
        ("npoints.x,nx", po::value<int>()->default_value(281), "set array size in the x direction")
        ("npoints.y,ny", po::value<int>()->default_value(281), "set array size in the y direction")
        ("npoints.eta,neta", po::value<int>()->default_value(121), "set array size in the eta_s direction");

    po::options_description side_size("Physical size of the grid");
    side_size.add_options()
        ("side_size.x,Lx", po::value<double>()->default_value(2), "set physical size in the x direction (in fm)")
        ("side_size.y,Ly", po::value<double>()->default_value(2), "set physical size in the y direction (in fm)")
        ("side_size.eta,Leta", po::value<double>()->default_value(2), "set physical size in the eta_s direction");

    po::options_description smearing("Smearing parton parameters");
    smearing.add_options()
        ("smearing.K,K", po::value<double>()->default_value(1.), "a normalization value of the generated energy density")
        ("smearing.sigma_r,sigma_r", po::value<double>()->default_value(.6), "radius in transverse direction over which the energy-momentum of the parton will be scattered.")
        ("smearing.sigma_eta,sigma_eta", po::value<double>()->default_value(.6), "length in longitudinal direction over which the energy-momentum of the parton will be scattered.")
        ("smearing.tau0,tau0", po::value<double>()->default_value(.4), "the time where we will intercept AMPT parton evolution and smear parton positions")
    ;

    po::options_description paths("Smearing parton parameters");
    paths.add_options()
        //("paths.parton_history,ph", po::value<std::string>()->default_value("ana/parton-collisionsHistory.dat"),
        //                     "path to AMPT parton collision history")
        //("paths.parton_init,pi", po::value<std::string>()->default_value("ana/parton-initial-afterPropagation.dat"),
        //                     "path to AMPT initial parton positions")
        ("paths.input,input",po::value<std::string>()->default_value("/storage/home/kpala/usphydro_analysis/sources/AMPT/ana"),
                             "path where AMPT stored its results")
        ("paths.output,output",po::value<std::string>()->default_value("./AMPT_smeared_ic.csv"),
                             "path to store IC")
    ;

    po::options_description sample_radius("Multiple of the smearing radius over which we will loop. Decrease for faster execution.");
    paths.add_options()
        ("sample_radius.xy", po::value<double>()->default_value(2.), "Scan size in transverse direction (in multiple of sigma_r)")
        ("sample_radius.eta,", po::value<double>()->default_value(2.), "Scan size in longitudinal direction (in multiple of sigma_eta)")
    ;

    po::options_description cmdline_options;        //List of inputs acceptable in the comand line
    po::options_description config_file_options;    //List of inputs acceptable in the config file
    po::options_description visible;                //List of inputs visible in the help menu
    po::positional_options_description pos_args;    //Configuration of positional arguments

    cmdline_options.add(generic).add(npoints).add(side_size).add(smearing).add(paths).add(sample_radius);
    config_file_options.add(npoints).add(side_size).add(smearing).add(paths).add(sample_radius);
    visible.add(generic).add(npoints).add(side_size).add(smearing).add(paths).add(sample_radius);

    pos_args.add("genesis-config",1);
    pos_args.add("output_path",1);

    po::variables_map vm; //Keep input values

    store(po::command_line_parser(ac, av).
            options(cmdline_options).positional(pos_args).run(), vm);
    notify(vm);

    //Read command-line inputs
    if (vm.count("help")) {
        std::cout << visible << "\n";
        exit(0);
    }

    std::ifstream ifs(conf_file.c_str());
    if (!ifs){
        std::cout << "\u001b[31m [ERROR]:\u001b[0m Can not open config file: " << conf_file << "\n";
        exit( -1);
    } else {
        store(po::parse_config_file(ifs, config_file_options), vm);
        notify(vm);
    }

    return vm;
}




int main(int argc, char** argv){
  std::unique_ptr<AMPTGenesis> genesis_ptr = std::make_unique<AMPTGenesis>();
  po::variables_map vm = get_input_parameters(argc, argv);
  genesis_ptr->output_file_path = vm["paths.output"].as<std::string>();
  genesis_ptr->input_folder_path = vm["paths.input"].as<std::string>() ;
  genesis_ptr->tau0 = vm["smearing.tau0"].as<double>();
  genesis_ptr->smearing_k = vm["smearing.K"].as<double>();
  genesis_ptr->nx = vm["npoints.x"].as<int>(); 
  genesis_ptr->ny = vm["npoints.y"].as<int>(); 
  genesis_ptr->neta = vm["npoints.eta"].as<int>(); 
  genesis_ptr->Lx = vm["side_size.x"].as<double>();
  genesis_ptr->Ly = vm["side_size.y"].as<double>(); 
  genesis_ptr->Leta = vm["side_size.eta"].as<double>(); 
  genesis_ptr->sigma_r = vm["smearing.sigma_r"].as<double>(); 
  genesis_ptr->sigma_eta = vm["smearing.sigma_eta"].as<double>(); 
  genesis_ptr->rxy = vm["sample_radius.xy"].as<double>(); 
  genesis_ptr->reta = vm["sample_radius.eta"].as<double>();

  genesis_ptr->run_genesis();
  std::cout << "Saving output to file...\n";
  genesis_ptr->output_to_file();
  std::cout << "Output stored in: " << genesis_ptr->output_file_path << "\n";
  return 0;
  

}