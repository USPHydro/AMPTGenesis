#include "./AMPTGenesis.cpp"


#include <string>
#include <iostream>
#include <fstream>
#include <sstream>
#include <vector>

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
        ("smearing.backpropagate", po::value<bool>()->default_value(false), "if true, partons that form after tau0 are free-streamed BACKWARD onto the tau0 surface (energy-conserving, original behaviour); if false they are dropped")
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
    po::options_description coordinates("Coordinates of the grid");
    coordinates.add_options()
        ("coordinates.system", po::value<std::string>()->default_value("hyperbolic"), "coordinate system to use (cartesian or hyperbolic)");

    po::options_description input_opts("Input options");
    input_opts.add_options()
        ("input.format", po::value<std::string>()->default_value("ampt"),
            "initial-condition input format: 'ampt' (parton collision history) or 'smash' (OSCAR2013 SMASH_IC particle list)");

    po::options_description output_opts("Output options");
    output_opts.add_options()
        ("output.output_diffusion", po::value<bool>()->default_value(false), "if true, append diffusion currents (qB, qS, qQ) to each output line")
        ("output.raw_tmunu", po::value<bool>()->default_value(false), "if true, write the raw (un-diagonalized) contravariant T^{mu nu} components + lab-frame currents instead of the Landau-matched eps/u/pi output (cuts on T^{tau tau}; energy_density_cutoff defaults to 1e-6 in this mode unless set)")
        ("output.energy_density_cutoff", po::value<double>()->default_value(0.15), "minimum energy density threshold for writing cells to output")
        ("output.sigma_scan", po::value<std::string>()->default_value(""), "comma-separated sigma list; if set, run eccentricity-ONLY scan (read partons once, loop sigmas, no IC written)")
        ("output.ecc_out", po::value<std::string>()->default_value("ecc_scan.txt"), "output file for the eccentricity scan")
        ("output.ecc_analytic", po::value<bool>()->default_value(false), "if true (with sigma_scan), use the GRID-FREE analytic eccentricity (much faster)")
        ("output.ecc_momentum", po::value<bool>()->default_value(false), "if true, write spatial eps (T^tt & Landau) and momentum anisotropy eps_p (ideal/full, mid/int) at the config sigma");


    po::options_description cmdline_options;        //List of inputs acceptable in the comand line
    po::options_description config_file_options;    //List of inputs acceptable in the config file
    po::options_description visible;                //List of inputs visible in the help menu
    po::positional_options_description pos_args;    //Configuration of positional arguments

    cmdline_options.add(generic).add(npoints).add(side_size).add(smearing).add(paths).add(sample_radius).add(coordinates).add(input_opts).add(output_opts);
    config_file_options.add(npoints).add(side_size).add(smearing).add(paths).add(sample_radius).add(coordinates).add(input_opts).add(output_opts);
    visible.add(generic).add(npoints).add(side_size).add(smearing).add(paths).add(sample_radius).add(coordinates).add(input_opts).add(output_opts);

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
  genesis_ptr->backpropagate = vm["smearing.backpropagate"].as<bool>();
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
  genesis_ptr->coordinate_system = vm["coordinates.system"].as<std::string>();
  genesis_ptr->input_format = vm["input.format"].as<std::string>();
  genesis_ptr->output_diffusion = vm["output.output_diffusion"].as<bool>();
  genesis_ptr->output_raw_tmunu = vm["output.raw_tmunu"].as<bool>();
  genesis_ptr->energy_density_cutoff = vm["output.energy_density_cutoff"].as<double>();
  // Raw-T^{mu nu} mode keys the cutoff on T^{tau tau} only to drop empty cells; the
  // Landau e/p/pi filtering does not apply, so use a small value (keep every non-empty
  // cell) unless the user set energy_density_cutoff explicitly.
  if (genesis_ptr->output_raw_tmunu && vm["output.energy_density_cutoff"].defaulted()) {
    genesis_ptr->energy_density_cutoff = 1e-6;
    std::cout << "[INFO]: raw_tmunu mode: energy_density_cutoff defaulted to "
              << genesis_ptr->energy_density_cutoff << " (drop empty cells only)\n";
  }

  // Spatial-vs-momentum anisotropy comparison (grid, single sigma)
  if (vm["output.ecc_momentum"].as<bool>()) {
    genesis_ptr->run_ecc_momentum(vm["output.ecc_out"].as<std::string>());
    std::cout << "Momentum-anisotropy comparison stored in: " << vm["output.ecc_out"].as<std::string>() << "\n";
    return 0;
  }

  // Fast eccentricity-only sigma scan: if output.sigma_scan is set, read partons once
  // and loop over the comma-separated sigma list, writing only eccentricities.
  std::string sigma_scan = vm["output.sigma_scan"].as<std::string>();
  if (!sigma_scan.empty()) {
    std::vector<double> sigmas;
    std::stringstream ss(sigma_scan); std::string tok;
    while (std::getline(ss, tok, ',')) { if (!tok.empty()) sigmas.push_back(std::stod(tok)); }
    std::string ecc_out = vm["output.ecc_out"].as<std::string>();
    if (vm["output.ecc_analytic"].as<bool>())
      genesis_ptr->run_ecc_analytic(sigmas, ecc_out);   // grid-free
    else
      genesis_ptr->run_ecc_scan(sigmas, ecc_out);       // grid-based
    std::cout << "Eccentricity scan stored in: " << ecc_out << "\n";
    return 0;
  }

  genesis_ptr->run_genesis();
  // In raw-T^{mu nu} mode run_genesis already wrote the file (it holds no
  // diagonalized vectors to re-dump), so only re-dump in the standard path.
  if (!genesis_ptr->output_raw_tmunu) {
    std::cout << "Saving output to file...\n";
    genesis_ptr->output_to_file();
  }
  std::cout << "Output stored in: " << genesis_ptr->output_file_path << "\n";
  return 0;
  

}