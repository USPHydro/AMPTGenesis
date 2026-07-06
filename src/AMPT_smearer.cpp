
#include "boost/multi_array.hpp"

#include "TParticlePDG.h"
#include "TDatabasePDG.h"

#include "TFile.h"
#include "TH1.h"
#include "TH2.h"

#include <string>
#include <vector>
#include <sstream>
#include <fstream>
#include <iostream>
#include "./parton_collision.cpp"


//#include "./progress_bar.cpp"

typedef std::array<std::array<double,4>,4> Mat4x4;
typedef std::vector<std::string> svec;

///\class AMPTSmearer
///\brief Tools to read the AMPT text files and creates the appropriate IC
class AMPTSmearer
{
private:
    std::string coordinates;
    std::string results_path;   ///< input folder (AMPT) or SMASH OSCAR file/dir
    svec cols_hist;
    svec init_parton;

    //Helpers for dealing with strings in a similar way as python
    svec readlines(std::string path) const;
    svec split(std::string input_string, char delimiter) const;
    void get_dnde(std::string results_path);
    PartonThermalized free_streamer(PartonCollision last_collision, double tau_f) const;
    double get_frac(double lower_bound, double upper_bound,
                                 double central_eta, double bin_width);
    int npartons;       ///< Number of partons in the event
    int ncollisions;    ///< Number of parton-parton collision (do not confuse with Ncoll, the number of nucleon-nucleon collisions)
    std::vector<std::vector<PartonCollision>> parton_histories; ///< Collision history of each parton
    const double e_charge = 0.30282212077; ///< Elementary charge in natural units (sqrt(4*pi*alpha))

public:
    std::vector<PartonThermalized> thermalized_partons;  // public: read by analytic ecc
    double impact_parameter;  ///< Impact parameter of the event
                             /// See https://arxiv.org/pdf/1910.08004.pdf
                             /// These values are outputted by AMPT
    double refmult1;          ///< Nch in |eta| < 0.5
    double refmult2;          ///< Nch in 0.5 < |eta| < 1.0
    double refmult3;          ///< Nch in |eta| < 1
    double FwdAll;            ///< Nch in 2.1 < |eta| < 5.1
    double Fwd1;              ///< Nch in 2.1 < |eta| < 3.0
    double Fwd2;              ///< Nch in 3.0 < |eta| < 4.0
    double Fwd3;              ///< Nch in 4.0 < |eta| < 5.0
    double Npart;             ///< Number of partipating nucleons
    double NpartTarg;         ///< Number of partipating nucleons in target
    double NpartProj;         ///< Number of partipating nucleons in projectile
    double NpartTargElastic;  ///< Number of elastically partipating nucleons in projectile
    double NpartProjElastic;
    int nx; int ny; int neta; double Lx;double Ly; double Leta; double dx; double dy; double deta;
                                double sigma_r; double sigma_eta; double tau0; double rxy; double reta;
    boost::multi_array<Mat4x4,3> Tmunu;
    boost::multi_array<double,3> rhob;
    boost::multi_array<double,3> j0;
    boost::multi_array<double,3> j1;
    boost::multi_array<double,3> j2;
    boost::multi_array<double,3> j3;
    boost::multi_array<double,3> j0e;
    boost::multi_array<double,3> j1e;
    boost::multi_array<double,3> j2e;
    boost::multi_array<double,3> j3e;
    boost::multi_array<double,3> j0s;
    boost::multi_array<double,3> j1s;
    boost::multi_array<double,3> j2s;
    boost::multi_array<double,3> j3s;


    AMPTSmearer(std::string results_path, double K,int nx_, int ny_, int neta_, double Lx_,double Ly_, double Leta_,
                                double sigma_r_, double sigma_eta_, double tau0_, double rxy_, double reta_, std::string coordinate_system);
    double K;
    ~AMPTSmearer();

    void parse_history();
    // Read a SMASH OSCAR2013 "SMASH_IC" particle list and fill thermalized_partons
    // directly. The particles are already on the IC hypersurface, so there is no
    // free-streaming step (no parse_history()/propagate()).
    void parse_smash();
    void propagate(double tau_f);
    void fill_Tmunu(double sr,double seta);

    std::array<double,4> net_p;

};

///\brief read a file and organizes its lines in a vector
///\param path The path where AMPT results are stored
svec AMPTSmearer::readlines(std::string path) const
{
    svec line_vecs;
    //Adapted from https://stackoverflow.com/a/7868998
    std::ifstream infile(path);
    std::string line;
    while (std::getline(infile, line))
        line_vecs.push_back(line);

    infile.close();
    return line_vecs;
}

///\brief Split a string using a delimiter
///\param input_string The string that we desire to split
///\param delimiter The delimiter that will be used to split the string
svec AMPTSmearer::split(std::string input_string, char delimiter) const
{
    svec split_string;
    //Adapted from https://stackoverflow.com/a/5167799
    std::istringstream stream_in(input_string);
    std::string line;
    while (std::getline(stream_in, line,delimiter))
        if(line.size() != 0)
            split_string.push_back(line);

    return split_string;
}


AMPTSmearer::AMPTSmearer(std::string results_path,double Kin,int nx_, int ny_, int neta_, double Lx_,double Ly_, double Leta_,
                                double sigma_r_, double sigma_eta_, double tau0_, double rxy_, double reta_, std::string coordinate_system):
refmult1(.0),
refmult2(.0),
refmult3(.0),
FwdAll(.0),
Fwd1(.0),
Fwd2(.0),
Fwd3(.0),
Tmunu(boost::extents[1][1][1]),
rhob(boost::extents[1][1][1]),
j0(boost::extents[1][1][1]),
j1(boost::extents[1][1][1]),
j2(boost::extents[1][1][1]),
j3(boost::extents[1][1][1]),
j0e(boost::extents[1][1][1]),
j1e(boost::extents[1][1][1]),
j2e(boost::extents[1][1][1]),
j3e(boost::extents[1][1][1]),
j0s(boost::extents[1][1][1]),
j1s(boost::extents[1][1][1]),
j2s(boost::extents[1][1][1]),
j3s(boost::extents[1][1][1]),
net_p({0, 0, 0, 0}),
coordinates(coordinate_system)
{
    nx = nx_; ny=ny_; neta=neta_; Lx=Lx_; Ly=Ly_; Leta=Leta_; 
                                sigma_r=sigma_r_; sigma_eta = sigma_eta_; tau0=tau0_; rxy=rxy_; reta = reta_;
    K = Kin;

    // The AMPT text files are read in parse_history(); the SMASH input
    // path (parse_smash) never touches them. Storing the path lets one
    // constructor serve both input formats.
    this->results_path = results_path;
}

void AMPTSmearer::get_dnde(std::string results_path){
    double const bin_width = .4;
    svec dn_de_data = readlines(results_path+"/dnde_ch.dat");
    for (std::string data_str : dn_de_data){
        svec data = split(data_str, ' ');
        if(data.size() == 2){
            double eta = atof(data[0].data());
            double n =  atof(data[1].data())*bin_width;

            // Central region
            if ( (eta - bin_width*.5 < 1.) && (eta + bin_width*.5 > -1.) ){ // |eta| < 1
                refmult3+= get_frac(-1.,1.,eta,bin_width)*n;
                if ( (eta - bin_width*.5 < -.5) && (eta + bin_width*.5 > -1) ) // -1 < eta < -.5
                    refmult2 += get_frac(-1.,-.5,eta,bin_width)*n;

                if ( (eta - bin_width*.5 < .5) && (eta + bin_width*.5 > -.5) ) // -.5 < eta < .5
                    refmult1 += get_frac(-.5,.5,eta,bin_width)*n;

                if ( (eta - bin_width*.5 < 1.) && (eta + bin_width*.5 > .5) ) // .5 < eta < 1.
                    refmult2 += get_frac(.5,1.,eta,bin_width)*n;
            }

            if ( (eta - bin_width*.5 < -2.1) && ((eta + bin_width*.5 > -5.1)) ){ // -5.1 < eta < -2.1
                FwdAll += get_frac(-5.1,-2.1,eta,bin_width)*n;
                if ( (eta - bin_width*.5 < -2.1) && ((eta + bin_width*.5 > -3.0)) ) // -3. < eta < -2.1
                    Fwd1 += get_frac(-3.,-2.1,eta,bin_width)*n;
                if ( (eta - bin_width*.5 < -3.0) && ((eta + bin_width*.5 > -4.0)) ) // -4.0 < eta < -3
                    Fwd2 += get_frac(-4.,-3.,eta,bin_width)*n;
                if ( (eta - bin_width*.5 < -4.0) && ((eta + bin_width*.5 > -5.0)) ) // -5.0 < eta < -4
                    Fwd3 += get_frac(-5.,-4.,eta,bin_width)*n;
            }

            if ( (eta - bin_width*.5 < 5.1) && ((eta + bin_width*.5 > 2.1)) ){ // 2.1 < eta < 5.1
                FwdAll += get_frac(2.1,5.1,eta,bin_width)*n;
                if ( (eta - bin_width*.5 < 3.) && ((eta + bin_width*.5 > 2.1)) ) // 2.1 < eta < 3
                    Fwd1 += get_frac(2.1,3.,eta,bin_width)*n;
                if ( (eta - bin_width*.5 < 4) && ((eta + bin_width*.5 > 3)) ) // 3 < eta < 4
                    Fwd2 += get_frac(3.,4.,eta,bin_width)*n;
                if ( (eta - bin_width*.5 < 5) && ((eta + bin_width*.5 > 4)) ) // -5.0 < eta < -4
                    Fwd3 += get_frac(4.,5.,eta,bin_width)*n;
            }

        } else break;
    }
}

double AMPTSmearer::get_frac(double lower_bound, double upper_bound, double central_eta, double bin_width){
    double eta_frac = std::min( bin_width*.5+central_eta - lower_bound, upper_bound - (central_eta - bin_width*.5))/bin_width;
    double bin_frac = std::min((upper_bound-lower_bound)/bin_width,1.);
    return std::min(eta_frac,bin_frac);
}

AMPTSmearer::~AMPTSmearer(){}

///\brief Creates the history of collision of each parton
void AMPTSmearer::parse_history(){
    // Read the AMPT text output (deferred from the constructor so the SMASH path
    // can construct the smearer without these files being present).
    cols_hist   = readlines(results_path+"/parton-collisionsHistory.dat");
    init_parton = readlines(results_path+"/parton-initial-afterPropagation.dat");

    std::string ampt_header = readlines(results_path+"/ampt.dat")[0];
    impact_parameter = atof(split(ampt_header,' ')[3].data());
    NpartTarg = atof(split(ampt_header,' ')[4].data());
    NpartProj = atof(split(ampt_header,' ')[5].data());
    Npart = NpartTarg+NpartProj;
    NpartTargElastic = atof(split(ampt_header,' ')[6].data());
    NpartProjElastic = atof(split(ampt_header,' ')[8].data());

    std::cout << "b = " << impact_parameter << " fm, Npart = " << Npart
              << " (targ " << NpartTarg << ", proj " << NpartProj << ")" << std::endl;

    npartons = atoi( split(cols_hist[0],' ')[1].data() );
    ncollisions = (cols_hist.size()-1)/5;

    //Parse information of hadron creation
    this->parton_histories.resize(npartons);
    for (int iparton = 1; iparton< npartons+1; ++iparton){
        svec parton_info = split(init_parton[iparton],' ');

        Vec3 pos({(double) atof(parton_info[5].data()),
                  (double) atof(parton_info[6].data()),
                  (double) atof(parton_info[7].data())});
        Vec3 mom({(double) atof(parton_info[1].data()),
                  (double) atof(parton_info[2].data()),
                  (double) atof(parton_info[3].data())});

        this->parton_histories[iparton-1].push_back(PartonCollision(
            atoi(parton_info[0].data()), //pid
            atof(parton_info[4].data()), //mass
            atof(parton_info[8].data()), //time
            pos, Vec3({0, 0, 0}), mom ));
            net_p[1] += mom[0];
            net_p[2] += mom[1];
            net_p[3] += mom[2];
            net_p[0] += sqrt(pow(atof(parton_info[4].data()),2)
                        + pow(mom[0],2) + pow(mom[1],2) + pow(mom[2],2) );
    }

    for (int icol = 0; icol < ncollisions; ++icol){
        int start_entry = icol*5+1;
        svec col_header = split(cols_hist[start_entry],' ');
        int iparton1 = atoi(col_header[3].data())-1;
        int iparton2 = atoi(col_header[4].data())-1;

        //Insert parton 1
        svec p1_incoming = split(cols_hist[start_entry+1], ' ');
        svec p2_incoming = split(cols_hist[start_entry+2], ' ');
        svec p1_outgoing = split(cols_hist[start_entry+3], ' ');
        svec p2_outgoing = split(cols_hist[start_entry+4], ' ');

        double col_t = atof(p1_outgoing[8].data());
        Vec3 col_pos({(double) atof(p1_outgoing[5].data()),
                      (double) atof(p1_outgoing[6].data()),
                      (double) atof(p1_outgoing[7].data())});

        Vec3 mom_in1({(double) atof(p1_incoming[1].data()),
                      (double) atof(p1_incoming[2].data()),
                      (double) atof(p1_incoming[3].data())});
        Vec3 mom_in2({(double) atof(p2_incoming[1].data()),
                      (double) atof(p2_incoming[2].data()),
                      (double) atof(p2_incoming[3].data())});

        Vec3 mom_out1({(double) atof(p1_outgoing[1].data()),
                       (double) atof(p1_outgoing[2].data()),
                       (double) atof(p1_outgoing[3].data())});
        Vec3 mom_out2({(double) atof(p2_outgoing[1].data()),
                       (double) atof(p2_outgoing[2].data()),
                       (double) atof(p2_outgoing[3].data())});


        this->parton_histories[iparton1].push_back(PartonCollision(
            atoi(p1_outgoing[0].data()),
            atof(p1_outgoing[4].data()),
            col_t, col_pos,
            mom_in1, mom_out1));

        this->parton_histories[iparton2].push_back(PartonCollision(
            atoi(p2_outgoing[0].data()),
            atof(p2_outgoing[4].data()),
            col_t, col_pos,
            mom_in2, mom_out2));

        //Check energy and momentum conservation
        if (fabs(mom_in1[0]+mom_in2[0]-mom_out1[0]-mom_out2[0]) > 1.1E-3)
            std::cout << "Net px = " << mom_in1[0]+mom_in2[0]-mom_out1[0]-mom_out2[0] << std::endl;
        if (fabs(mom_in1[1]+mom_in2[1]-mom_out1[1]-mom_out2[1]) > 1.1E-3)
            std::cout << "Net py = " << mom_in1[1]+mom_in2[1]-mom_out1[1]-mom_out2[1] << std::endl;
        if (fabs(mom_in1[2]+mom_in2[2]-mom_out1[2]-mom_out2[2]) > 1.1E-3)
            std::cout << "Net pz = " << mom_in1[2]+mom_in2[2]-mom_out1[2]-mom_out2[2] << std::endl;


    }

}

///\brief Read a SMASH OSCAR2013 "SMASH_IC" particle list into parton_histories.
///
/// The OSCAR2013Extended SMASH_IC format lists one hadron per data line at its
/// formation/last-interaction point:
///   t x y z mass p0 px py pz pdg ID charge ncoll form_time xsecfac
///   proc_id_origin proc_type_origin time_last_coll pdg_mother1 pdg_mother2
///   baryon_number strangeness
/// (units: fm for t,x,y,z; GeV for mass,p0..pz). Lines beginning with '#' are
/// comments / event markers.
///
/// Each hadron becomes a
/// single-entry collision history (its formation point). Everything downstream
/// is IDENTICAL to AMPT — propagate(tau0) free-streams every particle formed
/// before tau0 up to tau0 (and drops those formed later), then fill_Tmunu
/// deposits it. The only SMASH-specific bit is that the conserved charges
/// (B, Q, S) are taken verbatim from the file instead of derived from the PID.
void AMPTSmearer::parse_smash(){
    // Accept either a direct ".oscar" file or a directory holding SMASH_IC.oscar.
    std::string path = results_path;
    if (path.size() < 6 || path.substr(path.size()-6) != ".oscar")
        path += "/SMASH_IC.oscar";

    svec lines = readlines(path);
    if (lines.empty()){
        std::cerr << "[ERROR]: SMASH IC file empty or not found: " << path << std::endl;
        exit(1);
    }

    // Charges live in the last three named columns; require a full row.
    const size_t NCOL = 22;
    int nskip = 0;
    for (const std::string& line : lines){
        if (line.empty() || line[0] == '#') continue;
        svec f = split(line, ' ');
        if (f.size() < NCOL){ ++nskip; continue; }

        double t    = atof(f[0].data());
        Vec3   pos({ atof(f[1].data()), atof(f[2].data()), atof(f[3].data()) });
        double mass = atof(f[4].data());
        Vec3   mom({ atof(f[6].data()), atof(f[7].data()), atof(f[8].data()) });
        int    pid  = atoi(f[9].data());

        // One-entry collision history = the formation point. incoming_mom is
        // unused by free_streamer; outgoing_mom carries the velocity.
        PartonCollision pc(pid, mass, t, pos, Vec3({0,0,0}), mom);
        pc.use_stored_charges = true;
        pc.echarge = atof(f[11].data());  // electric charge
        pc.bcharge = atof(f[20].data());  // baryon number
        pc.scharge = atof(f[21].data());  // strangeness
        parton_histories.push_back(std::vector<PartonCollision>{pc});

        net_p[1] += mom[0];
        net_p[2] += mom[1];
        net_p[3] += mom[2];
        net_p[0] += sqrt(mass*mass + mom[0]*mom[0] + mom[1]*mom[1] + mom[2]*mom[2]);
    }
    npartons = (int)parton_histories.size();

    // These AMPT-event quantities are undefined for a SMASH IC; keep them at 0.
    impact_parameter = 0.; Npart = 0.; NpartTarg = 0.; NpartProj = 0.;

    std::cout << "[INFO]: Read " << npartons << " SMASH IC particles from " << path
              << " (skipped " << nskip << " malformed lines)" << std::endl;
}

///\brief propagates a parton to a time tau_f
///\param last_collision the last collision of the particle
///\param tau_f the time to which we desire to free-stream the particle
///\return A ThermalizedParton object
PartonThermalized AMPTSmearer::free_streamer(PartonCollision last_collision, double tau_f) const
{

    ///Aliases for particle properties (position, momenttum etc)
    Vec3 x0 = last_collision.x;
    double t0 = last_collision.t;
    Vec3 p = last_collision.outgoing_mom;
    double mass = last_collision.mass;

    //Compute Lorentz factor
    double gamma_v = sqrt(p[0]*p[0] + p[1]*p[1] + p[2]*p[2])/mass;
    double gamma_l = sqrt(1.+gamma_v*gamma_v); // = 1/sqrt(1-v^2)

    Vec3 vel = Vec3({(double) (p[0]/mass/gamma_l),
                     (double) (p[1]/mass/gamma_l),
                     (double) (p[2]/mass/gamma_l)});

    // For ultra-relativistic partons (gamma >> 1) floating-point cancellation in
    // 1 - v^2 makes |v|^2 tip just above 1. Clamp |v|^2 < 1 before the check.
    // A tiny excess is harmless rounding; a large one means the input mass or
    // momentum is bad, so warn only for those really problematic cases.
    double v2 = pow(vel[0],2) + pow(vel[1],2) + pow(vel[2],2);
    if (v2 >= 1.0) {
        const double v2_tol = 1e-9;  // rounding noise stays below this
        if (v2 - 1.0 > v2_tol)
            std::cout << "[WARN] superluminal parton (pid " << last_collision.pid
                      << ", mass " << mass << "): |v|^2 = " << v2
                      << " clamped to 1" << std::endl;
        double inv_v = 1.0/sqrt(v2);
        vel[0] *= inv_v; vel[1] *= inv_v; vel[2] *= inv_v;
        v2 = 1.0 - 1e-15;
    }

    double pos_init = x0[2]-vel[2]*t0; //Position back-propagated to t = 0 - This quantity appears many times
                                      //thus we compute it only once here
                     
    //Final time in cartesian coordinates
    double tf;
    if(coordinates == "cartesian"){
        tf = tau_f;
        //tf = tf/(1-pow(vel[2],2));
    }
    else if(coordinates == "hyperbolic"){
        double Delta= sqrt( pow(pos_init,2) + pow(tau_f,2)*(1-pow(vel[2],2)) );
        tf = vel[2]*pos_init + Delta;
        tf = tf/(1-pow(vel[2],2));
    }



    //Final position in cartesian coordinates

    Vec3 final_pos;
    final_pos[0] = x0[0] + vel[0]*(tf-t0);
    final_pos[1] = x0[1] + vel[1]*(tf-t0);
    final_pos[2] = x0[2] + vel[2]*(tf-t0);

    PartonThermalized result(last_collision.pid, mass, tf, final_pos, p);
    // Carry the conserved charges (set by the SMASH path) through free-streaming.
    result.use_stored_charges = last_collision.use_stored_charges;
    result.bcharge = last_collision.bcharge;
    result.echarge = last_collision.echarge;
    result.scharge = last_collision.scharge;
    return result;
}

///\brief Propagates all partons to
void AMPTSmearer::propagate(double tau_f){

    //ProgressBar pb(npartons, "Free-Streaming:");
    int nform_below = 0;
    int nform_0coll = 0;
    int nform_coll = 0;
    double max_eta_s = 0.;
    double max_x = 0.;
    double max_y = 0.;
    double max_z = 0.;
    double max_p = 0.;
    double max_px = 0.;
    double max_py = 0.;
    double max_pz = 0.;
    double t_form = 0.;
    std::cout << "Number of partons: " << npartons << std::endl;
    for (auto parton_cols : this->parton_histories){
        int ncols = parton_cols.size();
        double t_p1;
        double t_m1;
        if (coordinates == "cartesian"){
            t_m1 = parton_cols[0].t;
        } else if (coordinates == "hyperbolic"){
            t_m1 = parton_cols[0].tau;
        } else {
            std::cout << "Unknown coordinate system" << std::endl;
            exit(1);
        }
        double t_before;
        if(t_m1 < tau_f){
            nform_below++;
        }
        bool crossed = false;
        //formation time
        if ((ncols == 1) && (t_m1 <= tau_f)){
            nform_0coll++;
            thermalized_partons.push_back( free_streamer(parton_cols[ncols-1],tau_f) );
            crossed = true;
        } 
        else {

            for (int icol=0; icol<ncols-1; ++icol){
                if (coordinates == "cartesian"){
                    t_p1 = parton_cols[icol+1].t;
                    t_before = parton_cols[icol].t;
                    t_form = parton_cols[icol].t;
                } else if (coordinates == "hyperbolic"){
                    t_p1 = parton_cols[icol+1].tau;
                    t_before = parton_cols[icol].tau;
                    t_form = parton_cols[icol].tau;
                } else {
                    std::cout << "Unknown coordinate system" << std::endl;
                    exit(1);
                }
                if (t_p1 >= tau_f && t_form <= tau_f){
                //if (t_p1 >= tau_f ){
                    nform_coll++;
                    thermalized_partons.push_back( free_streamer(parton_cols[icol],tau_f) );
                    crossed = true;
                    break;
                }
            }
            double t_last = (coordinates=="hyperbolic") ? parton_cols[ncols-1].tau
                                                        : parton_cols[ncols-1].t;
            //if not crossed, free stream the last collision
            if (!crossed && t_last < tau_f && ncols != 1){
                // std::cout << "Last collision time: " << t_last << std::endl;  // silenced (per-parton)
                thermalized_partons.push_back( free_streamer(parton_cols[ncols-1],tau_f) );
                crossed = true;
            }
        }
        //#ifdef PROGRESSBAR
        //pb.step();
        //#endif
    }
    std::cout << "Number of partons below tau_f: " << nform_below << std::endl;
    std::cout << "Number of partons with 0 collisions accepted: " << nform_0coll << std::endl;
    std::cout << "Number of partons with collisions accepted: " << nform_coll << std::endl;
    std::cout << "tau_f = " << tau_f << std::endl;
    TFile* fdebug = new TFile("debug.root","recreate");
    TH1D* heta =  new TH1D("heta","Parton dN/deta",100,-5,5);
    TH1D* hY =  new TH1D("hY","Parton dN/deta",100,-5,5);
    TH1D* h_peta =  new TH1D("h_peta","Parton dN/deta",100,-5,5);
    TH2D* h_peta_vs_eta_s =  new TH2D("h_peta_vs_eta_s","Parton dN/deta",200,-10,10,200,-10,10);
    TH2D* hY_vs_eta_s =  new TH2D("hY_vs_eta_s","Parton dN/deta",200,-10,10,200,-10,10);
    TH2D* ht_vs_z =  new TH2D("hz_vs_t","Parton dN/deta",2000,-100,100,2000,0,100);
    std::cout << "Number of thermalized partons: " << thermalized_partons.size() << std::endl;
    for (auto p : thermalized_partons){

        double abs_p = sqrt( pow(p.outgoing_mom[0],2) + pow(p.outgoing_mom[1],2) + pow(p.outgoing_mom[2],2));
        double pz = p.outgoing_mom[2];
        double px = p.outgoing_mom[0];
        double py = p.outgoing_mom[1];
        double eta = 0.5*log( (abs_p + pz)/(abs_p - pz) );
        heta->Fill(eta);
        hY_vs_eta_s->Fill(p.eta_s,p.Y);
        hY->Fill(p.Y);
        h_peta->Fill(p.peta);
        h_peta_vs_eta_s->Fill(p.eta_s,p.peta);
        ht_vs_z->Fill(p.x[2],p.t);
        if ( fabs(p.eta_s) > max_eta_s )
            max_eta_s = fabs(p.eta_s);
        if ( fabs(p.x[0]) > max_x )
            max_x = fabs(p.x[0]);
        if ( fabs(p.x[1]) > max_y )
            max_y = fabs(p.x[1]);
        if ( fabs(abs_p) > max_p )
            max_p = fabs(abs_p);
        if ( fabs(p.x[2]) > max_z )
            max_z = fabs(p.x[2]);
        if ( fabs(px) > max_px )
            max_px = fabs(px);
        if ( fabs(py) > max_py )
            max_py = fabs(py);
        if ( fabs(pz) > max_pz )
            max_pz = fabs(pz);


    }
    fdebug->Write();
    fdebug->Close();
    // silenced per-event debug maxima:
    // std::cout<<"Largest x = "<< max_x << std::endl;  ... (x,y,eta_s,p,z,px,py,pz)
    return;

}

void AMPTSmearer::fill_Tmunu(double sr,double seta){

    double _sigma_r = sr;
    double _sigma_eta = seta;
    //Norm that will accompany the smearing
    double norm;
    //jacobian tau
    if (coordinates == "cartesian")
        norm = K/2./M_PI/pow(_sigma_r,2)/sqrt(2*M_PI)/_sigma_eta;
    else if (coordinates == "hyperbolic")
        norm = K/2./M_PI/pow(_sigma_r,2)/sqrt(2*M_PI)/_sigma_eta/tau0;
    
    
    double up = 0.;
    double dw =0.;
    double st = 0.;
    double ch =0.;
    double top =0.;
    double bot=0.;
    Tmunu.resize(boost::extents[nx][ny][neta]);
    rhob.resize(boost::extents[nx][ny][neta]);
    j0.resize(boost::extents[nx][ny][neta]);
    j1.resize(boost::extents[nx][ny][neta]);
    j2.resize(boost::extents[nx][ny][neta]);
    j3.resize(boost::extents[nx][ny][neta]);
    j0e.resize(boost::extents[nx][ny][neta]);
    j1e.resize(boost::extents[nx][ny][neta]);
    j2e.resize(boost::extents[nx][ny][neta]);
    j3e.resize(boost::extents[nx][ny][neta]);
    j0s.resize(boost::extents[nx][ny][neta]);
    j1s.resize(boost::extents[nx][ny][neta]);
    j2s.resize(boost::extents[nx][ny][neta]);
    j3s.resize(boost::extents[nx][ny][neta]);

    double dx = Lx/(nx-1);
    double dy = Ly/(ny-1);
    double deta = Leta/(neta-1);
    //Zero the arrays
    for (int ix=0; ix<nx; ++ix)
    for (int iy=0; iy<ny; ++iy)
    for (int ieta=0; ieta<neta; ++ieta){
        j0[ix][iy][ieta] = 0.;
        j1[ix][iy][ieta] = 0.;
        j2[ix][iy][ieta] = 0.;
        j3[ix][iy][ieta] = 0.;
        j0e[ix][iy][ieta] = 0.;
        j1e[ix][iy][ieta] = 0.;
        j2e[ix][iy][ieta] = 0.;
        j3e[ix][iy][ieta] = 0.;
        j0s[ix][iy][ieta] = 0.;
        j1s[ix][iy][ieta] = 0.;
        j2s[ix][iy][ieta] = 0.;
        j3s[ix][iy][ieta] = 0.;
        rhob[ix][iy][ieta] = 0.;
        for(int mu=0; mu<4; ++mu)
        for(int nu=mu; nu<4; ++nu)
            Tmunu[ix][iy][ieta][mu][nu] = 0.; 
    }

    auto smearing_func = [&_sigma_r, &_sigma_eta](Vec3 x0, Vec3 x){
        double arg = (pow(x[0]-x0[0],2) + pow(x[1]-x0[1],2))/pow(_sigma_r,2)
                     + pow(x[2] - x0[2],2)/pow(_sigma_eta,2);
        return std::exp(-arg*.5);
    };

    auto smearing_func_spline = [this](Vec3 x0, Vec3 x){
        double spline_norm_2d =  5./(14.*M_PI *this->sigma_r*this->sigma_r);
        double spline_norm_1d = 1./(6.*this->sigma_eta);
        double kernel_1d = 0.0;
        double kernel_2d = 0.0;
        double q_2d = sqrt(pow(x[0]-x0[0],2) + pow(x[1]-x0[1],2))/this->sigma_r;
        double q_1d =  sqrt(pow(x[2] - x0[2],2))/this->sigma_eta;

        if (q_2d <=1.){
            kernel_2d = spline_norm_2d*(std::pow(2.-q_2d,3) - 4.*pow(1.-q_2d,3));
        }
        else if( q_2d <= 2.){
            kernel_2d = spline_norm_2d*(std::pow(2.-q_2d,3));
        }
        else{
            kernel_2d = 0.0;
        }
        

        if (q_1d <= 1.){
            kernel_1d = spline_norm_1d*(std::pow(2.-q_1d,3) - 4.*pow(1.-q_1d,3));
        }
        else if( q_1d <= 2.){
            kernel_1d = spline_norm_1d*(std::pow(2.-q_1d,3));
        }
        else{
            kernel_1d = 0.0;
        }
        //jacobian tau
        double arg =  this->K*kernel_1d*kernel_2d;
        if (this->coordinates == "hyperbolic")
            arg = this->K*kernel_1d*kernel_2d/this->tau0;
        return arg;
    };

    //db to identify quarks
    TDatabasePDG db = TDatabasePDG();
    //#pragma omp parallel
    //#pragma omp for reduction(+:Tmunu)
    //ProgressBar pb(thermalized_partons.size(),"EMT building");
    //for (int iparton=0;iparton<npartons;++iparton){
    //    auto parton = thermalized_partons[iparton];
    for (auto parton : thermalized_partons){

        //Creates the tensor associated to the parton
        std::array<double,4> mom;

        //check for coordinate system
        if(coordinates == "cartesian"){
            mom[0] = std::sqrt(pow(parton.mass,2.) + pow(parton.outgoing_mom[0],2) + pow(parton.outgoing_mom[1],2) + pow(parton.outgoing_mom[2],2));
            mom[1] = parton.outgoing_mom[0];
            mom[2] = parton.outgoing_mom[1];
            mom[3] = parton.outgoing_mom[2];
        } else if (coordinates == "hyperbolic"){
            mom[0] = parton.ptau;
            mom[1] = parton.outgoing_mom[0];
            mom[2] = parton.outgoing_mom[1];
            mom[3] = parton.peta;
        } else {
            std::cout << "Unknown coordinate system" << std::endl;
        }

        Mat4x4 parton_Tmunu;
        for (int mu=0; mu<4; ++mu)
        for (int nu=mu; nu<4; ++nu)
            parton_Tmunu[mu][nu] = mom[mu]*mom[nu]/mom[0];


        //Conserved charges of the parton: Q = baryon number, Qe = electric, Qs = strangeness
        double Q = 0.;
        double Qe = 0.;
        double Qs = 0.;
        if (parton.use_stored_charges) {
            // SMASH hadrons carry B, Q, S explicitly in the OSCAR file.
            Q  = parton.bcharge;
            Qe = parton.echarge;
            Qs = parton.scharge;
        } else {
        TParticlePDG* particle = db.GetParticle(parton.pid);
        if (!std::string(particle->ParticleClass()).compare("Quark")) {
            if (parton.pid == 1) {          // up quark
                Q = 1. / 3.;
                Qe = 2. / 3. ;
                Qs = 0.0;
            } else if (parton.pid == -1) {  // anti-up quark
                Q = -1. / 3.;
                Qe = -2. / 3. ;
                Qs = 0.0;
            } else if (parton.pid == 2) {   // down quark
                Q = 1. / 3.;
                Qe = -1. / 3. ;
                Qs = 0.0;
            } else if (parton.pid == -2) {  // anti-down quark
                Q = -1. / 3.;
                Qe = 1. / 3. ;
                Qs = 0.0;
            } else if (parton.pid == 3) {   // strange quark
                Q = 1. / 3.;
                Qe = -1. / 3. ;
                Qs = -1.0;
            } else if (parton.pid == -3) {  // anti-strange quark
                Q = -1. / 3.;
                Qe = 1. / 3. ;
                Qs = 1.0;
            } 
            // sometimes one or two heavier quarks appear in the initial condition
            else if (parton.pid == 4) {   // charm quark
                Q = 1. / 3.;
                Qe = 2. / 3. ;
                Qs = 0.0;
            } else if (parton.pid == -4) {  // anti-charm quark
                Q = -1. / 3.;
                Qe = -2. / 3. ;
                Qs = 0.0;
            } else if (parton.pid == 5) {   // bottom quark
                Q = 1. / 3.;
                Qe = -1. / 3. ;
                Qs = 0.0;
            } else if (parton.pid == -5) {  // anti-bottom quark
                Q = -1. / 3.;
                Qe = 1. / 3. ;
                Qs = 0.0;
            } else if (parton.pid == 6) {   // top quark
                Q = 1. / 3.;
                Qe = 2. / 3. ;
                Qs = 0.0;
            } else if (parton.pid == -6) {  // anti-top quark
                Q = -1. / 3.;
                Qe = -2. / 3. ;
                Qs = 0.0;
            }

            // Add more else-if clauses for other quarks and antiquarks as needed
        }
        }  // end AMPT quark-charge branch (use_stored_charges == false)



        if(abs(parton.pid) == 1){dw +=1.;};
        if(abs(parton.pid) == 2){up +=1.;};
        if(abs(parton.pid) == 3){st +=1.;};
        if(abs(parton.pid) == 4){ch +=1.;};
        if(abs(parton.pid) == 5){bot +=1.;};
        if(abs(parton.pid) == 6){top +=1.;};

        double x0 = parton.x[0];
        double y0 = parton.x[1];
        double eta0;
        if(coordinates == "hyperbolic"){
            eta0 = parton.eta_s;
        } else if (coordinates == "cartesian"){
            eta0 = parton.x[2];
        }
        Vec3 pos0({x0,y0,eta0});

        //Constraint search for a cube in a range 8*sigma_r and 8*sigma_eta
        double min_ix = std::max<double>(floor( (x0-rxy*2.*_sigma_r+Lx/2)/dx ), .0);
        double min_iy = std::max<double>(floor( (y0-rxy*2.*_sigma_r+Ly/2)/dy ), .0);
        double min_ieta = std::max<double>(floor( (eta0-reta*2.*_sigma_eta+Leta/2)/deta ), .0);

        double max_ix = std::min<double>(ceil( (x0+rxy*2.*_sigma_r+Lx/2)/dx ), nx-1);
        double max_iy = std::min<double>(ceil( (y0+rxy*2.*_sigma_r+Ly/2)/dy ), ny-1);
        double max_ieta = std::min<double>(ceil( (eta0+reta*2.*_sigma_eta+Leta/2)/deta), neta-1);

        for (int ix=min_ix; ix<=max_ix; ++ix){
            double x = ix*dx - Lx*.5;
            for (int iy=min_iy; iy<=max_iy; ++iy){
                double y = iy*dy - Ly*.5;
                for (int ieta=min_ieta; ieta<=max_ieta; ++ieta){
                    double eta = ieta*deta - Leta*.5;
                    Vec3 pos({x, y, eta});
                    double smearing_factor_spline = smearing_func_spline(pos0,pos)/K; //
                    double smearing_factor_gaussian  = norm*smearing_func(pos0,pos)/K;
                    j0[ix][iy][ieta] += 1.*Q*smearing_factor_spline*mom[0]/(mom[0]);
                    j1[ix][iy][ieta] += 1.*Q*smearing_factor_spline*mom[1]/(mom[0]);
                    j2[ix][iy][ieta] += 1.*Q*smearing_factor_spline*mom[2]/(mom[0]);
                    j3[ix][iy][ieta] += 1.*Q*smearing_factor_spline*mom[3]/(mom[0]);
                    rhob[ix][iy][ieta] += 1.*Q*smearing_factor_spline;
                    j0e[ix][iy][ieta] += 1.*Qe*smearing_factor_spline*mom[0]/(mom[0]);
                    j1e[ix][iy][ieta] += 1.*Qe*smearing_factor_spline*mom[1]/(mom[0]);
                    j2e[ix][iy][ieta] += 1.*Qe*smearing_factor_spline*mom[2]/(mom[0]);
                    j3e[ix][iy][ieta] += 1.*Qe*smearing_factor_spline*mom[3]/(mom[0]);
                    j0s[ix][iy][ieta] += 1.*Qs*smearing_factor_spline*mom[0]/(mom[0]);
                    j1s[ix][iy][ieta] += 1.*Qs*smearing_factor_spline*mom[1]/(mom[0]);
                    j2s[ix][iy][ieta] += 1.*Qs*smearing_factor_spline*mom[2]/(mom[0]);
                    j3s[ix][iy][ieta] += 1.*Qs*smearing_factor_spline*mom[3]/(mom[0]);
                    for(int mu=0; mu<4; ++mu)
                    for(int nu=mu; nu<4; ++nu)
                        Tmunu[ix][iy][ieta][mu][nu] +=  K*parton_Tmunu[mu][nu]*smearing_factor_spline;
                }
            }
        }
        //#ifdef PROGRESSBAR
        //pb.step();
        //#endif
    }
    std::cout << "N_up: " << up << std::endl;
    std::cout << "N_down: "<< dw << std::endl;
    std::cout << "N_strange: "<< st << std::endl;
    std::cout << "N_charm: "<< ch << std::endl;
    std::cout << "N_bottom: "<< bot <<std::endl;
    std::cout << "N_top: "<< top <<std::endl;

    //Symmetrizes the tensor
    for (int ix=0; ix<nx; ++ix)
    for (int iy=0; iy<ny; ++iy)
    for (int ieta=0; ieta<neta; ++ieta)
    for(int mu=0; mu<4; ++mu)
    for(int nu=mu; nu<4; ++nu)
        Tmunu[ix][iy][ieta][nu][mu] = Tmunu[ix][iy][ieta][mu][nu];

}
