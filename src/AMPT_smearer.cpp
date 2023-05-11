
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
    std::vector<PartonThermalized> thermalized_partons;

public:
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

    AMPTSmearer(std::string results_path, double K,int nx_, int ny_, int neta_, double Lx_,double Ly_, double Leta_,
                                double sigma_r_, double sigma_eta_, double tau0_, double rxy_, double reta_);
    double K;
    ~AMPTSmearer();

    void parse_history();
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
                                double sigma_r_, double sigma_eta_, double tau0_, double rxy_, double reta_):
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
net_p({0, 0, 0, 0})




{
    nx = nx_; ny=ny_; neta=neta_; Lx=Lx_; Ly=Ly_; Leta=Leta_; 
                                sigma_r=sigma_r_; sigma_eta = sigma_eta_; tau0=tau0_; rxy=rxy_; reta = reta_;
    K = Kin;

    cols_hist = readlines(results_path+"/parton-collisionsHistory.dat");
    init_parton = readlines(results_path+"/parton-initial-afterPropagation.dat");

    std::string ampt_header = readlines(results_path+"/ampt.dat")[0];
    impact_parameter = atof(split(ampt_header,' ')[3].data());
    NpartTarg = atof(split(ampt_header,' ')[4].data());
    NpartProj = atof(split(ampt_header,' ')[5].data());
    Npart = NpartTarg+NpartProj;
    NpartTargElastic = atof(split(ampt_header,' ')[6].data());
    NpartProjElastic = atof(split(ampt_header,' ')[8].data());

    get_dnde(results_path);


    std::cout << "b..................: "<< impact_parameter << " fm"<<std::endl;
    std::cout << "Npart..............: "<< Npart << std::endl;
    std::cout << "NpartTarg..........: "<< NpartTarg << std::endl;
    std::cout << "NpartProj..........: "<< NpartProj << std::endl;
    std::cout << "NpartTargElastic...: "<< NpartTargElastic << std::endl;
    std::cout << "NpartProjElastic...: "<< NpartProjElastic << std::endl;
    std::cout << "refmult1...........: "<< refmult1 << std::endl;
    std::cout << "refmult2...........: "<< refmult2 << std::endl;
    std::cout << "refmult3...........: "<< refmult3 << std::endl;
    std::cout << "Fwd1...............: "<< Fwd1 << std::endl;
    std::cout << "Fwd2...............: "<< Fwd2 << std::endl;
    std::cout << "Fwd3...............: "<< Fwd3 << std::endl;
    std::cout << "FwdAll.............: "<< FwdAll << std::endl;
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


    //Check consistency on computting velocity
    if ( fabs(1./sqrt(1-pow(vel[0],2) - pow(vel[1],2) - pow(vel[2],2) ) - gamma_l )/gamma_l > 1.E-7){
        std::cout <<"Inconsistency in velocity: " << 1./sqrt(1-pow(vel[0],2) - pow(vel[1],2) - pow(vel[2],2) ) - gamma_l << std::endl;
        std::cout << "v = "<< pow(vel[0],2) + pow(vel[1],2) + pow(vel[2],2)
                  << ", gamma = "<< gamma_l << std::endl;
        std::cout << "Expected gamma = " << 1./sqrt(1-pow(vel[0],2) - pow(vel[1],2) - pow(vel[2],2) ) << std::endl;

    }

    double pos_init = x0[2]-vel[2]*t0; //Position back-propagated to t = 0 - This quantity appears many times
                                      //thus we compute it only once here

    //Final time in cartesian coordinates
    double Delta = sqrt( pow(pos_init,2) + pow(tau_f,2)*(1-pow(vel[2],2)) );
    double tf = vel[2]*pos_init + Delta;
    tf = tf/(1-pow(vel[2],2));

    //Final position in cartesian coordinates

    Vec3 final_pos;
    final_pos[0] = x0[0] + vel[0]*(tf-t0);
    final_pos[1] = x0[1] + vel[1]*(tf-t0);
    final_pos[2] = x0[2] + vel[2]*(tf-t0);

    return PartonThermalized(last_collision.pid, mass, tf, final_pos, p);
}

///\brief Propagates all partons to
void AMPTSmearer::propagate(double tau_f){

    //ProgressBar pb(npartons, "Free-Streaming:");
    double max_eta_s = 0.;
    double max_x = 0.;
    double max_y = 0.;
    double max_p = 0.;
    for (auto parton_cols : this->parton_histories){
        int ncols = parton_cols.size();

        if ((ncols == 1) || (parton_cols[ncols-1].tau < tau_f)){
            thermalized_partons.push_back( free_streamer(parton_cols[ncols-1],tau_f) );
        } else {
            for (int icol=0; icol<ncols-1; ++icol){
                if (parton_cols[icol+1].tau > tau_f){
                    thermalized_partons.push_back( free_streamer(parton_cols[icol],tau_f) );
                    break;
                }
            }
        }
        //#ifdef PROGRESSBAR
        //pb.step();
        //#endif
    }
    TFile* fdebug = new TFile("debug.root","recreate");
    TH1D* heta =  new TH1D("heta","Parton dN/deta",100,-5,5);
    TH1D* hY =  new TH1D("hY","Parton dN/deta",100,-5,5);
    TH1D* h_peta =  new TH1D("h_peta","Parton dN/deta",100,-5,5);
    TH2D* h_peta_vs_eta_s =  new TH2D("h_peta_vs_eta_s","Parton dN/deta",200,-10,10,200,-10,10);
    TH2D* hY_vs_eta_s =  new TH2D("hY_vs_eta_s","Parton dN/deta",200,-10,10,200,-10,10);
    TH2D* ht_vs_z =  new TH2D("hz_vs_t","Parton dN/deta",2000,-100,100,2000,0,100);

    for (auto p : thermalized_partons){

        double abs_p = sqrt( pow(p.outgoing_mom[0],2) + pow(p.outgoing_mom[1],2) + pow(p.outgoing_mom[2],2));
        double pz = p.outgoing_mom[2];
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


    }
    fdebug->Write();
    fdebug->Close();
    std::cout<<"Largest x = "<< max_x << std::endl;
    std::cout<<"Largest y = "<< max_y << std::endl;
    std::cout<<"Largest eta_s = "<< max_eta_s << std::endl;
    std::cout<<"Largest abs_p = "<< max_p << std::endl;
    return;

}

void AMPTSmearer::fill_Tmunu(double sr,double seta){

    double _sigma_r = sr;
    double _sigma_eta = seta;
    //Norm that will accompany the smearing
    double norm = K/2./M_PI/pow(_sigma_r,2)/sqrt(2*M_PI)/_sigma_eta/tau0;
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

        auto smearing_func_spline = [&_sigma_r, &_sigma_eta](Vec3 x0, Vec3 x){
        double spline_norm_2d =  15./(14.*_sigma_r*_sigma_r);
        double spline_norm_1d = 1./(6.*_sigma_eta);
        double kernel_1d = 0.0;
        double kernel_2d = 0.0;
        double q_2d = sqrt(pow(x[0]-x0[0],2) + pow(x[1]-x0[1],2))/_sigma_r;
        double q_1d =  sqrt(pow(x[2] - x0[2],2))/_sigma_eta;

        if (q_2d >= 2.){kernel_2d += 0.0;}
        if (q_2d<2. && q_2d >= 1. ){kernel_2d += spline_norm_2d*pow(2.-q_2d,3.);}
        if (q_2d>=0 && q_2d <1.){kernel_2d += spline_norm_2d*(pow(2.-q_2d,3.)-4.*pow(1.-q_2d,3.));}

        if (q_1d >= 2.){kernel_1d += 0.0;}
        if (q_1d<2. && q_1d >= 1. ){kernel_1d += spline_norm_1d*(2.-pow(2.-q_1d,3.));}
        if (q_1d>=0 && q_1d <1.){kernel_1d += spline_norm_1d*(pow(2.-q_1d,3.)-4.*pow(1.-q_1d,3.));}
        double arg = kernel_1d*kernel_2d;
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
        std::array<double,4> mom({parton.ptau, parton.outgoing_mom[0],
                                 parton.outgoing_mom[1],
                                 parton.peta});
        Mat4x4 parton_Tmunu;
        for (int mu=0; mu<4; ++mu)
        for (int nu=mu; nu<4; ++nu)
            parton_Tmunu[mu][nu] = mom[mu]*mom[nu]/parton.ptau;


        //Creates the baryon current associated to the parton
        TParticlePDG* particle = db.GetParticle(parton.pid);

        double Q = 0.;
        if(!std::string(particle->ParticleClass()).compare("Quark") )
            Q = parton.pid < 0 ? -1./3. : 1./3.;
        
        if(abs(parton.pid) == 1){dw +=1.;};
        if(abs(parton.pid) == 2){up +=1.;};
        if(abs(parton.pid) == 3){st +=1.;};
        if(abs(parton.pid) == 4){ch +=1.;};
        if(abs(parton.pid) == 5){bot +=1.;};
        if(abs(parton.pid) == 6){top +=1.;};

        double x0 = parton.x[0];
        double y0 = parton.x[1];
        double eta0 = parton.eta_s;
        Vec3 pos0({x0,y0,eta0});

        //Constraint search for a cube in a range 8*sigma_r and 8*sigma_eta
        double min_ix = std::max<double>(floor( (x0-rxy*100.*_sigma_r+Lx/2)/dx ), .0);
        double min_iy = std::max<double>(floor( (y0-rxy*100.*_sigma_r+Ly/2)/dy ), .0);
        double min_ieta = std::max<double>(floor( (eta0-reta*100.*_sigma_eta+Leta/2)/deta ), .0);

        double max_ix = std::min<double>(ceil( (x0+rxy*100.*_sigma_r+Lx/2)/dx ), nx-1);
        double max_iy = std::min<double>(ceil( (y0+rxy*100.*_sigma_r+Ly/2)/dy ), ny-1);
        double max_ieta = std::min<double>(ceil( (eta0+reta*100.*_sigma_eta+Leta/2)/deta), neta-1);

        for (int ix=min_ix; ix<=max_ix; ++ix){
            double x = ix*dx - Lx*.5;
            for (int iy=min_iy; iy<=max_iy; ++iy){
                double y = iy*dy - Ly*.5;
                for (int ieta=min_ieta; ieta<=max_ieta; ++ieta){
                    double eta = ieta*deta - Leta*.5;
                    Vec3 pos({x, y, eta});
                    double smearing_factor = smearing_func_spline(pos0,pos)/tau0; //if using gaussian multiply by norm and divide barionic sector by  k
                    double smearing_factor_gaussian  = norm*smearing_func(pos0,pos)/K;
                    j0[ix][iy][ieta] += 1.*Q*smearing_factor_gaussian*mom[0]/(parton.ptau);
                    j1[ix][iy][ieta] += 1.*Q*smearing_factor_gaussian*mom[1]/(parton.ptau);
                    j2[ix][iy][ieta] += 1.*Q*smearing_factor_gaussian*mom[2]/(parton.ptau);
                    j3[ix][iy][ieta] += 1.*Q*smearing_factor_gaussian*mom[3]/(parton.ptau);
                    rhob[ix][iy][ieta] += 1.*Q*smearing_factor_gaussian;
                    for(int mu=0; mu<4; ++mu)
                    for(int nu=mu; nu<4; ++nu)
                        Tmunu[ix][iy][ieta][mu][nu] +=  K*parton_Tmunu[mu][nu]*smearing_factor_gaussian;
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
