
#include <string>
#include <iostream>
#include <fstream>



//#define PROGRESSBAR

#include "./AMPT_smearer.cpp"
#include "./hydrodynamizer.cpp"

//#include "/storage/home/kpala/JETSCAPE/external_packages/music/src/eos.h"

// #include "/storage/home/kpala/JETSCAPE/external_packages/music/EOS"
class AMPTGenesis
{
    private:
    int dada;
    public:
        AMPTGenesis();
        ~AMPTGenesis();
        // Smearing time
        double tau_LandauMatch;

        std::string input_folder_path;
        std::string output_file_path;

        std::vector<double> final_energy_density;
        std::vector<double> final_pressure;      
        std::vector<double> final_ut;                        
        std::vector<double> final_ux;             
        std::vector<double> final_uy;            
        std::vector<double> final_un;               
        std::vector<double> final_pitt;               
        std::vector<double> final_pitx;           
        std::vector<double> final_pity;               
        std::vector<double> final_pitn;                 
        std::vector<double> final_pixx;                
        std::vector<double> final_pixy;                
        std::vector<double> final_pixn;                
        std::vector<double> final_piyy;                
        std::vector<double> final_piyn;                
        std::vector<double> final_pinn;                
        std::vector<double> final_Pi;              
        std::vector<double> final_rhob;       
        std::vector<double> final_q0;
        std::vector<double> final_q1;
        std::vector<double> final_q2;
        std::vector<double> final_q3;

        
        double smearing_k;
        int nx; int ny; int neta; double Lx;double Ly; double Leta; 
                        double sigma_r; double sigma_eta; double tau0; double rxy; double reta;

        void run_genesis();

        double global_conservation_law(Hydrodynamizer hydro);

        void write_vectors(Hydrodynamizer hydro, AMPTSmearer smearer);

        void output_to_vectors(std::vector<double>&, //e
                            std::vector<double>&, //p
                            std::vector<double>&, //ut
                            std::vector<double>&, //ux
                            std::vector<double>&, //uy
                            std::vector<double>&, //un
                            std::vector<double>&, //pitt
                            std::vector<double>&, //pitx
                            std::vector<double>&, //pity
                            std::vector<double>&, //pitn
                            std::vector<double>&, //pixx
                            std::vector<double>&, //pixy
                            std::vector<double>&, //pixn
                            std::vector<double>&, //piyy
                            std::vector<double>&, //piyn
                            std::vector<double>&, //pinn
                            std::vector<double>&, //Pi
                            double &, //tau
                            std::vector<double>&, //rho_b
                            std::vector<double>&, //q_0
                            std::vector<double>&, //q_1
                            std::vector<double>&, //q2
                            std::vector<double>&);//q3
        void output_to_file();
        void output_to_file_center();
};

AMPTGenesis::AMPTGenesis(){
    dada = 0;
}

AMPTGenesis::~AMPTGenesis(){
}


//use this function to return final hydro variables as vectors within JETSCAPE
void AMPTGenesis::output_to_vectors(std::vector<double> &energy_density_out,
                                        std::vector<double> &pressure_out,
                                        std::vector<double> &ut_out,
                                        std::vector<double> &ux_out,
                                        std::vector<double> &uy_out,
                                        std::vector<double> &un_out,
                                        std::vector<double> &pitt_out,
                                        std::vector<double> &pitx_out,
                                        std::vector<double> &pity_out,
                                        std::vector<double> &pitn_out,
                                        std::vector<double> &pixx_out,
                                        std::vector<double> &pixy_out,
                                        std::vector<double> &pixn_out,
                                        std::vector<double> &piyy_out,
                                        std::vector<double> &piyn_out,
                                        std::vector<double> &pinn_out,
                                        std::vector<double> &Pi_out,
                                        double &tau_hydro,
                                        std::vector<double> &rho_b_out,
                                        std::vector<double> &q0_out,
                                        std::vector<double> &q1_out,
                                        std::vector<double> &q2_out,
                                        std::vector<double> &q3_out) {
    energy_density_out = final_energy_density;
    pressure_out = final_pressure;
    ut_out = final_ut;
    ux_out = final_ux;
    uy_out = final_uy;
    un_out = final_un;
    pitt_out = final_pitt;
    pitx_out = final_pitx;
    pity_out = final_pity;
    pitn_out = final_pitn;
    pixx_out = final_pixx;
    pixy_out = final_pixy;
    pixn_out = final_pixn;
    piyy_out = final_piyy;
    piyn_out = final_piyn;
    pinn_out = final_pinn;
    Pi_out = final_Pi;
    tau_hydro = tau0;
    rho_b_out = final_rhob;
    q0_out = final_q0;
    q1_out = final_q1;
    q2_out = final_q2;
    q3_out = final_q3;
}


void AMPTGenesis::write_vectors(Hydrodynamizer hydro, AMPTSmearer smearer){

    double dx = Lx/(nx-1.);
    double dy = Ly/(ny-1.);
    double deta = Leta/(neta-1.);
    double Nbar = 0;

        for(int ix=0;ix<nx;++ix){
        double x = ix*dx - Lx*.5;
        for(int iy=0;iy<ny;++iy){
            double y = iy*dy - Ly*.5;
            for(int ieta=0;ieta<neta;++ieta){
                double eta = ieta*deta - Leta*.5;

                //std::ios::fmtflags bckp_flags = fout.flags();

                if (hydro.TmunuOut[ix][iy][ieta].eps < 0.1){
                final_energy_density.push_back(0.);
                final_pressure.push_back(0./3.); //ideal eos
                final_ut.push_back(1.);
                final_ux.push_back(0.);
                final_uy.push_back(0.);
                final_un.push_back(0.);
                final_pitt.push_back(0.); 
                final_pitx.push_back(0.);
                final_pity.push_back(0.);
                final_pitn.push_back(0.); 
                final_pixx.push_back(0.);
                final_pixy.push_back(0.);
                final_pixn.push_back(0.); 
                final_piyy.push_back(0.);
                final_piyn.push_back(0.);
                final_pinn.push_back(0.);
                final_Pi.push_back(0.);
                final_rhob.push_back(0.);
                final_q0.push_back(0.);
                final_q1.push_back(0.);
                final_q2.push_back(0.);
                final_q3.push_back(0.);
                    continue;
                }
                double u0 = sqrt(1.+pow(hydro.TmunuOut[ix][iy][ieta].u[0],2.)+pow(hydro.TmunuOut[ix][iy][ieta].u[1],2.)+pow(tau0*hydro.TmunuOut[ix][iy][ieta].u[2],2.));
                double rhob_ = smearer.j0[ix][iy][ieta]*u0 - smearer.j1[ix][iy][ieta]*hydro.TmunuOut[ix][iy][ieta].u[0]- smearer.j2[ix][iy][ieta]*hydro.TmunuOut[ix][iy][ieta].u[1]-tau0*tau0*smearer.j3[ix][iy][ieta]*hydro.TmunuOut[ix][iy][ieta].u[2];
                final_energy_density.push_back(hydro.TmunuOut[ix][iy][ieta].eps);
                final_pressure.push_back((hydro.TmunuOut[ix][iy][ieta].eps)/3.); //ideal eos
                final_ut.push_back(sqrt(1.+pow(hydro.TmunuOut[ix][iy][ieta].u[0],2.)+pow(hydro.TmunuOut[ix][iy][ieta].u[1],2.)+pow(tau0*hydro.TmunuOut[ix][iy][ieta].u[2],2.)));
                final_ux.push_back(hydro.TmunuOut[ix][iy][ieta].u[0]);
                final_uy.push_back(hydro.TmunuOut[ix][iy][ieta].u[1]);
                final_un.push_back(hydro.TmunuOut[ix][iy][ieta].u[2]);
                final_pitt.push_back(hydro.TmunuOut[ix][iy][ieta].pitautau); 
                final_pitx.push_back(hydro.TmunuOut[ix][iy][ieta].pitaux);
                final_pity.push_back(hydro.TmunuOut[ix][iy][ieta].pitauy);
                final_pitn.push_back(hydro.TmunuOut[ix][iy][ieta].pitaueta); 
                final_pixx.push_back(hydro.TmunuOut[ix][iy][ieta].pixx);
                final_pixy.push_back(hydro.TmunuOut[ix][iy][ieta].pixy);
                final_pixn.push_back(hydro.TmunuOut[ix][iy][ieta].pixeta); 
                final_piyy.push_back(hydro.TmunuOut[ix][iy][ieta].piyy);    
                final_piyn.push_back(hydro.TmunuOut[ix][iy][ieta].piyeta);
                final_pinn.push_back(hydro.TmunuOut[ix][iy][ieta].pietaeta);
                final_Pi.push_back(-hydro.TmunuOut[ix][iy][ieta].Tr/3.);
                
                Nbar += smearer.rhob[ix][iy][ieta]*dx*dy*tau0*tau0*deta;
                //final_rhob.push_back(smearer.rhob[ix][iy][ieta]);
                final_rhob.push_back(smearer.j0[ix][iy][ieta]*u0 - smearer.j1[ix][iy][ieta]*hydro.TmunuOut[ix][iy][ieta].u[0]- smearer.j2[ix][iy][ieta]*hydro.TmunuOut[ix][iy][ieta].u[1]-tau0*tau0*smearer.j3[ix][iy][ieta]*hydro.TmunuOut[ix][iy][ieta].u[2]);
                final_q0.push_back(smearer.j0[ix][iy][ieta]-rhob_*u0);
                final_q1.push_back(smearer.j1[ix][iy][ieta]-rhob_*hydro.TmunuOut[ix][iy][ieta].u[0]);
                final_q2.push_back(smearer.j2[ix][iy][ieta]-rhob_*hydro.TmunuOut[ix][iy][ieta].u[1]);
                final_q3.push_back(smearer.j3[ix][iy][ieta]-rhob_*hydro.TmunuOut[ix][iy][ieta].u[2]);
                //std::cout << smearer.j1[ix][iy][ieta]-rhob_*hydro.TmunuOut[ix][iy][ieta].u[0] <<std::endl;
                //std::cout <<"dif" <<rhob_-smearer.rhob[ix][iy][ieta]<<std::endl;

                //fout.flags(bckp_flags);
            }
        }
    }
    //std::cout << "Net Baryon Number" << Nbar << std::endl;
}

void AMPTGenesis::output_to_file(){
    std::string path_out = output_file_path;
    std::ofstream fout;
    fout.open(path_out, std::ofstream::out );
    double dx = Lx/(nx-1.);
    double dy = Ly/(ny-1.);
    double deta = Leta/(neta-1.);

    //fout << "# b = " << smearer.impact_parameter << std::endl;
    //fout << "# Npart = "<< smearer.Npart << std::endl;
    //fout << "# NpartTarg = "<< smearer.NpartTarg << std::endl;
    //fout << "# NpartProj = "<< smearer.NpartProj << std::endl;
    //fout << "# NpartTargElastic = "<< smearer.NpartTargElastic << std::endl;
    //fout << "# NpartProjElastic = "<< smearer.NpartProjElastic << std::endl;
    //fout << "# refmult1 = "<< smearer.refmult1 << std::endl;
    //fout << "# refmult2 = "<< smearer.refmult2 << std::endl;
    //fout << "# refmult3 = "<< smearer.refmult3 << std::endl;
    //fout << "# Fwd1 = "<< smearer.Fwd1 << std::endl;
    //fout << "# Fwd2 = "<< smearer.Fwd2 << std::endl;
    //fout << "# Fwd3 = "<< smearer.Fwd3 << std::endl;
    //fout << "# FwdAll = " << smearer.FwdAll << std::endl;
    fout << "# nx = " << nx << std::endl;
    fout << "# ny = " << ny << std::endl;
    fout << "# neta = " << neta << std::endl;
    fout << "# Lx = " << Lx << std::endl;
    fout << "# Ly = " << Ly << std::endl;
    fout << "# Leta = " << Leta << std::endl;
    fout << "#x y eta epsilon ux uy ueta trace pitautau pitaux pitauy pitaueta pixx pixy pixeta piyy piyeta pietaeta rhob qt qx qy qeta" <<std::endl;


        for(int ix=0;ix<nx;++ix){
        double x = ix*dx - Lx*.5;
        for(int iy=0;iy<ny;++iy){
            double y = iy*dy - Ly*.5;
            for(int ieta=0;ieta<neta;++ieta){
                double eta = ieta*deta - Leta*.5;
                const int idx = (ny*neta)*ix + neta*iy + ieta;
                //std::ios::fmtflags bckp_flags = fout.flags();
                
                fout << x << " " << y << " " << eta << " "
                //<< std::scientific << std::setprecision(16)
                << final_energy_density[idx] << " "
                << final_ux[idx] << " "
                << final_uy[idx] << " "
                << final_un[idx] << " "
                << -3.*final_Pi[idx] << " "
                << final_pitt[idx]<< " "
                << final_pitx[idx] << " "
                << final_pity[idx] << " "
                << final_pitn[idx] << " "
                <<final_pixx[idx] << " "
                << final_pixy[idx] << " "
                <<final_pixn[idx] << " "
                << final_piyy[idx] << " "
                << final_piyn[idx] << " "
                << final_pinn[idx] << " "
                << final_rhob[idx] << " "
                << final_q0[idx] << " "
                << final_q1[idx] << " "
                << final_q2[idx] << " "
                << final_q3[idx] << std::endl;

                //fout.flags(bckp_flags);
            }
        }

    }

    fout.flush();
    fout.close();

}


void AMPTGenesis::output_to_file_center(){
    std::string path_out = output_file_path + "_center";
    std::ofstream fout;
    fout.open(path_out, std::ofstream::out );
    double dx = Lx/(nx-1.);
    double dy = Ly/(ny-1.);
    double deta = Leta/(neta-1.);

    //fout << "# b = " << smearer.impact_parameter << std::endl;
    //fout << "# Npart = "<< smearer.Npart << std::endl;
    //fout << "# NpartTarg = "<< smearer.NpartTarg << std::endl;
    //fout << "# NpartProj = "<< smearer.NpartProj << std::endl;
    //fout << "# NpartTargElastic = "<< smearer.NpartTargElastic << std::endl;
    //fout << "# NpartProjElastic = "<< smearer.NpartProjElastic << std::endl;
    //fout << "# refmult1 = "<< smearer.refmult1 << std::endl;
    //fout << "# refmult2 = "<< smearer.refmult2 << std::endl;
    //fout << "# refmult3 = "<< smearer.refmult3 << std::endl;
    //fout << "# Fwd1 = "<< smearer.Fwd1 << std::endl;
    //fout << "# Fwd2 = "<< smearer.Fwd2 << std::endl;
    //fout << "# Fwd3 = "<< smearer.Fwd3 << std::endl;
    //fout << "# FwdAll = " << smearer.FwdAll << std::endl;
    fout << "# nx = " << nx << std::endl;
    fout << "# ny = " << ny << std::endl;
    fout << "# neta = " << neta << std::endl;
    fout << "# Lx = " << Lx << std::endl;
    fout << "# Ly = " << Ly << std::endl;
    fout << "# Leta = " << Leta << std::endl;
    fout << "#x y eta epsilon ux uy ueta trace pitautau pitaux pitauy pitaueta pixx pixy pixeta piyy piyeta pietaeta rhob" <<std::endl;


        for(int ix=0;ix<nx;++ix){
        double x = ix*dx - Lx*.5;
        for(int iy=0;iy<ny;++iy){
            double y = iy*dy - Ly*.5;
            for(int ieta=0;ieta<neta;++ieta){
                double eta = ieta*deta - Leta*.5;
                const int idx = (ny*neta)*ix + neta*iy + ieta;
                if(eta = 0){
                //std::ios::fmtflags bckp_flags = fout.flags();
                
                fout << x << " " << y << " " << eta << " "
                //<< std::scientific << std::setprecision(16)
                << final_energy_density[idx] << " "
                << final_ux[idx] << " "
                << final_uy[idx] << " "
                << final_un[idx] << " "
                << -3.*final_Pi[idx] << " "
                << final_pitt[idx]<< " "
                << final_pitx[idx] << " "
                << final_pity[idx] << " "
                << final_pitn[idx] << " "
                <<final_pixx[idx] << " "
                << final_pixy[idx] << " "
                <<final_pixn[idx] << " "
                << final_piyy[idx] << " "
                << final_piyn[idx] << " "
                << final_pinn[idx] << " "
                << final_rhob[idx] << std::endl;
                }
                //fout.flags(bckp_flags);
            }
        }

    }

    fout.flush();
    fout.close();

}



double AMPTGenesis::global_conservation_law(Hydrodynamizer hydro){

    double dx = Lx/(nx-1.);
    double dy = Ly/(ny-1.);
    double deta = Leta/(neta-1.);

    double Ttaut = 0;
    for(int ix=0;ix<nx;++ix)
        for(int iy=0;iy<ny;++iy)
            for(int ieta=0;ieta<neta;++ieta){
                double eta = ieta*deta - Leta*.5;
                double cosheta = cosh(eta); double sinheta = sinh(eta);

                double eps = hydro.TmunuOut[ix][iy][ieta].eps;
                double utau = sqrt(pow(hydro.TmunuOut[ix][iy][ieta].u[0],2)
                                  +pow(hydro.TmunuOut[ix][iy][ieta].u[1],2)
                                  +pow(tau0*hydro.TmunuOut[ix][iy][ieta].u[2],2)+1.);
                double ueta = hydro.TmunuOut[ix][iy][ieta].u[2];
                double Tr = hydro.TmunuOut[ix][iy][ieta].Tr;
                double pitautau = hydro.TmunuOut[ix][iy][ieta].pitautau;
                double pitaueta = hydro.TmunuOut[ix][iy][ieta].pitaueta;

                double Ttautau = (4.*eps-Tr)*utau*utau/3. - (eps-Tr)/3. + pitautau;
                double Ttaueta = (4.*eps-Tr)*utau*ueta/3. + pitaueta;

                if (eps < 0) std::cout<<"Negative eigenvalue found" << std::endl;
                if (utau < 0) std::cout<<"Negative utau found" << std::endl;
                //if (3*eps-Tr < 0) std::cout<<"Negative term found" << std::endl;
                //if (pitautau < 0) std::cout<<"Negative shear term found" << std::endl;


                Ttaut += (Ttautau*cosheta + tau0*Ttaueta*sinheta)*dx*dy*deta*tau0; //GeV
    }

    return Ttaut;
}

void AMPTGenesis::run_genesis(){

    std::cout<< "AMPT input folder:  "<< input_folder_path << std::endl;

    AMPTSmearer smearer(input_folder_path,smearing_k,nx,
        ny,neta,Lx,Ly,Leta,sigma_r,sigma_eta,tau0,rxy,reta);
    smearer.parse_history();
    std::cout << "[INFO]: History parsed" << std::endl << std::flush;
    smearer.propagate(tau0);
    std::cout << "[INFO]: Finished free-streaming step" << std::endl << std::flush;
    smearer.fill_Tmunu(sigma_r,sigma_eta);
    std::cout << "[INFO]: Finished filling tensor" << std::endl << std::flush;

    Hydrodynamizer hydro(smearer.Tmunu,
                        tau0,
                        nx,
                        ny,
                        neta);

    hydro.diagonalize();
    std::cout << "[INFO]: Finished diagonalization." << std::endl << std::flush;

    double Target_E = smearer.net_p[0];
    double Actual_E = global_conservation_law(hydro);
    std::cout << "Target total energy: " << Target_E <<std::endl;
    std::cout << "Actual total energy: " << Actual_E <<std::endl;



    write_vectors(hydro, smearer);
    output_to_file();
    //EOS eos(14);
    //std::cout << eos.get_pressure(0.2,0.1) << std::endl;
    //output_to_file_center();
    std::cout << "[INFO]: Finished output" << std::endl << std::flush;



}



//#endif