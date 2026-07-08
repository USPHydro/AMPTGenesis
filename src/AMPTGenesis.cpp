
#include <string>
#include <iostream>
#include <iomanip>
#include <array>
#include <cmath>

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
        std::string coordinate_system;
        std::string input_format;   ///< "ampt" (parton history) or "smash" (OSCAR2013 IC)

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
        std::vector<double> final_rhoe;
        std::vector<double> final_q0e;
        std::vector<double> final_q1e;
        std::vector<double> final_q2e;
        std::vector<double> final_q3e;
        std::vector<double> final_rhos;
        std::vector<double> final_q0s;
        std::vector<double> final_q1s;
        std::vector<double> final_q2s;
        std::vector<double> final_q3s;


        
        double smearing_k;
        bool output_diffusion;
        bool output_raw_tmunu;
        double energy_density_cutoff;
        int nx; int ny; int neta; double Lx;double Ly; double Leta;
                        double sigma_r; double sigma_eta; double tau0; double rxy; double reta;
        // back-propagate late-forming partons (formation tau > tau0) onto the tau0 surface
        bool backpropagate = false;

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
                            std::vector<double>&, //q3
                            std::vector<double>&, //rhoe
                            std::vector<double>&, //q0e
                            std::vector<double>&, //q1e
                            std::vector<double>&, //q2e
                            std::vector<double>&, //q3e
                            std::vector<double>&, //rhos
                            std::vector<double>&, //q0s
                            std::vector<double>&, //q1s
                            std::vector<double>&, //q2s
                            std::vector<double>&); //q3s
        void output_to_file();
        // Write the raw (un-diagonalized) contravariant T^{mu nu} components plus the
        // lab-frame conserved currents, skipping the Landau matching done by output_to_file.
        void output_tmunu_to_file(AMPTSmearer& smearer);
        void output_to_file_center();
        // Fast eccentricity-only sigma scan: read partons + free-stream ONCE, then
        // loop fill_Tmunu+diagonalize over the sigma list, writing only eccentricities.
        void run_ecc_scan(const std::vector<double>& sigmas, const std::string& ecc_out);
        // GRID-FREE analytic eccentricity scan: sum over partons, kernel handled
        // analytically (numerator sigma-independent, denominator from kernel radial moments).
        void run_ecc_analytic(const std::vector<double>& sigmas, const std::string& ecc_out);
        // Compare (at eta_s=0 slice) spatial eps2/eps3 weighted by T^tt vs Landau e,
        // and momentum anisotropy eps_p from the ideal vs the full T^ij.
        void run_ecc_momentum(const std::string& ecc_out);
};

AMPTGenesis::AMPTGenesis(){
    dada = 0;
    output_diffusion = false;
    output_raw_tmunu = false;
    energy_density_cutoff = 0.15;
    input_format = "ampt";
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
                                        std::vector<double> &q3_out,
                                        std::vector<double> &rhoe_out,
                                        std::vector<double> &q0e_out,
                                        std::vector<double> &q1e_out,
                                        std::vector<double> &q2e_out,
                                        std::vector<double> &q3e_out,
                                        std::vector<double> &rhos_out,
                                        std::vector<double> &q0s_out,
                                        std::vector<double> &q1s_out,
                                        std::vector<double> &q2s_out,
                                        std::vector<double> &q3s_out) {
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
    rhoe_out = final_rhoe;
    q0e_out = final_q0e;
    q1e_out = final_q1e;
    q2e_out = final_q2e;
    q3e_out = final_q3e;
    rhos_out = final_rhos;
    q0s_out = final_q0s;
    q1s_out = final_q1s;
    q2s_out = final_q2s;
    q3s_out = final_q3s;
}


void AMPTGenesis::write_vectors(Hydrodynamizer hydro, AMPTSmearer smearer){

    double dx = Lx/(nx-1.);
    double dy = Ly/(ny-1.);
    double sqrt_mg;
    if (coordinate_system == "cartesian"){
        sqrt_mg = 1.;
    } else if (coordinate_system == "hyperbolic"){
        sqrt_mg = tau0;
    }
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
                final_rhoe.push_back(0.);
                final_q0e.push_back(0.);
                final_q1e.push_back(0.);
                final_q2e.push_back(0.);
                final_q3e.push_back(0.);
                final_rhos.push_back(0.);
                final_q0s.push_back(0.);
                final_q1s.push_back(0.);
                final_q2s.push_back(0.);
                final_q3s.push_back(0.);
                    continue;
                }
                double u0 = sqrt(1.+pow(hydro.TmunuOut[ix][iy][ieta].u[0],2.)+pow(hydro.TmunuOut[ix][iy][ieta].u[1],2.)+pow(sqrt_mg*hydro.TmunuOut[ix][iy][ieta].u[2],2.));
                double rhob_ = smearer.j0[ix][iy][ieta]*u0 - smearer.j1[ix][iy][ieta]*hydro.TmunuOut[ix][iy][ieta].u[0]- smearer.j2[ix][iy][ieta]*hydro.TmunuOut[ix][iy][ieta].u[1]-sqrt_mg*sqrt_mg*smearer.j3[ix][iy][ieta]*hydro.TmunuOut[ix][iy][ieta].u[2];
                double rhoe_ = smearer.j0e[ix][iy][ieta]*u0-smearer.j1e[ix][iy][ieta]*hydro.TmunuOut[ix][iy][ieta].u[0]-smearer.j2e[ix][iy][ieta]*hydro.TmunuOut[ix][iy][ieta].u[1]-sqrt_mg*sqrt_mg*smearer.j3e[ix][iy][ieta]*hydro.TmunuOut[ix][iy][ieta].u[2];
                double rhos_ = smearer.j0s[ix][iy][ieta]*u0-smearer.j1s[ix][iy][ieta]*hydro.TmunuOut[ix][iy][ieta].u[0]-smearer.j2s[ix][iy][ieta]*hydro.TmunuOut[ix][iy][ieta].u[1]-sqrt_mg*sqrt_mg*smearer.j3s[ix][iy][ieta]*hydro.TmunuOut[ix][iy][ieta].u[2];
                final_energy_density.push_back(hydro.TmunuOut[ix][iy][ieta].eps);
                final_pressure.push_back((hydro.TmunuOut[ix][iy][ieta].eps)/3.); //ideal eos
                final_ut.push_back(sqrt(1.+pow(hydro.TmunuOut[ix][iy][ieta].u[0],2.)+pow(hydro.TmunuOut[ix][iy][ieta].u[1],2.)+pow(sqrt_mg*hydro.TmunuOut[ix][iy][ieta].u[2],2.)));
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
                final_rhob.push_back(smearer.j0[ix][iy][ieta]*u0 - smearer.j1[ix][iy][ieta]*hydro.TmunuOut[ix][iy][ieta].u[0]- smearer.j2[ix][iy][ieta]*hydro.TmunuOut[ix][iy][ieta].u[1]-sqrt_mg*sqrt_mg*smearer.j3[ix][iy][ieta]*hydro.TmunuOut[ix][iy][ieta].u[2]);
                final_q0.push_back(smearer.j0[ix][iy][ieta]-rhob_*u0);
                final_q1.push_back(smearer.j1[ix][iy][ieta]-rhob_*hydro.TmunuOut[ix][iy][ieta].u[0]);
                final_q2.push_back(smearer.j2[ix][iy][ieta]-rhob_*hydro.TmunuOut[ix][iy][ieta].u[1]);
                final_q3.push_back(smearer.j3[ix][iy][ieta]-rhob_*hydro.TmunuOut[ix][iy][ieta].u[2]);
                final_rhoe.push_back(smearer.j0e[ix][iy][ieta]*u0-smearer.j1e[ix][iy][ieta]*hydro.TmunuOut[ix][iy][ieta].u[0]-smearer.j2e[ix][iy][ieta]*hydro.TmunuOut[ix][iy][ieta].u[1]-sqrt_mg*sqrt_mg*smearer.j3e[ix][iy][ieta]*hydro.TmunuOut[ix][iy][ieta].u[2]);
                final_q0e.push_back(smearer.j0e[ix][iy][ieta]-rhoe_*u0);
                final_q1e.push_back(smearer.j1e[ix][iy][ieta]-rhoe_*hydro.TmunuOut[ix][iy][ieta].u[0]);
                final_q2e.push_back(smearer.j2e[ix][iy][ieta]-rhoe_*hydro.TmunuOut[ix][iy][ieta].u[1]);
                final_q3e.push_back(smearer.j3e[ix][iy][ieta]-rhoe_*hydro.TmunuOut[ix][iy][ieta].u[2]);
                final_rhos.push_back(smearer.j0s[ix][iy][ieta]*u0-smearer.j1s[ix][iy][ieta]*hydro.TmunuOut[ix][iy][ieta].u[0]-smearer.j2s[ix][iy][ieta]*hydro.TmunuOut[ix][iy][ieta].u[1]-sqrt_mg*sqrt_mg*smearer.j3s[ix][iy][ieta]*hydro.TmunuOut[ix][iy][ieta].u[2]);
                final_q0s.push_back(smearer.j0s[ix][iy][ieta]-rhos_*u0);
                final_q1s.push_back(smearer.j1s[ix][iy][ieta]-rhos_*hydro.TmunuOut[ix][iy][ieta].u[0]);
                final_q2s.push_back(smearer.j2s[ix][iy][ieta]-rhos_*hydro.TmunuOut[ix][iy][ieta].u[1]);
                final_q3s.push_back(smearer.j3s[ix][iy][ieta]-rhos_*hydro.TmunuOut[ix][iy][ieta].u[2]);
                //std::cout << smearer.j1[ix][iy][ieta]-rhob_*hydro.TmunuOut[ix][iy][ieta].u[0] <<std::endl;
                //std::cout <<"dif" <<rhob_-smearer.rhob[ix][iy][ieta]<<std::endl;

                //fout.flags(bckp_flags);
            }
        }
    }
    //std::cout << "Net Baryon Number" << Nbar << std::endl;
}


void AMPTGenesis::output_to_file() {

    std::string path_out = output_file_path;
    std::ofstream fout(path_out);
    if (!fout) {
        std::cerr << "Error: Unable to open file " << path_out << std::endl;
        return;
    }

    double dx = Lx / (nx - 1.0);
    double dy = Ly / (ny - 1.0);
    double deta = Leta / (neta - 1.0);

    double xmin = -Lx / 2.0;
    double ymin = -Ly / 2.0;
    double etamin = -Leta / 2.0;

    // Counters for cells relative to the energy_density_cutoff
    uint64_t cells_above_energy_density_cutoff = 0;
    uint64_t cells_below_or_equal_energy_density_cutoff = 0;

    // We will buffer the body first, but we can already write a header.
    // The first header line keeps the original Python-script format.
    fout << "#0 " << dx << " " << dy << " " << deta
         << " 0 " << xmin << " " << ymin << " " << etamin << "\n";

    // --- initial-state eccentricities from the FULL energy-density grid ---
    //     Energy-density weighted, recentered, eps_n = |<r^n e^{i n phi}>| / <r^n>.
    //     Uses EVERY cell, including epsilon <= energy_density_cutoff (no threshold).
    {
        const int NMAX = 5;
        int ieta0 = (int)((0.0 - etamin) / deta + 0.5);   // grid index nearest eta_s = 0
        if (ieta0 < 0) ieta0 = 0;
        if (ieta0 >= neta) ieta0 = neta - 1;
        // pass 1: energy-weighted transverse centroid (eta_s-integrated and midrapidity)
        double Wi = 0, Wix = 0, Wiy = 0, Wm = 0, Wmx = 0, Wmy = 0;
        for (int ix = 0; ix < nx; ++ix) { double x = ix * dx + xmin;
            for (int iy = 0; iy < ny; ++iy) { double y = iy * dy + ymin;
                for (int ie = 0; ie < neta; ++ie) {
                    double e = final_energy_density[(ny * neta) * ix + neta * iy + ie];
                    if (e <= 0.0) continue;
                    Wi += e; Wix += e * x; Wiy += e * y;
                    if (ie == ieta0) { Wm += e; Wmx += e * x; Wmy += e * y; }
                } } }
        double x0i = (Wi > 0) ? Wix / Wi : 0.0, y0i = (Wi > 0) ? Wiy / Wi : 0.0;
        double x0m = (Wm > 0) ? Wmx / Wm : 0.0, y0m = (Wm > 0) ? Wmy / Wm : 0.0;
        // pass 2: r^n-weighted complex moments
        double ci[NMAX + 1] = {0}, si[NMAX + 1] = {0}, di[NMAX + 1] = {0};
        double cm[NMAX + 1] = {0}, sm[NMAX + 1] = {0}, dm[NMAX + 1] = {0};
        for (int ix = 0; ix < nx; ++ix) { double x = ix * dx + xmin;
            for (int iy = 0; iy < ny; ++iy) { double y = iy * dy + ymin;
                for (int ie = 0; ie < neta; ++ie) {
                    double e = final_energy_density[(ny * neta) * ix + neta * iy + ie];
                    if (e <= 0.0) continue;
                    { double xr = x - x0i, yr = y - y0i;
                      double r = std::sqrt(xr * xr + yr * yr), phi = std::atan2(yr, xr);
                      for (int n = 1; n <= NMAX; ++n) { double rn = std::pow(r, n);
                          di[n] += e * rn; ci[n] += e * rn * std::cos(n * phi); si[n] += e * rn * std::sin(n * phi); } }
                    if (ie == ieta0) { double xr = x - x0m, yr = y - y0m;
                      double r = std::sqrt(xr * xr + yr * yr), phi = std::atan2(yr, xr);
                      for (int n = 1; n <= NMAX; ++n) { double rn = std::pow(r, n);
                          dm[n] += e * rn; cm[n] += e * rn * std::cos(n * phi); sm[n] += e * rn * std::sin(n * phi); } }
                } } }
        fout << "# eccentricities (energy-density weighted, ALL cells incl. epsilon<=cutoff, recentered)\n";
        fout << "# eps_n_eta_integrated";
        for (int n = 2; n <= NMAX; ++n) { double v = (di[n] > 0) ? std::sqrt(ci[n]*ci[n] + si[n]*si[n]) / di[n] : 0.0; fout << " eps" << n << "=" << v; }
        fout << "\n# eps_n_midrapidity";
        for (int n = 2; n <= NMAX; ++n) { double v = (dm[n] > 0) ? std::sqrt(cm[n]*cm[n] + sm[n]*sm[n]) / dm[n] : 0.0; fout << " eps" << n << "=" << v; }
        fout << "\n";
        double e2i=(di[2]>0)?std::sqrt(ci[2]*ci[2]+si[2]*si[2])/di[2]:0.0, e3i=(di[3]>0)?std::sqrt(ci[3]*ci[3]+si[3]*si[3])/di[3]:0.0;
        double e2m=(dm[2]>0)?std::sqrt(cm[2]*cm[2]+sm[2]*sm[2])/dm[2]:0.0, e3m=(dm[3]>0)?std::sqrt(cm[3]*cm[3]+sm[3]*sm[3])/dm[3]:0.0;
        std::cout << "[AMPTGenesis] eccentricities (all cells): eps2_int=" << e2i << " eps3_int=" << e3i
                  << " | eps2_mid=" << e2m << " eps3_mid=" << e3m << std::endl;
    }

    std::ostringstream buffer;

    for (int ix = 0; ix < nx; ++ix) {
        double x = ix * dx + xmin;
        for (int iy = 0; iy < ny; ++iy) {
            double y = iy * dy + ymin;
            for (int ieta = 0; ieta < neta; ++ieta) {
                double eta = ieta * deta + etamin;
                int idx = (ny * neta) * ix + neta * iy + ieta;

                double epsilon = final_energy_density[idx];

                if (epsilon > energy_density_cutoff) {
                    ++cells_above_energy_density_cutoff;

                    double ux = final_ux[idx];
                    double uy = final_uy[idx];
                    double un = final_un[idx];
                    double trace = -3.0 * final_Pi[idx];

                    double pixx = final_pixx[idx];
                    double pixy = final_pixy[idx];
                    double pixn = final_pixn[idx];
                    double piyy = final_piyy[idx];
                    double piyn = final_piyn[idx];
                    double pinn = final_pinn[idx];

                    double rhob = final_rhob[idx];
                    double rhoe = final_rhoe[idx];
                    double rhos = final_rhos[idx];

                    buffer << std::fixed << std::setprecision(8)
                           << x << " " << y << " " << eta << " "
                           << epsilon << " "
                           << rhob << " " << rhos << " " << rhoe << " "
                           << ux << " " << uy << " " << un << " "
                           << trace << " "
                           << pixx << " " << pixy << " " << pixn << " "
                           << piyy << " " << piyn << " " << pinn;
                    if (output_diffusion) {
                        buffer << " "
                               << final_q0[idx]  << " " << final_q1[idx]  << " " << final_q2[idx]  << " " << final_q3[idx]  << " "
                               << final_q0s[idx] << " " << final_q1s[idx] << " " << final_q2s[idx] << " " << final_q3s[idx] << " "
                               << final_q0e[idx] << " " << final_q1e[idx] << " " << final_q2e[idx] << " " << final_q3e[idx];
                    }
                    buffer << "\n";
                } else {
                    ++cells_below_or_equal_energy_density_cutoff;
                }
            }
        }
    }

    // Write counting info as comment lines so downstream readers can ignore if needed.
    fout << "# cells_above_energy_density_cutoff " << energy_density_cutoff << " : " << cells_above_energy_density_cutoff << "\n";
    fout << "# cells_below_or_equal_energy_density_cutoff " << energy_density_cutoff << " : " << cells_below_or_equal_energy_density_cutoff << "\n";

    // Now dump the buffered field lines
    fout << buffer.str();
    fout.close();
}

// ---- raw (un-diagonalized) T^{mu nu} output -----------------------------------
// Writes the 10 independent contravariant components of the deposited stress-energy
// tensor T^{mu nu} (indices 0=tau,1=x,2=y,3=eta) directly, WITHOUT the Landau
// decomposition (eps, u, shear pi, bulk Pi) performed in output_to_file(). The
// lab-frame conserved currents j^mu for baryon (B), strangeness (S) and electric
// charge (Q) are appended so the downstream code can do its own frame matching.
// Cell selection uses T^{tau tau} against energy_density_cutoff (T^{tau tau} >= eps,
// so this is the natural, slightly looser proxy for the Landau energy density).
void AMPTGenesis::output_tmunu_to_file(AMPTSmearer& smearer) {

    std::ofstream fout(output_file_path);
    if (!fout) {
        std::cerr << "Error: Unable to open file " << output_file_path << std::endl;
        return;
    }

    double dx = Lx / (nx - 1.0);
    double dy = Ly / (ny - 1.0);
    double deta = Leta / (neta - 1.0);

    double xmin = -Lx / 2.0;
    double ymin = -Ly / 2.0;
    double etamin = -Leta / 2.0;

    uint64_t cells_above_energy_density_cutoff = 0;
    uint64_t cells_below_or_equal_energy_density_cutoff = 0;

    // Same grid-metadata first line as output_to_file() for downstream compatibility.
    fout << "#0 " << dx << " " << dy << " " << deta
         << " 0 " << xmin << " " << ymin << " " << etamin << "\n";
    fout << "# raw (un-diagonalized) contravariant T^{mu nu} in "
         << coordinate_system << " coordinates; indices 0=tau,1=x,2=y,3=eta\n";
    fout << "# columns: x y eta"
         << " Ttt Ttx Tty Ttn Txx Txy Txn Tyy Tyn Tnn"
         << " jB0 jB1 jB2 jB3 jS0 jS1 jS2 jS3 jQ0 jQ1 jQ2 jQ3\n";

    std::ostringstream buffer;

    for (int ix = 0; ix < nx; ++ix) {
        double x = ix * dx + xmin;
        for (int iy = 0; iy < ny; ++iy) {
            double y = iy * dy + ymin;
            for (int ieta = 0; ieta < neta; ++ieta) {
                double eta = ieta * deta + etamin;
                const Mat4x4& T = smearer.Tmunu[ix][iy][ieta];

                if (T[0][0] <= energy_density_cutoff) {
                    ++cells_below_or_equal_energy_density_cutoff;
                    continue;
                }
                ++cells_above_energy_density_cutoff;

                buffer << std::fixed << std::setprecision(8)
                       << x << " " << y << " " << eta << " "
                       << T[0][0] << " " << T[0][1] << " " << T[0][2] << " " << T[0][3] << " "
                       << T[1][1] << " " << T[1][2] << " " << T[1][3] << " "
                       << T[2][2] << " " << T[2][3] << " "
                       << T[3][3] << " "
                       << smearer.j0[ix][iy][ieta]  << " " << smearer.j1[ix][iy][ieta]  << " " << smearer.j2[ix][iy][ieta]  << " " << smearer.j3[ix][iy][ieta]  << " "
                       << smearer.j0s[ix][iy][ieta] << " " << smearer.j1s[ix][iy][ieta] << " " << smearer.j2s[ix][iy][ieta] << " " << smearer.j3s[ix][iy][ieta] << " "
                       << smearer.j0e[ix][iy][ieta] << " " << smearer.j1e[ix][iy][ieta] << " " << smearer.j2e[ix][iy][ieta] << " " << smearer.j3e[ix][iy][ieta]
                       << "\n";
            }
        }
    }

    fout << "# cells_above_energy_density_cutoff " << energy_density_cutoff << " : " << cells_above_energy_density_cutoff << "\n";
    fout << "# cells_below_or_equal_energy_density_cutoff " << energy_density_cutoff << " : " << cells_below_or_equal_energy_density_cutoff << "\n";
    fout << buffer.str();
    fout.close();
}

// ---- Fast eccentricity-only sigma scan ----------------------------------------
// Reads partons + free-streams ONCE; for each sigma it refills T^{mu nu}, diagonalizes,
// and computes the recentered energy-density eccentricities (eps2..5, eta-integrated
// and midrapidity) directly from hydro.TmunuOut[..].eps. No grid/IC is written.
void AMPTGenesis::run_ecc_scan(const std::vector<double>& sigmas, const std::string& ecc_out){
    std::cout << "[ecc_scan] AMPT input folder: " << input_folder_path << std::endl;
    AMPTSmearer smearer(input_folder_path, smearing_k, nx, ny, neta, Lx, Ly, Leta,
                        sigma_r, sigma_eta, tau0, rxy, reta, coordinate_system);
    smearer.parse_history();
    smearer.propagate(tau0);          // sigma-independent: done ONCE
    std::cout << "[ecc_scan] free-streaming done; scanning " << sigmas.size() << " sigma values" << std::endl;

    std::ofstream fout(ecc_out);
    fout << "# sigma eps2_int eps3_int eps4_int eps5_int eps2_mid eps3_mid eps4_mid eps5_mid\n";

    const int NMAX = 5;
    double dx = Lx/(nx-1.0), dy = Ly/(ny-1.0), deta = Leta/(neta-1.0);
    double xmin = -Lx/2.0, ymin = -Ly/2.0, etamin = -Leta/2.0;
    int ieta0 = (int)((0.0 - etamin)/deta + 0.5);
    if (ieta0 < 0) ieta0 = 0; if (ieta0 >= neta) ieta0 = neta-1;

    for (double s : sigmas) {
        smearer.sigma_r = s; smearer.sigma_eta = s;   // spline kernel reads the MEMBERS
        smearer.fill_Tmunu(s, s);
        Hydrodynamizer hydro(smearer.Tmunu, tau0, nx, ny, neta, coordinate_system);
        hydro.diagonalize();
        // pass 1: energy-weighted transverse centroid
        double Wi=0,Wix=0,Wiy=0,Wm=0,Wmx=0,Wmy=0;
        for (int ix=0; ix<nx; ++ix){ double x=ix*dx+xmin;
          for (int iy=0; iy<ny; ++iy){ double y=iy*dy+ymin;
            for (int ie=0; ie<neta; ++ie){
              double e = hydro.TmunuOut[ix][iy][ie].eps;
              if (e <= 0.0) continue;   // no energy floor: use ALL cells
              Wi+=e; Wix+=e*x; Wiy+=e*y;
              if (ie==ieta0){ Wm+=e; Wmx+=e*x; Wmy+=e*y; }
            } } }
        double x0i=(Wi>0)?Wix/Wi:0.0, y0i=(Wi>0)?Wiy/Wi:0.0;
        double x0m=(Wm>0)?Wmx/Wm:0.0, y0m=(Wm>0)?Wmy/Wm:0.0;
        // pass 2: r^n-weighted complex moments
        double ci[NMAX+1]={0},si[NMAX+1]={0},di[NMAX+1]={0};
        double cm[NMAX+1]={0},sm[NMAX+1]={0},dm[NMAX+1]={0};
        for (int ix=0; ix<nx; ++ix){ double x=ix*dx+xmin;
          for (int iy=0; iy<ny; ++iy){ double y=iy*dy+ymin;
            for (int ie=0; ie<neta; ++ie){
              double e = hydro.TmunuOut[ix][iy][ie].eps;
              if (e <= 0.0) continue;   // no energy floor: use ALL cells
              { double xr=x-x0i, yr=y-y0i; double r=std::sqrt(xr*xr+yr*yr), phi=std::atan2(yr,xr);
                for (int n=1;n<=NMAX;++n){ double rn=std::pow(r,n); di[n]+=e*rn; ci[n]+=e*rn*std::cos(n*phi); si[n]+=e*rn*std::sin(n*phi); } }
              if (ie==ieta0){ double xr=x-x0m, yr=y-y0m; double r=std::sqrt(xr*xr+yr*yr), phi=std::atan2(yr,xr);
                for (int n=1;n<=NMAX;++n){ double rn=std::pow(r,n); dm[n]+=e*rn; cm[n]+=e*rn*std::cos(n*phi); sm[n]+=e*rn*std::sin(n*phi); } }
            } } }
        fout << s;
        for (int n=2;n<=NMAX;++n){ double v=(di[n]>0)?std::sqrt(ci[n]*ci[n]+si[n]*si[n])/di[n]:0.0; fout << " " << v; }
        for (int n=2;n<=NMAX;++n){ double v=(dm[n]>0)?std::sqrt(cm[n]*cm[n]+sm[n]*sm[n])/dm[n]:0.0; fout << " " << v; }
        fout << "\n" << std::flush;
        std::cout << "[ecc_scan] sigma=" << s << " done" << std::endl;
    }
    fout.close();
    std::cout << "[ecc_scan] wrote " << sigmas.size() << " sigma rows to " << ecc_out << std::endl;
}

// ---- GRID-FREE analytic eccentricity scan -------------------------------------
// Sum over partons. Numerator <r^n e^{i n phi}> = sum_i w_i z_i^n is sigma-INDEPENDENT
// (the isotropic kernel's mean is the parton position). Denominator <r^n> = sum_i w_i
// sigma^n G_n(r_i/sigma), where G_n(u) = int |u xhat + q|^n W2d(q) d^2q is the cubic-spline
// kernel's radial moment (precomputed by direct integration of the SAME kernel as fill_Tmunu).
// w_i = p^tau (T^{tau tau} weight); midrapidity uses w_i *= s1(|eta_s,i|/sigma_eta).
void AMPTGenesis::run_ecc_analytic(const std::vector<double>& sigmas, const std::string& ecc_out){
    std::cout << "[ecc_analytic] AMPT input folder: " << input_folder_path << std::endl;
    AMPTSmearer smearer(input_folder_path, smearing_k, nx, ny, neta, Lx, Ly, Leta,
                        sigma_r, sigma_eta, tau0, rxy, reta, coordinate_system);
    smearer.parse_history();
    smearer.propagate(tau0);
    const int Np = (int)smearer.thermalized_partons.size();
    std::vector<double> px(Np), py(Np), pe(Np), pw(Np);
    for (int i=0;i<Np;++i){ auto& P = smearer.thermalized_partons[i];
        px[i]=P.x[0]; py[i]=P.x[1]; pe[i]=P.eta_s; pw[i]=P.ptau; }
    std::cout << "[ecc_analytic] " << Np << " partons; scanning " << sigmas.size() << " sigma" << std::endl;

    // cubic-spline shapes (M4): 2D normalization 5/(14 pi), 1D 1/6 (norm cancels in ratios)
    auto s2 = [](double q){ q=std::fabs(q); if(q<=1.0) return std::pow(2-q,3)-4*std::pow(1-q,3);
                            if(q<=2.0) return std::pow(2-q,3); return 0.0; };
    auto s1 = [](double p){ p=std::fabs(p); if(p<=1.0) return std::pow(2-p,3)-4*std::pow(1-p,3);
                            if(p<=2.0) return std::pow(2-p,3); return 0.0; };
    const double N2 = 5.0/(14.0*M_PI);

    // precompute G_n(u), n=2..5, on a uniform u grid by integrating the 2D kernel
    const int NU=600; const double UMAX=40.0; const double du=UMAX/(NU-1);
    std::vector<std::array<double,6>> G(NU);
    const int NQ=400, NTH=160; const double dq=2.0/NQ, dth=2.0*M_PI/NTH;
    for (int iu=0; iu<NU; ++iu){ double u=iu*du; std::array<double,6> acc={0,0,0,0,0,0};
        for (int iq=0; iq<NQ; ++iq){ double q=(iq+0.5)*dq; double wq=N2*s2(q)*q*dq;
            for (int it=0; it<NTH; ++it){ double th=(it+0.5)*dth;
                double d2=u*u+q*q+2*u*q*std::cos(th); double d=std::sqrt(d2); double w=wq*dth;
                acc[2]+=w*d2; acc[3]+=w*d2*d; acc[4]+=w*d2*d2; acc[5]+=w*d2*d2*d; } }
        G[iu]=acc; }
    auto Gn = [&](int n,double u){ if(u<0)u=0; double f=u/du; int i=(int)f;
        if(i>=NU-1) return G[NU-1][n]; double t=f-i; return G[i][n]*(1-t)+G[i+1][n]*t; };

    std::ofstream fout(ecc_out);
    { double ICE=0; for(int i=0;i<Np;++i) ICE+=pw[i];   // IC transverse-integrated energy (sum p^tau)
      fout << "# IC_E_total " << ICE << "\n"; }
    fout << "# sigma eps2_int eps3_int eps4_int eps5_int eps2_mid eps3_mid eps4_mid eps5_mid\n";
    for (double s : sigmas){
        double out_i[6], out_m[6];
        // two passes: pass=0 eta-integrated (w=ptau), pass=1 midrapidity (w=ptau*s1)
        for (int pass=0; pass<2; ++pass){
            std::vector<double> w(Np);
            double W=0, Xc=0, Yc=0;
            for (int i=0;i<Np;++i){ w[i] = (pass==0) ? pw[i] : pw[i]*s1(std::fabs(pe[i])/s);
                W+=w[i]; Xc+=w[i]*px[i]; Yc+=w[i]*py[i]; }
            if (W<=0){ for(int n=2;n<=5;++n){ (pass?out_m:out_i)[n]=0.0; } continue; }
            Xc/=W; Yc/=W;
            double Nre[6]={0},Nim[6]={0},D[6]={0};
            for (int i=0;i<Np;++i){ double xr=px[i]-Xc, yr=py[i]-Yc;
                double r=std::sqrt(xr*xr+yr*yr), phi=std::atan2(yr,xr);
                for (int n=2;n<=5;++n){ double rn=std::pow(r,n);
                    Nre[n]+=w[i]*rn*std::cos(n*phi); Nim[n]+=w[i]*rn*std::sin(n*phi);
                    D[n]+=w[i]*std::pow(s,n)*Gn(n, r/s); } }
            for (int n=2;n<=5;++n){ double v=(D[n]>0)?std::sqrt(Nre[n]*Nre[n]+Nim[n]*Nim[n])/D[n]:0.0;
                (pass?out_m:out_i)[n]=v; }
        }
        fout << s;
        for (int n=2;n<=5;++n) fout << " " << out_i[n];
        for (int n=2;n<=5;++n) fout << " " << out_m[n];
        fout << "\n" << std::flush;
    }
    fout.close();
    std::cout << "[ecc_analytic] wrote " << sigmas.size() << " rows to " << ecc_out << std::endl;
}

// ---- spatial eccentricity vs momentum anisotropy (grid, single sigma) ---------
// Writes one line:
//   eps2_Ttt eps3_Ttt  eps2_Lan eps3_Lan  ep_ideal_mid ep_full_mid  ep_ideal_int ep_full_int
// spatial eps_n: energy-weighted recentered |<r^n e^{i n phi}>|/<r^n> at eta_s=0,
//   with weight = deposited T^{tau tau} ("Ttt") or Landau energy density ("Lan").
// ep = |int (Txx-Tyy+2i Txy) / int (Txx+Tyy)|, T^ij ideal=(e+p)u^i u^j + p (p=e/3)
//   or full = deposited T^ij; "mid"=eta_s=0 slice, "int"=all eta_s (tau=tau0 cancels).
void AMPTGenesis::run_ecc_momentum(const std::string& ecc_out){
    std::cout << "[ecc_mom] " << input_folder_path << " (sigma_r=" << sigma_r << ")" << std::endl;
    AMPTSmearer smearer(input_folder_path, smearing_k, nx, ny, neta, Lx, Ly, Leta,
                        sigma_r, sigma_eta, tau0, rxy, reta, coordinate_system);
    smearer.parse_history();
    smearer.propagate(tau0);
    smearer.fill_Tmunu(sigma_r, sigma_eta);
    Hydrodynamizer hydro(smearer.Tmunu, tau0, nx, ny, neta, coordinate_system);
    hydro.diagonalize();

    double dx = Lx/(nx-1.0), dy = Ly/(ny-1.0), deta = Leta/(neta-1.0);
    double xmin = -Lx/2.0, ymin = -Ly/2.0, etamin = -Leta/2.0;
    int ie0 = (int)((0.0 - etamin)/deta + 0.5); if(ie0<0)ie0=0; if(ie0>=neta)ie0=neta-1;

    // ---- spatial eccentricity at eta_s=0 (two-pass) for two weights ----
    auto spatial = [&](bool useLandau, double& e2, double& e3){
        double W=0,Wx=0,Wy=0;
        for(int ix=0;ix<nx;++ix){ double x=ix*dx+xmin;
          for(int iy=0;iy<ny;++iy){ double y=iy*dy+ymin;
            double w = useLandau ? hydro.TmunuOut[ix][iy][ie0].eps : smearer.Tmunu[ix][iy][ie0][0][0];
            if(w<=0) continue; W+=w; Wx+=w*x; Wy+=w*y; } }
        if(W<=0){ e2=e3=0; return; }
        double x0=Wx/W, y0=Wy/W;
        double c2=0,s2=0,d2=0,c3=0,s3=0,d3=0;
        for(int ix=0;ix<nx;++ix){ double x=ix*dx+xmin;
          for(int iy=0;iy<ny;++iy){ double y=iy*dy+ymin;
            double w = useLandau ? hydro.TmunuOut[ix][iy][ie0].eps : smearer.Tmunu[ix][iy][ie0][0][0];
            if(w<=0) continue;
            double xr=x-x0, yr=y-y0, r=std::sqrt(xr*xr+yr*yr), phi=std::atan2(yr,xr);
            double r2=r*r, r3=r2*r;
            d2+=w*r2; c2+=w*r2*std::cos(2*phi); s2+=w*r2*std::sin(2*phi);
            d3+=w*r3; c3+=w*r3*std::cos(3*phi); s3+=w*r3*std::sin(3*phi); } }
        e2=(d2>0)?std::sqrt(c2*c2+s2*s2)/d2:0; e3=(d3>0)?std::sqrt(c3*c3+s3*s3)/d3:0;
    };
    double e2T,e3T,e2L,e3L; spatial(false,e2T,e3T); spatial(true,e2L,e3L);

    // ---- momentum anisotropy ep (ideal & full) x (midrapidity & integrated) ----
    auto epval=[&](double num_re,double num_im,double den){ return (den>0)?std::sqrt(num_re*num_re+num_im*num_im)/den:0.0; };
    double Fre_m=0,Fim_m=0,Fd_m=0, Ire_m=0,Iim_m=0,Id_m=0;   // midrapidity (ie0)
    double Fre_i=0,Fim_i=0,Fd_i=0, Ire_i=0,Iim_i=0,Id_i=0;   // integrated (all eta)
    for(int ix=0;ix<nx;++ix)
      for(int iy=0;iy<ny;++iy)
        for(int ie=0;ie<neta;++ie){
            double Txx=smearer.Tmunu[ix][iy][ie][1][1];
            double Tyy=smearer.Tmunu[ix][iy][ie][2][2];
            double Txy=smearer.Tmunu[ix][iy][ie][1][2];
            double e=hydro.TmunuOut[ix][iy][ie].eps;
            double p=e/3.0, ux=hydro.TmunuOut[ix][iy][ie].u[0], uy=hydro.TmunuOut[ix][iy][ie].u[1];
            double Txxi=(e+p)*ux*ux+p, Tyyi=(e+p)*uy*uy+p, Txyi=(e+p)*ux*uy;
            // integrated (tau=tau0 cancels)
            Fre_i+=(Txx-Tyy); Fim_i+=2*Txy; Fd_i+=(Txx+Tyy);
            Ire_i+=(Txxi-Tyyi); Iim_i+=2*Txyi; Id_i+=(Txxi+Tyyi);
            if(ie==ie0){ Fre_m+=(Txx-Tyy); Fim_m+=2*Txy; Fd_m+=(Txx+Tyy);
                         Ire_m+=(Txxi-Tyyi); Iim_m+=2*Txyi; Id_m+=(Txxi+Tyyi); }
        }
    double ep_id_mid=epval(Ire_m,Iim_m,Id_m), ep_full_mid=epval(Fre_m,Fim_m,Fd_m);
    double ep_id_int=epval(Ire_i,Iim_i,Id_i), ep_full_int=epval(Fre_i,Fim_i,Fd_i);

    std::ofstream fout(ecc_out);
    fout << "# eps2_Ttt eps3_Ttt eps2_Lan eps3_Lan ep_ideal_mid ep_full_mid ep_ideal_int ep_full_int\n";
    fout << e2T<<" "<<e3T<<" "<<e2L<<" "<<e3L<<" "
         << ep_id_mid<<" "<<ep_full_mid<<" "<<ep_id_int<<" "<<ep_full_int<<"\n";
    fout.close();
    std::cout << "[ecc_mom] wrote " << ecc_out << std::endl;
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
                double utau;

                if (coordinate_system == "cartesian"){
                    utau =  sqrt(pow(hydro.TmunuOut[ix][iy][ieta].u[0],2)
                                  +pow(hydro.TmunuOut[ix][iy][ieta].u[1],2)
                                  +pow(hydro.TmunuOut[ix][iy][ieta].u[2],2)+1.);
                }
                else if (coordinate_system == "hyperbolic"){
                    utau = sqrt(pow(hydro.TmunuOut[ix][iy][ieta].u[0],2)
                                  +pow(hydro.TmunuOut[ix][iy][ieta].u[1],2)
                                  +pow(tau0*hydro.TmunuOut[ix][iy][ieta].u[2],2)+1.);
                }
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

                if (coordinate_system == "cartesian"){
                    Ttaut += (Ttautau)*dx*dy*deta;
                }
                else if (coordinate_system == "hyperbolic"){
                    Ttaut += (Ttautau*cosheta + tau0*Ttaueta*sinheta)*dx*dy*deta*tau0; //GeV
                }
    }

    return Ttaut;
}

void AMPTGenesis::run_genesis(){

    std::cout<< "[" << input_format << "] input folder:  "<< input_folder_path << std::endl;
    AMPTSmearer smearer(input_folder_path,smearing_k,nx,
        ny,neta,Lx,Ly,Leta,sigma_r,sigma_eta,tau0,rxy,reta,coordinate_system);
    // Only the parsing step differs between AMPT and SMASH; everything after
    // (free-streaming to tau0, tensor build, diagonalization, output) is identical.
    if (input_format == "smash") {
        smearer.parse_smash();
    } else if (input_format == "ampt") {
        smearer.parse_history();
    } else {
        std::cerr << "[ERROR]: unknown input.format '" << input_format
                  << "' (expected 'ampt' or 'smash')" << std::endl;
        exit(1);
    }
    std::cout << "[INFO]: Input parsed (" << input_format << ")" << std::endl << std::flush;
    smearer.backpropagate = backpropagate;
    smearer.propagate(tau0);
    std::cout << "[INFO]: Finished free-streaming step" << std::endl << std::flush;
    smearer.fill_Tmunu(sigma_r,sigma_eta);
    std::cout << "[INFO]: Finished filling tensor" << std::endl << std::flush;

    // Raw-tensor mode: dump the independent T^{mu nu} components directly and skip
    // the Landau diagonalization entirely (the downstream code does its own matching).
    if (output_raw_tmunu) {
        std::cout << "[INFO]: Raw T^{mu nu} output mode (no Landau diagonalization)." << std::endl << std::flush;
        // Sanity check: total energy reconstructed from the raw tensor.
        double dx_ = Lx/(nx-1.), dy_ = Ly/(ny-1.), deta_ = Leta/(neta-1.);
        double Etot = 0.;
        for (int ix=0; ix<nx; ++ix)
        for (int iy=0; iy<ny; ++iy)
        for (int ieta=0; ieta<neta; ++ieta){
            double eta = ieta*deta_ - Leta*.5;
            const Mat4x4& T = smearer.Tmunu[ix][iy][ieta];
            if (coordinate_system == "hyperbolic")
                Etot += (T[0][0]*cosh(eta) + tau0*T[0][3]*sinh(eta))*dx_*dy_*deta_*tau0;
            else
                Etot += T[0][0]*dx_*dy_*deta_;
        }
        std::cout << "Target total energy: " << smearer.net_p[0] << std::endl;
        std::cout << "Actual total energy (from raw T^{mu nu}): " << Etot << std::endl;
        output_tmunu_to_file(smearer);
        std::cout << "[INFO]: Finished output" << std::endl << std::flush;
        return;
    }

    Hydrodynamizer hydro(smearer.Tmunu,
                        tau0,
                        nx,
                        ny,
                        neta,
                        coordinate_system);

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