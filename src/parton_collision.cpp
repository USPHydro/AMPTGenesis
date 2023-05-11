#include <array>
#include <cmath>

///\class PartonCollision
///\brief Stores information of a parton that collides in a given position


//typedef std::array<double,4> Vec4;
typedef std::array<double,3> Vec3;

class PartonCollision
{
private:

    static constexpr double T0 = 0.E-3;

public:
    int pid;            ///< PDG PID of parton
    double mass;         ///< Parton mass
    double t;            ///< time of the collision
    Vec3 x;             ///< position of the collision
    Vec3 incoming_mom;  ///< Incoming parton momentum
    Vec3 outgoing_mom;  ///< Outgoing parton momentum
    double tau;
    double eta_s;


    PartonCollision(int pid, double mass, double time, Vec3 position,
                    Vec3 incoming_momentum, Vec3 outgoing_momentum);
    ~PartonCollision() = default;
};

///\brief Constructor of a parton collision
PartonCollision::PartonCollision(int pid, double mass, double time, Vec3 position,
                                 Vec3 incoming_momentum, Vec3 outgoing_momentum){
    this->pid = pid;
    this->mass = mass;
    this->t = time + T0;
    this->x = position;
    this->incoming_mom = incoming_momentum;
    this->outgoing_mom = outgoing_momentum;
    this->eta_s = .5*log((this->t+this->x[2])/(this->t-this->x[2]));
    this->tau = sqrt(pow(this->t,2) - pow(this->x[2],2));
}


///\class PartonThermalized
///\brief Once it is determined where the parton will be thermalized, we store
///its position and momentum in this class
class PartonThermalized : public PartonCollision
{
public:
    double E;    ///< Particle momentum
    double mT;   ///< Particle transverse mass
    double Y;    ///< Particle rapidity
    double ptau; ///< zeroth component of the momentum in hyperbolic coordinates
    double peta; ///< third component of the momentum in hyperbolic coordinates
    PartonThermalized(int pid, double mass, double time, Vec3 position,
                        Vec3 momentum);
    ~PartonThermalized() =  default;
};

PartonThermalized::PartonThermalized(int pid, double mass, double time,
                                     Vec3 position, Vec3 momentum) :
    PartonCollision(pid,mass,time,position,momentum,momentum) //Initializes base class
{
    Vec3 p = this->outgoing_mom;
    this->E = sqrt(pow(this->mass,2) + pow(p[0],2) + pow(p[1],2) + pow(p[2],2));
    this->mT = sqrt(pow(this->E,2) - pow(p[2],2));
    this->Y = .5*log((this->E+p[2])/(this->E-p[2]));

    //See arXiv:2011.03740, second column of page 2
    this->ptau = this->mT*cosh(this->Y-this->eta_s);
    this->peta = this->mT*sinh(this->Y-this->eta_s)/this->tau;
}

