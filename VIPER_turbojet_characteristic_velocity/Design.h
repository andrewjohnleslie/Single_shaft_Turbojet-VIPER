//
// Created by Andrew on 05/06/2023.
//

#ifndef VIPER_TURBOJET_DESIGN_DESIGN_H
#define VIPER_TURBOJET_DESIGN_DESIGN_H

#include <cmath>

class Design {

    // Public inheritance: Public member variables in the base class become protected in the derived class; therefore, not
    //                      accessible.
public:
// INTAKE
// =========================
// =========================

    //Member variables
    double gama = 1.4;
    double gamma_e = 1.3;
    double cp = 1005;
    double R = 287;
    double cpe = (gamma_e * R)/(gamma_e - 1);          // cpe = 1244 J/kgK

// Ambient conditions - SEA-LEVEL STATIC TEST BED (0 km, 0 m/s)
    /*
    double h = 10000;
    double T = 15 - (0.0065 * h);
    double T_a = T + 273.15;
    double T_sls = 288.15;
    double P_sls = 101300;

    double P_a = P_sls * (pow( (1 - (0.0065 * (h/T_sls))), 5.2561));

    double P_trop = 22625.04;
    double g = 9.81;

    //double T_a = 216.65;
    double ind = ((-1 * g)/(R * T_a)) * (h - 11000);
    //double P_a = P_trop * exp(ind);
    */

    double T_a = 288.15;                                               // At h=0 this is equivalent to T_sls_std
    double P_a = 101300;                                               // At h=0 this is equivalent to P_sls_std

    double M_f = 0;
    double a = pow((gama * R * T_a), 0.5);
    double V_f = a * M_f;

    //Member variables
    double T_02 = T_a * (1 + ((gama - 1) / 2) * pow(M_f, 2));
    double P_02 = P_a * pow((1 + ((gama - 1) / 2) * pow(M_f, 2)),gama/(gama - 1));
    double rho = P_02 / (R * T_02);


// COMPRESSOR
// =======================
// =======================

    double pi_c = 5.5; //4.7194693; //h = 7.8981219; // OPR = P_03/P_02
    double e_c = 0.9;

    double T_03 = T_02 * pow(pi_c, (gama - 1)/(e_c * gama));
    double P_03 = pi_c * P_02;
    double tau_c = T_03 / T_02;

    double DelT_c = T_03 - T_02;
    double DelT_t = (cp/cpe) * DelT_c;


// TURBINE
// =======================
// =======================

    double e_t = 0.85;
    double mdot_a = 23.81; //31.143808; //h = 8.9202948;
    double mdot_f = 0.4267; //0.53067478;  //H = 0.17371925;
    double mdot_tot = mdot_a + mdot_f;
    double LCV = 43000000;

    double energy_rel = mdot_f * LCV;               //[Energy_rel = Qr]
    // double sp_energy_rel = energy_rel/mdot_a                             [per kg of air-flow]

    double T_04 = ((energy_rel + (mdot_a * cp * (T_03 - 298))) / ((mdot_a + mdot_f) * cpe)) + 298;     // T04 = 1063.15;
    double P_04 = P_03;                          // Pressure loss in the combustor is neglected: P_04 = (1 - delP_c) * P_03

    double T_05 = T_04 - DelT_t;
    double P_05 = P_04 * pow((T_05/T_04), (gamma_e / (e_t * (gamma_e-1))));

    double x = P_05/P_04;                       // Turbine pressure ratio
    double y = P_05/P_a;                        // Nozzle pressure ratio

    //double Pcrit = 1 / pow((1 - ((gamma_e - 1) / (gamma_e + 1))), (gamma_e / (gamma_e - 1)));
    double Pcrit = pow(((gamma_e + 1)/2), ((gamma_e)/(gamma_e - 1)));                  // Equivalent to P05/Pcrit

    double m_norm_ge = (gamma_e / pow((gamma_e - 1), 0.5)) * pow(((gamma_e + 1) / 2), ((-1 * (gamma_e + 1)) / (2 * (gamma_e - 1))));      // norm_mf = 1.389 for gamma = 1.3
        // m_norm_4 == m_norm_9

    double A_turb =  ((mdot_a * pow((cpe * T_04), 0.5)) / (P_04 * m_norm_ge));
    double A_noz =  ((mdot_a * pow((cpe * T_05), 0.5)) / (P_05 * m_norm_ge));     // THIS COULD BE USED TO WORK OUT K_H
    double A_ratio = A_turb / A_noz;

    double k_H = (T_04 - T_05) / T_04;

    //double choke_t = P_05/Pcrit;        // If P_04/P* >= 1.832, the engine is choked for gamma values = 1.3!  pg 177

    double SigT = (2 * e_t * (gamma_e - 1)) / ((2 * gamma_e) - (e_t * (gamma_e - 1)));        // pg 180
    double SigP = (2 * gamma_e) / ((2 * gamma_e) - (e_t * (gamma_e - 1)));                    // pg 180

    double pi_t = pow((A_turb/A_noz), SigP);
    double tau_t = pow((A_turb/A_noz), SigT);

    double P_4 = P_04 / pow((1 + (0.5 * (gamma_e-1)) * pow(1,2)), gamma_e / (gamma_e-1));
    double T_4 = T_04 / (1 + (0.5 * (gamma_e-1)) * pow(1,2));       // The 1^2 represents choked conditions!
    double V_4 = pow(gamma_e * R * T_4, 0.5);
    double rho_4 = P_4 / (R * T_4);
    double A_4 = mdot_tot / (rho_4 * V_4);


// NOZZLE
// =======================
// =======================

    double T_09 = T_05;

    void myCustomFunction() {
        if ( y < Pcrit){                   // If P05/Pa < P05/Pcrit (Nozzle PR < Critical PR)
            std::cout << std::endl;
            std::cout << "Nozzle is not choked! \n" << std::endl;
        }else{
            std::cout << std::endl;
            std::cout << "Nozzle is choked! \n " << std::endl;
        }
    }

    void callmyCustomFunction(){
        myCustomFunction();
    }

    // The assumption of fully expanded flow to ambient conditions is used despite the equations bypassing the
    //      consideration for choked/unchoked flow. Therefore, the characteristic velocity is calculated to maintain
    //      consistency with the book!
    double P_9 = P_a;
    double T_9 = T_05 * pow((P_9/P_05), ((gamma_e - 1)/gamma_e));
    double V_9 = pow((2 * cpe * (T_05 - T_9)), 0.5);

    // The star variables represents actual static values when choked!
    double P_9_star = P_05 / pow((1 + (0.5 * (gamma_e-1)) * pow(1,2)), gamma_e / (gamma_e-1));
    double T_9_star = T_05 / (1 + (0.5 * (gamma_e-1)) * pow(1,2));
    double V_9_choke = pow(gamma_e * R * T_9_star, 0.5);
    double rho_9 = P_9 / (R * T_9_star);
    double A_9 = mdot_tot / (rho_9 * V_9);
    double A_ratio_fuel = A_4 / A_9;

    double Ft = (mdot_tot * V_9) - (mdot_a * V_f);
    double sfc = mdot_f / Ft;     // [kg/h/N]
    double Fs = Ft / mdot_a;

    /*
    double F_g = V_9 * mdot_tot;

        void  myCustomFunction1() {
            if (V_f == 0) {
                double Ft = F_g;
                double sfc = mdot_f / Ft;     // [kg/h/N]
                double Fs = Ft / mdot_a;

                std::cout << "Net thrust (F_T) [GROSS] = " << Ft << std::endl;
                std::cout << "Specific thrust (F_S) = " << Fs << std::endl;
                std::cout << "SFC = " << sfc << std::endl;

            } else {
                double Ft = (mdot_tot * V_9) - (mdot_a * V_f) + (A_noz * (P_05 - P_a));
                double sfc = mdot_f / Ft;     // [kg/h/N]
                double Fs = Ft / mdot_a;

                std::cout << "Net thrust (F_T) [NET] = " << Ft << std::endl;
                std::cout << "Specific thrust (F_S) = " << Fs << std::endl;
                std::cout << "SFC = " << sfc << std::endl;
            }
        }
    */

    double m2_bar = (mdot_a * pow(cp * T_02, 0.5)) / (P_02);
    double m2_bar_norm = m2_bar / m2_bar;

    /*
    void callmyCustomFunction1(){
        myCustomFunction1();
    }
    */

    Design() = default;          //Constructor
    ~Design() = default;         //Destructor
};
#endif //VIPER_TURBOJET_DESIGN_DESIGN_H
