//
// Created by Andrew on 06/06/2023.
//

#ifndef VIPER_TURBOJET_DESIGN_OFF_DESIGN_H
#define VIPER_TURBOJET_DESIGN_OFF_DESIGN_H

#include <cmath>
#include <fstream>
#include "Design.h"

class Off_design : public Design {

public:
    // Set off-design conditions (... km, ... m/s)
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


    double M_f = 0.8;
    double a = pow((gama * R * T_a), 0.5);
    double V_f = a * M_f;

    //Member variables
    double T02_OD = T_a * (1 + ((gama - 1) / 2) * pow(M_f, 2));
    double P02_OD = P_a * pow((1 + ((gama - 1) / 2) * pow(M_f, 2)),gama/(gama - 1));

    double T04_OD = 1063.150858;                         // T_04 = 1063.15 k - MARK'S VARIATION VALUE 
    double T03_OD = T02_OD + ((cpe/cp) * k_H * T04_OD);
    double T05_OD = T04_OD * (1 - k_H);

    double piC_OD = pow((T03_OD / T02_OD), (e_c * (gama/(gama - 1))));
    //double piC_OD = 5.17;
    double z_OD = T04_OD/T02_OD;
    double P03_OD = piC_OD * P02_OD;
    double P04_OD = P03_OD;
    double P05_OD = P04_OD * pow((T05_OD / T04_OD), (gamma_e / (e_t * (gamma_e - 1))));

    double f = ((cpe * (T04_OD - 298)) - (cp * (T03_OD - 298))) / (LCV - (cpe * T04_OD) + (298 * cpe));

    double T09_OD = T05_OD;

    // The assumption of fully expanded flow to ambient conditions is used despite the equations bypassing the
    //      consideration for choked/unchoked flow. Therefore, the characteristic velocity is calculated to maintain
    //      consistency with the book!
    double P9_OD = P_a;
    double T9_OD = T05_OD * pow((P9_OD/P05_OD), ((gamma_e - 1)/gamma_e));
    double V9_OD = pow((2 * cpe * (T05_OD - T9_OD)), 0.5);

    //double T9_OD = T05_OD * pow((P_a /P05_OD), ((gamma_e - 1)/ gamma_e));
    //double V9_OD = 1 * pow((gamma_e * R * T9_OD), 0.5);

    double x_OD = P05_OD/P04_OD;
    double y_OD = P05_OD/P_a;

    double mdot_a_OD = ((mdot_a * pow(cp * (T_03 - T_02), 0.5)) / P_03) * (P03_OD / pow(cp * (T03_OD - T02_OD), 0.5));

    double f_OD = ((cpe * (T04_OD - 298)) - (cp * (T03_OD - 298))) / (LCV - (cpe * T04_OD) + (298 * cpe));

    //double mdot_f_OD = (mdot_a_OD * ((cpe * (T04_OD - 298)) - (cp * (T03_OD -298)))) / (LCV - (cp * (T04_OD - 298)));
    double mdot_f_OD = mdot_a_OD * f_OD;

    double mdot_tot_OD = mdot_a_OD + mdot_f_OD;

    //double Ft_OD = (mdot_tot_OD * V9_OD);
    //double Ft_OD = (mdot_tot_OD * V9_OD) - (mdot_a_OD * V_f) + (A_noz * (P05_OD - P_a));
    double Ft_OD = (mdot_tot_OD * V9_OD) - (mdot_a_OD * V_f);               // [N]
    double Fg_OD = (mdot_tot_OD * V9_OD);
    double alpha = Ft_OD / Ft;
    double sfc_OD = mdot_f_OD / Ft_OD;
    double Fs_OD = Ft_OD / mdot_a_OD;

    double m2_bar_OD = (mdot_a_OD * pow(cp * T02_OD, 0.5)) / (P02_OD);
    double m2_bar_norm = m2_bar_OD / m2_bar;

    void Output() {
        std::ofstream outputFile("output.csv");
        if (outputFile.is_open()) {
            //outputFile << std::fixed << std::setprecision(10);
            // Export variable values to the CSV file
            outputFile << h << "," << Ft_OD << "," << Fg_OD << "," << std::setprecision(10) << sfc_OD << "," << Fs_OD
                       << std::setprecision(8) << ","
                       << T03_OD << "," << mdot_a_OD << "," << mdot_f_OD << "," << f_OD << "," << z_OD << ","
                       << piC_OD << "," << m2_bar_norm << "," << alpha << std::endl;

            // Close the file
            outputFile.close();
            std::cout << "\nData exported to output.csv" << std::endl;
        } else {
            std::cerr << "\nError opening file." << std::endl;
        }


    }

    void callOutput(){
        Output();
    }

    Off_design() = default;
    ~Off_design() = default;

};
#endif //VIPER_TURBOJET_DESIGN_OFF_DESIGN_H
