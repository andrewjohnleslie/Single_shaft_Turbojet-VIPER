// This test case is the preliminary test case for a simple turbojet (single-spool) engine.
// The engine being considered and modelled is the Rolls-Royce Viper [Mark 601] turbojet engine.
// This program places the previous code into an inherited class structure.

// The data provided for this engine is at a DESIGN POINT for a sea-level test bed.
// F_T = F_G = 15167 (N) [Static]
// sfc = 0.993 (kg/h/kg)
// m_a = 23.81 (kg/s)
// m_f = 0.4267 (kg/s)
// OPR [P_03/P_02] = 5.5

#include <iostream>
#include <iomanip>
#include <fstream>
#include <vector>
#include <string>
#include "Design.h"
#include "Off_design.h"

int main() {
    std::cout << "DESIGN POINT CALCULATION " << std::endl;
    std::cout << "============================== " << std::endl;
    std::cout << "      " << std::endl;

    std::cout << "Ta = " << Design().T_a << std::endl;
    std::cout << "Pa = " << Design().P_a << std::endl;
    //std::cout << "Rho_a = " << Design().rho << std::endl;
    std::cout << "T_02 = " << Design().T_02 << std::endl;
    std::cout << "P_02 = " << Design().P_02 << std::endl;
    std::cout << "T_03 = " << Design().T_03 << std::endl;
    std::cout << "P_03 = " << Design().P_03 << std::endl;
    std::cout << "Compressor TR = " << Design().tau_c << std::endl;
    std::cout << "Compressor PR = " << Design().pi_c << std::endl;
    std::cout << "T_04 = " << Design().T_04 << std::setprecision(5) << std::endl;
    std::cout << "P_04 = " << Design().P_04 << std::endl;
    std::cout << "T_05 = " << Design().T_05 << std::endl;
    std::cout << "P_05 = " << Design().P_05 << std::endl;
    std::cout << "Turbine TR = " << Design().tau_t << std::endl;
    std::cout << "Turbine PR = " << Design().x << std::endl;
    std::cout << "T_09 = " << Design().T_09 << std::endl;
    std::cout << "Nozzle PR = " << Design().y << std::endl;
    std::cout << "T_9 = " << Design().T_9 << std::endl;
    std::cout << "P_9 = " << Design().P_9 << std::endl;
    std::cout << "V_9 = " << Design().V_9 << std::endl;

    std::cout << "\nNet thrust (F_T) [GROSS] = " << Design().Ft << std::endl;
    std::cout << "Specific thrust (F_S) = " << Design().Fs << std::endl;
    std::cout << "SFC = " << Design().sfc << std::endl;

    std::cout << "\nmbar_choke_ge = " << Design().m_norm_ge << std::endl;
    std::cout << "K_H = " << Design().k_H << std::endl;
    std::cout << "Turbine nozzle area - NO FUEL = " << Design().A_turb << std::endl;
    std::cout << "Propulsive nozzle area - NO FUEL = " << Design().A_noz << std::endl;
    std::cout << "Area ratio (A4/A9) = " << Design().A_ratio << std::endl;
    std::cout << "Turbine nozzle area - FUEL = " << Design().A_4 << std::endl;
    std::cout << "Nozzle exit area - FUEL = " << Design().A_9 << std::endl;
    std::cout << "Area ratio (A4/A9) = " << Design().A_ratio_fuel << std::endl;
    //std::cout << "Critical pressure ratio = " << Design().Pcrit << std::endl;
    //std::cout << "Critical pressure = " << Design().Pc << std::endl;
    std::cout << "mbar_2 = " << Design().m2_bar << std::endl;

    Design nozzle;
    nozzle.callmyCustomFunction();

    //Design preform;
    //preform.callmyCustomFunction1();

    std::cout << "\nOFF-DESIGN CALCULATION " << std::endl;
    std::cout << "============================== \n \n";

    std::cout << "Flight speed = " << Off_design().V_f << std::endl;
    std::cout << "Off-design: Ta = " << Off_design().T_a << std::endl;
    std::cout << "Off-design: Pa = " << Off_design().P_a << std::endl;
    std::cout << "Off-design: T_02 = " << Off_design().T02_OD << std::endl;
    std::cout << "Off-design: P_02 = " << Off_design().P02_OD << std::endl;
    std::cout << "Off-design: T_03 = " << Off_design().T03_OD << std::endl;
    std::cout << "Off-design: P_03 = " << Off_design().P03_OD << std::endl;
    std::cout << "Off-design: T_04 = " << Off_design().T04_OD << std::endl;
    std::cout << "Off-design: P_04 = " << Off_design().P03_OD << std::endl;
    std::cout << "Off-design: T_05 = " << Off_design().T05_OD << std::endl;
    std::cout << "Off-design: P_05 = " << Off_design().P05_OD << std::endl;
    std::cout << "Off-design: T_09 = " << Off_design().T09_OD << std::endl;
    std::cout << "Off-design: P_9 = " << Off_design().P9_OD << std::endl;
    std::cout << "Off-design: T_9 = " << Off_design().T9_OD << std::endl;
    std::cout << "Off-design: V_9 = " << Off_design().V9_OD << std::endl;

    std::cout << "\nOff-design: Turbine PR = " << Off_design().x_OD << std::endl;
    std::cout << "Off-design: Compressor PR = " << Off_design().piC_OD << std::endl;
    std::cout << "Off-design: T_04/T_02 = " << Off_design().z_OD << std::endl;
    std::cout << "Off-design: Nozzle PR = " << Off_design().y_OD << std::endl;

    std::cout << "\n========================================================== " << std::endl;
    std::cout << "Off-design: Intake mass flow = " << Off_design().mdot_a_OD << std::endl;
    std::cout << "Off-design: fuel-air ratio = " << Off_design().f_OD<< std::endl;
    std::cout << "Off-design: Fuel mass flow (sheet) = " << Off_design().mdot_f_OD << std::endl;
    std::cout << "Off-design: Total mass flow = " << Off_design().mdot_tot_OD << std::endl;
    std::cout << "Off-design: mass flow [m2_bar] = " << Off_design().m2_bar_OD << std::endl;
    std::cout << "Off-design: Corrected mass flow = " << Off_design().m2_bar_norm << std::endl;
    std::cout << "========================================================== " << std::endl;
    std::cout << "Off-design: Thrust (F_T) = " << Off_design().Ft_OD << std::endl;
    std::cout << "Thrust Lapse = " << Off_design().alpha << std::endl;
    std::cout << "Off-design: Specific thrust (F_S) = " << Off_design().Fs_OD << std::endl;
    std::cout << "Off-design: SFC = " << Off_design().sfc_OD << std::endl;

    Off_design OD;
    OD.callOutput();
}
