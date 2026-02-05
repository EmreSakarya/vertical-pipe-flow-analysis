import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from iapws import IAPWS97

# --- 1. INPUT PARAMETERS ---
L = 17.2                # Length [m]
D = 0.32                # Diameter [m]
Area = np.pi * (D**2) / 4.0
m_dot = 121.0           # Mass flow rate [kg/s]
G = m_dot / Area        # Mass flux [kg/m^2-s]
q_prime = 1.57e6        # Linear Heat Rate [W/m]
P_in = 1.5              # Inlet Pressure [MPa]
T_in = 420.0            # Inlet Temperature [K]
g = 9.81                # Gravity [m/s^2]

# --- 2. ANALYTICAL CALCULATION ---
def solve_analytical():
    state_in = IAPWS97(P=P_in, T=T_in)
    h_in = state_in.h
    
    # Energy Balance
    total_heat = (q_prime * L) / 1000.0  # kW
    h_out = h_in + (total_heat / m_dot)
    
    # Average Properties estimation
    state_out_est = IAPWS97(P=P_in, h=h_out)
    rho_avg = (state_in.rho + state_out_est.rho) / 2
    mu_avg = (state_in.mu + state_out_est.mu) / 2
    
    # Pressure Drop Calculation
    Re = (G * D) / mu_avg
    f = 0.3164 * (Re**-0.25)
    
    dP_grav = rho_avg * g * L
    dP_fric = f * (L/D) * (G**2) / (2 * rho_avg)
    
    # Downward Flow: Gravity increases P, Friction decreases P
    P_out = P_in + (dP_grav - dP_fric) * 1e-6
    T_out = state_out_est.T
    
    return P_out, T_out

# --- 3. NUMERICAL CALCULATION (HEM) ---
def solve_hem(step_size):
    z = 0.0
    P_curr = P_in
    h_curr = IAPWS97(P=P_in, T=T_in).h
    
    results = []
    
    # Save inlet point
    inlet = IAPWS97(P=P_in, h=h_curr)
    results.append({
        'z': 0.0, 
        'P': P_curr, 
        'T': inlet.T, 
        'rho': inlet.rho, 
        'cp': inlet.cp, 
        'x': 0.0
    })
    
    while z < L - 1e-6:
        dz = min(step_size, L - z)
        
        # Energy Equation
        dh = (q_prime * dz / 1000.0) / m_dot
        h_next = h_curr + dh
        
        # Momentum Equation
        st = IAPWS97(P=P_curr, h=h_curr)
        rho = st.rho
        mu = st.mu
        
        Re = (G * D) / mu
        if Re > 0:
            f = 0.3164 * (Re**-0.25)
        else:
            f = 0
        
        dP_grav = rho * g * dz
        dP_fric = f * (dz/D) * (G**2) / (2*rho)
        
        # Update Pressure (Pa -> MPa conversion with 1e-6)
        P_next = P_curr + (dP_grav - dP_fric) * 1e-6
        
        z += dz
        P_curr = P_next
        h_curr = h_next
        
        # Data storage
        try:
            st_new = IAPWS97(P=P_curr, h=h_curr)
            x_val = st_new.x if st_new.x is not None else 0.0
            
            # Cp is problematic in two-phase region for some libraries, handle gracefully
            try:
                cp_val = st_new.cp
            except:
                cp_val = None 
            
            # Saturation Temp for reference
            try:
                tsat = IAPWS97(P=P_curr, x=0).T
            except:
                tsat = None

            results.append({
                'z': z, 
                'P': P_curr, 
                'T': st_new.T, 
                'rho': st_new.rho, 
                'cp': cp_val, 
                'x': x_val, 
                'T_sat': tsat
            })
        except:
            pass

    return pd.DataFrame(results)

# --- 4. MAIN EXECUTION AND PLOTTING ---
if __name__ == "__main__":
    # A. Analytical Solution
    P_anal, T_anal = solve_analytical()
    
    print(f"{'Method':<25} {'P_out (MPa)':<15} {'T_out (K)':<15} {'Quality'}")
    print("-" * 65)
    print(f"{'Analytical':<25} {P_anal:<15.4f} {T_anal:<15.2f} {0.0:<15}")
    
    # B. Numerical Solution Loop
    steps = [1.0, 0.1, 0.01]
    hem_data = {}
    
    for s in steps:
        df = solve_hem(s)
        hem_data[s] = df
        last = df.iloc[-1]
        print(f"{f'Numerical (dz={s})':<25} {last['P']:<15.4f} {last['T']:<15.2f} {last['x']:<15.5f}")

    # C. Specific Cp Distribution Table (From your Appendix)
    print("\n\n--- CP DISTRIBUTION (Detailed Profile, dz = 0.01 m) ---")
    print(f"{'Z (m)':<10} | {'cp (kJ/kg.K)':<15}")
    print("-" * 30)
    
    df_fine = hem_data[0.01] # Use the finest mesh result
    targets = [0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 17.2]
    
    for t in targets:
        # Find the row closest to target z
        idx = (df_fine['z'] - t).abs().idxmin()
        row = df_fine.loc[idx]
        print(f"{row['z']:<10.2f} | {row['cp']:<15.4f}")

    # D. Plotting
    fig, axes = plt.subplots(2, 2, figsize=(12, 10))
    
    # Density
    ax = axes[0,0]
    for s in steps: ax.plot(hem_data[s]['z'], hem_data[s]['rho'], label=f'dz={s}')
    ax.set_title('(a) Density Profile')
    ax.set_ylabel('Density (kg/m^3)')
    ax.set_xlabel('Distance (m)')
    ax.legend()
    ax.grid(True)
    
    # Cp
    ax = axes[0,1]
    for s in steps: ax.plot(hem_data[s]['z'], hem_data[s]['cp'], label=f'dz={s}')
    ax.set_title('(b) Specific Heat (Cp)')
    ax.set_ylabel('Cp (kJ/kg.K)')
    ax.set_xlabel('Distance (m)')
    ax.grid(True)
    
    # Temperature
    ax = axes[1,0]
    df_final = hem_data[0.01]
    if 'T_sat' in df_final.columns and not df_final['T_sat'].isna().all():
        ax.plot(df_final['z'], df_final['T_sat'], 'k--', label='T_sat', alpha=0.5)
        
    for s in steps: ax.plot(hem_data[s]['z'], hem_data[s]['T'], label=f'dz={s}')
    ax.set_title('(c) Temperature Profile')
    ax.set_ylabel('Temperature (K)')
    ax.set_xlabel('Distance (m)')
    ax.legend()
    ax.grid(True)
    
    # Pressure
    ax = axes[1,1]
    for s in steps: ax.plot(hem_data[s]['z'], hem_data[s]['P'], label=f'dz={s}')
    ax.set_title('(d) Pressure Profile')
    ax.set_ylabel('Pressure (MPa)')
    ax.set_xlabel('Distance (m)')
    ax.grid(True)
    
    plt.tight_layout()
    plt.savefig('hem_results.png') # Saves the graph
    print("\nGraphs saved as 'hem_results.png'")
