import numpy as np
import matplotlib.pyplot as plt
import matplotlib.animation as ani
from scipy.integrate import odeint

def info_of_reaction(num, T, time):
    reactants_dict = {1:['Cyclopropane'], 2:['CH3COOC2H5', 'OH-'], 3:['NO2'], 4:['H2', 'O2']}
    products_dict = {1:['Propene'], 2:['CH3COO-', 'C2H5OH'], 3:['NO', 'O2'], 4:['H2O']}
    coefficient_dict = {1:[-1, 1], 2:[-1, -1, 1, 1], 3:[-2, 2, 1], 4:[-2, -1, 2]}
    Ea_dict = {1:272, 2:48.3, 3:114, 4:165} # kJ
    A_dict = {1:1.48e15, 2:3.2e7, 3:4.0e9, 4:1.0e13}

    reactant = np.array(reactants_dict[num])
    product  = np.array(products_dict[num])
    coefficient = np.array(coefficient_dict[num])
    Ea = Ea_dict[num]*1000
    A = A_dict[num]
    R = 8.3145
    k = A * np.exp(-Ea/(R*(T+273.15)))
    A0 = []
    for i in reactants_dict[num]:
        a = float(input(f'Enter the initial concentration of {i} (M): '))
        A0.append(a)
    species = np.array(list(reactant) + list(product))
    return reactant, product, species, coefficient, k, A0, time

def info_of_series_reaction(num, T, time):
    reactants_dict = {5:['Ethanol'], 6:['C6H6']}
    inter_dict = {5:['Acetaldehyde'], 6:['C6H5NO2']}
    products_dict = {5:['Acetate'], 6:['C6H4(NO2)2']} 
    Ea_dict = {5:[45, 42], 6:[100, 120]} # kJ
    A_dict = {5:[5.0e6, 8.0e6], 6:[1.0e13, 1.0e13]}
    reactant = np.array(reactants_dict[num])
    intermediate  = np.array(inter_dict[num])
    product = np.array(products_dict[num])
    Ea = np.array(Ea_dict[num])*1000
    A = np.array(A_dict[num])
    R = 8.3145
    k = A * np.exp(-Ea/(R*(T+273.15)))
    A0 = []
    for i in reactants_dict[num]:
        a = float(input(f'Enter the initial concentration of {i} (M): '))
        A0.append(a)
    species = np.array(list(reactant) + list(intermediate) + list(product))
    return reactant, intermediate, product, species, k, A0, time

def reaction_model(C, t, k, num_of_reactants, coefficient):
    rate = k
    for i in range(num_of_reactants):
            rate *= (max(0, C[i])) ** abs(coefficient[i])
    dCdt = []
    for i in coefficient:
        dCdt.append(i * rate)#各物質濃度變化率
    return dCdt

def series_reaction_model(C, t, k):
    rate = k.copy()
    for i in range(len(rate)):
        rate[i] *= (max(0, C[i]))
    dCdt = []
    dCdt.append(-rate[0])
    dCdt.append(rate[0] - rate[1])
    dCdt.append(rate[1])
    return dCdt



def concentration_update(reactant, product, coefficient, k, A0, time):
    num_of_reactats = len(reactant)
    list_C0 = list(A0)
    for i in range(len(product)):
        list_C0.append(0)#各化合物初濃度
    C0 = np.array(list_C0)
    t  = np.linspace(0, time, 1000)
    C = odeint(reaction_model, C0, t, args=(k, num_of_reactats, coefficient))

    C[C < 0] = 0
    return t, C
    
def concentration_update_series(intermediate, product, k, A0, time):
    list_C0 = list(A0)
    for i in range(len(intermediate)+len(product)):
        list_C0.append(0)
    C0 = np.array(list_C0)
    t = np.linspace(0, time, 1000)
    C = odeint(series_reaction_model, C0, t, args = (k, ))

    C[C < 0] = 0
    return t, C



def format_reaction_str(names, coefs):
    parts = []
    for name, coef in zip(names, coefs):
        c = abs(coef)
        if c == 1:
            parts.append(name)
        else:
            parts.append(f"{int(c)} {name}")
    return " + ".join(parts)

def format_series_reaction_str(names):
    parts = []
    for name in names:
        parts.append(name)

    return " + ".join(parts)



print(  "1: Cyclopropane --> Propene \n" \
                "2: CH3COOC2H5 + OH- --> CH3COO- + C2H5OH \n" \
                "3: 2 NO2 --> 2 NO + O2  \n" \
                "4: 2 H2 + O2 --> 2 H2O \n"\
                "5: 乙醇代謝模擬: Ethanol --> Acetaldehyde --> Acetate\n"\
                "6: 苯的連續硝化反應: C6H6 --> C6H5NO2 --> C6H4(NO2)2"

                )
num = int(input("Choose a reaction (number): "))
T = float(input("Enter the reacting temperature(degree Celsius): "))
time = float(input("Enter how long the reaction last (s): "))

if num >=1 and num <=4: 
    reactant, product, species, coefficient, k, A0, time= info_of_reaction(num, T, time)
    print('k =', f'{k:.3e}')

    num_of_subs = len(reactant) + len(product)
    t, C = concentration_update(reactant, product, coefficient, k, A0, time)

    fig, ax = plt.subplots(figsize=(8, 6))

    species = list(reactant) + list(product)
    lines = []
    line_styles = ['-', '--', '-.', ':']
    for i in range(num_of_subs):
        line, = ax.plot([], [], lw=2, label=species[i], linestyle=line_styles[i % 4])
        lines.append(line)

    ax.set_xlim(0, time)
    ax.set_ylim(0, max(A0)*1.2)
    ax.set_xlabel('time (s)', fontsize=12)
    ax.set_ylabel('concentration (M)', fontsize=12)
    ax.grid(True, linestyle='--')
    ax.legend()

    coef_reactant = coefficient[:len(reactant)]
    coef_product = coefficient[len(reactant):]
    str_reactants = format_reaction_str(reactant, coef_reactant)
    str_products = format_reaction_str(product, coef_product)
    fig.suptitle(f"{str_reactants} --> {str_products}", fontsize=14, fontweight='bold')
    ax.set_title('Reaction Kinetics Simulation', fontsize=12)

    if k >= 100:
        ax.text(0.05, 0.95, f"k = {k:.2e}, Warning:k is big, the reaction may finish in a short time.", transform=ax.transAxes,
            fontsize=12, color="red", verticalalignment="top")
    elif k < 0.001:
        ax.text(0.05, 0.95, f"k = {k:.3e}, Warning:k is small, concentration may change very slowly.", transform=ax.transAxes,
            fontsize=12, color="red", verticalalignment="top")
    else:
        ax.text(0.05, 0.95, f"k = {k:.2f}", transform=ax.transAxes,
            fontsize=12, color="blue", verticalalignment="top")
    for i in range(num_of_subs):
        ax.set_title('reaction kinetics simulation', fontsize = 14)


    def init():
        for line in lines:
            line.set_data([], [])
        return lines

    def update(frames):
        t_frame = t[:frames]
        for i in range(num_of_subs):
            lines[i].set_data(t_frame, C[:frames, i])
        return lines

    draw = ani.FuncAnimation(fig , update, frames = len(t), init_func = init, interval = 20, blit = False)

    plt.tight_layout()
    plt.show()
    for i in range(len(species)):
        concentration = C[-1, i]
        if (concentration >= 1000 or concentration <= 0.01) and concentration != 0:
            print(f'[{species[i]}] = {concentration:.3e}')
        else:
            print(f'[{species[i]}] = {concentration:.3f}')
    print("conversion:")
    for i in range(len(reactant)):
        if C[0, i] != 0:
            conversion =(C[0, i] - C[-1, i]) / C[0, i] * 100
            print(f'{reactant[i]}: {conversion:.2f}%')

elif num >=5 and num <=6:
    reactant, intermediate, product, species, k, A0, time = info_of_series_reaction(num, T, time)
    for i in range(len(k)):
        print(f'k{i+1} =', f'{k[i]:.3e}')


    num_of_subs = len(reactant) + len(product) + len(intermediate)
    t, C = concentration_update_series(intermediate, product, k, A0, time)
    fig, ax = plt.subplots(figsize=(8, 6))

    species = list(reactant) + list(intermediate) + list(product)  
    lines = []
    line_styles = ['-', '--', '-.', ':']
    for i in range(num_of_subs):
        line, = ax.plot([], [], lw=2, label=species[i], linestyle=line_styles[i % 4])
        lines.append(line)

    ax.set_xlim(0, time)
    ax.set_ylim(0, max(A0)*1.2)
    ax.set_xlabel('time (s)', fontsize=12)
    ax.set_ylabel('concentration (M)', fontsize=12)
    ax.grid(True, linestyle='--')
    ax.legend()

    
    str_reactants = format_series_reaction_str(reactant)
    str_intermediate = format_series_reaction_str(intermediate)
    str_products = format_series_reaction_str(product)
    fig.suptitle(f"{str_reactants} --> {str_intermediate} --> {str_products}", fontsize=14, fontweight='bold')
    ax.set_title('Reaction Kinetics Simulation', fontsize=12)
    if np.any(k >= 100) :
        ax.text(0.05, 0.95, f"k1 = {k[0]:.2e}, k2 = {k[1]:.2e} \nWarning:k is big, the reaction may finish in a short time.", transform=ax.transAxes,
            fontsize=12, color="red", verticalalignment="top")
    elif np.any(k <= 0.001) :
        ax.text(0.05, 0.95, f"k1 = {k[0]:.3e}, k2 = {k[1]:.3e} \nWarning:k is small, concentration may change very slowly.", transform=ax.transAxes,
            fontsize=12, color="red", verticalalignment="top")
    else:
        ax.text(0.05, 0.95, f"k1 = {k[0]:.2f}, k2 = {k[1]:.2f}", transform=ax.transAxes,
            fontsize=12, color="blue", verticalalignment="top")
    for i in range(num_of_subs):
        ax.set_title('reaction kinetics simulation', fontsize = 14)


    def init():
        for line in lines:
            line.set_data([], [])
        return lines

    def update(frames):
        t_frame = t[:frames]
        for i in range(num_of_subs):
            lines[i].set_data(t_frame, C[:frames, i])
        return lines

    draw = ani.FuncAnimation(fig , update, frames = len(t), init_func = init, interval = 20, blit = False)

    plt.tight_layout()
    plt.show()

    print('Final concentration :')
    for i in range(len(species)):
        concentration = C[-1, i]
        if (concentration >= 1000 or concentration <= 0.01) and concentration != 0:
            print(f'[{species[i]}] = {concentration:.3e}')
        else:
            print(f'[{species[i]}] = {concentration:.3f}')

    max_idx = np.argmax(C[:, 1])
    t_max = t[max_idx]
    max_conc = C[max_idx, 1]
    print(f'The concentration of {species[1]} has maximum {max_conc:.3e} M at {t_max:.3f} second.')
else:
    print("Can't find the reaction")
