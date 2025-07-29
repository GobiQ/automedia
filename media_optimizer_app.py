#!/usr/bin/env python3
"""
Tissue Culture Media Optimizer - Streamlit Web App
Search for optimal gram-per-litre macro-salt recipes that satisfy
user-defined element ranges and nutrient-ratio windows.

Dependencies: streamlit, numpy, scipy, pandas
"""

import streamlit as st
import numpy as np
import pandas as pd
import random

# ─────────────────────────── Stoichiometry Database
STOICH_DATABASE = {
    # Common tissue culture salts
    'Ca(NO3)2·4H2O': {'Ca': 169.7, 'N_NO3': 118.8, 'category': 'Primary'},
    'KNO3': {'K': 386.7, 'N_NO3': 138.6, 'category': 'Primary'},
    'NH4NO3': {'N_NO3': 175.0, 'N_NH4': 175.0, 'category': 'Primary'},
    'NH4H2PO4': {'N_NH4': 121.7, 'P': 269.4, 'category': 'Primary'},
    'KH2PO4': {'K': 287.3, 'P': 227.7, 'category': 'Primary'},
    'K2SO4': {'K': 449.0, 'S': 184.0, 'category': 'Primary'},
    'MgSO4·7H2O': {'Mg': 98.6, 'S': 130.2, 'category': 'Primary'},
    'CaSO4·2H2O': {'Ca': 232.9, 'S': 186.3, 'category': 'Primary'},
    
    # Alternative salts
    'CaCl2·2H2O': {'Ca': 272.6, 'Cl': 482.9, 'category': 'Alternative'},
    'MgCl2·6H2O': {'Mg': 119.5, 'Cl': 348.8, 'category': 'Alternative'},
    'K2HPO4': {'K': 449.1, 'P': 178.2, 'category': 'Alternative'},
    'NaH2PO4·H2O': {'Na': 166.7, 'P': 224.6, 'category': 'Alternative'},
    '(NH4)2SO4': {'N_NH4': 212.1, 'S': 242.7, 'category': 'Alternative'},
    'NaNO3': {'Na': 270.6, 'N_NO3': 164.7, 'category': 'Alternative'},
    
    # Micronutrient salts
    'H3BO3': {'B': 174.8, 'category': 'Micronutrient'},
    'MnSO4·H2O': {'Mn': 363.0, 'S': 188.0, 'category': 'Micronutrient'},
    'ZnSO4·7H2O': {'Zn': 227.8, 'S': 111.5, 'category': 'Micronutrient'},
    'CuSO4·5H2O': {'Cu': 254.5, 'S': 128.3, 'category': 'Micronutrient'},
    'Na2MoO4·2H2O': {'Mo': 395.9, 'Na': 126.4, 'category': 'Micronutrient'},
    'FeSO4·7H2O': {'Fe': 201.5, 'S': 115.4, 'category': 'Micronutrient'},
    'Na2EDTA·2H2O': {'Na': 84.7, 'EDTA': 372.2, 'category': 'Micronutrient'},
    'MnCl2': {'Mn': 278.5, 'Cl': 504.8, 'category': 'Micronutrient'},
    'CuCl2': {'Cu': 370.5, 'Cl': 528.9, 'category': 'Micronutrient'},
}

# Default presets
DEFAULT_PRESETS = {
    "MS Medium": {
        'elements': {
            'N': (700, 1100), 'K': (700, 1000), 'P': (35, 60),
            'Ca': (200, 350), 'Mg': (80, 130), 'S': (60, 150)
        },
        'ratios': {
            'N:K': (0.9, 1.2), 'Ca:Mg': (2.0, 4.0), 'P:K': (0.05, 0.067)
        },
        'salts': ['Ca(NO3)2·4H2O', 'KNO3', 'NH4NO3', 'NH4H2PO4', 'KH2PO4', 'K2SO4', 'MgSO4·7H2O']
    },
    "Low Salt": {
        'elements': {
            'N': (400, 600), 'K': (400, 600), 'P': (20, 40),
            'Ca': (100, 200), 'Mg': (40, 80), 'S': (30, 80)
        },
        'ratios': {
            'N:K': (0.8, 1.3), 'Ca:Mg': (1.5, 3.5), 'P:K': (0.04, 0.08)
        },
        'salts': ['Ca(NO3)2·4H2O', 'KNO3', 'NH4H2PO4', 'MgSO4·7H2O']
    },
    "High Calcium": {
        'elements': {
            'N': (600, 900), 'K': (500, 800), 'P': (30, 50),
            'Ca': (300, 500), 'Mg': (60, 100), 'S': (50, 120)
        },
        'ratios': {
            'N:K': (1.0, 1.4), 'Ca:Mg': (3.0, 6.0), 'P:K': (0.04, 0.07)
        },
        'salts': ['Ca(NO3)2·4H2O', 'KNO3', 'NH4H2PO4', 'K2SO4', 'MgSO4·7H2O', 'CaSO4·2H2O']
    }
}

# ─────────────────────────── Helper Functions
def elemental_totals(g, selected_salts):
    """Calculate elemental totals from salt concentrations"""
    e = {k: 0 for k in ['N_NO3', 'N_NH4', 'N', 'K', 'Ca', 'Mg', 'P', 'S', 'Na', 'Cl', 'B', 'Mn', 'Zn', 'Cu', 'Mo', 'Fe', 'EDTA']}
    
    for grams, salt in zip(g, selected_salts):
        if salt in STOICH_DATABASE:
            for ion, mg_per_g in STOICH_DATABASE[salt].items():
                if ion in e:
                    e[ion] += grams * mg_per_g
    
    e['N'] = e['N_NO3'] + e['N_NH4']
    return e

def calculate_ratios(e):
    """Calculate key elemental ratios"""
    ratios = {}
    
    if e['K'] > 0:
        ratios['N:K'] = e['N'] / e['K']
        ratios['P:K'] = e['P'] / e['K']
    
    if e['Mg'] > 0:
        ratios['Ca:Mg'] = e['Ca'] / e['Mg']
    
    return ratios

def penalty_function(g, selected_salts, elem_bounds, ratio_bounds):
    """Penalty function for optimization"""
    e = elemental_totals(g, selected_salts)
    r = calculate_ratios(e)
    
    penalty = 0
    
    # Element constraints
    for el, (lo, hi) in elem_bounds.items():
        if el in e:
            if e[el] < lo:
                penalty += 1e6 + (lo - e[el]) ** 2
            elif e[el] > hi:
                penalty += 1e6 + (e[el] - hi) ** 2
    
    # Ratio constraints
    for ratio_name, (lo, hi) in ratio_bounds.items():
        if ratio_name in r:
            if r[ratio_name] < lo:
                penalty += 1e6 + (lo - r[ratio_name]) ** 2
            elif r[ratio_name] > hi:
                penalty += 1e6 + (r[ratio_name] - hi) ** 2
    
    # Minimize total salt concentration if feasible
    if penalty < 1e5:
        penalty += sum(g) * 0.01
    
    return penalty

def generate_salt_bounds(selected_salts):
    """Generate reasonable bounds for salt concentrations"""
    bounds_dict = {
        'Ca(NO3)2·4H2O': (0, 3.0),
        'KNO3': (0, 3.0),
        'NH4NO3': (0, 2.0),
        'NH4H2PO4': (0, 0.5),
        'KH2PO4': (0, 0.5),
        'K2SO4': (0, 1.0),
        'MgSO4·7H2O': (0, 1.5),
        'CaSO4·2H2O': (0, 0.5),
        'CaCl2·2H2O': (0, 2.0),
        'MgCl2·6H2O': (0, 1.0),
        'K2HPO4': (0, 0.5),
        'NaH2PO4·H2O': (0, 0.5),
        '(NH4)2SO4': (0, 1.0),
        'NaNO3': (0, 2.0),
        # Micronutrient bounds (much smaller concentrations)
        'H3BO3': (0, 0.01),
        'MnSO4·H2O': (0, 0.01),
        'ZnSO4·7H2O': (0, 0.01),
        'CuSO4·5H2O': (0, 0.01),
        'Na2MoO4·2H2O': (0, 0.001),
        'FeSO4·7H2O': (0, 0.01),
        'Na2EDTA·2H2O': (0, 0.01),
        'MnCl2': (0, 0.01),
        'CuCl2': (0, 0.01),
    }
    
    return [bounds_dict.get(salt, (0, 1.0)) for salt in selected_salts]

def differential_evolution_optimizer(objective_func, bounds, args, maxiter=1000, popsize=30, seed=42):
    """Custom implementation of differential evolution optimization"""
    random.seed(seed)
    np.random.seed(seed)
    
    # Initialize population
    population = []
    for _ in range(popsize):
        individual = np.array([random.uniform(lo, hi) for lo, hi in bounds])
        population.append(individual)
    
    best_individual = None
    best_fitness = np.inf
    
    # Evolution loop
    for generation in range(maxiter):
        for i in range(popsize):
            # Select three random individuals
            candidates = [j for j in range(popsize) if j != i]
            a, b, c = random.sample(candidates, 3)
            
            # Differential mutation
            F = 0.8  # Differential weight
            mutant = population[a] + F * (population[b] - population[c])
            
            # Crossover
            CR = 0.9  # Crossover probability
            trial = population[i].copy()
            for j in range(len(trial)):
                if random.random() < CR:
                    trial[j] = mutant[j]
            
            # Ensure bounds
            for j, (lo, hi) in enumerate(bounds):
                trial[j] = max(lo, min(hi, trial[j]))
            
            # Selection
            trial_fitness = objective_func(trial, *args)
            if trial_fitness < objective_func(population[i], *args):
                population[i] = trial
                if trial_fitness < best_fitness:
                    best_fitness = trial_fitness
                    best_individual = trial.copy()
    
    # Create result object similar to scipy
    class Result:
        def __init__(self, x, fun):
            self.x = x
            self.fun = fun
    
    return Result(best_individual, best_fitness)

def gradient_descent_optimizer(objective_func, x0, args, bounds, maxiter=1000):
    """Custom implementation of gradient descent optimization"""
    x = np.array(x0)
    learning_rate = 0.01
    
    for iteration in range(maxiter):
        # Simple finite difference gradient
        grad = np.zeros_like(x)
        h = 1e-6
        
        for i in range(len(x)):
            x_plus = x.copy()
            x_plus[i] += h
            x_minus = x.copy()
            x_minus[i] -= h
            
            grad[i] = (objective_func(x_plus, *args) - objective_func(x_minus, *args)) / (2 * h)
        
        # Update
        x_new = x - learning_rate * grad
        
        # Apply bounds
        for i, (lo, hi) in enumerate(bounds):
            x_new[i] = max(lo, min(hi, x_new[i]))
        
        # Check convergence
        if np.linalg.norm(x_new - x) < 1e-6:
            break
        
        x = x_new
    
    # Create result object similar to scipy
    class Result:
        def __init__(self, x, fun):
            self.x = x
            self.fun = fun
    
    return Result(x, objective_func(x, *args))

def optimize_media(selected_salts, elem_bounds, ratio_bounds, algorithm='DE', n_trials=1):
    """Optimize media composition using specified algorithm"""
    bounds = generate_salt_bounds(selected_salts)
    
    best_result = None
    best_penalty = np.inf
    
    for trial in range(n_trials):
        # Seeded Monte Carlo pre-search
        random.seed(42 + trial)
        np.random.seed(42 + trial)
        seed_pool = []
        
        for _ in range(5000):
            guess = np.array([random.uniform(lo, hi) for lo, hi in bounds])
            pen = penalty_function(guess, selected_salts, elem_bounds, ratio_bounds)
            if pen < 1e5:
                seed_pool.append(guess)
        
        if algorithm == 'DE':
            # Custom Differential Evolution
            result = differential_evolution_optimizer(
                penalty_function,
                bounds,
                args=(selected_salts, elem_bounds, ratio_bounds),
                maxiter=1000,
                popsize=30,
                seed=42 + trial
            )
        
        else:  # Gradient Descent (replaces SLSQP)
            if seed_pool:
                x0 = seed_pool[0]
            else:
                x0 = np.array([(lo + hi) / 2 for lo, hi in bounds])
            
            result = gradient_descent_optimizer(
                penalty_function,
                x0,
                args=(selected_salts, elem_bounds, ratio_bounds),
                bounds=bounds,
                maxiter=1000
            )
        
        if hasattr(result, 'fun') and result.fun < best_penalty:
            best_penalty = result.fun
            best_result = result
    
    return best_result

# ─────────────────────────── Streamlit App
def main():
    st.set_page_config(
        page_title="Tissue Culture Media Optimizer",
        page_icon="🧪",
        layout="wide"
    )
    
    st.title("🧪 Tissue Culture Media Optimizer")
    st.markdown("Optimize macro-salt recipes for tissue culture media using Monte Carlo seeding and evolutionary algorithms.")
    
    # Sidebar for main controls
    st.sidebar.header("Optimization Settings")
    
    # Preset selection
    preset_name = st.sidebar.selectbox(
        "Load Preset",
        ["Custom"] + list(DEFAULT_PRESETS.keys())
    )
    
    # Algorithm selection
    algorithm = st.sidebar.selectbox(
        "Optimization Algorithm",
        ["DE", "GD"],
        help="DE: Differential Evolution (global), GD: Gradient Descent (local)"
    )
    
    n_trials = st.sidebar.slider("Number of optimization trials", 1, 5, 1)
    
    # Main content area
    col1, col2 = st.columns([1, 1])
    
    with col1:
        st.header("Salt Selection")
        
        # Available salts grouped by category
        primary_salts = [salt for salt, data in STOICH_DATABASE.items() if data['category'] == 'Primary']
        alternative_salts = [salt for salt, data in STOICH_DATABASE.items() if data['category'] == 'Alternative']
        
        st.subheader("Primary Salts")
        selected_primary = []
        for salt in primary_salts:
            if preset_name != "Custom":
                default_val = salt in DEFAULT_PRESETS[preset_name]['salts']
            else:
                default_val = True  # Default all primary salts to checked
            
            if st.checkbox(salt, value=default_val, key=f"primary_{salt}"):
                selected_primary.append(salt)
        
        st.subheader("Alternative Salts")
        selected_alternative = []
        for salt in alternative_salts:
            if preset_name != "Custom":
                default_val = salt in DEFAULT_PRESETS[preset_name]['salts']
            else:
                default_val = True  # Default all alternative salts to checked
            
            if st.checkbox(salt, value=default_val, key=f"alt_{salt}"):
                selected_alternative.append(salt)
        
        st.subheader("Micronutrient Salts")
        micronutrient_salts = [salt for salt, data in STOICH_DATABASE.items() if data['category'] == 'Micronutrient']
        selected_micronutrients = []
        for salt in micronutrient_salts:
            if preset_name != "Custom":
                default_val = salt in DEFAULT_PRESETS[preset_name]['salts']
            else:
                default_val = True  # Default all micronutrients to checked
            
            if st.checkbox(salt, value=default_val, key=f"micronutrient_{salt}"):
                selected_micronutrients.append(salt)
        
        selected_salts = selected_primary + selected_alternative + selected_micronutrients
        
        if not selected_salts:
            st.warning("Please select at least one salt.")
            return
    
    with col2:
        st.header("Target Ranges")
        
        # Load preset values if selected
        if preset_name != "Custom":
            preset_data = DEFAULT_PRESETS[preset_name]
            default_elements = preset_data['elements']
            default_ratios = preset_data['ratios']
        else:
            default_elements = {
                'N': (700, 1100), 'K': (700, 1000), 'P': (35, 60),
                'Ca': (200, 350), 'Mg': (80, 130), 'S': (60, 150)
            }
            default_ratios = {
                'N:K': (0.9, 1.2), 'Ca:Mg': (2.0, 4.0), 'P:K': (0.05, 0.067)
            }
        
        st.subheader("Macronutrient Concentrations (mg/L)")
        elem_bounds = {}
        
        for element in ['N', 'K', 'P', 'Ca', 'Mg', 'S']:
            if element in default_elements:
                min_val, max_val = default_elements[element]
            else:
                min_val, max_val = 0, 1000
            
            col_min, col_max = st.columns(2)
            with col_min:
                min_bound = st.number_input(f"{element} min", value=float(min_val), key=f"{element}_min")
            with col_max:
                max_bound = st.number_input(f"{element} max", value=float(max_val), key=f"{element}_max")
            
            elem_bounds[element] = (min_bound, max_bound)
        
        st.subheader("Micronutrient Concentrations (mg/L)")
        micronutrient_defaults = {
            'B': (0.5, 3.0), 'Mn': (2.0, 10.0), 'Zn': (0.5, 2.0),
            'Cu': (0.01, 0.1), 'Mo': (0.01, 0.1), 'Fe': (2.0, 10.0)
        }
        
        for element in ['B', 'Mn', 'Zn', 'Cu', 'Mo', 'Fe']:
            if element in micronutrient_defaults:
                min_val, max_val = micronutrient_defaults[element]
            else:
                min_val, max_val = 0, 10
            
            col_min, col_max = st.columns(2)
            with col_min:
                min_bound = st.number_input(f"{element} min", value=float(min_val), key=f"{element}_min")
            with col_max:
                max_bound = st.number_input(f"{element} max", value=float(max_val), key=f"{element}_max")
            
            elem_bounds[element] = (min_bound, max_bound)
        
        st.subheader("Element Ratios")
        ratio_bounds = {}
        
        for ratio in ['N:K', 'Ca:Mg', 'P:K']:
            if ratio in default_ratios:
                min_val, max_val = default_ratios[ratio]
            else:
                min_val, max_val = 0.5, 2.0
            
            col_min, col_max = st.columns(2)
            with col_min:
                min_bound = st.number_input(f"{ratio} min", value=float(min_val), key=f"{ratio}_min", format="%.3f")
            with col_max:
                max_bound = st.number_input(f"{ratio} max", value=float(max_val), key=f"{ratio}_max", format="%.3f")
            
            ratio_bounds[ratio] = (min_bound, max_bound)
    
    # Optimization button and results
    st.header("Optimization Results")
    
    if st.button("🚀 Optimize Media Recipe", type="primary"):
        with st.spinner("Optimizing... This may take a few moments."):
            try:
                result = optimize_media(selected_salts, elem_bounds, ratio_bounds, algorithm, n_trials)
                
                if result is not None and hasattr(result, 'x'):
                    g_opt = np.round(result.x, 3)
                    e_opt = elemental_totals(g_opt, selected_salts)
                    r_opt = calculate_ratios(e_opt)
                    
                    # Display results in columns
                    col1, col2, col3 = st.columns(3)
                    
                    with col1:
                        st.subheader("📋 Recipe (g/L)")
                        recipe_data = []
                        for salt, grams in zip(selected_salts, g_opt):
                            if grams > 0.001:  # Only show salts with meaningful concentrations
                                recipe_data.append({'Salt': salt, 'Concentration (g/L)': f"{grams:.3f}"})
                        
                        if recipe_data:
                            recipe_df = pd.DataFrame(recipe_data)
                            st.dataframe(recipe_df, hide_index=True)
                        else:
                            st.warning("No viable recipe found with current constraints.")
                    
                    with col2:
                        st.subheader("🧮 Elemental Totals (mg/L)")
                        element_data = []
                        # Macronutrients
                        for element in ['N', 'K', 'P', 'Ca', 'Mg', 'S']:
                            if element in e_opt:
                                target_min, target_max = elem_bounds[element]
                                actual = e_opt[element]
                                status = "✅" if target_min <= actual <= target_max else "❌"
                                element_data.append({
                                    'Element': element,
                                    'Actual': f"{actual:.1f}",
                                    'Target': f"{target_min:.0f}-{target_max:.0f}",
                                    'Status': status
                                })
                        
                        # Micronutrients
                        for element in ['B', 'Mn', 'Zn', 'Cu', 'Mo', 'Fe']:
                            if element in e_opt:
                                target_min, target_max = elem_bounds[element]
                                actual = e_opt[element]
                                status = "✅" if target_min <= actual <= target_max else "❌"
                                element_data.append({
                                    'Element': element,
                                    'Actual': f"{actual:.3f}",
                                    'Target': f"{target_min:.3f}-{target_max:.3f}",
                                    'Status': status
                                })
                        
                        element_df = pd.DataFrame(element_data)
                        st.dataframe(element_df, hide_index=True)
                    
                    with col3:
                        st.subheader("📊 Key Ratios")
                        ratio_data = []
                        for ratio_name in ['N:K', 'Ca:Mg', 'P:K']:
                            if ratio_name in r_opt and ratio_name in ratio_bounds:
                                target_min, target_max = ratio_bounds[ratio_name]
                                actual = r_opt[ratio_name]
                                status = "✅" if target_min <= actual <= target_max else "❌"
                                ratio_data.append({
                                    'Ratio': ratio_name,
                                    'Actual': f"{actual:.3f}",
                                    'Target': f"{target_min:.3f}-{target_max:.3f}",
                                    'Status': status
                                })
                        
                        ratio_df = pd.DataFrame(ratio_data)
                        st.dataframe(ratio_df, hide_index=True)
                    
                    # Optimization info
                    if hasattr(result, 'fun'):
                        feasible = result.fun < 1e5
                        st.info(f"**Optimization Status:** {'✅ Feasible' if feasible else '❌ Infeasible'} | "
                               f"**Final Penalty:** {result.fun:.2e} | **Algorithm:** {algorithm}")
                    
                    # Visualization
                    if len([g for g in g_opt if g > 0.001]) > 0:
                        st.subheader("📈 Recipe Visualization")
                        
                        # Recipe bar chart using Streamlit
                        nonzero_salts = [(salt, grams) for salt, grams in zip(selected_salts, g_opt) if grams > 0.001]
                        if nonzero_salts:
                            salt_names, concentrations = zip(*nonzero_salts)
                            
                            # Create a simple bar chart using Streamlit
                            chart_data = pd.DataFrame({
                                'Salt': salt_names,
                                'Concentration (g/L)': concentrations
                            })
                            
                            st.bar_chart(chart_data.set_index('Salt'))
                            
                            # Also show as a table for precise values
                            st.subheader("📊 Detailed Recipe Breakdown")
                            st.dataframe(chart_data, hide_index=True)
                
                else:
                    st.error("Optimization failed. Try adjusting constraints or selecting different salts.")
            
            except Exception as e:
                st.error(f"An error occurred during optimization: {str(e)}")
    
    # Information section
    with st.expander("ℹ️ How to Use This Tool"):
        st.markdown("""
        1. **Select a preset** or use custom settings
        2. **Choose salts** from the available options (primary salts are commonly used)
        3. **Set target ranges** for elements (mg/L) and ratios
        4. **Choose algorithm**: DE for global optimization, SLSQP for local optimization
        5. **Click optimize** to find the best recipe
        
        **Tips:**
        - Start with fewer salts for simpler optimization
        - DE algorithm is more robust but slower
        - Increase trials for better solutions
        - Check that constraints are realistic and not conflicting
        """)
    
    with st.expander("🧪 Salt Information"):
        # Display stoichiometry table
        stoich_data = []
        for salt, data in STOICH_DATABASE.items():
            elements = ", ".join([f"{elem}: {conc:.1f}" for elem, conc in data.items() if elem != 'category'])
            stoich_data.append({
                'Salt': salt,
                'Category': data['category'],
                'Elements (mg/g salt)': elements
            })
        
        stoich_df = pd.DataFrame(stoich_data)
        st.dataframe(stoich_df, hide_index=True)

if __name__ == "__main__":
    main()
