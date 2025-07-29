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
    'MnCl2': {'Mn': 436.7, 'Cl': 563.3, 'category': 'Micronutrient'},  # Anhydrous
    'CuCl2': {'Cu': 472.8, 'Cl': 527.2, 'category': 'Micronutrient'},  # Anhydrous
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

def sq_violation(x, lo, hi):
    """Scaled squared violation for smoother penalty landscape"""
    if x < lo:
        return ((lo - x) / (hi - lo)) ** 2
    if x > hi:
        return ((x - hi) / (hi - lo)) ** 2
    return 0.0

def penalty_function(g, selected_salts, elem_bounds, ratio_bounds):
    """Penalty function for optimization with smoother landscape"""
    e = elemental_totals(g, selected_salts)
    r = calculate_ratios(e)
    
    penalty = 0
    
    # Check if micronutrient salts are being used
    micronutrient_salts_used = False
    micronutrient_salt_indices = [i for i, salt in enumerate(selected_salts) 
                                 if STOICH_DATABASE[salt]['category'] == 'Micronutrient']
    
    if micronutrient_salt_indices:
        micronutrient_concentrations = [g[i] for i in micronutrient_salt_indices]
        if any(conc > 1e-6 for conc in micronutrient_concentrations):
            micronutrient_salts_used = True
    
    # Element constraints with scaled violations
    for el, (lo, hi) in elem_bounds.items():
        if el in e:
            penalty += sq_violation(e[el], lo, hi) * 1000  # Weighted penalty
    
    # Ratio constraints with scaled violations
    for ratio_name, (lo, hi) in ratio_bounds.items():
        if ratio_name in r:
            penalty += sq_violation(r[ratio_name], lo, hi) * 1000  # Weighted penalty
    
    # Micronutrient inclusion penalty (softer)
    if not micronutrient_salts_used:
        penalty += 1000  # High but not impossible penalty
    
    # Additional penalty for micronutrient violations
    micronutrients = ['Cu', 'Mo', 'B', 'Mn', 'Zn', 'Fe']
    for micronutrient in micronutrients:
        if micronutrient in elem_bounds and micronutrient in e:
            target_min, target_max = elem_bounds[micronutrient]
            actual = e[micronutrient]
            penalty += sq_violation(actual, target_min, target_max) * 2000  # Higher weight for micronutrients
    
    # Minimize total salt concentration if feasible
    if penalty < 100:  # Much lower threshold
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
        # Micronutrient bounds (calculated for target ranges)
        'H3BO3': (0, 0.02),  # For 0.5-3.0 mg/L B
        'MnSO4·H2O': (0, 0.03),  # For 2.0-10.0 mg/L Mn
        'ZnSO4·7H2O': (0, 0.01),  # For 0.5-2.0 mg/L Zn
        'CuSO4·5H2O': (0, 0.001),  # For 0.01-0.1 mg/L Cu
        'Na2MoO4·2H2O': (0, 0.001),  # For 0.01-0.1 mg/L Mo
        'FeSO4·7H2O': (0, 0.05),  # For 2.0-10.0 mg/L Fe
        'Na2EDTA·2H2O': (0, 0.03),  # For Fe chelation
        'MnCl2': (0, 0.04),  # For 2.0-10.0 mg/L Mn
        'CuCl2': (0, 0.001),  # For 0.01-0.1 mg/L Cu
    }
    
    bounds = []
    for salt in selected_salts:
        if salt in bounds_dict:
            lo, hi = bounds_dict[salt]
            # Force minimum concentration for micronutrient salts
            if STOICH_DATABASE[salt]['category'] == 'Micronutrient':
                # Set minimum to a small but non-zero value
                lo = 0.0001  # 0.1 mg/L minimum
            bounds.append((lo, hi))
        else:
            bounds.append((0, 1.0))
    
    return bounds

def differential_evolution_optimizer(objective_func, bounds, args, maxiter=1000, popsize=30, seed=42, seed_pool=None):
    """Custom implementation of differential evolution optimization"""
    random.seed(seed)
    np.random.seed(seed)
    
    # Initialize population with seed pool if available
    population = []
    
    # Use seed pool for 30% of population if available
    if seed_pool and len(seed_pool) > 0:
        n_seeds = min(len(seed_pool), int(popsize * 0.3))
        for i in range(n_seeds):
            seed_individual = seed_pool[i % len(seed_pool)].copy()
            # Ensure bounds
            for j, (lo, hi) in enumerate(bounds):
                seed_individual[j] = max(lo, min(hi, seed_individual[j]))
            population.append(seed_individual)
    
    # Fill rest with random individuals
    while len(population) < popsize:
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
            
            # Differential mutation with adaptive F
            F = 0.8 if generation < maxiter // 2 else 0.5  # Adaptive differential weight
            mutant = population[a] + F * (population[b] - population[c])
            
            # Crossover with adaptive CR
            CR = 0.9 if generation < maxiter // 2 else 0.7  # Adaptive crossover probability
            trial = population[i].copy()
            for j in range(len(trial)):
                if random.random() < CR:
                    trial[j] = mutant[j]
            
            # Ensure bounds
            for j, (lo, hi) in enumerate(bounds):
                trial[j] = max(lo, min(hi, trial[j]))
                # Ensure micronutrient salts never go to zero
                if lo > 0:  # This is a micronutrient salt
                    trial[j] = max(lo, trial[j])  # Force minimum concentration
            
            # Selection
            trial_fitness = objective_func(trial, *args)
            current_fitness = objective_func(population[i], *args)
            if trial_fitness < current_fitness:
                population[i] = trial
                if trial_fitness < best_fitness:
                    best_fitness = trial_fitness
                    best_individual = trial.copy()
        
        # Periodically inject diversity for micronutrients
        if generation % 100 == 0 and generation > 0:
            # Force some individuals to include micronutrients
            for i in range(min(5, popsize // 2)):
                individual = population[i]
                # Ensure micronutrient salts have some concentration
                for j, (lo, hi) in enumerate(bounds):
                    if lo > 0:  # This is a micronutrient salt
                        individual[j] = random.uniform(lo, hi)
                population[i] = individual
    
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
        
        # Pre-seed with micronutrient salts
        micronutrient_indices = [i for i, salt in enumerate(selected_salts) 
                               if STOICH_DATABASE[salt]['category'] == 'Micronutrient']
        
        # Create multiple seeds with micronutrients included
        if micronutrient_indices:
            for seed_trial in range(10):  # Create 10 different seeds
                seed_with_micronutrients = np.array([(lo + hi) / 2 for lo, hi in bounds])
                # Set micronutrient salts to reasonable values
                for idx in micronutrient_indices:
                    lo, hi = bounds[idx]
                    # Use different values for each seed
                    seed_with_micronutrients[idx] = lo + (hi - lo) * (seed_trial / 10)
                pen = penalty_function(seed_with_micronutrients, selected_salts, elem_bounds, ratio_bounds)
                if pen < 1e8:  # Allow higher penalty for initial seeds
                    seed_pool.append(seed_with_micronutrients)
        
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
                seed=42 + trial,
                seed_pool=seed_pool
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
                    g_best = result.x.copy()  # Use unrounded for calculations
                    e_opt = elemental_totals(g_best, selected_salts)
                    r_opt = calculate_ratios(e_opt)
                    
                    # Display results in columns
                    col1, col2, col3 = st.columns(3)
                    
                    with col1:
                        st.subheader("📋 Recipe (g/L)")
                        recipe_data = []
                        g_disp = np.round(g_best, 6)  # Round only for display
                        for salt, grams in zip(selected_salts, g_disp):
                            if grams > 1e-6:  # Show micronutrient salts with very small concentrations
                                recipe_data.append({'Salt': salt, 'Concentration (g/L)': f"{grams:.6f}"})
                        
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
                            if element in e_opt and element in elem_bounds:
                                target_min, target_max = elem_bounds[element]
                                actual = e_opt[element]
                                status = "✅" if target_min <= actual <= target_max else "❌"
                                element_data.append({
                                    'Element': element,
                                    'Actual': f"{actual:.3f}",
                                    'Target': f"{target_min:.3f}-{target_max:.3f}",
                                    'Status': status
                                })
                            elif element in e_opt:
                                # Element calculated but not in bounds (shouldn't happen)
                                actual = e_opt[element]
                                element_data.append({
                                    'Element': element,
                                    'Actual': f"{actual:.3f}",
                                    'Target': "Not specified",
                                    'Status': "⚠️"
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
                    
                    # Debug information
                    with st.expander("🔍 Debug Information"):
                        st.write("**Elements being optimized:**")
                        st.write(list(elem_bounds.keys()))
                        st.write("**Elements calculated:**")
                        st.write(list(e_opt.keys()))
                        st.write("**Selected salts:**")
                        st.write(selected_salts)
                        st.write("**Salt concentrations (g/L):**")
                        for salt, conc in zip(selected_salts, g_best):
                            if conc > 1e-6:  # Show even very small concentrations
                                st.write(f"{salt}: {conc:.6f}")
                        
                        st.write("**Micronutrient Analysis:**")
                        for element in ['Cu', 'Mo', 'B', 'Mn', 'Zn', 'Fe']:
                            if element in e_opt and element in elem_bounds:
                                target_min, target_max = elem_bounds[element]
                                actual = e_opt[element]
                                st.write(f"{element}: {actual:.6f} mg/L (target: {target_min:.3f}-{target_max:.3f})")
                            elif element in e_opt:
                                actual = e_opt[element]
                                st.write(f"{element}: {actual:.6f} mg/L (no target set)")
                        
                        st.write("**Penalty breakdown:**")
                        e_debug = elemental_totals(g_best, selected_salts)
                        r_debug = calculate_ratios(e_debug)
                        penalty_debug = penalty_function(g_best, selected_salts, elem_bounds, ratio_bounds)
                        st.write(f"Total penalty: {penalty_debug:.2e}")
                        
                        # Show individual constraint violations
                        for el, (lo, hi) in elem_bounds.items():
                            if el in e_debug:
                                actual = e_debug[el]
                                if actual < lo:
                                    st.write(f"❌ {el}: {actual:.6f} < {lo:.6f}")
                                elif actual > hi:
                                    st.write(f"❌ {el}: {actual:.6f} > {hi:.6f}")
                                else:
                                    st.write(f"✅ {el}: {actual:.6f} in range [{lo:.6f}, {hi:.6f}]")
                        
                        # Check if micronutrient salts are in selected_salts
                        micronutrient_salts_in_selected = [salt for salt in selected_salts if STOICH_DATABASE[salt]['category'] == 'Micronutrient']
                        st.write("**Micronutrient salts in selection:**")
                        st.write(micronutrient_salts_in_selected)
                        
                        # Check if Cu and Mo salts are present
                        cu_salts = [salt for salt in selected_salts if 'Cu' in STOICH_DATABASE[salt]]
                        mo_salts = [salt for salt in selected_salts if 'Mo' in STOICH_DATABASE[salt]]
                        st.write("**Copper salts selected:**", cu_salts)
                        st.write("**Molybdenum salts selected:**", mo_salts)
                        
                        # Show what each salt contributes
                        st.write("**Salt contributions to Cu and Mo:**")
                        for salt, conc in zip(selected_salts, g_best):
                            if conc > 1e-6 and salt in STOICH_DATABASE:
                                salt_data = STOICH_DATABASE[salt]
                                if 'Cu' in salt_data:
                                    cu_contribution = conc * salt_data['Cu']
                                    st.write(f"{salt}: {cu_contribution:.6f} mg/L Cu")
                                if 'Mo' in salt_data:
                                    mo_contribution = conc * salt_data['Mo']
                                    st.write(f"{salt}: {mo_contribution:.6f} mg/L Mo")
                        
                        # Test penalty function with different scenarios
                        st.write("**Penalty Function Testing:**")
                        # Test with zero micronutrients
                        test_zero = np.zeros_like(g_best)
                        penalty_zero = penalty_function(test_zero, selected_salts, elem_bounds, ratio_bounds)
                        st.write(f"Penalty with zero micronutrients: {penalty_zero:.2e}")
                        
                        # Test with reasonable micronutrient values
                        test_micro = g_best.copy()
                        for i, salt in enumerate(selected_salts):
                            if STOICH_DATABASE[salt]['category'] == 'Micronutrient':
                                test_micro[i] = 0.001  # 1 mg/L equivalent
                        penalty_micro = penalty_function(test_micro, selected_salts, elem_bounds, ratio_bounds)
                        st.write(f"Penalty with micronutrients: {penalty_micro:.2e}")
                        
                        st.write(f"**Best solution penalty: {penalty_debug:.2e}**")
                        if penalty_debug > 1e6:
                            st.error("❌ Solution has high penalty - constraints not met!")
                        else:
                            st.success("✅ Solution has low penalty - constraints met!")
                    
                    # Visualization
                    if len([g for g in g_best if g > 1e-6]) > 0:
                        st.subheader("📈 Recipe Visualization")
                        
                        # Recipe bar chart using Streamlit
                        nonzero_salts = [(salt, grams) for salt, grams in zip(selected_salts, g_disp) if grams > 1e-6]
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
