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
import uuid

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
    
    # ADDED: Penalty for using multiple salts for the same micronutrient
    micronutrient_providers = {}  # element -> list of (salt_index, concentration)
    for i, (conc, salt) in enumerate(zip(g, selected_salts)):
        if conc > 1e-6 and salt in STOICH_DATABASE:
            salt_data = STOICH_DATABASE[salt]
            for element in ['Cu', 'Mo', 'B', 'Mn', 'Zn', 'Fe']:
                if element in salt_data and element in elem_bounds:
                    if element not in micronutrient_providers:
                        micronutrient_providers[element] = []
                    micronutrient_providers[element].append((i, conc))
    
    # Heavy penalty for using multiple salts for the same micronutrient
    for element, providers in micronutrient_providers.items():
        if len(providers) > 1:
            penalty += 50000  # Very high penalty for double-counting
    
    # Check if micronutrient salts are being used
    micronutrient_salts_used = False
    micronutrient_salt_indices = [i for i, salt in enumerate(selected_salts) 
                                 if STOICH_DATABASE[salt]['category'] == 'Micronutrient']
    
    if micronutrient_salt_indices:
        micronutrient_concentrations = [g[i] for i in micronutrient_salt_indices]
        if any(conc > 1e-6 for conc in micronutrient_concentrations):
            micronutrient_salts_used = True
    
    # Ratio constraints with scaled violations
    for ratio_name, (lo, hi) in ratio_bounds.items():
        if ratio_name in r:
            penalty += sq_violation(r[ratio_name], lo, hi) * 1000  # Weighted penalty
    
    # Micronutrient inclusion penalty (softer)
    if not micronutrient_salts_used:
        penalty += 1000  # High but not impossible penalty
    
    # FIXED: Check micronutrient coverage more thoroughly
    micronutrient_coverage = {'Cu': False, 'Mo': False, 'B': False, 'Mn': False, 'Zn': False, 'Fe': False}
    
    for element in micronutrient_coverage.keys():
        if element in elem_bounds:
            min_target, max_target = elem_bounds[element]
            actual = e.get(element, 0)
            if actual >= min_target:
                micronutrient_coverage[element] = True
    
    # Element constraints with scaled violations (FIXED: Only one calculation with proper Cu/Mo weighting)
    for el, (lo, hi) in elem_bounds.items():
        if el in e:
            violation = sq_violation(e[el], lo, hi)
            # FIXED: Higher penalty for Cu and Mo violations
            weight = 5000 if el in ['Cu', 'Mo'] else 1000
            penalty += violation * weight
    
    # FIXED: Micronutrient coverage penalty - especially for Cu and Mo
    uncovered_micronutrients = [el for el, covered in micronutrient_coverage.items() if not covered and el in elem_bounds]
    if uncovered_micronutrients:
        # Heavy penalty for missing Cu and Mo specifically
        cu_mo_missing = any(el in ['Cu', 'Mo'] for el in uncovered_micronutrients)
        base_penalty = 10000 if cu_mo_missing else 5000
        penalty += base_penalty * len(uncovered_micronutrients)
    
    # Minimize total salt concentration if feasible
    if penalty < 1000:  # Only if constraints are mostly satisfied
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
        'CuSO4·5H2O': (0.0001, 0.005),  # MINIMUM BOUND ADDED for Cu
        'Na2MoO4·2H2O': (0.0001, 0.005),  # MINIMUM BOUND ADDED for Mo
        'FeSO4·7H2O': (0, 0.05),  # For 2.0-10.0 mg/L Fe
        'Na2EDTA·2H2O': (0, 0.03),  # For Fe chelation
        'MnCl2': (0, 0.04),  # For 2.0-10.0 mg/L Mn
        'CuCl2': (0.0001, 0.005),  # MINIMUM BOUND ADDED for Cu
    }
    
    bounds = []
    for salt in selected_salts:
        if salt in bounds_dict:
            lo, hi = bounds_dict[salt]
            bounds.append((lo, hi))
        else:
            bounds.append((0, 1.0))
    
    return bounds

# Remove cache decorator completely to force fresh execution
def differential_evolution_optimizer(objective_func, bounds, args, maxiter=1000, popsize=30, seed=42, seed_pool=None):
    """Custom implementation of differential evolution optimization"""
    random.seed(seed)
    np.random.seed(seed)
    
    # Extract selected_salts from args for micronutrient identification
    selected_salts = args[0]  # First argument is selected_salts
    elem_bounds = args[1]     # Second argument is elem_bounds
    
    # FIXED: Pre-calculate optimal micronutrient concentrations with one-salt-per-element
    optimal_micronutrients = {}
    element_to_salt_map = {}  # Track which salt is assigned to each element
    
    for j, salt in enumerate(selected_salts):
        if STOICH_DATABASE[salt]['category'] == 'Micronutrient':
            # Find which element this salt provides (prioritize Cu, Mo first)
            for element in ['Cu', 'Mo', 'B', 'Mn', 'Zn', 'Fe']:
                if (element in elem_bounds and 
                    element in STOICH_DATABASE[salt] and
                    element not in element_to_salt_map):  # Only if element not already assigned
                    
                    min_val, max_val = elem_bounds[element]
                    # CHANGED: Target minimum + 20% instead of 40%
                    target = min_val + (max_val - min_val) * 0.2
                    mg_per_g = STOICH_DATABASE[salt][element]
                    optimal_g_per_l = target / mg_per_g
                    
                    # Ensure it fits within bounds
                    lo, hi = bounds[j]
                    if optimal_g_per_l <= hi and optimal_g_per_l >= lo:
                        optimal_micronutrients[j] = optimal_g_per_l
                        element_to_salt_map[element] = j  # Record which salt handles this element
                        break
    
    # CRITICAL: Zero out salts that provide elements already handled by other salts
    for j, salt in enumerate(selected_salts):
        if STOICH_DATABASE[salt]['category'] == 'Micronutrient' and j not in optimal_micronutrients:
            # This salt is not the primary provider for any element
            optimal_micronutrients[j] = 0.0  # Force to zero
    
    # Initialize population with seed pool if available
    population = []
    
    # Use seed pool for 50% of population if available (increased from 30%)
    if seed_pool and len(seed_pool) > 0:
        n_seeds = min(len(seed_pool), int(popsize * 0.5))
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
    
    # IMPROVED: Inject optimal micronutrients into ALL initial population
    for individual in population:
        for j, optimal_conc in optimal_micronutrients.items():
            individual[j] = optimal_conc
    
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
            
            # IMPROVED: Force micronutrients into trial solution
            for j, optimal_conc in optimal_micronutrients.items():
                if random.random() < 0.3:  # 30% chance to force optimal micronutrient
                    trial[j] = optimal_conc
            
            # FORCE: Ensure only one salt per element is used
            used_elements = set()
            for j, salt in enumerate(selected_salts):
                if salt in STOICH_DATABASE:
                    salt_data = STOICH_DATABASE[salt]
                    for element in ['Cu', 'Mo', 'Mn']:
                        if element in salt_data:
                            if element not in used_elements:
                                # First salt for this element - ensure minimum inclusion
                                bounds = generate_salt_bounds([salt])
                                lo, hi = bounds[0]
                                if trial[j] < lo:
                                    trial[j] = lo + (hi - lo) * 0.1
                                used_elements.add(element)
                            else:
                                # Second salt for this element - zero it out
                                trial[j] = 0.0
                            break  # Only process each salt once
            
            # Selection
            trial_fitness = objective_func(trial, *args)
            current_fitness = objective_func(population[i], *args)
            if trial_fitness < current_fitness:
                population[i] = trial
                if trial_fitness < best_fitness:
                    best_fitness = trial_fitness
                    best_individual = trial.copy()
        
        # IMPROVED: Inject optimal micronutrient concentrations more frequently
        if generation % 25 == 0:  # Every 25 generations instead of 100
            # Inject optimal micronutrient values into some individuals
            for i in range(min(10, popsize // 2)):  # More individuals
                individual = population[i]
                # Set micronutrient salts to their optimal values
                for j, optimal_conc in optimal_micronutrients.items():
                    individual[j] = optimal_conc
                
                # FIXED: Ensure only one salt per element is used
                used_elements = set()
                for j, salt in enumerate(selected_salts):
                    if salt in STOICH_DATABASE:
                        salt_data = STOICH_DATABASE[salt]
                        for element in ['Cu', 'Mo', 'Mn']:
                            if element in salt_data:
                                if element not in used_elements:
                                    used_elements.add(element)
                                else:
                                    # Second salt for this element - zero it out
                                    individual[j] = 0.0
                                break  # Only process each salt once
                
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

def calculate_micronutrient_seeds(selected_salts, elem_bounds):
    """Pre-calculate optimal micronutrient salt concentrations - FIXED for Cu/Mo"""
    micronutrient_seeds = {}
    
    # FIXED: Target the MINIMUM of each micronutrient range, not middle
    micronutrient_targets = {}
    for element in ['B', 'Mn', 'Zn', 'Cu', 'Mo', 'Fe']:
        if element in elem_bounds:
            min_val, max_val = elem_bounds[element]
            # CHANGED: Target minimum + 20% of range instead of middle
            micronutrient_targets[element] = min_val + (max_val - min_val) * 0.2
    
    # FIXED: Prioritize Cu and Mo salts first
    priority_elements = ['Cu', 'Mo', 'B', 'Mn', 'Zn', 'Fe']
    
    for element in priority_elements:
        if element in micronutrient_targets:
            target_concentration = micronutrient_targets[element]
            
            # Find the best salt for this element
            best_salt = None
            best_efficiency = 0
            
            for salt in selected_salts:
                if (salt in STOICH_DATABASE and 
                    STOICH_DATABASE[salt]['category'] == 'Micronutrient' and
                    element in STOICH_DATABASE[salt]):
                    
                    mg_per_g = STOICH_DATABASE[salt][element]
                    required_g_per_l = target_concentration / mg_per_g
                    
                    # Check if it fits within bounds
                    bounds = generate_salt_bounds([salt])
                    lo, hi = bounds[0]
                    
                    # FIXED: Use efficiency metric (mg element per g salt) and bounds check
                    if required_g_per_l <= hi:
                        efficiency = mg_per_g  # Higher efficiency = more element per gram
                        if efficiency > best_efficiency:
                            best_efficiency = efficiency
                            best_salt = salt
            
            # Set ONLY the best salt for this element
            if best_salt is not None:
                mg_per_g = STOICH_DATABASE[best_salt][element]
                required_g_per_l = target_concentration / mg_per_g
                micronutrient_seeds[best_salt] = required_g_per_l
                
                # FIXED: Zero out other salts providing the same element
                for salt in selected_salts:
                    if (salt != best_salt and 
                        salt in STOICH_DATABASE and 
                        STOICH_DATABASE[salt]['category'] == 'Micronutrient' and
                        element in STOICH_DATABASE[salt]):
                        micronutrient_seeds[salt] = 0.0  # Zero out competing salts
    
    return micronutrient_seeds

# Remove cache decorator completely to force fresh execution
def force_micronutrients_in_solution_v2_2_9_final_cache_buster(g_best, selected_salts, elem_bounds, cache_buster=None):
    """Force micronutrients to meet minimum targets if they're missing - VERSION 2.2.9 FINAL"""
    g_forced = g_best.copy()
    
    # Store forcing info for display (print doesn't show in Streamlit)
    if not hasattr(force_micronutrients_in_solution_v2_2_9_final_cache_buster, 'forcing_log'):
        force_micronutrients_in_solution_v2_2_9_final_cache_buster.forcing_log = []
    
    force_micronutrients_in_solution_v2_2_9_final_cache_buster.forcing_log.append("=== FORCING FUNCTION STARTED ===")
    
    # Check each micronutrient
    for element in ['Cu', 'Mo', 'B', 'Mn', 'Zn', 'Fe']:
        if element in elem_bounds:
            min_target, max_target = elem_bounds[element]
            current_total = 0
            
            # Calculate current contribution
            for i, salt in enumerate(selected_salts):
                if salt in STOICH_DATABASE:
                    salt_data = STOICH_DATABASE[salt]
                    if element in salt_data:
                        current_total += g_forced[i] * salt_data[element]
            
            force_micronutrients_in_solution_v2_2_9_final_cache_buster.forcing_log.append(
                f"Checking {element}: current={current_total:.8f} mg/L, target={min_target:.6f} mg/L"
            )
            
            # If below minimum, force the appropriate salt
            if current_total < min_target:
                # Find the best salt to add for this micronutrient
                best_salt_index = None
                best_salt_name = None
                
                # Log available salts for this element
                available_salts = []
                for i, salt in enumerate(selected_salts):
                    if salt in STOICH_DATABASE:
                        salt_data = STOICH_DATABASE[salt]
                        if element in salt_data:
                            available_salts.append((i, salt, salt_data[element]))
                
                force_micronutrients_in_solution_v2_2_9_final_cache_buster.forcing_log.append(
                    f"Available {element} salts: {[salt for _, salt, _ in available_salts]}"
                )
                
                for i, salt in enumerate(selected_salts):
                    if salt in STOICH_DATABASE:
                        salt_data = STOICH_DATABASE[salt]
                        if element in salt_data:
                            mg_per_g = salt_data[element]
                            
                            # FIXED: Calculate TOTAL needed, not just deficit
                            total_needed_g_per_l = min_target / mg_per_g
                            
                            # FIXED: Get bounds for this specific salt
                            all_bounds = generate_salt_bounds(selected_salts)
                            lo, hi = all_bounds[i]  # Correct index from selected_salts
                            
                            force_micronutrients_in_solution_v2_2_9_final_cache_buster.forcing_log.append(
                                f"  {salt}: {mg_per_g:.1f} mg/g, need {total_needed_g_per_l:.8f} g/L, bounds [{lo:.6f}, {hi:.6f}]"
                            )
                            
                            if total_needed_g_per_l <= hi:
                                best_salt_index = i
                                best_salt_name = salt
                                force_micronutrients_in_solution_v2_2_9_final_cache_buster.forcing_log.append(
                                    f"  ✅ Selected {salt} for {element}"
                                )
                                break
                            else:
                                force_micronutrients_in_solution_v2_2_9_final_cache_buster.forcing_log.append(
                                    f"  ❌ {salt} exceeds bounds (need {total_needed_g_per_l:.8f}, max {hi:.6f})"
                                )
                
                # Force the micronutrient salt to the total needed amount
                if best_salt_index is not None:
                    # FIXED: Set to total needed, and zero out other salts providing this element
                    total_needed_g_per_l = min_target / STOICH_DATABASE[best_salt_name][element]
                    
                    # Zero out other salts that provide this element to avoid double-counting
                    for i, salt in enumerate(selected_salts):
                        if (i != best_salt_index and salt in STOICH_DATABASE and 
                            element in STOICH_DATABASE[salt]):
                            g_forced[i] = 0.0
                    
                    # Set the best salt to provide exactly the minimum needed
                    g_forced[best_salt_index] = total_needed_g_per_l
                    
                    force_micronutrients_in_solution_v2_2_9_final_cache_buster.forcing_log.append(
                        f"FORCED {element}: Set {best_salt_name} to {total_needed_g_per_l:.8f} g/L to meet minimum {min_target:.6f} mg/L"
                    )
                else:
                    force_micronutrients_in_solution_v2_2_9_final_cache_buster.forcing_log.append(
                        f"❌ FAILED to force {element}: No suitable salt found within bounds!"
                    )
    
    force_micronutrients_in_solution_v2_2_9_final_cache_buster.forcing_log.append("=== FORCING FUNCTION COMPLETED ===")
    return g_forced

def calculate_cross_contributions(micronutrient_seeds, selected_salts):
    """Calculate what micronutrients contribute to macronutrients"""
    cross_contributions = {'N': 0, 'K': 0, 'P': 0, 'Ca': 0, 'Mg': 0, 'S': 0, 'Na': 0, 'Cl': 0}
    
    for salt, concentration in micronutrient_seeds.items():
        if salt in STOICH_DATABASE:
            salt_data = STOICH_DATABASE[salt]
            for element, mg_per_g in salt_data.items():
                if element in cross_contributions:
                    cross_contributions[element] += concentration * mg_per_g
    
    return cross_contributions

def generate_smart_seeds(selected_salts, elem_bounds, ratio_bounds, n_seeds=20):
    """Generate smart seeds using analytical chemistry knowledge"""
    seeds = []
    
    # Get micronutrient anchor points
    micronutrient_seeds = calculate_micronutrient_seeds(selected_salts, elem_bounds)
    cross_contributions = calculate_cross_contributions(micronutrient_seeds, selected_salts)
    
    # Create base seed with micronutrients locked
    base_seed = np.zeros(len(selected_salts))
    for i, salt in enumerate(selected_salts):
        if salt in micronutrient_seeds:
            base_seed[i] = micronutrient_seeds[salt]
    
    # FORCE: Ensure only one salt per element is used
    used_elements = set()
    for i, salt in enumerate(selected_salts):
        if salt in STOICH_DATABASE:
            salt_data = STOICH_DATABASE[salt]
            for element in ['Cu', 'Mo', 'Mn']:
                if element in salt_data and element not in used_elements:
                    if base_seed[i] == 0:  # If not already set by micronutrient_seeds
                        # Set to minimum bound to ensure inclusion
                        bounds = generate_salt_bounds([salt])
                        lo, hi = bounds[0]
                        base_seed[i] = lo + (hi - lo) * 0.1  # 10% of range
                        used_elements.add(element)  # Mark this element as used
                        break  # Only use one salt per element
    
    # Strategy variations for macronutrients
    strategies = [
        "high_ca",    # Use Ca(NO3)2 + CaSO4 heavily
        "high_k",     # Use KNO3 + K2SO4 heavily  
        "balanced",   # Standard ratios
        "low_salt",   # Minimal concentrations
        "high_p"      # Use phosphate salts heavily
    ]
    
    for seed_idx in range(n_seeds):
        seed = base_seed.copy()
        strategy = strategies[seed_idx % len(strategies)]
        
        # Adjust macronutrient salts based on strategy
        for i, salt in enumerate(selected_salts):
            if STOICH_DATABASE[salt]['category'] != 'Micronutrient':
                lo, hi = generate_salt_bounds([salt])[0]
                
                if strategy == "high_ca" and "Ca" in salt:
                    seed[i] = (lo + hi) * 0.8  # High Ca strategy
                elif strategy == "high_k" and "K" in salt:
                    seed[i] = (lo + hi) * 0.8  # High K strategy
                elif strategy == "balanced":
                    seed[i] = (lo + hi) * 0.5  # Balanced strategy
                elif strategy == "low_salt":
                    seed[i] = lo + (hi - lo) * 0.2  # Low salt strategy
                elif strategy == "high_p" and "P" in salt:
                    seed[i] = (lo + hi) * 0.8  # High P strategy
                else:
                    seed[i] = (lo + hi) * 0.5  # Default balanced
        
        # FIXED: Only perturb macronutrient salts, leave micronutrients at optimal values
        for i in range(len(seed)):
            salt = selected_salts[i]
            if seed[i] > 0 and STOICH_DATABASE[salt]['category'] != 'Micronutrient':
                # Add ±15% perturbations only to macronutrient salts
                seed[i] *= (0.85 + 0.3 * random.random())
            # Micronutrient salts remain at their analytically calculated optimal values
        
        # Ensure bounds
        bounds = generate_salt_bounds(selected_salts)
        for i, (lo, hi) in enumerate(bounds):
            seed[i] = max(lo, min(hi, seed[i]))
        
        seeds.append(seed)
    
    return seeds

def optimize_media(selected_salts, elem_bounds, ratio_bounds, algorithm='DE', n_trials=1):
    """Optimize media composition using specified algorithm"""
    bounds = generate_salt_bounds(selected_salts)
    
    best_result = None
    best_penalty = np.inf
    
    for trial in range(n_trials):
        # Smart seeding instead of random Monte Carlo
        random.seed(42 + trial)
        np.random.seed(42 + trial)
        seed_pool = []
        
        # Generate smart seeds (95%+ should be feasible)
        smart_seeds = generate_smart_seeds(selected_salts, elem_bounds, ratio_bounds, n_seeds=20)
        seed_pool.extend(smart_seeds)
        
        # Add a few random seeds for diversity (but much fewer)
        for _ in range(100):  # Reduced from 5000 to 100
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
    # Version check - this should appear if new code is running
    import time
    timestamp = int(time.time())
    
    st.set_page_config(
        page_title="Tissue Culture Media Optimizer",
        page_icon="🧪",
        layout="wide"
    )
    
    # Version tracking and cache clearing
    st.sidebar.markdown("---")
    st.sidebar.subheader("🔄 Cache & Version")
    
    # Add a button to clear cache
    if st.sidebar.button("🗑️ Clear Cache", help="Clear Streamlit cache if updates aren't showing"):
        st.cache_data.clear()
        st.cache_resource.clear()
        st.rerun()
    
    # Force version check
    if st.sidebar.button("🔄 Force Refresh", help="Force refresh to ensure latest code is running"):
        st.rerun()
    
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
    
    n_trials = st.sidebar.slider("Number of optimization trials", 1, 3, 1, 
                                 help="Smart seeding reduces trials needed from 5+ to 2-3")
    
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
                # Default most micronutrients to checked, but MnSO4·H2O and CuSO4·5H2O unchecked
                if salt in ['MnSO4·H2O', 'CuSO4·5H2O']:
                    default_val = False
                else:
                    default_val = True
            
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
            'Cu': (0.05, 0.3), 'Mo': (0.05, 0.2), 'Fe': (2.0, 10.0)  # INCREASED Cu max to 0.3
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
        with st.spinner("🧠 Using smart seeding with analytical chemistry knowledge..."):
            try:
                result = optimize_media(selected_salts, elem_bounds, ratio_bounds, algorithm, n_trials)
                
                if result is not None and hasattr(result, 'x'):
                    g_best = result.x.copy()  # Use unrounded for calculations
                    
                    st.write(f"**Pre-forcing penalty:** {result.fun:.2e}")
                    
                    # ALWAYS force micronutrients to meet minimum targets, regardless of penalty
                    st.write("**🔧 Calling forcing function...**")
                    st.write("**🔧 CACHE BUSTER: This should show if new code is running**")
                    st.write("**🔧 VERSION 2.2.9: Fixed duplicate penalty calculation**")
                    
                    # Show pre-forcing Cu and Mo totals
                    pre_cu_total = 0
                    pre_mo_total = 0
                    for i, salt in enumerate(selected_salts):
                        if salt in STOICH_DATABASE:
                            if 'Cu' in STOICH_DATABASE[salt]:
                                pre_cu_total += g_best[i] * STOICH_DATABASE[salt]['Cu']
                            if 'Mo' in STOICH_DATABASE[salt]:
                                pre_mo_total += g_best[i] * STOICH_DATABASE[salt]['Mo']
                    st.write(f"**Pre-forcing totals:** Cu={pre_cu_total:.8f} mg/L, Mo={pre_mo_total:.8f} mg/L")
                    
                    # Force cache refresh by adding a unique parameter
                    import time
                    cache_buster = f"{int(time.time())}_{timestamp}"
                    st.write(f"**🔧 Cache buster:** {cache_buster}")
                    g_best = force_micronutrients_in_solution_v2_2_9_final_cache_buster(g_best, selected_salts, elem_bounds, cache_buster)
                    
                    # Show post-forcing Cu and Mo totals
                    post_cu_total = 0
                    post_mo_total = 0
                    for i, salt in enumerate(selected_salts):
                        if salt in STOICH_DATABASE:
                            if 'Cu' in STOICH_DATABASE[salt]:
                                post_cu_total += g_best[i] * STOICH_DATABASE[salt]['Cu']
                            if 'Mo' in STOICH_DATABASE[salt]:
                                post_mo_total += g_best[i] * STOICH_DATABASE[salt]['Mo']
                    st.write(f"**Post-forcing totals:** Cu={post_cu_total:.8f} mg/L, Mo={post_mo_total:.8f} mg/L")
                    
                    st.write("**🔧 Forcing function completed**")
                    
                    # Recalculate penalty after forcing
                    final_penalty = penalty_function(g_best, selected_salts, elem_bounds, ratio_bounds)
                    
                    st.write(f"**Post-forcing penalty:** {final_penalty:.2e}")
                    
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
                        feasible = final_penalty < 1e5
                        st.info(f"**Optimization Status:** {'✅ Feasible' if feasible else '❌ Infeasible'} | "
                               f"**Original Penalty:** {result.fun:.2e} | **Final Penalty:** {final_penalty:.2e} | **Algorithm:** {algorithm}")
                    
                    # Debug information
                    with st.expander("🔍 Debug Information"):
                        # Version and cache info
                        st.write("**🔄 Version Info:**")
                        st.write("Version: 2.2.5 (Fixed smart seeding to target 40% of range)")
                        st.write("Cache TTL: Removed cache decorators")
                        st.write("Last Update: Fixed smart seeding to target 40% of range")
                        
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
                        
                        # Show smart seeding information
                        st.write("**Smart Seeding Analysis:**")
                        micronutrient_seeds = calculate_micronutrient_seeds(selected_salts, elem_bounds)
                        st.write("**Target micronutrient concentrations:**")
                        for element in ['B', 'Mn', 'Zn', 'Cu', 'Mo', 'Fe']:
                            if element in elem_bounds:
                                min_val, max_val = elem_bounds[element]
                                target = (min_val + max_val) / 2
                                st.write(f"{element}: {target:.3f} mg/L (target: {min_val:.3f}-{max_val:.3f})")
                        
                        st.write("**Calculated salt concentrations for micronutrients:**")
                        for salt, conc in micronutrient_seeds.items():
                            st.write(f"{salt}: {conc:.8f} g/L")
                        
                        # Check if Cu and Mo salts are in the smart seeds
                        cu_salts_in_seeds = [salt for salt in micronutrient_seeds.keys() if 'Cu' in STOICH_DATABASE[salt]]
                        mo_salts_in_seeds = [salt for salt in micronutrient_seeds.keys() if 'Mo' in STOICH_DATABASE[salt]]
                        st.write(f"**Cu salts in smart seeds:** {cu_salts_in_seeds}")
                        st.write(f"**Mo salts in smart seeds:** {mo_salts_in_seeds}")
                        
                        # Verify smart seed values are optimal
                        st.write("**Smart Seed Verification:**")
                        for salt, conc in micronutrient_seeds.items():
                            if 'Cu' in STOICH_DATABASE[salt]:
                                expected_cu = conc * STOICH_DATABASE[salt]['Cu']
                                if 'Cu' in elem_bounds:
                                    min_cu, max_cu = elem_bounds['Cu']
                                    target_cu = (min_cu + max_cu) / 2
                                    st.write(f"  {salt}: {conc:.8f} g/L → {expected_cu:.6f} mg/L Cu (target: {target_cu:.6f})")
                            elif 'Mo' in STOICH_DATABASE[salt]:
                                expected_mo = conc * STOICH_DATABASE[salt]['Mo']
                                if 'Mo' in elem_bounds:
                                    min_mo, max_mo = elem_bounds['Mo']
                                    target_mo = (min_mo + max_mo) / 2
                                    st.write(f"  {salt}: {conc:.8f} g/L → {expected_mo:.6f} mg/L Mo (target: {target_mo:.6f})")
                        
                        # Show what Cu and Mo targets should be
                        if 'Cu' in elem_bounds:
                            min_cu, max_cu = elem_bounds['Cu']
                            target_cu = (min_cu + max_cu) / 2
                            st.write(f"**Cu target in smart seeds:** {target_cu:.6f} mg/L")
                            # Calculate what the salt concentration should be
                            cu_salt = None
                            for salt in selected_salts:
                                if salt in STOICH_DATABASE and 'Cu' in STOICH_DATABASE[salt]:
                                    cu_salt = salt
                                    break
                            if cu_salt:
                                required_g_per_l = target_cu / STOICH_DATABASE[cu_salt]['Cu']
                                st.write(f"**Required {cu_salt} concentration:** {required_g_per_l:.8f} g/L")
                        if 'Mo' in elem_bounds:
                            min_mo, max_mo = elem_bounds['Mo']
                            target_mo = (min_mo + max_mo) / 2
                            st.write(f"**Mo target in smart seeds:** {target_mo:.6f} mg/L")
                            # Calculate what the salt concentration should be
                            mo_salt = None
                            for salt in selected_salts:
                                if salt in STOICH_DATABASE and 'Mo' in STOICH_DATABASE[salt]:
                                    mo_salt = salt
                                    break
                            if mo_salt:
                                required_g_per_l = target_mo / STOICH_DATABASE[mo_salt]['Mo']
                                st.write(f"**Required {mo_salt} concentration:** {required_g_per_l:.8f} g/L")
                        
                        cross_contributions = calculate_cross_contributions(micronutrient_seeds, selected_salts)
                        st.write("**Cross-contributions from micronutrients:**")
                        for element, contribution in cross_contributions.items():
                            if contribution > 0:
                                st.write(f"{element}: {contribution:.3f} mg/L")
                        
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
                        
                        # Show ALL salt concentrations (including zeros)
                        st.write("**ALL Salt Concentrations (including zeros):**")
                        for salt, conc in zip(selected_salts, g_best):
                            st.write(f"{salt}: {conc:.8f} g/L")
                        
                        # Show raw elemental totals
                        st.write("**Raw Elemental Totals (mg/L):**")
                        for element in ['Cu', 'Mo', 'B', 'Mn', 'Zn', 'Fe']:
                            if element in e_opt:
                                st.write(f"{element}: {e_opt[element]:.8f} mg/L")
                        
                        # Show raw ratios
                        st.write("**Raw Ratios:**")
                        for ratio_name in ['N:K', 'Ca:Mg', 'P:K']:
                            if ratio_name in r_opt:
                                st.write(f"{ratio_name}: {r_opt[ratio_name]:.8f}")
                        
                        # Detailed Cu and Mo analysis
                        st.write("**🔍 Detailed Cu and Mo Analysis:**")
                        
                        # Check Cu
                        cu_salts = [salt for salt in selected_salts if 'Cu' in STOICH_DATABASE[salt]]
                        st.write(f"**Copper salts available:** {cu_salts}")
                        cu_total = 0
                        for salt, conc in zip(selected_salts, g_best):
                            if salt in STOICH_DATABASE and 'Cu' in STOICH_DATABASE[salt]:
                                cu_contribution = conc * STOICH_DATABASE[salt]['Cu']
                                cu_total += cu_contribution
                                st.write(f"  {salt}: {conc:.8f} g/L × {STOICH_DATABASE[salt]['Cu']:.1f} mg/g = {cu_contribution:.8f} mg/L Cu")
                        st.write(f"**Total Cu: {cu_total:.8f} mg/L**")
                        if 'Cu' in elem_bounds:
                            min_cu, max_cu = elem_bounds['Cu']
                            st.write(f"**Cu target: {min_cu:.6f}-{max_cu:.6f} mg/L**")
                            if cu_total < min_cu:
                                st.error(f"❌ Cu below minimum! Need {min_cu - cu_total:.8f} mg/L more")
                            elif cu_total > max_cu:
                                st.error(f"❌ Cu above maximum! {cu_total - max_cu:.8f} mg/L too much")
                            else:
                                st.success(f"✅ Cu in range!")
                        
                        # Check Mo
                        mo_salts = [salt for salt in selected_salts if 'Mo' in STOICH_DATABASE[salt]]
                        st.write(f"**Molybdenum salts available:** {mo_salts}")
                        mo_total = 0
                        for salt, conc in zip(selected_salts, g_best):
                            if salt in STOICH_DATABASE and 'Mo' in STOICH_DATABASE[salt]:
                                mo_contribution = conc * STOICH_DATABASE[salt]['Mo']
                                mo_total += mo_contribution
                                st.write(f"  {salt}: {conc:.8f} g/L × {STOICH_DATABASE[salt]['Mo']:.1f} mg/g = {mo_contribution:.8f} mg/L Mo")
                        st.write(f"**Total Mo: {mo_total:.8f} mg/L**")
                        if 'Mo' in elem_bounds:
                            min_mo, max_mo = elem_bounds['Mo']
                            st.write(f"**Mo target: {min_mo:.6f}-{max_mo:.6f} mg/L**")
                            if mo_total < min_mo:
                                st.error(f"❌ Mo below minimum! Need {min_mo - mo_total:.8f} mg/L more")
                            elif mo_total > max_mo:
                                st.error(f"❌ Mo above maximum! {mo_total - max_mo:.8f} mg/L too much")
                            else:
                                st.success(f"✅ Mo in range!")
                        
                        # Check if forcing was applied
                        st.write("**🔧 Forcing Analysis:**")
                        original_g = result.x.copy()
                        forced_g = g_best.copy()
                        forcing_applied = False
                        for i, (orig, forced) in enumerate(zip(original_g, forced_g)):
                            if abs(orig - forced) > 1e-10:
                                st.write(f"  {selected_salts[i]}: {orig:.8f} → {forced:.8f} g/L (FORCED)")
                                forcing_applied = True
                        if not forcing_applied:
                            st.write("  No forcing applied - solution already met targets")
                        
                        # Show forcing calculations for Cu and Mo
                        if forcing_applied:
                            st.write("**🔧 Forcing Calculations:**")
                            for element in ['Cu', 'Mo']:
                                if element in elem_bounds:
                                    min_target, max_target = elem_bounds[element]
                                    
                                    # Calculate original total
                                    original_total = 0
                                    for i, salt in enumerate(selected_salts):
                                        if salt in STOICH_DATABASE and element in STOICH_DATABASE[salt]:
                                            original_total += original_g[i] * STOICH_DATABASE[salt][element]
                                    
                                    # Calculate forced total
                                    forced_total = 0
                                    for i, salt in enumerate(selected_salts):
                                        if salt in STOICH_DATABASE and element in STOICH_DATABASE[salt]:
                                            forced_total += forced_g[i] * STOICH_DATABASE[salt][element]
                                    
                                    st.write(f"  {element}: {original_total:.8f} → {forced_total:.8f} mg/L (target: {min_target:.6f}-{max_target:.6f})")
                                    
                                    if original_total < min_target:
                                        deficit = min_target - original_total
                                        st.write(f"    Deficit was: {deficit:.8f} mg/L")
                                        st.write(f"    Total needed: {min_target:.8f} mg/L")
                        
                        # Show forcing function logs
                        if hasattr(force_micronutrients_in_solution_v2_2_9_final_cache_buster, 'forcing_log'):
                            st.write("**🔧 Forcing Function Logs:**")
                            for log_entry in force_micronutrients_in_solution_v2_2_9_final_cache_buster.forcing_log:
                                st.write(f"  {log_entry}")
                            # Clear the log for next run
                            force_micronutrients_in_solution_v2_2_9_final_cache_buster.forcing_log = []
                        else:
                            st.write("**🔧 Forcing Function Logs:** No forcing applied")
                        
                        # Cu/Mo Feasibility Check
                        st.write("**🔍 Cu/Mo Feasibility Check:**")
                        for element in ['Cu', 'Mo']:
                            if element in elem_bounds:
                                min_target, max_target = elem_bounds[element]
                                st.write(f"**{element} target:** {min_target:.6f}-{max_target:.6f} mg/L")
                                
                                # Find all salts that provide this element
                                providing_salts = []
                                for salt in selected_salts:
                                    if salt in STOICH_DATABASE and element in STOICH_DATABASE[salt]:
                                        mg_per_g = STOICH_DATABASE[salt][element]
                                        bounds = generate_salt_bounds([salt])
                                        lo, hi = bounds[0]
                                        
                                        # Calculate achievable range
                                        min_achievable = lo * mg_per_g
                                        max_achievable = hi * mg_per_g
                                        
                                        providing_salts.append({
                                            'salt': salt,
                                            'mg_per_g': mg_per_g,
                                            'bounds': (lo, hi),
                                            'achievable_range': (min_achievable, max_achievable),
                                            'can_meet_min': max_achievable >= min_target,
                                            'can_meet_max': min_achievable <= max_target
                                        })
                                
                                st.write(f"**Salts providing {element}:**")
                                for salt_info in providing_salts:
                                    st.write(f"  {salt_info['salt']}:")
                                    st.write(f"    {salt_info['mg_per_g']:.1f} mg/g")
                                    st.write(f"    Bounds: {salt_info['bounds'][0]:.6f}-{salt_info['bounds'][1]:.6f} g/L")
                                    st.write(f"    Achievable: {salt_info['achievable_range'][0]:.6f}-{salt_info['achievable_range'][1]:.6f} mg/L")
                                    st.write(f"    Can meet minimum: {salt_info['can_meet_min']}")
                                
                                # Check if any salt can meet the minimum
                                can_meet_target = any(salt['can_meet_min'] for salt in providing_salts)
                                st.write(f"**{element} target achievable:** {can_meet_target}")
                                
                                if not can_meet_target:
                                    st.error(f"❌ NO SALT CAN MEET {element} MINIMUM TARGET!")
                                    st.write("   Consider:")
                                    st.write("   1. Lowering the minimum target")
                                    st.write("   2. Increasing salt bounds")
                                    st.write("   3. Adding different Cu/Mo salts")
                        
                        # Detailed forcing analysis
                        st.write("**🔧 Detailed Forcing Analysis:**")
                        for element in ['Cu', 'Mo']:
                            if element in elem_bounds:
                                min_target, max_target = elem_bounds[element]
                                st.write(f"**{element} Analysis:**")
                                
                                # Show all salts that provide this element
                                element_salts = []
                                for i, salt in enumerate(selected_salts):
                                    if salt in STOICH_DATABASE and element in STOICH_DATABASE[salt]:
                                        element_salts.append((i, salt, STOICH_DATABASE[salt][element]))
                                
                                st.write(f"  Salts providing {element}: {[salt for _, salt, _ in element_salts]}")
                                
                                # Show original vs forced concentrations
                                for i, salt, mg_per_g in element_salts:
                                    original_conc = original_g[i]
                                    forced_conc = forced_g[i]
                                    original_contribution = original_conc * mg_per_g
                                    forced_contribution = forced_conc * mg_per_g
                                    
                                    st.write(f"    {salt}: {original_conc:.8f} → {forced_conc:.8f} g/L")
                                    st.write(f"      Contribution: {original_contribution:.8f} → {forced_contribution:.8f} mg/L {element}")
                                
                                # Show totals
                                original_total = sum(original_g[i] * mg_per_g for i, _, mg_per_g in element_salts)
                                forced_total = sum(forced_g[i] * mg_per_g for i, _, mg_per_g in element_salts)
                                st.write(f"  Total {element}: {original_total:.8f} → {forced_total:.8f} mg/L (target: {min_target:.6f}-{max_target:.6f})")
                                
                                if original_total < min_target:
                                    st.write(f"  ✅ Forcing was needed and applied")
                                else:
                                    st.write(f"  ℹ️ No forcing needed - already met target")
                        
                        # Show DE optimal micronutrient injection values
                        st.write("**🧬 DE Optimal Micronutrient Injection:**")
                        if algorithm == 'DE':
                            # Recalculate what DE should be injecting
                            optimal_micronutrients = {}
                            for j, salt in enumerate(selected_salts):
                                if STOICH_DATABASE[salt]['category'] == 'Micronutrient':
                                    for element in ['B', 'Mn', 'Zn', 'Cu', 'Mo', 'Fe']:
                                        if element in elem_bounds and element in STOICH_DATABASE[salt]:
                                            min_val, max_val = elem_bounds[element]
                                            target = (min_val + max_val) / 2
                                            mg_per_g = STOICH_DATABASE[salt][element]
                                            optimal_g_per_l = target / mg_per_g
                                            lo, hi = generate_salt_bounds([salt])[0]
                                            optimal_micronutrients[j] = max(lo, min(hi, optimal_g_per_l))
                                            st.write(f"  {salt}: {optimal_g_per_l:.8f} g/L → {target:.6f} mg/L {element}")
                                            break
                        
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
