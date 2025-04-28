import streamlit as st
import numpy as np
import pandas as pd
import math

def calculate_fault_current(base_kva, voltage_lt, impedance_percent):
    """Calculate fault current based on transformer parameters
    
    Modified to use standard formula:
    If = (kVA * 1000) / (√3 * V * Z%)
    where If is in Amps, kVA is transformer rating, V is in volts, Z% is in percent
    """
    # Return fault current in Amps
    return (base_kva * 1000) / (math.sqrt(3) * voltage_lt * 1000 * (impedance_percent / 100))

def calculate_conductor_area_simple(fault_current, operation_time, k_factor):
    """Calculate required cross-sectional area of conductor using simple formula"""
    return (fault_current * math.sqrt(operation_time)) / k_factor

def calculate_conductor_area_ieee(fault_current_ka, tc, tcap, rho, k0, tm, ta):
    """Calculate required cross-sectional area of conductor using IEEE-80 formula
    
    Formula: A = I * √(tc * TCAP * ρr * 10⁻⁴) / (log₁₀(K₀ + Tm) - log₁₀(K₀ + Ta))
    where:
    - A = conductor cross-sectional area in mm²
    - I = fault current in kA (will be multiplied by 1000 for calculation)
    - tc = fault duration in seconds
    - TCAP = thermal capacity of material in J/(cm³·°C)
    - ρr = resistivity of material at 20°C in μΩ-cm 
    - K₀ = material constant
    - Tm = maximum allowable temperature (°C)
    - Ta = ambient temperature (°C)
    """
    # Convert fault current from kA to A for calculation
    fault_current = fault_current_ka * 1000
    
    # Factor from IEEE-80 - includes unit conversion factors for SI units to get mm²
    # Note the 10⁻⁴ factor in the formula is for unit conversions within IEEE-80
    numerator = fault_current * math.sqrt(tc * tcap * rho * 1e-4)
    denominator = math.log10(k0 + tm) - math.log10(k0 + ta)
    if denominator == 0:
        raise ValueError("Denominator in IEEE-80 formula is zero. Ensure max temp > ambient temp.")
    
    return numerator / denominator

def calculate_strip_resistance(soil_resistivity, length_cm, width_cm, thickness_cm, num_runs):
    """Calculate resistance of earthing strip"""
    total_length_cm = length_cm * num_runs
    return (100 * soil_resistivity * math.log(4 * total_length_cm / thickness_cm)) / (2 * math.pi * total_length_cm)

def calculate_rod_resistance(soil_resistivity, length_cm, diameter_cm, num_electrodes):
    """Calculate resistance of rod electrodes"""
    single_rod = (100 * soil_resistivity * math.log(2 * length_cm / diameter_cm)) / (2 * math.pi * length_cm)
    return single_rod / num_electrodes

def calculate_net_resistance(strip_resistance, rod_resistance):
    """Calculate combined resistance of the earthing system"""
    return 1 / ((1 / strip_resistance) + (1 / rod_resistance))

# Set up the Streamlit app
st.set_page_config(page_title="Solar PV AC Earthing Calculator (IEEE-80)", layout="wide")

st.title("Solar PV AC Earthing System Calculator (IEEE-80)")
st.markdown("""
This application calculates the required earthing system parameters for a solar PV installation's AC side 
based on IEEE-80 standards. Enter your system parameters below to get the complete earthing design.
""")

# Create sidebar for inputs
st.sidebar.header("System Parameters")

# Transformer parameters
st.sidebar.subheader("Transformer Parameters")
base_kva = st.sidebar.number_input("Transformer Capacity (kVA)", min_value=100, value=2000, step=100)
voltage_lt = st.sidebar.number_input("Secondary Voltage (kV)", min_value=0.1, value=0.433, step=0.001, format="%.3f")
impedance_percent = st.sidebar.number_input("Transformer Impedance (%)", min_value=1.0, value=6.0, step=0.1)

# Material selection
st.sidebar.subheader("Material Selection")
material_categories = ["Copper", "Aluminum", "Steel"]
material_category = st.sidebar.selectbox("Material Category", material_categories)

# IEEE-80 material constants
ieee_materials = {
    "Copper": {
        "Copper, annealed soft-drawn": {
            "conductivity": 100.0,
            "resistivity": 1.72,
            "tcap": 3.4,
            "fusing_temp": 1083,
            "alpha_r": 0.00393,
            "k0": 234
        },
        "Copper, commercial hard-drawn": {
            "conductivity": 97.0,
            "resistivity": 1.78,
            "tcap": 3.4,
            "fusing_temp": 1084,
            "alpha_r": 0.00381,
            "k0": 242
        }
    },
    "Aluminum": {
        "Aluminum": {
            "conductivity": 61.0,
            "resistivity": 2.8,
            "tcap": 2.7,
            "fusing_temp": 660,
            "alpha_r": 0.00403,
            "k0": 228
        },
        "Aluminum-clad steel wire": {
            "conductivity": 20.3,
            "resistivity": 8.48,
            "tcap": 3.561,
            "fusing_temp": 657,
            "alpha_r": 0.00360,
            "k0": 258
        }
    },
    "Steel": {
        "Steel, 1020": {
            "conductivity": 10.8,
            "resistivity": 15.90,
            "tcap": 3.8,
            "fusing_temp": 1510,
            "alpha_r": 0.00377,
            "k0": 245
        },
        "Zinc-coated steel rod": {
            "conductivity": 8.6,
            "resistivity": 20.10,
            "tcap": 3.9,
            "fusing_temp": 419,
            "alpha_r": 0.00320,
            "k0": 293
        },
        "Stainless steel, 304": {
            "conductivity": 2.4,
            "resistivity": 72.00,
            "tcap": 4.0,
            "fusing_temp": 1400,
            "alpha_r": 0.00130,
            "k0": 749
        }
    }
}

material_options = list(ieee_materials[material_category].keys())
material = st.sidebar.selectbox(f"Select {material_category} Type", material_options)

# Get material constants
material_constants = ieee_materials[material_category][material]

# Temperature parameters
st.sidebar.subheader("Temperature Parameters")
ambient_temp = st.sidebar.number_input("Ambient Temperature (°C)", min_value=0, value=40, step=1)
# Ensure max_temp is always at least 1°C greater than ambient_temp
min_max_temp = ambient_temp + 1
max_temp = st.sidebar.number_input(
    "Maximum Allowable Temperature (°C)", 
    min_value=min_max_temp, 
    value=min(material_constants["fusing_temp"] - 100, 300) if min(material_constants["fusing_temp"] - 100, 300) > min_max_temp else min_max_temp,
    step=10,
    help=f"Max fusing temperature for {material}: {material_constants['fusing_temp']}°C"
)

# Fault calculation parameters
st.sidebar.subheader("Fault Calculation Parameters")
fault_duration_options = ["0.5 seconds", "1 second", "3 seconds"]
fault_duration = st.sidebar.radio("Fault Duration", fault_duration_options)

# Set operation time from fault duration
operation_time_map = {"0.5 seconds": 0.5, "1 second": 1.0, "3 seconds": 3.0}
operation_time = operation_time_map[fault_duration]

corrosion_allowance = st.sidebar.number_input("Corrosion Allowance (%)", min_value=0, value=30, step=5)

# Set resistance threshold
resistance_threshold = st.sidebar.number_input("Maximum Allowable Resistance (Ω)", min_value=0.1, value=1.0, step=0.1, 
                                               help="Solar PV systems typically require < 1 Ω")

# Soil and installation parameters
st.sidebar.subheader("Soil & Installation Parameters")
soil_resistivity = st.sidebar.number_input("Soil Resistivity (Ω-m)", min_value=1, value=300, step=1)
strip_length = st.sidebar.number_input("Length of Earth Strip (m)", min_value=1, value=205, step=1)
num_runs = st.sidebar.number_input("Number of Strip Runs", min_value=1, value=2, step=1)

# Electrode parameters
st.sidebar.subheader("Electrode Parameters")
num_electrodes = st.sidebar.number_input("Number of Electrodes", min_value=1, value=4, step=1)
electrode_length = st.sidebar.number_input("Electrode Length (m)", min_value=1.0, value=3.0, step=0.5)
electrode_diameter = st.sidebar.number_input("Electrode Diameter (mm)", min_value=10.0, value=17.2, step=0.1)

# Main calculations
fault_current = calculate_fault_current(base_kva, voltage_lt, impedance_percent)
fault_current_ka = fault_current / 1000

# Required conductor area using IEEE-80 formula
required_area_ieee = calculate_conductor_area_ieee(
    fault_current_ka,
    operation_time,
    material_constants["tcap"],
    material_constants["resistivity"],
    material_constants["k0"],
    max_temp,
    ambient_temp
)

# Apply corrosion allowance
total_required_area = required_area_ieee * (1 + corrosion_allowance / 100)

# Strip selection based on material
strip_dimension_options = {
    "Copper": {
        "width": [20, 25, 30, 40, 50],
        "thickness": [3, 5, 6, 8]
    },
    "Aluminum": {
        "width": [25, 32, 40, 50, 63],
        "thickness": [3, 6, 8, 10]
    },
    "Steel": {
        "width": [25, 30, 40, 50, 75, 100],
        "thickness": [3, 5, 6, 8, 10]
    }
}

# Function to automatically select optimal strip dimensions
def select_optimal_strip_dimensions(required_area, material_category, num_runs):
    """
    Automatically select the most economical strip dimensions that meet or exceed the required area.
    Returns: (width_mm, thickness_mm, actual_area_mm2)
    """
    width_options = strip_dimension_options[material_category]["width"]
    thickness_options = strip_dimension_options[material_category]["thickness"]
    
    # Calculate the area needed per single strip run
    area_per_run = required_area / num_runs
    
    # Sort options by area to ensure we find the smallest viable combination
    combinations = []
    for width in width_options:
        for thickness in thickness_options:
            area = width * thickness
            combinations.append((width, thickness, area))
    
    # Sort by area (smallest first)
    combinations.sort(key=lambda x: x[2])
    
    # Find the first combination that meets or exceeds the required area per run
    for width, thickness, area in combinations:
        if area >= area_per_run:
            total_area = area * num_runs
            return width, thickness, total_area
    
    # If no combination is sufficient, return the largest available
    width, thickness, area = combinations[-1]
    total_area = area * num_runs
    return width, thickness, total_area

# Automatically select optimal strip dimensions based on required area
selected_width, selected_thickness, actual_area = select_optimal_strip_dimensions(
    total_required_area,
    material_category,
    num_runs
)

# Calculate resistances
strip_resistance = calculate_strip_resistance(
    soil_resistivity, 
    strip_length * 100,  # Convert to cm
    selected_width / 10,    # Convert to cm
    selected_thickness / 10,  # Convert to cm
    num_runs
)

rod_resistance = calculate_rod_resistance(
    soil_resistivity,
    electrode_length * 100,  # Convert to cm
    electrode_diameter / 10,  # Convert to cm
    num_electrodes
)

net_resistance = calculate_net_resistance(strip_resistance, rod_resistance)

# Display results
col1, col2 = st.columns(2)

with col1:
    st.header("Fault Current Calculation")
    fault_data = {
        "Parameter": ["Base kVA", "Secondary Voltage (kV)", "Impedance (%)", "Fault Current (A)", "Fault Current (kA)"],
        "Value": [base_kva, voltage_lt, impedance_percent, f"{fault_current:.2f}", f"{fault_current_ka:.2f}"]
    }
    st.table(pd.DataFrame(fault_data))
    
    st.header("IEEE-80 Material Properties")
    material_data = {
        "Property": [
            "Material Type",
            "Conductivity (% IACS)",
            "Resistivity at 20°C (μΩ-cm)",
            "Thermal Capacity (J/cm³·°C)",
            "Fusing Temperature (°C)",
            "Temperature Coefficient αr (1/°C)",
            "K₀ at 0°C"
        ],
        "Value": [
            material,
            f"{material_constants['conductivity']}",
            f"{material_constants['resistivity']}",
            f"{material_constants['tcap']}",
            f"{material_constants['fusing_temp']}",
            f"{material_constants['alpha_r']}",
            f"{material_constants['k0']}"
        ]
    }
    st.table(pd.DataFrame(material_data))
    
with col2:
    st.header("Conductor Sizing")
    st.write(f"**Fault Duration:** {fault_duration}")
    
    sizing_data = {
        "Parameter": [
            "Required Area per IEEE-80 (mm²)",
            "With Corrosion Allowance (mm²)",
            "Automatically Selected Strip Width (mm)",
            "Automatically Selected Strip Thickness (mm)",
            "Actual Total Area (mm²)",
            "Is Selected Size Adequate?",
        ],
        "Value": [
            f"{required_area_ieee:.2f}",
            f"{total_required_area:.2f}",
            f"{selected_width}",
            f"{selected_thickness}",
            f"{actual_area:.2f}",
            "✓ Yes" if actual_area >= total_required_area else "✗ No"
        ]
    }
    st.table(pd.DataFrame(sizing_data))
    
    st.header("Earthing Resistance Calculation")
    resistance_data = {
        "Parameter": [
            "Soil Resistivity (Ω-m)",
            "Strip Length (m)",
            "Total Strip Length (m)",
            "Strip Resistance (Ω)",
            "Number of Electrodes",
            "Electrode Length (m)",
            "Electrode Diameter (mm)",
            "Rod Electrode Resistance (Ω)",
            "Net System Resistance (Ω)",
            "Maximum Allowable Resistance (Ω)",
            "Status"
        ],
        "Value": [
            f"{soil_resistivity}",
            f"{strip_length}",
            f"{strip_length * num_runs}",
            f"{strip_resistance:.3f}",
            f"{num_electrodes}",
            f"{electrode_length}",
            f"{electrode_diameter}",
            f"{rod_resistance:.3f}",
            f"{net_resistance:.3f}",
            f"{resistance_threshold}",
            "✅ Acceptable" if net_resistance < resistance_threshold else "❌ Too High"
        ]
    }
    st.table(pd.DataFrame(resistance_data))

# Reference: Show IEEE-80 material table
st.header("Reference: IEEE-80-2013 Material Constants")
ieee_material_data = []

for category, materials in ieee_materials.items():
    for mat_name, constants in materials.items():
        ieee_material_data.append({
            "Material": mat_name,
            "Conductivity (% IACS)": constants["conductivity"],
            "Resistivity (μΩ-cm)": constants["resistivity"],
            "Thermal Capacity (J/cm³·°C)": constants["tcap"],
            "Fusing Temp. (°C)": constants["fusing_temp"],
            "αr (1/°C)": constants["alpha_r"],
            "K₀ at 0°C": constants["k0"]
        })

st.table(pd.DataFrame(ieee_material_data))

# Visualization
st.header("System Adequacy Visualization")
col1, col2, col3 = st.columns(3)

with col1:
    area_adequacy = actual_area / total_required_area * 100
    area_diff = actual_area - total_required_area
    st.metric("Area Adequacy", f"{area_adequacy:.1f}%", 
              f"{area_diff:.1f} mm²")
    
    # Progress bar for area adequacy - ensure non-negative value
    st.progress(max(min(area_adequacy/100, 1.0), 0.0))

with col2:
    # Fix the resistance adequacy calculation to avoid negative values
    if net_resistance < resistance_threshold:
        resistance_adequacy = (resistance_threshold - net_resistance) / resistance_threshold * 100
        resistance_diff = resistance_threshold - net_resistance
        resistance_status = f"{resistance_diff:.2f} Ω below max"
    else:
        resistance_adequacy = 0
        resistance_diff = net_resistance - resistance_threshold
        resistance_status = f"{resistance_diff:.2f} Ω above max"
    
    st.metric("Resistance Adequacy", f"{resistance_adequacy:.1f}%", resistance_status)
    
    # Progress bar - ensure non-negative value
    st.progress(max(min(resistance_adequacy/100, 1.0), 0.0))

with col3:
    # Safety factor calculation - ensure non-negative
    safety_factor = min(max(area_adequacy/100, 0.0), max(resistance_adequacy/100, 0.0)) * 100
    st.metric("Overall Safety Factor", f"{safety_factor:.1f}%")
    
    # Progress bar - ensure non-negative value
    st.progress(max(min(safety_factor/100, 1.0), 0.0))

# Market Availability Warning
if actual_area < total_required_area:
    st.warning(f"""
    ⚠️ **Market Availability Consideration**
    
    The current selection provides {actual_area} mm², which is {total_required_area - actual_area:.2f} mm² 
    less than the calculated requirement of {total_required_area:.2f} mm².
    
    This may be due to market availability constraints. If this is an intentional compromise, please:
    1. Ensure this deviation is properly documented
    2. Consider additional safety measures
    3. Obtain proper engineering approval for this deviation
    """)

# Additional design recommendations
st.header("Earthing System Design Recommendations")

if actual_area < total_required_area:
    st.warning(f"""
    ⚠️ The current strip size ({selected_width} × {selected_thickness} mm) with {num_runs} runs 
    is insufficient for the calculated fault current.
    
    Consider one of the following options:
    - Increase strip width or thickness
    - Add more parallel runs
    - Use a material with higher current carrying capacity
    """)
else:
    st.success(f"""
    ✅ The selected strip size ({selected_width} × {selected_thickness} mm) with {num_runs} runs 
    is adequate for the calculated fault current.
    """)

if net_resistance >= resistance_threshold:
    st.warning(f"""
    ⚠️ The calculated system resistance ({net_resistance:.2f} Ω) exceeds the maximum threshold of {resistance_threshold} Ω.
    
    Consider one of the following options:
    - Add more electrodes
    - Increase electrode length
    - Add more earthing strips
    - Implement soil treatment to reduce resistivity
    """)
else:
    st.success(f"""
    ✅ The calculated system resistance ({net_resistance:.2f} Ω) is below the maximum threshold of {resistance_threshold} Ω.
    """)

# IEEE Formula explanation
st.header("IEEE-80 Formula Explanation")
st.markdown(f"""
The IEEE-80 standard uses a comprehensive formula that accounts for thermal properties and temperature rise:

```
A = I × √(tc × TCAP × ρr × 10⁻⁴) / (log₁₀(K₀ + Tm) - log₁₀(K₀ + Ta))
```

Where:
- A = conductor cross-sectional area in mm²
- I = fault current in kA (multiplied by 1000 in the calculation)
- tc = fault duration in seconds ({operation_time} s)
- TCAP = thermal capacity of material in J/(cm³·°C) ({material_constants['tcap']} for {material})
- ρr = resistivity of material at 20°C in μΩ-cm ({material_constants['resistivity']} for {material})
- K₀ = material constant ({material_constants['k0']} for {material})
- Tm = maximum allowable temperature ({max_temp}°C)
- Ta = ambient temperature ({ambient_temp}°C)

This produces a more precise requirement than the simpler k-factor method.
""")

# Display design summary
st.header("Earthing System Design Summary")
st.markdown(f"""
### AC Earthing System Design (IEEE-80)

1. **Fault Current Calculation**:
   - System Fault Current: **{fault_current_ka:.2f} kA**
   - Based on {base_kva} kVA transformer with {impedance_percent}% impedance
   - Fault Duration: {fault_duration}

2. **Earthing Conductor**:
   - Material: **{material_category} ({material})**
   - Required: **{total_required_area:.2f} mm²** (including {corrosion_allowance}% corrosion allowance)
   - Selected: **{selected_width} × {selected_thickness} mm Strip** with **{num_runs}** runs
   - Total Area: **{actual_area} mm²**
   - Adequacy: **{"Adequate" if actual_area >= total_required_area else "Inadequate"}**

3. **Earthing Electrodes**:
   - **{num_electrodes}** electrodes
   - Each **{electrode_length} m** long with **{electrode_diameter} mm** diameter

4. **System Resistance**:
   - Net Resistance: **{net_resistance:.3f} Ω**
   - Standard Requirement for Solar PV: < {resistance_threshold} Ω
   - Status: **{"Acceptable" if net_resistance < resistance_threshold else "Unacceptable"}**

This design {"complies with" if (actual_area >= total_required_area and net_resistance < resistance_threshold) else "does not fully comply with"} IEEE-80 standards for earthing systems.
""")

# Note about earthing design
st.info("""
**Note**: This calculation is based on IEEE-80 standards for earthing system design. 
The actual implementation should consider local regulations, specific site conditions, 
and be approved by a qualified electrical engineer.

Solar PV systems typically require a resistance threshold of <1 ohm for proper operation and safety.
Any deviations from calculated requirements should be properly documented and approved.
""")

# IEEE Standard Footer
st.markdown("""
---
### IEEE-80-2013 Standard Information
IEEE Std 80-2013, "IEEE Guide for Safety in AC Substation Grounding," provides comprehensive guidance for 
designing safe grounding systems. The methods have been adapted here for solar PV applications.
""")