# MadGraph5+MadSpin+Py8+DELPHES 
## Dark Matter Simulation Full Workflow
### Model by Benjamin Fuks: [DMSimpt-GitHub](https://github.com/BFuks/DMSimpt)

This document describes how to set up and run MG5 v3.6.3, and analyze a DM process in a Conda environment, including showering with Pythia8, MadSpin decays, and a **separate** detector simulation with Delphes (necessary because the built-in detector FAST SIM option is not available for this aMC@NLO generation).


## We are going to run it inside the same conda environment mg5_py39
### MG Version used is also the same: 3_6_3

# 1. Download & extract MG5_aMC v3.6.3
```bash=
wget https://launchpad.net/mg5amcnlo/3.0/3.6.x/+download/MG5_aMC_v3.6.3.tar.gz
tar -xvzf MG5_aMC_v3.6.3.tar.gz
cd MG5_aMC_v3_6_3
git clone https://github.com/BFuks/DMSimpt
```
# 2. Activate environment & launch MadGraph5_aMC@NLO
```bash=
conda activate mg5_py39
export MAKEFLAGS="-j$(nproc)"
./bin/MG5_aMC
```

# 3. Inside MG5_aMC prompt
### Install dependencies
We won't need `MadSTR` for this, but will do for other processes soon to come by in the tutorial.
```bash=
MG5_aMC> install lhapdf6
MG5_aMC> install pythia8
MG5_aMC> install Delphes
MG5_aMC> install MadSTR
MG5_aMC> install pythia8_hepmc3
```

# 4. Import model & define particles
```bash=
MG5_aMC> import model ./DMSimpt/DMSimpt_v2_0-F3S_cr --modelname
MG5_aMC> define yy2=yf3u2
MG5_aMC> define yy2~=yf3u2~
MG5_aMC> define dm=xs
```
# 5. Generate the process & output
```bash=
MG5_aMC> generate p p > yy2 yy2~ / a z yf3qu1 yf3qu2 yf3qu3 yf3qd1 yf3qd2  \
    yf3qd3 ys3u2 yf3u1 yf3u3 yf3d1 yf3d2 yf3d3 t yf3qu1~ yf3qu2~ yf3qu3~  \
    yf3qd1~ yf3qd2~ yf3qd3~ ys3u2~ yf3u1~ yf3u3~ yf3d1~ yf3d2~ yf3d3~ ys3qu1  \
    ys3qu2 ys3qu3 ys3qd1 ys3qd2 ys3qd3 ys3u1 ys3u3 ys3d1 ys3d2 ys3d3 t~ ys3qu1~  \
    ys3qu2~ ys3qu3~ ys3qd1~ ys3qd2~ ys3qd3~ ys3u1~ ys3u3~ ys3d1~ ys3d2~ ys3d3~  \
    xc xc~ xm xd xd~ xv xw xw~ [QCD]
MG5_aMC> output yy_qcd
```
At this point, MG5 will prompt you to install a bunch of tools. Hit Enter. Let it install those automatically for loop calculations.

# 6. Launch the run
```bash=
MG5_aMC> launch yy_qcd
```

Once the prompt to change shower and madspin comes up, type:
```bash=
> shower=PYTHIA8    (Press Enter)
> madspin=ON        (Press Enter)
#Press Enter once again
# on being prompted to edit MadSpin card type in
> 3
```
and then 

```vim=
vim> Delete existing content using ESC :1,%d Enter
vim> Paste the MadSpin Card below without redundancy; save and exit (ESC :wq)
``` 

# 7. Edit the MadSpin card 

##    The `import model` directive must be at the top to preserve definitions
Replace my username with yours; MS can't resolve $USER
MadSpin can't identify the model, so we have to force feed the model into it.
```bash=
set max_weight_ps_point 400
import model /home/your_username/mg5-tutorial/madgraph_tutorial/MG5_aMC_v3_6_3/DMSimpt/DMSimpt_v2_0-F3S_cr --bypass_check
# specify the decay for the final state particles
decay t    > w+  b, w+  > all all
decay t~   > w-  b~, w-  > all all
decay w+   > all all
decay w-   > all all
decay z    > all all
decay yf3u2  > xs c
decay yf3u2~ > xs c~
# running the actual code
launch
```
Let that thing cook, bake, whatever.

# 8. Separate Delphes detector simulation 
### Detector option is not available in aMC@NLO for this process
This needs troubleshooting; but we can do it outside of MG5.

ab refers to run number like 01, 07 etc. Also need to unzip the HEP MC GZ file before running DELPHES fast sim
```bash=
cd /home/$USER/mg5-tutorial/madgraph_tutorial/MG5_aMC_v3_6_3/Delphes
gunzip /home/$USER/mg5-tutorial/madgraph_tutorial/MG5_aMC_v3_6_3/yy_qcd/Events/run_ab_decayed_1/events_PYTHIA8_0.hepmc.gz
# Have to use HEPMC2 because HEPMC3 gives "Invalid Vertex Format"
./DelphesHepMC2 cards/delphes_card_ATLAS.tcl delphes.root /home/$USER/mg5-tutorial/madgraph_tutorial/MG5_aMC_v3_6_3/yy_qcd/Events/run_xy_decayed_1/events_PYTHIA8_0.hepmc
```

# 9. Analyze "delphes.root" using Python scripts
```bash=
cd /home/$USER/mg5-tutorial/madgraph_tutorial/MG5_aMC_v3_6_3/Delphes/
```
## Prepare a script `analysis.py` or open a Jupyter  to run the cell(s) below
```python=
import os
import uproot
import numpy as np
import matplotlib.pyplot as plt

# 1. Locate delphes.root automatically
delphes_file = None
for root, dirs, files in os.walk(os.getcwd()):
    if 'delphes.root' in files:
        delphes_file = os.path.join(root, 'delphes.root')
        break

if delphes_file is None:
    raise FileNotFoundError("delphes.root not found in current directory tree.")
print(f"Using Delphes file: {delphes_file}")

# 2. Open ROOT file and extract JET and MET branches
with uproot.open(delphes_file) as f:
    tree = f['Delphes']
    jet_pt = tree['Jet.PT'].array(library="np")
    met_pt = tree['MissingET.MET'].array(library="np")

# 3. Flatten arrays and compute basic statistics
jet_pt_flat = np.concatenate(jet_pt)
met_pt_flat = np.concatenate(met_pt)

print(f"Jet pT: mean = {np.mean(jet_pt_flat):.2f} GeV, std = {np.std(jet_pt_flat):.2f} GeV")
print(f"MET   : mean = {np.mean(met_pt_flat):.2f} GeV, std = {np.std(met_pt_flat):.2f} GeV")

# 4. Plot distributions
plt.figure()
plt.hist(jet_pt_flat, bins=50)
plt.title("Jet pT Distribution")
plt.xlabel("pT [GeV]")
plt.ylabel("Entries")
plt.show()

plt.figure()
plt.hist(met_pt_flat, bins=50)
plt.title("Missing ET Distribution")
plt.xlabel("MET [GeV]")
plt.ylabel("Entries")
plt.show()
```

# To Plot the mass of the xs DM scalar run
```python=
import os
import uproot
import awkward as ak
import numpy as np
import matplotlib.pyplot as plt

# 1. Locate delphes.root
delphes_file = None
for root, dirs, files in os.walk(os.getcwd()):
    if "delphes.root" in files:
        delphes_file = os.path.join(root, "delphes.root")
        break
if delphes_file is None:
    raise FileNotFoundError("delphes.root not found")
print(f"Using: {delphes_file}")

# 2. Open tree and read PID + mass (or compute it)
with uproot.open(delphes_file) as f:
    tree   = f["Delphes"]
    pid_ak = tree["Particle/Particle.PID"].array(library="ak")
    if "Particle/Particle.Mass" in tree.keys():
        mass_ak = tree["Particle/Particle.Mass"].array(library="ak")
    else:
        E_ak  = tree["Particle/Particle.E"].array(library="ak")
        px_ak = tree["Particle/Particle.Px"].array(library="ak")
        py_ak = tree["Particle/Particle.Py"].array(library="ak")
        pz_ak = tree["Particle/Particle.Pz"].array(library="ak")
        mass_ak = ak.sqrt(ak.clip(E_ak**2 - (px_ak**2 + py_ak**2 + pz_ak**2), 0, None))

# 3. Select PDGID = ±51
mask     = (pid_ak == 51) | (pid_ak == -51)
selected = mass_ak[mask]

# 4. Flatten and convert to a pure Python list, then to NumPy
flat_mass = ak.flatten(selected)
mass_list = ak.to_list(flat_mass)          # list of floats
mass_np   = np.array(mass_list, dtype=float)

# 5. Ensure that particle 51 (or –51) exists
if mass_np.size == 0:
    print("No entries found for PDGID = ±51.")
else:
    print(f"Found {mass_np.size} entries for PDGID = ±51")

    # 6. Plot mass distribution
    plt.figure()
    plt.hist(mass_np, bins=50)
    plt.title("Mass Distribution for PDGID = ±51")
    plt.xlabel("Mass [GeV]")
    plt.ylabel("Entries")
    plt.tight_layout()
    plt.show()
```

# Fancy Gemini Style
```python=
import os
import uproot
import awkward as ak
import numpy as np
import matplotlib.pyplot as plt

def analyze_delphes_data(pdgid_target=51, bins=50):
    """
    Analyzes Delphes.root file to plot the mass distribution of a specific particle PDGID.

    Args:
        pdgid_target (int): The absolute PDGID of the particle to analyze (e.g., 51 for ±51).
        bins (int): The number of bins for the histogram.
    """
    # --- 1. Locate delphes.root ---
    delphes_file = None
    try:
        print("Searching for 'delphes.root' in current and subdirectories...")
        for root, dirs, files in os.walk(os.getcwd()):
            if "delphes.root" in files:
                delphes_file = os.path.join(root, "delphes.root")
                break
        if delphes_file is None:
            raise FileNotFoundError("The 'delphes.root' file was not found in the current directory or any subdirectories.")
        print(f"Successfully found Delphes file: {delphes_file}")

    except FileNotFoundError as e:
        print(f"Error: {e}")
        print("Please ensure 'delphes.root' is accessible from the script's execution path.")
        return # Exit the function if file not found
    except Exception as e:
        print(f"An unexpected error occurred during file search: {e}")
        return

    # --- 2. Open tree and read PID + mass (or compute it) ---
    pid_ak = ak.Array([]) # Initialize to empty awkward arrays
    mass_ak = ak.Array([]) # to ensure they are defined in case of early exit

    try:
        with uproot.open(delphes_file) as f:
            if "Delphes" not in f:
                raise KeyError("The 'Delphes' TTree was not found in the specified 'delphes.root' file.")
            tree = f["Delphes"]

            # Define expected branch names based on common Delphes output
            # If these are still incorrect for your file, you'll need to adjust them.
            # (Based on your previous output, these appear to be correct for the ROOT file structure)
            PID_BRANCH_NAME = "Particle/Particle.PID"
            MASS_BRANCH_NAME = "Particle/Particle.Mass"
            E_BRANCH_NAME = "Particle/Particle.E"
            PX_BRANCH_NAME = "Particle/Particle.Px"
            PY_BRANCH_NAME = "Particle/Particle.Py"
            PZ_BRANCH_NAME = "Particle/Particle.Pz"

            # Check for PID branch existence
            if PID_BRANCH_NAME not in tree.keys():
                raise KeyError(f"Particle ID branch '{PID_BRANCH_NAME}' not found in 'Delphes' tree. "
                               "Please verify the exact branch name in your Delphes.root file.")

            # Read PID directly
            pid_ak = tree[PID_BRANCH_NAME].array(library="ak")

            # Check for Mass branch existence or prepare for calculation
            if MASS_BRANCH_NAME in tree.keys():
                mass_ak = tree[MASS_BRANCH_NAME].array(library="ak")
            else:
                # Check if all 4-momentum components are present for mass calculation
                momentum_branches = [E_BRANCH_NAME, PX_BRANCH_NAME, PY_BRANCH_NAME, PZ_BRANCH_NAME]
                if not all(b in tree.keys() for b in momentum_branches):
                    missing_momentum_branches = [b for b in momentum_branches if b not in tree.keys()]
                    raise KeyError(f"Required 4-momentum branches for mass calculation are missing: {missing_momentum_branches}. "
                                   "Cannot calculate particle mass without them.")

                # Read 4-momentum components
                E_ak = tree[E_BRANCH_NAME].array(library="ak")
                px_ak = tree[PX_BRANCH_NAME].array(library="ak")
                py_ak = tree[PY_BRANCH_NAME].array(library="ak")
                pz_ak = tree[PZ_BRANCH_NAME].array(library="ak")

                # Compute invariant mass, ensuring non-negative argument to sqrt
                mass_squared = E_ak**2 - (px_ak**2 + py_ak**2 + pz_ak**2)
                mass_ak = ak.sqrt(ak.clip(mass_squared, 0, None))
                print("Particle mass calculated from 4-momentum (E, Px, Py, Pz).")

    except uproot.exceptions.KeyError as e: # Catch specific uproot KeyErrors
        print(f"Data Access Error: {e}")
        print("This likely means a TTree or branch name is incorrect. Please double-check.")
        return
    except Exception as e:
        print(f"An unexpected error occurred during data reading or initial processing: {e}")
        return

    # --- 3. Select PDGID ---
    # Using abs() for cleanliness when comparing with pdgid_target
    mask = (abs(pid_ak) == pdgid_target)
    selected_masses = mass_ak[mask]

    # --- 4. Flatten and convert to NumPy ---
    # More robust conversion for plotting
    try:
        flat_mass = ak.flatten(selected_masses)
        mass_np = np.asarray(flat_mass, dtype=float)
        if mass_np.ndim != 1: # Ensure it's a 1D array
            raise ValueError(f"Flattened mass array is not 1-dimensional. Shape: {mass_np.shape}")
    except Exception as e:
        print(f"Error during data flattening or conversion to NumPy: {e}")
        print("This might indicate an issue with the structure of the selected awkward array.")
        return

    # --- 5. Ensure that selected particles exist ---
    if mass_np.size == 0:
        print(f"No entries found for PDGID = ±{pdgid_target}. No plot will be generated.")
        return # Exit if no data to plot
    else:
        print(f"Found {mass_np.size} entries for PDGID = ±{pdgid_target}")

        # --- 6. Plot mass distribution ---
        try:
            plt.figure(figsize=(10, 6))

            # Handle cases for histogram bins
            if mass_np.min() != mass_np.max():
                plt.hist(mass_np, bins=bins, edgecolor='black', alpha=0.7)
            else:
                # If all masses are identical, a single bin is appropriate
                plt.xlim(9, 11) # Extend x-axis slightly
                plt.hist(mass_np, bins=20, edgecolor='black', alpha=0.7)
                print(f"Warning: All selected masses are identical ({mass_np[0]:.2f} GeV). Plotting with 1 bin.")

            plt.title(f"Mass Distribution for PDGID = ±{pdgid_target} (N={mass_np.size})", fontsize=14)
            plt.xlabel("Mass [GeV]", fontsize=12)
            plt.ylabel("Entries", fontsize=12)
            #plt.grid(axis='y', alpha=0.75)

            # Apply log scale if data is not all zero and has some spread
            if mass_np.size > 0 and np.all(mass_np >= 0): # Check for non-negative values for log scale
                if np.unique(mass_np).size > 1: # Only use log scale if there's more than one unique value
                    plt.yscale('log')

            plt.tight_layout()
            plt.show()

        except Exception as e:
            print(f"An error occurred during plotting: {e}")
            print("Please check your matplotlib installation and data validity for plotting.")
            return

    print("\nDelphes analysis script finished.")

# --- How to run the analysis ---
if __name__ == "__main__":
    # You can change the PDGID here, e.g., to 6 for top quarks
    analyze_delphes_data(pdgid_target=51, bins=50)
    # analyze_delphes_data(pdgid_target=6, bins=100) # Example for top quarks
```
