import nglview as nv
import random
import time
import MDAnalysis as mda
import numpy as np
import py3Dmol
from MDAnalysis.coordinates.PDB import PDBWriter
import tempfile
import os


class SimulationManager:
    def __init__(self, protein_name):
        self.protein_name = protein_name
        self.steps = [
            "pdb_to_gro",
            "define_box",
            "solvate",
            "equilibration",
            "production_md"
        ]
        self.completed_steps = set()
        self.last_view = None  # Store the last view for visualization

    def clear_view(self):
        """Clears all representations from the last visualization."""
        if self.last_view:
            self.last_view.clear_representations()
            print("🔄 Cleared previous visual representations.")

    def execute_step(self, step_name):
        """Executes a step only if the previous step is completed and returns the view."""
        if step_name not in self.steps:
            print(f"❌ {step_name} is not a valid step.")
            return None

        step_index = self.steps.index(step_name)
        previous_step = self.steps[step_index - 1] if step_index > 0 else None

        if previous_step and previous_step not in self.completed_steps:
            print(f"⚠️ You must complete '{previous_step}' before running '{step_name}'.")
            return None

        # Clear previous visualization before executing
        #self.clear_view()

        # Execute the step function
        print(f"\n✅ Executing: {step_name}")
        u, view = globals()[step_name](self.protein_name)
        self.completed_steps.add(step_name)
        self.last_view = view  # Store the latest view for visualization
        return view  # Return the view for interactive use

    def show_progress(self):
        """Displays which steps have been completed and which are pending."""
        print("\n📌 Simulation Progress:")
        for step in self.steps:
            status = "✅ Done" if step in self.completed_steps else "⏳ Pending"
            print(f"- {step}: {status}")

    def get_last_view(self):
        """Returns the last viewed structure for visualization."""
        if self.last_view:
            return self.last_view
        else:
            print("⚠️ No structures have been loaded yet.")
            return None


def fake_log_message(step_name):
    """Generates a unique fake log message to simulate MD output."""
    print(f"\n=== Running {step_name} ===")
    time.sleep(1)

    def loading_animation(duration=3, message="Processing"):
        """Simulates a loading animation."""
        print(message, end="", flush=True)
        for _ in range(duration):
            time.sleep(1)
            print(".", end="", flush=True)
        print(" Done.")

    if step_name == "pdb2gmx":
        print("Converting PDB to GROMACS format...")
        loading_animation(3, "Assigning atom types")
        print("\nAtom types assigned.")
        print("\nConverted PDB successfully to GRO.")

    elif step_name == "Defining PBC Box":
        print("Setting periodic boundary conditions...")
        loading_animation(2, "Calculating box dimensions")
        print("\nBox dimensions set successfully.")
        print("\nChecking molecular placement in the box...")
        loading_animation(2, "Ensuring no overlaps detected")
        print("\nSetup of the PBC box was successful.")

    elif step_name == "Solvating the System":
        print("Adding water molecules to system...")
        loading_animation(3, "Calculating solvation shell")
        print("\nWater molecules added. System fully solvated.")
        time.sleep(2)
        print("\nAdding Ions to neutralize charge in system...")
        loading_animation(3, "Placing counterions")
        print(f"\nAdded {random.randint(10, 100)} Cl atoms and {random.randint(10, 100)} Na atoms")
        print("\nSystem neutralized.")

    elif step_name == "Energy Minimization":
        print("Starting Energy minimization process...")
        loading_animation(5, "Optimizing molecular geometry")
        time.sleep(3)
        print("\nRunning steepest descent algorithm...")
        loading_animation(5, "Minimizing energy")
        print(f"\nEnergy minimized successfully.")

    elif step_name == "Equilibration":
        print("\nInitializing NVT equilibration...")
        loading_animation(4, "Scaling velocities to 300K")
        print("\nTemperature stabilized at 300K. System ready for NPT.")
        print("\nStarting NPT equilibration...")
        loading_animation(4, "Adjusting pressure and density")
        print("\nPressure stabilized. Density converged.")

    elif step_name == "Production MD Run":
        print("Starting production MD simulation...")
        for step in range(0, 101, 10):  # Goes from 0 to 100 ns in 10 seconds
            time.sleep(1)
            print(f"Simulation progress: {step} ns / 100 ns", end="\r", flush=True)
        print("\nProduction MD run completed successfully.")
        print("\nRemoving water molecules from system...")
        loading_animation(3, "Extracting protein structure")
        print("\nSolvent molecules removed. Showing protein-only structure.")
        print("\nCorrecting periodic boundary condition")
        loading_animation(3, "Centering protein complex")
        print("\nPeriodic boundary condition sucessfully removed")
    
    else:
        print("Step completed successfully.")

def load_structure(file_path):
    """Loads a structure file into an MDAnalysis Universe and displays it in NGLView."""
    u = mda.Universe(file_path)
    #zn = u.select_atoms("resname ZN")
    view = nv.show_mdanalysis(u, default_representation=False)
    
    # Clear previous representations
    #view.clear_representations()
    
    view.add_cartoon("protein")
    #view.add_surface(zn.residues, color="yellow")  # Highlights Zn atoms
    return u, view  # Returns Universe & NGLView instance

def pdb_to_gro(protein_name):
    fake_log_message("pdb2gmx")
    u, view = load_structure(f"{protein_name}/{protein_name}_pdb2gmx.gro")
    zn = u.select_atoms("resname ZN")
    view.add_surface(zn.residues, color="yellow")
    return u, view

def define_box(protein_name):
    fake_log_message("Defining PBC Box")
    u, view = load_structure(f"{protein_name}/{protein_name}_box.pdb")
    
    # Ensure the box is displayed correctly
    box_vectors = u.dimensions[:3]
    print(f"PBC Box Dimensions: {box_vectors}")
    zn = u.select_atoms("resname ZN")
    view.add_surface(zn.residues, color="yellow")
    view.add_unitcell()
    view.center()
    return u, view

def solvate(protein_name):
    fake_log_message("Solvating the System")
    u = mda.Universe(f"{protein_name}/{protein_name}_sol.pdb")
    zn = u.select_atoms("resname ZN")
    sol = u.select_atoms("resname SOL")
    view = nv.show_mdanalysis(u, default_representation=False)
    # Clear previous representations
    view.clear_representations()
    view.add_cartoon("protein")
    view.add_point(sol.residues, color="cyan")
    view.add_surface(zn.residues, color="yellow")  # Highlights Zn atoms
    view.center()
    return u, view


def equilibration(protein_name):
    fake_log_message("Energy Minimization")
    fake_log_message("Equilibration")
    u, v = load_structure(f"{protein_name}/{protein_name}_steep1_only_prot.pdb")
    zn = u.select_atoms("resname ZN")
    v.add_surface(zn.residues, color="yellow")
    v.center()
    return u, v

def production_md(protein_name):
    """Loads an MD trajectory with optimized memory usage and stable playback."""
    fake_log_message("Production MD Run")

    
    xtc_file = f"{protein_name}/{protein_name}_rep1_dt100.xtc"
    pdb_file = f"{protein_name}/{protein_name}_steep1_only_prot.pdb"

    u = mda.Universe(pdb_file, xtc_file)

    view = py3Dmol.view(width=800, height=600)

    frames = []

    # temporary file reused each frame
    with tempfile.NamedTemporaryFile(suffix=".pdb", delete=False) as tmp:
        tmp_name = tmp.name

    for ts in u.trajectory[::10]:
        with PDBWriter(tmp_name) as W:
            W.write(u.atoms)

        with open(tmp_name) as f:
            frames.append(f.read())

    os.remove(tmp_name)

    #  REAL trajectory
    view.addModelsAsFrames("".join(frames), "pdb")

    view.setStyle({}, {})  # clear default style

    # Protein cartoon (blue)
    view.setStyle(
        {"cartoon": {"color": "blue"}}
    )
    # Zinc residues (yellow surface)
    view.addStyle(
        {"resn": "ZN"},
        {
            "sphere": {"color": "yellow", "scale": 1.5},
            "surface": {"opacity": 1.0, "color": "yellow"}
        }
    )

    view.animate({
            "loop": "forward",
            "interval": 50
        })

    view.zoomTo()

    return u, view
