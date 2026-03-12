import io
from pdbfixer import PDBFixer
from openmm.app import PDBFile, ForceField, Simulation
from openmm import LangevinIntegrator
from openmm import unit


class HydrogenationService:
    def __init__(self):
        """
        INITIALIZES THE SERVICE FOR ADDING HYDROGEN ATOMS AND OPTIMIZING MOLECULAR GEOMETRY.
        """
        pass

    def add_hydrogen_atoms(self, pdb_content: str, ph: float = 7.0, optimize: bool = False) -> str:
        """
        ADDS MISSING HYDROGEN ATOMS AND OPTIMIZES THEM.
        INCLUDES HETEROGEN REMOVAL TO PREVENT FORCEFIELD ERRORS.
        """
        input_stream = io.StringIO(pdb_content)
        fixer = PDBFixer(pdbfile=input_stream)

        # 1. Basic cleanup - VERY IMPORTANT
        # This removes Glycerol (GOL) and other small molecules that
        # lack templates in the standard Amber14 force field.
        if optimize:
            fixer.removeHeterogens(False)  # False keeps water, but we usually want to remove it for simple min

        fixer.findMissingResidues()
        fixer.findMissingAtoms()
        fixer.addMissingHydrogens(ph)

        if optimize:
            try:
                # Forcefield setup
                forcefield = ForceField('amber14-all.xml', 'amber14/tip3p.xml')

                # We use the fixer's topology which now only contains recognized residues
                system = forcefield.createSystem(fixer.topology)

                integrator = LangevinIntegrator(
                    300 * unit.kelvin,
                    1 / unit.picosecond,
                    0.002 * unit.picoseconds
                )

                simulation = Simulation(fixer.topology, system, integrator)
                simulation.context.setPositions(fixer.positions)
                simulation.minimizeEnergy(maxIterations=100)

                state = simulation.context.getState(getPositions=True)
                fixer.positions = state.getPositions()
            except Exception as e:
                # Fallback: If optimization still fails, we at least have the unoptimized hydrogens
                print(f"Optimization failed, returning unoptimized structure: {e}")

        output_stream = io.StringIO()
        PDBFile.writeFile(fixer.topology, fixer.positions, output_stream, keepIds=True)

        return output_stream.getvalue()