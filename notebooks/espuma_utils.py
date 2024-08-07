from espuma import Case_Directory
from pyvista import Line
from numpy import argmin
from re import findall

def get_profile_over_depth(of: Case_Directory, t_target_days:float = -1):
    """
    Return a profile over depth at the latest time.
    """
    reader = of.get_vtk_reader()
    if t_target_days == -1:
        time_of_interest_index = -1
    else:
        time_of_interest_index = argmin([abs(t_target_days - t/86400) for t in reader.time_values])
    
    time_of_interest = reader.time_values[time_of_interest_index]
    reader.set_active_time_value(time_of_interest)
    mesh = reader.read()  #<- Read the data
    internalMesh = mesh["internalMesh"]

    (xi, yi, zi, xf ,yf, zf) = internalMesh.bounds

    _ = Line(
        a := [xi, yi, zi],
        b := [xi, yi, zf]
    )

    line_probe = internalMesh.sample_over_line(a, b)
    
    return line_probe


def get_stacked_biomass(of: Case_Directory):
    reader = of.get_vtk_reader()

    stacked = {}
    for t in reader.time_values:
        reader.set_active_time_value(t)
        mesh = reader.read()  #<- Read the data
        internalMesh = mesh["internalMesh"]
        
        stacked[f"{t:.2f}"] = internalMesh.integrate_data()
    
    return stacked

def get_timestep_integrated(of: Case_Directory, t_target_days:float = -1):
    """
    Return an integrated array for a single timestep.
    """

    reader = of.get_vtk_reader()
    time_of_interest_index = argmin([abs(t_target_days - t/86400) for t in reader.time_values])
    time_of_interest = reader.time_values[time_of_interest_index]
    reader.set_active_time_value(time_of_interest)

    mesh = reader.read()  #<- Read the data
    internalMesh = mesh["internalMesh"]
    
    return internalMesh.integrate_data()

def to_numeric(s:str):
    '''Quick regex to extract the values of an openfoam dict entry'''
    numbers = findall(r"[-]?\d+[.]?\d*", s)
    return [float(f) for f in numbers] if len(numbers) > 1 else float(numbers[0])

def main():
    return 0

if __name__ == "__main__":
    exit(main())