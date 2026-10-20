# augmenta_config.yaml e' il dump dei parametri effettivi (i default) scritto
# da AugmentA in <cartella mesh>/_configs/, tenuto come registrazione. Per
# rifare la run con quei parametri esatti si passa --config_dir.
cd /headless/augmenta
python main.py --mesh /headless/data/ssm_mean_LA.vtk --closed_surface 0 --open_orifices 0 --atrium LA --SSM_fitting 0 --resample_input 0 --normals_outside 1 --add_bridges 0 --MRI 0 --ofmt vtk
