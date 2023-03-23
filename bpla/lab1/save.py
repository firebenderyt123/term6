import os
import pandas as pd


def save_array(data, titles=[], path="out/states.xlsx"):

    if not os.path.exists(os.path.dirname(path)):
        os.makedirs(os.path.dirname(path))

    default_titles = [
        'x',
        'y',
        'z',
        'Vx',
        'Vy',
        'Vz',
        'psi',
        'teta',
        'gamma',
        'Om1',
        'Om2',
        'Om3',
        'omega0',
        'omega1',
        'omega2',
        'omega3',
        'omega4',
        'omega5',
        'eps0',
        'eps1',
        'eps2',
        'eps3',
        'eps4',
        'eps5',
        'F0_x',
        'F0_y',
        'F0_z',
        'F1_x',
        'F1_y',
        'F1_z',
        'F2_x',
        'F2_y',
        'F2_z',
        'F3_x',
        'F3_y',
        'F3_z',
        'F4_x',
        'F4_y',
        'F4_z',
        'F5_x',
        'F5_y',
        'F5_z',
        'Lam_w',
        'Lam_x',
        'Lam_y',
        'Lam_z',
        'Fe_x',
        'Fe_y',
        'Fe_z',
        'v_zsk_x',
        'v_zsk_y',
        'v_zsk_z',
        'Fs_zsk_x',
        'Fs_zsk_y',
        'Fs_zsk_z',
        'Fs_x',
        'Fs_y',
        'Fs_z',
        'Mg_zsk_x',
        'Mg_zsk_y',
        'Mg_zsk_z',
        'dH_x',
        'dH_y',
        'dH_z',
        'MFe_x',
        'MFe_y',
        'MFe_z',
        'MAe_x',
        'MAe_y',
        'MAe_z'
    ]

    if len(titles) != len(default_titles):
        titles = default_titles

    df = pd.DataFrame(data, columns=titles)

    df.to_excel(path, index=False)
