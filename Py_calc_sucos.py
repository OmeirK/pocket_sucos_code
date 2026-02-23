## Adapt the sucos calculation code from the Runs N Poses publication

import argparse
import os
import pandas as pd
import numpy as np
from rdkit import Chem, DataStructs, RDConfig
from rdkit.Chem import AllChem, rdShapeAlign, rdShapeHelpers
from rdkit.Chem.FeatMaps import FeatMaps
from plinder.core import get_config
from plinder.core.scores import query_index
from plinder.core import PlinderSystem

parser = argparse.ArgumentParser()

parser.add_argument('--qcov_file', '-q', help='tsv file with pocket_qcov values')
parser.add_argument('--lig_sdf', '-sdf', help='.sdf file with the query ligand')
parser.add_argument('--outfile', '-o', help='Name of the output .tsv file')

args = parser.parse_args()

# Adapted from https://github.com/susanhleung/SuCOS
# Initialize feature factory for pharmacophore scoring
FDEF = AllChem.BuildFeatureFactory(
    os.path.join(RDConfig.RDDataDir, "BaseFeatures.fdef")
)

# Feature map parameters
FEAT_MAP_PARAMS = {k: FeatMaps.FeatMapParams() for k in FDEF.GetFeatureFamilies()}

# Feature types to keep for pharmacophore scoring
PHARMACOPHORE_FEATURES = (
    "Donor",
    "Acceptor",
    "NegIonizable",
    "PosIonizable",
    "ZnBinder",
    "Aromatic",
    "Hydrophobe",
    "LumpedHydrophobe",
)

def get_feature_map_score(
    mol_1: Chem.Mol,
    mol_2: Chem.Mol,
    score_mode: FeatMaps.FeatMapScoreMode = FeatMaps.FeatMapScoreMode.All,
) -> float:
    feat_lists = []
    for molecule in [mol_1, mol_2]:
        raw_feats = FDEF.GetFeaturesForMol(molecule)
        feat_lists.append([
            f for f in raw_feats if f.GetFamily() in PHARMACOPHORE_FEATURES
        ])

    feat_maps = [
        FeatMaps.FeatMap(feats=x, weights=[1] * len(x), params=FEAT_MAP_PARAMS)
        for x in feat_lists
    ]
    feat_maps[0].scoreMode = score_mode

    score = feat_maps[0].ScoreFeats(feat_lists[1])
    return score / min(feat_maps[0].GetNumFeatures(), len(feat_lists[1]))

def get_sucos_score(
    mol_1: Chem.Mol,
    mol_2: Chem.Mol,
    score_mode: FeatMaps.FeatMapScoreMode = FeatMaps.FeatMapScoreMode.All,
) -> float:
    fm_score = get_feature_map_score(mol_1, mol_2, score_mode)
    fm_score = np.clip(fm_score, 0, 1)

    protrude_dist = rdShapeHelpers.ShapeProtrudeDist(
        mol_1, mol_2, allowReordering=False
    )
    protrude_dist = np.clip(protrude_dist, 0, 1)

    return 0.5 * fm_score + 0.5 * (1 - protrude_dist)

def main():
    df = pd.read_csv(args.qcov_file, delimiter='\t')
    err_log = []

    outlines = ['query\ttarget\ttarget_rec_chain\tligand_chain\trelease_date\taln_evalue\tpocket_qcov\tsucos\tsucos_pocket']

    q_mol = Chem.MolFromMolFile(args.lig_sdf)

    for i, target in enumerate(df['target']):
        query = df['query'].iloc[i]
        rec_ch = df['target_rec_chain'].iloc[i]
        lig_ch = df['ligand_chain'].iloc[i]
        pocket_qcov = df['pocket_qcov'].iloc[i]
        rls_date = df['release_date'].iloc[i]
        evalue = df['aln_evalue'].iloc[i]
        u_mtx = df['rot_mtx'].iloc[i]
        t_vec = df['trans_vec'].iloc[i]

        #print(i, target, rls_date, pocket_qcov)

        #if pocket_qcov == 0:
        #    sucos = 0
        #    sucos_pocket = 0
        #    outlines.append(f'{query}\t{target}\t{rec_ch}\t{lig_ch}\t{rls_date}\t{evalue}\t{pocket_qcov}\t{sucos}\t{sucos_pocket}')
        #    continue

        plinder_system = PlinderSystem(system_id=target)
        target_sdfs = plinder_system.ligand_sdfs

        target_sdf = target_sdfs[lig_ch]
        t_mol = Chem.MolFromMolFile(target_sdf)
        #print(lig_ch, target_sdf, Chem.MolToSmiles(t_mol))

        rotation = np.array(list(map(float, u_mtx.split(','))))
        translation = np.array(list(map(float, t_vec.split(','))))

        conf = t_mol.GetConformer()
        coords = np.array([
            list(conf.GetAtomPosition(i))
            for i in range(t_mol.GetNumAtoms())
        ])
        rotated_coords = coords @ rotation.reshape(3, 3).T + translation
        for i in range(t_mol.GetNumAtoms()):
            conf.SetAtomPosition(i, rotated_coords[i])
        
        try:
            sucos = get_sucos_score(q_mol, t_mol)
            sucos_pocket = sucos*pocket_qcov
            outlines.append(f'{query}\t{target}\t{rec_ch}\t{lig_ch}\t{rls_date}\t{evalue}\t{pocket_qcov}\t{sucos}\t{sucos_pocket}')
        except Exception as e:
            err_log.append(f'SuCOS error for: {target} {lig_ch} {Chem.MolToSmiles(t_mol)}:\n')
            err_log.append(str(e)+ '\n')

        #print('\t', pocket_qcov, sucos, sucos)

        with open(args.outfile, 'w') as fo:
            fo.write('\n'.join(outlines))

        err_out = os.path.basename(args.outfile).split('.')
        err_out = '.'.join(err_out[:-1]) + '.err'

        with open(f'{os.path.dirname(os.path.abspath(args.outfile))}/{err_out}', 'w') as fo:
            fo.write(''.join(err_log))



if __name__=='__main__':
    main()
