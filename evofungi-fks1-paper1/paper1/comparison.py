"""Fixed matched-mask, position-held-out comparisons; no tuning or imputation."""
import math
import statistics

from .esm_contract import CHECKPOINTS


def ensemble_mean(members):
    if set(members) != set(CHECKPOINTS):
        raise ValueError('All five fixed ESM-1v checkpoints are required')
    keys = set(members[CHECKPOINTS[0]])
    if not keys or any(set(members[name]) != keys for name in CHECKPOINTS):
        raise ValueError('Checkpoint AA masks differ')
    result = {}
    for key in sorted(keys):
        values = [float(members[name][key]) for name in CHECKPOINTS]
        if not all(math.isfinite(x) for x in values):
            raise ValueError('Nonfinite checkpoint score; no partial ensemble')
        result[key] = statistics.mean(values)
    return result


def matched_lopo(frame, endpoint):
    """Fixed OLS: GEMME+ESM versus GEMME+ESM+Evo2, one shared AA mask.

    Returns OOF predictions and train-fold preprocessing/rank audit, without
    calculating inferential intervals or claiming zero-shot supervised results.
    """
    import numpy as np
    required = ['aa_variant_id', 'aa_position', 'gemme_primary', 'esm1v_ensemble', 'evo2_primary_AA_score', endpoint]
    if not set(required).issubset(frame.columns):
        raise ValueError('Explicit AA/position/predictor/endpoint columns required')
    if frame['aa_variant_id'].duplicated().any():
        raise ValueError('Duplicate AA phenotype; codon/replicate rows are not independent variants')
    features = ['gemme_primary', 'esm1v_ensemble', 'evo2_primary_AA_score']
    values = frame[features + [endpoint]].to_numpy(dtype=float)
    keep = np.isfinite(values).all(axis=1)
    data = frame.loc[keep, required].copy().reset_index(drop=True)
    dropped = frame.loc[~keep, 'aa_variant_id'].tolist()
    groups = data['aa_position'].to_numpy()
    positions = sorted(set(groups))
    if len(positions) < 3:
        raise ValueError('At least three physical positions required for grouped comparison')
    audits = []
    for label, columns in [('GEMME_ESM', features[:2]), ('GEMME_ESM_Evo2', features)]:
        predictions = np.full(len(data), np.nan)
        baseline = np.full(len(data), np.nan)
        for position in positions:
            train = groups != position
            test = ~train
            xtrain = data.loc[train, columns].to_numpy(dtype=float)
            xtest = data.loc[test, columns].to_numpy(dtype=float)
            ytrain = data.loc[train, endpoint].to_numpy(dtype=float)
            mean = xtrain.mean(axis=0)
            sd = xtrain.std(axis=0, ddof=0)
            scale = np.where(sd == 0, 1.0, sd)
            design = np.column_stack([np.ones(len(xtrain)), (xtrain-mean)/scale])
            heldout = np.column_stack([np.ones(len(xtest)), (xtest-mean)/scale])
            coefficients, _, rank, singular = np.linalg.lstsq(design, ytrain, rcond=None)
            predictions[test] = heldout @ coefficients
            baseline[test] = ytrain.mean()
            condition = float(singular[0]/singular[-1]) if singular[-1]>0 else None
            audits.append({'model':label,'heldout_position':int(position),'train_positions':[int(p) for p in positions if p != position],'train_AA_ids':data.loc[train,'aa_variant_id'].tolist(),'test_AA_ids':data.loc[test,'aa_variant_id'].tolist(),'features':columns,'train_means':mean.tolist(),'train_scales':scale.tolist(),'constant_training_features':[c for c,s in zip(columns,sd) if s==0], 'singular_values':singular.tolist(),'condition_number':condition,'rank':int(rank),'rank_deficient':bool(rank < design.shape[1]),'tuning_performed':False})
        if not np.isfinite(predictions).all():
            raise ValueError('Missing/nonfinite held-out predictions')
        data[label + '_OOF'] = predictions
        data[label + '_train_mean_baseline'] = baseline
    return data, audits, dropped
