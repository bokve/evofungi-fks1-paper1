import importlib.util
import unittest

from paper1.comparison import ensemble_mean, matched_lopo
from paper1.esm_contract import CHECKPOINTS


class EnsembleTests(unittest.TestCase):
    def test_missing_checkpoint_and_different_mask_are_rejected(self):
        members = {name:{'AA':float(i)} for i,name in enumerate(CHECKPOINTS)}
        self.assertEqual(ensemble_mean(members), {'AA':2.0})
        with self.assertRaises(ValueError):
            ensemble_mean({name:members[name] for name in CHECKPOINTS[:-1]})
        members[CHECKPOINTS[0]] = {'other_AA':1.0}
        with self.assertRaises(ValueError):
            ensemble_mean(members)


@unittest.skipUnless(importlib.util.find_spec('numpy') and importlib.util.find_spec('pandas'), 'Comparison CPU dependencies remain on Colab')
class GroupedComparisonTests(unittest.TestCase):
    def test_heldout_y_does_not_enter_fit_and_masks_are_shared(self):
        import numpy as np
        import pandas as pd
        rows = []
        for pos in range(1,5):
            for alt in range(3):
                rows.append({'aa_variant_id':f'test:{pos}:{alt}','aa_position':pos,'gemme_primary':pos+alt/3,'esm1v_ensemble':alt-pos/4,'evo2_primary_AA_score':pos*alt+1,'y':pos+alt})
        frame = pd.DataFrame(rows)
        frame.loc[0,'evo2_primary_AA_score'] = np.nan
        a,audit,dropped = matched_lopo(frame,'y')
        changed = frame.copy()
        changed.loc[changed.aa_position==1,'y'] = 999999
        b,_,_ = matched_lopo(changed,'y')
        self.assertEqual(dropped,['test:1:0'])
        for model in ['GEMME_ESM','GEMME_ESM_Evo2']:
            np.testing.assert_allclose(a.loc[a.aa_position==1,model+'_OOF'],b.loc[b.aa_position==1,model+'_OOF'])
        for record in audit:
            self.assertNotIn(record['heldout_position'],record['train_positions'])
            self.assertFalse(set(record['train_AA_ids']) & set(record['test_AA_ids']))
        self.assertEqual([r['test_AA_ids'] for r in audit[:4]], [r['test_AA_ids'] for r in audit[4:]])


if __name__=='__main__':unittest.main()
