import importlib.util
import unittest

@unittest.skipUnless(all(importlib.util.find_spec(p) for p in ['numpy','pandas','scipy']), 'Scientific numeric dependencies run on Colab')
class ExtensionTests(unittest.TestCase):
    def test_replacement_control_changes_difference_without_affecting_other_AA(self):
        import pandas as pd
        import numpy as np
        from paper1.panels import replace_panel,DRUGS
        ids=['FKS1:V641W','FKS1:L642K','FKS1:F639A']
        data=pd.DataFrame({'aa_variant_id':ids,'fitness_none':[1.,2.,3.]})
        for d in DRUGS:data['fitness_'+d]=[3.,4.,5.];data['drug_specific_'+d]=2.
        data=data.set_index('aa_variant_id',drop=False)
        repl=pd.DataFrame([{'aa_variant_id':a,'compound':c,'s':10. if c=='none' else 12.} for a in ids[:2] for c in ['none',*DRUGS]])
        result=replace_panel(data,repl)
        np.testing.assert_equal(result.loc[ids[2],data.columns].to_numpy(),data.loc[ids[2]].to_numpy())
        self.assertEqual(result.loc[ids[0],'drug_specific_micafungin'],2.)
        with self.assertRaises(ValueError):replace_panel(data,pd.concat([repl,repl.iloc[:1]]))

    def test_conditional_draws_preserve_position_and_hotspot_and_undefined_is_NA(self):
        import numpy as np
        from paper1.metrics import signed_rho,draw_positions
        self.assertTrue(np.isnan(signed_rho([1,1,1],[1,2,3])[0]))
        self.assertAlmostEqual(signed_rho([1,2,2,3],[4,3,3,1])[0],-1.)
        hs={1:'HS1',2:'HS1',3:'HS2'};a=draw_positions([1,2,3],hs,B=25);b=draw_positions([1,2,3],hs,B=25)
        np.testing.assert_array_equal(a,b)
        self.assertTrue((a[:,-1]==3).all());self.assertTrue(np.isin(a[:,:2],[1,2]).all())

    def test_Q2_uses_training_mean_and_retains_negative(self):
        import pandas as pd
        from paper1.metrics import prediction_metrics
        f=pd.DataFrame({'aa_position':[1,1,2],'y':[1.,2.,3.],'m_OOF':[9.,9.,9.],'m_train_mean_baseline':[0.,0.,0.]})
        r=prediction_metrics(f,'y','m')
        self.assertAlmostEqual(r['pooled_Q2'],1-(64+49+36)/14)
        self.assertLess(r['equal_position_Q2'],0)
        f['m_train_mean_baseline']=f.y
        import math
        self.assertTrue(math.isnan(prediction_metrics(f,'y','m')['pooled_Q2']))

    def test_constant_and_rank_deficient_features_are_retained_and_audited(self):
        import pandas as pd
        from paper1.comparison import matched_lopo
        rows=[dict(aa_variant_id=f'{p}:{i}',aa_position=p,gemme_primary=1.,esm1v_ensemble=float(i),evo2_primary_AA_score=2.*i,y=p+i) for p in range(4) for i in range(3)]
        _,audit,_=matched_lopo(pd.DataFrame(rows),'y')
        self.assertTrue(all(x['rank_deficient'] for x in audit))
        self.assertTrue(all('gemme_primary' in x['constant_training_features'] for x in audit))
        self.assertTrue(all(len(x['singular_values'])==len(x['features'])+1 for x in audit))

    def test_paired_effect_is_median_of_differences_and_estimability_is_shared(self):
        import numpy as np
        import pandas as pd
        from paper1.metrics import rank_summaries,draw_positions,paired_positions
        records=pd.DataFrame([{'aa_position':p,'model':m,'rho':v} for p,a,b in [(1,.8,.7),(2,.9,-.8),(3,-.2,-.1),(4,np.nan,.5)] for m,v in [('A',a),('B',b)]])
        pos=paired_positions(records,['A','B']);self.assertEqual(pos.tolist(),[1,2,3])
        draws=draw_positions(pos,{1:'HS1',2:'HS1',3:'HS2',4:'HS2'},B=50)
        summaries,_=rank_summaries(records,['A','B'],[('A','B')],draws)
        difference=next(x['point'] for x in summaries if x['contrast_or_model']=='A minus B')
        self.assertAlmostEqual(difference,.1)
        self.assertNotAlmostEqual(difference,np.median([.8,.9,-.2])-np.median([.7,-.8,-.1]))

    def test_position_deletion_and_HS2_removal_rebuild_training(self):
        import numpy as np
        import pandas as pd
        from paper1.analysis import analyse
        frame=pd.DataFrame([dict(aa_variant_id=f'{p}:{i}',aa_position=p,hotspot='HS2' if p==4 else 'HS1',
              gemme_primary=float(i+p),esm1v_ensemble=float(i*i-p),evo2_primary_AA_score=float(i*p),
              y=float(i+p) if p!=4 else float(1000+100*i)) for p in range(5) for i in range(4)])
        full=analyse(frame,'y',uncertainty=False)
        reduced=analyse(frame[frame.hotspot!='HS2'],'y',uncertainty=False)
        for audit in reduced['audit']:
            self.assertNotIn(4,audit['train_positions']);self.assertNotEqual(4,audit['heldout_position'])
        a=full['OOF'].query('aa_position==1').GEMME_ESM_OOF.to_numpy()
        b=reduced['OOF'].query('aa_position==1').GEMME_ESM_OOF.to_numpy()
        self.assertTrue(np.any(np.abs(a-b)>1e-6))

if __name__=='__main__':unittest.main()
