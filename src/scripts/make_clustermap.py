import sys
sys.path.append('/lustre/groups/ml01/workspace/samantha.bening/Bachelor/')
from importlib import reload
import genereporter.sample_pipeline as spModule
reload(spModule)
sp = spModule.SamplePipeline(wdir="/lustre/groups/ml01/workspace/samantha.bening/Bachelor/", adata="data2/veo_ibd_balanced.h5ad")

GOI = 'CASP8' # set the GOI!!!

sp.pl_sample_celltype(GOI, z_score=True)

print("Done!")