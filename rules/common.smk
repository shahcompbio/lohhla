import os
import numpy as np
import pandas as pd
import warnings

class WgsRunInfo:
    def __init__(self, config):
        self.metadata = pd.read_csv(config['cluster_paths_sheet'])
        self.metadata = self.metadata[self.metadata['isabl_application'].isin(['WGS-ALIGNMENT', 'WGS-POLYSOLVER'])]
        sample_ids = set(self.metadata.query('sample_category == "TUMOR"')['isabl_sample_id'].unique())
        self.paths = self.metadata.set_index(['isabl_sample_id', 'isabl_application', 'result_type'])['result_filepath'].to_dict()
        self.paired_dict = self.metadata.query('isabl_sample_id != \'Unspecified Aliquot\'').set_index(['isabl_patient_id', 'sample_category'])['isabl_sample_id'].to_dict()

wgsruninfo = WgsRunInfo(config)
