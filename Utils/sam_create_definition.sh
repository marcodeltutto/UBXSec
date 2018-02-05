#! /bin/bash

# Create 'group1' definition with files from prod_v04_26_04_05 and version v04%
samweb create-definition mdeltutt_prod_reco2_extbnb_v9_mcc8_group1 "defname: prod_reco2_extbnb_v9_mcc8 and isdescendantof: ( file_type data and data_tier raw and ub_project.name swizzle_trigger_streams and ub_project.stage mergeext_bnb and ub_project.version prod_v04_26_04_05 and version v04%)"

# Create 'group2' definition with files from prod_v04_26_04_05 and version v05%
samweb create-definition mdeltutt_prod_reco2_extbnb_v9_mcc8_group2 "defname: prod_reco2_extbnb_v9_mcc8 and isdescendantof: ( file_type data and data_tier raw and ub_project.name swizzle_trigger_streams and ub_project.stage mergeext_bnb and ub_project.version prod_v04_26_04_05 and version v05%)"



-bash-4.1$ samweb list-files --summary  "defname:prod_reco_optfilter_bnb_v11_unblind_mcc8a"
File count: 258
Total size: 16701709083
Event count:  912
-bash-4.1$ samweb list-files --summary  "defname:prod_reco_optfilter_bnb_v11_unblind_mcc8b"
File count: 3852
Total size: 3453288850193
Event count:  191131
-bash-4.1$ samweb list-files --summary  "defname:prod_reco_optfilter_bnb_v11_unblind_mcc8c"
File count: 0
Total size: None
Event count:  None
-bash-4.1$ samweb list-files --summary  "defname:prod_reco_optfilter_bnb_v11_unblind_mcc8"
File count: 4110
Total size: 3469990559276
Event count:  192043
-bash-4.1$ samweb list-files --summary  "defname:prod_reco_optfilter_extbnb_v11_mcc8a"
File count: 6969
Total size: 9412118136113
Event count:  521271
-bash-4.1$ samweb list-files --summary  "defname:prod_reco_optfilter_extbnb_v11_mcc8b"
File count: 11312
Total size: 8420718763202
Event count:  463273
-bash-4.1$ samweb list-files --summary  "defname:prod_reco_optfilter_extbnb_v11_mcc8"
File count: 18281
Total size: 17832836899315
Event count:  984544
