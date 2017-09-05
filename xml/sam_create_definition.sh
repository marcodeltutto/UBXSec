#! /bin/bash

# Create 'group1' definition with files from prod_v04_26_04_05 and version v04%
samweb create-definition mdeltutt_prod_reco2_extbnb_v8_mcc8_group1 "defname: prod_reco2_extbnb_v8_mcc8 and isdescendantof: ( file_type data and data_tier raw and ub_project.name swizzle_trigger_streams and ub_project.stage mergeext_bnb and ub_project.version prod_v04_26_04_05 and version v04%)"

# Create 'group2' definition with files from prod_v04_26_04_05 and version v05%
samweb create-definition mdeltutt_prod_reco2_extbnb_v8_mcc8_group2 "defname: prod_reco2_extbnb_v8_mcc8 and isdescendantof: ( file_type data and data_tier raw and ub_project.name swizzle_trigger_streams and ub_project.stage mergeext_bnb and ub_project.version prod_v04_26_04_05 and version v05%)"
