
from Boolean2Mechanism import importModel, importLibrary, importLabels, \
    MotifBuilder, Combine_and_Build


# ----- Simple example (A = B and C) -----------------------------------------------

model_nodes = importModel('Example_model.py')
library = importLibrary('B2M.mrl')
importLabels('Example_labels.csv', library, model_nodes)
MotifBuilder(model_nodes, library)

# print
# print '--------------------------------------'
# for each in model_nodes:
#     print
#     print each
#     for item in model_nodes[each].motifs:
#         print
#         for every in item:
#             print every
#
# print
# print '--------------------------------------'

Combine_and_Build(model_nodes, library, top='all')
