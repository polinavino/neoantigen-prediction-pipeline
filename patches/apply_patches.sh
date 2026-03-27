#!/bin/bash
PVAC_LIB=$(python3 -c "import pvactools; import os; print(os.path.dirname(pvactools.__file__))")/lib

sed -i "" "s/asparagine_proline_bond_count/aspartate_proline_bond_count/g" $PVAC_LIB/calculate_manufacturability.py
sed -i "" "s|http://tools-cluster-interface.iedb.org/tools_api/mhci/|https://tools-cluster-interface.iedb.org/tools_api/mhci/|g" $PVAC_LIB/prediction_class.py
sed -i "" "s|http://tools-cluster-interface.iedb.org/tools_api/mhcii/|https://tools-cluster-interface.iedb.org/tools_api/mhcii/|g" $PVAC_LIB/prediction_class.py
sed -i "" "s/requests.post(self.url, data=data)/requests.post(self.url, data=data, timeout=300)/g" $PVAC_LIB/prediction_class.py
sed -i "" 's/arguments = \["python", script/arguments = ["python3", script/' $PVAC_LIB/prediction_class.py

echo "Patches applied successfully"
