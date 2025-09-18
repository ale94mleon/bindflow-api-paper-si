from toff import Parameterize
import logging
logging.basicConfig(level=logging.INFO, format='%(asctime)s - %(levelname)s - %(name)s - %(message)s')

ff = {
    "gaff": "gaff-2.11",
    "openff": "openff_unconstrained-2.0.0.offxml"
}

for force_field_type, force_field_code in ff.items():
    print(force_field_type)

    parameterizer = Parameterize(force_field_type=force_field_type,
                                force_field_code=force_field_code,
                                ext_types=['gro', 'top'],
                                overwrite=True, safe_naming_prefix='y',
                                out_dir=force_field_code)

    parameterizer("OA.sdf", mol_resi_name="HOST")
