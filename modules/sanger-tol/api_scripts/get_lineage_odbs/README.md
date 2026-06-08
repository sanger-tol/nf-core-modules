# get_odbs.py

Originally written by Matthieu Muffato (mm49)

Modified by Damon-Lee Pointon (dp24)

## Folder Structure

```
.
├── resources
│   └── usr
│       └── assets
│       │   ├── mapping_taxids-busco_dataset_name.eukaryota_odb10.2019-12-16.txt
│       │   └── mapping_taxids-busco_dataset_name.eukaryota_odb12.2025-01-15.txt
│       └── bin
│           └── get_odbs.py
├── tests
│   ├── main.nf.test
│   ├── main.nf.test.snap
│   └── nextflow.config
│
├── environment.yml
├── README.md
├── main.nf
└── meta.yml
```

## Sources

ODB mapping files are provided by the [BUSCO](https://busco.ezlab.org/) project at the `busco-data.ezlab.org` data repo:

ODB files must be the merged collection of eukaryota, bacteria and archaea. For example:

To make the `odb10` merged mapping file, the below files must be merged into a single file such as `mapping_files/odb10_mapping_file.txt`:

- https://busco-data.s3.amazonaws.com/placement_files/mapping_taxids-busco_dataset_name.archaea_odb10.2019-12-16.txt.tar.gz"
- https://busco-data.s3.amazonaws.com/placement_files/mapping_taxids-busco_dataset_name.bacteria_odb10.2019-12-16.txt.tar.gz"
- https://busco-data.s3.amazonaws.com/placement_files/mapping_taxids-busco_dataset_name.eukaryota_odb10.2019-12-16.txt.tar.gz"

To make the `odb12` merged mapping file, the below files must be merged into a single file such as `mapping_files/odb12_mapping_file.txt`:

- https://busco-data.s3.amazonaws.com/placement_files/mapping_taxids-busco_dataset_name.archaea_odb12.2024-11-15.txt.tar.gz"
- https://busco-data.s3.amazonaws.com/placement_files/mapping_taxids-busco_dataset_name.bacteria_odb12.2024-11-15.txt.tar.gz"
- https://busco-data.s3.amazonaws.com/placement_files/mapping_taxids-busco_dataset_name.eukaryota_odb12.2024-11-15.txt.tar.gz"
