By default, simulation results are saved in this folder.

Generated output files will have the following folder structure:

output   
├── sim_name   
│   ├── drs_<drug_name>   
│   │   ├── drs_<drug_name>_<replicate_id>   
│   │   |    ├── <drug_name>_<dose>   
│   │   |    |  ├── output_g1.pkl   
│   │   |    |  ├── output_g1.pkl   
│   │   |    |  ├── output_g1.pkl   
│   │   |    |  ├── output_g1.pkl   
│   │   |    |  ├── output_g1.pkl   

RootFolder   
├── Folder1   
│   ├── SubFolder1   
│   │   ├── Level1   
│   │   │   ├── Level2   
│   │   │   │   ├── File1.txt   
│   │   │   │   ├── File2.txt   
│   │   │   │   ├── File3.txt   
│   │   │   │   ├── File4.txt   
│   │   │   │   └── File5.txt   
