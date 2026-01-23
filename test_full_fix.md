## Notes for Running the Full Test Workflow

1. **MultiQC container**
   - `nextflow.config` has been updated with a valid MultiQC image for all processes that require it.

2. **Test samplesheet**
   - `test/samplesheet_test.csv` has been updated so that folder paths point to the current test data location.

3. **GTDB database configuration**
   - After downloading the GTDB database, update the `gtdb_db` parameter in:
     - `conf/test_full.config`
     - `nextflow.config`
     - or set it as an argument `--gtdb_db`
   - Set the parameter to the correct path of the downloaded GTDB database (the **extracted** directory, not the original compressed .gz file).
   - There is a subfoler `skani` in the extracted directory. You will need to change it to `fastani`  
   ```bash 
   mv skani fastani
   ```

4. **GTDBTK_CLASSIFYWF**
   - GTDB-Tk uses TMPDIR by default, previously it ended up on NFS, now we force it to local disk → no more tempdir crash

5. **Running the full test workflow**
   - Use the `singularity` profile.
   - Disable the `midas2` process.
   - The following command runs a full test of all tools in the pipeline using the provided test dataset, assuming you're in the same folder where you cloned the repo. Leave `--gtdb_db` alone if it's the first time you run the pipeline:

   ```bash
   nextflow run main.nf (or uel3/nf-uno) \
     -profile singularity,test_full \
     --outdir <output_folder> \
     --skip_midas2 \
     --input test/samplesheet_test.csv \
     --gtdb_db <path to the extracted database> 
