- [ ] Remove all camel casing function names/vars, change to underscores and lowercase.
- [X] Reorder amino acids alphabetically.
- [X] Add normalize to moran and geary auto
- [X] For quasi seq order you should be able to pass in name of distance matrix file with or without .json
- [X] For each descriptor, check valid amino acids in seq, if not then raise custom error.
- [ ] Add dimensions of each descriptor to readme and docs.
- [ ] Round SOCN to 3 d.p
- [X] Unit test the data type for each column in all descriptors.
- [X] Remove .empty tests from unit tests as validating shape of DF will test for emptiness.
- [X] Add 0 to singualr descriptor columns, e.g polarizability_CTD_C_1 -> polarizability_CTD_C_01 
- [X] Change sec_struct to secondary_struct.
- [ ] Reduce number of tests by iterating over list of protein seqs.
- [X] Test dtypes of output dataframe -> test_autocorrelation 
- [X] Add shape to comment on testing shape unittests.
- [X] Mention lag is similar to gap between 2 amino acids.
- [X] Go through test_quasi file, double checking correct values.
- [X] Append distance matrix to SOCN & Quasi columns, SW or G.
- [X] Change quasi sequence order -> sequence_order.
- [X] Calculate all SOCN, for both matrices, append to single output df.
- [X] SOCN done, quasi done.
- [ ] Reread descriptor comments and explanations.
- [X] Change SOCNUm to SOCN.
- [X] Pseudo AAC has to explicitly use hydrophoobicity, hydrophilicity and side-chain/residue mass values. Can't find corresponding values in aaindex so just hard code them in.
- [X] If no properties input to pseudo or amp comp funcs then use hydo, hydrophi, residue by default. Accept list of aaindex1 codes, if str input then cast to list.
- [X] Uppercase sequence on input, remove whitespace. 
- [X] Move references to top of each module.
- [ ] Add the equations to comments in some of the descriptors? (https://bmcbioinformatics.biomedcentral.com/articles/10.1186/s12859-015-0554-8#Sec10)
- [X] Input property in CTD funcs can be used with closeness function.
- [X] Double check functions that use aa_composition values, aa_comp func returns series rather than dict.
- [X] Rather than iterate over range of lags, use different lag in each sequence test.
- [X] Change max_lag to lag 
- [ ] Create demo on Notebook.
- [X] Add descriptor abbreviations to each functiosn comments, change abbreviations of Pseudo AAComp -> PAAComp.
- [X] Add references to readme text.
- [X] In readme, add output of each function below its usage.
- [ ] Add reference numbers to comments in descriptor functions - double check existing ones are correct.
- [X] Add lag and weight param validation to sequence order module.
- [X] Change QSOrder to QSO.
- [ ] Rewrite APAAComp descriptor comments to mention its dimensions change with lamda. 
- [X] For all functions that have lag in them:
        #raise value error if int cant be parsed from input lag
        try:
            lag = int(lag)
        except:
            raise ValueError("Invalid lag value input, integer cannot be parsed from {}".format(lag))
- [ ] Add logo/image to main readme.
- [ ] Add emojis to readme.
- [ ] Add releases.
