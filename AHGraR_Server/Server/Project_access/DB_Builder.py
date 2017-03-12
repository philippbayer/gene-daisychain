# Provides project build functions to AHGraR-Server
# All functions return either a String or Null
# Functions directly accessible by user query always return a string via socket connection
import os
import subprocess
from itertools import islice
from CSV_creator.annotation_to_csv import AnnoToCSV
from CSV_creator.protein_cluster_to_csv import ClusterToCSV
from Parser.FASTA_parser import FastaParser
#import Parser.GFF3_parser_gffutils
from Parser.GFF3_parser_gffutils_v2 import GFF3Parser_v2
from random import choice
from neo4j.v1 import GraphDatabase, basic_auth
import time
import json

class DBBuilder:
    def __init__(self, main_db_connection, task_manager, send_data, ahgrar_config):
        self.main_db_conn = main_db_connection
        self.task_mngr = task_manager
        self.send_data = send_data
        self.ahgrar_config = ahgrar_config

    # Close connection to main-DB
    def close_connection(self):
        self.main_db_conn.close()

    # Reply to request send from a user app
    # User_request is a list produced by the "_" split command
    # e.g. [STAT, ProjectID, TaskID1, TaskID2]
    def evaluate_user_request(self, user_request):
        # Set GFF3 parser for some or all GFF3 files in a project
        print(user_request)
        print(len(user_request))
        if user_request[0] == "GFF3" and len(user_request) == 7 and user_request[1].isdigit():
            # Call format: ProjectID, parent_feat, sub_feature, name_attr, descr_attr
            self.set_gff3_parser(user_request[1],user_request[2], user_request[3],
                                 user_request[4], user_request[5], user_request[6])
        # Build the neo4j-based project database from the previously added files
        elif user_request[0] == "DB" and len(user_request) == 2 and user_request[1].isdigit():
            self.build_db(user_request[1])
        else:
            self.send_data("-3")

    # Build a neo4j graph database of the genes and proteins in the project data files
    # Two types of input for each species are possible:
    # (1.) Protein-fasta plus Gene annotation in GFF3 or CVS format
    # (2.) Nucleotide-fasta plus Gene annotation in GFF3 format
    # For (2.) a Protein-fasta is generated from the nucleotide sequence and then GFF3 gene annotation
    # Afterwards, the approach is the same as for (1.)
    # The build_db functions consists of multiple parts
    # After finishing each part, the task status is updated.
    # If one part fails, an error message is stored in the task result output.
    # (1.) Check for valid file combinations/input
    # (2.)
    def build_db(self, proj_id):
        # First, define a new task for the db build and return task-id to user
        task_id = self.task_mngr.define_task(proj_id, "Building project DB")
        # Send task-id to user
        self.send_data(task_id)
        # Then, retrieve the list of files in this project (excluding hidden files)
        # Only GFF3 files have the anno_mapping and feat_hierarchie field
        # For all other files, the return value for this field is None
        self.task_mngr.set_task_status(proj_id, task_id, "Collecting files")
        file_list = list(self.main_db_conn.run("MATCH(proj:Project)-[:has_files]->(:File_Manager)-[:file]->(file:File) "
                              "WHERE ID(proj)={proj_id} AND file.hidden = 'False' "
                              "RETURN file.filename, file.filetype, file.species, file.variant, file.parent_feat, "
                                            "file.sub_features, file.name_attr, file.desc_attr ORDER BY file.filename",
                          {"proj_id":int(proj_id)}))
        # Convert file_list into a dictionary:
        file_dict = {}
        # Keys are (Species, Variant) and entries are a list of files
        for file in file_list:
            file_dict[(file["file.species"],file["file.variant"])] = []
        for file in file_list:
            file_dict[(file["file.species"],file["file.variant"])].append(
                (file["file.filename"],file["file.filetype"],file["file.parent_feat"],file["file.sub_features"],
                 file["file.name_attr"],file["file.desc_attr"]))

        # Check if each entry in the database consists of exactly two files, one fasta and one annotation file
        # If not, remove that entry from the database
        for species in file_dict.keys():
            file_types = [item[1] for item in file_dict[species]]
            if len(file_types) != 2 or "annotation" not in file_types or "genome" not in file_types:
                del file_dict[species]

        # Initialize the annotation to csv format parser
        self.task_mngr.set_task_status(proj_id, task_id, "Parsing annotation data")
        # Initialize CSV parser:  CSV output directory, gene transcript and translation output files
        # All genes from all files are combined into one set of output files
        anno_to_csv_parser = AnnoToCSV(os.path.join("Projects", proj_id, "CSV"),
                                       os.path.join("Projects", proj_id, "BlastDB", "transcripts.faa"),
                                       os.path.join("Projects", proj_id, "BlastDB", "translations.faa"))
        # Then convert every annotation file into a Neo4j-specific CSV file format
        for species in file_dict.keys():
            try:
                # Retrieve name of annotation and genome file. Sort file list alphanumerical. Since
                # annotation < genome the annotation file is nr. 0, the genome file nr. 1
                anno_file = sorted(file_dict[species], key=lambda x: x[1])[0]
                genome_file = sorted(file_dict[species], key=lambda x: x[1])[1]
                print(anno_file)
                print(genome_file)
                anno_to_csv_parser.create_csv(os.path.join("Projects", proj_id, "Files", anno_file[0]),
                                              os.path.join("Projects", proj_id, "Files", genome_file[0]),
                                              anno_file[2], anno_file[3],anno_file[4],anno_file[5])
            except (IndexError, KeyError):
                self.task_mngr.set_task_status(proj_id, task_id, "Failed")
                self.task_mngr.add_task_results(proj_id, task_id, "Failed: Annotation parsing")
                return
        # Build the BLAST databases using BLAST+ makeblastdb
        # Define File folder path:
        self.task_mngr.set_task_status(proj_id, task_id, "Building Blast+ DB")
        BlastDB_path = os.path.join("Projects", str(proj_id), "BlastDB")
        makeblastdb_path = os.path.join(self.ahgrar_config["AHGraR_Server"]["blast+_path"], "makeblastdb")
        # Create transcript blast DB
        subprocess.run(
            [makeblastdb_path, "-dbtype", "nucl", "-in", os.path.join(BlastDB_path, "transcripts.faa"),
             "-parse_seqids", "-hash_index", "-out", os.path.join(BlastDB_path, "transcript_db")], check=True)
        subprocess.run(
            [makeblastdb_path, "-dbtype", "prot", "-in", os.path.join(BlastDB_path, "translations.faa"),
             "-parse_seqids", "-hash_index", "-out", os.path.join(BlastDB_path, "translation_db")], check=True)
        # Perform an all vs all blastn search
        self.task_mngr.set_task_status(proj_id, task_id, "All vs. all BlastN")
        blastn_path = os.path.join(self.ahgrar_config["AHGraR_Server"]["blast+_path"], "blastn")
        cpu_cores = self.ahgrar_config["AHGraR_Server"]["cpu_cores"]
        subprocess.run(
            [blastn_path, "-query", os.path.join(BlastDB_path, "transcripts.faa"), "-db",
             os.path.join(BlastDB_path, "transcript_db"), "-outfmt", "6 qseqid sseqid evalue pident",
                                         "-out", os.path.join(BlastDB_path, "transcripts.blastn"), "-evalue", "0.05",
                            "-num_threads", cpu_cores, "-parse_deflines"])
        # Perform an all vs all blastp search
        self.task_mngr.set_task_status(proj_id, task_id, "All vs. all BlastP")
        blastp_path = os.path.join(self.ahgrar_config["AHGraR_Server"]["blast+_path"], "blastp")
        subprocess.run(
            [blastp_path, "-query", os.path.join(BlastDB_path, "translations.faa"), "-db",
             os.path.join(BlastDB_path, "translation_db"), "-outfmt", "6 qseqid sseqid evalue pident",
             "-out", os.path.join(BlastDB_path, "translations.blastp"), "-evalue", "0.05",
             "-num_threads", cpu_cores, "-parse_deflines"])
        # Extract sequence match identity from blast result files
        # Create new blastn/blastp result files lacking the percent match ID column (ABC files)
        # Dump dict with geneID/geneID/PercentMatch and protID/protID/PercentMatch as json
        gene_gene_percentID = {}
        with open(os.path.join(BlastDB_path, "transcripts.blastn"), "r") as nt_blast_file:
            with open(os.path.join(BlastDB_path, "transcripts.abc"), "w") as  nt_blast_abc_file:
             for line in nt_blast_file:
                    line = line.split("\t")
                    gene_gene_percentID[(line[0],line[1])]=line[3]
                    nt_blast_abc_file.write("\t".join(line[:3]))
        with open(os.path.join(BlastDB_path, "transcripts_pid.json"), 'w') as dict_dump:
            json.dump(gene_gene_percentID, dict_dump)
        prot_prot_percentID = {}
        with open(os.path.join(BlastDB_path, "translations.blastp"), "r") as prot_blast_file:
            with open(os.path.join(BlastDB_path, "translations.abc"), "w") as  prot_blast_abc_file:
             for line in prot_blast_file:
                    line = line.split("\t")
                    prot_prot_percentID[(line[0],line[1])]=line[3]
                    prot_blast_abc_file.write("\t".join(line[:3]))
        with open(os.path.join(BlastDB_path, "translations_pid.json"), 'w') as dict_dump:
            json.dump(prot_prot_percentID, dict_dump)
        return

        # Cluster all-vs.-all BlastP results into protein homology groups
        self.task_mngr.set_task_status(proj_id, task_id, "Cluster BlastP results")
        # 1. Concatenate all BlastP Results into one "ABC" file
        print("1. Concatenate all BlastP Results into one ABC file")
        with open(os.path.join(BlastDB_path, "blastp.abc"), 'w') as combined_file:
            for species in file_dict:
                prot_fasta_file_name = sorted(file_dict[species], key=lambda x: x[1])[1][0]
                prot_blastp_file_name = prot_fasta_file_name[:prot_fasta_file_name.rfind(".")]+".blastp"
                with open(os.path.join(BlastDB_path, prot_blastp_file_name),"r") as blastp_file:
                    for line in blastp_file:
                        combined_file.write(line)
        # 2. Convert ABC file into a network and dictionary file.
        print("2. Convert ABC file into a network and dictionary file.")
        subprocess.run(["mcxload", "-abc", os.path.join(BlastDB_path, "blastp.abc"), "--stream-mirror", "--stream-neg-log10", "-stream-tf",
                        "ceil(200)", "-o", os.path.join(BlastDB_path, "blastp.mci"), "-write-tab", os.path.join(BlastDB_path, "blastp.tab")], check=True)
        # 3.Cluster results
        print("3.Cluster results")
        subprocess.run(["mcl", os.path.join(BlastDB_path, "blastp.mci"), "-te", "8", "-I", "1.4", "-use-tab",
                        os.path.join(BlastDB_path, "blastp.tab"), "-o",
                        os.path.join(BlastDB_path, "protein_cluster_1.4.clstr")], check=True)
        subprocess.run(["mcl", os.path.join(BlastDB_path, "blastp.mci"), "-te", "8", "-I", "2.0", "-use-tab",
                        os.path.join(BlastDB_path, "blastp.tab"), "-o",
                        os.path.join(BlastDB_path, "protein_cluster_2.0.clstr")], check=True)
        subprocess.run(["mcl", os.path.join(BlastDB_path, "blastp.mci"), "-te", "8", "-I", "4.0", "-use-tab",
                        os.path.join(BlastDB_path, "blastp.tab"), "-o",
                        os.path.join(BlastDB_path, "protein_cluster_4.0.clstr")],check=True)
        subprocess.run(["mcl", os.path.join(BlastDB_path, "blastp.mci"), "-te", "8", "-I", "6.0", "-use-tab",
                        os.path.join(BlastDB_path, "blastp.tab"), "-o",
                        os.path.join(BlastDB_path, "protein_cluster_6.0.clstr")],check=True)
        subprocess.run(["mcl", os.path.join(BlastDB_path, "blastp.mci"), "-te", "8", "-I", "8.0", "-use-tab",
                        os.path.join(BlastDB_path, "blastp.tab"), "-o",
                        os.path.join(BlastDB_path, "protein_cluster_8.0.clstr")], check=True)
        subprocess.run(["mcl", os.path.join(BlastDB_path, "blastp.mci"), "-te", "8", "-I", "10.0", "-use-tab",
                        os.path.join(BlastDB_path, "blastp.tab"), "-o",
                        os.path.join(BlastDB_path, "protein_cluster_10.0.clstr")], check=True)
        # Parse MCL cluster files and create CSV files describing the homology relationships between gene nodes
        self.task_mngr.set_task_status(proj_id, task_id, "Write CSV files for clusters")
        cluster_to_csv_parser = ClusterToCSV(proj_id)
        # Load all species-specific protein to gene dicts into cluster2csv parser
        for species in file_dict:
            cluster_to_csv_parser.add_protein_dict("_".join([species[0],species[1]]))
        # Create CSV file for all MCL-generated clusters
        cluster_to_csv_parser.create_csv()
        # Create a new Neo4j graph database from node and relationship CSV files
        print("4. Create DB")
        self.task_mngr.set_task_status(proj_id, task_id, "Building database from CSV files")
        # Use neo4j-admin to create a database from the CSV files
        # The database is created within the projects neo4j folder
        try:
            subprocess.run([
                os.path.join("Projects", str(proj_id), "proj_graph_db", "bin", "neo4j-admin"),
                "import","--id-type","STRING",
                "--nodes:Gene", os.path.join("Projects", str(proj_id),"CSV", "gene_nodes.csv"),
                "--relationships:5_NB", os.path.join("Projects",str(proj_id), "CSV", "gene_5nb.csv"),
                "--relationships:3_NB", os.path.join("Projects", str(proj_id), "CSV", "gene_3nb.csv"),
                "--nodes:Protein", os.path.join("Projects", str(proj_id), "CSV", "protein_nodes.csv"),
                "--relationships:CODING", os.path.join("Projects", str(proj_id), "CSV", "gene_protein_coding.csv"),
                "--relationships:HOMOLOG", os.path.join("Projects", str(proj_id),"CSV", "protein_hmlg.csv")],
                check=True, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
            # Change project status to DB_BUILD
            self.main_db_conn.run(
                "MATCH (proj:Project) WHERE ID(proj) = {proj_id} SET proj.status = {new_status}"
                , {"proj_id": int(proj_id), "new_status": "DB_BUILD"})
        except subprocess.CalledProcessError as err:
            # Change project status to DB_BUILD_FAILED in case build failed
            self.main_db_conn.run(
                "MATCH (proj:Project) WHERE ID(proj) = {proj_id} SET proj.status = {new_status}"
                , {"proj_id": int(proj_id), "new_status": "DB_BUILD_FAILED"})
            print(err.stdout)
            print(err.stderr)
        # Set the Neo4j admin password for the project database
        chars = "ABCDEFGHJKLMNPQRSTUVWXYZ23456789"
        neo4j_pw = ''.join(choice(chars) for _ in range(50))
        # Write password to project folder
        with open(os.path.join("Projects", str(proj_id), "access"), "w") as file:
            file.write(neo4j_pw)
        # Start up the database
        try:
            subprocess.run([os.path.join("Projects", str(proj_id), "proj_graph_db", "bin", "neo4j-admin"),
                    "set-initial-password", neo4j_pw], check=True, stdout=subprocess.PIPE, stderr =subprocess.PIPE)
            subprocess.run([os.path.join("Projects", str(proj_id), "proj_graph_db", "bin", "neo4j"),
                            "start"], check=True, stdout=subprocess.PIPE,stderr=subprocess.PIPE)
            # Wait for database to startup
            while True:
                print("Check db status")
                time.sleep(60)
                status = subprocess.run([os.path.join("Projects", str(proj_id), "proj_graph_db", "bin", "neo4j"),
                                "status"], check=True, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
                if status.returncode == 0:
                    break
                else:
                    print("Wait for DB")
            # Change project status to DB_RUNNING
            self.main_db_conn.run(
                "MATCH (proj:Project) WHERE ID(proj) = {proj_id} SET proj.status = {new_status}"
                , {"proj_id": int(proj_id), "new_status": "DB_RUNNING"})
        except subprocess.CalledProcessError as err:
            # Change project status to DB_START FAILED in case the build database could not be started up
            self.main_db_conn.run(
                "MATCH (proj:Project) WHERE ID(proj) = {proj_id} SET proj.status = {new_status}"
                , {"proj_id": int(proj_id), "new_status": "DB_START_FAILED"})
            print(err.stdout)
            print(err.stderr)
        # Build indices on node properties
        # First connect to the newly build database
        # Retrieve the bolt port number
        print("5. Create Indices")
        self.task_mngr.set_task_status(proj_id, task_id, "Start building indices on project db")
        bolt_port = self.main_db_conn.run("MATCH(proj:Project) WHERE ID(proj)={proj_id} "
                                          "RETURN proj.bolt_port", {"proj_id": int(proj_id)}).single()[0]
        print(bolt_port)
        # Connect to the project DB
        project_db_driver = GraphDatabase.driver("bolt://localhost:" + str(bolt_port),
                                                 auth=basic_auth("neo4j", neo4j_pw), encrypted=False)
        project_db_conn = project_db_driver.session()
        # Build indices
        project_db_conn.run("CREATE INDEX ON :Gene(gene_name)")
        project_db_conn.run("CREATE INDEX ON :Gene(geneId)")
        project_db_conn.run("CREATE INDEX ON :Gene(species)")
        project_db_conn.run("CREATE INDEX ON :Gene(chromosome)")
        project_db_conn.run("CREATE INDEX ON :Gene(gene_descr)")
        project_db_conn.run("CREATE INDEX ON :Protein(protein_descr)")
        project_db_conn.run("CREATE INDEX ON :Protein(proteinId)")
        project_db_conn.run("CREATE INDEX ON :Protein(protein_name)")
        project_db_conn.close()
        self.task_mngr.set_task_status(proj_id, task_id, "Finished")
        print("Finished")



    # For one GFF3 file (or all GFF3 files) in a project, set the annotation mapper and the feature hierarchy
    # Function initializes an instance of the GFF3-parser to check the validity of the annotation mapper string
    # and the feature hierarchy string and then uses the GFF3-parser to parse the beginning of one GFF3 file.
    # Result returned by this function is the first gene and protein node retrieved by the parsing test.
    def set_gff3_parser(self, proj_id, parent_feat, sub_features, name_attr, desc_attr, file_name):
        # Restore function parameters by replacing "\t" back to "_"
        proj_id = proj_id.replace("\t", "_")
        parent_feat = parent_feat.replace("\t", "_")
        sub_features = sub_features.replace("\t", "_")
        name_attr = name_attr.replace("\t", "_")
        desc_attr = desc_attr.replace("\t", "_")
        file_name = file_name.replace("\t", "_")
        # Create a new task and return task-id to user
        task_id = self.task_mngr.define_task(proj_id, "Configure GFF3 parser")
        # # Send task-id to user
        # self.send_data(task_id)
        # # If no file names were specified, set the GFF3 parser for all GFF3 files that are not hidden
        # if not file_names:
        #     file_list = list(
        #         self.main_db_conn.run("MATCH(proj:Project)-[:has_files]->(:File_Manager)-[:file]->(file:File) "
        #                               "WHERE ID(proj)={proj_id} AND file.filetype = 'gff3' AND file.hidden = 'False' "
        #                               "RETURN (file.filename)",
        #                               {"proj_id": int(proj_id)}))
        # If file names were specified, check if they point to existing, non-hidden gff3 files
        # else:
        #     file_list = self.main_db_conn.run(
        #         "MATCH(proj:Project)-[:has_files]->(fileMngr:File_Manager) WHERE ID(proj)={proj_id} "
        #         "MATCH (fileMngr)-[:file]->(file:File) WHERE file.filename IN {file_list} "
        #         "AND file.filetype = 'gff3' AND file.hidden = 'False' RETURN file.filename",
        #         {"proj_id": int(proj_id), "file_list": file_names})
        # Add relative path to file_list
        file_path = os.path.join("Projects", str(proj_id), "Files", file_name)
        # First, test whether the annotation_mapping and the feature_hierarchy fulfill or formal requirements
        # For this, initialize an GFF3_parser instance
        # Do not provide a file path here as this class is used only to verify
        # the correctness of the annotation mapping and the feature hierarchy
        # gff3_parser = Parser.GFF3_parser_gffutils.GFF3Parser("",0,0)
        # valid_annotation_mapper = gff3_parser.set_annotation_mapper(annotation_mapping)
        # valid_feature_hierarchy = gff3_parser.set_feature_hierarchy(feature_hierarchy)
        # # If the gff3 parser rejects one or both of the strings, put these failure into the task status / results
        # # and exit here
        # if not (valid_annotation_mapper and valid_feature_hierarchy):
        #     self.task_mngr.set_task_status(proj_id, task_id, "failed: invalid syntax")
        #     self.task_mngr.add_task_results(proj_id, task_id, "Correct Annotation: " + str(valid_annotation_mapper)+ " Correct Hierarchy: "+ str(valid_feature_hierarchy))
        #     return
        # Else add them to the main-db
        self.main_db_conn.run("MATCH(proj:Project)-[:has_files]->(fileMngr:File_Manager) WHERE ID(proj)={proj_id} "
                    "MATCH (fileMngr)-[:file]->(file:File) WHERE file.filename = {file_name} "
                    "AND file.filetype = 'annotation' AND file.hidden = 'False' "
                    "SET file.parent_feat = {parent_feat} SET file.sub_features = {sub_features} "
                              "SET file.name_attr = {name_attr} SET file.desc_attr = {desc_attr} ",
                    {"proj_id": int(proj_id), "file_name": file_name,
                     "parent_feat": parent_feat, "sub_features": sub_features,
                     "name_attr":name_attr, "desc_attr":desc_attr})
        self.task_mngr.set_task_status(proj_id, task_id, "Added annotation to main-db")
        # Test the parsing of the GFF3 file
        # Make a copy of the gff3-fie, copying only the first 100 lines
        # Try to retrieve one gene node
        # Return this node
        with open(file_path, "r") as original_gff3_file:
            head = list(islice(original_gff3_file, 100))
        with open(os.path.join("Projects", str(proj_id), "Files", "tmp.gff3"), "w") as head_gff3_file:
            for line in head:
                head_gff3_file.write(line)
        # Retrieve the name of the corresponding genome sequence
        genome_file = file_path[:file_path.rfind(".")]+".faa"
        print(genome_file)
        gff3_parser_v2 = GFF3Parser_v2(os.path.join("Projects", str(proj_id), "Files", "tmp.gff3"), genome_file, True, 0,
                                       parent_feat, sub_features, name_attr, desc_attr)
        gene_list = gff3_parser_v2.parse_gff3_file()
        # Return (at max) the first three genes in the gene list
        return_gene_list = []
        for gene in gene_list:
            # Stop when there are already three genes in the return list
            if len(return_gene_list) >= 3:
                break
            # Convert every item of a gene into string format
            gene = [str(item) for item in gene]
            # Get transcript
            #gene_transcript = gff3_parser_v2.get_nt_sequence(gene[0])
            # Get translation
            #gene_translation = gff3_parser_v2.get_prot_sequence(gene[0])
            # Add both to gene list
            #gene.append(gene_transcript)
            #gene.append(gene_translation)
            # Add gene to return list, removing the node ID and species name
            return_gene_list.append(gene[2:])
        # Delete temporary GFF3 file
        os.remove(os.path.join("Projects", str(proj_id), "Files", "tmp.gff3"))
        # Delete temporary transcript and translation files
        gff3_parser_v2.delete_transcripts_translations()
        self.task_mngr.set_task_status(proj_id, task_id, "Finished")
        self.send_data("\n".join(["\t".join(item) for item in return_gene_list]))
        if len(gene_list) >= 1:
            self.task_mngr.add_task_results(proj_id, task_id, "Success")
        else:
            self.task_mngr.add_task_results(proj_id, task_id, "Failed")




