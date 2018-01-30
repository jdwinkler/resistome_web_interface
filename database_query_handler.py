import psycopg2
import psycopg2.extras
import os
from collections import defaultdict


class ResistomeDBHandler:

    def __init__(self):

        try:
            user_name = os.environ['user_name'.upper()]
            password = os.environ['sql_password'.upper()]
        except:
            user_name = 'james'
            password = 'winkler'

        connection = psycopg2.connect(
            "dbname='resistome' user='%s' host='localhost' password='%s'" % (user_name, password))
        connection.set_session(readonly=True)
        cursor = connection.cursor(cursor_factory=psycopg2.extras.RealDictCursor)

        self.connection = connection
        self.cursor = cursor

    def obtain_mutants_mutations_annotations(self, std_gene_names):

        '''
        :param std_gene_names: 
        :return: 
        '''

        query_genes = [x[1] for x in std_gene_names]

        # get mutant id, name pairs for provided gene list
        self.cursor.execute('select resistome.mutations.mutant_id,'
                            'resistome.mutations.name '
                            'from resistome.mutations where resistome_mutations.name = ANY(%s) ', (query_genes,))

        id_name_pairs = defaultdict(list)
        for record in self.cursor:
            id_name_pairs[record['name']].append(record['mutant_id'])

        self.cursor.execute('select resistome.paper.title, '
                            'resistome.paper.doi,'
                            'resistome.phenotypes.phenotype,'
                            'resistome.phenotypes.phenotype_class, '
                            'resistome.phenotypes.ontology_root, '
                            'resistome.phenotypes.phenotype_type, '
                            'resistome.mutations.name,'
                            'resistome.annotation.mutation_type,'
                            'resistome.annotation.annotation,'
                            'resistome.mutants.mutant_id,'
                            'where resistome.mutants.mutant_id = ANY(%s) '
                            'inner join resistome.papers on (resistome.papers.paper_id = resistome.mutants.mutant_id)'
                            'inner join resistome.mutations on (resistome.mutations.mutant_id = resistome.mutants.mutant_id)'
                            'inner join resistome.annotation on (resistome.annotation.gene_id = resistome.mutations.gene_id)'
                            'inner join resistome.phenotypes on (resistome.phenotypes.mutant_id = resistome.mutants.mutant_id)')



        return None

    def standardize_input_gene_names(self, gene_names):

        """
        
        :param gene_names: list of strings represent E. coli genes
        :return: converted_gene_names that are standardized into bnumbers, standard names
        
        """

        self.cursor.execute('select name, strain, species_accession, mg1655_accession '
                            'from resistome.gene_standardization where name = ANY(%s)', gene_names)

        requested_genes = set(gene_names)
        found_genes = set()

        output = []

        for record in self.cursor:
            output.append((record['name'], record['mg1655_accession']))
            found_genes.add(record['name'])

        remaining_genes = requested_genes - found_genes
        for gene in remaining_genes:
            output.append((gene, gene))

        return output

