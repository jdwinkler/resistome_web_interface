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

    def get_mutant_output(self, gene_names, specific_flag=False):

        """
        
        :param gene_names: 
        :param specific_flag: 
        :return: 
        """

        if not isinstance(gene_names, list):
            gene_names = [gene_names]

        std_gene_names, display_names = self.standardize_input_gene_names(gene_names)

        query_genes = [x[1] for x in std_gene_names]

        # get mutant id, name pairs for provided gene list
        self.cursor.execute('select resistome.mutations.mutant_id,'
                            'resistome.mutations.name '
                            'from resistome.mutations where resistome.mutations.name = ANY(%s) ', (query_genes,))

        results = self.cursor.fetchall()
        mutant_ids = []
        mutant_dict = defaultdict(set)
        for result in results:
            mutant_ids.append(result['mutant_id'])
            mutant_dict[result['mutant_id']].add(result['name'])

        gene_names = [x[0] + ' (%s)' % x[1] for x in std_gene_names]

        return gene_names, self.prep_mutant_output(mutant_dict,
                                       self.query_mutant_genotypes(mutant_ids),
                                       only_affected_genes=specific_flag,
                                        display_converter=display_names)

    def query_mutant_genotypes(self, mutant_ids):

        self.cursor.execute('select array_agg(distinct resistome.papers.title) as title, '
                       'array_agg(distinct resistome.papers.doi) as doi,'
                        'array_agg(resistome.phenotypes.phenotype) as phenotype,'
                       'array_agg(distinct resistome.phenotypes.phenotype_class) as phenotype_class, '
                       'array_agg(distinct resistome.phenotypes.ontology_root) as ontology_root, '
                       'array_agg(resistome.phenotypes.phenotype_type) as phenotype_type, '
                       'array_agg(resistome.mutations.name) as genes,'
                       'array_agg(resistome.annotations.mutation_type) as mutation_type,'
                       'array_agg(resistome.annotations.annotation) as annotation,'
                       'array_agg(resistome.annotations.gene_id) as gene_ids,'
                       'resistome.mutants.mutant_id '
                       'from resistome.mutants '
                       'inner join resistome.papers on (resistome.papers.paper_id = resistome.mutants.paper_id) '
                       'inner join resistome.mutations on (resistome.mutations.mutant_id = resistome.mutants.mutant_id) '
                       'inner join resistome.annotations on (resistome.annotations.gene_id = resistome.mutations.gene_id) '
                       'inner join resistome.phenotypes on (resistome.phenotypes.mutant_id = resistome.mutants.mutant_id) '
                       'where resistome.mutants.mutant_id = ANY(%s) '
                       'group by resistome.mutants.mutant_id ',
                       (mutant_ids,))

        mutants = self.cursor.fetchall()

        if mutants is None or len(mutants) == 0:
            return []
        else:
            return mutants

    def standardize_input_gene_names(self, gene_names):

        """
        
        :param gene_names: list of strings represent E. coli genes
        :return: converted_gene_names that are standardized into bnumbers, standard names
        
        """

        gene_names = [x.upper() for x in gene_names]

        self.cursor.execute('select name, mg1655_accession, species_accession '
                            'from resistome.gene_standardization where name = ANY(%s) '
                            'OR mg1655_accession = ANY(%s) '
                            'OR species_accession = ANY(%s)', (gene_names, gene_names, gene_names))

        requested_genes = set(gene_names)
        found_genes = set()

        output = []

        display_name = dict()

        for record in self.cursor:
            output.append((record['name'], record['species_accession']))
            found_genes.add(record['name'])
            found_genes.add(record['mg1655_accession'])
            found_genes.add(record['species_accession'])

            display_name[record['mg1655_accession']] = record['mg1655_accession']
            display_name[record['species_accession']] = record['species_accession']
            display_name[record['name']] = record['mg1655_accession']

        remaining_genes = requested_genes - found_genes
        for gene in remaining_genes:
            output.append((gene, gene))

        return output, display_name

    @staticmethod
    def prep_mutant_output(mutant_to_queried_genes, records, only_affected_genes=False, display_converter=None):

        output_text_lines = []

        if display_converter is None:
            display_converter = dict()

        for record in records:

            title = record['title'][0]
            doi = record['doi'][0]
            phenotype = record['phenotype']
            phenotype_type = record['phenotype_type']

            filter_phenotypes = set()

            for p, t in zip(phenotype, phenotype_type):
                filter_phenotypes.add((p, t))

            filter_phenotypes = list(filter_phenotypes)
            phenotype = [x[0] for x in filter_phenotypes]
            phenotype_type = [x[1] for x in filter_phenotypes]

            root = record['phenotype_class']

            genes_to_highlight = mutant_to_queried_genes[record['mutant_id']]

            mutated_genes = record['genes']
            mutation_types = record['mutation_type']
            annotations = record['annotation']
            gene_ids = record['gene_ids']

            affected_genes = sorted(list(set(display_converter.get(x, x) for x in set(genes_to_highlight) & set(mutated_genes))))

            gene_annotation_output = set()

            for gene, m_type, annotation, gene_id in zip(mutated_genes, mutation_types, annotations, gene_ids):

                if only_affected_genes and gene not in genes_to_highlight:
                    continue

                if 'large_' in m_type:
                    gene_annotation_output.add(ResistomeDBHandler.standard_mutation_formatting(m_type, annotation))
                else:
                    gene_annotation_output.add(display_converter.get(gene, gene) + ': ' + ResistomeDBHandler.standard_mutation_formatting(m_type, annotation))

            gene_annotation_output = sorted(list(gene_annotation_output))

            output_text_lines.append(
                (affected_genes, doi, phenotype, phenotype_type, root, gene_annotation_output))

        return output_text_lines

    @staticmethod
    def standard_mutation_formatting(mutation_type, annotation):

        '''

        :param mutation_type: 
        :param annotation: 
        :return: 

        '''

        if mutation_type == 'aa_snps':
            str_output = []
            for aa_array in annotation['aa_snps']:
                str_output.append(aa_array[1] + str(aa_array[0]) + aa_array[2])
            return ', '.join(str_output)
        elif mutation_type == 'nuc_snps':
            str_output = []
            for nuc_array in annotation['nuc_snps']:
                str_output.append(nuc_array[1] + str(nuc_array[0]) + nuc_array[2])
            return 'SNP(s) ' + ', '.join(str_output)
        elif mutation_type == 'indel':
            str_output = []
            for indel_array in annotation['indel']:
                prefix = '+' if indel_array[1] > 0 else '-'
                try:
                    str_output.append(
                        'indel ' + str(indel_array[0]) + '|' + prefix + str(abs(indel_array[1])) + ' bp|' +
                        indel_array[2])
                except:
                    str_output.append('indel, but location/size unclear')
            return ','.join(str_output)
        elif mutation_type == 'is_insertion':
            return 'IS insertion ' + annotation['is_insertion'][0]
        elif mutation_type == 'del':
            return '(deleted)'
        elif mutation_type == 'intergenic':
            return 'intergenic ' + annotation['intergenic'][0] + '/' + annotation['intergenic'][1]
        elif mutation_type == 'amplified':
            return '(amplified %iX)' % annotation['amplified']
        elif mutation_type == 'duplication':
            return annotation['duplication'][0] + '|' + str(annotation['duplication'][1]) + 'bp|' + \
                   annotation['duplication'][2]
        elif mutation_type == 'large_amplification':
            return annotation[mutation_type][0] + '-' + annotation[mutation_type][1] + ' (amplified %iX)' % \
                                                                                       annotation[mutation_type][2]
        elif mutation_type == 'large_deletion':
            return annotation[mutation_type][0] + '-' + annotation[mutation_type][1] + ' (deleted)'
        elif mutation_type == 'larger_inversion':
            return annotation[mutation_type][0] + ' <=> ' + annotation[mutation_type][1] + ' (inverted)'
        elif mutation_type == 'mutated':
            return '(mutated)'
        elif mutation_type == 'oe':
            return 'overexpressed'
        elif mutation_type == 'plasmid':
            return 'plasmid expression'
        elif mutation_type == 'truncated':
            return 'truncated'
        elif mutation_type == 'rep':
            return 'repressed'
        elif mutation_type == 'integrated':
            return 'integrated'
        else:
            raise ValueError('Unknown mutation type: %s' % mutation_type)

