import psycopg2
import psycopg2.extras
import os
import urlparse
import cPickle
import recommendation_system as rs
from collections import defaultdict


class ResistomeDBHandler:

    def __init__(self, path_to_serialized_vectors):

        try:
            # heroku hosted version
            urlparse.uses_netloc.append("postgres")
            url = urlparse.urlparse(os.environ["DATABASE_URL"])

            connection = psycopg2.connect(
                database=url.path[1:],
                user=url.username,
                password=url.password,
                host=url.hostname,
                port=url.port
            )

        except:
            # local version of the database
            user_name = 'james'
            password = 'winkler'

            connection = psycopg2.connect("dbname='resistome' user='%s' host='localhost' password='%s'"
                                          % (user_name, password))

        # remove possibility of a query editing data
        connection.set_session(readonly=True)
        cursor = connection.cursor(cursor_factory=psycopg2.extras.RealDictCursor)

        self.connection = connection
        self.cursor = cursor

        # load data used for recommendation system
        self.vector_db = cPickle.load(open(path_to_serialized_vectors, 'rU'))

    def __del__(self):

        self.connection.close()

    def get_phenotype_listing(self):

        self.cursor.execute('SELECT resistome.phenotype_standardization.standard_name,'
                            'resistome.phenotype_standardization.phenotype_type '
                            'from resistome.phenotype_standardization')

        output = set()

        for record in self.cursor:

            output.add(record['standard_name'])
            output.add(record['phenotype_type'])

        return sorted(list(output))

    def find_similar_genotypes(self, gene_names, feature_types):

        if feature_types != 'gene' and feature_types != 'go':
            raise ValueError('Unknown feature type: %s' % feature_types)

        std_gene_names, gene_display_names = self.standardize_input_gene_names(gene_names)

        std_genes = [x[1] for x in std_gene_names]

        query_vector, query_features = rs.build_proposed_vector(self.cursor, std_genes, feature_types)

        pw_distances = sorted(rs.pairwise_distance_vector(query_vector,
                                                          self.vector_db[feature_types],
                                                          method='jaccard'),
                              key=lambda x: x[1])

        min_mutants = 25
        mutant_ids = []

        scores_dict = dict()

        for (m_id, score) in pw_distances:

            # pull mutant ids if minimum count hasn't been met yet
            if m_id <= self.vector_db[feature_types][1] or len(mutant_ids) < min_mutants:
                mutant_ids.append(m_id)
                scores_dict[m_id] = '%0.3f' % (1.0 - score)

        gene_text_output = self.prep_gene_mutant_output(defaultdict(list),
                                                        self.query_mutant_genotypes(mutant_ids, ge_flag=False),
                                                        only_affected_genes=False,
                                                        display_converter=gene_display_names)

        sorted_text_output = []

        for m_dict in gene_text_output:
            sorted_text_output.append((scores_dict[m_dict['id']], m_dict))

        sorted_text_output = sorted(sorted_text_output, key=lambda x: x[0], reverse=True)
        sorted_text_output = [x[1] for x in sorted_text_output]

        gene_names = [x[0] + ' (%s)' % x[1] for x in std_gene_names]

        return gene_names, query_features, scores_dict, sorted_text_output

    def get_mutant_output(self, gene_names, phenotype_names, specific_flag=False, ge_flag=False):

        """
        
        Extract mutant output related to gene names, phenotype names.
        
        :param gene_names: 
        :param phenotype_names
        :param specific_flag: 
        :param ge_flag
        :return: 
        """

        if not isinstance(gene_names, list):
            gene_names = [gene_names]

        std_gene_names, gene_display_names = self.standardize_input_gene_names(gene_names)
        std_phenotype_names, _ = self.standardize_input_phenotype_names(phenotype_names)

        query_genes = [x[1] for x in std_gene_names]
        query_phenotypes = [x[1] for x in std_phenotype_names]

        # get mutant id, name pairs for provided gene list
        self.cursor.execute('select resistome.mutations.mutant_id,'
                            'resistome.mutations.name '
                            'from resistome.mutations where resistome.mutations.name = ANY(%s) ', (query_genes,))

        gene_results = self.cursor.fetchall()
        gene_mutant_ids = []
        gene_mutant_dict = defaultdict(set)
        for result in gene_results:
            gene_mutant_ids.append(result['mutant_id'])
            gene_mutant_dict[result['mutant_id']].add(result['name'])

        self.cursor.execute('select resistome.mutants.mutant_id,'
                            'resistome.phenotypes.phenotype '
                            'from resistome.mutants '
                            'inner join resistome.phenotypes on (resistome.mutants.mutant_id = resistome.phenotypes.mutant_id)'
                            'where resistome.phenotypes.phenotype = ANY(%s) ', (query_phenotypes,))

        pheno_results = self.cursor.fetchall()
        pheno_mutant_ids = []
        pheno_mutant_dict = defaultdict(set)
        for result in pheno_results:
            pheno_mutant_ids.append(result['mutant_id'])
            pheno_mutant_dict[result['mutant_id']].add(result['phenotype'])

        gene_names = [x[0] + ' (%s)' % x[1] for x in std_gene_names]
        pheno_names = [x[0] + ' (%s)' % x[1] for x in std_phenotype_names]

        gene_text_output = self.prep_gene_mutant_output(gene_mutant_dict,
                                                        self.query_mutant_genotypes(gene_mutant_ids, ge_flag=ge_flag),
                                                        only_affected_genes=specific_flag,
                                                        display_converter=gene_display_names)

        pheno_text_output = self.prep_pheno_mutant_output(pheno_mutant_dict,
                                                        self.query_mutant_genotypes(pheno_mutant_ids, ge_flag=ge_flag),
                                                        only_affected_phenotypes=specific_flag,
                                                        display_converter=None)

        gene_text_output.extend(pheno_text_output)

        return gene_names, pheno_names, gene_text_output

    def query_mutant_genotypes(self, mutant_ids, ge_flag=False):

        if ge_flag:
            self.cursor.execute('select array_agg(distinct resistome.papers.title) as title, '
                                'array_agg(distinct resistome.papers.doi) as doi,'
                                'array_agg(resistome.phenotypes.phenotype) as phenotype,'
                                'array_agg(distinct resistome.phenotypes.phenotype_class) as phenotype_class, '
                                'array_agg(distinct resistome.phenotypes.ontology_root) as ontology_root, '
                                'array_agg(resistome.phenotypes.phenotype_type) as phenotype_type, '
                                'array_agg(resistome.mutations.name) as genes,'
                                'array_agg(resistome.annotations.mutation_type) as mutation_type,'
                                'array_agg(resistome.annotations.annotation) as annotation,'
                                'array_agg(resistome.annotations.gene_id) as gene_ids, '
                                'array_agg(resistome.expressions.name) as de_genes, '
                                'array_agg(resistome.expressions.fold_change) as fold_changes, '
                                'array_agg(resistome.expressions.status) as status,'
                                'resistome.mutants.mutant_id, '
                                'resistome.mutants.species, '
                                'resistome.mutants.strain '
                                'from resistome.mutants '
                                'inner join resistome.papers on (resistome.papers.paper_id = resistome.mutants.paper_id) '
                                'inner join resistome.mutations on (resistome.mutations.mutant_id = resistome.mutants.mutant_id) '
                                'inner join resistome.annotations on (resistome.annotations.gene_id = resistome.mutations.gene_id) '
                                'inner join resistome.phenotypes on (resistome.phenotypes.mutant_id = resistome.mutants.mutant_id) '
                                'left join resistome.expressions on (resistome.expressions.mutant_id = resistome.mutants.mutant_id) '
                                'where resistome.mutants.mutant_id = ANY(%s) '
                                'group by resistome.mutants.mutant_id ',
                                (mutant_ids,))
        else:
            self.cursor.execute('select array_agg(distinct resistome.papers.title) as title, '
                                'array_agg(distinct resistome.papers.doi) as doi,'
                                'array_agg(resistome.phenotypes.phenotype) as phenotype,'
                                'array_agg(distinct resistome.phenotypes.phenotype_class) as phenotype_class, '
                                'array_agg(distinct resistome.phenotypes.ontology_root) as ontology_root, '
                                'array_agg(resistome.phenotypes.phenotype_type) as phenotype_type, '
                                'array_agg(resistome.mutations.name) as genes,'
                                'array_agg(resistome.annotations.mutation_type) as mutation_type,'
                                'array_agg(resistome.annotations.annotation) as annotation,'
                                'array_agg(resistome.annotations.gene_id) as gene_ids, '
                                'resistome.mutants.mutant_id, '
                                'resistome.mutants.species, '
                                'resistome.mutants.strain '
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

    def standardize_input_phenotype_names(self, phenotype_names):

        """

        :param phenotype_names: list of strings represent E. coli genes
        :return: converted_gene_names that are standardized into bnumbers, standard names

        """

        phenotype_names = [x.upper() for x in phenotype_names]

        self.cursor.execute('select name, standard_name from resistome.phenotype_standardization '
                            'where resistome.phenotype_standardization.name = ANY(%s) OR '
                            'resistome.phenotype_standardization.standard_name = ANY(%s) OR '
                            'resistome.phenotype_standardization.phenotype_type = ANY(%s)', (phenotype_names,
                                                                                             phenotype_names,
                                                                                             [x.lower() for x in phenotype_names]))

        requested_phenotypes = set(phenotype_names)
        found_phenotypes = set()

        output = []

        display_name = dict()

        for record in self.cursor:
            output.append((record['name'], record['standard_name']))
            found_phenotypes.add(record['name'])
            found_phenotypes.add(record['standard_name'])

            display_name[record['name']] = record['name']
            display_name[record['standard_name']] = record['standard_name']

        remaining_phenotypes = requested_phenotypes - found_phenotypes
        for pheno in remaining_phenotypes:
            output.append((pheno, pheno))

        return output, display_name

    @staticmethod
    def prep_pheno_mutant_output(mutant_to_queried_phenotypes, records, only_affected_phenotypes=False,
                                 display_converter=None):

        output_text_lines = []

        if display_converter is None:
            display_converter = dict()

        for record in records:

            title = record['title'][0]
            mutant_id = record['mutant_id']
            species = record['species']
            strain = record['strain']
            doi = record['doi'][0]
            phenotype = record['phenotype']
            phenotype_type = record['phenotype_type']

            phenotype_type = ['Sensitive' if x == 'S' else 'Resistant' for x in phenotype_type]

            phenotypes_to_highlight = mutant_to_queried_phenotypes[record['mutant_id']]

            filter_phenotypes = set()

            for p, t in zip(phenotype, phenotype_type):

                if only_affected_phenotypes and p not in phenotypes_to_highlight:
                    pass
                else:
                    filter_phenotypes.add((p, t))

            filter_phenotypes = list(filter_phenotypes)
            phenotype = [x[0] for x in filter_phenotypes]
            phenotype_type = [x[1] for x in filter_phenotypes]

            root = record['phenotype_class']
            mutated_genes = record['genes']
            mutation_types = record['mutation_type']
            annotations = record['annotation']
            gene_ids = record['gene_ids']

            affected_phenotypes = sorted(list(set(phenotypes_to_highlight) & set(phenotype)))

            gene_annotation_output = set()

            for gene, m_type, annotation in zip(mutated_genes, mutation_types, annotations):

                if 'large_' in m_type:
                    gene_annotation_output.add(ResistomeDBHandler.standard_mutation_formatting(m_type, annotation))
                else:
                    gene_annotation_output.add(gene + ' ' + ResistomeDBHandler.standard_mutation_formatting(m_type,
                                                                                                            annotation))

            expression_gene_names = record.get('de_genes', [])
            expression_fold_change = record.get('fold_changes', [])
            directionality = record.get('status', [])

            expression_annotation_output = set()

            for gene, fold_change, status in zip(expression_gene_names, expression_fold_change, directionality):

                if fold_change is not None:
                    expression_annotation_output.add(gene + ' (%f)' % str(fold_change))
                elif status is not None:
                    expression_annotation_output.add(gene + ' (%s)' % status)

            gene_annotation_output = sorted(list(gene_annotation_output))
            expression_annotation_output = sorted(list(expression_annotation_output))

            output_dict = {'title': title,
                           'doi': doi,
                           'id': mutant_id,
                           'species': species,
                           'strain': strain,
                           'phenotypes': phenotype,
                           'phenotype_types': phenotype_type,
                           'root': root,
                           'annotations': gene_annotation_output,
                           'expressions': expression_annotation_output,
                           'affected_genes' : ['N/A'],
                           'affected_phenotypes': affected_phenotypes}

            output_text_lines.append(output_dict)

        return output_text_lines

    @staticmethod
    def prep_gene_mutant_output(mutant_to_queried_genes, records, only_affected_genes=False, display_converter=None):

        output_mutant_dicts = []

        if display_converter is None:
            display_converter = dict()

        for record in records:

            title = record['title'][0]
            mutant_id = record['mutant_id']
            species = record['species']
            strain = record['strain']
            doi = record['doi'][0]
            phenotype = record['phenotype']
            phenotype_type = record['phenotype_type']

            phenotype_type = ['Sensitive' if x == 'S' else 'Resistant' for x in phenotype_type]

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
                    gene_annotation_output.add(display_converter.get(gene, gene) + ' ' +
                                               ResistomeDBHandler.standard_mutation_formatting(m_type, annotation))

            expression_gene_names = record.get('de_genes', [])
            expression_fold_change = record.get('fold_changes', [])
            directionality = record.get('status', [])

            expression_annotation_output = set()

            for gene, fold_change, status in zip(expression_gene_names, expression_fold_change, directionality):

                if fold_change is not None:
                    expression_annotation_output.add(gene + ' (%f)' % fold_change)
                elif status is not None:
                    expression_annotation_output.add(gene + ' (%s)' % status)

            gene_annotation_output = sorted(list(gene_annotation_output))
            expression_annotation_output = sorted(list(expression_annotation_output))

            output_dict = {'title': title,
                           'doi': doi,
                           'id': mutant_id,
                           'species': species,
                           'strain': strain,
                           'phenotypes': phenotype,
                           'phenotype_types': phenotype_type,
                           'root': root,
                           'annotations': gene_annotation_output,
                           'expressions': expression_annotation_output,
                           'affected_genes' : affected_genes,
                           'affected_phenotypes': ['N/A']}

            output_mutant_dicts.append(output_dict)

        return output_mutant_dicts

    @staticmethod
    def standard_mutation_formatting(mutation_type, annotation):

        """
        
        Generates a 'visually pleasing' version of internal Resistome mutation data depending on the mutation type
        and annotation data available.
        
        :param mutation_type: 
        :param annotation: 
        :return: 
        """

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
            prefix = '+' if annotation['duplication'][1] > 0 else '-'
            return annotation['duplication'][0] + '|' + prefix + str(annotation['duplication'][1]) + ' bp|' + \
                   annotation['duplication'][2]
        elif mutation_type == 'large_amplification':
            return annotation[mutation_type][0] + '-' + annotation[mutation_type][1] + ' (amplified %iX)' % \
                                                                                       annotation[mutation_type][2]
        elif mutation_type == 'large_deletion':
            return annotation[mutation_type][0] + '-' + annotation[mutation_type][1] + ' (deleted)'
        elif mutation_type == 'large_inversion':
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
        elif mutation_type == 'frameshift':
            return 'frameshift'
        elif mutation_type == 'rbs_tuned':
            return '(engineered RBS)'
        else:
            return '(%s)' % mutation_type
