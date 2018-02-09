import itertools
import math
import os


class Feature(object):
    # encapsulates feature details

    def __init__(self, feature_name, feature_details):

        self.name = feature_name.upper()
        self.details = set(feature_details)

    def __eq__(self, other):

        if (not isinstance(other, Feature)):
            return False
        else:
            return self.name == other.name and self.details == other.details

    def __hash__(self):

        return hash(self.name)

    def __str__(self):

        return self.name


class Vector:
    # feature name converter is standardized gene name to a list of features of some sort (string to list map)
    # examples: seed ontology, go ontology

    def __init__(self, internal_id, raw_features):

        self.id = internal_id
        self.feature_set = self.build_internal_features(raw_features)

    def build_internal_features(self, feature_set):

        output_set = set()
        output_set.update(feature_set)

        return output_set

    def __len__(self):

        return len(self.feature_set)

    def combine(self, v2):

        return Vector(self.id + v2.id, self.feature_set.union(v2.feature_set))

    def union(self, v2):

        return len(self.feature_set.union(v2.feature_set))

    def intersection(self, v2):

        return len(self.feature_set.intersection(v2.feature_set))

    def difference(self, v2):

        return len(self.feature_set.difference(v2.feature_set))

    def distance(self, v2, method='euclidean'):

        if method == 'euclidean':
            return math.sqrt(self.difference(v2)) / float(len(self))

        if method == 'similarity':
            return float(self.intersection(v2)) / float(len(self))

        if method == 'cosine':
            return 1 - float(self.intersection(v2)) / (math.sqrt(len(self)) * math.sqrt(len(v2)))

        if method == 'ochiai':
            product = set()
            for combo in itertools.product(self.feature_set, v2.feature_set):
                product.add(combo)
            return 1 - float(self.intersection(v2)) / float(len(product))

        if method == 'jaccard':
            return 1.0 - float(self.intersection(v2)) / (float(len(self) + len(v2)) - float(self.intersection(v2)))

        raise AssertionError('Requested method %s not yet implemented' % method)


def build_proposed_vector(cursor, genes, type_of_feature):

    # assumes genes is already standardized

    features = genes

    if type_of_feature == 'go':

        cursor.execute('select go_term as go,'
                       'from resistome.gene_ontology '
                       'where resistome.gene_ontology.accession = ANY(%s)', (genes,))

        features = []
        for record in cursor:
            features.append(record['go'])

    elif type_of_feature == 'metabolite':

        cursor.execute('select metabolite as metabolite,'
                       'from resistome.metabolomics '
                       'where resistome.metabolomics.accession = ANY(%s)', (genes,))

        features = []
        for record in cursor:
            features.append(record['metabolite'])

    return Vector('test', [Feature(x, []) for x in features])


def pairwise_distance_vector(proposed_vector, vector_set, method='euclidean'):

    distances = []

    for vector in vector_set:
        distances.append((vector.id, proposed_vector.distance(vector, method=method)))

    return distances