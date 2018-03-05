# todo: add tornado skeleton

import tornado.ioloop
import tornado.web
import tornado.httpserver
import os.path
from database_query_handler import ResistomeDBHandler
import uuid
import datetime

MAIN_DIR = ''
QUERY_DIR = ''
TEMP_DIR = ''
resistome_handler = None
phenotype_listing = None

# display order for web tables

display_order = ['doi',
                 'id',
                 'phenotypes',
                 'phenotype_types',
                 'annotations']

# column order for dicts representing output genotypes in downloadable tables

serialize_order = ['title',
                   'doi',
                   'id',
                   'affected_genes',
                   'affected_phenotypes',
                   'species',
                   'strain',
                   'phenotypes',
                   'phenotype_types',
                   'root',
                   'annotations',
                   'expressions']


class MainHandler(tornado.web.RequestHandler):
    def get(self):

        """
        
        Main page rendering.
        
        :return: 
        """

        if 'summary' in self.request.uri:
            self.render(os.path.join(MAIN_DIR, 'summary.html'))
        else:
            self.render(os.path.join(MAIN_DIR, 'main.html'))

    def post(self):

        """
        
        Renders a search page based the value provided in the 'search' parameter when post request is received.
        If the value is missing or does not fall under the categories below, the main page is rendered.
        
        :return: 
        """

        query_type = self.get_argument('search', default=None)

        if 'Gene' in query_type:

            self.render(os.path.join(QUERY_DIR, 'gene_search.html'))

        elif 'Phenotype' in query_type:

            self.render(os.path.join(QUERY_DIR, 'phenotype_search.html'), phenotype_list = phenotype_listing)

        elif 'Genotype' in query_type:

            self.render(os.path.join(QUERY_DIR, 'rs_search.html'))

        else:

            self.render(os.path.join(MAIN_DIR, 'main.html'))


class FileHandler(tornado.web.RequestHandler):

    def get(self, fname):

        """
        
        Serves up temporary files representing data requested during searches.
        
        :param fname: 
        :return: 
        """

        file_name = os.path.join(TEMP_DIR, fname)

        if file_name is None:
            self.finish()

        self.set_header('Content-Type', 'application/text')
        self.set_header('Content-Disposition', 'attachment; filename=resistome_data.txt')
        with open(file_name, 'rU') as f:
            for line in f:
                self.write(line)
        self.finish()


class QueryHandler(tornado.web.RequestHandler):

    def mutant_dict_to_web_tuple(self, mutant_dict):

        output = []

        for x in display_order:

            if isinstance(mutant_dict[x], list):
                if len(mutant_dict[x]) > 0:
                    output.append(mutant_dict[x])
                else:
                    output.append('N/A')
            else:
                output.append(str(mutant_dict[x]))

        return tuple(output)

    def mutant_dict_to_serialized_tuples(self, mutant_dict):

        output = []

        for x in serialize_order:

            if isinstance(mutant_dict[x], list):
                output.append(','.join(mutant_dict[x]))
            else:
                output.append(str(mutant_dict[x]))

        return tuple(output)

    def output_to_file(self, queried_gene_names, queried_phenotype_names, mutant_text_data):

        now = datetime.datetime.now()

        file_name = str(uuid.uuid4())

        fhandle = open(os.path.join(TEMP_DIR, file_name), 'w+')

        fhandle.write('# Queried genes: %s' % ','.join(queried_gene_names) + '\n')
        fhandle.write('# Queried phenotypes: %s' % ','.join(queried_phenotype_names) + '\n')
        fhandle.write('# Date of search: %s' % now.strftime("%Y-%m-%d %H:%M") + '\n')

        fhandle.write('\t'.join(serialize_order) + '\n')

        for x in mutant_text_data:

            fhandle.write('\t'.join(x) + '\n')

        fhandle.close()

        return file_name

    def post(self):

        recommendation_request = self.get_argument('recommendation', False)
        recommendation_request = True if recommendation_request == 'true' else False

        if recommendation_request:

            feature_type = self.get_argument('ftype', 'gene')

            multi_gene_query = self.get_arguments('gene_list')
            multi_gene_query = multi_gene_query[0].split('\n')
            multi_gene_query = [x.strip() for x in multi_gene_query]

            g_names, query_names, scores_dict, mutant_text_array = resistome_handler.find_similar_genotypes(multi_gene_query,
                                                                                                        feature_type)

            serialized_output = [self.mutant_dict_to_serialized_tuples(x) for x in mutant_text_array]

            fname = self.output_to_file(g_names, query_names, serialized_output)

            rs_display_order = ['doi',
                             'id',
                             'score',
                             'phenotypes',
                             'phenotype_types',
                             'annotations']

            self.render(os.path.join(QUERY_DIR, 'rs_results.html'),
                        records=mutant_text_array,
                        converted_gene_names=g_names,
                        converted_feature_names=query_names,
                        scores=scores_dict,
                        temp_file_name=fname,
                        field_order=rs_display_order)

        else:

            specific_flag = self.get_argument('specific', True)
            specific_flag = True if specific_flag == 'true' else False
            ge_flag = True if self.get_argument('include_ge', None) is not None else False

            multi_gene_query = self.get_arguments('gene_list')
            multi_phenotype_query = self.get_arguments('phenotype_list')

            if len(multi_gene_query) == 1:
                multi_gene_query = multi_gene_query[0].split('\n')
                multi_gene_query = [x.strip() for x in multi_gene_query]

            if len(multi_phenotype_query) == 1:
                multi_phenotype_query = multi_phenotype_query[0].split('\n')
                multi_phenotype_query = [x.strip() for x in multi_phenotype_query]

            g_names, p_names, mutant_text_array = resistome_handler.get_mutant_output(multi_gene_query,
                                                                           multi_phenotype_query,
                                                                           specific_flag=specific_flag,
                                                                                      ge_flag = ge_flag)

            serialized_output = [self.mutant_dict_to_serialized_tuples(x) for x in mutant_text_array]

            fname = self.output_to_file(g_names, p_names, serialized_output)

            self.render(os.path.join(QUERY_DIR, 'results.html'),
                        records=mutant_text_array,
                        converted_gene_names=g_names,
                        converted_phenotype_names=p_names,
                        temp_file_name=fname,
                        field_order=display_order)


def main():

    global MAIN_DIR, QUERY_DIR, TEMP_DIR, phenotype_listing, resistome_handler

    current_location = os.path.realpath(__file__)
    path = os.path.split(current_location)[0]

    STATIC_DIR = os.path.join(path, 'static')
    TEMPLATE_DIR = os.path.join(path, 'templates')
    MAIN_DIR = os.path.join(path, 'main')
    TEMP_DIR = os.path.join(path, 'tmp')
    QUERY_DIR = os.path.join(path, 'search')

    resistome_handler = ResistomeDBHandler(os.path.join(STATIC_DIR, 'Serialized Vector Sets.obj'))
    phenotype_listing = resistome_handler.get_phenotype_listing()
    temp = []

    for phenotype in phenotype_listing:

        if phenotype.islower():
            temp.append((phenotype, phenotype + ' (category)'))
        else:
            temp.append((phenotype, phenotype))

    phenotype_listing = temp

    handlers = [
        (r"/", MainHandler),
        (r"/summary", MainHandler),
        (r"/query/", QueryHandler),
        (r"/temp/(.*)", FileHandler),
        (r'/static/(.*)', tornado.web.StaticFileHandler, {'path': STATIC_DIR})]

    application = tornado.web.Application(handlers, debug=True, template_path=TEMPLATE_DIR)

    http_server = tornado.httpserver.HTTPServer(application)

    port = int(os.environ.get("PORT", 5000))
    http_server.listen(port)
    http_server.start()
    tornado.ioloop.IOLoop.instance().start()

if __name__ == "__main__":
    main()

