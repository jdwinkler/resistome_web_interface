# todo: add tornado skeleton

import tornado.ioloop
import tornado.web
import os.path
from email.mime.multipart import MIMEMultipart
from email.mime.text import MIMEText
from collections import defaultdict
from threading import Thread
from database_query_handler import ResistomeDBHandler

MAIN_DIR = ''
QUERY_DIR = ''
resistome_handler = None

class MainHandler(tornado.web.RequestHandler):
    def get(self):

        self.render(os.path.join(MAIN_DIR, 'main.html'))

    def post(self):

        query_type = self.get_argument('search', default=None)

        if 'Gene' in query_type:

            self.render(os.path.join(QUERY_DIR, 'gene_search.html'))

        elif 'Phenotype' in query_type:

            self.render(os.path.join(QUERY_DIR, 'phenotype_search.html'))

        else:

            pass


class QueryHandler(tornado.web.RequestHandler):

    def get(self):

        query_type = self.get_argument('search', default='')

        if 'Gene' in query_type:

            self.render(os.path.join(QUERY_DIR, 'gene_search.html'))

    def post(self):

        single_gene_query = self.get_argument('single_gene', default=None)

        specific_flag = self.get_argument('specific', True)

        specific_flag = True if specific_flag == 'true' else False

        if single_gene_query is None:
            single_gene_query = []
        else:
            single_gene_query = [single_gene_query]

        multi_gene_query = self.get_arguments('gene_list')

        if len(multi_gene_query) == 1:
            multi_gene_query = multi_gene_query[0].split('\n')
            multi_gene_query = [x.strip() for x in multi_gene_query]

        records = []
        records.extend(single_gene_query)
        records.extend(multi_gene_query)

        names, mutant_text_array = resistome_handler.get_mutant_output(records, specific_flag=specific_flag)

        self.render(os.path.join(QUERY_DIR, 'results.html'), records=mutant_text_array, converted_gene_names = names)

def main():

    global MAIN_DIR, QUERY_DIR, resistome_handler

    resistome_handler = ResistomeDBHandler()

    current_location = os.path.realpath(__file__)
    path = os.path.split(current_location)[0]

    STATIC_DIR = os.path.join(path, 'static')
    TEMPLATE_DIR = os.path.join(path, 'templates')
    MAIN_DIR = os.path.join(path, 'main')
    QUERY_DIR = os.path.join(path, 'search')

    handlers = [
        (r"/", MainHandler),
        (r"/query/", QueryHandler),
        (r'/static/(.*)', tornado.web.StaticFileHandler, {'path': STATIC_DIR})]

    application = tornado.web.Application(handlers, debug=True, template_path=TEMPLATE_DIR)

    application.listen(8080)
    tornado.ioloop.IOLoop.instance().start()

if __name__ == "__main__":
    main()

