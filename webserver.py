# todo: add tornado skeleton

import tornado.ioloop
import tornado.web
import os.path
from email.mime.multipart import MIMEMultipart
from email.mime.text import MIMEText
from collections import defaultdict
from threading import Thread

MAIN_DIR = ''
QUERY_DIR = ''


class MainHandler(tornado.web.RequestHandler):
    def get(self):
        self.render(os.path.join(MAIN_DIR, 'main.html'))


class QueryHandler(tornado.web.RequestHandler):

    def post(self):

        single_gene_query = self.get_argument('single_gene', default=None)

        if single_gene_query is None:
            single_gene_query = []
        else:
            single_gene_query = [single_gene_query]

        multi_gene_query = self.get_arguments('gene_list')

        records = []
        records.extend(single_gene_query)
        records.extend(multi_gene_query)

        self.render(os.path.join(QUERY_DIR, 'results.html'), gene_query=records)

def main():

    global MAIN_DIR, QUERY_DIR

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

