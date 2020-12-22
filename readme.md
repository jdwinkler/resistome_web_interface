Update instructions:

1. Deposit the SQL dump of the database onto publicly accessible storage, e.g. a presigned AWS S3 link or a random 
file hosting website. Must be accessible to unauthenticated HTTPS requests.
2. Log into the Heroku CLI.
3. Log into the Heroku web interface. Drop the existing postgres database (make sure you have a backup!)
4. Began the restoration process outlined [here](https://devcenter.heroku.com/articles/heroku-postgres-import-export).
5. Make sure the website is accessible and works for simple queries after the restore is complete.

For the restore process, check the Postgres web add-on for the DATABASE_URL to use (should be 
postgres-random name-random number).

Note: code lives [here](https://bitbucket.org/jdwinkler/genome-engineering/src/master/).