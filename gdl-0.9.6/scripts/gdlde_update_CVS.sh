#
# AC 01/08/2015
# update an already existing CVS with sub-directories 
# should be run at top of the root of the CVS (e.g. gdlde/)
#
echo 'This is for the GDLDE part : IDE for GDL'
echo ''
#
echo 'just press enter for the passwd !'
echo ''
#
CVS_PATH=/cvsroot/gnudatalanguage
CVS_SITE=anonymous@gnudatalanguage.cvs.sourceforge.net
#
cvs -d:pserver:$CVS_SITE:$CVS_PATH login
cvs -z3 -d:pserver:$CVS_SITE:$CVS_PATH update -d gdlde
#
