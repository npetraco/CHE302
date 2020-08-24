# Handy Makefile tocleanup and push changes out to MAT301 website 

default:
	cp pchem2.html index.html
	git add --all
	git commit -m "update"
	git push -u origin master

copy:
	cp pchem2.html index.html
