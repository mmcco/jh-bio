# A Bash two-liner to defy gc's attempt to force me to use tabs by removing the -tabs option from gofmt.

gofmt -s=true -w=true repeatgenome/*.go
perl -pi -e 's/\t/    /g' repeatgenome/*.go
