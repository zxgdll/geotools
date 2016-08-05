#!/bin/sh

echo "Running tests..."

for f in test_*.py;
do
	python $f
	if [ $? -ne 0 ];
	then
		return 1
	fi
done

echo "Done"

return 0
