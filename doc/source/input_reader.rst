.. module:: source.input_reader

Input Reader
============

Class for parsing and storing the values in the input file.

All keywords are listed in 'dict' along with default values.

Keywords that take a tuple instead take a string of space
separated numbers.

the ideal line in the input file is:

keyword = value

though a '=' between the keyword and the value is the only
required part.

The class is initialized without arguments.  The dictionary
is populated by calling the read_input function with the name
and path of the input file.

____________

.. autoclass:: ReadInput
   :members: read_input

