# 5G-NR-LDPC
5G NR LDPC C++ implementation

C++ 11 or newer version needed to run the examples and use the code. No additional package is used, only standard library.

A good starting point is to run `example_run.m`

This implementation uses min-sum offset decoding algorithms, read `LDPC_decoder_help_doc.pdf` for detail.

Any questions, eamil juquan.justin.mao@gmail.com 

### Compling issue

* sprintf_s was not declared in this scope
   * please make sure c++ 11 or newer version used
   * refer to https://stackoverflow.com/questions/4828228/sprintf-s-was-not-declared-in-this-scope for more information
