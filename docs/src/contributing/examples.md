# Adding examples

Examples are a great way to showcase how to make use of a specific feature.
We consider two types of examples:

1. a simple example (running in < 30 seconds)
2. a complex example (anything that's not in the first category)

The simple examples are usually to be preferred since they can also be used as
*tests*.
For a big simulation (e.g.: use of LBPS for Probabilistic Matrix Factorisation)
we recommend the use of separate Github repos including notebooks analysing the
results.

A simple example can often be considered as a form of test by itself except with
more documentation in order to clarify what is being done.
The process in this project is therefore

1. write the example as a test (see [this example](https://github.com/alan-turing-institute/PDMP.jl/blob/master/test/ex_gbps1.jl) for reference)
2. run `readexamples.jl` in `doc/src` which will generate a markdown file which
can be served by the documentation website.

## Syntax of the test file

The part of the test you want to appear in the documentation should be between
two tags `#@startexample` and `#@endexample`:

```julia
using Base.Test
#@startexample
...
#@endexample
@test ...
```

Typically the part between the tags will contain no tests and will just showcase
how to use a specific set of features while the part outside of the tags may
ensure that the results obtained are correct.

The text that will be in the documentation should be written as multiline
comments i.e.:

```julia
#=
Everything here will be written as Markdown in the doc.
=#
```

The rest of the code and comments that are in between the tags `#@startexample`
and `#@endexample` will be copied in the doc as code blocks.
A very simple example would be:

```julia
using Base.Test
#@startexample
#=
# Name of the Test

In this example we will show how to find the maximum over a collection of random
numbers in Julia for all the numbers that are negative.
=#
a = randn(500)
# we use a comprehension
m = maximum(a[i] for i in 1:length(a) if a[i]<0)
#=
That's it!
=#
#@endexample
@test m < 0.0
```

This will be converted to Markdown:

```markdown
# Name of the Test

In this example we will show how to find the maximum over a collection of random
numbers in Julia for all the numbers that are negative.
 ```julia
a = randn(500)
# we use a comprehension
m = maximum(a[i] for i in 1:length(a) if a[i]<0)
 ```
That's it!
```

*Remark*, the spaces in front of the triple backticks in the snippet above are
not actually generated when you use `readexamples.jl`.
The spaces are used here in order to escape these triple backticks so that the
snippet does not end up being unduly fragmented in three pieces.

## Full workflow

1. Write the example in the `PDMP/test/` directory following the examples
present there and using the syntax mentioned above.
Use the name convention `ex_NAME.jl`.
2. Make sure this runs and tests pass.
3. Add the test at the bottom of `PDMP/test/runtests.jl` following code already
present there.
4. Navigate to `PDMP/docs` in your terminal and `julia readexamples`. This will
generate the `.md` file corresponding to your example and place it in
`PDMP/docs/src/examples/`.
It will also do the same with all other examples (effectively refreshing those).
5. In `PDMP/docs/make.jl` add the page name `ex_NAME.md` in the list under
`Examples` following code already present there.
6. Git add, commit and push.
