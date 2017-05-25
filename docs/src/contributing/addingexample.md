# Adding an example

Examples are a great way to showcase how to make use of a specific feature. We consider two types of examples:

1. a simple example (running in < 10 seconds)
2. a complex example (the complement of the first category)

The first category is great as they can be used as *tests* and as documentation.

The process is somewhat automated here and essentially all you have to do is write the example in the test directory and comment it accordingly, we show this below.

### Syntax for the example

Let's say you have an example which can be run in a few seconds and uses a new feature. You can effectively use it as a unit test by itself.
To respect conventions, please name your example `ex_NAME.jl` and put it in `test/`.

Start your code by

```julia
using Base.Test
```

Then a few markers should be considered:

* `#@startexample NAME_OF_EXAMPLE` indicates that you start the code of the example proper, everything between that mark and the end flag will appear in the doc.
* Encapsulates all the explanations you want to appear in markdown between `#=` and `=#`, all the other comments will be taken as part of the julia code and shown in code blocks.
* `#@endexample` to indicate that the example is finished
* write a few tests that check that the example produces the right answer (unit test)

So it should look like (a full example can be seen [here](https://github.com/alan-turing-institute/PDMP.jl/blob/master/test/ex_gbps1.jl))

```julia
using Base.Test
#@startexample A simple example
#=
In this example we will show how to find the maximum over a collection of
random numbers in Julia for all the numbers that are *negative*.
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

Make sure the tests pass! This will generate the following markdown (see next point for the command that generates it):

```markdown
# Name of the Test

In this example we will show how to find the maximum over a collection of
random numbers in Julia for all the numbers that are negative.
 ```julia
a = randn(500)
# we use a comprehension
m = maximum(a[i] for i in 1:length(a) if a[i]<0)
 ```
That's it!
```

*Remark*, the spaces in front of the triple backticks in the markdown snippet above are not actually generated when you use `readexamples.jl`.
The spaces are used here in order to escape these triple backticks so that the
snippet does not end up being unduly fragmented in three pieces.

### Declaring your example

You have to mention your example in a few spots:

* in `test/runtests.jl`, add a line at the bottom following the examples already present i.e. something like `@testset "ex_NAME"    begin include("ex_NAME.jl") end` (make sure this passes!)
* in `docs/src/make.jl`, add a line under `"Examples"` following the syntax of examples already present so something like: `"Name of your expl" => "examples/ex_name_of_example.jl"`

Finally, to generate the markdown run the following command (this will act for all the examples at once so it will also refresh any other modification you may have added to other examples):

```
julia docs/readexamples.jl
```
