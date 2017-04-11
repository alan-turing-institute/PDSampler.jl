# Memo

Change of delimiters

```bash
a=ratings.dat; sed s/::/,/g $a > _tmp; mv _tmp $a
```

```julia
a = readdlm("ratings.dat",',',Int64)
```
