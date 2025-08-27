# Overview of Commands

There are two sets of options:

1. options for mpirun
    - these directly follow the mpirun command
2. options for the R script
    - these follow the script's name: [smmala.R](./smmala.R)
## Perform MCMC Sampling

Here is a command that can be used to perform sampling, with explanations:

```sh
mpirun -H localhost:8 -N 8 ./smmala.R -N 10000 -c 6 --prefix test | tee smmala.log
```

Of course mpirun is one of several alternative MPI implementations,
like OpenMPI. The option `-H` declares which machines to use, with
`:8` setting how many MPI workers the machine can run (max), and `-N`
sets the number of worker-slots to actually use (per-machine).

The `smmala.R` script itself takes several arguments:

- `-N` size of sample, the script will take several smaller samples first; the smaller samples are used to:
    + tweak the acceptance rate
    + perform a burn-in to allow the Markov chains to converge to th target distribution
- `-c` how many bunches of samples to take, i.e. cycles of MCMC to perform
    + the first half will be used to tune the acceptance rate
    + the last of these will have the requested sample size
    + the first of these will be very small
    + sample sizes increase on a logarithmic scale between the first and last sample
- `--prefix` this string will be used in all file-names

## Reasonable Choices

A big `-c` value can help on machines with small random access memory. A big overall sample can be achieved in two ways:

Many small samples:

```sh
mpirun -H localhost:8 -N 8 ./smmala.R -N 1000 -c 60 --prefix many-files | tee smmala.log
```

Few, but big samples:

```sh
mpirun -H localhost:8 -N 8 ./smmala.R -N 60000 -c 4 --prefix few-files | tee smmala.log
```

## Find Number of Cores

On GNU/Linux, the following command can be used to find the number of cpu-cores:

```sh
N=$(grep -c processor /proc/cpuinfo)
```

So, the above command changes to the more automatic:

```sh
N=$(grep -c processor /proc/cpuinfo)
mpirun -H localhost:$((N)) -N $((N)) ./smmala.R -N 10000 -c 8 --prefix test | tee smmala.log
```

## TEE

The `tee` command will write the output of `smmala.R` to the screen,
but also to the given file (`smmala.log`), so that this information is
not lost afterwards.


# After Sampling

In plots, we typically don't want to draw more lines than necessary to
represent the distribution. For this reason we want to determine the
auto-correlation length of the sample, and thin it out to its effective
size.

We simulate that smaller sample.

For convenience, the sample will be ordered by log-likelihood
(decreasing), so that the maximum-likelihood estimate is the first
row, and the first simulation:

```R
y[[i]]$func[,,1]
```

is the maximum-likelihood output of the simulation of experiment _i_.

## One Big Sample File

The very last cycle of MCMC (it has the maximal, requested sample
size) is saved to a dedicated file, using this pattern:

```
PREFIX-final-sample-RANK-MPI_COMM_SIZE-lb100-LOGBETA.RDS`
```

This can be loaded to postprocess the sample:

```R
source("post.R")
x <- readRDS("test-final-sample-0-8-lb100-0.RDS")
res <- sampleQuality(x)
```

Where we had 8 MPI workers (8 cores) and LOG-BETA is 0 (\(\beta=1\)).

## Many Small Sample Files

If the number of cycles `-c` is very big, then it is better to use _all_
or _most_ generated samples. For this purpose we call `manyFilesSample()`.

The script [post.R](./post.R) will load a function into R which loads
the generated sample for plotting. It determines the auto-correlation
length via [hadron](https://github.com/HISKP-LQCD/hadron) if it is
installed and cruedly guessed if it isn't. (TODO: alternative method
to determine auto-correltation).

It will also simulate the model, using that sample and make a few
plots regarding the sample's quality:

```R
source("post.R")
files <- dir(pattern="^test-main-branch-Sample") # or whatever prefix you chose
res <- manyFilesSample(files,beta=1.0)       # or pick a different beta

print(comment(res))
print(res$tau)
print(names(res))

print(dim(res$sample))
```
