To build, standard ISGR

```
$ $ISDC_ENV/ac_stuff/configure && make
```

You can also install it in your `ISDC_ENV`, if it works for you:

```
$ make install
```


To run:

```
./mimosa
```

Alternatively:

```
rm out_cat.fits out_mosa_*; mimosa inCat="$ISDC_REF_CAT[ISGRI_FLAG==1]"
```
