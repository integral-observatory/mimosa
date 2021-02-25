| | |
| :-- | :-- |
| Author |  Aleksandra Gros, CEA |
| Adopted for OSA and maintained by | Volodymyr Savchenko, ISDC |


To build, follow [standard OSA procedure](https://www.isdc.unige.ch/integral/download/osa/doc/11.1/osa_inst_guide.pdf). Note that OSA is has to be installed and initialized prior to building the component. ISDC_ENV is set by OSA init. 

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
