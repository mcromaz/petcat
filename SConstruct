#dbg = Environment(CCFLAGS = "-pedantic -g -D_FILE_OFFSET_BITS=64 -DTIMESTAMP -DONEINT", LIBS=["m",])
dbg = Environment(CCFLAGS = "-g -D_FILE_OFFSET_BITS=64 -DTIMESTAMP -DSAMPLE25", LIBS=["m",])
#dbg = Environment(CCFLAGS = "-g -D_FILE_OFFSET_BITS=64 -DTIMESTAMP", LIBS=["m",])
dbg.Program(["petcat.c", "preproc.c", "decompose.c", "config.c", "mode3io.c", "fal.c", "talign.c", "utils.c", "grid.c", "fitter.c", "eval.c", "matinv.c", "interpolate.c", "read_basis.c", "diagcnt.c", "log3.c", "log2.c", "log.c", "wrap.c", "lh.c", "oneInt.c"])
#dbg.Program(["sim.c", "decompose.c", "grid.c", "fitter.c", "eval.c", "matinv.c", "interpolate.c", "read_basis.c", "diagcnt.c", "log2.c", "log.c", "wrap.c", "lh.c", "oneInt.c"])
dbg.Program(["vegcat.c", "preproc.c", "talign.c", "sint.c", "decompose.c", "grid.c", "fitter.c", "eval.c", "interpolate.c", "matinv.c", "read_basis.c", "utils.c", "config.c", "fal.c", "diagcnt.c", "log.c", "lh.c", "wrap.c"])
dbg.Program(["cevt.c",])
dbg.Program(["basisPts.c"])
