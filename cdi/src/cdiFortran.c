/* Automatically generated by make_fint.c, don't edit! */

#if defined (HAVE_CONFIG_H)
#  include "config.h"
#endif

#if ! defined (CDI_H_)
#  include "cdi.h"
#endif

#if defined (HAVE_CF_INTERFACE)

#if ! defined (__CFORTRAN_LOADED)
#  include "cfortran.h"
#endif


/*  Byte order  */


/*  Error identifier  */


/*  File types  */


/*  Compress types  */


/*  external data types  */


/*  internal data types  */


/*  Chunks  */


/*  GRID types  */


/*  ZAXIS types  */


/*  SUBTYPE types  */


/*  Data structure defining a key-value search, possibly with multiple
   key-value pairs in combination.

   Currently, only multiple pairs combined by AND are supported.  */


/*  TIME types  */


/*  TSTEP types  */


/*  TAXIS types  */


/*  TUNIT types  */


/*  CALENDAR types  */


/*  number of unsigned char needed to store UUID  */


/*  Structs that are used to return data to the user  */


/*  Opaque types  */


/*  CDI control routines  */

FCALLSCSUB0 (cdiReset, CDIRESET, cdireset)
FCALLSCFUN1 (STRING, cdiStringError, CDISTRINGERROR, cdistringerror, INT)
FCALLSCSUB1 (cdiDebug, CDIDEBUG, cdidebug, INT)
FCALLSCFUN0 (STRING, cdiLibraryVersion, CDILIBRARYVERSION, cdilibraryversion)
FCALLSCSUB0 (cdiPrintVersion, CDIPRINTVERSION, cdiprintversion)
FCALLSCFUN1 (INT, cdiHaveFiletype, CDIHAVEFILETYPE, cdihavefiletype, INT)
FCALLSCSUB1 (cdiDefMissval, CDIDEFMISSVAL, cdidefmissval, DOUBLE)
FCALLSCFUN0 (DOUBLE, cdiInqMissval, CDIINQMISSVAL, cdiinqmissval)
FCALLSCSUB2 (cdiDefGlobal, CDIDEFGLOBAL, cdidefglobal, STRING, INT)
FCALLSCFUN0 (INT, namespaceNew, NAMESPACENEW, namespacenew)
FCALLSCSUB1 (namespaceSetActive, NAMESPACESETACTIVE, namespacesetactive, INT)
FCALLSCFUN0 (INT, namespaceGetActive, NAMESPACEGETACTIVE, namespacegetactive)
FCALLSCSUB1 (namespaceDelete, NAMESPACEDELETE, namespacedelete, INT)

/*  CDI converter routines  */


/*  parameter  */

FCALLSCSUB3 (cdiParamToString, CDIPARAMTOSTRING, cdiparamtostring, INT, PSTRING, INT)
FCALLSCSUB4 (cdiDecodeParam, CDIDECODEPARAM, cdidecodeparam, INT, PINT, PINT, PINT)
FCALLSCFUN3 (INT, cdiEncodeParam, CDIENCODEPARAM, cdiencodeparam, INT, INT, INT)

/*  date format:  YYYYMMDD  */


/*  time format:    hhmmss  */

FCALLSCSUB4 (cdiDecodeDate, CDIDECODEDATE, cdidecodedate, INT, PINT, PINT, PINT)
FCALLSCFUN3 (INT, cdiEncodeDate, CDIENCODEDATE, cdiencodedate, INT, INT, INT)
FCALLSCSUB4 (cdiDecodeTime, CDIDECODETIME, cdidecodetime, INT, PINT, PINT, PINT)
FCALLSCFUN3 (INT, cdiEncodeTime, CDIENCODETIME, cdiencodetime, INT, INT, INT)

/*  STREAM control routines  */

FCALLSCFUN2 (INT, cdiGetFiletype, CDIGETFILETYPE, cdigetfiletype, STRING, PINT)
FCALLSCFUN1 (INT, streamOpenRead, STREAMOPENREAD, streamopenread, STRING)
FCALLSCFUN2 (INT, streamOpenWrite, STREAMOPENWRITE, streamopenwrite, STRING, INT)
FCALLSCFUN1 (INT, streamOpenAppend, STREAMOPENAPPEND, streamopenappend, STRING)
FCALLSCSUB1 (streamClose, STREAMCLOSE, streamclose, INT)
FCALLSCSUB1 (streamSync, STREAMSYNC, streamsync, INT)
FCALLSCSUB2 (streamDefVlist, STREAMDEFVLIST, streamdefvlist, INT, INT)
FCALLSCFUN1 (INT, streamInqVlist, STREAMINQVLIST, streaminqvlist, INT)
FCALLSCFUN1 (INT, streamInqFiletype, STREAMINQFILETYPE, streaminqfiletype, INT)
FCALLSCSUB2 (streamDefByteorder, STREAMDEFBYTEORDER, streamdefbyteorder, INT, INT)
FCALLSCFUN1 (INT, streamInqByteorder, STREAMINQBYTEORDER, streaminqbyteorder, INT)
FCALLSCSUB2 (streamDefCompType, STREAMDEFCOMPTYPE, streamdefcomptype, INT, INT)
FCALLSCFUN1 (INT, streamInqCompType, STREAMINQCOMPTYPE, streaminqcomptype, INT)
FCALLSCSUB2 (streamDefCompLevel, STREAMDEFCOMPLEVEL, streamdefcomplevel, INT, INT)
FCALLSCFUN1 (INT, streamInqCompLevel, STREAMINQCOMPLEVEL, streaminqcomplevel, INT)
FCALLSCFUN2 (INT, streamDefTimestep, STREAMDEFTIMESTEP, streamdeftimestep, INT, INT)
FCALLSCFUN2 (INT, streamInqTimestep, STREAMINQTIMESTEP, streaminqtimestep, INT, INT)
FCALLSCFUN1 (INT, streamInqCurTimestepID, STREAMINQCURTIMESTEPID, streaminqcurtimestepid, INT)
FCALLSCFUN1 (STRING, streamFilename, STREAMFILENAME, streamfilename, INT)
FCALLSCFUN1 (STRING, streamFilesuffix, STREAMFILESUFFIX, streamfilesuffix, INT)
FCALLSCFUN1 (INT, streamInqNvars, STREAMINQNVARS, streaminqnvars, INT)

/*  STREAM var I/O routines  */

FCALLSCSUB4 (streamWriteVar, STREAMWRITEVAR, streamwritevar, INT, INT, DOUBLEV, INT)
FCALLSCSUB4 (streamWriteVarF, STREAMWRITEVARF, streamwritevarf, INT, INT, FLOATV, INT)
FCALLSCSUB4 (streamReadVar, STREAMREADVAR, streamreadvar, INT, INT, DOUBLEV, PINT)
FCALLSCSUB4 (streamReadVarF, STREAMREADVARF, streamreadvarf, INT, INT, FLOATV, PINT)
FCALLSCSUB5 (streamWriteVarSlice, STREAMWRITEVARSLICE, streamwritevarslice, INT, INT, INT, DOUBLEV, INT)
FCALLSCSUB5 (streamWriteVarSliceF, STREAMWRITEVARSLICEF, streamwritevarslicef, INT, INT, INT, FLOATV, INT)
FCALLSCSUB5 (streamReadVarSlice, STREAMREADVARSLICE, streamreadvarslice, INT, INT, INT, DOUBLEV, PINT)
FCALLSCSUB5 (streamReadVarSliceF, STREAMREADVARSLICEF, streamreadvarslicef, INT, INT, INT, FLOATV, PINT)
FCALLSCSUB5 (streamWriteVarChunk, STREAMWRITEVARCHUNK, streamwritevarchunk, INT, INT, INTVV, DOUBLEV, INT)

/*  STREAM record I/O routines  */

FCALLSCSUB3 (streamDefRecord, STREAMDEFRECORD, streamdefrecord, INT, INT, INT)
FCALLSCSUB3 (streamInqRecord, STREAMINQRECORD, streaminqrecord, INT, PINT, PINT)
FCALLSCSUB3 (streamWriteRecord, STREAMWRITERECORD, streamwriterecord, INT, DOUBLEV, INT)
FCALLSCSUB3 (streamWriteRecordF, STREAMWRITERECORDF, streamwriterecordf, INT, FLOATV, INT)
FCALLSCSUB3 (streamReadRecord, STREAMREADRECORD, streamreadrecord, INT, DOUBLEV, PINT)
FCALLSCSUB2 (streamCopyRecord, STREAMCOPYRECORD, streamcopyrecord, INT, INT)

/*  File driven I/O (may yield better performance than using the streamXXX functions)  */


/*  VLIST routines  */

FCALLSCFUN0 (INT, vlistCreate, VLISTCREATE, vlistcreate)
FCALLSCSUB1 (vlistDestroy, VLISTDESTROY, vlistdestroy, INT)
FCALLSCFUN1 (INT, vlistDuplicate, VLISTDUPLICATE, vlistduplicate, INT)
FCALLSCSUB2 (vlistCopy, VLISTCOPY, vlistcopy, INT, INT)
FCALLSCSUB2 (vlistCopyFlag, VLISTCOPYFLAG, vlistcopyflag, INT, INT)
FCALLSCSUB1 (vlistClearFlag, VLISTCLEARFLAG, vlistclearflag, INT)
FCALLSCSUB2 (vlistCat, VLISTCAT, vlistcat, INT, INT)
FCALLSCSUB2 (vlistMerge, VLISTMERGE, vlistmerge, INT, INT)
FCALLSCSUB1 (vlistPrint, VLISTPRINT, vlistprint, INT)
FCALLSCFUN1 (INT, vlistNumber, VLISTNUMBER, vlistnumber, INT)
FCALLSCFUN1 (INT, vlistNvars, VLISTNVARS, vlistnvars, INT)
FCALLSCFUN1 (INT, vlistNgrids, VLISTNGRIDS, vlistngrids, INT)
FCALLSCFUN1 (INT, vlistNzaxis, VLISTNZAXIS, vlistnzaxis, INT)
FCALLSCFUN1 (INT, vlistNsubtypes, VLISTNSUBTYPES, vlistnsubtypes, INT)
FCALLSCSUB2 (vlistDefNtsteps, VLISTDEFNTSTEPS, vlistdefntsteps, INT, INT)
FCALLSCFUN1 (INT, vlistNtsteps, VLISTNTSTEPS, vlistntsteps, INT)
FCALLSCFUN1 (INT, vlistGridsizeMax, VLISTGRIDSIZEMAX, vlistgridsizemax, INT)
FCALLSCFUN2 (INT, vlistGrid, VLISTGRID, vlistgrid, INT, INT)
FCALLSCFUN2 (INT, vlistGridIndex, VLISTGRIDINDEX, vlistgridindex, INT, INT)
FCALLSCSUB3 (vlistChangeGridIndex, VLISTCHANGEGRIDINDEX, vlistchangegridindex, INT, INT, INT)
FCALLSCSUB3 (vlistChangeGrid, VLISTCHANGEGRID, vlistchangegrid, INT, INT, INT)
FCALLSCFUN2 (INT, vlistZaxis, VLISTZAXIS, vlistzaxis, INT, INT)
FCALLSCFUN2 (INT, vlistZaxisIndex, VLISTZAXISINDEX, vlistzaxisindex, INT, INT)
FCALLSCSUB3 (vlistChangeZaxisIndex, VLISTCHANGEZAXISINDEX, vlistchangezaxisindex, INT, INT, INT)
FCALLSCSUB3 (vlistChangeZaxis, VLISTCHANGEZAXIS, vlistchangezaxis, INT, INT, INT)
FCALLSCFUN1 (INT, vlistNrecs, VLISTNRECS, vlistnrecs, INT)
FCALLSCFUN2 (INT, vlistSubtype, VLISTSUBTYPE, vlistsubtype, INT, INT)
FCALLSCFUN2 (INT, vlistSubtypeIndex, VLISTSUBTYPEINDEX, vlistsubtypeindex, INT, INT)
FCALLSCSUB2 (vlistDefTaxis, VLISTDEFTAXIS, vlistdeftaxis, INT, INT)
FCALLSCFUN1 (INT, vlistInqTaxis, VLISTINQTAXIS, vlistinqtaxis, INT)
FCALLSCSUB2 (vlistDefTable, VLISTDEFTABLE, vlistdeftable, INT, INT)
FCALLSCFUN1 (INT, vlistInqTable, VLISTINQTABLE, vlistinqtable, INT)
FCALLSCSUB2 (vlistDefInstitut, VLISTDEFINSTITUT, vlistdefinstitut, INT, INT)
FCALLSCFUN1 (INT, vlistInqInstitut, VLISTINQINSTITUT, vlistinqinstitut, INT)
FCALLSCSUB2 (vlistDefModel, VLISTDEFMODEL, vlistdefmodel, INT, INT)
FCALLSCFUN1 (INT, vlistInqModel, VLISTINQMODEL, vlistinqmodel, INT)

/*  VLIST VAR routines  */

FCALLSCFUN5 (INT, vlistDefVarTiles, VLISTDEFVARTILES, vlistdefvartiles, INT, INT, INT, INT, INT)
FCALLSCFUN4 (INT, vlistDefVar, VLISTDEFVAR, vlistdefvar, INT, INT, INT, INT)
FCALLSCSUB3 (vlistChangeVarGrid, VLISTCHANGEVARGRID, vlistchangevargrid, INT, INT, INT)
FCALLSCSUB3 (vlistChangeVarZaxis, VLISTCHANGEVARZAXIS, vlistchangevarzaxis, INT, INT, INT)
FCALLSCSUB5 (vlistInqVar, VLISTINQVAR, vlistinqvar, INT, INT, PINT, PINT, PINT)
FCALLSCFUN2 (INT, vlistInqVarGrid, VLISTINQVARGRID, vlistinqvargrid, INT, INT)
FCALLSCFUN2 (INT, vlistInqVarZaxis, VLISTINQVARZAXIS, vlistinqvarzaxis, INT, INT)

/*  used in MPIOM  */

FCALLSCFUN2 (INT, vlistInqVarID, VLISTINQVARID, vlistinqvarid, INT, INT)
FCALLSCSUB3 (vlistDefVarTsteptype, VLISTDEFVARTSTEPTYPE, vlistdefvartsteptype, INT, INT, INT)
FCALLSCFUN2 (INT, vlistInqVarTsteptype, VLISTINQVARTSTEPTYPE, vlistinqvartsteptype, INT, INT)
FCALLSCSUB3 (vlistDefVarCompType, VLISTDEFVARCOMPTYPE, vlistdefvarcomptype, INT, INT, INT)
FCALLSCFUN2 (INT, vlistInqVarCompType, VLISTINQVARCOMPTYPE, vlistinqvarcomptype, INT, INT)
FCALLSCSUB3 (vlistDefVarCompLevel, VLISTDEFVARCOMPLEVEL, vlistdefvarcomplevel, INT, INT, INT)
FCALLSCFUN2 (INT, vlistInqVarCompLevel, VLISTINQVARCOMPLEVEL, vlistinqvarcomplevel, INT, INT)
FCALLSCSUB3 (vlistDefVarParam, VLISTDEFVARPARAM, vlistdefvarparam, INT, INT, INT)
FCALLSCFUN2 (INT, vlistInqVarParam, VLISTINQVARPARAM, vlistinqvarparam, INT, INT)
FCALLSCSUB3 (vlistDefVarCode, VLISTDEFVARCODE, vlistdefvarcode, INT, INT, INT)
FCALLSCFUN2 (INT, vlistInqVarCode, VLISTINQVARCODE, vlistinqvarcode, INT, INT)
FCALLSCSUB3 (vlistDefVarDatatype, VLISTDEFVARDATATYPE, vlistdefvardatatype, INT, INT, INT)
FCALLSCFUN2 (INT, vlistInqVarDatatype, VLISTINQVARDATATYPE, vlistinqvardatatype, INT, INT)
FCALLSCSUB3 (vlistDefVarChunkType, VLISTDEFVARCHUNKTYPE, vlistdefvarchunktype, INT, INT, INT)
FCALLSCFUN2 (INT, vlistInqVarChunkType, VLISTINQVARCHUNKTYPE, vlistinqvarchunktype, INT, INT)
FCALLSCSUB3 (vlistDefVarXYZ, VLISTDEFVARXYZ, vlistdefvarxyz, INT, INT, INT)
FCALLSCFUN2 (INT, vlistInqVarXYZ, VLISTINQVARXYZ, vlistinqvarxyz, INT, INT)
FCALLSCFUN2 (INT, vlistInqVarNumber, VLISTINQVARNUMBER, vlistinqvarnumber, INT, INT)
FCALLSCSUB3 (vlistDefVarInstitut, VLISTDEFVARINSTITUT, vlistdefvarinstitut, INT, INT, INT)
FCALLSCFUN2 (INT, vlistInqVarInstitut, VLISTINQVARINSTITUT, vlistinqvarinstitut, INT, INT)
FCALLSCSUB3 (vlistDefVarModel, VLISTDEFVARMODEL, vlistdefvarmodel, INT, INT, INT)
FCALLSCFUN2 (INT, vlistInqVarModel, VLISTINQVARMODEL, vlistinqvarmodel, INT, INT)
FCALLSCSUB3 (vlistDefVarTable, VLISTDEFVARTABLE, vlistdefvartable, INT, INT, INT)
FCALLSCFUN2 (INT, vlistInqVarTable, VLISTINQVARTABLE, vlistinqvartable, INT, INT)
FCALLSCSUB3 (vlistDefVarName, VLISTDEFVARNAME, vlistdefvarname, INT, INT, STRING)
FCALLSCSUB3 (vlistInqVarName, VLISTINQVARNAME, vlistinqvarname, INT, INT, PSTRING)
FCALLSCFUN2 (STRING, vlistCopyVarName, VLISTCOPYVARNAME, vlistcopyvarname, INT, INT)
FCALLSCSUB3 (vlistDefVarStdname, VLISTDEFVARSTDNAME, vlistdefvarstdname, INT, INT, STRING)
FCALLSCSUB3 (vlistInqVarStdname, VLISTINQVARSTDNAME, vlistinqvarstdname, INT, INT, PSTRING)
FCALLSCSUB3 (vlistDefVarLongname, VLISTDEFVARLONGNAME, vlistdefvarlongname, INT, INT, STRING)
FCALLSCSUB3 (vlistInqVarLongname, VLISTINQVARLONGNAME, vlistinqvarlongname, INT, INT, PSTRING)
FCALLSCSUB3 (vlistDefVarUnits, VLISTDEFVARUNITS, vlistdefvarunits, INT, INT, STRING)
FCALLSCSUB3 (vlistInqVarUnits, VLISTINQVARUNITS, vlistinqvarunits, INT, INT, PSTRING)
FCALLSCSUB3 (vlistDefVarMissval, VLISTDEFVARMISSVAL, vlistdefvarmissval, INT, INT, DOUBLE)
FCALLSCFUN2 (DOUBLE, vlistInqVarMissval, VLISTINQVARMISSVAL, vlistinqvarmissval, INT, INT)
FCALLSCSUB3 (vlistDefVarExtra, VLISTDEFVAREXTRA, vlistdefvarextra, INT, INT, STRING)
FCALLSCSUB3 (vlistInqVarExtra, VLISTINQVAREXTRA, vlistinqvarextra, INT, INT, PSTRING)
FCALLSCSUB3 (vlistDefVarScalefactor, VLISTDEFVARSCALEFACTOR, vlistdefvarscalefactor, INT, INT, DOUBLE)
FCALLSCFUN2 (DOUBLE, vlistInqVarScalefactor, VLISTINQVARSCALEFACTOR, vlistinqvarscalefactor, INT, INT)
FCALLSCSUB3 (vlistDefVarAddoffset, VLISTDEFVARADDOFFSET, vlistdefvaraddoffset, INT, INT, DOUBLE)
FCALLSCFUN2 (DOUBLE, vlistInqVarAddoffset, VLISTINQVARADDOFFSET, vlistinqvaraddoffset, INT, INT)
FCALLSCSUB3 (vlistDefVarTimave, VLISTDEFVARTIMAVE, vlistdefvartimave, INT, INT, INT)
FCALLSCFUN2 (INT, vlistInqVarTimave, VLISTINQVARTIMAVE, vlistinqvartimave, INT, INT)
FCALLSCSUB3 (vlistDefVarTimaccu, VLISTDEFVARTIMACCU, vlistdefvartimaccu, INT, INT, INT)
FCALLSCFUN2 (INT, vlistInqVarTimaccu, VLISTINQVARTIMACCU, vlistinqvartimaccu, INT, INT)
FCALLSCSUB3 (vlistDefVarTypeOfGeneratingProcess, VLISTDEFVARTYPEOFGENERATINGPROCESS, vlistdefvartypeofgeneratingprocess, INT, INT, INT)
FCALLSCFUN2 (INT, vlistInqVarTypeOfGeneratingProcess, VLISTINQVARTYPEOFGENERATINGPROCESS, vlistinqvartypeofgeneratingprocess, INT, INT)
FCALLSCSUB3 (vlistDefVarProductDefinitionTemplate, VLISTDEFVARPRODUCTDEFINITIONTEMPLATE, vlistdefvarproductdefinitiontemplate, INT, INT, INT)
FCALLSCFUN2 (INT, vlistInqVarProductDefinitionTemplate, VLISTINQVARPRODUCTDEFINITIONTEMPLATE, vlistinqvarproductdefinitiontemplate, INT, INT)
FCALLSCFUN2 (INT, vlistInqVarSize, VLISTINQVARSIZE, vlistinqvarsize, INT, INT)
FCALLSCSUB4 (vlistDefIndex, VLISTDEFINDEX, vlistdefindex, INT, INT, INT, INT)
FCALLSCFUN3 (INT, vlistInqIndex, VLISTINQINDEX, vlistinqindex, INT, INT, INT)
FCALLSCSUB4 (vlistDefFlag, VLISTDEFFLAG, vlistdefflag, INT, INT, INT, INT)
FCALLSCFUN3 (INT, vlistInqFlag, VLISTINQFLAG, vlistinqflag, INT, INT, INT)
FCALLSCFUN2 (INT, vlistFindVar, VLISTFINDVAR, vlistfindvar, INT, INT)
FCALLSCFUN3 (INT, vlistFindLevel, VLISTFINDLEVEL, vlistfindlevel, INT, INT, INT)
FCALLSCFUN2 (INT, vlistMergedVar, VLISTMERGEDVAR, vlistmergedvar, INT, INT)
FCALLSCFUN3 (INT, vlistMergedLevel, VLISTMERGEDLEVEL, vlistmergedlevel, INT, INT, INT)

/*  Ensemble info routines  */

FCALLSCSUB5 (vlistDefVarEnsemble, VLISTDEFVARENSEMBLE, vlistdefvarensemble, INT, INT, INT, INT, INT)
FCALLSCFUN5 (INT, vlistInqVarEnsemble, VLISTINQVARENSEMBLE, vlistinqvarensemble, INT, INT, PINT, PINT, PINT)
FCALLSCSUB0 (cdiClearAdditionalKeys, CDICLEARADDITIONALKEYS, cdiclearadditionalkeys)
FCALLSCSUB1 (cdiDefAdditionalKey, CDIDEFADDITIONALKEY, cdidefadditionalkey, STRING)
FCALLSCSUB4 (vlistDefVarIntKey, VLISTDEFVARINTKEY, vlistdefvarintkey, INT, INT, STRING, INT)
FCALLSCSUB4 (vlistDefVarDblKey, VLISTDEFVARDBLKEY, vlistdefvardblkey, INT, INT, STRING, DOUBLE)
FCALLSCFUN3 (INT, vlistHasVarKey, VLISTHASVARKEY, vlisthasvarkey, INT, INT, STRING)
FCALLSCFUN3 (DOUBLE, vlistInqVarDblKey, VLISTINQVARDBLKEY, vlistinqvardblkey, INT, INT, STRING)
FCALLSCFUN3 (INT, vlistInqVarIntKey, VLISTINQVARINTKEY, vlistinqvarintkey, INT, INT, STRING)

/*  VLIST attributes  */

FCALLSCFUN3 (INT, vlistInqNatts, VLISTINQNATTS, vlistinqnatts, INT, INT, PINT)
FCALLSCFUN6 (INT, vlistInqAtt, VLISTINQATT, vlistinqatt, INT, INT, INT, PSTRING, PINT, PINT)
FCALLSCFUN3 (INT, vlistDelAtt, VLISTDELATT, vlistdelatt, INT, INT, STRING)
FCALLSCFUN6 (INT, vlistDefAttInt, VLISTDEFATTINT, vlistdefattint, INT, INT, STRING, INT, INT, INTV)
FCALLSCFUN6 (INT, vlistDefAttFlt, VLISTDEFATTFLT, vlistdefattflt, INT, INT, STRING, INT, INT, DOUBLEV)
FCALLSCFUN5 (INT, vlistDefAttTxt, VLISTDEFATTTXT, vlistdefatttxt, INT, INT, STRING, INT, PPSTRING)
FCALLSCFUN5 (INT, vlistInqAttInt, VLISTINQATTINT, vlistinqattint, INT, INT, STRING, INT, INTV)
FCALLSCFUN5 (INT, vlistInqAttFlt, VLISTINQATTFLT, vlistinqattflt, INT, INT, STRING, INT, DOUBLEV)
FCALLSCFUN5 (INT, vlistInqAttTxt, VLISTINQATTTXT, vlistinqatttxt, INT, INT, STRING, INT, PPSTRING)

/*  GRID routines  */

FCALLSCSUB2 (gridName, GRIDNAME, gridname, INT, PSTRING)
FCALLSCFUN1 (STRING, gridNamePtr, GRIDNAMEPTR, gridnameptr, INT)
FCALLSCSUB1 (gridCompress, GRIDCOMPRESS, gridcompress, INT)
FCALLSCSUB2 (gridDefMaskGME, GRIDDEFMASKGME, griddefmaskgme, INT, INTV)
FCALLSCFUN2 (INT, gridInqMaskGME, GRIDINQMASKGME, gridinqmaskgme, INT, INTV)
FCALLSCSUB2 (gridDefMask, GRIDDEFMASK, griddefmask, INT, INTV)
FCALLSCFUN2 (INT, gridInqMask, GRIDINQMASK, gridinqmask, INT, INTV)
FCALLSCSUB3 (gridPrint, GRIDPRINT, gridprint, INT, INT, INT)
FCALLSCFUN2 (INT, gridCreate, GRIDCREATE, gridcreate, INT, INT)
FCALLSCSUB1 (gridDestroy, GRIDDESTROY, griddestroy, INT)
FCALLSCFUN1 (INT, gridDuplicate, GRIDDUPLICATE, gridduplicate, INT)
FCALLSCFUN1 (INT, gridInqType, GRIDINQTYPE, gridinqtype, INT)
FCALLSCFUN1 (INT, gridInqSize, GRIDINQSIZE, gridinqsize, INT)
FCALLSCSUB2 (gridDefXsize, GRIDDEFXSIZE, griddefxsize, INT, INT)
FCALLSCFUN1 (INT, gridInqXsize, GRIDINQXSIZE, gridinqxsize, INT)
FCALLSCSUB2 (gridDefYsize, GRIDDEFYSIZE, griddefysize, INT, INT)
FCALLSCFUN1 (INT, gridInqYsize, GRIDINQYSIZE, gridinqysize, INT)
FCALLSCSUB2 (gridDefNP, GRIDDEFNP, griddefnp, INT, INT)
FCALLSCFUN1 (INT, gridInqNP, GRIDINQNP, gridinqnp, INT)
FCALLSCSUB2 (gridDefXvals, GRIDDEFXVALS, griddefxvals, INT, DOUBLEV)
FCALLSCFUN2 (INT, gridInqXvals, GRIDINQXVALS, gridinqxvals, INT, DOUBLEV)
FCALLSCSUB2 (gridDefYvals, GRIDDEFYVALS, griddefyvals, INT, DOUBLEV)
FCALLSCFUN2 (INT, gridInqYvals, GRIDINQYVALS, gridinqyvals, INT, DOUBLEV)
FCALLSCSUB2 (gridDefXname, GRIDDEFXNAME, griddefxname, INT, STRING)
FCALLSCSUB2 (gridInqXname, GRIDINQXNAME, gridinqxname, INT, PSTRING)
FCALLSCSUB2 (gridDefXlongname, GRIDDEFXLONGNAME, griddefxlongname, INT, STRING)
FCALLSCSUB2 (gridInqXlongname, GRIDINQXLONGNAME, gridinqxlongname, INT, PSTRING)
FCALLSCSUB2 (gridDefXunits, GRIDDEFXUNITS, griddefxunits, INT, STRING)
FCALLSCSUB2 (gridInqXunits, GRIDINQXUNITS, gridinqxunits, INT, PSTRING)
FCALLSCSUB2 (gridDefYname, GRIDDEFYNAME, griddefyname, INT, STRING)
FCALLSCSUB2 (gridInqYname, GRIDINQYNAME, gridinqyname, INT, PSTRING)
FCALLSCSUB2 (gridDefYlongname, GRIDDEFYLONGNAME, griddefylongname, INT, STRING)
FCALLSCSUB2 (gridInqYlongname, GRIDINQYLONGNAME, gridinqylongname, INT, PSTRING)
FCALLSCSUB2 (gridDefYunits, GRIDDEFYUNITS, griddefyunits, INT, STRING)
FCALLSCSUB2 (gridInqYunits, GRIDINQYUNITS, gridinqyunits, INT, PSTRING)
FCALLSCSUB2 (gridInqXstdname, GRIDINQXSTDNAME, gridinqxstdname, INT, PSTRING)
FCALLSCSUB2 (gridInqYstdname, GRIDINQYSTDNAME, gridinqystdname, INT, PSTRING)
FCALLSCSUB2 (gridDefPrec, GRIDDEFPREC, griddefprec, INT, INT)
FCALLSCFUN1 (INT, gridInqPrec, GRIDINQPREC, gridinqprec, INT)
FCALLSCFUN2 (DOUBLE, gridInqXval, GRIDINQXVAL, gridinqxval, INT, INT)
FCALLSCFUN2 (DOUBLE, gridInqYval, GRIDINQYVAL, gridinqyval, INT, INT)
FCALLSCFUN1 (DOUBLE, gridInqXinc, GRIDINQXINC, gridinqxinc, INT)
FCALLSCFUN1 (DOUBLE, gridInqYinc, GRIDINQYINC, gridinqyinc, INT)
FCALLSCFUN1 (INT, gridIsCircular, GRIDISCIRCULAR, gridiscircular, INT)
FCALLSCFUN1 (INT, gridIsRotated, GRIDISROTATED, gridisrotated, INT)
FCALLSCSUB2 (gridDefXpole, GRIDDEFXPOLE, griddefxpole, INT, DOUBLE)
FCALLSCFUN1 (DOUBLE, gridInqXpole, GRIDINQXPOLE, gridinqxpole, INT)
FCALLSCSUB2 (gridDefYpole, GRIDDEFYPOLE, griddefypole, INT, DOUBLE)
FCALLSCFUN1 (DOUBLE, gridInqYpole, GRIDINQYPOLE, gridinqypole, INT)
FCALLSCSUB2 (gridDefAngle, GRIDDEFANGLE, griddefangle, INT, DOUBLE)
FCALLSCFUN1 (DOUBLE, gridInqAngle, GRIDINQANGLE, gridinqangle, INT)
FCALLSCFUN1 (INT, gridInqTrunc, GRIDINQTRUNC, gridinqtrunc, INT)
FCALLSCSUB2 (gridDefTrunc, GRIDDEFTRUNC, griddeftrunc, INT, INT)

/*  Hexagonal GME grid  */

FCALLSCSUB2 (gridDefGMEnd, GRIDDEFGMEND, griddefgmend, INT, INT)
FCALLSCFUN1 (INT, gridInqGMEnd, GRIDINQGMEND, gridinqgmend, INT)
FCALLSCSUB2 (gridDefGMEni, GRIDDEFGMENI, griddefgmeni, INT, INT)
FCALLSCFUN1 (INT, gridInqGMEni, GRIDINQGMENI, gridinqgmeni, INT)
FCALLSCSUB2 (gridDefGMEni2, GRIDDEFGMENI2, griddefgmeni2, INT, INT)
FCALLSCFUN1 (INT, gridInqGMEni2, GRIDINQGMENI2, gridinqgmeni2, INT)
FCALLSCSUB2 (gridDefGMEni3, GRIDDEFGMENI3, griddefgmeni3, INT, INT)
FCALLSCFUN1 (INT, gridInqGMEni3, GRIDINQGMENI3, gridinqgmeni3, INT)

/*  Reference of an unstructured grid  */

FCALLSCSUB2 (gridDefNumber, GRIDDEFNUMBER, griddefnumber, INT, INT)
FCALLSCFUN1 (INT, gridInqNumber, GRIDINQNUMBER, gridinqnumber, INT)
FCALLSCSUB2 (gridDefPosition, GRIDDEFPOSITION, griddefposition, INT, INT)
FCALLSCFUN1 (INT, gridInqPosition, GRIDINQPOSITION, gridinqposition, INT)
FCALLSCSUB2 (gridDefReference, GRIDDEFREFERENCE, griddefreference, INT, STRING)
FCALLSCFUN2 (INT, gridInqReference, GRIDINQREFERENCE, gridinqreference, INT, PSTRING)
FCALLSCSUB2 (gridDefUUID, GRIDDEFUUID, griddefuuid, INT, PVOID)
FCALLSCSUB2 (gridInqUUID, GRIDINQUUID, gridinquuid, INT, PVOID)

/*  Lambert Conformal Conic grid (GRIB version)  */

FCALLSCSUB10 (gridDefLCC, GRIDDEFLCC, griddeflcc, INT, DOUBLE, DOUBLE, DOUBLE, DOUBLE, DOUBLE, DOUBLE, DOUBLE, INT, INT)
FCALLSCSUB10 (gridInqLCC, GRIDINQLCC, gridinqlcc, INT, PDOUBLE, PDOUBLE, PDOUBLE, PDOUBLE, PDOUBLE, PDOUBLE, PDOUBLE, PINT, PINT)

/*  Lambert Conformal Conic 2 grid (PROJ version)  */

FCALLSCSUB6 (gridDefLcc2, GRIDDEFLCC2, griddeflcc2, INT, DOUBLE, DOUBLE, DOUBLE, DOUBLE, DOUBLE)
FCALLSCSUB6 (gridInqLcc2, GRIDINQLCC2, gridinqlcc2, INT, PDOUBLE, PDOUBLE, PDOUBLE, PDOUBLE, PDOUBLE)

/*  Lambert Azimuthal Equal Area grid  */

FCALLSCSUB4 (gridDefLaea, GRIDDEFLAEA, griddeflaea, INT, DOUBLE, DOUBLE, DOUBLE)
FCALLSCSUB4 (gridInqLaea, GRIDINQLAEA, gridinqlaea, INT, PDOUBLE, PDOUBLE, PDOUBLE)
FCALLSCSUB2 (gridDefArea, GRIDDEFAREA, griddefarea, INT, DOUBLEV)
FCALLSCSUB2 (gridInqArea, GRIDINQAREA, gridinqarea, INT, DOUBLEV)
FCALLSCFUN1 (INT, gridHasArea, GRIDHASAREA, gridhasarea, INT)
FCALLSCSUB2 (gridDefNvertex, GRIDDEFNVERTEX, griddefnvertex, INT, INT)
FCALLSCFUN1 (INT, gridInqNvertex, GRIDINQNVERTEX, gridinqnvertex, INT)
FCALLSCSUB2 (gridDefXbounds, GRIDDEFXBOUNDS, griddefxbounds, INT, DOUBLEV)
FCALLSCFUN2 (INT, gridInqXbounds, GRIDINQXBOUNDS, gridinqxbounds, INT, DOUBLEV)
FCALLSCSUB2 (gridDefYbounds, GRIDDEFYBOUNDS, griddefybounds, INT, DOUBLEV)
FCALLSCFUN2 (INT, gridInqYbounds, GRIDINQYBOUNDS, gridinqybounds, INT, DOUBLEV)
FCALLSCSUB3 (gridDefRowlon, GRIDDEFROWLON, griddefrowlon, INT, INT, INTV)
FCALLSCSUB2 (gridInqRowlon, GRIDINQROWLON, gridinqrowlon, INT, INTV)
FCALLSCSUB2 (gridChangeType, GRIDCHANGETYPE, gridchangetype, INT, INT)
FCALLSCSUB2 (gridDefComplexPacking, GRIDDEFCOMPLEXPACKING, griddefcomplexpacking, INT, INT)
FCALLSCFUN1 (INT, gridInqComplexPacking, GRIDINQCOMPLEXPACKING, gridinqcomplexpacking, INT)

/*  ZAXIS routines  */

FCALLSCSUB2 (zaxisName, ZAXISNAME, zaxisname, INT, PSTRING)
FCALLSCFUN2 (INT, zaxisCreate, ZAXISCREATE, zaxiscreate, INT, INT)
FCALLSCSUB1 (zaxisDestroy, ZAXISDESTROY, zaxisdestroy, INT)
FCALLSCFUN1 (INT, zaxisInqType, ZAXISINQTYPE, zaxisinqtype, INT)
FCALLSCFUN1 (INT, zaxisInqSize, ZAXISINQSIZE, zaxisinqsize, INT)
FCALLSCFUN1 (INT, zaxisDuplicate, ZAXISDUPLICATE, zaxisduplicate, INT)
FCALLSCSUB2 (zaxisResize, ZAXISRESIZE, zaxisresize, INT, INT)
FCALLSCSUB2 (zaxisPrint, ZAXISPRINT, zaxisprint, INT, INT)
FCALLSCSUB2 (zaxisDefLevels, ZAXISDEFLEVELS, zaxisdeflevels, INT, DOUBLEV)
FCALLSCSUB2 (zaxisInqLevels, ZAXISINQLEVELS, zaxisinqlevels, INT, DOUBLEV)
FCALLSCSUB3 (zaxisDefLevel, ZAXISDEFLEVEL, zaxisdeflevel, INT, INT, DOUBLE)
FCALLSCFUN2 (DOUBLE, zaxisInqLevel, ZAXISINQLEVEL, zaxisinqlevel, INT, INT)
FCALLSCSUB2 (zaxisDefNlevRef, ZAXISDEFNLEVREF, zaxisdefnlevref, INT, INT)
FCALLSCFUN1 (INT, zaxisInqNlevRef, ZAXISINQNLEVREF, zaxisinqnlevref, INT)
FCALLSCSUB2 (zaxisDefNumber, ZAXISDEFNUMBER, zaxisdefnumber, INT, INT)
FCALLSCFUN1 (INT, zaxisInqNumber, ZAXISINQNUMBER, zaxisinqnumber, INT)
FCALLSCSUB2 (zaxisDefUUID, ZAXISDEFUUID, zaxisdefuuid, INT, PVOID)
FCALLSCSUB2 (zaxisInqUUID, ZAXISINQUUID, zaxisinquuid, INT, PVOID)
FCALLSCSUB2 (zaxisDefName, ZAXISDEFNAME, zaxisdefname, INT, STRING)
FCALLSCSUB2 (zaxisInqName, ZAXISINQNAME, zaxisinqname, INT, PSTRING)
FCALLSCSUB2 (zaxisDefLongname, ZAXISDEFLONGNAME, zaxisdeflongname, INT, STRING)
FCALLSCSUB2 (zaxisInqLongname, ZAXISINQLONGNAME, zaxisinqlongname, INT, PSTRING)
FCALLSCSUB2 (zaxisDefUnits, ZAXISDEFUNITS, zaxisdefunits, INT, STRING)
FCALLSCSUB2 (zaxisInqUnits, ZAXISINQUNITS, zaxisinqunits, INT, PSTRING)
FCALLSCSUB2 (zaxisInqStdname, ZAXISINQSTDNAME, zaxisinqstdname, INT, PSTRING)
FCALLSCSUB2 (zaxisDefPsName, ZAXISDEFPSNAME, zaxisdefpsname, INT, STRING)
FCALLSCSUB2 (zaxisInqPsName, ZAXISINQPSNAME, zaxisinqpsname, INT, PSTRING)
FCALLSCSUB2 (zaxisDefPrec, ZAXISDEFPREC, zaxisdefprec, INT, INT)
FCALLSCFUN1 (INT, zaxisInqPrec, ZAXISINQPREC, zaxisinqprec, INT)
FCALLSCSUB2 (zaxisDefPositive, ZAXISDEFPOSITIVE, zaxisdefpositive, INT, INT)
FCALLSCFUN1 (INT, zaxisInqPositive, ZAXISINQPOSITIVE, zaxisinqpositive, INT)
FCALLSCSUB1 (zaxisDefScalar, ZAXISDEFSCALAR, zaxisdefscalar, INT)
FCALLSCFUN1 (INT, zaxisInqScalar, ZAXISINQSCALAR, zaxisinqscalar, INT)
FCALLSCSUB2 (zaxisDefLtype, ZAXISDEFLTYPE, zaxisdefltype, INT, INT)
FCALLSCFUN1 (INT, zaxisInqLtype, ZAXISINQLTYPE, zaxisinqltype, INT)
FCALLSCSUB3 (zaxisDefVct, ZAXISDEFVCT, zaxisdefvct, INT, INT, DOUBLEV)
FCALLSCSUB2 (zaxisInqVct, ZAXISINQVCT, zaxisinqvct, INT, DOUBLEV)
FCALLSCFUN1 (INT, zaxisInqVctSize, ZAXISINQVCTSIZE, zaxisinqvctsize, INT)
FCALLSCSUB2 (zaxisDefLbounds, ZAXISDEFLBOUNDS, zaxisdeflbounds, INT, DOUBLEV)
FCALLSCFUN2 (INT, zaxisInqLbounds, ZAXISINQLBOUNDS, zaxisinqlbounds, INT, DOUBLEV)
FCALLSCFUN2 (DOUBLE, zaxisInqLbound, ZAXISINQLBOUND, zaxisinqlbound, INT, INT)
FCALLSCSUB2 (zaxisDefUbounds, ZAXISDEFUBOUNDS, zaxisdefubounds, INT, DOUBLEV)
FCALLSCFUN2 (INT, zaxisInqUbounds, ZAXISINQUBOUNDS, zaxisinqubounds, INT, DOUBLEV)
FCALLSCFUN2 (DOUBLE, zaxisInqUbound, ZAXISINQUBOUND, zaxisinqubound, INT, INT)
FCALLSCSUB2 (zaxisDefWeights, ZAXISDEFWEIGHTS, zaxisdefweights, INT, DOUBLEV)
FCALLSCFUN2 (INT, zaxisInqWeights, ZAXISINQWEIGHTS, zaxisinqweights, INT, DOUBLEV)
FCALLSCSUB2 (zaxisChangeType, ZAXISCHANGETYPE, zaxischangetype, INT, INT)

/*  TAXIS routines  */

FCALLSCFUN1 (INT, taxisCreate, TAXISCREATE, taxiscreate, INT)
FCALLSCSUB1 (taxisDestroy, TAXISDESTROY, taxisdestroy, INT)
FCALLSCFUN1 (INT, taxisDuplicate, TAXISDUPLICATE, taxisduplicate, INT)
FCALLSCSUB2 (taxisCopyTimestep, TAXISCOPYTIMESTEP, taxiscopytimestep, INT, INT)
FCALLSCSUB2 (taxisDefType, TAXISDEFTYPE, taxisdeftype, INT, INT)
FCALLSCSUB2 (taxisDefVdate, TAXISDEFVDATE, taxisdefvdate, INT, INT)
FCALLSCSUB2 (taxisDefVtime, TAXISDEFVTIME, taxisdefvtime, INT, INT)
FCALLSCFUN1 (INT, taxisInqVdate, TAXISINQVDATE, taxisinqvdate, INT)
FCALLSCFUN1 (INT, taxisInqVtime, TAXISINQVTIME, taxisinqvtime, INT)
FCALLSCSUB2 (taxisDefRdate, TAXISDEFRDATE, taxisdefrdate, INT, INT)
FCALLSCSUB2 (taxisDefRtime, TAXISDEFRTIME, taxisdefrtime, INT, INT)
FCALLSCFUN1 (INT, taxisInqRdate, TAXISINQRDATE, taxisinqrdate, INT)
FCALLSCFUN1 (INT, taxisInqRtime, TAXISINQRTIME, taxisinqrtime, INT)
FCALLSCSUB2 (taxisDefFdate, TAXISDEFFDATE, taxisdeffdate, INT, INT)
FCALLSCSUB2 (taxisDefFtime, TAXISDEFFTIME, taxisdefftime, INT, INT)
FCALLSCFUN1 (INT, taxisInqFdate, TAXISINQFDATE, taxisinqfdate, INT)
FCALLSCFUN1 (INT, taxisInqFtime, TAXISINQFTIME, taxisinqftime, INT)
FCALLSCFUN1 (INT, taxisHasBounds, TAXISHASBOUNDS, taxishasbounds, INT)
FCALLSCSUB1 (taxisDeleteBounds, TAXISDELETEBOUNDS, taxisdeletebounds, INT)
FCALLSCSUB3 (taxisDefVdateBounds, TAXISDEFVDATEBOUNDS, taxisdefvdatebounds, INT, INT, INT)
FCALLSCSUB3 (taxisDefVtimeBounds, TAXISDEFVTIMEBOUNDS, taxisdefvtimebounds, INT, INT, INT)
FCALLSCSUB3 (taxisInqVdateBounds, TAXISINQVDATEBOUNDS, taxisinqvdatebounds, INT, PINT, PINT)
FCALLSCSUB3 (taxisInqVtimeBounds, TAXISINQVTIMEBOUNDS, taxisinqvtimebounds, INT, PINT, PINT)
FCALLSCSUB2 (taxisDefCalendar, TAXISDEFCALENDAR, taxisdefcalendar, INT, INT)
FCALLSCFUN1 (INT, taxisInqCalendar, TAXISINQCALENDAR, taxisinqcalendar, INT)
FCALLSCSUB2 (taxisDefTunit, TAXISDEFTUNIT, taxisdeftunit, INT, INT)
FCALLSCFUN1 (INT, taxisInqTunit, TAXISINQTUNIT, taxisinqtunit, INT)
FCALLSCSUB2 (taxisDefForecastTunit, TAXISDEFFORECASTTUNIT, taxisdefforecasttunit, INT, INT)
FCALLSCFUN1 (INT, taxisInqForecastTunit, TAXISINQFORECASTTUNIT, taxisinqforecasttunit, INT)
FCALLSCSUB2 (taxisDefForecastPeriod, TAXISDEFFORECASTPERIOD, taxisdefforecastperiod, INT, DOUBLE)
FCALLSCFUN1 (DOUBLE, taxisInqForecastPeriod, TAXISINQFORECASTPERIOD, taxisinqforecastperiod, INT)
FCALLSCSUB2 (taxisDefNumavg, TAXISDEFNUMAVG, taxisdefnumavg, INT, INT)
FCALLSCFUN1 (INT, taxisInqType, TAXISINQTYPE, taxisinqtype, INT)
FCALLSCFUN1 (INT, taxisInqNumavg, TAXISINQNUMAVG, taxisinqnumavg, INT)
FCALLSCFUN1 (STRING, tunitNamePtr, TUNITNAMEPTR, tunitnameptr, INT)

/*  Institut routines  */

FCALLSCFUN4 (INT, institutDef, INSTITUTDEF, institutdef, INT, INT, STRING, STRING)
FCALLSCFUN4 (INT, institutInq, INSTITUTINQ, institutinq, INT, INT, STRING, STRING)
FCALLSCFUN0 (INT, institutInqNumber, INSTITUTINQNUMBER, institutinqnumber)
FCALLSCFUN1 (INT, institutInqCenter, INSTITUTINQCENTER, institutinqcenter, INT)
FCALLSCFUN1 (INT, institutInqSubcenter, INSTITUTINQSUBCENTER, institutinqsubcenter, INT)
FCALLSCFUN1 (STRING, institutInqNamePtr, INSTITUTINQNAMEPTR, institutinqnameptr, INT)
FCALLSCFUN1 (STRING, institutInqLongnamePtr, INSTITUTINQLONGNAMEPTR, institutinqlongnameptr, INT)

/*  Model routines  */

FCALLSCFUN3 (INT, modelDef, MODELDEF, modeldef, INT, INT, STRING)
FCALLSCFUN3 (INT, modelInq, MODELINQ, modelinq, INT, INT, STRING)
FCALLSCFUN1 (INT, modelInqInstitut, MODELINQINSTITUT, modelinqinstitut, INT)
FCALLSCFUN1 (INT, modelInqGribID, MODELINQGRIBID, modelinqgribid, INT)
FCALLSCFUN1 (STRING, modelInqNamePtr, MODELINQNAMEPTR, modelinqnameptr, INT)

/*  Table routines  */

FCALLSCSUB2 (tableWriteC, TABLEWRITEC, tablewritec, STRING, INT)
FCALLSCSUB2 (tableWrite, TABLEWRITE, tablewrite, STRING, INT)
FCALLSCFUN1 (INT, tableRead, TABLEREAD, tableread, STRING)
FCALLSCFUN3 (INT, tableDef, TABLEDEF, tabledef, INT, INT, STRING)
FCALLSCFUN1 (STRING, tableInqNamePtr, TABLEINQNAMEPTR, tableinqnameptr, INT)
FCALLSCSUB5 (tableDefEntry, TABLEDEFENTRY, tabledefentry, INT, INT, STRING, STRING, STRING)
FCALLSCFUN3 (INT, tableInq, TABLEINQ, tableinq, INT, INT, STRING)
FCALLSCFUN0 (INT, tableInqNumber, TABLEINQNUMBER, tableinqnumber)
FCALLSCFUN1 (INT, tableInqNum, TABLEINQNUM, tableinqnum, INT)
FCALLSCFUN1 (INT, tableInqModel, TABLEINQMODEL, tableinqmodel, INT)
FCALLSCSUB5 (tableInqPar, TABLEINQPAR, tableinqpar, INT, INT, PSTRING, PSTRING, PSTRING)
FCALLSCFUN3 (INT, tableInqParCode, TABLEINQPARCODE, tableinqparcode, INT, PSTRING, PINT)
FCALLSCFUN3 (INT, tableInqParName, TABLEINQPARNAME, tableinqparname, INT, INT, PSTRING)
FCALLSCFUN3 (INT, tableInqParLongname, TABLEINQPARLONGNAME, tableinqparlongname, INT, INT, PSTRING)
FCALLSCFUN3 (INT, tableInqParUnits, TABLEINQPARUNITS, tableinqparunits, INT, INT, PSTRING)
FCALLSCFUN2 (STRING, tableInqParNamePtr, TABLEINQPARNAMEPTR, tableinqparnameptr, INT, INT)
FCALLSCFUN2 (STRING, tableInqParLongnamePtr, TABLEINQPARLONGNAMEPTR, tableinqparlongnameptr, INT, INT)
FCALLSCFUN2 (STRING, tableInqParUnitsPtr, TABLEINQPARUNITSPTR, tableinqparunitsptr, INT, INT)

/*  History routines  */

FCALLSCSUB3 (streamDefHistory, STREAMDEFHISTORY, streamdefhistory, INT, INT, STRING)
FCALLSCFUN1 (INT, streamInqHistorySize, STREAMINQHISTORYSIZE, streaminqhistorysize, INT)
FCALLSCSUB2 (streamInqHistoryString, STREAMINQHISTORYSTRING, streaminqhistorystring, INT, PSTRING)

/*  Subtype routines  */

FCALLSCFUN1 (INT, subtypeCreate, SUBTYPECREATE, subtypecreate, INT)

/*  Gives a textual summary of the variable subtype  */

FCALLSCSUB1 (subtypePrint, SUBTYPEPRINT, subtypeprint, INT)

/*  Compares two subtype data structures.  */

FCALLSCFUN2 (INT, subtypeCompare, SUBTYPECOMPARE, subtypecompare, INT, INT)
FCALLSCFUN1 (INT, subtypeInqSize, SUBTYPEINQSIZE, subtypeinqsize, INT)
FCALLSCFUN1 (INT, subtypeInqActiveIndex, SUBTYPEINQACTIVEINDEX, subtypeinqactiveindex, INT)
FCALLSCSUB2 (subtypeDefActiveIndex, SUBTYPEDEFACTIVEINDEX, subtypedefactiveindex, INT, INT)

/*  Generate a "query object" out of a key-value pair.  */


/*  Generate an AND-combined "query object" out of two previous
         query objects.  */

FCALLSCFUN3 (INT, subtypeInqTile, SUBTYPEINQTILE, subtypeinqtile, INT, INT, INT)
FCALLSCFUN2 (INT, vlistInqVarSubtype, VLISTINQVARSUBTYPE, vlistinqvarsubtype, INT, INT)
FCALLSCSUB3 (gribapiLibraryVersion, GRIBAPILIBRARYVERSION, gribapilibraryversion, PINT, PINT, PINT)

#endif
