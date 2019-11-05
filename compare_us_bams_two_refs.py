import sys
import os
from optparse import OptionParser
from collections import defaultdict

# import tempfile

import pysam

def processSingly( r1, r2, flagRequire, minq, minqdiff, minnmdiff, nmtag, nmsign=1 ):
    ismapped = (r1.flag&0x4==0, r2.flag&0x4==0)

    if flagRequire!=0:
        f1reqPass,f2reqPass=(r1.flag&flagRequire)!=0,(r2.flag&flagRequire)!=0
    else:
        f1reqPass,f2reqPass=True,True

    if ismapped[0] and ismapped[1] and (f1reqPass or f2reqPass):
        mq1='mq1pass' if r1.mapq>=minq else 'mq1fail'
        mq2='mq2pass' if r2.mapq>=minq else 'mq2fail'
        
        if r1.mapq-r2.mapq >= minqdiff:
            md = 'mq1diffoverThresh'
        elif r2.mapq-r1.mapq >= minqdiff:
            md = 'mq2diffoverThresh'
        else:
            md = 'mqWithinThresh'

        nm1 = dict(r1.tags)[nmtag]
        nm2 = dict(r2.tags)[nmtag]

        # NM - lower is better
        if (nm1 - nm2)*nmsign >= minnmdiff:
            nmd = 'nm1diffoverThresh'
        elif (nm2 - nm1)*nmsign >= minnmdiff:
            nmd = 'nm2diffoverThresh'
        else:
            nmd = 'nmWithinThresh'

        s = ('isMappedBoth', mq1, mq2, md, nmd)

    elif ismapped[0] and not ismapped[1]:
        if r1.mapq>=minq and f1reqPass:
            s = ('isMapped1', 'mq1pass', 'NA', 'NA', 'NA')
        else:
            s = ('isMapped1', 'mq1fail', 'NA', 'NA', 'NA')
    elif ismapped[1] and not ismapped[0]:
        if r2.mapq>=minq and f2reqPass:
            s = ('isMapped2', 'NA', 'mq2pass', 'NA', 'NA')
        else:
            s = ('isMapped2', 'NA', 'mq2fail', 'NA', 'NA')
    else:
        s = ('isMappedNone','NA','NA','NA', 'NA')

    return s


if __name__=='__main__':
    opts=OptionParser()
    opts.add_option('','--inBam1',dest='inBam1')
    opts.add_option('','--inBam2',dest='inBam2')

    opts.add_option('','--outBam1',default=None,dest='outBam1')
    opts.add_option('','--outBam2',default=None,dest='outBam2')

    opts.add_option('','--refname1',default=None,dest='refname1')
    opts.add_option('','--refname2',default=None,dest='refname2')

    opts.add_option('','--flagRequireEither',default=0,type=int,dest='flagRequireEither')

    opts.add_option('','--mapqMinThresh',default=20,type=int,dest='mapqMinThresh')
    opts.add_option('','--mapqDiffThresh',default=20,type=int,dest='mapqDiffThresh')
    opts.add_option('','--minNumMismatchDiffThresh',default=1,type=int,dest='minNumMismatchDiffThresh')

    opts.add_option('','--useASinsteadofNM',default=False,action='store_true',dest='useASinsteadofNM')

    opts.add_option('','--outSummary',dest='outSummary')

    opts.add_option('','--kvPassthru',default=None,dest='kvPassthru')

    opts.add_option('','--analyzeMates',default=False,action='store_true',dest='analyzeMates')

    opts.add_option('','--outMatedOnlySummary',default=None,dest='outMatedOnlySummary')
    opts.add_option('','--outSingletonOnlySummary',default=None,dest='outSingletonOnlySummary')

    (o,args)=opts.parse_args()

    kvPassthru = [] if o.kvPassthru is None else [ (kv.split(':')[0],kv.split(':')[1]) for kv in o.kvPassthru.split(',')]
 
    bi1=pysam.Samfile(o.inBam1,'rb')
    bi2=pysam.Samfile(o.inBam2,'rb')

    if o.outBam1 is not None:
        assert o.outBam1 is not None and o.outBam2 is not None
        bo1=pysam.Samfile(o.outBam1,'wb',template=bi1)
        bo2=pysam.Samfile(o.outBam2,'wb',template=bi2)
    else:
        bo1=bo2=None

    # mapped: isMappedBoth, isMapped1, isMapped2, isMappedNone
    # mq >= minthresh: mqAboveMin1, mqAboveMin2
    # mqdifference: mq1diffOverThresh, mq2diffOverThresh, mqWithinThresh
    
    nmtag = 'AS' if o.useASinsteadofNM else 'NM'
    nmsign = 1 if o.useASinsteadofNM else -1

    if not o.analyzeMates:
        mStatusCount = defaultdict( int )

        mapqMinThresh, mapqDiffThresh, minNumMismatchDiffThresh = o.mapqMinThresh, o.mapqDiffThresh, o.minNumMismatchDiffThresh
        flagRequireEither = o.flagRequireEither

        ctr = 0

        if bo1 is None:
            for r1 in bi1:
                r2=bi2.next()
                assert r1.qname == r2.qname, ctr

                ctr+=1
                if ctr % 50000 == 0 :
                    print '%d..'%ctr

                s = processSingly(r1, r2, flagRequireEither, mapqMinThresh, mapqDiffThresh, minNumMismatchDiffThresh, nmtag, nmsign)
                mStatusCount[s]+=1
        else:
            for r1 in bi1:
                r2=bi2.next()
                assert r1.qname == r2.qname, ctr

                ctr+=1
                if ctr % 50000 == 0 :
                    print '%d..'%ctr

                s = processSingly(r1, r2, flagRequireEither, mapqMinThresh, mapqDiffThresh, minNumMismatchDiffThresh, nmtag, nmsign)
                mStatusCount[s]+=1

                r1.tags+=[('Z1','_'.join(s))]
                r2.tags+=[('Z1','_'.join(s))]
                bo1.write(r1)
                bo2.write(r2)

    else:

        mSingleonlyStatusCount = defaultdict( int )
        mPairStatusCount = defaultdict( int )

        mapqMinThresh, mapqDiffThresh, minNumMismatchDiffThresh = o.mapqMinThresh, o.mapqDiffThresh, o.minNumMismatchDiffThresh
        flagRequireEither = o.flagRequireEither

        r1prev,r2prev=None,None
        sprev=None

        ctr = 0

        if bo1 is None:
            for r1 in bi1:
                r2=bi2.next()
                assert r1.qname == r2.qname

                ctr+=1
                if ctr % 50000 == 0 :
                    print '%d..'%ctr

                s = processSingly(r1, r2, flagRequireEither, mapqMinThresh, mapqDiffThresh, minNumMismatchDiffThresh, nmtag, nmsign)

                if r1prev is not None:
                    if r1.qname == r1prev.qname:
                        # prev read was paired.  we are looking at the second mate.
                        mPairStatusCount[ (sprev,s) ] += 1
                        sprev,r1prev,r2prev=None,None,None
                    else:
                        # prev read was singleton
                        mSingleonlyStatusCount[ sprev ] += 1
                        sprev,r1prev,r2prev = s,r1,r2                
                else:
                    sprev,r1prev,r2prev = s,r1,r2                

            if r1prev is not None:
                mSingleonlyStatusCount[ sprev ] += 1

        else:

            for r1 in bi1:
                r2=bi2.next()
                assert r1.qname == r2.qname

                ctr+=1
                if ctr % 50000 == 0 :
                    print '%d..'%ctr
                    # break

                s = processSingly(r1, r2, flagRequireEither, mapqMinThresh, mapqDiffThresh, minNumMismatchDiffThresh, nmtag, nmsign)

                if r1prev is not None:
                    if r1.qname == r1prev.qname:
                        # prev read was paired.  we are looking at the second mate.
                        mPairStatusCount[ (s,sprev) ] += 1

                        r1.tags+=[('Z1','_'.join(s))]
                        r2.tags+=[('Z1','_'.join(s))]
                        r1.tags+=[('Z2','_'.join(sprev))]
                        r2.tags+=[('Z2','_'.join(sprev))]

                        r1prev.tags+=[('Z1','_'.join(sprev))]
                        r2prev.tags+=[('Z1','_'.join(sprev))]
                        r1prev.tags+=[('Z2','_'.join(s))]
                        r2prev.tags+=[('Z2','_'.join(s))]

                        bo1.write(r1prev)
                        bo2.write(r2prev)
                        bo1.write(r1)
                        bo2.write(r2)

                        sprev,r1prev,r2prev=None,None,None
                    else:
                        # prev read was singleton
                        mSingleonlyStatusCount[ sprev ] += 1

                        r1prev.tags+=[('Z1','_'.join(sprev))]
                        r2prev.tags+=[('Z1','_'.join(sprev))]

                        bo1.write(r1prev)
                        bo2.write(r2prev)

                        sprev,r1prev,r2prev = s,r1,r2                
                else:
                    sprev,r1prev,r2prev = s,r1,r2                

            if r1prev is not None:
                mSingleonlyStatusCount[ sprev ] += 1

                r1prev.tags+=[('Z1','_'.join(sprev))]
                r2prev.tags+=[('Z1','_'.join(sprev))]

                bo1.write(r1prev)
                bo2.write(r2prev)

            bo1.close()
            bo2.close()

        # sum the mated regions and the singletons to get the marginal counts for each status
        mStatusCount = defaultdict(int)
        for k in mSingleonlyStatusCount:
            mStatusCount[k]+=mSingleonlyStatusCount[k]
        for kk in mPairStatusCount:
            mStatusCount[kk[0]]+=mPairStatusCount[kk]
            mStatusCount[kk[1]]+=mPairStatusCount[kk]

    filOut = open(o.outSummary,'w')

    hdrOut = [ 'numReads', 'mapStatus', 'mq1status', 'mq2status', 'mqdiffStatus', 'nmdiffStatus' ] + [ kv[0] for kv in kvPassthru ]
    lv = [ v for v in hdrOut ]
    for iv in xrange(len(lv)):
        if o.refname1 is not None:
            lv[iv] = lv[iv].replace('1',o.refname1)
        if o.refname2 is not None:
            lv[iv] = lv[iv].replace('2',o.refname2)
    filOut.write( tabnl(hdrOut) )

    for k in sorted( mStatusCount.keys(), key=lambda j:-mStatusCount[j] ):
        lv = [ v for v in k ]
        for iv in xrange(len(lv)):
            if o.refname1 is not None:
                lv[iv] = lv[iv].replace('1',o.refname1)
            if o.refname2 is not None:
                lv[iv] = lv[iv].replace('2',o.refname2)

        lout=[ mStatusCount[k] ]+lv+[ kv[1] for kv in kvPassthru ]
        filOut.write(tabnl(lout))

    filOut.close()

    if o.analyzeMates:
        if o.outSingletonOnlySummary is not None:
            filOut = open(o.outSingletonOnlySummary,'w')

            hdrOut = [ 'numReads', 'mapStatus', 'mq1status', 'mq2status', 'mqdiffStatus', 'nmdiffStatus' ] + [ kv[0] for kv in kvPassthru ]
            lv = [ v for v in hdrOut ]
            for iv in xrange(len(lv)):
                if o.refname1 is not None:
                    lv[iv] = lv[iv].replace('1',o.refname1)
                if o.refname2 is not None:
                    lv[iv] = lv[iv].replace('2',o.refname2)
            filOut.write( tabnl(hdrOut) )

            for k in sorted( mSingleonlyStatusCount.keys(), key=lambda j:-mSingleonlyStatusCount[j] ):
                lv = [ v for v in k ]
                for iv in xrange(len(lv)):
                    if o.refname1 is not None:
                        lv[iv] = lv[iv].replace('1',o.refname1)
                    if o.refname2 is not None:
                        lv[iv] = lv[iv].replace('2',o.refname2)

                lout=[ mSingleonlyStatusCount[k] ]+lv+[ kv[1] for kv in kvPassthru ]
                filOut.write(tabnl(lout))

            filOut.close()

        if o.outMatedOnlySummary is not None:
            filOut = open(o.outMatedOnlySummary,'w')

            hdrOut = [ 'numReads', 
                       'firstRd_mapStatus','secondRd_mapStatus',
                       'firstRd_mq1status', 'secondRd_mq1status', 
                       'firstRd_mq2status', 'secondRd_mq2status', 
                       'firstRd_mqdiffStatus', 'secondRd_mqdiffStatus', 
                       'firstRd_nmdiffStatus', 'secondRd_nmdiffStatus' ] + [ kv[0] for kv in kvPassthru ]
            lv = [ v for v in hdrOut ]
            for iv in xrange(len(lv)):
                if o.refname1 is not None:
                    lv[iv] = lv[iv].replace('1',o.refname1)
                if o.refname2 is not None:
                    lv[iv] = lv[iv].replace('2',o.refname2)
            filOut.write( tabnl(hdrOut) )

            for kk in sorted( mPairStatusCount.keys(), key=lambda j:-mPairStatusCount[j] ):
                # lv = [ v for v in k ]
                lv = [ kk[0][0], kk[1][0], kk[0][1], kk[1][1], kk[0][2],kk[1][2], kk[0][3],kk[1][3], kk[0][4],kk[1][4] ]
                for iv in xrange(len(lv)):
                    if o.refname1 is not None:
                        lv[iv] = lv[iv].replace('1',o.refname1)
                    if o.refname2 is not None:
                        lv[iv] = lv[iv].replace('2',o.refname2)

                lout=[ mPairStatusCount[kk] ]+lv+[ kv[1] for kv in kvPassthru ]
                filOut.write(tabnl(lout))

            filOut.close()
                


