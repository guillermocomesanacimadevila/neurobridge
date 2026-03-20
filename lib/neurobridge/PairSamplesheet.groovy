class PairSamplesheet {

    static Map buildMeta(row) {
        def meta = [:]
        meta.trait1     = row.trait1.toString().trim()
        meta.trait2     = row.trait2.toString().trim()
        meta.cases1     = row.cases1
        meta.controls1  = row.controls1
        meta.cases2     = row.cases2
        meta.controls2  = row.controls2
        meta.pop_prev1  = row.pop_prev1
        meta.pop_prev2  = row.pop_prev2
        return meta
    }
}