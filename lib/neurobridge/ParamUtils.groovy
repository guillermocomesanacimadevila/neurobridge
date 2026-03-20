class ParamUtils {
    static void requireParam(value, String name) {
        if( !value ) {
            throw new IllegalArgumentException("Missing --${name}")
        }
    }
}