# import helper
# import bwa
# import bcftools

if __name__ == "__main__":
    import helper
    args = helper.get_variant_args()

    import bwa
    bwa.bwa_pipeline(args)
    helper.run(args)

    import bcftools
    bcftools.bcftools_pipeline(args)
    helper.run(args)
