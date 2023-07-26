# import helper
# import bwa
# import bcftools

if __name__ == "__main__":
    import helper
    args = helper.get_variant_args()

    import bwa1
    bwa1.helper.run(args)

    import bcftools1
    bcftools1.helper.run(args)
