from modules.logger_factory import LoggerFactory

from .re_functions import REbase_org, insert_site_CDS, sites_in_cds, remove_site_from_plasmid

logger = LoggerFactory.create_logger("RE")


def parse_inp1(usr_inp1):
    'separate the input into 2 lists- list of optimized organism names'
    'and deoptimized organism names'
    optimized_org_names = []
    deoptimized_org_names = []
    for org_name, org_dict in usr_inp1['organisms'].items():
        try:
            if org_dict['optimized']:
                optimized_org_names.append(org_name)
            else:
                deoptimized_org_names.append(org_name)
        except:
            continue
    return optimized_org_names, deoptimized_org_names


def multi_organism_RE_dict(org_list, cds_aa):
    RE_dict = {}
    for org in org_list:
        RE_dict[org] = REbase_org(org, cds_aa)
        logger.info(f'Number of {org} restriction enzymes found in sequence: {len(RE_dict[org])}')
    return RE_dict



def multi_org_remove_site(optimized_RE_dict, cds_nt):
    RE_dict = {}
    for org_RE_dict in optimized_RE_dict.values():
        RE_dict = {**RE_dict, **org_RE_dict}

    try:
        cds_nt = remove_site_from_plasmid(RE_dict, cds_nt)
    except:
        logger.info('Not all restriction requirements are applicable')
    return cds_nt


def multi_org_final_found_sites(RE_dict, final_cds_nt):
    for org, org_enzyme_dict in RE_dict.items():
        found_sites_dict = sites_in_cds(org_enzyme_dict, final_cds_nt)
        logger.info(f'For {org}, {len(found_sites_dict)} sites were found in the final coding sequence,'
                    f'belonging to the following enzymes')
        for enzyme_name in found_sites_dict.keys():
            logger.info(f'{enzyme_name} : {found_sites_dict[enzyme_name]}')
            # logger.info(found_sites_dict)
#            ambiguous_site = found_sites_dict[enzyme_name]['ambiguous_site']
#            logger.info(f'{enzyme_name}, which has the ambiguous site {ambiguous_site}')


def total_sequence_analysis(optimized_RE_dict, deoptimized_RE_dict, final_cds_nt):
    logger.info('Sites from optimized organisms: ')
    multi_org_final_found_sites(optimized_RE_dict, final_cds_nt)
    logger.info('Sites from deoptimized organisms: ')
    multi_org_final_found_sites(deoptimized_RE_dict, final_cds_nt)
