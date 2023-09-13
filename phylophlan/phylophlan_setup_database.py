#!/usr/bin/env python


__author__ = ('Francesco Asnicar (f.asnicar@unitn.it), '
              'Francesco Beghini (francesco.beghini@unitn.it), '
              'Claudia Mengoni (claudia.mengoni@studenti.unitn.it), '
              'Mattia Bolzan (mattia.bolzan@unitn.it), '
              'Nicola Segata (nicola.segata@unitn.it)')
__version__ = '3.0.25'
__date__ = '27 March 2023'


import sys
import bz2
import glob
from Bio import SeqIO
from Bio.SeqRecord import SeqRecord
import os
import argparse as ap
import re
import json
import zlib
from xml.etree import ElementTree
from urllib.request import urlretrieve
from urllib.parse import urlencode
from urllib.parse import urlparse
from urllib.parse import parse_qs
from urllib.request import Request
from urllib.request import urlopen
import requests
from requests.adapters import HTTPAdapter
from requests.adapters import Retry
import time


if sys.version_info[0] < 3:
    raise Exception("PhyloPhlAn {} requires Python 3, your current Python version is {}.{}.{}"
                    .format(__version__, sys.version_info[0], sys.version_info[1], sys.version_info[2]))

DB_TYPE_CHOICES = ['n', 'a']
GENOME_EXTENSION = '.fna'
PROTEOME_EXTENSION = '.faa'
#DOWNLOAD_URL = "https://www.dropbox.com/s/jrituki21vx6p30/taxa2core.txt?dl=1"
DOWNLOAD_URL = "http://cmprod1.cibio.unitn.it/databases/PhyloPhlAn/taxa2core.txt"


def info(s, init_new_line=False, exit=False, exit_value=0):
    if init_new_line:
        sys.stdout.write('\n')

    sys.stdout.write('{}'.format(s))
    sys.stdout.flush()

    if exit:
        sys.exit(exit_value)


def error(s, init_new_line=False, exit=False, exit_value=1):
    if init_new_line:
        sys.stderr.write('\n')

    sys.stderr.write('[e] {}\n'.format(s))
    sys.stderr.flush()

    if exit:
        sys.exit(exit_value)


def read_params():
    p = ap.ArgumentParser(description=("The phylophlan_setup_database.py script can be used to either format an input folder or "
                                       "multi-fasta file to be used as database in phylophlan.py, or automatically download a "
                                       "pre-identified set of core UniRef90 proteins for the taxonomic label of a given species"),
                          formatter_class=ap.ArgumentDefaultsHelpFormatter)

    group = p.add_mutually_exclusive_group()
    group.add_argument('-i', '--input', type=str,
                       help=("Specify the path to either the folder containing the marker files or the file of markers, in "
                             "(multi-)fasta format"))
    group.add_argument('-g', '--get_core_proteins', type=str, default=None,
                       help=('Specify the taxonomic label for which download the set of core proteins. The label must represent a '
                             'species: "--get_core_proteins s__Escherichia_coli"'))

    p.add_argument('--database_update', action='store_true', default=False, help="Update the databases file")
    p.add_argument('-o', '--output', type=str, default=None, help="Specify path to the output folder where to save the database")
    p.add_argument('-d', '--db_name', type=str, help="Specify the name of the output database")
    p.add_argument('-e', '--input_extension', type=str, default=None,
                   help="Specify the extension of the input file(s) specified via -i/--input")
    p.add_argument('-t', '--db_type', default=None, choices=DB_TYPE_CHOICES,
                   help='Specify the type of the database, where "n" stands for nucleotides and "a" for amino acids')
    p.add_argument('-x', '--output_extension', type=str, default=None, help="Set the database output extension")
    p.add_argument('--overwrite', action='store_true', default=False, help="If specified and the output file exists it will be overwritten")
    p.add_argument('--citation', action='version',
                   version=('Asnicar, F., Thomas, A.M., Beghini, F. et al. '
                            'Precise phylogenetic analysis of microbial isolates and genomes from metagenomes using PhyloPhlAn 3.0. '
                            'Nat Commun 11, 2500 (2020). '
                            'https://doi.org/10.1038/s41467-020-16366-7'),
                   help="Show citation")
    p.add_argument('--verbose', action='store_true', default=False, help="Prints more stuff")
    p.add_argument('-v', '--version', action='version',
                   version='phylophlan_setup_database.py version {} ({})'.format(__version__, __date__),
                   help="Prints the current phylophlan_setup_database.py version and exit")
    return p.parse_args()


def check_params(args, verbose=False):
    if args.database_update:
        database_update(update=args.database_update, verbose=args.verbose)
        args.database_update = False

    if not args.input and not args.get_core_proteins:
        error('either -i/--input or -g/--get_core_proteins must be specified', init_new_line=True, exit=True)

    if args.input:
        if (not os.path.isdir(args.input)) and (not os.path.isfile(args.input)):
            error('input must be either a folder or a file', exit=True)

        if not args.db_name:
            error('-d/--db_name must be specified', exit=True)

        if os.path.isdir(args.input):
            if not args.input_extension:
                error('input is a folder, hence --input_extension must be specified', exit=True)

    if args.get_core_proteins:
        # not handling keys different than species (Tue 19 Jun 2018)
        if args.get_core_proteins[:3] != 's__':
            error('The taxonomic label provided "{}" does not starts with "s__"'.format(args.get_core_proteins),
                  exit=True)

        if not args.output:
            output_tmp = os.path.dirname(args.get_core_proteins)

            if output_tmp:
                args.output = output_tmp
            else:
                args.output = args.get_core_proteins

            if verbose:
                info('Setting output folder "{}"\n'.format(args.output))
        elif os.path.isdir(args.output) and (not args.output.endswith(args.get_core_proteins)):
            args.output = os.path.join(args.output, args.get_core_proteins)

            if verbose:
                info('Setting output folder "{}"\n'.format(args.output))

        args.input = args.output
        args.input_extension = PROTEOME_EXTENSION
        args.output_extension = PROTEOME_EXTENSION

        if not args.db_name:
            args.db_name = args.get_core_proteins

            if verbose:
                info('Setting db_name to "{}"\n'.format(args.db_name))

    if not args.db_type:
        if not args.output_extension:
            error('either -t/--db_type or -x/--output_extension must be specified', exit=True)
        else:
            if not args.output_extension.startswith('.'):
                args.output_extension = '.' + args.output_extension

            if args.output_extension.endswith('.'):
                args.output_extension = args.output_extension[:-1]
    else:
        if not args.output_extension:
            if 'n' == args.db_type:
                args.output_extension = GENOME_EXTENSION
            if 'a' == args.db_type:
                args.output_extension = PROTEOME_EXTENSION

            if verbose:
                info('Setting output extension to "{}"\n'.format(args.output_extension))
        else:
            error("both -t/--db_type and -x/--output_extension were specified, don't know which one to use!",
                  exit=True)

    if not args.output:
        args.output = args.input if os.path.isdir(args.input) else os.path.dirname(args.input)

        if verbose:
            info('Output folder not specified, setting to "{}"\n'.format(args.output))

    if verbose:
        info('Arguments: {}\n'.format(vars(args)))


def check_response(response):
    try:
        response.raise_for_status()
    except requests.HTTPError:
        error(response.json(), init_new_line=True, exit=True)


def submit_id_mapping(api_URL, IDs_list, from_db, to_db):
    # valid from_db to_db pairs at https://rest.uniprot.org/configure/idmapping/fields
    request = requests.post(f"{api_URL}/idmapping/run", data={"from": from_db, "to": to_db, "ids": ",".join(IDs_list)})
    check_response(request)
    return request.json()["jobId"]


def get_next_link(headers):
    re_next_link = re.compile(r'<(.+)>; rel="next"')

    if "Link" in headers:
        match = re_next_link.match(headers["Link"])

        if match:
            return match.group(1)


def check_id_mapping_results_ready(session, api_URL, job_ID, polling_interval, verbose=False):
    while True:
        request = session.get(f"{api_URL}/idmapping/status/{job_ID}")
        check_response(request)
        j = request.json()

        if "jobStatus" in j:
            if j["jobStatus"] == "RUNNING":
                if verbose:
                    info(f"Check Retrying in {polling_interval}s\n")

                time.sleep(polling_interval)
            else:
                error("check_id_mapping_results_ready(): " + request["jobStatus"], init_new_line=True, exit=True)
        else:
            return bool(j["results"] or j["failedIds"])


def get_batch(session, batch_response, file_format, compressed):
    batch_url = get_next_link(batch_response.headers)

    while batch_url:
        batch_response = session.get(batch_url)
        batch_response.raise_for_status()
        yield decode_results(batch_response, file_format, compressed)
        batch_url = get_next_link(batch_response.headers)


def combine_batches(all_results, batch_results, file_format):
    if file_format == "json":
        for key in ("results", "failedIds"):
            if key in batch_results and batch_results[key]:
                all_results[key] += batch_results[key]
    elif file_format == "tsv":
        return all_results + batch_results[1:]
    else:
        return all_results + batch_results
    return all_results


def get_id_mapping_results_link(session, api_URL, job_ID):
    url = f"{api_URL}/idmapping/details/{job_ID}"
    request = session.get(url)
    check_response(request)
    return request.json()["redirectURL"]


def decode_results(response, file_format, compressed):
    if compressed:
        decompressed = zlib.decompress(response.content, 16 + zlib.MAX_WBITS)

        if file_format == "json":
            j = json.loads(decompressed.decode("utf-8"))
            return j
        elif file_format == "tsv":
            return [line for line in decompressed.decode("utf-8").split("\n") if line]
        elif file_format == "xlsx":
            return [decompressed]
        elif file_format == "xml":
            return [decompressed.decode("utf-8")]
        else:
            return decompressed.decode("utf-8")
    elif file_format == "json":
        return response.json()
    elif file_format == "tsv":
        return [line for line in response.text.split("\n") if line]
    elif file_format == "xlsx":
        return [response.content]
    elif file_format == "xml":
        return [response.text]

    return response.text


def get_xml_namespace(element):
    m = re.match(r"\{(.*)\}", element.tag)
    return m.groups()[0] if m else ""


def merge_xml_results(xml_results):
    merged_root = ElementTree.fromstring(xml_results[0])

    for result in xml_results[1:]:
        root = ElementTree.fromstring(result)

        for child in root.findall("{http://uniprot.org/uniprot}entry"):
            merged_root.insert(-1, child)

    ElementTree.register_namespace("", get_xml_namespace(merged_root[0]))
    return ElementTree.tostring(merged_root, encoding="utf-8", xml_declaration=True)


def get_id_mapping_results_search(session, url):
    parsed = urlparse(url)
    query = parse_qs(parsed.query)
    file_format = query["format"][0] if "format" in query else "json"
    size = int(query["size"][0]) if "size" in query else 500
    query["size"] = size
    compressed = query["compressed"][0].lower() == "true" if "compressed" in query else False
    parsed = parsed._replace(query=urlencode(query, doseq=True))
    url = parsed.geturl()
    request = session.get(url)
    check_response(request)
    results = decode_results(request, file_format, compressed)

    for batch in get_batch(session, request, file_format, compressed):
        results = combine_batches(results, batch, file_format)

    if file_format == "xml":
        return merge_xml_results(results)

    return results


def database_update(update=False, verbose=False):
    taxa2core_file_latest = None
    taxa2core_file = os.path.basename(DOWNLOAD_URL).replace('?dl=1', '')
    download(DOWNLOAD_URL, taxa2core_file, overwrite=update, verbose=verbose)

    with open(taxa2core_file) as f:
        for r in f:
            if not r.startswith('#'):
                taxa2core_file_latest, taxa2core_file_latest_url = r.strip().split('\t')
                break  # file should contains only one line, i.e., the name of the latest taxa2core file

    taxa2core_file_latest = taxa2core_file_latest.replace('?dl=1', '')
    download(taxa2core_file_latest_url, taxa2core_file_latest, overwrite=update, verbose=verbose)

    return taxa2core_file_latest


def byte_to_megabyte(byte):
    """
    Convert byte value to megabyte
    """

    return (byte / 1048576)


class ReportHook():

    def __init__(self):
        self.start_time = time.time()

    def report(self, blocknum, block_size, total_size):
        """
        Print download progress message
        """

        if blocknum == 0:
            self.start_time = time.time()

            if total_size > 0:
                info("Downloading file of size: {:.2f} MB\n".format(byte_to_megabyte(total_size)))
        else:
            total_downloaded = blocknum * block_size
            status = "{:3.2f} MB ".format(byte_to_megabyte(total_downloaded))

            if total_size > 0:
                percent_downloaded = total_downloaded * 100.0 / total_size
                # use carriage return plus sys.stderr to overwrite stderr
                download_rate = total_downloaded / (time.time() - self.start_time)
                estimated_time = (total_size - total_downloaded) / download_rate
                estimated_minutes = int(estimated_time / 60.0)
                estimated_seconds = estimated_time - estimated_minutes * 60.0
                status += ("{:3.2f} %  {:5.2f} MB/sec {:2.0f} min {:2.0f} sec "
                           .format(percent_downloaded, byte_to_megabyte(download_rate),
                                   estimated_minutes, estimated_seconds))

            status += "        \r"
            info(status)


def download(url, download_file, overwrite=False, verbose=False):
    """
    Download a file from a url
    """

    if overwrite or not os.path.isfile(download_file):
        try:
            if verbose:
                info('Downloading "{}" to "{}"\n'.format(url, download_file))

            urlretrieve(url, download_file, reporthook=ReportHook().report)
            info('\n')
        except EnvironmentError:
            error('unable to download "{}"'.format(url))
    elif verbose:
        info('File "{}" present\n'.format(download_file))


def create_folder(output, verbose=False):
    if not os.path.exists(output):
        if verbose:
            info('Creating output folder "{}"\n'.format(output))

        os.mkdir(output, mode=0o775)
    elif verbose:
        info('Output "{}" folder present\n'.format(output))


def get_core_proteins(taxa2core_file, taxa_label, output, output_extension, verbose=False):
    core_proteins = {}
    url = None
    retry2download = []
    not_mapped = []
    not_mapped_again = []

    # taxa2core format
    #
    #   # taxid   taxonomy      UniRefURL       core_list
    #   12345     k__|...|s__   http://..{}..   UP123;UP456;...

    for r in bz2.open(taxa2core_file, 'rt'):
        if r.startswith('#'):
            continue

        r_clean = r.strip().split('\t')

        if taxa_label in r_clean[1]:
            if taxa_label == r_clean[1].split('|')[-1]:
                url = r_clean[2]
                core_proteins[r_clean[1]] = r_clean[3].split(';')
                break
            else:
                core_proteins[r_clean[1]] = None

    if not len(core_proteins):
        error('no entry found for "{}", please check the taxonomic label provided'.format(taxa_label),
              exit=True)
    elif len([k for k, v in core_proteins.items() if v is not None]) > 1:
        error('{} entries found for "{}":\n{}    please check the taxonomic label provided'
              .format(len(core_proteins), taxa_label, '    - {}\n'.join(core_proteins.keys())), exit=True)

    for lbl, core_prots in core_proteins.items():
        if core_prots is None:
            continue

        if verbose:
            info('Downloading {} core proteins for {}\n'.format(len(core_prots), lbl))

        for core_prot in core_prots:
            local_prot = os.path.join(output, core_prot + output_extension)
            download(url.format(core_prot), local_prot, verbose=verbose)

            if not os.path.exists(local_prot):
                retry2download.append(core_prot)

    # try to re-map the ids in case "not mapped" store in not_mapped
    if retry2download:
        resolved_uniref90s = resolve_IDs(retry2download, verbose=verbose)

        for uniref90_id in (x.replace('UniRef90_', '') for x in resolved_uniref90s.values()):
            local_prot = os.path.join(output, uniref90_id + output_extension)
            download(url.format(uniref90_id), local_prot, verbose=verbose)

            if not os.path.exists(local_prot):
                not_mapped.append(uniref90_id)

        if len(resolved_uniref90s) != len(retry2download):
            # probably deleted proteins in the Uniprot versions, try to download their latest version in any case
            for ur90 in set(resolved_uniref90s.keys()) - set(retry2download):
                local_prot = os.path.join(output, resolved_uniref90s[ur90] + output_extension)
                download('https://www.uniprot.org/uniprot/{}.fasta?version=*'.format(ur90), local_prot, verbose=verbose)

                if not os.path.exists(local_prot):
                    not_mapped_again.append(core_prot)

        # really don't know what else to try... I'm sorry!
        if not_mapped_again:
            nd_out = os.path.join(output, taxa_label + '_core_proteins_not_mapped.txt')

            if verbose:
                info('There are {} core proteins that could not be downloaded, writing thier IDs to "{}"\n'
                     .format(len(not_mapped_again), nd_out))

            with open(nd_out, 'w') as f:
                f.write('\n'.join(not_mapped_again) + '\n')


def resolve_IDs(IDs_list, verbose=False):
    if verbose:
        info("Resolving {} UniRef90 ID, please wait this might take some time\n".format(len(IDs_list)))

    polling_interval = 5  # in seconds
    api_URL = "https://rest.uniprot.org"
    retries = Retry(total=5, backoff_factor=0.25, status_forcelist=[500, 502, 503, 504])
    session = requests.Session()
    session.mount("https://", HTTPAdapter(max_retries=retries))
    job_ID = submit_id_mapping(api_URL, IDs_list, from_db="UniProtKB_AC-ID", to_db="UniRef90")

    if check_id_mapping_results_ready(session, api_URL, job_ID, polling_interval, verbose=verbose):
        link = get_id_mapping_results_link(session, api_URL, job_ID)
        results = get_id_mapping_results_search(session, link)

    return dict([(i['from'], i['to']['id']) for i in results['results']])


def create_database(db_name, inputt, input_ext, output, overwrite, verbose=False):
    seqs = []

    if os.path.exists(output) and (not overwrite):
        error('output file exists and --overwrite not specified', exit=True)

    if os.path.isdir(inputt):
        for marker in glob.iglob(os.path.join(inputt, '*' + input_ext + '*')):
            seqs += [SeqRecord(record.seq,
                               id='_'.join([db_name.replace('_', '-').replace(',', '-').replace(':', '').replace('|', '-'),
                                            record.id.replace('_', '-').replace(',', '-').replace(':', '').replace('|', '-'),
                                            str(count)]),
                               description='')
                     for count, record in enumerate(SeqIO.parse(bz2.open(marker, 'rt') if marker.endswith('.bz2')
                                                                else open(marker), "fasta"))]
    else:
        seqs = [SeqRecord(record.seq,
                          id='_'.join([db_name.replace('_', '-').replace(',', '-').replace(':', '').replace('|', '-'),
                                       record.id.replace('_', '-').replace(',', '-').replace(':', '').split('|')[0].replace('|', '-'),
                                       str(count)]),
                          description='')
                for count, record in enumerate(SeqIO.parse(bz2.open(inputt, 'rt') if inputt.endswith('.bz2')
                                                           else open(inputt), "fasta"))]

    if not seqs:
        error('no sequences found, make sure the input folder/file provided is not empty', exit=True)

    if verbose:
        info('Writing output database "{}"\n'.format(output))

    with open(output, 'w') as f:
        SeqIO.write(seqs, f, "fasta")


def phylophlan_setup_database():
    args = read_params()

    if args.verbose:
        info('phylophlan_setup_database.py version {} ({})\n'.format(__version__, __date__))
        info('Command line: {}\n\n'.format(' '.join(sys.argv)), init_new_line=True)

    check_params(args, verbose=args.verbose)
    create_folder(args.output, verbose=args.verbose)

    if args.get_core_proteins:
        taxa2core_file_latest = database_update(update=args.database_update, verbose=args.verbose)
        get_core_proteins(taxa2core_file_latest, args.get_core_proteins, args.output, args.output_extension, verbose=args.verbose)

    create_database(args.db_name, args.input, args.input_extension, os.path.join(args.output, args.db_name + args.output_extension),
                    args.overwrite, verbose=args.verbose)


if __name__ == '__main__':
    t0 = time.time()
    phylophlan_setup_database()
    t1 = time.time()
    info('Total elapsed time {}s\n'.format(int(t1 - t0)), init_new_line=True)
    sys.exit(0)
