#!/usr/bin/env python

""" MultiQC module to parse output from somalier """

from __future__ import print_function
from collections import OrderedDict
import logging
import csv

import random
from collections import defaultdict
from math import isnan, isinf
from multiqc import config
from multiqc.plots import bargraph
from multiqc.plots import heatmap
from multiqc.plots import scatter
from multiqc.modules.base_module import BaseMultiqcModule

# Initialise the logger
log = logging.getLogger(__name__)


class MultiqcModule(BaseMultiqcModule):
    """
    somalier module class, parses stderr logs.
    """

    def __init__(self):

        # Initialise the parent object
        super(MultiqcModule, self).__init__(name='somalier', anchor='somalier',
        href='https://github.com/brentp/somalier',
        info="calculates genotype :: pedigree correspondence checks from sketches derived from BAM/CRAM or VCF")

        # Find and load any somalier reports
        self.somalier_data = dict()
        self.somalier_ancestry_cats = list()
        self.somalier_length_counts = dict()
        self.somalier_length_exp = dict()
        self.somalier_length_obsexp = dict()

        # parse somalier summary file
        for f in self.find_log_files('somalier/samples'):
            parsed_data = self.parse_somalier_summary(f)
            if parsed_data is not None:
                for s_name in parsed_data:
                    s_name = self.clean_s_name(s_name, f['root'])
                    if s_name in self.somalier_data.keys():
                            log.debug("Duplicate sample name found! Overwriting: {}".format(s_name))
                    self.add_data_source(f, s_name)
                    self.somalier_data[s_name] = parsed_data[s_name]
        # parse somalier CSV files # TODO: CSV or TSV?
        for pattern in ['pairs']:
            sp_key = 'somalier/{}'.format(pattern)
            for f in self.find_log_files(sp_key):
                # some columns have the same name in het_check and sex_check (median_depth)
                # pass pattern to parse_somalier_pairs_csv so the column names can include pattern to
                # avoid being overwritten
                parsed_data = self.parse_somalier_pairs_csv(f, pattern)
                if parsed_data is not None:
                    for s_name in parsed_data:
                        s_name = self.clean_s_name(s_name, f['root'])
                        if s_name in self.somalier_data.keys():
                            log.debug("Duplicate sample name found! Overwriting: {}".format(s_name))
                        self.add_data_source(f, s_name)
                        self.somalier_data[s_name] = parsed_data[s_name]

        # parse somalier ancestry files
        for f in self.find_log_files('somalier/ancestry_prediction', filehandles=True):
            self.parse_somalier_ancestry_prediction(f)
            
        # parse somalier pc files
        for f in self.find_log_files('somalier/ancestry_pcs', filehandles=True):
            self.parse_somalier_ancestry_pcs(f)

        # Filter to strip out ignored sample names
        self.somalier_data = self.ignore_samples(self.somalier_data)

        if len(self.somalier_data) == 0:
            raise UserWarning

        log.info("Found {} reports".format(len(self.somalier_data)))

        # Write parsed report data to a file
        self.write_data_file(self.somalier_data, 'multiqc_somalier')

        # Basic Stats Table
        self.somalier_general_stats_table()

        # Relatedness plot
        self.somalier_relatedness_plot()

        self.somalier_relatedness_heatmap_plot()

        # hetcheck plot
        self.somalier_het_check_plot()

        self.somalier_sex_check_plot()

        # ancestry plots
        self.somalier_ancestry_barplot()

        self.somalier_ancestry_pca_plot()

    def parse_somalier_summary(self, f):
        """ Go through log file looking for somalier output """
        parsed_data = dict()
        headers = None
        sample_i = -100
        for l in f['f'].splitlines():
            s = l.split("\t")
            if headers is None:
                s[0] = s[0].lstrip('#')
                headers = s
                sample_i = headers.index("sample")
            else:
                parsed_data[s[sample_i]] = dict()
                for i, v in enumerate(s):
                    if i != sample_i:
                        try:
                            parsed_data[s[sample_i]][headers[i]] = float(v)
                        except ValueError:
                            parsed_data[s[sample_i]][headers[i]] = v
        if len(parsed_data) == 0:
            return None
        return parsed_data

    def parse_somalier_pairs_csv(self, f, pattern):
        """ Parse csv output from somalier """
        parsed_data = dict()
        headers = None
        s_name_idx = None
        for l in f['f'].splitlines():
            s = l.lstrip('#').split("\t")
            if headers is None:
                headers = s
                try:
                    s_name_idx = [headers.index("sample")]
                except ValueError:
                    try:
                        s_name_idx = [headers.index("sample_a"), headers.index("sample_b")]
                    except ValueError:
                        log.warn("Could not find sample name in somalier output: {}".format(f['fn']))
                        return None
            else:
                s_name = '*'.join([s[idx] for idx in s_name_idx]) # TODO: not safe, but works
                parsed_data[s_name] = dict()
                for i, v in enumerate(s):
                    if i not in s_name_idx: #skip i=0,1, i.e. sample_a, sample_b
                        if (isnan(float(v)) or isinf(float(v))):
                            log.info("Found Inf or NaN value. Overwriting with -2.") # TODO: find better solution
                            v = -2
                        try:
                            # add the pattern as a suffix to key
                            parsed_data[s_name][headers[i] + "_" + pattern] = float(v)
                        except ValueError:
                            # add the pattern as a suffix to key
                            parsed_data[s_name][headers[i] + "_" + pattern] = v
                            

        if len(parsed_data) == 0:
            return None
        return parsed_data

    def parse_somalier_ancestry_prediction(self, f):
        parsed_data = dict()
        reader = csv.DictReader(f['f'])
        idx = "#sample"

        if (reader.fieldnames is not None):
            for c in reader.fieldnames:
                if c != idx:
                    self.somalier_ancestry_cats.append(c)
            for row in reader:
                parsed_data[row[idx]] = {k:v for k,v in row.items() if k != idx}

            if len(parsed_data) > 0:
                for s_name in parsed_data:
                    s_name = self.clean_s_name(s_name, f['root'])
                    if s_name in self.somalier_data.keys():
                        intersect_keys = parsed_data[s_name].keys() & self.somalier_data.keys()
                        if len(intersect_keys) > 0:
                            log.debug("Duplicate sample name found! Overwriting: {} : {}".format(s_name, intersect_keys))
                    self.add_data_source(f, s_name)
                    try:
                        self.somalier_data[s_name].update(parsed_data[s_name])
                    except KeyError:
                        self.somalier_data[s_name] = parsed_data[s_name]
        else:
            log.warn("Detected empty file: {}".format(f['fn']))

    def parse_somalier_ancestry_pcs(self, f):
        bg_pc1 = []
        bg_pc2 = []
        bg_ancestry = []
        parsed_data = dict()

        reader = csv.DictReader(f['f'])
        idx = "ancestry"
        
        if (reader.fieldnames is not None):
            for row in reader:
                if (row[idx] in self.somalier_ancestry_cats):
                    bg_pc1.append(float(row["PC1"]))
                    bg_pc2.append(float(row["PC2"]))
                    bg_ancestry.append(row["ancestry"])
                else:
                    d = {k:float(v) for k,v in row.items() if k != idx}
                    parsed_data[row[idx]] = d
            
            if len(parsed_data) > 0:
                self.somalier_data["background_pcs"] = {"PC1":bg_pc1, "PC2":bg_pc2, "ancestry":bg_ancestry}

                for s_name in parsed_data:
                    s_name = self.clean_s_name(s_name, f['root'])
                    if s_name in self.somalier_data.keys():
                        intersect_keys = parsed_data[s_name].keys() & self.somalier_data.keys()
                        if len(intersect_keys) > 0:
                            log.debug("Duplicate sample name found! Overwriting: {} : {}".format(s_name, intersect_keys))
                    self.add_data_source(f, s_name)
                    try:
                        self.somalier_data[s_name].update(parsed_data[s_name])
                    except KeyError:
                        self.somalier_data[s_name] = parsed_data[s_name]
        else:
            log.warn("Detected empty file: {}".format(f['fn']))

    def somalier_general_stats_table(self):
        """ Take the parsed stats from the somalier report and add it to the
        basic stats table at the top of the report """

        headers = OrderedDict()
        headers['X_het'] = {
            'title': 'heterozygous variants on X chromosome',
        }
        headers['X_hom_alt'] = {
            'title': 'homozygous alternate variants on X chromosome',
        }
        self.general_stats_addcols(self.somalier_data, headers)

    def somalier_relatedness_plot(self):
        data = dict()
        for s_name, d in self.somalier_data.items():
            if 'ibs0_pairs' in d and 'ibs2_pairs' in d:
                data[s_name] = {
                    'x': d['ibs0_pairs'],
                    'y': d['ibs2_pairs']
                }
            if 'relatedness_pairs' in d:
                if d['expected_relatedness_pairs'] < 0.125:
                    data[s_name]['color'] = 'rgba(109, 164, 202, 0.9)'
                elif d['expected_relatedness_pairs'] < 0.5:
                    data[s_name]['color'] = 'rgba(250, 160, 81, 0.8)'
                else:
                    data[s_name]['color'] = 'rgba(43, 159, 43, 0.8)'

        pconfig = {
            'id': 'somalier_relatedness_plot',
            'title': 'somalier: Relatedness Plot',
            'xlab': 'IBS0 (no alleles shared)',
            'ylab': 'IBS2 (both alleles shared)',
        }

        if len(data) > 0:
            self.add_section (
                name = 'Relatedness',
                anchor = 'somalier-relatedness-plot',
                description = """Shared allele rates between sample pairs.  Points are coloured by degree of expected-relatedness:
                <span style="color: #6DA4CA;">less than 0.25</span>,
                <span style="color: #FAA051;">0.25 - 0.5</span>,
                <span style="color: #2B9F2B;">greather than 0.5</span>.""",
                plot = scatter.plot(data, pconfig)
            )

    def somalier_relatedness_heatmap_plot(self):
        # inspiration: MultiQC/modules/vcftools/relatedness2.py
        
        data = []
        labels = set()
        rels = defaultdict(dict)
        for s_name, d in self.somalier_data.items():
            if 'relatedness_pairs' in d:
                a, b = s_name.split('*')
                labels.add(a)
                labels.add(b)
                rels[a][b] = float(d["relatedness_pairs"])

        # impose alphabetical order and avoid json serialisation errors in utils.report
        labels = sorted(labels)
        
        for x in labels:
            line = []
            for y in labels:
                try: # 
                    line.append(rels[x][y])
                except KeyError:
                    try:
                        line.append(rels[y][x])
                    except KeyError:
                        try:
                            line.append(rels[x][x])
                        except KeyError:
                            line.append(-2)
            data.append(line)

        pconfig = {
            'id': 'somalier_relatedness_heatmap_plot',
            'title': 'somalier: Relatedness Heatmap Plot',
            'xlab': 'Sample A',
            'ylab': 'Sample B',
        }

        if len(data) > 0:
            self.add_section (
                name = 'Relatedness Heatmap',
                anchor = 'somalier-relatedness-heatmap-plot',
                description = """Heatmap displaying relatedness of sample pairs.""",
                # plot = scatter.plot(data, pconfig)
                plot = heatmap.plot(
                    data = data,
                    xcats = labels,
                    ycats = labels,
                    pconfig = pconfig,
                )
            )
    
    def somalier_het_check_plot(self):
        """plot the het_check scatter plot"""
        # empty dictionary to add sample names, and dictionary of values
        data = {}

        # for each sample, and list in self.somalier_data
        for s_name, d in self.somalier_data.items():
            # check the sample contains the required columns
            if 'gt_depth_mean' in d and 'ab_std' in d:
                # add sample to dictionary with value as a dictionary of points to plot
                data[s_name] = {
                    'x': d['gt_depth_mean'],
                    'y': d['ab_std']
                }

        pconfig = {
            'id': 'somalier_het_check_plot',
            'title': 'somalier: Het Check',
            'xlab': 'mean depth',
            'ylab': 'standard deviation of allele-balance',
        }
        
        if len(data) > 0:
            self.add_section (
                name = 'Het Check',
                description = "Std devation of heterozygous allele balance against mean depth.",
                helptext = """A high standard deviation in allele balance suggests contamination.
                """,
                anchor = 'somalier-hetcheck-plot',
                plot = scatter.plot(data, pconfig)
            )

    def somalier_sex_check_plot(self):
        data = {}
        sex_index = {"female": 0, "male": 1, "unknown": 2}

        for s_name, d in self.somalier_data.items():
            if 'X_depth_mean' in d and 'pedigree_sex' in d:
                data[s_name] = {
                    'x': (random.random() - 0.5) * 0.1 + sex_index.get(d['pedigree_sex'], 2),
                    'y': d["X_depth_mean"]
                }

        pconfig = {
            'id': 'somalier_sex_check_plot',
            'title': 'somalier: Sex Check',
            'xlab': 'Sex From Ped',
            'ylab': 'Scaled mean depth on X',
            'categories': ["Female", "Male", "Unknown"]
        }
        if len(data) > 0:
            self.add_section(
                name = 'Sex Check',
                description = "Predicted sex against scaled depth on X",
                helptext = """
                Higher values of depth, low values suggest male.
                """,
                anchor='somalier-sexcheck-plot',
                plot=scatter.plot(data, pconfig)
            )

    def somalier_ancestry_barplot(self):
        data = dict()
        
        for s_name, d in self.somalier_data.items():
            # ensure that only relevant items are added, 
            # i.e. only ancestry category values
            data[s_name] = {k:v for k,v in d.items() if k in self.somalier_ancestry_cats}

        pconfig = {
            'id' : 'somalier_ancestry_barplot',
            'cpswitch_c_active' : False,
            'title' : 'somalier: Ancestry Prediction, Barplot',
            'ylab' : 'Predicted Ancestry'
        }

        if len(data) > 0:
            self.add_section(
                name = "Ancestry Barplot",
                description = "Predicted ancestry of samples.",
                helptext = """Only non-zero prediction categories are shown, e.g.
                if SAS ancestry prediction is 0 for all samples, the category is 
                not shown""",
                anchor = "somalier-ancestry-barplot",
                plot = bargraph.plot(data=data, pconfig=pconfig)
            )

    def somalier_ancestry_pca_plot(self):
        ancestry_colors = {
            'SAS': 'rgb(68,1,81,1)',
            'EAS': 'rgb(59,81,139,1)',
            'AMR': 'rgb(33,144,141,1)',
            'AFR': 'rgb(92,200,99,1)',
            'EUR': 'rgb(253,231,37,1)'
        }
        background_ancestry_colors = {
            'SAS': 'rgb(68,1,81,0.1)',
            'EAS': 'rgb(59,81,139,0.1)',
            'AMR': 'rgb(33,144,141,0.1)',
            'AFR': 'rgb(92,200,99,0.1)',
            'EUR': 'rgb(253,231,37,0.1)'
        }
        default_color = '#000000'
        default_background_color = 'rgb(211,211,211,0.05)'
        data = OrderedDict()

        d = self.somalier_data.pop("background_pcs", {})
        if d:
            background = [{'x': pc1,
                        'y': pc2,
                        'color': default_background_color,
                        'name': ancestry,
                        'marker_size': 1}
                        for pc1, pc2, ancestry in zip(d['PC1'], d['PC2'], d['ancestry'])]
            data["background"] = background
        
        for s_name, d in self.somalier_data.items():
            if 'PC1' in d and 'PC2' in d:
                data[s_name] = {
                    'x': d['PC1'],
                    'y': d['PC2'],
                    # 'color': ancestry_colors.get(d['ancestry-prediction'], default_color)
                    'color' : '#000000'
                }
        pconfig = {
            'id' : 'somalier_ancestry_pca_plot',
            'title' : 'somalier: Ancestry Prediction, PCA',
            'xlab': 'PC1',
            'ylab': 'PC2',
            'marker_size': 5,
            'marker_line_width': 0
        }

        if len(data) > 0:
            self.add_section(
                name = "Ancestry PCA Plot",
                description = "TBA",
                helptext = """TBA""",
                anchor = "somalier-ancestry-pca-plot",
                plot = scatter.plot(data, pconfig)
            )


