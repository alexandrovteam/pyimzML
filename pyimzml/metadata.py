from warnings import warn

from pyimzml.ontology.ontology import lookup_and_convert_cv_param, convert_xml_value, convert_term_name

XMLNS_PREFIX = "{http://psi.hupo.org/ms/mzml}"

class _ParseUtils:
    """
    Utility class for common parsing patterns and tracking created param groups so that
    their refs to other param groups can later be linked up.
    """
    def __init__(self):
        self.created_param_groups = []

    def param_group(self, node, **extra_fields):
        pg = ParamGroup(node, **extra_fields)
        self.created_param_groups.append(pg)
        return pg

    def optional_param_group(self, parent_node, xpath, **extra_fields):
        node = parent_node.find(xpath.format(XMLNS_PREFIX))
        return self.param_group(node, **extra_fields) if node is not None else None

    def param_groups_by_id(self, parent_node, xpath):
        return dict(
            (n.attrib.get('id', idx), self.param_group(n))
            for idx, n in enumerate(parent_node.findall(xpath.format(XMLNS_PREFIX)))
        )

    def param_groups_list(self, parent_node, xpath):
        return [self.param_group(n) for n in parent_node.findall(xpath.format(XMLNS_PREFIX))]

    def refs_list(self, parent_node, xpath):
        return [n.attrib.get('ref') for n in parent_node.findall(xpath.format(XMLNS_PREFIX))]


class Metadata:
    def __init__(self, root):
        pu = _ParseUtils()

        fd_node = root.find('{0}fileDescription'.format(XMLNS_PREFIX))
        self.file_description = pu.param_group(
            fd_node.find('{0}fileContent'.format(XMLNS_PREFIX)),
            source_files=pu.param_groups_by_id(fd_node, '{0}sourceFileList/{0}sourceFile'),
            contacts=pu.param_groups_list(fd_node, '{0}contact'),
        )

        self.referenceable_param_groups = pu.param_groups_by_id(
            root,
            '{0}referenceableParamGroupList/{0}referenceableParamGroup'
        )
        self.samples = pu.param_groups_by_id(root, '{0}sampleList/{0}sample')
        self.softwares = pu.param_groups_by_id(root, '{0}softwareList/{0}software')

        self.scan_settings = {}
        for node in root.findall('{0}scanSettingsList/{0}scanSettings'.format(XMLNS_PREFIX)):
            self.scan_settings[node.attrib.get('id')] = pu.param_group(
                node,
                source_file_refs=pu.refs_list(node, '{0}sourceFileRefList/{0}sourceFileRef'),
                targets=pu.param_groups_by_id(node, '{0}targetList/{0}target'),
            )

        self.instrument_configurations = {}
        for node in root.findall('{0}instrumentConfigurationList/{0}instrumentConfiguration'.format(XMLNS_PREFIX)):
            self.instrument_configurations[node.attrib.get('id')] = pu.param_group(
                node,
                components=pu.param_groups_list(node, '{0}componentList/*'),
                software_ref=next(iter(pu.refs_list(node, '{0}softwareRef')), None),
            )

        self.data_processings = {}
        for node in root.findall('{0}dataProcessingList/{0}dataProcessing'.format(XMLNS_PREFIX)):
            self.data_processings[node.attrib.get('id')] = pu.param_group(
                node,
                methods=pu.param_groups_list(node, '{0}processingMethod')
            )

        # Apply referenceable_param_groups
        for pg in pu.created_param_groups:
            pg.apply_referenceable_param_groups(self.referenceable_param_groups)


class ParamGroup:
    def __init__(self, elem, **extra_attrs):
        """
        Parses an XML element representing a group of controlled vocabulary parameters.

        :param elem:              an XML element containing cvParam children
        :param extra_attrs:       extra attributes to assign to the class instance
        """
        self._param_group_refs = [
            ref.get('ref')
            for ref in elem.findall('{0}referenceableParamGroupRef'.format(XMLNS_PREFIX))
        ]
        self.type = elem.tag.replace(XMLNS_PREFIX, '')

        # Tuples of (name, accession, parsed_value, raw_value, unit_name, unit_accession)
        # These are kept in an array as the imzML spec allows multiple uses of accession numbers
        # in the same block
        self.cv_params = []
        for node in elem.findall('{0}cvParam'.format(XMLNS_PREFIX)):
            accession = node.attrib.get('accession')
            raw_value = node.attrib.get('value')
            unit_accession = node.attrib.get('unitAccession')
            name, parsed_value, unit_name = lookup_and_convert_cv_param(accession, raw_value, unit_accession)
            self.cv_params.append(
                (name, accession, parsed_value, raw_value, unit_name, unit_accession)
            )

        # Tuples of (name, type, parsed_value, raw_value, unit_name, unit_accession)
        self.user_params = []
        for node in elem.findall('{0}userParam'.format(XMLNS_PREFIX)):
            name = node.attrib.get('name')
            dtype = node.attrib.get('dtype')
            raw_value = node.attrib.get('value')
            parsed_value = convert_xml_value(raw_value, dtype)
            unit_accession = node.attrib.get('unitAccession')
            unit_name = convert_term_name(unit_accession)
            self.user_params.append(
                (name, dtype, parsed_value, raw_value, unit_name, unit_accession)
            )

        # Mapping of CV param name to parsed value
        self.param_by_name = {}
        self.param_by_name.update((param[0], param[2]) for param in self.user_params)
        self.param_by_name.update((param[0], param[2]) for param in self.cv_params)
        # Mapping of CV param accession to parsed value
        self.param_by_accession = {
            param[1]: param[2] for param in self.cv_params
        }

        self.attrs = elem.attrib

        self.extra_fields = extra_attrs
        for k, v in extra_attrs.items():
            setattr(self, k, v)

    def __getitem__(self, key):
        try:
            return self.param_by_accession[key]
        except KeyError:
            return self.param_by_name[key]

    def __contains__(self, key):
        return key in self.param_by_accession or key in self.param_by_name

    def apply_referenceable_param_groups(self, rpgs):
        for ref in self._param_group_refs[::-1]:
            rpg = rpgs.get(ref)
            if rpg:
                self.cv_params = rpg.cv_params + self.cv_params

                for name, accession, parsed_value, *_ in rpg.cv_params:
                    if name is not None and name != accession:
                        self.param_by_name.setdefault(name, parsed_value)
                    self.param_by_accession.setdefault(accession, parsed_value)
            else:
                warn('ReferenceableParamGroup "%s" not found' % ref)

    def pretty(self):
        """
        Flattens attributes, params and extra fields into a single dict keyed by name.
        This function is intended to help human inspection. For programmatic access to specific fields,
        always use the `attrs`, `param_by_name`, `param_by_accession`, etc. instance attributes instead.

        :return:
        """
        def deep_pretty(obj):
            if isinstance(obj, ParamGroup):
                return obj.pretty()
            if isinstance(obj, list):
                return [deep_pretty(item) for item in obj]
            if isinstance(obj, dict):
                return {k: deep_pretty(v) for k, v in obj.items()}
            return obj

        result = {
            'type': self.type,
        }
        result.update(self.attrs)
        result.update(self.param_by_name)
        result.update(deep_pretty(self.extra_fields))

        return result


class SpectrumData(ParamGroup):
    def __init__(self, root, referenceable_param_groups):
        pu = _ParseUtils()

        scan_list_params = pu.optional_param_group(root, '{0}scanList')
        scans = []
        for node in root.findall('{0}scanList/{0}scan'.format(XMLNS_PREFIX)):
            scans.append(
                pu.param_group(
                    node,
                    scan_windows=pu.param_groups_list(node, '{0}scanWindowList/{0}scanWindow')
                )
            )

        precursors = []
        for node in root.findall('{0}precursorList/{0}precursor'.format(XMLNS_PREFIX)):
            precursors.append(
                pu.param_group(
                    node,
                    isolation_window=pu.optional_param_group(node, '{0}isolationWindow'),
                    selected_ions=pu.param_groups_list(node, '{0}selectedIonList/{0}selectedIon'),
                    activation=pu.optional_param_group(node, '{0}activation'),
                )
            )

        products = []
        for node in root.findall('{0}productList/{0}product'.format(XMLNS_PREFIX)):
            products.append(
                pu.param_group(
                    node,
                    isolation_window=pu.optional_param_group(node, '{0}isolationWindow'),
                )
            )

        binary_data_arrays = pu.param_groups_list(root, '{0}binaryDataArrayList/{0}binaryDataArray')

        super().__init__(
            root,
            scan_list_params=scan_list_params,
            scans=scans,
            precursors=precursors,
            products=products,
            binary_data_arrays=binary_data_arrays,
        )

        for pg in pu.created_param_groups:
            pg.apply_referenceable_param_groups(referenceable_param_groups)

        self.apply_referenceable_param_groups(referenceable_param_groups)
