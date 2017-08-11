# standard libraries
import numpy as np
import logging

# third party libraries
import uuid

# local libraries
from . import sampler

class DoubleGaussianFilterAMOperationDelegate(object):

    def __init__(self, api):
        self.__api = api
        self.panel_id = 'Sampler-Panel'
        self.panel_name = 'Sampler'
        self.panel_positions = ['left', 'right']
        self.panel_position = 'right'
        self.fft_data_item = None
        self.small_ellipse_region = None
        self.big_ellipse_region = None
        self.line_profile_data_item = None
        self.interval_region = None
        parameters = {'base_vec_1': [0,0], 'base_vec_2': [0,0], 'sample_rate': 1, 'offset': [0,0],
                           'average_unitcell_shape': [64,64], 'periodic_repeats': 1}
        self.source_data_item = None
        self.sampler = sampler.Sampler(**parameters)

    def create_panel_widget(self, ui, document_controller):

        def base_vec_1_finished(text, update_sampling=False):
            if len(text) > 0:
                try:
                    b_value = float(base_vec_1_field_b.text)
                    a_value = float(base_vec_1_field_a.text)
                except ValueError:
                    pass
                else:
                    if (np.array([b_value, a_value]) != self.sampler.base_vec_1).any():
                        self.sampler.base_vec_1 = np.array([b_value, a_value])
                        if update_sampling:
                            self.update_calculation()

            base_vec_1_field_b.text = '{:g}'.format(self.sampler.base_vec_1[0])
            base_vec_1_field_a.text = '{:g}'.format(self.sampler.base_vec_1[1])

        def base_vec_2_finished(text, update_sampling=False):
            if len(text) > 0:
                try:
                    b_value = float(base_vec_2_field_b.text)
                    a_value = float(base_vec_2_field_a.text)
                except ValueError:
                    pass
                else:
                    if (np.array([b_value, a_value]) != self.sampler.base_vec_2).any():
                        self.sampler.base_vec_2 = np.array([b_value, a_value])
                        if update_sampling:
                            self.update_calculation()

            base_vec_2_field_b.text = '{:g}'.format(self.sampler.base_vec_2[0])
            base_vec_2_field_a.text = '{:g}'.format(self.sampler.base_vec_2[1])

        def sample_rate_finished(text, update_sampling=False):
            if len(text) > 0:
                try:
                    value = int(text)
                except ValueError:
                    pass
                else:
                    if value != self.sampler.sample_rate:
                        self.sampler.sample_rate = value
                        if update_sampling:
                            self.update_calculation()

            sample_rate_field.text = '{:g}'.format(self.sampler.sample_rate)

        def offset_finished(text, update_sampling=False):
            if len(text) > 0:
                try:
                    y_value = float(offset_field_y.text)
                    x_value = float(offset_field_x.text)
                except ValueError:
                    pass
                else:
                    if (np.array([y_value, x_value]) != self.sampler.offset).any():
                        self.sampler.offset = np.array([y_value, x_value])
                        if update_sampling:
                            self.update_calculation()

            offset_field_y.text = '{:g}'.format(self.sampler.offset[0])
            offset_field_x.text = '{:g}'.format(self.sampler.offset[1])

        def average_shape_finished(text, update_sampling=False):
            if len(text) > 0:
                try:
                    y_value = int(average_unitcell_shape_field_y.text)
                    x_value = int(average_unitcell_shape_field_x.text)
                except ValueError:
                    pass
                else:
                    if (np.array([y_value, x_value]) != self.sampler.average_unitcell_shape).any():
                        self.sampler.average_unitcell_shape = np.array([y_value, x_value])
                        if update_sampling:
                            self.update_calculation()

            average_unitcell_shape_field_y.text = '{:g}'.format(self.sampler.average_unitcell_shape[0])
            average_unitcell_shape_field_x.text = '{:g}'.format(self.sampler.average_unitcell_shape[1])

        def periodic_repeats_finished(text):
            if len(text) > 0:
                try:
                    value = float(text)
                except ValueError:
                    pass
                else:
                    if value != self.sampler.periodic_repeats:
                        self.sampler.periodic_repeats = value

            periodic_repeats_field.text = '{:g}'.format(self.sampler.periodic_repeats)

        def run_button_clicked():
            self.source_data_item = document_controller.target_data_item
            data = self.source_data_item.xdata.data
            if not data.dtype == np.uint16:
                data = sampler.calculate_counts(data)
            self.sampler.image = data
            self.update_calculation()

        def get_start_values_from_fft_button_clicked():
            self.fft_data_item = document_controller.target_data_item
            if self.fft_data_item is not None:
                if self.source_data_item is not None:
                    self.update_metadata(self.source_data_item, 'sampler.fft_uuid', self.fft_data_item.uuid.hex)
                base_vecs = []
                for graphic in self.fft_data_item.graphics:
                    if graphic.type == 'line-region':
                        base_vecs.append(graphic)
                    if len(base_vecs) == 2:
                        break
                if len(base_vecs) < 2:
                    logging.info('You have to add two lines to the FFT in order to read start vectors.')
                    return
                v1 = np.array(base_vecs[0].end)*np.array(self.fft_data_item.xdata.data_shape) - np.array(self.fft_data_item.xdata.data_shape)/2
                v2 = np.array(base_vecs[1].end)*np.array(self.fft_data_item.xdata.data_shape) - np.array(self.fft_data_item.xdata.data_shape)/2
                base = self.sampler.calculate_base_from_fft(v1, v2)
                self.sampler.base_vec_1=base[0]*np.array(self.fft_data_item.xdata.data_shape)
                self.sampler.base_vec_2=base[1]*np.array(self.fft_data_item.xdata.data_shape)
                base_vec_1_finished('')
                base_vec_2_finished('')


        column = ui.create_column_widget()
        base_vec_1_label = ui.create_label_widget('Base Vector 1 ')
        base_vec_1_field_a = ui.create_line_edit_widget()
        base_vec_1_field_b = ui.create_line_edit_widget()
        base_vec_2_label = ui.create_label_widget('Base Vector 2 ')
        base_vec_2_field_a = ui.create_line_edit_widget()
        base_vec_2_field_b = ui.create_line_edit_widget()
        offset_label = ui.create_label_widget('Offset ')
        offset_field_x = ui.create_line_edit_widget()
        offset_field_y = ui.create_line_edit_widget()
        sample_rate_label = ui.create_label_widget('Sample rate ')
        sample_rate_field = ui.create_line_edit_widget()
        average_unitcell_shape_label = ui.create_label_widget('Avg uc shape ')
        average_unitcell_shape_field_x = ui.create_line_edit_widget()
        average_unitcell_shape_field_y = ui.create_line_edit_widget()
        periodic_repeats_label = ui.create_label_widget('Periodically repeat result ')
        periodic_repeats_field = ui.create_line_edit_widget()

        run_button = ui.create_push_button_widget('Run')
        get_start_values_from_fft_button = ui.create_push_button_widget('Get vectors from FFT')

        split_row = ui.create_row_widget()
        base_vec_column = ui.create_column_widget()
        offset_column = ui.create_column_widget()
        base_vec_row_1 = ui.create_row_widget()
        base_vec_row_2 = ui.create_row_widget()
        offset_row = ui.create_row_widget()
        average_unitcell_shape_row = ui.create_row_widget()
        sample_rate_row = ui.create_row_widget()
        run_row = ui.create_row_widget()


        split_row.add_spacing(5)
        split_row.add(base_vec_column)
        split_row.add_spacing(5)
        split_row.add(offset_column)
        split_row.add_stretch()
        split_row.add_spacing(5)

        base_vec_column.add(base_vec_row_1)
        base_vec_column.add(base_vec_row_2)

        offset_column.add(offset_row)

        base_vec_row_1.add(base_vec_1_label)
        base_vec_row_1.add(base_vec_1_field_b)
        base_vec_row_1.add_spacing(5)
        base_vec_row_1.add(base_vec_1_field_a)
        base_vec_row_1.add_stretch()

        base_vec_row_2.add(base_vec_2_label)
        base_vec_row_2.add(base_vec_2_field_b)
        base_vec_row_2.add_spacing(5)
        base_vec_row_2.add(base_vec_2_field_a)
        base_vec_row_2.add_stretch()

        offset_row.add(offset_label)
        offset_row.add(offset_field_y)
        offset_row.add_spacing(5)
        offset_row.add(offset_field_x)
        offset_row.add_stretch()

        average_unitcell_shape_row.add_spacing(5)
        average_unitcell_shape_row.add(average_unitcell_shape_label)
        average_unitcell_shape_row.add(average_unitcell_shape_field_y)
        average_unitcell_shape_row.add_spacing(5)
        average_unitcell_shape_row.add(average_unitcell_shape_field_x)
        average_unitcell_shape_row.add_spacing(5)
        average_unitcell_shape_row.add_stretch()

        sample_rate_row.add_spacing(5)
        sample_rate_row.add(sample_rate_label)
        sample_rate_row.add(sample_rate_field)
        sample_rate_row.add_spacing(10)
        sample_rate_row.add(periodic_repeats_label)
        sample_rate_row.add(periodic_repeats_field)
        sample_rate_row.add_spacing(5)
        sample_rate_row.add_stretch()

        run_row.add_spacing(5)
        run_row.add(run_button)
        run_row.add_stretch()
        run_row.add(get_start_values_from_fft_button)
        run_row.add_spacing(5)

        column.add_spacing(5)
        column.add(split_row)
        column.add(average_unitcell_shape_row)
        column.add(sample_rate_row)
        column.add(run_row)
        column.add_spacing(5)
        column.add_stretch()

        base_vec_1_field_a._widget.on_return_pressed = lambda: base_vec_1_finished(base_vec_1_field_a.text, update_sampling=True)
        base_vec_1_field_b._widget.on_return_pressed = lambda: base_vec_1_finished(base_vec_1_field_b.text, update_sampling=True)
        base_vec_2_field_a._widget.on_return_pressed = lambda: base_vec_2_finished(base_vec_2_field_a.text, update_sampling=True)
        base_vec_2_field_b._widget.on_return_pressed = lambda: base_vec_2_finished(base_vec_2_field_b.text, update_sampling=True)
        offset_field_x._widget.on_return_pressed = lambda: offset_finished(offset_field_x.text, update_sampling=True)
        offset_field_y._widget.on_return_pressed = lambda: offset_finished(offset_field_y.text, update_sampling=True)
        sample_rate_field._widget.on_return_pressed = lambda: sample_rate_finished(sample_rate_field.text, update_sampling=True)
        average_unitcell_shape_field_x._widget.on_return_pressed = lambda: average_shape_finished(average_unitcell_shape_field_x.text, update_sampling=True)
        average_unitcell_shape_field_y._widget.on_return_pressed = lambda: average_shape_finished(average_unitcell_shape_field_y.text, update_sampling=True)

        base_vec_1_field_a.on_editing_finished = base_vec_1_finished
        base_vec_1_field_b.on_editing_finished = base_vec_1_finished
        base_vec_2_field_a.on_editing_finished = base_vec_2_finished
        base_vec_2_field_b.on_editing_finished = base_vec_2_finished
        offset_field_x.on_editing_finished = offset_finished
        offset_field_y.on_editing_finished = offset_finished
        sample_rate_field.on_editing_finished = sample_rate_finished
        average_unitcell_shape_field_x.on_editing_finished = average_shape_finished
        average_unitcell_shape_field_y.on_editing_finished = average_shape_finished
        periodic_repeats_field.on_editing_finished = periodic_repeats_finished

        run_button.on_clicked = run_button_clicked
        get_start_values_from_fft_button.on_clicked = get_start_values_from_fft_button_clicked

        base_vec_1_finished('')
        base_vec_2_finished('')
        offset_finished('')
        sample_rate_finished('')
        average_shape_finished('')
        periodic_repeats_finished('')

        return column


    def update_calculation(self):
        self.sampler.load_c_sampler()
        result_data_item = self.get_result_data_item()
        if result_data_item is not None:
            print(self.sampler.sample_image())
            self.update_metadata(self.source_data_item, 'sampler.result_uuid', result_data_item.uuid.hex)
            result_data_item.set_data(self.sampler.average_unitcell)
            pretty_data_item = self.get_pretty_data_item()
            if pretty_data_item is not None:
                self.update_metadata(self.source_data_item, 'sampler.pretty_uuid', pretty_data_item.uuid.hex)
                self.sampler.make_pretty_output()
                pretty_data_item.set_data(self.sampler.pretty_unitcell)

    def update_metadata(self, data_item, key, value):
        metadata = data_item.metadata
        metadata[key] = value
        data_item.set_metadata(metadata)

    def get_metadata_value(self, data_item, key):
        return data_item.metadata.get(key)

    def update_session_metadata(self, key, value):
        metadata = self.__api.library._document_model.session_metadata.copy()
        metadata[key] = value
        self.__api.library._document_model.session_metadata = metadata

    def get_session_metadata_value(self, key):
        try:
            value = self.__api.library._document_model.session_metadata[key]
        except KeyError:
            return None
        else:
            return value

    def get_fft_data_item(self):
        fft_uuid = self.get_session_metadata_value('sampler.fft_uuid')
        fft_data_item = None
        if fft_uuid is not None:
            fft_data_item = self.__api.library.get_data_item_by_uuid(uuid.UUID(fft_uuid))
        if fft_data_item is None:
            fft_data_item = self.__api.library.create_data_item("Filtered FFT")
            self.update_session_metadata('sampler.fft_uuid', fft_data_item.uuid.hex)
        return fft_data_item

    def get_pretty_data_item(self):
        if self.source_data_item is None:
            return

        pretty_uuid = self.get_metadata_value(self.source_data_item, 'sampler.pretty_uuid')

        if pretty_uuid is not None:
            pretty_data_item = self.__api.library.get_data_item_by_uuid(uuid.UUID(pretty_uuid))
        else:
            pretty_data_item = None

        if pretty_data_item is None:
            pretty_data_item = self.__api.library.create_data_item('Pretty sampling result of ' + self.source_data_item.title)

        return pretty_data_item

    def get_result_data_item(self):
        if self.source_data_item is None:
            return

        result_uuid = self.get_metadata_value(self.source_data_item, 'sampler.result_uuid')

        if result_uuid is not None:
            result_data_item = self.__api.library.get_data_item_by_uuid(uuid.UUID(result_uuid))
        else:
            result_data_item = None

        if result_data_item is None:
            result_data_item = self.__api.library.create_data_item('Sampling result of ' + self.source_data_item.title)

        return result_data_item


class SamplerExtension(object):

    # required for Swift to recognize this as an extension class.
    extension_id = "univie.extensions.sampler"

    def __init__(self, api_broker):
        # grab the api object.
        api = api_broker.get_api(version="1", ui_version="1")
        # be sure to keep a reference or it will be closed immediately.
        self.__panel_ref = api.create_panel(DoubleGaussianFilterAMOperationDelegate(api))

    def close(self):
        # close will be called when the extension is unloaded. in turn, close any references so they get closed. this
        # is not strictly necessary since the references will be deleted naturally when this object is deleted.
        self.__panel_ref.close()
        self.__panel_ref = None
