import h5py as h5
import numpy as np
import matplotlib.pyplot as plt
from obspy.core import UTCDateTime as udt

# TEST:

def load_data(filename):
    data = {}
    with h5.File(filename, mode='r') as f:
        data['metadata'] = {}
        data['metadata']['stations'] = f['stations'][()].astype('U').tolist()
        data['metadata']['components'] = f['components'][()].astype('U').tolist()
        data['metadata']['date'] = udt(f['date'][()])
        data['metadata']['sampling_rate'] = f['sampling_rate'][()]
        data['waveforms'] = f['waveforms'][()]
    return data

def load_template(filename):
    template = {}
    with h5.File(filename, mode='r') as f:
        for key in f.keys():
            template[key] = f[key][()]
    return template

def load_cc(filename):
    with h5.File(filename, mode='r') as f:
        cc_sum = f['cc_sum'][()]
    return cc_sum

def load_detections(filename, tid,):
    meta = filename + 'meta.h5'
    wav = filename + 'wav.h5'
    detections = {}
    with h5.File(meta, mode='r') as f:
        detections['metadata'] = {}
        for key in f[str(tid)].keys():
            detections['metadata'][key] = f[str(tid)][key][()]
    with h5.File(wav, mode='r') as f:
        detections['waveforms'] = f[str(tid)]['waveforms'][()]
    detections['metadata']['stations'] = detections['metadata']['stations'].astype('U')
    detections['metadata']['components'] = detections['metadata']['components'].astype('U')
    return detections

def extract_new_detections(data,
                           cc_sums,
                           moveout_array,
                           n_mad=10.,
                           template_duration=8.,
                           step=1,
                           extracted_duration=60.):

    n_templates = cc_sums.shape[0] #! number of templates
    n_stations = moveout_array.shape[1] #! number of stations
    n_components = moveout_array.shape[2] #! number of components
    n_extracted_samples = np.int32(extracted_duration * data['metadata']['sampling_rate']) #! number of extracted samples
    buffer_extracted_events = 10. 

    list_metadata = []
    list_waveforms = []
    for i in range(n_templates):
        cc_sum = cc_sums[i, :]

        #hugo adds this to deal with zero values
        nonzero_val = np.where(cc_sum != 0.)
        nonzero_val = list(nonzero_val)[0]
        cc_sum -= np.median(cc_sum[nonzero_val])
        threshold = n_mad * np.median(np.abs(cc_sum[nonzero_val]))
        ##cc_sum -= np.median(cc_sum))
        ##threshold = n_mad * np.median(np.abs(cc_sum))
        # ------------------
        cc_idx = np.argwhere(cc_sum > threshold)
        detections = cc_idx * step

        # only keep highest correlation coefficient for grouped detections
        # we assume the LAST COMPONENT IS THE VERTICAL COMPONENT
        d_mv = moveout_array[i, :, 0] - moveout_array[i, :, -1] #! vector with length of stations
        # fix the maximum window size to 3 times the template duration
        # fix the minimum window size to 1 time the templare duration
        # in between: choose an adaptive size based on the median
        # P-S time
        search_win = min(np.int32(3. * template_duration * 
                                  data['metadata']['sampling_rate'] / step),
                         max(np.int32(1. * np.median(d_mv[d_mv != 0]) / step),
                             np.int32(template_duration *
                                      data['metadata']['sampling_rate'] / step))) #! scalar at the end
        
        for j in range(cc_idx.size):
            idx = np.arange(max(0, cc_idx[j] - search_win // 2),
                            min(cc_sum.size-1, cc_idx[j] + search_win // 2),
                            dtype=np.int32)
            idx_to_update = np.where(cc_idx == cc_idx[j])[0]
            cc_idx[idx_to_update] = np.argmax(cc_sum[idx]) + idx[0]

        cc_idx = np.unique(cc_idx)
        detections = cc_idx * step

        # after this step, we can have detections closest than search_win / 2
        cc_idx = list(cc_idx)
        n_removed = 0
        for j in range(1, detections.size):
            if (cc_idx[j-n_removed] - cc_idx[j-n_removed-1]) < search_win // 2:
                if cc_sum[cc_idx[j-n_removed]] > cc_sum[cc_idx[j-n_removed-1]]:
                    cc_idx.remove(cc_idx[j-n_removed-1])
                else:
                    cc_idx.remove(cc_idx[j-n_removed])
                n_removed += 1
        cc_idx = np.asarray(cc_idx)
        detections = cc_idx * step

        n_multiplets = len(detections)
        # ------------------------------------------------------
        metadata_events = {}
        waveforms_events = {}
        origin_times = np.zeros(n_multiplets, dtype=np.float64)
        correlation_coefficients = np.zeros(n_multiplets, dtype=np.float32)
        waveforms = np.zeros((n_multiplets, n_stations,
                              n_components, n_extracted_samples), dtype=np.float32)
        idx_min = 0  # can't extract continuous data before index 0
        idx_max = data['waveforms'].shape[-1]  # can't extract continuous data after
        #                                        the last sample of the day
        for d in range(n_multiplets):
            origin_time = udt(data['metadata']['date']) \
                          + detections[d] / data['metadata']['sampling_rate']
            origin_times[d] = origin_time.timestamp \
                - buffer_extracted_events
            correlation_coefficients[d] = cc_sum[cc_idx[d]]
            #print(origin_time) # that was testing with HUGO
            # -----------------------------------------
            # take care of not selecting out-of-bound indexes:
            id1 = detections[d] - np.int32(buffer_extracted_events
                                           * data['metadata']['sampling_rate'])
            if id1 < idx_min:
                # will have to zero-pad the beginning of the extracted sequence
                dn_b = idx_min - id1
                id2 = np.int32(id1 + n_extracted_samples)
                id1 = np.int32(idx_min)
            else:
                dn_b = 0
                id2 = id1 + n_extracted_samples
            if id2 > idx_max:
                # will have to zero-pad the end of the extracted sequence
                dn_e = id2 - idx_max
                id2 = np.int32(idx_max)
            else:
                dn_e = 0
            #print(id1, id2, d) # by David
            waveforms[d, :, :, :] = np.concatenate((np.zeros((n_stations,
                                                              n_components,
                                                              dn_b),
                                                             dtype=np.float32),
                                                    data['waveforms'][:,
                                                                     :,
                                                                      id1:id2],
                                                    np.zeros((n_stations,
                                                              n_components,
                                                              dn_e),
                                                             dtype=np.float32)),
                                                   axis=-1)
            # -----------------------------------------
        metadata_events.update({'template_id'                :   np.array([i])})
        metadata_events.update({'stations'                   :   np.asarray(data['metadata']['stations']).astype('S')})
        metadata_events.update({'components'                 :   np.asarray(data['metadata']['components']).astype('S')})
        metadata_events.update({'origin_times'               :   origin_times})
        metadata_events.update({'correlation_coefficients'   :   correlation_coefficients})
        waveforms_events.update({'waveforms'                 :   waveforms})

        list_metadata.append(metadata_events)
        list_waveforms.append(waveforms_events)
    return list_metadata, list_waveforms


def write_new_detections(filename, metadata, waveforms, db_path=''):
    filename_meta = db_path + filename + 'meta.h5'
    filename_wave = db_path + filename + 'wav.h5'
    n_templates = len(metadata)
    with h5.File(filename_meta, mode='w') as f:
        for t in range(n_templates):
            if len(metadata[t]['origin_times']) == 0:
                # no detection
                continue
            f.create_group('{:d}'.format(metadata[t]['template_id'][0]))
            for key in metadata[t].keys():
                f['{:d}'.format(metadata[t]['template_id'][0])].create_dataset(key, data=metadata[t][key], compression='gzip')
    with h5.File(filename_wave, mode='w') as f:
        for t in range(n_templates):
            if len(metadata[t]['origin_times']) == 0:
                # no detection
                continue
            f.create_group('{:d}'.format(metadata[t]['template_id'][0]))
            f['{:d}'.format(metadata[t]['template_id'][0])].create_dataset('waveforms', data=waveforms[t]['waveforms'], compression='lzf')
            print('{:d} events detected with Template {:d}'.format(waveforms[t]['waveforms'].shape[0], metadata[t]['template_id'][0]))


def plot_n_detections(detections, n_best, template,path_to_plot,center_event, station_idx=0):
    template_waveforms = template['waveforms'][station_idx, :, :]
    detection_waveforms = detections['waveforms'][:, station_idx, :, :]
    # only select the n_best best detections
    I = np.argsort(detections['metadata']['correlation_coefficients'])[::-1]
    I = I[:n_best]
    chronological_order = np.argsort(detections['metadata']['origin_times'][I])
    I = I[chronological_order]
    # select subsets
    detection_waveforms = detection_waveforms[I, :, :]
    OT = detections['metadata']['origin_times'][I]
    CC = detections['metadata']['correlation_coefficients'][I]
    n_components = detection_waveforms.shape[1]
    # define some useful variables before starting plotting
    buffer_extracted_events = np.int32(10. * template['sampling_rate'])
    moveouts = np.array([template['moveouts_S'][station_idx],
                         template['moveouts_S'][station_idx],
                         template['moveouts_P'][station_idx]])
    duration = template['waveforms'].shape[-1]
    components = ['N', 'E', 'Z']
    # start plotting
    time = np.linspace(0., 8., duration)
    figsize = (28, 12)
    plt.figure('detections', figsize=figsize)
    plt.suptitle('Station {}'.format(template['stations'][station_idx].decode('utf-8')))
    for c in range(n_components):
        plt.subplot(n_best + 1, n_components, 1 + c)
        plt.title(components[c])
        plt.plot(time, template_waveforms[c, :], lw=3, color='C3', label='Template')
        plt.xlim(time.min(), time.max())
        if c == 2:
            plt.legend(bbox_to_anchor=(1.05, 1), loc='upper left', borderaxespad=0., handlelength=0.1)
        for n in range(n_best):
            plt.subplot(n_best + 1, n_components, (1 + n)*n_components + c + 1)
            idx_start = int(buffer_extracted_events + moveouts[c])
            idx_end = int(idx_start + duration)
            plt.plot(time, detection_waveforms[n, c, idx_start:idx_end], lw=3., color='C0',
                     label=udt(OT[n])\
                             .strftime('%Y,%m,%d -- %H:%M:%S'))
            plt.xlim(time.min(), time.max())
            if c == 2:
                plt.legend(bbox_to_anchor=(1.05, 1), loc='upper left', borderaxespad=0.)
            if c == 0:
                plt.text(0.05, 0.5, 'CC = {:.2f}'.format(CC[n]), transform=plt.gca().transAxes,
                        bbox=dict(facecolor='white', alpha=0.5))
            if n == n_best - 1:
                plt.xlabel('Time (s)')
    plt.subplots_adjust(top=0.91,
            bottom=0.075,
            left=0.06,
            right=0.885,
            hspace=0.2,
            wspace=0.2)
    plt.savefig(path_to_plot+str(n_best)+'_'+center_event+'_new_detections.png')
    plt.close()
