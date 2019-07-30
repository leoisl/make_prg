#!/usr/bin/env python3
from Bio import AlignIO
import logging
from sklearn.cluster import KMeans
import numpy as np
import sys

def remove_duplicates(seqs):
    seen = set()
    for x in seqs:
        if x in seen:
            continue
        seen.add(x)
        yield x


def cluster(msa_file, min_match_length):
    interval_alignment = AlignIO.read(msa_file, "fasta") #get the MSA of the (sub)-alignment
    interval_seq_dict = {} #alignment seq -> sequence id if alignment seq >= self.min_match_length
    small_interval_seq_dict = {} #alignment seq -> sequence id if alignment seq < self.min_match_length
    seq_dict_keys = [] #all keys (alignment seqs)

    for record in interval_alignment: #add the alignments to the inverval_seq_dict or to the small_interval_seq_dict, depending on their size without -
        seq = str(record.seq).replace('-', '') #remove "-" from the alignment
        if seq in list(interval_seq_dict.keys()):
            interval_seq_dict[seq].append(record.id) #just add a new id to the list
        elif seq in list(small_interval_seq_dict.keys()):
            small_interval_seq_dict[seq].append(record.id) #just add a new id to the list
        elif len(seq) >= min_match_length:
            interval_seq_dict[seq] = [record.id]
            seq_dict_keys.append(seq)
        else:
            small_interval_seq_dict[seq] = [record.id]
            seq_dict_keys.append(seq)

    assert len(seq_dict_keys) == len(
        list(remove_duplicates(seq_dict_keys))), "error, have duplicate dictionary keys"
    assert len([key for key in list(interval_seq_dict.keys()) if
                key in list(small_interval_seq_dict.keys())]) == 0, "error, should have no overlap of keys"
    assert len([key for key in list(small_interval_seq_dict.keys()) if
                key in list(interval_seq_dict.keys())]) == 0, "error, should have no overlap of keys"

    logging.debug("Add classes corresponding to %d small sequences" % len(list(small_interval_seq_dict.keys())))

    logging.debug("Now add classes corresponding to %d longer sequences" % len(list(interval_seq_dict.keys())))
    interval_seqs = list(interval_seq_dict.keys())
    big_return_id_lists = []
    if len(interval_seqs) > 1:
        # first transform sequences into kmer occurance vectors using a dict
        logging.debug("First transform sequences into kmer occurance vectors")

        # make dict based on all kmers in all sequences
        kmer_dict = {} #associate each kmer to an ID
        n = 0 #n = kmer ID
        for j, seq in enumerate(interval_seqs): #goes through all large enough seqs
            for i in range(len(seq) - min_match_length + 1): #goes through all kmers
                kmer = seq[i:i + min_match_length] #get the kmer
                if kmer not in kmer_dict.keys(): #is this kmer new?
                    kmer_dict[kmer] = n #associate the kmer to an id
                    n += 1
        logging.debug("self.kmer_dict = %s"%str(kmer_dict))
        logging.debug("These vectors have length %d" % n)


        # transform to vectors using dict
        # describes the kmers of each sequence as a kmer spectrum (each kmer has an ID, and seq_kmer_counts denotes the occurence of each kmer in each sequence)
        seq_kmer_counts = np.zeros(shape=(len(interval_seqs), n))
        for j, seq in enumerate(interval_seqs):
            counts = np.zeros(n)
            for i in range(len(seq) - min_match_length + 1):
                kmer = seq[i:i + min_match_length]
                counts[kmer_dict[kmer]] += 1
            seq_kmer_counts[j] = counts

        logging.debug("seq_kmer_counts = %s" % str(seq_kmer_counts))

        # cluster sequences using kmeans
        logging.debug("Now cluster:")
        kmeans = KMeans(n_clusters=1, random_state=2).fit(seq_kmer_counts)
        pre_cluster_inertia = kmeans.inertia_

        if pre_cluster_inertia == 0:
            logging.debug("pre_cluster_intertia is 0!")
            for key in list(interval_seq_dict.keys()):
                logging.debug("seq: %s, num_seqs with this seq: %d", key, len(interval_seq_dict[key]))

        cluster_inertia = pre_cluster_inertia
        number_of_clusters = 1
        logging.debug("number of clusters: %d, inertia: %f", number_of_clusters, cluster_inertia)
        while (cluster_inertia > 0
               and cluster_inertia > pre_cluster_inertia / 2 #we cluster until we reach less than half of the initial cluster inertia
               and number_of_clusters <= len(interval_seqs)):
            number_of_clusters += 1
            kmeans = KMeans(n_clusters=number_of_clusters, random_state=2).fit(seq_kmer_counts)
            cluster_inertia = kmeans.inertia_
            logging.debug("number of clusters: %d, inertia: %f", number_of_clusters, cluster_inertia)

        # now extract the equivalence class details from this partition and return
        logging.debug("Extract equivalence classes from this partition")
        if pre_cluster_inertia > 0:
            equiv_class_ids = list(kmeans.predict(seq_kmer_counts))
            for i in range(max(equiv_class_ids) + 1):
                big_return_id_lists.append([])
            for i, val in enumerate(equiv_class_ids):
                big_return_id_lists[val].extend(interval_seq_dict[interval_seqs[i]])
        else:
            logging.debug("default to not clustering")
            big_return_id_lists = [interval_seq_dict[key] for key in interval_seq_dict.keys()]
    elif len(interval_seqs) == 1:
        big_return_id_lists = [interval_seq_dict[interval_seqs[0]]]

    # now merge big and small return_id_lists so as to maintain the order of seqs before
    logging.debug("Merge return id lists for the partitions")
    return_id_lists = []
    added_ids = []
    big_keys = list(interval_seq_dict.keys())
    small_keys = list(small_interval_seq_dict.keys())
    for seq in seq_dict_keys:
        if seq in small_keys:
            logging.debug("add (small) return ids: %s" % small_interval_seq_dict[seq])
            return_id_lists.append(small_interval_seq_dict[seq])
        elif seq in big_keys:
            not_added = [nid for nid in interval_seq_dict[seq] if nid not in added_ids]
            if len(not_added) == len(interval_seq_dict[seq]):
                logging.debug("want to add (big) return ids: %s" % interval_seq_dict[seq])
                for i in range(len(big_return_id_lists)):
                    if interval_seq_dict[seq][0] in big_return_id_lists[i]:
                        logging.debug("add (big) return ids %d: %s" % (i, big_return_id_lists[i]))
                        return_id_lists.append(big_return_id_lists[i])
                        added_ids.extend(return_id_lists[-1])
                        break
            else:
                assert len(
                    not_added) == 0, "Equivalent sequences should be in same part of partition and are not"
        else:
            logging.warning("Key %s doesn't seem to be in either big keys or small keys")
    assert len(interval_alignment) == sum([len(i) for i in return_id_lists]), \
        "I seem to have lost (or gained?) some sequences in the process of clustering"
    assert len(return_id_lists) > 1, \
        "should have some alternate alleles, not only one sequence, this is a non-match interval"


    for cluster in return_id_lists:
        for sequence_index in cluster:
            sequence_index = int(sequence_index)
            print(interval_alignment[sequence_index])

cluster(sys.argv[1], int(sys.argv[2]))