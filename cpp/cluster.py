import sys
from sklearn.cluster import KMeans

def cluster(seq_kmer_counts):
    kmeans = KMeans(n_clusters=1, random_state=2).fit(seq_kmer_counts)
    pre_cluster_inertia = kmeans.inertia_

    '''
    if pre_cluster_inertia == 0:
        logging.debug("pre_cluster_intertia is 0!")
        for key in list(interval_seq_dict.keys()):
            logging.debug("seq: %s, num_seqs with this seq: %d", key, len(interval_seq_dict[key]))
    '''

    cluster_inertia = pre_cluster_inertia
    number_of_clusters = 1
    #logging.debug("number of clusters: %d, inertia: %f", number_of_clusters, cluster_inertia)
    while (cluster_inertia > 0
           and cluster_inertia > pre_cluster_inertia / 2 #we cluster until we reach less than half of the initial cluster inertia
           and number_of_clusters <= len(seq_kmer_counts)):
        number_of_clusters += 1
        kmeans = KMeans(n_clusters=number_of_clusters, random_state=2).fit(seq_kmer_counts)
        cluster_inertia = kmeans.inertia_
        #logging.debug("number of clusters: %d, inertia: %f", number_of_clusters, cluster_inertia)

    # now extract the equivalence class details from this partition and return
    #logging.debug("Extract equivalence classes from this partition")
    big_return_id_lists = []
    if pre_cluster_inertia > 0:
        equiv_class_ids = list(kmeans.predict(seq_kmer_counts))
        for i in range(max(equiv_class_ids) + 1):
            big_return_id_lists.append([])
        for i, val in enumerate(equiv_class_ids):
            big_return_id_lists[val].append(i)
    else:
        #logging.debug("default to not clustering")
        big_return_id_lists = [[i] for i in range(len(seq_kmer_counts))]

    return big_return_id_lists

with open(sys.argv[1], "r") as clusterInputFile:
    seq_kmer_counts = []
    for line in clusterInputFile:
        seq_kmer_counts.append([int(x) for x in line.split()])

    clusters = cluster(seq_kmer_counts)

    with open(sys.argv[2], "w") as clusterOutputFile:
        clusterOutputFile.write(clusters)