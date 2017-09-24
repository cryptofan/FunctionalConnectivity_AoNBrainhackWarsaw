Social Computing Data Repository - Basic Information
==========================================================================
Dataset Name: YouTube 
Abstract: YouTube is a multi-dimensional network consists of various type of interactions.
Number of Nodes: 15,088
Number of Interactions: 5
Maximum Number of Edges: 5,574,249
Missing Values: No

Source:
==========================================================================
Lei Tang, Huan Liu
Computer Science and Engineering, 
Arizona State University. 
Email: L.Tang@asu.edu, Huan.Liu@asu.edu


Data Set Information:
==========================================================================
[I]. Brief description
This is the data set crawled on Dec, 2008 from YouTube. (http://www.youtube.com/). YouTube is a video sharing site  where various interactions occur between users.  In particular, we crawled 30, 522 user profiles.  For each user, we crawl his/her contacts, subscriptions and favorite videos. To avoid sample selection bias, we choose authors of 100 recently uploaded videos as seed set. This crawling reaches in total 848, 003 users and 1,299,642 videos. However, not all users sharing all kinds of information. After removing those users, we have 15, 088 active user profiles. 

Based on the crawled information, we construct 5 different interactions between the 15, 088 users.  Specifically, they are:
1. the contact network between the 15,088 users; 
2. the number of shared friends between two users in the 848, 003 (excluding the 15,088) contacts;
3. the number of shared subscriptions between two users;
4. the number of shared subscribers between two users;
5. the number of shared favoriate videos.

Details can be found in the related reference.

[II]. Basic statistics
Number of users : 15,088
Number of Interaction Types: 5

[III]. The data format

6 files are included:

1. [1-5]-edges.csv
-- they are the csv format of interactions. Each csv file represents one type of interaction. It is composed of three columns, with the first two representing the user ids, and the last representing the intensity of interaction. Here is an example:

   1,58,3

The interaction intensity between users 1 and 58 is 3.

Our network is symmetric, so we only show the interaction once. That is,  58,1,3 will not show up if 1,58,3 is already there.

2. nodes.csv
-- it's the file of all the users. This file works as a dictionary of all the users in this data set. It's useful for fast reference. It contains
all the node ids used in the dataset

Relevant Papers:
==========================================================================

[1]. Lei Tang, Xufei Wang, and Huan Liu. "Uncovering Groups via Heterogeneous Interaction Analysis", IEEE International Conference on Data Mining (ICDM09), Dec. 6-9, 2009. Miami Florida.
