% Zeitableitung der translatorischen Teilmatrix der analytischen Jacobi-Matrix für Segment Nr. 6 (0=Basis) von
% S6RPRRRR12
% Use Code from Maple symbolic Code Generation
%
% analytische Jacobi-Matrix: Differentieller Zusammenhang zwischen
% Endeffektorposition und verallgemeinerten Koordinaten.
%
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% qJD [6x1]
%   Generalized joint velocities
% r_i_i_C [3x1]
%   Ortsvektor vom KörperKS-Ursprung zum gesuchten Punkt
% pkin [14x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,alpha2,alpha3,alpha4,d1,d3,d4,d5,d6,theta2]';
%
% Output:
% JaD_transl [3x6]
%   Zeitableitung der translatorischen Teilmatrix der analytischen Jacobi-Matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-02-26 21:21
% Revision: d75aae1ac561373cd3be920984c3df29a1c2ecc8 (2019-02-26)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function JaD_transl = S6RPRRRR12_jacobiaD_transl_6_sym_varpar(qJ, qJD, r_i_i_C, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(3,1),zeros(14,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPRRRR12_jacobiaD_transl_6_sym_varpar: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RPRRRR12_jacobiaD_transl_6_sym_varpar: qJD has to be [6x1] (double)');
assert(isa(r_i_i_C,'double') && isreal(r_i_i_C) && all(size(r_i_i_C) == [3 1]), ...
	'S6RPRRRR12_jacobiaD_transl_6_sym_varpar: Position vector r_i_i_C has to be [3x1] double');
assert(isreal(pkin) && all(size(pkin) == [14 1]), ...
  'S6RPRRRR12_jacobiaD_transl_6_sym_varpar: pkin has to be [14x1] (double)');

%% Symbolic Calculation
% From jacobiaD_transl_6_floatb_twist_matlab.m
% OptimizationMode: 2
% StartTime: 2019-02-26 21:21:20
% EndTime: 2019-02-26 21:21:23
% DurationCPUTime: 3.77s
% Computational Cost: add. (5036->259), mult. (16095->427), div. (0->0), fcn. (18656->18), ass. (0->156)
t876 = sin(qJ(3));
t880 = cos(qJ(3));
t877 = sin(qJ(1));
t868 = sin(pkin(14));
t971 = cos(pkin(6));
t944 = t868 * t971;
t970 = cos(pkin(14));
t975 = cos(qJ(1));
t910 = -t877 * t970 - t975 * t944;
t931 = t971 * t970;
t860 = t877 * t868 - t975 * t931;
t872 = cos(pkin(7));
t870 = sin(pkin(6));
t969 = sin(pkin(7));
t943 = t870 * t969;
t914 = t860 * t872 + t975 * t943;
t839 = -t876 * t910 + t914 * t880;
t840 = t914 * t876 + t910 * t880;
t871 = cos(pkin(8));
t875 = sin(qJ(4));
t869 = sin(pkin(8));
t898 = t860 * t969;
t949 = t870 * t975;
t915 = -t872 * t949 + t898;
t908 = t915 * t869;
t974 = cos(qJ(4));
t813 = t840 * t974 + (t839 * t871 - t908) * t875;
t829 = t839 * t869 + t915 * t871;
t874 = sin(qJ(5));
t879 = cos(qJ(5));
t796 = t813 * t879 - t829 * t874;
t873 = sin(qJ(6));
t1005 = t796 * t873;
t878 = cos(qJ(6));
t1004 = t796 * t878;
t911 = t975 * t868 + t877 * t931;
t858 = t911 * qJD(1);
t945 = t858 * t969;
t957 = qJD(1) * t877;
t947 = t870 * t957;
t916 = t872 * t947 + t945;
t912 = -t877 * t944 + t975 * t970;
t859 = t912 * qJD(1);
t937 = t877 * t943;
t924 = qJD(1) * t937;
t913 = -t858 * t872 + t924;
t986 = t840 * qJD(3) - t859 * t876 + t913 * t880;
t805 = t869 * t986 - t916 * t871;
t1003 = t796 * qJD(5) - t805 * t879;
t929 = t813 * t874 + t829 * t879;
t1002 = t929 * qJD(5) - t805 * t874;
t999 = t813 * qJD(4);
t909 = t916 * t869;
t995 = (t871 * t986 + t909) * t875;
t948 = t871 * t974;
t992 = t839 * t948 - t840 * t875;
t976 = r_i_i_C(3) + pkin(13);
t984 = pkin(10) * t872 + qJ(2);
t933 = qJD(3) * t943;
t979 = t914 * qJD(1) - qJD(3) * t912;
t904 = qJD(1) * t910;
t906 = t911 * t872;
t981 = -qJD(3) * t906 + t904;
t885 = -t979 * t880 + (t877 * t933 + t981) * t876;
t895 = -t906 + t937;
t841 = -t912 * t876 + t895 * t880;
t951 = t841 * t974;
t983 = qJD(4) * t951 - t885 * t875;
t892 = qJD(1) * t915;
t881 = t885 * t869 - t871 * t892;
t842 = t895 * t876 + t912 * t880;
t959 = t870 * t877;
t894 = t872 * t959 + t911 * t969;
t893 = t894 * t869;
t815 = t842 * t974 + (t841 * t871 + t893) * t875;
t831 = -t841 * t869 + t894 * t871;
t980 = -t815 * t874 + t831 * t879;
t921 = qJD(6) * (t873 * r_i_i_C(1) + t878 * r_i_i_C(2));
t930 = t971 * t969;
t942 = t872 * t970;
t853 = (t868 * t880 + t876 * t942) * t870 + t876 * t930;
t972 = pkin(11) * t869;
t961 = t869 * t874;
t960 = t869 * t879;
t958 = t871 * t875;
t956 = qJD(3) * t860;
t955 = qJD(3) * t876;
t954 = qJD(4) * t875;
t953 = qJD(6) * t873;
t952 = qJD(6) * t878;
t950 = t869 * t974;
t946 = t842 * t954;
t941 = t842 * t948;
t939 = qJD(1) * t949;
t935 = t880 * t942;
t798 = t815 * t879 + t831 * t874;
t922 = t880 * t930;
t852 = t922 + (-t868 * t876 + t935) * t870;
t907 = t971 * t872 - t970 * t943;
t902 = t907 * t869;
t828 = t853 * t974 + (t852 * t871 + t902) * t875;
t836 = -t852 * t869 + t907 * t871;
t807 = t828 * t879 + t836 * t874;
t928 = -t828 * t874 + t836 * t879;
t927 = t878 * r_i_i_C(1) - t873 * r_i_i_C(2) + pkin(5);
t925 = t880 * t933;
t823 = -t839 * t974 + t840 * t958;
t799 = t823 * t879 - t840 * t961;
t825 = -t842 * t958 + t951;
t800 = t825 * t879 + t842 * t961;
t833 = t852 * t974 - t853 * t958;
t826 = t833 * t879 + t853 * t961;
t919 = t839 * t875 + t840 * t948;
t918 = -t852 * t875 - t853 * t948;
t905 = -t976 * t874 - t927 * t879 - pkin(4);
t901 = -t915 * t950 + t992;
t810 = -t974 * t908 + t992;
t891 = t974 * t893;
t890 = t869 * t892;
t888 = t879 * t921 + (t927 * t874 - t976 * t879) * qJD(5);
t827 = -t852 * t948 + t853 * t875 - t974 * t902;
t886 = t986 * t974;
t882 = t885 * t974;
t863 = t975 * t925;
t850 = t853 * qJD(3);
t849 = t870 * t868 * t955 + (-t870 * t935 - t922) * qJD(3);
t824 = t841 * t875 + t941;
t821 = t863 + (t872 * t956 - t859) * t880 + (-qJD(3) * t910 - t913) * t876;
t819 = t876 * t924 + t859 * t880 + t910 * t955 - t863 + (-t858 * t876 - t880 * t956) * t872;
t818 = t979 * t876 + t877 * t925 + t981 * t880;
t814 = -t841 * t948 + t842 * t875 - t891;
t809 = t918 * qJD(4) + t849 * t958 - t850 * t974;
t808 = t833 * qJD(4) - t849 * t948 - t850 * t875;
t802 = -t827 * qJD(4) - t849 * t974 - t850 * t958;
t801 = t828 * qJD(4) - t849 * t875 + t850 * t948;
t793 = -t849 * t961 + t809 * t879 + (-t833 * t874 + t853 * t960) * qJD(5);
t791 = t919 * qJD(4) - t819 * t958 + t886;
t790 = t823 * qJD(4) + t819 * t948 + t875 * t986;
t789 = -qJD(4) * t941 - t818 * t958 - t841 * t954 - t882;
t788 = t818 * t948 - t871 * t946 + t983;
t787 = t928 * qJD(5) + t802 * t879 + t850 * t961;
t785 = t901 * qJD(4) + t821 * t974 - t995;
t784 = t821 * t875 + t916 * t950 + t948 * t986 + t999;
t783 = -t810 * qJD(4) + t819 * t974 + t995;
t782 = t819 * t875 - t871 * t886 - t974 * t909 - t999;
t781 = qJD(4) * t891 + t818 * t974 + t983 * t871 - t875 * t890 - t946;
t780 = t815 * qJD(4) + t818 * t875 + t871 * t882 + t974 * t890;
t779 = t819 * t961 + t791 * t879 + (-t823 * t874 - t840 * t960) * qJD(5);
t777 = t818 * t961 + t789 * t879 + (-t825 * t874 + t842 * t960) * qJD(5);
t775 = t785 * t879 - t1002;
t773 = t783 * t879 + t1002;
t771 = t980 * qJD(5) + t781 * t879 + t881 * t874;
t770 = t798 * qJD(5) + t781 * t874 - t881 * t879;
t769 = t771 * t878 + t780 * t873 + (-t798 * t873 + t814 * t878) * qJD(6);
t768 = -t771 * t873 + t780 * t878 + (-t798 * t878 - t814 * t873) * qJD(6);
t1 = [(t775 * t878 + t784 * t873) * r_i_i_C(1) + (-t775 * t873 + t784 * t878) * r_i_i_C(2) + t775 * pkin(5) + t785 * pkin(4) + t784 * pkin(12) + t821 * pkin(3) - t859 * pkin(2) - pkin(10) * t945 + qJD(2) * t949 + t976 * (t785 * t874 + t1003) + ((-t878 * t901 - t1005) * r_i_i_C(1) + (t873 * t901 - t1004) * r_i_i_C(2)) * qJD(6) + t805 * pkin(11) + (-t975 * pkin(1) - t984 * t959) * qJD(1), t939 (t777 * t878 + t788 * t873 + (-t800 * t873 + t824 * t878) * qJD(6)) * r_i_i_C(1) + (-t777 * t873 + t788 * t878 + (-t800 * t878 - t824 * t873) * qJD(6)) * r_i_i_C(2) + t777 * pkin(5) + t789 * pkin(4) + t788 * pkin(12) - t885 * pkin(3) + t818 * t972 + t976 * (t800 * qJD(5) + t789 * t874 - t818 * t960) (t781 * t873 + t815 * t952) * r_i_i_C(1) + (t781 * t878 - t815 * t953) * r_i_i_C(2) + t781 * pkin(12) + t905 * t780 + t888 * t814, -t927 * t770 + t976 * t771 - t980 * t921, t768 * r_i_i_C(1) - t769 * r_i_i_C(2); -qJD(1) * pkin(10) * t898 - pkin(1) * t957 + pkin(2) * t904 + t818 * pkin(3) + t781 * pkin(4) + t771 * pkin(5) + t881 * pkin(11) + t780 * pkin(12) + t769 * r_i_i_C(1) + t768 * r_i_i_C(2) + qJD(2) * t959 + t976 * t770 + t984 * t939, t947 (t779 * t878 + t790 * t873 + (-t799 * t873 - t878 * t919) * qJD(6)) * r_i_i_C(1) + (-t779 * t873 + t790 * t878 + (-t799 * t878 + t873 * t919) * qJD(6)) * r_i_i_C(2) + t779 * pkin(5) + t791 * pkin(4) + t790 * pkin(12) + t986 * pkin(3) + t819 * t972 + t976 * (t799 * qJD(5) + t791 * t874 - t819 * t960) (t783 * t873 - t813 * t952) * r_i_i_C(1) + (t783 * t878 + t813 * t953) * r_i_i_C(2) + t783 * pkin(12) + t905 * t782 + t888 * t810, t976 * t773 - t929 * t921 + t927 * (-t783 * t874 + t1003) (-t773 * t873 + t782 * t878) * r_i_i_C(1) + (-t773 * t878 - t782 * t873) * r_i_i_C(2) + ((-t810 * t873 + t1004) * r_i_i_C(1) + (-t810 * t878 - t1005) * r_i_i_C(2)) * qJD(6); 0, 0 (t793 * t878 + t808 * t873) * r_i_i_C(1) + (-t793 * t873 + t808 * t878) * r_i_i_C(2) + t793 * pkin(5) + t809 * pkin(4) + t808 * pkin(12) - t850 * pkin(3) - t849 * t972 + t976 * (t826 * qJD(5) + t809 * t874 + t849 * t960) + ((-t826 * t873 - t878 * t918) * r_i_i_C(1) + (-t826 * t878 + t873 * t918) * r_i_i_C(2)) * qJD(6) (t802 * t873 + t828 * t952) * r_i_i_C(1) + (t802 * t878 - t828 * t953) * r_i_i_C(2) + t802 * pkin(12) + t905 * t801 + t888 * t827, t976 * t787 - t928 * t921 + t927 * (-t807 * qJD(5) - t802 * t874 + t850 * t960) (-t787 * t873 + t801 * t878) * r_i_i_C(1) + (-t787 * t878 - t801 * t873) * r_i_i_C(2) + ((-t807 * t878 - t827 * t873) * r_i_i_C(1) + (t807 * t873 - t827 * t878) * r_i_i_C(2)) * qJD(6);];
JaD_transl  = t1;
