% Zeitableitung der translatorischen Teilmatrix der analytischen Jacobi-Matrix für Segment Nr. 6 (0=Basis) von
% S6PRRRRR6
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
%   pkin=[a2,a3,a4,a5,a6,alpha2,alpha3,alpha4,d2,d3,d4,d5,d6,theta1]';
%
% Output:
% JaD_transl [3x6]
%   Zeitableitung der translatorischen Teilmatrix der analytischen Jacobi-Matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-02-26 20:22
% Revision: d75aae1ac561373cd3be920984c3df29a1c2ecc8 (2019-02-26)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function JaD_transl = S6PRRRRR6_jacobiaD_transl_6_sym_varpar(qJ, qJD, r_i_i_C, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(3,1),zeros(14,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6PRRRRR6_jacobiaD_transl_6_sym_varpar: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6PRRRRR6_jacobiaD_transl_6_sym_varpar: qJD has to be [6x1] (double)');
assert(isa(r_i_i_C,'double') && isreal(r_i_i_C) && all(size(r_i_i_C) == [3 1]), ...
	'S6PRRRRR6_jacobiaD_transl_6_sym_varpar: Position vector r_i_i_C has to be [3x1] double');
assert(isreal(pkin) && all(size(pkin) == [14 1]), ...
  'S6PRRRRR6_jacobiaD_transl_6_sym_varpar: pkin has to be [14x1] (double)');

%% Symbolic Calculation
% From jacobiaD_transl_6_floatb_twist_matlab.m
% OptimizationMode: 2
% StartTime: 2019-02-26 20:22:00
% EndTime: 2019-02-26 20:22:03
% DurationCPUTime: 2.95s
% Computational Cost: add. (4964->335), mult. (16032->574), div. (0->0), fcn. (18521->18), ass. (0->192)
t1001 = r_i_i_C(3) + pkin(13);
t893 = sin(pkin(6));
t895 = cos(pkin(7));
t901 = sin(qJ(2));
t905 = cos(qJ(3));
t971 = t901 * t905;
t900 = sin(qJ(3));
t906 = cos(qJ(2));
t972 = t900 * t906;
t930 = t895 * t972 + t971;
t919 = t930 * qJD(3);
t931 = -t895 * t971 - t972;
t910 = (t931 * qJD(2) - t919) * t893;
t892 = sin(pkin(7));
t896 = cos(pkin(6));
t980 = t892 * t896;
t961 = t900 * t980;
t952 = qJD(3) * t961;
t1007 = t910 - t952;
t998 = cos(pkin(14));
t955 = t998 * t906;
t890 = sin(pkin(14));
t986 = t890 * t901;
t879 = t896 * t955 - t986;
t875 = t879 * qJD(2);
t956 = t998 * t901;
t985 = t890 * t906;
t922 = -t896 * t956 - t985;
t876 = t922 * qJD(2);
t957 = t893 * t998;
t951 = t892 * t957;
t923 = -t879 * t895 + t951;
t974 = t895 * t905;
t989 = t922 * t905;
t909 = (t923 * t900 + t989) * qJD(3) - t875 * t900 + t876 * t974;
t925 = t896 * t986 - t955;
t924 = t896 * t985 + t956;
t987 = t890 * t893;
t934 = t892 * t987 - t895 * t924;
t857 = t934 * t900 - t905 * t925;
t877 = t924 * qJD(2);
t878 = t925 * qJD(2);
t1006 = -t857 * qJD(3) + t877 * t900 + t878 * t974;
t897 = sin(qJ(6));
t902 = cos(qJ(6));
t928 = qJD(6) * (r_i_i_C(1) * t897 + r_i_i_C(2) * t902);
t854 = t900 * t922 - t923 * t905;
t975 = t895 * t900;
t855 = t879 * t975 - t900 * t951 - t989;
t899 = sin(qJ(4));
t904 = cos(qJ(4));
t891 = sin(pkin(8));
t926 = -t879 * t892 - t895 * t957;
t917 = t926 * t891;
t894 = cos(pkin(8));
t976 = t894 * t904;
t1005 = t854 * t976 - t855 * t899 + t904 * t917;
t973 = t900 * t901;
t962 = t895 * t973;
t970 = t905 * t906;
t929 = t962 - t970;
t872 = t929 * t893;
t871 = t931 * t893;
t979 = t892 * t901;
t960 = t891 * t979;
t954 = t893 * t960;
t927 = t871 * t894 + t954;
t1004 = t872 * t899 + t927 * t904;
t936 = t905 * t924 - t925 * t975;
t862 = t900 * t924 + t925 * t974;
t984 = t891 * t892;
t940 = t862 * t894 - t925 * t984;
t1003 = t899 * t936 + t940 * t904;
t937 = -t879 * t905 - t922 * t975;
t860 = -t879 * t900 + t922 * t974;
t941 = t860 * t894 - t922 * t984;
t1002 = t899 * t937 + t941 * t904;
t1000 = pkin(10) * t892;
t999 = pkin(11) * t891;
t996 = t857 * t899;
t932 = t895 * t970 - t973;
t858 = (-t932 * qJD(2) + t929 * qJD(3)) * t893;
t995 = t858 * t891;
t869 = t930 * t893 + t961;
t992 = t869 * t899;
t990 = t876 * t900;
t898 = sin(qJ(5));
t983 = t891 * t898;
t903 = cos(qJ(5));
t982 = t891 * t903;
t981 = t892 * t894;
t978 = t893 * t906;
t977 = t894 * t899;
t969 = qJD(2) * t893;
t968 = qJD(3) * t900;
t967 = qJD(3) * t905;
t966 = qJD(6) * t897;
t965 = qJD(6) * t902;
t963 = t904 * t984;
t959 = t894 * t979;
t958 = t906 * t969;
t953 = t892 * t958;
t948 = t899 * t952;
t812 = t855 * t904 + (t854 * t894 + t917) * t899;
t839 = -t854 * t891 + t926 * t894;
t796 = t812 * t903 + t839 * t898;
t947 = -t812 * t898 + t839 * t903;
t856 = t900 * t925 + t934 * t905;
t935 = t892 * t924 + t895 * t987;
t921 = t935 * t891;
t915 = t856 * t894 + t921;
t814 = t857 * t904 + t915 * t899;
t840 = -t856 * t891 + t935 * t894;
t798 = t814 * t903 + t840 * t898;
t946 = -t814 * t898 + t840 * t903;
t820 = t941 * t899 - t904 * t937;
t844 = -t860 * t891 - t922 * t981;
t803 = t820 * t903 + t844 * t898;
t822 = t940 * t899 - t904 * t936;
t845 = -t862 * t891 - t925 * t981;
t804 = t822 * t903 + t845 * t898;
t868 = t932 * t893 + t905 * t980;
t933 = -t892 * t978 + t896 * t895;
t920 = t933 * t891;
t914 = t868 * t894 + t920;
t834 = t869 * t904 + t914 * t899;
t853 = -t868 * t891 + t933 * t894;
t810 = t834 * t903 + t853 * t898;
t945 = -t834 * t898 + t853 * t903;
t847 = -t872 * t904 + t927 * t899;
t864 = -t871 * t891 + t893 * t959;
t827 = t847 * t903 + t864 * t898;
t944 = r_i_i_C(1) * t902 - r_i_i_C(2) * t897 + pkin(5);
t824 = t854 * t904 - t855 * t977;
t807 = t824 * t903 + t855 * t983;
t826 = t856 * t904 - t857 * t977;
t808 = t826 * t903 + t857 * t983;
t835 = t937 * qJD(3) - t875 * t974 - t990;
t817 = -t835 * t891 + t875 * t981;
t943 = -t835 * t894 - t875 * t984;
t837 = t936 * qJD(3) + t877 * t974 - t878 * t900;
t818 = -t837 * t891 - t877 * t981;
t942 = -t837 * t894 + t877 * t984;
t843 = t868 * t904 - t869 * t977;
t830 = t843 * t903 + t869 * t983;
t823 = t854 * t899 + t855 * t976;
t825 = t856 * t899 + t857 * t976;
t842 = t868 * t899 + t869 * t976;
t918 = t858 * t894 + t891 * t953;
t913 = -t1001 * t898 - t944 * t903 - pkin(4);
t908 = t903 * t928 + (-t1001 * t903 + t944 * t898) * qJD(5);
t907 = t894 * t909;
t859 = (-t930 * qJD(2) + t931 * qJD(3)) * t893;
t852 = -t962 * t969 - t893 * t901 * t968 + (t958 + (t895 * t978 + t980) * qJD(3)) * t905;
t848 = t894 * t953 - t995;
t841 = t891 * t952 + (t891 * t919 + (-t931 * t891 + t959) * qJD(2)) * t893;
t838 = t862 * qJD(3) + t877 * t975 + t878 * t905;
t836 = t860 * qJD(3) - t875 * t975 + t876 * t905;
t833 = -t868 * t976 - t904 * t920 + t992;
t832 = t878 * t975 + t925 * t968 + (t934 * qJD(3) - t877) * t905;
t831 = t875 * t905 + t922 * t968 - t951 * t967 + (t879 * t967 + t990) * t895;
t816 = -t1006 * t891 - t878 * t981;
t815 = -t876 * t981 - t891 * t909;
t813 = -t856 * t976 - t904 * t921 + t996;
t806 = qJD(4) * t1004 + t859 * t904 + t918 * t899;
t805 = t847 * qJD(4) + t859 * t899 - t918 * t904;
t802 = -t842 * qJD(4) + t1007 * t904 - t852 * t977;
t801 = t843 * qJD(4) + t852 * t976 + t899 * t910 - t948;
t800 = -t894 * t948 + t852 * t904 + (-t894 * t919 + (t931 * t894 + t960) * qJD(2)) * t899 * t893 + (t914 * t904 - t992) * qJD(4);
t799 = -qJD(2) * t904 * t954 + t834 * qJD(4) - t1007 * t976 + t852 * t899;
t794 = -t825 * qJD(4) + t1006 * t904 - t832 * t977;
t793 = t826 * qJD(4) + t1006 * t899 + t832 * t976;
t792 = -t823 * qJD(4) - t831 * t977 + t904 * t909;
t791 = t824 * qJD(4) + t831 * t976 + t899 * t909;
t790 = qJD(4) * t1003 + t838 * t904 - t942 * t899;
t789 = t822 * qJD(4) + t838 * t899 + t942 * t904;
t788 = qJD(4) * t1002 + t836 * t904 - t943 * t899;
t787 = t820 * qJD(4) + t836 * t899 + t943 * t904;
t786 = t806 * t903 + t848 * t898 + (-t847 * t898 + t864 * t903) * qJD(5);
t784 = t852 * t983 + t802 * t903 + (-t843 * t898 + t869 * t982) * qJD(5);
t782 = t832 * t904 + (t1006 * t894 - t878 * t984) * t899 + (t915 * t904 - t996) * qJD(4);
t781 = t814 * qJD(4) - t1006 * t976 + t832 * t899 + t878 * t963;
t780 = t831 * t904 + (-t876 * t984 + t907) * t899 + t1005 * qJD(4);
t779 = t812 * qJD(4) + t831 * t899 + t876 * t963 - t904 * t907;
t778 = t945 * qJD(5) + t800 * t903 + t841 * t898;
t776 = t832 * t983 + t794 * t903 + (-t826 * t898 + t857 * t982) * qJD(5);
t774 = t831 * t983 + t792 * t903 + (-t824 * t898 + t855 * t982) * qJD(5);
t772 = t790 * t903 + t818 * t898 + (-t822 * t898 + t845 * t903) * qJD(5);
t770 = t788 * t903 + t817 * t898 + (-t820 * t898 + t844 * t903) * qJD(5);
t768 = t946 * qJD(5) + t782 * t903 + t816 * t898;
t766 = t947 * qJD(5) + t780 * t903 + t815 * t898;
t1 = [0 (t772 * t902 + t789 * t897) * r_i_i_C(1) + (-t772 * t897 + t789 * t902) * r_i_i_C(2) + t772 * pkin(5) + t790 * pkin(4) + t789 * pkin(12) + t838 * pkin(3) + t878 * pkin(2) - t877 * t1000 + t1001 * (t804 * qJD(5) + t790 * t898 - t818 * t903) + ((-t1003 * t902 - t804 * t897) * r_i_i_C(1) + (t1003 * t897 - t804 * t902) * r_i_i_C(2)) * qJD(6) + t818 * pkin(11) (t776 * t902 + t793 * t897) * r_i_i_C(1) + (-t776 * t897 + t793 * t902) * r_i_i_C(2) + t776 * pkin(5) + t794 * pkin(4) + t793 * pkin(12) + t832 * t999 + t1001 * (t808 * qJD(5) + t794 * t898 - t832 * t982) + ((-t808 * t897 + t825 * t902) * r_i_i_C(1) + (-t808 * t902 - t825 * t897) * r_i_i_C(2)) * qJD(6) + t1006 * pkin(3) (t782 * t897 + t814 * t965) * r_i_i_C(1) + (t782 * t902 - t814 * t966) * r_i_i_C(2) + t782 * pkin(12) + t913 * t781 + t908 * t813, t1001 * t768 - t946 * t928 + t944 * (-qJD(5) * t798 - t782 * t898 + t816 * t903) (-t768 * t897 + t781 * t902) * r_i_i_C(1) + (-t768 * t902 - t781 * t897) * r_i_i_C(2) + ((-t798 * t902 - t813 * t897) * r_i_i_C(1) + (t798 * t897 - t813 * t902) * r_i_i_C(2)) * qJD(6); 0 (t770 * t902 + t787 * t897) * r_i_i_C(1) + (-t770 * t897 + t787 * t902) * r_i_i_C(2) + t770 * pkin(5) + t788 * pkin(4) + t787 * pkin(12) + t836 * pkin(3) + t876 * pkin(2) + t875 * t1000 + t1001 * (t803 * qJD(5) + t788 * t898 - t817 * t903) + ((-t1002 * t902 - t803 * t897) * r_i_i_C(1) + (t1002 * t897 - t803 * t902) * r_i_i_C(2)) * qJD(6) + t817 * pkin(11) (t774 * t902 + t791 * t897 + (-t807 * t897 + t823 * t902) * qJD(6)) * r_i_i_C(1) + (-t774 * t897 + t791 * t902 + (-t807 * t902 - t823 * t897) * qJD(6)) * r_i_i_C(2) + t774 * pkin(5) + t792 * pkin(4) + t791 * pkin(12) + t909 * pkin(3) + t831 * t999 + t1001 * (t807 * qJD(5) + t792 * t898 - t831 * t982) (t780 * t897 + t812 * t965) * r_i_i_C(1) + (t780 * t902 - t812 * t966) * r_i_i_C(2) + t780 * pkin(12) + t913 * t779 - t908 * t1005, t1001 * t766 - t947 * t928 + t944 * (-qJD(5) * t796 - t780 * t898 + t815 * t903) (-t766 * t897 + t779 * t902) * r_i_i_C(1) + (-t766 * t902 - t779 * t897) * r_i_i_C(2) + ((t1005 * t897 - t796 * t902) * r_i_i_C(1) + (t1005 * t902 + t796 * t897) * r_i_i_C(2)) * qJD(6); 0 (t786 * t902 + t805 * t897) * r_i_i_C(1) + (-t786 * t897 + t805 * t902) * r_i_i_C(2) + t786 * pkin(5) + t806 * pkin(4) + t805 * pkin(12) + t859 * pkin(3) - pkin(11) * t995 + t1001 * (t827 * qJD(5) + t806 * t898 - t848 * t903) + ((-t1004 * t902 - t827 * t897) * r_i_i_C(1) + (t1004 * t897 - t827 * t902) * r_i_i_C(2)) * qJD(6) + (-pkin(2) * t901 + (pkin(11) * t894 + pkin(10)) * t906 * t892) * t969 (t784 * t902 + t801 * t897) * r_i_i_C(1) + (-t784 * t897 + t801 * t902) * r_i_i_C(2) + t784 * pkin(5) + t802 * pkin(4) + t801 * pkin(12) + t852 * t999 + t1001 * (t830 * qJD(5) + t802 * t898 - t852 * t982) + ((-t830 * t897 + t842 * t902) * r_i_i_C(1) + (-t830 * t902 - t842 * t897) * r_i_i_C(2)) * qJD(6) + t1007 * pkin(3) (t800 * t897 + t834 * t965) * r_i_i_C(1) + (t800 * t902 - t834 * t966) * r_i_i_C(2) + t800 * pkin(12) + t913 * t799 + t908 * t833, t1001 * t778 - t945 * t928 + t944 * (-qJD(5) * t810 - t800 * t898 + t841 * t903) (-t778 * t897 + t799 * t902) * r_i_i_C(1) + (-t778 * t902 - t799 * t897) * r_i_i_C(2) + ((-t810 * t902 - t833 * t897) * r_i_i_C(1) + (t810 * t897 - t833 * t902) * r_i_i_C(2)) * qJD(6);];
JaD_transl  = t1;
