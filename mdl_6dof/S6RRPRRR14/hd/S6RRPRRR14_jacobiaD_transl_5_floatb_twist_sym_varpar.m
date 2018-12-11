% Zeitableitung der translatorischen Teilmatrix der analytischen Jacobi-Matrix für Segment Nr. 5 (0=Basis) von
% S6RRPRRR14
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
%   pkin=[a2,a3,a4,a5,a6,alpha2,alpha3,alpha4,d1,d2,d4,d5,d6,theta3]';
%
% Output:
% JaD_transl [3x6]
%   Zeitableitung der translatorischen Teilmatrix der analytischen Jacobi-Matrix

% Quelle: HybrDyn-Toolbox (ehem. IRT-Maple-Toolbox)
% Datum: 2018-12-10 18:38
% Revision: bb42a8b95257d9bc83910d26e849f5825122f662 (2018-12-05)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für mechatronische Systeme, Universität Hannover

function JaD_transl = S6RRPRRR14_jacobiaD_transl_5_floatb_twist_sym_varpar(qJ, qJD, r_i_i_C, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(3,1),zeros(14,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRPRRR14_jacobiaD_transl_5_floatb_twist_sym_varpar: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RRPRRR14_jacobiaD_transl_5_floatb_twist_sym_varpar: qJD has to be [6x1] (double)');
assert(isa(r_i_i_C,'double') && isreal(r_i_i_C) && all(size(r_i_i_C) == [3 1]), ...
	'S6RRPRRR14_jacobiaD_transl_5_floatb_twist_sym_varpar: Position vector r_i_i_C has to be [3x1] double');
assert(isreal(pkin) && all(size(pkin) == [14 1]), ...
  'S6RRPRRR14_jacobiaD_transl_5_floatb_twist_sym_varpar: pkin has to be [14x1] (double)');

%% Symbolic Calculation
% From jacobiaD_transl_5_floatb_twist_matlab.m
% OptimizationMode: 2
% StartTime: 2018-12-10 18:38:30
% EndTime: 2018-12-10 18:38:32
% DurationCPUTime: 2.29s
% Computational Cost: add. (6780->262), mult. (8221->405), div. (0->0), fcn. (6786->28), ass. (0->150)
t918 = cos(qJ(2));
t919 = cos(qJ(1));
t945 = pkin(6) + qJ(2);
t929 = sin(t945) / 0.2e1;
t946 = pkin(6) - qJ(2);
t936 = sin(t946);
t878 = t929 - t936 / 0.2e1;
t920 = qJD(2) * t878;
t914 = sin(qJ(2));
t915 = sin(qJ(1));
t932 = cos(t945) / 0.2e1;
t938 = cos(t946);
t925 = t938 / 0.2e1 + t932;
t921 = t919 * t914 + t915 * t925;
t950 = qJD(2) * t915;
t842 = qJD(1) * t921 + t918 * t950 + t919 * t920;
t868 = t925 * qJD(2);
t889 = t914 * t950;
t951 = qJD(1) * t915;
t843 = -t878 * t951 - t889 + (qJD(1) * t918 + t868) * t919;
t902 = pkin(7) + pkin(14);
t890 = sin(t902) / 0.2e1;
t903 = pkin(7) - pkin(14);
t900 = sin(t903);
t871 = t890 + t900 / 0.2e1;
t891 = cos(t903) / 0.2e1;
t901 = cos(t902);
t873 = t891 + t901 / 0.2e1;
t904 = sin(pkin(14));
t907 = sin(pkin(6));
t942 = t907 * t951;
t803 = -t842 * t873 - t843 * t904 + t871 * t942;
t872 = t890 - t900 / 0.2e1;
t874 = t891 - t901 / 0.2e1;
t908 = cos(pkin(14));
t804 = -t842 * t872 + t843 * t908 + t874 * t942;
t922 = t919 * t925;
t855 = t915 * t914 - t922;
t856 = t878 * t919 + t915 * t918;
t952 = t907 * t919;
t821 = t855 * t873 + t856 * t904 + t871 * t952;
t822 = t855 * t872 - t856 * t908 + t874 * t952;
t906 = sin(pkin(7));
t910 = cos(pkin(7));
t830 = t842 * t906 + t910 * t942;
t851 = -t855 * t906 + t910 * t952;
t943 = pkin(8) + qJ(4);
t934 = sin(t943);
t927 = t934 / 0.2e1;
t944 = pkin(8) - qJ(4);
t935 = sin(t944);
t928 = t935 / 0.2e1;
t875 = t927 + t928;
t863 = t875 * qJD(4);
t930 = cos(t943) / 0.2e1;
t937 = cos(t944);
t880 = t937 / 0.2e1 + t930;
t865 = t880 * qJD(4);
t876 = t927 - t935 / 0.2e1;
t879 = t930 - t937 / 0.2e1;
t917 = cos(qJ(4));
t913 = sin(qJ(4));
t948 = qJD(4) * t913;
t773 = t803 * t876 + t804 * t917 - t821 * t865 + t822 * t948 - t830 * t879 - t851 * t863;
t905 = sin(pkin(8));
t909 = cos(pkin(8));
t793 = t803 * t905 - t830 * t909;
t912 = sin(qJ(5));
t916 = cos(qJ(5));
t966 = -t773 * t916 + t793 * t912;
t965 = t773 * t912 + t793 * t916;
t788 = t821 * t876 + t822 * t917 - t851 * t879;
t808 = t821 * t905 - t851 * t909;
t964 = t788 * t916 - t808 * t912;
t963 = -t788 * t912 - t808 * t916;
t864 = (t928 - t934 / 0.2e1) * qJD(4);
t866 = t879 * qJD(4);
t947 = qJD(4) * t917;
t962 = t803 * t880 - t804 * t913 - t821 * t864 + t822 * t947 + t830 * t875 - t851 * t866;
t881 = t932 - t938 / 0.2e1;
t949 = qJD(2) * t919;
t840 = qJD(1) * t856 + t915 * t868 + t914 * t949;
t961 = r_i_i_C(3) + pkin(12);
t839 = -qJD(1) * t922 + t914 * t951 + t915 * t920 - t918 * t949;
t811 = -t839 * t904 + t840 * t873;
t960 = t811 * t905;
t860 = t915 * t878 - t919 * t918;
t844 = t860 * qJD(1) - t868 * t919 + t889;
t813 = t842 * t904 + t844 * t873;
t959 = t813 * t905;
t877 = t929 + t936 / 0.2e1;
t867 = t877 * qJD(2);
t869 = t881 * qJD(2);
t847 = -t867 * t873 - t869 * t904;
t958 = t847 * t905;
t957 = t869 * t906;
t956 = t879 * t906;
t955 = t906 * qJ(3);
t954 = t906 * t909;
t953 = t907 * t915;
t941 = qJD(1) * t952;
t940 = qJ(3) * t910 + pkin(10);
t939 = -pkin(11) * t909 - qJ(3);
t801 = t839 * t873 + t840 * t904 + t871 * t941;
t828 = -t839 * t906 + t910 * t941;
t791 = -t801 * t905 + t828 * t909;
t926 = t916 * r_i_i_C(1) - t912 * r_i_i_C(2) + pkin(4);
t853 = t906 * t921 + t910 * t953;
t924 = qJD(5) * (-t912 * r_i_i_C(1) - t916 * r_i_i_C(2));
t802 = t839 * t872 - t840 * t908 + t874 * t941;
t823 = t860 * t904 + t871 * t953 - t873 * t921;
t824 = -t860 * t908 - t872 * t921 + t874 * t953;
t770 = t801 * t876 + t802 * t917 + t823 * t865 - t824 * t948 - t828 * t879 + t853 * t863;
t911 = cos(pkin(6));
t836 = t871 * t911 + t873 * t877 + t881 * t904;
t837 = t872 * t877 + t874 * t911 - t881 * t908;
t845 = -t867 * t904 + t869 * t873;
t846 = t867 * t908 + t869 * t872;
t854 = -t877 * t906 + t910 * t911;
t782 = t836 * t865 - t837 * t948 + t845 * t876 + t846 * t917 + t854 * t863 + t869 * t956;
t850 = t872 * t881 + t877 * t908;
t849 = t873 * t881 - t877 * t904;
t848 = -t867 * t872 + t869 * t908;
t835 = t860 * t872 - t908 * t921;
t834 = t860 * t873 + t904 * t921;
t833 = -t855 * t908 - t856 * t872;
t832 = t855 * t904 - t856 * t873;
t831 = -t849 * t905 - t881 * t954;
t827 = t867 * t954 - t958;
t826 = -t845 * t905 - t869 * t954;
t818 = -t836 * t905 + t854 * t909;
t817 = -t834 * t905 - t860 * t954;
t816 = -t832 * t905 + t856 * t954;
t814 = -t842 * t908 + t844 * t872;
t812 = t839 * t908 + t840 * t872;
t810 = -t823 * t905 + t853 * t909;
t807 = t849 * t876 + t850 * t917 + t881 * t956;
t799 = t836 * t876 + t837 * t917 - t854 * t879;
t797 = t834 * t876 + t835 * t917 + t860 * t956;
t796 = t832 * t876 + t833 * t917 - t856 * t956;
t795 = -t844 * t954 - t959;
t794 = -t840 * t954 - t960;
t790 = t823 * t876 + t824 * t917 - t853 * t879;
t785 = -t850 * t948 + t847 * t876 + t848 * t917 + t849 * t865 + (-t863 * t881 - t867 * t879) * t906;
t780 = -t833 * t948 + t813 * t876 + t814 * t917 + t832 * t865 + (t844 * t879 + t856 * t863) * t906;
t778 = -t835 * t948 + t811 * t876 + t812 * t917 + t834 * t865 + (t840 * t879 - t860 * t863) * t906;
t769 = -t801 * t880 + t802 * t913 - t823 * t864 + t824 * t947 - t828 * t875 - t853 * t866;
t768 = t770 * t916 + t791 * t912 + (-t790 * t912 + t810 * t916) * qJD(5);
t767 = -t770 * t912 + t791 * t916 + (-t790 * t916 - t810 * t912) * qJD(5);
t1 = [t966 * r_i_i_C(1) + t965 * r_i_i_C(2) - t773 * pkin(4) - t804 * pkin(3) - t843 * pkin(2) - t842 * t955 + t961 * t962 + (t963 * r_i_i_C(1) - t964 * r_i_i_C(2)) * qJD(5) + t851 * qJD(3) + t793 * pkin(11) + (-t919 * pkin(1) - t940 * t953) * qJD(1) (t778 * t916 + t794 * t912) * r_i_i_C(1) + (-t778 * t912 + t794 * t916) * r_i_i_C(2) + t778 * pkin(4) + t812 * pkin(3) - pkin(11) * t960 + t839 * pkin(2) + t961 * (t835 * t947 - t811 * t880 + t812 * t913 - t834 * t864 + (t840 * t875 + t860 * t866) * t906) + ((-t797 * t912 + t817 * t916) * r_i_i_C(1) + (-t797 * t916 - t817 * t912) * r_i_i_C(2)) * qJD(5) + (-t860 * qJD(3) + t840 * t939) * t906, t828, t961 * t770 + (t823 * t880 - t824 * t913 + t853 * t875) * t924 - t926 * t769, r_i_i_C(1) * t767 - r_i_i_C(2) * t768, 0; -t839 * t955 - t840 * pkin(2) + t802 * pkin(3) + t770 * pkin(4) + t768 * r_i_i_C(1) + t767 * r_i_i_C(2) + t961 * t769 + t853 * qJD(3) + t791 * pkin(11) + (-pkin(1) * t915 + t940 * t952) * qJD(1) (t780 * t916 + t795 * t912) * r_i_i_C(1) + (-t780 * t912 + t795 * t916) * r_i_i_C(2) + t780 * pkin(4) + t814 * pkin(3) - pkin(11) * t959 - t842 * pkin(2) + t961 * (t833 * t947 - t813 * t880 + t814 * t913 - t832 * t864 + (t844 * t875 - t856 * t866) * t906) + ((-t796 * t912 + t816 * t916) * r_i_i_C(1) + (-t796 * t916 - t816 * t912) * r_i_i_C(2)) * qJD(5) + (qJD(3) * t856 + t844 * t939) * t906, t830, t961 * t773 + (-t821 * t880 + t822 * t913 - t851 * t875) * t924 + t926 * t962, -t965 * r_i_i_C(1) + t966 * r_i_i_C(2) + (t964 * r_i_i_C(1) + t963 * r_i_i_C(2)) * qJD(5), 0; 0 (t785 * t916 + t827 * t912) * r_i_i_C(1) + (-t785 * t912 + t827 * t916) * r_i_i_C(2) + t785 * pkin(4) + t848 * pkin(3) - pkin(11) * t958 + t869 * pkin(2) + t961 * (t850 * t947 - t847 * t880 + t848 * t913 - t849 * t864 + (t866 * t881 - t867 * t875) * t906) + ((-t807 * t912 + t831 * t916) * r_i_i_C(1) + (-t807 * t916 - t831 * t912) * r_i_i_C(2)) * qJD(5) + (-t881 * qJD(3) - t867 * t939) * t906, -t957, t961 * t782 + (t836 * t880 - t837 * t913 + t854 * t875) * t924 + t926 * (t836 * t864 - t837 * t947 + t845 * t880 - t846 * t913 + t854 * t866 - t875 * t957) (-t782 * t912 + t826 * t916) * r_i_i_C(1) + (-t782 * t916 - t826 * t912) * r_i_i_C(2) + ((-t799 * t916 - t818 * t912) * r_i_i_C(1) + (t799 * t912 - t818 * t916) * r_i_i_C(2)) * qJD(5), 0;];
JaD_transl  = t1;
