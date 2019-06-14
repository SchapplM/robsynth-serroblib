% Calculate vector of inverse dynamics joint torques and base forces with Newton-Euler
% S6PRPRRR7
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% qJD [6x1]
%   Generalized joint velocities
% qJDD [6x1]
%   Generalized joint accelerations
% g [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [14x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,alpha2,alpha3,alpha4,d2,d4,d5,d6,theta1,theta3]';
% m_mdh [7x1]
%   mass of all robot links (including the base)
% mrSges [7x3]
%  first moment of all robot links (mass times center of mass in body frames)
%  rows: links of the robot (starting with base)
%  columns: x-, y-, z-coordinates
% Ifges [7x6]
%   inertia of all robot links about their respective body frame origins, in body frames
%   rows: links of the robot (starting with base)
%   columns: xx, yy, zz, xy, xz, yz (see inertial_parameters_convert_par1_par2.m)
% 
% Output:
% tauJB [(6+6)x1]
%   joint torques and base forces of inverse dynamics (contains inertial, gravitational coriolis and centrifugal forces)

% Quelle: HybrDyn-Toolbox
% Datum: 2019-05-05 02:02
% Revision: 89c353f7eff3bd693eda4e29f35b2761dbc3ada0 (2019-05-03)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauJB = S6PRPRRR7_invdynJB_fixb_snew_vp2(qJ, qJD, qJDD, g, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(6,1),zeros(3,1),zeros(14,1),zeros(7,1),zeros(7,3),zeros(7,6)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6PRPRRR7_invdynJB_fixb_snew_vp2: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6PRPRRR7_invdynJB_fixb_snew_vp2: qJD has to be [6x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [6 1]), ...
  'S6PRPRRR7_invdynJB_fixb_snew_vp2: qJDD has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6PRPRRR7_invdynJB_fixb_snew_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [14 1]), ...
  'S6PRPRRR7_invdynJB_fixb_snew_vp2: pkin has to be [14x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6PRPRRR7_invdynJB_fixb_snew_vp2: m has to be [7x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [7,3]), ...
  'S6PRPRRR7_invdynJB_fixb_snew_vp2: mrSges has to be [7x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [7 6]), ...
  'S6PRPRRR7_invdynJB_fixb_snew_vp2: Ifges has to be [7x6] (double)'); 

%% Symbolic Calculation
% From invdyn_fixb_NewtonEuler_linkframe_tauJB_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-05-05 01:54:00
% EndTime: 2019-05-05 01:55:38
% DurationCPUTime: 100.73s
% Computational Cost: add. (1663980->356), mult. (4698873->499), div. (0->0), fcn. (4009726->18), ass. (0->180)
t1001 = cos(qJ(4));
t988 = sin(pkin(7));
t992 = cos(pkin(8));
t1041 = t988 * t992;
t987 = sin(pkin(8));
t990 = cos(pkin(14));
t993 = cos(pkin(7));
t1017 = t1041 * t990 + t987 * t993;
t985 = sin(pkin(14));
t1045 = t985 * t988;
t997 = sin(qJ(4));
t1009 = t1001 * t1017 - t1045 * t997;
t956 = t1009 * qJD(2);
t1038 = t992 * t997;
t1043 = t987 * t997;
t1007 = t993 * t1043 + (t1001 * t985 + t1038 * t990) * t988;
t957 = t1007 * qJD(2);
t940 = -t957 * qJD(4) + qJDD(2) * t1009;
t1035 = qJD(2) * t988;
t1027 = qJD(3) * t1035;
t1039 = t990 * t993;
t1042 = t988 * t990;
t1003 = qJD(2) ^ 2;
t1046 = qJ(3) * t988;
t1002 = cos(qJ(2));
t994 = cos(pkin(6));
t1032 = t1002 * t994;
t989 = sin(pkin(6));
t1033 = t1002 * t989;
t986 = sin(pkin(13));
t991 = cos(pkin(13));
t979 = g(1) * t986 - g(2) * t991;
t980 = -g(1) * t991 - g(2) * t986;
t984 = -g(3) + qJDD(1);
t998 = sin(qJ(2));
t952 = t1032 * t979 + t1033 * t984 - t980 * t998;
t950 = qJDD(2) * pkin(2) + t1003 * t1046 + t952;
t969 = -t979 * t989 + t984 * t994;
t1028 = -0.2e1 * t1027 * t985 + t1039 * t950 + t1042 * t969;
t1037 = t994 * t998;
t1040 = t989 * t998;
t953 = t1002 * t980 + t1037 * t979 + t1040 * t984;
t951 = -pkin(2) * t1003 + qJDD(2) * t1046 + t953;
t1048 = pkin(10) * t985;
t1020 = -pkin(3) * t990 - t1048 * t987;
t966 = t1020 * t1035;
t1015 = pkin(10) * t1017;
t967 = qJD(2) * t1015;
t909 = (pkin(3) * qJDD(2) + qJD(2) * t967) * t993 + (-t951 + (-pkin(10) * qJDD(2) * t992 - qJD(2) * t966) * t988) * t985 + t1028;
t1044 = t985 * t993;
t925 = t950 * t1044 + t969 * t1045 + (0.2e1 * t1027 + t951) * t990;
t971 = (pkin(3) * t993 - t1041 * t1048) * qJD(2);
t910 = (t1042 * t966 - t971 * t993) * qJD(2) + qJDD(2) * t1015 + t925;
t1030 = t969 * t993 + qJDD(3);
t923 = (-t950 + t1020 * qJDD(2) + (-t967 * t990 + t971 * t985) * qJD(2)) * t988 + t1030;
t895 = (t909 * t992 + t923 * t987) * t1001 - t997 * t910;
t1047 = Ifges(4,3) * t993;
t1019 = mrSges(4,1) * t993 - mrSges(4,3) * t1045;
t1000 = cos(qJ(5));
t896 = t1001 * t910 + t1038 * t909 + t1043 * t923;
t939 = -pkin(4) * t956 - pkin(11) * t957;
t1016 = -t1042 * t987 + t992 * t993;
t968 = qJD(2) * t1016 + qJD(4);
t964 = t968 ^ 2;
t965 = qJDD(2) * t1016 + qJDD(4);
t892 = -pkin(4) * t964 + pkin(11) * t965 + t939 * t956 + t896;
t897 = -t987 * t909 + t923 * t992;
t941 = t956 * qJD(4) + qJDD(2) * t1007;
t894 = (-t956 * t968 - t941) * pkin(11) + (t957 * t968 - t940) * pkin(4) + t897;
t996 = sin(qJ(5));
t888 = t1000 * t892 + t894 * t996;
t943 = t1000 * t968 - t957 * t996;
t944 = t1000 * t957 + t968 * t996;
t927 = -pkin(5) * t943 - pkin(12) * t944;
t937 = qJDD(5) - t940;
t955 = qJD(5) - t956;
t954 = t955 ^ 2;
t886 = -pkin(5) * t954 + pkin(12) * t937 + t927 * t943 + t888;
t891 = -t965 * pkin(4) - t964 * pkin(11) + t939 * t957 - t895;
t919 = -qJD(5) * t944 + t1000 * t965 - t941 * t996;
t920 = qJD(5) * t943 + t1000 * t941 + t965 * t996;
t889 = (-t943 * t955 - t920) * pkin(12) + (t944 * t955 - t919) * pkin(5) + t891;
t995 = sin(qJ(6));
t999 = cos(qJ(6));
t882 = -t886 * t995 + t889 * t999;
t929 = -t944 * t995 + t955 * t999;
t900 = qJD(6) * t929 + t920 * t999 + t937 * t995;
t930 = t944 * t999 + t955 * t995;
t908 = -mrSges(7,1) * t929 + mrSges(7,2) * t930;
t942 = qJD(6) - t943;
t911 = -mrSges(7,2) * t942 + mrSges(7,3) * t929;
t917 = qJDD(6) - t919;
t880 = m(7) * t882 + mrSges(7,1) * t917 - mrSges(7,3) * t900 - t908 * t930 + t911 * t942;
t883 = t886 * t999 + t889 * t995;
t899 = -qJD(6) * t930 - t920 * t995 + t937 * t999;
t912 = mrSges(7,1) * t942 - mrSges(7,3) * t930;
t881 = m(7) * t883 - mrSges(7,2) * t917 + mrSges(7,3) * t899 + t908 * t929 - t912 * t942;
t873 = t880 * t999 + t881 * t995;
t931 = -mrSges(6,2) * t955 + mrSges(6,3) * t943;
t932 = mrSges(6,1) * t955 - mrSges(6,3) * t944;
t1006 = -m(6) * t891 + mrSges(6,1) * t919 - mrSges(6,2) * t920 + t931 * t943 - t932 * t944 - t873;
t938 = -mrSges(5,1) * t956 + mrSges(5,2) * t957;
t945 = -mrSges(5,2) * t968 + mrSges(5,3) * t956;
t869 = m(5) * t895 + mrSges(5,1) * t965 - mrSges(5,3) * t941 - t938 * t957 + t945 * t968 + t1006;
t1034 = t1001 * t869;
t874 = -t880 * t995 + t881 * t999;
t926 = -mrSges(6,1) * t943 + mrSges(6,2) * t944;
t872 = m(6) * t888 - mrSges(6,2) * t937 + mrSges(6,3) * t919 + t926 * t943 - t932 * t955 + t874;
t887 = t1000 * t894 - t892 * t996;
t885 = -pkin(5) * t937 - pkin(12) * t954 + t927 * t944 - t887;
t884 = -m(7) * t885 + mrSges(7,1) * t899 - mrSges(7,2) * t900 + t911 * t929 - t912 * t930;
t878 = m(6) * t887 + mrSges(6,1) * t937 - mrSges(6,3) * t920 - t926 * t944 + t931 * t955 + t884;
t1025 = t1000 * t872 - t878 * t996;
t946 = mrSges(5,1) * t968 - mrSges(5,3) * t957;
t863 = m(5) * t896 - mrSges(5,2) * t965 + mrSges(5,3) * t940 + t938 * t956 - t946 * t968 + t1025;
t866 = t1000 * t878 + t872 * t996;
t865 = m(5) * t897 - mrSges(5,1) * t940 + mrSges(5,2) * t941 - t945 * t956 + t946 * t957 + t866;
t852 = t1034 * t992 + t1038 * t863 - t987 * t865;
t924 = -t985 * t951 + t1028;
t1023 = -mrSges(4,1) * t990 + mrSges(4,2) * t985;
t970 = t1023 * t1035;
t1018 = -mrSges(4,2) * t993 + mrSges(4,3) * t1042;
t975 = t1018 * qJD(2);
t848 = m(4) * t924 + t1019 * qJDD(2) + (-t1045 * t970 + t975 * t993) * qJD(2) + t852;
t851 = t1034 * t987 + t1043 * t863 + t865 * t992;
t936 = -t988 * t950 + t1030;
t974 = t1019 * qJD(2);
t850 = m(4) * t936 + (t1023 * qJDD(2) + (t974 * t985 - t975 * t990) * qJD(2)) * t988 + t851;
t857 = t1001 * t863 - t997 * t869;
t856 = m(4) * t925 + t1018 * qJDD(2) + (t1042 * t970 - t974 * t993) * qJD(2) + t857;
t837 = t1039 * t848 + t1044 * t856 - t850 * t988;
t833 = m(3) * t952 + qJDD(2) * mrSges(3,1) - mrSges(3,2) * t1003 + t837;
t836 = t1042 * t848 + t1045 * t856 + t850 * t993;
t835 = m(3) * t969 + t836;
t842 = -t848 * t985 + t856 * t990;
t841 = m(3) * t953 - mrSges(3,1) * t1003 - qJDD(2) * mrSges(3,2) + t842;
t823 = t1032 * t833 + t1037 * t841 - t835 * t989;
t821 = m(2) * t979 + t823;
t830 = t1002 * t841 - t833 * t998;
t829 = m(2) * t980 + t830;
t1036 = t821 * t991 + t829 * t986;
t822 = t1033 * t833 + t1040 * t841 + t835 * t994;
t1026 = -t821 * t986 + t829 * t991;
t1024 = m(2) * t984 + t822;
t1022 = Ifges(4,5) * t985 + Ifges(4,6) * t990;
t1010 = Ifges(4,6) * t993 + (Ifges(4,4) * t985 + Ifges(4,2) * t990) * t988;
t901 = Ifges(7,5) * t930 + Ifges(7,6) * t929 + Ifges(7,3) * t942;
t903 = Ifges(7,1) * t930 + Ifges(7,4) * t929 + Ifges(7,5) * t942;
t875 = -mrSges(7,1) * t885 + mrSges(7,3) * t883 + Ifges(7,4) * t900 + Ifges(7,2) * t899 + Ifges(7,6) * t917 - t901 * t930 + t903 * t942;
t902 = Ifges(7,4) * t930 + Ifges(7,2) * t929 + Ifges(7,6) * t942;
t876 = mrSges(7,2) * t885 - mrSges(7,3) * t882 + Ifges(7,1) * t900 + Ifges(7,4) * t899 + Ifges(7,5) * t917 + t901 * t929 - t902 * t942;
t913 = Ifges(6,5) * t944 + Ifges(6,6) * t943 + Ifges(6,3) * t955;
t914 = Ifges(6,4) * t944 + Ifges(6,2) * t943 + Ifges(6,6) * t955;
t858 = mrSges(6,2) * t891 - mrSges(6,3) * t887 + Ifges(6,1) * t920 + Ifges(6,4) * t919 + Ifges(6,5) * t937 - pkin(12) * t873 - t875 * t995 + t876 * t999 + t913 * t943 - t914 * t955;
t1005 = mrSges(7,1) * t882 - mrSges(7,2) * t883 + Ifges(7,5) * t900 + Ifges(7,6) * t899 + Ifges(7,3) * t917 + t902 * t930 - t903 * t929;
t915 = Ifges(6,1) * t944 + Ifges(6,4) * t943 + Ifges(6,5) * t955;
t859 = -mrSges(6,1) * t891 + mrSges(6,3) * t888 + Ifges(6,4) * t920 + Ifges(6,2) * t919 + Ifges(6,6) * t937 - pkin(5) * t873 - t913 * t944 + t915 * t955 - t1005;
t933 = Ifges(5,5) * t957 + Ifges(5,6) * t956 + Ifges(5,3) * t968;
t934 = Ifges(5,4) * t957 + Ifges(5,2) * t956 + Ifges(5,6) * t968;
t844 = mrSges(5,2) * t897 - mrSges(5,3) * t895 + Ifges(5,1) * t941 + Ifges(5,4) * t940 + Ifges(5,5) * t965 - pkin(11) * t866 + t1000 * t858 - t859 * t996 + t933 * t956 - t934 * t968;
t1004 = mrSges(6,1) * t887 - mrSges(6,2) * t888 + Ifges(6,5) * t920 + Ifges(6,6) * t919 + Ifges(6,3) * t937 + pkin(5) * t884 + pkin(12) * t874 + t875 * t999 + t876 * t995 + t914 * t944 - t915 * t943;
t935 = Ifges(5,1) * t957 + Ifges(5,4) * t956 + Ifges(5,5) * t968;
t845 = -mrSges(5,1) * t897 + mrSges(5,3) * t896 + Ifges(5,4) * t941 + Ifges(5,2) * t940 + Ifges(5,6) * t965 - pkin(4) * t866 - t933 * t957 + t935 * t968 - t1004;
t1013 = pkin(10) * t857 + t1001 * t845 + t844 * t997;
t843 = mrSges(5,1) * t895 - mrSges(5,2) * t896 + Ifges(5,5) * t941 + Ifges(5,6) * t940 + Ifges(5,3) * t965 + pkin(4) * t1006 + pkin(11) * t1025 + t1000 * t859 + t996 * t858 + t957 * t934 - t956 * t935;
t960 = (t1022 * t988 + t1047) * qJD(2);
t1011 = Ifges(4,5) * t993 + (Ifges(4,1) * t985 + Ifges(4,4) * t990) * t988;
t962 = t1011 * qJD(2);
t825 = -mrSges(4,1) * t936 + mrSges(4,3) * t925 - pkin(3) * t851 - t987 * t843 + (-t1045 * t960 + t962 * t993) * qJD(2) + t1013 * t992 + t1010 * qJDD(2);
t961 = t1010 * qJD(2);
t826 = mrSges(4,2) * t936 - mrSges(4,3) * t924 + t1001 * t844 - t997 * t845 + (t1042 * t960 - t961 * t993) * qJD(2) + (-t851 * t987 - t852 * t992) * pkin(10) + t1011 * qJDD(2);
t1012 = qJ(3) * t842 + t825 * t990 + t826 * t985;
t824 = qJDD(2) * t1047 + mrSges(4,1) * t924 - mrSges(4,2) * t925 + pkin(3) * t852 + t992 * t843 + t1013 * t987 + (t1022 * qJDD(2) + (t961 * t985 - t962 * t990) * qJD(2)) * t988;
t818 = -mrSges(3,1) * t969 + mrSges(3,3) * t953 + t1003 * Ifges(3,5) + Ifges(3,6) * qJDD(2) - pkin(2) * t836 + t1012 * t993 - t988 * t824;
t819 = mrSges(3,2) * t969 - mrSges(3,3) * t952 + Ifges(3,5) * qJDD(2) - t1003 * Ifges(3,6) - t985 * t825 + t990 * t826 + (-t836 * t988 - t837 * t993) * qJ(3);
t1014 = pkin(9) * t830 + t1002 * t818 + t819 * t998;
t817 = mrSges(3,1) * t952 - mrSges(3,2) * t953 + Ifges(3,3) * qJDD(2) + pkin(2) * t837 + t1012 * t988 + t993 * t824;
t816 = mrSges(2,2) * t984 - mrSges(2,3) * t979 + t1002 * t819 - t998 * t818 + (-t822 * t989 - t823 * t994) * pkin(9);
t815 = -mrSges(2,1) * t984 + mrSges(2,3) * t980 - pkin(1) * t822 + t1014 * t994 - t989 * t817;
t1 = [-m(1) * g(1) + t1026; -m(1) * g(2) + t1036; -m(1) * g(3) + t1024; -mrSges(1,2) * g(3) + mrSges(1,3) * g(2) - qJ(1) * t1036 - t815 * t986 + t816 * t991; mrSges(1,1) * g(3) - mrSges(1,3) * g(1) + qJ(1) * t1026 + t991 * t815 + t986 * t816; -mrSges(1,1) * g(2) + mrSges(2,1) * t979 + mrSges(1,2) * g(1) - mrSges(2,2) * t980 + pkin(1) * t823 + t1014 * t989 + t994 * t817; t1024; t817; t850; t843; t1004; t1005;];
tauJB  = t1;
