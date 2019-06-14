% Calculate vector of inverse dynamics joint torques and base forces with Newton-Euler
% S6RRPRPR6
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
% pkin [11x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,alpha2,d1,d2,d4,d6,theta3]';
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
% Datum: 2019-05-06 14:30
% Revision: 89c353f7eff3bd693eda4e29f35b2761dbc3ada0 (2019-05-03)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauJB = S6RRPRPR6_invdynJB_fixb_snew_vp2(qJ, qJD, qJDD, g, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(6,1),zeros(3,1),zeros(11,1),zeros(7,1),zeros(7,3),zeros(7,6)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRPRPR6_invdynJB_fixb_snew_vp2: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RRPRPR6_invdynJB_fixb_snew_vp2: qJD has to be [6x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [6 1]), ...
  'S6RRPRPR6_invdynJB_fixb_snew_vp2: qJDD has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6RRPRPR6_invdynJB_fixb_snew_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [11 1]), ...
  'S6RRPRPR6_invdynJB_fixb_snew_vp2: pkin has to be [11x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6RRPRPR6_invdynJB_fixb_snew_vp2: m has to be [7x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [7,3]), ...
  'S6RRPRPR6_invdynJB_fixb_snew_vp2: mrSges has to be [7x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [7 6]), ...
  'S6RRPRPR6_invdynJB_fixb_snew_vp2: Ifges has to be [7x6] (double)'); 

%% Symbolic Calculation
% From invdyn_fixb_NewtonEuler_linkframe_tauJB_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-05-06 14:22:03
% EndTime: 2019-05-06 14:22:37
% DurationCPUTime: 22.39s
% Computational Cost: add. (315743->381), mult. (829884->476), div. (0->0), fcn. (645900->12), ass. (0->160)
t1048 = Ifges(5,1) + Ifges(6,2);
t1041 = Ifges(5,4) + Ifges(6,6);
t1040 = Ifges(5,5) - Ifges(6,4);
t1047 = -Ifges(5,2) - Ifges(6,3);
t1039 = Ifges(5,6) - Ifges(6,5);
t1046 = Ifges(5,3) + Ifges(6,1);
t1002 = cos(qJ(2));
t996 = cos(pkin(6));
t1025 = t1002 * t996;
t1004 = qJD(1) ^ 2;
t994 = sin(pkin(6));
t1042 = pkin(8) * t994;
t1000 = sin(qJ(1));
t1003 = cos(qJ(1));
t984 = t1000 * g(1) - g(2) * t1003;
t975 = qJDD(1) * pkin(1) + t1004 * t1042 + t984;
t985 = -g(1) * t1003 - g(2) * t1000;
t976 = -pkin(1) * t1004 + qJDD(1) * t1042 + t985;
t999 = sin(qJ(2));
t1016 = t975 * t1025 - t999 * t976;
t1024 = t1004 * t994 ^ 2;
t1023 = qJD(1) * t1002;
t978 = (qJD(2) * t1023 + qJDD(1) * t999) * t994;
t987 = qJDD(1) * t996 + qJDD(2);
t988 = qJD(1) * t996 + qJD(2);
t904 = t987 * pkin(2) - t978 * qJ(3) + (pkin(2) * t999 * t1024 + (qJ(3) * qJD(1) * t988 - g(3)) * t994) * t1002 + t1016;
t1021 = t1002 ^ 2 * t1024;
t1035 = t996 * t999;
t1036 = t994 * t999;
t943 = -g(3) * t1036 + t1002 * t976 + t975 * t1035;
t1028 = qJD(1) * t999;
t1022 = t994 * t1028;
t972 = pkin(2) * t988 - qJ(3) * t1022;
t979 = (-qJD(2) * t1028 + qJDD(1) * t1002) * t994;
t907 = -pkin(2) * t1021 + qJ(3) * t979 - t972 * t988 + t943;
t1029 = qJD(1) * t994;
t993 = sin(pkin(11));
t995 = cos(pkin(11));
t970 = (t1002 * t993 + t995 * t999) * t1029;
t884 = -0.2e1 * qJD(3) * t970 + t995 * t904 - t993 * t907;
t1001 = cos(qJ(6));
t1043 = cos(qJ(4));
t1020 = t994 * t1023;
t969 = t995 * t1020 - t1022 * t993;
t885 = 0.2e1 * qJD(3) * t969 + t993 * t904 + t995 * t907;
t945 = -pkin(3) * t969 - pkin(9) * t970;
t986 = t988 ^ 2;
t883 = -pkin(3) * t986 + pkin(9) * t987 + t945 * t969 + t885;
t961 = -t996 * g(3) - t994 * t975;
t925 = -t979 * pkin(2) - qJ(3) * t1021 + t972 * t1022 + qJDD(3) + t961;
t949 = -t978 * t993 + t979 * t995;
t950 = t978 * t995 + t979 * t993;
t888 = (-t969 * t988 - t950) * pkin(9) + (t970 * t988 - t949) * pkin(3) + t925;
t998 = sin(qJ(4));
t879 = t1043 * t883 + t998 * t888;
t954 = -t1043 * t988 + t998 * t970;
t955 = t1043 * t970 + t998 * t988;
t926 = pkin(4) * t954 - qJ(5) * t955;
t948 = qJDD(4) - t949;
t968 = qJD(4) - t969;
t967 = t968 ^ 2;
t1012 = -t967 * pkin(4) + t948 * qJ(5) - t954 * t926 + t879;
t921 = t955 * qJD(4) - t1043 * t987 + t998 * t950;
t936 = pkin(5) * t955 - pkin(10) * t968;
t953 = t954 ^ 2;
t873 = -t921 * pkin(5) - t953 * pkin(10) + ((2 * qJD(5)) + t936) * t968 + t1012;
t997 = sin(qJ(6));
t931 = t1001 * t968 + t954 * t997;
t892 = -qJD(6) * t931 + t1001 * t921 - t948 * t997;
t930 = t1001 * t954 - t968 * t997;
t893 = qJD(6) * t930 + t1001 * t948 + t921 * t997;
t952 = qJD(6) + t955;
t905 = -mrSges(7,2) * t952 + mrSges(7,3) * t930;
t906 = mrSges(7,1) * t952 - mrSges(7,3) * t931;
t1013 = -m(7) * t873 + t892 * mrSges(7,1) - t893 * mrSges(7,2) + t930 * t905 - t931 * t906;
t1044 = -2 * qJD(5);
t875 = t1044 * t968 - t1012;
t933 = mrSges(6,1) * t955 + mrSges(6,2) * t968;
t1008 = -m(6) * t875 + t948 * mrSges(6,3) + t968 * t933 - t1013;
t1030 = t1040 * t968 - t1041 * t954 + t1048 * t955;
t1031 = t1039 * t968 + t1041 * t955 + t1047 * t954;
t1037 = t954 * t968;
t878 = t1043 * t888 - t998 * t883;
t876 = -t948 * pkin(4) - t967 * qJ(5) + t955 * t926 + qJDD(5) - t878;
t922 = -t954 * qJD(4) + t1043 * t950 + t998 * t987;
t871 = (t954 * t955 - t948) * pkin(10) + (t922 + t1037) * pkin(5) + t876;
t882 = -t987 * pkin(3) - t986 * pkin(9) + t970 * t945 - t884;
t1007 = (-t922 + t1037) * qJ(5) + t882 + (pkin(4) * t968 + t1044) * t955;
t874 = -t953 * pkin(5) - t955 * t936 + (pkin(4) + pkin(10)) * t921 + t1007;
t869 = t1001 * t871 - t874 * t997;
t898 = -mrSges(7,1) * t930 + mrSges(7,2) * t931;
t920 = qJDD(6) + t922;
t866 = m(7) * t869 + mrSges(7,1) * t920 - t893 * mrSges(7,3) - t898 * t931 + t905 * t952;
t870 = t1001 * t874 + t871 * t997;
t867 = m(7) * t870 - mrSges(7,2) * t920 + t892 * mrSges(7,3) + t898 * t930 - t906 * t952;
t857 = t1001 * t866 + t997 * t867;
t928 = -mrSges(6,2) * t954 - mrSges(6,3) * t955;
t1010 = -m(6) * t876 - t922 * mrSges(6,1) - t955 * t928 - t857;
t932 = mrSges(6,1) * t954 - mrSges(6,3) * t968;
t855 = t948 * mrSges(6,2) + t968 * t932 - t1010;
t894 = Ifges(7,5) * t931 + Ifges(7,6) * t930 + Ifges(7,3) * t952;
t896 = Ifges(7,1) * t931 + Ifges(7,4) * t930 + Ifges(7,5) * t952;
t858 = -mrSges(7,1) * t873 + mrSges(7,3) * t870 + Ifges(7,4) * t893 + Ifges(7,2) * t892 + Ifges(7,6) * t920 - t894 * t931 + t896 * t952;
t895 = Ifges(7,4) * t931 + Ifges(7,2) * t930 + Ifges(7,6) * t952;
t859 = mrSges(7,2) * t873 - mrSges(7,3) * t869 + Ifges(7,1) * t893 + Ifges(7,4) * t892 + Ifges(7,5) * t920 + t894 * t930 - t895 * t952;
t1045 = t1030 * t954 + t1031 * t955 + t1046 * t948 - t1039 * t921 + t1040 * t922 + mrSges(5,1) * t878 - mrSges(5,2) * t879 + mrSges(6,2) * t876 - mrSges(6,3) * t875 - pkin(4) * t855 - pkin(10) * t857 + qJ(5) * (-t921 * mrSges(6,1) - t954 * t928 + t1008) + t1001 * t859 - t997 * t858;
t927 = mrSges(5,1) * t954 + mrSges(5,2) * t955;
t934 = -mrSges(5,2) * t968 - mrSges(5,3) * t954;
t854 = m(5) * t878 - t922 * mrSges(5,3) - t955 * t927 + (-t932 + t934) * t968 + (mrSges(5,1) - mrSges(6,2)) * t948 + t1010;
t935 = mrSges(5,1) * t968 - mrSges(5,3) * t955;
t862 = m(5) * t879 - t948 * mrSges(5,2) - t968 * t935 + (-t927 - t928) * t954 + (-mrSges(5,3) - mrSges(6,1)) * t921 + t1008;
t1018 = t1043 * t862 - t854 * t998;
t944 = -mrSges(4,1) * t969 + mrSges(4,2) * t970;
t957 = mrSges(4,1) * t988 - mrSges(4,3) * t970;
t846 = m(4) * t885 - mrSges(4,2) * t987 + mrSges(4,3) * t949 + t944 * t969 - t957 * t988 + t1018;
t1033 = t1001 * t867 - t997 * t866;
t877 = t921 * pkin(4) + t1007;
t1014 = -m(6) * t877 + t921 * mrSges(6,2) + t954 * t932 - t1033;
t1006 = -m(5) * t882 - t954 * t934 - t921 * mrSges(5,1) + (t933 - t935) * t955 + (-mrSges(5,2) + mrSges(6,3)) * t922 + t1014;
t956 = -mrSges(4,2) * t988 + mrSges(4,3) * t969;
t852 = m(4) * t884 + t987 * mrSges(4,1) - t950 * mrSges(4,3) - t970 * t944 + t988 * t956 + t1006;
t843 = t993 * t846 + t995 * t852;
t1026 = t1002 * t994;
t942 = -g(3) * t1026 + t1016;
t974 = -mrSges(3,2) * t988 + mrSges(3,3) * t1020;
t977 = (-mrSges(3,1) * t1002 + mrSges(3,2) * t999) * t1029;
t841 = m(3) * t942 + mrSges(3,1) * t987 - mrSges(3,3) * t978 - t1022 * t977 + t974 * t988 + t843;
t1019 = t995 * t846 - t852 * t993;
t973 = mrSges(3,1) * t988 - mrSges(3,3) * t1022;
t842 = m(3) * t943 - mrSges(3,2) * t987 + mrSges(3,3) * t979 + t1020 * t977 - t973 * t988 + t1019;
t850 = t1043 * t854 + t998 * t862;
t849 = m(4) * t925 - t949 * mrSges(4,1) + t950 * mrSges(4,2) - t969 * t956 + t970 * t957 + t850;
t848 = m(3) * t961 - t979 * mrSges(3,1) + t978 * mrSges(3,2) + (-t1002 * t974 + t973 * t999) * t1029 + t849;
t827 = t841 * t1025 + t842 * t1035 - t848 * t994;
t824 = m(2) * t984 + qJDD(1) * mrSges(2,1) - mrSges(2,2) * t1004 + t827;
t832 = t1002 * t842 - t841 * t999;
t830 = m(2) * t985 - mrSges(2,1) * t1004 - qJDD(1) * mrSges(2,2) + t832;
t1034 = t1000 * t830 + t1003 * t824;
t1032 = t1039 * t954 - t1040 * t955 - t1046 * t968;
t826 = t841 * t1026 + t842 * t1036 + t996 * t848;
t1015 = -t1000 * t824 + t1003 * t830;
t856 = -t922 * mrSges(6,3) - t955 * t933 - t1014;
t834 = -mrSges(5,1) * t882 - mrSges(6,1) * t875 + mrSges(6,2) * t877 + mrSges(5,3) * t879 - pkin(4) * t856 - pkin(5) * t1013 - pkin(10) * t1033 - t1001 * t858 + t1030 * t968 + t1032 * t955 + t1039 * t948 + t1041 * t922 + t1047 * t921 - t997 * t859;
t1009 = mrSges(7,1) * t869 - mrSges(7,2) * t870 + Ifges(7,5) * t893 + Ifges(7,6) * t892 + Ifges(7,3) * t920 + t931 * t895 - t930 * t896;
t835 = mrSges(6,1) * t876 + mrSges(5,2) * t882 - mrSges(5,3) * t878 - mrSges(6,3) * t877 + pkin(5) * t857 - qJ(5) * t856 - t1031 * t968 + t1032 * t954 + t1040 * t948 - t1041 * t921 + t1048 * t922 + t1009;
t938 = Ifges(4,5) * t970 + Ifges(4,6) * t969 + Ifges(4,3) * t988;
t939 = Ifges(4,4) * t970 + Ifges(4,2) * t969 + Ifges(4,6) * t988;
t822 = mrSges(4,2) * t925 - mrSges(4,3) * t884 + Ifges(4,1) * t950 + Ifges(4,4) * t949 + Ifges(4,5) * t987 - pkin(9) * t850 + t1043 * t835 - t998 * t834 + t969 * t938 - t988 * t939;
t940 = Ifges(4,1) * t970 + Ifges(4,4) * t969 + Ifges(4,5) * t988;
t833 = -mrSges(4,1) * t925 + mrSges(4,3) * t885 + Ifges(4,4) * t950 + Ifges(4,2) * t949 + Ifges(4,6) * t987 - pkin(3) * t850 - t970 * t938 + t988 * t940 - t1045;
t958 = Ifges(3,3) * t988 + (Ifges(3,5) * t999 + Ifges(3,6) * t1002) * t1029;
t960 = Ifges(3,5) * t988 + (Ifges(3,1) * t999 + Ifges(3,4) * t1002) * t1029;
t817 = -mrSges(3,1) * t961 + mrSges(3,3) * t943 + Ifges(3,4) * t978 + Ifges(3,2) * t979 + Ifges(3,6) * t987 - pkin(2) * t849 + qJ(3) * t1019 - t1022 * t958 + t993 * t822 + t995 * t833 + t988 * t960;
t959 = Ifges(3,6) * t988 + (Ifges(3,4) * t999 + Ifges(3,2) * t1002) * t1029;
t819 = mrSges(3,2) * t961 - mrSges(3,3) * t942 + Ifges(3,1) * t978 + Ifges(3,4) * t979 + Ifges(3,5) * t987 - qJ(3) * t843 + t1020 * t958 + t822 * t995 - t833 * t993 - t959 * t988;
t821 = Ifges(3,5) * t978 + Ifges(3,6) * t979 + mrSges(3,1) * t942 - mrSges(3,2) * t943 + Ifges(4,5) * t950 + Ifges(4,6) * t949 + t970 * t939 - t969 * t940 + mrSges(4,1) * t884 - mrSges(4,2) * t885 + t998 * t835 + t1043 * t834 + pkin(3) * t1006 + pkin(9) * t1018 + pkin(2) * t843 + (Ifges(3,3) + Ifges(4,3)) * t987 + (-t1002 * t960 + t959 * t999) * t1029;
t1011 = mrSges(2,1) * t984 - mrSges(2,2) * t985 + Ifges(2,3) * qJDD(1) + pkin(1) * t827 + t817 * t1026 + t819 * t1036 + t832 * t1042 + t996 * t821;
t815 = -mrSges(2,2) * g(3) - mrSges(2,3) * t984 + Ifges(2,5) * qJDD(1) - t1004 * Ifges(2,6) + t1002 * t819 - t999 * t817 + (-t826 * t994 - t827 * t996) * pkin(8);
t814 = mrSges(2,1) * g(3) + mrSges(2,3) * t985 + t1004 * Ifges(2,5) + Ifges(2,6) * qJDD(1) - pkin(1) * t826 - t994 * t821 + (pkin(8) * t832 + t1002 * t817 + t819 * t999) * t996;
t1 = [-m(1) * g(1) + t1015; -m(1) * g(2) + t1034; (-m(1) - m(2)) * g(3) + t826; -mrSges(1,2) * g(3) + mrSges(1,3) * g(2) - pkin(7) * t1034 - t1000 * t814 + t1003 * t815; mrSges(1,1) * g(3) - mrSges(1,3) * g(1) + pkin(7) * t1015 + t1000 * t815 + t1003 * t814; -mrSges(1,1) * g(2) + mrSges(1,2) * g(1) + t1011; t1011; t821; t849; t1045; t855; t1009;];
tauJB  = t1;
