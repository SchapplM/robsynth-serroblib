% Calculate vector of inverse dynamics joint torques and base forces with Newton-Euler
% S6RRPPRR4
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
%   pkin=[a2,a3,a4,a5,a6,alpha2,d1,d2,d5,d6,theta3]';
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
% Datum: 2019-05-06 10:27
% Revision: 89c353f7eff3bd693eda4e29f35b2761dbc3ada0 (2019-05-03)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauJB = S6RRPPRR4_invdynJB_fixb_snew_vp2(qJ, qJD, qJDD, g, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(6,1),zeros(3,1),zeros(11,1),zeros(7,1),zeros(7,3),zeros(7,6)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRPPRR4_invdynJB_fixb_snew_vp2: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RRPPRR4_invdynJB_fixb_snew_vp2: qJD has to be [6x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [6 1]), ...
  'S6RRPPRR4_invdynJB_fixb_snew_vp2: qJDD has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6RRPPRR4_invdynJB_fixb_snew_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [11 1]), ...
  'S6RRPPRR4_invdynJB_fixb_snew_vp2: pkin has to be [11x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6RRPPRR4_invdynJB_fixb_snew_vp2: m has to be [7x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [7,3]), ...
  'S6RRPPRR4_invdynJB_fixb_snew_vp2: mrSges has to be [7x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [7 6]), ...
  'S6RRPPRR4_invdynJB_fixb_snew_vp2: Ifges has to be [7x6] (double)'); 

%% Symbolic Calculation
% From invdyn_fixb_NewtonEuler_linkframe_tauJB_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-05-06 10:20:57
% EndTime: 2019-05-06 10:21:17
% DurationCPUTime: 19.79s
% Computational Cost: add. (281205->380), mult. (756345->475), div. (0->0), fcn. (578517->12), ass. (0->163)
t1079 = -2 * qJD(4);
t1078 = Ifges(4,1) + Ifges(5,2);
t1077 = -Ifges(5,1) - Ifges(4,3);
t1073 = Ifges(4,5) - Ifges(5,4);
t1076 = -Ifges(4,2) - Ifges(5,3);
t1072 = Ifges(4,6) - Ifges(5,5);
t1071 = -Ifges(5,6) - Ifges(4,4);
t1023 = cos(pkin(6));
t1016 = qJD(1) * t1023 + qJD(2);
t1014 = t1016 ^ 2;
t1015 = qJDD(1) * t1023 + qJDD(2);
t1021 = sin(pkin(11));
t1069 = cos(pkin(11));
t1022 = sin(pkin(6));
t1026 = sin(qJ(2));
t1030 = cos(qJ(2));
t1059 = qJD(1) * qJD(2);
t1005 = (qJDD(1) * t1026 + t1030 * t1059) * t1022;
t1027 = sin(qJ(1));
t1031 = cos(qJ(1));
t1012 = t1027 * g(1) - g(2) * t1031;
t1032 = qJD(1) ^ 2;
t1070 = pkin(8) * t1022;
t1002 = qJDD(1) * pkin(1) + t1032 * t1070 + t1012;
t1013 = -g(1) * t1031 - g(2) * t1027;
t1003 = -pkin(1) * t1032 + qJDD(1) * t1070 + t1013;
t1052 = t1023 * t1030;
t1044 = t1002 * t1052 - t1026 * t1003;
t1056 = t1022 ^ 2 * t1032;
t926 = t1015 * pkin(2) - t1005 * qJ(3) + (pkin(2) * t1026 * t1056 + (qJ(3) * qJD(1) * t1016 - g(3)) * t1022) * t1030 + t1044;
t1006 = (qJDD(1) * t1030 - t1026 * t1059) * t1022;
t1049 = t1030 ^ 2 * t1056;
t1053 = t1023 * t1026;
t1055 = t1022 * t1026;
t960 = -g(3) * t1055 + t1002 * t1053 + t1030 * t1003;
t1051 = qJD(1) * t1055;
t999 = pkin(2) * t1016 - qJ(3) * t1051;
t929 = -pkin(2) * t1049 + qJ(3) * t1006 - t1016 * t999 + t960;
t1067 = t1021 * t926 + t1069 * t929;
t1054 = t1022 * t1030;
t1050 = qJD(1) * t1054;
t995 = t1021 * t1051 - t1050 * t1069;
t1058 = qJD(1) * t1022;
t996 = (t1021 * t1030 + t1026 * t1069) * t1058;
t961 = pkin(3) * t995 - qJ(4) * t996;
t1075 = t1014 * pkin(3) - t1015 * qJ(4) + t1016 * t1079 + t995 * t961 - t1067;
t912 = -0.2e1 * qJD(3) * t996 - t1021 * t929 + t1069 * t926;
t1074 = mrSges(4,1) - mrSges(5,2);
t1001 = -mrSges(3,2) * t1016 + mrSges(3,3) * t1050;
t1004 = (-mrSges(3,1) * t1030 + mrSges(3,2) * t1026) * t1058;
t1025 = sin(qJ(5));
t1029 = cos(qJ(5));
t1024 = sin(qJ(6));
t1028 = cos(qJ(6));
t1060 = t1016 * t995;
t909 = -t1015 * pkin(3) - t1014 * qJ(4) + t996 * t961 + qJDD(4) - t912;
t970 = t1005 * t1069 + t1021 * t1006;
t903 = (t995 * t996 - t1015) * pkin(9) + (t970 + t1060) * pkin(4) + t909;
t985 = -t1023 * g(3) - t1022 * t1002;
t941 = -t1006 * pkin(2) - qJ(3) * t1049 + t999 * t1051 + qJDD(3) + t985;
t1033 = (-t970 + t1060) * qJ(4) + t941 + (pkin(3) * t1016 + t1079) * t996;
t969 = t1021 * t1005 - t1006 * t1069;
t980 = pkin(4) * t996 - pkin(9) * t1016;
t994 = t995 ^ 2;
t907 = -t994 * pkin(4) - t996 * t980 + (pkin(3) + pkin(9)) * t969 + t1033;
t900 = t1025 * t903 + t1029 * t907;
t974 = -t1016 * t1025 + t1029 * t995;
t975 = t1016 * t1029 + t1025 * t995;
t943 = -pkin(5) * t974 - pkin(10) * t975;
t968 = qJDD(5) + t970;
t993 = qJD(5) + t996;
t992 = t993 ^ 2;
t897 = -pkin(5) * t992 + pkin(10) * t968 + t943 * t974 + t900;
t1062 = qJD(3) * t995;
t988 = -0.2e1 * t1062;
t905 = -t969 * pkin(4) - t994 * pkin(9) + t1016 * t980 - t1075 + t988;
t938 = -qJD(5) * t975 - t1015 * t1025 + t1029 * t969;
t939 = qJD(5) * t974 + t1015 * t1029 + t1025 * t969;
t901 = (-t974 * t993 - t939) * pkin(10) + (t975 * t993 - t938) * pkin(5) + t905;
t894 = -t1024 * t897 + t1028 * t901;
t947 = -t1024 * t975 + t1028 * t993;
t916 = qJD(6) * t947 + t1024 * t968 + t1028 * t939;
t948 = t1024 * t993 + t1028 * t975;
t922 = -mrSges(7,1) * t947 + mrSges(7,2) * t948;
t973 = qJD(6) - t974;
t927 = -mrSges(7,2) * t973 + mrSges(7,3) * t947;
t937 = qJDD(6) - t938;
t891 = m(7) * t894 + mrSges(7,1) * t937 - t916 * mrSges(7,3) - t922 * t948 + t927 * t973;
t895 = t1024 * t901 + t1028 * t897;
t915 = -qJD(6) * t948 - t1024 * t939 + t1028 * t968;
t928 = mrSges(7,1) * t973 - mrSges(7,3) * t948;
t892 = m(7) * t895 - mrSges(7,2) * t937 + t915 * mrSges(7,3) + t922 * t947 - t928 * t973;
t1047 = -t1024 * t891 + t1028 * t892;
t942 = -mrSges(6,1) * t974 + mrSges(6,2) * t975;
t950 = mrSges(6,1) * t993 - mrSges(6,3) * t975;
t880 = m(6) * t900 - mrSges(6,2) * t968 + mrSges(6,3) * t938 + t942 * t974 - t950 * t993 + t1047;
t899 = -t1025 * t907 + t1029 * t903;
t896 = -pkin(5) * t968 - pkin(10) * t992 + t943 * t975 - t899;
t1040 = -m(7) * t896 + t915 * mrSges(7,1) - t916 * mrSges(7,2) + t947 * t927 - t928 * t948;
t949 = -mrSges(6,2) * t993 + mrSges(6,3) * t974;
t887 = m(6) * t899 + mrSges(6,1) * t968 - mrSges(6,3) * t939 - t942 * t975 + t949 * t993 + t1040;
t875 = t1025 * t880 + t1029 * t887;
t963 = -mrSges(5,2) * t995 - mrSges(5,3) * t996;
t1039 = -m(5) * t909 - t970 * mrSges(5,1) - t996 * t963 - t875;
t978 = mrSges(5,1) * t995 - mrSges(5,3) * t1016;
t1063 = -mrSges(4,2) * t1016 - mrSges(4,3) * t995 - t978;
t962 = mrSges(4,1) * t995 + mrSges(4,2) * t996;
t869 = m(4) * t912 - t970 * mrSges(4,3) + t1015 * t1074 + t1063 * t1016 - t996 * t962 + t1039;
t882 = t1024 * t892 + t1028 * t891;
t1037 = -m(6) * t905 + t938 * mrSges(6,1) - t939 * mrSges(6,2) + t974 * t949 - t975 * t950 - t882;
t908 = 0.2e1 * t1062 + t1075;
t979 = mrSges(5,1) * t996 + mrSges(5,2) * t1016;
t1035 = -m(5) * t908 + t1015 * mrSges(5,3) + t1016 * t979 - t1037;
t913 = t988 + t1067;
t977 = mrSges(4,1) * t1016 - mrSges(4,3) * t996;
t878 = m(4) * t913 - t1015 * mrSges(4,2) - t1016 * t977 + (-mrSges(4,3) - mrSges(5,1)) * t969 + (-t962 - t963) * t995 + t1035;
t865 = t1021 * t878 + t1069 * t869;
t959 = -g(3) * t1054 + t1044;
t863 = m(3) * t959 + mrSges(3,1) * t1015 - mrSges(3,3) * t1005 + t1001 * t1016 - t1004 * t1051 + t865;
t1000 = mrSges(3,1) * t1016 - mrSges(3,3) * t1051;
t1048 = -t1021 * t869 + t1069 * t878;
t864 = m(3) * t960 - mrSges(3,2) * t1015 + mrSges(3,3) * t1006 - t1000 * t1016 + t1004 * t1050 + t1048;
t1046 = -t1025 * t887 + t1029 * t880;
t911 = t969 * pkin(3) + t1033;
t1041 = m(5) * t911 - t970 * mrSges(5,3) - t996 * t979 + t1046;
t872 = m(4) * t941 + t970 * mrSges(4,2) + t1063 * t995 + t1074 * t969 + t996 * t977 + t1041;
t871 = t872 + m(3) * t985 + t1005 * mrSges(3,2) - t1006 * mrSges(3,1) + (t1000 * t1026 - t1001 * t1030) * t1058;
t851 = -t1022 * t871 + t863 * t1052 + t864 * t1053;
t848 = m(2) * t1012 + qJDD(1) * mrSges(2,1) - mrSges(2,2) * t1032 + t851;
t857 = -t1026 * t863 + t1030 * t864;
t855 = m(2) * t1013 - mrSges(2,1) * t1032 - qJDD(1) * mrSges(2,2) + t857;
t1068 = t1027 * t855 + t1031 * t848;
t1066 = t1016 * t1077 + t1072 * t995 - t1073 * t996;
t1065 = t1016 * t1072 - t1071 * t996 + t1076 * t995;
t1064 = t1016 * t1073 + t1071 * t995 + t1078 * t996;
t850 = t1023 * t871 + t863 * t1054 + t864 * t1055;
t1045 = -t1027 * t848 + t1031 * t855;
t917 = Ifges(7,5) * t948 + Ifges(7,6) * t947 + Ifges(7,3) * t973;
t919 = Ifges(7,1) * t948 + Ifges(7,4) * t947 + Ifges(7,5) * t973;
t885 = -mrSges(7,1) * t896 + mrSges(7,3) * t895 + Ifges(7,4) * t916 + Ifges(7,2) * t915 + Ifges(7,6) * t937 - t917 * t948 + t919 * t973;
t918 = Ifges(7,4) * t948 + Ifges(7,2) * t947 + Ifges(7,6) * t973;
t886 = mrSges(7,2) * t896 - mrSges(7,3) * t894 + Ifges(7,1) * t916 + Ifges(7,4) * t915 + Ifges(7,5) * t937 + t917 * t947 - t918 * t973;
t930 = Ifges(6,5) * t975 + Ifges(6,6) * t974 + Ifges(6,3) * t993;
t931 = Ifges(6,4) * t975 + Ifges(6,2) * t974 + Ifges(6,6) * t993;
t866 = mrSges(6,2) * t905 - mrSges(6,3) * t899 + Ifges(6,1) * t939 + Ifges(6,4) * t938 + Ifges(6,5) * t968 - pkin(10) * t882 - t1024 * t885 + t1028 * t886 + t930 * t974 - t931 * t993;
t1034 = mrSges(7,1) * t894 - mrSges(7,2) * t895 + Ifges(7,5) * t916 + Ifges(7,6) * t915 + Ifges(7,3) * t937 + t918 * t948 - t919 * t947;
t932 = Ifges(6,1) * t975 + Ifges(6,4) * t974 + Ifges(6,5) * t993;
t867 = -mrSges(6,1) * t905 + mrSges(6,3) * t900 + Ifges(6,4) * t939 + Ifges(6,2) * t938 + Ifges(6,6) * t968 - pkin(5) * t882 - t930 * t975 + t932 * t993 - t1034;
t874 = -t969 * mrSges(5,2) - t995 * t978 + t1041;
t846 = -mrSges(4,1) * t941 - mrSges(5,1) * t908 + mrSges(5,2) * t911 + mrSges(4,3) * t913 - pkin(3) * t874 - pkin(4) * t1037 - pkin(9) * t1046 + t1072 * t1015 + t1064 * t1016 - t1025 * t866 - t1029 * t867 + t1066 * t996 - t1071 * t970 + t1076 * t969;
t1036 = mrSges(6,1) * t899 - mrSges(6,2) * t900 + Ifges(6,5) * t939 + Ifges(6,6) * t938 + Ifges(6,3) * t968 + pkin(5) * t1040 + pkin(10) * t1047 + t1024 * t886 + t1028 * t885 + t975 * t931 - t974 * t932;
t852 = mrSges(5,1) * t909 + mrSges(4,2) * t941 - mrSges(4,3) * t912 - mrSges(5,3) * t911 + pkin(4) * t875 - qJ(4) * t874 + t1073 * t1015 - t1065 * t1016 + t1066 * t995 + t1071 * t969 + t1078 * t970 + t1036;
t982 = Ifges(3,3) * t1016 + (Ifges(3,5) * t1026 + Ifges(3,6) * t1030) * t1058;
t984 = Ifges(3,5) * t1016 + (Ifges(3,1) * t1026 + Ifges(3,4) * t1030) * t1058;
t841 = -mrSges(3,1) * t985 + mrSges(3,3) * t960 + Ifges(3,4) * t1005 + Ifges(3,2) * t1006 + Ifges(3,6) * t1015 - pkin(2) * t872 + qJ(3) * t1048 + t1016 * t984 + t1021 * t852 - t1051 * t982 + t1069 * t846;
t983 = Ifges(3,6) * t1016 + (Ifges(3,4) * t1026 + Ifges(3,2) * t1030) * t1058;
t843 = mrSges(3,2) * t985 - mrSges(3,3) * t959 + Ifges(3,1) * t1005 + Ifges(3,4) * t1006 + Ifges(3,5) * t1015 - qJ(3) * t865 - t1016 * t983 - t1021 * t846 + t1050 * t982 + t1069 * t852;
t873 = t1015 * mrSges(5,2) + t1016 * t978 - t1039;
t845 = mrSges(4,1) * t912 - mrSges(4,2) * t913 - t1025 * t867 + pkin(2) * t865 - mrSges(5,3) * t908 + mrSges(5,2) * t909 + t1029 * t866 + Ifges(3,5) * t1005 + mrSges(3,1) * t959 - mrSges(3,2) * t960 - pkin(9) * t875 + Ifges(3,6) * t1006 + qJ(4) * t1035 - pkin(3) * t873 + t1065 * t996 + (-qJ(4) * t963 + t1064) * t995 + t1073 * t970 + (-mrSges(5,1) * qJ(4) - t1072) * t969 + (t1026 * t983 - t1030 * t984) * t1058 + (Ifges(3,3) - t1077) * t1015;
t1038 = mrSges(2,1) * t1012 - mrSges(2,2) * t1013 + Ifges(2,3) * qJDD(1) + pkin(1) * t851 + t1023 * t845 + t841 * t1054 + t843 * t1055 + t857 * t1070;
t839 = -mrSges(2,2) * g(3) - mrSges(2,3) * t1012 + Ifges(2,5) * qJDD(1) - t1032 * Ifges(2,6) - t1026 * t841 + t1030 * t843 + (-t1022 * t850 - t1023 * t851) * pkin(8);
t838 = mrSges(2,1) * g(3) + mrSges(2,3) * t1013 + t1032 * Ifges(2,5) + Ifges(2,6) * qJDD(1) - pkin(1) * t850 - t1022 * t845 + (pkin(8) * t857 + t1026 * t843 + t1030 * t841) * t1023;
t1 = [-m(1) * g(1) + t1045; -m(1) * g(2) + t1068; (-m(1) - m(2)) * g(3) + t850; -mrSges(1,2) * g(3) + mrSges(1,3) * g(2) - pkin(7) * t1068 - t1027 * t838 + t1031 * t839; mrSges(1,1) * g(3) - mrSges(1,3) * g(1) + pkin(7) * t1045 + t1027 * t839 + t1031 * t838; -mrSges(1,1) * g(2) + mrSges(1,2) * g(1) + t1038; t1038; t845; t872; t873; t1036; t1034;];
tauJB  = t1;
