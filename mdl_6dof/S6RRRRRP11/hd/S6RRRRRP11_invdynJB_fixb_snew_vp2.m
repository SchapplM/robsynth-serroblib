% Calculate vector of inverse dynamics joint torques and base forces with Newton-Euler
% S6RRRRRP11
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
% pkin [12x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,alpha2,alpha3,d1,d2,d3,d4,d5]';
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
% Datum: 2019-05-08 07:15
% Revision: 89c353f7eff3bd693eda4e29f35b2761dbc3ada0 (2019-05-03)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauJB = S6RRRRRP11_invdynJB_fixb_snew_vp2(qJ, qJD, qJDD, g, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(6,1),zeros(3,1),zeros(12,1),zeros(7,1),zeros(7,3),zeros(7,6)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRRRRP11_invdynJB_fixb_snew_vp2: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RRRRRP11_invdynJB_fixb_snew_vp2: qJD has to be [6x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [6 1]), ...
  'S6RRRRRP11_invdynJB_fixb_snew_vp2: qJDD has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6RRRRRP11_invdynJB_fixb_snew_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [12 1]), ...
  'S6RRRRRP11_invdynJB_fixb_snew_vp2: pkin has to be [12x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6RRRRRP11_invdynJB_fixb_snew_vp2: m has to be [7x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [7,3]), ...
  'S6RRRRRP11_invdynJB_fixb_snew_vp2: mrSges has to be [7x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [7 6]), ...
  'S6RRRRRP11_invdynJB_fixb_snew_vp2: Ifges has to be [7x6] (double)'); 

%% Symbolic Calculation
% From invdyn_fixb_NewtonEuler_linkframe_tauJB_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-05-08 06:47:17
% EndTime: 2019-05-08 06:49:10
% DurationCPUTime: 83.18s
% Computational Cost: add. (1389807->397), mult. (3428212->514), div. (0->0), fcn. (2885023->14), ass. (0->174)
t1083 = Ifges(6,1) + Ifges(7,1);
t1078 = Ifges(6,4) + Ifges(7,4);
t1077 = Ifges(6,5) + Ifges(7,5);
t1082 = Ifges(6,2) + Ifges(7,2);
t1076 = Ifges(6,6) + Ifges(7,6);
t1081 = Ifges(6,3) + Ifges(7,3);
t1021 = sin(pkin(7));
t1023 = cos(pkin(7));
t1027 = sin(qJ(3));
t1032 = cos(qJ(3));
t1024 = cos(pkin(6));
t1018 = qJD(1) * t1024 + qJD(2);
t1022 = sin(pkin(6));
t1033 = cos(qJ(2));
t1058 = t1022 * t1033;
t1050 = qJD(1) * t1058;
t1002 = (t1018 * t1021 + t1023 * t1050) * pkin(10);
t1028 = sin(qJ(2));
t1064 = qJD(1) * t1022;
t1073 = pkin(10) * t1021;
t1006 = (-pkin(2) * t1033 - t1028 * t1073) * t1064;
t1062 = qJD(1) * t1033;
t1012 = (qJD(2) * t1062 + qJDD(1) * t1028) * t1022;
t1017 = qJDD(1) * t1024 + qJDD(2);
t1029 = sin(qJ(1));
t1034 = cos(qJ(1));
t1015 = t1029 * g(1) - g(2) * t1034;
t1035 = qJD(1) ^ 2;
t1074 = pkin(9) * t1022;
t1009 = qJDD(1) * pkin(1) + t1035 * t1074 + t1015;
t1016 = -g(1) * t1034 - g(2) * t1029;
t1010 = -pkin(1) * t1035 + qJDD(1) * t1074 + t1016;
t1054 = t1024 * t1033;
t1046 = t1009 * t1054 - t1028 * t1010;
t1063 = qJD(1) * t1028;
t1072 = pkin(10) * t1023;
t962 = -t1012 * t1072 + t1017 * pkin(2) + t1018 * t1002 + (-g(3) * t1033 - t1006 * t1063) * t1022 + t1046;
t1059 = t1022 * t1028;
t1051 = qJD(1) * t1059;
t1005 = pkin(2) * t1018 - t1051 * t1072;
t1013 = (-qJD(2) * t1063 + qJDD(1) * t1033) * t1022;
t1042 = t1013 * t1023 + t1017 * t1021;
t1055 = t1024 * t1028;
t1065 = t1009 * t1055 + t1033 * t1010;
t963 = -t1018 * t1005 + (-g(3) * t1028 + t1006 * t1062) * t1022 + t1042 * pkin(10) + t1065;
t1071 = t1024 * g(3);
t968 = -t1012 * t1073 - t1013 * pkin(2) - t1071 + (-t1009 + (-t1002 * t1033 + t1005 * t1028) * qJD(1)) * t1022;
t929 = -t1027 * t963 + (t1021 * t968 + t1023 * t962) * t1032;
t1056 = t1023 * t1033;
t1061 = t1021 * t1027;
t993 = t1018 * t1061 + (t1027 * t1056 + t1028 * t1032) * t1064;
t979 = -t993 * qJD(3) - t1027 * t1012 + t1032 * t1042;
t1060 = t1021 * t1032;
t992 = (-t1027 * t1028 + t1032 * t1056) * t1064 + t1018 * t1060;
t1025 = sin(qJ(5));
t1030 = cos(qJ(5));
t1003 = t1018 * t1023 - t1021 * t1050 + qJD(3);
t1026 = sin(qJ(4));
t1031 = cos(qJ(4));
t985 = t1003 * t1026 + t1031 * t993;
t991 = qJD(4) - t992;
t971 = -t1025 * t985 + t1030 * t991;
t972 = t1025 * t991 + t1030 * t985;
t984 = t1003 * t1031 - t1026 * t993;
t983 = qJD(5) - t984;
t1067 = t1077 * t983 + t1078 * t971 + t1083 * t972;
t1068 = -t1076 * t983 - t1078 * t972 - t1082 * t971;
t1001 = t1003 ^ 2;
t1057 = t1023 * t1027;
t930 = t1032 * t963 + t962 * t1057 + t968 * t1061;
t982 = -pkin(3) * t992 - pkin(11) * t993;
t994 = -t1013 * t1021 + t1017 * t1023 + qJDD(3);
t921 = -pkin(3) * t1001 + pkin(11) * t994 + t982 * t992 + t930;
t940 = -t1021 * t962 + t1023 * t968;
t980 = t992 * qJD(3) + t1032 * t1012 + t1027 * t1042;
t923 = (-t1003 * t992 - t980) * pkin(11) + (t1003 * t993 - t979) * pkin(3) + t940;
t917 = t1026 * t923 + t1031 * t921;
t965 = -pkin(4) * t984 - pkin(12) * t985;
t978 = qJDD(4) - t979;
t990 = t991 ^ 2;
t912 = -pkin(4) * t990 + pkin(12) * t978 + t965 * t984 + t917;
t920 = -t994 * pkin(3) - t1001 * pkin(11) + t993 * t982 - t929;
t947 = -qJD(4) * t985 - t1026 * t980 + t1031 * t994;
t948 = qJD(4) * t984 + t1026 * t994 + t1031 * t980;
t915 = (-t984 * t991 - t948) * pkin(12) + (t985 * t991 - t947) * pkin(4) + t920;
t907 = -t1025 * t912 + t1030 * t915;
t928 = qJD(5) * t971 + t1025 * t978 + t1030 * t948;
t946 = qJDD(5) - t947;
t903 = -0.2e1 * qJD(6) * t972 + (t971 * t983 - t928) * qJ(6) + (t971 * t972 + t946) * pkin(5) + t907;
t950 = -mrSges(7,2) * t983 + mrSges(7,3) * t971;
t1053 = m(7) * t903 + t946 * mrSges(7,1) + t983 * t950;
t942 = -mrSges(7,1) * t971 + mrSges(7,2) * t972;
t901 = -t928 * mrSges(7,3) - t972 * t942 + t1053;
t908 = t1025 * t915 + t1030 * t912;
t927 = -qJD(5) * t972 - t1025 * t948 + t1030 * t978;
t952 = pkin(5) * t983 - qJ(6) * t972;
t970 = t971 ^ 2;
t906 = -pkin(5) * t970 + qJ(6) * t927 + 0.2e1 * qJD(6) * t971 - t952 * t983 + t908;
t1080 = mrSges(6,1) * t907 + mrSges(7,1) * t903 - mrSges(6,2) * t908 - mrSges(7,2) * t906 + pkin(5) * t901 - t1067 * t971 - t1068 * t972 + t1076 * t927 + t1077 * t928 + t1081 * t946;
t1079 = -mrSges(6,2) - mrSges(7,2);
t1008 = -mrSges(3,2) * t1018 + mrSges(3,3) * t1050;
t1011 = (-mrSges(3,1) * t1033 + mrSges(3,2) * t1028) * t1064;
t943 = -mrSges(6,1) * t971 + mrSges(6,2) * t972;
t951 = -mrSges(6,2) * t983 + mrSges(6,3) * t971;
t895 = m(6) * t907 + t946 * mrSges(6,1) + t983 * t951 + (-t942 - t943) * t972 + (-mrSges(6,3) - mrSges(7,3)) * t928 + t1053;
t1052 = m(7) * t906 + t927 * mrSges(7,3) + t971 * t942;
t953 = mrSges(7,1) * t983 - mrSges(7,3) * t972;
t1066 = -mrSges(6,1) * t983 + mrSges(6,3) * t972 - t953;
t897 = m(6) * t908 + t927 * mrSges(6,3) + t1066 * t983 + t1079 * t946 + t971 * t943 + t1052;
t894 = -t1025 * t895 + t1030 * t897;
t964 = -mrSges(5,1) * t984 + mrSges(5,2) * t985;
t974 = mrSges(5,1) * t991 - mrSges(5,3) * t985;
t891 = m(5) * t917 - mrSges(5,2) * t978 + mrSges(5,3) * t947 + t964 * t984 - t974 * t991 + t894;
t916 = -t1026 * t921 + t1031 * t923;
t911 = -pkin(4) * t978 - pkin(12) * t990 + t985 * t965 - t916;
t909 = -pkin(5) * t927 - qJ(6) * t970 + t952 * t972 + qJDD(6) + t911;
t1045 = -m(7) * t909 + t927 * mrSges(7,1) + t971 * t950;
t900 = -m(6) * t911 + t927 * mrSges(6,1) + t1066 * t972 + t1079 * t928 + t971 * t951 + t1045;
t973 = -mrSges(5,2) * t991 + mrSges(5,3) * t984;
t899 = m(5) * t916 + t978 * mrSges(5,1) - t948 * mrSges(5,3) - t985 * t964 + t991 * t973 + t900;
t1048 = -t1026 * t899 + t1031 * t891;
t981 = -mrSges(4,1) * t992 + mrSges(4,2) * t993;
t987 = mrSges(4,1) * t1003 - mrSges(4,3) * t993;
t881 = m(4) * t930 - mrSges(4,2) * t994 + mrSges(4,3) * t979 - t1003 * t987 + t981 * t992 + t1048;
t884 = t1026 * t891 + t1031 * t899;
t986 = -mrSges(4,2) * t1003 + mrSges(4,3) * t992;
t883 = m(4) * t940 - mrSges(4,1) * t979 + mrSges(4,2) * t980 - t986 * t992 + t987 * t993 + t884;
t893 = t1025 * t897 + t1030 * t895;
t1037 = -m(5) * t920 + t947 * mrSges(5,1) - mrSges(5,2) * t948 + t984 * t973 - t974 * t985 - t893;
t888 = m(4) * t929 + mrSges(4,1) * t994 - mrSges(4,3) * t980 + t1003 * t986 - t981 * t993 + t1037;
t870 = t1023 * t1032 * t888 - t1021 * t883 + t881 * t1057;
t988 = -g(3) * t1058 + t1046;
t866 = m(3) * t988 + mrSges(3,1) * t1017 - mrSges(3,3) * t1012 + t1008 * t1018 - t1011 * t1051 + t870;
t1007 = mrSges(3,1) * t1018 - mrSges(3,3) * t1051;
t869 = t1023 * t883 + t888 * t1060 + t881 * t1061;
t998 = -t1022 * t1009 - t1071;
t868 = m(3) * t998 - t1013 * mrSges(3,1) + t1012 * mrSges(3,2) + (t1007 * t1028 - t1008 * t1033) * t1064 + t869;
t876 = -t1027 * t888 + t1032 * t881;
t989 = -g(3) * t1059 + t1065;
t875 = m(3) * t989 - mrSges(3,2) * t1017 + mrSges(3,3) * t1013 - t1007 * t1018 + t1011 * t1050 + t876;
t855 = -t1022 * t868 + t866 * t1054 + t875 * t1055;
t852 = m(2) * t1015 + qJDD(1) * mrSges(2,1) - mrSges(2,2) * t1035 + t855;
t862 = -t1028 * t866 + t1033 * t875;
t860 = m(2) * t1016 - mrSges(2,1) * t1035 - qJDD(1) * mrSges(2,2) + t862;
t1070 = t1029 * t860 + t1034 * t852;
t1069 = -t1076 * t971 - t1077 * t972 - t1081 * t983;
t854 = t1024 * t868 + t866 * t1058 + t875 * t1059;
t1047 = -t1029 * t852 + t1034 * t860;
t904 = t928 * mrSges(7,2) + t972 * t953 - t1045;
t885 = -mrSges(6,1) * t911 + mrSges(6,3) * t908 - mrSges(7,1) * t909 + mrSges(7,3) * t906 - pkin(5) * t904 + qJ(6) * t1052 + (-qJ(6) * t953 + t1067) * t983 + t1069 * t972 + (-mrSges(7,2) * qJ(6) + t1076) * t946 + t1078 * t928 + t1082 * t927;
t892 = mrSges(6,2) * t911 + mrSges(7,2) * t909 - mrSges(6,3) * t907 - mrSges(7,3) * t903 - qJ(6) * t901 + t1068 * t983 - t1069 * t971 + t1077 * t946 + t1078 * t927 + t1083 * t928;
t955 = Ifges(5,5) * t985 + Ifges(5,6) * t984 + Ifges(5,3) * t991;
t956 = Ifges(5,4) * t985 + Ifges(5,2) * t984 + Ifges(5,6) * t991;
t871 = mrSges(5,2) * t920 - mrSges(5,3) * t916 + Ifges(5,1) * t948 + Ifges(5,4) * t947 + Ifges(5,5) * t978 - pkin(12) * t893 - t1025 * t885 + t1030 * t892 + t955 * t984 - t956 * t991;
t957 = Ifges(5,1) * t985 + Ifges(5,4) * t984 + Ifges(5,5) * t991;
t877 = -mrSges(5,1) * t920 + mrSges(5,3) * t917 + Ifges(5,4) * t948 + Ifges(5,2) * t947 + Ifges(5,6) * t978 - pkin(4) * t893 - t985 * t955 + t991 * t957 - t1080;
t975 = Ifges(4,5) * t993 + Ifges(4,6) * t992 + Ifges(4,3) * t1003;
t976 = Ifges(4,4) * t993 + Ifges(4,2) * t992 + Ifges(4,6) * t1003;
t857 = mrSges(4,2) * t940 - mrSges(4,3) * t929 + Ifges(4,1) * t980 + Ifges(4,4) * t979 + Ifges(4,5) * t994 - pkin(11) * t884 - t1003 * t976 - t1026 * t877 + t1031 * t871 + t975 * t992;
t1036 = mrSges(5,1) * t916 - mrSges(5,2) * t917 + Ifges(5,5) * t948 + Ifges(5,6) * t947 + Ifges(5,3) * t978 + pkin(4) * t900 + pkin(12) * t894 + t1025 * t892 + t1030 * t885 + t985 * t956 - t984 * t957;
t977 = Ifges(4,1) * t993 + Ifges(4,4) * t992 + Ifges(4,5) * t1003;
t863 = -mrSges(4,1) * t940 + mrSges(4,3) * t930 + Ifges(4,4) * t980 + Ifges(4,2) * t979 + Ifges(4,6) * t994 - pkin(3) * t884 + t1003 * t977 - t993 * t975 - t1036;
t1040 = pkin(10) * t876 + t1027 * t857 + t1032 * t863;
t856 = mrSges(4,1) * t929 - mrSges(4,2) * t930 + Ifges(4,5) * t980 + Ifges(4,6) * t979 + Ifges(4,3) * t994 + pkin(3) * t1037 + pkin(11) * t1048 + t1026 * t871 + t1031 * t877 + t993 * t976 - t992 * t977;
t996 = Ifges(3,6) * t1018 + (Ifges(3,4) * t1028 + Ifges(3,2) * t1033) * t1064;
t997 = Ifges(3,5) * t1018 + (Ifges(3,1) * t1028 + Ifges(3,4) * t1033) * t1064;
t846 = mrSges(3,1) * t988 - mrSges(3,2) * t989 + Ifges(3,5) * t1012 + Ifges(3,6) * t1013 + Ifges(3,3) * t1017 + pkin(2) * t870 + t1023 * t856 + (t1028 * t996 - t1033 * t997) * t1064 + t1040 * t1021;
t995 = Ifges(3,3) * t1018 + (Ifges(3,5) * t1028 + Ifges(3,6) * t1033) * t1064;
t848 = -mrSges(3,1) * t998 + mrSges(3,3) * t989 + Ifges(3,4) * t1012 + Ifges(3,2) * t1013 + Ifges(3,6) * t1017 - pkin(2) * t869 + t1018 * t997 - t1021 * t856 + t1023 * t1040 - t1051 * t995;
t850 = t995 * t1050 + mrSges(3,2) * t998 - mrSges(3,3) * t988 + Ifges(3,1) * t1012 + Ifges(3,4) * t1013 + Ifges(3,5) * t1017 - t1018 * t996 - t1027 * t863 + t1032 * t857 + (-t1021 * t869 - t1023 * t870) * pkin(10);
t1039 = mrSges(2,1) * t1015 - mrSges(2,2) * t1016 + Ifges(2,3) * qJDD(1) + pkin(1) * t855 + t1024 * t846 + t848 * t1058 + t850 * t1059 + t862 * t1074;
t844 = -mrSges(2,2) * g(3) - mrSges(2,3) * t1015 + Ifges(2,5) * qJDD(1) - t1035 * Ifges(2,6) - t1028 * t848 + t1033 * t850 + (-t1022 * t854 - t1024 * t855) * pkin(9);
t843 = mrSges(2,1) * g(3) + mrSges(2,3) * t1016 + t1035 * Ifges(2,5) + Ifges(2,6) * qJDD(1) - pkin(1) * t854 - t1022 * t846 + (pkin(9) * t862 + t1028 * t850 + t1033 * t848) * t1024;
t1 = [-m(1) * g(1) + t1047; -m(1) * g(2) + t1070; (-m(1) - m(2)) * g(3) + t854; -mrSges(1,2) * g(3) + mrSges(1,3) * g(2) - pkin(8) * t1070 - t1029 * t843 + t1034 * t844; mrSges(1,1) * g(3) - mrSges(1,3) * g(1) + pkin(8) * t1047 + t1029 * t844 + t1034 * t843; -mrSges(1,1) * g(2) + mrSges(1,2) * g(1) + t1039; t1039; t846; t856; t1036; t1080; t904;];
tauJB  = t1;
