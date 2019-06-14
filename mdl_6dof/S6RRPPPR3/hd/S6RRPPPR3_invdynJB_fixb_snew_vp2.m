% Calculate vector of inverse dynamics joint torques and base forces with Newton-Euler
% S6RRPPPR3
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
% pkin [9x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d2,d6,theta5]';
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
% Datum: 2019-05-06 08:36
% Revision: 89c353f7eff3bd693eda4e29f35b2761dbc3ada0 (2019-05-03)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauJB = S6RRPPPR3_invdynJB_fixb_snew_vp2(qJ, qJD, qJDD, g, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(6,1),zeros(3,1),zeros(9,1),zeros(7,1),zeros(7,3),zeros(7,6)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRPPPR3_invdynJB_fixb_snew_vp2: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RRPPPR3_invdynJB_fixb_snew_vp2: qJD has to be [6x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [6 1]), ...
  'S6RRPPPR3_invdynJB_fixb_snew_vp2: qJDD has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6RRPPPR3_invdynJB_fixb_snew_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S6RRPPPR3_invdynJB_fixb_snew_vp2: pkin has to be [9x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6RRPPPR3_invdynJB_fixb_snew_vp2: m has to be [7x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [7,3]), ...
  'S6RRPPPR3_invdynJB_fixb_snew_vp2: mrSges has to be [7x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [7 6]), ...
  'S6RRPPPR3_invdynJB_fixb_snew_vp2: Ifges has to be [7x6] (double)'); 

%% Symbolic Calculation
% From invdyn_fixb_NewtonEuler_linkframe_tauJB_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-05-06 08:34:00
% EndTime: 2019-05-06 08:34:07
% DurationCPUTime: 5.34s
% Computational Cost: add. (54428->346), mult. (118777->411), div. (0->0), fcn. (63287->8), ass. (0->138)
t1088 = -Ifges(3,1) - Ifges(4,1) - Ifges(5,2);
t1065 = Ifges(3,4) - Ifges(4,5) + Ifges(5,4);
t1064 = Ifges(3,5) + Ifges(4,4) + Ifges(5,6);
t1063 = Ifges(3,6) - Ifges(4,6) + Ifges(5,5);
t1087 = Ifges(3,3) + Ifges(4,2) + Ifges(5,3);
t1086 = Ifges(4,3) + Ifges(5,1) + Ifges(3,2);
t1029 = sin(pkin(9));
t1030 = cos(pkin(9));
t1033 = sin(qJ(2));
t1036 = cos(qJ(2));
t1070 = qJD(1) * t1033;
t1004 = -qJD(2) * mrSges(4,1) + mrSges(4,2) * t1070;
t1002 = qJD(2) * mrSges(5,2) - mrSges(5,3) * t1070;
t1068 = qJD(3) * qJD(2);
t1017 = 0.2e1 * t1068;
t1038 = qJD(2) ^ 2;
t1069 = qJD(1) * t1036;
t1075 = -pkin(2) - qJ(5);
t1001 = -qJD(2) * pkin(3) - qJ(4) * t1070;
t1039 = qJD(1) ^ 2;
t1067 = t1036 ^ 2 * t1039;
t1034 = sin(qJ(1));
t1037 = cos(qJ(1));
t1009 = -g(1) * t1037 - g(2) * t1034;
t969 = -pkin(1) * t1039 + qJDD(1) * pkin(7) + t1009;
t943 = -t1033 * g(3) + t1036 * t969;
t989 = (-pkin(2) * t1036 - qJ(3) * t1033) * qJD(1);
t1083 = qJDD(2) * qJ(3) + t989 * t1069 + t943;
t1056 = qJD(2) * t1070;
t995 = qJDD(1) * t1036 - t1056;
t1084 = pkin(3) * t1067 + t995 * qJ(4) - qJD(2) * t1001 + 0.2e1 * qJD(4) * t1069 - t1083;
t992 = (pkin(4) * t1033 + qJ(5) * t1036) * qJD(1);
t910 = qJDD(2) * pkin(4) + t1075 * t1038 - t992 * t1069 + qJDD(5) + t1017 - t1084;
t948 = -qJDD(2) * t1030 + t1029 * t995;
t981 = qJD(2) * t1029 + t1030 * t1069;
t950 = pkin(5) * t1070 + pkin(8) * t981;
t980 = -qJD(2) * t1030 + t1029 * t1069;
t977 = t980 ^ 2;
t903 = -t948 * pkin(5) - t977 * pkin(8) - t981 * t950 + t910;
t1032 = sin(qJ(6));
t1035 = cos(qJ(6));
t940 = t1032 * t980 - t1035 * t981;
t949 = -qJDD(2) * t1029 - t1030 * t995;
t920 = -t940 * qJD(6) - t1032 * t949 + t1035 * t948;
t939 = t1032 * t981 + t1035 * t980;
t921 = t939 * qJD(6) + t1032 * t948 + t1035 * t949;
t1012 = qJD(6) + t1070;
t932 = -mrSges(7,2) * t1012 + t939 * mrSges(7,3);
t933 = mrSges(7,1) * t1012 - t940 * mrSges(7,3);
t1051 = m(7) * t903 - t920 * mrSges(7,1) + t921 * mrSges(7,2) - t939 * t932 + t940 * t933;
t946 = -mrSges(6,2) * t1070 + mrSges(6,3) * t980;
t947 = mrSges(6,1) * t1070 + mrSges(6,3) * t981;
t894 = m(6) * t910 - t948 * mrSges(6,1) + t949 * mrSges(6,2) - t980 * t946 - t981 * t947 + t1051;
t1073 = pkin(2) * t1038;
t918 = -0.2e1 * t1068 + t1073 + t1084;
t1043 = -m(5) * t918 + qJDD(2) * mrSges(5,1) - t995 * mrSges(5,3) + qJD(2) * t1002 + t894;
t930 = t1017 - t1073 + t1083;
t990 = (-mrSges(4,1) * t1036 - mrSges(4,3) * t1033) * qJD(1);
t1042 = m(4) * t930 + qJDD(2) * mrSges(4,3) + qJD(2) * t1004 + t990 * t1069 + t1043;
t1059 = -t1064 * qJD(2) + (t1033 * t1088 - t1065 * t1036) * qJD(1);
t1060 = -t1063 * qJD(2) + (-t1033 * t1065 - t1036 * t1086) * qJD(1);
t1078 = 2 * qJD(5);
t1058 = qJD(2) * t1069;
t1008 = t1034 * g(1) - t1037 * g(2);
t968 = -qJDD(1) * pkin(1) - t1039 * pkin(7) - t1008;
t994 = qJDD(1) * t1033 + t1058;
t1049 = -t995 * pkin(2) + t968 + (-t1058 - t994) * qJ(3);
t1044 = -qJ(4) * t1067 + qJDD(4) - t1049 + (0.2e1 * qJD(3) + t1001) * t1070;
t1074 = pkin(3) + qJ(5);
t906 = t1074 * t995 + (pkin(4) * t1036 + t1075 * t1033) * qJD(2) * qJD(1) + t1044 + t994 * pkin(4);
t942 = -t1036 * g(3) - t1033 * t969;
t1052 = t989 * t1070 + qJDD(3) - t942;
t1066 = t1036 * t1039;
t1081 = -0.2e1 * qJD(4) * t1070 + (t1058 - t994) * qJ(4);
t911 = (-pkin(4) - qJ(3)) * t1038 + (-pkin(3) * t1066 - qJD(1) * t992) * t1033 + (-pkin(2) - t1074) * qJDD(2) + t1052 + t1081;
t900 = -t1029 * t911 + t1030 * t906 + t981 * t1078;
t898 = (t980 * t1070 - t949) * pkin(8) + (-t980 * t981 + t994) * pkin(5) + t900;
t901 = t1029 * t906 + t1030 * t911 + t980 * t1078;
t899 = -pkin(5) * t977 + t948 * pkin(8) - t950 * t1070 + t901;
t896 = -t1032 * t899 + t1035 * t898;
t927 = -mrSges(7,1) * t939 + mrSges(7,2) * t940;
t987 = qJDD(6) + t994;
t892 = m(7) * t896 + mrSges(7,1) * t987 - t921 * mrSges(7,3) + t1012 * t932 - t940 * t927;
t897 = t1032 * t898 + t1035 * t899;
t893 = m(7) * t897 - mrSges(7,2) * t987 + t920 * mrSges(7,3) - t1012 * t933 + t939 * t927;
t1055 = -t1032 * t892 + t1035 * t893;
t922 = Ifges(7,5) * t940 + Ifges(7,6) * t939 + Ifges(7,3) * t1012;
t924 = Ifges(7,1) * t940 + Ifges(7,4) * t939 + Ifges(7,5) * t1012;
t883 = -mrSges(7,1) * t903 + mrSges(7,3) * t897 + Ifges(7,4) * t921 + Ifges(7,2) * t920 + Ifges(7,6) * t987 + t1012 * t924 - t940 * t922;
t923 = Ifges(7,4) * t940 + Ifges(7,2) * t939 + Ifges(7,6) * t1012;
t884 = mrSges(7,2) * t903 - mrSges(7,3) * t896 + Ifges(7,1) * t921 + Ifges(7,4) * t920 + Ifges(7,5) * t987 - t1012 * t923 + t939 * t922;
t934 = -Ifges(6,5) * t981 + Ifges(6,6) * t980 + Ifges(6,3) * t1070;
t936 = -Ifges(6,1) * t981 + Ifges(6,4) * t980 + Ifges(6,5) * t1070;
t864 = -mrSges(6,1) * t910 + mrSges(6,3) * t901 + Ifges(6,4) * t949 + Ifges(6,2) * t948 + Ifges(6,6) * t994 - pkin(5) * t1051 + pkin(8) * t1055 + t1032 * t884 + t1035 * t883 + t936 * t1070 + t981 * t934;
t882 = t1032 * t893 + t1035 * t892;
t935 = -Ifges(6,4) * t981 + Ifges(6,2) * t980 + Ifges(6,6) * t1070;
t865 = mrSges(6,2) * t910 - mrSges(6,3) * t900 + Ifges(6,1) * t949 + Ifges(6,4) * t948 + Ifges(6,5) * t994 - pkin(8) * t882 - t1032 * t883 + t1035 * t884 - t935 * t1070 + t934 * t980;
t1005 = -qJD(2) * mrSges(5,1) + mrSges(5,3) * t1069;
t941 = -mrSges(6,1) * t980 - mrSges(6,2) * t981;
t880 = m(6) * t900 + mrSges(6,1) * t994 - t949 * mrSges(6,3) + t946 * t1070 + t941 * t981 + t882;
t881 = m(6) * t901 - mrSges(6,2) * t994 + t948 * mrSges(6,3) - t947 * t1070 + t941 * t980 + t1055;
t876 = -t1029 * t880 + t1030 * t881;
t931 = -qJDD(2) * pkin(2) - t1038 * qJ(3) + t1052;
t919 = (-t1033 * t1066 - qJDD(2)) * pkin(3) + t931 + t1081;
t993 = (mrSges(5,1) * t1033 - mrSges(5,2) * t1036) * qJD(1);
t1050 = -m(5) * t919 + t993 * t1070 - t876;
t1045 = qJDD(2) * mrSges(5,2) + qJD(2) * t1005 - t1050;
t1007 = mrSges(4,2) * t1069 + qJD(2) * mrSges(4,3);
t1082 = m(4) * t931 - qJDD(2) * mrSges(4,1) - qJD(2) * t1007;
t872 = t990 * t1070 + (mrSges(4,2) - mrSges(5,3)) * t994 + t1045 + t1082;
t874 = -t994 * mrSges(5,3) + t1045;
t1085 = -(t1060 * t1033 - t1059 * t1036) * qJD(1) + t1087 * qJDD(2) + t1063 * t995 + t1064 * t994 + mrSges(3,1) * t942 - mrSges(4,1) * t931 - mrSges(5,1) * t918 - mrSges(3,2) * t943 + mrSges(5,2) * t919 + mrSges(4,3) * t930 - pkin(2) * t872 - pkin(3) * t874 + pkin(4) * t894 + qJ(3) * (t995 * mrSges(4,2) - t993 * t1069 + t1042) - qJ(5) * t876 - t1029 * t865 - t1030 * t864;
t1076 = mrSges(3,3) + mrSges(4,2);
t1006 = -qJD(2) * mrSges(3,2) + mrSges(3,3) * t1069;
t991 = (-mrSges(3,1) * t1036 + mrSges(3,2) * t1033) * qJD(1);
t870 = m(3) * t942 + (mrSges(3,1) - mrSges(5,2)) * qJDD(2) + (-t1005 + t1006) * qJD(2) + (-t990 - t991) * t1070 + (mrSges(5,3) - t1076) * t994 + t1050 - t1082;
t1003 = qJD(2) * mrSges(3,1) - mrSges(3,3) * t1070;
t887 = (t991 - t993) * t1069 + t1042 + t1076 * t995 - qJDD(2) * mrSges(3,2) - qJD(2) * t1003 + m(3) * t943;
t1054 = -t1033 * t870 + t1036 * t887;
t861 = m(2) * t1009 - mrSges(2,1) * t1039 - qJDD(1) * mrSges(2,2) + t1054;
t875 = t1029 * t881 + t1030 * t880;
t913 = -pkin(2) * t1056 + t995 * pkin(3) + t1044;
t873 = m(5) * t913 + t994 * mrSges(5,1) - t995 * mrSges(5,2) + t1002 * t1070 - t1005 * t1069 + t875;
t928 = (pkin(2) * qJD(2) - 0.2e1 * qJD(3)) * t1070 + t1049;
t871 = m(4) * t928 - mrSges(4,1) * t995 - t994 * mrSges(4,3) - t1004 * t1070 - t1007 * t1069 - t873;
t1041 = -m(3) * t968 + t995 * mrSges(3,1) - mrSges(3,2) * t994 - t1003 * t1070 + t1006 * t1069 - t871;
t867 = m(2) * t1008 + qJDD(1) * mrSges(2,1) - mrSges(2,2) * t1039 + t1041;
t1071 = t1034 * t861 + t1037 * t867;
t863 = t1033 * t887 + t1036 * t870;
t1061 = -t1087 * qJD(2) + (-t1033 * t1064 - t1036 * t1063) * qJD(1);
t1053 = -t1034 * t867 + t1037 * t861;
t856 = -pkin(2) * t871 - qJ(4) * t1043 + qJ(5) * t875 - mrSges(4,1) * t928 + mrSges(4,2) * t930 - mrSges(5,2) * t913 + mrSges(5,3) * t918 - mrSges(3,1) * t968 + mrSges(3,3) * t943 + t1029 * t864 - t1030 * t865 + pkin(3) * t873 + t1086 * t995 + t1065 * t994 + t1063 * qJDD(2) - t1059 * qJD(2) + (qJ(4) * t1036 * t993 + t1061 * t1033) * qJD(1);
t1047 = mrSges(7,1) * t896 - mrSges(7,2) * t897 + Ifges(7,5) * t921 + Ifges(7,6) * t920 + Ifges(7,3) * t987 + t940 * t923 - t939 * t924;
t858 = -t1061 * t1069 + t1065 * t995 + (Ifges(6,3) - t1088) * t994 - qJ(3) * t871 + t1060 * qJD(2) + t1064 * qJDD(2) + pkin(4) * t875 - mrSges(4,3) * t928 + mrSges(4,2) * t931 + mrSges(5,1) * t913 - mrSges(5,3) * t919 - qJ(4) * t874 + Ifges(6,6) * t948 + Ifges(6,5) * t949 + mrSges(3,2) * t968 - mrSges(3,3) * t942 - mrSges(6,2) * t901 + pkin(5) * t882 - t980 * t936 - t981 * t935 + t1047 + mrSges(6,1) * t900;
t1048 = mrSges(2,1) * t1008 - mrSges(2,2) * t1009 + Ifges(2,3) * qJDD(1) + pkin(1) * t1041 + pkin(7) * t1054 + t1033 * t858 + t1036 * t856;
t854 = mrSges(2,1) * g(3) + mrSges(2,3) * t1009 + t1039 * Ifges(2,5) + Ifges(2,6) * qJDD(1) - pkin(1) * t863 - t1085;
t853 = -mrSges(2,2) * g(3) - mrSges(2,3) * t1008 + Ifges(2,5) * qJDD(1) - Ifges(2,6) * t1039 - pkin(7) * t863 - t1033 * t856 + t1036 * t858;
t1 = [-m(1) * g(1) + t1053; -m(1) * g(2) + t1071; (-m(1) - m(2)) * g(3) + t863; -mrSges(1,2) * g(3) + mrSges(1,3) * g(2) - pkin(6) * t1071 - t1034 * t854 + t1037 * t853; mrSges(1,1) * g(3) - mrSges(1,3) * g(1) + pkin(6) * t1053 + t1034 * t853 + t1037 * t854; -mrSges(1,1) * g(2) + mrSges(1,2) * g(1) + t1048; t1048; t1085; t872; t873; t894; t1047;];
tauJB  = t1;
