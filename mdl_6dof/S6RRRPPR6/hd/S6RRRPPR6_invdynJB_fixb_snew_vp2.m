% Calculate vector of inverse dynamics joint torques and base forces with Newton-Euler
% S6RRRPPR6
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
%   pkin=[a2,a3,a4,a5,a6,alpha2,d1,d2,d3,d6,theta4]';
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
% Datum: 2019-05-07 05:37
% Revision: 89c353f7eff3bd693eda4e29f35b2761dbc3ada0 (2019-05-03)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauJB = S6RRRPPR6_invdynJB_fixb_snew_vp2(qJ, qJD, qJDD, g, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(6,1),zeros(3,1),zeros(11,1),zeros(7,1),zeros(7,3),zeros(7,6)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRRPPR6_invdynJB_fixb_snew_vp2: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RRRPPR6_invdynJB_fixb_snew_vp2: qJD has to be [6x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [6 1]), ...
  'S6RRRPPR6_invdynJB_fixb_snew_vp2: qJDD has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6RRRPPR6_invdynJB_fixb_snew_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [11 1]), ...
  'S6RRRPPR6_invdynJB_fixb_snew_vp2: pkin has to be [11x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6RRRPPR6_invdynJB_fixb_snew_vp2: m has to be [7x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [7,3]), ...
  'S6RRRPPR6_invdynJB_fixb_snew_vp2: mrSges has to be [7x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [7 6]), ...
  'S6RRRPPR6_invdynJB_fixb_snew_vp2: Ifges has to be [7x6] (double)'); 

%% Symbolic Calculation
% From invdyn_fixb_NewtonEuler_linkframe_tauJB_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-05-07 05:28:26
% EndTime: 2019-05-07 05:28:53
% DurationCPUTime: 25.21s
% Computational Cost: add. (408849->383), mult. (906373->477), div. (0->0), fcn. (717778->12), ass. (0->161)
t1065 = -2 * qJD(4);
t1064 = Ifges(5,1) + Ifges(6,2);
t1063 = Ifges(6,1) + Ifges(5,3);
t1059 = Ifges(5,4) + Ifges(6,6);
t1058 = Ifges(5,5) - Ifges(6,4);
t1062 = -Ifges(5,2) - Ifges(6,3);
t1057 = Ifges(5,6) - Ifges(6,5);
t1011 = sin(pkin(6));
t1015 = sin(qJ(2));
t1019 = cos(qJ(2));
t997 = (qJD(1) * qJD(2) * t1015 - qJDD(1) * t1019) * t1011;
t1010 = sin(pkin(11));
t1054 = cos(pkin(11));
t1041 = t1011 * t1019;
t1036 = qJD(1) * t1041;
t1000 = -qJD(3) + t1036;
t1014 = sin(qJ(3));
t1018 = cos(qJ(3));
t1012 = cos(pkin(6));
t1006 = qJD(1) * t1012 + qJD(2);
t1004 = t1006 ^ 2;
t1005 = qJDD(1) * t1012 + qJDD(2);
t1043 = qJD(1) * t1019;
t1040 = t1012 * t1015;
t1016 = sin(qJ(1));
t1020 = cos(qJ(1));
t1002 = t1016 * g(1) - g(2) * t1020;
t1021 = qJD(1) ^ 2;
t1056 = pkin(8) * t1011;
t992 = qJDD(1) * pkin(1) + t1021 * t1056 + t1002;
t1003 = -g(1) * t1020 - g(2) * t1016;
t993 = -pkin(1) * t1021 + qJDD(1) * t1056 + t1003;
t1047 = t1019 * t993 + t992 * t1040;
t1044 = qJD(1) * t1011;
t995 = (-pkin(2) * t1019 - pkin(9) * t1015) * t1044;
t943 = -t1004 * pkin(2) + t1005 * pkin(9) + (-g(3) * t1015 + t995 * t1043) * t1011 + t1047;
t1055 = t1012 * g(3);
t996 = (qJD(2) * t1043 + qJDD(1) * t1015) * t1011;
t944 = t997 * pkin(2) - t996 * pkin(9) - t1055 + (-t992 + (pkin(2) * t1015 - pkin(9) * t1019) * t1006 * qJD(1)) * t1011;
t908 = -t1014 * t943 + t1018 * t944;
t1042 = t1011 * t1015;
t1037 = qJD(1) * t1042;
t984 = t1006 * t1018 - t1014 * t1037;
t962 = qJD(3) * t984 + t1005 * t1014 + t1018 * t996;
t985 = t1006 * t1014 + t1018 * t1037;
t989 = qJDD(3) + t997;
t897 = (-t1000 * t984 - t962) * qJ(4) + (t984 * t985 + t989) * pkin(3) + t908;
t909 = t1014 * t944 + t1018 * t943;
t961 = -qJD(3) * t985 + t1005 * t1018 - t1014 * t996;
t974 = -pkin(3) * t1000 - qJ(4) * t985;
t983 = t984 ^ 2;
t900 = -pkin(3) * t983 + qJ(4) * t961 + t1000 * t974 + t909;
t971 = t1010 * t984 + t1054 * t985;
t892 = -t1010 * t900 + t1054 * t897 + t971 * t1065;
t1013 = sin(qJ(6));
t1017 = cos(qJ(6));
t1051 = t1010 * t897 + t1054 * t900;
t999 = t1000 ^ 2;
t1031 = t999 * pkin(4) - t989 * qJ(5) - t1051;
t1060 = -2 * qJD(5);
t925 = t1010 * t962 - t1054 * t961;
t970 = t1010 * t985 - t1054 * t984;
t936 = pkin(4) * t970 - qJ(5) * t971;
t952 = pkin(5) * t971 + pkin(10) * t1000;
t966 = t970 * t1065;
t969 = t970 ^ 2;
t888 = -t925 * pkin(5) - t969 * pkin(10) - t970 * t936 + t966 + (t1060 - t952) * t1000 - t1031;
t947 = -t1000 * t1017 + t1013 * t970;
t906 = -qJD(6) * t947 - t1013 * t989 + t1017 * t925;
t946 = t1000 * t1013 + t1017 * t970;
t907 = qJD(6) * t946 + t1013 * t925 + t1017 * t989;
t968 = qJD(6) + t971;
t917 = -mrSges(7,2) * t968 + mrSges(7,3) * t946;
t918 = mrSges(7,1) * t968 - mrSges(7,3) * t947;
t1029 = -m(7) * t888 + t906 * mrSges(7,1) - t907 * mrSges(7,2) + t946 * t917 - t947 * t918;
t890 = 0.2e1 * qJD(5) * t1000 + ((2 * qJD(4)) + t936) * t970 + t1031;
t949 = mrSges(6,1) * t971 - mrSges(6,2) * t1000;
t1025 = -m(6) * t890 + t989 * mrSges(6,3) - t1000 * t949 - t1029;
t1048 = t1058 * t1000 + t1059 * t970 - t1064 * t971;
t1049 = -t1057 * t1000 + t1059 * t971 + t1062 * t970;
t1045 = t1000 * t970;
t891 = -t989 * pkin(4) - t999 * qJ(5) + t971 * t936 + qJDD(5) - t892;
t926 = t1010 * t961 + t1054 * t962;
t886 = (t970 * t971 - t989) * pkin(10) + (t926 - t1045) * pkin(5) + t891;
t1039 = t1012 * t1019;
t963 = -g(3) * t1041 - t1015 * t993 + t992 * t1039;
t942 = -t1005 * pkin(2) - t1004 * pkin(9) + t995 * t1037 - t963;
t902 = -t961 * pkin(3) - t983 * qJ(4) + t985 * t974 + qJDD(4) + t942;
t1023 = (-t926 - t1045) * qJ(5) + t902 + (-t1000 * pkin(4) + t1060) * t971;
t889 = t1023 + (pkin(4) + pkin(10)) * t925 - t969 * pkin(5) - t971 * t952;
t884 = -t1013 * t889 + t1017 * t886;
t916 = -mrSges(7,1) * t946 + mrSges(7,2) * t947;
t924 = qJDD(6) + t926;
t881 = m(7) * t884 + mrSges(7,1) * t924 - mrSges(7,3) * t907 - t916 * t947 + t917 * t968;
t885 = t1013 * t886 + t1017 * t889;
t882 = m(7) * t885 - mrSges(7,2) * t924 + mrSges(7,3) * t906 + t916 * t946 - t918 * t968;
t872 = t1013 * t882 + t1017 * t881;
t938 = -mrSges(6,2) * t970 - mrSges(6,3) * t971;
t1028 = -m(6) * t891 - t926 * mrSges(6,1) - t971 * t938 - t872;
t937 = mrSges(5,1) * t970 + mrSges(5,2) * t971;
t948 = mrSges(6,1) * t970 + mrSges(6,3) * t1000;
t950 = mrSges(5,2) * t1000 - mrSges(5,3) * t970;
t865 = m(5) * t892 - t926 * mrSges(5,3) - t971 * t937 + (mrSges(5,1) - mrSges(6,2)) * t989 + (t948 - t950) * t1000 + t1028;
t893 = t966 + t1051;
t951 = -mrSges(5,1) * t1000 - mrSges(5,3) * t971;
t877 = m(5) * t893 - t989 * mrSges(5,2) + t1000 * t951 + (-t937 - t938) * t970 + (-mrSges(5,3) - mrSges(6,1)) * t925 + t1025;
t863 = t1010 * t877 + t1054 * t865;
t870 = t989 * mrSges(6,2) - t1000 * t948 - t1028;
t910 = Ifges(7,5) * t947 + Ifges(7,6) * t946 + Ifges(7,3) * t968;
t912 = Ifges(7,1) * t947 + Ifges(7,4) * t946 + Ifges(7,5) * t968;
t875 = -mrSges(7,1) * t888 + mrSges(7,3) * t885 + Ifges(7,4) * t907 + Ifges(7,2) * t906 + Ifges(7,6) * t924 - t910 * t947 + t912 * t968;
t911 = Ifges(7,4) * t947 + Ifges(7,2) * t946 + Ifges(7,6) * t968;
t876 = mrSges(7,2) * t888 - mrSges(7,3) * t884 + Ifges(7,1) * t907 + Ifges(7,4) * t906 + Ifges(7,5) * t924 + t910 * t946 - t911 * t968;
t956 = Ifges(4,4) * t985 + Ifges(4,2) * t984 - Ifges(4,6) * t1000;
t957 = Ifges(4,1) * t985 + Ifges(4,4) * t984 - Ifges(4,5) * t1000;
t1061 = (Ifges(4,3) + t1063) * t989 - t1048 * t970 + t1049 * t971 - t1057 * t925 + t1058 * t926 + mrSges(4,1) * t908 + mrSges(5,1) * t892 + mrSges(6,2) * t891 + Ifges(4,5) * t962 + Ifges(4,6) * t961 + pkin(3) * t863 + qJ(5) * (-t925 * mrSges(6,1) - t970 * t938 + t1025) + t1017 * t876 + t985 * t956 - mrSges(4,2) * t909 - mrSges(5,2) * t893 - mrSges(6,3) * t890 - pkin(4) * t870 - pkin(10) * t872 - t1013 * t875 - t984 * t957;
t972 = -mrSges(4,1) * t984 + mrSges(4,2) * t985;
t973 = mrSges(4,2) * t1000 + mrSges(4,3) * t984;
t861 = m(4) * t908 + mrSges(4,1) * t989 - mrSges(4,3) * t962 - t1000 * t973 - t972 * t985 + t863;
t1035 = -t1010 * t865 + t1054 * t877;
t975 = -mrSges(4,1) * t1000 - mrSges(4,3) * t985;
t862 = m(4) * t909 - mrSges(4,2) * t989 + mrSges(4,3) * t961 + t1000 * t975 + t972 * t984 + t1035;
t1034 = -t1014 * t861 + t1018 * t862;
t964 = -g(3) * t1042 + t1047;
t990 = mrSges(3,1) * t1006 - mrSges(3,3) * t1037;
t994 = (-mrSges(3,1) * t1019 + mrSges(3,2) * t1015) * t1044;
t853 = m(3) * t964 - mrSges(3,2) * t1005 - mrSges(3,3) * t997 - t1006 * t990 + t994 * t1036 + t1034;
t856 = t1014 * t862 + t1018 * t861;
t979 = -t1011 * t992 - t1055;
t991 = -mrSges(3,2) * t1006 + mrSges(3,3) * t1036;
t855 = m(3) * t979 + t997 * mrSges(3,1) + t996 * mrSges(3,2) + (t1015 * t990 - t1019 * t991) * t1044 + t856;
t1052 = -t1013 * t881 + t1017 * t882;
t895 = t925 * pkin(4) + t1023;
t871 = m(6) * t895 - t925 * mrSges(6,2) - t926 * mrSges(6,3) - t970 * t948 - t971 * t949 + t1052;
t869 = m(5) * t902 + t925 * mrSges(5,1) + t926 * mrSges(5,2) + t970 * t950 + t971 * t951 + t871;
t1024 = -m(4) * t942 + t961 * mrSges(4,1) - t962 * mrSges(4,2) + t984 * t973 - t985 * t975 - t869;
t868 = m(3) * t963 + t1005 * mrSges(3,1) - t996 * mrSges(3,3) + t1006 * t991 - t994 * t1037 + t1024;
t843 = -t1011 * t855 + t868 * t1039 + t853 * t1040;
t840 = m(2) * t1002 + qJDD(1) * mrSges(2,1) - mrSges(2,2) * t1021 + t843;
t848 = -t1015 * t868 + t1019 * t853;
t846 = m(2) * t1003 - mrSges(2,1) * t1021 - qJDD(1) * mrSges(2,2) + t848;
t1053 = t1016 * t846 + t1020 * t840;
t1050 = t1063 * t1000 + t1057 * t970 - t1058 * t971;
t842 = t1012 * t855 + t868 * t1041 + t853 * t1042;
t1033 = -t1016 * t840 + t1020 * t846;
t849 = -mrSges(5,1) * t902 - mrSges(6,1) * t890 + mrSges(6,2) * t895 + mrSges(5,3) * t893 - pkin(4) * t871 - pkin(5) * t1029 - pkin(10) * t1052 + t1048 * t1000 - t1013 * t876 - t1017 * t875 + t1050 * t971 + t1057 * t989 + t1059 * t926 + t1062 * t925;
t1026 = mrSges(7,1) * t884 - mrSges(7,2) * t885 + Ifges(7,5) * t907 + Ifges(7,6) * t906 + Ifges(7,3) * t924 + t947 * t911 - t946 * t912;
t857 = mrSges(6,1) * t891 + mrSges(5,2) * t902 - mrSges(5,3) * t892 - mrSges(6,3) * t895 + pkin(5) * t872 - qJ(5) * t871 + t1049 * t1000 + t1050 * t970 + t1058 * t989 - t1059 * t925 + t1064 * t926 + t1026;
t955 = Ifges(4,5) * t985 + Ifges(4,6) * t984 - Ifges(4,3) * t1000;
t835 = -mrSges(4,1) * t942 + mrSges(4,3) * t909 + Ifges(4,4) * t962 + Ifges(4,2) * t961 + Ifges(4,6) * t989 - pkin(3) * t869 + qJ(4) * t1035 - t1000 * t957 + t1010 * t857 + t1054 * t849 - t985 * t955;
t838 = mrSges(4,2) * t942 - mrSges(4,3) * t908 + Ifges(4,1) * t962 + Ifges(4,4) * t961 + Ifges(4,5) * t989 - qJ(4) * t863 + t1000 * t956 - t1010 * t849 + t1054 * t857 + t984 * t955;
t977 = Ifges(3,6) * t1006 + (Ifges(3,4) * t1015 + Ifges(3,2) * t1019) * t1044;
t978 = Ifges(3,5) * t1006 + (Ifges(3,1) * t1015 + Ifges(3,4) * t1019) * t1044;
t832 = Ifges(3,5) * t996 - Ifges(3,6) * t997 + Ifges(3,3) * t1005 + mrSges(3,1) * t963 - mrSges(3,2) * t964 + t1014 * t838 + t1018 * t835 + pkin(2) * t1024 + pkin(9) * t1034 + (t1015 * t977 - t1019 * t978) * t1044;
t976 = Ifges(3,3) * t1006 + (Ifges(3,5) * t1015 + Ifges(3,6) * t1019) * t1044;
t834 = mrSges(3,2) * t979 - mrSges(3,3) * t963 + Ifges(3,1) * t996 - Ifges(3,4) * t997 + Ifges(3,5) * t1005 - pkin(9) * t856 - t1006 * t977 - t1014 * t835 + t1018 * t838 + t976 * t1036;
t837 = -mrSges(3,1) * t979 + mrSges(3,3) * t964 + Ifges(3,4) * t996 - Ifges(3,2) * t997 + Ifges(3,6) * t1005 - pkin(2) * t856 + t1006 * t978 - t976 * t1037 - t1061;
t1027 = mrSges(2,1) * t1002 - mrSges(2,2) * t1003 + Ifges(2,3) * qJDD(1) + pkin(1) * t843 + t1012 * t832 + t837 * t1041 + t834 * t1042 + t848 * t1056;
t830 = -mrSges(2,2) * g(3) - mrSges(2,3) * t1002 + Ifges(2,5) * qJDD(1) - t1021 * Ifges(2,6) - t1015 * t837 + t1019 * t834 + (-t1011 * t842 - t1012 * t843) * pkin(8);
t829 = mrSges(2,1) * g(3) + mrSges(2,3) * t1003 + t1021 * Ifges(2,5) + Ifges(2,6) * qJDD(1) - pkin(1) * t842 - t1011 * t832 + (pkin(8) * t848 + t1015 * t834 + t1019 * t837) * t1012;
t1 = [-m(1) * g(1) + t1033; -m(1) * g(2) + t1053; (-m(1) - m(2)) * g(3) + t842; -mrSges(1,2) * g(3) + mrSges(1,3) * g(2) - pkin(7) * t1053 - t1016 * t829 + t1020 * t830; mrSges(1,1) * g(3) - mrSges(1,3) * g(1) + pkin(7) * t1033 + t1016 * t830 + t1020 * t829; -mrSges(1,1) * g(2) + mrSges(1,2) * g(1) + t1027; t1027; t832; t1061; t869; t870; t1026;];
tauJB  = t1;
