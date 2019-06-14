% Calculate vector of inverse dynamics joint torques and base forces with Newton-Euler
% S6RRRRPR10
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
%   pkin=[a2,a3,a4,a5,a6,alpha2,d1,d2,d3,d4,d6]';
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
% Datum: 2019-05-07 23:04
% Revision: 89c353f7eff3bd693eda4e29f35b2761dbc3ada0 (2019-05-03)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauJB = S6RRRRPR10_invdynJB_fixb_snew_vp2(qJ, qJD, qJDD, g, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(6,1),zeros(3,1),zeros(11,1),zeros(7,1),zeros(7,3),zeros(7,6)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRRRPR10_invdynJB_fixb_snew_vp2: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RRRRPR10_invdynJB_fixb_snew_vp2: qJD has to be [6x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [6 1]), ...
  'S6RRRRPR10_invdynJB_fixb_snew_vp2: qJDD has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6RRRRPR10_invdynJB_fixb_snew_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [11 1]), ...
  'S6RRRRPR10_invdynJB_fixb_snew_vp2: pkin has to be [11x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6RRRRPR10_invdynJB_fixb_snew_vp2: m has to be [7x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [7,3]), ...
  'S6RRRRPR10_invdynJB_fixb_snew_vp2: mrSges has to be [7x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [7 6]), ...
  'S6RRRRPR10_invdynJB_fixb_snew_vp2: Ifges has to be [7x6] (double)'); 

%% Symbolic Calculation
% From invdyn_fixb_NewtonEuler_linkframe_tauJB_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-05-07 22:51:55
% EndTime: 2019-05-07 22:52:31
% DurationCPUTime: 27.50s
% Computational Cost: add. (450222->381), mult. (965379->474), div. (0->0), fcn. (772833->12), ass. (0->162)
t1065 = Ifges(5,4) + Ifges(6,6);
t1075 = -Ifges(5,2) - Ifges(6,3);
t1071 = Ifges(5,6) - Ifges(6,5);
t1023 = sin(qJ(6));
t1028 = cos(qJ(6));
t1021 = sin(pkin(6));
t1030 = cos(qJ(2));
t1053 = t1021 * t1030;
t1048 = qJD(1) * t1053;
t1011 = qJD(3) - t1048;
t1008 = -qJD(4) - t1011;
t1007 = t1008 ^ 2;
t1024 = sin(qJ(4));
t1066 = cos(qJ(4));
t1025 = sin(qJ(3));
t1029 = cos(qJ(3));
t1026 = sin(qJ(2));
t1056 = qJD(1) * t1021;
t1004 = (-pkin(2) * t1030 - pkin(9) * t1026) * t1056;
t1022 = cos(pkin(6));
t1017 = qJD(1) * t1022 + qJD(2);
t1015 = t1017 ^ 2;
t1016 = qJDD(1) * t1022 + qJDD(2);
t1055 = qJD(1) * t1030;
t1027 = sin(qJ(1));
t1031 = cos(qJ(1));
t1012 = t1027 * g(1) - g(2) * t1031;
t1032 = qJD(1) ^ 2;
t1064 = pkin(8) * t1021;
t1001 = qJDD(1) * pkin(1) + t1032 * t1064 + t1012;
t1013 = -g(1) * t1031 - g(2) * t1027;
t1050 = qJDD(1) * t1021;
t1002 = -pkin(1) * t1032 + pkin(8) * t1050 + t1013;
t1052 = t1022 * t1026;
t1058 = t1001 * t1052 + t1030 * t1002;
t951 = -t1015 * pkin(2) + t1016 * pkin(9) + (-g(3) * t1026 + t1004 * t1055) * t1021 + t1058;
t1005 = (qJD(2) * t1055 + qJDD(1) * t1026) * t1021;
t1054 = t1021 * t1026;
t1049 = qJD(1) * t1054;
t1006 = -qJD(2) * t1049 + t1030 * t1050;
t1063 = t1022 * g(3);
t952 = -t1006 * pkin(2) - t1005 * pkin(9) - t1063 + (-t1001 + (pkin(2) * t1026 - pkin(9) * t1030) * t1017 * qJD(1)) * t1021;
t910 = -t1025 * t951 + t1029 * t952;
t993 = t1017 * t1029 - t1025 * t1049;
t970 = qJD(3) * t993 + t1005 * t1029 + t1016 * t1025;
t994 = t1017 * t1025 + t1029 * t1049;
t998 = qJDD(3) - t1006;
t899 = (t1011 * t993 - t970) * pkin(10) + (t993 * t994 + t998) * pkin(3) + t910;
t911 = t1025 * t952 + t1029 * t951;
t969 = -qJD(3) * t994 - t1005 * t1025 + t1016 * t1029;
t981 = pkin(3) * t1011 - pkin(10) * t994;
t992 = t993 ^ 2;
t902 = -pkin(3) * t992 + pkin(10) * t969 - t1011 * t981 + t911;
t897 = t1024 * t899 + t1066 * t902;
t976 = t1024 * t994 - t1066 * t993;
t977 = t1024 * t993 + t1066 * t994;
t944 = pkin(4) * t976 - qJ(5) * t977;
t997 = qJDD(4) + t998;
t1042 = -t1007 * pkin(4) + t997 * qJ(5) - t976 * t944 + t897;
t1067 = -2 * qJD(5);
t926 = t977 * qJD(4) + t1024 * t970 - t1066 * t969;
t960 = pkin(5) * t977 + pkin(11) * t1008;
t975 = t976 ^ 2;
t888 = -t926 * pkin(5) - t975 * pkin(11) + (t1067 - t960) * t1008 + t1042;
t955 = -t1008 * t1028 + t1023 * t976;
t906 = -qJD(6) * t955 - t1023 * t997 + t1028 * t926;
t954 = t1008 * t1023 + t1028 * t976;
t907 = qJD(6) * t954 + t1023 * t926 + t1028 * t997;
t974 = qJD(6) + t977;
t933 = -mrSges(7,2) * t974 + mrSges(7,3) * t954;
t934 = mrSges(7,1) * t974 - mrSges(7,3) * t955;
t1043 = -m(7) * t888 + t906 * mrSges(7,1) - t907 * mrSges(7,2) + t954 * t933 - t955 * t934;
t890 = 0.2e1 * qJD(5) * t1008 - t1042;
t957 = mrSges(6,1) * t977 - mrSges(6,2) * t1008;
t1038 = -m(6) * t890 + t997 * mrSges(6,3) - t1008 * t957 - t1043;
t1072 = Ifges(5,5) - Ifges(6,4);
t1073 = Ifges(5,1) + Ifges(6,2);
t1059 = t1072 * t1008 + t1065 * t976 - t1073 * t977;
t1069 = -t1071 * t1008 + t1065 * t977 + t1075 * t976;
t1070 = Ifges(5,3) + Ifges(6,1);
t1057 = t1008 * t976;
t896 = -t1024 * t902 + t1066 * t899;
t892 = -t997 * pkin(4) - t1007 * qJ(5) + t977 * t944 + qJDD(5) - t896;
t927 = -t976 * qJD(4) + t1024 * t969 + t1066 * t970;
t886 = (t976 * t977 - t997) * pkin(11) + (t927 - t1057) * pkin(5) + t892;
t1051 = t1022 * t1030;
t971 = -g(3) * t1053 + t1001 * t1051 - t1026 * t1002;
t950 = -t1016 * pkin(2) - t1015 * pkin(9) + t1004 * t1049 - t971;
t909 = -t969 * pkin(3) - t992 * pkin(10) + t994 * t981 + t950;
t1035 = (-t927 - t1057) * qJ(5) + t909 + (-pkin(4) * t1008 + t1067) * t977;
t889 = t1035 + (pkin(4) + pkin(11)) * t926 - t977 * t960 - t975 * pkin(5);
t884 = -t1023 * t889 + t1028 * t886;
t925 = qJDD(6) + t927;
t932 = -mrSges(7,1) * t954 + mrSges(7,2) * t955;
t881 = m(7) * t884 + mrSges(7,1) * t925 - mrSges(7,3) * t907 - t932 * t955 + t933 * t974;
t885 = t1023 * t886 + t1028 * t889;
t882 = m(7) * t885 - mrSges(7,2) * t925 + mrSges(7,3) * t906 + t932 * t954 - t934 * t974;
t870 = t1023 * t882 + t1028 * t881;
t946 = -mrSges(6,2) * t976 - mrSges(6,3) * t977;
t1041 = -m(6) * t892 - t927 * mrSges(6,1) - t977 * t946 - t870;
t956 = mrSges(6,1) * t976 + mrSges(6,3) * t1008;
t869 = t997 * mrSges(6,2) - t1008 * t956 - t1041;
t912 = Ifges(7,5) * t955 + Ifges(7,6) * t954 + Ifges(7,3) * t974;
t914 = Ifges(7,1) * t955 + Ifges(7,4) * t954 + Ifges(7,5) * t974;
t872 = -mrSges(7,1) * t888 + mrSges(7,3) * t885 + Ifges(7,4) * t907 + Ifges(7,2) * t906 + Ifges(7,6) * t925 - t912 * t955 + t914 * t974;
t913 = Ifges(7,4) * t955 + Ifges(7,2) * t954 + Ifges(7,6) * t974;
t873 = mrSges(7,2) * t888 - mrSges(7,3) * t884 + Ifges(7,1) * t907 + Ifges(7,4) * t906 + Ifges(7,5) * t925 + t912 * t954 - t913 * t974;
t1074 = -mrSges(5,2) * t897 - mrSges(6,3) * t890 - pkin(4) * t869 - pkin(11) * t870 - t1023 * t872 + t1028 * t873 + qJ(5) * (-t976 * t946 + t1038) + mrSges(6,2) * t892 + mrSges(5,1) * t896 + t1070 * t997 + t1069 * t977 + t1072 * t927 + (-mrSges(6,1) * qJ(5) - t1071) * t926 - t1059 * t976;
t945 = mrSges(5,1) * t976 + mrSges(5,2) * t977;
t958 = mrSges(5,2) * t1008 - mrSges(5,3) * t976;
t866 = m(5) * t896 - t927 * mrSges(5,3) - t977 * t945 + (mrSges(5,1) - mrSges(6,2)) * t997 + (t956 - t958) * t1008 + t1041;
t959 = -mrSges(5,1) * t1008 - mrSges(5,3) * t977;
t876 = m(5) * t897 - t997 * mrSges(5,2) + t1008 * t959 + (-t945 - t946) * t976 + (-mrSges(5,3) - mrSges(6,1)) * t926 + t1038;
t861 = t1024 * t876 + t1066 * t866;
t964 = Ifges(4,4) * t994 + Ifges(4,2) * t993 + Ifges(4,6) * t1011;
t965 = Ifges(4,1) * t994 + Ifges(4,4) * t993 + Ifges(4,5) * t1011;
t1068 = mrSges(4,1) * t910 - mrSges(4,2) * t911 + Ifges(4,5) * t970 + Ifges(4,6) * t969 + Ifges(4,3) * t998 + pkin(3) * t861 + t994 * t964 - t993 * t965 + t1074;
t1003 = (-mrSges(3,1) * t1030 + mrSges(3,2) * t1026) * t1056;
t978 = -mrSges(4,1) * t993 + mrSges(4,2) * t994;
t979 = -mrSges(4,2) * t1011 + mrSges(4,3) * t993;
t859 = m(4) * t910 + mrSges(4,1) * t998 - mrSges(4,3) * t970 + t1011 * t979 - t978 * t994 + t861;
t1046 = -t1024 * t866 + t1066 * t876;
t980 = mrSges(4,1) * t1011 - mrSges(4,3) * t994;
t860 = m(4) * t911 - mrSges(4,2) * t998 + mrSges(4,3) * t969 - t1011 * t980 + t978 * t993 + t1046;
t1045 = -t1025 * t859 + t1029 * t860;
t972 = -g(3) * t1054 + t1058;
t999 = mrSges(3,1) * t1017 - mrSges(3,3) * t1049;
t851 = m(3) * t972 - mrSges(3,2) * t1016 + mrSges(3,3) * t1006 + t1003 * t1048 - t1017 * t999 + t1045;
t1000 = -mrSges(3,2) * t1017 + mrSges(3,3) * t1048;
t854 = t1025 * t860 + t1029 * t859;
t985 = -t1021 * t1001 - t1063;
t853 = m(3) * t985 - t1006 * mrSges(3,1) + t1005 * mrSges(3,2) + (-t1000 * t1030 + t1026 * t999) * t1056 + t854;
t1061 = -t1023 * t881 + t1028 * t882;
t894 = t926 * pkin(4) + t1035;
t867 = m(6) * t894 - t926 * mrSges(6,2) - t927 * mrSges(6,3) - t976 * t956 - t977 * t957 + t1061;
t1037 = m(5) * t909 + t926 * mrSges(5,1) + t927 * mrSges(5,2) + t976 * t958 + t977 * t959 + t867;
t1034 = -m(4) * t950 + t969 * mrSges(4,1) - t970 * mrSges(4,2) + t993 * t979 - t994 * t980 - t1037;
t864 = m(3) * t971 + t1016 * mrSges(3,1) - t1005 * mrSges(3,3) + t1017 * t1000 - t1003 * t1049 + t1034;
t841 = -t1021 * t853 + t864 * t1051 + t851 * t1052;
t838 = m(2) * t1012 + qJDD(1) * mrSges(2,1) - mrSges(2,2) * t1032 + t841;
t846 = -t1026 * t864 + t1030 * t851;
t844 = m(2) * t1013 - mrSges(2,1) * t1032 - qJDD(1) * mrSges(2,2) + t846;
t1062 = t1027 * t844 + t1031 * t838;
t1060 = t1008 * t1070 + t1071 * t976 - t1072 * t977;
t840 = t1022 * t853 + t864 * t1053 + t851 * t1054;
t1044 = -t1027 * t838 + t1031 * t844;
t847 = -mrSges(5,1) * t909 - mrSges(6,1) * t890 + mrSges(6,2) * t894 + mrSges(5,3) * t897 - pkin(4) * t867 - pkin(5) * t1043 - pkin(11) * t1061 + t1059 * t1008 - t1023 * t873 - t1028 * t872 + t1060 * t977 + t1065 * t927 + t1071 * t997 + t1075 * t926;
t1039 = mrSges(7,1) * t884 - mrSges(7,2) * t885 + Ifges(7,5) * t907 + Ifges(7,6) * t906 + Ifges(7,3) * t925 + t955 * t913 - t954 * t914;
t855 = mrSges(6,1) * t892 + mrSges(5,2) * t909 - mrSges(5,3) * t896 - mrSges(6,3) * t894 + pkin(5) * t870 - qJ(5) * t867 + t1069 * t1008 + t1060 * t976 - t1065 * t926 + t1072 * t997 + t1073 * t927 + t1039;
t963 = Ifges(4,5) * t994 + Ifges(4,6) * t993 + Ifges(4,3) * t1011;
t833 = -mrSges(4,1) * t950 + mrSges(4,3) * t911 + Ifges(4,4) * t970 + Ifges(4,2) * t969 + Ifges(4,6) * t998 - pkin(3) * t1037 + pkin(10) * t1046 + t1011 * t965 + t1024 * t855 + t1066 * t847 - t994 * t963;
t836 = mrSges(4,2) * t950 - mrSges(4,3) * t910 + Ifges(4,1) * t970 + Ifges(4,4) * t969 + Ifges(4,5) * t998 - pkin(10) * t861 - t1011 * t964 - t1024 * t847 + t1066 * t855 + t993 * t963;
t983 = Ifges(3,6) * t1017 + (Ifges(3,4) * t1026 + Ifges(3,2) * t1030) * t1056;
t984 = Ifges(3,5) * t1017 + (Ifges(3,1) * t1026 + Ifges(3,4) * t1030) * t1056;
t830 = Ifges(3,5) * t1005 + Ifges(3,6) * t1006 + Ifges(3,3) * t1016 + mrSges(3,1) * t971 - mrSges(3,2) * t972 + t1025 * t836 + t1029 * t833 + pkin(2) * t1034 + pkin(9) * t1045 + (t1026 * t983 - t1030 * t984) * t1056;
t982 = Ifges(3,3) * t1017 + (Ifges(3,5) * t1026 + Ifges(3,6) * t1030) * t1056;
t832 = mrSges(3,2) * t985 - mrSges(3,3) * t971 + Ifges(3,1) * t1005 + Ifges(3,4) * t1006 + Ifges(3,5) * t1016 - pkin(9) * t854 - t1017 * t983 - t1025 * t833 + t1029 * t836 + t1048 * t982;
t835 = -mrSges(3,1) * t985 + mrSges(3,3) * t972 + Ifges(3,4) * t1005 + Ifges(3,2) * t1006 + Ifges(3,6) * t1016 - pkin(2) * t854 + t1017 * t984 - t982 * t1049 - t1068;
t1040 = mrSges(2,1) * t1012 - mrSges(2,2) * t1013 + Ifges(2,3) * qJDD(1) + pkin(1) * t841 + t1022 * t830 + t835 * t1053 + t832 * t1054 + t846 * t1064;
t828 = -mrSges(2,2) * g(3) - mrSges(2,3) * t1012 + Ifges(2,5) * qJDD(1) - t1032 * Ifges(2,6) - t1026 * t835 + t1030 * t832 + (-t1021 * t840 - t1022 * t841) * pkin(8);
t827 = mrSges(2,1) * g(3) + mrSges(2,3) * t1013 + t1032 * Ifges(2,5) + Ifges(2,6) * qJDD(1) - pkin(1) * t840 - t1021 * t830 + (pkin(8) * t846 + t1026 * t832 + t1030 * t835) * t1022;
t1 = [-m(1) * g(1) + t1044; -m(1) * g(2) + t1062; (-m(1) - m(2)) * g(3) + t840; -mrSges(1,2) * g(3) + mrSges(1,3) * g(2) - pkin(7) * t1062 - t1027 * t827 + t1031 * t828; mrSges(1,1) * g(3) - mrSges(1,3) * g(1) + pkin(7) * t1044 + t1027 * t828 + t1031 * t827; -mrSges(1,1) * g(2) + mrSges(1,2) * g(1) + t1040; t1040; t830; t1068; t1074; t869; t1039;];
tauJB  = t1;
