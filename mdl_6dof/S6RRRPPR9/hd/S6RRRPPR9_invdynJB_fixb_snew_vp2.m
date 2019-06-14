% Calculate vector of inverse dynamics joint torques and base forces with Newton-Euler
% S6RRRPPR9
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
% Datum: 2019-05-07 06:43
% Revision: 89c353f7eff3bd693eda4e29f35b2761dbc3ada0 (2019-05-03)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauJB = S6RRRPPR9_invdynJB_fixb_snew_vp2(qJ, qJD, qJDD, g, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(6,1),zeros(3,1),zeros(11,1),zeros(7,1),zeros(7,3),zeros(7,6)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRRPPR9_invdynJB_fixb_snew_vp2: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RRRPPR9_invdynJB_fixb_snew_vp2: qJD has to be [6x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [6 1]), ...
  'S6RRRPPR9_invdynJB_fixb_snew_vp2: qJDD has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6RRRPPR9_invdynJB_fixb_snew_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [11 1]), ...
  'S6RRRPPR9_invdynJB_fixb_snew_vp2: pkin has to be [11x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6RRRPPR9_invdynJB_fixb_snew_vp2: m has to be [7x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [7,3]), ...
  'S6RRRPPR9_invdynJB_fixb_snew_vp2: mrSges has to be [7x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [7 6]), ...
  'S6RRRPPR9_invdynJB_fixb_snew_vp2: Ifges has to be [7x6] (double)'); 

%% Symbolic Calculation
% From invdyn_fixb_NewtonEuler_linkframe_tauJB_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-05-07 06:33:35
% EndTime: 2019-05-07 06:34:01
% DurationCPUTime: 23.53s
% Computational Cost: add. (382330->379), mult. (836497->476), div. (0->0), fcn. (655577->12), ass. (0->159)
t1060 = -2 * qJD(4);
t1059 = Ifges(5,1) + Ifges(6,1);
t1052 = Ifges(5,4) - Ifges(6,5);
t1051 = -Ifges(5,5) - Ifges(6,4);
t1058 = Ifges(5,2) + Ifges(6,3);
t1057 = -Ifges(6,2) - Ifges(5,3);
t1050 = Ifges(5,6) - Ifges(6,6);
t1006 = sin(pkin(6));
t1010 = sin(qJ(2));
t1013 = cos(qJ(2));
t990 = (qJD(1) * qJD(2) * t1010 - qJDD(1) * t1013) * t1006;
t1005 = sin(pkin(11));
t1047 = cos(pkin(11));
t1009 = sin(qJ(3));
t1054 = cos(qJ(3));
t1007 = cos(pkin(6));
t1000 = qJDD(1) * t1007 + qJDD(2);
t1037 = qJD(1) * t1013;
t1034 = t1007 * t1010;
t1015 = qJD(1) ^ 2;
t1049 = pkin(8) * t1006;
t1011 = sin(qJ(1));
t1014 = cos(qJ(1));
t997 = t1011 * g(1) - g(2) * t1014;
t985 = qJDD(1) * pkin(1) + t1015 * t1049 + t997;
t998 = -g(1) * t1014 - g(2) * t1011;
t986 = -pkin(1) * t1015 + qJDD(1) * t1049 + t998;
t1040 = t1013 * t986 + t985 * t1034;
t1038 = qJD(1) * t1006;
t988 = (-pkin(2) * t1013 - pkin(9) * t1010) * t1038;
t1001 = qJD(1) * t1007 + qJD(2);
t999 = t1001 ^ 2;
t927 = -t999 * pkin(2) + t1000 * pkin(9) + (-g(3) * t1010 + t1037 * t988) * t1006 + t1040;
t1048 = t1007 * g(3);
t989 = (qJD(2) * t1037 + qJDD(1) * t1010) * t1006;
t928 = t990 * pkin(2) - t989 * pkin(9) - t1048 + (-t985 + (pkin(2) * t1010 - pkin(9) * t1013) * t1001 * qJD(1)) * t1006;
t904 = t1009 * t928 + t1054 * t927;
t1036 = t1006 * t1010;
t1032 = qJD(1) * t1036;
t978 = -t1054 * t1001 + t1009 * t1032;
t979 = t1009 * t1001 + t1032 * t1054;
t959 = pkin(3) * t978 - qJ(4) * t979;
t982 = qJDD(3) + t990;
t1035 = t1006 * t1013;
t1031 = qJD(1) * t1035;
t996 = qJD(3) - t1031;
t995 = t996 ^ 2;
t894 = -pkin(3) * t995 + qJ(4) * t982 - t959 * t978 + t904;
t1033 = t1007 * t1013;
t957 = -g(3) * t1035 - t1010 * t986 + t1033 * t985;
t926 = -t1000 * pkin(2) - t999 * pkin(9) + t988 * t1032 - t957;
t955 = qJD(3) * t979 - t1054 * t1000 + t1009 * t989;
t956 = -t978 * qJD(3) + t1009 * t1000 + t1054 * t989;
t897 = (t978 * t996 - t956) * qJ(4) + (t979 * t996 + t955) * pkin(3) + t926;
t965 = t1005 * t996 + t1047 * t979;
t889 = -t1005 * t894 + t1047 * t897 + t965 * t1060;
t903 = -t1009 * t927 + t1054 * t928;
t1021 = t982 * pkin(3) + t995 * qJ(4) - t979 * t959 - qJDD(4) + t903;
t964 = t1005 * t979 - t1047 * t996;
t1046 = t964 * t978;
t936 = t1005 * t982 + t1047 * t956;
t1056 = (-t936 + t1046) * qJ(5) - t1021;
t1055 = 2 * qJD(5);
t1053 = -mrSges(5,3) - mrSges(6,2);
t1008 = sin(qJ(6));
t1012 = cos(qJ(6));
t937 = pkin(4) * t964 - qJ(5) * t965;
t977 = t978 ^ 2;
t887 = -t955 * pkin(4) - t977 * qJ(5) + t965 * t937 + qJDD(5) - t889;
t881 = (-t936 - t1046) * pkin(10) + (t964 * t965 - t955) * pkin(5) + t887;
t890 = t1005 * t897 + t1047 * t894 + t964 * t1060;
t886 = -pkin(4) * t977 + t955 * qJ(5) + t978 * t1055 - t964 * t937 + t890;
t935 = t1005 * t956 - t1047 * t982;
t944 = -pkin(5) * t978 - pkin(10) * t965;
t963 = t964 ^ 2;
t882 = -pkin(5) * t963 + pkin(10) * t935 + t944 * t978 + t886;
t879 = -t1008 * t882 + t1012 * t881;
t931 = -t1008 * t965 + t1012 * t964;
t902 = qJD(6) * t931 + t1008 * t935 + t1012 * t936;
t932 = t1008 * t964 + t1012 * t965;
t909 = -mrSges(7,1) * t931 + mrSges(7,2) * t932;
t975 = qJD(6) - t978;
t912 = -mrSges(7,2) * t975 + mrSges(7,3) * t931;
t953 = qJDD(6) - t955;
t875 = m(7) * t879 + mrSges(7,1) * t953 - mrSges(7,3) * t902 - t909 * t932 + t912 * t975;
t880 = t1008 * t881 + t1012 * t882;
t901 = -qJD(6) * t932 - t1008 * t936 + t1012 * t935;
t913 = mrSges(7,1) * t975 - mrSges(7,3) * t932;
t876 = m(7) * t880 - mrSges(7,2) * t953 + mrSges(7,3) * t901 + t909 * t931 - t913 * t975;
t1029 = -t1008 * t875 + t1012 * t876;
t943 = -mrSges(6,1) * t978 + mrSges(6,2) * t965;
t1024 = m(6) * t886 + t955 * mrSges(6,3) + t978 * t943 + t1029;
t938 = mrSges(6,1) * t964 - mrSges(6,3) * t965;
t1041 = -mrSges(5,1) * t964 - mrSges(5,2) * t965 - t938;
t942 = mrSges(5,1) * t978 - mrSges(5,3) * t965;
t865 = m(5) * t890 - t955 * mrSges(5,2) + t1041 * t964 + t1053 * t935 - t978 * t942 + t1024;
t868 = t1008 * t876 + t1012 * t875;
t941 = -mrSges(6,2) * t964 + mrSges(6,3) * t978;
t1019 = -m(6) * t887 + t955 * mrSges(6,1) + t978 * t941 - t868;
t1026 = -mrSges(5,2) * t978 - mrSges(5,3) * t964;
t866 = m(5) * t889 + t955 * mrSges(5,1) + t1026 * t978 + t1041 * t965 + t1053 * t936 + t1019;
t863 = -t1005 * t866 + t1047 * t865;
t960 = mrSges(4,1) * t978 + mrSges(4,2) * t979;
t967 = mrSges(4,1) * t996 - mrSges(4,3) * t979;
t861 = m(4) * t904 - mrSges(4,2) * t982 - mrSges(4,3) * t955 - t960 * t978 - t967 * t996 + t863;
t884 = -t963 * pkin(10) + (-pkin(4) - pkin(5)) * t935 + (-pkin(4) * t978 + t1055 + t944) * t965 - t1056;
t1022 = -m(7) * t884 + t901 * mrSges(7,1) - t902 * mrSges(7,2) + t931 * t912 - t932 * t913;
t888 = -0.2e1 * qJD(5) * t965 + (t965 * t978 + t935) * pkin(4) + t1056;
t877 = m(6) * t888 + t935 * mrSges(6,1) - t936 * mrSges(6,3) + t964 * t941 - t965 * t943 + t1022;
t873 = -m(5) * t1021 + t935 * mrSges(5,1) + t936 * mrSges(5,2) + t964 * t1026 + t965 * t942 + t877;
t966 = -mrSges(4,2) * t996 - mrSges(4,3) * t978;
t872 = m(4) * t903 + t982 * mrSges(4,1) - t956 * mrSges(4,3) - t979 * t960 + t996 * t966 - t873;
t1028 = -t1009 * t872 + t1054 * t861;
t958 = -g(3) * t1036 + t1040;
t983 = mrSges(3,1) * t1001 - mrSges(3,3) * t1032;
t987 = (-mrSges(3,1) * t1013 + mrSges(3,2) * t1010) * t1038;
t851 = m(3) * t958 - mrSges(3,2) * t1000 - mrSges(3,3) * t990 - t1001 * t983 + t1031 * t987 + t1028;
t855 = t1009 * t861 + t1054 * t872;
t971 = -t1006 * t985 - t1048;
t984 = -mrSges(3,2) * t1001 + mrSges(3,3) * t1031;
t853 = m(3) * t971 + t990 * mrSges(3,1) + t989 * mrSges(3,2) + (t1010 * t983 - t1013 * t984) * t1038 + t855;
t862 = t1005 * t865 + t1047 * t866;
t1017 = -m(4) * t926 - t955 * mrSges(4,1) - t956 * mrSges(4,2) - t978 * t966 - t979 * t967 - t862;
t858 = m(3) * t957 + t1000 * mrSges(3,1) - t989 * mrSges(3,3) + t1001 * t984 - t1032 * t987 + t1017;
t840 = -t1006 * t853 + t858 * t1033 + t851 * t1034;
t837 = m(2) * t997 + qJDD(1) * mrSges(2,1) - mrSges(2,2) * t1015 + t840;
t846 = -t1010 * t858 + t1013 * t851;
t844 = m(2) * t998 - mrSges(2,1) * t1015 - qJDD(1) * mrSges(2,2) + t846;
t1045 = t1011 * t844 + t1014 * t837;
t1044 = -t1050 * t978 - t1052 * t965 + t1058 * t964;
t1043 = t1050 * t964 + t1051 * t965 + t1057 * t978;
t1042 = -t1051 * t978 - t1052 * t964 + t1059 * t965;
t839 = t1007 * t853 + t858 * t1035 + t851 * t1036;
t1027 = -t1011 * t837 + t1014 * t844;
t905 = Ifges(7,5) * t932 + Ifges(7,6) * t931 + Ifges(7,3) * t975;
t907 = Ifges(7,1) * t932 + Ifges(7,4) * t931 + Ifges(7,5) * t975;
t869 = -mrSges(7,1) * t884 + mrSges(7,3) * t880 + Ifges(7,4) * t902 + Ifges(7,2) * t901 + Ifges(7,6) * t953 - t905 * t932 + t907 * t975;
t906 = Ifges(7,4) * t932 + Ifges(7,2) * t931 + Ifges(7,6) * t975;
t870 = mrSges(7,2) * t884 - mrSges(7,3) * t879 + Ifges(7,1) * t902 + Ifges(7,4) * t901 + Ifges(7,5) * t953 + t905 * t931 - t906 * t975;
t847 = mrSges(5,1) * t1021 - mrSges(6,1) * t888 + mrSges(6,2) * t886 + mrSges(5,3) * t890 - pkin(4) * t877 - pkin(5) * t1022 - pkin(10) * t1029 - t1008 * t870 - t1012 * t869 + t1042 * t978 + t1043 * t965 + t1050 * t955 + t1052 * t936 - t1058 * t935;
t854 = -mrSges(5,2) * t1021 + mrSges(6,2) * t887 - mrSges(5,3) * t889 - mrSges(6,3) * t888 - pkin(10) * t868 - qJ(5) * t877 - t1008 * t869 + t1012 * t870 + t1043 * t964 + t1044 * t978 - t1051 * t955 - t1052 * t935 + t1059 * t936;
t946 = Ifges(4,5) * t979 - Ifges(4,6) * t978 + Ifges(4,3) * t996;
t947 = Ifges(4,4) * t979 - Ifges(4,2) * t978 + Ifges(4,6) * t996;
t835 = mrSges(4,2) * t926 - mrSges(4,3) * t903 + Ifges(4,1) * t956 - Ifges(4,4) * t955 + Ifges(4,5) * t982 - qJ(4) * t862 - t1005 * t847 + t1047 * t854 - t978 * t946 - t996 * t947;
t1018 = mrSges(7,1) * t879 - mrSges(7,2) * t880 + Ifges(7,5) * t902 + Ifges(7,6) * t901 + Ifges(7,3) * t953 + t932 * t906 - t931 * t907;
t867 = t936 * mrSges(6,2) + t965 * t938 - t1019;
t948 = Ifges(4,1) * t979 - Ifges(4,4) * t978 + Ifges(4,5) * t996;
t841 = t1044 * t965 + (qJ(5) * t938 - t1042) * t964 + (-Ifges(4,2) + t1057) * t955 + t1051 * t936 + (mrSges(6,2) * qJ(5) + t1050) * t935 + t1018 - qJ(5) * t1024 + t996 * t948 + Ifges(4,6) * t982 - t979 * t946 + Ifges(4,4) * t956 - mrSges(4,1) * t926 + mrSges(4,3) * t904 - mrSges(5,1) * t889 + mrSges(5,2) * t890 - mrSges(6,3) * t886 + mrSges(6,1) * t887 + pkin(5) * t868 + pkin(4) * t867 - pkin(3) * t862;
t969 = Ifges(3,6) * t1001 + (Ifges(3,4) * t1010 + Ifges(3,2) * t1013) * t1038;
t970 = Ifges(3,5) * t1001 + (Ifges(3,1) * t1010 + Ifges(3,4) * t1013) * t1038;
t830 = Ifges(3,5) * t989 - Ifges(3,6) * t990 + Ifges(3,3) * t1000 + mrSges(3,1) * t957 - mrSges(3,2) * t958 + t1009 * t835 + t1054 * t841 + pkin(2) * t1017 + pkin(9) * t1028 + (t1010 * t969 - t1013 * t970) * t1038;
t968 = Ifges(3,3) * t1001 + (Ifges(3,5) * t1010 + Ifges(3,6) * t1013) * t1038;
t832 = mrSges(3,2) * t971 - mrSges(3,3) * t957 + Ifges(3,1) * t989 - Ifges(3,4) * t990 + Ifges(3,5) * t1000 - pkin(9) * t855 - t1001 * t969 - t1009 * t841 + t1031 * t968 + t1054 * t835;
t1016 = mrSges(4,1) * t903 - mrSges(4,2) * t904 + Ifges(4,5) * t956 - Ifges(4,6) * t955 + Ifges(4,3) * t982 - pkin(3) * t873 + qJ(4) * t863 + t1005 * t854 + t1047 * t847 + t979 * t947 + t978 * t948;
t834 = -mrSges(3,1) * t971 + mrSges(3,3) * t958 + Ifges(3,4) * t989 - Ifges(3,2) * t990 + Ifges(3,6) * t1000 - pkin(2) * t855 + t1001 * t970 - t1032 * t968 - t1016;
t1020 = mrSges(2,1) * t997 - mrSges(2,2) * t998 + Ifges(2,3) * qJDD(1) + pkin(1) * t840 + t1007 * t830 + t834 * t1035 + t832 * t1036 + t846 * t1049;
t828 = -mrSges(2,2) * g(3) - mrSges(2,3) * t997 + Ifges(2,5) * qJDD(1) - t1015 * Ifges(2,6) - t1010 * t834 + t1013 * t832 + (-t1006 * t839 - t1007 * t840) * pkin(8);
t827 = mrSges(2,1) * g(3) + mrSges(2,3) * t998 + t1015 * Ifges(2,5) + Ifges(2,6) * qJDD(1) - pkin(1) * t839 - t1006 * t830 + (pkin(8) * t846 + t1010 * t832 + t1013 * t834) * t1007;
t1 = [-m(1) * g(1) + t1027; -m(1) * g(2) + t1045; (-m(1) - m(2)) * g(3) + t839; -mrSges(1,2) * g(3) + mrSges(1,3) * g(2) - pkin(7) * t1045 - t1011 * t827 + t1014 * t828; mrSges(1,1) * g(3) - mrSges(1,3) * g(1) + pkin(7) * t1027 + t1011 * t828 + t1014 * t827; -mrSges(1,1) * g(2) + mrSges(1,2) * g(1) + t1020; t1020; t830; t1016; t873; t867; t1018;];
tauJB  = t1;
