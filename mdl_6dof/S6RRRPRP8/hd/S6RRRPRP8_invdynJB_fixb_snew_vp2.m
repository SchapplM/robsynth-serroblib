% Calculate vector of inverse dynamics joint torques and base forces with Newton-Euler
% S6RRRPRP8
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
%   pkin=[a2,a3,a4,a5,a6,d1,d2,d3,d5]';
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
% Datum: 2019-05-07 08:37
% Revision: 89c353f7eff3bd693eda4e29f35b2761dbc3ada0 (2019-05-03)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauJB = S6RRRPRP8_invdynJB_fixb_snew_vp2(qJ, qJD, qJDD, g, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(6,1),zeros(3,1),zeros(9,1),zeros(7,1),zeros(7,3),zeros(7,6)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRRPRP8_invdynJB_fixb_snew_vp2: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RRRPRP8_invdynJB_fixb_snew_vp2: qJD has to be [6x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [6 1]), ...
  'S6RRRPRP8_invdynJB_fixb_snew_vp2: qJDD has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6RRRPRP8_invdynJB_fixb_snew_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S6RRRPRP8_invdynJB_fixb_snew_vp2: pkin has to be [9x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6RRRPRP8_invdynJB_fixb_snew_vp2: m has to be [7x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [7,3]), ...
  'S6RRRPRP8_invdynJB_fixb_snew_vp2: mrSges has to be [7x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [7 6]), ...
  'S6RRRPRP8_invdynJB_fixb_snew_vp2: Ifges has to be [7x6] (double)'); 

%% Symbolic Calculation
% From invdyn_fixb_NewtonEuler_linkframe_tauJB_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-05-07 08:29:58
% EndTime: 2019-05-07 08:30:14
% DurationCPUTime: 8.07s
% Computational Cost: add. (88585->342), mult. (174832->398), div. (0->0), fcn. (115432->8), ass. (0->139)
t1020 = Ifges(6,4) + Ifges(7,4);
t1038 = Ifges(6,2) + Ifges(7,2);
t1032 = Ifges(6,6) + Ifges(7,6);
t1034 = Ifges(6,5) + Ifges(7,5);
t1035 = Ifges(6,1) + Ifges(7,1);
t981 = sin(qJ(2));
t1008 = qJD(1) * t981;
t1024 = cos(qJ(3));
t980 = sin(qJ(3));
t956 = -qJD(2) * t1024 + t1008 * t980;
t957 = t980 * qJD(2) + t1008 * t1024;
t979 = sin(qJ(5));
t983 = cos(qJ(5));
t921 = t956 * t983 - t957 * t979;
t922 = t956 * t979 + t957 * t983;
t984 = cos(qJ(2));
t1007 = qJD(1) * t984;
t970 = qJD(3) - t1007;
t968 = qJD(5) - t970;
t1013 = t1020 * t921 + t1034 * t968 + t1035 * t922;
t1029 = t1020 * t922 + t1032 * t968 + t1038 * t921;
t1030 = Ifges(6,3) + Ifges(7,3);
t1016 = t956 * t970;
t1006 = qJD(1) * qJD(2);
t1002 = t984 * t1006;
t1003 = t981 * t1006;
t982 = sin(qJ(1));
t985 = cos(qJ(1));
t966 = t982 * g(1) - t985 * g(2);
t987 = qJD(1) ^ 2;
t946 = -qJDD(1) * pkin(1) - t987 * pkin(7) - t966;
t960 = qJDD(1) * t981 + t1002;
t961 = t984 * qJDD(1) - t1003;
t892 = (-t960 - t1002) * pkin(8) + (-t961 + t1003) * pkin(2) + t946;
t967 = -g(1) * t985 - g(2) * t982;
t947 = -pkin(1) * t987 + qJDD(1) * pkin(7) + t967;
t933 = -g(3) * t981 + t984 * t947;
t959 = (-pkin(2) * t984 - pkin(8) * t981) * qJD(1);
t986 = qJD(2) ^ 2;
t897 = -pkin(2) * t986 + qJDD(2) * pkin(8) + t1007 * t959 + t933;
t876 = t1024 * t892 - t980 * t897;
t925 = pkin(3) * t956 - qJ(4) * t957;
t955 = qJDD(3) - t961;
t969 = t970 ^ 2;
t862 = -t955 * pkin(3) - t969 * qJ(4) + t957 * t925 + qJDD(4) - t876;
t919 = -t956 * qJD(3) + t980 * qJDD(2) + t1024 * t960;
t854 = (-t919 - t1016) * pkin(9) + (t956 * t957 - t955) * pkin(4) + t862;
t1025 = 2 * qJD(4);
t877 = t1024 * t897 + t980 * t892;
t860 = -pkin(3) * t969 + t955 * qJ(4) + t1025 * t970 - t956 * t925 + t877;
t918 = t957 * qJD(3) - qJDD(2) * t1024 + t980 * t960;
t934 = -pkin(4) * t970 - pkin(9) * t957;
t954 = t956 ^ 2;
t857 = -pkin(4) * t954 + pkin(9) * t918 + t934 * t970 + t860;
t848 = t983 * t854 - t979 * t857;
t875 = qJD(5) * t921 + t918 * t979 + t919 * t983;
t951 = qJDD(5) - t955;
t844 = -0.2e1 * qJD(6) * t922 + (t921 * t968 - t875) * qJ(6) + (t921 * t922 + t951) * pkin(5) + t848;
t898 = -mrSges(7,2) * t968 + mrSges(7,3) * t921;
t1005 = m(7) * t844 + t951 * mrSges(7,1) + t968 * t898;
t889 = -mrSges(7,1) * t921 + mrSges(7,2) * t922;
t840 = -t875 * mrSges(7,3) - t922 * t889 + t1005;
t849 = t979 * t854 + t983 * t857;
t874 = -qJD(5) * t922 + t918 * t983 - t919 * t979;
t900 = pkin(5) * t968 - qJ(6) * t922;
t920 = t921 ^ 2;
t846 = -pkin(5) * t920 + qJ(6) * t874 + 0.2e1 * qJD(6) * t921 - t900 * t968 + t849;
t1037 = mrSges(6,1) * t848 + mrSges(7,1) * t844 - mrSges(6,2) * t849 - mrSges(7,2) * t846 + pkin(5) * t840 - t1013 * t921 + t1029 * t922 + t1030 * t951 + t1032 * t874 + t1034 * t875;
t1036 = Ifges(4,1) + Ifges(5,1);
t1021 = Ifges(4,4) - Ifges(5,5);
t1019 = Ifges(4,5) + Ifges(5,4);
t1033 = Ifges(4,2) + Ifges(5,3);
t1018 = Ifges(4,6) - Ifges(5,6);
t1031 = Ifges(4,3) + Ifges(5,2);
t1010 = t1019 * t970 - t1021 * t956 + t1036 * t957;
t1012 = -t1018 * t970 - t1021 * t957 + t1033 * t956;
t926 = mrSges(5,1) * t956 - mrSges(5,3) * t957;
t890 = -mrSges(6,1) * t921 + mrSges(6,2) * t922;
t899 = -mrSges(6,2) * t968 + mrSges(6,3) * t921;
t835 = m(6) * t848 + t951 * mrSges(6,1) + t968 * t899 + (-t889 - t890) * t922 + (-mrSges(6,3) - mrSges(7,3)) * t875 + t1005;
t1004 = m(7) * t846 + t874 * mrSges(7,3) + t921 * t889;
t901 = mrSges(7,1) * t968 - mrSges(7,3) * t922;
t902 = mrSges(6,1) * t968 - mrSges(6,3) * t922;
t837 = m(6) * t849 + t874 * mrSges(6,3) + t921 * t890 + (-t901 - t902) * t968 + (-mrSges(6,2) - mrSges(7,2)) * t951 + t1004;
t831 = t983 * t835 + t979 * t837;
t931 = -mrSges(5,2) * t956 + mrSges(5,3) * t970;
t993 = -m(5) * t862 + t955 * mrSges(5,1) + t970 * t931 - t831;
t829 = t919 * mrSges(5,2) + t957 * t926 - t993;
t930 = -mrSges(5,1) * t970 + mrSges(5,2) * t957;
t998 = -t979 * t835 + t983 * t837;
t996 = m(5) * t860 + t955 * mrSges(5,3) + t970 * t930 + t998;
t1028 = t1010 * t956 - t1012 * t957 + t1031 * t955 - t1018 * t918 + t1019 * t919 + mrSges(4,1) * t876 - mrSges(5,1) * t862 - mrSges(4,2) * t877 + mrSges(5,3) * t860 - pkin(3) * t829 - pkin(4) * t831 + qJ(4) * (-t918 * mrSges(5,2) - t956 * t926 + t996) - t1037;
t1011 = t1018 * t956 - t1019 * t957 - t1031 * t970;
t1014 = -t1030 * t968 - t1032 * t921 - t1034 * t922;
t1023 = pkin(3) * t970;
t932 = -t984 * g(3) - t981 * t947;
t896 = -qJDD(2) * pkin(2) - t986 * pkin(8) + t959 * t1008 - t932;
t995 = t918 * pkin(3) + t896 + (t1016 - t919) * qJ(4);
t858 = -pkin(4) * t918 - pkin(9) * t954 - t995 + (-t1023 + t1025 + t934) * t957;
t851 = -pkin(5) * t874 - qJ(6) * t920 + t900 * t922 + qJDD(6) + t858;
t841 = m(7) * t851 - t874 * mrSges(7,1) + t875 * mrSges(7,2) - t921 * t898 + t922 * t901;
t824 = -mrSges(6,1) * t858 + mrSges(6,3) * t849 - mrSges(7,1) * t851 + mrSges(7,3) * t846 - pkin(5) * t841 + qJ(6) * t1004 + (-qJ(6) * t901 + t1013) * t968 + (-mrSges(7,2) * qJ(6) + t1032) * t951 + t1014 * t922 + t1020 * t875 + t1038 * t874;
t830 = mrSges(6,2) * t858 + mrSges(7,2) * t851 - mrSges(6,3) * t848 - mrSges(7,3) * t844 - qJ(6) * t840 - t1014 * t921 + t1020 * t874 - t1029 * t968 + t1034 * t951 + t1035 * t875;
t861 = (-(2 * qJD(4)) + t1023) * t957 + t995;
t992 = -m(6) * t858 + t874 * mrSges(6,1) - t875 * mrSges(6,2) + t921 * t899 - t922 * t902 - t841;
t838 = m(5) * t861 + t918 * mrSges(5,1) - t919 * mrSges(5,3) - t957 * t930 + t956 * t931 + t992;
t809 = -mrSges(4,1) * t896 - mrSges(5,1) * t861 + mrSges(5,2) * t860 + mrSges(4,3) * t877 - pkin(3) * t838 - pkin(4) * t992 - pkin(9) * t998 + t1010 * t970 + t1011 * t957 + t1018 * t955 + t1021 * t919 - t1033 * t918 - t983 * t824 - t979 * t830;
t810 = mrSges(4,2) * t896 + mrSges(5,2) * t862 - mrSges(4,3) * t876 - mrSges(5,3) * t861 - pkin(9) * t831 - qJ(4) * t838 + t1011 * t956 + t1012 * t970 + t1019 * t955 - t1021 * t918 + t1036 * t919 - t979 * t824 + t983 * t830;
t1009 = -mrSges(4,1) * t956 - mrSges(4,2) * t957 - t926;
t1022 = -mrSges(4,3) - mrSges(5,2);
t929 = mrSges(4,1) * t970 - mrSges(4,3) * t957;
t826 = m(4) * t877 - t955 * mrSges(4,2) + t1009 * t956 + t1022 * t918 - t970 * t929 + t996;
t928 = -mrSges(4,2) * t970 - mrSges(4,3) * t956;
t827 = m(4) * t876 + t955 * mrSges(4,1) + t1009 * t957 + t1022 * t919 + t970 * t928 + t993;
t823 = t1024 * t826 - t827 * t980;
t834 = -m(4) * t896 - t918 * mrSges(4,1) - t919 * mrSges(4,2) - t956 * t928 - t957 * t929 - t838;
t944 = Ifges(3,6) * qJD(2) + (Ifges(3,4) * t981 + Ifges(3,2) * t984) * qJD(1);
t945 = Ifges(3,5) * qJD(2) + (Ifges(3,1) * t981 + Ifges(3,4) * t984) * qJD(1);
t1027 = mrSges(3,1) * t932 - mrSges(3,2) * t933 + Ifges(3,5) * t960 + Ifges(3,6) * t961 + Ifges(3,3) * qJDD(2) + pkin(2) * t834 + pkin(8) * t823 + (t944 * t981 - t945 * t984) * qJD(1) + t1024 * t809 + t980 * t810;
t958 = (-mrSges(3,1) * t984 + mrSges(3,2) * t981) * qJD(1);
t963 = qJD(2) * mrSges(3,1) - mrSges(3,3) * t1008;
t821 = m(3) * t933 - qJDD(2) * mrSges(3,2) + mrSges(3,3) * t961 - qJD(2) * t963 + t1007 * t958 + t823;
t964 = -qJD(2) * mrSges(3,2) + mrSges(3,3) * t1007;
t833 = m(3) * t932 + qJDD(2) * mrSges(3,1) - t960 * mrSges(3,3) + qJD(2) * t964 - t1008 * t958 + t834;
t999 = t984 * t821 - t833 * t981;
t813 = m(2) * t967 - mrSges(2,1) * t987 - qJDD(1) * mrSges(2,2) + t999;
t822 = t1024 * t827 + t980 * t826;
t990 = -m(3) * t946 + t961 * mrSges(3,1) - t960 * mrSges(3,2) + t964 * t1007 - t1008 * t963 - t822;
t817 = m(2) * t966 + qJDD(1) * mrSges(2,1) - t987 * mrSges(2,2) + t990;
t1015 = t982 * t813 + t985 * t817;
t815 = t981 * t821 + t984 * t833;
t1000 = t985 * t813 - t817 * t982;
t943 = Ifges(3,3) * qJD(2) + (Ifges(3,5) * t981 + Ifges(3,6) * t984) * qJD(1);
t806 = mrSges(3,2) * t946 - mrSges(3,3) * t932 + Ifges(3,1) * t960 + Ifges(3,4) * t961 + Ifges(3,5) * qJDD(2) - pkin(8) * t822 - qJD(2) * t944 + t1007 * t943 + t1024 * t810 - t980 * t809;
t808 = -mrSges(3,1) * t946 + mrSges(3,3) * t933 + Ifges(3,4) * t960 + Ifges(3,2) * t961 + Ifges(3,6) * qJDD(2) - pkin(2) * t822 + qJD(2) * t945 - t1008 * t943 - t1028;
t994 = mrSges(2,1) * t966 - mrSges(2,2) * t967 + Ifges(2,3) * qJDD(1) + pkin(1) * t990 + pkin(7) * t999 + t981 * t806 + t984 * t808;
t804 = mrSges(2,1) * g(3) + mrSges(2,3) * t967 + t987 * Ifges(2,5) + Ifges(2,6) * qJDD(1) - pkin(1) * t815 - t1027;
t803 = -mrSges(2,2) * g(3) - mrSges(2,3) * t966 + Ifges(2,5) * qJDD(1) - Ifges(2,6) * t987 - pkin(7) * t815 + t806 * t984 - t808 * t981;
t1 = [-m(1) * g(1) + t1000; -m(1) * g(2) + t1015; (-m(1) - m(2)) * g(3) + t815; -mrSges(1,2) * g(3) + mrSges(1,3) * g(2) - pkin(6) * t1015 + t985 * t803 - t982 * t804; mrSges(1,1) * g(3) - mrSges(1,3) * g(1) + pkin(6) * t1000 + t982 * t803 + t985 * t804; -mrSges(1,1) * g(2) + mrSges(1,2) * g(1) + t994; t994; t1027; t1028; t829; t1037; t841;];
tauJB  = t1;
