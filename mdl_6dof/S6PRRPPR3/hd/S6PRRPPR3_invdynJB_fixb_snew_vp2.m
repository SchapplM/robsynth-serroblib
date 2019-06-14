% Calculate vector of inverse dynamics joint torques and base forces with Newton-Euler
% S6PRRPPR3
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
% pkin [10x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,alpha2,d2,d3,d6,theta1]';
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
% Datum: 2019-05-05 03:01
% Revision: 89c353f7eff3bd693eda4e29f35b2761dbc3ada0 (2019-05-03)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauJB = S6PRRPPR3_invdynJB_fixb_snew_vp2(qJ, qJD, qJDD, g, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(6,1),zeros(3,1),zeros(10,1),zeros(7,1),zeros(7,3),zeros(7,6)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6PRRPPR3_invdynJB_fixb_snew_vp2: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6PRRPPR3_invdynJB_fixb_snew_vp2: qJD has to be [6x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [6 1]), ...
  'S6PRRPPR3_invdynJB_fixb_snew_vp2: qJDD has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6PRRPPR3_invdynJB_fixb_snew_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S6PRRPPR3_invdynJB_fixb_snew_vp2: pkin has to be [10x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6PRRPPR3_invdynJB_fixb_snew_vp2: m has to be [7x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [7,3]), ...
  'S6PRRPPR3_invdynJB_fixb_snew_vp2: mrSges has to be [7x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [7 6]), ...
  'S6PRRPPR3_invdynJB_fixb_snew_vp2: Ifges has to be [7x6] (double)'); 

%% Symbolic Calculation
% From invdyn_fixb_NewtonEuler_linkframe_tauJB_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-05-05 02:57:56
% EndTime: 2019-05-05 02:58:02
% DurationCPUTime: 4.14s
% Computational Cost: add. (43533->304), mult. (86337->366), div. (0->0), fcn. (47256->10), ass. (0->136)
t1014 = Ifges(4,1) + Ifges(5,1) + Ifges(6,2);
t986 = Ifges(4,4) - Ifges(5,5) + Ifges(6,4);
t985 = Ifges(4,5) + Ifges(5,4) + Ifges(6,6);
t1013 = Ifges(4,2) + Ifges(5,3) + Ifges(6,1);
t984 = Ifges(4,6) - Ifges(5,6) + Ifges(6,5);
t1012 = Ifges(4,3) + Ifges(5,2) + Ifges(6,3);
t961 = qJD(3) ^ 2;
t950 = sin(pkin(10));
t952 = cos(pkin(10));
t924 = g(1) * t950 - g(2) * t952;
t925 = -g(1) * t952 - g(2) * t950;
t948 = -g(3) + qJDD(1);
t960 = cos(qJ(2));
t953 = cos(pkin(6));
t957 = sin(qJ(2));
t994 = t953 * t957;
t951 = sin(pkin(6));
t996 = t951 * t957;
t865 = t924 * t994 + t960 * t925 + t948 * t996;
t962 = qJD(2) ^ 2;
t863 = -pkin(2) * t962 + qJDD(2) * pkin(8) + t865;
t879 = -t924 * t951 + t948 * t953;
t956 = sin(qJ(3));
t959 = cos(qJ(3));
t857 = -t956 * t863 + t959 * t879;
t914 = (-pkin(3) * t959 - qJ(4) * t956) * qJD(2);
t990 = qJD(2) * t956;
t974 = t914 * t990 + qJDD(4) - t857;
t855 = -qJDD(3) * pkin(3) - t961 * qJ(4) + t974;
t989 = qJD(2) * t959;
t932 = mrSges(5,2) * t989 + qJD(3) * mrSges(5,3);
t1008 = m(5) * t855 - qJDD(3) * mrSges(5,1) - qJD(3) * t932;
t915 = (-mrSges(5,1) * t959 - mrSges(5,3) * t956) * qJD(2);
t988 = qJD(2) * qJD(3);
t979 = t959 * t988;
t919 = qJDD(2) * t956 + t979;
t930 = -qJD(3) * mrSges(6,1) + mrSges(6,3) * t989;
t1002 = pkin(4) + pkin(9);
t1003 = -pkin(3) - pkin(9);
t978 = t956 * t988;
t920 = qJDD(2) * t959 - t978;
t1004 = 2 * qJD(4);
t926 = -qJD(3) * pkin(4) - qJ(5) * t990;
t993 = t953 * t960;
t995 = t951 * t960;
t864 = t924 * t993 - t957 * t925 + t948 * t995;
t862 = -qJDD(2) * pkin(2) - t962 * pkin(8) - t864;
t971 = -t920 * pkin(3) + t862 + (-t919 - t979) * qJ(4);
t997 = t959 ^ 2 * t962;
t966 = -qJ(5) * t997 + qJDD(5) - t971 + (t1004 + t926) * t990;
t845 = t966 + (pkin(5) * t959 + t1003 * t956) * t988 + t1002 * t920 + t919 * pkin(5);
t987 = qJD(2) * qJD(5);
t1009 = -0.2e1 * t956 * t987 + (-t919 + t979) * qJ(5);
t918 = (pkin(5) * t956 + pkin(9) * t959) * qJD(2);
t992 = t959 * t962;
t848 = (-pkin(5) - qJ(4)) * t961 + (-pkin(4) * t992 - qJD(2) * t918) * t956 + (-pkin(3) - t1002) * qJDD(3) + t974 + t1009;
t955 = sin(qJ(6));
t958 = cos(qJ(6));
t842 = t845 * t958 - t848 * t955;
t912 = -qJD(3) * t958 + t955 * t989;
t874 = qJD(6) * t912 - qJDD(3) * t955 - t920 * t958;
t913 = -qJD(3) * t955 - t958 * t989;
t875 = -mrSges(7,1) * t912 + mrSges(7,2) * t913;
t937 = qJD(6) + t990;
t877 = -mrSges(7,2) * t937 + mrSges(7,3) * t912;
t908 = qJDD(6) + t919;
t839 = m(7) * t842 + mrSges(7,1) * t908 - mrSges(7,3) * t874 - t875 * t913 + t877 * t937;
t843 = t845 * t955 + t848 * t958;
t873 = -qJD(6) * t913 - qJDD(3) * t958 + t920 * t955;
t878 = mrSges(7,1) * t937 - mrSges(7,3) * t913;
t840 = m(7) * t843 - mrSges(7,2) * t908 + mrSges(7,3) * t873 + t875 * t912 - t878 * t937;
t829 = -t955 * t839 + t958 * t840;
t850 = (-t956 * t992 - qJDD(3)) * pkin(4) + t855 + t1009;
t917 = (mrSges(6,1) * t956 - mrSges(6,2) * t959) * qJD(2);
t973 = -m(6) * t850 + t917 * t990 - t829;
t968 = qJDD(3) * mrSges(6,2) + qJD(3) * t930 - t973;
t825 = t915 * t990 + (mrSges(5,2) - mrSges(6,3)) * t919 + t968 + t1008;
t827 = -t919 * mrSges(6,3) + t968;
t942 = qJD(3) * t1004;
t858 = t959 * t863 + t956 * t879;
t1010 = qJDD(3) * qJ(4) + t914 * t989 + t858;
t970 = pkin(4) * t997 + t920 * qJ(5) - t1010;
t847 = qJDD(3) * pkin(5) + qJD(3) * t926 + t942 + t1003 * t961 + (-0.2e1 * qJD(5) - t918) * t989 - t970;
t866 = Ifges(7,5) * t913 + Ifges(7,6) * t912 + Ifges(7,3) * t937;
t868 = Ifges(7,1) * t913 + Ifges(7,4) * t912 + Ifges(7,5) * t937;
t833 = -mrSges(7,1) * t847 + mrSges(7,3) * t843 + Ifges(7,4) * t874 + Ifges(7,2) * t873 + Ifges(7,6) * t908 - t866 * t913 + t868 * t937;
t867 = Ifges(7,4) * t913 + Ifges(7,2) * t912 + Ifges(7,6) * t937;
t834 = mrSges(7,2) * t847 - mrSges(7,3) * t842 + Ifges(7,1) * t874 + Ifges(7,4) * t873 + Ifges(7,5) * t908 + t866 * t912 - t867 * t937;
t844 = -m(7) * t847 + t873 * mrSges(7,1) - t874 * mrSges(7,2) + t912 * t877 - t913 * t878;
t1000 = t961 * pkin(3);
t1005 = -2 * qJD(4);
t849 = 0.2e1 * t959 * t987 + t1000 + (t1005 - t926) * qJD(3) + t970;
t854 = t942 - t1000 + t1010;
t929 = -qJD(3) * mrSges(5,1) + mrSges(5,2) * t990;
t927 = qJD(3) * mrSges(6,2) - mrSges(6,3) * t990;
t967 = -m(6) * t849 + qJDD(3) * mrSges(6,1) - t920 * mrSges(6,3) + qJD(3) * t927 - t844;
t965 = m(5) * t854 + qJDD(3) * mrSges(5,3) + qJD(3) * t929 + t915 * t989 + t967;
t980 = -t985 * qJD(3) + (-t1014 * t956 - t986 * t959) * qJD(2);
t981 = -t984 * qJD(3) + (-t1013 * t959 - t986 * t956) * qJD(2);
t1011 = -(t981 * t956 - t980 * t959) * qJD(2) + t1012 * qJDD(3) + t985 * t919 + t984 * t920 + mrSges(4,1) * t857 - mrSges(5,1) * t855 - mrSges(6,1) * t849 - mrSges(4,2) * t858 + mrSges(6,2) * t850 + mrSges(5,3) * t854 - pkin(3) * t825 - pkin(4) * t827 - pkin(5) * t844 - pkin(9) * t829 + qJ(4) * (t920 * mrSges(5,2) - t917 * t989 + t965) - t958 * t833 - t955 * t834;
t999 = mrSges(4,3) + mrSges(5,2);
t916 = (-mrSges(4,1) * t959 + mrSges(4,2) * t956) * qJD(2);
t931 = -qJD(3) * mrSges(4,2) + mrSges(4,3) * t989;
t823 = m(4) * t857 + (mrSges(4,1) - mrSges(6,2)) * qJDD(3) + (-t930 + t931) * qJD(3) + (-t915 - t916) * t990 + (mrSges(6,3) - t999) * t919 + t973 - t1008;
t928 = qJD(3) * mrSges(4,1) - mrSges(4,3) * t990;
t832 = t999 * t920 + (t916 - t917) * t989 + t965 - qJDD(3) * mrSges(4,2) + m(4) * t858 - qJD(3) * t928;
t976 = -t823 * t956 + t959 * t832;
t815 = m(3) * t865 - mrSges(3,1) * t962 - qJDD(2) * mrSges(3,2) + t976;
t818 = t959 * t823 + t956 * t832;
t817 = m(3) * t879 + t818;
t828 = t958 * t839 + t955 * t840;
t852 = -pkin(3) * t978 + t920 * pkin(4) + t966;
t826 = m(6) * t852 + t919 * mrSges(6,1) - t920 * mrSges(6,2) + t927 * t990 - t930 * t989 + t828;
t856 = (pkin(3) * qJD(3) + t1005) * t990 + t971;
t824 = m(5) * t856 - mrSges(5,1) * t920 - t919 * mrSges(5,3) - t929 * t990 - t932 * t989 - t826;
t964 = -m(4) * t862 + t920 * mrSges(4,1) - mrSges(4,2) * t919 - t928 * t990 + t931 * t989 - t824;
t821 = m(3) * t864 + qJDD(2) * mrSges(3,1) - mrSges(3,2) * t962 + t964;
t806 = t815 * t994 - t817 * t951 + t821 * t993;
t804 = m(2) * t924 + t806;
t811 = t960 * t815 - t821 * t957;
t810 = m(2) * t925 + t811;
t991 = t952 * t804 + t950 * t810;
t805 = t815 * t996 + t953 * t817 + t821 * t995;
t982 = -t1012 * qJD(3) + (-t985 * t956 - t984 * t959) * qJD(2);
t977 = -t804 * t950 + t952 * t810;
t975 = m(2) * t948 + t805;
t802 = pkin(4) * t826 - pkin(3) * t824 + pkin(9) * t828 + mrSges(6,3) * t849 - mrSges(6,2) * t852 + mrSges(5,2) * t854 - mrSges(5,1) * t856 + mrSges(4,3) * t858 - mrSges(4,1) * t862 + t955 * t833 - t958 * t834 - qJ(5) * t967 + t1013 * t920 + t986 * t919 + t984 * qJDD(3) - t980 * qJD(3) + (qJ(5) * t917 * t959 + t982 * t956) * qJD(2);
t969 = mrSges(7,1) * t842 - mrSges(7,2) * t843 + Ifges(7,5) * t874 + Ifges(7,6) * t873 + Ifges(7,3) * t908 + t913 * t867 - t912 * t868;
t807 = mrSges(6,1) * t852 + mrSges(4,2) * t862 + mrSges(5,2) * t855 - mrSges(4,3) * t857 - mrSges(5,3) * t856 - mrSges(6,3) * t850 + pkin(5) * t828 - qJ(4) * t824 - qJ(5) * t827 + t981 * qJD(3) + t985 * qJDD(3) + t1014 * t919 + t986 * t920 - t982 * t989 + t969;
t800 = mrSges(3,2) * t879 - mrSges(3,3) * t864 + Ifges(3,5) * qJDD(2) - Ifges(3,6) * t962 - pkin(8) * t818 - t802 * t956 + t807 * t959;
t801 = -mrSges(3,1) * t879 + mrSges(3,3) * t865 + t962 * Ifges(3,5) + Ifges(3,6) * qJDD(2) - pkin(2) * t818 - t1011;
t972 = pkin(7) * t811 + t800 * t957 + t801 * t960;
t799 = mrSges(3,1) * t864 - mrSges(3,2) * t865 + Ifges(3,3) * qJDD(2) + pkin(2) * t964 + pkin(8) * t976 + t959 * t802 + t956 * t807;
t798 = mrSges(2,2) * t948 - mrSges(2,3) * t924 + t960 * t800 - t957 * t801 + (-t805 * t951 - t806 * t953) * pkin(7);
t797 = -mrSges(2,1) * t948 + mrSges(2,3) * t925 - pkin(1) * t805 - t951 * t799 + t972 * t953;
t1 = [-m(1) * g(1) + t977; -m(1) * g(2) + t991; -m(1) * g(3) + t975; -mrSges(1,2) * g(3) + mrSges(1,3) * g(2) - qJ(1) * t991 - t950 * t797 + t952 * t798; mrSges(1,1) * g(3) - mrSges(1,3) * g(1) + qJ(1) * t977 + t952 * t797 + t950 * t798; -mrSges(1,1) * g(2) + mrSges(2,1) * t924 + mrSges(1,2) * g(1) - mrSges(2,2) * t925 + pkin(1) * t806 + t953 * t799 + t972 * t951; t975; t799; t1011; t825; t826; t969;];
tauJB  = t1;
