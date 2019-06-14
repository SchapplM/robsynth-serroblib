% Calculate vector of inverse dynamics joint torques and base forces with Newton-Euler
% S6RRPPRR6
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
%   pkin=[a2,a3,a4,a5,a6,d1,d2,d5,d6,theta4]';
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
% Datum: 2019-05-06 10:56
% Revision: 89c353f7eff3bd693eda4e29f35b2761dbc3ada0 (2019-05-03)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauJB = S6RRPPRR6_invdynJB_fixb_snew_vp2(qJ, qJD, qJDD, g, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(6,1),zeros(3,1),zeros(10,1),zeros(7,1),zeros(7,3),zeros(7,6)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRPPRR6_invdynJB_fixb_snew_vp2: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RRPPRR6_invdynJB_fixb_snew_vp2: qJD has to be [6x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [6 1]), ...
  'S6RRPPRR6_invdynJB_fixb_snew_vp2: qJDD has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6RRPPRR6_invdynJB_fixb_snew_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S6RRPPRR6_invdynJB_fixb_snew_vp2: pkin has to be [10x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6RRPPRR6_invdynJB_fixb_snew_vp2: m has to be [7x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [7,3]), ...
  'S6RRPPRR6_invdynJB_fixb_snew_vp2: mrSges has to be [7x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [7 6]), ...
  'S6RRPPRR6_invdynJB_fixb_snew_vp2: Ifges has to be [7x6] (double)'); 

%% Symbolic Calculation
% From invdyn_fixb_NewtonEuler_linkframe_tauJB_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-05-06 10:53:16
% EndTime: 2019-05-06 10:53:27
% DurationCPUTime: 10.25s
% Computational Cost: add. (149199->364), mult. (340780->447), div. (0->0), fcn. (225279->10), ass. (0->145)
t1023 = Ifges(3,1) + Ifges(4,1);
t1015 = Ifges(3,4) - Ifges(4,5);
t1014 = Ifges(3,5) + Ifges(4,4);
t1022 = Ifges(3,2) + Ifges(4,3);
t1013 = Ifges(3,6) - Ifges(4,6);
t1021 = Ifges(3,3) + Ifges(4,2);
t979 = sin(qJ(2));
t983 = cos(qJ(2));
t1008 = t1014 * qJD(2) + (t1015 * t983 + t1023 * t979) * qJD(1);
t1009 = -t1013 * qJD(2) + (-t1015 * t979 - t1022 * t983) * qJD(1);
t1007 = qJD(1) * t979;
t943 = (-mrSges(4,1) * t983 - mrSges(4,3) * t979) * qJD(1);
t1005 = qJD(1) * qJD(2);
t1002 = t983 * t1005;
t945 = t979 * qJDD(1) + t1002;
t986 = qJD(1) ^ 2;
t1012 = t983 ^ 2 * t986;
t1006 = qJD(1) * t983;
t1017 = 2 * qJD(3);
t980 = sin(qJ(1));
t984 = cos(qJ(1));
t956 = -t984 * g(1) - t980 * g(2);
t933 = -t986 * pkin(1) + qJDD(1) * pkin(7) + t956;
t912 = -t979 * g(3) + t983 * t933;
t942 = (-pkin(2) * t983 - qJ(3) * t979) * qJD(1);
t985 = qJD(2) ^ 2;
t891 = -t985 * pkin(2) + qJDD(2) * qJ(3) + qJD(2) * t1017 + t942 * t1006 + t912;
t1003 = t979 * t1005;
t946 = t983 * qJDD(1) - t1003;
t950 = -qJD(2) * pkin(3) - qJ(4) * t1007;
t887 = -pkin(3) * t1012 - t946 * qJ(4) + qJD(2) * t950 + t891;
t911 = -t983 * g(3) - t979 * t933;
t894 = -qJDD(2) * pkin(2) - t985 * qJ(3) + t942 * t1007 + qJDD(3) - t911;
t888 = (-t945 + t1002) * qJ(4) + (-t979 * t983 * t986 - qJDD(2)) * pkin(3) + t894;
t974 = sin(pkin(10));
t975 = cos(pkin(10));
t925 = (-t974 * t983 + t975 * t979) * qJD(1);
t853 = -0.2e1 * qJD(4) * t925 - t974 * t887 + t975 * t888;
t910 = t975 * t945 - t974 * t946;
t924 = (-t974 * t979 - t975 * t983) * qJD(1);
t850 = (-qJD(2) * t924 - t910) * pkin(8) + (t924 * t925 - qJDD(2)) * pkin(4) + t853;
t854 = 0.2e1 * qJD(4) * t924 + t975 * t887 + t974 * t888;
t909 = -t974 * t945 - t975 * t946;
t915 = -qJD(2) * pkin(4) - pkin(8) * t925;
t923 = t924 ^ 2;
t852 = -pkin(4) * t923 + pkin(8) * t909 + qJD(2) * t915 + t854;
t978 = sin(qJ(5));
t982 = cos(qJ(5));
t847 = t978 * t850 + t982 * t852;
t902 = t978 * t924 + t982 * t925;
t870 = -t902 * qJD(5) + t982 * t909 - t978 * t910;
t901 = t982 * t924 - t978 * t925;
t885 = -mrSges(6,1) * t901 + mrSges(6,2) * t902;
t967 = -qJD(2) + qJD(5);
t896 = t967 * mrSges(6,1) - t902 * mrSges(6,3);
t966 = -qJDD(2) + qJDD(5);
t886 = -pkin(5) * t901 - pkin(9) * t902;
t965 = t967 ^ 2;
t844 = -t965 * pkin(5) + t966 * pkin(9) + t901 * t886 + t847;
t955 = t980 * g(1) - t984 * g(2);
t932 = -qJDD(1) * pkin(1) - t986 * pkin(7) - t955;
t995 = -t946 * pkin(2) + t932 + (-t1002 - t945) * qJ(3);
t874 = -pkin(2) * t1003 + t946 * pkin(3) - qJ(4) * t1012 + qJDD(4) - t995 + (t1017 + t950) * t1007;
t856 = -t909 * pkin(4) - t923 * pkin(8) + t925 * t915 + t874;
t871 = t901 * qJD(5) + t978 * t909 + t982 * t910;
t848 = t856 + (t902 * t967 - t870) * pkin(5) + (-t901 * t967 - t871) * pkin(9);
t977 = sin(qJ(6));
t981 = cos(qJ(6));
t841 = -t977 * t844 + t981 * t848;
t892 = -t977 * t902 + t981 * t967;
t859 = t892 * qJD(6) + t981 * t871 + t977 * t966;
t869 = qJDD(6) - t870;
t893 = t981 * t902 + t977 * t967;
t872 = -mrSges(7,1) * t892 + mrSges(7,2) * t893;
t897 = qJD(6) - t901;
t875 = -mrSges(7,2) * t897 + mrSges(7,3) * t892;
t837 = m(7) * t841 + mrSges(7,1) * t869 - mrSges(7,3) * t859 - t872 * t893 + t875 * t897;
t842 = t981 * t844 + t977 * t848;
t858 = -t893 * qJD(6) - t977 * t871 + t981 * t966;
t876 = mrSges(7,1) * t897 - mrSges(7,3) * t893;
t838 = m(7) * t842 - mrSges(7,2) * t869 + mrSges(7,3) * t858 + t872 * t892 - t876 * t897;
t997 = -t977 * t837 + t981 * t838;
t823 = m(6) * t847 - t966 * mrSges(6,2) + t870 * mrSges(6,3) + t901 * t885 - t967 * t896 + t997;
t846 = t982 * t850 - t978 * t852;
t895 = -t967 * mrSges(6,2) + t901 * mrSges(6,3);
t843 = -t966 * pkin(5) - t965 * pkin(9) + t902 * t886 - t846;
t992 = -m(7) * t843 + t858 * mrSges(7,1) - t859 * mrSges(7,2) + t892 * t875 - t893 * t876;
t833 = m(6) * t846 + t966 * mrSges(6,1) - t871 * mrSges(6,3) - t902 * t885 + t967 * t895 + t992;
t816 = t978 * t823 + t982 * t833;
t906 = -mrSges(5,1) * t924 + mrSges(5,2) * t925;
t913 = qJD(2) * mrSges(5,2) + mrSges(5,3) * t924;
t814 = m(5) * t853 - qJDD(2) * mrSges(5,1) - mrSges(5,3) * t910 - qJD(2) * t913 - t906 * t925 + t816;
t914 = -qJD(2) * mrSges(5,1) - mrSges(5,3) * t925;
t998 = t982 * t823 - t978 * t833;
t815 = m(5) * t854 + qJDD(2) * mrSges(5,2) + t909 * mrSges(5,3) + qJD(2) * t914 + t924 * t906 + t998;
t810 = t975 * t814 + t974 * t815;
t954 = mrSges(4,2) * t1006 + qJD(2) * mrSges(4,3);
t991 = -m(4) * t894 + qJDD(2) * mrSges(4,1) + qJD(2) * t954 - t810;
t809 = t945 * mrSges(4,2) + t943 * t1007 - t991;
t899 = Ifges(5,4) * t925 + Ifges(5,2) * t924 - Ifges(5,6) * qJD(2);
t900 = Ifges(5,1) * t925 + Ifges(5,4) * t924 - Ifges(5,5) * qJD(2);
t860 = Ifges(7,5) * t893 + Ifges(7,6) * t892 + Ifges(7,3) * t897;
t862 = Ifges(7,1) * t893 + Ifges(7,4) * t892 + Ifges(7,5) * t897;
t830 = -mrSges(7,1) * t843 + mrSges(7,3) * t842 + Ifges(7,4) * t859 + Ifges(7,2) * t858 + Ifges(7,6) * t869 - t860 * t893 + t862 * t897;
t861 = Ifges(7,4) * t893 + Ifges(7,2) * t892 + Ifges(7,6) * t897;
t831 = mrSges(7,2) * t843 - mrSges(7,3) * t841 + Ifges(7,1) * t859 + Ifges(7,4) * t858 + Ifges(7,5) * t869 + t860 * t892 - t861 * t897;
t878 = Ifges(6,4) * t902 + Ifges(6,2) * t901 + Ifges(6,6) * t967;
t879 = Ifges(6,1) * t902 + Ifges(6,4) * t901 + Ifges(6,5) * t967;
t990 = -mrSges(6,1) * t846 + mrSges(6,2) * t847 - Ifges(6,5) * t871 - Ifges(6,6) * t870 - Ifges(6,3) * t966 - pkin(5) * t992 - pkin(9) * t997 - t981 * t830 - t977 * t831 - t902 * t878 + t901 * t879;
t952 = -qJD(2) * mrSges(4,1) + mrSges(4,2) * t1007;
t999 = -t974 * t814 + t975 * t815;
t994 = m(4) * t891 + qJDD(2) * mrSges(4,3) + qJD(2) * t952 + t943 * t1006 + t999;
t1020 = -(t1008 * t983 + t1009 * t979) * qJD(1) + (Ifges(5,3) + t1021) * qJDD(2) + t1013 * t946 + t1014 * t945 + mrSges(3,1) * t911 - mrSges(4,1) * t894 - mrSges(5,1) * t853 - mrSges(3,2) * t912 + mrSges(5,2) * t854 + mrSges(4,3) * t891 - Ifges(5,5) * t910 - Ifges(5,6) * t909 - pkin(2) * t809 - pkin(3) * t810 - pkin(4) * t816 + qJ(3) * (t946 * mrSges(4,2) + t994) - t925 * t899 + t924 * t900 + t990;
t1016 = mrSges(3,3) + mrSges(4,2);
t944 = (-mrSges(3,1) * t983 + mrSges(3,2) * t979) * qJD(1);
t951 = qJD(2) * mrSges(3,1) - mrSges(3,3) * t1007;
t806 = m(3) * t912 - qJDD(2) * mrSges(3,2) - qJD(2) * t951 + t944 * t1006 + t1016 * t946 + t994;
t953 = -qJD(2) * mrSges(3,2) + mrSges(3,3) * t1006;
t807 = m(3) * t911 + qJDD(2) * mrSges(3,1) + qJD(2) * t953 - t1016 * t945 + (-t943 - t944) * t1007 + t991;
t1000 = t983 * t806 - t979 * t807;
t799 = m(2) * t956 - t986 * mrSges(2,1) - qJDD(1) * mrSges(2,2) + t1000;
t826 = t981 * t837 + t977 * t838;
t996 = m(6) * t856 - t870 * mrSges(6,1) + t871 * mrSges(6,2) - t901 * t895 + t902 * t896 + t826;
t824 = m(5) * t874 - t909 * mrSges(5,1) + t910 * mrSges(5,2) - t924 * t913 + t925 * t914 + t996;
t889 = (pkin(2) * qJD(2) - (2 * qJD(3))) * t1007 + t995;
t820 = m(4) * t889 - t946 * mrSges(4,1) - t945 * mrSges(4,3) - t954 * t1006 - t952 * t1007 - t824;
t988 = -m(3) * t932 + t946 * mrSges(3,1) - t945 * mrSges(3,2) + t953 * t1006 - t951 * t1007 - t820;
t818 = m(2) * t955 + qJDD(1) * mrSges(2,1) - t986 * mrSges(2,2) + t988;
t1011 = t980 * t799 + t984 * t818;
t801 = t979 * t806 + t983 * t807;
t1010 = t1021 * qJD(2) + (t1013 * t983 + t1014 * t979) * qJD(1);
t1001 = t984 * t799 - t980 * t818;
t877 = Ifges(6,5) * t902 + Ifges(6,6) * t901 + Ifges(6,3) * t967;
t811 = mrSges(6,2) * t856 - mrSges(6,3) * t846 + Ifges(6,1) * t871 + Ifges(6,4) * t870 + Ifges(6,5) * t966 - pkin(9) * t826 - t977 * t830 + t981 * t831 + t901 * t877 - t967 * t878;
t989 = mrSges(7,1) * t841 - mrSges(7,2) * t842 + Ifges(7,5) * t859 + Ifges(7,6) * t858 + Ifges(7,3) * t869 + t893 * t861 - t892 * t862;
t812 = -mrSges(6,1) * t856 + mrSges(6,3) * t847 + Ifges(6,4) * t871 + Ifges(6,2) * t870 + Ifges(6,6) * t966 - pkin(5) * t826 - t902 * t877 + t967 * t879 - t989;
t898 = Ifges(5,5) * t925 + Ifges(5,6) * t924 - Ifges(5,3) * qJD(2);
t796 = -mrSges(5,1) * t874 + mrSges(5,3) * t854 + Ifges(5,4) * t910 + Ifges(5,2) * t909 - Ifges(5,6) * qJDD(2) - pkin(4) * t996 + pkin(8) * t998 - qJD(2) * t900 + t978 * t811 + t982 * t812 - t925 * t898;
t802 = mrSges(5,2) * t874 - mrSges(5,3) * t853 + Ifges(5,1) * t910 + Ifges(5,4) * t909 - Ifges(5,5) * qJDD(2) - pkin(8) * t816 + qJD(2) * t899 + t982 * t811 - t978 * t812 + t924 * t898;
t793 = -mrSges(3,1) * t932 - mrSges(4,1) * t889 + mrSges(4,2) * t891 + mrSges(3,3) * t912 - pkin(2) * t820 + pkin(3) * t824 - qJ(4) * t999 + t1008 * qJD(2) + t1013 * qJDD(2) - t1010 * t1007 + t1015 * t945 + t1022 * t946 - t975 * t796 - t974 * t802;
t795 = mrSges(3,2) * t932 + mrSges(4,2) * t894 - mrSges(3,3) * t911 - mrSges(4,3) * t889 - qJ(3) * t820 - qJ(4) * t810 + t1009 * qJD(2) + t1014 * qJDD(2) + t1010 * t1006 + t1015 * t946 + t1023 * t945 - t974 * t796 + t975 * t802;
t993 = mrSges(2,1) * t955 - mrSges(2,2) * t956 + Ifges(2,3) * qJDD(1) + pkin(1) * t988 + pkin(7) * t1000 + t983 * t793 + t979 * t795;
t791 = mrSges(2,1) * g(3) + mrSges(2,3) * t956 + t986 * Ifges(2,5) + Ifges(2,6) * qJDD(1) - pkin(1) * t801 - t1020;
t790 = -mrSges(2,2) * g(3) - mrSges(2,3) * t955 + Ifges(2,5) * qJDD(1) - t986 * Ifges(2,6) - pkin(7) * t801 - t979 * t793 + t983 * t795;
t1 = [-m(1) * g(1) + t1001; -m(1) * g(2) + t1011; (-m(1) - m(2)) * g(3) + t801; -mrSges(1,2) * g(3) + mrSges(1,3) * g(2) - pkin(6) * t1011 + t984 * t790 - t980 * t791; mrSges(1,1) * g(3) - mrSges(1,3) * g(1) + pkin(6) * t1001 + t980 * t790 + t984 * t791; -mrSges(1,1) * g(2) + mrSges(1,2) * g(1) + t993; t993; t1020; t809; t824; -t990; t989;];
tauJB  = t1;
