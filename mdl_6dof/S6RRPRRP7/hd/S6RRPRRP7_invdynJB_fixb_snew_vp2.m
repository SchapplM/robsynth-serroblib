% Calculate vector of inverse dynamics joint torques and base forces with Newton-Euler
% S6RRPRRP7
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
%   pkin=[a2,a3,a4,a5,a6,d1,d2,d4,d5]';
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
% Datum: 2019-05-06 18:17
% Revision: 89c353f7eff3bd693eda4e29f35b2761dbc3ada0 (2019-05-03)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauJB = S6RRPRRP7_invdynJB_fixb_snew_vp2(qJ, qJD, qJDD, g, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(6,1),zeros(3,1),zeros(9,1),zeros(7,1),zeros(7,3),zeros(7,6)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRPRRP7_invdynJB_fixb_snew_vp2: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RRPRRP7_invdynJB_fixb_snew_vp2: qJD has to be [6x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [6 1]), ...
  'S6RRPRRP7_invdynJB_fixb_snew_vp2: qJDD has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6RRPRRP7_invdynJB_fixb_snew_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S6RRPRRP7_invdynJB_fixb_snew_vp2: pkin has to be [9x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6RRPRRP7_invdynJB_fixb_snew_vp2: m has to be [7x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [7,3]), ...
  'S6RRPRRP7_invdynJB_fixb_snew_vp2: mrSges has to be [7x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [7 6]), ...
  'S6RRPRRP7_invdynJB_fixb_snew_vp2: Ifges has to be [7x6] (double)'); 

%% Symbolic Calculation
% From invdyn_fixb_NewtonEuler_linkframe_tauJB_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-05-06 18:13:09
% EndTime: 2019-05-06 18:13:20
% DurationCPUTime: 6.82s
% Computational Cost: add. (74675->341), mult. (152003->404), div. (0->0), fcn. (94246->8), ass. (0->138)
t1027 = Ifges(3,1) + Ifges(4,1);
t1026 = Ifges(6,1) + Ifges(7,1);
t1013 = Ifges(3,4) - Ifges(4,5);
t1012 = Ifges(6,4) - Ifges(7,5);
t1011 = Ifges(3,5) + Ifges(4,4);
t1010 = -Ifges(6,5) - Ifges(7,4);
t1025 = Ifges(3,2) + Ifges(4,3);
t1024 = Ifges(6,2) + Ifges(7,3);
t1009 = Ifges(3,6) - Ifges(4,6);
t1008 = Ifges(6,6) - Ifges(7,6);
t1023 = Ifges(3,3) + Ifges(4,2);
t1022 = -Ifges(6,3) - Ifges(7,2);
t968 = sin(qJ(2));
t971 = cos(qJ(2));
t932 = (-mrSges(4,1) * t971 - mrSges(4,3) * t968) * qJD(1);
t994 = qJD(1) * qJD(2);
t992 = t971 * t994;
t934 = qJDD(1) * t968 + t992;
t974 = qJD(1) ^ 2;
t1005 = t971 ^ 2 * t974;
t1017 = 2 * qJD(3);
t969 = sin(qJ(1));
t972 = cos(qJ(1));
t945 = -g(1) * t972 - g(2) * t969;
t922 = -pkin(1) * t974 + qJDD(1) * pkin(7) + t945;
t903 = -g(3) * t968 + t971 * t922;
t931 = (-pkin(2) * t971 - qJ(3) * t968) * qJD(1);
t973 = qJD(2) ^ 2;
t995 = qJD(1) * t971;
t880 = -pkin(2) * t973 + qJDD(2) * qJ(3) + qJD(2) * t1017 + t931 * t995 + t903;
t991 = t968 * t994;
t935 = qJDD(1) * t971 - t991;
t996 = qJD(1) * t968;
t943 = -qJD(2) * pkin(3) - pkin(8) * t996;
t859 = -pkin(3) * t1005 - pkin(8) * t935 + qJD(2) * t943 + t880;
t902 = -t971 * g(3) - t968 * t922;
t886 = -qJDD(2) * pkin(2) - t973 * qJ(3) + t931 * t996 + qJDD(3) - t902;
t860 = (-t934 + t992) * pkin(8) + (-t968 * t971 * t974 - qJDD(2)) * pkin(3) + t886;
t967 = sin(qJ(4));
t970 = cos(qJ(4));
t849 = t970 * t859 + t967 * t860;
t920 = (-t967 * t971 + t968 * t970) * qJD(1);
t887 = -qJD(4) * t920 - t934 * t967 - t935 * t970;
t919 = (t967 * t968 + t970 * t971) * qJD(1);
t898 = mrSges(5,1) * t919 + mrSges(5,2) * t920;
t958 = -qJD(2) + qJD(4);
t905 = mrSges(5,1) * t958 - mrSges(5,3) * t920;
t957 = -qJDD(2) + qJDD(4);
t1016 = cos(qJ(5));
t966 = sin(qJ(5));
t900 = -t1016 * t958 + t966 * t920;
t901 = t1016 * t920 + t966 * t958;
t873 = mrSges(7,1) * t900 - mrSges(7,3) * t901;
t1000 = -mrSges(6,1) * t900 - mrSges(6,2) * t901 - t873;
t1014 = -mrSges(6,3) - mrSges(7,2);
t944 = t969 * g(1) - t972 * g(2);
t921 = -qJDD(1) * pkin(1) - t974 * pkin(7) - t944;
t983 = -t935 * pkin(2) + t921 + (-t934 - t992) * qJ(3);
t851 = -pkin(2) * t991 + t935 * pkin(3) - pkin(8) * t1005 - t983 + (t1017 + t943) * t996;
t888 = -qJD(4) * t919 + t934 * t970 - t935 * t967;
t843 = t851 + (t919 * t958 - t888) * pkin(9) + (t920 * t958 - t887) * pkin(4);
t899 = pkin(4) * t919 - pkin(9) * t920;
t956 = t958 ^ 2;
t846 = -pkin(4) * t956 + pkin(9) * t957 - t899 * t919 + t849;
t841 = t1016 * t846 + t966 * t843;
t854 = t901 * qJD(5) - t1016 * t957 + t966 * t888;
t885 = qJDD(5) - t887;
t912 = qJD(5) + t919;
t891 = mrSges(6,1) * t912 - mrSges(6,3) * t901;
t872 = pkin(5) * t900 - qJ(6) * t901;
t908 = t912 ^ 2;
t837 = -pkin(5) * t908 + t885 * qJ(6) + 0.2e1 * qJD(6) * t912 - t872 * t900 + t841;
t892 = -mrSges(7,1) * t912 + mrSges(7,2) * t901;
t993 = m(7) * t837 + t885 * mrSges(7,3) + t912 * t892;
t828 = m(6) * t841 - t885 * mrSges(6,2) + t1000 * t900 + t1014 * t854 - t912 * t891 + t993;
t840 = t1016 * t843 - t966 * t846;
t855 = -t900 * qJD(5) + t1016 * t888 + t966 * t957;
t890 = -mrSges(6,2) * t912 - mrSges(6,3) * t900;
t838 = -t885 * pkin(5) - t908 * qJ(6) + t901 * t872 + qJDD(6) - t840;
t889 = -mrSges(7,2) * t900 + mrSges(7,3) * t912;
t986 = -m(7) * t838 + t885 * mrSges(7,1) + t912 * t889;
t830 = m(6) * t840 + t885 * mrSges(6,1) + t1000 * t901 + t1014 * t855 + t912 * t890 + t986;
t987 = t1016 * t828 - t830 * t966;
t817 = m(5) * t849 - mrSges(5,2) * t957 + mrSges(5,3) * t887 - t898 * t919 - t905 * t958 + t987;
t848 = -t967 * t859 + t970 * t860;
t904 = -mrSges(5,2) * t958 - mrSges(5,3) * t919;
t845 = -t957 * pkin(4) - t956 * pkin(9) + t920 * t899 - t848;
t839 = -0.2e1 * qJD(6) * t901 + (t900 * t912 - t855) * qJ(6) + (t901 * t912 + t854) * pkin(5) + t845;
t835 = m(7) * t839 + t854 * mrSges(7,1) - t855 * mrSges(7,3) + t889 * t900 - t901 * t892;
t977 = -m(6) * t845 - t854 * mrSges(6,1) - t855 * mrSges(6,2) - t900 * t890 - t891 * t901 - t835;
t825 = m(5) * t848 + mrSges(5,1) * t957 - mrSges(5,3) * t888 - t898 * t920 + t904 * t958 + t977;
t811 = t967 * t817 + t970 * t825;
t942 = mrSges(4,2) * t995 + qJD(2) * mrSges(4,3);
t980 = -m(4) * t886 + qJDD(2) * mrSges(4,1) + qJD(2) * t942 - t811;
t810 = t934 * mrSges(4,2) + t932 * t996 - t980;
t1001 = t1010 * t912 + t1012 * t900 - t1026 * t901;
t1003 = t1008 * t900 + t1010 * t901 + t1022 * t912;
t819 = -mrSges(6,1) * t845 - mrSges(7,1) * t839 + mrSges(7,2) * t837 + mrSges(6,3) * t841 - pkin(5) * t835 - t1001 * t912 + t1003 * t901 + t1008 * t885 + t1012 * t855 - t1024 * t854;
t1002 = -t1008 * t912 - t1012 * t901 + t1024 * t900;
t821 = mrSges(6,2) * t845 + mrSges(7,2) * t838 - mrSges(6,3) * t840 - mrSges(7,3) * t839 - qJ(6) * t835 + t1002 * t912 + t1003 * t900 - t1010 * t885 - t1012 * t854 + t1026 * t855;
t894 = Ifges(5,4) * t920 - Ifges(5,2) * t919 + Ifges(5,6) * t958;
t895 = Ifges(5,1) * t920 - Ifges(5,4) * t919 + Ifges(5,5) * t958;
t979 = -mrSges(5,1) * t848 + mrSges(5,2) * t849 - Ifges(5,5) * t888 - Ifges(5,6) * t887 - Ifges(5,3) * t957 - pkin(4) * t977 - pkin(9) * t987 - t1016 * t819 - t966 * t821 - t920 * t894 - t919 * t895;
t940 = -qJD(2) * mrSges(4,1) + mrSges(4,2) * t996;
t988 = t970 * t817 - t967 * t825;
t982 = m(4) * t880 + qJDD(2) * mrSges(4,3) + qJD(2) * t940 + t932 * t995 + t988;
t997 = t1011 * qJD(2) + (t1013 * t971 + t1027 * t968) * qJD(1);
t998 = -t1009 * qJD(2) + (-t1013 * t968 - t1025 * t971) * qJD(1);
t1021 = -(t968 * t998 + t971 * t997) * qJD(1) + t1023 * qJDD(2) + t1009 * t935 + t1011 * t934 + mrSges(3,1) * t902 - mrSges(4,1) * t886 - mrSges(3,2) * t903 + mrSges(4,3) * t880 - pkin(2) * t810 - pkin(3) * t811 + qJ(3) * (mrSges(4,2) * t935 + t982) + t979;
t834 = t855 * mrSges(7,2) + t901 * t873 - t986;
t1018 = -t1001 * t900 - t1002 * t901 - t1022 * t885 - t1008 * t854 - t1010 * t855 + mrSges(6,1) * t840 - mrSges(7,1) * t838 - mrSges(6,2) * t841 + mrSges(7,3) * t837 - pkin(5) * t834 + qJ(6) * (-t854 * mrSges(7,2) - t900 * t873 + t993);
t1015 = mrSges(3,3) + mrSges(4,2);
t933 = (-mrSges(3,1) * t971 + mrSges(3,2) * t968) * qJD(1);
t939 = qJD(2) * mrSges(3,1) - mrSges(3,3) * t996;
t807 = m(3) * t903 - qJDD(2) * mrSges(3,2) - qJD(2) * t939 + t1015 * t935 + t933 * t995 + t982;
t941 = -qJD(2) * mrSges(3,2) + mrSges(3,3) * t995;
t808 = m(3) * t902 + qJDD(2) * mrSges(3,1) + qJD(2) * t941 - t1015 * t934 + (-t932 - t933) * t996 + t980;
t989 = t971 * t807 - t808 * t968;
t799 = m(2) * t945 - mrSges(2,1) * t974 - qJDD(1) * mrSges(2,2) + t989;
t871 = (pkin(2) * qJD(2) - (2 * qJD(3))) * t996 + t983;
t823 = t1016 * t830 + t966 * t828;
t984 = -m(5) * t851 + t887 * mrSges(5,1) - t888 * mrSges(5,2) - t919 * t904 - t920 * t905 - t823;
t815 = m(4) * t871 - mrSges(4,1) * t935 - t934 * mrSges(4,3) - t940 * t996 - t942 * t995 + t984;
t976 = -m(3) * t921 + t935 * mrSges(3,1) - mrSges(3,2) * t934 - t939 * t996 + t941 * t995 - t815;
t813 = m(2) * t944 + qJDD(1) * mrSges(2,1) - mrSges(2,2) * t974 + t976;
t1004 = t969 * t799 + t972 * t813;
t801 = t968 * t807 + t971 * t808;
t999 = t1023 * qJD(2) + (t1009 * t971 + t1011 * t968) * qJD(1);
t990 = t972 * t799 - t813 * t969;
t893 = Ifges(5,5) * t920 - Ifges(5,6) * t919 + Ifges(5,3) * t958;
t802 = mrSges(5,2) * t851 - mrSges(5,3) * t848 + Ifges(5,1) * t888 + Ifges(5,4) * t887 + Ifges(5,5) * t957 - pkin(9) * t823 + t1016 * t821 - t966 * t819 - t919 * t893 - t958 * t894;
t803 = -mrSges(5,1) * t851 + mrSges(5,3) * t849 + Ifges(5,4) * t888 + Ifges(5,2) * t887 + Ifges(5,6) * t957 - pkin(4) * t823 - t920 * t893 + t958 * t895 - t1018;
t794 = -mrSges(3,1) * t921 - mrSges(4,1) * t871 + mrSges(4,2) * t880 + mrSges(3,3) * t903 - pkin(2) * t815 - pkin(3) * t984 - pkin(8) * t988 + t997 * qJD(2) + t1009 * qJDD(2) + t1013 * t934 + t1025 * t935 - t967 * t802 - t970 * t803 - t999 * t996;
t796 = mrSges(3,2) * t921 + mrSges(4,2) * t886 - mrSges(3,3) * t902 - mrSges(4,3) * t871 - pkin(8) * t811 - qJ(3) * t815 + t998 * qJD(2) + t1011 * qJDD(2) + t1013 * t935 + t1027 * t934 + t970 * t802 - t967 * t803 + t999 * t995;
t981 = mrSges(2,1) * t944 - mrSges(2,2) * t945 + Ifges(2,3) * qJDD(1) + pkin(1) * t976 + pkin(7) * t989 + t971 * t794 + t968 * t796;
t792 = mrSges(2,1) * g(3) + mrSges(2,3) * t945 + t974 * Ifges(2,5) + Ifges(2,6) * qJDD(1) - pkin(1) * t801 - t1021;
t791 = -mrSges(2,2) * g(3) - mrSges(2,3) * t944 + Ifges(2,5) * qJDD(1) - Ifges(2,6) * t974 - pkin(7) * t801 - t794 * t968 + t796 * t971;
t1 = [-m(1) * g(1) + t990; -m(1) * g(2) + t1004; (-m(1) - m(2)) * g(3) + t801; -mrSges(1,2) * g(3) + mrSges(1,3) * g(2) - pkin(6) * t1004 + t972 * t791 - t969 * t792; mrSges(1,1) * g(3) - mrSges(1,3) * g(1) + pkin(6) * t990 + t969 * t791 + t972 * t792; -mrSges(1,1) * g(2) + mrSges(1,2) * g(1) + t981; t981; t1021; t810; -t979; t1018; t834;];
tauJB  = t1;
