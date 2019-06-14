% Calculate vector of inverse dynamics joint torques and base forces with Newton-Euler
% S6RRRPPR7
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
%   pkin=[a2,a3,a4,a5,a6,d1,d2,d3,d6,theta5]';
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
% Datum: 2019-05-07 05:59
% Revision: 89c353f7eff3bd693eda4e29f35b2761dbc3ada0 (2019-05-03)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauJB = S6RRRPPR7_invdynJB_fixb_snew_vp2(qJ, qJD, qJDD, g, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(6,1),zeros(3,1),zeros(10,1),zeros(7,1),zeros(7,3),zeros(7,6)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRRPPR7_invdynJB_fixb_snew_vp2: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RRRPPR7_invdynJB_fixb_snew_vp2: qJD has to be [6x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [6 1]), ...
  'S6RRRPPR7_invdynJB_fixb_snew_vp2: qJDD has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6RRRPPR7_invdynJB_fixb_snew_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S6RRRPPR7_invdynJB_fixb_snew_vp2: pkin has to be [10x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6RRRPPR7_invdynJB_fixb_snew_vp2: m has to be [7x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [7,3]), ...
  'S6RRRPPR7_invdynJB_fixb_snew_vp2: mrSges has to be [7x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [7 6]), ...
  'S6RRRPPR7_invdynJB_fixb_snew_vp2: Ifges has to be [7x6] (double)'); 

%% Symbolic Calculation
% From invdyn_fixb_NewtonEuler_linkframe_tauJB_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-05-07 05:52:57
% EndTime: 2019-05-07 05:53:14
% DurationCPUTime: 12.63s
% Computational Cost: add. (195117->364), mult. (397542->441), div. (0->0), fcn. (270417->10), ass. (0->147)
t1027 = Ifges(4,1) + Ifges(5,1);
t1018 = Ifges(4,4) - Ifges(5,5);
t1017 = Ifges(4,5) + Ifges(5,4);
t1026 = Ifges(4,2) + Ifges(5,3);
t1025 = Ifges(5,2) + Ifges(4,3);
t1016 = -Ifges(5,6) + Ifges(4,6);
t1021 = cos(qJ(3));
t984 = sin(qJ(2));
t1009 = qJD(1) * t984;
t983 = sin(qJ(3));
t956 = -t1021 * qJD(2) + t983 * t1009;
t987 = cos(qJ(2));
t1008 = t987 * qJD(1);
t970 = -qJD(3) + t1008;
t1015 = t956 * t970;
t1007 = qJD(1) * qJD(2);
t1005 = t987 * t1007;
t985 = sin(qJ(1));
t988 = cos(qJ(1));
t966 = t985 * g(1) - t988 * g(2);
t990 = qJD(1) ^ 2;
t946 = -qJDD(1) * pkin(1) - t990 * pkin(7) - t966;
t960 = t984 * qJDD(1) + t1005;
t972 = t984 * t1007;
t961 = t987 * qJDD(1) - t972;
t895 = (-t960 - t1005) * pkin(8) + (-t961 + t972) * pkin(2) + t946;
t967 = -t988 * g(1) - t985 * g(2);
t947 = -t990 * pkin(1) + qJDD(1) * pkin(7) + t967;
t936 = -t984 * g(3) + t987 * t947;
t959 = (-pkin(2) * t987 - pkin(8) * t984) * qJD(1);
t989 = qJD(2) ^ 2;
t899 = -t989 * pkin(2) + qJDD(2) * pkin(8) + t959 * t1008 + t936;
t878 = t1021 * t895 - t983 * t899;
t957 = t983 * qJD(2) + t1021 * t1009;
t927 = t956 * pkin(3) - t957 * qJ(4);
t955 = -qJDD(3) + t961;
t969 = t970 ^ 2;
t871 = t955 * pkin(3) - t969 * qJ(4) + t957 * t927 + qJDD(4) - t878;
t921 = -t956 * qJD(3) + t983 * qJDD(2) + t1021 * t960;
t856 = (-t921 + t1015) * qJ(5) + (t956 * t957 + t955) * pkin(4) + t871;
t1022 = -2 * qJD(4);
t879 = t1021 * t899 + t983 * t895;
t869 = -t969 * pkin(3) - t955 * qJ(4) + t970 * t1022 - t956 * t927 + t879;
t920 = t957 * qJD(3) - t1021 * qJDD(2) + t983 * t960;
t931 = t970 * pkin(4) - t957 * qJ(5);
t954 = t956 ^ 2;
t860 = -t954 * pkin(4) + t920 * qJ(5) - t970 * t931 + t869;
t979 = sin(pkin(10));
t980 = cos(pkin(10));
t924 = t979 * t956 + t980 * t957;
t850 = -0.2e1 * qJD(5) * t924 + t980 * t856 - t979 * t860;
t888 = t979 * t920 + t980 * t921;
t923 = t980 * t956 - t979 * t957;
t848 = (t923 * t970 - t888) * pkin(9) + (t923 * t924 + t955) * pkin(5) + t850;
t851 = 0.2e1 * qJD(5) * t923 + t979 * t856 + t980 * t860;
t887 = t980 * t920 - t979 * t921;
t902 = t970 * pkin(5) - t924 * pkin(9);
t922 = t923 ^ 2;
t849 = -t922 * pkin(5) + t887 * pkin(9) - t970 * t902 + t851;
t982 = sin(qJ(6));
t986 = cos(qJ(6));
t846 = t986 * t848 - t982 * t849;
t891 = t986 * t923 - t982 * t924;
t866 = t891 * qJD(6) + t982 * t887 + t986 * t888;
t892 = t982 * t923 + t986 * t924;
t877 = -t891 * mrSges(7,1) + t892 * mrSges(7,2);
t968 = qJD(6) + t970;
t880 = -t968 * mrSges(7,2) + t891 * mrSges(7,3);
t951 = qJDD(6) + t955;
t842 = m(7) * t846 + t951 * mrSges(7,1) - t866 * mrSges(7,3) - t892 * t877 + t968 * t880;
t847 = t982 * t848 + t986 * t849;
t865 = -t892 * qJD(6) + t986 * t887 - t982 * t888;
t881 = t968 * mrSges(7,1) - t892 * mrSges(7,3);
t843 = m(7) * t847 - t951 * mrSges(7,2) + t865 * mrSges(7,3) + t891 * t877 - t968 * t881;
t832 = t986 * t842 + t982 * t843;
t893 = -t923 * mrSges(6,1) + t924 * mrSges(6,2);
t900 = -t970 * mrSges(6,2) + t923 * mrSges(6,3);
t830 = m(6) * t850 + t955 * mrSges(6,1) - t888 * mrSges(6,3) - t924 * t893 + t970 * t900 + t832;
t1001 = -t982 * t842 + t986 * t843;
t901 = t970 * mrSges(6,1) - t924 * mrSges(6,3);
t831 = m(6) * t851 - t955 * mrSges(6,2) + t887 * mrSges(6,3) + t923 * t893 - t970 * t901 + t1001;
t1002 = -t979 * t830 + t980 * t831;
t1011 = -t1017 * t970 - t1018 * t956 + t1027 * t957;
t1012 = t1016 * t956 - t1017 * t957 + t1025 * t970;
t1020 = pkin(3) * t970;
t935 = -t987 * g(3) - t984 * t947;
t898 = -qJDD(2) * pkin(2) - t989 * pkin(8) + t959 * t1009 - t935;
t997 = t920 * pkin(3) + t898 + (-t1015 - t921) * qJ(4);
t859 = -t920 * pkin(4) - t954 * qJ(5) + qJDD(5) - t997 + ((2 * qJD(4)) + t1020 + t931) * t957;
t853 = -t887 * pkin(5) - t922 * pkin(9) + t924 * t902 + t859;
t1000 = m(7) * t853 - t865 * mrSges(7,1) + t866 * mrSges(7,2) - t891 * t880 + t892 * t881;
t872 = Ifges(7,5) * t892 + Ifges(7,6) * t891 + Ifges(7,3) * t968;
t874 = Ifges(7,1) * t892 + Ifges(7,4) * t891 + Ifges(7,5) * t968;
t833 = -mrSges(7,1) * t853 + mrSges(7,3) * t847 + Ifges(7,4) * t866 + Ifges(7,2) * t865 + Ifges(7,6) * t951 - t892 * t872 + t968 * t874;
t873 = Ifges(7,4) * t892 + Ifges(7,2) * t891 + Ifges(7,6) * t968;
t834 = mrSges(7,2) * t853 - mrSges(7,3) * t846 + Ifges(7,1) * t866 + Ifges(7,4) * t865 + Ifges(7,5) * t951 + t891 * t872 - t968 * t873;
t884 = Ifges(6,5) * t924 + Ifges(6,6) * t923 + Ifges(6,3) * t970;
t886 = Ifges(6,1) * t924 + Ifges(6,4) * t923 + Ifges(6,5) * t970;
t821 = -mrSges(6,1) * t859 + mrSges(6,3) * t851 + Ifges(6,4) * t888 + Ifges(6,2) * t887 + Ifges(6,6) * t955 - pkin(5) * t1000 + pkin(9) * t1001 + t986 * t833 + t982 * t834 - t924 * t884 + t970 * t886;
t885 = Ifges(6,4) * t924 + Ifges(6,2) * t923 + Ifges(6,6) * t970;
t822 = mrSges(6,2) * t859 - mrSges(6,3) * t850 + Ifges(6,1) * t888 + Ifges(6,4) * t887 + Ifges(6,5) * t955 - pkin(9) * t832 - t982 * t833 + t986 * t834 + t923 * t884 - t970 * t885;
t844 = m(6) * t859 - t887 * mrSges(6,1) + t888 * mrSges(6,2) - t923 * t900 + t924 * t901 + t1000;
t870 = (t1022 - t1020) * t957 + t997;
t933 = t970 * mrSges(5,1) + t957 * mrSges(5,2);
t934 = -t956 * mrSges(5,2) - t970 * mrSges(5,3);
t838 = m(5) * t870 + t920 * mrSges(5,1) - t921 * mrSges(5,3) - t957 * t933 + t956 * t934 - t844;
t806 = -mrSges(4,1) * t898 - mrSges(5,1) * t870 + mrSges(5,2) * t869 + mrSges(4,3) * t879 - pkin(3) * t838 + pkin(4) * t844 - qJ(5) * t1002 - t1011 * t970 + t1012 * t957 - t1016 * t955 + t1018 * t921 - t1026 * t920 - t980 * t821 - t979 * t822;
t1013 = t1016 * t970 - t1018 * t957 + t1026 * t956;
t828 = t980 * t830 + t979 * t831;
t807 = mrSges(4,2) * t898 + mrSges(5,2) * t871 - mrSges(4,3) * t878 - mrSges(5,3) * t870 - qJ(4) * t838 - qJ(5) * t828 + t1012 * t956 - t1013 * t970 - t1017 * t955 - t1018 * t920 + t1027 * t921 - t979 * t821 + t980 * t822;
t928 = t956 * mrSges(5,1) - t957 * mrSges(5,3);
t1010 = -t956 * mrSges(4,1) - t957 * mrSges(4,2) - t928;
t1019 = -mrSges(4,3) - mrSges(5,2);
t932 = -t970 * mrSges(4,1) - t957 * mrSges(4,3);
t998 = m(5) * t869 - t955 * mrSges(5,3) - t970 * t933 + t1002;
t824 = m(4) * t879 + t955 * mrSges(4,2) + t1010 * t956 + t1019 * t920 + t970 * t932 + t998;
t930 = t970 * mrSges(4,2) - t956 * mrSges(4,3);
t995 = -m(5) * t871 - t955 * mrSges(5,1) - t970 * t934 - t828;
t825 = m(4) * t878 - t955 * mrSges(4,1) + t1010 * t957 + t1019 * t921 - t970 * t930 + t995;
t820 = t1021 * t824 - t983 * t825;
t837 = -m(4) * t898 - t920 * mrSges(4,1) - t921 * mrSges(4,2) - t956 * t930 - t957 * t932 - t838;
t944 = Ifges(3,6) * qJD(2) + (Ifges(3,4) * t984 + Ifges(3,2) * t987) * qJD(1);
t945 = Ifges(3,5) * qJD(2) + (Ifges(3,1) * t984 + Ifges(3,4) * t987) * qJD(1);
t1024 = mrSges(3,1) * t935 - mrSges(3,2) * t936 + Ifges(3,5) * t960 + Ifges(3,6) * t961 + Ifges(3,3) * qJDD(2) + pkin(2) * t837 + pkin(8) * t820 + (t944 * t984 - t945 * t987) * qJD(1) + t1021 * t806 + t983 * t807;
t827 = t921 * mrSges(5,2) + t957 * t928 - t995;
t994 = -mrSges(7,1) * t846 + mrSges(7,2) * t847 - Ifges(7,5) * t866 - Ifges(7,6) * t865 - Ifges(7,3) * t951 - t892 * t873 + t891 * t874;
t1023 = -(Ifges(6,3) + t1025) * t955 + t1011 * t956 - t1013 * t957 - t1016 * t920 + t1017 * t921 + mrSges(4,1) * t878 - mrSges(5,1) * t871 - mrSges(6,1) * t850 - mrSges(4,2) * t879 + mrSges(6,2) * t851 + mrSges(5,3) * t869 - Ifges(6,5) * t888 - Ifges(6,6) * t887 - pkin(3) * t827 - pkin(4) * t828 - pkin(5) * t832 + qJ(4) * (-t920 * mrSges(5,2) - t956 * t928 + t998) - t924 * t885 + t923 * t886 + t994;
t958 = (-mrSges(3,1) * t987 + mrSges(3,2) * t984) * qJD(1);
t963 = qJD(2) * mrSges(3,1) - mrSges(3,3) * t1009;
t818 = m(3) * t936 - qJDD(2) * mrSges(3,2) + t961 * mrSges(3,3) - qJD(2) * t963 + t958 * t1008 + t820;
t964 = -qJD(2) * mrSges(3,2) + mrSges(3,3) * t1008;
t836 = m(3) * t935 + qJDD(2) * mrSges(3,1) - t960 * mrSges(3,3) + qJD(2) * t964 - t958 * t1009 + t837;
t1003 = t987 * t818 - t984 * t836;
t810 = m(2) * t967 - t990 * mrSges(2,1) - qJDD(1) * mrSges(2,2) + t1003;
t819 = t1021 * t825 + t983 * t824;
t993 = -m(3) * t946 + t961 * mrSges(3,1) - t960 * mrSges(3,2) + t964 * t1008 - t963 * t1009 - t819;
t814 = m(2) * t966 + qJDD(1) * mrSges(2,1) - t990 * mrSges(2,2) + t993;
t1014 = t985 * t810 + t988 * t814;
t812 = t984 * t818 + t987 * t836;
t1004 = t988 * t810 - t985 * t814;
t943 = Ifges(3,3) * qJD(2) + (Ifges(3,5) * t984 + Ifges(3,6) * t987) * qJD(1);
t803 = mrSges(3,2) * t946 - mrSges(3,3) * t935 + Ifges(3,1) * t960 + Ifges(3,4) * t961 + Ifges(3,5) * qJDD(2) - pkin(8) * t819 - qJD(2) * t944 + t943 * t1008 + t1021 * t807 - t983 * t806;
t805 = -mrSges(3,1) * t946 + mrSges(3,3) * t936 + Ifges(3,4) * t960 + Ifges(3,2) * t961 + Ifges(3,6) * qJDD(2) - pkin(2) * t819 + qJD(2) * t945 - t943 * t1009 - t1023;
t996 = mrSges(2,1) * t966 - mrSges(2,2) * t967 + Ifges(2,3) * qJDD(1) + pkin(1) * t993 + pkin(7) * t1003 + t984 * t803 + t987 * t805;
t801 = mrSges(2,1) * g(3) + mrSges(2,3) * t967 + t990 * Ifges(2,5) + Ifges(2,6) * qJDD(1) - pkin(1) * t812 - t1024;
t800 = -mrSges(2,2) * g(3) - mrSges(2,3) * t966 + Ifges(2,5) * qJDD(1) - t990 * Ifges(2,6) - pkin(7) * t812 + t987 * t803 - t984 * t805;
t1 = [-m(1) * g(1) + t1004; -m(1) * g(2) + t1014; (-m(1) - m(2)) * g(3) + t812; -mrSges(1,2) * g(3) + mrSges(1,3) * g(2) - pkin(6) * t1014 + t988 * t800 - t985 * t801; mrSges(1,1) * g(3) - mrSges(1,3) * g(1) + pkin(6) * t1004 + t985 * t800 + t988 * t801; -mrSges(1,1) * g(2) + mrSges(1,2) * g(1) + t996; t996; t1024; t1023; t827; t844; -t994;];
tauJB  = t1;
