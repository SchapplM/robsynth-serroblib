% Calculate vector of inverse dynamics joint torques and base forces with Newton-Euler
% S6RRRRPP8
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
%   pkin=[a2,a3,a4,a5,a6,alpha2,d1,d2,d3,d4]';
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
% Datum: 2019-05-07 19:10
% Revision: 89c353f7eff3bd693eda4e29f35b2761dbc3ada0 (2019-05-03)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauJB = S6RRRRPP8_invdynJB_fixb_snew_vp2(qJ, qJD, qJDD, g, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(6,1),zeros(3,1),zeros(10,1),zeros(7,1),zeros(7,3),zeros(7,6)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRRRPP8_invdynJB_fixb_snew_vp2: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RRRRPP8_invdynJB_fixb_snew_vp2: qJD has to be [6x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [6 1]), ...
  'S6RRRRPP8_invdynJB_fixb_snew_vp2: qJDD has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6RRRRPP8_invdynJB_fixb_snew_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S6RRRRPP8_invdynJB_fixb_snew_vp2: pkin has to be [10x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6RRRRPP8_invdynJB_fixb_snew_vp2: m has to be [7x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [7,3]), ...
  'S6RRRRPP8_invdynJB_fixb_snew_vp2: mrSges has to be [7x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [7 6]), ...
  'S6RRRRPP8_invdynJB_fixb_snew_vp2: Ifges has to be [7x6] (double)'); 

%% Symbolic Calculation
% From invdyn_fixb_NewtonEuler_linkframe_tauJB_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-05-07 18:58:32
% EndTime: 2019-05-07 18:58:49
% DurationCPUTime: 12.94s
% Computational Cost: add. (210079->355), mult. (445784->435), div. (0->0), fcn. (344646->10), ass. (0->146)
t1030 = Ifges(5,1) + Ifges(6,1) + Ifges(7,1);
t1005 = Ifges(5,4) - Ifges(6,5) - Ifges(7,4);
t1004 = Ifges(5,5) + Ifges(6,4) - Ifges(7,5);
t1029 = -Ifges(5,2) - Ifges(6,3) - Ifges(7,2);
t1003 = -Ifges(6,6) + Ifges(7,6) + Ifges(5,6);
t1028 = Ifges(5,3) + Ifges(6,2) + Ifges(7,3);
t1006 = qJD(1) * qJD(2);
t975 = sin(pkin(6));
t979 = sin(qJ(2));
t981 = cos(qJ(2));
t961 = (-qJDD(1) * t981 + t1006 * t979) * t975;
t1021 = cos(qJ(4));
t1022 = cos(qJ(3));
t1007 = qJD(1) * t981;
t976 = cos(pkin(6));
t1013 = t976 * t979;
t1019 = pkin(8) * t975;
t980 = sin(qJ(1));
t982 = cos(qJ(1));
t968 = t980 * g(1) - g(2) * t982;
t983 = qJD(1) ^ 2;
t956 = qJDD(1) * pkin(1) + t1019 * t983 + t968;
t969 = -g(1) * t982 - g(2) * t980;
t957 = -pkin(1) * t983 + qJDD(1) * t1019 + t969;
t1009 = t956 * t1013 + t981 * t957;
t1008 = qJD(1) * t975;
t959 = (-pkin(2) * t981 - pkin(9) * t979) * t1008;
t971 = qJD(1) * t976 + qJD(2);
t970 = t971 ^ 2;
t992 = qJDD(1) * t976 + qJDD(2);
t903 = t992 * pkin(9) - t970 * pkin(2) + (-g(3) * t979 + t1007 * t959) * t975 + t1009;
t1018 = t976 * g(3);
t960 = (qJDD(1) * t979 + t1006 * t981) * t975;
t904 = t961 * pkin(2) - t960 * pkin(9) - t1018 + (-t956 + (pkin(2) * t979 - pkin(9) * t981) * t971 * qJD(1)) * t975;
t978 = sin(qJ(3));
t872 = t1022 * t903 + t978 * t904;
t997 = t979 * t1008;
t949 = t1022 * t971 - t978 * t997;
t950 = t1022 * t997 + t978 * t971;
t933 = -pkin(3) * t949 - pkin(10) * t950;
t953 = qJDD(3) + t961;
t996 = t975 * t1007;
t967 = qJD(3) - t996;
t966 = t967 ^ 2;
t868 = -pkin(3) * t966 + pkin(10) * t953 + t933 * t949 + t872;
t1012 = t976 * t981;
t1014 = t975 * t981;
t930 = -g(3) * t1014 + t1012 * t956 - t979 * t957;
t902 = -pkin(2) * t992 - t970 * pkin(9) + t959 * t997 - t930;
t928 = -qJD(3) * t950 + t1022 * t992 - t960 * t978;
t929 = t949 * qJD(3) + t1022 * t960 + t978 * t992;
t870 = (-t949 * t967 - t929) * pkin(10) + (t950 * t967 - t928) * pkin(3) + t902;
t977 = sin(qJ(4));
t863 = t1021 * t870 - t977 * t868;
t935 = -t1021 * t967 + t950 * t977;
t936 = t1021 * t950 + t977 * t967;
t907 = pkin(4) * t935 - qJ(5) * t936;
t926 = qJDD(4) - t928;
t947 = qJD(4) - t949;
t946 = t947 ^ 2;
t861 = -t926 * pkin(4) - t946 * qJ(5) + t936 * t907 + qJDD(5) - t863;
t912 = -mrSges(6,2) * t935 + mrSges(6,3) * t947;
t1027 = -m(6) * t861 + t926 * mrSges(6,1) + t947 * t912;
t1016 = t935 * t947;
t880 = -t935 * qJD(4) + t1021 * t929 + t977 * t953;
t871 = t1022 * t904 - t978 * t903;
t988 = t953 * pkin(3) + t966 * pkin(10) - t950 * t933 + t871;
t1026 = (-t880 + t1016) * qJ(5) - t988;
t913 = mrSges(7,2) * t947 + mrSges(7,3) * t935;
t1024 = -0.2e1 * t936;
t854 = qJD(6) * t1024 + (-t880 - t1016) * qJ(6) + (t935 * t936 - t926) * pkin(5) + t861;
t909 = -mrSges(7,1) * t935 + mrSges(7,2) * t936;
t991 = -m(7) * t854 + t880 * mrSges(7,3) + t936 * t909;
t852 = -t926 * mrSges(7,1) - t947 * t913 - t991;
t908 = mrSges(6,1) * t935 - mrSges(6,3) * t936;
t849 = t880 * mrSges(6,2) + t936 * t908 - t1027 + t852;
t1023 = 2 * qJD(5);
t864 = t1021 * t868 + t977 * t870;
t860 = -pkin(4) * t946 + t926 * qJ(5) + t947 * t1023 - t935 * t907 + t864;
t879 = qJD(4) * t936 - t1021 * t953 + t929 * t977;
t915 = -pkin(5) * t947 - qJ(6) * t936;
t934 = t935 ^ 2;
t856 = -pkin(5) * t934 + qJ(6) * t879 + 0.2e1 * qJD(6) * t935 + t915 * t947 + t860;
t916 = -mrSges(7,1) * t947 - mrSges(7,3) * t936;
t1001 = m(7) * t856 + t879 * mrSges(7,3) + t935 * t909;
t918 = -mrSges(6,1) * t947 + mrSges(6,2) * t936;
t990 = m(6) * t860 + t926 * mrSges(6,3) + t947 * t918 + t1001;
t998 = t1004 * t947 - t1005 * t935 + t1030 * t936;
t999 = t1003 * t947 + t1005 * t936 + t1029 * t935;
t1025 = t1028 * t926 - t1003 * t879 + t1004 * t880 + t935 * t998 + t936 * t999 + mrSges(5,1) * t863 - mrSges(6,1) * t861 - mrSges(7,1) * t854 - mrSges(5,2) * t864 + mrSges(7,2) * t856 + mrSges(6,3) * t860 - pkin(4) * t849 - pkin(5) * t852 + qJ(5) * (-t879 * mrSges(6,2) + t926 * mrSges(7,2) - t935 * t908 + t947 * t916 + t990);
t1017 = -mrSges(5,3) - mrSges(6,2);
t1015 = t975 * t979;
t931 = -g(3) * t1015 + t1009;
t954 = mrSges(3,1) * t971 - mrSges(3,3) * t997;
t958 = (-mrSges(3,1) * t981 + mrSges(3,2) * t979) * t1008;
t1010 = -mrSges(5,1) * t935 - mrSges(5,2) * t936 - t908;
t914 = -mrSges(5,2) * t947 - mrSges(5,3) * t935;
t845 = m(5) * t863 + (t913 + t914) * t947 + t1010 * t936 + (mrSges(5,1) + mrSges(7,1)) * t926 + t1017 * t880 + t991 + t1027;
t917 = mrSges(5,1) * t947 - mrSges(5,3) * t936;
t847 = m(5) * t864 + (t916 - t917) * t947 + t1010 * t935 + (-mrSges(5,2) + mrSges(7,2)) * t926 + t1017 * t879 + t990;
t842 = t1021 * t847 - t845 * t977;
t932 = -mrSges(4,1) * t949 + mrSges(4,2) * t950;
t938 = mrSges(4,1) * t967 - mrSges(4,3) * t950;
t840 = m(4) * t872 - mrSges(4,2) * t953 + mrSges(4,3) * t928 + t932 * t949 - t938 * t967 + t842;
t858 = -t934 * qJ(6) + qJDD(6) + (-pkin(4) - pkin(5)) * t879 + (-pkin(4) * t947 + t1023 + t915) * t936 - t1026;
t853 = m(7) * t858 - t879 * mrSges(7,1) + t880 * mrSges(7,2) - t935 * t913 + t936 * t916;
t862 = qJD(5) * t1024 + (t936 * t947 + t879) * pkin(4) + t1026;
t851 = m(6) * t862 + mrSges(6,1) * t879 - t880 * mrSges(6,3) + t912 * t935 - t936 * t918 - t853;
t848 = m(5) * t988 - t879 * mrSges(5,1) - mrSges(5,2) * t880 - t935 * t914 - t917 * t936 - t851;
t937 = -mrSges(4,2) * t967 + mrSges(4,3) * t949;
t844 = m(4) * t871 + mrSges(4,1) * t953 - mrSges(4,3) * t929 - t932 * t950 + t937 * t967 + t848;
t994 = t1022 * t840 - t978 * t844;
t830 = m(3) * t931 - mrSges(3,2) * t992 - t961 * mrSges(3,3) - t971 * t954 + t958 * t996 + t994;
t833 = t1022 * t844 + t978 * t840;
t942 = -t975 * t956 - t1018;
t955 = -mrSges(3,2) * t971 + mrSges(3,3) * t996;
t832 = m(3) * t942 + t961 * mrSges(3,1) + t960 * mrSges(3,2) + (t954 * t979 - t955 * t981) * t1008 + t833;
t841 = t1021 * t845 + t977 * t847;
t986 = -m(4) * t902 + t928 * mrSges(4,1) - t929 * mrSges(4,2) + t949 * t937 - t950 * t938 - t841;
t837 = m(3) * t930 + mrSges(3,1) * t992 - t960 * mrSges(3,3) + t971 * t955 - t958 * t997 + t986;
t818 = t837 * t1012 + t830 * t1013 - t832 * t975;
t815 = m(2) * t968 + qJDD(1) * mrSges(2,1) - mrSges(2,2) * t983 + t818;
t825 = t981 * t830 - t837 * t979;
t823 = m(2) * t969 - mrSges(2,1) * t983 - qJDD(1) * mrSges(2,2) + t825;
t1011 = t982 * t815 + t980 * t823;
t817 = t837 * t1014 + t830 * t1015 + t976 * t832;
t1000 = t1003 * t935 - t1004 * t936 - t1028 * t947;
t995 = -t815 * t980 + t982 * t823;
t826 = mrSges(5,1) * t988 + mrSges(5,3) * t864 - mrSges(6,1) * t862 + mrSges(6,2) * t860 + mrSges(7,1) * t858 - mrSges(7,3) * t856 + pkin(5) * t853 - qJ(6) * t1001 - pkin(4) * t851 + (-qJ(6) * t916 + t998) * t947 + t1000 * t936 + (-mrSges(7,2) * qJ(6) + t1003) * t926 + t1005 * t880 + t1029 * t879;
t834 = -mrSges(5,2) * t988 + mrSges(6,2) * t861 + mrSges(7,2) * t858 - mrSges(5,3) * t863 - mrSges(6,3) * t862 - mrSges(7,3) * t854 - qJ(5) * t851 - qJ(6) * t852 + t1000 * t935 + t1004 * t926 - t1005 * t879 + t1030 * t880 - t999 * t947;
t922 = Ifges(4,5) * t950 + Ifges(4,6) * t949 + Ifges(4,3) * t967;
t923 = Ifges(4,4) * t950 + Ifges(4,2) * t949 + Ifges(4,6) * t967;
t819 = mrSges(4,2) * t902 - mrSges(4,3) * t871 + Ifges(4,1) * t929 + Ifges(4,4) * t928 + Ifges(4,5) * t953 - pkin(10) * t841 + t1021 * t834 - t977 * t826 + t949 * t922 - t967 * t923;
t924 = Ifges(4,1) * t950 + Ifges(4,4) * t949 + Ifges(4,5) * t967;
t820 = -mrSges(4,1) * t902 + mrSges(4,3) * t872 + Ifges(4,4) * t929 + Ifges(4,2) * t928 + Ifges(4,6) * t953 - pkin(3) * t841 - t950 * t922 + t967 * t924 - t1025;
t940 = Ifges(3,6) * t971 + (Ifges(3,4) * t979 + Ifges(3,2) * t981) * t1008;
t941 = Ifges(3,5) * t971 + (Ifges(3,1) * t979 + Ifges(3,4) * t981) * t1008;
t809 = Ifges(3,5) * t960 - Ifges(3,6) * t961 + Ifges(3,3) * t992 + mrSges(3,1) * t930 - mrSges(3,2) * t931 + t978 * t819 + t1022 * t820 + pkin(2) * t986 + pkin(9) * t994 + (t940 * t979 - t941 * t981) * t1008;
t939 = Ifges(3,3) * t971 + (Ifges(3,5) * t979 + Ifges(3,6) * t981) * t1008;
t811 = mrSges(3,2) * t942 - mrSges(3,3) * t930 + Ifges(3,1) * t960 - Ifges(3,4) * t961 + Ifges(3,5) * t992 - pkin(9) * t833 + t1022 * t819 - t978 * t820 + t939 * t996 - t971 * t940;
t984 = mrSges(4,1) * t871 - mrSges(4,2) * t872 + Ifges(4,5) * t929 + Ifges(4,6) * t928 + Ifges(4,3) * t953 + pkin(3) * t848 + pkin(10) * t842 + t1021 * t826 + t977 * t834 + t950 * t923 - t949 * t924;
t813 = -mrSges(3,1) * t942 + mrSges(3,3) * t931 + Ifges(3,4) * t960 - Ifges(3,2) * t961 + Ifges(3,6) * t992 - pkin(2) * t833 - t939 * t997 + t971 * t941 - t984;
t987 = mrSges(2,1) * t968 - mrSges(2,2) * t969 + Ifges(2,3) * qJDD(1) + pkin(1) * t818 + t813 * t1014 + t811 * t1015 + t825 * t1019 + t976 * t809;
t807 = -mrSges(2,2) * g(3) - mrSges(2,3) * t968 + Ifges(2,5) * qJDD(1) - t983 * Ifges(2,6) + t981 * t811 - t979 * t813 + (-t817 * t975 - t818 * t976) * pkin(8);
t806 = mrSges(2,1) * g(3) + mrSges(2,3) * t969 + t983 * Ifges(2,5) + Ifges(2,6) * qJDD(1) - pkin(1) * t817 - t975 * t809 + (pkin(8) * t825 + t811 * t979 + t813 * t981) * t976;
t1 = [-m(1) * g(1) + t995; -m(1) * g(2) + t1011; (-m(1) - m(2)) * g(3) + t817; -mrSges(1,2) * g(3) + mrSges(1,3) * g(2) - pkin(7) * t1011 - t980 * t806 + t982 * t807; mrSges(1,1) * g(3) - mrSges(1,3) * g(1) + pkin(7) * t995 + t982 * t806 + t980 * t807; -mrSges(1,1) * g(2) + mrSges(1,2) * g(1) + t987; t987; t809; t984; t1025; t849; t853;];
tauJB  = t1;
