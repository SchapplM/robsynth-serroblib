% Calculate vector of inverse dynamics joint torques and base forces with Newton-Euler
% S6RPRPPR4
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
%   pkin=[a2,a3,a4,a5,a6,d1,d3,d6,theta2,theta4]';
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
% Datum: 2019-05-05 16:51
% Revision: 89c353f7eff3bd693eda4e29f35b2761dbc3ada0 (2019-05-03)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauJB = S6RPRPPR4_invdynJB_fixb_snew_vp2(qJ, qJD, qJDD, g, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(6,1),zeros(3,1),zeros(10,1),zeros(7,1),zeros(7,3),zeros(7,6)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPRPPR4_invdynJB_fixb_snew_vp2: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RPRPPR4_invdynJB_fixb_snew_vp2: qJD has to be [6x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [6 1]), ...
  'S6RPRPPR4_invdynJB_fixb_snew_vp2: qJDD has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6RPRPPR4_invdynJB_fixb_snew_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S6RPRPPR4_invdynJB_fixb_snew_vp2: pkin has to be [10x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6RPRPPR4_invdynJB_fixb_snew_vp2: m has to be [7x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [7,3]), ...
  'S6RPRPPR4_invdynJB_fixb_snew_vp2: mrSges has to be [7x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [7 6]), ...
  'S6RPRPPR4_invdynJB_fixb_snew_vp2: Ifges has to be [7x6] (double)'); 

%% Symbolic Calculation
% From invdyn_fixb_NewtonEuler_linkframe_tauJB_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-05-05 16:48:29
% EndTime: 2019-05-05 16:48:39
% DurationCPUTime: 9.78s
% Computational Cost: add. (131506->344), mult. (325623->417), div. (0->0), fcn. (238882->10), ass. (0->150)
t1025 = -2 * qJD(4);
t1024 = Ifges(5,1) + Ifges(6,1);
t1016 = Ifges(5,4) - Ifges(6,5);
t1015 = -Ifges(5,5) - Ifges(6,4);
t1023 = Ifges(5,2) + Ifges(6,3);
t1022 = -Ifges(6,2) - Ifges(5,3);
t1014 = Ifges(5,6) - Ifges(6,6);
t973 = qJD(1) ^ 2;
t1012 = cos(pkin(10));
t1019 = cos(qJ(3));
t966 = cos(pkin(9));
t1018 = pkin(2) * t966;
t969 = sin(qJ(1));
t971 = cos(qJ(1));
t950 = -g(1) * t971 - g(2) * t969;
t945 = -pkin(1) * t973 + qJDD(1) * qJ(2) + t950;
t965 = sin(pkin(9));
t1000 = qJD(1) * qJD(2);
t996 = -t966 * g(3) - 0.2e1 * t965 * t1000;
t904 = (-pkin(7) * qJDD(1) + t1018 * t973 - t945) * t965 + t996;
t962 = t966 ^ 2;
t1010 = t962 * t973;
t930 = -g(3) * t965 + (0.2e1 * t1000 + t945) * t966;
t998 = qJDD(1) * t966;
t913 = -pkin(2) * t1010 + pkin(7) * t998 + t930;
t968 = sin(qJ(3));
t879 = t1019 * t913 + t968 * t904;
t1004 = qJD(1) * t965;
t997 = t966 * t1019;
t943 = -qJD(1) * t997 + t1004 * t968;
t983 = t1019 * t965 + t966 * t968;
t944 = t983 * qJD(1);
t920 = pkin(3) * t943 - qJ(4) * t944;
t972 = qJD(3) ^ 2;
t870 = -pkin(3) * t972 + qJDD(3) * qJ(4) - t920 * t943 + t879;
t1001 = t943 * qJD(3);
t1003 = qJD(3) * t944;
t961 = t965 ^ 2;
t949 = t969 * g(1) - t971 * g(2);
t990 = qJDD(2) - t949;
t926 = (-pkin(1) - t1018) * qJDD(1) + (-qJ(2) + (-t961 - t962) * pkin(7)) * t973 + t990;
t999 = qJDD(1) * t965;
t927 = -qJDD(1) * t997 + t968 * t999 + t1003;
t928 = qJDD(1) * t983 - t1001;
t872 = (-t928 + t1001) * qJ(4) + (t927 + t1003) * pkin(3) + t926;
t964 = sin(pkin(10));
t935 = t964 * qJD(3) + t1012 * t944;
t859 = t1012 * t872 + t935 * t1025 - t964 * t870;
t934 = -qJD(3) * t1012 + t964 * t944;
t1011 = t934 * t943;
t912 = t964 * qJDD(3) + t1012 * t928;
t878 = t1019 * t904 - t968 * t913;
t978 = qJDD(3) * pkin(3) + t972 * qJ(4) - t944 * t920 - qJDD(4) + t878;
t1021 = (-t912 + t1011) * qJ(5) - t978;
t1020 = 2 * qJD(5);
t1017 = -mrSges(5,3) - mrSges(6,2);
t1013 = mrSges(3,2) * t965;
t897 = mrSges(6,1) * t934 - mrSges(6,3) * t935;
t1005 = -mrSges(5,1) * t934 - mrSges(5,2) * t935 - t897;
t860 = t1012 * t870 + t934 * t1025 + t964 * t872;
t908 = mrSges(5,1) * t943 - mrSges(5,3) * t935;
t911 = -qJDD(3) * t1012 + t964 * t928;
t896 = pkin(4) * t934 - qJ(5) * t935;
t942 = t943 ^ 2;
t856 = -pkin(4) * t942 + t927 * qJ(5) + t943 * t1020 - t934 * t896 + t860;
t909 = -mrSges(6,1) * t943 + mrSges(6,2) * t935;
t857 = -t927 * pkin(4) - t942 * qJ(5) + t935 * t896 + qJDD(5) - t859;
t851 = (-t912 - t1011) * pkin(8) + (t934 * t935 - t927) * pkin(5) + t857;
t910 = -pkin(5) * t943 - pkin(8) * t935;
t933 = t934 ^ 2;
t852 = -pkin(5) * t933 + pkin(8) * t911 + t910 * t943 + t856;
t967 = sin(qJ(6));
t970 = cos(qJ(6));
t849 = t851 * t970 - t852 * t967;
t894 = t934 * t970 - t935 * t967;
t869 = qJD(6) * t894 + t911 * t967 + t912 * t970;
t895 = t934 * t967 + t935 * t970;
t877 = -mrSges(7,1) * t894 + mrSges(7,2) * t895;
t940 = qJD(6) - t943;
t882 = -mrSges(7,2) * t940 + mrSges(7,3) * t894;
t925 = qJDD(6) - t927;
t845 = m(7) * t849 + mrSges(7,1) * t925 - mrSges(7,3) * t869 - t877 * t895 + t882 * t940;
t850 = t851 * t967 + t852 * t970;
t868 = -qJD(6) * t895 + t911 * t970 - t912 * t967;
t883 = mrSges(7,1) * t940 - mrSges(7,3) * t895;
t846 = m(7) * t850 - mrSges(7,2) * t925 + mrSges(7,3) * t868 + t877 * t894 - t883 * t940;
t992 = -t967 * t845 + t970 * t846;
t982 = m(6) * t856 + t927 * mrSges(6,3) + t943 * t909 + t992;
t834 = m(5) * t860 - t927 * mrSges(5,2) + t1005 * t934 + t1017 * t911 - t943 * t908 + t982;
t838 = t970 * t845 + t967 * t846;
t907 = -mrSges(6,2) * t934 + mrSges(6,3) * t943;
t979 = -m(6) * t857 + t927 * mrSges(6,1) + t943 * t907 - t838;
t989 = -mrSges(5,2) * t943 - mrSges(5,3) * t934;
t836 = m(5) * t859 + t927 * mrSges(5,1) + t1005 * t935 + t1017 * t912 + t943 * t989 + t979;
t831 = t1012 * t834 - t836 * t964;
t921 = mrSges(4,1) * t943 + mrSges(4,2) * t944;
t937 = qJD(3) * mrSges(4,1) - mrSges(4,3) * t944;
t829 = m(4) * t879 - qJDD(3) * mrSges(4,2) - mrSges(4,3) * t927 - qJD(3) * t937 - t921 * t943 + t831;
t858 = -0.2e1 * qJD(5) * t935 + (t935 * t943 + t911) * pkin(4) + t1021;
t854 = -t933 * pkin(8) + (-pkin(4) - pkin(5)) * t911 + (-pkin(4) * t943 + t1020 + t910) * t935 - t1021;
t981 = -m(7) * t854 + t868 * mrSges(7,1) - t869 * mrSges(7,2) + t894 * t882 - t895 * t883;
t847 = m(6) * t858 + t911 * mrSges(6,1) - t912 * mrSges(6,3) + t934 * t907 - t935 * t909 + t981;
t843 = -m(5) * t978 + t911 * mrSges(5,1) + t912 * mrSges(5,2) + t935 * t908 + t934 * t989 + t847;
t936 = -qJD(3) * mrSges(4,2) - mrSges(4,3) * t943;
t842 = m(4) * t878 + qJDD(3) * mrSges(4,1) - t928 * mrSges(4,3) + qJD(3) * t936 - t944 * t921 - t843;
t822 = t1019 * t842 + t968 * t829;
t929 = -t965 * t945 + t996;
t984 = mrSges(3,3) * qJDD(1) + t973 * (-mrSges(3,1) * t966 + t1013);
t820 = m(3) * t929 - t965 * t984 + t822;
t993 = t1019 * t829 - t968 * t842;
t821 = m(3) * t930 + t966 * t984 + t993;
t994 = -t820 * t965 + t966 * t821;
t812 = m(2) * t950 - mrSges(2,1) * t973 - qJDD(1) * mrSges(2,2) + t994;
t941 = -qJDD(1) * pkin(1) - t973 * qJ(2) + t990;
t830 = t1012 * t836 + t964 * t834;
t976 = m(4) * t926 + t927 * mrSges(4,1) + t928 * mrSges(4,2) + t943 * t936 + t944 * t937 + t830;
t975 = -m(3) * t941 + mrSges(3,1) * t998 - t976 + (t961 * t973 + t1010) * mrSges(3,3);
t824 = t975 + m(2) * t949 + (mrSges(2,1) - t1013) * qJDD(1) - t973 * mrSges(2,2);
t1009 = t969 * t812 + t971 * t824;
t814 = t966 * t820 + t965 * t821;
t1008 = -t1014 * t943 - t1016 * t935 + t1023 * t934;
t1007 = t1014 * t934 + t1015 * t935 + t1022 * t943;
t1006 = -t1015 * t943 - t1016 * t934 + t1024 * t935;
t995 = t971 * t812 - t824 * t969;
t988 = Ifges(3,1) * t965 + Ifges(3,4) * t966;
t987 = Ifges(3,4) * t965 + Ifges(3,2) * t966;
t986 = Ifges(3,5) * t965 + Ifges(3,6) * t966;
t873 = Ifges(7,5) * t895 + Ifges(7,6) * t894 + Ifges(7,3) * t940;
t875 = Ifges(7,1) * t895 + Ifges(7,4) * t894 + Ifges(7,5) * t940;
t839 = -mrSges(7,1) * t854 + mrSges(7,3) * t850 + Ifges(7,4) * t869 + Ifges(7,2) * t868 + Ifges(7,6) * t925 - t873 * t895 + t875 * t940;
t874 = Ifges(7,4) * t895 + Ifges(7,2) * t894 + Ifges(7,6) * t940;
t840 = mrSges(7,2) * t854 - mrSges(7,3) * t849 + Ifges(7,1) * t869 + Ifges(7,4) * t868 + Ifges(7,5) * t925 + t873 * t894 - t874 * t940;
t815 = mrSges(5,1) * t978 - mrSges(6,1) * t858 + mrSges(6,2) * t856 + mrSges(5,3) * t860 - pkin(4) * t847 - pkin(5) * t981 - pkin(8) * t992 + t1006 * t943 + t1007 * t935 + t1014 * t927 + t1016 * t912 - t1023 * t911 - t970 * t839 - t967 * t840;
t816 = -mrSges(5,2) * t978 + mrSges(6,2) * t857 - mrSges(5,3) * t859 - mrSges(6,3) * t858 - pkin(8) * t838 - qJ(5) * t847 + t1007 * t934 + t1008 * t943 - t1015 * t927 - t1016 * t911 + t1024 * t912 - t967 * t839 + t970 * t840;
t914 = Ifges(4,5) * t944 - Ifges(4,6) * t943 + Ifges(4,3) * qJD(3);
t915 = Ifges(4,4) * t944 - Ifges(4,2) * t943 + Ifges(4,6) * qJD(3);
t808 = mrSges(4,2) * t926 - mrSges(4,3) * t878 + Ifges(4,1) * t928 - Ifges(4,4) * t927 + Ifges(4,5) * qJDD(3) - qJ(4) * t830 - qJD(3) * t915 + t1012 * t816 - t964 * t815 - t943 * t914;
t837 = t912 * mrSges(6,2) + t935 * t897 - t979;
t916 = Ifges(4,1) * t944 - Ifges(4,4) * t943 + Ifges(4,5) * qJD(3);
t977 = mrSges(7,1) * t849 - mrSges(7,2) * t850 + Ifges(7,5) * t869 + Ifges(7,6) * t868 + Ifges(7,3) * t925 + t895 * t874 - t894 * t875;
t809 = t977 - t944 * t914 + Ifges(4,4) * t928 + qJD(3) * t916 - mrSges(4,1) * t926 + mrSges(4,3) * t879 + mrSges(5,2) * t860 - mrSges(6,3) * t856 + mrSges(6,1) * t857 - mrSges(5,1) * t859 + pkin(5) * t838 + pkin(4) * t837 + Ifges(4,6) * qJDD(3) - pkin(3) * t830 - qJ(5) * t982 + t1008 * t935 + (qJ(5) * t897 - t1006) * t934 + (-Ifges(4,2) + t1022) * t927 + t1015 * t912 + (mrSges(6,2) * qJ(5) + t1014) * t911;
t947 = t986 * qJD(1);
t804 = -mrSges(3,1) * t941 + mrSges(3,3) * t930 - pkin(2) * t976 + pkin(7) * t993 + qJDD(1) * t987 - t1004 * t947 + t1019 * t809 + t968 * t808;
t806 = t966 * qJD(1) * t947 + mrSges(3,2) * t941 - mrSges(3,3) * t929 - pkin(7) * t822 + qJDD(1) * t988 + t1019 * t808 - t968 * t809;
t826 = mrSges(3,2) * t999 - t975;
t980 = mrSges(2,1) * t949 - mrSges(2,2) * t950 + Ifges(2,3) * qJDD(1) - pkin(1) * t826 + qJ(2) * t994 + t966 * t804 + t965 * t806;
t974 = mrSges(4,1) * t878 - mrSges(4,2) * t879 + Ifges(4,5) * t928 - Ifges(4,6) * t927 + Ifges(4,3) * qJDD(3) - pkin(3) * t843 + qJ(4) * t831 + t1012 * t815 + t964 * t816 + t944 * t915 + t943 * t916;
t807 = -t974 + mrSges(2,3) * t950 - mrSges(3,1) * t929 + mrSges(3,2) * t930 + mrSges(2,1) * g(3) - pkin(2) * t822 - pkin(1) * t814 + (Ifges(2,6) - t986) * qJDD(1) + (-t965 * t987 + t966 * t988 + Ifges(2,5)) * t973;
t802 = -mrSges(2,2) * g(3) - mrSges(2,3) * t949 + Ifges(2,5) * qJDD(1) - Ifges(2,6) * t973 - qJ(2) * t814 - t804 * t965 + t806 * t966;
t1 = [-m(1) * g(1) + t995; -m(1) * g(2) + t1009; (-m(1) - m(2)) * g(3) + t814; -mrSges(1,2) * g(3) + mrSges(1,3) * g(2) - pkin(6) * t1009 + t971 * t802 - t969 * t807; mrSges(1,1) * g(3) - mrSges(1,3) * g(1) + pkin(6) * t995 + t969 * t802 + t971 * t807; -mrSges(1,1) * g(2) + mrSges(1,2) * g(1) + t980; t980; t826; t974; t843; t837; t977;];
tauJB  = t1;
