% Calculate vector of inverse dynamics joint torques and base forces with Newton-Euler
% S6RRRPRP10
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
%   pkin=[a2,a3,a4,a5,a6,alpha2,d1,d2,d3,d5,theta4]';
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
% Datum: 2019-05-07 09:06
% Revision: 89c353f7eff3bd693eda4e29f35b2761dbc3ada0 (2019-05-03)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauJB = S6RRRPRP10_invdynJB_fixb_snew_vp2(qJ, qJD, qJDD, g, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(6,1),zeros(3,1),zeros(11,1),zeros(7,1),zeros(7,3),zeros(7,6)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRRPRP10_invdynJB_fixb_snew_vp2: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RRRPRP10_invdynJB_fixb_snew_vp2: qJD has to be [6x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [6 1]), ...
  'S6RRRPRP10_invdynJB_fixb_snew_vp2: qJDD has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6RRRPRP10_invdynJB_fixb_snew_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [11 1]), ...
  'S6RRRPRP10_invdynJB_fixb_snew_vp2: pkin has to be [11x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6RRRPRP10_invdynJB_fixb_snew_vp2: m has to be [7x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [7,3]), ...
  'S6RRRPRP10_invdynJB_fixb_snew_vp2: mrSges has to be [7x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [7 6]), ...
  'S6RRRPRP10_invdynJB_fixb_snew_vp2: Ifges has to be [7x6] (double)'); 

%% Symbolic Calculation
% From invdyn_fixb_NewtonEuler_linkframe_tauJB_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-05-07 08:53:24
% EndTime: 2019-05-07 08:54:02
% DurationCPUTime: 33.06s
% Computational Cost: add. (540796->377), mult. (1179232->477), div. (0->0), fcn. (934956->12), ass. (0->156)
t1008 = Ifges(6,1) + Ifges(7,1);
t999 = Ifges(6,4) - Ifges(7,5);
t998 = -Ifges(6,5) - Ifges(7,4);
t1007 = Ifges(6,2) + Ifges(7,3);
t997 = Ifges(6,6) - Ifges(7,6);
t1006 = -Ifges(6,3) - Ifges(7,2);
t959 = sin(pkin(6));
t964 = sin(qJ(2));
t966 = cos(qJ(2));
t983 = qJD(1) * qJD(2);
t944 = (-qJDD(1) * t966 + t964 * t983) * t959;
t1003 = cos(qJ(5));
t1004 = cos(qJ(3));
t961 = cos(pkin(6));
t954 = qJD(1) * t961 + qJD(2);
t963 = sin(qJ(3));
t985 = qJD(1) * t959;
t981 = t964 * t985;
t933 = t1004 * t981 + t963 * t954;
t984 = qJD(1) * t966;
t980 = t959 * t984;
t949 = qJD(3) - t980;
t958 = sin(pkin(11));
t960 = cos(pkin(11));
t920 = -t933 * t958 + t949 * t960;
t921 = t933 * t960 + t949 * t958;
t962 = sin(qJ(5));
t893 = -t1003 * t920 + t921 * t962;
t932 = -t1004 * t954 + t963 * t981;
t943 = (qJDD(1) * t964 + t966 * t983) * t959;
t953 = qJDD(1) * t961 + qJDD(2);
t911 = -t932 * qJD(3) + t1004 * t943 + t963 * t953;
t936 = qJDD(3) + t944;
t896 = -t911 * t958 + t936 * t960;
t897 = t911 * t960 + t936 * t958;
t861 = -t893 * qJD(5) + t1003 * t897 + t962 * t896;
t894 = t1003 * t921 + t962 * t920;
t873 = mrSges(7,1) * t893 - mrSges(7,3) * t894;
t942 = (-pkin(2) * t966 - pkin(9) * t964) * t985;
t952 = t954 ^ 2;
t1002 = pkin(8) * t959;
t965 = sin(qJ(1));
t967 = cos(qJ(1));
t950 = t965 * g(1) - g(2) * t967;
t968 = qJD(1) ^ 2;
t939 = qJDD(1) * pkin(1) + t1002 * t968 + t950;
t951 = -g(1) * t967 - g(2) * t965;
t940 = -pkin(1) * t968 + qJDD(1) * t1002 + t951;
t993 = t961 * t964;
t986 = t939 * t993 + t966 * t940;
t890 = -t952 * pkin(2) + t953 * pkin(9) + (-g(3) * t964 + t942 * t984) * t959 + t986;
t1001 = t961 * g(3);
t891 = t944 * pkin(2) - t943 * pkin(9) - t1001 + (-t939 + (pkin(2) * t964 - pkin(9) * t966) * t954 * qJD(1)) * t959;
t863 = t1004 * t890 + t963 * t891;
t914 = pkin(3) * t932 - qJ(4) * t933;
t948 = t949 ^ 2;
t854 = -pkin(3) * t948 + qJ(4) * t936 - t914 * t932 + t863;
t992 = t961 * t966;
t994 = t959 * t966;
t912 = -g(3) * t994 + t939 * t992 - t964 * t940;
t889 = -t953 * pkin(2) - t952 * pkin(9) + t942 * t981 - t912;
t910 = qJD(3) * t933 - t1004 * t953 + t943 * t963;
t857 = (t932 * t949 - t911) * qJ(4) + (t933 * t949 + t910) * pkin(3) + t889;
t849 = -0.2e1 * qJD(4) * t921 - t958 * t854 + t960 * t857;
t846 = (t920 * t932 - t897) * pkin(10) + (t920 * t921 + t910) * pkin(4) + t849;
t850 = 0.2e1 * qJD(4) * t920 + t960 * t854 + t958 * t857;
t901 = pkin(4) * t932 - pkin(10) * t921;
t919 = t920 ^ 2;
t848 = -pkin(4) * t919 + pkin(10) * t896 - t901 * t932 + t850;
t842 = t1003 * t846 - t962 * t848;
t872 = pkin(5) * t893 - qJ(6) * t894;
t908 = qJDD(5) + t910;
t931 = qJD(5) + t932;
t930 = t931 ^ 2;
t841 = -t908 * pkin(5) - t930 * qJ(6) + t894 * t872 + qJDD(6) - t842;
t877 = -mrSges(7,2) * t893 + mrSges(7,3) * t931;
t976 = -m(7) * t841 + t908 * mrSges(7,1) + t931 * t877;
t837 = t861 * mrSges(7,2) + t894 * t873 - t976;
t843 = t1003 * t848 + t962 * t846;
t840 = -pkin(5) * t930 + qJ(6) * t908 + 0.2e1 * qJD(6) * t931 - t872 * t893 + t843;
t860 = t894 * qJD(5) - t1003 * t896 + t962 * t897;
t880 = -mrSges(7,1) * t931 + mrSges(7,2) * t894;
t982 = m(7) * t840 + t908 * mrSges(7,3) + t931 * t880;
t988 = t1008 * t894 - t999 * t893 - t998 * t931;
t990 = t1007 * t893 - t894 * t999 - t931 * t997;
t1005 = t997 * t860 + t998 * t861 - t988 * t893 + t990 * t894 + t1006 * t908 - mrSges(6,1) * t842 + mrSges(7,1) * t841 + mrSges(6,2) * t843 - mrSges(7,3) * t840 + pkin(5) * t837 - qJ(6) * (-t860 * mrSges(7,2) - t893 * t873 + t982);
t1000 = -mrSges(6,3) - mrSges(7,2);
t995 = t959 * t964;
t913 = -g(3) * t995 + t986;
t937 = mrSges(3,1) * t954 - mrSges(3,3) * t981;
t941 = (-mrSges(3,1) * t966 + mrSges(3,2) * t964) * t985;
t879 = mrSges(6,1) * t931 - mrSges(6,3) * t894;
t987 = -mrSges(6,1) * t893 - mrSges(6,2) * t894 - t873;
t830 = m(6) * t843 - t908 * mrSges(6,2) + t1000 * t860 - t931 * t879 + t893 * t987 + t982;
t878 = -mrSges(6,2) * t931 - mrSges(6,3) * t893;
t832 = m(6) * t842 + t908 * mrSges(6,1) + t1000 * t861 + t931 * t878 + t894 * t987 + t976;
t825 = t1003 * t832 + t962 * t830;
t898 = -mrSges(5,1) * t920 + mrSges(5,2) * t921;
t975 = -mrSges(5,2) * t932 + mrSges(5,3) * t920;
t823 = m(5) * t849 + t910 * mrSges(5,1) - t897 * mrSges(5,3) - t921 * t898 + t932 * t975 + t825;
t900 = mrSges(5,1) * t932 - mrSges(5,3) * t921;
t977 = t1003 * t830 - t832 * t962;
t824 = m(5) * t850 - mrSges(5,2) * t910 + mrSges(5,3) * t896 + t898 * t920 - t900 * t932 + t977;
t821 = -t823 * t958 + t960 * t824;
t915 = mrSges(4,1) * t932 + mrSges(4,2) * t933;
t923 = mrSges(4,1) * t949 - mrSges(4,3) * t933;
t819 = m(4) * t863 - mrSges(4,2) * t936 - mrSges(4,3) * t910 - t915 * t932 - t923 * t949 + t821;
t862 = t1004 * t891 - t963 * t890;
t853 = -t936 * pkin(3) - t948 * qJ(4) + t933 * t914 + qJDD(4) - t862;
t851 = -t896 * pkin(4) - t919 * pkin(10) + t921 * t901 + t853;
t844 = -0.2e1 * qJD(6) * t894 + (t893 * t931 - t861) * qJ(6) + (t894 * t931 + t860) * pkin(5) + t851;
t838 = m(7) * t844 + t860 * mrSges(7,1) - t861 * mrSges(7,3) + t893 * t877 - t894 * t880;
t970 = m(6) * t851 + t860 * mrSges(6,1) + t861 * mrSges(6,2) + t893 * t878 + t894 * t879 + t838;
t835 = m(5) * t853 - t896 * mrSges(5,1) + t897 * mrSges(5,2) + t921 * t900 - t920 * t975 + t970;
t922 = -mrSges(4,2) * t949 - mrSges(4,3) * t932;
t834 = m(4) * t862 + t936 * mrSges(4,1) - t911 * mrSges(4,3) - t933 * t915 + t949 * t922 - t835;
t978 = t1004 * t819 - t834 * t963;
t810 = m(3) * t913 - mrSges(3,2) * t953 - mrSges(3,3) * t944 - t937 * t954 + t941 * t980 + t978;
t813 = t1004 * t834 + t963 * t819;
t927 = -t959 * t939 - t1001;
t938 = -mrSges(3,2) * t954 + mrSges(3,3) * t980;
t812 = m(3) * t927 + t944 * mrSges(3,1) + t943 * mrSges(3,2) + (t937 * t964 - t938 * t966) * t985 + t813;
t820 = t823 * t960 + t824 * t958;
t972 = -m(4) * t889 - t910 * mrSges(4,1) - mrSges(4,2) * t911 - t932 * t922 - t923 * t933 - t820;
t816 = m(3) * t912 + mrSges(3,1) * t953 - mrSges(3,3) * t943 + t938 * t954 - t941 * t981 + t972;
t798 = t810 * t993 - t812 * t959 + t816 * t992;
t795 = m(2) * t950 + qJDD(1) * mrSges(2,1) - mrSges(2,2) * t968 + t798;
t803 = t966 * t810 - t816 * t964;
t801 = m(2) * t951 - mrSges(2,1) * t968 - qJDD(1) * mrSges(2,2) + t803;
t991 = t967 * t795 + t965 * t801;
t989 = t1006 * t931 + t893 * t997 + t894 * t998;
t797 = t810 * t995 + t961 * t812 + t816 * t994;
t979 = -t795 * t965 + t967 * t801;
t826 = -mrSges(6,1) * t851 - mrSges(7,1) * t844 + mrSges(7,2) * t840 + mrSges(6,3) * t843 - pkin(5) * t838 - t1007 * t860 + t999 * t861 + t989 * t894 + t997 * t908 + t988 * t931;
t827 = mrSges(6,2) * t851 + mrSges(7,2) * t841 - mrSges(6,3) * t842 - mrSges(7,3) * t844 - qJ(6) * t838 + t1008 * t861 - t999 * t860 + t989 * t893 - t998 * t908 + t990 * t931;
t881 = Ifges(5,5) * t921 + Ifges(5,6) * t920 + Ifges(5,3) * t932;
t883 = Ifges(5,1) * t921 + Ifges(5,4) * t920 + Ifges(5,5) * t932;
t805 = -mrSges(5,1) * t853 + mrSges(5,3) * t850 + Ifges(5,4) * t897 + Ifges(5,2) * t896 + Ifges(5,6) * t910 - pkin(4) * t970 + pkin(10) * t977 + t1003 * t826 + t962 * t827 - t921 * t881 + t932 * t883;
t882 = Ifges(5,4) * t921 + Ifges(5,2) * t920 + Ifges(5,6) * t932;
t806 = mrSges(5,2) * t853 - mrSges(5,3) * t849 + Ifges(5,1) * t897 + Ifges(5,4) * t896 + Ifges(5,5) * t910 - pkin(10) * t825 + t1003 * t827 - t962 * t826 + t920 * t881 - t932 * t882;
t904 = Ifges(4,5) * t933 - Ifges(4,6) * t932 + Ifges(4,3) * t949;
t905 = Ifges(4,4) * t933 - Ifges(4,2) * t932 + Ifges(4,6) * t949;
t793 = mrSges(4,2) * t889 - mrSges(4,3) * t862 + Ifges(4,1) * t911 - Ifges(4,4) * t910 + Ifges(4,5) * t936 - qJ(4) * t820 - t805 * t958 + t806 * t960 - t904 * t932 - t905 * t949;
t906 = Ifges(4,1) * t933 - Ifges(4,4) * t932 + Ifges(4,5) * t949;
t804 = t1005 - pkin(3) * t820 - mrSges(4,1) * t889 + t949 * t906 + mrSges(4,3) * t863 - mrSges(5,1) * t849 + mrSges(5,2) * t850 - Ifges(5,6) * t896 - t933 * t904 + Ifges(4,6) * t936 + t920 * t883 - t921 * t882 - pkin(4) * t825 - Ifges(5,5) * t897 + (-Ifges(5,3) - Ifges(4,2)) * t910 + Ifges(4,4) * t911;
t925 = Ifges(3,6) * t954 + (Ifges(3,4) * t964 + Ifges(3,2) * t966) * t985;
t926 = Ifges(3,5) * t954 + (Ifges(3,1) * t964 + Ifges(3,4) * t966) * t985;
t788 = Ifges(3,5) * t943 - Ifges(3,6) * t944 + Ifges(3,3) * t953 + mrSges(3,1) * t912 - mrSges(3,2) * t913 + t963 * t793 + t1004 * t804 + pkin(2) * t972 + pkin(9) * t978 + (t925 * t964 - t926 * t966) * t985;
t924 = Ifges(3,3) * t954 + (Ifges(3,5) * t964 + Ifges(3,6) * t966) * t985;
t790 = mrSges(3,2) * t927 - mrSges(3,3) * t912 + Ifges(3,1) * t943 - Ifges(3,4) * t944 + Ifges(3,5) * t953 - pkin(9) * t813 + t1004 * t793 - t963 * t804 + t924 * t980 - t954 * t925;
t969 = mrSges(4,1) * t862 - mrSges(4,2) * t863 + Ifges(4,5) * t911 - Ifges(4,6) * t910 + Ifges(4,3) * t936 - pkin(3) * t835 + qJ(4) * t821 + t960 * t805 + t958 * t806 + t933 * t905 + t932 * t906;
t792 = -mrSges(3,1) * t927 + mrSges(3,3) * t913 + Ifges(3,4) * t943 - Ifges(3,2) * t944 + Ifges(3,6) * t953 - pkin(2) * t813 - t924 * t981 + t954 * t926 - t969;
t973 = mrSges(2,1) * t950 - mrSges(2,2) * t951 + Ifges(2,3) * qJDD(1) + pkin(1) * t798 + t803 * t1002 + t961 * t788 + t790 * t995 + t792 * t994;
t786 = -mrSges(2,2) * g(3) - mrSges(2,3) * t950 + Ifges(2,5) * qJDD(1) - t968 * Ifges(2,6) + t966 * t790 - t964 * t792 + (-t797 * t959 - t798 * t961) * pkin(8);
t785 = mrSges(2,1) * g(3) + mrSges(2,3) * t951 + t968 * Ifges(2,5) + Ifges(2,6) * qJDD(1) - pkin(1) * t797 - t959 * t788 + (pkin(8) * t803 + t790 * t964 + t792 * t966) * t961;
t1 = [-m(1) * g(1) + t979; -m(1) * g(2) + t991; (-m(1) - m(2)) * g(3) + t797; -mrSges(1,2) * g(3) + mrSges(1,3) * g(2) - pkin(7) * t991 - t965 * t785 + t967 * t786; mrSges(1,1) * g(3) - mrSges(1,3) * g(1) + pkin(7) * t979 + t967 * t785 + t965 * t786; -mrSges(1,1) * g(2) + mrSges(1,2) * g(1) + t973; t973; t788; t969; t835; -t1005; t837;];
tauJB  = t1;
