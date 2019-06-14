% Calculate vector of inverse dynamics joint torques and base forces with Newton-Euler
% S6RRPPPR1
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
%   pkin=[a2,a3,a4,a5,a6,d1,d2,d6,theta3,theta4]';
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
% Datum: 2019-05-06 08:15
% Revision: 89c353f7eff3bd693eda4e29f35b2761dbc3ada0 (2019-05-03)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauJB = S6RRPPPR1_invdynJB_fixb_snew_vp2(qJ, qJD, qJDD, g, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(6,1),zeros(3,1),zeros(10,1),zeros(7,1),zeros(7,3),zeros(7,6)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRPPPR1_invdynJB_fixb_snew_vp2: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RRPPPR1_invdynJB_fixb_snew_vp2: qJD has to be [6x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [6 1]), ...
  'S6RRPPPR1_invdynJB_fixb_snew_vp2: qJDD has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6RRPPPR1_invdynJB_fixb_snew_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S6RRPPPR1_invdynJB_fixb_snew_vp2: pkin has to be [10x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6RRPPPR1_invdynJB_fixb_snew_vp2: m has to be [7x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [7,3]), ...
  'S6RRPPPR1_invdynJB_fixb_snew_vp2: mrSges has to be [7x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [7 6]), ...
  'S6RRPPPR1_invdynJB_fixb_snew_vp2: Ifges has to be [7x6] (double)'); 

%% Symbolic Calculation
% From invdyn_fixb_NewtonEuler_linkframe_tauJB_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-05-06 08:11:14
% EndTime: 2019-05-06 08:11:25
% DurationCPUTime: 10.46s
% Computational Cost: add. (144460->368), mult. (343289->448), div. (0->0), fcn. (237599->10), ass. (0->149)
t1020 = -2 * qJD(4);
t1019 = Ifges(5,1) + Ifges(6,1);
t1011 = Ifges(5,4) - Ifges(6,5);
t1010 = -Ifges(5,5) - Ifges(6,4);
t1018 = Ifges(5,2) + Ifges(6,3);
t1017 = -Ifges(6,2) - Ifges(5,3);
t1009 = Ifges(5,6) - Ifges(6,6);
t1006 = cos(pkin(10));
t1007 = cos(pkin(9));
t967 = sin(pkin(9));
t969 = sin(qJ(2));
t972 = cos(qJ(2));
t941 = (t1007 * t969 + t967 * t972) * qJD(1);
t966 = sin(pkin(10));
t927 = -qJD(2) * t1006 + t941 * t966;
t996 = qJD(1) * t972;
t997 = qJD(1) * t969;
t940 = -t1007 * t996 + t967 * t997;
t1005 = t927 * t940;
t993 = qJD(1) * qJD(2);
t952 = qJDD(1) * t969 + t972 * t993;
t953 = qJDD(1) * t972 - t969 * t993;
t923 = t1007 * t952 + t967 * t953;
t909 = t966 * qJDD(2) + t1006 * t923;
t1016 = (-t909 + t1005) * qJ(5);
t970 = sin(qJ(1));
t973 = cos(qJ(1));
t959 = -g(1) * t973 - g(2) * t970;
t975 = qJD(1) ^ 2;
t947 = -pkin(1) * t975 + qJDD(1) * pkin(7) + t959;
t1004 = t969 * t947;
t1013 = pkin(2) * t975;
t895 = qJDD(2) * pkin(2) - t952 * qJ(3) - t1004 + (qJ(3) * t993 + t1013 * t969 - g(3)) * t972;
t930 = -g(3) * t969 + t972 * t947;
t955 = qJD(2) * pkin(2) - qJ(3) * t997;
t965 = t972 ^ 2;
t899 = qJ(3) * t953 - qJD(2) * t955 - t1013 * t965 + t930;
t873 = -0.2e1 * qJD(3) * t940 + t1007 * t899 + t967 * t895;
t914 = pkin(3) * t940 - qJ(4) * t941;
t974 = qJD(2) ^ 2;
t860 = -pkin(3) * t974 + qJDD(2) * qJ(4) - t914 * t940 + t873;
t958 = t970 * g(1) - t973 * g(2);
t984 = -qJDD(1) * pkin(1) - t958;
t901 = -t953 * pkin(2) + qJDD(3) + t955 * t997 + (-qJ(3) * t965 - pkin(7)) * t975 + t984;
t922 = -t1007 * t953 + t952 * t967;
t862 = (qJD(2) * t940 - t923) * qJ(4) + (qJD(2) * t941 + t922) * pkin(3) + t901;
t928 = t966 * qJD(2) + t1006 * t941;
t854 = t1006 * t862 + t1020 * t928 - t966 * t860;
t1000 = -t1010 * t940 - t1011 * t927 + t1019 * t928;
t1001 = t1009 * t927 + t1010 * t928 + t1017 * t940;
t896 = pkin(4) * t927 - qJ(5) * t928;
t939 = t940 ^ 2;
t852 = -t922 * pkin(4) - t939 * qJ(5) + t928 * t896 + qJDD(5) - t854;
t846 = (-t909 - t1005) * pkin(8) + (t927 * t928 - t922) * pkin(5) + t852;
t1014 = 2 * qJD(5);
t855 = t1006 * t860 + t1020 * t927 + t966 * t862;
t851 = -pkin(4) * t939 + t922 * qJ(5) + t1014 * t940 - t927 * t896 + t855;
t905 = -pkin(5) * t940 - pkin(8) * t928;
t908 = -qJDD(2) * t1006 + t923 * t966;
t926 = t927 ^ 2;
t847 = -pkin(5) * t926 + pkin(8) * t908 + t905 * t940 + t851;
t968 = sin(qJ(6));
t971 = cos(qJ(6));
t845 = t846 * t968 + t847 * t971;
t995 = qJD(3) * t941;
t935 = -0.2e1 * t995;
t999 = t1007 * t895 - t967 * t899;
t979 = qJDD(2) * pkin(3) + t974 * qJ(4) - t941 * t914 - qJDD(4) + t999;
t850 = -t926 * pkin(8) + t935 + (-pkin(4) - pkin(5)) * t908 - t1016 + (-pkin(4) * t940 + t1014 + t905) * t928 + t979;
t892 = t927 * t968 + t928 * t971;
t866 = -qJD(6) * t892 + t908 * t971 - t909 * t968;
t891 = t927 * t971 - t928 * t968;
t867 = qJD(6) * t891 + t908 * t968 + t909 * t971;
t938 = qJD(6) - t940;
t868 = Ifges(7,5) * t892 + Ifges(7,6) * t891 + Ifges(7,3) * t938;
t870 = Ifges(7,1) * t892 + Ifges(7,4) * t891 + Ifges(7,5) * t938;
t921 = qJDD(6) - t922;
t834 = -mrSges(7,1) * t850 + mrSges(7,3) * t845 + Ifges(7,4) * t867 + Ifges(7,2) * t866 + Ifges(7,6) * t921 - t868 * t892 + t870 * t938;
t844 = t846 * t971 - t847 * t968;
t869 = Ifges(7,4) * t892 + Ifges(7,2) * t891 + Ifges(7,6) * t938;
t835 = mrSges(7,2) * t850 - mrSges(7,3) * t844 + Ifges(7,1) * t867 + Ifges(7,4) * t866 + Ifges(7,5) * t921 + t868 * t891 - t869 * t938;
t859 = 0.2e1 * t995 - t979;
t853 = -0.2e1 * qJD(5) * t928 + t1016 + (t928 * t940 + t908) * pkin(4) + t859;
t902 = -mrSges(6,2) * t927 + mrSges(6,3) * t940;
t904 = -mrSges(6,1) * t940 + mrSges(6,2) * t928;
t877 = -mrSges(7,2) * t938 + mrSges(7,3) * t891;
t878 = mrSges(7,1) * t938 - mrSges(7,3) * t892;
t982 = -m(7) * t850 + t866 * mrSges(7,1) - t867 * mrSges(7,2) + t891 * t877 - t892 * t878;
t842 = m(6) * t853 + t908 * mrSges(6,1) - t909 * mrSges(6,3) + t927 * t902 - t928 * t904 + t982;
t874 = -mrSges(7,1) * t891 + mrSges(7,2) * t892;
t840 = m(7) * t844 + mrSges(7,1) * t921 - mrSges(7,3) * t867 - t874 * t892 + t877 * t938;
t841 = m(7) * t845 - mrSges(7,2) * t921 + mrSges(7,3) * t866 + t874 * t891 - t878 * t938;
t989 = -t968 * t840 + t971 * t841;
t810 = -mrSges(5,1) * t859 - mrSges(6,1) * t853 + mrSges(6,2) * t851 + mrSges(5,3) * t855 - pkin(4) * t842 - pkin(5) * t982 - pkin(8) * t989 + t1000 * t940 + t1001 * t928 + t1009 * t922 + t1011 * t909 - t1018 * t908 - t971 * t834 - t968 * t835;
t1002 = -t1009 * t940 - t1011 * t928 + t1018 * t927;
t833 = t971 * t840 + t968 * t841;
t816 = mrSges(5,2) * t859 + mrSges(6,2) * t852 - mrSges(5,3) * t854 - mrSges(6,3) * t853 - pkin(8) * t833 - qJ(5) * t842 + t1001 * t927 + t1002 * t940 - t1010 * t922 - t1011 * t908 + t1019 * t909 - t968 * t834 + t971 * t835;
t1012 = -mrSges(5,3) - mrSges(6,2);
t903 = mrSges(5,1) * t940 - mrSges(5,3) * t928;
t983 = m(6) * t851 + t922 * mrSges(6,3) + t940 * t904 + t989;
t897 = mrSges(6,1) * t927 - mrSges(6,3) * t928;
t998 = -mrSges(5,1) * t927 - mrSges(5,2) * t928 - t897;
t829 = m(5) * t855 - t922 * mrSges(5,2) + t1012 * t908 - t940 * t903 + t927 * t998 + t983;
t980 = -m(6) * t852 + t922 * mrSges(6,1) + t940 * t902 - t833;
t987 = -mrSges(5,2) * t940 - mrSges(5,3) * t927;
t831 = m(5) * t854 + t922 * mrSges(5,1) + t1012 * t909 + t928 * t998 + t940 * t987 + t980;
t826 = t1006 * t829 - t831 * t966;
t915 = mrSges(4,1) * t940 + mrSges(4,2) * t941;
t932 = qJD(2) * mrSges(4,1) - mrSges(4,3) * t941;
t823 = m(4) * t873 - qJDD(2) * mrSges(4,2) - mrSges(4,3) * t922 - qJD(2) * t932 - t915 * t940 + t826;
t838 = m(5) * t859 + t908 * mrSges(5,1) + t909 * mrSges(5,2) + t928 * t903 + t927 * t987 + t842;
t872 = t935 + t999;
t931 = -qJD(2) * mrSges(4,2) - mrSges(4,3) * t940;
t837 = m(4) * t872 + qJDD(2) * mrSges(4,1) - t923 * mrSges(4,3) + qJD(2) * t931 - t941 * t915 - t838;
t817 = t1007 * t837 + t967 * t823;
t911 = Ifges(4,4) * t941 - Ifges(4,2) * t940 + Ifges(4,6) * qJD(2);
t912 = Ifges(4,1) * t941 - Ifges(4,4) * t940 + Ifges(4,5) * qJD(2);
t929 = -t972 * g(3) - t1004;
t943 = Ifges(3,6) * qJD(2) + (Ifges(3,4) * t969 + Ifges(3,2) * t972) * qJD(1);
t944 = Ifges(3,5) * qJD(2) + (Ifges(3,1) * t969 + Ifges(3,4) * t972) * qJD(1);
t1015 = (t943 * t969 - t944 * t972) * qJD(1) + (Ifges(3,3) + Ifges(4,3)) * qJDD(2) + mrSges(3,1) * t929 + mrSges(4,1) * t872 - mrSges(3,2) * t930 - mrSges(4,2) * t873 + Ifges(3,5) * t952 + Ifges(4,5) * t923 + Ifges(3,6) * t953 - Ifges(4,6) * t922 + pkin(2) * t817 - pkin(3) * t838 + qJ(4) * t826 + t1006 * t810 + t966 * t816 + t941 * t911 + t940 * t912;
t951 = (-mrSges(3,1) * t972 + mrSges(3,2) * t969) * qJD(1);
t957 = -qJD(2) * mrSges(3,2) + mrSges(3,3) * t996;
t814 = m(3) * t929 + qJDD(2) * mrSges(3,1) - mrSges(3,3) * t952 + qJD(2) * t957 - t951 * t997 + t817;
t956 = qJD(2) * mrSges(3,1) - mrSges(3,3) * t997;
t990 = t1007 * t823 - t837 * t967;
t815 = m(3) * t930 - qJDD(2) * mrSges(3,2) + mrSges(3,3) * t953 - qJD(2) * t956 + t951 * t996 + t990;
t991 = -t814 * t969 + t972 * t815;
t807 = m(2) * t959 - mrSges(2,1) * t975 - qJDD(1) * mrSges(2,2) + t991;
t825 = t1006 * t831 + t966 * t829;
t824 = m(4) * t901 + t922 * mrSges(4,1) + mrSges(4,2) * t923 + t940 * t931 + t932 * t941 + t825;
t946 = -t975 * pkin(7) + t984;
t977 = -m(3) * t946 + t953 * mrSges(3,1) - mrSges(3,2) * t952 - t956 * t997 + t957 * t996 - t824;
t819 = m(2) * t958 + qJDD(1) * mrSges(2,1) - mrSges(2,2) * t975 + t977;
t1003 = t970 * t807 + t973 * t819;
t809 = t972 * t814 + t969 * t815;
t992 = t973 * t807 - t819 * t970;
t910 = Ifges(4,5) * t941 - Ifges(4,6) * t940 + Ifges(4,3) * qJD(2);
t803 = mrSges(4,2) * t901 - mrSges(4,3) * t872 + Ifges(4,1) * t923 - Ifges(4,4) * t922 + Ifges(4,5) * qJDD(2) - qJ(4) * t825 - qJD(2) * t911 + t1006 * t816 - t966 * t810 - t940 * t910;
t832 = t909 * mrSges(6,2) + t928 * t897 - t980;
t978 = mrSges(7,1) * t844 - mrSges(7,2) * t845 + Ifges(7,5) * t867 + Ifges(7,6) * t866 + Ifges(7,3) * t921 + t892 * t869 - t891 * t870;
t804 = Ifges(4,6) * qJDD(2) + pkin(4) * t832 + pkin(5) * t833 - qJ(5) * t983 + t1002 * t928 + (qJ(5) * t897 - t1000) * t927 + (-Ifges(4,2) + t1017) * t922 + t1010 * t909 + (mrSges(6,2) * qJ(5) + t1009) * t908 - pkin(3) * t825 + t978 - mrSges(6,3) * t851 + mrSges(6,1) * t852 - mrSges(5,1) * t854 + mrSges(5,2) * t855 + mrSges(4,3) * t873 - mrSges(4,1) * t901 + qJD(2) * t912 + Ifges(4,4) * t923 - t941 * t910;
t942 = Ifges(3,3) * qJD(2) + (Ifges(3,5) * t969 + Ifges(3,6) * t972) * qJD(1);
t799 = -mrSges(3,1) * t946 + mrSges(3,3) * t930 + Ifges(3,4) * t952 + Ifges(3,2) * t953 + Ifges(3,6) * qJDD(2) - pkin(2) * t824 + qJ(3) * t990 + qJD(2) * t944 + t1007 * t804 + t967 * t803 - t942 * t997;
t801 = mrSges(3,2) * t946 - mrSges(3,3) * t929 + Ifges(3,1) * t952 + Ifges(3,4) * t953 + Ifges(3,5) * qJDD(2) - qJ(3) * t817 - qJD(2) * t943 + t1007 * t803 - t967 * t804 + t942 * t996;
t981 = mrSges(2,1) * t958 - mrSges(2,2) * t959 + Ifges(2,3) * qJDD(1) + pkin(1) * t977 + pkin(7) * t991 + t972 * t799 + t969 * t801;
t802 = mrSges(2,1) * g(3) + mrSges(2,3) * t959 + t975 * Ifges(2,5) + Ifges(2,6) * qJDD(1) - pkin(1) * t809 - t1015;
t797 = -mrSges(2,2) * g(3) - mrSges(2,3) * t958 + Ifges(2,5) * qJDD(1) - Ifges(2,6) * t975 - pkin(7) * t809 - t799 * t969 + t801 * t972;
t1 = [-m(1) * g(1) + t992; -m(1) * g(2) + t1003; (-m(1) - m(2)) * g(3) + t809; -mrSges(1,2) * g(3) + mrSges(1,3) * g(2) - pkin(6) * t1003 + t973 * t797 - t970 * t802; mrSges(1,1) * g(3) - mrSges(1,3) * g(1) + pkin(6) * t992 + t970 * t797 + t973 * t802; -mrSges(1,1) * g(2) + mrSges(1,2) * g(1) + t981; t981; t1015; t824; t838; t832; t978;];
tauJB  = t1;
