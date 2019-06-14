% Calculate vector of inverse dynamics joint torques and base forces with Newton-Euler
% S6RRPPPR4
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
%   pkin=[a2,a3,a4,a5,a6,d1,d2,d6,theta4]';
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
% Datum: 2019-05-06 08:45
% Revision: 89c353f7eff3bd693eda4e29f35b2761dbc3ada0 (2019-05-03)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauJB = S6RRPPPR4_invdynJB_fixb_snew_vp2(qJ, qJD, qJDD, g, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(6,1),zeros(3,1),zeros(9,1),zeros(7,1),zeros(7,3),zeros(7,6)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRPPPR4_invdynJB_fixb_snew_vp2: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RRPPPR4_invdynJB_fixb_snew_vp2: qJD has to be [6x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [6 1]), ...
  'S6RRPPPR4_invdynJB_fixb_snew_vp2: qJDD has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6RRPPPR4_invdynJB_fixb_snew_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S6RRPPPR4_invdynJB_fixb_snew_vp2: pkin has to be [9x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6RRPPPR4_invdynJB_fixb_snew_vp2: m has to be [7x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [7,3]), ...
  'S6RRPPPR4_invdynJB_fixb_snew_vp2: mrSges has to be [7x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [7 6]), ...
  'S6RRPPPR4_invdynJB_fixb_snew_vp2: Ifges has to be [7x6] (double)'); 

%% Symbolic Calculation
% From invdyn_fixb_NewtonEuler_linkframe_tauJB_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-05-06 08:42:05
% EndTime: 2019-05-06 08:42:13
% DurationCPUTime: 5.36s
% Computational Cost: add. (53151->350), mult. (117445->406), div. (0->0), fcn. (66602->8), ass. (0->144)
t1027 = -2 * qJD(4);
t1026 = Ifges(5,1) + Ifges(6,1);
t1013 = Ifges(3,4) + Ifges(4,6);
t1012 = Ifges(5,4) - Ifges(6,5);
t1011 = Ifges(3,5) - Ifges(4,4);
t1010 = Ifges(5,5) + Ifges(6,4);
t1025 = Ifges(3,2) + Ifges(4,3);
t1024 = Ifges(4,2) + Ifges(3,1);
t1023 = Ifges(5,2) + Ifges(6,3);
t1022 = -Ifges(6,2) - Ifges(5,3);
t1009 = Ifges(3,6) - Ifges(4,5);
t1008 = Ifges(5,6) - Ifges(6,6);
t1021 = Ifges(3,3) + Ifges(4,1);
t959 = sin(pkin(9));
t960 = cos(pkin(9));
t966 = cos(qJ(2));
t996 = qJD(1) * t966;
t920 = t959 * qJD(2) + t960 * t996;
t921 = t960 * qJD(2) - t959 * t996;
t963 = sin(qJ(2));
t994 = t963 * qJD(1);
t1001 = t1010 * t994 - t1012 * t920 + t1026 * t921;
t1002 = t1008 * t920 - t1010 * t921 + t1022 * t994;
t969 = qJD(1) ^ 2;
t1006 = t963 ^ 2 * t969;
t993 = qJD(1) * qJD(2);
t990 = t963 * t993;
t932 = t966 * qJDD(1) - t990;
t938 = pkin(3) * t994 - qJD(2) * qJ(4);
t958 = t966 ^ 2;
t989 = t966 * t993;
t931 = t963 * qJDD(1) + t989;
t964 = sin(qJ(1));
t967 = cos(qJ(1));
t941 = t964 * g(1) - t967 * g(2);
t982 = -qJDD(1) * pkin(1) - t941;
t974 = pkin(2) * t990 - 0.2e1 * qJD(3) * t994 + (-t931 - t989) * qJ(3) + t982;
t854 = -t938 * t994 + (-pkin(3) * t958 - pkin(7)) * t969 + (-pkin(2) - qJ(4)) * t932 + t974;
t942 = -t967 * g(1) - t964 * g(2);
t915 = -t969 * pkin(1) + qJDD(1) * pkin(7) + t942;
t889 = -t966 * g(3) - t963 * t915;
t928 = (-pkin(2) * t966 - qJ(3) * t963) * qJD(1);
t968 = qJD(2) ^ 2;
t869 = -qJDD(2) * pkin(2) - t968 * qJ(3) + t928 * t994 + qJDD(3) - t889;
t864 = (-t963 * t966 * t969 - qJDD(2)) * qJ(4) + (t931 - t989) * pkin(3) + t869;
t848 = t921 * t1027 - t959 * t854 + t960 * t864;
t884 = t920 * pkin(4) - t921 * qJ(5);
t846 = -t931 * pkin(4) - qJ(5) * t1006 + t921 * t884 + qJDD(5) - t848;
t898 = t960 * qJDD(2) - t959 * t932;
t991 = t920 * t994;
t842 = (-t898 - t991) * pkin(8) + (t920 * t921 - t931) * pkin(5) + t846;
t1016 = 2 * qJD(5);
t849 = t920 * t1027 + t960 * t854 + t959 * t864;
t845 = -pkin(4) * t1006 + t931 * qJ(5) + t994 * t1016 - t920 * t884 + t849;
t897 = t959 * qJDD(2) + t960 * t932;
t899 = -pkin(5) * t994 - t921 * pkin(8);
t918 = t920 ^ 2;
t843 = -t918 * pkin(5) + t897 * pkin(8) + t899 * t994 + t845;
t962 = sin(qJ(6));
t965 = cos(qJ(6));
t841 = t962 * t842 + t965 * t843;
t1019 = (-t898 + t991) * qJ(5);
t992 = qJD(3) * qJD(2);
t950 = -0.2e1 * t992;
t890 = -t963 * g(3) + t966 * t915;
t983 = t968 * pkin(2) - qJDD(2) * qJ(3) - t928 * t996 - t890;
t975 = t958 * t969 * qJ(4) - t932 * pkin(3) - qJD(2) * t938 - qJDD(4) + t983;
t847 = -t918 * pkin(8) + t950 + (-pkin(4) - pkin(5)) * t897 - t1019 + (-pkin(4) * t994 + t1016 + t899) * t921 + t975;
t883 = t962 * t920 + t965 * t921;
t856 = -t883 * qJD(6) + t965 * t897 - t962 * t898;
t882 = t965 * t920 - t962 * t921;
t857 = t882 * qJD(6) + t962 * t897 + t965 * t898;
t945 = qJD(6) - t994;
t861 = Ifges(7,5) * t883 + Ifges(7,6) * t882 + Ifges(7,3) * t945;
t863 = Ifges(7,1) * t883 + Ifges(7,4) * t882 + Ifges(7,5) * t945;
t927 = qJDD(6) - t931;
t829 = -mrSges(7,1) * t847 + mrSges(7,3) * t841 + Ifges(7,4) * t857 + Ifges(7,2) * t856 + Ifges(7,6) * t927 - t883 * t861 + t945 * t863;
t840 = t965 * t842 - t962 * t843;
t862 = Ifges(7,4) * t883 + Ifges(7,2) * t882 + Ifges(7,6) * t945;
t830 = mrSges(7,2) * t847 - mrSges(7,3) * t840 + Ifges(7,1) * t857 + Ifges(7,4) * t856 + Ifges(7,5) * t927 + t882 * t861 - t945 * t862;
t860 = 0.2e1 * t992 - t975;
t851 = -0.2e1 * qJD(5) * t921 + t1019 + (t921 * t994 + t897) * pkin(4) + t860;
t893 = -t920 * mrSges(6,2) + mrSges(6,3) * t994;
t896 = -mrSges(6,1) * t994 + t921 * mrSges(6,2);
t870 = -t945 * mrSges(7,2) + t882 * mrSges(7,3);
t871 = t945 * mrSges(7,1) - t883 * mrSges(7,3);
t977 = -m(7) * t847 + t856 * mrSges(7,1) - t857 * mrSges(7,2) + t882 * t870 - t883 * t871;
t839 = m(6) * t851 + t897 * mrSges(6,1) - t898 * mrSges(6,3) + t920 * t893 - t921 * t896 + t977;
t866 = -t882 * mrSges(7,1) + t883 * mrSges(7,2);
t837 = m(7) * t840 + t927 * mrSges(7,1) - t857 * mrSges(7,3) - t883 * t866 + t945 * t870;
t838 = m(7) * t841 - t927 * mrSges(7,2) + t856 * mrSges(7,3) + t882 * t866 - t945 * t871;
t986 = -t962 * t837 + t965 * t838;
t808 = -mrSges(5,1) * t860 - mrSges(6,1) * t851 + mrSges(6,2) * t845 + mrSges(5,3) * t849 - pkin(4) * t839 - pkin(5) * t977 - pkin(8) * t986 + t1001 * t994 + t1002 * t921 + t1008 * t931 + t1012 * t898 - t1023 * t897 - t965 * t829 - t962 * t830;
t1003 = -t1008 * t994 - t1012 * t921 + t1023 * t920;
t828 = t965 * t837 + t962 * t838;
t814 = mrSges(5,2) * t860 + mrSges(6,2) * t846 - mrSges(5,3) * t848 - mrSges(6,3) * t851 - pkin(8) * t828 - qJ(5) * t839 + t1002 * t920 + t1003 * t994 + t1010 * t931 - t1012 * t897 + t1026 * t898 - t962 * t829 + t965 * t830;
t929 = (mrSges(4,2) * t966 - mrSges(4,3) * t963) * qJD(1);
t939 = -mrSges(4,1) * t996 - qJD(2) * mrSges(4,3);
t885 = t920 * mrSges(6,1) - t921 * mrSges(6,3);
t1000 = -t920 * mrSges(5,1) - t921 * mrSges(5,2) - t885;
t1014 = -mrSges(5,3) - mrSges(6,2);
t895 = mrSges(5,1) * t994 - t921 * mrSges(5,3);
t980 = m(6) * t845 + t931 * mrSges(6,3) + t896 * t994 + t986;
t824 = m(5) * t849 - t931 * mrSges(5,2) + t1000 * t920 + t1014 * t897 - t895 * t994 + t980;
t894 = -mrSges(5,2) * t994 - t920 * mrSges(5,3);
t976 = -m(6) * t846 + t931 * mrSges(6,1) + t893 * t994 - t828;
t826 = m(5) * t848 + t931 * mrSges(5,1) + t1000 * t921 + t1014 * t898 + t894 * t994 + t976;
t822 = t959 * t824 + t960 * t826;
t979 = -m(4) * t869 - t931 * mrSges(4,1) - t822;
t821 = qJDD(2) * mrSges(4,2) + qJD(2) * t939 + t929 * t994 - t979;
t868 = t950 + t983;
t835 = m(5) * t860 + t897 * mrSges(5,1) + t898 * mrSges(5,2) + t920 * t894 + t921 * t895 + t839;
t940 = mrSges(4,1) * t994 + qJD(2) * mrSges(4,2);
t970 = -m(4) * t868 + qJDD(2) * mrSges(4,3) + qJD(2) * t940 + t929 * t996 + t835;
t997 = t1011 * qJD(2) + (t1013 * t966 + t1024 * t963) * qJD(1);
t998 = t1009 * qJD(2) + (t1013 * t963 + t1025 * t966) * qJD(1);
t1020 = (t998 * t963 - t997 * t966) * qJD(1) + t1021 * qJDD(2) + t1009 * t932 + t1011 * t931 + mrSges(3,1) * t889 - mrSges(3,2) * t890 + mrSges(4,2) * t869 - mrSges(4,3) * t868 - pkin(2) * t821 + qJ(3) * (t932 * mrSges(4,1) + t970) - qJ(4) * t822 - t959 * t808 + t960 * t814;
t1015 = t969 * pkin(7);
t930 = (-mrSges(3,1) * t966 + mrSges(3,2) * t963) * qJD(1);
t937 = -qJD(2) * mrSges(3,2) + mrSges(3,3) * t996;
t819 = m(3) * t889 - t931 * mrSges(3,3) + (mrSges(3,1) - mrSges(4,2)) * qJDD(2) + (t937 - t939) * qJD(2) + (-t929 - t930) * t994 + t979;
t936 = qJD(2) * mrSges(3,1) - mrSges(3,3) * t994;
t833 = (mrSges(3,3) + mrSges(4,1)) * t932 + m(3) * t890 - qJD(2) * t936 - qJDD(2) * mrSges(3,2) + t930 * t996 + t970;
t987 = -t963 * t819 + t966 * t833;
t811 = m(2) * t942 - t969 * mrSges(2,1) - qJDD(1) * mrSges(2,2) + t987;
t914 = t982 - t1015;
t1004 = t960 * t824 - t959 * t826;
t867 = -t932 * pkin(2) - t1015 + t974;
t981 = -m(4) * t867 - t932 * mrSges(4,2) + t940 * t994 - t1004;
t972 = -m(3) * t914 + t937 * t996 + t932 * mrSges(3,1) + (-mrSges(3,2) + mrSges(4,3)) * t931 + (-t936 * t963 - t939 * t966) * qJD(1) + t981;
t816 = m(2) * t941 + qJDD(1) * mrSges(2,1) - t969 * mrSges(2,2) + t972;
t1005 = t964 * t811 + t967 * t816;
t813 = t966 * t819 + t963 * t833;
t999 = t1021 * qJD(2) + (t1009 * t966 + t1011 * t963) * qJD(1);
t988 = t967 * t811 - t964 * t816;
t820 = -t931 * mrSges(4,3) + t939 * t996 - t981;
t805 = -mrSges(3,1) * t914 - mrSges(4,1) * t868 + mrSges(4,2) * t867 + mrSges(3,3) * t890 - pkin(2) * t820 + pkin(3) * t835 - qJ(4) * t1004 + t997 * qJD(2) + t1009 * qJDD(2) + t1013 * t931 + t1025 * t932 - t960 * t808 - t959 * t814 - t999 * t994;
t827 = t898 * mrSges(6,2) + t921 * t885 - t976;
t973 = mrSges(7,1) * t840 - mrSges(7,2) * t841 + Ifges(7,5) * t857 + Ifges(7,6) * t856 + Ifges(7,3) * t927 + t883 * t862 - t882 * t863;
t807 = -t973 + t1010 * t898 + (-qJ(5) * mrSges(6,2) - t1008) * t897 - t998 * qJD(2) - t1003 * t921 + (-qJ(5) * t885 + t1001) * t920 + t999 * t996 - pkin(4) * t827 + t1013 * t932 + (-t1022 + t1024) * t931 + mrSges(3,2) * t914 - pkin(5) * t828 - qJ(3) * t820 - mrSges(4,3) * t867 + mrSges(4,1) * t869 - mrSges(3,3) * t889 + t1011 * qJDD(2) + pkin(3) * t822 + mrSges(6,3) * t845 - mrSges(6,1) * t846 + mrSges(5,1) * t848 - mrSges(5,2) * t849 + qJ(5) * t980;
t978 = mrSges(2,1) * t941 - mrSges(2,2) * t942 + Ifges(2,3) * qJDD(1) + pkin(1) * t972 + pkin(7) * t987 + t966 * t805 + t963 * t807;
t803 = mrSges(2,1) * g(3) + mrSges(2,3) * t942 + t969 * Ifges(2,5) + Ifges(2,6) * qJDD(1) - pkin(1) * t813 - t1020;
t802 = -mrSges(2,2) * g(3) - mrSges(2,3) * t941 + Ifges(2,5) * qJDD(1) - t969 * Ifges(2,6) - pkin(7) * t813 - t963 * t805 + t966 * t807;
t1 = [-m(1) * g(1) + t988; -m(1) * g(2) + t1005; (-m(1) - m(2)) * g(3) + t813; -mrSges(1,2) * g(3) + mrSges(1,3) * g(2) - pkin(6) * t1005 + t967 * t802 - t964 * t803; mrSges(1,1) * g(3) - mrSges(1,3) * g(1) + pkin(6) * t988 + t964 * t802 + t967 * t803; -mrSges(1,1) * g(2) + mrSges(1,2) * g(1) + t978; t978; t1020; t821; t835; t827; t973;];
tauJB  = t1;
