% Calculate vector of inverse dynamics joint torques and base forces with Newton-Euler
% S6RRRRPR5
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
%   pkin=[a2,a3,a4,a5,a6,d1,d2,d3,d4,d6]';
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
% Datum: 2019-05-07 20:33
% Revision: 89c353f7eff3bd693eda4e29f35b2761dbc3ada0 (2019-05-03)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauJB = S6RRRRPR5_invdynJB_fixb_snew_vp2(qJ, qJD, qJDD, g, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(6,1),zeros(3,1),zeros(10,1),zeros(7,1),zeros(7,3),zeros(7,6)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRRRPR5_invdynJB_fixb_snew_vp2: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RRRRPR5_invdynJB_fixb_snew_vp2: qJD has to be [6x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [6 1]), ...
  'S6RRRRPR5_invdynJB_fixb_snew_vp2: qJDD has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6RRRRPR5_invdynJB_fixb_snew_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S6RRRRPR5_invdynJB_fixb_snew_vp2: pkin has to be [10x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6RRRRPR5_invdynJB_fixb_snew_vp2: m has to be [7x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [7,3]), ...
  'S6RRRRPR5_invdynJB_fixb_snew_vp2: mrSges has to be [7x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [7 6]), ...
  'S6RRRRPR5_invdynJB_fixb_snew_vp2: Ifges has to be [7x6] (double)'); 

%% Symbolic Calculation
% From invdyn_fixb_NewtonEuler_linkframe_tauJB_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-05-07 20:26:57
% EndTime: 2019-05-07 20:27:19
% DurationCPUTime: 12.84s
% Computational Cost: add. (200147->369), mult. (399619->446), div. (0->0), fcn. (284345->10), ass. (0->150)
t1023 = Ifges(5,1) + Ifges(6,1);
t1013 = Ifges(5,4) - Ifges(6,5);
t1012 = Ifges(5,5) + Ifges(6,4);
t1022 = -Ifges(5,2) - Ifges(6,3);
t1011 = Ifges(5,6) - Ifges(6,6);
t1021 = Ifges(5,3) + Ifges(6,2);
t1016 = cos(qJ(4));
t972 = sin(qJ(3));
t973 = sin(qJ(2));
t976 = cos(qJ(3));
t977 = cos(qJ(2));
t944 = (t972 * t977 + t973 * t976) * qJD(1);
t967 = qJD(2) + qJD(3);
t971 = sin(qJ(4));
t929 = -t1016 * t967 + t971 * t944;
t1001 = qJD(1) * t977;
t1002 = qJD(1) * t973;
t943 = t976 * t1001 - t1002 * t972;
t939 = qJD(4) - t943;
t1009 = t929 * t939;
t1000 = qJD(1) * qJD(2);
t952 = qJDD(1) * t973 + t1000 * t977;
t953 = qJDD(1) * t977 - t1000 * t973;
t915 = qJD(3) * t943 + t952 * t976 + t953 * t972;
t966 = qJDD(2) + qJDD(3);
t877 = -t929 * qJD(4) + t1016 * t915 + t971 * t966;
t974 = sin(qJ(1));
t978 = cos(qJ(1));
t959 = -g(1) * t978 - g(2) * t974;
t979 = qJD(1) ^ 2;
t946 = -pkin(1) * t979 + qJDD(1) * pkin(7) + t959;
t1008 = t973 * t946;
t1015 = pkin(2) * t979;
t902 = qJDD(2) * pkin(2) - t952 * pkin(8) - t1008 + (pkin(8) * t1000 + t1015 * t973 - g(3)) * t977;
t932 = -g(3) * t973 + t977 * t946;
t957 = qJD(2) * pkin(2) - pkin(8) * t1002;
t969 = t977 ^ 2;
t903 = pkin(8) * t953 - qJD(2) * t957 - t1015 * t969 + t932;
t872 = t976 * t902 - t972 * t903;
t927 = -pkin(3) * t943 - pkin(9) * t944;
t965 = t967 ^ 2;
t989 = t966 * pkin(3) + t965 * pkin(9) - t944 * t927 + t872;
t1020 = (-t877 + t1009) * qJ(5) - t989;
t873 = t972 * t902 + t976 * t903;
t914 = -qJD(3) * t944 - t952 * t972 + t976 * t953;
t926 = -mrSges(4,1) * t943 + mrSges(4,2) * t944;
t934 = mrSges(4,1) * t967 - mrSges(4,3) * t944;
t930 = t1016 * t944 + t971 * t967;
t899 = mrSges(6,1) * t929 - mrSges(6,3) * t930;
t1003 = -mrSges(5,1) * t929 - mrSges(5,2) * t930 - t899;
t1014 = -mrSges(5,3) - mrSges(6,2);
t958 = t974 * g(1) - t978 * g(2);
t991 = -qJDD(1) * pkin(1) - t958;
t916 = -t953 * pkin(2) + t957 * t1002 + (-pkin(8) * t969 - pkin(7)) * t979 + t991;
t860 = (-t943 * t967 - t915) * pkin(9) + (t944 * t967 - t914) * pkin(3) + t916;
t864 = -pkin(3) * t965 + pkin(9) * t966 + t927 * t943 + t873;
t851 = t1016 * t864 + t971 * t860;
t876 = t930 * qJD(4) - t1016 * t966 + t971 * t915;
t913 = qJDD(4) - t914;
t919 = mrSges(5,1) * t939 - mrSges(5,3) * t930;
t1017 = 2 * qJD(5);
t898 = pkin(4) * t929 - qJ(5) * t930;
t938 = t939 ^ 2;
t847 = -pkin(4) * t938 + t913 * qJ(5) + t1017 * t939 - t929 * t898 + t851;
t920 = -mrSges(6,1) * t939 + mrSges(6,2) * t930;
t850 = t1016 * t860 - t971 * t864;
t848 = -t913 * pkin(4) - t938 * qJ(5) + t930 * t898 + qJDD(5) - t850;
t842 = (-t877 - t1009) * pkin(10) + (t929 * t930 - t913) * pkin(5) + t848;
t921 = -pkin(5) * t939 - pkin(10) * t930;
t928 = t929 ^ 2;
t843 = -pkin(5) * t928 + pkin(10) * t876 + t921 * t939 + t847;
t970 = sin(qJ(6));
t975 = cos(qJ(6));
t840 = t842 * t975 - t843 * t970;
t892 = t929 * t975 - t930 * t970;
t857 = qJD(6) * t892 + t876 * t970 + t877 * t975;
t893 = t929 * t970 + t930 * t975;
t870 = -mrSges(7,1) * t892 + mrSges(7,2) * t893;
t937 = qJD(6) - t939;
t880 = -mrSges(7,2) * t937 + mrSges(7,3) * t892;
t908 = qJDD(6) - t913;
t837 = m(7) * t840 + mrSges(7,1) * t908 - t857 * mrSges(7,3) - t870 * t893 + t880 * t937;
t841 = t842 * t970 + t843 * t975;
t856 = -qJD(6) * t893 + t876 * t975 - t877 * t970;
t881 = mrSges(7,1) * t937 - mrSges(7,3) * t893;
t838 = m(7) * t841 - mrSges(7,2) * t908 + t856 * mrSges(7,3) + t870 * t892 - t881 * t937;
t995 = -t970 * t837 + t975 * t838;
t990 = m(6) * t847 + t913 * mrSges(6,3) + t939 * t920 + t995;
t823 = m(5) * t851 - t913 * mrSges(5,2) + t1003 * t929 + t1014 * t876 - t939 * t919 + t990;
t918 = -mrSges(5,2) * t939 - mrSges(5,3) * t929;
t828 = t975 * t837 + t970 * t838;
t917 = -mrSges(6,2) * t929 + mrSges(6,3) * t939;
t987 = -m(6) * t848 + t913 * mrSges(6,1) + t939 * t917 - t828;
t825 = m(5) * t850 + t913 * mrSges(5,1) + t1003 * t930 + t1014 * t877 + t939 * t918 + t987;
t996 = t1016 * t823 - t825 * t971;
t817 = m(4) * t873 - mrSges(4,2) * t966 + mrSges(4,3) * t914 + t926 * t943 - t934 * t967 + t996;
t933 = -mrSges(4,2) * t967 + mrSges(4,3) * t943;
t849 = -0.2e1 * qJD(5) * t930 + (t930 * t939 + t876) * pkin(4) + t1020;
t845 = -t928 * pkin(10) + (-pkin(4) - pkin(5)) * t876 + (-pkin(4) * t939 + t1017 + t921) * t930 - t1020;
t993 = -m(7) * t845 + t856 * mrSges(7,1) - t857 * mrSges(7,2) + t892 * t880 - t893 * t881;
t835 = m(6) * t849 + mrSges(6,1) * t876 - t877 * mrSges(6,3) + t917 * t929 - t930 * t920 + t993;
t982 = m(5) * t989 - t876 * mrSges(5,1) - mrSges(5,2) * t877 - t929 * t918 - t919 * t930 - t835;
t832 = m(4) * t872 + mrSges(4,1) * t966 - mrSges(4,3) * t915 - t926 * t944 + t933 * t967 + t982;
t811 = t972 * t817 + t976 * t832;
t931 = -t977 * g(3) - t1008;
t941 = Ifges(3,6) * qJD(2) + (Ifges(3,4) * t973 + Ifges(3,2) * t977) * qJD(1);
t942 = Ifges(3,5) * qJD(2) + (Ifges(3,1) * t973 + Ifges(3,4) * t977) * qJD(1);
t1004 = t1012 * t939 - t1013 * t929 + t1023 * t930;
t1006 = t1011 * t929 - t1012 * t930 - t1021 * t939;
t865 = Ifges(7,5) * t893 + Ifges(7,6) * t892 + Ifges(7,3) * t937;
t867 = Ifges(7,1) * t893 + Ifges(7,4) * t892 + Ifges(7,5) * t937;
t829 = -mrSges(7,1) * t845 + mrSges(7,3) * t841 + Ifges(7,4) * t857 + Ifges(7,2) * t856 + Ifges(7,6) * t908 - t865 * t893 + t867 * t937;
t866 = Ifges(7,4) * t893 + Ifges(7,2) * t892 + Ifges(7,6) * t937;
t830 = mrSges(7,2) * t845 - mrSges(7,3) * t840 + Ifges(7,1) * t857 + Ifges(7,4) * t856 + Ifges(7,5) * t908 + t865 * t892 - t866 * t937;
t803 = mrSges(5,1) * t989 - mrSges(6,1) * t849 + mrSges(6,2) * t847 + mrSges(5,3) * t851 - pkin(4) * t835 - pkin(5) * t993 - pkin(10) * t995 + t1004 * t939 + t1006 * t930 + t1011 * t913 + t1013 * t877 + t1022 * t876 - t975 * t829 - t970 * t830;
t1005 = t1011 * t939 + t1013 * t930 + t1022 * t929;
t805 = -mrSges(5,2) * t989 + mrSges(6,2) * t848 - mrSges(5,3) * t850 - mrSges(6,3) * t849 - pkin(10) * t828 - qJ(5) * t835 - t1005 * t939 + t1006 * t929 + t1012 * t913 - t1013 * t876 + t1023 * t877 - t970 * t829 + t975 * t830;
t923 = Ifges(4,4) * t944 + Ifges(4,2) * t943 + Ifges(4,6) * t967;
t924 = Ifges(4,1) * t944 + Ifges(4,4) * t943 + Ifges(4,5) * t967;
t984 = -mrSges(4,1) * t872 + mrSges(4,2) * t873 - Ifges(4,5) * t915 - Ifges(4,6) * t914 - Ifges(4,3) * t966 - pkin(3) * t982 - pkin(9) * t996 - t1016 * t803 - t971 * t805 - t944 * t923 + t943 * t924;
t1019 = mrSges(3,1) * t931 - mrSges(3,2) * t932 + Ifges(3,5) * t952 + Ifges(3,6) * t953 + Ifges(3,3) * qJDD(2) + pkin(2) * t811 + (t941 * t973 - t942 * t977) * qJD(1) - t984;
t827 = t877 * mrSges(6,2) + t930 * t899 - t987;
t986 = -mrSges(7,1) * t840 + mrSges(7,2) * t841 - Ifges(7,5) * t857 - Ifges(7,6) * t856 - Ifges(7,3) * t908 - t893 * t866 + t892 * t867;
t1018 = t1004 * t929 + t1005 * t930 + t1021 * t913 - t1011 * t876 + t1012 * t877 + mrSges(5,1) * t850 - mrSges(6,1) * t848 - mrSges(5,2) * t851 + mrSges(6,3) * t847 - pkin(4) * t827 - pkin(5) * t828 + qJ(5) * (-t876 * mrSges(6,2) - t929 * t899 + t990) + t986;
t951 = (-mrSges(3,1) * t977 + mrSges(3,2) * t973) * qJD(1);
t956 = -qJD(2) * mrSges(3,2) + mrSges(3,3) * t1001;
t809 = m(3) * t931 + qJDD(2) * mrSges(3,1) - mrSges(3,3) * t952 + qJD(2) * t956 - t1002 * t951 + t811;
t955 = qJD(2) * mrSges(3,1) - mrSges(3,3) * t1002;
t997 = t976 * t817 - t832 * t972;
t810 = m(3) * t932 - qJDD(2) * mrSges(3,2) + mrSges(3,3) * t953 - qJD(2) * t955 + t1001 * t951 + t997;
t998 = -t809 * t973 + t977 * t810;
t799 = m(2) * t959 - mrSges(2,1) * t979 - qJDD(1) * mrSges(2,2) + t998;
t945 = -t979 * pkin(7) + t991;
t819 = t1016 * t825 + t971 * t823;
t985 = m(4) * t916 - t914 * mrSges(4,1) + mrSges(4,2) * t915 - t943 * t933 + t934 * t944 + t819;
t983 = -m(3) * t945 + t953 * mrSges(3,1) - mrSges(3,2) * t952 + t956 * t1001 - t1002 * t955 - t985;
t813 = m(2) * t958 + qJDD(1) * mrSges(2,1) - mrSges(2,2) * t979 + t983;
t1007 = t974 * t799 + t978 * t813;
t801 = t977 * t809 + t973 * t810;
t999 = t978 * t799 - t813 * t974;
t922 = Ifges(4,5) * t944 + Ifges(4,6) * t943 + Ifges(4,3) * t967;
t795 = mrSges(4,2) * t916 - mrSges(4,3) * t872 + Ifges(4,1) * t915 + Ifges(4,4) * t914 + Ifges(4,5) * t966 - pkin(9) * t819 + t1016 * t805 - t971 * t803 + t943 * t922 - t967 * t923;
t796 = -mrSges(4,1) * t916 + mrSges(4,3) * t873 + Ifges(4,4) * t915 + Ifges(4,2) * t914 + Ifges(4,6) * t966 - pkin(3) * t819 - t944 * t922 + t967 * t924 - t1018;
t940 = Ifges(3,3) * qJD(2) + (Ifges(3,5) * t973 + Ifges(3,6) * t977) * qJD(1);
t791 = -mrSges(3,1) * t945 + mrSges(3,3) * t932 + Ifges(3,4) * t952 + Ifges(3,2) * t953 + Ifges(3,6) * qJDD(2) - pkin(2) * t985 + pkin(8) * t997 + qJD(2) * t942 - t1002 * t940 + t972 * t795 + t976 * t796;
t793 = mrSges(3,2) * t945 - mrSges(3,3) * t931 + Ifges(3,1) * t952 + Ifges(3,4) * t953 + Ifges(3,5) * qJDD(2) - pkin(8) * t811 - qJD(2) * t941 + t1001 * t940 + t795 * t976 - t796 * t972;
t988 = mrSges(2,1) * t958 - mrSges(2,2) * t959 + Ifges(2,3) * qJDD(1) + pkin(1) * t983 + pkin(7) * t998 + t977 * t791 + t973 * t793;
t794 = mrSges(2,1) * g(3) + mrSges(2,3) * t959 + t979 * Ifges(2,5) + Ifges(2,6) * qJDD(1) - pkin(1) * t801 - t1019;
t789 = -mrSges(2,2) * g(3) - mrSges(2,3) * t958 + Ifges(2,5) * qJDD(1) - Ifges(2,6) * t979 - pkin(7) * t801 - t791 * t973 + t793 * t977;
t1 = [-m(1) * g(1) + t999; -m(1) * g(2) + t1007; (-m(1) - m(2)) * g(3) + t801; -mrSges(1,2) * g(3) + mrSges(1,3) * g(2) - pkin(6) * t1007 + t978 * t789 - t974 * t794; mrSges(1,1) * g(3) - mrSges(1,3) * g(1) + pkin(6) * t999 + t974 * t789 + t978 * t794; -mrSges(1,1) * g(2) + mrSges(1,2) * g(1) + t988; t988; t1019; -t984; t1018; t827; -t986;];
tauJB  = t1;
