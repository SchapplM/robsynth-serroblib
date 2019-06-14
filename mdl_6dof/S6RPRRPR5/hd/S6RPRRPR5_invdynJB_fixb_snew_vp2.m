% Calculate vector of inverse dynamics joint torques and base forces with Newton-Euler
% S6RPRRPR5
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
%   pkin=[a2,a3,a4,a5,a6,d1,d3,d4,d6,theta2]';
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
% Datum: 2019-05-05 22:42
% Revision: 89c353f7eff3bd693eda4e29f35b2761dbc3ada0 (2019-05-03)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauJB = S6RPRRPR5_invdynJB_fixb_snew_vp2(qJ, qJD, qJDD, g, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(6,1),zeros(3,1),zeros(10,1),zeros(7,1),zeros(7,3),zeros(7,6)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPRRPR5_invdynJB_fixb_snew_vp2: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RPRRPR5_invdynJB_fixb_snew_vp2: qJD has to be [6x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [6 1]), ...
  'S6RPRRPR5_invdynJB_fixb_snew_vp2: qJDD has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6RPRRPR5_invdynJB_fixb_snew_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S6RPRRPR5_invdynJB_fixb_snew_vp2: pkin has to be [10x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6RPRRPR5_invdynJB_fixb_snew_vp2: m has to be [7x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [7,3]), ...
  'S6RPRRPR5_invdynJB_fixb_snew_vp2: mrSges has to be [7x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [7 6]), ...
  'S6RPRRPR5_invdynJB_fixb_snew_vp2: Ifges has to be [7x6] (double)'); 

%% Symbolic Calculation
% From invdyn_fixb_NewtonEuler_linkframe_tauJB_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-05-05 22:39:18
% EndTime: 2019-05-05 22:39:31
% DurationCPUTime: 12.63s
% Computational Cost: add. (173956->345), mult. (424264->413), div. (0->0), fcn. (327463->10), ass. (0->149)
t1005 = Ifges(5,4) + Ifges(6,6);
t1016 = -Ifges(5,2) - Ifges(6,3);
t1012 = Ifges(5,6) - Ifges(6,5);
t1007 = cos(qJ(4));
t959 = sin(pkin(10));
t960 = cos(pkin(10));
t963 = sin(qJ(3));
t966 = cos(qJ(3));
t982 = -t959 * t963 + t960 * t966;
t936 = t982 * qJD(1);
t983 = t959 * t966 + t960 * t963;
t937 = t983 * qJD(1);
t962 = sin(qJ(4));
t917 = -t1007 * t936 + t937 * t962;
t918 = t1007 * t937 + t962 * t936;
t957 = qJD(3) + qJD(4);
t1010 = t1005 * t918 + t1012 * t957 + t1016 * t917;
t1011 = Ifges(5,3) + Ifges(6,1);
t1013 = Ifges(5,5) - Ifges(6,4);
t904 = mrSges(6,1) * t917 - mrSges(6,3) * t957;
t954 = qJDD(3) + qJDD(4);
t1003 = t917 * t957;
t1006 = pkin(2) * t960;
t964 = sin(qJ(1));
t967 = cos(qJ(1));
t943 = -g(1) * t967 - g(2) * t964;
t968 = qJD(1) ^ 2;
t938 = -pkin(1) * t968 + qJDD(1) * qJ(2) + t943;
t995 = qJD(1) * qJD(2);
t993 = -t960 * g(3) - 0.2e1 * t959 * t995;
t908 = (-pkin(7) * qJDD(1) + t1006 * t968 - t938) * t959 + t993;
t956 = t960 ^ 2;
t1002 = t956 * t968;
t928 = -g(3) * t959 + (t938 + 0.2e1 * t995) * t960;
t994 = qJDD(1) * t960;
t910 = -pkin(2) * t1002 + pkin(7) * t994 + t928;
t880 = t966 * t908 - t963 * t910;
t996 = t936 * qJD(3);
t926 = qJDD(1) * t983 + t996;
t847 = (-t926 + t996) * pkin(8) + (t936 * t937 + qJDD(3)) * pkin(3) + t880;
t881 = t963 * t908 + t966 * t910;
t925 = -t937 * qJD(3) + qJDD(1) * t982;
t931 = qJD(3) * pkin(3) - pkin(8) * t937;
t935 = t936 ^ 2;
t855 = -pkin(3) * t935 + pkin(8) * t925 - qJD(3) * t931 + t881;
t843 = t1007 * t847 - t962 * t855;
t894 = pkin(4) * t917 - qJ(5) * t918;
t953 = t957 ^ 2;
t839 = -t954 * pkin(4) - t953 * qJ(5) + t918 * t894 + qJDD(5) - t843;
t878 = -t917 * qJD(4) + t1007 * t926 + t962 * t925;
t833 = (t917 * t918 - t954) * pkin(9) + (t878 + t1003) * pkin(5) + t839;
t877 = qJD(4) * t918 - t1007 * t925 + t926 * t962;
t906 = pkin(5) * t918 - pkin(9) * t957;
t913 = t917 ^ 2;
t1008 = -2 * qJD(5);
t955 = t959 ^ 2;
t942 = t964 * g(1) - t967 * g(2);
t987 = qJDD(2) - t942;
t924 = (-pkin(1) - t1006) * qJDD(1) + (-qJ(2) + (-t955 - t956) * pkin(7)) * t968 + t987;
t863 = -t925 * pkin(3) - t935 * pkin(8) + t937 * t931 + t924;
t970 = (-t878 + t1003) * qJ(5) + t863 + (t957 * pkin(4) + t1008) * t918;
t834 = t970 - t918 * t906 - t913 * pkin(5) + (pkin(4) + pkin(9)) * t877;
t961 = sin(qJ(6));
t965 = cos(qJ(6));
t831 = t833 * t965 - t834 * t961;
t898 = t917 * t965 - t957 * t961;
t852 = qJD(6) * t898 + t877 * t961 + t954 * t965;
t876 = qJDD(6) + t878;
t899 = t917 * t961 + t957 * t965;
t879 = -mrSges(7,1) * t898 + mrSges(7,2) * t899;
t912 = qJD(6) + t918;
t882 = -mrSges(7,2) * t912 + mrSges(7,3) * t898;
t828 = m(7) * t831 + mrSges(7,1) * t876 - t852 * mrSges(7,3) - t879 * t899 + t882 * t912;
t832 = t833 * t961 + t834 * t965;
t851 = -qJD(6) * t899 + t877 * t965 - t954 * t961;
t883 = mrSges(7,1) * t912 - mrSges(7,3) * t899;
t829 = m(7) * t832 - mrSges(7,2) * t876 + t851 * mrSges(7,3) + t879 * t898 - t883 * t912;
t817 = t965 * t828 + t961 * t829;
t896 = -mrSges(6,2) * t917 - mrSges(6,3) * t918;
t978 = -m(6) * t839 - t878 * mrSges(6,1) - t918 * t896 - t817;
t816 = t954 * mrSges(6,2) + t957 * t904 - t978;
t844 = t1007 * t855 + t962 * t847;
t977 = -t953 * pkin(4) + t954 * qJ(5) - t917 * t894 + t844;
t836 = -t877 * pkin(5) - t913 * pkin(9) + ((2 * qJD(5)) + t906) * t957 + t977;
t856 = Ifges(7,5) * t899 + Ifges(7,6) * t898 + Ifges(7,3) * t912;
t858 = Ifges(7,1) * t899 + Ifges(7,4) * t898 + Ifges(7,5) * t912;
t819 = -mrSges(7,1) * t836 + mrSges(7,3) * t832 + Ifges(7,4) * t852 + Ifges(7,2) * t851 + Ifges(7,6) * t876 - t856 * t899 + t858 * t912;
t857 = Ifges(7,4) * t899 + Ifges(7,2) * t898 + Ifges(7,6) * t912;
t820 = mrSges(7,2) * t836 - mrSges(7,3) * t831 + Ifges(7,1) * t852 + Ifges(7,4) * t851 + Ifges(7,5) * t876 + t856 * t898 - t857 * t912;
t837 = t1008 * t957 - t977;
t905 = mrSges(6,1) * t918 + mrSges(6,2) * t957;
t980 = -m(7) * t836 + t851 * mrSges(7,1) - t852 * mrSges(7,2) + t898 * t882 - t899 * t883;
t974 = -m(6) * t837 + t954 * mrSges(6,3) + t957 * t905 - t980;
t1014 = Ifges(5,1) + Ifges(6,2);
t998 = -t1005 * t917 + t1013 * t957 + t1014 * t918;
t1015 = -mrSges(5,2) * t844 - mrSges(6,3) * t837 - pkin(4) * t816 - pkin(9) * t817 - t961 * t819 + t965 * t820 + qJ(5) * (-t917 * t896 + t974) + mrSges(6,2) * t839 + mrSges(5,1) * t843 + t1011 * t954 + t1010 * t918 + t1013 * t878 + (-qJ(5) * mrSges(6,1) - t1012) * t877 + t998 * t917;
t895 = mrSges(5,1) * t917 + mrSges(5,2) * t918;
t902 = -mrSges(5,2) * t957 - mrSges(5,3) * t917;
t813 = m(5) * t843 - t878 * mrSges(5,3) - t918 * t895 + (t902 - t904) * t957 + (mrSges(5,1) - mrSges(6,2)) * t954 + t978;
t903 = mrSges(5,1) * t957 - mrSges(5,3) * t918;
t823 = m(5) * t844 - t954 * mrSges(5,2) - t957 * t903 + (-t895 - t896) * t917 + (-mrSges(5,3) - mrSges(6,1)) * t877 + t974;
t807 = t1007 * t813 + t962 * t823;
t915 = Ifges(4,4) * t937 + Ifges(4,2) * t936 + Ifges(4,6) * qJD(3);
t916 = Ifges(4,1) * t937 + Ifges(4,4) * t936 + Ifges(4,5) * qJD(3);
t1009 = mrSges(4,1) * t880 - mrSges(4,2) * t881 + Ifges(4,5) * t926 + Ifges(4,6) * t925 + Ifges(4,3) * qJDD(3) + pkin(3) * t807 + t937 * t915 - t936 * t916 + t1015;
t1004 = mrSges(3,2) * t959;
t922 = -mrSges(4,1) * t936 + mrSges(4,2) * t937;
t929 = -qJD(3) * mrSges(4,2) + mrSges(4,3) * t936;
t805 = m(4) * t880 + qJDD(3) * mrSges(4,1) - mrSges(4,3) * t926 + qJD(3) * t929 - t922 * t937 + t807;
t930 = qJD(3) * mrSges(4,1) - mrSges(4,3) * t937;
t988 = t1007 * t823 - t813 * t962;
t806 = m(4) * t881 - qJDD(3) * mrSges(4,2) + mrSges(4,3) * t925 - qJD(3) * t930 + t922 * t936 + t988;
t800 = t966 * t805 + t963 * t806;
t927 = -t959 * t938 + t993;
t981 = mrSges(3,3) * qJDD(1) + t968 * (-mrSges(3,1) * t960 + t1004);
t798 = m(3) * t927 - t959 * t981 + t800;
t989 = -t963 * t805 + t966 * t806;
t799 = m(3) * t928 + t960 * t981 + t989;
t990 = -t798 * t959 + t960 * t799;
t791 = m(2) * t943 - mrSges(2,1) * t968 - qJDD(1) * mrSges(2,2) + t990;
t934 = -qJDD(1) * pkin(1) - t968 * qJ(2) + t987;
t1000 = -t961 * t828 + t965 * t829;
t841 = t877 * pkin(4) + t970;
t814 = m(6) * t841 - t877 * mrSges(6,2) - t878 * mrSges(6,3) - t917 * t904 - t918 * t905 + t1000;
t975 = m(5) * t863 + t877 * mrSges(5,1) + t878 * mrSges(5,2) + t917 * t902 + t918 * t903 + t814;
t973 = m(4) * t924 - t925 * mrSges(4,1) + t926 * mrSges(4,2) - t936 * t929 + t937 * t930 + t975;
t971 = -m(3) * t934 + mrSges(3,1) * t994 - t973 + (t955 * t968 + t1002) * mrSges(3,3);
t809 = t971 + (mrSges(2,1) - t1004) * qJDD(1) - t968 * mrSges(2,2) + m(2) * t942;
t1001 = t964 * t791 + t967 * t809;
t793 = t960 * t798 + t959 * t799;
t999 = -t1011 * t957 + t1012 * t917 - t1013 * t918;
t984 = Ifges(3,5) * t959 + Ifges(3,6) * t960;
t997 = t968 * t984;
t991 = t967 * t791 - t809 * t964;
t986 = Ifges(3,1) * t959 + Ifges(3,4) * t960;
t985 = Ifges(3,4) * t959 + Ifges(3,2) * t960;
t794 = -mrSges(5,1) * t863 - mrSges(6,1) * t837 + mrSges(6,2) * t841 + mrSges(5,3) * t844 - pkin(4) * t814 - pkin(5) * t980 - pkin(9) * t1000 + t1005 * t878 + t1012 * t954 + t1016 * t877 - t965 * t819 - t961 * t820 + t999 * t918 + t998 * t957;
t976 = mrSges(7,1) * t831 - mrSges(7,2) * t832 + Ifges(7,5) * t852 + Ifges(7,6) * t851 + Ifges(7,3) * t876 + t899 * t857 - t898 * t858;
t801 = mrSges(6,1) * t839 + mrSges(5,2) * t863 - mrSges(5,3) * t843 - mrSges(6,3) * t841 + pkin(5) * t817 - qJ(5) * t814 - t1005 * t877 - t1010 * t957 + t1013 * t954 + t1014 * t878 + t999 * t917 + t976;
t914 = Ifges(4,5) * t937 + Ifges(4,6) * t936 + Ifges(4,3) * qJD(3);
t787 = -mrSges(4,1) * t924 + mrSges(4,3) * t881 + Ifges(4,4) * t926 + Ifges(4,2) * t925 + Ifges(4,6) * qJDD(3) - pkin(3) * t975 + pkin(8) * t988 + qJD(3) * t916 + t1007 * t794 + t962 * t801 - t937 * t914;
t788 = mrSges(4,2) * t924 - mrSges(4,3) * t880 + Ifges(4,1) * t926 + Ifges(4,4) * t925 + Ifges(4,5) * qJDD(3) - pkin(8) * t807 - qJD(3) * t915 + t1007 * t801 - t962 * t794 + t936 * t914;
t783 = -mrSges(3,1) * t934 + mrSges(3,3) * t928 - pkin(2) * t973 + pkin(7) * t989 + qJDD(1) * t985 + t966 * t787 + t963 * t788 - t959 * t997;
t785 = mrSges(3,2) * t934 - mrSges(3,3) * t927 - pkin(7) * t800 + qJDD(1) * t986 - t963 * t787 + t966 * t788 + t960 * t997;
t811 = qJDD(1) * t1004 - t971;
t979 = mrSges(2,1) * t942 - mrSges(2,2) * t943 + Ifges(2,3) * qJDD(1) - pkin(1) * t811 + qJ(2) * t990 + t960 * t783 + t959 * t785;
t786 = -pkin(2) * t800 + mrSges(2,1) * g(3) - pkin(1) * t793 + (Ifges(2,6) - t984) * qJDD(1) + mrSges(2,3) * t943 - mrSges(3,1) * t927 + mrSges(3,2) * t928 + (-t959 * t985 + t960 * t986 + Ifges(2,5)) * t968 - t1009;
t781 = -mrSges(2,2) * g(3) - mrSges(2,3) * t942 + Ifges(2,5) * qJDD(1) - Ifges(2,6) * t968 - qJ(2) * t793 - t783 * t959 + t785 * t960;
t1 = [-m(1) * g(1) + t991; -m(1) * g(2) + t1001; (-m(1) - m(2)) * g(3) + t793; -mrSges(1,2) * g(3) + mrSges(1,3) * g(2) - pkin(6) * t1001 + t967 * t781 - t964 * t786; mrSges(1,1) * g(3) - mrSges(1,3) * g(1) + pkin(6) * t991 + t964 * t781 + t967 * t786; -mrSges(1,1) * g(2) + mrSges(1,2) * g(1) + t979; t979; t811; t1009; t1015; t816; t976;];
tauJB  = t1;
