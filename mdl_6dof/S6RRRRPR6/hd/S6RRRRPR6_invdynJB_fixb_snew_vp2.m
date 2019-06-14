% Calculate vector of inverse dynamics joint torques and base forces with Newton-Euler
% S6RRRRPR6
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
%   pkin=[a2,a3,a4,a5,a6,d1,d2,d3,d4,d6,theta5]';
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
% Datum: 2019-05-07 20:55
% Revision: 89c353f7eff3bd693eda4e29f35b2761dbc3ada0 (2019-05-03)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauJB = S6RRRRPR6_invdynJB_fixb_snew_vp2(qJ, qJD, qJDD, g, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(6,1),zeros(3,1),zeros(11,1),zeros(7,1),zeros(7,3),zeros(7,6)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRRRPR6_invdynJB_fixb_snew_vp2: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RRRRPR6_invdynJB_fixb_snew_vp2: qJD has to be [6x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [6 1]), ...
  'S6RRRRPR6_invdynJB_fixb_snew_vp2: qJDD has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6RRRRPR6_invdynJB_fixb_snew_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [11 1]), ...
  'S6RRRRPR6_invdynJB_fixb_snew_vp2: pkin has to be [11x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6RRRRPR6_invdynJB_fixb_snew_vp2: m has to be [7x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [7,3]), ...
  'S6RRRRPR6_invdynJB_fixb_snew_vp2: mrSges has to be [7x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [7 6]), ...
  'S6RRRRPR6_invdynJB_fixb_snew_vp2: Ifges has to be [7x6] (double)'); 

%% Symbolic Calculation
% From invdyn_fixb_NewtonEuler_linkframe_tauJB_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-05-07 20:48:45
% EndTime: 2019-05-07 20:49:54
% DurationCPUTime: 50.30s
% Computational Cost: add. (824314->387), mult. (1723623->484), div. (0->0), fcn. (1279920->12), ass. (0->153)
t954 = sin(qJ(1));
t959 = cos(qJ(1));
t940 = t954 * g(1) - t959 * g(2);
t961 = qJD(1) ^ 2;
t923 = -qJDD(1) * pkin(1) - t961 * pkin(7) - t940;
t953 = sin(qJ(2));
t958 = cos(qJ(2));
t977 = qJD(1) * qJD(2);
t976 = t958 * t977;
t934 = qJDD(1) * t953 + t976;
t944 = t953 * t977;
t935 = qJDD(1) * t958 - t944;
t889 = (-t934 - t976) * pkin(8) + (-t935 + t944) * pkin(2) + t923;
t941 = -g(1) * t959 - g(2) * t954;
t924 = -pkin(1) * t961 + qJDD(1) * pkin(7) + t941;
t911 = -g(3) * t953 + t958 * t924;
t933 = (-pkin(2) * t958 - pkin(8) * t953) * qJD(1);
t960 = qJD(2) ^ 2;
t978 = qJD(1) * t958;
t892 = -pkin(2) * t960 + qJDD(2) * pkin(8) + t933 * t978 + t911;
t952 = sin(qJ(3));
t957 = cos(qJ(3));
t871 = t957 * t889 - t952 * t892;
t979 = qJD(1) * t953;
t930 = qJD(2) * t957 - t952 * t979;
t903 = qJD(3) * t930 + qJDD(2) * t952 + t934 * t957;
t929 = qJDD(3) - t935;
t931 = qJD(2) * t952 + t957 * t979;
t943 = qJD(3) - t978;
t851 = (t930 * t943 - t903) * pkin(9) + (t930 * t931 + t929) * pkin(3) + t871;
t872 = t952 * t889 + t957 * t892;
t902 = -qJD(3) * t931 + qJDD(2) * t957 - t934 * t952;
t912 = pkin(3) * t943 - pkin(9) * t931;
t928 = t930 ^ 2;
t853 = -pkin(3) * t928 + pkin(9) * t902 - t912 * t943 + t872;
t951 = sin(qJ(4));
t956 = cos(qJ(4));
t832 = t956 * t851 - t951 * t853;
t905 = t930 * t956 - t931 * t951;
t870 = qJD(4) * t905 + t902 * t951 + t903 * t956;
t906 = t930 * t951 + t931 * t956;
t925 = qJDD(4) + t929;
t942 = qJD(4) + t943;
t821 = (t905 * t942 - t870) * qJ(5) + (t905 * t906 + t925) * pkin(4) + t832;
t833 = t951 * t851 + t956 * t853;
t869 = -qJD(4) * t906 + t902 * t956 - t903 * t951;
t894 = pkin(4) * t942 - qJ(5) * t906;
t904 = t905 ^ 2;
t829 = -pkin(4) * t904 + qJ(5) * t869 - t894 * t942 + t833;
t948 = sin(pkin(11));
t949 = cos(pkin(11));
t885 = t905 * t948 + t906 * t949;
t815 = -0.2e1 * qJD(5) * t885 + t949 * t821 - t948 * t829;
t848 = t869 * t948 + t870 * t949;
t884 = t905 * t949 - t906 * t948;
t812 = (t884 * t942 - t848) * pkin(10) + (t884 * t885 + t925) * pkin(5) + t815;
t816 = 0.2e1 * qJD(5) * t884 + t948 * t821 + t949 * t829;
t847 = t869 * t949 - t870 * t948;
t875 = pkin(5) * t942 - pkin(10) * t885;
t883 = t884 ^ 2;
t813 = -pkin(5) * t883 + pkin(10) * t847 - t875 * t942 + t816;
t950 = sin(qJ(6));
t955 = cos(qJ(6));
t811 = t812 * t950 + t813 * t955;
t910 = -t958 * g(3) - t953 * t924;
t891 = -qJDD(2) * pkin(2) - pkin(8) * t960 + t933 * t979 - t910;
t864 = -pkin(3) * t902 - pkin(9) * t928 + t931 * t912 + t891;
t835 = -pkin(4) * t869 - qJ(5) * t904 + t906 * t894 + qJDD(5) + t864;
t818 = -pkin(5) * t847 - pkin(10) * t883 + t875 * t885 + t835;
t862 = t884 * t950 + t885 * t955;
t826 = -qJD(6) * t862 + t847 * t955 - t848 * t950;
t861 = t884 * t955 - t885 * t950;
t827 = qJD(6) * t861 + t847 * t950 + t848 * t955;
t937 = qJD(6) + t942;
t836 = Ifges(7,5) * t862 + Ifges(7,6) * t861 + Ifges(7,3) * t937;
t838 = Ifges(7,1) * t862 + Ifges(7,4) * t861 + Ifges(7,5) * t937;
t919 = qJDD(6) + t925;
t798 = -mrSges(7,1) * t818 + mrSges(7,3) * t811 + Ifges(7,4) * t827 + Ifges(7,2) * t826 + Ifges(7,6) * t919 - t836 * t862 + t838 * t937;
t810 = t812 * t955 - t813 * t950;
t837 = Ifges(7,4) * t862 + Ifges(7,2) * t861 + Ifges(7,6) * t937;
t799 = mrSges(7,2) * t818 - mrSges(7,3) * t810 + Ifges(7,1) * t827 + Ifges(7,4) * t826 + Ifges(7,5) * t919 + t836 * t861 - t837 * t937;
t856 = Ifges(6,5) * t885 + Ifges(6,6) * t884 + Ifges(6,3) * t942;
t858 = Ifges(6,1) * t885 + Ifges(6,4) * t884 + Ifges(6,5) * t942;
t854 = -mrSges(7,2) * t937 + mrSges(7,3) * t861;
t855 = mrSges(7,1) * t937 - mrSges(7,3) * t862;
t970 = m(7) * t818 - t826 * mrSges(7,1) + t827 * mrSges(7,2) - t861 * t854 + t862 * t855;
t841 = -mrSges(7,1) * t861 + mrSges(7,2) * t862;
t803 = m(7) * t810 + mrSges(7,1) * t919 - mrSges(7,3) * t827 - t841 * t862 + t854 * t937;
t804 = m(7) * t811 - mrSges(7,2) * t919 + mrSges(7,3) * t826 + t841 * t861 - t855 * t937;
t971 = -t803 * t950 + t955 * t804;
t782 = -mrSges(6,1) * t835 + mrSges(6,3) * t816 + Ifges(6,4) * t848 + Ifges(6,2) * t847 + Ifges(6,6) * t925 - pkin(5) * t970 + pkin(10) * t971 + t955 * t798 + t950 * t799 - t885 * t856 + t942 * t858;
t797 = t955 * t803 + t950 * t804;
t857 = Ifges(6,4) * t885 + Ifges(6,2) * t884 + Ifges(6,6) * t942;
t783 = mrSges(6,2) * t835 - mrSges(6,3) * t815 + Ifges(6,1) * t848 + Ifges(6,4) * t847 + Ifges(6,5) * t925 - pkin(10) * t797 - t798 * t950 + t799 * t955 + t856 * t884 - t857 * t942;
t873 = -mrSges(6,2) * t942 + mrSges(6,3) * t884;
t874 = mrSges(6,1) * t942 - mrSges(6,3) * t885;
t808 = m(6) * t835 - t847 * mrSges(6,1) + t848 * mrSges(6,2) - t884 * t873 + t885 * t874 + t970;
t876 = Ifges(5,5) * t906 + Ifges(5,6) * t905 + Ifges(5,3) * t942;
t878 = Ifges(5,1) * t906 + Ifges(5,4) * t905 + Ifges(5,5) * t942;
t863 = -mrSges(6,1) * t884 + mrSges(6,2) * t885;
t794 = m(6) * t815 + mrSges(6,1) * t925 - mrSges(6,3) * t848 - t863 * t885 + t873 * t942 + t797;
t795 = m(6) * t816 - mrSges(6,2) * t925 + mrSges(6,3) * t847 + t863 * t884 - t874 * t942 + t971;
t972 = -t794 * t948 + t949 * t795;
t776 = -mrSges(5,1) * t864 + mrSges(5,3) * t833 + Ifges(5,4) * t870 + Ifges(5,2) * t869 + Ifges(5,6) * t925 - pkin(4) * t808 + qJ(5) * t972 + t949 * t782 + t948 * t783 - t906 * t876 + t942 * t878;
t790 = t949 * t794 + t948 * t795;
t877 = Ifges(5,4) * t906 + Ifges(5,2) * t905 + Ifges(5,6) * t942;
t777 = mrSges(5,2) * t864 - mrSges(5,3) * t832 + Ifges(5,1) * t870 + Ifges(5,4) * t869 + Ifges(5,5) * t925 - qJ(5) * t790 - t782 * t948 + t783 * t949 + t876 * t905 - t877 * t942;
t896 = Ifges(4,5) * t931 + Ifges(4,6) * t930 + Ifges(4,3) * t943;
t898 = Ifges(4,1) * t931 + Ifges(4,4) * t930 + Ifges(4,5) * t943;
t893 = -mrSges(5,2) * t942 + mrSges(5,3) * t905;
t895 = mrSges(5,1) * t942 - mrSges(5,3) * t906;
t965 = m(5) * t864 - t869 * mrSges(5,1) + t870 * mrSges(5,2) - t905 * t893 + t906 * t895 + t808;
t886 = -mrSges(5,1) * t905 + mrSges(5,2) * t906;
t787 = m(5) * t832 + mrSges(5,1) * t925 - mrSges(5,3) * t870 - t886 * t906 + t893 * t942 + t790;
t788 = m(5) * t833 - mrSges(5,2) * t925 + mrSges(5,3) * t869 + t886 * t905 - t895 * t942 + t972;
t973 = -t787 * t951 + t956 * t788;
t759 = -mrSges(4,1) * t891 + mrSges(4,3) * t872 + Ifges(4,4) * t903 + Ifges(4,2) * t902 + Ifges(4,6) * t929 - pkin(3) * t965 + pkin(9) * t973 + t956 * t776 + t951 * t777 - t931 * t896 + t943 * t898;
t781 = t956 * t787 + t951 * t788;
t897 = Ifges(4,4) * t931 + Ifges(4,2) * t930 + Ifges(4,6) * t943;
t760 = mrSges(4,2) * t891 - mrSges(4,3) * t871 + Ifges(4,1) * t903 + Ifges(4,4) * t902 + Ifges(4,5) * t929 - pkin(9) * t781 - t776 * t951 + t777 * t956 + t896 * t930 - t897 * t943;
t907 = -mrSges(4,1) * t930 + mrSges(4,2) * t931;
t908 = -mrSges(4,2) * t943 + mrSges(4,3) * t930;
t779 = m(4) * t871 + mrSges(4,1) * t929 - mrSges(4,3) * t903 - t907 * t931 + t908 * t943 + t781;
t909 = mrSges(4,1) * t943 - mrSges(4,3) * t931;
t780 = m(4) * t872 - mrSges(4,2) * t929 + mrSges(4,3) * t902 + t907 * t930 - t909 * t943 + t973;
t775 = -t779 * t952 + t957 * t780;
t807 = -m(4) * t891 + t902 * mrSges(4,1) - t903 * mrSges(4,2) + t930 * t908 - t931 * t909 - t965;
t921 = Ifges(3,6) * qJD(2) + (Ifges(3,4) * t953 + Ifges(3,2) * t958) * qJD(1);
t922 = Ifges(3,5) * qJD(2) + (Ifges(3,1) * t953 + Ifges(3,4) * t958) * qJD(1);
t981 = mrSges(3,1) * t910 - mrSges(3,2) * t911 + Ifges(3,5) * t934 + Ifges(3,6) * t935 + Ifges(3,3) * qJDD(2) + pkin(2) * t807 + pkin(8) * t775 + t957 * t759 + t952 * t760 + (t921 * t953 - t922 * t958) * qJD(1);
t932 = (-mrSges(3,1) * t958 + mrSges(3,2) * t953) * qJD(1);
t938 = qJD(2) * mrSges(3,1) - mrSges(3,3) * t979;
t773 = m(3) * t911 - qJDD(2) * mrSges(3,2) + mrSges(3,3) * t935 - qJD(2) * t938 + t932 * t978 + t775;
t939 = -qJD(2) * mrSges(3,2) + mrSges(3,3) * t978;
t806 = m(3) * t910 + qJDD(2) * mrSges(3,1) - t934 * mrSges(3,3) + qJD(2) * t939 - t932 * t979 + t807;
t974 = t958 * t773 - t806 * t953;
t765 = m(2) * t941 - mrSges(2,1) * t961 - qJDD(1) * mrSges(2,2) + t974;
t774 = t779 * t957 + t780 * t952;
t966 = -m(3) * t923 + t935 * mrSges(3,1) - mrSges(3,2) * t934 - t938 * t979 + t939 * t978 - t774;
t769 = m(2) * t940 + qJDD(1) * mrSges(2,1) - mrSges(2,2) * t961 + t966;
t980 = t954 * t765 + t959 * t769;
t767 = t953 * t773 + t958 * t806;
t975 = t959 * t765 - t769 * t954;
t920 = Ifges(3,3) * qJD(2) + (Ifges(3,5) * t953 + Ifges(3,6) * t958) * qJD(1);
t758 = mrSges(3,2) * t923 - mrSges(3,3) * t910 + Ifges(3,1) * t934 + Ifges(3,4) * t935 + Ifges(3,5) * qJDD(2) - pkin(8) * t774 - qJD(2) * t921 - t759 * t952 + t760 * t957 + t920 * t978;
t967 = -mrSges(7,1) * t810 + mrSges(7,2) * t811 - Ifges(7,5) * t827 - Ifges(7,6) * t826 - Ifges(7,3) * t919 - t862 * t837 + t861 * t838;
t963 = -mrSges(5,1) * t832 - mrSges(6,1) * t815 + mrSges(5,2) * t833 + mrSges(6,2) * t816 - Ifges(5,5) * t870 - Ifges(6,5) * t848 - Ifges(5,6) * t869 - Ifges(6,6) * t847 - pkin(4) * t790 - pkin(5) * t797 - t885 * t857 + t884 * t858 - t906 * t877 + t905 * t878 + t967 + (-Ifges(6,3) - Ifges(5,3)) * t925;
t962 = mrSges(4,1) * t871 - mrSges(4,2) * t872 + Ifges(4,5) * t903 + Ifges(4,6) * t902 + Ifges(4,3) * t929 + pkin(3) * t781 + t931 * t897 - t930 * t898 - t963;
t762 = -mrSges(3,1) * t923 + mrSges(3,3) * t911 + Ifges(3,4) * t934 + Ifges(3,2) * t935 + Ifges(3,6) * qJDD(2) - pkin(2) * t774 + qJD(2) * t922 - t920 * t979 - t962;
t968 = mrSges(2,1) * t940 - mrSges(2,2) * t941 + Ifges(2,3) * qJDD(1) + pkin(1) * t966 + pkin(7) * t974 + t953 * t758 + t958 * t762;
t756 = mrSges(2,1) * g(3) + mrSges(2,3) * t941 + t961 * Ifges(2,5) + Ifges(2,6) * qJDD(1) - pkin(1) * t767 - t981;
t755 = -mrSges(2,2) * g(3) - mrSges(2,3) * t940 + Ifges(2,5) * qJDD(1) - Ifges(2,6) * t961 - pkin(7) * t767 + t758 * t958 - t762 * t953;
t1 = [-m(1) * g(1) + t975; -m(1) * g(2) + t980; (-m(1) - m(2)) * g(3) + t767; -mrSges(1,2) * g(3) + mrSges(1,3) * g(2) - pkin(6) * t980 + t959 * t755 - t954 * t756; mrSges(1,1) * g(3) - mrSges(1,3) * g(1) + pkin(6) * t975 + t954 * t755 + t959 * t756; -mrSges(1,1) * g(2) + mrSges(1,2) * g(1) + t968; t968; t981; t962; -t963; t808; -t967;];
tauJB  = t1;
