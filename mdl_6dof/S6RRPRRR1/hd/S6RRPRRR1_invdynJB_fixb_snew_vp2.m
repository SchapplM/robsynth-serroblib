% Calculate vector of inverse dynamics joint torques and base forces with Newton-Euler
% S6RRPRRR1
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
%   pkin=[a2,a3,a4,a5,a6,d1,d2,d4,d5,d6,theta3]';
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
% Datum: 2019-05-06 19:35
% Revision: 89c353f7eff3bd693eda4e29f35b2761dbc3ada0 (2019-05-03)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauJB = S6RRPRRR1_invdynJB_fixb_snew_vp2(qJ, qJD, qJDD, g, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(6,1),zeros(3,1),zeros(11,1),zeros(7,1),zeros(7,3),zeros(7,6)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRPRRR1_invdynJB_fixb_snew_vp2: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RRPRRR1_invdynJB_fixb_snew_vp2: qJD has to be [6x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [6 1]), ...
  'S6RRPRRR1_invdynJB_fixb_snew_vp2: qJDD has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6RRPRRR1_invdynJB_fixb_snew_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [11 1]), ...
  'S6RRPRRR1_invdynJB_fixb_snew_vp2: pkin has to be [11x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6RRPRRR1_invdynJB_fixb_snew_vp2: m has to be [7x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [7,3]), ...
  'S6RRPRRR1_invdynJB_fixb_snew_vp2: mrSges has to be [7x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [7 6]), ...
  'S6RRPRRR1_invdynJB_fixb_snew_vp2: Ifges has to be [7x6] (double)'); 

%% Symbolic Calculation
% From invdyn_fixb_NewtonEuler_linkframe_tauJB_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-05-06 19:31:19
% EndTime: 2019-05-06 19:31:59
% DurationCPUTime: 39.70s
% Computational Cost: add. (608159->389), mult. (1465920->490), div. (0->0), fcn. (1123477->12), ass. (0->154)
t948 = sin(qJ(2));
t953 = cos(qJ(2));
t973 = qJD(1) * qJD(2);
t925 = qJDD(1) * t948 + t953 * t973;
t949 = sin(qJ(1));
t954 = cos(qJ(1));
t932 = -g(1) * t954 - g(2) * t949;
t955 = qJD(1) ^ 2;
t920 = -pkin(1) * t955 + qJDD(1) * pkin(7) + t932;
t977 = t948 * t920;
t979 = pkin(2) * t955;
t887 = qJDD(2) * pkin(2) - t925 * qJ(3) - t977 + (qJ(3) * t973 + t948 * t979 - g(3)) * t953;
t906 = -g(3) * t948 + t953 * t920;
t926 = qJDD(1) * t953 - t948 * t973;
t975 = qJD(1) * t948;
t928 = qJD(2) * pkin(2) - qJ(3) * t975;
t942 = t953 ^ 2;
t888 = qJ(3) * t926 - qJD(2) * t928 - t942 * t979 + t906;
t943 = sin(pkin(11));
t944 = cos(pkin(11));
t915 = (t943 * t953 + t944 * t948) * qJD(1);
t861 = -0.2e1 * qJD(3) * t915 + t944 * t887 - t943 * t888;
t904 = t925 * t944 + t926 * t943;
t914 = (-t943 * t948 + t944 * t953) * qJD(1);
t846 = (qJD(2) * t914 - t904) * pkin(8) + (t914 * t915 + qJDD(2)) * pkin(3) + t861;
t862 = 0.2e1 * qJD(3) * t914 + t943 * t887 + t944 * t888;
t903 = -t925 * t943 + t926 * t944;
t909 = qJD(2) * pkin(3) - pkin(8) * t915;
t913 = t914 ^ 2;
t848 = -pkin(3) * t913 + pkin(8) * t903 - qJD(2) * t909 + t862;
t947 = sin(qJ(4));
t952 = cos(qJ(4));
t825 = t952 * t846 - t947 * t848;
t897 = t914 * t952 - t915 * t947;
t869 = qJD(4) * t897 + t903 * t947 + t904 * t952;
t898 = t914 * t947 + t915 * t952;
t939 = qJDD(2) + qJDD(4);
t940 = qJD(2) + qJD(4);
t821 = (t897 * t940 - t869) * pkin(9) + (t897 * t898 + t939) * pkin(4) + t825;
t826 = t947 * t846 + t952 * t848;
t868 = -qJD(4) * t898 + t903 * t952 - t904 * t947;
t892 = pkin(4) * t940 - pkin(9) * t898;
t893 = t897 ^ 2;
t823 = -pkin(4) * t893 + pkin(9) * t868 - t892 * t940 + t826;
t946 = sin(qJ(5));
t951 = cos(qJ(5));
t818 = t946 * t821 + t951 * t823;
t882 = t897 * t946 + t898 * t951;
t837 = -qJD(5) * t882 + t868 * t951 - t869 * t946;
t881 = t897 * t951 - t898 * t946;
t857 = -mrSges(6,1) * t881 + mrSges(6,2) * t882;
t937 = qJD(5) + t940;
t873 = mrSges(6,1) * t937 - mrSges(6,3) * t882;
t936 = qJDD(5) + t939;
t858 = -pkin(5) * t881 - pkin(10) * t882;
t935 = t937 ^ 2;
t815 = -pkin(5) * t935 + pkin(10) * t936 + t858 * t881 + t818;
t931 = t949 * g(1) - t954 * g(2);
t965 = -qJDD(1) * pkin(1) - t931;
t889 = -t926 * pkin(2) + qJDD(3) + t928 * t975 + (-qJ(3) * t942 - pkin(7)) * t955 + t965;
t860 = -t903 * pkin(3) - t913 * pkin(8) + t915 * t909 + t889;
t831 = -t868 * pkin(4) - t893 * pkin(9) + t898 * t892 + t860;
t838 = qJD(5) * t881 + t868 * t946 + t869 * t951;
t819 = t831 + (-t881 * t937 - t838) * pkin(10) + (t882 * t937 - t837) * pkin(5);
t945 = sin(qJ(6));
t950 = cos(qJ(6));
t812 = -t815 * t945 + t819 * t950;
t870 = -t882 * t945 + t937 * t950;
t829 = qJD(6) * t870 + t838 * t950 + t936 * t945;
t836 = qJDD(6) - t837;
t871 = t882 * t950 + t937 * t945;
t849 = -mrSges(7,1) * t870 + mrSges(7,2) * t871;
t877 = qJD(6) - t881;
t850 = -mrSges(7,2) * t877 + mrSges(7,3) * t870;
t808 = m(7) * t812 + mrSges(7,1) * t836 - mrSges(7,3) * t829 - t849 * t871 + t850 * t877;
t813 = t815 * t950 + t819 * t945;
t828 = -qJD(6) * t871 - t838 * t945 + t936 * t950;
t851 = mrSges(7,1) * t877 - mrSges(7,3) * t871;
t809 = m(7) * t813 - mrSges(7,2) * t836 + mrSges(7,3) * t828 + t849 * t870 - t851 * t877;
t967 = -t808 * t945 + t950 * t809;
t795 = m(6) * t818 - mrSges(6,2) * t936 + mrSges(6,3) * t837 + t857 * t881 - t873 * t937 + t967;
t817 = t821 * t951 - t823 * t946;
t872 = -mrSges(6,2) * t937 + mrSges(6,3) * t881;
t814 = -pkin(5) * t936 - pkin(10) * t935 + t858 * t882 - t817;
t962 = -m(7) * t814 + t828 * mrSges(7,1) - mrSges(7,2) * t829 + t870 * t850 - t851 * t871;
t804 = m(6) * t817 + mrSges(6,1) * t936 - mrSges(6,3) * t838 - t857 * t882 + t872 * t937 + t962;
t788 = t946 * t795 + t951 * t804;
t883 = -mrSges(5,1) * t897 + mrSges(5,2) * t898;
t890 = -mrSges(5,2) * t940 + mrSges(5,3) * t897;
t785 = m(5) * t825 + mrSges(5,1) * t939 - mrSges(5,3) * t869 - t883 * t898 + t890 * t940 + t788;
t891 = mrSges(5,1) * t940 - mrSges(5,3) * t898;
t968 = t951 * t795 - t804 * t946;
t786 = m(5) * t826 - mrSges(5,2) * t939 + mrSges(5,3) * t868 + t883 * t897 - t891 * t940 + t968;
t779 = t952 * t785 + t947 * t786;
t901 = -mrSges(4,1) * t914 + mrSges(4,2) * t915;
t907 = -qJD(2) * mrSges(4,2) + mrSges(4,3) * t914;
t777 = m(4) * t861 + qJDD(2) * mrSges(4,1) - mrSges(4,3) * t904 + qJD(2) * t907 - t901 * t915 + t779;
t908 = qJD(2) * mrSges(4,1) - mrSges(4,3) * t915;
t969 = -t785 * t947 + t952 * t786;
t778 = m(4) * t862 - qJDD(2) * mrSges(4,2) + mrSges(4,3) * t903 - qJD(2) * t908 + t901 * t914 + t969;
t772 = t944 * t777 + t943 * t778;
t895 = Ifges(4,4) * t915 + Ifges(4,2) * t914 + Ifges(4,6) * qJD(2);
t896 = Ifges(4,1) * t915 + Ifges(4,4) * t914 + Ifges(4,5) * qJD(2);
t905 = -t953 * g(3) - t977;
t917 = Ifges(3,6) * qJD(2) + (Ifges(3,4) * t948 + Ifges(3,2) * t953) * qJD(1);
t918 = Ifges(3,5) * qJD(2) + (Ifges(3,1) * t948 + Ifges(3,4) * t953) * qJD(1);
t875 = Ifges(5,4) * t898 + Ifges(5,2) * t897 + Ifges(5,6) * t940;
t876 = Ifges(5,1) * t898 + Ifges(5,4) * t897 + Ifges(5,5) * t940;
t839 = Ifges(7,5) * t871 + Ifges(7,6) * t870 + Ifges(7,3) * t877;
t841 = Ifges(7,1) * t871 + Ifges(7,4) * t870 + Ifges(7,5) * t877;
t801 = -mrSges(7,1) * t814 + mrSges(7,3) * t813 + Ifges(7,4) * t829 + Ifges(7,2) * t828 + Ifges(7,6) * t836 - t839 * t871 + t841 * t877;
t840 = Ifges(7,4) * t871 + Ifges(7,2) * t870 + Ifges(7,6) * t877;
t802 = mrSges(7,2) * t814 - mrSges(7,3) * t812 + Ifges(7,1) * t829 + Ifges(7,4) * t828 + Ifges(7,5) * t836 + t839 * t870 - t840 * t877;
t853 = Ifges(6,4) * t882 + Ifges(6,2) * t881 + Ifges(6,6) * t937;
t854 = Ifges(6,1) * t882 + Ifges(6,4) * t881 + Ifges(6,5) * t937;
t960 = -mrSges(6,1) * t817 + mrSges(6,2) * t818 - Ifges(6,5) * t838 - Ifges(6,6) * t837 - Ifges(6,3) * t936 - pkin(5) * t962 - pkin(10) * t967 - t950 * t801 - t945 * t802 - t882 * t853 + t881 * t854;
t958 = -mrSges(5,1) * t825 + mrSges(5,2) * t826 - Ifges(5,5) * t869 - Ifges(5,6) * t868 - Ifges(5,3) * t939 - pkin(4) * t788 - t898 * t875 + t897 * t876 + t960;
t980 = mrSges(3,1) * t905 + mrSges(4,1) * t861 - mrSges(3,2) * t906 - mrSges(4,2) * t862 + Ifges(3,5) * t925 + Ifges(4,5) * t904 + Ifges(3,6) * t926 + Ifges(4,6) * t903 + pkin(2) * t772 + pkin(3) * t779 + (t917 * t948 - t918 * t953) * qJD(1) + (Ifges(3,3) + Ifges(4,3)) * qJDD(2) + t915 * t895 - t914 * t896 - t958;
t924 = (-mrSges(3,1) * t953 + mrSges(3,2) * t948) * qJD(1);
t974 = qJD(1) * t953;
t930 = -qJD(2) * mrSges(3,2) + mrSges(3,3) * t974;
t770 = m(3) * t905 + qJDD(2) * mrSges(3,1) - mrSges(3,3) * t925 + qJD(2) * t930 - t924 * t975 + t772;
t929 = qJD(2) * mrSges(3,1) - mrSges(3,3) * t975;
t970 = -t777 * t943 + t944 * t778;
t771 = m(3) * t906 - qJDD(2) * mrSges(3,2) + mrSges(3,3) * t926 - qJD(2) * t929 + t924 * t974 + t970;
t971 = -t770 * t948 + t953 * t771;
t763 = m(2) * t932 - mrSges(2,1) * t955 - qJDD(1) * mrSges(2,2) + t971;
t797 = t950 * t808 + t945 * t809;
t964 = m(6) * t831 - t837 * mrSges(6,1) + t838 * mrSges(6,2) - t881 * t872 + t882 * t873 + t797;
t961 = m(5) * t860 - t868 * mrSges(5,1) + t869 * mrSges(5,2) - t897 * t890 + t898 * t891 + t964;
t792 = m(4) * t889 - t903 * mrSges(4,1) + t904 * mrSges(4,2) - t914 * t907 + t915 * t908 + t961;
t919 = -t955 * pkin(7) + t965;
t957 = -m(3) * t919 + t926 * mrSges(3,1) - t925 * mrSges(3,2) - t929 * t975 + t930 * t974 - t792;
t790 = m(2) * t931 + qJDD(1) * mrSges(2,1) - t955 * mrSges(2,2) + t957;
t976 = t949 * t763 + t954 * t790;
t765 = t953 * t770 + t948 * t771;
t972 = t954 * t763 - t790 * t949;
t852 = Ifges(6,5) * t882 + Ifges(6,6) * t881 + Ifges(6,3) * t937;
t780 = mrSges(6,2) * t831 - mrSges(6,3) * t817 + Ifges(6,1) * t838 + Ifges(6,4) * t837 + Ifges(6,5) * t936 - pkin(10) * t797 - t801 * t945 + t802 * t950 + t852 * t881 - t853 * t937;
t959 = mrSges(7,1) * t812 - mrSges(7,2) * t813 + Ifges(7,5) * t829 + Ifges(7,6) * t828 + Ifges(7,3) * t836 + t840 * t871 - t841 * t870;
t781 = -mrSges(6,1) * t831 + mrSges(6,3) * t818 + Ifges(6,4) * t838 + Ifges(6,2) * t837 + Ifges(6,6) * t936 - pkin(5) * t797 - t852 * t882 + t854 * t937 - t959;
t874 = Ifges(5,5) * t898 + Ifges(5,6) * t897 + Ifges(5,3) * t940;
t766 = -mrSges(5,1) * t860 + mrSges(5,3) * t826 + Ifges(5,4) * t869 + Ifges(5,2) * t868 + Ifges(5,6) * t939 - pkin(4) * t964 + pkin(9) * t968 + t946 * t780 + t951 * t781 - t898 * t874 + t940 * t876;
t773 = mrSges(5,2) * t860 - mrSges(5,3) * t825 + Ifges(5,1) * t869 + Ifges(5,4) * t868 + Ifges(5,5) * t939 - pkin(9) * t788 + t780 * t951 - t781 * t946 + t874 * t897 - t875 * t940;
t894 = Ifges(4,5) * t915 + Ifges(4,6) * t914 + Ifges(4,3) * qJD(2);
t759 = -mrSges(4,1) * t889 + mrSges(4,3) * t862 + Ifges(4,4) * t904 + Ifges(4,2) * t903 + Ifges(4,6) * qJDD(2) - pkin(3) * t961 + pkin(8) * t969 + qJD(2) * t896 + t952 * t766 + t947 * t773 - t915 * t894;
t760 = mrSges(4,2) * t889 - mrSges(4,3) * t861 + Ifges(4,1) * t904 + Ifges(4,4) * t903 + Ifges(4,5) * qJDD(2) - pkin(8) * t779 - qJD(2) * t895 - t766 * t947 + t773 * t952 + t894 * t914;
t916 = Ifges(3,3) * qJD(2) + (Ifges(3,5) * t948 + Ifges(3,6) * t953) * qJD(1);
t755 = -mrSges(3,1) * t919 + mrSges(3,3) * t906 + Ifges(3,4) * t925 + Ifges(3,2) * t926 + Ifges(3,6) * qJDD(2) - pkin(2) * t792 + qJ(3) * t970 + qJD(2) * t918 + t944 * t759 + t943 * t760 - t916 * t975;
t757 = mrSges(3,2) * t919 - mrSges(3,3) * t905 + Ifges(3,1) * t925 + Ifges(3,4) * t926 + Ifges(3,5) * qJDD(2) - qJ(3) * t772 - qJD(2) * t917 - t759 * t943 + t760 * t944 + t916 * t974;
t963 = mrSges(2,1) * t931 - mrSges(2,2) * t932 + Ifges(2,3) * qJDD(1) + pkin(1) * t957 + pkin(7) * t971 + t953 * t755 + t948 * t757;
t758 = mrSges(2,1) * g(3) + mrSges(2,3) * t932 + t955 * Ifges(2,5) + Ifges(2,6) * qJDD(1) - pkin(1) * t765 - t980;
t753 = -mrSges(2,2) * g(3) - mrSges(2,3) * t931 + Ifges(2,5) * qJDD(1) - Ifges(2,6) * t955 - pkin(7) * t765 - t755 * t948 + t757 * t953;
t1 = [-m(1) * g(1) + t972; -m(1) * g(2) + t976; (-m(1) - m(2)) * g(3) + t765; -mrSges(1,2) * g(3) + mrSges(1,3) * g(2) - pkin(6) * t976 + t954 * t753 - t949 * t758; mrSges(1,1) * g(3) - mrSges(1,3) * g(1) + pkin(6) * t972 + t949 * t753 + t954 * t758; -mrSges(1,1) * g(2) + mrSges(1,2) * g(1) + t963; t963; t980; t792; -t958; -t960; t959;];
tauJB  = t1;
