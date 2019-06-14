% Calculate vector of inverse dynamics joint torques and base forces with Newton-Euler
% S6PRRRRR4
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
% pkin [13x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,alpha2,alpha3,d2,d3,d4,d5,d6,theta1]';
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
% Datum: 2019-05-05 11:33
% Revision: 89c353f7eff3bd693eda4e29f35b2761dbc3ada0 (2019-05-03)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauJB = S6PRRRRR4_invdynJB_fixb_snew_vp2(qJ, qJD, qJDD, g, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(6,1),zeros(3,1),zeros(13,1),zeros(7,1),zeros(7,3),zeros(7,6)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6PRRRRR4_invdynJB_fixb_snew_vp2: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6PRRRRR4_invdynJB_fixb_snew_vp2: qJD has to be [6x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [6 1]), ...
  'S6PRRRRR4_invdynJB_fixb_snew_vp2: qJDD has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6PRRRRR4_invdynJB_fixb_snew_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [13 1]), ...
  'S6PRRRRR4_invdynJB_fixb_snew_vp2: pkin has to be [13x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6PRRRRR4_invdynJB_fixb_snew_vp2: m has to be [7x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [7,3]), ...
  'S6PRRRRR4_invdynJB_fixb_snew_vp2: mrSges has to be [7x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [7 6]), ...
  'S6PRRRRR4_invdynJB_fixb_snew_vp2: Ifges has to be [7x6] (double)'); 

%% Symbolic Calculation
% From invdyn_fixb_NewtonEuler_linkframe_tauJB_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-05-05 11:26:26
% EndTime: 2019-05-05 11:27:17
% DurationCPUTime: 50.40s
% Computational Cost: add. (890819->354), mult. (1839535->467), div. (0->0), fcn. (1498359->16), ass. (0->160)
t918 = sin(pkin(13));
t921 = cos(pkin(13));
t908 = g(1) * t918 - g(2) * t921;
t909 = -g(1) * t921 - g(2) * t918;
t917 = -g(3) + qJDD(1);
t928 = sin(qJ(2));
t923 = cos(pkin(6));
t933 = cos(qJ(2));
t955 = t923 * t933;
t920 = sin(pkin(6));
t958 = t920 * t933;
t878 = t908 * t955 - t909 * t928 + t917 * t958;
t919 = sin(pkin(7));
t934 = qJD(2) ^ 2;
t875 = pkin(9) * t919 * t934 + qJDD(2) * pkin(2) + t878;
t956 = t923 * t928;
t959 = t920 * t928;
t879 = t908 * t956 + t933 * t909 + t917 * t959;
t951 = qJDD(2) * t919;
t876 = -pkin(2) * t934 + pkin(9) * t951 + t879;
t893 = -t908 * t920 + t917 * t923;
t922 = cos(pkin(7));
t927 = sin(qJ(3));
t932 = cos(qJ(3));
t846 = -t927 * t876 + (t875 * t922 + t893 * t919) * t932;
t915 = qJD(2) * t922 + qJD(3);
t952 = qJD(2) * t932;
t949 = t919 * t952;
t897 = -mrSges(4,2) * t915 + mrSges(4,3) * t949;
t953 = qJD(2) * t919;
t898 = (-mrSges(4,1) * t932 + mrSges(4,2) * t927) * t953;
t900 = (qJD(3) * t952 + qJDD(2) * t927) * t919;
t914 = qJDD(2) * t922 + qJDD(3);
t899 = (-pkin(3) * t932 - pkin(10) * t927) * t953;
t913 = t915 ^ 2;
t950 = t927 * t953;
t833 = -t914 * pkin(3) - t913 * pkin(10) + t899 * t950 - t846;
t926 = sin(qJ(4));
t931 = cos(qJ(4));
t892 = t915 * t926 + t931 * t950;
t867 = -qJD(4) * t892 - t900 * t926 + t914 * t931;
t891 = t915 * t931 - t926 * t950;
t868 = qJD(4) * t891 + t900 * t931 + t914 * t926;
t907 = qJD(4) - t949;
t880 = -mrSges(5,2) * t907 + mrSges(5,3) * t891;
t881 = mrSges(5,1) * t907 - mrSges(5,3) * t892;
t957 = t922 * t927;
t960 = t919 * t927;
t847 = t875 * t957 + t932 * t876 + t893 * t960;
t834 = -pkin(3) * t913 + pkin(10) * t914 + t899 * t949 + t847;
t888 = t922 * t893;
t901 = -qJD(3) * t950 + t932 * t951;
t844 = -t901 * pkin(3) - t900 * pkin(10) + t888 + (-t875 + (pkin(3) * t927 - pkin(10) * t932) * t915 * qJD(2)) * t919;
t822 = -t926 * t834 + t931 * t844;
t895 = qJDD(4) - t901;
t819 = (t891 * t907 - t868) * pkin(11) + (t891 * t892 + t895) * pkin(4) + t822;
t823 = t931 * t834 + t926 * t844;
t882 = pkin(4) * t907 - pkin(11) * t892;
t890 = t891 ^ 2;
t821 = -pkin(4) * t890 + pkin(11) * t867 - t882 * t907 + t823;
t925 = sin(qJ(5));
t930 = cos(qJ(5));
t816 = t925 * t819 + t930 * t821;
t873 = t891 * t930 - t892 * t925;
t874 = t891 * t925 + t892 * t930;
t855 = -pkin(5) * t873 - pkin(12) * t874;
t894 = qJDD(5) + t895;
t906 = qJD(5) + t907;
t905 = t906 ^ 2;
t813 = -pkin(5) * t905 + pkin(12) * t894 + t855 * t873 + t816;
t824 = -t867 * pkin(4) - t890 * pkin(11) + t892 * t882 + t833;
t839 = -qJD(5) * t874 + t867 * t930 - t868 * t925;
t840 = qJD(5) * t873 + t867 * t925 + t868 * t930;
t817 = (-t873 * t906 - t840) * pkin(12) + (t874 * t906 - t839) * pkin(5) + t824;
t924 = sin(qJ(6));
t929 = cos(qJ(6));
t810 = -t813 * t924 + t817 * t929;
t857 = -t874 * t924 + t906 * t929;
t827 = qJD(6) * t857 + t840 * t929 + t894 * t924;
t838 = qJDD(6) - t839;
t858 = t874 * t929 + t906 * t924;
t845 = -mrSges(7,1) * t857 + mrSges(7,2) * t858;
t870 = qJD(6) - t873;
t848 = -mrSges(7,2) * t870 + mrSges(7,3) * t857;
t806 = m(7) * t810 + mrSges(7,1) * t838 - mrSges(7,3) * t827 - t845 * t858 + t848 * t870;
t811 = t813 * t929 + t817 * t924;
t826 = -qJD(6) * t858 - t840 * t924 + t894 * t929;
t849 = mrSges(7,1) * t870 - mrSges(7,3) * t858;
t807 = m(7) * t811 - mrSges(7,2) * t838 + mrSges(7,3) * t826 + t845 * t857 - t849 * t870;
t795 = t929 * t806 + t924 * t807;
t859 = -mrSges(6,2) * t906 + mrSges(6,3) * t873;
t860 = mrSges(6,1) * t906 - mrSges(6,3) * t874;
t939 = m(6) * t824 - t839 * mrSges(6,1) + mrSges(6,2) * t840 - t873 * t859 + t860 * t874 + t795;
t936 = -m(5) * t833 + t867 * mrSges(5,1) - mrSges(5,2) * t868 + t891 * t880 - t881 * t892 - t939;
t790 = m(4) * t846 + mrSges(4,1) * t914 - mrSges(4,3) * t900 + t897 * t915 - t898 * t950 + t936;
t961 = t790 * t932;
t896 = mrSges(4,1) * t915 - mrSges(4,3) * t950;
t854 = -mrSges(6,1) * t873 + mrSges(6,2) * t874;
t945 = -t806 * t924 + t929 * t807;
t793 = m(6) * t816 - mrSges(6,2) * t894 + mrSges(6,3) * t839 + t854 * t873 - t860 * t906 + t945;
t815 = t819 * t930 - t821 * t925;
t812 = -pkin(5) * t894 - pkin(12) * t905 + t855 * t874 - t815;
t940 = -m(7) * t812 + t826 * mrSges(7,1) - mrSges(7,2) * t827 + t857 * t848 - t849 * t858;
t802 = m(6) * t815 + mrSges(6,1) * t894 - mrSges(6,3) * t840 - t854 * t874 + t859 * t906 + t940;
t787 = t925 * t793 + t930 * t802;
t877 = -mrSges(5,1) * t891 + mrSges(5,2) * t892;
t785 = m(5) * t822 + mrSges(5,1) * t895 - mrSges(5,3) * t868 - t877 * t892 + t880 * t907 + t787;
t946 = t930 * t793 - t802 * t925;
t786 = m(5) * t823 - mrSges(5,2) * t895 + mrSges(5,3) * t867 + t877 * t891 - t881 * t907 + t946;
t947 = -t785 * t926 + t931 * t786;
t776 = m(4) * t847 - mrSges(4,2) * t914 + mrSges(4,3) * t901 - t896 * t915 + t898 * t949 + t947;
t779 = t931 * t785 + t926 * t786;
t856 = -t919 * t875 + t888;
t778 = m(4) * t856 - t901 * mrSges(4,1) + t900 * mrSges(4,2) + (t896 * t927 - t897 * t932) * t953 + t779;
t765 = t776 * t957 - t778 * t919 + t922 * t961;
t761 = m(3) * t878 + qJDD(2) * mrSges(3,1) - mrSges(3,2) * t934 + t765;
t764 = t776 * t960 + t922 * t778 + t919 * t961;
t763 = m(3) * t893 + t764;
t772 = t932 * t776 - t790 * t927;
t771 = m(3) * t879 - mrSges(3,1) * t934 - qJDD(2) * mrSges(3,2) + t772;
t751 = t761 * t955 - t763 * t920 + t771 * t956;
t749 = m(2) * t908 + t751;
t757 = -t761 * t928 + t933 * t771;
t756 = m(2) * t909 + t757;
t954 = t921 * t749 + t918 * t756;
t750 = t761 * t958 + t923 * t763 + t771 * t959;
t948 = -t749 * t918 + t921 * t756;
t944 = m(2) * t917 + t750;
t828 = Ifges(7,5) * t858 + Ifges(7,6) * t857 + Ifges(7,3) * t870;
t830 = Ifges(7,1) * t858 + Ifges(7,4) * t857 + Ifges(7,5) * t870;
t799 = -mrSges(7,1) * t812 + mrSges(7,3) * t811 + Ifges(7,4) * t827 + Ifges(7,2) * t826 + Ifges(7,6) * t838 - t828 * t858 + t830 * t870;
t829 = Ifges(7,4) * t858 + Ifges(7,2) * t857 + Ifges(7,6) * t870;
t800 = mrSges(7,2) * t812 - mrSges(7,3) * t810 + Ifges(7,1) * t827 + Ifges(7,4) * t826 + Ifges(7,5) * t838 + t828 * t857 - t829 * t870;
t850 = Ifges(6,5) * t874 + Ifges(6,6) * t873 + Ifges(6,3) * t906;
t851 = Ifges(6,4) * t874 + Ifges(6,2) * t873 + Ifges(6,6) * t906;
t780 = mrSges(6,2) * t824 - mrSges(6,3) * t815 + Ifges(6,1) * t840 + Ifges(6,4) * t839 + Ifges(6,5) * t894 - pkin(12) * t795 - t799 * t924 + t800 * t929 + t850 * t873 - t851 * t906;
t852 = Ifges(6,1) * t874 + Ifges(6,4) * t873 + Ifges(6,5) * t906;
t937 = mrSges(7,1) * t810 - mrSges(7,2) * t811 + Ifges(7,5) * t827 + Ifges(7,6) * t826 + Ifges(7,3) * t838 + t829 * t858 - t830 * t857;
t781 = -mrSges(6,1) * t824 + mrSges(6,3) * t816 + Ifges(6,4) * t840 + Ifges(6,2) * t839 + Ifges(6,6) * t894 - pkin(5) * t795 - t850 * t874 + t852 * t906 - t937;
t861 = Ifges(5,5) * t892 + Ifges(5,6) * t891 + Ifges(5,3) * t907;
t863 = Ifges(5,1) * t892 + Ifges(5,4) * t891 + Ifges(5,5) * t907;
t766 = -mrSges(5,1) * t833 + mrSges(5,3) * t823 + Ifges(5,4) * t868 + Ifges(5,2) * t867 + Ifges(5,6) * t895 - pkin(4) * t939 + pkin(11) * t946 + t925 * t780 + t930 * t781 - t892 * t861 + t907 * t863;
t862 = Ifges(5,4) * t892 + Ifges(5,2) * t891 + Ifges(5,6) * t907;
t767 = mrSges(5,2) * t833 - mrSges(5,3) * t822 + Ifges(5,1) * t868 + Ifges(5,4) * t867 + Ifges(5,5) * t895 - pkin(11) * t787 + t780 * t930 - t781 * t925 + t861 * t891 - t862 * t907;
t885 = Ifges(4,6) * t915 + (Ifges(4,4) * t927 + Ifges(4,2) * t932) * t953;
t886 = Ifges(4,5) * t915 + (Ifges(4,1) * t927 + Ifges(4,4) * t932) * t953;
t752 = Ifges(4,5) * t900 + Ifges(4,6) * t901 + Ifges(4,3) * t914 + mrSges(4,1) * t846 - mrSges(4,2) * t847 + t926 * t767 + t931 * t766 + pkin(3) * t936 + pkin(10) * t947 + (t885 * t927 - t886 * t932) * t953;
t884 = Ifges(4,3) * t915 + (Ifges(4,5) * t927 + Ifges(4,6) * t932) * t953;
t753 = mrSges(4,2) * t856 - mrSges(4,3) * t846 + Ifges(4,1) * t900 + Ifges(4,4) * t901 + Ifges(4,5) * t914 - pkin(10) * t779 - t766 * t926 + t767 * t931 + t884 * t949 - t885 * t915;
t938 = -mrSges(6,1) * t815 + mrSges(6,2) * t816 - Ifges(6,5) * t840 - Ifges(6,6) * t839 - Ifges(6,3) * t894 - pkin(5) * t940 - pkin(12) * t945 - t929 * t799 - t924 * t800 - t874 * t851 + t873 * t852;
t935 = mrSges(5,1) * t822 - mrSges(5,2) * t823 + Ifges(5,5) * t868 + Ifges(5,6) * t867 + Ifges(5,3) * t895 + pkin(4) * t787 + t892 * t862 - t891 * t863 - t938;
t758 = -mrSges(4,1) * t856 + mrSges(4,3) * t847 + Ifges(4,4) * t900 + Ifges(4,2) * t901 + Ifges(4,6) * t914 - pkin(3) * t779 - t884 * t950 + t915 * t886 - t935;
t941 = pkin(9) * t772 + t753 * t927 + t758 * t932;
t746 = -mrSges(3,1) * t893 + mrSges(3,3) * t879 + t934 * Ifges(3,5) + Ifges(3,6) * qJDD(2) - pkin(2) * t764 - t919 * t752 + t941 * t922;
t747 = mrSges(3,2) * t893 - mrSges(3,3) * t878 + Ifges(3,5) * qJDD(2) - t934 * Ifges(3,6) + t932 * t753 - t927 * t758 + (-t764 * t919 - t765 * t922) * pkin(9);
t942 = pkin(8) * t757 + t746 * t933 + t747 * t928;
t745 = mrSges(3,1) * t878 - mrSges(3,2) * t879 + Ifges(3,3) * qJDD(2) + pkin(2) * t765 + t922 * t752 + t941 * t919;
t744 = mrSges(2,2) * t917 - mrSges(2,3) * t908 - t928 * t746 + t933 * t747 + (-t750 * t920 - t751 * t923) * pkin(8);
t743 = -mrSges(2,1) * t917 + mrSges(2,3) * t909 - pkin(1) * t750 - t920 * t745 + t942 * t923;
t1 = [-m(1) * g(1) + t948; -m(1) * g(2) + t954; -m(1) * g(3) + t944; -mrSges(1,2) * g(3) + mrSges(1,3) * g(2) - qJ(1) * t954 - t918 * t743 + t921 * t744; mrSges(1,1) * g(3) - mrSges(1,3) * g(1) + qJ(1) * t948 + t921 * t743 + t918 * t744; -mrSges(1,1) * g(2) + mrSges(2,1) * t908 + mrSges(1,2) * g(1) - mrSges(2,2) * t909 + pkin(1) * t751 + t923 * t745 + t942 * t920; t944; t745; t752; t935; -t938; t937;];
tauJB  = t1;
