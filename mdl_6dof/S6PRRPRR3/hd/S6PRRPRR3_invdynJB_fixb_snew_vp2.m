% Calculate vector of inverse dynamics joint torques and base forces with Newton-Euler
% S6PRRPRR3
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
%   pkin=[a2,a3,a4,a5,a6,alpha2,alpha3,d2,d3,d5,d6,theta1,theta4]';
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
% Datum: 2019-05-05 04:54
% Revision: 89c353f7eff3bd693eda4e29f35b2761dbc3ada0 (2019-05-03)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauJB = S6PRRPRR3_invdynJB_fixb_snew_vp2(qJ, qJD, qJDD, g, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(6,1),zeros(3,1),zeros(13,1),zeros(7,1),zeros(7,3),zeros(7,6)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6PRRPRR3_invdynJB_fixb_snew_vp2: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6PRRPRR3_invdynJB_fixb_snew_vp2: qJD has to be [6x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [6 1]), ...
  'S6PRRPRR3_invdynJB_fixb_snew_vp2: qJDD has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6PRRPRR3_invdynJB_fixb_snew_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [13 1]), ...
  'S6PRRPRR3_invdynJB_fixb_snew_vp2: pkin has to be [13x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6PRRPRR3_invdynJB_fixb_snew_vp2: m has to be [7x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [7,3]), ...
  'S6PRRPRR3_invdynJB_fixb_snew_vp2: mrSges has to be [7x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [7 6]), ...
  'S6PRRPRR3_invdynJB_fixb_snew_vp2: Ifges has to be [7x6] (double)'); 

%% Symbolic Calculation
% From invdyn_fixb_NewtonEuler_linkframe_tauJB_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-05-05 04:47:08
% EndTime: 2019-05-05 04:47:52
% DurationCPUTime: 42.80s
% Computational Cost: add. (697838->351), mult. (1695088->467), div. (0->0), fcn. (1351240->16), ass. (0->159)
t963 = -2 * qJD(4);
t918 = sin(pkin(13));
t922 = cos(pkin(13));
t928 = sin(qJ(3));
t932 = cos(qJ(3));
t920 = sin(pkin(7));
t951 = qJD(2) * t920;
t895 = (t918 * t928 - t922 * t932) * t951;
t919 = sin(pkin(12));
t923 = cos(pkin(12));
t907 = t919 * g(1) - t923 * g(2);
t908 = -t923 * g(1) - t919 * g(2);
t917 = -g(3) + qJDD(1);
t929 = sin(qJ(2));
t925 = cos(pkin(6));
t933 = cos(qJ(2));
t953 = t925 * t933;
t921 = sin(pkin(6));
t957 = t921 * t933;
t877 = t907 * t953 - t929 * t908 + t917 * t957;
t934 = qJD(2) ^ 2;
t962 = pkin(9) * t920;
t869 = qJDD(2) * pkin(2) + t934 * t962 + t877;
t954 = t925 * t929;
t958 = t921 * t929;
t878 = t907 * t954 + t933 * t908 + t917 * t958;
t870 = -t934 * pkin(2) + qJDD(2) * t962 + t878;
t897 = -t921 * t907 + t925 * t917;
t924 = cos(pkin(7));
t955 = t924 * t932;
t959 = t920 * t932;
t841 = t869 * t955 - t928 * t870 + t897 * t959;
t949 = qJD(2) * qJD(3);
t902 = (qJDD(2) * t928 + t932 * t949) * t920;
t912 = t924 * qJDD(2) + qJDD(3);
t913 = t924 * qJD(2) + qJD(3);
t946 = t932 * t951;
t961 = t920 ^ 2 * t934;
t833 = (t913 * t946 - t902) * qJ(4) + (t928 * t932 * t961 + t912) * pkin(3) + t841;
t956 = t924 * t928;
t960 = t920 * t928;
t842 = t869 * t956 + t932 * t870 + t897 * t960;
t947 = t928 * t951;
t898 = t913 * pkin(3) - qJ(4) * t947;
t903 = (qJDD(2) * t932 - t928 * t949) * t920;
t948 = t932 ^ 2 * t961;
t834 = -pkin(3) * t948 + t903 * qJ(4) - t913 * t898 + t842;
t896 = (t918 * t932 + t922 * t928) * t951;
t823 = t922 * t833 - t918 * t834 + t896 * t963;
t824 = t918 * t833 + t922 * t834 + t895 * t963;
t871 = t895 * mrSges(5,1) + t896 * mrSges(5,2);
t875 = -t918 * t902 + t922 * t903;
t883 = t913 * mrSges(5,1) - t896 * mrSges(5,3);
t872 = t895 * pkin(4) - t896 * pkin(10);
t911 = t913 ^ 2;
t822 = -t911 * pkin(4) + t912 * pkin(10) - t895 * t872 + t824;
t854 = -t920 * t869 + t924 * t897;
t843 = -t903 * pkin(3) - qJ(4) * t948 + t898 * t947 + qJDD(4) + t854;
t876 = t922 * t902 + t918 * t903;
t826 = (t895 * t913 - t876) * pkin(10) + (t896 * t913 - t875) * pkin(4) + t843;
t927 = sin(qJ(5));
t931 = cos(qJ(5));
t819 = t931 * t822 + t927 * t826;
t880 = -t927 * t896 + t931 * t913;
t881 = t931 * t896 + t927 * t913;
t856 = -t880 * pkin(5) - t881 * pkin(11);
t874 = qJDD(5) - t875;
t894 = qJD(5) + t895;
t893 = t894 ^ 2;
t816 = -t893 * pkin(5) + t874 * pkin(11) + t880 * t856 + t819;
t821 = -t912 * pkin(4) - t911 * pkin(10) + t896 * t872 - t823;
t852 = -t881 * qJD(5) - t927 * t876 + t931 * t912;
t853 = t880 * qJD(5) + t931 * t876 + t927 * t912;
t817 = (-t880 * t894 - t853) * pkin(11) + (t881 * t894 - t852) * pkin(5) + t821;
t926 = sin(qJ(6));
t930 = cos(qJ(6));
t812 = -t926 * t816 + t930 * t817;
t858 = -t926 * t881 + t930 * t894;
t829 = t858 * qJD(6) + t930 * t853 + t926 * t874;
t859 = t930 * t881 + t926 * t894;
t839 = -t858 * mrSges(7,1) + t859 * mrSges(7,2);
t879 = qJD(6) - t880;
t844 = -t879 * mrSges(7,2) + t858 * mrSges(7,3);
t851 = qJDD(6) - t852;
t810 = m(7) * t812 + t851 * mrSges(7,1) - t829 * mrSges(7,3) - t859 * t839 + t879 * t844;
t813 = t930 * t816 + t926 * t817;
t828 = -t859 * qJD(6) - t926 * t853 + t930 * t874;
t845 = t879 * mrSges(7,1) - t859 * mrSges(7,3);
t811 = m(7) * t813 - t851 * mrSges(7,2) + t828 * mrSges(7,3) + t858 * t839 - t879 * t845;
t804 = -t926 * t810 + t930 * t811;
t855 = -t880 * mrSges(6,1) + t881 * mrSges(6,2);
t861 = t894 * mrSges(6,1) - t881 * mrSges(6,3);
t802 = m(6) * t819 - t874 * mrSges(6,2) + t852 * mrSges(6,3) + t880 * t855 - t894 * t861 + t804;
t818 = -t927 * t822 + t931 * t826;
t815 = -t874 * pkin(5) - t893 * pkin(11) + t881 * t856 - t818;
t814 = -m(7) * t815 + t828 * mrSges(7,1) - t829 * mrSges(7,2) + t858 * t844 - t859 * t845;
t860 = -t894 * mrSges(6,2) + t880 * mrSges(6,3);
t808 = m(6) * t818 + t874 * mrSges(6,1) - t853 * mrSges(6,3) - t881 * t855 + t894 * t860 + t814;
t943 = t931 * t802 - t927 * t808;
t793 = m(5) * t824 - t912 * mrSges(5,2) + t875 * mrSges(5,3) - t895 * t871 - t913 * t883 + t943;
t882 = -t913 * mrSges(5,2) - t895 * mrSges(5,3);
t803 = t930 * t810 + t926 * t811;
t937 = -m(6) * t821 + t852 * mrSges(6,1) - t853 * mrSges(6,2) + t880 * t860 - t881 * t861 - t803;
t799 = m(5) * t823 + t912 * mrSges(5,1) - t876 * mrSges(5,3) - t896 * t871 + t913 * t882 + t937;
t788 = t918 * t793 + t922 * t799;
t900 = -t913 * mrSges(4,2) + mrSges(4,3) * t946;
t901 = (-mrSges(4,1) * t932 + mrSges(4,2) * t928) * t951;
t786 = m(4) * t841 + t912 * mrSges(4,1) - t902 * mrSges(4,3) + t913 * t900 - t901 * t947 + t788;
t899 = t913 * mrSges(4,1) - mrSges(4,3) * t947;
t944 = t922 * t793 - t918 * t799;
t787 = m(4) * t842 - t912 * mrSges(4,2) + t903 * mrSges(4,3) - t913 * t899 + t901 * t946 + t944;
t797 = t927 * t802 + t931 * t808;
t796 = m(5) * t843 - t875 * mrSges(5,1) + t876 * mrSges(5,2) + t895 * t882 + t896 * t883 + t797;
t795 = m(4) * t854 - t903 * mrSges(4,1) + t902 * mrSges(4,2) + (t899 * t928 - t900 * t932) * t951 + t796;
t773 = t786 * t955 + t787 * t956 - t920 * t795;
t769 = m(3) * t877 + qJDD(2) * mrSges(3,1) - t934 * mrSges(3,2) + t773;
t772 = t786 * t959 + t787 * t960 + t924 * t795;
t771 = m(3) * t897 + t772;
t779 = -t928 * t786 + t932 * t787;
t778 = m(3) * t878 - t934 * mrSges(3,1) - qJDD(2) * mrSges(3,2) + t779;
t759 = t769 * t953 - t921 * t771 + t778 * t954;
t757 = m(2) * t907 + t759;
t765 = -t929 * t769 + t933 * t778;
t764 = m(2) * t908 + t765;
t952 = t923 * t757 + t919 * t764;
t758 = t769 * t957 + t925 * t771 + t778 * t958;
t945 = -t919 * t757 + t923 * t764;
t941 = m(2) * t917 + t758;
t835 = Ifges(7,5) * t859 + Ifges(7,6) * t858 + Ifges(7,3) * t879;
t837 = Ifges(7,1) * t859 + Ifges(7,4) * t858 + Ifges(7,5) * t879;
t805 = -mrSges(7,1) * t815 + mrSges(7,3) * t813 + Ifges(7,4) * t829 + Ifges(7,2) * t828 + Ifges(7,6) * t851 - t859 * t835 + t879 * t837;
t836 = Ifges(7,4) * t859 + Ifges(7,2) * t858 + Ifges(7,6) * t879;
t806 = mrSges(7,2) * t815 - mrSges(7,3) * t812 + Ifges(7,1) * t829 + Ifges(7,4) * t828 + Ifges(7,5) * t851 + t858 * t835 - t879 * t836;
t846 = Ifges(6,5) * t881 + Ifges(6,6) * t880 + Ifges(6,3) * t894;
t847 = Ifges(6,4) * t881 + Ifges(6,2) * t880 + Ifges(6,6) * t894;
t789 = mrSges(6,2) * t821 - mrSges(6,3) * t818 + Ifges(6,1) * t853 + Ifges(6,4) * t852 + Ifges(6,5) * t874 - pkin(11) * t803 - t926 * t805 + t930 * t806 + t880 * t846 - t894 * t847;
t848 = Ifges(6,1) * t881 + Ifges(6,4) * t880 + Ifges(6,5) * t894;
t936 = mrSges(7,1) * t812 - mrSges(7,2) * t813 + Ifges(7,5) * t829 + Ifges(7,6) * t828 + Ifges(7,3) * t851 + t859 * t836 - t858 * t837;
t790 = -mrSges(6,1) * t821 + mrSges(6,3) * t819 + Ifges(6,4) * t853 + Ifges(6,2) * t852 + Ifges(6,6) * t874 - pkin(5) * t803 - t881 * t846 + t894 * t848 - t936;
t863 = Ifges(5,4) * t896 - Ifges(5,2) * t895 + Ifges(5,6) * t913;
t864 = Ifges(5,1) * t896 - Ifges(5,4) * t895 + Ifges(5,5) * t913;
t887 = Ifges(4,6) * t913 + (Ifges(4,4) * t928 + Ifges(4,2) * t932) * t951;
t888 = Ifges(4,5) * t913 + (Ifges(4,1) * t928 + Ifges(4,4) * t932) * t951;
t766 = Ifges(4,5) * t902 + Ifges(4,6) * t903 + mrSges(4,1) * t841 - mrSges(4,2) * t842 + Ifges(5,5) * t876 + Ifges(5,6) * t875 + t896 * t863 + t895 * t864 + mrSges(5,1) * t823 - mrSges(5,2) * t824 + t927 * t789 + t931 * t790 + pkin(4) * t937 + pkin(10) * t943 + pkin(3) * t788 + (Ifges(4,3) + Ifges(5,3)) * t912 + (t887 * t928 - t888 * t932) * t951;
t862 = Ifges(5,5) * t896 - Ifges(5,6) * t895 + Ifges(5,3) * t913;
t774 = mrSges(5,2) * t843 - mrSges(5,3) * t823 + Ifges(5,1) * t876 + Ifges(5,4) * t875 + Ifges(5,5) * t912 - pkin(10) * t797 + t931 * t789 - t927 * t790 - t895 * t862 - t913 * t863;
t935 = mrSges(6,1) * t818 - mrSges(6,2) * t819 + Ifges(6,5) * t853 + Ifges(6,6) * t852 + Ifges(6,3) * t874 + pkin(5) * t814 + pkin(11) * t804 + t930 * t805 + t926 * t806 + t881 * t847 - t880 * t848;
t780 = -mrSges(5,1) * t843 + mrSges(5,3) * t824 + Ifges(5,4) * t876 + Ifges(5,2) * t875 + Ifges(5,6) * t912 - pkin(4) * t797 - t896 * t862 + t913 * t864 - t935;
t886 = Ifges(4,3) * t913 + (Ifges(4,5) * t928 + Ifges(4,6) * t932) * t951;
t760 = -mrSges(4,1) * t854 + mrSges(4,3) * t842 + Ifges(4,4) * t902 + Ifges(4,2) * t903 + Ifges(4,6) * t912 - pkin(3) * t796 + qJ(4) * t944 + t918 * t774 + t922 * t780 - t886 * t947 + t913 * t888;
t761 = mrSges(4,2) * t854 - mrSges(4,3) * t841 + Ifges(4,1) * t902 + Ifges(4,4) * t903 + Ifges(4,5) * t912 - qJ(4) * t788 + t922 * t774 - t918 * t780 + t886 * t946 - t913 * t887;
t938 = pkin(9) * t779 + t760 * t932 + t761 * t928;
t754 = -mrSges(3,1) * t897 + mrSges(3,3) * t878 + t934 * Ifges(3,5) + Ifges(3,6) * qJDD(2) - pkin(2) * t772 - t920 * t766 + t938 * t924;
t755 = mrSges(3,2) * t897 - mrSges(3,3) * t877 + Ifges(3,5) * qJDD(2) - t934 * Ifges(3,6) - t928 * t760 + t932 * t761 + (-t772 * t920 - t773 * t924) * pkin(9);
t939 = pkin(8) * t765 + t754 * t933 + t755 * t929;
t753 = mrSges(3,1) * t877 - mrSges(3,2) * t878 + Ifges(3,3) * qJDD(2) + pkin(2) * t773 + t924 * t766 + t938 * t920;
t752 = mrSges(2,2) * t917 - mrSges(2,3) * t907 - t929 * t754 + t933 * t755 + (-t758 * t921 - t759 * t925) * pkin(8);
t751 = -mrSges(2,1) * t917 + mrSges(2,3) * t908 - pkin(1) * t758 - t921 * t753 + t939 * t925;
t1 = [-m(1) * g(1) + t945; -m(1) * g(2) + t952; -m(1) * g(3) + t941; -mrSges(1,2) * g(3) + mrSges(1,3) * g(2) - qJ(1) * t952 - t919 * t751 + t923 * t752; mrSges(1,1) * g(3) - mrSges(1,3) * g(1) + qJ(1) * t945 + t923 * t751 + t919 * t752; -mrSges(1,1) * g(2) + mrSges(2,1) * t907 + mrSges(1,2) * g(1) - mrSges(2,2) * t908 + pkin(1) * t759 + t925 * t753 + t939 * t921; t941; t753; t766; t796; t935; t936;];
tauJB  = t1;
