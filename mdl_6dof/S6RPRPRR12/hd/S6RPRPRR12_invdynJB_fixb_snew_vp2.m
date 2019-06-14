% Calculate vector of inverse dynamics joint torques and base forces with Newton-Euler
% S6RPRPRR12
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
%   pkin=[a2,a3,a4,a5,a6,d1,d3,d5,d6]';
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
% Datum: 2019-05-05 20:47
% Revision: 89c353f7eff3bd693eda4e29f35b2761dbc3ada0 (2019-05-03)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauJB = S6RPRPRR12_invdynJB_fixb_snew_vp2(qJ, qJD, qJDD, g, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(6,1),zeros(3,1),zeros(9,1),zeros(7,1),zeros(7,3),zeros(7,6)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPRPRR12_invdynJB_fixb_snew_vp2: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RPRPRR12_invdynJB_fixb_snew_vp2: qJD has to be [6x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [6 1]), ...
  'S6RPRPRR12_invdynJB_fixb_snew_vp2: qJDD has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6RPRPRR12_invdynJB_fixb_snew_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S6RPRPRR12_invdynJB_fixb_snew_vp2: pkin has to be [9x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6RPRPRR12_invdynJB_fixb_snew_vp2: m has to be [7x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [7,3]), ...
  'S6RPRPRR12_invdynJB_fixb_snew_vp2: mrSges has to be [7x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [7 6]), ...
  'S6RPRPRR12_invdynJB_fixb_snew_vp2: Ifges has to be [7x6] (double)'); 

%% Symbolic Calculation
% From invdyn_fixb_NewtonEuler_linkframe_tauJB_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-05-05 20:45:09
% EndTime: 2019-05-05 20:45:16
% DurationCPUTime: 4.60s
% Computational Cost: add. (50833->329), mult. (101504->387), div. (0->0), fcn. (56432->8), ass. (0->138)
t935 = -2 * qJD(4);
t934 = Ifges(4,5) - Ifges(5,4);
t933 = Ifges(4,6) - Ifges(5,5);
t932 = (Ifges(4,3) + Ifges(5,1));
t881 = sin(qJ(1));
t885 = cos(qJ(1));
t856 = -t885 * g(1) - t881 * g(2);
t903 = qJDD(1) * qJ(2) + (2 * qJD(2) * qJD(1)) + t856;
t887 = qJD(1) ^ 2;
t929 = (-pkin(1) - pkin(7));
t914 = t929 * t887;
t812 = t914 + t903;
t880 = sin(qJ(3));
t884 = cos(qJ(3));
t917 = qJD(1) * qJD(3);
t913 = t884 * t917;
t845 = qJDD(1) * t880 + t913;
t862 = t880 * t917;
t846 = qJDD(1) * t884 - t862;
t918 = qJD(1) * t880;
t850 = -(qJD(3) * mrSges(4,2)) - mrSges(4,3) * t918;
t863 = t884 * qJD(1);
t851 = (qJD(3) * mrSges(4,1)) - mrSges(4,3) * t863;
t895 = pkin(3) * t913 + t863 * t935 + t903 + (-t846 + t862) * qJ(4);
t785 = t845 * pkin(3) + t895 + t914;
t852 = mrSges(5,1) * t918 - (qJD(3) * mrSges(5,3));
t853 = mrSges(5,1) * t863 + (qJD(3) * mrSges(5,2));
t854 = pkin(4) * t863 - (qJD(3) * pkin(8));
t875 = t880 ^ 2;
t928 = pkin(3) + pkin(8);
t773 = -t854 * t863 + t928 * t845 + (-pkin(4) * t875 + t929) * t887 + t895;
t842 = (pkin(3) * t880 - qJ(4) * t884) * qJD(1);
t886 = qJD(3) ^ 2;
t855 = t881 * g(1) - t885 * g(2);
t902 = -t887 * qJ(2) + qJDD(2) - t855;
t813 = qJDD(1) * t929 + t902;
t920 = t884 * t813;
t901 = -t886 * qJ(4) + t842 * t863 + qJDD(4) - t920;
t927 = pkin(8) * t887;
t778 = t846 * pkin(4) - t928 * qJDD(3) + (pkin(4) * t917 + t884 * t927 - g(3)) * t880 + t901;
t879 = sin(qJ(5));
t883 = cos(qJ(5));
t762 = -t879 * t773 + t883 * t778;
t840 = -qJD(3) * t879 + t883 * t918;
t800 = qJD(5) * t840 + qJDD(3) * t883 + t845 * t879;
t839 = qJDD(5) + t846;
t841 = qJD(3) * t883 + t879 * t918;
t860 = t863 + qJD(5);
t759 = (t840 * t860 - t800) * pkin(9) + (t840 * t841 + t839) * pkin(5) + t762;
t763 = t883 * t773 + t879 * t778;
t799 = -qJD(5) * t841 - qJDD(3) * t879 + t845 * t883;
t811 = pkin(5) * t860 - pkin(9) * t841;
t838 = t840 ^ 2;
t760 = -pkin(5) * t838 + t799 * pkin(9) - t811 * t860 + t763;
t878 = sin(qJ(6));
t882 = cos(qJ(6));
t757 = t759 * t882 - t760 * t878;
t801 = t840 * t882 - t841 * t878;
t771 = t801 * qJD(6) + t799 * t878 + t800 * t882;
t802 = t840 * t878 + t841 * t882;
t784 = -mrSges(7,1) * t801 + mrSges(7,2) * t802;
t857 = qJD(6) + t860;
t790 = -mrSges(7,2) * t857 + t801 * mrSges(7,3);
t828 = qJDD(6) + t839;
t754 = m(7) * t757 + t828 * mrSges(7,1) - t771 * mrSges(7,3) - t802 * t784 + t790 * t857;
t758 = t759 * t878 + t760 * t882;
t770 = -t802 * qJD(6) + t799 * t882 - t800 * t878;
t791 = mrSges(7,1) * t857 - t802 * mrSges(7,3);
t755 = m(7) * t758 - t828 * mrSges(7,2) + t770 * mrSges(7,3) + t801 * t784 - t791 * t857;
t745 = t882 * t754 + t878 * t755;
t804 = -mrSges(6,1) * t840 + mrSges(6,2) * t841;
t807 = -mrSges(6,2) * t860 + mrSges(6,3) * t840;
t742 = m(6) * t762 + mrSges(6,1) * t839 - t800 * mrSges(6,3) - t804 * t841 + t807 * t860 + t745;
t808 = mrSges(6,1) * t860 - mrSges(6,3) * t841;
t909 = -t754 * t878 + t882 * t755;
t743 = m(6) * t763 - mrSges(6,2) * t839 + t799 * mrSges(6,3) + t804 * t840 - t808 * t860 + t909;
t910 = -t879 * t742 + t883 * t743;
t894 = m(5) * t785 - t846 * mrSges(5,3) - (t852 * t880 + t853 * t884) * qJD(1) + t910;
t924 = mrSges(4,1) - mrSges(5,2);
t931 = -m(4) * t812 - t846 * mrSges(4,2) - t845 * t924 - t850 * t918 - t851 * t863 - t894;
t806 = -g(3) * t884 + t880 * t813;
t786 = pkin(3) * t886 - qJDD(3) * qJ(4) + (qJD(3) * t935) + t842 * t918 - t806;
t926 = t880 * g(3);
t925 = mrSges(2,1) - mrSges(3,2);
t923 = Ifges(2,5) - Ifges(3,4);
t922 = (-Ifges(2,6) + Ifges(3,5));
t921 = -Ifges(5,6) - Ifges(4,4);
t805 = t920 + t926;
t740 = t883 * t742 + t879 * t743;
t788 = -qJDD(3) * pkin(3) + t901 - t926;
t898 = -m(5) * t788 - t846 * mrSges(5,1) - t740;
t843 = (-mrSges(5,2) * t880 - mrSges(5,3) * t884) * qJD(1);
t907 = qJD(1) * (-t843 - (mrSges(4,1) * t880 + mrSges(4,2) * t884) * qJD(1));
t736 = m(4) * t805 - t846 * mrSges(4,3) + t924 * qJDD(3) + (t850 - t852) * qJD(3) + t884 * t907 + t898;
t777 = -pkin(4) * t845 + qJD(3) * t854 - t875 * t927 - t786;
t765 = -t799 * pkin(5) - pkin(9) * t838 + t811 * t841 + t777;
t900 = m(7) * t765 - t770 * mrSges(7,1) + t771 * mrSges(7,2) - t801 * t790 + t802 * t791;
t893 = -m(6) * t777 + t799 * mrSges(6,1) - t800 * mrSges(6,2) + t840 * t807 - t841 * t808 - t900;
t891 = -m(5) * t786 + qJDD(3) * mrSges(5,3) + qJD(3) * t853 - t893;
t749 = t891 - qJDD(3) * mrSges(4,2) + m(4) * t806 - qJD(3) * t851 + t880 * t907 + (-mrSges(4,3) - mrSges(5,1)) * t845;
t728 = t884 * t736 + t880 * t749;
t818 = -qJDD(1) * pkin(1) + t902;
t899 = -m(3) * t818 + (t887 * mrSges(3,3)) - t728;
t724 = m(2) * t855 - (t887 * mrSges(2,2)) + qJDD(1) * t925 + t899;
t816 = t887 * pkin(1) - t903;
t889 = -m(3) * t816 + (t887 * mrSges(3,2)) + qJDD(1) * mrSges(3,3) - t931;
t734 = m(2) * t856 - (t887 * mrSges(2,1)) - qJDD(1) * mrSges(2,2) + t889;
t919 = t885 * t724 + t881 * t734;
t912 = -t724 * t881 + t885 * t734;
t911 = -t880 * t736 + t884 * t749;
t908 = qJD(1) * (-(t932 * qJD(3)) + (t880 * t933 - t884 * t934) * qJD(1));
t780 = Ifges(7,4) * t802 + Ifges(7,2) * t801 + Ifges(7,6) * t857;
t781 = Ifges(7,1) * t802 + Ifges(7,4) * t801 + Ifges(7,5) * t857;
t897 = mrSges(7,1) * t757 - mrSges(7,2) * t758 + Ifges(7,5) * t771 + Ifges(7,6) * t770 + Ifges(7,3) * t828 + t802 * t780 - t801 * t781;
t779 = Ifges(7,5) * t802 + Ifges(7,6) * t801 + Ifges(7,3) * t857;
t746 = -mrSges(7,1) * t765 + mrSges(7,3) * t758 + Ifges(7,4) * t771 + Ifges(7,2) * t770 + Ifges(7,6) * t828 - t802 * t779 + t781 * t857;
t747 = mrSges(7,2) * t765 - mrSges(7,3) * t757 + Ifges(7,1) * t771 + Ifges(7,4) * t770 + Ifges(7,5) * t828 + t801 * t779 - t780 * t857;
t792 = Ifges(6,5) * t841 + Ifges(6,6) * t840 + Ifges(6,3) * t860;
t794 = Ifges(6,1) * t841 + Ifges(6,4) * t840 + Ifges(6,5) * t860;
t729 = -mrSges(6,1) * t777 + mrSges(6,3) * t763 + Ifges(6,4) * t800 + Ifges(6,2) * t799 + Ifges(6,6) * t839 - pkin(5) * t900 + pkin(9) * t909 + t882 * t746 + t878 * t747 - t841 * t792 + t860 * t794;
t793 = Ifges(6,4) * t841 + Ifges(6,2) * t840 + Ifges(6,6) * t860;
t731 = mrSges(6,2) * t777 - mrSges(6,3) * t762 + Ifges(6,1) * t800 + Ifges(6,4) * t799 + Ifges(6,5) * t839 - pkin(9) * t745 - t746 * t878 + t747 * t882 + t792 * t840 - t793 * t860;
t737 = -t845 * mrSges(5,2) + t894;
t823 = (Ifges(4,5) * qJD(3)) + (Ifges(4,1) * t884 - Ifges(4,4) * t880) * qJD(1);
t825 = (Ifges(5,4) * qJD(3)) + (-Ifges(5,2) * t884 + Ifges(5,6) * t880) * qJD(1);
t720 = -mrSges(4,1) * t812 + mrSges(4,3) * t806 - mrSges(5,1) * t786 + mrSges(5,2) * t785 - t879 * t731 - t883 * t729 - pkin(4) * t893 - pkin(8) * t910 - pkin(3) * t737 - t921 * t846 + (-Ifges(4,2) - Ifges(5,3)) * t845 + t933 * qJDD(3) + (t823 - t825) * qJD(3) + t884 * t908;
t822 = (Ifges(4,6) * qJD(3)) + (Ifges(4,4) * t884 - Ifges(4,2) * t880) * qJD(1);
t824 = (Ifges(5,5) * qJD(3)) + (-Ifges(5,6) * t884 + Ifges(5,3) * t880) * qJD(1);
t890 = mrSges(6,1) * t762 - mrSges(6,2) * t763 + Ifges(6,5) * t800 + Ifges(6,6) * t799 + Ifges(6,3) * t839 + pkin(5) * t745 + t841 * t793 - t840 * t794 + t897;
t722 = t921 * t845 + (-t822 + t824) * qJD(3) + t890 + (Ifges(4,1) + Ifges(5,2)) * t846 + t880 * t908 + mrSges(4,2) * t812 + pkin(4) * t740 - mrSges(5,3) * t785 + mrSges(5,1) * t788 - mrSges(4,3) * t805 + t934 * qJDD(3) - qJ(4) * t737;
t726 = qJDD(1) * mrSges(3,2) - t899;
t892 = mrSges(2,1) * t855 - mrSges(2,2) * t856 + mrSges(3,2) * t818 - mrSges(3,3) * t816 - pkin(1) * t726 - pkin(7) * t728 + qJ(2) * t889 - t720 * t880 + t884 * t722 + (Ifges(3,1) + Ifges(2,3)) * qJDD(1);
t739 = qJDD(3) * mrSges(5,2) + qJD(3) * t852 + t843 * t863 - t898;
t888 = -mrSges(4,2) * t806 - mrSges(5,3) * t786 - pkin(8) * t740 - t879 * t729 + t883 * t731 - pkin(3) * t739 + qJ(4) * (-t843 * t918 + t891) + mrSges(5,2) * t788 + mrSges(4,1) * t805 + t823 * t918 + t822 * t863 + (-t824 * t884 - t825 * t880) * qJD(1) + t934 * t846 + (-qJ(4) * mrSges(5,1) - t933) * t845 + t932 * qJDD(3);
t727 = -m(3) * g(3) + t911;
t719 = pkin(2) * t728 - qJ(2) * t727 - mrSges(2,3) * t855 + (-mrSges(2,2) + mrSges(3,3)) * g(3) + mrSges(3,1) * t818 + t923 * qJDD(1) + (t922 * t887) + t888;
t718 = -mrSges(3,1) * t816 + mrSges(2,3) * t856 - pkin(1) * t727 - pkin(2) * t931 - pkin(7) * t911 + t925 * g(3) - t922 * qJDD(1) - t884 * t720 - t880 * t722 + t923 * t887;
t1 = [-m(1) * g(1) + t912; -m(1) * g(2) + t919; (-m(1) - m(2) - m(3)) * g(3) + t911; -mrSges(1,2) * g(3) + mrSges(1,3) * g(2) - pkin(6) * t919 - t881 * t718 + t885 * t719; mrSges(1,1) * g(3) - mrSges(1,3) * g(1) + pkin(6) * t912 + t885 * t718 + t881 * t719; -mrSges(1,1) * g(2) + mrSges(1,2) * g(1) + t892; t892; t726; t888; t739; t890; t897;];
tauJB  = t1;
