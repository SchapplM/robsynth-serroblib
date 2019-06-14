% Calculate vector of inverse dynamics joint torques and base forces with Newton-Euler
% S6PRRRPR1
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
% pkin [12x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,alpha2,d2,d3,d4,d6,theta1,theta5]';
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
% Datum: 2019-05-05 07:03
% Revision: 89c353f7eff3bd693eda4e29f35b2761dbc3ada0 (2019-05-03)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauJB = S6PRRRPR1_invdynJB_fixb_snew_vp2(qJ, qJD, qJDD, g, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(6,1),zeros(3,1),zeros(12,1),zeros(7,1),zeros(7,3),zeros(7,6)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6PRRRPR1_invdynJB_fixb_snew_vp2: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6PRRRPR1_invdynJB_fixb_snew_vp2: qJD has to be [6x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [6 1]), ...
  'S6PRRRPR1_invdynJB_fixb_snew_vp2: qJDD has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6PRRRPR1_invdynJB_fixb_snew_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [12 1]), ...
  'S6PRRRPR1_invdynJB_fixb_snew_vp2: pkin has to be [12x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6PRRRPR1_invdynJB_fixb_snew_vp2: m has to be [7x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [7,3]), ...
  'S6PRRRPR1_invdynJB_fixb_snew_vp2: mrSges has to be [7x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [7 6]), ...
  'S6PRRRPR1_invdynJB_fixb_snew_vp2: Ifges has to be [7x6] (double)'); 

%% Symbolic Calculation
% From invdyn_fixb_NewtonEuler_linkframe_tauJB_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-05-05 06:59:21
% EndTime: 2019-05-05 06:59:45
% DurationCPUTime: 23.78s
% Computational Cost: add. (411071->343), mult. (877776->443), div. (0->0), fcn. (652633->14), ass. (0->147)
t888 = sin(pkin(11));
t891 = cos(pkin(11));
t873 = g(1) * t888 - g(2) * t891;
t874 = -g(1) * t891 - g(2) * t888;
t886 = -g(3) + qJDD(1);
t889 = sin(pkin(6));
t892 = cos(pkin(6));
t896 = sin(qJ(2));
t900 = cos(qJ(2));
t845 = -t896 * t874 + (t873 * t892 + t886 * t889) * t900;
t924 = t892 * t896;
t925 = t889 * t896;
t846 = t873 * t924 + t900 * t874 + t886 * t925;
t901 = qJD(2) ^ 2;
t839 = -pkin(2) * t901 + qJDD(2) * pkin(8) + t846;
t857 = -t873 * t889 + t886 * t892;
t895 = sin(qJ(3));
t899 = cos(qJ(3));
t825 = -t895 * t839 + t899 * t857;
t920 = qJD(2) * qJD(3);
t919 = t899 * t920;
t871 = qJDD(2) * t895 + t919;
t810 = (-t871 + t919) * pkin(9) + (t895 * t899 * t901 + qJDD(3)) * pkin(3) + t825;
t826 = t899 * t839 + t895 * t857;
t872 = qJDD(2) * t899 - t895 * t920;
t922 = qJD(2) * t895;
t878 = qJD(3) * pkin(3) - pkin(9) * t922;
t885 = t899 ^ 2;
t813 = -pkin(3) * t885 * t901 + pkin(9) * t872 - qJD(3) * t878 + t826;
t894 = sin(qJ(4));
t898 = cos(qJ(4));
t790 = t898 * t810 - t894 * t813;
t863 = (-t894 * t895 + t898 * t899) * qJD(2);
t835 = qJD(4) * t863 + t871 * t898 + t872 * t894;
t864 = (t894 * t899 + t895 * t898) * qJD(2);
t883 = qJDD(3) + qJDD(4);
t884 = qJD(3) + qJD(4);
t786 = (t863 * t884 - t835) * qJ(5) + (t863 * t864 + t883) * pkin(4) + t790;
t791 = t894 * t810 + t898 * t813;
t834 = -qJD(4) * t864 - t871 * t894 + t872 * t898;
t855 = pkin(4) * t884 - qJ(5) * t864;
t859 = t863 ^ 2;
t788 = -pkin(4) * t859 + qJ(5) * t834 - t855 * t884 + t791;
t887 = sin(pkin(12));
t890 = cos(pkin(12));
t849 = t863 * t890 - t864 * t887;
t927 = 2 * qJD(5);
t783 = t887 * t786 + t890 * t788 + t849 * t927;
t811 = t834 * t890 - t835 * t887;
t850 = t863 * t887 + t864 * t890;
t822 = -mrSges(6,1) * t849 + mrSges(6,2) * t850;
t837 = mrSges(6,1) * t884 - mrSges(6,3) * t850;
t823 = -pkin(5) * t849 - pkin(10) * t850;
t882 = t884 ^ 2;
t780 = -pkin(5) * t882 + pkin(10) * t883 + t823 * t849 + t783;
t907 = -qJDD(2) * pkin(2) - t845;
t824 = -t872 * pkin(3) + t878 * t922 + (-pkin(9) * t885 - pkin(8)) * t901 + t907;
t793 = -t834 * pkin(4) - t859 * qJ(5) + t864 * t855 + qJDD(5) + t824;
t812 = t834 * t887 + t835 * t890;
t784 = (-t849 * t884 - t812) * pkin(10) + (t850 * t884 - t811) * pkin(5) + t793;
t893 = sin(qJ(6));
t897 = cos(qJ(6));
t777 = -t780 * t893 + t784 * t897;
t831 = -t850 * t893 + t884 * t897;
t796 = qJD(6) * t831 + t812 * t897 + t883 * t893;
t809 = qJDD(6) - t811;
t832 = t850 * t897 + t884 * t893;
t814 = -mrSges(7,1) * t831 + mrSges(7,2) * t832;
t841 = qJD(6) - t849;
t815 = -mrSges(7,2) * t841 + mrSges(7,3) * t831;
t773 = m(7) * t777 + mrSges(7,1) * t809 - mrSges(7,3) * t796 - t814 * t832 + t815 * t841;
t778 = t780 * t897 + t784 * t893;
t795 = -qJD(6) * t832 - t812 * t893 + t883 * t897;
t816 = mrSges(7,1) * t841 - mrSges(7,3) * t832;
t774 = m(7) * t778 - mrSges(7,2) * t809 + mrSges(7,3) * t795 + t814 * t831 - t816 * t841;
t914 = -t773 * t893 + t897 * t774;
t759 = m(6) * t783 - mrSges(6,2) * t883 + mrSges(6,3) * t811 + t822 * t849 - t837 * t884 + t914;
t912 = -t890 * t786 + t887 * t788;
t782 = -0.2e1 * qJD(5) * t850 - t912;
t836 = -mrSges(6,2) * t884 + mrSges(6,3) * t849;
t779 = -t883 * pkin(5) - t882 * pkin(10) + (t927 + t823) * t850 + t912;
t908 = -m(7) * t779 + t795 * mrSges(7,1) - mrSges(7,2) * t796 + t831 * t815 - t816 * t832;
t769 = m(6) * t782 + mrSges(6,1) * t883 - mrSges(6,3) * t812 - t822 * t850 + t836 * t884 + t908;
t753 = t887 * t759 + t890 * t769;
t851 = -mrSges(5,1) * t863 + mrSges(5,2) * t864;
t854 = -mrSges(5,2) * t884 + mrSges(5,3) * t863;
t750 = m(5) * t790 + mrSges(5,1) * t883 - mrSges(5,3) * t835 - t851 * t864 + t854 * t884 + t753;
t856 = mrSges(5,1) * t884 - mrSges(5,3) * t864;
t915 = t890 * t759 - t769 * t887;
t751 = m(5) * t791 - mrSges(5,2) * t883 + mrSges(5,3) * t834 + t851 * t863 - t856 * t884 + t915;
t744 = t898 * t750 + t894 * t751;
t861 = Ifges(4,6) * qJD(3) + (Ifges(4,4) * t895 + Ifges(4,2) * t899) * qJD(2);
t862 = Ifges(4,5) * qJD(3) + (Ifges(4,1) * t895 + Ifges(4,4) * t899) * qJD(2);
t797 = Ifges(7,5) * t832 + Ifges(7,6) * t831 + Ifges(7,3) * t841;
t799 = Ifges(7,1) * t832 + Ifges(7,4) * t831 + Ifges(7,5) * t841;
t766 = -mrSges(7,1) * t779 + mrSges(7,3) * t778 + Ifges(7,4) * t796 + Ifges(7,2) * t795 + Ifges(7,6) * t809 - t797 * t832 + t799 * t841;
t798 = Ifges(7,4) * t832 + Ifges(7,2) * t831 + Ifges(7,6) * t841;
t767 = mrSges(7,2) * t779 - mrSges(7,3) * t777 + Ifges(7,1) * t796 + Ifges(7,4) * t795 + Ifges(7,5) * t809 + t797 * t831 - t798 * t841;
t818 = Ifges(6,4) * t850 + Ifges(6,2) * t849 + Ifges(6,6) * t884;
t819 = Ifges(6,1) * t850 + Ifges(6,4) * t849 + Ifges(6,5) * t884;
t843 = Ifges(5,4) * t864 + Ifges(5,2) * t863 + Ifges(5,6) * t884;
t844 = Ifges(5,1) * t864 + Ifges(5,4) * t863 + Ifges(5,5) * t884;
t904 = -mrSges(5,1) * t790 - mrSges(6,1) * t782 + mrSges(5,2) * t791 + mrSges(6,2) * t783 - pkin(4) * t753 - pkin(5) * t908 - pkin(10) * t914 - t897 * t766 - t893 * t767 - t850 * t818 + t849 * t819 + t863 * t844 - Ifges(6,6) * t811 - Ifges(6,5) * t812 - t864 * t843 - Ifges(5,6) * t834 - Ifges(5,5) * t835 + (-Ifges(5,3) - Ifges(6,3)) * t883;
t928 = mrSges(4,1) * t825 - mrSges(4,2) * t826 + Ifges(4,5) * t871 + Ifges(4,6) * t872 + Ifges(4,3) * qJDD(3) + pkin(3) * t744 + (t861 * t895 - t862 * t899) * qJD(2) - t904;
t838 = -t901 * pkin(8) + t907;
t875 = qJD(3) * mrSges(4,1) - mrSges(4,3) * t922;
t921 = qJD(2) * t899;
t876 = -qJD(3) * mrSges(4,2) + mrSges(4,3) * t921;
t762 = t897 * t773 + t893 * t774;
t760 = m(6) * t793 - t811 * mrSges(6,1) + t812 * mrSges(6,2) - t849 * t836 + t850 * t837 + t762;
t906 = m(5) * t824 - t834 * mrSges(5,1) + mrSges(5,2) * t835 - t863 * t854 + t856 * t864 + t760;
t903 = -m(4) * t838 + t872 * mrSges(4,1) - mrSges(4,2) * t871 - t875 * t922 + t876 * t921 - t906;
t756 = m(3) * t845 + qJDD(2) * mrSges(3,1) - mrSges(3,2) * t901 + t903;
t926 = t756 * t900;
t870 = (-mrSges(4,1) * t899 + mrSges(4,2) * t895) * qJD(2);
t742 = m(4) * t825 + qJDD(3) * mrSges(4,1) - mrSges(4,3) * t871 + qJD(3) * t876 - t870 * t922 + t744;
t916 = -t750 * t894 + t898 * t751;
t743 = m(4) * t826 - qJDD(3) * mrSges(4,2) + mrSges(4,3) * t872 - qJD(3) * t875 + t870 * t921 + t916;
t917 = -t742 * t895 + t899 * t743;
t734 = m(3) * t846 - mrSges(3,1) * t901 - qJDD(2) * mrSges(3,2) + t917;
t737 = t899 * t742 + t895 * t743;
t736 = m(3) * t857 + t737;
t725 = t734 * t924 - t736 * t889 + t892 * t926;
t723 = m(2) * t873 + t725;
t729 = t900 * t734 - t756 * t896;
t728 = m(2) * t874 + t729;
t923 = t891 * t723 + t888 * t728;
t724 = t734 * t925 + t892 * t736 + t889 * t926;
t918 = -t723 * t888 + t891 * t728;
t913 = m(2) * t886 + t724;
t817 = Ifges(6,5) * t850 + Ifges(6,6) * t849 + Ifges(6,3) * t884;
t745 = mrSges(6,2) * t793 - mrSges(6,3) * t782 + Ifges(6,1) * t812 + Ifges(6,4) * t811 + Ifges(6,5) * t883 - pkin(10) * t762 - t766 * t893 + t767 * t897 + t817 * t849 - t818 * t884;
t905 = mrSges(7,1) * t777 - mrSges(7,2) * t778 + Ifges(7,5) * t796 + Ifges(7,6) * t795 + Ifges(7,3) * t809 + t798 * t832 - t799 * t831;
t746 = -mrSges(6,1) * t793 + mrSges(6,3) * t783 + Ifges(6,4) * t812 + Ifges(6,2) * t811 + Ifges(6,6) * t883 - pkin(5) * t762 - t817 * t850 + t819 * t884 - t905;
t842 = Ifges(5,5) * t864 + Ifges(5,6) * t863 + Ifges(5,3) * t884;
t730 = -mrSges(5,1) * t824 + mrSges(5,3) * t791 + Ifges(5,4) * t835 + Ifges(5,2) * t834 + Ifges(5,6) * t883 - pkin(4) * t760 + qJ(5) * t915 + t887 * t745 + t890 * t746 - t864 * t842 + t884 * t844;
t738 = mrSges(5,2) * t824 - mrSges(5,3) * t790 + Ifges(5,1) * t835 + Ifges(5,4) * t834 + Ifges(5,5) * t883 - qJ(5) * t753 + t745 * t890 - t746 * t887 + t842 * t863 - t843 * t884;
t860 = Ifges(4,3) * qJD(3) + (Ifges(4,5) * t895 + Ifges(4,6) * t899) * qJD(2);
t719 = -mrSges(4,1) * t838 + mrSges(4,3) * t826 + Ifges(4,4) * t871 + Ifges(4,2) * t872 + Ifges(4,6) * qJDD(3) - pkin(3) * t906 + pkin(9) * t916 + qJD(3) * t862 + t898 * t730 + t894 * t738 - t860 * t922;
t720 = mrSges(4,2) * t838 - mrSges(4,3) * t825 + Ifges(4,1) * t871 + Ifges(4,4) * t872 + Ifges(4,5) * qJDD(3) - pkin(9) * t744 - qJD(3) * t861 - t730 * t894 + t738 * t898 + t860 * t921;
t718 = mrSges(3,2) * t857 - mrSges(3,3) * t845 + Ifges(3,5) * qJDD(2) - Ifges(3,6) * t901 - pkin(8) * t737 - t719 * t895 + t720 * t899;
t721 = -mrSges(3,1) * t857 + mrSges(3,3) * t846 + t901 * Ifges(3,5) + Ifges(3,6) * qJDD(2) - pkin(2) * t737 - t928;
t909 = pkin(7) * t729 + t718 * t896 + t721 * t900;
t717 = mrSges(3,1) * t845 - mrSges(3,2) * t846 + Ifges(3,3) * qJDD(2) + pkin(2) * t903 + pkin(8) * t917 + t899 * t719 + t895 * t720;
t716 = mrSges(2,2) * t886 - mrSges(2,3) * t873 + t900 * t718 - t896 * t721 + (-t724 * t889 - t725 * t892) * pkin(7);
t715 = -mrSges(2,1) * t886 + mrSges(2,3) * t874 - pkin(1) * t724 - t889 * t717 + t909 * t892;
t1 = [-m(1) * g(1) + t918; -m(1) * g(2) + t923; -m(1) * g(3) + t913; -mrSges(1,2) * g(3) + mrSges(1,3) * g(2) - qJ(1) * t923 - t888 * t715 + t891 * t716; mrSges(1,1) * g(3) - mrSges(1,3) * g(1) + qJ(1) * t918 + t891 * t715 + t888 * t716; -mrSges(1,1) * g(2) + mrSges(2,1) * t873 + mrSges(1,2) * g(1) - mrSges(2,2) * t874 + pkin(1) * t725 + t892 * t717 + t909 * t889; t913; t717; t928; -t904; t760; t905;];
tauJB  = t1;
