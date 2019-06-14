% Calculate vector of inverse dynamics joint torques and base forces with Newton-Euler
% S6PRRRPP1
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
%   pkin=[a2,a3,a4,a5,a6,alpha2,d2,d3,d4,theta1,theta5]';
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
% Datum: 2019-05-05 06:40
% Revision: 89c353f7eff3bd693eda4e29f35b2761dbc3ada0 (2019-05-03)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauJB = S6PRRRPP1_invdynJB_fixb_snew_vp2(qJ, qJD, qJDD, g, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(6,1),zeros(3,1),zeros(11,1),zeros(7,1),zeros(7,3),zeros(7,6)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6PRRRPP1_invdynJB_fixb_snew_vp2: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6PRRRPP1_invdynJB_fixb_snew_vp2: qJD has to be [6x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [6 1]), ...
  'S6PRRRPP1_invdynJB_fixb_snew_vp2: qJDD has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6PRRRPP1_invdynJB_fixb_snew_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [11 1]), ...
  'S6PRRRPP1_invdynJB_fixb_snew_vp2: pkin has to be [11x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6PRRRPP1_invdynJB_fixb_snew_vp2: m has to be [7x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [7,3]), ...
  'S6PRRRPP1_invdynJB_fixb_snew_vp2: mrSges has to be [7x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [7 6]), ...
  'S6PRRRPP1_invdynJB_fixb_snew_vp2: Ifges has to be [7x6] (double)'); 

%% Symbolic Calculation
% From invdyn_fixb_NewtonEuler_linkframe_tauJB_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-05-05 06:34:16
% EndTime: 2019-05-05 06:34:27
% DurationCPUTime: 11.18s
% Computational Cost: add. (185477->317), mult. (367161->396), div. (0->0), fcn. (257908->12), ass. (0->137)
t892 = Ifges(6,1) + Ifges(7,1);
t885 = Ifges(6,4) - Ifges(7,5);
t884 = Ifges(6,5) + Ifges(7,4);
t891 = -Ifges(6,2) - Ifges(7,3);
t890 = -Ifges(7,2) - Ifges(6,3);
t883 = Ifges(6,6) - Ifges(7,6);
t843 = sin(pkin(10));
t845 = cos(pkin(10));
t833 = t843 * g(1) - t845 * g(2);
t834 = -t845 * g(1) - t843 * g(2);
t841 = -g(3) + qJDD(1);
t844 = sin(pkin(6));
t846 = cos(pkin(6));
t849 = sin(qJ(2));
t852 = cos(qJ(2));
t793 = -t849 * t834 + (t833 * t846 + t841 * t844) * t852;
t879 = t846 * t849;
t880 = t844 * t849;
t794 = t833 * t879 + t852 * t834 + t841 * t880;
t854 = qJD(2) ^ 2;
t788 = -t854 * pkin(2) + qJDD(2) * pkin(8) + t794;
t812 = -t844 * t833 + t846 * t841;
t848 = sin(qJ(3));
t851 = cos(qJ(3));
t778 = -t848 * t788 + t851 * t812;
t830 = (-pkin(3) * t851 - pkin(9) * t848) * qJD(2);
t853 = qJD(3) ^ 2;
t873 = qJD(2) * t848;
t759 = -qJDD(3) * pkin(3) - t853 * pkin(9) + t830 * t873 - t778;
t847 = sin(qJ(4));
t850 = cos(qJ(4));
t828 = t847 * qJD(3) + t850 * t873;
t871 = qJD(2) * qJD(3);
t867 = t851 * t871;
t831 = t848 * qJDD(2) + t867;
t803 = -t828 * qJD(4) + t850 * qJDD(3) - t847 * t831;
t872 = t851 * qJD(2);
t839 = qJD(4) - t872;
t810 = t839 * pkin(4) - t828 * qJ(5);
t827 = t850 * qJD(3) - t847 * t873;
t823 = t827 ^ 2;
t757 = -t803 * pkin(4) - t823 * qJ(5) + t828 * t810 + qJDD(5) + t759;
t804 = t827 * qJD(4) + t847 * qJDD(3) + t850 * t831;
t842 = sin(pkin(11));
t882 = cos(pkin(11));
t773 = -t882 * t803 + t842 * t804;
t774 = t842 * t803 + t882 * t804;
t805 = -t882 * t827 + t842 * t828;
t806 = t842 * t827 + t882 * t828;
t750 = -0.2e1 * qJD(6) * t806 + (t805 * t839 - t774) * qJ(6) + (t806 * t839 + t773) * pkin(5) + t757;
t791 = -t839 * mrSges(7,1) + t806 * mrSges(7,2);
t792 = -t805 * mrSges(7,2) + t839 * mrSges(7,3);
t743 = m(7) * t750 + t773 * mrSges(7,1) - t774 * mrSges(7,3) - t806 * t791 + t805 * t792;
t779 = t851 * t788 + t848 * t812;
t760 = -t853 * pkin(3) + qJDD(3) * pkin(9) + t830 * t872 + t779;
t787 = -qJDD(2) * pkin(2) - t854 * pkin(8) - t793;
t868 = t848 * t871;
t832 = t851 * qJDD(2) - t868;
t763 = (-t831 - t867) * pkin(9) + (-t832 + t868) * pkin(3) + t787;
t755 = -t847 * t760 + t850 * t763;
t824 = qJDD(4) - t832;
t752 = (t827 * t839 - t804) * qJ(5) + (t827 * t828 + t824) * pkin(4) + t755;
t756 = t850 * t760 + t847 * t763;
t754 = -t823 * pkin(4) + t803 * qJ(5) - t839 * t810 + t756;
t887 = -2 * qJD(5);
t748 = t842 * t752 + t882 * t754 + t805 * t887;
t780 = t805 * pkin(5) - t806 * qJ(6);
t838 = t839 ^ 2;
t745 = -t838 * pkin(5) + t824 * qJ(6) + 0.2e1 * qJD(6) * t839 - t805 * t780 + t748;
t875 = -t885 * t805 + t806 * t892 + t884 * t839;
t877 = t883 * t805 - t884 * t806 + t890 * t839;
t730 = -mrSges(6,1) * t757 - mrSges(7,1) * t750 + mrSges(7,2) * t745 + mrSges(6,3) * t748 - pkin(5) * t743 + t773 * t891 + t885 * t774 + t877 * t806 + t883 * t824 + t875 * t839;
t859 = t882 * t752 - t842 * t754;
t746 = -t824 * pkin(5) - t838 * qJ(6) + qJDD(6) + ((2 * qJD(5)) + t780) * t806 - t859;
t747 = t806 * t887 + t859;
t876 = t891 * t805 + t885 * t806 + t883 * t839;
t731 = mrSges(6,2) * t757 + mrSges(7,2) * t746 - mrSges(6,3) * t747 - mrSges(7,3) * t750 - qJ(6) * t743 - t885 * t773 + t774 * t892 + t877 * t805 + t884 * t824 - t876 * t839;
t789 = -t839 * mrSges(6,2) - t805 * mrSges(6,3);
t790 = t839 * mrSges(6,1) - t806 * mrSges(6,3);
t740 = m(6) * t757 + t773 * mrSges(6,1) + t774 * mrSges(6,2) + t805 * t789 + t806 * t790 + t743;
t796 = Ifges(5,5) * t828 + Ifges(5,6) * t827 + Ifges(5,3) * t839;
t798 = Ifges(5,1) * t828 + Ifges(5,4) * t827 + Ifges(5,5) * t839;
t869 = m(7) * t745 + t824 * mrSges(7,3) + t839 * t791;
t781 = t805 * mrSges(7,1) - t806 * mrSges(7,3);
t874 = -t805 * mrSges(6,1) - t806 * mrSges(6,2) - t781;
t886 = -mrSges(6,3) - mrSges(7,2);
t734 = m(6) * t748 - t824 * mrSges(6,2) + t886 * t773 - t839 * t790 + t874 * t805 + t869;
t862 = -m(7) * t746 + t824 * mrSges(7,1) + t839 * t792;
t736 = m(6) * t747 + t824 * mrSges(6,1) + t886 * t774 + t839 * t789 + t874 * t806 + t862;
t864 = t882 * t734 - t842 * t736;
t709 = -mrSges(5,1) * t759 + mrSges(5,3) * t756 + Ifges(5,4) * t804 + Ifges(5,2) * t803 + Ifges(5,6) * t824 - pkin(4) * t740 + qJ(5) * t864 + t882 * t730 + t842 * t731 - t828 * t796 + t839 * t798;
t729 = t842 * t734 + t882 * t736;
t797 = Ifges(5,4) * t828 + Ifges(5,2) * t827 + Ifges(5,6) * t839;
t710 = mrSges(5,2) * t759 - mrSges(5,3) * t755 + Ifges(5,1) * t804 + Ifges(5,4) * t803 + Ifges(5,5) * t824 - qJ(5) * t729 - t842 * t730 + t882 * t731 + t827 * t796 - t839 * t797;
t807 = -t827 * mrSges(5,1) + t828 * mrSges(5,2);
t809 = -t839 * mrSges(5,2) + t827 * mrSges(5,3);
t727 = m(5) * t755 + t824 * mrSges(5,1) - t804 * mrSges(5,3) - t828 * t807 + t839 * t809 + t729;
t811 = t839 * mrSges(5,1) - t828 * mrSges(5,3);
t728 = m(5) * t756 - t824 * mrSges(5,2) + t803 * mrSges(5,3) + t827 * t807 - t839 * t811 + t864;
t725 = -t847 * t727 + t850 * t728;
t739 = -m(5) * t759 + t803 * mrSges(5,1) - t804 * mrSges(5,2) + t827 * t809 - t828 * t811 - t740;
t817 = Ifges(4,6) * qJD(3) + (Ifges(4,4) * t848 + Ifges(4,2) * t851) * qJD(2);
t818 = Ifges(4,5) * qJD(3) + (Ifges(4,1) * t848 + Ifges(4,4) * t851) * qJD(2);
t889 = mrSges(4,1) * t778 - mrSges(4,2) * t779 + Ifges(4,5) * t831 + Ifges(4,6) * t832 + Ifges(4,3) * qJDD(3) + pkin(3) * t739 + pkin(9) * t725 + t850 * t709 + t847 * t710 + (t848 * t817 - t851 * t818) * qJD(2);
t742 = t774 * mrSges(7,2) + t806 * t781 - t862;
t888 = -t883 * t773 + t884 * t774 + t875 * t805 + t876 * t806 + (Ifges(5,3) - t890) * t824 + mrSges(5,1) * t755 + mrSges(6,1) * t747 - mrSges(7,1) * t746 - mrSges(5,2) * t756 - mrSges(6,2) * t748 + mrSges(7,3) * t745 + Ifges(5,5) * t804 + Ifges(5,6) * t803 + pkin(4) * t729 - pkin(5) * t742 + qJ(6) * (-t773 * mrSges(7,2) - t805 * t781 + t869) + t828 * t797 - t827 * t798;
t724 = t850 * t727 + t847 * t728;
t835 = qJD(3) * mrSges(4,1) - mrSges(4,3) * t873;
t836 = -qJD(3) * mrSges(4,2) + mrSges(4,3) * t872;
t857 = -m(4) * t787 + t832 * mrSges(4,1) - t831 * mrSges(4,2) - t835 * t873 + t836 * t872 - t724;
t720 = m(3) * t793 + qJDD(2) * mrSges(3,1) - t854 * mrSges(3,2) + t857;
t881 = t720 * t852;
t829 = (-mrSges(4,1) * t851 + mrSges(4,2) * t848) * qJD(2);
t723 = m(4) * t779 - qJDD(3) * mrSges(4,2) + t832 * mrSges(4,3) - qJD(3) * t835 + t829 * t872 + t725;
t738 = m(4) * t778 + qJDD(3) * mrSges(4,1) - t831 * mrSges(4,3) + qJD(3) * t836 - t829 * t873 + t739;
t865 = t851 * t723 - t848 * t738;
t714 = m(3) * t794 - t854 * mrSges(3,1) - qJDD(2) * mrSges(3,2) + t865;
t717 = t848 * t723 + t851 * t738;
t716 = m(3) * t812 + t717;
t703 = t714 * t879 - t844 * t716 + t846 * t881;
t701 = m(2) * t833 + t703;
t707 = t852 * t714 - t849 * t720;
t706 = m(2) * t834 + t707;
t878 = t845 * t701 + t843 * t706;
t702 = t714 * t880 + t846 * t716 + t844 * t881;
t866 = -t843 * t701 + t845 * t706;
t863 = m(2) * t841 + t702;
t816 = Ifges(4,3) * qJD(3) + (Ifges(4,5) * t848 + Ifges(4,6) * t851) * qJD(2);
t699 = mrSges(4,2) * t787 - mrSges(4,3) * t778 + Ifges(4,1) * t831 + Ifges(4,4) * t832 + Ifges(4,5) * qJDD(3) - pkin(9) * t724 - qJD(3) * t817 - t847 * t709 + t850 * t710 + t816 * t872;
t708 = -mrSges(4,1) * t787 + mrSges(4,3) * t779 + Ifges(4,4) * t831 + Ifges(4,2) * t832 + Ifges(4,6) * qJDD(3) - pkin(3) * t724 + qJD(3) * t818 - t816 * t873 - t888;
t697 = mrSges(3,2) * t812 - mrSges(3,3) * t793 + Ifges(3,5) * qJDD(2) - t854 * Ifges(3,6) - pkin(8) * t717 + t851 * t699 - t848 * t708;
t698 = -mrSges(3,1) * t812 + mrSges(3,3) * t794 + t854 * Ifges(3,5) + Ifges(3,6) * qJDD(2) - pkin(2) * t717 - t889;
t858 = pkin(7) * t707 + t697 * t849 + t698 * t852;
t696 = mrSges(3,1) * t793 - mrSges(3,2) * t794 + Ifges(3,3) * qJDD(2) + pkin(2) * t857 + pkin(8) * t865 + t848 * t699 + t851 * t708;
t695 = mrSges(2,2) * t841 - mrSges(2,3) * t833 + t852 * t697 - t849 * t698 + (-t702 * t844 - t703 * t846) * pkin(7);
t694 = -mrSges(2,1) * t841 + mrSges(2,3) * t834 - pkin(1) * t702 - t844 * t696 + t858 * t846;
t1 = [-m(1) * g(1) + t866; -m(1) * g(2) + t878; -m(1) * g(3) + t863; -mrSges(1,2) * g(3) + mrSges(1,3) * g(2) - qJ(1) * t878 - t843 * t694 + t845 * t695; mrSges(1,1) * g(3) - mrSges(1,3) * g(1) + qJ(1) * t866 + t845 * t694 + t843 * t695; -mrSges(1,1) * g(2) + mrSges(2,1) * t833 + mrSges(1,2) * g(1) - mrSges(2,2) * t834 + pkin(1) * t703 + t846 * t696 + t858 * t844; t863; t696; t889; t888; t740; t742;];
tauJB  = t1;
