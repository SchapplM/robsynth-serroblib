% Calculate vector of inverse dynamics joint torques and base forces with Newton-Euler
% S5RRRRR10
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% qJD [5x1]
%   Generalized joint velocities
% qJDD [5x1]
%   Generalized joint accelerations
% g [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [10x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,alpha2,d1,d2,d3,d4,d5]';
% m_mdh [6x1]
%   mass of all robot links (including the base)
% mrSges [6x3]
%  first moment of all robot links (mass times center of mass in body frames)
%  rows: links of the robot (starting with base)
%  columns: x-, y-, z-coordinates
% Ifges [6x6]
%   inertia of all robot links about their respective body frame origins, in body frames
%   rows: links of the robot (starting with base)
%   columns: xx, yy, zz, xy, xz, yz (see inertial_parameters_convert_par1_par2.m)
% 
% Output:
% tauJB [(6+5)x1]
%   joint torques and base forces of inverse dynamics (contains inertial, gravitational coriolis and centrifugal forces)

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-31 22:37
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauJB = S5RRRRR10_invdynJB_fixb_snew_vp2(qJ, qJD, qJDD, g, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(5,1),zeros(3,1),zeros(10,1),zeros(6,1),zeros(6,3),zeros(6,6)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RRRRR10_invdynJB_fixb_snew_vp2: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RRRRR10_invdynJB_fixb_snew_vp2: qJD has to be [5x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [5 1]), ...
  'S5RRRRR10_invdynJB_fixb_snew_vp2: qJDD has to be [5x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S5RRRRR10_invdynJB_fixb_snew_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S5RRRRR10_invdynJB_fixb_snew_vp2: pkin has to be [10x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5RRRRR10_invdynJB_fixb_snew_vp2: m has to be [6x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [6,3]), ...
  'S5RRRRR10_invdynJB_fixb_snew_vp2: mrSges has to be [6x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [6 6]), ...
  'S5RRRRR10_invdynJB_fixb_snew_vp2: Ifges has to be [6x6] (double)'); 

%% Symbolic Calculation
% From invdyn_fixb_NewtonEuler_linkframe_tauJB_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 22:33:16
% EndTime: 2019-12-31 22:33:33
% DurationCPUTime: 16.55s
% Computational Cost: add. (269671->328), mult. (578983->427), div. (0->0), fcn. (456982->12), ass. (0->140)
t798 = sin(pkin(5));
t834 = pkin(7) * t798;
t799 = cos(pkin(5));
t833 = t799 * g(3);
t803 = sin(qJ(2));
t832 = t798 * t803;
t808 = cos(qJ(2));
t831 = t798 * t808;
t830 = t799 * t803;
t829 = t799 * t808;
t804 = sin(qJ(1));
t809 = cos(qJ(1));
t789 = t804 * g(1) - t809 * g(2);
t810 = qJD(1) ^ 2;
t779 = qJDD(1) * pkin(1) + t810 * t834 + t789;
t790 = -t809 * g(1) - t804 * g(2);
t824 = qJDD(1) * t798;
t780 = -t810 * pkin(1) + pkin(7) * t824 + t790;
t827 = t779 * t830 + t808 * t780;
t755 = -g(3) * t832 + t827;
t794 = t799 * qJD(1) + qJD(2);
t826 = qJD(1) * t798;
t823 = t803 * t826;
t777 = t794 * mrSges(3,1) - mrSges(3,3) * t823;
t781 = (-mrSges(3,1) * t808 + mrSges(3,2) * t803) * t826;
t784 = -qJD(2) * t823 + t808 * t824;
t793 = t799 * qJDD(1) + qJDD(2);
t782 = (-pkin(2) * t808 - pkin(8) * t803) * t826;
t792 = t794 ^ 2;
t825 = qJD(1) * t808;
t740 = -t792 * pkin(2) + t793 * pkin(8) + (-g(3) * t803 + t782 * t825) * t798 + t827;
t783 = (qJD(2) * t825 + qJDD(1) * t803) * t798;
t741 = -t784 * pkin(2) - t783 * pkin(8) - t833 + (-t779 + (pkin(2) * t803 - pkin(8) * t808) * t794 * qJD(1)) * t798;
t802 = sin(qJ(3));
t807 = cos(qJ(3));
t714 = -t802 * t740 + t807 * t741;
t771 = t807 * t794 - t802 * t823;
t753 = t771 * qJD(3) + t807 * t783 + t802 * t793;
t772 = t802 * t794 + t807 * t823;
t776 = qJDD(3) - t784;
t822 = t798 * t825;
t788 = qJD(3) - t822;
t707 = (t771 * t788 - t753) * pkin(9) + (t771 * t772 + t776) * pkin(3) + t714;
t715 = t807 * t740 + t802 * t741;
t752 = -t772 * qJD(3) - t802 * t783 + t807 * t793;
t762 = t788 * pkin(3) - t772 * pkin(9);
t770 = t771 ^ 2;
t709 = -t770 * pkin(3) + t752 * pkin(9) - t788 * t762 + t715;
t801 = sin(qJ(4));
t806 = cos(qJ(4));
t705 = t801 * t707 + t806 * t709;
t758 = t801 * t771 + t806 * t772;
t724 = -t758 * qJD(4) + t806 * t752 - t801 * t753;
t757 = t806 * t771 - t801 * t772;
t734 = -t757 * mrSges(5,1) + t758 * mrSges(5,2);
t786 = qJD(4) + t788;
t745 = t786 * mrSges(5,1) - t758 * mrSges(5,3);
t775 = qJDD(4) + t776;
t735 = -t757 * pkin(4) - t758 * pkin(10);
t785 = t786 ^ 2;
t701 = -t785 * pkin(4) + t775 * pkin(10) + t757 * t735 + t705;
t754 = -g(3) * t831 + t779 * t829 - t803 * t780;
t739 = -t793 * pkin(2) - t792 * pkin(8) + t782 * t823 - t754;
t713 = -t752 * pkin(3) - t770 * pkin(9) + t772 * t762 + t739;
t725 = t757 * qJD(4) + t801 * t752 + t806 * t753;
t702 = (-t757 * t786 - t725) * pkin(10) + (t758 * t786 - t724) * pkin(4) + t713;
t800 = sin(qJ(5));
t805 = cos(qJ(5));
t698 = -t800 * t701 + t805 * t702;
t742 = -t800 * t758 + t805 * t786;
t712 = t742 * qJD(5) + t805 * t725 + t800 * t775;
t723 = qJDD(5) - t724;
t743 = t805 * t758 + t800 * t786;
t727 = -t742 * mrSges(6,1) + t743 * mrSges(6,2);
t756 = qJD(5) - t757;
t728 = -t756 * mrSges(6,2) + t742 * mrSges(6,3);
t694 = m(6) * t698 + t723 * mrSges(6,1) - t712 * mrSges(6,3) - t743 * t727 + t756 * t728;
t699 = t805 * t701 + t800 * t702;
t711 = -t743 * qJD(5) - t800 * t725 + t805 * t775;
t729 = t756 * mrSges(6,1) - t743 * mrSges(6,3);
t695 = m(6) * t699 - t723 * mrSges(6,2) + t711 * mrSges(6,3) + t742 * t727 - t756 * t729;
t818 = -t800 * t694 + t805 * t695;
t681 = m(5) * t705 - t775 * mrSges(5,2) + t724 * mrSges(5,3) + t757 * t734 - t786 * t745 + t818;
t704 = t806 * t707 - t801 * t709;
t744 = -t786 * mrSges(5,2) + t757 * mrSges(5,3);
t700 = -t775 * pkin(4) - t785 * pkin(10) + t758 * t735 - t704;
t817 = -m(6) * t700 + t711 * mrSges(6,1) - t712 * mrSges(6,2) + t742 * t728 - t743 * t729;
t690 = m(5) * t704 + t775 * mrSges(5,1) - t725 * mrSges(5,3) - t758 * t734 + t786 * t744 + t817;
t675 = t801 * t681 + t806 * t690;
t759 = -t771 * mrSges(4,1) + t772 * mrSges(4,2);
t760 = -t788 * mrSges(4,2) + t771 * mrSges(4,3);
t673 = m(4) * t714 + t776 * mrSges(4,1) - t753 * mrSges(4,3) - t772 * t759 + t788 * t760 + t675;
t761 = t788 * mrSges(4,1) - t772 * mrSges(4,3);
t819 = t806 * t681 - t801 * t690;
t674 = m(4) * t715 - t776 * mrSges(4,2) + t752 * mrSges(4,3) + t771 * t759 - t788 * t761 + t819;
t820 = -t802 * t673 + t807 * t674;
t664 = m(3) * t755 - t793 * mrSges(3,2) + t784 * mrSges(3,3) - t794 * t777 + t781 * t822 + t820;
t667 = t807 * t673 + t802 * t674;
t766 = -t798 * t779 - t833;
t778 = -t794 * mrSges(3,2) + mrSges(3,3) * t822;
t666 = m(3) * t766 - t784 * mrSges(3,1) + t783 * mrSges(3,2) + (t777 * t803 - t778 * t808) * t826 + t667;
t683 = t805 * t694 + t800 * t695;
t815 = m(5) * t713 - t724 * mrSges(5,1) + t725 * mrSges(5,2) - t757 * t744 + t758 * t745 + t683;
t812 = -m(4) * t739 + t752 * mrSges(4,1) - t753 * mrSges(4,2) + t771 * t760 - t772 * t761 - t815;
t678 = m(3) * t754 + t793 * mrSges(3,1) - t783 * mrSges(3,3) + t794 * t778 - t781 * t823 + t812;
t653 = t664 * t830 - t798 * t666 + t678 * t829;
t650 = m(2) * t789 + qJDD(1) * mrSges(2,1) - t810 * mrSges(2,2) + t653;
t660 = t808 * t664 - t803 * t678;
t658 = m(2) * t790 - t810 * mrSges(2,1) - qJDD(1) * mrSges(2,2) + t660;
t828 = t809 * t650 + t804 * t658;
t652 = t664 * t832 + t799 * t666 + t678 * t831;
t821 = -t804 * t650 + t809 * t658;
t716 = Ifges(6,5) * t743 + Ifges(6,6) * t742 + Ifges(6,3) * t756;
t718 = Ifges(6,1) * t743 + Ifges(6,4) * t742 + Ifges(6,5) * t756;
t687 = -mrSges(6,1) * t700 + mrSges(6,3) * t699 + Ifges(6,4) * t712 + Ifges(6,2) * t711 + Ifges(6,6) * t723 - t743 * t716 + t756 * t718;
t717 = Ifges(6,4) * t743 + Ifges(6,2) * t742 + Ifges(6,6) * t756;
t688 = mrSges(6,2) * t700 - mrSges(6,3) * t698 + Ifges(6,1) * t712 + Ifges(6,4) * t711 + Ifges(6,5) * t723 + t742 * t716 - t756 * t717;
t730 = Ifges(5,5) * t758 + Ifges(5,6) * t757 + Ifges(5,3) * t786;
t731 = Ifges(5,4) * t758 + Ifges(5,2) * t757 + Ifges(5,6) * t786;
t668 = mrSges(5,2) * t713 - mrSges(5,3) * t704 + Ifges(5,1) * t725 + Ifges(5,4) * t724 + Ifges(5,5) * t775 - pkin(10) * t683 - t800 * t687 + t805 * t688 + t757 * t730 - t786 * t731;
t732 = Ifges(5,1) * t758 + Ifges(5,4) * t757 + Ifges(5,5) * t786;
t813 = mrSges(6,1) * t698 - mrSges(6,2) * t699 + Ifges(6,5) * t712 + Ifges(6,6) * t711 + Ifges(6,3) * t723 + t743 * t717 - t742 * t718;
t669 = -mrSges(5,1) * t713 + mrSges(5,3) * t705 + Ifges(5,4) * t725 + Ifges(5,2) * t724 + Ifges(5,6) * t775 - pkin(4) * t683 - t758 * t730 + t786 * t732 - t813;
t746 = Ifges(4,5) * t772 + Ifges(4,6) * t771 + Ifges(4,3) * t788;
t748 = Ifges(4,1) * t772 + Ifges(4,4) * t771 + Ifges(4,5) * t788;
t654 = -mrSges(4,1) * t739 + mrSges(4,3) * t715 + Ifges(4,4) * t753 + Ifges(4,2) * t752 + Ifges(4,6) * t776 - pkin(3) * t815 + pkin(9) * t819 + t801 * t668 + t806 * t669 - t772 * t746 + t788 * t748;
t747 = Ifges(4,4) * t772 + Ifges(4,2) * t771 + Ifges(4,6) * t788;
t655 = mrSges(4,2) * t739 - mrSges(4,3) * t714 + Ifges(4,1) * t753 + Ifges(4,4) * t752 + Ifges(4,5) * t776 - pkin(9) * t675 + t806 * t668 - t801 * t669 + t771 * t746 - t788 * t747;
t764 = Ifges(3,6) * t794 + (Ifges(3,4) * t803 + Ifges(3,2) * t808) * t826;
t765 = Ifges(3,5) * t794 + (Ifges(3,1) * t803 + Ifges(3,4) * t808) * t826;
t644 = Ifges(3,5) * t783 + Ifges(3,6) * t784 + Ifges(3,3) * t793 + mrSges(3,1) * t754 - mrSges(3,2) * t755 + t802 * t655 + t807 * t654 + pkin(2) * t812 + pkin(8) * t820 + (t764 * t803 - t765 * t808) * t826;
t763 = Ifges(3,3) * t794 + (Ifges(3,5) * t803 + Ifges(3,6) * t808) * t826;
t646 = mrSges(3,2) * t766 - mrSges(3,3) * t754 + Ifges(3,1) * t783 + Ifges(3,4) * t784 + Ifges(3,5) * t793 - pkin(8) * t667 - t802 * t654 + t807 * t655 + t763 * t822 - t794 * t764;
t814 = -mrSges(5,1) * t704 + mrSges(5,2) * t705 - Ifges(5,5) * t725 - Ifges(5,6) * t724 - Ifges(5,3) * t775 - pkin(4) * t817 - pkin(10) * t818 - t805 * t687 - t800 * t688 - t758 * t731 + t757 * t732;
t811 = mrSges(4,1) * t714 - mrSges(4,2) * t715 + Ifges(4,5) * t753 + Ifges(4,6) * t752 + Ifges(4,3) * t776 + pkin(3) * t675 + t772 * t747 - t771 * t748 - t814;
t648 = -mrSges(3,1) * t766 + mrSges(3,3) * t755 + Ifges(3,4) * t783 + Ifges(3,2) * t784 + Ifges(3,6) * t793 - pkin(2) * t667 - t763 * t823 + t794 * t765 - t811;
t816 = mrSges(2,1) * t789 - mrSges(2,2) * t790 + Ifges(2,3) * qJDD(1) + pkin(1) * t653 + t799 * t644 + t646 * t832 + t648 * t831 + t660 * t834;
t642 = -mrSges(2,2) * g(3) - mrSges(2,3) * t789 + Ifges(2,5) * qJDD(1) - t810 * Ifges(2,6) + t808 * t646 - t803 * t648 + (-t652 * t798 - t653 * t799) * pkin(7);
t641 = mrSges(2,1) * g(3) + mrSges(2,3) * t790 + t810 * Ifges(2,5) + Ifges(2,6) * qJDD(1) - pkin(1) * t652 - t798 * t644 + (pkin(7) * t660 + t646 * t803 + t648 * t808) * t799;
t1 = [-m(1) * g(1) + t821; -m(1) * g(2) + t828; (-m(1) - m(2)) * g(3) + t652; -mrSges(1,2) * g(3) + mrSges(1,3) * g(2) - pkin(6) * t828 - t804 * t641 + t809 * t642; mrSges(1,1) * g(3) - mrSges(1,3) * g(1) + pkin(6) * t821 + t809 * t641 + t804 * t642; -mrSges(1,1) * g(2) + mrSges(1,2) * g(1) + t816; t816; t644; t811; -t814; t813;];
tauJB = t1;
