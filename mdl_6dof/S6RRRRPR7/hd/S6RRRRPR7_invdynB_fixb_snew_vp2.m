% Calculate vector of inverse dynamics base forces with Newton-Euler for
% S6RRRRPR7
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
%   pkin=[a2,a3,a4,a5,a6,alpha2,d1,d2,d3,d4,d6,theta5]';
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
% tauB [6x1]
%   base forces of inverse dynamics (contains inertial, gravitational coriolis and centrifugal forces)

% Quelle: HybrDyn-Toolbox
% Datum: 2019-05-07 21:26
% Revision: 89c353f7eff3bd693eda4e29f35b2761dbc3ada0 (2019-05-03)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauB = S6RRRRPR7_invdynB_fixb_snew_vp2(qJ, qJD, qJDD, g, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(6,1),zeros(3,1),zeros(12,1),zeros(7,1),zeros(7,3),zeros(7,6)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRRRPR7_invdynB_fixb_snew_vp2: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RRRRPR7_invdynB_fixb_snew_vp2: qJD has to be [6x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [6 1]), ...
  'S6RRRRPR7_invdynB_fixb_snew_vp2: qJDD has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6RRRRPR7_invdynB_fixb_snew_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [12 1]), ...
  'S6RRRRPR7_invdynB_fixb_snew_vp2: pkin has to be [12x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6RRRRPR7_invdynB_fixb_snew_vp2: m has to be [7x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [7,3]), ...
  'S6RRRRPR7_invdynB_fixb_snew_vp2: mrSges has to be [7x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [7 6]), ...
  'S6RRRRPR7_invdynB_fixb_snew_vp2: Ifges has to be [7x6] (double)'); 

%% Symbolic Calculation
% From invdyn_fixb_NewtonEuler_linkframe_tauB_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-05-07 21:13:41
% EndTime: 2019-05-07 21:15:06
% DurationCPUTime: 62.63s
% Computational Cost: add. (1043554->399), mult. (2287071->516), div. (0->0), fcn. (1875061->14), ass. (0->164)
t805 = 2 * qJD(5);
t770 = cos(pkin(6));
t804 = t770 * g(3);
t768 = sin(pkin(6));
t774 = sin(qJ(2));
t803 = t768 * t774;
t779 = cos(qJ(2));
t802 = t768 * t779;
t801 = t770 * t774;
t800 = t770 * t779;
t775 = sin(qJ(1));
t780 = cos(qJ(1));
t759 = t775 * g(1) - g(2) * t780;
t781 = qJD(1) ^ 2;
t750 = pkin(8) * t768 * t781 + qJDD(1) * pkin(1) + t759;
t760 = -g(1) * t780 - g(2) * t775;
t795 = qJDD(1) * t768;
t751 = -pkin(1) * t781 + pkin(8) * t795 + t760;
t798 = t750 * t801 + t779 * t751;
t727 = -g(3) * t803 + t798;
t764 = qJD(1) * t770 + qJD(2);
t797 = qJD(1) * t768;
t794 = t774 * t797;
t748 = mrSges(3,1) * t764 - mrSges(3,3) * t794;
t752 = (-mrSges(3,1) * t779 + mrSges(3,2) * t774) * t797;
t755 = -qJD(2) * t794 + t779 * t795;
t763 = qJDD(1) * t770 + qJDD(2);
t753 = (-pkin(2) * t779 - pkin(9) * t774) * t797;
t762 = t764 ^ 2;
t796 = qJD(1) * t779;
t713 = -t762 * pkin(2) + t763 * pkin(9) + (-g(3) * t774 + t753 * t796) * t768 + t798;
t754 = (qJD(2) * t796 + qJDD(1) * t774) * t768;
t714 = -t755 * pkin(2) - t754 * pkin(9) - t804 + (-t750 + (pkin(2) * t774 - pkin(9) * t779) * t764 * qJD(1)) * t768;
t773 = sin(qJ(3));
t778 = cos(qJ(3));
t690 = -t773 * t713 + t778 * t714;
t742 = t764 * t778 - t773 * t794;
t725 = qJD(3) * t742 + t754 * t778 + t763 * t773;
t743 = t764 * t773 + t778 * t794;
t747 = qJDD(3) - t755;
t793 = t768 * t796;
t758 = qJD(3) - t793;
t671 = (t742 * t758 - t725) * pkin(10) + (t742 * t743 + t747) * pkin(3) + t690;
t691 = t778 * t713 + t773 * t714;
t724 = -qJD(3) * t743 - t754 * t773 + t763 * t778;
t734 = pkin(3) * t758 - pkin(10) * t743;
t741 = t742 ^ 2;
t678 = -pkin(3) * t741 + pkin(10) * t724 - t734 * t758 + t691;
t772 = sin(qJ(4));
t777 = cos(qJ(4));
t658 = t777 * t671 - t772 * t678;
t729 = t742 * t777 - t743 * t772;
t694 = qJD(4) * t729 + t724 * t772 + t725 * t777;
t730 = t742 * t772 + t743 * t777;
t746 = qJDD(4) + t747;
t757 = qJD(4) + t758;
t655 = (t729 * t757 - t694) * qJ(5) + (t729 * t730 + t746) * pkin(4) + t658;
t659 = t772 * t671 + t777 * t678;
t693 = -qJD(4) * t730 + t724 * t777 - t725 * t772;
t716 = pkin(4) * t757 - qJ(5) * t730;
t728 = t729 ^ 2;
t657 = -pkin(4) * t728 + qJ(5) * t693 - t716 * t757 + t659;
t767 = sin(pkin(12));
t769 = cos(pkin(12));
t706 = t729 * t769 - t730 * t767;
t652 = t767 * t655 + t769 * t657 + t706 * t805;
t675 = t693 * t769 - t694 * t767;
t707 = t729 * t767 + t730 * t769;
t688 = -mrSges(6,1) * t706 + mrSges(6,2) * t707;
t698 = mrSges(6,1) * t757 - mrSges(6,3) * t707;
t689 = -pkin(5) * t706 - pkin(11) * t707;
t756 = t757 ^ 2;
t650 = -pkin(5) * t756 + pkin(11) * t746 + t689 * t706 + t652;
t726 = -g(3) * t802 + t750 * t800 - t774 * t751;
t712 = -t763 * pkin(2) - t762 * pkin(9) + t753 * t794 - t726;
t682 = -t724 * pkin(3) - t741 * pkin(10) + t743 * t734 + t712;
t661 = -t693 * pkin(4) - t728 * qJ(5) + t730 * t716 + qJDD(5) + t682;
t676 = t693 * t767 + t694 * t769;
t653 = t661 + (-t706 * t757 - t676) * pkin(11) + (t707 * t757 - t675) * pkin(5);
t771 = sin(qJ(6));
t776 = cos(qJ(6));
t647 = -t650 * t771 + t653 * t776;
t695 = -t707 * t771 + t757 * t776;
t664 = qJD(6) * t695 + t676 * t776 + t746 * t771;
t674 = qJDD(6) - t675;
t696 = t707 * t776 + t757 * t771;
t679 = -mrSges(7,1) * t695 + mrSges(7,2) * t696;
t705 = qJD(6) - t706;
t680 = -mrSges(7,2) * t705 + mrSges(7,3) * t695;
t645 = m(7) * t647 + mrSges(7,1) * t674 - mrSges(7,3) * t664 - t679 * t696 + t680 * t705;
t648 = t650 * t776 + t653 * t771;
t663 = -qJD(6) * t696 - t676 * t771 + t746 * t776;
t681 = mrSges(7,1) * t705 - mrSges(7,3) * t696;
t646 = m(7) * t648 - mrSges(7,2) * t674 + mrSges(7,3) * t663 + t679 * t695 - t681 * t705;
t788 = -t645 * t771 + t776 * t646;
t636 = m(6) * t652 - mrSges(6,2) * t746 + mrSges(6,3) * t675 + t688 * t706 - t698 * t757 + t788;
t787 = -t769 * t655 + t767 * t657;
t651 = -0.2e1 * qJD(5) * t707 - t787;
t697 = -mrSges(6,2) * t757 + mrSges(6,3) * t706;
t649 = -t746 * pkin(5) - t756 * pkin(11) + (t805 + t689) * t707 + t787;
t784 = -m(7) * t649 + t663 * mrSges(7,1) - mrSges(7,2) * t664 + t695 * t680 - t681 * t696;
t641 = m(6) * t651 + mrSges(6,1) * t746 - mrSges(6,3) * t676 - t688 * t707 + t697 * t757 + t784;
t630 = t767 * t636 + t769 * t641;
t708 = -mrSges(5,1) * t729 + mrSges(5,2) * t730;
t715 = -mrSges(5,2) * t757 + mrSges(5,3) * t729;
t628 = m(5) * t658 + mrSges(5,1) * t746 - mrSges(5,3) * t694 - t708 * t730 + t715 * t757 + t630;
t717 = mrSges(5,1) * t757 - mrSges(5,3) * t730;
t789 = t769 * t636 - t641 * t767;
t629 = m(5) * t659 - mrSges(5,2) * t746 + mrSges(5,3) * t693 + t708 * t729 - t717 * t757 + t789;
t622 = t777 * t628 + t772 * t629;
t731 = -mrSges(4,1) * t742 + mrSges(4,2) * t743;
t732 = -mrSges(4,2) * t758 + mrSges(4,3) * t742;
t620 = m(4) * t690 + mrSges(4,1) * t747 - mrSges(4,3) * t725 - t731 * t743 + t732 * t758 + t622;
t733 = mrSges(4,1) * t758 - mrSges(4,3) * t743;
t790 = -t628 * t772 + t777 * t629;
t621 = m(4) * t691 - mrSges(4,2) * t747 + mrSges(4,3) * t724 + t731 * t742 - t733 * t758 + t790;
t791 = -t620 * t773 + t778 * t621;
t612 = m(3) * t727 - mrSges(3,2) * t763 + mrSges(3,3) * t755 - t748 * t764 + t752 * t793 + t791;
t615 = t778 * t620 + t773 * t621;
t738 = -t768 * t750 - t804;
t749 = -mrSges(3,2) * t764 + mrSges(3,3) * t793;
t614 = m(3) * t738 - t755 * mrSges(3,1) + t754 * mrSges(3,2) + (t748 * t774 - t749 * t779) * t797 + t615;
t637 = t776 * t645 + t771 * t646;
t786 = m(6) * t661 - t675 * mrSges(6,1) + t676 * mrSges(6,2) - t706 * t697 + t707 * t698 + t637;
t783 = m(5) * t682 - t693 * mrSges(5,1) + t694 * mrSges(5,2) - t729 * t715 + t730 * t717 + t786;
t782 = -m(4) * t712 + t724 * mrSges(4,1) - t725 * mrSges(4,2) + t742 * t732 - t743 * t733 - t783;
t633 = m(3) * t726 + t763 * mrSges(3,1) - t754 * mrSges(3,3) + t764 * t749 - t752 * t794 + t782;
t603 = t612 * t801 - t614 * t768 + t633 * t800;
t601 = m(2) * t759 + qJDD(1) * mrSges(2,1) - mrSges(2,2) * t781 + t603;
t607 = t779 * t612 - t633 * t774;
t606 = m(2) * t760 - mrSges(2,1) * t781 - qJDD(1) * mrSges(2,2) + t607;
t799 = t780 * t601 + t775 * t606;
t602 = t612 * t803 + t770 * t614 + t633 * t802;
t792 = -t601 * t775 + t780 * t606;
t665 = Ifges(7,5) * t696 + Ifges(7,6) * t695 + Ifges(7,3) * t705;
t667 = Ifges(7,1) * t696 + Ifges(7,4) * t695 + Ifges(7,5) * t705;
t638 = -mrSges(7,1) * t649 + mrSges(7,3) * t648 + Ifges(7,4) * t664 + Ifges(7,2) * t663 + Ifges(7,6) * t674 - t665 * t696 + t667 * t705;
t666 = Ifges(7,4) * t696 + Ifges(7,2) * t695 + Ifges(7,6) * t705;
t639 = mrSges(7,2) * t649 - mrSges(7,3) * t647 + Ifges(7,1) * t664 + Ifges(7,4) * t663 + Ifges(7,5) * t674 + t665 * t695 - t666 * t705;
t683 = Ifges(6,5) * t707 + Ifges(6,6) * t706 + Ifges(6,3) * t757;
t684 = Ifges(6,4) * t707 + Ifges(6,2) * t706 + Ifges(6,6) * t757;
t623 = mrSges(6,2) * t661 - mrSges(6,3) * t651 + Ifges(6,1) * t676 + Ifges(6,4) * t675 + Ifges(6,5) * t746 - pkin(11) * t637 - t638 * t771 + t639 * t776 + t683 * t706 - t684 * t757;
t685 = Ifges(6,1) * t707 + Ifges(6,4) * t706 + Ifges(6,5) * t757;
t624 = -mrSges(6,1) * t661 - mrSges(7,1) * t647 + mrSges(7,2) * t648 + mrSges(6,3) * t652 + Ifges(6,4) * t676 - Ifges(7,5) * t664 + Ifges(6,2) * t675 + Ifges(6,6) * t746 - Ifges(7,6) * t663 - Ifges(7,3) * t674 - pkin(5) * t637 - t666 * t696 + t667 * t695 - t683 * t707 + t685 * t757;
t699 = Ifges(5,5) * t730 + Ifges(5,6) * t729 + Ifges(5,3) * t757;
t701 = Ifges(5,1) * t730 + Ifges(5,4) * t729 + Ifges(5,5) * t757;
t608 = -mrSges(5,1) * t682 + mrSges(5,3) * t659 + Ifges(5,4) * t694 + Ifges(5,2) * t693 + Ifges(5,6) * t746 - pkin(4) * t786 + qJ(5) * t789 + t767 * t623 + t769 * t624 - t730 * t699 + t757 * t701;
t700 = Ifges(5,4) * t730 + Ifges(5,2) * t729 + Ifges(5,6) * t757;
t616 = mrSges(5,2) * t682 - mrSges(5,3) * t658 + Ifges(5,1) * t694 + Ifges(5,4) * t693 + Ifges(5,5) * t746 - qJ(5) * t630 + t623 * t769 - t624 * t767 + t699 * t729 - t700 * t757;
t718 = Ifges(4,5) * t743 + Ifges(4,6) * t742 + Ifges(4,3) * t758;
t720 = Ifges(4,1) * t743 + Ifges(4,4) * t742 + Ifges(4,5) * t758;
t597 = -mrSges(4,1) * t712 + mrSges(4,3) * t691 + Ifges(4,4) * t725 + Ifges(4,2) * t724 + Ifges(4,6) * t747 - pkin(3) * t783 + pkin(10) * t790 + t777 * t608 + t772 * t616 - t743 * t718 + t758 * t720;
t719 = Ifges(4,4) * t743 + Ifges(4,2) * t742 + Ifges(4,6) * t758;
t598 = mrSges(4,2) * t712 - mrSges(4,3) * t690 + Ifges(4,1) * t725 + Ifges(4,4) * t724 + Ifges(4,5) * t747 - pkin(10) * t622 - t608 * t772 + t616 * t777 + t718 * t742 - t719 * t758;
t735 = Ifges(3,3) * t764 + (Ifges(3,5) * t774 + Ifges(3,6) * t779) * t797;
t736 = Ifges(3,6) * t764 + (Ifges(3,4) * t774 + Ifges(3,2) * t779) * t797;
t596 = mrSges(3,2) * t738 - mrSges(3,3) * t726 + Ifges(3,1) * t754 + Ifges(3,4) * t755 + Ifges(3,5) * t763 - pkin(9) * t615 - t597 * t773 + t598 * t778 + t735 * t793 - t736 * t764;
t737 = Ifges(3,5) * t764 + (Ifges(3,1) * t774 + Ifges(3,4) * t779) * t797;
t599 = Ifges(3,6) * t763 + t764 * t737 - pkin(11) * t788 - t735 * t794 - Ifges(4,3) * t747 + Ifges(3,4) * t754 + Ifges(3,2) * t755 - mrSges(3,1) * t738 + t742 * t720 - t743 * t719 + t729 * t701 - t730 * t700 - Ifges(4,6) * t724 - Ifges(4,5) * t725 + mrSges(3,3) * t727 + t706 * t685 - t707 * t684 - Ifges(5,6) * t693 - Ifges(5,5) * t694 + mrSges(4,2) * t691 - mrSges(4,1) * t690 - Ifges(6,5) * t676 - Ifges(6,6) * t675 - mrSges(5,1) * t658 + mrSges(5,2) * t659 - mrSges(6,1) * t651 + mrSges(6,2) * t652 - t771 * t639 - pkin(2) * t615 - pkin(4) * t630 + (-Ifges(6,3) - Ifges(5,3)) * t746 - pkin(3) * t622 - pkin(5) * t784 - t776 * t638;
t785 = pkin(8) * t607 + t596 * t774 + t599 * t779;
t595 = Ifges(3,5) * t754 + Ifges(3,6) * t755 + Ifges(3,3) * t763 + mrSges(3,1) * t726 - mrSges(3,2) * t727 + t773 * t598 + t778 * t597 + pkin(2) * t782 + pkin(9) * t791 + (t736 * t774 - t737 * t779) * t797;
t594 = -mrSges(2,2) * g(3) - mrSges(2,3) * t759 + Ifges(2,5) * qJDD(1) - t781 * Ifges(2,6) + t779 * t596 - t774 * t599 + (-t602 * t768 - t603 * t770) * pkin(8);
t593 = mrSges(2,1) * g(3) + mrSges(2,3) * t760 + t781 * Ifges(2,5) + Ifges(2,6) * qJDD(1) - pkin(1) * t602 - t768 * t595 + t770 * t785;
t1 = [-m(1) * g(1) + t792; -m(1) * g(2) + t799; (-m(1) - m(2)) * g(3) + t602; -mrSges(1,2) * g(3) + mrSges(1,3) * g(2) - pkin(7) * t799 - t775 * t593 + t780 * t594; mrSges(1,1) * g(3) - mrSges(1,3) * g(1) + pkin(7) * t792 + t780 * t593 + t775 * t594; -mrSges(1,1) * g(2) + mrSges(2,1) * t759 + mrSges(1,2) * g(1) - mrSges(2,2) * t760 + Ifges(2,3) * qJDD(1) + pkin(1) * t603 + t770 * t595 + t768 * t785;];
tauB  = t1;
