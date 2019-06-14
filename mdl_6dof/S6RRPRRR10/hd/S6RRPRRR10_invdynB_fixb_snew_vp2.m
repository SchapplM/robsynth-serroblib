% Calculate vector of inverse dynamics base forces with Newton-Euler for
% S6RRPRRR10
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
%   pkin=[a2,a3,a4,a5,a6,alpha2,d1,d2,d4,d5,d6,theta3]';
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
% Datum: 2019-05-06 23:51
% Revision: 89c353f7eff3bd693eda4e29f35b2761dbc3ada0 (2019-05-03)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauB = S6RRPRRR10_invdynB_fixb_snew_vp2(qJ, qJD, qJDD, g, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(6,1),zeros(3,1),zeros(12,1),zeros(7,1),zeros(7,3),zeros(7,6)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRPRRR10_invdynB_fixb_snew_vp2: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RRPRRR10_invdynB_fixb_snew_vp2: qJD has to be [6x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [6 1]), ...
  'S6RRPRRR10_invdynB_fixb_snew_vp2: qJDD has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6RRPRRR10_invdynB_fixb_snew_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [12 1]), ...
  'S6RRPRRR10_invdynB_fixb_snew_vp2: pkin has to be [12x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6RRPRRR10_invdynB_fixb_snew_vp2: m has to be [7x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [7,3]), ...
  'S6RRPRRR10_invdynB_fixb_snew_vp2: mrSges has to be [7x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [7 6]), ...
  'S6RRPRRR10_invdynB_fixb_snew_vp2: Ifges has to be [7x6] (double)'); 

%% Symbolic Calculation
% From invdyn_fixb_NewtonEuler_linkframe_tauB_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-05-06 23:39:50
% EndTime: 2019-05-06 23:40:56
% DurationCPUTime: 53.21s
% Computational Cost: add. (859766->398), mult. (1943252->516), div. (0->0), fcn. (1588952->14), ass. (0->162)
t762 = cos(pkin(6));
t795 = t762 * g(3);
t760 = sin(pkin(6));
t766 = sin(qJ(2));
t794 = t760 * t766;
t771 = cos(qJ(2));
t793 = t760 * t771;
t792 = t762 * t766;
t791 = t762 * t771;
t767 = sin(qJ(1));
t772 = cos(qJ(1));
t751 = t767 * g(1) - g(2) * t772;
t773 = qJD(1) ^ 2;
t743 = pkin(8) * t760 * t773 + qJDD(1) * pkin(1) + t751;
t752 = -g(1) * t772 - g(2) * t767;
t786 = qJDD(1) * t760;
t744 = -pkin(1) * t773 + pkin(8) * t786 + t752;
t789 = t743 * t792 + t771 * t744;
t713 = -g(3) * t794 + t789;
t756 = qJD(1) * t762 + qJD(2);
t788 = qJD(1) * t760;
t785 = t766 * t788;
t741 = mrSges(3,1) * t756 - mrSges(3,3) * t785;
t746 = (-mrSges(3,1) * t771 + mrSges(3,2) * t766) * t788;
t748 = -qJD(2) * t785 + t771 * t786;
t755 = qJDD(1) * t762 + qJDD(2);
t745 = (-pkin(2) * t771 - qJ(3) * t766) * t788;
t754 = t756 ^ 2;
t787 = qJD(1) * t771;
t699 = -t754 * pkin(2) + t755 * qJ(3) + (-g(3) * t766 + t745 * t787) * t760 + t789;
t747 = (qJD(2) * t787 + qJDD(1) * t766) * t760;
t700 = -t748 * pkin(2) - t795 - t747 * qJ(3) + (-t743 + (pkin(2) * t766 - qJ(3) * t771) * t756 * qJD(1)) * t760;
t759 = sin(pkin(12));
t761 = cos(pkin(12));
t737 = t756 * t759 + t761 * t785;
t665 = -0.2e1 * qJD(3) * t737 - t759 * t699 + t761 * t700;
t724 = t747 * t761 + t755 * t759;
t736 = t756 * t761 - t759 * t785;
t784 = t760 * t787;
t656 = (-t736 * t784 - t724) * pkin(9) + (t736 * t737 - t748) * pkin(3) + t665;
t666 = 0.2e1 * qJD(3) * t736 + t761 * t699 + t759 * t700;
t723 = -t747 * t759 + t755 * t761;
t725 = -pkin(3) * t784 - pkin(9) * t737;
t735 = t736 ^ 2;
t663 = -pkin(3) * t735 + pkin(9) * t723 + t725 * t784 + t666;
t765 = sin(qJ(4));
t770 = cos(qJ(4));
t648 = t765 * t656 + t770 * t663;
t717 = t736 * t765 + t737 * t770;
t684 = -t717 * qJD(4) + t723 * t770 - t765 * t724;
t716 = t736 * t770 - t765 * t737;
t693 = -mrSges(5,1) * t716 + mrSges(5,2) * t717;
t750 = qJD(4) - t784;
t705 = mrSges(5,1) * t750 - mrSges(5,3) * t717;
t740 = qJDD(4) - t748;
t694 = -pkin(4) * t716 - pkin(10) * t717;
t749 = t750 ^ 2;
t646 = -pkin(4) * t749 + pkin(10) * t740 + t694 * t716 + t648;
t712 = -g(3) * t793 + t743 * t791 - t766 * t744;
t698 = -t755 * pkin(2) - t754 * qJ(3) + t745 * t785 + qJDD(3) - t712;
t670 = -t723 * pkin(3) - t735 * pkin(9) + t737 * t725 + t698;
t685 = qJD(4) * t716 + t723 * t765 + t724 * t770;
t654 = (-t716 * t750 - t685) * pkin(10) + (t717 * t750 - t684) * pkin(4) + t670;
t764 = sin(qJ(5));
t769 = cos(qJ(5));
t641 = -t764 * t646 + t769 * t654;
t702 = -t717 * t764 + t750 * t769;
t669 = qJD(5) * t702 + t685 * t769 + t740 * t764;
t683 = qJDD(5) - t684;
t703 = t717 * t769 + t750 * t764;
t715 = qJD(5) - t716;
t639 = (t702 * t715 - t669) * pkin(11) + (t702 * t703 + t683) * pkin(5) + t641;
t642 = t769 * t646 + t764 * t654;
t668 = -qJD(5) * t703 - t685 * t764 + t740 * t769;
t688 = pkin(5) * t715 - pkin(11) * t703;
t701 = t702 ^ 2;
t640 = -pkin(5) * t701 + pkin(11) * t668 - t688 * t715 + t642;
t763 = sin(qJ(6));
t768 = cos(qJ(6));
t637 = t639 * t768 - t640 * t763;
t678 = t702 * t768 - t703 * t763;
t651 = qJD(6) * t678 + t668 * t763 + t669 * t768;
t679 = t702 * t763 + t703 * t768;
t664 = -mrSges(7,1) * t678 + mrSges(7,2) * t679;
t711 = qJD(6) + t715;
t671 = -mrSges(7,2) * t711 + mrSges(7,3) * t678;
t681 = qJDD(6) + t683;
t635 = m(7) * t637 + mrSges(7,1) * t681 - mrSges(7,3) * t651 - t664 * t679 + t671 * t711;
t638 = t639 * t763 + t640 * t768;
t650 = -qJD(6) * t679 + t668 * t768 - t669 * t763;
t672 = mrSges(7,1) * t711 - mrSges(7,3) * t679;
t636 = m(7) * t638 - mrSges(7,2) * t681 + mrSges(7,3) * t650 + t664 * t678 - t672 * t711;
t627 = t768 * t635 + t763 * t636;
t680 = -mrSges(6,1) * t702 + mrSges(6,2) * t703;
t686 = -mrSges(6,2) * t715 + mrSges(6,3) * t702;
t625 = m(6) * t641 + mrSges(6,1) * t683 - mrSges(6,3) * t669 - t680 * t703 + t686 * t715 + t627;
t687 = mrSges(6,1) * t715 - mrSges(6,3) * t703;
t779 = -t635 * t763 + t768 * t636;
t626 = m(6) * t642 - mrSges(6,2) * t683 + mrSges(6,3) * t668 + t680 * t702 - t687 * t715 + t779;
t780 = -t625 * t764 + t769 * t626;
t620 = m(5) * t648 - mrSges(5,2) * t740 + mrSges(5,3) * t684 + t693 * t716 - t705 * t750 + t780;
t647 = t656 * t770 - t765 * t663;
t704 = -mrSges(5,2) * t750 + mrSges(5,3) * t716;
t645 = -pkin(4) * t740 - pkin(10) * t749 + t717 * t694 - t647;
t643 = -pkin(5) * t668 - pkin(11) * t701 + t688 * t703 + t645;
t777 = m(7) * t643 - t650 * mrSges(7,1) + mrSges(7,2) * t651 - t678 * t671 + t672 * t679;
t775 = -m(6) * t645 + t668 * mrSges(6,1) - mrSges(6,2) * t669 + t702 * t686 - t687 * t703 - t777;
t631 = m(5) * t647 + mrSges(5,1) * t740 - mrSges(5,3) * t685 - t693 * t717 + t704 * t750 + t775;
t612 = t765 * t620 + t770 * t631;
t718 = -mrSges(4,1) * t736 + mrSges(4,2) * t737;
t721 = mrSges(4,2) * t784 + mrSges(4,3) * t736;
t610 = m(4) * t665 - mrSges(4,1) * t748 - mrSges(4,3) * t724 - t718 * t737 - t721 * t784 + t612;
t722 = -mrSges(4,1) * t784 - mrSges(4,3) * t737;
t781 = t770 * t620 - t631 * t765;
t611 = m(4) * t666 + mrSges(4,2) * t748 + mrSges(4,3) * t723 + t718 * t736 + t722 * t784 + t781;
t782 = -t610 * t759 + t761 * t611;
t602 = m(3) * t713 - mrSges(3,2) * t755 + mrSges(3,3) * t748 - t741 * t756 + t746 * t784 + t782;
t605 = t761 * t610 + t759 * t611;
t729 = -t760 * t743 - t795;
t742 = -mrSges(3,2) * t756 + mrSges(3,3) * t784;
t604 = m(3) * t729 - t748 * mrSges(3,1) + t747 * mrSges(3,2) + (t741 * t766 - t742 * t771) * t788 + t605;
t621 = t769 * t625 + t764 * t626;
t776 = m(5) * t670 - t684 * mrSges(5,1) + mrSges(5,2) * t685 - t716 * t704 + t705 * t717 + t621;
t774 = -m(4) * t698 + t723 * mrSges(4,1) - mrSges(4,2) * t724 + t736 * t721 - t722 * t737 - t776;
t617 = m(3) * t712 + mrSges(3,1) * t755 - mrSges(3,3) * t747 + t742 * t756 - t746 * t785 + t774;
t593 = t602 * t792 - t604 * t760 + t617 * t791;
t591 = m(2) * t751 + qJDD(1) * mrSges(2,1) - mrSges(2,2) * t773 + t593;
t597 = t771 * t602 - t617 * t766;
t596 = m(2) * t752 - mrSges(2,1) * t773 - qJDD(1) * mrSges(2,2) + t597;
t790 = t772 * t591 + t767 * t596;
t592 = t602 * t794 + t762 * t604 + t617 * t793;
t783 = -t591 * t767 + t772 * t596;
t657 = Ifges(7,5) * t679 + Ifges(7,6) * t678 + Ifges(7,3) * t711;
t659 = Ifges(7,1) * t679 + Ifges(7,4) * t678 + Ifges(7,5) * t711;
t628 = -mrSges(7,1) * t643 + mrSges(7,3) * t638 + Ifges(7,4) * t651 + Ifges(7,2) * t650 + Ifges(7,6) * t681 - t657 * t679 + t659 * t711;
t658 = Ifges(7,4) * t679 + Ifges(7,2) * t678 + Ifges(7,6) * t711;
t629 = mrSges(7,2) * t643 - mrSges(7,3) * t637 + Ifges(7,1) * t651 + Ifges(7,4) * t650 + Ifges(7,5) * t681 + t657 * t678 - t658 * t711;
t673 = Ifges(6,5) * t703 + Ifges(6,6) * t702 + Ifges(6,3) * t715;
t675 = Ifges(6,1) * t703 + Ifges(6,4) * t702 + Ifges(6,5) * t715;
t613 = -mrSges(6,1) * t645 + mrSges(6,3) * t642 + Ifges(6,4) * t669 + Ifges(6,2) * t668 + Ifges(6,6) * t683 - pkin(5) * t777 + pkin(11) * t779 + t768 * t628 + t763 * t629 - t703 * t673 + t715 * t675;
t674 = Ifges(6,4) * t703 + Ifges(6,2) * t702 + Ifges(6,6) * t715;
t614 = mrSges(6,2) * t645 - mrSges(6,3) * t641 + Ifges(6,1) * t669 + Ifges(6,4) * t668 + Ifges(6,5) * t683 - pkin(11) * t627 - t628 * t763 + t629 * t768 + t673 * t702 - t674 * t715;
t689 = Ifges(5,5) * t717 + Ifges(5,6) * t716 + Ifges(5,3) * t750;
t690 = Ifges(5,4) * t717 + Ifges(5,2) * t716 + Ifges(5,6) * t750;
t598 = mrSges(5,2) * t670 - mrSges(5,3) * t647 + Ifges(5,1) * t685 + Ifges(5,4) * t684 + Ifges(5,5) * t740 - pkin(10) * t621 - t613 * t764 + t614 * t769 + t689 * t716 - t690 * t750;
t691 = Ifges(5,1) * t717 + Ifges(5,4) * t716 + Ifges(5,5) * t750;
t606 = Ifges(5,4) * t685 + Ifges(5,2) * t684 + Ifges(5,6) * t740 - t717 * t689 + t750 * t691 - mrSges(5,1) * t670 + mrSges(5,3) * t648 - Ifges(6,5) * t669 - Ifges(6,6) * t668 - Ifges(6,3) * t683 - t703 * t674 + t702 * t675 - mrSges(6,1) * t641 + mrSges(6,2) * t642 - Ifges(7,5) * t651 - Ifges(7,6) * t650 - Ifges(7,3) * t681 - t679 * t658 + t678 * t659 - mrSges(7,1) * t637 + mrSges(7,2) * t638 - pkin(5) * t627 - pkin(4) * t621;
t706 = Ifges(4,5) * t737 + Ifges(4,6) * t736 - Ifges(4,3) * t784;
t708 = Ifges(4,1) * t737 + Ifges(4,4) * t736 - Ifges(4,5) * t784;
t587 = -mrSges(4,1) * t698 + mrSges(4,3) * t666 + Ifges(4,4) * t724 + Ifges(4,2) * t723 - Ifges(4,6) * t748 - pkin(3) * t776 + pkin(9) * t781 + t765 * t598 + t770 * t606 - t737 * t706 - t708 * t784;
t707 = Ifges(4,4) * t737 + Ifges(4,2) * t736 - Ifges(4,6) * t784;
t589 = mrSges(4,2) * t698 - mrSges(4,3) * t665 + Ifges(4,1) * t724 + Ifges(4,4) * t723 - Ifges(4,5) * t748 - pkin(9) * t612 + t598 * t770 - t606 * t765 + t706 * t736 + t707 * t784;
t726 = Ifges(3,3) * t756 + (Ifges(3,5) * t766 + Ifges(3,6) * t771) * t788;
t727 = Ifges(3,6) * t756 + (Ifges(3,4) * t766 + Ifges(3,2) * t771) * t788;
t586 = mrSges(3,2) * t729 - mrSges(3,3) * t712 + Ifges(3,1) * t747 + Ifges(3,4) * t748 + Ifges(3,5) * t755 - qJ(3) * t605 - t587 * t759 + t589 * t761 + t726 * t784 - t727 * t756;
t728 = Ifges(3,5) * t756 + (Ifges(3,1) * t766 + Ifges(3,4) * t771) * t788;
t588 = -mrSges(5,1) * t647 + mrSges(5,2) * t648 - pkin(2) * t605 + (Ifges(3,2) + Ifges(4,3)) * t748 - t769 * t613 - t764 * t614 + Ifges(3,6) * t755 + t756 * t728 + Ifges(3,4) * t747 + t736 * t708 - t737 * t707 - Ifges(5,3) * t740 - Ifges(4,5) * t724 - mrSges(3,1) * t729 + t716 * t691 - t717 * t690 - Ifges(4,6) * t723 + mrSges(3,3) * t713 - Ifges(5,6) * t684 - Ifges(5,5) * t685 - mrSges(4,1) * t665 + mrSges(4,2) * t666 - pkin(4) * t775 - pkin(3) * t612 - pkin(10) * t780 - t726 * t785;
t778 = pkin(8) * t597 + t586 * t766 + t588 * t771;
t585 = Ifges(3,5) * t747 + Ifges(3,6) * t748 + Ifges(3,3) * t755 + mrSges(3,1) * t712 - mrSges(3,2) * t713 + t759 * t589 + t761 * t587 + pkin(2) * t774 + qJ(3) * t782 + (t727 * t766 - t728 * t771) * t788;
t584 = -mrSges(2,2) * g(3) - mrSges(2,3) * t751 + Ifges(2,5) * qJDD(1) - t773 * Ifges(2,6) + t771 * t586 - t766 * t588 + (-t592 * t760 - t593 * t762) * pkin(8);
t583 = mrSges(2,1) * g(3) + mrSges(2,3) * t752 + t773 * Ifges(2,5) + Ifges(2,6) * qJDD(1) - pkin(1) * t592 - t760 * t585 + t762 * t778;
t1 = [-m(1) * g(1) + t783; -m(1) * g(2) + t790; (-m(1) - m(2)) * g(3) + t592; -mrSges(1,2) * g(3) + mrSges(1,3) * g(2) - pkin(7) * t790 - t767 * t583 + t772 * t584; mrSges(1,1) * g(3) - mrSges(1,3) * g(1) + pkin(7) * t783 + t772 * t583 + t767 * t584; -mrSges(1,1) * g(2) + mrSges(2,1) * t751 + mrSges(1,2) * g(1) - mrSges(2,2) * t752 + Ifges(2,3) * qJDD(1) + pkin(1) * t593 + t762 * t585 + t760 * t778;];
tauB  = t1;
