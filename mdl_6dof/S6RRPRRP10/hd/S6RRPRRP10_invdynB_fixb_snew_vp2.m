% Calculate vector of inverse dynamics base forces with Newton-Euler for
% S6RRPRRP10
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
%   pkin=[a2,a3,a4,a5,a6,alpha2,d1,d2,d4,d5,theta3]';
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
% Datum: 2019-05-06 18:51
% Revision: 89c353f7eff3bd693eda4e29f35b2761dbc3ada0 (2019-05-03)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauB = S6RRPRRP10_invdynB_fixb_snew_vp2(qJ, qJD, qJDD, g, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(6,1),zeros(3,1),zeros(11,1),zeros(7,1),zeros(7,3),zeros(7,6)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRPRRP10_invdynB_fixb_snew_vp2: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RRPRRP10_invdynB_fixb_snew_vp2: qJD has to be [6x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [6 1]), ...
  'S6RRPRRP10_invdynB_fixb_snew_vp2: qJDD has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6RRPRRP10_invdynB_fixb_snew_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [11 1]), ...
  'S6RRPRRP10_invdynB_fixb_snew_vp2: pkin has to be [11x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6RRPRRP10_invdynB_fixb_snew_vp2: m has to be [7x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [7,3]), ...
  'S6RRPRRP10_invdynB_fixb_snew_vp2: mrSges has to be [7x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [7 6]), ...
  'S6RRPRRP10_invdynB_fixb_snew_vp2: Ifges has to be [7x6] (double)'); 

%% Symbolic Calculation
% From invdyn_fixb_NewtonEuler_linkframe_tauB_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-05-06 18:40:47
% EndTime: 2019-05-06 18:41:14
% DurationCPUTime: 24.34s
% Computational Cost: add. (384704->374), mult. (871766->474), div. (0->0), fcn. (698135->12), ass. (0->152)
t785 = Ifges(6,1) + Ifges(7,1);
t778 = Ifges(6,4) - Ifges(7,5);
t784 = -Ifges(6,5) - Ifges(7,4);
t783 = Ifges(6,2) + Ifges(7,3);
t776 = Ifges(6,6) - Ifges(7,6);
t782 = -Ifges(6,3) - Ifges(7,2);
t781 = cos(qJ(5));
t742 = cos(pkin(6));
t780 = t742 * g(3);
t779 = -mrSges(6,3) - mrSges(7,2);
t740 = sin(pkin(6));
t745 = sin(qJ(2));
t775 = t740 * t745;
t748 = cos(qJ(2));
t774 = t740 * t748;
t773 = t742 * t745;
t772 = t742 * t748;
t746 = sin(qJ(1));
t749 = cos(qJ(1));
t731 = t746 * g(1) - t749 * g(2);
t750 = qJD(1) ^ 2;
t722 = t750 * t740 * pkin(8) + qJDD(1) * pkin(1) + t731;
t732 = -t749 * g(1) - t746 * g(2);
t763 = qJDD(1) * t740;
t723 = -t750 * pkin(1) + pkin(8) * t763 + t732;
t766 = t722 * t773 + t748 * t723;
t693 = -g(3) * t775 + t766;
t736 = t742 * qJD(1) + qJD(2);
t765 = qJD(1) * t740;
t761 = t745 * t765;
t720 = t736 * mrSges(3,1) - mrSges(3,3) * t761;
t725 = (-mrSges(3,1) * t748 + mrSges(3,2) * t745) * t765;
t727 = -qJD(2) * t761 + t748 * t763;
t735 = t742 * qJDD(1) + qJDD(2);
t724 = (-pkin(2) * t748 - qJ(3) * t745) * t765;
t734 = t736 ^ 2;
t764 = qJD(1) * t748;
t681 = -t734 * pkin(2) + t735 * qJ(3) + (-g(3) * t745 + t724 * t764) * t740 + t766;
t726 = (qJD(2) * t764 + qJDD(1) * t745) * t740;
t682 = -t727 * pkin(2) - t780 - t726 * qJ(3) + (-t722 + (pkin(2) * t745 - qJ(3) * t748) * t736 * qJD(1)) * t740;
t739 = sin(pkin(11));
t741 = cos(pkin(11));
t716 = t739 * t736 + t741 * t761;
t640 = -0.2e1 * qJD(3) * t716 - t739 * t681 + t741 * t682;
t703 = t741 * t726 + t739 * t735;
t715 = t741 * t736 - t739 * t761;
t760 = t740 * t764;
t636 = (-t715 * t760 - t703) * pkin(9) + (t715 * t716 - t727) * pkin(3) + t640;
t641 = 0.2e1 * qJD(3) * t715 + t741 * t681 + t739 * t682;
t702 = -t739 * t726 + t741 * t735;
t704 = -pkin(3) * t760 - t716 * pkin(9);
t713 = t715 ^ 2;
t639 = -t713 * pkin(3) + t702 * pkin(9) + t704 * t760 + t641;
t744 = sin(qJ(4));
t747 = cos(qJ(4));
t632 = t744 * t636 + t747 * t639;
t697 = t744 * t715 + t747 * t716;
t665 = -t697 * qJD(4) + t747 * t702 - t744 * t703;
t696 = t747 * t715 - t744 * t716;
t675 = -t696 * mrSges(5,1) + t697 * mrSges(5,2);
t730 = qJD(4) - t760;
t686 = t730 * mrSges(5,1) - t697 * mrSges(5,3);
t719 = qJDD(4) - t727;
t676 = -t696 * pkin(4) - t697 * pkin(10);
t729 = t730 ^ 2;
t630 = -t729 * pkin(4) + t719 * pkin(10) + t696 * t676 + t632;
t692 = -g(3) * t774 + t722 * t772 - t745 * t723;
t680 = -t735 * pkin(2) - t734 * qJ(3) + t724 * t761 + qJDD(3) - t692;
t646 = -t702 * pkin(3) - t713 * pkin(9) + t716 * t704 + t680;
t666 = t696 * qJD(4) + t744 * t702 + t747 * t703;
t634 = (-t696 * t730 - t666) * pkin(10) + (t697 * t730 - t665) * pkin(4) + t646;
t743 = sin(qJ(5));
t627 = t781 * t630 + t743 * t634;
t684 = t781 * t697 + t743 * t730;
t644 = t684 * qJD(5) + t743 * t666 - t781 * t719;
t664 = qJDD(5) - t665;
t695 = qJD(5) - t696;
t669 = t695 * mrSges(6,1) - t684 * mrSges(6,3);
t683 = t743 * t697 - t781 * t730;
t658 = t683 * pkin(5) - t684 * qJ(6);
t694 = t695 ^ 2;
t623 = -t694 * pkin(5) + t664 * qJ(6) + 0.2e1 * qJD(6) * t695 - t683 * t658 + t627;
t670 = -t695 * mrSges(7,1) + t684 * mrSges(7,2);
t762 = m(7) * t623 + t664 * mrSges(7,3) + t695 * t670;
t659 = t683 * mrSges(7,1) - t684 * mrSges(7,3);
t767 = -t683 * mrSges(6,1) - t684 * mrSges(6,2) - t659;
t618 = m(6) * t627 - t664 * mrSges(6,2) + t779 * t644 - t695 * t669 + t767 * t683 + t762;
t626 = -t743 * t630 + t781 * t634;
t645 = -t683 * qJD(5) + t781 * t666 + t743 * t719;
t668 = -t695 * mrSges(6,2) - t683 * mrSges(6,3);
t624 = -t664 * pkin(5) - t694 * qJ(6) + t684 * t658 + qJDD(6) - t626;
t667 = -t683 * mrSges(7,2) + t695 * mrSges(7,3);
t755 = -m(7) * t624 + t664 * mrSges(7,1) + t695 * t667;
t620 = m(6) * t626 + t664 * mrSges(6,1) + t779 * t645 + t695 * t668 + t767 * t684 + t755;
t756 = t781 * t618 - t743 * t620;
t610 = m(5) * t632 - t719 * mrSges(5,2) + t665 * mrSges(5,3) + t696 * t675 - t730 * t686 + t756;
t631 = t747 * t636 - t744 * t639;
t685 = -t730 * mrSges(5,2) + t696 * mrSges(5,3);
t629 = -t719 * pkin(4) - t729 * pkin(10) + t697 * t676 - t631;
t625 = -0.2e1 * qJD(6) * t684 + (t683 * t695 - t645) * qJ(6) + (t684 * t695 + t644) * pkin(5) + t629;
t621 = m(7) * t625 + t644 * mrSges(7,1) - t645 * mrSges(7,3) + t683 * t667 - t684 * t670;
t752 = -m(6) * t629 - t644 * mrSges(6,1) - t645 * mrSges(6,2) - t683 * t668 - t684 * t669 - t621;
t615 = m(5) * t631 + t719 * mrSges(5,1) - t666 * mrSges(5,3) - t697 * t675 + t730 * t685 + t752;
t604 = t744 * t610 + t747 * t615;
t698 = -t715 * mrSges(4,1) + t716 * mrSges(4,2);
t700 = mrSges(4,2) * t760 + t715 * mrSges(4,3);
t602 = m(4) * t640 - t727 * mrSges(4,1) - t703 * mrSges(4,3) - t716 * t698 - t700 * t760 + t604;
t701 = -mrSges(4,1) * t760 - t716 * mrSges(4,3);
t757 = t747 * t610 - t744 * t615;
t603 = m(4) * t641 + t727 * mrSges(4,2) + t702 * mrSges(4,3) + t715 * t698 + t701 * t760 + t757;
t758 = -t739 * t602 + t741 * t603;
t593 = m(3) * t693 - t735 * mrSges(3,2) + t727 * mrSges(3,3) - t736 * t720 + t725 * t760 + t758;
t596 = t741 * t602 + t739 * t603;
t708 = -t740 * t722 - t780;
t721 = -t736 * mrSges(3,2) + mrSges(3,3) * t760;
t595 = m(3) * t708 - t727 * mrSges(3,1) + t726 * mrSges(3,2) + (t720 * t745 - t721 * t748) * t765 + t596;
t613 = t743 * t618 + t781 * t620;
t753 = m(5) * t646 - t665 * mrSges(5,1) + t666 * mrSges(5,2) - t696 * t685 + t697 * t686 + t613;
t751 = -m(4) * t680 + t702 * mrSges(4,1) - t703 * mrSges(4,2) + t715 * t700 - t716 * t701 - t753;
t607 = m(3) * t692 + t735 * mrSges(3,1) - t726 * mrSges(3,3) + t736 * t721 - t725 * t761 + t751;
t584 = t593 * t773 - t740 * t595 + t607 * t772;
t582 = m(2) * t731 + qJDD(1) * mrSges(2,1) - t750 * mrSges(2,2) + t584;
t589 = t748 * t593 - t745 * t607;
t588 = m(2) * t732 - t750 * mrSges(2,1) - qJDD(1) * mrSges(2,2) + t589;
t771 = t749 * t582 + t746 * t588;
t770 = t683 * t783 - t684 * t778 - t695 * t776;
t769 = t683 * t776 + t684 * t784 + t695 * t782;
t768 = -t778 * t683 + t684 * t785 - t784 * t695;
t583 = t593 * t775 + t742 * t595 + t607 * t774;
t759 = -t746 * t582 + t749 * t588;
t611 = -mrSges(6,1) * t629 - mrSges(7,1) * t625 + mrSges(7,2) * t623 + mrSges(6,3) * t627 - pkin(5) * t621 - t644 * t783 + t778 * t645 + t776 * t664 + t769 * t684 + t768 * t695;
t612 = mrSges(6,2) * t629 + mrSges(7,2) * t624 - mrSges(6,3) * t626 - mrSges(7,3) * t625 - qJ(6) * t621 - t778 * t644 + t645 * t785 - t664 * t784 + t769 * t683 + t770 * t695;
t671 = Ifges(5,5) * t697 + Ifges(5,6) * t696 + Ifges(5,3) * t730;
t672 = Ifges(5,4) * t697 + Ifges(5,2) * t696 + Ifges(5,6) * t730;
t597 = mrSges(5,2) * t646 - mrSges(5,3) * t631 + Ifges(5,1) * t666 + Ifges(5,4) * t665 + Ifges(5,5) * t719 - pkin(10) * t613 - t743 * t611 + t781 * t612 + t696 * t671 - t730 * t672;
t673 = Ifges(5,1) * t697 + Ifges(5,4) * t696 + Ifges(5,5) * t730;
t598 = Ifges(5,4) * t666 + Ifges(5,2) * t665 + Ifges(5,6) * t719 - t697 * t671 + t730 * t673 - mrSges(5,1) * t646 + mrSges(5,3) * t632 - mrSges(6,1) * t626 + mrSges(6,2) * t627 + mrSges(7,1) * t624 - mrSges(7,3) * t623 - pkin(5) * t755 - qJ(6) * t762 - pkin(4) * t613 + (pkin(5) * t659 + t770) * t684 + (qJ(6) * t659 - t768) * t683 + t782 * t664 + (pkin(5) * mrSges(7,2) + t784) * t645 + (qJ(6) * mrSges(7,2) + t776) * t644;
t687 = Ifges(4,5) * t716 + Ifges(4,6) * t715 - Ifges(4,3) * t760;
t689 = Ifges(4,1) * t716 + Ifges(4,4) * t715 - Ifges(4,5) * t760;
t580 = -mrSges(4,1) * t680 + mrSges(4,3) * t641 + Ifges(4,4) * t703 + Ifges(4,2) * t702 - Ifges(4,6) * t727 - pkin(3) * t753 + pkin(9) * t757 + t744 * t597 + t747 * t598 - t716 * t687 - t689 * t760;
t688 = Ifges(4,4) * t716 + Ifges(4,2) * t715 - Ifges(4,6) * t760;
t585 = mrSges(4,2) * t680 - mrSges(4,3) * t640 + Ifges(4,1) * t703 + Ifges(4,4) * t702 - Ifges(4,5) * t727 - pkin(9) * t604 + t747 * t597 - t744 * t598 + t715 * t687 + t688 * t760;
t705 = Ifges(3,3) * t736 + (Ifges(3,5) * t745 + Ifges(3,6) * t748) * t765;
t706 = Ifges(3,6) * t736 + (Ifges(3,4) * t745 + Ifges(3,2) * t748) * t765;
t578 = mrSges(3,2) * t708 - mrSges(3,3) * t692 + Ifges(3,1) * t726 + Ifges(3,4) * t727 + Ifges(3,5) * t735 - qJ(3) * t596 - t739 * t580 + t741 * t585 + t705 * t760 - t736 * t706;
t707 = Ifges(3,5) * t736 + (Ifges(3,1) * t745 + Ifges(3,4) * t748) * t765;
t579 = -pkin(2) * t596 - mrSges(4,1) * t640 + mrSges(4,2) * t641 - pkin(3) * t604 - t781 * t611 + (Ifges(3,2) + Ifges(4,3)) * t727 - t705 * t761 - pkin(10) * t756 - pkin(4) * t752 - mrSges(5,1) * t631 + mrSges(5,2) * t632 - Ifges(5,6) * t665 - Ifges(5,5) * t666 + mrSges(3,3) * t693 + t696 * t673 - t697 * t672 - Ifges(4,6) * t702 - Ifges(4,5) * t703 - mrSges(3,1) * t708 + t715 * t689 - t716 * t688 - Ifges(5,3) * t719 + Ifges(3,4) * t726 + Ifges(3,6) * t735 + t736 * t707 - t743 * t612;
t754 = pkin(8) * t589 + t578 * t745 + t579 * t748;
t577 = Ifges(3,5) * t726 + Ifges(3,6) * t727 + Ifges(3,3) * t735 + mrSges(3,1) * t692 - mrSges(3,2) * t693 + t739 * t585 + t741 * t580 + pkin(2) * t751 + qJ(3) * t758 + (t706 * t745 - t707 * t748) * t765;
t576 = -mrSges(2,2) * g(3) - mrSges(2,3) * t731 + Ifges(2,5) * qJDD(1) - t750 * Ifges(2,6) + t748 * t578 - t745 * t579 + (-t583 * t740 - t584 * t742) * pkin(8);
t575 = mrSges(2,1) * g(3) + mrSges(2,3) * t732 + t750 * Ifges(2,5) + Ifges(2,6) * qJDD(1) - pkin(1) * t583 - t740 * t577 + t754 * t742;
t1 = [-m(1) * g(1) + t759; -m(1) * g(2) + t771; (-m(1) - m(2)) * g(3) + t583; -mrSges(1,2) * g(3) + mrSges(1,3) * g(2) - pkin(7) * t771 - t746 * t575 + t749 * t576; mrSges(1,1) * g(3) - mrSges(1,3) * g(1) + pkin(7) * t759 + t749 * t575 + t746 * t576; -mrSges(1,1) * g(2) + mrSges(2,1) * t731 + mrSges(1,2) * g(1) - mrSges(2,2) * t732 + Ifges(2,3) * qJDD(1) + pkin(1) * t584 + t742 * t577 + t754 * t740;];
tauB  = t1;
