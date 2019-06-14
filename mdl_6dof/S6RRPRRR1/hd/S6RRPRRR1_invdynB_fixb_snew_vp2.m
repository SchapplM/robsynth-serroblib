% Calculate vector of inverse dynamics base forces with Newton-Euler for
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
% tauB [6x1]
%   base forces of inverse dynamics (contains inertial, gravitational coriolis and centrifugal forces)

% Quelle: HybrDyn-Toolbox
% Datum: 2019-05-06 19:35
% Revision: 89c353f7eff3bd693eda4e29f35b2761dbc3ada0 (2019-05-03)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauB = S6RRPRRR1_invdynB_fixb_snew_vp2(qJ, qJD, qJDD, g, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(6,1),zeros(3,1),zeros(11,1),zeros(7,1),zeros(7,3),zeros(7,6)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRPRRR1_invdynB_fixb_snew_vp2: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RRPRRR1_invdynB_fixb_snew_vp2: qJD has to be [6x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [6 1]), ...
  'S6RRPRRR1_invdynB_fixb_snew_vp2: qJDD has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6RRPRRR1_invdynB_fixb_snew_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [11 1]), ...
  'S6RRPRRR1_invdynB_fixb_snew_vp2: pkin has to be [11x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6RRPRRR1_invdynB_fixb_snew_vp2: m has to be [7x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [7,3]), ...
  'S6RRPRRR1_invdynB_fixb_snew_vp2: mrSges has to be [7x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [7 6]), ...
  'S6RRPRRR1_invdynB_fixb_snew_vp2: Ifges has to be [7x6] (double)'); 

%% Symbolic Calculation
% From invdyn_fixb_NewtonEuler_linkframe_tauB_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-05-06 19:30:46
% EndTime: 2019-05-06 19:31:19
% DurationCPUTime: 31.62s
% Computational Cost: add. (471365->387), mult. (1135998->490), div. (0->0), fcn. (870270->12), ass. (0->149)
t734 = qJD(1) ^ 2;
t752 = pkin(2) * t734;
t728 = sin(qJ(1));
t733 = cos(qJ(1));
t714 = -g(1) * t733 - g(2) * t728;
t703 = -pkin(1) * t734 + qJDD(1) * pkin(7) + t714;
t727 = sin(qJ(2));
t751 = t727 * t703;
t732 = cos(qJ(2));
t747 = qJD(1) * qJD(2);
t708 = qJDD(1) * t727 + t732 * t747;
t670 = qJDD(2) * pkin(2) - t708 * qJ(3) - t751 + (qJ(3) * t747 + t727 * t752 - g(3)) * t732;
t689 = -g(3) * t727 + t732 * t703;
t709 = qJDD(1) * t732 - t727 * t747;
t749 = qJD(1) * t727;
t710 = qJD(2) * pkin(2) - qJ(3) * t749;
t721 = t732 ^ 2;
t671 = qJ(3) * t709 - qJD(2) * t710 - t721 * t752 + t689;
t722 = sin(pkin(11));
t723 = cos(pkin(11));
t698 = (t722 * t732 + t723 * t727) * qJD(1);
t647 = -0.2e1 * qJD(3) * t698 + t723 * t670 - t722 * t671;
t687 = t708 * t723 + t709 * t722;
t697 = (-t722 * t727 + t723 * t732) * qJD(1);
t632 = (qJD(2) * t697 - t687) * pkin(8) + (t697 * t698 + qJDD(2)) * pkin(3) + t647;
t648 = 0.2e1 * qJD(3) * t697 + t722 * t670 + t723 * t671;
t686 = -t708 * t722 + t709 * t723;
t692 = qJD(2) * pkin(3) - pkin(8) * t698;
t696 = t697 ^ 2;
t634 = -pkin(3) * t696 + pkin(8) * t686 - qJD(2) * t692 + t648;
t726 = sin(qJ(4));
t731 = cos(qJ(4));
t614 = t731 * t632 - t726 * t634;
t680 = t697 * t731 - t698 * t726;
t652 = qJD(4) * t680 + t686 * t726 + t687 * t731;
t681 = t697 * t726 + t698 * t731;
t719 = qJDD(2) + qJDD(4);
t720 = qJD(2) + qJD(4);
t611 = (t680 * t720 - t652) * pkin(9) + (t680 * t681 + t719) * pkin(4) + t614;
t615 = t726 * t632 + t731 * t634;
t651 = -qJD(4) * t681 + t686 * t731 - t687 * t726;
t675 = pkin(4) * t720 - pkin(9) * t681;
t676 = t680 ^ 2;
t613 = -pkin(4) * t676 + pkin(9) * t651 - t675 * t720 + t615;
t725 = sin(qJ(5));
t730 = cos(qJ(5));
t608 = t725 * t611 + t730 * t613;
t665 = t680 * t725 + t681 * t730;
t624 = -qJD(5) * t665 + t651 * t730 - t652 * t725;
t664 = t680 * t730 - t681 * t725;
t643 = -mrSges(6,1) * t664 + mrSges(6,2) * t665;
t717 = qJD(5) + t720;
t656 = mrSges(6,1) * t717 - mrSges(6,3) * t665;
t716 = qJDD(5) + t719;
t644 = -pkin(5) * t664 - pkin(10) * t665;
t715 = t717 ^ 2;
t606 = -pkin(5) * t715 + pkin(10) * t716 + t644 * t664 + t608;
t713 = t728 * g(1) - t733 * g(2);
t740 = -qJDD(1) * pkin(1) - t713;
t672 = -t709 * pkin(2) + qJDD(3) + t710 * t749 + (-qJ(3) * t721 - pkin(7)) * t734 + t740;
t646 = -t686 * pkin(3) - t696 * pkin(8) + t698 * t692 + t672;
t620 = -t651 * pkin(4) - t676 * pkin(9) + t681 * t675 + t646;
t625 = qJD(5) * t664 + t651 * t725 + t652 * t730;
t609 = t620 + (-t664 * t717 - t625) * pkin(10) + (t665 * t717 - t624) * pkin(5);
t724 = sin(qJ(6));
t729 = cos(qJ(6));
t603 = -t606 * t724 + t609 * t729;
t653 = -t665 * t724 + t717 * t729;
t618 = qJD(6) * t653 + t625 * t729 + t716 * t724;
t623 = qJDD(6) - t624;
t654 = t665 * t729 + t717 * t724;
t635 = -mrSges(7,1) * t653 + mrSges(7,2) * t654;
t660 = qJD(6) - t664;
t636 = -mrSges(7,2) * t660 + mrSges(7,3) * t653;
t601 = m(7) * t603 + mrSges(7,1) * t623 - mrSges(7,3) * t618 - t635 * t654 + t636 * t660;
t604 = t606 * t729 + t609 * t724;
t617 = -qJD(6) * t654 - t625 * t724 + t716 * t729;
t637 = mrSges(7,1) * t660 - mrSges(7,3) * t654;
t602 = m(7) * t604 - mrSges(7,2) * t623 + mrSges(7,3) * t617 + t635 * t653 - t637 * t660;
t741 = -t601 * t724 + t729 * t602;
t592 = m(6) * t608 - mrSges(6,2) * t716 + mrSges(6,3) * t624 + t643 * t664 - t656 * t717 + t741;
t607 = t611 * t730 - t613 * t725;
t655 = -mrSges(6,2) * t717 + mrSges(6,3) * t664;
t605 = -pkin(5) * t716 - pkin(10) * t715 + t644 * t665 - t607;
t738 = -m(7) * t605 + t617 * mrSges(7,1) - mrSges(7,2) * t618 + t653 * t636 - t637 * t654;
t597 = m(6) * t607 + mrSges(6,1) * t716 - mrSges(6,3) * t625 - t643 * t665 + t655 * t717 + t738;
t587 = t725 * t592 + t730 * t597;
t666 = -mrSges(5,1) * t680 + mrSges(5,2) * t681;
t673 = -mrSges(5,2) * t720 + mrSges(5,3) * t680;
t585 = m(5) * t614 + mrSges(5,1) * t719 - mrSges(5,3) * t652 - t666 * t681 + t673 * t720 + t587;
t674 = mrSges(5,1) * t720 - mrSges(5,3) * t681;
t742 = t730 * t592 - t597 * t725;
t586 = m(5) * t615 - mrSges(5,2) * t719 + mrSges(5,3) * t651 + t666 * t680 - t674 * t720 + t742;
t579 = t731 * t585 + t726 * t586;
t684 = -mrSges(4,1) * t697 + mrSges(4,2) * t698;
t690 = -qJD(2) * mrSges(4,2) + mrSges(4,3) * t697;
t577 = m(4) * t647 + qJDD(2) * mrSges(4,1) - mrSges(4,3) * t687 + qJD(2) * t690 - t684 * t698 + t579;
t691 = qJD(2) * mrSges(4,1) - mrSges(4,3) * t698;
t743 = -t585 * t726 + t731 * t586;
t578 = m(4) * t648 - qJDD(2) * mrSges(4,2) + mrSges(4,3) * t686 - qJD(2) * t691 + t684 * t697 + t743;
t572 = t723 * t577 + t722 * t578;
t688 = -t732 * g(3) - t751;
t707 = (-mrSges(3,1) * t732 + mrSges(3,2) * t727) * qJD(1);
t748 = qJD(1) * t732;
t712 = -qJD(2) * mrSges(3,2) + mrSges(3,3) * t748;
t570 = m(3) * t688 + qJDD(2) * mrSges(3,1) - mrSges(3,3) * t708 + qJD(2) * t712 - t707 * t749 + t572;
t711 = qJD(2) * mrSges(3,1) - mrSges(3,3) * t749;
t744 = -t577 * t722 + t723 * t578;
t571 = m(3) * t689 - qJDD(2) * mrSges(3,2) + mrSges(3,3) * t709 - qJD(2) * t711 + t707 * t748 + t744;
t745 = -t570 * t727 + t732 * t571;
t564 = m(2) * t714 - mrSges(2,1) * t734 - qJDD(1) * mrSges(2,2) + t745;
t702 = -t734 * pkin(7) + t740;
t593 = t729 * t601 + t724 * t602;
t739 = m(6) * t620 - t624 * mrSges(6,1) + t625 * mrSges(6,2) - t664 * t655 + t665 * t656 + t593;
t737 = m(5) * t646 - t651 * mrSges(5,1) + t652 * mrSges(5,2) - t680 * t673 + t681 * t674 + t739;
t736 = m(4) * t672 - t686 * mrSges(4,1) + t687 * mrSges(4,2) - t697 * t690 + t698 * t691 + t737;
t735 = -m(3) * t702 + t709 * mrSges(3,1) - t708 * mrSges(3,2) - t711 * t749 + t712 * t748 - t736;
t589 = m(2) * t713 + qJDD(1) * mrSges(2,1) - t734 * mrSges(2,2) + t735;
t750 = t728 * t564 + t733 * t589;
t565 = t732 * t570 + t727 * t571;
t746 = t733 * t564 - t589 * t728;
t701 = Ifges(3,5) * qJD(2) + (Ifges(3,1) * t727 + Ifges(3,4) * t732) * qJD(1);
t700 = Ifges(3,6) * qJD(2) + (Ifges(3,4) * t727 + Ifges(3,2) * t732) * qJD(1);
t699 = Ifges(3,3) * qJD(2) + (Ifges(3,5) * t727 + Ifges(3,6) * t732) * qJD(1);
t679 = Ifges(4,1) * t698 + Ifges(4,4) * t697 + Ifges(4,5) * qJD(2);
t678 = Ifges(4,4) * t698 + Ifges(4,2) * t697 + Ifges(4,6) * qJD(2);
t677 = Ifges(4,5) * t698 + Ifges(4,6) * t697 + Ifges(4,3) * qJD(2);
t659 = Ifges(5,1) * t681 + Ifges(5,4) * t680 + Ifges(5,5) * t720;
t658 = Ifges(5,4) * t681 + Ifges(5,2) * t680 + Ifges(5,6) * t720;
t657 = Ifges(5,5) * t681 + Ifges(5,6) * t680 + Ifges(5,3) * t720;
t640 = Ifges(6,1) * t665 + Ifges(6,4) * t664 + Ifges(6,5) * t717;
t639 = Ifges(6,4) * t665 + Ifges(6,2) * t664 + Ifges(6,6) * t717;
t638 = Ifges(6,5) * t665 + Ifges(6,6) * t664 + Ifges(6,3) * t717;
t628 = Ifges(7,1) * t654 + Ifges(7,4) * t653 + Ifges(7,5) * t660;
t627 = Ifges(7,4) * t654 + Ifges(7,2) * t653 + Ifges(7,6) * t660;
t626 = Ifges(7,5) * t654 + Ifges(7,6) * t653 + Ifges(7,3) * t660;
t595 = mrSges(7,2) * t605 - mrSges(7,3) * t603 + Ifges(7,1) * t618 + Ifges(7,4) * t617 + Ifges(7,5) * t623 + t626 * t653 - t627 * t660;
t594 = -mrSges(7,1) * t605 + mrSges(7,3) * t604 + Ifges(7,4) * t618 + Ifges(7,2) * t617 + Ifges(7,6) * t623 - t626 * t654 + t628 * t660;
t581 = -mrSges(6,1) * t620 - mrSges(7,1) * t603 + mrSges(7,2) * t604 + mrSges(6,3) * t608 + Ifges(6,4) * t625 - Ifges(7,5) * t618 + Ifges(6,2) * t624 + Ifges(6,6) * t716 - Ifges(7,6) * t617 - Ifges(7,3) * t623 - pkin(5) * t593 - t627 * t654 + t628 * t653 - t638 * t665 + t640 * t717;
t580 = mrSges(6,2) * t620 - mrSges(6,3) * t607 + Ifges(6,1) * t625 + Ifges(6,4) * t624 + Ifges(6,5) * t716 - pkin(10) * t593 - t594 * t724 + t595 * t729 + t638 * t664 - t639 * t717;
t573 = mrSges(5,2) * t646 - mrSges(5,3) * t614 + Ifges(5,1) * t652 + Ifges(5,4) * t651 + Ifges(5,5) * t719 - pkin(9) * t587 + t580 * t730 - t581 * t725 + t657 * t680 - t658 * t720;
t566 = -mrSges(5,1) * t646 + mrSges(5,3) * t615 + Ifges(5,4) * t652 + Ifges(5,2) * t651 + Ifges(5,6) * t719 - pkin(4) * t739 + pkin(9) * t742 + t725 * t580 + t730 * t581 - t681 * t657 + t720 * t659;
t561 = mrSges(4,2) * t672 - mrSges(4,3) * t647 + Ifges(4,1) * t687 + Ifges(4,4) * t686 + Ifges(4,5) * qJDD(2) - pkin(8) * t579 - qJD(2) * t678 - t566 * t726 + t573 * t731 + t677 * t697;
t560 = -mrSges(4,1) * t672 + mrSges(4,3) * t648 + Ifges(4,4) * t687 + Ifges(4,2) * t686 + Ifges(4,6) * qJDD(2) - pkin(3) * t737 + pkin(8) * t743 + qJD(2) * t679 + t731 * t566 + t726 * t573 - t698 * t677;
t559 = -pkin(10) * t741 + t734 * Ifges(2,5) - t729 * t594 - t724 * t595 - Ifges(6,3) * t716 - Ifges(5,3) * t719 - Ifges(3,5) * t708 - Ifges(3,6) * t709 + mrSges(2,3) * t714 + t697 * t679 - t698 * t678 - t681 * t658 - Ifges(4,6) * t686 - Ifges(4,5) * t687 - mrSges(3,1) * t688 + mrSges(3,2) * t689 + t680 * t659 + t664 * t640 - t665 * t639 - Ifges(5,6) * t651 - Ifges(5,5) * t652 - mrSges(4,1) * t647 + mrSges(4,2) * t648 - Ifges(6,6) * t624 - Ifges(6,5) * t625 - mrSges(5,1) * t614 + mrSges(5,2) * t615 + mrSges(6,2) * t608 - mrSges(6,1) * t607 + mrSges(2,1) * g(3) - pkin(5) * t738 + (-t700 * t727 + t701 * t732) * qJD(1) - pkin(4) * t587 + Ifges(2,6) * qJDD(1) + (-Ifges(3,3) - Ifges(4,3)) * qJDD(2) - pkin(1) * t565 - pkin(3) * t579 - pkin(2) * t572;
t558 = mrSges(3,2) * t702 - mrSges(3,3) * t688 + Ifges(3,1) * t708 + Ifges(3,4) * t709 + Ifges(3,5) * qJDD(2) - qJ(3) * t572 - qJD(2) * t700 - t560 * t722 + t561 * t723 + t699 * t748;
t557 = -mrSges(3,1) * t702 + mrSges(3,3) * t689 + Ifges(3,4) * t708 + Ifges(3,2) * t709 + Ifges(3,6) * qJDD(2) - pkin(2) * t736 + qJ(3) * t744 + qJD(2) * t701 + t723 * t560 + t722 * t561 - t699 * t749;
t556 = -mrSges(2,2) * g(3) - mrSges(2,3) * t713 + Ifges(2,5) * qJDD(1) - Ifges(2,6) * t734 - pkin(7) * t565 - t557 * t727 + t558 * t732;
t1 = [-m(1) * g(1) + t746; -m(1) * g(2) + t750; (-m(1) - m(2)) * g(3) + t565; -mrSges(1,2) * g(3) + mrSges(1,3) * g(2) - pkin(6) * t750 + t733 * t556 - t728 * t559; mrSges(1,1) * g(3) - mrSges(1,3) * g(1) + pkin(6) * t746 + t728 * t556 + t733 * t559; -mrSges(1,1) * g(2) + mrSges(2,1) * t713 + mrSges(1,2) * g(1) - mrSges(2,2) * t714 + Ifges(2,3) * qJDD(1) + pkin(1) * t735 + pkin(7) * t745 + t732 * t557 + t727 * t558;];
tauB  = t1;
