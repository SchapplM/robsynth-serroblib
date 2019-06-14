% Calculate vector of inverse dynamics base forces with Newton-Euler for
% S6PRPRRR3
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
%   pkin=[a2,a3,a4,a5,a6,alpha2,d2,d4,d5,d6,theta1,theta3]';
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
% Datum: 2019-05-05 00:44
% Revision: 89c353f7eff3bd693eda4e29f35b2761dbc3ada0 (2019-05-03)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauB = S6PRPRRR3_invdynB_fixb_snew_vp2(qJ, qJD, qJDD, g, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(6,1),zeros(3,1),zeros(12,1),zeros(7,1),zeros(7,3),zeros(7,6)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6PRPRRR3_invdynB_fixb_snew_vp2: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6PRPRRR3_invdynB_fixb_snew_vp2: qJD has to be [6x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [6 1]), ...
  'S6PRPRRR3_invdynB_fixb_snew_vp2: qJDD has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6PRPRRR3_invdynB_fixb_snew_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [12 1]), ...
  'S6PRPRRR3_invdynB_fixb_snew_vp2: pkin has to be [12x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6PRPRRR3_invdynB_fixb_snew_vp2: m has to be [7x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [7,3]), ...
  'S6PRPRRR3_invdynB_fixb_snew_vp2: mrSges has to be [7x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [7 6]), ...
  'S6PRPRRR3_invdynB_fixb_snew_vp2: Ifges has to be [7x6] (double)'); 

%% Symbolic Calculation
% From invdyn_fixb_NewtonEuler_linkframe_tauB_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-05-05 00:40:46
% EndTime: 2019-05-05 00:41:07
% DurationCPUTime: 20.45s
% Computational Cost: add. (325817->320), mult. (735743->410), div. (0->0), fcn. (577982->14), ass. (0->143)
t686 = qJD(2) ^ 2;
t673 = sin(pkin(11));
t676 = cos(pkin(11));
t658 = g(1) * t673 - g(2) * t676;
t659 = -g(1) * t676 - g(2) * t673;
t671 = -g(3) + qJDD(1);
t674 = sin(pkin(6));
t677 = cos(pkin(6));
t681 = sin(qJ(2));
t685 = cos(qJ(2));
t636 = -t681 * t659 + (t658 * t677 + t671 * t674) * t685;
t675 = cos(pkin(12));
t716 = pkin(3) * t675;
t672 = sin(pkin(12));
t715 = mrSges(4,2) * t672;
t689 = qJDD(3) - t636;
t626 = -qJDD(2) * pkin(2) - t686 * qJ(3) + t689;
t668 = t672 ^ 2;
t669 = t675 ^ 2;
t622 = (-pkin(2) - t716) * qJDD(2) + (-qJ(3) + (-t668 - t669) * pkin(8)) * t686 + t689;
t680 = sin(qJ(4));
t684 = cos(qJ(4));
t695 = t672 * t684 + t675 * t680;
t652 = t695 * qJD(2);
t694 = -t672 * t680 + t675 * t684;
t642 = -t652 * qJD(4) + qJDD(2) * t694;
t651 = t694 * qJD(2);
t707 = t651 * qJD(4);
t643 = qJDD(2) * t695 + t707;
t646 = -qJD(4) * mrSges(5,2) + mrSges(5,3) * t651;
t647 = qJD(4) * mrSges(5,1) - mrSges(5,3) * t652;
t711 = t677 * t681;
t712 = t674 * t681;
t637 = t658 * t711 + t685 * t659 + t671 * t712;
t629 = -pkin(2) * t686 + qJDD(2) * qJ(3) + t637;
t649 = -t658 * t674 + t671 * t677;
t706 = qJD(2) * qJD(3);
t709 = t675 * t649 - 0.2e1 * t672 * t706;
t609 = (-pkin(8) * qJDD(2) + t686 * t716 - t629) * t672 + t709;
t621 = t672 * t649 + (t629 + 0.2e1 * t706) * t675;
t705 = qJDD(2) * t675;
t713 = t669 * t686;
t612 = -pkin(3) * t713 + pkin(8) * t705 + t621;
t589 = t684 * t609 - t680 * t612;
t586 = (-t643 + t707) * pkin(9) + (t651 * t652 + qJDD(4)) * pkin(4) + t589;
t590 = t680 * t609 + t684 * t612;
t648 = qJD(4) * pkin(4) - pkin(9) * t652;
t650 = t651 ^ 2;
t588 = -pkin(4) * t650 + pkin(9) * t642 - qJD(4) * t648 + t590;
t679 = sin(qJ(5));
t683 = cos(qJ(5));
t583 = t679 * t586 + t683 * t588;
t634 = t651 * t683 - t652 * t679;
t635 = t651 * t679 + t652 * t683;
t619 = -pkin(5) * t634 - pkin(10) * t635;
t670 = qJD(4) + qJD(5);
t666 = t670 ^ 2;
t667 = qJDD(4) + qJDD(5);
t581 = -pkin(5) * t666 + pkin(10) * t667 + t619 * t634 + t583;
t595 = -t642 * pkin(4) - t650 * pkin(9) + t652 * t648 + t622;
t603 = -qJD(5) * t635 + t642 * t683 - t643 * t679;
t604 = qJD(5) * t634 + t642 * t679 + t643 * t683;
t584 = (-t634 * t670 - t604) * pkin(10) + (t635 * t670 - t603) * pkin(5) + t595;
t678 = sin(qJ(6));
t682 = cos(qJ(6));
t578 = -t581 * t678 + t584 * t682;
t623 = -t635 * t678 + t670 * t682;
t593 = qJD(6) * t623 + t604 * t682 + t667 * t678;
t602 = qJDD(6) - t603;
t624 = t635 * t682 + t670 * t678;
t605 = -mrSges(7,1) * t623 + mrSges(7,2) * t624;
t630 = qJD(6) - t634;
t610 = -mrSges(7,2) * t630 + mrSges(7,3) * t623;
t576 = m(7) * t578 + mrSges(7,1) * t602 - mrSges(7,3) * t593 - t605 * t624 + t610 * t630;
t579 = t581 * t682 + t584 * t678;
t592 = -qJD(6) * t624 - t604 * t678 + t667 * t682;
t611 = mrSges(7,1) * t630 - mrSges(7,3) * t624;
t577 = m(7) * t579 - mrSges(7,2) * t602 + mrSges(7,3) * t592 + t605 * t623 - t611 * t630;
t568 = t682 * t576 + t678 * t577;
t627 = -mrSges(6,2) * t670 + mrSges(6,3) * t634;
t628 = mrSges(6,1) * t670 - mrSges(6,3) * t635;
t692 = m(6) * t595 - t603 * mrSges(6,1) + t604 * mrSges(6,2) - t634 * t627 + t635 * t628 + t568;
t688 = m(5) * t622 - t642 * mrSges(5,1) + t643 * mrSges(5,2) - t651 * t646 + t652 * t647 + t692;
t687 = -m(4) * t626 + mrSges(4,1) * t705 - t688 + (t668 * t686 + t713) * mrSges(4,3);
t564 = t687 + (mrSges(3,1) - t715) * qJDD(2) - t686 * mrSges(3,2) + m(3) * t636;
t714 = t564 * t685;
t618 = -mrSges(6,1) * t634 + mrSges(6,2) * t635;
t700 = -t576 * t678 + t682 * t577;
t567 = m(6) * t583 - mrSges(6,2) * t667 + mrSges(6,3) * t603 + t618 * t634 - t628 * t670 + t700;
t582 = t586 * t683 - t588 * t679;
t580 = -pkin(5) * t667 - pkin(10) * t666 + t619 * t635 - t582;
t690 = -m(7) * t580 + t592 * mrSges(7,1) - mrSges(7,2) * t593 + t623 * t610 - t611 * t624;
t572 = m(6) * t582 + mrSges(6,1) * t667 - mrSges(6,3) * t604 - t618 * t635 + t627 * t670 + t690;
t561 = t679 * t567 + t683 * t572;
t640 = -mrSges(5,1) * t651 + mrSges(5,2) * t652;
t559 = m(5) * t589 + qJDD(4) * mrSges(5,1) - mrSges(5,3) * t643 + qJD(4) * t646 - t640 * t652 + t561;
t701 = t683 * t567 - t572 * t679;
t560 = m(5) * t590 - qJDD(4) * mrSges(5,2) + mrSges(5,3) * t642 - qJD(4) * t647 + t640 * t651 + t701;
t553 = t684 * t559 + t680 * t560;
t620 = -t629 * t672 + t709;
t693 = mrSges(4,3) * qJDD(2) + t686 * (-mrSges(4,1) * t675 + t715);
t551 = m(4) * t620 - t672 * t693 + t553;
t702 = -t680 * t559 + t684 * t560;
t552 = m(4) * t621 + t675 * t693 + t702;
t703 = -t551 * t672 + t675 * t552;
t543 = m(3) * t637 - mrSges(3,1) * t686 - qJDD(2) * mrSges(3,2) + t703;
t546 = t675 * t551 + t672 * t552;
t545 = m(3) * t649 + t546;
t534 = t543 * t711 - t545 * t674 + t677 * t714;
t532 = m(2) * t658 + t534;
t538 = t685 * t543 - t564 * t681;
t537 = m(2) * t659 + t538;
t710 = t676 * t532 + t673 * t537;
t697 = Ifges(4,5) * t672 + Ifges(4,6) * t675;
t708 = t686 * t697;
t533 = t543 * t712 + t677 * t545 + t674 * t714;
t704 = -t532 * t673 + t676 * t537;
t699 = Ifges(4,1) * t672 + Ifges(4,4) * t675;
t698 = Ifges(4,4) * t672 + Ifges(4,2) * t675;
t596 = Ifges(7,5) * t624 + Ifges(7,6) * t623 + Ifges(7,3) * t630;
t598 = Ifges(7,1) * t624 + Ifges(7,4) * t623 + Ifges(7,5) * t630;
t569 = -mrSges(7,1) * t580 + mrSges(7,3) * t579 + Ifges(7,4) * t593 + Ifges(7,2) * t592 + Ifges(7,6) * t602 - t596 * t624 + t598 * t630;
t597 = Ifges(7,4) * t624 + Ifges(7,2) * t623 + Ifges(7,6) * t630;
t570 = mrSges(7,2) * t580 - mrSges(7,3) * t578 + Ifges(7,1) * t593 + Ifges(7,4) * t592 + Ifges(7,5) * t602 + t596 * t623 - t597 * t630;
t613 = Ifges(6,5) * t635 + Ifges(6,6) * t634 + Ifges(6,3) * t670;
t614 = Ifges(6,4) * t635 + Ifges(6,2) * t634 + Ifges(6,6) * t670;
t554 = mrSges(6,2) * t595 - mrSges(6,3) * t582 + Ifges(6,1) * t604 + Ifges(6,4) * t603 + Ifges(6,5) * t667 - pkin(10) * t568 - t569 * t678 + t570 * t682 + t613 * t634 - t614 * t670;
t615 = Ifges(6,1) * t635 + Ifges(6,4) * t634 + Ifges(6,5) * t670;
t555 = -mrSges(6,1) * t595 - mrSges(7,1) * t578 + mrSges(7,2) * t579 + mrSges(6,3) * t583 + Ifges(6,4) * t604 - Ifges(7,5) * t593 + Ifges(6,2) * t603 + Ifges(6,6) * t667 - Ifges(7,6) * t592 - Ifges(7,3) * t602 - pkin(5) * t568 - t597 * t624 + t598 * t623 - t613 * t635 + t615 * t670;
t631 = Ifges(5,5) * t652 + Ifges(5,6) * t651 + Ifges(5,3) * qJD(4);
t633 = Ifges(5,1) * t652 + Ifges(5,4) * t651 + Ifges(5,5) * qJD(4);
t539 = -mrSges(5,1) * t622 + mrSges(5,3) * t590 + Ifges(5,4) * t643 + Ifges(5,2) * t642 + Ifges(5,6) * qJDD(4) - pkin(4) * t692 + pkin(9) * t701 + qJD(4) * t633 + t679 * t554 + t683 * t555 - t652 * t631;
t632 = Ifges(5,4) * t652 + Ifges(5,2) * t651 + Ifges(5,6) * qJD(4);
t547 = mrSges(5,2) * t622 - mrSges(5,3) * t589 + Ifges(5,1) * t643 + Ifges(5,4) * t642 + Ifges(5,5) * qJDD(4) - pkin(9) * t561 - qJD(4) * t632 + t554 * t683 - t555 * t679 + t631 * t651;
t528 = -mrSges(4,1) * t626 + mrSges(4,3) * t621 - pkin(3) * t688 + pkin(8) * t702 + qJDD(2) * t698 + t684 * t539 + t680 * t547 - t672 * t708;
t529 = mrSges(4,2) * t626 - mrSges(4,3) * t620 - pkin(8) * t553 + qJDD(2) * t699 - t680 * t539 + t684 * t547 + t675 * t708;
t527 = mrSges(3,2) * t649 - mrSges(3,3) * t636 + Ifges(3,5) * qJDD(2) - Ifges(3,6) * t686 - qJ(3) * t546 - t528 * t672 + t529 * t675;
t530 = -Ifges(5,3) * qJDD(4) - pkin(3) * t553 - pkin(2) * t546 + (-t672 * t698 + t675 * t699 + Ifges(3,5)) * t686 - pkin(10) * t700 + (Ifges(3,6) - t697) * qJDD(2) - pkin(5) * t690 - t682 * t569 - t678 * t570 - Ifges(6,3) * t667 - t652 * t632 - Ifges(5,5) * t643 - mrSges(3,1) * t649 + t651 * t633 + t634 * t615 - t635 * t614 + mrSges(3,3) * t637 - Ifges(5,6) * t642 - mrSges(4,1) * t620 + mrSges(4,2) * t621 - Ifges(6,6) * t603 - Ifges(6,5) * t604 - mrSges(5,1) * t589 + mrSges(5,2) * t590 - mrSges(6,1) * t582 + mrSges(6,2) * t583 - pkin(4) * t561;
t691 = pkin(7) * t538 + t527 * t681 + t530 * t685;
t526 = Ifges(3,3) * qJDD(2) + mrSges(3,1) * t636 - mrSges(3,2) * t637 + t672 * t529 + t675 * t528 + pkin(2) * (-qJDD(2) * t715 + t687) + qJ(3) * t703;
t525 = mrSges(2,2) * t671 - mrSges(2,3) * t658 + t685 * t527 - t681 * t530 + (-t533 * t674 - t534 * t677) * pkin(7);
t524 = -mrSges(2,1) * t671 + mrSges(2,3) * t659 - pkin(1) * t533 - t674 * t526 + t677 * t691;
t1 = [-m(1) * g(1) + t704; -m(1) * g(2) + t710; -m(1) * g(3) + m(2) * t671 + t533; -mrSges(1,2) * g(3) + mrSges(1,3) * g(2) - qJ(1) * t710 - t673 * t524 + t676 * t525; mrSges(1,1) * g(3) - mrSges(1,3) * g(1) + qJ(1) * t704 + t676 * t524 + t673 * t525; -mrSges(1,1) * g(2) + mrSges(2,1) * t658 + mrSges(1,2) * g(1) - mrSges(2,2) * t659 + pkin(1) * t534 + t677 * t526 + t674 * t691;];
tauB  = t1;
