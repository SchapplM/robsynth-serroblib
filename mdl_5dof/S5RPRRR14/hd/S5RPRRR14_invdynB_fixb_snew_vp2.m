% Calculate vector of inverse dynamics base forces with Newton-Euler for
% S5RPRRR14
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
% pkin [11x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,alpha2,alpha3,d1,d3,d4,d5,theta2]';
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
% tauB [6x1]
%   base forces of inverse dynamics (contains inertial, gravitational coriolis and centrifugal forces)

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-31 19:22
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauB = S5RPRRR14_invdynB_fixb_snew_vp2(qJ, qJD, qJDD, g, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(5,1),zeros(3,1),zeros(11,1),zeros(6,1),zeros(6,3),zeros(6,6)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RPRRR14_invdynB_fixb_snew_vp2: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RPRRR14_invdynB_fixb_snew_vp2: qJD has to be [5x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [5 1]), ...
  'S5RPRRR14_invdynB_fixb_snew_vp2: qJDD has to be [5x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S5RPRRR14_invdynB_fixb_snew_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [11 1]), ...
  'S5RPRRR14_invdynB_fixb_snew_vp2: pkin has to be [11x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5RPRRR14_invdynB_fixb_snew_vp2: m has to be [6x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [6,3]), ...
  'S5RPRRR14_invdynB_fixb_snew_vp2: mrSges has to be [6x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [6 6]), ...
  'S5RPRRR14_invdynB_fixb_snew_vp2: Ifges has to be [6x6] (double)'); 

%% Symbolic Calculation
% From invdyn_fixb_NewtonEuler_linkframe_tauB_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 19:18:18
% EndTime: 2019-12-31 19:18:48
% DurationCPUTime: 28.48s
% Computational Cost: add. (389417->329), mult. (1213610->460), div. (0->0), fcn. (1013197->14), ass. (0->156)
t646 = sin(pkin(11));
t648 = sin(pkin(5));
t649 = cos(pkin(11));
t651 = cos(pkin(5));
t654 = sin(qJ(3));
t650 = cos(pkin(6));
t658 = cos(qJ(3));
t690 = t650 * t658;
t647 = sin(pkin(6));
t695 = t647 * t658;
t662 = (-t646 * t654 + t649 * t690) * t648 + t651 * t695;
t621 = t662 * qJD(1);
t691 = t650 * t654;
t696 = t647 * t654;
t664 = t651 * t696 + (t646 * t658 + t649 * t691) * t648;
t622 = t664 * qJD(1);
t610 = -t622 * qJD(3) + t662 * qJDD(1);
t693 = t648 * t650;
t633 = (t647 * t651 + t649 * t693) * qJD(1) * pkin(8);
t655 = sin(qJ(1));
t659 = cos(qJ(1));
t644 = -t659 * g(1) - t655 * g(2);
t660 = qJD(1) ^ 2;
t699 = qJ(2) * t648;
t637 = -t660 * pkin(1) + qJDD(1) * t699 + t644;
t702 = pkin(8) * t646;
t675 = -pkin(2) * t649 - t647 * t702;
t688 = qJD(1) * t648;
t700 = pkin(8) * qJDD(1);
t670 = qJD(1) * t675 * t688 + t650 * t700;
t643 = t655 * g(1) - t659 * g(2);
t636 = qJDD(1) * pkin(1) + t660 * t699 + t643;
t684 = qJD(2) * t688;
t692 = t649 * t651;
t694 = t648 * t649;
t676 = -g(3) * t694 + t636 * t692 - 0.2e1 * t646 * t684;
t592 = (pkin(2) * qJDD(1) + qJD(1) * t633) * t651 + (-t670 * t648 - t637) * t646 + t676;
t638 = (pkin(2) * t651 - t693 * t702) * qJD(1);
t697 = t646 * t651;
t685 = t636 * t697 + (t637 + 0.2e1 * t684) * t649;
t593 = (-qJD(1) * t638 + t647 * t700) * t651 + (-g(3) * t646 + t670 * t649) * t648 + t685;
t683 = -t651 * g(3) + qJDD(2);
t601 = (-t636 + t675 * qJDD(1) + (-t633 * t649 + t638 * t646) * qJD(1)) * t648 + t683;
t568 = -t654 * t593 + (t592 * t650 + t601 * t647) * t658;
t701 = Ifges(3,3) * t651;
t698 = t646 * t648;
t569 = t592 * t691 + t658 * t593 + t601 * t696;
t608 = -t621 * mrSges(4,1) + t622 * mrSges(4,2);
t671 = -t647 * t694 + t650 * t651;
t634 = t671 * qJD(1) + qJD(3);
t618 = t634 * mrSges(4,1) - t622 * mrSges(4,3);
t631 = t671 * qJDD(1) + qJDD(3);
t609 = -t621 * pkin(3) - t622 * pkin(9);
t630 = t634 ^ 2;
t565 = -t630 * pkin(3) + t631 * pkin(9) + t621 * t609 + t569;
t577 = -t647 * t592 + t650 * t601;
t611 = t621 * qJD(3) + t664 * qJDD(1);
t567 = (-t621 * t634 - t611) * pkin(9) + (t622 * t634 - t610) * pkin(3) + t577;
t653 = sin(qJ(4));
t657 = cos(qJ(4));
t562 = t657 * t565 + t653 * t567;
t616 = t657 * t622 + t653 * t634;
t587 = -t616 * qJD(4) - t653 * t611 + t657 * t631;
t615 = -t653 * t622 + t657 * t634;
t594 = -t615 * mrSges(5,1) + t616 * mrSges(5,2);
t620 = qJD(4) - t621;
t603 = t620 * mrSges(5,1) - t616 * mrSges(5,3);
t607 = qJDD(4) - t610;
t595 = -t615 * pkin(4) - t616 * pkin(10);
t619 = t620 ^ 2;
t559 = -t619 * pkin(4) + t607 * pkin(10) + t615 * t595 + t562;
t564 = -t631 * pkin(3) - t630 * pkin(9) + t622 * t609 - t568;
t588 = t615 * qJD(4) + t657 * t611 + t653 * t631;
t560 = (-t615 * t620 - t588) * pkin(10) + (t616 * t620 - t587) * pkin(4) + t564;
t652 = sin(qJ(5));
t656 = cos(qJ(5));
t556 = -t652 * t559 + t656 * t560;
t599 = -t652 * t616 + t656 * t620;
t572 = t599 * qJD(5) + t656 * t588 + t652 * t607;
t600 = t656 * t616 + t652 * t620;
t578 = -t599 * mrSges(6,1) + t600 * mrSges(6,2);
t612 = qJD(5) - t615;
t579 = -t612 * mrSges(6,2) + t599 * mrSges(6,3);
t585 = qJDD(5) - t587;
t554 = m(6) * t556 + t585 * mrSges(6,1) - t572 * mrSges(6,3) - t600 * t578 + t612 * t579;
t557 = t656 * t559 + t652 * t560;
t571 = -t600 * qJD(5) - t652 * t588 + t656 * t607;
t580 = t612 * mrSges(6,1) - t600 * mrSges(6,3);
t555 = m(6) * t557 - t585 * mrSges(6,2) + t571 * mrSges(6,3) + t599 * t578 - t612 * t580;
t680 = -t652 * t554 + t656 * t555;
t547 = m(5) * t562 - t607 * mrSges(5,2) + t587 * mrSges(5,3) + t615 * t594 - t620 * t603 + t680;
t561 = -t653 * t565 + t657 * t567;
t602 = -t620 * mrSges(5,2) + t615 * mrSges(5,3);
t558 = -t607 * pkin(4) - t619 * pkin(10) + t616 * t595 - t561;
t667 = -m(6) * t558 + t571 * mrSges(6,1) - t572 * mrSges(6,2) + t599 * t579 - t600 * t580;
t552 = m(5) * t561 + t607 * mrSges(5,1) - t588 * mrSges(5,3) - t616 * t594 + t620 * t602 + t667;
t681 = t657 * t547 - t653 * t552;
t538 = m(4) * t569 - t631 * mrSges(4,2) + t610 * mrSges(4,3) + t621 * t608 - t634 * t618 + t681;
t541 = t653 * t547 + t657 * t552;
t617 = -t634 * mrSges(4,2) + t621 * mrSges(4,3);
t540 = m(4) * t577 - t610 * mrSges(4,1) + t611 * mrSges(4,2) - t621 * t617 + t622 * t618 + t541;
t548 = t656 * t554 + t652 * t555;
t661 = -m(5) * t564 + t587 * mrSges(5,1) - t588 * mrSges(5,2) + t615 * t602 - t616 * t603 - t548;
t544 = m(4) * t568 + t631 * mrSges(4,1) - t611 * mrSges(4,3) - t622 * t608 + t634 * t617 + t661;
t527 = t538 * t691 - t647 * t540 + t544 * t690;
t613 = -t646 * t637 + t676;
t679 = -mrSges(3,1) * t649 + mrSges(3,2) * t646;
t635 = t679 * t688;
t673 = -mrSges(3,2) * t651 + mrSges(3,3) * t694;
t640 = t673 * qJD(1);
t674 = mrSges(3,1) * t651 - mrSges(3,3) * t698;
t523 = m(3) * t613 + t674 * qJDD(1) + (-t635 * t698 + t640 * t651) * qJD(1) + t527;
t526 = t538 * t696 + t650 * t540 + t544 * t695;
t623 = -t648 * t636 + t683;
t639 = t674 * qJD(1);
t525 = m(3) * t623 + (t679 * qJDD(1) + (t639 * t646 - t640 * t649) * qJD(1)) * t648 + t526;
t532 = t658 * t538 - t654 * t544;
t614 = -g(3) * t698 + t685;
t531 = m(3) * t614 + t673 * qJDD(1) + (t635 * t694 - t639 * t651) * qJD(1) + t532;
t513 = t523 * t692 - t648 * t525 + t531 * t697;
t511 = m(2) * t643 + qJDD(1) * mrSges(2,1) - t660 * mrSges(2,2) + t513;
t517 = -t646 * t523 + t649 * t531;
t516 = m(2) * t644 - t660 * mrSges(2,1) - qJDD(1) * mrSges(2,2) + t517;
t689 = t659 * t511 + t655 * t516;
t512 = t523 * t694 + t651 * t525 + t531 * t698;
t682 = -t655 * t511 + t659 * t516;
t678 = Ifges(3,5) * t646 + Ifges(3,6) * t649;
t573 = Ifges(6,5) * t600 + Ifges(6,6) * t599 + Ifges(6,3) * t612;
t575 = Ifges(6,1) * t600 + Ifges(6,4) * t599 + Ifges(6,5) * t612;
t549 = -mrSges(6,1) * t558 + mrSges(6,3) * t557 + Ifges(6,4) * t572 + Ifges(6,2) * t571 + Ifges(6,6) * t585 - t600 * t573 + t612 * t575;
t574 = Ifges(6,4) * t600 + Ifges(6,2) * t599 + Ifges(6,6) * t612;
t550 = mrSges(6,2) * t558 - mrSges(6,3) * t556 + Ifges(6,1) * t572 + Ifges(6,4) * t571 + Ifges(6,5) * t585 + t599 * t573 - t612 * t574;
t581 = Ifges(5,5) * t616 + Ifges(5,6) * t615 + Ifges(5,3) * t620;
t582 = Ifges(5,4) * t616 + Ifges(5,2) * t615 + Ifges(5,6) * t620;
t533 = mrSges(5,2) * t564 - mrSges(5,3) * t561 + Ifges(5,1) * t588 + Ifges(5,4) * t587 + Ifges(5,5) * t607 - pkin(10) * t548 - t652 * t549 + t656 * t550 + t615 * t581 - t620 * t582;
t583 = Ifges(5,1) * t616 + Ifges(5,4) * t615 + Ifges(5,5) * t620;
t534 = -mrSges(5,1) * t564 - mrSges(6,1) * t556 + mrSges(6,2) * t557 + mrSges(5,3) * t562 + Ifges(5,4) * t588 - Ifges(6,5) * t572 + Ifges(5,2) * t587 + Ifges(5,6) * t607 - Ifges(6,6) * t571 - Ifges(6,3) * t585 - pkin(4) * t548 - t600 * t574 + t599 * t575 - t616 * t581 + t620 * t583;
t604 = Ifges(4,5) * t622 + Ifges(4,6) * t621 + Ifges(4,3) * t634;
t605 = Ifges(4,4) * t622 + Ifges(4,2) * t621 + Ifges(4,6) * t634;
t519 = mrSges(4,2) * t577 - mrSges(4,3) * t568 + Ifges(4,1) * t611 + Ifges(4,4) * t610 + Ifges(4,5) * t631 - pkin(9) * t541 + t657 * t533 - t653 * t534 + t621 * t604 - t634 * t605;
t606 = Ifges(4,1) * t622 + Ifges(4,4) * t621 + Ifges(4,5) * t634;
t520 = Ifges(4,4) * t611 + Ifges(4,2) * t610 + Ifges(4,6) * t631 - t622 * t604 + t634 * t606 - mrSges(4,1) * t577 + mrSges(4,3) * t569 - Ifges(5,5) * t588 - Ifges(5,6) * t587 - Ifges(5,3) * t607 - t616 * t582 + t615 * t583 - mrSges(5,1) * t561 + mrSges(5,2) * t562 - t652 * t550 - t656 * t549 - pkin(4) * t667 - pkin(10) * t680 - pkin(3) * t541;
t669 = pkin(8) * t532 + t519 * t654 + t520 * t658;
t518 = mrSges(4,1) * t568 - mrSges(4,2) * t569 + Ifges(4,5) * t611 + Ifges(4,6) * t610 + Ifges(4,3) * t631 + pkin(3) * t661 + pkin(9) * t681 + t653 * t533 + t657 * t534 + t622 * t605 - t621 * t606;
t626 = (t678 * t648 + t701) * qJD(1);
t666 = Ifges(3,5) * t651 + (Ifges(3,1) * t646 + Ifges(3,4) * t649) * t648;
t628 = t666 * qJD(1);
t665 = Ifges(3,6) * t651 + (Ifges(3,4) * t646 + Ifges(3,2) * t649) * t648;
t508 = -mrSges(3,1) * t623 + mrSges(3,3) * t614 - pkin(2) * t526 - t647 * t518 + (-t626 * t698 + t628 * t651) * qJD(1) + t669 * t650 + t665 * qJDD(1);
t627 = t665 * qJD(1);
t509 = mrSges(3,2) * t623 - mrSges(3,3) * t613 + t658 * t519 - t654 * t520 + (t626 * t694 - t627 * t651) * qJD(1) + (-t526 * t647 - t527 * t650) * pkin(8) + t666 * qJDD(1);
t668 = qJ(2) * t517 + t508 * t649 + t509 * t646;
t507 = qJDD(1) * t701 + mrSges(3,1) * t613 - mrSges(3,2) * t614 + pkin(2) * t527 + t650 * t518 + t669 * t647 + (t678 * qJDD(1) + (t627 * t646 - t628 * t649) * qJD(1)) * t648;
t506 = -mrSges(2,2) * g(3) - mrSges(2,3) * t643 + Ifges(2,5) * qJDD(1) - t660 * Ifges(2,6) - t646 * t508 + t649 * t509 + (-t512 * t648 - t513 * t651) * qJ(2);
t505 = mrSges(2,1) * g(3) + mrSges(2,3) * t644 + t660 * Ifges(2,5) + Ifges(2,6) * qJDD(1) - pkin(1) * t512 - t648 * t507 + t668 * t651;
t1 = [-m(1) * g(1) + t682; -m(1) * g(2) + t689; (-m(1) - m(2)) * g(3) + t512; -mrSges(1,2) * g(3) + mrSges(1,3) * g(2) - pkin(7) * t689 - t655 * t505 + t659 * t506; mrSges(1,1) * g(3) - mrSges(1,3) * g(1) + pkin(7) * t682 + t659 * t505 + t655 * t506; -mrSges(1,1) * g(2) + mrSges(2,1) * t643 + mrSges(1,2) * g(1) - mrSges(2,2) * t644 + Ifges(2,3) * qJDD(1) + pkin(1) * t513 + t651 * t507 + t668 * t648;];
tauB = t1;
