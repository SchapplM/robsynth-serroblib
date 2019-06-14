% Calculate vector of inverse dynamics base forces with Newton-Euler for
% S6RPPRPR7
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
% pkin [10x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d4,d6,theta3,theta5]';
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
% Datum: 2019-05-05 14:33
% Revision: 89c353f7eff3bd693eda4e29f35b2761dbc3ada0 (2019-05-03)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauB = S6RPPRPR7_invdynB_fixb_snew_vp2(qJ, qJD, qJDD, g, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(6,1),zeros(3,1),zeros(10,1),zeros(7,1),zeros(7,3),zeros(7,6)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPPRPR7_invdynB_fixb_snew_vp2: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RPPRPR7_invdynB_fixb_snew_vp2: qJD has to be [6x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [6 1]), ...
  'S6RPPRPR7_invdynB_fixb_snew_vp2: qJDD has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6RPPRPR7_invdynB_fixb_snew_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S6RPPRPR7_invdynB_fixb_snew_vp2: pkin has to be [10x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6RPPRPR7_invdynB_fixb_snew_vp2: m has to be [7x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [7,3]), ...
  'S6RPPRPR7_invdynB_fixb_snew_vp2: mrSges has to be [7x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [7 6]), ...
  'S6RPPRPR7_invdynB_fixb_snew_vp2: Ifges has to be [7x6] (double)'); 

%% Symbolic Calculation
% From invdyn_fixb_NewtonEuler_linkframe_tauB_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-05-05 14:31:40
% EndTime: 2019-05-05 14:31:47
% DurationCPUTime: 6.43s
% Computational Cost: add. (87404->313), mult. (205421->383), div. (0->0), fcn. (144602->10), ass. (0->130)
t642 = sin(qJ(1));
t645 = cos(qJ(1));
t620 = t642 * g(1) - t645 * g(2);
t647 = qJD(1) ^ 2;
t655 = -t647 * qJ(2) + qJDD(2) - t620;
t678 = -pkin(1) - qJ(3);
t684 = -(2 * qJD(1) * qJD(3)) + t678 * qJDD(1) + t655;
t637 = sin(pkin(9));
t630 = t637 ^ 2;
t639 = cos(pkin(9));
t675 = t639 ^ 2 + t630;
t668 = t675 * mrSges(4,3);
t621 = -t645 * g(1) - t642 * g(2);
t683 = qJDD(1) * qJ(2) + (2 * qJD(2) * qJD(1)) + t621;
t682 = pkin(3) * t647;
t681 = mrSges(2,1) - mrSges(3,2);
t680 = -Ifges(3,4) + Ifges(2,5);
t679 = -Ifges(2,6) + Ifges(3,5);
t677 = mrSges(4,2) * t639;
t597 = t637 * g(3) + t684 * t639;
t579 = (-pkin(7) * qJDD(1) - t637 * t682) * t639 + t597;
t598 = -t639 * g(3) + t684 * t637;
t670 = qJDD(1) * t637;
t580 = -pkin(7) * t670 - t630 * t682 + t598;
t641 = sin(qJ(4));
t644 = cos(qJ(4));
t565 = t641 * t579 + t644 * t580;
t659 = t637 * t644 + t639 * t641;
t616 = t659 * qJD(1);
t658 = -t637 * t641 + t639 * t644;
t617 = t658 * qJD(1);
t594 = t616 * mrSges(5,1) + t617 * mrSges(5,2);
t672 = t617 * qJD(4);
t599 = qJDD(1) * t659 + t672;
t611 = qJD(4) * mrSges(5,1) - t617 * mrSges(5,3);
t593 = t616 * pkin(4) - t617 * qJ(5);
t646 = qJD(4) ^ 2;
t552 = -t646 * pkin(4) + qJDD(4) * qJ(5) - t616 * t593 + t565;
t654 = qJDD(3) + t683;
t588 = pkin(3) * t670 + (-t675 * pkin(7) + t678) * t647 + t654;
t673 = t616 * qJD(4);
t600 = qJDD(1) * t658 - t673;
t558 = (-t600 + t673) * qJ(5) + (t599 + t672) * pkin(4) + t588;
t636 = sin(pkin(10));
t638 = cos(pkin(10));
t605 = t636 * qJD(4) + t638 * t617;
t547 = -0.2e1 * qJD(5) * t605 - t636 * t552 + t638 * t558;
t586 = t636 * qJDD(4) + t638 * t600;
t604 = t638 * qJD(4) - t636 * t617;
t545 = (t604 * t616 - t586) * pkin(8) + (t604 * t605 + t599) * pkin(5) + t547;
t548 = 0.2e1 * qJD(5) * t604 + t638 * t552 + t636 * t558;
t584 = t616 * pkin(5) - t605 * pkin(8);
t585 = t638 * qJDD(4) - t636 * t600;
t603 = t604 ^ 2;
t546 = -t603 * pkin(5) + t585 * pkin(8) - t616 * t584 + t548;
t640 = sin(qJ(6));
t643 = cos(qJ(6));
t543 = t643 * t545 - t640 * t546;
t572 = t643 * t604 - t640 * t605;
t557 = t572 * qJD(6) + t640 * t585 + t643 * t586;
t573 = t640 * t604 + t643 * t605;
t563 = -t572 * mrSges(7,1) + t573 * mrSges(7,2);
t614 = qJD(6) + t616;
t566 = -t614 * mrSges(7,2) + t572 * mrSges(7,3);
t596 = qJDD(6) + t599;
t541 = m(7) * t543 + t596 * mrSges(7,1) - t557 * mrSges(7,3) - t573 * t563 + t614 * t566;
t544 = t640 * t545 + t643 * t546;
t556 = -t573 * qJD(6) + t643 * t585 - t640 * t586;
t567 = t614 * mrSges(7,1) - t573 * mrSges(7,3);
t542 = m(7) * t544 - t596 * mrSges(7,2) + t556 * mrSges(7,3) + t572 * t563 - t614 * t567;
t533 = t643 * t541 + t640 * t542;
t574 = -t604 * mrSges(6,1) + t605 * mrSges(6,2);
t582 = -t616 * mrSges(6,2) + t604 * mrSges(6,3);
t531 = m(6) * t547 + t599 * mrSges(6,1) - t586 * mrSges(6,3) - t605 * t574 + t616 * t582 + t533;
t583 = t616 * mrSges(6,1) - t605 * mrSges(6,3);
t663 = -t640 * t541 + t643 * t542;
t532 = m(6) * t548 - t599 * mrSges(6,2) + t585 * mrSges(6,3) + t604 * t574 - t616 * t583 + t663;
t664 = -t636 * t531 + t638 * t532;
t526 = m(5) * t565 - qJDD(4) * mrSges(5,2) - t599 * mrSges(5,3) - qJD(4) * t611 - t616 * t594 + t664;
t564 = t644 * t579 - t641 * t580;
t610 = -qJD(4) * mrSges(5,2) - t616 * mrSges(5,3);
t551 = -qJDD(4) * pkin(4) - t646 * qJ(5) + t617 * t593 + qJDD(5) - t564;
t549 = -t585 * pkin(5) - t603 * pkin(8) + t605 * t584 + t551;
t652 = m(7) * t549 - t556 * mrSges(7,1) + t557 * mrSges(7,2) - t572 * t566 + t573 * t567;
t648 = -m(6) * t551 + t585 * mrSges(6,1) - t586 * mrSges(6,2) + t604 * t582 - t605 * t583 - t652;
t537 = m(5) * t564 + qJDD(4) * mrSges(5,1) - t600 * mrSges(5,3) + qJD(4) * t610 - t617 * t594 + t648;
t518 = t641 * t526 + t644 * t537;
t657 = -qJDD(1) * mrSges(4,3) - t647 * (mrSges(4,1) * t637 + t677);
t516 = m(4) * t597 + t639 * t657 + t518;
t665 = t644 * t526 - t641 * t537;
t517 = m(4) * t598 + t637 * t657 + t665;
t513 = t639 * t516 + t637 * t517;
t615 = -qJDD(1) * pkin(1) + t655;
t653 = -m(3) * t615 + t647 * mrSges(3,3) - t513;
t511 = m(2) * t620 - t647 * mrSges(2,2) + t681 * qJDD(1) + t653;
t613 = t647 * pkin(1) - t683;
t609 = t678 * t647 + t654;
t527 = t638 * t531 + t636 * t532;
t651 = m(5) * t588 + t599 * mrSges(5,1) + t600 * mrSges(5,2) + t616 * t610 + t617 * t611 + t527;
t650 = m(4) * t609 + mrSges(4,1) * t670 + qJDD(1) * t677 + t651;
t649 = -m(3) * t613 + t647 * mrSges(3,2) + qJDD(1) * mrSges(3,3) + t650;
t523 = t649 + (-mrSges(2,1) - t668) * t647 - qJDD(1) * mrSges(2,2) + m(2) * t621;
t676 = t645 * t511 + t642 * t523;
t660 = Ifges(4,5) * t639 - Ifges(4,6) * t637;
t674 = t647 * t660;
t667 = -t642 * t511 + t645 * t523;
t666 = -t637 * t516 + t639 * t517;
t662 = Ifges(4,1) * t639 - Ifges(4,4) * t637;
t661 = Ifges(4,4) * t639 - Ifges(4,2) * t637;
t591 = Ifges(5,1) * t617 - Ifges(5,4) * t616 + Ifges(5,5) * qJD(4);
t590 = Ifges(5,4) * t617 - Ifges(5,2) * t616 + Ifges(5,6) * qJD(4);
t589 = Ifges(5,5) * t617 - Ifges(5,6) * t616 + Ifges(5,3) * qJD(4);
t570 = Ifges(6,1) * t605 + Ifges(6,4) * t604 + Ifges(6,5) * t616;
t569 = Ifges(6,4) * t605 + Ifges(6,2) * t604 + Ifges(6,6) * t616;
t568 = Ifges(6,5) * t605 + Ifges(6,6) * t604 + Ifges(6,3) * t616;
t561 = Ifges(7,1) * t573 + Ifges(7,4) * t572 + Ifges(7,5) * t614;
t560 = Ifges(7,4) * t573 + Ifges(7,2) * t572 + Ifges(7,6) * t614;
t559 = Ifges(7,5) * t573 + Ifges(7,6) * t572 + Ifges(7,3) * t614;
t535 = mrSges(7,2) * t549 - mrSges(7,3) * t543 + Ifges(7,1) * t557 + Ifges(7,4) * t556 + Ifges(7,5) * t596 + t572 * t559 - t614 * t560;
t534 = -mrSges(7,1) * t549 + mrSges(7,3) * t544 + Ifges(7,4) * t557 + Ifges(7,2) * t556 + Ifges(7,6) * t596 - t573 * t559 + t614 * t561;
t520 = mrSges(6,2) * t551 - mrSges(6,3) * t547 + Ifges(6,1) * t586 + Ifges(6,4) * t585 + Ifges(6,5) * t599 - pkin(8) * t533 - t640 * t534 + t643 * t535 + t604 * t568 - t616 * t569;
t519 = -mrSges(6,1) * t551 + mrSges(6,3) * t548 + Ifges(6,4) * t586 + Ifges(6,2) * t585 + Ifges(6,6) * t599 - pkin(5) * t652 + pkin(8) * t663 + t643 * t534 + t640 * t535 - t605 * t568 + t616 * t570;
t514 = Ifges(5,4) * t600 + Ifges(5,6) * qJDD(4) - t617 * t589 + qJD(4) * t591 - mrSges(5,1) * t588 + mrSges(5,3) * t565 - Ifges(6,5) * t586 - Ifges(6,6) * t585 - t605 * t569 + t604 * t570 - mrSges(6,1) * t547 + mrSges(6,2) * t548 - Ifges(7,5) * t557 - Ifges(7,6) * t556 - Ifges(7,3) * t596 - t573 * t560 + t572 * t561 - mrSges(7,1) * t543 + mrSges(7,2) * t544 - pkin(5) * t533 - pkin(4) * t527 + (-Ifges(5,2) - Ifges(6,3)) * t599;
t512 = -m(3) * g(3) + t666;
t509 = mrSges(5,2) * t588 - mrSges(5,3) * t564 + Ifges(5,1) * t600 - Ifges(5,4) * t599 + Ifges(5,5) * qJDD(4) - qJ(5) * t527 - qJD(4) * t590 - t636 * t519 + t638 * t520 - t616 * t589;
t508 = mrSges(4,2) * t609 - mrSges(4,3) * t597 - pkin(7) * t518 + qJDD(1) * t662 + t644 * t509 - t641 * t514 - t637 * t674;
t507 = -mrSges(4,1) * t609 + mrSges(4,3) * t598 - pkin(3) * t651 + pkin(7) * t665 + t661 * qJDD(1) + t641 * t509 + t644 * t514 - t639 * t674;
t506 = pkin(2) * t513 - qJ(2) * t512 + Ifges(5,3) * qJDD(4) + pkin(3) * t518 + mrSges(5,1) * t564 - mrSges(5,2) * t565 + mrSges(4,1) * t597 - mrSges(4,2) * t598 - Ifges(5,6) * t599 + Ifges(5,5) * t600 + pkin(4) * t648 + mrSges(3,1) * t615 + t616 * t591 + t617 * t590 - mrSges(2,3) * t620 + qJ(5) * t664 + t636 * t520 + t638 * t519 + (-mrSges(2,2) + mrSges(3,3)) * g(3) + (t660 + t680) * qJDD(1) + (t637 * t662 + t639 * t661 + t679) * t647;
t505 = mrSges(2,3) * t621 - mrSges(3,1) * t613 - t637 * t508 - t639 * t507 + pkin(2) * t650 - qJ(3) * t666 - pkin(1) * t512 - t679 * qJDD(1) + t681 * g(3) + (-pkin(2) * t668 + t680) * t647;
t1 = [-m(1) * g(1) + t667; -m(1) * g(2) + t676; (-m(1) - m(2) - m(3)) * g(3) + t666; -mrSges(1,2) * g(3) + mrSges(1,3) * g(2) - pkin(6) * t676 - t642 * t505 + t645 * t506; mrSges(1,1) * g(3) - mrSges(1,3) * g(1) + pkin(6) * t667 + t645 * t505 + t642 * t506; pkin(1) * t653 + qJ(2) * (-t647 * t668 + t649) + t639 * t508 - t637 * t507 - qJ(3) * t513 + mrSges(2,1) * t620 - mrSges(2,2) * t621 + mrSges(3,2) * t615 - mrSges(3,3) * t613 - mrSges(1,1) * g(2) + mrSges(1,2) * g(1) + (-pkin(1) * mrSges(3,2) + Ifges(3,1) + Ifges(2,3)) * qJDD(1);];
tauB  = t1;
