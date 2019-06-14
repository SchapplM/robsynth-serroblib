% Calculate vector of inverse dynamics base forces with Newton-Euler for
% S6RPRPPR6
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
%   pkin=[a2,a3,a4,a5,a6,d1,d3,d6,theta4,theta5]';
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
% Datum: 2019-05-05 17:10
% Revision: 89c353f7eff3bd693eda4e29f35b2761dbc3ada0 (2019-05-03)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauB = S6RPRPPR6_invdynB_fixb_snew_vp2(qJ, qJD, qJDD, g, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(6,1),zeros(3,1),zeros(10,1),zeros(7,1),zeros(7,3),zeros(7,6)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPRPPR6_invdynB_fixb_snew_vp2: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RPRPPR6_invdynB_fixb_snew_vp2: qJD has to be [6x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [6 1]), ...
  'S6RPRPPR6_invdynB_fixb_snew_vp2: qJDD has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6RPRPPR6_invdynB_fixb_snew_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S6RPRPPR6_invdynB_fixb_snew_vp2: pkin has to be [10x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6RPRPPR6_invdynB_fixb_snew_vp2: m has to be [7x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [7,3]), ...
  'S6RPRPPR6_invdynB_fixb_snew_vp2: mrSges has to be [7x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [7 6]), ...
  'S6RPRPPR6_invdynB_fixb_snew_vp2: Ifges has to be [7x6] (double)'); 

%% Symbolic Calculation
% From invdyn_fixb_NewtonEuler_linkframe_tauB_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-05-05 17:07:54
% EndTime: 2019-05-05 17:08:00
% DurationCPUTime: 6.69s
% Computational Cost: add. (95555->337), mult. (216289->417), div. (0->0), fcn. (143658->10), ass. (0->129)
t677 = -2 * qJD(4);
t645 = sin(qJ(1));
t648 = cos(qJ(1));
t628 = t645 * g(1) - t648 * g(2);
t650 = qJD(1) ^ 2;
t657 = -t650 * qJ(2) + qJDD(2) - t628;
t676 = -pkin(1) - pkin(7);
t606 = qJDD(1) * t676 + t657;
t644 = sin(qJ(3));
t647 = cos(qJ(3));
t593 = t644 * g(3) + t647 * t606;
t668 = qJD(1) * qJD(3);
t666 = t644 * t668;
t624 = qJDD(1) * t647 - t666;
t569 = (-t624 - t666) * qJ(4) + (-t644 * t647 * t650 + qJDD(3)) * pkin(3) + t593;
t594 = -g(3) * t647 + t644 * t606;
t623 = -qJDD(1) * t644 - t647 * t668;
t670 = qJD(1) * t647;
t626 = qJD(3) * pkin(3) - qJ(4) * t670;
t636 = t644 ^ 2;
t570 = -pkin(3) * t636 * t650 + qJ(4) * t623 - qJD(3) * t626 + t594;
t640 = sin(pkin(9));
t642 = cos(pkin(9));
t671 = qJD(1) * t644;
t613 = -t640 * t671 + t642 * t670;
t553 = t569 * t642 - t640 * t570 + t613 * t677;
t629 = -t648 * g(1) - t645 * g(2);
t658 = qJDD(1) * qJ(2) + 0.2e1 * qJD(2) * qJD(1) + t629;
t675 = mrSges(2,1) - mrSges(3,2);
t674 = Ifges(2,5) - Ifges(3,4);
t673 = -Ifges(2,6) + Ifges(3,5);
t612 = (t640 * t647 + t642 * t644) * qJD(1);
t554 = t640 * t569 + t642 * t570 + t612 * t677;
t588 = mrSges(5,1) * t612 + mrSges(5,2) * t613;
t591 = -t642 * t623 + t624 * t640;
t605 = qJD(3) * mrSges(5,1) - mrSges(5,3) * t613;
t587 = pkin(4) * t612 - qJ(5) * t613;
t649 = qJD(3) ^ 2;
t546 = -pkin(4) * t649 + qJDD(3) * qJ(5) - t587 * t612 + t554;
t574 = -t623 * pkin(3) + qJDD(4) + t626 * t670 + (-qJ(4) * t636 + t676) * t650 + t658;
t592 = t623 * t640 + t624 * t642;
t549 = (qJD(3) * t612 - t592) * qJ(5) + (qJD(3) * t613 + t591) * pkin(4) + t574;
t639 = sin(pkin(10));
t641 = cos(pkin(10));
t599 = qJD(3) * t639 + t613 * t641;
t541 = -0.2e1 * qJD(5) * t599 - t639 * t546 + t641 * t549;
t582 = qJDD(3) * t639 + t592 * t641;
t598 = qJD(3) * t641 - t613 * t639;
t539 = (t598 * t612 - t582) * pkin(8) + (t598 * t599 + t591) * pkin(5) + t541;
t542 = 0.2e1 * qJD(5) * t598 + t641 * t546 + t639 * t549;
t579 = pkin(5) * t612 - pkin(8) * t599;
t581 = qJDD(3) * t641 - t592 * t639;
t597 = t598 ^ 2;
t540 = -pkin(5) * t597 + pkin(8) * t581 - t579 * t612 + t542;
t643 = sin(qJ(6));
t646 = cos(qJ(6));
t537 = t539 * t646 - t540 * t643;
t572 = t598 * t646 - t599 * t643;
t552 = qJD(6) * t572 + t581 * t643 + t582 * t646;
t573 = t598 * t643 + t599 * t646;
t559 = -mrSges(7,1) * t572 + mrSges(7,2) * t573;
t610 = qJD(6) + t612;
t560 = -mrSges(7,2) * t610 + mrSges(7,3) * t572;
t590 = qJDD(6) + t591;
t535 = m(7) * t537 + mrSges(7,1) * t590 - mrSges(7,3) * t552 - t559 * t573 + t560 * t610;
t538 = t539 * t643 + t540 * t646;
t551 = -qJD(6) * t573 + t581 * t646 - t582 * t643;
t561 = mrSges(7,1) * t610 - mrSges(7,3) * t573;
t536 = m(7) * t538 - mrSges(7,2) * t590 + mrSges(7,3) * t551 + t559 * t572 - t561 * t610;
t527 = t646 * t535 + t643 * t536;
t575 = -mrSges(6,1) * t598 + mrSges(6,2) * t599;
t577 = -mrSges(6,2) * t612 + mrSges(6,3) * t598;
t525 = m(6) * t541 + mrSges(6,1) * t591 - mrSges(6,3) * t582 - t575 * t599 + t577 * t612 + t527;
t578 = mrSges(6,1) * t612 - mrSges(6,3) * t599;
t661 = -t535 * t643 + t646 * t536;
t526 = m(6) * t542 - mrSges(6,2) * t591 + mrSges(6,3) * t581 + t575 * t598 - t578 * t612 + t661;
t662 = -t525 * t639 + t641 * t526;
t520 = m(5) * t554 - qJDD(3) * mrSges(5,2) - mrSges(5,3) * t591 - qJD(3) * t605 - t588 * t612 + t662;
t604 = -qJD(3) * mrSges(5,2) - mrSges(5,3) * t612;
t545 = -qJDD(3) * pkin(4) - qJ(5) * t649 + t613 * t587 + qJDD(5) - t553;
t543 = -pkin(5) * t581 - pkin(8) * t597 + t579 * t599 + t545;
t655 = m(7) * t543 - t551 * mrSges(7,1) + mrSges(7,2) * t552 - t572 * t560 + t561 * t573;
t652 = -m(6) * t545 + t581 * mrSges(6,1) - mrSges(6,2) * t582 + t598 * t577 - t578 * t599 - t655;
t531 = m(5) * t553 + qJDD(3) * mrSges(5,1) - mrSges(5,3) * t592 + qJD(3) * t604 - t588 * t613 + t652;
t512 = t640 * t520 + t642 * t531;
t622 = (mrSges(4,1) * t644 + mrSges(4,2) * t647) * qJD(1);
t625 = -qJD(3) * mrSges(4,2) - mrSges(4,3) * t671;
t510 = m(4) * t593 + qJDD(3) * mrSges(4,1) - mrSges(4,3) * t624 + qJD(3) * t625 - t622 * t670 + t512;
t627 = qJD(3) * mrSges(4,1) - mrSges(4,3) * t670;
t663 = t642 * t520 - t531 * t640;
t511 = m(4) * t594 - qJDD(3) * mrSges(4,2) + mrSges(4,3) * t623 - qJD(3) * t627 - t622 * t671 + t663;
t507 = t647 * t510 + t644 * t511;
t611 = -qJDD(1) * pkin(1) + t657;
t656 = -m(3) * t611 + t650 * mrSges(3,3) - t507;
t505 = m(2) * t628 - t650 * mrSges(2,2) + qJDD(1) * t675 + t656;
t607 = t650 * pkin(1) - t658;
t603 = t650 * t676 + t658;
t521 = t641 * t525 + t639 * t526;
t654 = m(5) * t574 + mrSges(5,1) * t591 + t592 * mrSges(5,2) + t604 * t612 + t613 * t605 + t521;
t653 = -m(4) * t603 + mrSges(4,1) * t623 - t624 * mrSges(4,2) - t625 * t671 - t627 * t670 - t654;
t651 = -m(3) * t607 + t650 * mrSges(3,2) + qJDD(1) * mrSges(3,3) - t653;
t517 = m(2) * t629 - mrSges(2,1) * t650 - qJDD(1) * mrSges(2,2) + t651;
t672 = t648 * t505 + t645 * t517;
t665 = -t505 * t645 + t648 * t517;
t664 = -t644 * t510 + t647 * t511;
t616 = Ifges(4,5) * qJD(3) + (Ifges(4,1) * t647 - Ifges(4,4) * t644) * qJD(1);
t615 = Ifges(4,6) * qJD(3) + (Ifges(4,4) * t647 - Ifges(4,2) * t644) * qJD(1);
t614 = Ifges(4,3) * qJD(3) + (Ifges(4,5) * t647 - Ifges(4,6) * t644) * qJD(1);
t585 = Ifges(5,1) * t613 - Ifges(5,4) * t612 + Ifges(5,5) * qJD(3);
t584 = Ifges(5,4) * t613 - Ifges(5,2) * t612 + Ifges(5,6) * qJD(3);
t583 = Ifges(5,5) * t613 - Ifges(5,6) * t612 + Ifges(5,3) * qJD(3);
t564 = Ifges(6,1) * t599 + Ifges(6,4) * t598 + Ifges(6,5) * t612;
t563 = Ifges(6,4) * t599 + Ifges(6,2) * t598 + Ifges(6,6) * t612;
t562 = Ifges(6,5) * t599 + Ifges(6,6) * t598 + Ifges(6,3) * t612;
t557 = Ifges(7,1) * t573 + Ifges(7,4) * t572 + Ifges(7,5) * t610;
t556 = Ifges(7,4) * t573 + Ifges(7,2) * t572 + Ifges(7,6) * t610;
t555 = Ifges(7,5) * t573 + Ifges(7,6) * t572 + Ifges(7,3) * t610;
t529 = mrSges(7,2) * t543 - mrSges(7,3) * t537 + Ifges(7,1) * t552 + Ifges(7,4) * t551 + Ifges(7,5) * t590 + t555 * t572 - t556 * t610;
t528 = -mrSges(7,1) * t543 + mrSges(7,3) * t538 + Ifges(7,4) * t552 + Ifges(7,2) * t551 + Ifges(7,6) * t590 - t555 * t573 + t557 * t610;
t514 = mrSges(6,2) * t545 - mrSges(6,3) * t541 + Ifges(6,1) * t582 + Ifges(6,4) * t581 + Ifges(6,5) * t591 - pkin(8) * t527 - t528 * t643 + t529 * t646 + t562 * t598 - t563 * t612;
t513 = -mrSges(6,1) * t545 + mrSges(6,3) * t542 + Ifges(6,4) * t582 + Ifges(6,2) * t581 + Ifges(6,6) * t591 - pkin(5) * t655 + pkin(8) * t661 + t646 * t528 + t643 * t529 - t599 * t562 + t612 * t564;
t508 = Ifges(5,4) * t592 + Ifges(5,6) * qJDD(3) - t613 * t583 + qJD(3) * t585 - mrSges(5,1) * t574 + mrSges(5,3) * t554 - Ifges(6,5) * t582 - Ifges(6,6) * t581 - t599 * t563 + t598 * t564 - mrSges(6,1) * t541 + mrSges(6,2) * t542 - Ifges(7,5) * t552 - Ifges(7,6) * t551 - Ifges(7,3) * t590 - t573 * t556 + t572 * t557 - mrSges(7,1) * t537 + mrSges(7,2) * t538 - pkin(5) * t527 - pkin(4) * t521 + (-Ifges(5,2) - Ifges(6,3)) * t591;
t506 = -m(3) * g(3) + t664;
t503 = mrSges(5,2) * t574 - mrSges(5,3) * t553 + Ifges(5,1) * t592 - Ifges(5,4) * t591 + Ifges(5,5) * qJDD(3) - qJ(5) * t521 - qJD(3) * t584 - t513 * t639 + t514 * t641 - t583 * t612;
t502 = mrSges(4,2) * t603 - mrSges(4,3) * t593 + Ifges(4,1) * t624 + Ifges(4,4) * t623 + Ifges(4,5) * qJDD(3) - qJ(4) * t512 - qJD(3) * t615 + t503 * t642 - t508 * t640 - t614 * t671;
t501 = -mrSges(4,1) * t603 + mrSges(4,3) * t594 + Ifges(4,4) * t624 + Ifges(4,2) * t623 + Ifges(4,6) * qJDD(3) - pkin(3) * t654 + qJ(4) * t663 + qJD(3) * t616 + t640 * t503 + t642 * t508 - t614 * t670;
t500 = pkin(2) * t507 - qJ(2) * t506 + t641 * t513 + qJ(5) * t662 + t639 * t514 + Ifges(4,6) * t623 + Ifges(4,5) * t624 - mrSges(2,3) * t628 + mrSges(3,1) * t611 + t612 * t585 + t613 * t584 + Ifges(5,5) * t592 + mrSges(4,1) * t593 - mrSges(4,2) * t594 + pkin(4) * t652 - Ifges(5,6) * t591 + mrSges(5,1) * t553 - mrSges(5,2) * t554 + pkin(3) * t512 + t673 * t650 + (Ifges(4,3) + Ifges(5,3)) * qJDD(3) + t674 * qJDD(1) + (t615 * t647 + t616 * t644) * qJD(1) + (-mrSges(2,2) + mrSges(3,3)) * g(3);
t499 = -mrSges(3,1) * t607 + mrSges(2,3) * t629 - pkin(1) * t506 - pkin(2) * t653 - pkin(7) * t664 + g(3) * t675 - qJDD(1) * t673 - t647 * t501 - t644 * t502 + t650 * t674;
t1 = [-m(1) * g(1) + t665; -m(1) * g(2) + t672; (-m(1) - m(2) - m(3)) * g(3) + t664; -mrSges(1,2) * g(3) + mrSges(1,3) * g(2) - pkin(6) * t672 - t645 * t499 + t648 * t500; mrSges(1,1) * g(3) - mrSges(1,3) * g(1) + pkin(6) * t665 + t648 * t499 + t645 * t500; pkin(1) * t656 + qJ(2) * t651 + t647 * t502 - t644 * t501 - pkin(7) * t507 + mrSges(2,1) * t628 - mrSges(2,2) * t629 + mrSges(3,2) * t611 - mrSges(3,3) * t607 - mrSges(1,1) * g(2) + mrSges(1,2) * g(1) + (-mrSges(3,2) * pkin(1) + Ifges(3,1) + Ifges(2,3)) * qJDD(1);];
tauB  = t1;
