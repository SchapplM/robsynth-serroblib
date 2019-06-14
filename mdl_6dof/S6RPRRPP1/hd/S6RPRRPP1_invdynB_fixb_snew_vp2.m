% Calculate vector of inverse dynamics base forces with Newton-Euler for
% S6RPRRPP1
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
%   pkin=[a2,a3,a4,a5,a6,d1,d3,d4,theta2,theta5]';
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
% Datum: 2019-05-05 21:19
% Revision: 89c353f7eff3bd693eda4e29f35b2761dbc3ada0 (2019-05-03)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauB = S6RPRRPP1_invdynB_fixb_snew_vp2(qJ, qJD, qJDD, g, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(6,1),zeros(3,1),zeros(10,1),zeros(7,1),zeros(7,3),zeros(7,6)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPRRPP1_invdynB_fixb_snew_vp2: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RPRRPP1_invdynB_fixb_snew_vp2: qJD has to be [6x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [6 1]), ...
  'S6RPRRPP1_invdynB_fixb_snew_vp2: qJDD has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6RPRRPP1_invdynB_fixb_snew_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S6RPRRPP1_invdynB_fixb_snew_vp2: pkin has to be [10x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6RPRRPP1_invdynB_fixb_snew_vp2: m has to be [7x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [7,3]), ...
  'S6RPRRPP1_invdynB_fixb_snew_vp2: mrSges has to be [7x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [7 6]), ...
  'S6RPRRPP1_invdynB_fixb_snew_vp2: Ifges has to be [7x6] (double)'); 

%% Symbolic Calculation
% From invdyn_fixb_NewtonEuler_linkframe_tauB_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-05-05 21:15:47
% EndTime: 2019-05-05 21:15:56
% DurationCPUTime: 6.69s
% Computational Cost: add. (101162->319), mult. (200161->387), div. (0->0), fcn. (129294->10), ass. (0->126)
t664 = Ifges(6,1) + Ifges(7,1);
t659 = Ifges(6,4) - Ifges(7,5);
t658 = Ifges(6,5) + Ifges(7,4);
t663 = Ifges(6,2) + Ifges(7,3);
t662 = -Ifges(7,2) - Ifges(6,3);
t657 = Ifges(6,6) - Ifges(7,6);
t661 = -2 * qJD(5);
t660 = -mrSges(6,3) - mrSges(7,2);
t656 = cos(pkin(10));
t628 = sin(qJ(1));
t631 = cos(qJ(1));
t615 = t628 * g(1) - g(2) * t631;
t607 = qJDD(1) * pkin(1) + t615;
t616 = -g(1) * t631 - g(2) * t628;
t633 = qJD(1) ^ 2;
t609 = -pkin(1) * t633 + t616;
t624 = sin(pkin(9));
t625 = cos(pkin(9));
t586 = t624 * t607 + t625 * t609;
t573 = -pkin(2) * t633 + qJDD(1) * pkin(7) + t586;
t622 = -g(3) + qJDD(2);
t627 = sin(qJ(3));
t630 = cos(qJ(3));
t564 = t630 * t573 + t627 * t622;
t608 = (-mrSges(4,1) * t630 + mrSges(4,2) * t627) * qJD(1);
t648 = qJD(1) * qJD(3);
t645 = t627 * t648;
t612 = qJDD(1) * t630 - t645;
t650 = qJD(1) * t627;
t613 = qJD(3) * mrSges(4,1) - mrSges(4,3) * t650;
t585 = t625 * t607 - t624 * t609;
t572 = -qJDD(1) * pkin(2) - t633 * pkin(7) - t585;
t644 = t630 * t648;
t611 = qJDD(1) * t627 + t644;
t549 = (-t611 - t644) * pkin(8) + (-t612 + t645) * pkin(3) + t572;
t610 = (-pkin(3) * t630 - pkin(8) * t627) * qJD(1);
t632 = qJD(3) ^ 2;
t649 = qJD(1) * t630;
t557 = -pkin(3) * t632 + qJDD(3) * pkin(8) + t610 * t649 + t564;
t626 = sin(qJ(4));
t629 = cos(qJ(4));
t536 = t629 * t549 - t626 * t557;
t605 = qJD(3) * t629 - t626 * t650;
t582 = qJD(4) * t605 + qJDD(3) * t626 + t611 * t629;
t604 = qJDD(4) - t612;
t606 = qJD(3) * t626 + t629 * t650;
t618 = qJD(4) - t649;
t532 = (t605 * t618 - t582) * qJ(5) + (t605 * t606 + t604) * pkin(4) + t536;
t537 = t626 * t549 + t629 * t557;
t581 = -qJD(4) * t606 + qJDD(3) * t629 - t611 * t626;
t589 = pkin(4) * t618 - qJ(5) * t606;
t603 = t605 ^ 2;
t534 = -pkin(4) * t603 + qJ(5) * t581 - t589 * t618 + t537;
t623 = sin(pkin(10));
t583 = -t656 * t605 + t623 * t606;
t528 = t623 * t532 + t656 * t534 + t583 * t661;
t550 = -t656 * t581 + t623 * t582;
t584 = t623 * t605 + t656 * t606;
t566 = mrSges(6,1) * t618 - mrSges(6,3) * t584;
t558 = pkin(5) * t583 - qJ(6) * t584;
t617 = t618 ^ 2;
t525 = -pkin(5) * t617 + qJ(6) * t604 + 0.2e1 * qJD(6) * t618 - t558 * t583 + t528;
t567 = -mrSges(7,1) * t618 + mrSges(7,2) * t584;
t647 = m(7) * t525 + t604 * mrSges(7,3) + t618 * t567;
t559 = mrSges(7,1) * t583 - mrSges(7,3) * t584;
t651 = -mrSges(6,1) * t583 - mrSges(6,2) * t584 - t559;
t518 = m(6) * t528 - t604 * mrSges(6,2) + t660 * t550 - t618 * t566 + t651 * t583 + t647;
t637 = t656 * t532 - t623 * t534;
t527 = t584 * t661 + t637;
t551 = t623 * t581 + t656 * t582;
t565 = -mrSges(6,2) * t618 - mrSges(6,3) * t583;
t526 = -t604 * pkin(5) - t617 * qJ(6) + qJDD(6) + ((2 * qJD(5)) + t558) * t584 - t637;
t568 = -mrSges(7,2) * t583 + mrSges(7,3) * t618;
t638 = -m(7) * t526 + t604 * mrSges(7,1) + t618 * t568;
t520 = m(6) * t527 + t604 * mrSges(6,1) + t660 * t551 + t618 * t565 + t651 * t584 + t638;
t513 = t623 * t518 + t656 * t520;
t587 = -mrSges(5,1) * t605 + mrSges(5,2) * t606;
t588 = -mrSges(5,2) * t618 + mrSges(5,3) * t605;
t511 = m(5) * t536 + mrSges(5,1) * t604 - mrSges(5,3) * t582 - t587 * t606 + t588 * t618 + t513;
t590 = mrSges(5,1) * t618 - mrSges(5,3) * t606;
t639 = t656 * t518 - t520 * t623;
t512 = m(5) * t537 - mrSges(5,2) * t604 + mrSges(5,3) * t581 + t587 * t605 - t590 * t618 + t639;
t640 = -t511 * t626 + t629 * t512;
t508 = m(4) * t564 - qJDD(3) * mrSges(4,2) + mrSges(4,3) * t612 - qJD(3) * t613 + t608 * t649 + t640;
t563 = -t627 * t573 + t630 * t622;
t614 = -qJD(3) * mrSges(4,2) + mrSges(4,3) * t649;
t556 = -qJDD(3) * pkin(3) - t632 * pkin(8) + t610 * t650 - t563;
t535 = -t581 * pkin(4) - t603 * qJ(5) + t606 * t589 + qJDD(5) + t556;
t530 = -0.2e1 * qJD(6) * t584 + (t583 * t618 - t551) * qJ(6) + (t584 * t618 + t550) * pkin(5) + t535;
t523 = m(7) * t530 + t550 * mrSges(7,1) - t551 * mrSges(7,3) - t584 * t567 + t583 * t568;
t636 = m(6) * t535 + t550 * mrSges(6,1) + t551 * mrSges(6,2) + t583 * t565 + t584 * t566 + t523;
t634 = -m(5) * t556 + t581 * mrSges(5,1) - t582 * mrSges(5,2) + t605 * t588 - t606 * t590 - t636;
t522 = m(4) * t563 + qJDD(3) * mrSges(4,1) - t611 * mrSges(4,3) + qJD(3) * t614 - t608 * t650 + t634;
t641 = t630 * t508 - t522 * t627;
t502 = m(3) * t586 - mrSges(3,1) * t633 - qJDD(1) * mrSges(3,2) + t641;
t509 = t511 * t629 + t512 * t626;
t635 = -m(4) * t572 + t612 * mrSges(4,1) - mrSges(4,2) * t611 - t613 * t650 + t614 * t649 - t509;
t505 = m(3) * t585 + qJDD(1) * mrSges(3,1) - mrSges(3,2) * t633 + t635;
t496 = t624 * t502 + t625 * t505;
t494 = m(2) * t615 + qJDD(1) * mrSges(2,1) - mrSges(2,2) * t633 + t496;
t642 = t625 * t502 - t505 * t624;
t495 = m(2) * t616 - mrSges(2,1) * t633 - qJDD(1) * mrSges(2,2) + t642;
t655 = t631 * t494 + t628 * t495;
t503 = t627 * t508 + t630 * t522;
t654 = t663 * t583 - t659 * t584 - t657 * t618;
t653 = t657 * t583 - t658 * t584 + t662 * t618;
t652 = -t659 * t583 + t664 * t584 + t658 * t618;
t646 = m(3) * t622 + t503;
t643 = -t494 * t628 + t631 * t495;
t598 = Ifges(4,5) * qJD(3) + (Ifges(4,1) * t627 + Ifges(4,4) * t630) * qJD(1);
t597 = Ifges(4,6) * qJD(3) + (Ifges(4,4) * t627 + Ifges(4,2) * t630) * qJD(1);
t596 = Ifges(4,3) * qJD(3) + (Ifges(4,5) * t627 + Ifges(4,6) * t630) * qJD(1);
t576 = Ifges(5,1) * t606 + Ifges(5,4) * t605 + Ifges(5,5) * t618;
t575 = Ifges(5,4) * t606 + Ifges(5,2) * t605 + Ifges(5,6) * t618;
t574 = Ifges(5,5) * t606 + Ifges(5,6) * t605 + Ifges(5,3) * t618;
t515 = mrSges(6,2) * t535 + mrSges(7,2) * t526 - mrSges(6,3) * t527 - mrSges(7,3) * t530 - qJ(6) * t523 - t659 * t550 + t664 * t551 + t653 * t583 + t658 * t604 + t654 * t618;
t514 = -mrSges(6,1) * t535 - mrSges(7,1) * t530 + mrSges(7,2) * t525 + mrSges(6,3) * t528 - pkin(5) * t523 - t663 * t550 + t659 * t551 + t653 * t584 + t657 * t604 + t652 * t618;
t499 = mrSges(5,2) * t556 - mrSges(5,3) * t536 + Ifges(5,1) * t582 + Ifges(5,4) * t581 + Ifges(5,5) * t604 - qJ(5) * t513 - t623 * t514 + t656 * t515 + t605 * t574 - t618 * t575;
t498 = -mrSges(5,1) * t556 + mrSges(5,3) * t537 + Ifges(5,4) * t582 + Ifges(5,2) * t581 + Ifges(5,6) * t604 - pkin(4) * t636 + qJ(5) * t639 + t656 * t514 + t623 * t515 - t606 * t574 + t618 * t576;
t497 = (-Ifges(5,3) + t662) * t604 + Ifges(4,6) * qJDD(3) - pkin(5) * t638 + Ifges(4,4) * t611 + Ifges(4,2) * t612 - t606 * t575 + t605 * t576 + qJD(3) * t598 - Ifges(5,6) * t581 - Ifges(5,5) * t582 + mrSges(4,3) * t564 - mrSges(4,1) * t572 - mrSges(5,1) * t536 + mrSges(5,2) * t537 - mrSges(7,3) * t525 + mrSges(7,1) * t526 - mrSges(6,1) * t527 + mrSges(6,2) * t528 - t596 * t650 + (qJ(6) * t559 - t652) * t583 + (pkin(5) * t559 + t654) * t584 - pkin(4) * t513 - pkin(3) * t509 + (mrSges(7,2) * qJ(6) + t657) * t550 + (mrSges(7,2) * pkin(5) - t658) * t551 - qJ(6) * t647;
t490 = mrSges(4,2) * t572 - mrSges(4,3) * t563 + Ifges(4,1) * t611 + Ifges(4,4) * t612 + Ifges(4,5) * qJDD(3) - pkin(8) * t509 - qJD(3) * t597 - t498 * t626 + t499 * t629 + t596 * t649;
t489 = Ifges(3,6) * qJDD(1) + t633 * Ifges(3,5) - mrSges(3,1) * t622 + mrSges(3,3) * t586 - Ifges(4,5) * t611 - Ifges(4,6) * t612 - Ifges(4,3) * qJDD(3) - mrSges(4,1) * t563 + mrSges(4,2) * t564 - t626 * t499 - t629 * t498 - pkin(3) * t634 - pkin(8) * t640 - pkin(2) * t503 + (-t597 * t627 + t598 * t630) * qJD(1);
t488 = mrSges(3,2) * t622 - mrSges(3,3) * t585 + Ifges(3,5) * qJDD(1) - Ifges(3,6) * t633 - pkin(7) * t503 + t490 * t630 - t497 * t627;
t487 = -mrSges(2,2) * g(3) - mrSges(2,3) * t615 + Ifges(2,5) * qJDD(1) - Ifges(2,6) * t633 - qJ(2) * t496 + t488 * t625 - t489 * t624;
t486 = mrSges(2,1) * g(3) + mrSges(2,3) * t616 + t633 * Ifges(2,5) + Ifges(2,6) * qJDD(1) - pkin(1) * t646 + qJ(2) * t642 + t624 * t488 + t625 * t489;
t1 = [-m(1) * g(1) + t643; -m(1) * g(2) + t655; (-m(1) - m(2)) * g(3) + t646; -mrSges(1,2) * g(3) + mrSges(1,3) * g(2) - pkin(6) * t655 - t628 * t486 + t631 * t487; mrSges(1,1) * g(3) - mrSges(1,3) * g(1) + pkin(6) * t643 + t631 * t486 + t628 * t487; pkin(1) * t496 + mrSges(2,1) * t615 - mrSges(2,2) * t616 + t627 * t490 + t630 * t497 + pkin(2) * t635 + pkin(7) * t641 + mrSges(3,1) * t585 - mrSges(3,2) * t586 - mrSges(1,1) * g(2) + mrSges(1,2) * g(1) + (Ifges(2,3) + Ifges(3,3)) * qJDD(1);];
tauB  = t1;
