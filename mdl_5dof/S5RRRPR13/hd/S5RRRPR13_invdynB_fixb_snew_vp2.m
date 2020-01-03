% Calculate vector of inverse dynamics base forces with Newton-Euler for
% S5RRRPR13
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
% pkin [9x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,alpha2,d1,d2,d3,d5]';
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
% Datum: 2019-12-31 21:48
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauB = S5RRRPR13_invdynB_fixb_snew_vp2(qJ, qJD, qJDD, g, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(5,1),zeros(3,1),zeros(9,1),zeros(6,1),zeros(6,3),zeros(6,6)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RRRPR13_invdynB_fixb_snew_vp2: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RRRPR13_invdynB_fixb_snew_vp2: qJD has to be [5x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [5 1]), ...
  'S5RRRPR13_invdynB_fixb_snew_vp2: qJDD has to be [5x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S5RRRPR13_invdynB_fixb_snew_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S5RRRPR13_invdynB_fixb_snew_vp2: pkin has to be [9x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5RRRPR13_invdynB_fixb_snew_vp2: m has to be [6x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [6,3]), ...
  'S5RRRPR13_invdynB_fixb_snew_vp2: mrSges has to be [6x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [6 6]), ...
  'S5RRRPR13_invdynB_fixb_snew_vp2: Ifges has to be [6x6] (double)'); 

%% Symbolic Calculation
% From invdyn_fixb_NewtonEuler_linkframe_tauB_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 21:43:49
% EndTime: 2019-12-31 21:44:00
% DurationCPUTime: 6.56s
% Computational Cost: add. (71814->307), mult. (155280->383), div. (0->0), fcn. (115420->10), ass. (0->131)
t653 = Ifges(4,1) + Ifges(5,2);
t646 = Ifges(4,4) + Ifges(5,6);
t645 = Ifges(4,5) - Ifges(5,4);
t652 = -Ifges(4,2) - Ifges(5,3);
t644 = Ifges(4,6) - Ifges(5,5);
t651 = Ifges(4,3) + Ifges(5,1);
t607 = sin(pkin(5));
t611 = sin(qJ(2));
t614 = cos(qJ(2));
t630 = qJD(1) * qJD(2);
t594 = (-qJDD(1) * t614 + t611 * t630) * t607;
t650 = -2 * qJD(4);
t649 = cos(qJ(3));
t648 = pkin(7) * t607;
t608 = cos(pkin(5));
t647 = t608 * g(3);
t604 = t608 * qJD(1) + qJD(2);
t610 = sin(qJ(3));
t632 = qJD(1) * t607;
t629 = t611 * t632;
t581 = -t649 * t604 + t610 * t629;
t631 = qJD(1) * t614;
t628 = t607 * t631;
t598 = -qJD(3) + t628;
t643 = t581 * t598;
t642 = t607 * t611;
t641 = t607 * t614;
t640 = t608 * t611;
t639 = t608 * t614;
t612 = sin(qJ(1));
t615 = cos(qJ(1));
t600 = t612 * g(1) - t615 * g(2);
t616 = qJD(1) ^ 2;
t589 = qJDD(1) * pkin(1) + t616 * t648 + t600;
t601 = -t615 * g(1) - t612 * g(2);
t590 = -t616 * pkin(1) + qJDD(1) * t648 + t601;
t633 = t589 * t640 + t614 * t590;
t559 = -g(3) * t642 + t633;
t587 = t604 * mrSges(3,1) - mrSges(3,3) * t629;
t591 = (-mrSges(3,1) * t614 + mrSges(3,2) * t611) * t632;
t603 = t608 * qJDD(1) + qJDD(2);
t592 = (-pkin(2) * t614 - pkin(8) * t611) * t632;
t602 = t604 ^ 2;
t537 = -t602 * pkin(2) + t603 * pkin(8) + (-g(3) * t611 + t592 * t631) * t607 + t633;
t593 = (qJDD(1) * t611 + t614 * t630) * t607;
t538 = t594 * pkin(2) - t593 * pkin(8) - t647 + (-t589 + (pkin(2) * t611 - pkin(8) * t614) * t604 * qJD(1)) * t607;
t524 = -t610 * t537 + t649 * t538;
t557 = -t581 * qJD(3) + t649 * t593 + t610 * t603;
t582 = t610 * t604 + t649 * t629;
t561 = t581 * mrSges(4,1) + t582 * mrSges(4,2);
t566 = t598 * mrSges(4,2) - t581 * mrSges(4,3);
t568 = t581 * mrSges(5,1) + t598 * mrSges(5,3);
t586 = qJDD(3) + t594;
t560 = t581 * pkin(3) - t582 * qJ(4);
t597 = t598 ^ 2;
t522 = -t586 * pkin(3) - t597 * qJ(4) + t582 * t560 + qJDD(4) - t524;
t517 = (t581 * t582 - t586) * pkin(9) + (t557 - t643) * pkin(4) + t522;
t556 = t582 * qJD(3) + t610 * t593 - t649 * t603;
t570 = t582 * pkin(4) + t598 * pkin(9);
t580 = t581 ^ 2;
t558 = -g(3) * t641 + t589 * t639 - t611 * t590;
t536 = -t603 * pkin(2) - t602 * pkin(8) + t592 * t629 - t558;
t617 = (-t557 - t643) * qJ(4) + t536 + (-t598 * pkin(3) + t650) * t582;
t520 = -t580 * pkin(4) - t582 * t570 + (pkin(3) + pkin(9)) * t556 + t617;
t609 = sin(qJ(5));
t613 = cos(qJ(5));
t515 = t613 * t517 - t609 * t520;
t564 = t613 * t581 + t609 * t598;
t528 = t564 * qJD(5) + t609 * t556 + t613 * t586;
t565 = t609 * t581 - t613 * t598;
t539 = -t564 * mrSges(6,1) + t565 * mrSges(6,2);
t579 = qJD(5) + t582;
t542 = -t579 * mrSges(6,2) + t564 * mrSges(6,3);
t553 = qJDD(5) + t557;
t513 = m(6) * t515 + t553 * mrSges(6,1) - t528 * mrSges(6,3) - t565 * t539 + t579 * t542;
t516 = t609 * t517 + t613 * t520;
t527 = -t565 * qJD(5) + t613 * t556 - t609 * t586;
t543 = t579 * mrSges(6,1) - t565 * mrSges(6,3);
t514 = m(6) * t516 - t553 * mrSges(6,2) + t527 * mrSges(6,3) + t564 * t539 - t579 * t543;
t505 = t613 * t513 + t609 * t514;
t562 = -t581 * mrSges(5,2) - t582 * mrSges(5,3);
t621 = -m(5) * t522 - t557 * mrSges(5,1) - t582 * t562 - t505;
t503 = m(4) * t524 - t557 * mrSges(4,3) - t582 * t561 + (-t566 + t568) * t598 + (mrSges(4,1) - mrSges(5,2)) * t586 + t621;
t525 = t649 * t537 + t610 * t538;
t567 = -t598 * mrSges(4,1) - t582 * mrSges(4,3);
t620 = -t597 * pkin(3) + t586 * qJ(4) - t581 * t560 + t525;
t521 = 0.2e1 * qJD(4) * t598 - t620;
t569 = t582 * mrSges(5,1) - t598 * mrSges(5,2);
t519 = -t556 * pkin(4) - t580 * pkin(9) + (t650 - t570) * t598 + t620;
t622 = -m(6) * t519 + t527 * mrSges(6,1) - t528 * mrSges(6,2) + t564 * t542 - t565 * t543;
t619 = -m(5) * t521 + t586 * mrSges(5,3) - t598 * t569 - t622;
t510 = m(4) * t525 - t586 * mrSges(4,2) + t598 * t567 + (-t561 - t562) * t581 + (-mrSges(4,3) - mrSges(5,1)) * t556 + t619;
t626 = -t610 * t503 + t649 * t510;
t495 = m(3) * t559 - t603 * mrSges(3,2) - t594 * mrSges(3,3) - t604 * t587 + t591 * t628 + t626;
t498 = t649 * t503 + t610 * t510;
t575 = -t607 * t589 - t647;
t588 = -t604 * mrSges(3,2) + mrSges(3,3) * t628;
t497 = m(3) * t575 + t594 * mrSges(3,1) + t593 * mrSges(3,2) + (t587 * t611 - t588 * t614) * t632 + t498;
t523 = t556 * pkin(3) + t617;
t637 = -t609 * t513 + t613 * t514;
t625 = -m(5) * t523 + t556 * mrSges(5,2) + t581 * t568 - t637;
t618 = -m(4) * t536 - t556 * mrSges(4,1) - t581 * t566 + (-t567 + t569) * t582 + (-mrSges(4,2) + mrSges(5,3)) * t557 + t625;
t501 = m(3) * t558 + t603 * mrSges(3,1) - t593 * mrSges(3,3) + t604 * t588 - t591 * t629 + t618;
t485 = t495 * t640 - t607 * t497 + t501 * t639;
t483 = m(2) * t600 + qJDD(1) * mrSges(2,1) - t616 * mrSges(2,2) + t485;
t491 = t614 * t495 - t611 * t501;
t490 = m(2) * t601 - t616 * mrSges(2,1) - qJDD(1) * mrSges(2,2) + t491;
t638 = t615 * t483 + t612 * t490;
t636 = t644 * t581 - t645 * t582 + t651 * t598;
t635 = t652 * t581 + t646 * t582 - t644 * t598;
t634 = t646 * t581 - t653 * t582 + t645 * t598;
t484 = t495 * t642 + t608 * t497 + t501 * t641;
t627 = -t612 * t483 + t615 * t490;
t504 = -t557 * mrSges(5,3) - t582 * t569 - t625;
t529 = Ifges(6,5) * t565 + Ifges(6,6) * t564 + Ifges(6,3) * t579;
t531 = Ifges(6,1) * t565 + Ifges(6,4) * t564 + Ifges(6,5) * t579;
t506 = -mrSges(6,1) * t519 + mrSges(6,3) * t516 + Ifges(6,4) * t528 + Ifges(6,2) * t527 + Ifges(6,6) * t553 - t565 * t529 + t579 * t531;
t530 = Ifges(6,4) * t565 + Ifges(6,2) * t564 + Ifges(6,6) * t579;
t507 = mrSges(6,2) * t519 - mrSges(6,3) * t515 + Ifges(6,1) * t528 + Ifges(6,4) * t527 + Ifges(6,5) * t553 + t564 * t529 - t579 * t530;
t486 = -mrSges(4,1) * t536 - mrSges(5,1) * t521 + mrSges(5,2) * t523 + mrSges(4,3) * t525 - pkin(3) * t504 - pkin(4) * t622 - pkin(9) * t637 - t613 * t506 - t609 * t507 + t652 * t556 + t646 * t557 + t636 * t582 + t644 * t586 + t634 * t598;
t487 = mrSges(5,1) * t522 + mrSges(6,1) * t515 + mrSges(4,2) * t536 - mrSges(6,2) * t516 - mrSges(4,3) * t524 - mrSges(5,3) * t523 + Ifges(6,5) * t528 + Ifges(6,6) * t527 + Ifges(6,3) * t553 + pkin(4) * t505 - qJ(4) * t504 + t565 * t530 - t564 * t531 + t635 * t598 + t645 * t586 + t636 * t581 + t653 * t557 - t646 * t556;
t572 = Ifges(3,3) * t604 + (Ifges(3,5) * t611 + Ifges(3,6) * t614) * t632;
t573 = Ifges(3,6) * t604 + (Ifges(3,4) * t611 + Ifges(3,2) * t614) * t632;
t480 = mrSges(3,2) * t575 - mrSges(3,3) * t558 + Ifges(3,1) * t593 - Ifges(3,4) * t594 + Ifges(3,5) * t603 - pkin(8) * t498 - t610 * t486 + t649 * t487 + t572 * t628 - t604 * t573;
t574 = Ifges(3,5) * t604 + (Ifges(3,1) * t611 + Ifges(3,4) * t614) * t632;
t481 = Ifges(3,4) * t593 - Ifges(3,2) * t594 + pkin(9) * t505 + t609 * t506 - pkin(3) * (t598 * t568 + t621) - t613 * t507 + mrSges(3,3) * t559 + Ifges(3,6) * t603 + t604 * t574 - mrSges(4,1) * t524 + mrSges(4,2) * t525 - mrSges(3,1) * t575 - qJ(4) * t619 + mrSges(5,3) * t521 - mrSges(5,2) * t522 - t572 * t629 - pkin(2) * t498 + (pkin(3) * mrSges(5,2) - t651) * t586 - t635 * t582 + (qJ(4) * t562 + t634) * t581 - t645 * t557 + (qJ(4) * mrSges(5,1) + t644) * t556;
t623 = pkin(7) * t491 + t480 * t611 + t481 * t614;
t479 = Ifges(3,5) * t593 - Ifges(3,6) * t594 + Ifges(3,3) * t603 + mrSges(3,1) * t558 - mrSges(3,2) * t559 + t610 * t487 + t649 * t486 + pkin(2) * t618 + pkin(8) * t626 + (t573 * t611 - t574 * t614) * t632;
t478 = -mrSges(2,2) * g(3) - mrSges(2,3) * t600 + Ifges(2,5) * qJDD(1) - t616 * Ifges(2,6) + t614 * t480 - t611 * t481 + (-t484 * t607 - t485 * t608) * pkin(7);
t477 = mrSges(2,1) * g(3) + mrSges(2,3) * t601 + t616 * Ifges(2,5) + Ifges(2,6) * qJDD(1) - pkin(1) * t484 - t607 * t479 + t623 * t608;
t1 = [-m(1) * g(1) + t627; -m(1) * g(2) + t638; (-m(1) - m(2)) * g(3) + t484; -mrSges(1,2) * g(3) + mrSges(1,3) * g(2) - pkin(6) * t638 - t612 * t477 + t615 * t478; mrSges(1,1) * g(3) - mrSges(1,3) * g(1) + pkin(6) * t627 + t615 * t477 + t612 * t478; -mrSges(1,1) * g(2) + mrSges(2,1) * t600 + mrSges(1,2) * g(1) - mrSges(2,2) * t601 + Ifges(2,3) * qJDD(1) + pkin(1) * t485 + t608 * t479 + t607 * t623;];
tauB = t1;
