% Calculate vector of inverse dynamics base forces with Newton-Euler for
% S5RRPRR10
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
% pkin [10x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,alpha2,d1,d2,d4,d5,theta3]';
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
% Datum: 2019-12-31 20:28
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauB = S5RRPRR10_invdynB_fixb_snew_vp2(qJ, qJD, qJDD, g, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(5,1),zeros(3,1),zeros(10,1),zeros(6,1),zeros(6,3),zeros(6,6)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RRPRR10_invdynB_fixb_snew_vp2: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RRPRR10_invdynB_fixb_snew_vp2: qJD has to be [5x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [5 1]), ...
  'S5RRPRR10_invdynB_fixb_snew_vp2: qJDD has to be [5x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S5RRPRR10_invdynB_fixb_snew_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S5RRPRR10_invdynB_fixb_snew_vp2: pkin has to be [10x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5RRPRR10_invdynB_fixb_snew_vp2: m has to be [6x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [6,3]), ...
  'S5RRPRR10_invdynB_fixb_snew_vp2: mrSges has to be [6x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [6 6]), ...
  'S5RRPRR10_invdynB_fixb_snew_vp2: Ifges has to be [6x6] (double)'); 

%% Symbolic Calculation
% From invdyn_fixb_NewtonEuler_linkframe_tauB_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 20:24:38
% EndTime: 2019-12-31 20:24:50
% DurationCPUTime: 11.47s
% Computational Cost: add. (161058->323), mult. (422320->426), div. (0->0), fcn. (323902->12), ass. (0->135)
t638 = -2 * qJD(3);
t601 = sin(pkin(10));
t603 = cos(pkin(10));
t607 = sin(qJ(2));
t611 = cos(qJ(2));
t602 = sin(pkin(5));
t630 = qJD(1) * t602;
t582 = (t601 * t607 - t603 * t611) * t630;
t628 = qJD(1) * qJD(2);
t591 = (qJDD(1) * t607 + t611 * t628) * t602;
t604 = cos(pkin(5));
t596 = qJDD(1) * t604 + qJDD(2);
t597 = qJD(1) * t604 + qJD(2);
t608 = sin(qJ(1));
t612 = cos(qJ(1));
t593 = t608 * g(1) - g(2) * t612;
t613 = qJD(1) ^ 2;
t637 = pkin(7) * t602;
t588 = qJDD(1) * pkin(1) + t613 * t637 + t593;
t594 = -g(1) * t612 - g(2) * t608;
t589 = -pkin(1) * t613 + qJDD(1) * t637 + t594;
t632 = t604 * t611;
t619 = t588 * t632 - t607 * t589;
t636 = t602 ^ 2 * t613;
t533 = t596 * pkin(2) - t591 * qJ(3) + (pkin(2) * t607 * t636 + (qJ(3) * qJD(1) * t597 - g(3)) * t602) * t611 + t619;
t633 = t604 * t607;
t635 = t602 * t607;
t559 = -g(3) * t635 + t588 * t633 + t611 * t589;
t626 = t607 * t630;
t585 = pkin(2) * t597 - qJ(3) * t626;
t592 = (qJDD(1) * t611 - t607 * t628) * t602;
t627 = t611 ^ 2 * t636;
t536 = -pkin(2) * t627 + qJ(3) * t592 - t585 * t597 + t559;
t583 = (t601 * t611 + t603 * t607) * t630;
t518 = t533 * t603 - t601 * t536 + t583 * t638;
t634 = t602 * t611;
t519 = t601 * t533 + t603 * t536 + t582 * t638;
t560 = mrSges(4,1) * t582 + mrSges(4,2) * t583;
t564 = -t591 * t601 + t592 * t603;
t570 = mrSges(4,1) * t597 - mrSges(4,3) * t583;
t561 = pkin(3) * t582 - pkin(8) * t583;
t595 = t597 ^ 2;
t517 = -pkin(3) * t595 + pkin(8) * t596 - t561 * t582 + t519;
t574 = -g(3) * t604 - t588 * t602;
t546 = -pkin(2) * t592 - qJ(3) * t627 + t585 * t626 + qJDD(3) + t574;
t565 = t591 * t603 + t592 * t601;
t521 = (t582 * t597 - t565) * pkin(8) + (t583 * t597 - t564) * pkin(3) + t546;
t606 = sin(qJ(4));
t610 = cos(qJ(4));
t514 = t610 * t517 + t606 * t521;
t568 = t583 * t610 + t597 * t606;
t543 = -qJD(4) * t568 - t565 * t606 + t596 * t610;
t567 = -t583 * t606 + t597 * t610;
t547 = -mrSges(5,1) * t567 + mrSges(5,2) * t568;
t581 = qJD(4) + t582;
t553 = mrSges(5,1) * t581 - mrSges(5,3) * t568;
t563 = qJDD(4) - t564;
t548 = -pkin(4) * t567 - pkin(9) * t568;
t580 = t581 ^ 2;
t511 = -pkin(4) * t580 + pkin(9) * t563 + t548 * t567 + t514;
t516 = -pkin(3) * t596 - pkin(8) * t595 + t583 * t561 - t518;
t544 = qJD(4) * t567 + t565 * t610 + t596 * t606;
t512 = (-t567 * t581 - t544) * pkin(9) + (t568 * t581 - t543) * pkin(4) + t516;
t605 = sin(qJ(5));
t609 = cos(qJ(5));
t508 = -t511 * t605 + t512 * t609;
t550 = -t568 * t605 + t581 * t609;
t524 = qJD(5) * t550 + t544 * t609 + t563 * t605;
t551 = t568 * t609 + t581 * t605;
t529 = -mrSges(6,1) * t550 + mrSges(6,2) * t551;
t566 = qJD(5) - t567;
t534 = -mrSges(6,2) * t566 + mrSges(6,3) * t550;
t542 = qJDD(5) - t543;
t506 = m(6) * t508 + mrSges(6,1) * t542 - mrSges(6,3) * t524 - t529 * t551 + t534 * t566;
t509 = t511 * t609 + t512 * t605;
t523 = -qJD(5) * t551 - t544 * t605 + t563 * t609;
t535 = mrSges(6,1) * t566 - mrSges(6,3) * t551;
t507 = m(6) * t509 - mrSges(6,2) * t542 + mrSges(6,3) * t523 + t529 * t550 - t535 * t566;
t621 = -t506 * t605 + t609 * t507;
t499 = m(5) * t514 - mrSges(5,2) * t563 + mrSges(5,3) * t543 + t547 * t567 - t553 * t581 + t621;
t513 = -t517 * t606 + t521 * t610;
t552 = -mrSges(5,2) * t581 + mrSges(5,3) * t567;
t510 = -pkin(4) * t563 - pkin(9) * t580 + t548 * t568 - t513;
t616 = -m(6) * t510 + t523 * mrSges(6,1) - mrSges(6,2) * t524 + t550 * t534 - t535 * t551;
t504 = m(5) * t513 + mrSges(5,1) * t563 - mrSges(5,3) * t544 - t547 * t568 + t552 * t581 + t616;
t622 = t610 * t499 - t504 * t606;
t491 = m(4) * t519 - mrSges(4,2) * t596 + mrSges(4,3) * t564 - t560 * t582 - t570 * t597 + t622;
t569 = -mrSges(4,2) * t597 - mrSges(4,3) * t582;
t500 = t506 * t609 + t507 * t605;
t614 = -m(5) * t516 + t543 * mrSges(5,1) - mrSges(5,2) * t544 + t567 * t552 - t553 * t568 - t500;
t496 = m(4) * t518 + mrSges(4,1) * t596 - mrSges(4,3) * t565 - t560 * t583 + t569 * t597 + t614;
t486 = t601 * t491 + t603 * t496;
t558 = -g(3) * t634 + t619;
t625 = t611 * t630;
t587 = -mrSges(3,2) * t597 + mrSges(3,3) * t625;
t590 = (-mrSges(3,1) * t611 + mrSges(3,2) * t607) * t630;
t484 = m(3) * t558 + mrSges(3,1) * t596 - mrSges(3,3) * t591 + t587 * t597 - t590 * t626 + t486;
t586 = mrSges(3,1) * t597 - mrSges(3,3) * t626;
t623 = t603 * t491 - t601 * t496;
t485 = m(3) * t559 - mrSges(3,2) * t596 + mrSges(3,3) * t592 - t586 * t597 + t590 * t625 + t623;
t494 = t606 * t499 + t610 * t504;
t615 = m(4) * t546 - mrSges(4,1) * t564 + t565 * mrSges(4,2) + t569 * t582 + t583 * t570 + t494;
t493 = m(3) * t574 - mrSges(3,1) * t592 + mrSges(3,2) * t591 + (t586 * t607 - t587 * t611) * t630 + t615;
t472 = t484 * t632 + t485 * t633 - t493 * t602;
t470 = m(2) * t593 + qJDD(1) * mrSges(2,1) - mrSges(2,2) * t613 + t472;
t477 = -t484 * t607 + t611 * t485;
t476 = m(2) * t594 - mrSges(2,1) * t613 - qJDD(1) * mrSges(2,2) + t477;
t631 = t612 * t470 + t608 * t476;
t471 = t484 * t634 + t485 * t635 + t604 * t493;
t624 = -t470 * t608 + t612 * t476;
t525 = Ifges(6,5) * t551 + Ifges(6,6) * t550 + Ifges(6,3) * t566;
t527 = Ifges(6,1) * t551 + Ifges(6,4) * t550 + Ifges(6,5) * t566;
t501 = -mrSges(6,1) * t510 + mrSges(6,3) * t509 + Ifges(6,4) * t524 + Ifges(6,2) * t523 + Ifges(6,6) * t542 - t525 * t551 + t527 * t566;
t526 = Ifges(6,4) * t551 + Ifges(6,2) * t550 + Ifges(6,6) * t566;
t502 = mrSges(6,2) * t510 - mrSges(6,3) * t508 + Ifges(6,1) * t524 + Ifges(6,4) * t523 + Ifges(6,5) * t542 + t525 * t550 - t526 * t566;
t537 = Ifges(5,5) * t568 + Ifges(5,6) * t567 + Ifges(5,3) * t581;
t538 = Ifges(5,4) * t568 + Ifges(5,2) * t567 + Ifges(5,6) * t581;
t487 = mrSges(5,2) * t516 - mrSges(5,3) * t513 + Ifges(5,1) * t544 + Ifges(5,4) * t543 + Ifges(5,5) * t563 - pkin(9) * t500 - t501 * t605 + t502 * t609 + t537 * t567 - t538 * t581;
t539 = Ifges(5,1) * t568 + Ifges(5,4) * t567 + Ifges(5,5) * t581;
t488 = -mrSges(5,1) * t516 - mrSges(6,1) * t508 + mrSges(6,2) * t509 + mrSges(5,3) * t514 + Ifges(5,4) * t544 - Ifges(6,5) * t524 + Ifges(5,2) * t543 + Ifges(5,6) * t563 - Ifges(6,6) * t523 - Ifges(6,3) * t542 - pkin(4) * t500 - t526 * t551 + t527 * t550 - t537 * t568 + t539 * t581;
t554 = Ifges(4,5) * t583 - Ifges(4,6) * t582 + Ifges(4,3) * t597;
t555 = Ifges(4,4) * t583 - Ifges(4,2) * t582 + Ifges(4,6) * t597;
t473 = mrSges(4,2) * t546 - mrSges(4,3) * t518 + Ifges(4,1) * t565 + Ifges(4,4) * t564 + Ifges(4,5) * t596 - pkin(8) * t494 + t487 * t610 - t488 * t606 - t554 * t582 - t555 * t597;
t556 = Ifges(4,1) * t583 - Ifges(4,4) * t582 + Ifges(4,5) * t597;
t478 = Ifges(4,4) * t565 + Ifges(4,2) * t564 + Ifges(4,6) * t596 - t583 * t554 + t597 * t556 - mrSges(4,1) * t546 + mrSges(4,3) * t519 - Ifges(5,5) * t544 - Ifges(5,6) * t543 - Ifges(5,3) * t563 - t568 * t538 + t567 * t539 - mrSges(5,1) * t513 + mrSges(5,2) * t514 - t605 * t502 - t609 * t501 - pkin(4) * t616 - pkin(9) * t621 - pkin(3) * t494;
t571 = Ifges(3,3) * t597 + (Ifges(3,5) * t607 + Ifges(3,6) * t611) * t630;
t573 = Ifges(3,5) * t597 + (Ifges(3,1) * t607 + Ifges(3,4) * t611) * t630;
t466 = -mrSges(3,1) * t574 + mrSges(3,3) * t559 + Ifges(3,4) * t591 + Ifges(3,2) * t592 + Ifges(3,6) * t596 - pkin(2) * t615 + qJ(3) * t623 + t601 * t473 + t603 * t478 - t571 * t626 + t597 * t573;
t572 = Ifges(3,6) * t597 + (Ifges(3,4) * t607 + Ifges(3,2) * t611) * t630;
t467 = mrSges(3,2) * t574 - mrSges(3,3) * t558 + Ifges(3,1) * t591 + Ifges(3,4) * t592 + Ifges(3,5) * t596 - qJ(3) * t486 + t473 * t603 - t478 * t601 + t571 * t625 - t572 * t597;
t617 = pkin(7) * t477 + t466 * t611 + t467 * t607;
t468 = Ifges(3,5) * t591 + Ifges(3,6) * t592 + mrSges(3,1) * t558 - mrSges(3,2) * t559 + Ifges(4,5) * t565 + Ifges(4,6) * t564 + t583 * t555 + t582 * t556 + mrSges(4,1) * t518 - mrSges(4,2) * t519 + t606 * t487 + t610 * t488 + pkin(3) * t614 + pkin(8) * t622 + pkin(2) * t486 + (Ifges(3,3) + Ifges(4,3)) * t596 + (t572 * t607 - t573 * t611) * t630;
t465 = -mrSges(2,2) * g(3) - mrSges(2,3) * t593 + Ifges(2,5) * qJDD(1) - t613 * Ifges(2,6) - t607 * t466 + t611 * t467 + (-t471 * t602 - t472 * t604) * pkin(7);
t464 = mrSges(2,1) * g(3) + mrSges(2,3) * t594 + t613 * Ifges(2,5) + Ifges(2,6) * qJDD(1) - pkin(1) * t471 - t602 * t468 + t604 * t617;
t1 = [-m(1) * g(1) + t624; -m(1) * g(2) + t631; (-m(1) - m(2)) * g(3) + t471; -mrSges(1,2) * g(3) + mrSges(1,3) * g(2) - pkin(6) * t631 - t608 * t464 + t612 * t465; mrSges(1,1) * g(3) - mrSges(1,3) * g(1) + pkin(6) * t624 + t612 * t464 + t608 * t465; -mrSges(1,1) * g(2) + mrSges(2,1) * t593 + mrSges(1,2) * g(1) - mrSges(2,2) * t594 + Ifges(2,3) * qJDD(1) + pkin(1) * t472 + t468 * t604 + t602 * t617;];
tauB = t1;
