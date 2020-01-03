% Calculate vector of inverse dynamics base forces with Newton-Euler for
% S5RRRRP10
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
%   pkin=[a2,a3,a4,a5,alpha2,d1,d2,d3,d4]';
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
% Datum: 2019-12-31 22:13
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauB = S5RRRRP10_invdynB_fixb_snew_vp2(qJ, qJD, qJDD, g, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(5,1),zeros(3,1),zeros(9,1),zeros(6,1),zeros(6,3),zeros(6,6)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RRRRP10_invdynB_fixb_snew_vp2: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RRRRP10_invdynB_fixb_snew_vp2: qJD has to be [5x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [5 1]), ...
  'S5RRRRP10_invdynB_fixb_snew_vp2: qJDD has to be [5x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S5RRRRP10_invdynB_fixb_snew_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S5RRRRP10_invdynB_fixb_snew_vp2: pkin has to be [9x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5RRRRP10_invdynB_fixb_snew_vp2: m has to be [6x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [6,3]), ...
  'S5RRRRP10_invdynB_fixb_snew_vp2: mrSges has to be [6x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [6 6]), ...
  'S5RRRRP10_invdynB_fixb_snew_vp2: Ifges has to be [6x6] (double)'); 

%% Symbolic Calculation
% From invdyn_fixb_NewtonEuler_linkframe_tauB_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 22:08:46
% EndTime: 2019-12-31 22:08:57
% DurationCPUTime: 6.81s
% Computational Cost: add. (94819->303), mult. (203158->384), div. (0->0), fcn. (154656->10), ass. (0->128)
t630 = Ifges(5,1) + Ifges(6,1);
t623 = Ifges(5,4) + Ifges(6,4);
t622 = Ifges(5,5) + Ifges(6,5);
t629 = Ifges(5,2) + Ifges(6,2);
t628 = Ifges(5,6) + Ifges(6,6);
t627 = Ifges(5,3) + Ifges(6,3);
t585 = sin(pkin(5));
t589 = sin(qJ(2));
t593 = cos(qJ(2));
t608 = qJD(1) * qJD(2);
t575 = (-qJDD(1) * t593 + t589 * t608) * t585;
t626 = pkin(7) * t585;
t586 = cos(pkin(5));
t625 = t586 * g(3);
t624 = -mrSges(5,2) - mrSges(6,2);
t620 = t585 * t589;
t619 = t585 * t593;
t618 = t586 * t589;
t617 = t586 * t593;
t590 = sin(qJ(1));
t594 = cos(qJ(1));
t578 = t590 * g(1) - t594 * g(2);
t595 = qJD(1) ^ 2;
t570 = qJDD(1) * pkin(1) + t595 * t626 + t578;
t579 = -t594 * g(1) - t590 * g(2);
t571 = -t595 * pkin(1) + qJDD(1) * t626 + t579;
t611 = t570 * t618 + t593 * t571;
t548 = -g(3) * t620 + t611;
t582 = t586 * qJD(1) + qJD(2);
t610 = qJD(1) * t585;
t605 = t589 * t610;
t568 = t582 * mrSges(3,1) - mrSges(3,3) * t605;
t572 = (-mrSges(3,1) * t593 + mrSges(3,2) * t589) * t610;
t581 = t586 * qJDD(1) + qJDD(2);
t573 = (-pkin(2) * t593 - pkin(8) * t589) * t610;
t580 = t582 ^ 2;
t609 = qJD(1) * t593;
t527 = -t580 * pkin(2) + t581 * pkin(8) + (-g(3) * t589 + t573 * t609) * t585 + t611;
t574 = (qJDD(1) * t589 + t593 * t608) * t585;
t528 = t575 * pkin(2) - t574 * pkin(8) - t625 + (-t570 + (pkin(2) * t589 - pkin(8) * t593) * t582 * qJD(1)) * t585;
t588 = sin(qJ(3));
t592 = cos(qJ(3));
t508 = t592 * t527 + t588 * t528;
t564 = t588 * t582 + t592 * t605;
t545 = -t564 * qJD(3) - t588 * t574 + t592 * t581;
t563 = t592 * t582 - t588 * t605;
t549 = -t563 * mrSges(4,1) + t564 * mrSges(4,2);
t604 = t585 * t609;
t577 = qJD(3) - t604;
t555 = t577 * mrSges(4,1) - t564 * mrSges(4,3);
t567 = qJDD(3) + t575;
t550 = -t563 * pkin(3) - t564 * pkin(9);
t576 = t577 ^ 2;
t503 = -t576 * pkin(3) + t567 * pkin(9) + t563 * t550 + t508;
t547 = -g(3) * t619 + t570 * t617 - t589 * t571;
t526 = -t581 * pkin(2) - t580 * pkin(8) + t573 * t605 - t547;
t546 = t563 * qJD(3) + t592 * t574 + t588 * t581;
t506 = (-t563 * t577 - t546) * pkin(9) + (t564 * t577 - t545) * pkin(3) + t526;
t587 = sin(qJ(4));
t591 = cos(qJ(4));
t498 = -t587 * t503 + t591 * t506;
t552 = -t587 * t564 + t591 * t577;
t513 = t552 * qJD(4) + t591 * t546 + t587 * t567;
t553 = t591 * t564 + t587 * t577;
t530 = -t552 * mrSges(6,1) + t553 * mrSges(6,2);
t531 = -t552 * mrSges(5,1) + t553 * mrSges(5,2);
t562 = qJD(4) - t563;
t534 = -t562 * mrSges(5,2) + t552 * mrSges(5,3);
t543 = qJDD(4) - t545;
t495 = -0.2e1 * qJD(5) * t553 + (t552 * t562 - t513) * qJ(5) + (t552 * t553 + t543) * pkin(4) + t498;
t533 = -t562 * mrSges(6,2) + t552 * mrSges(6,3);
t607 = m(6) * t495 + t543 * mrSges(6,1) + t562 * t533;
t488 = m(5) * t498 + t543 * mrSges(5,1) + t562 * t534 + (-t530 - t531) * t553 + (-mrSges(5,3) - mrSges(6,3)) * t513 + t607;
t499 = t591 * t503 + t587 * t506;
t512 = -t553 * qJD(4) - t587 * t546 + t591 * t567;
t535 = t562 * pkin(4) - t553 * qJ(5);
t551 = t552 ^ 2;
t497 = -t551 * pkin(4) + t512 * qJ(5) + 0.2e1 * qJD(5) * t552 - t562 * t535 + t499;
t606 = m(6) * t497 + t512 * mrSges(6,3) + t552 * t530;
t536 = t562 * mrSges(6,1) - t553 * mrSges(6,3);
t612 = -t562 * mrSges(5,1) + t553 * mrSges(5,3) - t536;
t490 = m(5) * t499 + t512 * mrSges(5,3) + t552 * t531 + t624 * t543 + t612 * t562 + t606;
t601 = -t587 * t488 + t591 * t490;
t485 = m(4) * t508 - t567 * mrSges(4,2) + t545 * mrSges(4,3) + t563 * t549 - t577 * t555 + t601;
t507 = -t588 * t527 + t592 * t528;
t554 = -t577 * mrSges(4,2) + t563 * mrSges(4,3);
t502 = -t567 * pkin(3) - t576 * pkin(9) + t564 * t550 - t507;
t500 = -t512 * pkin(4) - t551 * qJ(5) + t553 * t535 + qJDD(5) + t502;
t600 = m(6) * t500 - t512 * mrSges(6,1) - t552 * t533;
t596 = -m(5) * t502 + t512 * mrSges(5,1) + t624 * t513 + t552 * t534 + t612 * t553 - t600;
t492 = m(4) * t507 + t567 * mrSges(4,1) - t546 * mrSges(4,3) - t564 * t549 + t577 * t554 + t596;
t602 = t592 * t485 - t588 * t492;
t475 = m(3) * t548 - t581 * mrSges(3,2) - t575 * mrSges(3,3) - t582 * t568 + t572 * t604 + t602;
t478 = t588 * t485 + t592 * t492;
t559 = -t585 * t570 - t625;
t569 = -t582 * mrSges(3,2) + mrSges(3,3) * t604;
t477 = m(3) * t559 + t575 * mrSges(3,1) + t574 * mrSges(3,2) + (t568 * t589 - t569 * t593) * t610 + t478;
t487 = t591 * t488 + t587 * t490;
t597 = -m(4) * t526 + t545 * mrSges(4,1) - t546 * mrSges(4,2) + t563 * t554 - t564 * t555 - t487;
t482 = m(3) * t547 + t581 * mrSges(3,1) - t574 * mrSges(3,3) + t582 * t569 - t572 * t605 + t597;
t465 = t475 * t618 - t585 * t477 + t482 * t617;
t463 = m(2) * t578 + qJDD(1) * mrSges(2,1) - t595 * mrSges(2,2) + t465;
t470 = t593 * t475 - t589 * t482;
t469 = m(2) * t579 - t595 * mrSges(2,1) - qJDD(1) * mrSges(2,2) + t470;
t616 = t594 * t463 + t590 * t469;
t615 = t628 * t552 + t622 * t553 + t627 * t562;
t614 = -t629 * t552 - t623 * t553 - t628 * t562;
t613 = t623 * t552 + t630 * t553 + t622 * t562;
t464 = t475 * t620 + t586 * t477 + t482 * t619;
t603 = -t590 * t463 + t594 * t469;
t479 = -mrSges(5,1) * t502 + mrSges(5,3) * t499 - mrSges(6,1) * t500 + mrSges(6,3) * t497 - pkin(4) * t600 + qJ(5) * t606 + (-qJ(5) * t536 + t613) * t562 + (-pkin(4) * t536 - t615) * t553 + (-qJ(5) * mrSges(6,2) + t628) * t543 + (-pkin(4) * mrSges(6,2) + t623) * t513 + t629 * t512;
t493 = -t513 * mrSges(6,3) - t553 * t530 + t607;
t486 = mrSges(5,2) * t502 + mrSges(6,2) * t500 - mrSges(5,3) * t498 - mrSges(6,3) * t495 - qJ(5) * t493 + t623 * t512 + t630 * t513 + t622 * t543 + t615 * t552 + t614 * t562;
t539 = Ifges(4,5) * t564 + Ifges(4,6) * t563 + Ifges(4,3) * t577;
t540 = Ifges(4,4) * t564 + Ifges(4,2) * t563 + Ifges(4,6) * t577;
t466 = mrSges(4,2) * t526 - mrSges(4,3) * t507 + Ifges(4,1) * t546 + Ifges(4,4) * t545 + Ifges(4,5) * t567 - pkin(9) * t487 - t587 * t479 + t591 * t486 + t563 * t539 - t577 * t540;
t541 = Ifges(4,1) * t564 + Ifges(4,4) * t563 + Ifges(4,5) * t577;
t471 = -mrSges(4,1) * t526 - mrSges(5,1) * t498 - mrSges(6,1) * t495 + mrSges(5,2) * t499 + mrSges(6,2) * t497 + mrSges(4,3) * t508 + Ifges(4,4) * t546 + Ifges(4,2) * t545 + Ifges(4,6) * t567 - pkin(3) * t487 - pkin(4) * t493 - t564 * t539 + t577 * t541 + t614 * t553 + t613 * t552 - t627 * t543 - t622 * t513 - t628 * t512;
t556 = Ifges(3,3) * t582 + (Ifges(3,5) * t589 + Ifges(3,6) * t593) * t610;
t557 = Ifges(3,6) * t582 + (Ifges(3,4) * t589 + Ifges(3,2) * t593) * t610;
t460 = mrSges(3,2) * t559 - mrSges(3,3) * t547 + Ifges(3,1) * t574 - Ifges(3,4) * t575 + Ifges(3,5) * t581 - pkin(8) * t478 + t592 * t466 - t588 * t471 + t556 * t604 - t582 * t557;
t558 = Ifges(3,5) * t582 + (Ifges(3,1) * t589 + Ifges(3,4) * t593) * t610;
t461 = Ifges(3,4) * t574 - Ifges(3,2) * t575 + Ifges(3,6) * t581 - t556 * t605 + t582 * t558 - mrSges(3,1) * t559 + mrSges(3,3) * t548 - Ifges(4,5) * t546 - Ifges(4,6) * t545 - Ifges(4,3) * t567 - t564 * t540 + t563 * t541 - mrSges(4,1) * t507 + mrSges(4,2) * t508 - t587 * t486 - t591 * t479 - pkin(3) * t596 - pkin(9) * t601 - pkin(2) * t478;
t598 = pkin(7) * t470 + t460 * t589 + t461 * t593;
t459 = Ifges(3,5) * t574 - Ifges(3,6) * t575 + Ifges(3,3) * t581 + mrSges(3,1) * t547 - mrSges(3,2) * t548 + t588 * t466 + t592 * t471 + pkin(2) * t597 + pkin(8) * t602 + (t557 * t589 - t558 * t593) * t610;
t458 = -mrSges(2,2) * g(3) - mrSges(2,3) * t578 + Ifges(2,5) * qJDD(1) - t595 * Ifges(2,6) + t593 * t460 - t589 * t461 + (-t464 * t585 - t465 * t586) * pkin(7);
t457 = mrSges(2,1) * g(3) + mrSges(2,3) * t579 + t595 * Ifges(2,5) + Ifges(2,6) * qJDD(1) - pkin(1) * t464 - t585 * t459 + t598 * t586;
t1 = [-m(1) * g(1) + t603; -m(1) * g(2) + t616; (-m(1) - m(2)) * g(3) + t464; -mrSges(1,2) * g(3) + mrSges(1,3) * g(2) - pkin(6) * t616 - t590 * t457 + t594 * t458; mrSges(1,1) * g(3) - mrSges(1,3) * g(1) + pkin(6) * t603 + t594 * t457 + t590 * t458; -mrSges(1,1) * g(2) + mrSges(2,1) * t578 + mrSges(1,2) * g(1) - mrSges(2,2) * t579 + Ifges(2,3) * qJDD(1) + pkin(1) * t465 + t586 * t459 + t598 * t585;];
tauB = t1;
