% Calculate vector of inverse dynamics base forces with Newton-Euler for
% S5RRRRR14
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
%   pkin=[a2,a3,a4,a5,alpha3,d1,d2,d3,d4,d5]';
% m [6x1]
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
% Datum: 2024-09-27 18:44
% Revision: 4128c0dfad36ceb08eaa7a9fbc8f797c888f59be (2024-07-29)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauB = S5RRRRR14_invdynB_fixb_snew_vp2(qJ, qJD, qJDD, g, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(5,1),zeros(3,1),zeros(10,1),zeros(6,1),zeros(6,3),zeros(6,6)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RRRRR14_invdynB_fixb_snew_vp2: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RRRRR14_invdynB_fixb_snew_vp2: qJD has to be [5x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [5 1]), ...
  'S5RRRRR14_invdynB_fixb_snew_vp2: qJDD has to be [5x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S5RRRRR14_invdynB_fixb_snew_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S5RRRRR14_invdynB_fixb_snew_vp2: pkin has to be [10x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5RRRRR14_invdynB_fixb_snew_vp2: m has to be [6x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [6,3]), ...
  'S5RRRRR14_invdynB_fixb_snew_vp2: mrSges has to be [6x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [6 6]), ...
  'S5RRRRR14_invdynB_fixb_snew_vp2: Ifges has to be [6x6] (double)'); 

%% Symbolic Calculation
% From invdyn_fixb_NewtonEuler_linkframe_tauB_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2024-09-27 18:42:59
% EndTime: 2024-09-27 18:43:14
% DurationCPUTime: 6.20s
% Computational Cost: add. (226229->283), mult. (346876->373), div. (0->0), fcn. (247899->12), ass. (0->123)
t576 = sin(pkin(5));
t608 = pkin(8) * t576;
t574 = qJD(1) + qJD(2);
t571 = t574 ^ 2;
t607 = t571 * t576 ^ 2;
t606 = t574 * t576;
t580 = sin(qJ(3));
t605 = t576 * t580;
t585 = cos(qJ(3));
t604 = t576 * t585;
t577 = cos(pkin(5));
t603 = t577 * t580;
t602 = t577 * t585;
t572 = qJDD(1) + qJDD(2);
t600 = qJD(3) * t574;
t553 = (t572 * t580 + t585 * t600) * t576;
t564 = t577 * t572 + qJDD(3);
t565 = t577 * t574 + qJD(3);
t582 = sin(qJ(1));
t587 = cos(qJ(1));
t566 = t582 * g(1) - g(2) * t587;
t559 = qJDD(1) * pkin(1) + t566;
t567 = -g(1) * t587 - g(2) * t582;
t588 = qJD(1) ^ 2;
t561 = -pkin(1) * t588 + t567;
t581 = sin(qJ(2));
t586 = cos(qJ(2));
t544 = t586 * t559 - t561 * t581;
t535 = pkin(2) * t572 + t571 * t608 + t544;
t545 = t581 * t559 + t586 * t561;
t536 = -pkin(2) * t571 + t572 * t608 + t545;
t592 = t535 * t602 - t580 * t536;
t508 = t564 * pkin(3) - t553 * pkin(9) + (pkin(3) * t580 * t607 + (pkin(9) * t565 * t574 - g(3)) * t576) * t585 + t592;
t518 = -g(3) * t605 + t535 * t603 + t585 * t536;
t598 = t574 * t605;
t552 = pkin(3) * t565 - pkin(9) * t598;
t554 = (t572 * t585 - t580 * t600) * t576;
t599 = t585 ^ 2 * t607;
t509 = -pkin(3) * t599 + pkin(9) * t554 - t552 * t565 + t518;
t579 = sin(qJ(4));
t584 = cos(qJ(4));
t498 = t584 * t508 - t579 * t509;
t547 = (-t579 * t580 + t584 * t585) * t606;
t521 = qJD(4) * t547 + t553 * t584 + t554 * t579;
t548 = (t579 * t585 + t580 * t584) * t606;
t562 = qJDD(4) + t564;
t563 = qJD(4) + t565;
t496 = (t547 * t563 - t521) * pkin(10) + (t547 * t548 + t562) * pkin(4) + t498;
t499 = t579 * t508 + t584 * t509;
t520 = -qJD(4) * t548 - t553 * t579 + t554 * t584;
t539 = pkin(4) * t563 - pkin(10) * t548;
t546 = t547 ^ 2;
t497 = -pkin(4) * t546 + pkin(10) * t520 - t539 * t563 + t499;
t578 = sin(qJ(5));
t583 = cos(qJ(5));
t494 = t496 * t583 - t497 * t578;
t528 = t547 * t583 - t548 * t578;
t504 = qJD(5) * t528 + t520 * t578 + t521 * t583;
t529 = t547 * t578 + t548 * t583;
t515 = -mrSges(6,1) * t528 + mrSges(6,2) * t529;
t560 = qJD(5) + t563;
t522 = -mrSges(6,2) * t560 + mrSges(6,3) * t528;
t558 = qJDD(5) + t562;
t492 = m(6) * t494 + mrSges(6,1) * t558 - mrSges(6,3) * t504 - t515 * t529 + t522 * t560;
t495 = t496 * t578 + t497 * t583;
t503 = -qJD(5) * t529 + t520 * t583 - t521 * t578;
t523 = mrSges(6,1) * t560 - mrSges(6,3) * t529;
t493 = m(6) * t495 - mrSges(6,2) * t558 + mrSges(6,3) * t503 + t515 * t528 - t523 * t560;
t484 = t583 * t492 + t578 * t493;
t530 = -mrSges(5,1) * t547 + mrSges(5,2) * t548;
t537 = -mrSges(5,2) * t563 + mrSges(5,3) * t547;
t482 = m(5) * t498 + mrSges(5,1) * t562 - mrSges(5,3) * t521 - t530 * t548 + t537 * t563 + t484;
t538 = mrSges(5,1) * t563 - mrSges(5,3) * t548;
t593 = -t492 * t578 + t583 * t493;
t483 = m(5) * t499 - mrSges(5,2) * t562 + mrSges(5,3) * t520 + t530 * t547 - t538 * t563 + t593;
t478 = t584 * t482 + t579 * t483;
t517 = -g(3) * t604 + t592;
t597 = t574 * t604;
t550 = -mrSges(4,2) * t565 + mrSges(4,3) * t597;
t551 = (-mrSges(4,1) * t585 + mrSges(4,2) * t580) * t606;
t476 = m(4) * t517 + mrSges(4,1) * t564 - mrSges(4,3) * t553 + t550 * t565 - t551 * t598 + t478;
t549 = mrSges(4,1) * t565 - mrSges(4,3) * t598;
t594 = -t482 * t579 + t584 * t483;
t477 = m(4) * t518 - mrSges(4,2) * t564 + mrSges(4,3) * t554 - t549 * t565 + t551 * t597 + t594;
t531 = -g(3) * t577 - t535 * t576;
t516 = -pkin(3) * t554 - pkin(9) * t599 + t552 * t598 + t531;
t501 = -pkin(4) * t520 - pkin(10) * t546 + t539 * t548 + t516;
t590 = m(6) * t501 - t503 * mrSges(6,1) + t504 * mrSges(6,2) - t528 * t522 + t529 * t523;
t589 = m(5) * t516 - t520 * mrSges(5,1) + t521 * mrSges(5,2) - t547 * t537 + t548 * t538 + t590;
t488 = t589 + (t549 * t580 - t550 * t585) * t606 + t553 * mrSges(4,2) - t554 * mrSges(4,1) + m(4) * t531;
t463 = t476 * t602 + t477 * t603 - t488 * t576;
t461 = m(3) * t544 + mrSges(3,1) * t572 - mrSges(3,2) * t571 + t463;
t467 = -t476 * t580 + t477 * t585;
t466 = m(3) * t545 - mrSges(3,1) * t571 - mrSges(3,2) * t572 + t467;
t457 = t461 * t586 + t466 * t581;
t455 = m(2) * t566 + qJDD(1) * mrSges(2,1) - mrSges(2,2) * t588 + t457;
t595 = -t461 * t581 + t466 * t586;
t456 = m(2) * t567 - mrSges(2,1) * t588 - qJDD(1) * mrSges(2,2) + t595;
t601 = t455 * t587 + t456 * t582;
t462 = t476 * t604 + t477 * t605 + t577 * t488;
t596 = -t455 * t582 + t456 * t587;
t510 = Ifges(6,5) * t529 + Ifges(6,6) * t528 + Ifges(6,3) * t560;
t512 = Ifges(6,1) * t529 + Ifges(6,4) * t528 + Ifges(6,5) * t560;
t485 = -mrSges(6,1) * t501 + mrSges(6,3) * t495 + Ifges(6,4) * t504 + Ifges(6,2) * t503 + Ifges(6,6) * t558 - t510 * t529 + t512 * t560;
t511 = Ifges(6,4) * t529 + Ifges(6,2) * t528 + Ifges(6,6) * t560;
t486 = mrSges(6,2) * t501 - mrSges(6,3) * t494 + Ifges(6,1) * t504 + Ifges(6,4) * t503 + Ifges(6,5) * t558 + t510 * t528 - t511 * t560;
t524 = Ifges(5,5) * t548 + Ifges(5,6) * t547 + Ifges(5,3) * t563;
t526 = Ifges(5,1) * t548 + Ifges(5,4) * t547 + Ifges(5,5) * t563;
t469 = -mrSges(5,1) * t516 + mrSges(5,3) * t499 + Ifges(5,4) * t521 + Ifges(5,2) * t520 + Ifges(5,6) * t562 - pkin(4) * t590 + pkin(10) * t593 + t583 * t485 + t578 * t486 - t548 * t524 + t563 * t526;
t525 = Ifges(5,4) * t548 + Ifges(5,2) * t547 + Ifges(5,6) * t563;
t470 = mrSges(5,2) * t516 - mrSges(5,3) * t498 + Ifges(5,1) * t521 + Ifges(5,4) * t520 + Ifges(5,5) * t562 - pkin(10) * t484 - t485 * t578 + t486 * t583 + t524 * t547 - t525 * t563;
t540 = Ifges(4,3) * t565 + (Ifges(4,5) * t580 + Ifges(4,6) * t585) * t606;
t542 = Ifges(4,5) * t565 + (Ifges(4,1) * t580 + Ifges(4,4) * t585) * t606;
t458 = -mrSges(4,1) * t531 + mrSges(4,3) * t518 + Ifges(4,4) * t553 + Ifges(4,2) * t554 + Ifges(4,6) * t564 - pkin(3) * t589 + pkin(9) * t594 + t584 * t469 + t579 * t470 - t540 * t598 + t565 * t542;
t541 = Ifges(4,6) * t565 + (Ifges(4,4) * t580 + Ifges(4,2) * t585) * t606;
t459 = mrSges(4,2) * t531 - mrSges(4,3) * t517 + Ifges(4,1) * t553 + Ifges(4,4) * t554 + Ifges(4,5) * t564 - pkin(9) * t478 - t469 * t579 + t470 * t584 + t540 * t597 - t541 * t565;
t591 = pkin(8) * t467 + t458 * t585 + t459 * t580;
t468 = Ifges(4,5) * t553 + Ifges(4,6) * t554 + Ifges(4,3) * t564 + mrSges(4,1) * t517 - mrSges(4,2) * t518 + Ifges(5,5) * t521 + Ifges(5,6) * t520 + Ifges(5,3) * t562 + t548 * t525 - t547 * t526 + mrSges(5,1) * t498 - mrSges(5,2) * t499 + Ifges(6,5) * t504 + Ifges(6,6) * t503 + Ifges(6,3) * t558 + t529 * t511 - t528 * t512 + mrSges(6,1) * t494 - mrSges(6,2) * t495 + pkin(4) * t484 + pkin(3) * t478 + (t541 * t580 - t542 * t585) * t606;
t451 = -mrSges(3,2) * g(3) - mrSges(3,3) * t544 + Ifges(3,5) * t572 - t571 * Ifges(3,6) - t580 * t458 + t585 * t459 + (-t462 * t576 - t463 * t577) * pkin(8);
t450 = mrSges(3,1) * g(3) + mrSges(3,3) * t545 + t571 * Ifges(3,5) + Ifges(3,6) * t572 - pkin(2) * t462 - t576 * t468 + t577 * t591;
t449 = -mrSges(2,2) * g(3) - mrSges(2,3) * t566 + Ifges(2,5) * qJDD(1) - Ifges(2,6) * t588 - pkin(7) * t457 - t450 * t581 + t451 * t586;
t448 = Ifges(2,6) * qJDD(1) + t588 * Ifges(2,5) + mrSges(2,1) * g(3) + mrSges(2,3) * t567 + t581 * t451 + t586 * t450 - pkin(1) * (-m(3) * g(3) + t462) + pkin(7) * t595;
t1 = [-m(1) * g(1) + t596; -m(1) * g(2) + t601; (-m(1) - m(2) - m(3)) * g(3) + t462; -mrSges(1,2) * g(3) + mrSges(1,3) * g(2) - pkin(6) * t601 - t448 * t582 + t449 * t587; mrSges(1,1) * g(3) - mrSges(1,3) * g(1) + pkin(6) * t596 + t587 * t448 + t582 * t449; -mrSges(1,1) * g(2) + mrSges(2,1) * t566 + mrSges(3,1) * t544 + mrSges(1,2) * g(1) - mrSges(2,2) * t567 - mrSges(3,2) * t545 + Ifges(2,3) * qJDD(1) + Ifges(3,3) * t572 + pkin(1) * t457 + pkin(2) * t463 + t577 * t468 + t576 * t591;];
tauB = t1;
