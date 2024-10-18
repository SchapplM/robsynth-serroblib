% Calculate vector of inverse dynamics base forces with Newton-Euler for
% S5PRRRR11
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
%   pkin=[a2,a3,a4,a5,alpha3,d2,d3,d4,d5,theta1]';
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
% Datum: 2024-09-27 21:46
% Revision: 4128c0dfad36ceb08eaa7a9fbc8f797c888f59be (2024-07-29)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauB = S5PRRRR11_invdynB_fixb_snew_vp2(qJ, qJD, qJDD, g, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(5,1),zeros(3,1),zeros(10,1),zeros(6,1),zeros(6,3),zeros(6,6)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5PRRRR11_invdynB_fixb_snew_vp2: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5PRRRR11_invdynB_fixb_snew_vp2: qJD has to be [5x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [5 1]), ...
  'S5PRRRR11_invdynB_fixb_snew_vp2: qJDD has to be [5x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S5PRRRR11_invdynB_fixb_snew_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S5PRRRR11_invdynB_fixb_snew_vp2: pkin has to be [10x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5PRRRR11_invdynB_fixb_snew_vp2: m has to be [6x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [6,3]), ...
  'S5PRRRR11_invdynB_fixb_snew_vp2: mrSges has to be [6x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [6 6]), ...
  'S5PRRRR11_invdynB_fixb_snew_vp2: Ifges has to be [6x6] (double)'); 

%% Symbolic Calculation
% From invdyn_fixb_NewtonEuler_linkframe_tauB_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2024-09-27 21:45:35
% EndTime: 2024-09-27 21:45:41
% DurationCPUTime: 2.81s
% Computational Cost: add. (131329->269), mult. (326724->359), div. (0->0), fcn. (247899->12), ass. (0->119)
t558 = sin(pkin(5));
t589 = pkin(7) * t558;
t569 = qJD(2) ^ 2;
t588 = t558 ^ 2 * t569;
t563 = sin(qJ(3));
t587 = t558 * t563;
t567 = cos(qJ(3));
t586 = t558 * t567;
t560 = cos(pkin(5));
t585 = t560 * t563;
t584 = t560 * t567;
t557 = sin(pkin(10));
t559 = cos(pkin(10));
t542 = t557 * g(1) - t559 * g(2);
t543 = -t559 * g(1) - t557 * g(2);
t564 = sin(qJ(2));
t568 = cos(qJ(2));
t525 = t568 * t542 - t564 * t543;
t517 = qJDD(2) * pkin(2) + t569 * t589 + t525;
t526 = t564 * t542 + t568 * t543;
t518 = -t569 * pkin(2) + qJDD(2) * t589 + t526;
t556 = -g(3) + qJDD(1);
t499 = t517 * t584 - t563 * t518 + t556 * t586;
t581 = qJD(2) * qJD(3);
t535 = (qJDD(2) * t563 + t567 * t581) * t558;
t548 = t560 * qJDD(2) + qJDD(3);
t549 = t560 * qJD(2) + qJD(3);
t582 = qJD(2) * t558;
t578 = t567 * t582;
t490 = (t549 * t578 - t535) * pkin(8) + (t563 * t567 * t588 + t548) * pkin(3) + t499;
t500 = t517 * t585 + t567 * t518 + t556 * t587;
t579 = t563 * t582;
t534 = t549 * pkin(3) - pkin(8) * t579;
t536 = (qJDD(2) * t567 - t563 * t581) * t558;
t580 = t567 ^ 2 * t588;
t491 = -pkin(3) * t580 + t536 * pkin(8) - t549 * t534 + t500;
t562 = sin(qJ(4));
t566 = cos(qJ(4));
t480 = t566 * t490 - t562 * t491;
t529 = (-t562 * t563 + t566 * t567) * t582;
t503 = t529 * qJD(4) + t566 * t535 + t562 * t536;
t530 = (t562 * t567 + t563 * t566) * t582;
t546 = qJDD(4) + t548;
t547 = qJD(4) + t549;
t478 = (t529 * t547 - t503) * pkin(9) + (t529 * t530 + t546) * pkin(4) + t480;
t481 = t562 * t490 + t566 * t491;
t502 = -t530 * qJD(4) - t562 * t535 + t566 * t536;
t521 = t547 * pkin(4) - t530 * pkin(9);
t528 = t529 ^ 2;
t479 = -t528 * pkin(4) + t502 * pkin(9) - t547 * t521 + t481;
t561 = sin(qJ(5));
t565 = cos(qJ(5));
t476 = t565 * t478 - t561 * t479;
t510 = t565 * t529 - t561 * t530;
t486 = t510 * qJD(5) + t561 * t502 + t565 * t503;
t511 = t561 * t529 + t565 * t530;
t496 = -t510 * mrSges(6,1) + t511 * mrSges(6,2);
t541 = qJD(5) + t547;
t504 = -t541 * mrSges(6,2) + t510 * mrSges(6,3);
t540 = qJDD(5) + t546;
t474 = m(6) * t476 + t540 * mrSges(6,1) - t486 * mrSges(6,3) - t511 * t496 + t541 * t504;
t477 = t561 * t478 + t565 * t479;
t485 = -t511 * qJD(5) + t565 * t502 - t561 * t503;
t505 = t541 * mrSges(6,1) - t511 * mrSges(6,3);
t475 = m(6) * t477 - t540 * mrSges(6,2) + t485 * mrSges(6,3) + t510 * t496 - t541 * t505;
t466 = t565 * t474 + t561 * t475;
t513 = -t529 * mrSges(5,1) + t530 * mrSges(5,2);
t519 = -t547 * mrSges(5,2) + t529 * mrSges(5,3);
t464 = m(5) * t480 + t546 * mrSges(5,1) - t503 * mrSges(5,3) - t530 * t513 + t547 * t519 + t466;
t520 = t547 * mrSges(5,1) - t530 * mrSges(5,3);
t574 = -t561 * t474 + t565 * t475;
t465 = m(5) * t481 - t546 * mrSges(5,2) + t502 * mrSges(5,3) + t529 * t513 - t547 * t520 + t574;
t460 = t566 * t464 + t562 * t465;
t532 = -t549 * mrSges(4,2) + mrSges(4,3) * t578;
t533 = (-mrSges(4,1) * t567 + mrSges(4,2) * t563) * t582;
t458 = m(4) * t499 + t548 * mrSges(4,1) - t535 * mrSges(4,3) + t549 * t532 - t533 * t579 + t460;
t531 = t549 * mrSges(4,1) - mrSges(4,3) * t579;
t575 = -t562 * t464 + t566 * t465;
t459 = m(4) * t500 - t548 * mrSges(4,2) + t536 * mrSges(4,3) - t549 * t531 + t533 * t578 + t575;
t512 = -t558 * t517 + t560 * t556;
t498 = -t536 * pkin(3) - pkin(8) * t580 + t534 * t579 + t512;
t483 = -t502 * pkin(4) - t528 * pkin(9) + t530 * t521 + t498;
t571 = m(6) * t483 - t485 * mrSges(6,1) + t486 * mrSges(6,2) - t510 * t504 + t511 * t505;
t570 = m(5) * t498 - t502 * mrSges(5,1) + t503 * mrSges(5,2) - t529 * t519 + t530 * t520 + t571;
t470 = t570 + (t531 * t563 - t532 * t567) * t582 + t535 * mrSges(4,2) - t536 * mrSges(4,1) + m(4) * t512;
t445 = t458 * t584 + t459 * t585 - t558 * t470;
t443 = m(3) * t525 + qJDD(2) * mrSges(3,1) - t569 * mrSges(3,2) + t445;
t449 = -t563 * t458 + t567 * t459;
t448 = m(3) * t526 - t569 * mrSges(3,1) - qJDD(2) * mrSges(3,2) + t449;
t439 = t568 * t443 + t564 * t448;
t437 = m(2) * t542 + t439;
t576 = -t564 * t443 + t568 * t448;
t438 = m(2) * t543 + t576;
t583 = t559 * t437 + t557 * t438;
t444 = t458 * t586 + t459 * t587 + t560 * t470;
t577 = -t557 * t437 + t559 * t438;
t573 = m(3) * t556 + t444;
t492 = Ifges(6,5) * t511 + Ifges(6,6) * t510 + Ifges(6,3) * t541;
t494 = Ifges(6,1) * t511 + Ifges(6,4) * t510 + Ifges(6,5) * t541;
t467 = -mrSges(6,1) * t483 + mrSges(6,3) * t477 + Ifges(6,4) * t486 + Ifges(6,2) * t485 + Ifges(6,6) * t540 - t511 * t492 + t541 * t494;
t493 = Ifges(6,4) * t511 + Ifges(6,2) * t510 + Ifges(6,6) * t541;
t468 = mrSges(6,2) * t483 - mrSges(6,3) * t476 + Ifges(6,1) * t486 + Ifges(6,4) * t485 + Ifges(6,5) * t540 + t510 * t492 - t541 * t493;
t506 = Ifges(5,5) * t530 + Ifges(5,6) * t529 + Ifges(5,3) * t547;
t508 = Ifges(5,1) * t530 + Ifges(5,4) * t529 + Ifges(5,5) * t547;
t451 = -mrSges(5,1) * t498 + mrSges(5,3) * t481 + Ifges(5,4) * t503 + Ifges(5,2) * t502 + Ifges(5,6) * t546 - pkin(4) * t571 + pkin(9) * t574 + t565 * t467 + t561 * t468 - t530 * t506 + t547 * t508;
t507 = Ifges(5,4) * t530 + Ifges(5,2) * t529 + Ifges(5,6) * t547;
t452 = mrSges(5,2) * t498 - mrSges(5,3) * t480 + Ifges(5,1) * t503 + Ifges(5,4) * t502 + Ifges(5,5) * t546 - pkin(9) * t466 - t561 * t467 + t565 * t468 + t529 * t506 - t547 * t507;
t522 = Ifges(4,3) * t549 + (Ifges(4,5) * t563 + Ifges(4,6) * t567) * t582;
t524 = Ifges(4,5) * t549 + (Ifges(4,1) * t563 + Ifges(4,4) * t567) * t582;
t440 = -mrSges(4,1) * t512 + mrSges(4,3) * t500 + Ifges(4,4) * t535 + Ifges(4,2) * t536 + Ifges(4,6) * t548 - pkin(3) * t570 + pkin(8) * t575 + t566 * t451 + t562 * t452 - t522 * t579 + t549 * t524;
t523 = Ifges(4,6) * t549 + (Ifges(4,4) * t563 + Ifges(4,2) * t567) * t582;
t441 = mrSges(4,2) * t512 - mrSges(4,3) * t499 + Ifges(4,1) * t535 + Ifges(4,4) * t536 + Ifges(4,5) * t548 - pkin(8) * t460 - t562 * t451 + t566 * t452 + t522 * t578 - t549 * t523;
t572 = pkin(7) * t449 + t440 * t567 + t441 * t563;
t450 = Ifges(4,5) * t535 + Ifges(4,6) * t536 + Ifges(4,3) * t548 + mrSges(4,1) * t499 - mrSges(4,2) * t500 + Ifges(5,5) * t503 + Ifges(5,6) * t502 + Ifges(5,3) * t546 + t530 * t507 - t529 * t508 + mrSges(5,1) * t480 - mrSges(5,2) * t481 + Ifges(6,5) * t486 + Ifges(6,6) * t485 + Ifges(6,3) * t540 + t511 * t493 - t510 * t494 + mrSges(6,1) * t476 - mrSges(6,2) * t477 + pkin(4) * t466 + pkin(3) * t460 + (t523 * t563 - t524 * t567) * t582;
t433 = mrSges(3,2) * t556 - mrSges(3,3) * t525 + Ifges(3,5) * qJDD(2) - t569 * Ifges(3,6) - t563 * t440 + t567 * t441 + (-t444 * t558 - t445 * t560) * pkin(7);
t432 = -mrSges(3,1) * t556 + mrSges(3,3) * t526 + t569 * Ifges(3,5) + Ifges(3,6) * qJDD(2) - pkin(2) * t444 - t558 * t450 + t560 * t572;
t431 = mrSges(2,2) * t556 - mrSges(2,3) * t542 - pkin(6) * t439 - t564 * t432 + t568 * t433;
t430 = -mrSges(2,1) * t556 + mrSges(2,3) * t543 - pkin(1) * t573 + pkin(6) * t576 + t568 * t432 + t564 * t433;
t1 = [-m(1) * g(1) + t577; -m(1) * g(2) + t583; -m(1) * g(3) + m(2) * t556 + t573; -mrSges(1,2) * g(3) + mrSges(1,3) * g(2) - qJ(1) * t583 - t557 * t430 + t559 * t431; mrSges(1,1) * g(3) - mrSges(1,3) * g(1) + qJ(1) * t577 + t559 * t430 + t557 * t431; -mrSges(1,1) * g(2) + mrSges(2,1) * t542 + mrSges(3,1) * t525 + mrSges(1,2) * g(1) - mrSges(2,2) * t543 - mrSges(3,2) * t526 + Ifges(3,3) * qJDD(2) + pkin(1) * t439 + pkin(2) * t445 + t560 * t450 + t572 * t558;];
tauB = t1;
