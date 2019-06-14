% Calculate vector of inverse dynamics joint torques for with Newton-Euler
% S6RRRRRR5
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
%   pkin=[a2,a3,a4,a5,a6,alpha2,d1,d2,d3,d4,d5,d6]';
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
% tauJ [6x1]
%   joint torques of inverse dynamics (contains inertial, gravitational coriolis and centrifugal forces)

% Quelle: HybrDyn-Toolbox
% Datum: 2019-05-08 10:25
% Revision: 89c353f7eff3bd693eda4e29f35b2761dbc3ada0 (2019-05-03)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauJ = S6RRRRRR5_invdynJ_fixb_snew_vp2(qJ, qJD, qJDD, g, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(6,1),zeros(3,1),zeros(12,1),zeros(7,1),zeros(7,3),zeros(7,6)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRRRRR5_invdynJ_fixb_snew_vp2: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RRRRRR5_invdynJ_fixb_snew_vp2: qJD has to be [6x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [6 1]), ...
  'S6RRRRRR5_invdynJ_fixb_snew_vp2: qJDD has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6RRRRRR5_invdynJ_fixb_snew_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [12 1]), ...
  'S6RRRRRR5_invdynJ_fixb_snew_vp2: pkin has to be [12x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6RRRRRR5_invdynJ_fixb_snew_vp2: m has to be [7x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [7,3]), ...
  'S6RRRRRR5_invdynJ_fixb_snew_vp2: mrSges has to be [7x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [7 6]), ...
  'S6RRRRRR5_invdynJ_fixb_snew_vp2: Ifges has to be [7x6] (double)'); 

%% Symbolic Calculation
% From invdyn_fixb_NewtonEuler_linkframe_tauJ_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-05-08 10:07:45
% EndTime: 2019-05-08 10:08:12
% DurationCPUTime: 19.24s
% Computational Cost: add. (307137->364), mult. (665361->470), div. (0->0), fcn. (549900->14), ass. (0->155)
t552 = sin(pkin(6));
t591 = pkin(8) * t552;
t553 = cos(pkin(6));
t590 = g(3) * t553;
t566 = qJD(1) ^ 2;
t559 = sin(qJ(1));
t565 = cos(qJ(1));
t580 = t559 * g(1) - g(2) * t565;
t536 = qJDD(1) * pkin(1) + t566 * t591 + t580;
t589 = t536 * t553;
t558 = sin(qJ(2));
t588 = t552 * t558;
t564 = cos(qJ(2));
t587 = t552 * t564;
t585 = qJD(1) * t552;
t539 = (-t564 * pkin(2) - t558 * pkin(9)) * t585;
t549 = qJD(1) * t553 + qJD(2);
t547 = t549 ^ 2;
t548 = qJDD(1) * t553 + qJDD(2);
t584 = qJD(1) * t564;
t575 = -g(1) * t565 - g(2) * t559;
t583 = qJDD(1) * t552;
t537 = -pkin(1) * t566 + pkin(8) * t583 + t575;
t586 = t564 * t537 + t558 * t589;
t496 = -pkin(2) * t547 + pkin(9) * t548 + (-g(3) * t558 + t539 * t584) * t552 + t586;
t540 = (qJD(2) * t584 + qJDD(1) * t558) * t552;
t582 = t558 * t585;
t541 = -qJD(2) * t582 + t564 * t583;
t497 = -pkin(2) * t541 - pkin(9) * t540 - t590 + (-t536 + (pkin(2) * t558 - pkin(9) * t564) * t549 * qJD(1)) * t552;
t557 = sin(qJ(3));
t563 = cos(qJ(3));
t471 = -t496 * t557 + t563 * t497;
t528 = t549 * t563 - t557 * t582;
t508 = qJD(3) * t528 + t540 * t563 + t548 * t557;
t529 = t549 * t557 + t563 * t582;
t533 = qJDD(3) - t541;
t581 = t552 * t584;
t545 = qJD(3) - t581;
t456 = (t528 * t545 - t508) * pkin(10) + (t528 * t529 + t533) * pkin(3) + t471;
t472 = t563 * t496 + t557 * t497;
t507 = -qJD(3) * t529 - t540 * t557 + t548 * t563;
t517 = pkin(3) * t545 - pkin(10) * t529;
t527 = t528 ^ 2;
t459 = -pkin(3) * t527 + pkin(10) * t507 - t517 * t545 + t472;
t556 = sin(qJ(4));
t562 = cos(qJ(4));
t436 = t562 * t456 - t459 * t556;
t512 = t528 * t562 - t529 * t556;
t477 = qJD(4) * t512 + t507 * t556 + t508 * t562;
t513 = t528 * t556 + t529 * t562;
t532 = qJDD(4) + t533;
t544 = qJD(4) + t545;
t432 = (t512 * t544 - t477) * pkin(11) + (t512 * t513 + t532) * pkin(4) + t436;
t437 = t556 * t456 + t562 * t459;
t476 = -qJD(4) * t513 + t507 * t562 - t508 * t556;
t500 = pkin(4) * t544 - pkin(11) * t513;
t511 = t512 ^ 2;
t434 = -pkin(4) * t511 + pkin(11) * t476 - t500 * t544 + t437;
t555 = sin(qJ(5));
t561 = cos(qJ(5));
t429 = t555 * t432 + t561 * t434;
t490 = t512 * t555 + t513 * t561;
t448 = -qJD(5) * t490 + t476 * t561 - t477 * t555;
t489 = t512 * t561 - t513 * t555;
t469 = -mrSges(6,1) * t489 + mrSges(6,2) * t490;
t543 = qJD(5) + t544;
t482 = mrSges(6,1) * t543 - mrSges(6,3) * t490;
t526 = qJDD(5) + t532;
t470 = -pkin(5) * t489 - pkin(12) * t490;
t542 = t543 ^ 2;
t426 = -pkin(5) * t542 + pkin(12) * t526 + t470 * t489 + t429;
t509 = -g(3) * t587 - t558 * t537 + t564 * t589;
t495 = -pkin(2) * t548 - pkin(9) * t547 + t539 * t582 - t509;
t463 = -pkin(3) * t507 - pkin(10) * t527 + t529 * t517 + t495;
t442 = -pkin(4) * t476 - pkin(11) * t511 + t513 * t500 + t463;
t449 = qJD(5) * t489 + t476 * t555 + t477 * t561;
t430 = (t490 * t543 - t448) * pkin(5) + (-t489 * t543 - t449) * pkin(12) + t442;
t554 = sin(qJ(6));
t560 = cos(qJ(6));
t423 = -t426 * t554 + t430 * t560;
t479 = -t490 * t554 + t543 * t560;
t440 = qJD(6) * t479 + t449 * t560 + t526 * t554;
t447 = qJDD(6) - t448;
t480 = t490 * t560 + t543 * t554;
t460 = -mrSges(7,1) * t479 + mrSges(7,2) * t480;
t488 = qJD(6) - t489;
t461 = -mrSges(7,2) * t488 + mrSges(7,3) * t479;
t419 = m(7) * t423 + mrSges(7,1) * t447 - mrSges(7,3) * t440 - t460 * t480 + t461 * t488;
t424 = t426 * t560 + t430 * t554;
t439 = -qJD(6) * t480 - t449 * t554 + t526 * t560;
t462 = mrSges(7,1) * t488 - mrSges(7,3) * t480;
t420 = m(7) * t424 - mrSges(7,2) * t447 + mrSges(7,3) * t439 + t460 * t479 - t462 * t488;
t576 = -t419 * t554 + t560 * t420;
t406 = m(6) * t429 - mrSges(6,2) * t526 + mrSges(6,3) * t448 + t469 * t489 - t482 * t543 + t576;
t428 = t432 * t561 - t434 * t555;
t481 = -mrSges(6,2) * t543 + mrSges(6,3) * t489;
t425 = -pkin(5) * t526 - pkin(12) * t542 + t470 * t490 - t428;
t573 = -m(7) * t425 + t439 * mrSges(7,1) - mrSges(7,2) * t440 + t479 * t461 - t462 * t480;
t415 = m(6) * t428 + mrSges(6,1) * t526 - mrSges(6,3) * t449 - t469 * t490 + t481 * t543 + t573;
t402 = t555 * t406 + t561 * t415;
t491 = -mrSges(5,1) * t512 + mrSges(5,2) * t513;
t498 = -mrSges(5,2) * t544 + mrSges(5,3) * t512;
t399 = m(5) * t436 + mrSges(5,1) * t532 - mrSges(5,3) * t477 - t491 * t513 + t498 * t544 + t402;
t499 = mrSges(5,1) * t544 - mrSges(5,3) * t513;
t577 = t561 * t406 - t415 * t555;
t400 = m(5) * t437 - mrSges(5,2) * t532 + mrSges(5,3) * t476 + t491 * t512 - t499 * t544 + t577;
t393 = t562 * t399 + t556 * t400;
t514 = -mrSges(4,1) * t528 + mrSges(4,2) * t529;
t515 = -mrSges(4,2) * t545 + mrSges(4,3) * t528;
t391 = m(4) * t471 + mrSges(4,1) * t533 - mrSges(4,3) * t508 - t514 * t529 + t515 * t545 + t393;
t516 = mrSges(4,1) * t545 - mrSges(4,3) * t529;
t578 = -t399 * t556 + t562 * t400;
t392 = m(4) * t472 - mrSges(4,2) * t533 + mrSges(4,3) * t507 + t514 * t528 - t516 * t545 + t578;
t386 = t563 * t391 + t557 * t392;
t408 = t560 * t419 + t554 * t420;
t579 = -t391 * t557 + t563 * t392;
t574 = m(6) * t442 - t448 * mrSges(6,1) + t449 * mrSges(6,2) - t489 * t481 + t490 * t482 + t408;
t450 = Ifges(7,5) * t480 + Ifges(7,6) * t479 + Ifges(7,3) * t488;
t452 = Ifges(7,1) * t480 + Ifges(7,4) * t479 + Ifges(7,5) * t488;
t412 = -mrSges(7,1) * t425 + mrSges(7,3) * t424 + Ifges(7,4) * t440 + Ifges(7,2) * t439 + Ifges(7,6) * t447 - t450 * t480 + t452 * t488;
t451 = Ifges(7,4) * t480 + Ifges(7,2) * t479 + Ifges(7,6) * t488;
t413 = mrSges(7,2) * t425 - mrSges(7,3) * t423 + Ifges(7,1) * t440 + Ifges(7,4) * t439 + Ifges(7,5) * t447 + t450 * t479 - t451 * t488;
t465 = Ifges(6,4) * t490 + Ifges(6,2) * t489 + Ifges(6,6) * t543;
t466 = Ifges(6,1) * t490 + Ifges(6,4) * t489 + Ifges(6,5) * t543;
t572 = -mrSges(6,1) * t428 + mrSges(6,2) * t429 - Ifges(6,5) * t449 - Ifges(6,6) * t448 - Ifges(6,3) * t526 - pkin(5) * t573 - pkin(12) * t576 - t560 * t412 - t554 * t413 - t490 * t465 + t489 * t466;
t571 = m(5) * t463 - t476 * mrSges(5,1) + t477 * mrSges(5,2) - t512 * t498 + t513 * t499 + t574;
t570 = mrSges(7,1) * t423 - mrSges(7,2) * t424 + Ifges(7,5) * t440 + Ifges(7,6) * t439 + Ifges(7,3) * t447 + t451 * t480 - t452 * t479;
t484 = Ifges(5,4) * t513 + Ifges(5,2) * t512 + Ifges(5,6) * t544;
t485 = Ifges(5,1) * t513 + Ifges(5,4) * t512 + Ifges(5,5) * t544;
t569 = -mrSges(5,1) * t436 + mrSges(5,2) * t437 - Ifges(5,5) * t477 - Ifges(5,6) * t476 - Ifges(5,3) * t532 - pkin(4) * t402 - t513 * t484 + t512 * t485 + t572;
t568 = -m(4) * t495 + t507 * mrSges(4,1) - t508 * mrSges(4,2) + t528 * t515 - t529 * t516 - t571;
t502 = Ifges(4,4) * t529 + Ifges(4,2) * t528 + Ifges(4,6) * t545;
t503 = Ifges(4,1) * t529 + Ifges(4,4) * t528 + Ifges(4,5) * t545;
t567 = mrSges(4,1) * t471 - mrSges(4,2) * t472 + Ifges(4,5) * t508 + Ifges(4,6) * t507 + Ifges(4,3) * t533 + pkin(3) * t393 + t529 * t502 - t528 * t503 - t569;
t538 = (-t564 * mrSges(3,1) + t558 * mrSges(3,2)) * t585;
t535 = -mrSges(3,2) * t549 + mrSges(3,3) * t581;
t534 = mrSges(3,1) * t549 - mrSges(3,3) * t582;
t521 = -t536 * t552 - t590;
t520 = Ifges(3,5) * t549 + (t558 * Ifges(3,1) + t564 * Ifges(3,4)) * t585;
t519 = Ifges(3,6) * t549 + (t558 * Ifges(3,4) + t564 * Ifges(3,2)) * t585;
t518 = Ifges(3,3) * t549 + (t558 * Ifges(3,5) + t564 * Ifges(3,6)) * t585;
t510 = -g(3) * t588 + t586;
t501 = Ifges(4,5) * t529 + Ifges(4,6) * t528 + Ifges(4,3) * t545;
t483 = Ifges(5,5) * t513 + Ifges(5,6) * t512 + Ifges(5,3) * t544;
t464 = Ifges(6,5) * t490 + Ifges(6,6) * t489 + Ifges(6,3) * t543;
t403 = m(3) * t509 + t548 * mrSges(3,1) - t540 * mrSges(3,3) + t549 * t535 - t538 * t582 + t568;
t395 = -mrSges(6,1) * t442 + mrSges(6,3) * t429 + Ifges(6,4) * t449 + Ifges(6,2) * t448 + Ifges(6,6) * t526 - pkin(5) * t408 - t464 * t490 + t466 * t543 - t570;
t394 = mrSges(6,2) * t442 - mrSges(6,3) * t428 + Ifges(6,1) * t449 + Ifges(6,4) * t448 + Ifges(6,5) * t526 - pkin(12) * t408 - t412 * t554 + t413 * t560 + t464 * t489 - t465 * t543;
t387 = mrSges(5,2) * t463 - mrSges(5,3) * t436 + Ifges(5,1) * t477 + Ifges(5,4) * t476 + Ifges(5,5) * t532 - pkin(11) * t402 + t394 * t561 - t395 * t555 + t483 * t512 - t484 * t544;
t385 = m(3) * t510 - mrSges(3,2) * t548 + mrSges(3,3) * t541 - t534 * t549 + t538 * t581 + t579;
t384 = -mrSges(5,1) * t463 + mrSges(5,3) * t437 + Ifges(5,4) * t477 + Ifges(5,2) * t476 + Ifges(5,6) * t532 - pkin(4) * t574 + pkin(11) * t577 + t555 * t394 + t561 * t395 - t513 * t483 + t544 * t485;
t383 = mrSges(4,2) * t495 - mrSges(4,3) * t471 + Ifges(4,1) * t508 + Ifges(4,4) * t507 + Ifges(4,5) * t533 - pkin(10) * t393 - t384 * t556 + t387 * t562 + t501 * t528 - t502 * t545;
t382 = -mrSges(4,1) * t495 + mrSges(4,3) * t472 + Ifges(4,4) * t508 + Ifges(4,2) * t507 + Ifges(4,6) * t533 - pkin(3) * t571 + pkin(10) * t578 + t562 * t384 + t556 * t387 - t529 * t501 + t545 * t503;
t381 = Ifges(3,5) * t540 + Ifges(3,6) * t541 + Ifges(3,3) * t548 + mrSges(3,1) * t509 - mrSges(3,2) * t510 + t557 * t383 + t563 * t382 + pkin(2) * t568 + pkin(9) * t579 + (t558 * t519 - t564 * t520) * t585;
t1 = [Ifges(2,3) * qJDD(1) + mrSges(2,1) * t580 - mrSges(2,2) * t575 + (mrSges(3,2) * t521 - mrSges(3,3) * t509 + Ifges(3,1) * t540 + Ifges(3,4) * t541 + Ifges(3,5) * t548 - pkin(9) * t386 - t557 * t382 + t563 * t383 + t518 * t581 - t549 * t519) * t588 + (-mrSges(3,1) * t521 + mrSges(3,3) * t510 + Ifges(3,4) * t540 + Ifges(3,2) * t541 + Ifges(3,6) * t548 - pkin(2) * t386 - t518 * t582 + t549 * t520 - t567) * t587 + t553 * t381 + pkin(1) * ((t385 * t558 + t403 * t564) * t553 + (-m(3) * t521 + t541 * mrSges(3,1) - t540 * mrSges(3,2) + (-t534 * t558 + t535 * t564) * t585 - t386) * t552) + (t564 * t385 - t558 * t403) * t591; t381; t567; -t569; -t572; t570;];
tauJ  = t1;
