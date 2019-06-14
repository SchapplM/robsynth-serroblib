% Calculate vector of inverse dynamics joint torques for with Newton-Euler
% S6RRPRPR9
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
%   pkin=[a2,a3,a4,a5,a6,alpha2,d1,d2,d4,d6,theta3,theta5]';
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
% Datum: 2019-05-06 15:23
% Revision: 89c353f7eff3bd693eda4e29f35b2761dbc3ada0 (2019-05-03)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauJ = S6RRPRPR9_invdynJ_fixb_snew_vp2(qJ, qJD, qJDD, g, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(6,1),zeros(3,1),zeros(12,1),zeros(7,1),zeros(7,3),zeros(7,6)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRPRPR9_invdynJ_fixb_snew_vp2: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RRPRPR9_invdynJ_fixb_snew_vp2: qJD has to be [6x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [6 1]), ...
  'S6RRPRPR9_invdynJ_fixb_snew_vp2: qJDD has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6RRPRPR9_invdynJ_fixb_snew_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [12 1]), ...
  'S6RRPRPR9_invdynJ_fixb_snew_vp2: pkin has to be [12x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6RRPRPR9_invdynJ_fixb_snew_vp2: m has to be [7x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [7,3]), ...
  'S6RRPRPR9_invdynJ_fixb_snew_vp2: mrSges has to be [7x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [7 6]), ...
  'S6RRPRPR9_invdynJ_fixb_snew_vp2: Ifges has to be [7x6] (double)'); 

%% Symbolic Calculation
% From invdyn_fixb_NewtonEuler_linkframe_tauJ_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-05-06 15:14:35
% EndTime: 2019-05-06 15:14:51
% DurationCPUTime: 14.72s
% Computational Cost: add. (226601->360), mult. (518653->472), div. (0->0), fcn. (422125->14), ass. (0->149)
t573 = cos(qJ(4));
t540 = sin(pkin(6));
t572 = pkin(8) * t540;
t543 = cos(pkin(6));
t571 = t543 * g(3);
t551 = qJD(1) ^ 2;
t547 = sin(qJ(1));
t550 = cos(qJ(1));
t561 = t547 * g(1) - g(2) * t550;
t524 = qJDD(1) * pkin(1) + t551 * t572 + t561;
t570 = t524 * t543;
t546 = sin(qJ(2));
t569 = t540 * t546;
t549 = cos(qJ(2));
t568 = t540 * t549;
t566 = qJD(1) * t540;
t526 = (-pkin(2) * t549 - qJ(3) * t546) * t566;
t535 = qJD(1) * t543 + qJD(2);
t533 = t535 ^ 2;
t534 = qJDD(1) * t543 + qJDD(2);
t565 = qJD(1) * t549;
t557 = -g(1) * t550 - g(2) * t547;
t564 = qJDD(1) * t540;
t525 = -pkin(1) * t551 + pkin(8) * t564 + t557;
t567 = t549 * t525 + t546 * t570;
t480 = -t533 * pkin(2) + t534 * qJ(3) + (-g(3) * t546 + t526 * t565) * t540 + t567;
t528 = (qJD(2) * t565 + qJDD(1) * t546) * t540;
t563 = t546 * t566;
t529 = -qJD(2) * t563 + t549 * t564;
t481 = -t529 * pkin(2) - t571 - t528 * qJ(3) + (-t524 + (pkin(2) * t546 - qJ(3) * t549) * t535 * qJD(1)) * t540;
t539 = sin(pkin(11));
t542 = cos(pkin(11));
t518 = t535 * t539 + t542 * t563;
t448 = -0.2e1 * qJD(3) * t518 - t539 * t480 + t542 * t481;
t505 = t528 * t542 + t534 * t539;
t517 = t535 * t542 - t539 * t563;
t562 = t540 * t565;
t439 = (-t517 * t562 - t505) * pkin(9) + (t517 * t518 - t529) * pkin(3) + t448;
t449 = 0.2e1 * qJD(3) * t517 + t542 * t480 + t539 * t481;
t504 = -t528 * t539 + t534 * t542;
t506 = -pkin(3) * t562 - pkin(9) * t518;
t516 = t517 ^ 2;
t446 = -pkin(3) * t516 + pkin(9) * t504 + t506 * t562 + t449;
t545 = sin(qJ(4));
t431 = t545 * t439 + t573 * t446;
t497 = -t573 * t517 + t518 * t545;
t498 = t545 * t517 + t573 * t518;
t474 = pkin(4) * t497 - qJ(5) * t498;
t521 = qJDD(4) - t529;
t531 = qJD(4) - t562;
t530 = t531 ^ 2;
t429 = -pkin(4) * t530 + qJ(5) * t521 - t474 * t497 + t431;
t494 = -g(3) * t568 - t546 * t525 + t549 * t570;
t479 = -t534 * pkin(2) - t533 * qJ(3) + t526 * t563 + qJDD(3) - t494;
t450 = -t504 * pkin(3) - t516 * pkin(9) + t518 * t506 + t479;
t466 = qJD(4) * t498 - t573 * t504 + t505 * t545;
t467 = -t497 * qJD(4) + t545 * t504 + t573 * t505;
t434 = (t497 * t531 - t467) * qJ(5) + (t498 * t531 + t466) * pkin(4) + t450;
t538 = sin(pkin(12));
t541 = cos(pkin(12));
t486 = t498 * t541 + t531 * t538;
t424 = -0.2e1 * qJD(5) * t486 - t538 * t429 + t541 * t434;
t460 = t467 * t541 + t521 * t538;
t485 = -t498 * t538 + t531 * t541;
t422 = (t485 * t497 - t460) * pkin(10) + (t485 * t486 + t466) * pkin(5) + t424;
t425 = 0.2e1 * qJD(5) * t485 + t541 * t429 + t538 * t434;
t459 = -t467 * t538 + t521 * t541;
t469 = pkin(5) * t497 - pkin(10) * t486;
t484 = t485 ^ 2;
t423 = -pkin(5) * t484 + pkin(10) * t459 - t469 * t497 + t425;
t544 = sin(qJ(6));
t548 = cos(qJ(6));
t420 = t422 * t548 - t423 * t544;
t461 = t485 * t548 - t486 * t544;
t437 = qJD(6) * t461 + t459 * t544 + t460 * t548;
t462 = t485 * t544 + t486 * t548;
t447 = -mrSges(7,1) * t461 + mrSges(7,2) * t462;
t496 = qJD(6) + t497;
t451 = -mrSges(7,2) * t496 + mrSges(7,3) * t461;
t465 = qJDD(6) + t466;
t417 = m(7) * t420 + mrSges(7,1) * t465 - mrSges(7,3) * t437 - t447 * t462 + t451 * t496;
t421 = t422 * t544 + t423 * t548;
t436 = -qJD(6) * t462 + t459 * t548 - t460 * t544;
t452 = mrSges(7,1) * t496 - mrSges(7,3) * t462;
t418 = m(7) * t421 - mrSges(7,2) * t465 + mrSges(7,3) * t436 + t447 * t461 - t452 * t496;
t409 = t548 * t417 + t544 * t418;
t463 = -mrSges(6,1) * t485 + mrSges(6,2) * t486;
t556 = -mrSges(6,2) * t497 + mrSges(6,3) * t485;
t407 = m(6) * t424 + t466 * mrSges(6,1) - t460 * mrSges(6,3) - t486 * t463 + t497 * t556 + t409;
t468 = mrSges(6,1) * t497 - mrSges(6,3) * t486;
t558 = -t417 * t544 + t548 * t418;
t408 = m(6) * t425 - mrSges(6,2) * t466 + mrSges(6,3) * t459 + t463 * t485 - t468 * t497 + t558;
t403 = -t407 * t538 + t541 * t408;
t475 = mrSges(5,1) * t497 + mrSges(5,2) * t498;
t488 = mrSges(5,1) * t531 - mrSges(5,3) * t498;
t400 = m(5) * t431 - mrSges(5,2) * t521 - mrSges(5,3) * t466 - t475 * t497 - t488 * t531 + t403;
t430 = t573 * t439 - t545 * t446;
t428 = -t521 * pkin(4) - t530 * qJ(5) + t498 * t474 + qJDD(5) - t430;
t426 = -t459 * pkin(5) - t484 * pkin(10) + t486 * t469 + t428;
t555 = m(7) * t426 - t436 * mrSges(7,1) + mrSges(7,2) * t437 - t461 * t451 + t452 * t462;
t419 = m(6) * t428 - t459 * mrSges(6,1) + mrSges(6,2) * t460 + t468 * t486 - t485 * t556 + t555;
t487 = -mrSges(5,2) * t531 - mrSges(5,3) * t497;
t413 = m(5) * t430 + mrSges(5,1) * t521 - mrSges(5,3) * t467 - t475 * t498 + t487 * t531 - t419;
t394 = t545 * t400 + t573 * t413;
t499 = -mrSges(4,1) * t517 + mrSges(4,2) * t518;
t502 = mrSges(4,2) * t562 + mrSges(4,3) * t517;
t392 = m(4) * t448 - mrSges(4,1) * t529 - mrSges(4,3) * t505 - t499 * t518 - t502 * t562 + t394;
t503 = -mrSges(4,1) * t562 - mrSges(4,3) * t518;
t559 = t573 * t400 - t413 * t545;
t393 = m(4) * t449 + mrSges(4,2) * t529 + mrSges(4,3) * t504 + t499 * t517 + t503 * t562 + t559;
t387 = t542 * t392 + t539 * t393;
t402 = t541 * t407 + t538 * t408;
t560 = -t392 * t539 + t542 * t393;
t554 = m(5) * t450 + t466 * mrSges(5,1) + mrSges(5,2) * t467 + t497 * t487 + t488 * t498 + t402;
t441 = Ifges(7,4) * t462 + Ifges(7,2) * t461 + Ifges(7,6) * t496;
t442 = Ifges(7,1) * t462 + Ifges(7,4) * t461 + Ifges(7,5) * t496;
t553 = mrSges(7,1) * t420 - mrSges(7,2) * t421 + Ifges(7,5) * t437 + Ifges(7,6) * t436 + Ifges(7,3) * t465 + t462 * t441 - t461 * t442;
t401 = m(4) * t479 - t504 * mrSges(4,1) + mrSges(4,2) * t505 - t517 * t502 + t503 * t518 + t554;
t440 = Ifges(7,5) * t462 + Ifges(7,6) * t461 + Ifges(7,3) * t496;
t410 = -mrSges(7,1) * t426 + mrSges(7,3) * t421 + Ifges(7,4) * t437 + Ifges(7,2) * t436 + Ifges(7,6) * t465 - t440 * t462 + t442 * t496;
t411 = mrSges(7,2) * t426 - mrSges(7,3) * t420 + Ifges(7,1) * t437 + Ifges(7,4) * t436 + Ifges(7,5) * t465 + t440 * t461 - t441 * t496;
t453 = Ifges(6,5) * t486 + Ifges(6,6) * t485 + Ifges(6,3) * t497;
t455 = Ifges(6,1) * t486 + Ifges(6,4) * t485 + Ifges(6,5) * t497;
t395 = -mrSges(6,1) * t428 + mrSges(6,3) * t425 + Ifges(6,4) * t460 + Ifges(6,2) * t459 + Ifges(6,6) * t466 - pkin(5) * t555 + pkin(10) * t558 + t548 * t410 + t544 * t411 - t486 * t453 + t497 * t455;
t454 = Ifges(6,4) * t486 + Ifges(6,2) * t485 + Ifges(6,6) * t497;
t396 = mrSges(6,2) * t428 - mrSges(6,3) * t424 + Ifges(6,1) * t460 + Ifges(6,4) * t459 + Ifges(6,5) * t466 - pkin(10) * t409 - t410 * t544 + t411 * t548 + t453 * t485 - t454 * t497;
t471 = Ifges(5,4) * t498 - Ifges(5,2) * t497 + Ifges(5,6) * t531;
t472 = Ifges(5,1) * t498 - Ifges(5,4) * t497 + Ifges(5,5) * t531;
t552 = mrSges(5,1) * t430 - mrSges(5,2) * t431 + Ifges(5,5) * t467 - Ifges(5,6) * t466 + Ifges(5,3) * t521 - pkin(4) * t419 + qJ(5) * t403 + t541 * t395 + t538 * t396 + t498 * t471 + t497 * t472;
t527 = (-mrSges(3,1) * t549 + mrSges(3,2) * t546) * t566;
t523 = -mrSges(3,2) * t535 + mrSges(3,3) * t562;
t522 = mrSges(3,1) * t535 - mrSges(3,3) * t563;
t510 = -t540 * t524 - t571;
t509 = Ifges(3,5) * t535 + (Ifges(3,1) * t546 + Ifges(3,4) * t549) * t566;
t508 = Ifges(3,6) * t535 + (t546 * Ifges(3,4) + Ifges(3,2) * t549) * t566;
t507 = Ifges(3,3) * t535 + (Ifges(3,5) * t546 + Ifges(3,6) * t549) * t566;
t495 = -g(3) * t569 + t567;
t491 = Ifges(4,1) * t518 + Ifges(4,4) * t517 - Ifges(4,5) * t562;
t490 = Ifges(4,4) * t518 + Ifges(4,2) * t517 - Ifges(4,6) * t562;
t489 = Ifges(4,5) * t518 + Ifges(4,6) * t517 - Ifges(4,3) * t562;
t470 = Ifges(5,5) * t498 - Ifges(5,6) * t497 + Ifges(5,3) * t531;
t397 = m(3) * t494 + mrSges(3,1) * t534 - mrSges(3,3) * t528 + t523 * t535 - t527 * t563 - t401;
t388 = -t498 * t470 + t485 * t455 - t486 * t454 + Ifges(5,4) * t467 - Ifges(6,6) * t459 - Ifges(6,5) * t460 - mrSges(5,1) * t450 + mrSges(5,3) * t431 + mrSges(6,2) * t425 - mrSges(6,1) * t424 - pkin(5) * t409 - t553 + (-Ifges(5,2) - Ifges(6,3)) * t466 - pkin(4) * t402 + t531 * t472 + Ifges(5,6) * t521;
t386 = m(3) * t495 - mrSges(3,2) * t534 + mrSges(3,3) * t529 - t522 * t535 + t527 * t562 + t560;
t385 = mrSges(5,2) * t450 - mrSges(5,3) * t430 + Ifges(5,1) * t467 - Ifges(5,4) * t466 + Ifges(5,5) * t521 - qJ(5) * t402 - t395 * t538 + t396 * t541 - t470 * t497 - t471 * t531;
t384 = mrSges(4,2) * t479 - mrSges(4,3) * t448 + Ifges(4,1) * t505 + Ifges(4,4) * t504 - Ifges(4,5) * t529 - pkin(9) * t394 + t573 * t385 - t545 * t388 + t517 * t489 + t490 * t562;
t383 = -mrSges(4,1) * t479 + mrSges(4,3) * t449 + Ifges(4,4) * t505 + Ifges(4,2) * t504 - Ifges(4,6) * t529 - pkin(3) * t554 + pkin(9) * t559 + t545 * t385 + t573 * t388 - t518 * t489 - t491 * t562;
t382 = Ifges(3,5) * t528 + Ifges(3,6) * t529 + Ifges(3,3) * t534 + mrSges(3,1) * t494 - mrSges(3,2) * t495 + t539 * t384 + t542 * t383 - pkin(2) * t401 + qJ(3) * t560 + (t508 * t546 - t509 * t549) * t566;
t1 = [Ifges(2,3) * qJDD(1) + mrSges(2,1) * t561 - mrSges(2,2) * t557 + (mrSges(3,2) * t510 - mrSges(3,3) * t494 + Ifges(3,1) * t528 + Ifges(3,4) * t529 + Ifges(3,5) * t534 - qJ(3) * t387 - t383 * t539 + t384 * t542 + t507 * t562 - t508 * t535) * t569 + (-Ifges(4,5) * t505 - mrSges(3,1) * t510 - pkin(2) * t387 + t517 * t491 - Ifges(4,6) * t504 + mrSges(3,3) * t495 - mrSges(4,1) * t448 + mrSges(4,2) * t449 - t507 * t563 - t552 + (Ifges(3,2) + Ifges(4,3)) * t529 + t535 * t509 + Ifges(3,6) * t534 + Ifges(3,4) * t528 - t518 * t490 - pkin(3) * t394) * t568 + t543 * t382 + pkin(1) * ((t386 * t546 + t397 * t549) * t543 + (-m(3) * t510 + t529 * mrSges(3,1) - t528 * mrSges(3,2) + (-t522 * t546 + t523 * t549) * t566 - t387) * t540) + (t386 * t549 - t397 * t546) * t572; t382; t401; t552; t419; t553;];
tauJ  = t1;
