% Calculate vector of inverse dynamics joint torques for with Newton-Euler
% S6RRRRPP7
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
% pkin [11x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,alpha2,d1,d2,d3,d4,theta5]';
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
% Datum: 2019-05-07 18:56
% Revision: 89c353f7eff3bd693eda4e29f35b2761dbc3ada0 (2019-05-03)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauJ = S6RRRRPP7_invdynJ_fixb_snew_vp2(qJ, qJD, qJDD, g, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(6,1),zeros(3,1),zeros(11,1),zeros(7,1),zeros(7,3),zeros(7,6)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRRRPP7_invdynJ_fixb_snew_vp2: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RRRRPP7_invdynJ_fixb_snew_vp2: qJD has to be [6x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [6 1]), ...
  'S6RRRRPP7_invdynJ_fixb_snew_vp2: qJDD has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6RRRRPP7_invdynJ_fixb_snew_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [11 1]), ...
  'S6RRRRPP7_invdynJ_fixb_snew_vp2: pkin has to be [11x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6RRRRPP7_invdynJ_fixb_snew_vp2: m has to be [7x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [7,3]), ...
  'S6RRRRPP7_invdynJ_fixb_snew_vp2: mrSges has to be [7x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [7 6]), ...
  'S6RRRRPP7_invdynJ_fixb_snew_vp2: Ifges has to be [7x6] (double)'); 

%% Symbolic Calculation
% From invdyn_fixb_NewtonEuler_linkframe_tauJ_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-05-07 18:41:43
% EndTime: 2019-05-07 18:41:54
% DurationCPUTime: 8.27s
% Computational Cost: add. (123167->338), mult. (264394->430), div. (0->0), fcn. (210590->12), ass. (0->144)
t568 = Ifges(6,1) + Ifges(7,1);
t560 = Ifges(6,4) - Ifges(7,5);
t559 = Ifges(6,5) + Ifges(7,4);
t567 = -Ifges(6,2) - Ifges(7,3);
t566 = -Ifges(7,2) - Ifges(6,3);
t558 = Ifges(6,6) - Ifges(7,6);
t521 = sin(pkin(6));
t525 = sin(qJ(2));
t529 = cos(qJ(2));
t546 = qJD(1) * qJD(2);
t512 = (-qJDD(1) * t529 + t525 * t546) * t521;
t548 = qJD(1) * t521;
t510 = (-pkin(2) * t529 - pkin(9) * t525) * t548;
t522 = cos(pkin(6));
t517 = qJD(1) * t522 + qJD(2);
t515 = t517 ^ 2;
t516 = qJDD(1) * t522 + qJDD(2);
t547 = qJD(1) * t529;
t531 = qJD(1) ^ 2;
t526 = sin(qJ(1));
t530 = cos(qJ(1));
t537 = -g(1) * t530 - g(2) * t526;
t563 = pkin(8) * t521;
t508 = -pkin(1) * t531 + qJDD(1) * t563 + t537;
t541 = t526 * g(1) - g(2) * t530;
t507 = qJDD(1) * pkin(1) + t531 * t563 + t541;
t556 = t507 * t522;
t549 = t529 * t508 + t525 * t556;
t462 = -t515 * pkin(2) + t516 * pkin(9) + (-g(3) * t525 + t510 * t547) * t521 + t549;
t511 = (qJDD(1) * t525 + t529 * t546) * t521;
t562 = t522 * g(3);
t463 = t512 * pkin(2) - t511 * pkin(9) - t562 + (-t507 + (pkin(2) * t525 - pkin(9) * t529) * t517 * qJD(1)) * t521;
t524 = sin(qJ(3));
t528 = cos(qJ(3));
t430 = t528 * t462 + t524 * t463;
t543 = t525 * t548;
t500 = t517 * t528 - t524 * t543;
t501 = t517 * t524 + t528 * t543;
t485 = -pkin(3) * t500 - pkin(10) * t501;
t504 = qJDD(3) + t512;
t542 = t521 * t547;
t514 = qJD(3) - t542;
t513 = t514 ^ 2;
t420 = -pkin(3) * t513 + pkin(10) * t504 + t485 * t500 + t430;
t554 = t521 * t529;
t482 = -g(3) * t554 - t525 * t508 + t529 * t556;
t461 = -t516 * pkin(2) - t515 * pkin(9) + t510 * t543 - t482;
t480 = -qJD(3) * t501 - t511 * t524 + t516 * t528;
t481 = qJD(3) * t500 + t511 * t528 + t516 * t524;
t423 = (-t500 * t514 - t481) * pkin(10) + (t501 * t514 - t480) * pkin(3) + t461;
t523 = sin(qJ(4));
t527 = cos(qJ(4));
t415 = -t523 * t420 + t527 * t423;
t488 = -t501 * t523 + t514 * t527;
t448 = qJD(4) * t488 + t481 * t527 + t504 * t523;
t478 = qJDD(4) - t480;
t489 = t501 * t527 + t514 * t523;
t499 = qJD(4) - t500;
t412 = (t488 * t499 - t448) * qJ(5) + (t488 * t489 + t478) * pkin(4) + t415;
t416 = t527 * t420 + t523 * t423;
t447 = -qJD(4) * t489 - t481 * t523 + t504 * t527;
t470 = pkin(4) * t499 - qJ(5) * t489;
t487 = t488 ^ 2;
t414 = -pkin(4) * t487 + qJ(5) * t447 - t470 * t499 + t416;
t520 = sin(pkin(11));
t557 = cos(pkin(11));
t465 = -t488 * t557 + t489 * t520;
t564 = -2 * qJD(5);
t408 = t520 * t412 + t557 * t414 + t465 * t564;
t427 = -t447 * t557 + t448 * t520;
t466 = t520 * t488 + t489 * t557;
t451 = mrSges(6,1) * t499 - mrSges(6,3) * t466;
t440 = pkin(5) * t465 - qJ(6) * t466;
t498 = t499 ^ 2;
t405 = -pkin(5) * t498 + qJ(6) * t478 + 0.2e1 * qJD(6) * t499 - t440 * t465 + t408;
t452 = -mrSges(7,1) * t499 + mrSges(7,2) * t466;
t544 = m(7) * t405 + t478 * mrSges(7,3) + t499 * t452;
t441 = mrSges(7,1) * t465 - mrSges(7,3) * t466;
t550 = -mrSges(6,1) * t465 - mrSges(6,2) * t466 - t441;
t561 = -mrSges(6,3) - mrSges(7,2);
t394 = m(6) * t408 - t478 * mrSges(6,2) + t427 * t561 - t499 * t451 + t465 * t550 + t544;
t536 = t412 * t557 - t520 * t414;
t407 = t466 * t564 + t536;
t428 = t520 * t447 + t448 * t557;
t450 = -mrSges(6,2) * t499 - mrSges(6,3) * t465;
t406 = -t478 * pkin(5) - t498 * qJ(6) + qJDD(6) + ((2 * qJD(5)) + t440) * t466 - t536;
t449 = -mrSges(7,2) * t465 + mrSges(7,3) * t499;
t538 = -m(7) * t406 + t478 * mrSges(7,1) + t499 * t449;
t396 = m(6) * t407 + t478 * mrSges(6,1) + t428 * t561 + t499 * t450 + t466 * t550 + t538;
t389 = t520 * t394 + t557 * t396;
t402 = t428 * mrSges(7,2) + t466 * t441 - t538;
t454 = Ifges(5,4) * t489 + Ifges(5,2) * t488 + Ifges(5,6) * t499;
t455 = Ifges(5,1) * t489 + Ifges(5,4) * t488 + Ifges(5,5) * t499;
t551 = -t560 * t465 + t568 * t466 + t559 * t499;
t552 = t567 * t465 + t560 * t466 + t558 * t499;
t565 = -t427 * t558 + t428 * t559 + t465 * t551 + t466 * t552 + (Ifges(5,3) - t566) * t478 + mrSges(5,1) * t415 + mrSges(6,1) * t407 - mrSges(7,1) * t406 - mrSges(5,2) * t416 - mrSges(6,2) * t408 + mrSges(7,3) * t405 + Ifges(5,5) * t448 + Ifges(5,6) * t447 + pkin(4) * t389 - pkin(5) * t402 + qJ(6) * (-t427 * mrSges(7,2) - t465 * t441 + t544) + t489 * t454 - t488 * t455;
t555 = t521 * t525;
t467 = -mrSges(5,1) * t488 + mrSges(5,2) * t489;
t469 = -mrSges(5,2) * t499 + mrSges(5,3) * t488;
t387 = m(5) * t415 + mrSges(5,1) * t478 - mrSges(5,3) * t448 - t467 * t489 + t469 * t499 + t389;
t471 = mrSges(5,1) * t499 - mrSges(5,3) * t489;
t539 = t557 * t394 - t396 * t520;
t388 = m(5) * t416 - mrSges(5,2) * t478 + mrSges(5,3) * t447 + t467 * t488 - t471 * t499 + t539;
t385 = -t387 * t523 + t527 * t388;
t484 = -mrSges(4,1) * t500 + mrSges(4,2) * t501;
t491 = mrSges(4,1) * t514 - mrSges(4,3) * t501;
t383 = m(4) * t430 - mrSges(4,2) * t504 + mrSges(4,3) * t480 + t484 * t500 - t491 * t514 + t385;
t429 = -t524 * t462 + t528 * t463;
t419 = -t504 * pkin(3) - t513 * pkin(10) + t501 * t485 - t429;
t417 = -t447 * pkin(4) - t487 * qJ(5) + t489 * t470 + qJDD(5) + t419;
t410 = -0.2e1 * qJD(6) * t466 + (t465 * t499 - t428) * qJ(6) + (t466 * t499 + t427) * pkin(5) + t417;
t403 = m(7) * t410 + t427 * mrSges(7,1) - t428 * mrSges(7,3) + t465 * t449 - t466 * t452;
t400 = m(6) * t417 + t427 * mrSges(6,1) + mrSges(6,2) * t428 + t465 * t450 + t451 * t466 + t403;
t399 = -m(5) * t419 + t447 * mrSges(5,1) - mrSges(5,2) * t448 + t488 * t469 - t471 * t489 - t400;
t490 = -mrSges(4,2) * t514 + mrSges(4,3) * t500;
t398 = m(4) * t429 + mrSges(4,1) * t504 - mrSges(4,3) * t481 - t484 * t501 + t490 * t514 + t399;
t379 = t524 * t383 + t528 * t398;
t553 = t558 * t465 - t559 * t466 + t566 * t499;
t540 = t528 * t383 - t398 * t524;
t384 = t387 * t527 + t388 * t523;
t534 = -m(4) * t461 + t480 * mrSges(4,1) - mrSges(4,2) * t481 + t500 * t490 - t491 * t501 - t384;
t390 = -mrSges(6,1) * t417 - mrSges(7,1) * t410 + mrSges(7,2) * t405 + mrSges(6,3) * t408 - pkin(5) * t403 + t567 * t427 + t560 * t428 + t553 * t466 + t558 * t478 + t551 * t499;
t391 = mrSges(6,2) * t417 + mrSges(7,2) * t406 - mrSges(6,3) * t407 - mrSges(7,3) * t410 - qJ(6) * t403 - t560 * t427 + t568 * t428 + t553 * t465 + t559 * t478 - t552 * t499;
t453 = Ifges(5,5) * t489 + Ifges(5,6) * t488 + Ifges(5,3) * t499;
t376 = -mrSges(5,1) * t419 + mrSges(5,3) * t416 + Ifges(5,4) * t448 + Ifges(5,2) * t447 + Ifges(5,6) * t478 - pkin(4) * t400 + qJ(5) * t539 + t390 * t557 + t520 * t391 - t489 * t453 + t499 * t455;
t377 = mrSges(5,2) * t419 - mrSges(5,3) * t415 + Ifges(5,1) * t448 + Ifges(5,4) * t447 + Ifges(5,5) * t478 - qJ(5) * t389 - t520 * t390 + t391 * t557 + t488 * t453 - t499 * t454;
t475 = Ifges(4,4) * t501 + Ifges(4,2) * t500 + Ifges(4,6) * t514;
t476 = Ifges(4,1) * t501 + Ifges(4,4) * t500 + Ifges(4,5) * t514;
t533 = mrSges(4,1) * t429 - mrSges(4,2) * t430 + Ifges(4,5) * t481 + Ifges(4,6) * t480 + Ifges(4,3) * t504 + pkin(3) * t399 + pkin(10) * t385 + t527 * t376 + t523 * t377 + t501 * t475 - t500 * t476;
t509 = (-mrSges(3,1) * t529 + mrSges(3,2) * t525) * t548;
t506 = -mrSges(3,2) * t517 + mrSges(3,3) * t542;
t505 = mrSges(3,1) * t517 - mrSges(3,3) * t543;
t495 = -t521 * t507 - t562;
t494 = Ifges(3,5) * t517 + (Ifges(3,1) * t525 + Ifges(3,4) * t529) * t548;
t493 = Ifges(3,6) * t517 + (Ifges(3,4) * t525 + Ifges(3,2) * t529) * t548;
t492 = Ifges(3,3) * t517 + (Ifges(3,5) * t525 + Ifges(3,6) * t529) * t548;
t483 = -g(3) * t555 + t549;
t474 = Ifges(4,5) * t501 + Ifges(4,6) * t500 + Ifges(4,3) * t514;
t380 = m(3) * t482 + mrSges(3,1) * t516 - mrSges(3,3) * t511 + t506 * t517 - t509 * t543 + t534;
t378 = m(3) * t483 - mrSges(3,2) * t516 - mrSges(3,3) * t512 - t505 * t517 + t509 * t542 + t540;
t375 = -mrSges(4,1) * t461 + mrSges(4,3) * t430 + Ifges(4,4) * t481 + Ifges(4,2) * t480 + Ifges(4,6) * t504 - pkin(3) * t384 - t501 * t474 + t514 * t476 - t565;
t374 = mrSges(4,2) * t461 - mrSges(4,3) * t429 + Ifges(4,1) * t481 + Ifges(4,4) * t480 + Ifges(4,5) * t504 - pkin(10) * t384 - t376 * t523 + t377 * t527 + t474 * t500 - t475 * t514;
t373 = Ifges(3,5) * t511 - Ifges(3,6) * t512 + Ifges(3,3) * t516 + mrSges(3,1) * t482 - mrSges(3,2) * t483 + t524 * t374 + t528 * t375 + pkin(2) * t534 + pkin(9) * t540 + (t493 * t525 - t494 * t529) * t548;
t1 = [Ifges(2,3) * qJDD(1) + mrSges(2,1) * t541 - mrSges(2,2) * t537 + (mrSges(3,2) * t495 - mrSges(3,3) * t482 + Ifges(3,1) * t511 - Ifges(3,4) * t512 + Ifges(3,5) * t516 - pkin(9) * t379 + t374 * t528 - t375 * t524 + t492 * t542 - t493 * t517) * t555 + (-mrSges(3,1) * t495 + mrSges(3,3) * t483 + Ifges(3,4) * t511 - Ifges(3,2) * t512 + Ifges(3,6) * t516 - pkin(2) * t379 - t492 * t543 + t517 * t494 - t533) * t554 + t522 * t373 + pkin(1) * ((t378 * t525 + t380 * t529) * t522 + (-m(3) * t495 - t512 * mrSges(3,1) - t511 * mrSges(3,2) + (-t505 * t525 + t506 * t529) * t548 - t379) * t521) + (t378 * t529 - t380 * t525) * t563; t373; t533; t565; t400; t402;];
tauJ  = t1;
