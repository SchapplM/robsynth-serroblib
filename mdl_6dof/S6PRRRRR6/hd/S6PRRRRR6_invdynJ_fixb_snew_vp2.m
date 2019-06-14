% Calculate vector of inverse dynamics joint torques for with Newton-Euler
% S6PRRRRR6
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
% pkin [14x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,alpha2,alpha3,alpha4,d2,d3,d4,d5,d6,theta1]';
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
% Datum: 2019-05-05 12:50
% Revision: 89c353f7eff3bd693eda4e29f35b2761dbc3ada0 (2019-05-03)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauJ = S6PRRRRR6_invdynJ_fixb_snew_vp2(qJ, qJD, qJDD, g, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(6,1),zeros(3,1),zeros(14,1),zeros(7,1),zeros(7,3),zeros(7,6)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6PRRRRR6_invdynJ_fixb_snew_vp2: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6PRRRRR6_invdynJ_fixb_snew_vp2: qJD has to be [6x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [6 1]), ...
  'S6PRRRRR6_invdynJ_fixb_snew_vp2: qJDD has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6PRRRRR6_invdynJ_fixb_snew_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [14 1]), ...
  'S6PRRRRR6_invdynJ_fixb_snew_vp2: pkin has to be [14x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6PRRRRR6_invdynJ_fixb_snew_vp2: m has to be [7x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [7,3]), ...
  'S6PRRRRR6_invdynJ_fixb_snew_vp2: mrSges has to be [7x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [7 6]), ...
  'S6PRRRRR6_invdynJ_fixb_snew_vp2: Ifges has to be [7x6] (double)'); 

%% Symbolic Calculation
% From invdyn_fixb_NewtonEuler_linkframe_tauJ_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-05-05 12:35:18
% EndTime: 2019-05-05 12:35:31
% DurationCPUTime: 12.14s
% Computational Cost: add. (196459->319), mult. (462189->435), div. (0->0), fcn. (390390->18), ass. (0->150)
t514 = sin(pkin(14));
t518 = cos(pkin(14));
t507 = g(1) * t514 - g(2) * t518;
t513 = -g(3) + qJDD(1);
t517 = sin(pkin(6));
t521 = cos(pkin(6));
t563 = t507 * t521 + t513 * t517;
t508 = -g(1) * t518 - g(2) * t514;
t526 = sin(qJ(2));
t531 = cos(qJ(2));
t482 = t531 * t508 + t563 * t526;
t532 = qJD(2) ^ 2;
t516 = sin(pkin(7));
t559 = t516 * pkin(10);
t479 = -pkin(2) * t532 + qJDD(2) * t559 + t482;
t520 = cos(pkin(7));
t512 = qJD(2) * t520 + qJD(3);
t515 = sin(pkin(8));
t519 = cos(pkin(8));
t530 = cos(qJ(3));
t547 = qJD(2) * t516;
t543 = t530 * t547;
t494 = (t512 * t515 + t519 * t543) * pkin(11);
t525 = sin(qJ(3));
t560 = pkin(11) * t525;
t498 = (-t530 * pkin(3) - t515 * t560) * t547;
t546 = qJD(2) * qJD(3);
t502 = (qJDD(2) * t525 + t530 * t546) * t516;
t511 = qJDD(2) * t520 + qJDD(3);
t481 = -t508 * t526 + t563 * t531;
t478 = qJDD(2) * pkin(2) + t532 * t559 + t481;
t550 = t520 * t530;
t497 = -t507 * t517 + t513 * t521;
t558 = t497 * t516;
t548 = t478 * t550 + t530 * t558;
t561 = pkin(11) * t519;
t437 = -t502 * t561 + pkin(3) * t511 + t494 * t512 + (-t498 * t547 - t479) * t525 + t548;
t551 = t520 * t525;
t455 = t478 * t551 + t530 * t479 + t525 * t558;
t544 = t525 * t547;
t496 = pkin(3) * t512 - t544 * t561;
t503 = (qJDD(2) * t530 - t525 * t546) * t516;
t539 = t503 * t519 + t511 * t515;
t438 = pkin(11) * t539 - t496 * t512 + t498 * t543 + t455;
t492 = t520 * t497;
t562 = pkin(11) * t515;
t451 = -t502 * t562 - pkin(3) * t503 + t492 + (-t478 + (-t494 * t530 + t496 * t525) * qJD(2)) * t516;
t524 = sin(qJ(4));
t529 = cos(qJ(4));
t423 = -t524 * t438 + (t437 * t519 + t451 * t515) * t529;
t552 = t519 * t530;
t555 = t515 * t524;
t485 = t512 * t555 + (t524 * t552 + t525 * t529) * t547;
t465 = -qJD(4) * t485 - t502 * t524 + t529 * t539;
t554 = t515 * t529;
t484 = (-t525 * t524 + t529 * t552) * t547 + t512 * t554;
t553 = t519 * t524;
t424 = t437 * t553 + t529 * t438 + t451 * t555;
t467 = -mrSges(5,1) * t484 + mrSges(5,2) * t485;
t495 = t512 * t519 - t515 * t543 + qJD(4);
t474 = mrSges(5,1) * t495 - mrSges(5,3) * t485;
t486 = -t503 * t515 + t511 * t519 + qJDD(4);
t468 = -pkin(4) * t484 - pkin(12) * t485;
t493 = t495 ^ 2;
t420 = -pkin(4) * t493 + pkin(12) * t486 + t468 * t484 + t424;
t428 = -t437 * t515 + t519 * t451;
t466 = qJD(4) * t484 + t502 * t529 + t524 * t539;
t422 = (-t484 * t495 - t466) * pkin(12) + (t485 * t495 - t465) * pkin(4) + t428;
t523 = sin(qJ(5));
t528 = cos(qJ(5));
t416 = t528 * t420 + t523 * t422;
t471 = -t485 * t523 + t495 * t528;
t472 = t485 * t528 + t495 * t523;
t453 = -pkin(5) * t471 - pkin(13) * t472;
t464 = qJDD(5) - t465;
t483 = qJD(5) - t484;
t480 = t483 ^ 2;
t414 = -pkin(5) * t480 + pkin(13) * t464 + t453 * t471 + t416;
t419 = -pkin(4) * t486 - pkin(12) * t493 + t485 * t468 - t423;
t441 = -qJD(5) * t472 - t466 * t523 + t486 * t528;
t442 = qJD(5) * t471 + t466 * t528 + t486 * t523;
t417 = (-t471 * t483 - t442) * pkin(13) + (t472 * t483 - t441) * pkin(5) + t419;
t522 = sin(qJ(6));
t527 = cos(qJ(6));
t411 = -t414 * t522 + t417 * t527;
t457 = -t472 * t522 + t483 * t527;
t427 = qJD(6) * t457 + t442 * t527 + t464 * t522;
t458 = t472 * t527 + t483 * t522;
t433 = -mrSges(7,1) * t457 + mrSges(7,2) * t458;
t440 = qJDD(6) - t441;
t470 = qJD(6) - t471;
t443 = -mrSges(7,2) * t470 + mrSges(7,3) * t457;
t408 = m(7) * t411 + mrSges(7,1) * t440 - mrSges(7,3) * t427 - t433 * t458 + t443 * t470;
t412 = t414 * t527 + t417 * t522;
t426 = -qJD(6) * t458 - t442 * t522 + t464 * t527;
t444 = mrSges(7,1) * t470 - mrSges(7,3) * t458;
t409 = m(7) * t412 - mrSges(7,2) * t440 + mrSges(7,3) * t426 + t433 * t457 - t444 * t470;
t402 = -t408 * t522 + t527 * t409;
t452 = -mrSges(6,1) * t471 + mrSges(6,2) * t472;
t460 = mrSges(6,1) * t483 - mrSges(6,3) * t472;
t400 = m(6) * t416 - mrSges(6,2) * t464 + mrSges(6,3) * t441 + t452 * t471 - t460 * t483 + t402;
t415 = -t420 * t523 + t422 * t528;
t413 = -pkin(5) * t464 - pkin(13) * t480 + t453 * t472 - t415;
t410 = -m(7) * t413 + t426 * mrSges(7,1) - mrSges(7,2) * t427 + t457 * t443 - t444 * t458;
t459 = -mrSges(6,2) * t483 + mrSges(6,3) * t471;
t406 = m(6) * t415 + mrSges(6,1) * t464 - mrSges(6,3) * t442 - t452 * t472 + t459 * t483 + t410;
t541 = t528 * t400 - t406 * t523;
t391 = m(5) * t424 - mrSges(5,2) * t486 + mrSges(5,3) * t465 + t467 * t484 - t474 * t495 + t541;
t394 = t523 * t400 + t528 * t406;
t473 = -mrSges(5,2) * t495 + mrSges(5,3) * t484;
t393 = m(5) * t428 - mrSges(5,1) * t465 + mrSges(5,2) * t466 - t473 * t484 + t474 * t485 + t394;
t401 = t408 * t527 + t409 * t522;
t535 = -m(6) * t419 + t441 * mrSges(6,1) - mrSges(6,2) * t442 + t471 * t459 - t460 * t472 - t401;
t397 = m(5) * t423 + mrSges(5,1) * t486 - mrSges(5,3) * t466 - t467 * t485 + t473 * t495 + t535;
t381 = t519 * t529 * t397 + t391 * t553 - t393 * t515;
t454 = -t479 * t525 + t548;
t500 = -mrSges(4,2) * t512 + mrSges(4,3) * t543;
t501 = (-t530 * mrSges(4,1) + t525 * mrSges(4,2)) * t547;
t378 = m(4) * t454 + mrSges(4,1) * t511 - mrSges(4,3) * t502 + t500 * t512 - t501 * t544 + t381;
t385 = t529 * t391 - t397 * t524;
t499 = mrSges(4,1) * t512 - mrSges(4,3) * t544;
t384 = m(4) * t455 - mrSges(4,2) * t511 + mrSges(4,3) * t503 - t499 * t512 + t501 * t543 + t385;
t549 = t378 * t550 + t384 * t551;
t380 = t391 * t555 + t519 * t393 + t397 * t554;
t542 = -t378 * t525 + t530 * t384;
t429 = Ifges(7,5) * t458 + Ifges(7,6) * t457 + Ifges(7,3) * t470;
t431 = Ifges(7,1) * t458 + Ifges(7,4) * t457 + Ifges(7,5) * t470;
t403 = -mrSges(7,1) * t413 + mrSges(7,3) * t412 + Ifges(7,4) * t427 + Ifges(7,2) * t426 + Ifges(7,6) * t440 - t429 * t458 + t431 * t470;
t430 = Ifges(7,4) * t458 + Ifges(7,2) * t457 + Ifges(7,6) * t470;
t404 = mrSges(7,2) * t413 - mrSges(7,3) * t411 + Ifges(7,1) * t427 + Ifges(7,4) * t426 + Ifges(7,5) * t440 + t429 * t457 - t430 * t470;
t447 = Ifges(6,5) * t472 + Ifges(6,6) * t471 + Ifges(6,3) * t483;
t448 = Ifges(6,4) * t472 + Ifges(6,2) * t471 + Ifges(6,6) * t483;
t386 = mrSges(6,2) * t419 - mrSges(6,3) * t415 + Ifges(6,1) * t442 + Ifges(6,4) * t441 + Ifges(6,5) * t464 - pkin(13) * t401 - t403 * t522 + t404 * t527 + t447 * t471 - t448 * t483;
t449 = Ifges(6,1) * t472 + Ifges(6,4) * t471 + Ifges(6,5) * t483;
t534 = mrSges(7,1) * t411 - mrSges(7,2) * t412 + Ifges(7,5) * t427 + Ifges(7,6) * t426 + Ifges(7,3) * t440 + t430 * t458 - t457 * t431;
t387 = -mrSges(6,1) * t419 + mrSges(6,3) * t416 + Ifges(6,4) * t442 + Ifges(6,2) * t441 + Ifges(6,6) * t464 - pkin(5) * t401 - t447 * t472 + t449 * t483 - t534;
t461 = Ifges(5,5) * t485 + Ifges(5,6) * t484 + Ifges(5,3) * t495;
t462 = Ifges(5,4) * t485 + Ifges(5,2) * t484 + Ifges(5,6) * t495;
t375 = mrSges(5,2) * t428 - mrSges(5,3) * t423 + Ifges(5,1) * t466 + Ifges(5,4) * t465 + Ifges(5,5) * t486 - pkin(12) * t394 + t386 * t528 - t387 * t523 + t461 * t484 - t462 * t495;
t463 = Ifges(5,1) * t485 + Ifges(5,4) * t484 + Ifges(5,5) * t495;
t533 = mrSges(6,1) * t415 - mrSges(6,2) * t416 + Ifges(6,5) * t442 + Ifges(6,6) * t441 + Ifges(6,3) * t464 + pkin(5) * t410 + pkin(13) * t402 + t527 * t403 + t522 * t404 + t472 * t448 - t471 * t449;
t376 = -mrSges(5,1) * t428 + mrSges(5,3) * t424 + Ifges(5,4) * t466 + Ifges(5,2) * t465 + Ifges(5,6) * t486 - pkin(4) * t394 - t485 * t461 + t495 * t463 - t533;
t536 = pkin(11) * t385 + t375 * t524 + t376 * t529;
t491 = Ifges(4,5) * t512 + (t525 * Ifges(4,1) + t530 * Ifges(4,4)) * t547;
t490 = Ifges(4,6) * t512 + (t525 * Ifges(4,4) + t530 * Ifges(4,2)) * t547;
t469 = -t478 * t516 + t492;
t379 = m(4) * t469 - mrSges(4,1) * t503 + mrSges(4,2) * t502 + (t499 * t525 - t500 * t530) * t547 + t380;
t374 = mrSges(5,1) * t423 - mrSges(5,2) * t424 + Ifges(5,5) * t466 + Ifges(5,6) * t465 + Ifges(5,3) * t486 + pkin(4) * t535 + pkin(12) * t541 + t523 * t386 + t528 * t387 + t485 * t462 - t484 * t463;
t373 = mrSges(4,1) * t454 - mrSges(4,2) * t455 + Ifges(4,5) * t502 + Ifges(4,6) * t503 + Ifges(4,3) * t511 + pkin(3) * t381 + t374 * t519 + (t525 * t490 - t530 * t491) * t547 + t536 * t515;
t1 = [m(2) * t513 + t521 * (m(3) * t497 + t379 * t520 + (t378 * t530 + t384 * t525) * t516) + (t526 * (m(3) * t482 - mrSges(3,1) * t532 - qJDD(2) * mrSges(3,2) + t542) + t531 * (m(3) * t481 + qJDD(2) * mrSges(3,1) - mrSges(3,2) * t532 - t379 * t516 + t549)) * t517; Ifges(3,3) * qJDD(2) + mrSges(3,1) * t481 - mrSges(3,2) * t482 + t520 * t373 + pkin(2) * t549 + (t525 * (mrSges(4,2) * t469 - mrSges(4,3) * t454 + Ifges(4,1) * t502 + Ifges(4,4) * t503 + Ifges(4,5) * t511 + t375 * t529 - t376 * t524 - t380 * t562 - t490 * t512) + t530 * (-mrSges(4,1) * t469 + mrSges(4,3) * t455 + Ifges(4,4) * t502 + Ifges(4,2) * t503 + Ifges(4,6) * t511 - pkin(3) * t380 - t374 * t515 + t491 * t512) - pkin(2) * t379 + pkin(10) * t542 + (-t381 * t560 + t530 * t536) * t519) * t516; t373; t374; t533; t534;];
tauJ  = t1;
