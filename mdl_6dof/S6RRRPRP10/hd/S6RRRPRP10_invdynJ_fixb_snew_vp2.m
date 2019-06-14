% Calculate vector of inverse dynamics joint torques for with Newton-Euler
% S6RRRPRP10
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
%   pkin=[a2,a3,a4,a5,a6,alpha2,d1,d2,d3,d5,theta4]';
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
% Datum: 2019-05-07 09:06
% Revision: 89c353f7eff3bd693eda4e29f35b2761dbc3ada0 (2019-05-03)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauJ = S6RRRPRP10_invdynJ_fixb_snew_vp2(qJ, qJD, qJDD, g, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(6,1),zeros(3,1),zeros(11,1),zeros(7,1),zeros(7,3),zeros(7,6)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRRPRP10_invdynJ_fixb_snew_vp2: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RRRPRP10_invdynJ_fixb_snew_vp2: qJD has to be [6x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [6 1]), ...
  'S6RRRPRP10_invdynJ_fixb_snew_vp2: qJDD has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6RRRPRP10_invdynJ_fixb_snew_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [11 1]), ...
  'S6RRRPRP10_invdynJ_fixb_snew_vp2: pkin has to be [11x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6RRRPRP10_invdynJ_fixb_snew_vp2: m has to be [7x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [7,3]), ...
  'S6RRRPRP10_invdynJ_fixb_snew_vp2: mrSges has to be [7x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [7 6]), ...
  'S6RRRPRP10_invdynJ_fixb_snew_vp2: Ifges has to be [7x6] (double)'); 

%% Symbolic Calculation
% From invdyn_fixb_NewtonEuler_linkframe_tauJ_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-05-07 08:52:23
% EndTime: 2019-05-07 08:52:42
% DurationCPUTime: 8.51s
% Computational Cost: add. (117040->338), mult. (255146->430), div. (0->0), fcn. (202614->12), ass. (0->142)
t569 = Ifges(6,1) + Ifges(7,1);
t560 = Ifges(6,4) - Ifges(7,5);
t559 = -Ifges(6,5) - Ifges(7,4);
t568 = Ifges(6,2) + Ifges(7,3);
t558 = Ifges(6,6) - Ifges(7,6);
t567 = -Ifges(6,3) - Ifges(7,2);
t522 = sin(pkin(6));
t527 = sin(qJ(2));
t529 = cos(qJ(2));
t546 = qJD(1) * qJD(2);
t511 = (-qJDD(1) * t529 + t527 * t546) * t522;
t524 = cos(pkin(6));
t518 = qJD(1) * t524 + qJD(2);
t526 = sin(qJ(3));
t548 = qJD(1) * t522;
t544 = t527 * t548;
t565 = cos(qJ(3));
t500 = t518 * t526 + t544 * t565;
t547 = qJD(1) * t529;
t543 = t522 * t547;
t515 = qJD(3) - t543;
t521 = sin(pkin(11));
t523 = cos(pkin(11));
t487 = -t500 * t521 + t515 * t523;
t488 = t500 * t523 + t515 * t521;
t525 = sin(qJ(5));
t564 = cos(qJ(5));
t460 = -t487 * t564 + t488 * t525;
t499 = -t518 * t565 + t526 * t544;
t510 = (qJDD(1) * t527 + t529 * t546) * t522;
t517 = qJDD(1) * t524 + qJDD(2);
t478 = -qJD(3) * t499 + t510 * t565 + t517 * t526;
t503 = qJDD(3) + t511;
t463 = -t478 * t521 + t503 * t523;
t464 = t478 * t523 + t503 * t521;
t428 = -qJD(5) * t460 + t463 * t525 + t464 * t564;
t461 = t487 * t525 + t488 * t564;
t440 = mrSges(7,1) * t460 - mrSges(7,3) * t461;
t509 = (-pkin(2) * t529 - pkin(9) * t527) * t548;
t516 = t518 ^ 2;
t531 = qJD(1) ^ 2;
t528 = sin(qJ(1));
t530 = cos(qJ(1));
t538 = -g(1) * t530 - g(2) * t528;
t563 = pkin(8) * t522;
t507 = -pkin(1) * t531 + qJDD(1) * t563 + t538;
t542 = g(1) * t528 - g(2) * t530;
t506 = qJDD(1) * pkin(1) + t531 * t563 + t542;
t556 = t506 * t524;
t549 = t507 * t529 + t527 * t556;
t457 = -pkin(2) * t516 + pkin(9) * t517 + (-g(3) * t527 + t509 * t547) * t522 + t549;
t562 = g(3) * t524;
t458 = pkin(2) * t511 - pkin(9) * t510 - t562 + (-t506 + (pkin(2) * t527 - pkin(9) * t529) * t518 * qJD(1)) * t522;
t430 = t457 * t565 + t458 * t526;
t481 = pkin(3) * t499 - qJ(4) * t500;
t514 = t515 ^ 2;
t421 = -pkin(3) * t514 + qJ(4) * t503 - t481 * t499 + t430;
t554 = t522 * t529;
t479 = -g(3) * t554 - t507 * t527 + t529 * t556;
t456 = -pkin(2) * t517 - pkin(9) * t516 + t509 * t544 - t479;
t477 = qJD(3) * t500 + t510 * t526 - t517 * t565;
t424 = (t499 * t515 - t478) * qJ(4) + (t500 * t515 + t477) * pkin(3) + t456;
t416 = -0.2e1 * qJD(4) * t488 - t421 * t521 + t424 * t523;
t413 = (t487 * t499 - t464) * pkin(10) + (t487 * t488 + t477) * pkin(4) + t416;
t417 = 0.2e1 * qJD(4) * t487 + t421 * t523 + t424 * t521;
t468 = pkin(4) * t499 - pkin(10) * t488;
t486 = t487 ^ 2;
t415 = -pkin(4) * t486 + pkin(10) * t463 - t468 * t499 + t417;
t409 = t413 * t564 - t415 * t525;
t439 = pkin(5) * t460 - qJ(6) * t461;
t475 = qJDD(5) + t477;
t498 = qJD(5) + t499;
t497 = t498 ^ 2;
t408 = -pkin(5) * t475 - qJ(6) * t497 + t439 * t461 + qJDD(6) - t409;
t444 = -mrSges(7,2) * t460 + mrSges(7,3) * t498;
t539 = -m(7) * t408 + mrSges(7,1) * t475 + t444 * t498;
t404 = mrSges(7,2) * t428 + t440 * t461 - t539;
t410 = t413 * t525 + t415 * t564;
t407 = -pkin(5) * t497 + qJ(6) * t475 + 0.2e1 * qJD(6) * t498 - t439 * t460 + t410;
t427 = qJD(5) * t461 - t463 * t564 + t464 * t525;
t447 = -mrSges(7,1) * t498 + mrSges(7,2) * t461;
t545 = m(7) * t407 + mrSges(7,3) * t475 + t447 * t498;
t551 = t460 * t560 - t461 * t569 + t498 * t559;
t553 = t460 * t568 - t461 * t560 - t498 * t558;
t566 = t558 * t427 + t559 * t428 + t551 * t460 + t553 * t461 + t567 * t475 - mrSges(6,1) * t409 + mrSges(7,1) * t408 + mrSges(6,2) * t410 - mrSges(7,3) * t407 + pkin(5) * t404 - qJ(6) * (-mrSges(7,2) * t427 - t440 * t460 + t545);
t561 = -mrSges(6,3) - mrSges(7,2);
t555 = t522 * t527;
t446 = mrSges(6,1) * t498 - mrSges(6,3) * t461;
t550 = -mrSges(6,1) * t460 - mrSges(6,2) * t461 - t440;
t397 = m(6) * t410 - mrSges(6,2) * t475 + t427 * t561 - t446 * t498 + t460 * t550 + t545;
t445 = -mrSges(6,2) * t498 - mrSges(6,3) * t460;
t399 = m(6) * t409 + mrSges(6,1) * t475 + t428 * t561 + t445 * t498 + t461 * t550 + t539;
t392 = t397 * t525 + t399 * t564;
t465 = -mrSges(5,1) * t487 + mrSges(5,2) * t488;
t537 = -mrSges(5,2) * t499 + mrSges(5,3) * t487;
t390 = m(5) * t416 + t477 * mrSges(5,1) - t464 * mrSges(5,3) - t488 * t465 + t499 * t537 + t392;
t467 = mrSges(5,1) * t499 - mrSges(5,3) * t488;
t540 = t397 * t564 - t399 * t525;
t391 = m(5) * t417 - mrSges(5,2) * t477 + mrSges(5,3) * t463 + t465 * t487 - t467 * t499 + t540;
t388 = -t390 * t521 + t391 * t523;
t482 = mrSges(4,1) * t499 + mrSges(4,2) * t500;
t490 = mrSges(4,1) * t515 - mrSges(4,3) * t500;
t386 = m(4) * t430 - mrSges(4,2) * t503 - mrSges(4,3) * t477 - t482 * t499 - t490 * t515 + t388;
t429 = -t457 * t526 + t458 * t565;
t420 = -t503 * pkin(3) - t514 * qJ(4) + t481 * t500 + qJDD(4) - t429;
t418 = -t463 * pkin(4) - t486 * pkin(10) + t468 * t488 + t420;
t411 = -0.2e1 * qJD(6) * t461 + (t460 * t498 - t428) * qJ(6) + (t461 * t498 + t427) * pkin(5) + t418;
t405 = m(7) * t411 + mrSges(7,1) * t427 - mrSges(7,3) * t428 + t444 * t460 - t447 * t461;
t533 = m(6) * t418 + mrSges(6,1) * t427 + mrSges(6,2) * t428 + t445 * t460 + t446 * t461 + t405;
t402 = m(5) * t420 - mrSges(5,1) * t463 + mrSges(5,2) * t464 + t467 * t488 - t487 * t537 + t533;
t489 = -mrSges(4,2) * t515 - mrSges(4,3) * t499;
t401 = m(4) * t429 + mrSges(4,1) * t503 - mrSges(4,3) * t478 - t482 * t500 + t489 * t515 - t402;
t382 = t386 * t526 + t401 * t565;
t552 = t460 * t558 + t461 * t559 + t498 * t567;
t541 = t386 * t565 - t401 * t526;
t387 = t390 * t523 + t391 * t521;
t535 = -m(4) * t456 - mrSges(4,1) * t477 - mrSges(4,2) * t478 - t489 * t499 - t490 * t500 - t387;
t393 = -mrSges(6,1) * t418 - mrSges(7,1) * t411 + mrSges(7,2) * t407 + mrSges(6,3) * t410 - pkin(5) * t405 - t427 * t568 + t428 * t560 + t461 * t552 + t475 * t558 - t498 * t551;
t394 = mrSges(6,2) * t418 + mrSges(7,2) * t408 - mrSges(6,3) * t409 - mrSges(7,3) * t411 - qJ(6) * t405 - t427 * t560 + t428 * t569 + t460 * t552 - t475 * t559 + t498 * t553;
t448 = Ifges(5,5) * t488 + Ifges(5,6) * t487 + Ifges(5,3) * t499;
t450 = Ifges(5,1) * t488 + Ifges(5,4) * t487 + Ifges(5,5) * t499;
t379 = -mrSges(5,1) * t420 + mrSges(5,3) * t417 + Ifges(5,4) * t464 + Ifges(5,2) * t463 + Ifges(5,6) * t477 - pkin(4) * t533 + pkin(10) * t540 + t393 * t564 + t525 * t394 - t488 * t448 + t499 * t450;
t449 = Ifges(5,4) * t488 + Ifges(5,2) * t487 + Ifges(5,6) * t499;
t380 = mrSges(5,2) * t420 - mrSges(5,3) * t416 + Ifges(5,1) * t464 + Ifges(5,4) * t463 + Ifges(5,5) * t477 - pkin(10) * t392 - t393 * t525 + t394 * t564 + t448 * t487 - t449 * t499;
t472 = Ifges(4,4) * t500 - Ifges(4,2) * t499 + Ifges(4,6) * t515;
t473 = Ifges(4,1) * t500 - Ifges(4,4) * t499 + Ifges(4,5) * t515;
t532 = mrSges(4,1) * t429 - mrSges(4,2) * t430 + Ifges(4,5) * t478 - Ifges(4,6) * t477 + Ifges(4,3) * t503 - pkin(3) * t402 + qJ(4) * t388 + t523 * t379 + t521 * t380 + t500 * t472 + t499 * t473;
t508 = (-mrSges(3,1) * t529 + mrSges(3,2) * t527) * t548;
t505 = -mrSges(3,2) * t518 + mrSges(3,3) * t543;
t504 = mrSges(3,1) * t518 - mrSges(3,3) * t544;
t494 = -t506 * t522 - t562;
t493 = Ifges(3,5) * t518 + (Ifges(3,1) * t527 + Ifges(3,4) * t529) * t548;
t492 = Ifges(3,6) * t518 + (Ifges(3,4) * t527 + Ifges(3,2) * t529) * t548;
t491 = Ifges(3,3) * t518 + (Ifges(3,5) * t527 + Ifges(3,6) * t529) * t548;
t480 = -g(3) * t555 + t549;
t471 = Ifges(4,5) * t500 - Ifges(4,6) * t499 + Ifges(4,3) * t515;
t383 = m(3) * t479 + mrSges(3,1) * t517 - mrSges(3,3) * t510 + t505 * t518 - t508 * t544 + t535;
t381 = m(3) * t480 - mrSges(3,2) * t517 - mrSges(3,3) * t511 - t504 * t518 + t508 * t543 + t541;
t378 = t566 + t515 * t473 + Ifges(4,6) * t503 - t500 * t471 + t487 * t450 - t488 * t449 + Ifges(4,4) * t478 - Ifges(5,5) * t464 - Ifges(5,6) * t463 - mrSges(4,1) * t456 + mrSges(4,3) * t430 - mrSges(5,1) * t416 + mrSges(5,2) * t417 - pkin(4) * t392 + (-Ifges(5,3) - Ifges(4,2)) * t477 - pkin(3) * t387;
t377 = mrSges(4,2) * t456 - mrSges(4,3) * t429 + Ifges(4,1) * t478 - Ifges(4,4) * t477 + Ifges(4,5) * t503 - qJ(4) * t387 - t379 * t521 + t380 * t523 - t471 * t499 - t472 * t515;
t376 = Ifges(3,5) * t510 - Ifges(3,6) * t511 + Ifges(3,3) * t517 + mrSges(3,1) * t479 - mrSges(3,2) * t480 + t526 * t377 + t565 * t378 + pkin(2) * t535 + pkin(9) * t541 + (t492 * t527 - t493 * t529) * t548;
t1 = [Ifges(2,3) * qJDD(1) + mrSges(2,1) * t542 - mrSges(2,2) * t538 + (mrSges(3,2) * t494 - mrSges(3,3) * t479 + Ifges(3,1) * t510 - Ifges(3,4) * t511 + Ifges(3,5) * t517 - pkin(9) * t382 + t377 * t565 - t378 * t526 + t491 * t543 - t492 * t518) * t555 + (-mrSges(3,1) * t494 + mrSges(3,3) * t480 + Ifges(3,4) * t510 - Ifges(3,2) * t511 + Ifges(3,6) * t517 - pkin(2) * t382 - t491 * t544 + t518 * t493 - t532) * t554 + t524 * t376 + pkin(1) * ((t381 * t527 + t383 * t529) * t524 + (-m(3) * t494 - t511 * mrSges(3,1) - t510 * mrSges(3,2) + (-t504 * t527 + t505 * t529) * t548 - t382) * t522) + (t381 * t529 - t383 * t527) * t563; t376; t532; t402; -t566; t404;];
tauJ  = t1;
