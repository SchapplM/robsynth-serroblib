% Calculate vector of inverse dynamics joint torques for with Newton-Euler
% S6RRPRPR8
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
% pkin [10x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d2,d4,d6,theta3]';
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
% Datum: 2019-05-06 15:01
% Revision: 89c353f7eff3bd693eda4e29f35b2761dbc3ada0 (2019-05-03)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauJ = S6RRPRPR8_invdynJ_fixb_snew_vp2(qJ, qJD, qJDD, g, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(6,1),zeros(3,1),zeros(10,1),zeros(7,1),zeros(7,3),zeros(7,6)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRPRPR8_invdynJ_fixb_snew_vp2: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RRPRPR8_invdynJ_fixb_snew_vp2: qJD has to be [6x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [6 1]), ...
  'S6RRPRPR8_invdynJ_fixb_snew_vp2: qJDD has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6RRPRPR8_invdynJ_fixb_snew_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S6RRPRPR8_invdynJ_fixb_snew_vp2: pkin has to be [10x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6RRPRPR8_invdynJ_fixb_snew_vp2: m has to be [7x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [7,3]), ...
  'S6RRPRPR8_invdynJ_fixb_snew_vp2: mrSges has to be [7x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [7 6]), ...
  'S6RRPRPR8_invdynJ_fixb_snew_vp2: Ifges has to be [7x6] (double)'); 

%% Symbolic Calculation
% From invdyn_fixb_NewtonEuler_linkframe_tauJ_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-05-06 14:55:22
% EndTime: 2019-05-06 14:55:30
% DurationCPUTime: 5.13s
% Computational Cost: add. (46691->328), mult. (100935->398), div. (0->0), fcn. (71297->10), ass. (0->130)
t567 = Ifges(5,1) + Ifges(6,1);
t557 = Ifges(5,4) - Ifges(6,5);
t556 = Ifges(5,5) + Ifges(6,4);
t566 = -Ifges(5,2) - Ifges(6,3);
t555 = Ifges(5,6) - Ifges(6,6);
t564 = Ifges(5,3) + Ifges(6,2);
t521 = sin(pkin(10));
t522 = cos(pkin(10));
t525 = sin(qJ(2));
t549 = qJD(1) * t525;
t505 = qJD(2) * t522 - t521 * t549;
t506 = qJD(2) * t521 + t522 * t549;
t524 = sin(qJ(4));
t559 = cos(qJ(4));
t476 = -t505 * t559 + t506 * t524;
t528 = cos(qJ(2));
t547 = qJD(1) * qJD(2);
t546 = t528 * t547;
t510 = qJDD(1) * t525 + t546;
t485 = qJDD(2) * t522 - t510 * t521;
t486 = qJDD(2) * t521 + t510 * t522;
t442 = -t476 * qJD(4) + t524 * t485 + t486 * t559;
t477 = t524 * t505 + t506 * t559;
t455 = mrSges(6,1) * t476 - mrSges(6,3) * t477;
t531 = qJD(1) ^ 2;
t526 = sin(qJ(1));
t529 = cos(qJ(1));
t545 = g(1) * t526 - t529 * g(2);
t498 = -qJDD(1) * pkin(1) - pkin(7) * t531 - t545;
t518 = t525 * t547;
t511 = qJDD(1) * t528 - t518;
t459 = (-t510 - t546) * qJ(3) + (-t511 + t518) * pkin(2) + t498;
t540 = -g(1) * t529 - g(2) * t526;
t499 = -pkin(1) * t531 + qJDD(1) * pkin(7) + t540;
t481 = -g(3) * t525 + t528 * t499;
t508 = (-pkin(2) * t528 - qJ(3) * t525) * qJD(1);
t530 = qJD(2) ^ 2;
t548 = qJD(1) * t528;
t464 = -pkin(2) * t530 + qJDD(2) * qJ(3) + t508 * t548 + t481;
t431 = -0.2e1 * qJD(3) * t506 + t522 * t459 - t464 * t521;
t421 = (-t505 * t548 - t486) * pkin(8) + (t505 * t506 - t511) * pkin(3) + t431;
t432 = 0.2e1 * qJD(3) * t505 + t521 * t459 + t522 * t464;
t487 = -pkin(3) * t548 - pkin(8) * t506;
t504 = t505 ^ 2;
t424 = -pkin(3) * t504 + pkin(8) * t485 + t487 * t548 + t432;
t411 = t421 * t559 - t524 * t424;
t454 = pkin(4) * t476 - qJ(5) * t477;
t507 = qJDD(4) - t511;
t517 = qJD(4) - t548;
t516 = t517 ^ 2;
t408 = -t507 * pkin(4) - t516 * qJ(5) + t477 * t454 + qJDD(5) - t411;
t554 = t476 * t517;
t402 = (-t442 - t554) * pkin(9) + (t476 * t477 - t507) * pkin(5) + t408;
t412 = t524 * t421 + t559 * t424;
t560 = 2 * qJD(5);
t407 = -pkin(4) * t516 + t507 * qJ(5) - t476 * t454 + t517 * t560 + t412;
t441 = qJD(4) * t477 - t485 * t559 + t486 * t524;
t469 = -pkin(5) * t517 - pkin(9) * t477;
t475 = t476 ^ 2;
t403 = -pkin(5) * t475 + pkin(9) * t441 + t469 * t517 + t407;
t523 = sin(qJ(6));
t527 = cos(qJ(6));
t400 = t402 * t527 - t403 * t523;
t452 = t476 * t527 - t477 * t523;
t418 = qJD(6) * t452 + t441 * t523 + t442 * t527;
t453 = t476 * t523 + t477 * t527;
t430 = -mrSges(7,1) * t452 + mrSges(7,2) * t453;
t515 = qJD(6) - t517;
t436 = -mrSges(7,2) * t515 + mrSges(7,3) * t452;
t503 = qJDD(6) - t507;
t397 = m(7) * t400 + mrSges(7,1) * t503 - mrSges(7,3) * t418 - t430 * t453 + t436 * t515;
t401 = t402 * t523 + t403 * t527;
t417 = -qJD(6) * t453 + t441 * t527 - t442 * t523;
t437 = mrSges(7,1) * t515 - mrSges(7,3) * t453;
t398 = m(7) * t401 - mrSges(7,2) * t503 + mrSges(7,3) * t417 + t430 * t452 - t437 * t515;
t391 = t397 * t527 + t398 * t523;
t468 = -mrSges(6,2) * t476 + mrSges(6,3) * t517;
t537 = -m(6) * t408 + t507 * mrSges(6,1) + t517 * t468 - t391;
t390 = mrSges(6,2) * t442 + t455 * t477 - t537;
t426 = Ifges(7,4) * t453 + Ifges(7,2) * t452 + Ifges(7,6) * t515;
t427 = Ifges(7,1) * t453 + Ifges(7,4) * t452 + Ifges(7,5) * t515;
t535 = -mrSges(7,1) * t400 + mrSges(7,2) * t401 - Ifges(7,5) * t418 - Ifges(7,6) * t417 - Ifges(7,3) * t503 - t453 * t426 + t452 * t427;
t467 = -mrSges(6,1) * t517 + mrSges(6,2) * t477;
t542 = -t397 * t523 + t527 * t398;
t538 = m(6) * t407 + t507 * mrSges(6,3) + t517 * t467 + t542;
t551 = -t557 * t476 + t567 * t477 + t556 * t517;
t552 = t566 * t476 + t557 * t477 + t555 * t517;
t565 = t555 * t441 - t556 * t442 - t564 * t507 - mrSges(5,1) * t411 + mrSges(6,1) * t408 + mrSges(5,2) * t412 - mrSges(6,3) * t407 + pkin(4) * t390 + pkin(5) * t391 - qJ(5) * (-mrSges(6,2) * t441 - t455 * t476 + t538) - t535 - t552 * t477 - t551 * t476;
t480 = -t528 * g(3) - t525 * t499;
t536 = qJDD(2) * pkin(2) + qJ(3) * t530 - t508 * t549 - qJDD(3) + t480;
t534 = pkin(3) * t485 + pkin(8) * t504 - t506 * t487 + t536;
t561 = (-t442 + t554) * qJ(5) - t534;
t558 = -mrSges(5,3) - mrSges(6,2);
t466 = mrSges(5,1) * t517 - mrSges(5,3) * t477;
t550 = -mrSges(5,1) * t476 - mrSges(5,2) * t477 - t455;
t386 = m(5) * t412 - mrSges(5,2) * t507 + t441 * t558 - t466 * t517 + t476 * t550 + t538;
t465 = -mrSges(5,2) * t517 - mrSges(5,3) * t476;
t388 = m(5) * t411 + mrSges(5,1) * t507 + t442 * t558 + t465 * t517 + t477 * t550 + t537;
t383 = t524 * t386 + t559 * t388;
t553 = t555 * t476 - t477 * t556 - t564 * t517;
t478 = -mrSges(4,1) * t505 + mrSges(4,2) * t506;
t483 = mrSges(4,2) * t548 + mrSges(4,3) * t505;
t381 = m(4) * t431 - mrSges(4,1) * t511 - mrSges(4,3) * t486 - t478 * t506 - t483 * t548 + t383;
t484 = -mrSges(4,1) * t548 - mrSges(4,3) * t506;
t543 = t559 * t386 - t388 * t524;
t382 = m(4) * t432 + mrSges(4,2) * t511 + mrSges(4,3) * t485 + t478 * t505 + t484 * t548 + t543;
t544 = -t381 * t521 + t522 * t382;
t405 = -pkin(9) * t475 + (-pkin(4) - pkin(5)) * t441 + (-pkin(4) * t517 + t469 + t560) * t477 - t561;
t539 = -m(7) * t405 + t417 * mrSges(7,1) - t418 * mrSges(7,2) + t452 * t436 - t453 * t437;
t377 = t381 * t522 + t382 * t521;
t410 = -0.2e1 * qJD(5) * t477 + (t477 * t517 + t441) * pkin(4) + t561;
t395 = m(6) * t410 + t441 * mrSges(6,1) - t442 * mrSges(6,3) - t477 * t467 + t476 * t468 + t539;
t533 = -m(5) * t534 + t441 * mrSges(5,1) + t442 * mrSges(5,2) + t476 * t465 + t477 * t466 + t395;
t394 = -m(4) * t536 - t485 * mrSges(4,1) + t486 * mrSges(4,2) - t505 * t483 + t506 * t484 + t533;
t513 = -qJD(2) * mrSges(3,2) + mrSges(3,3) * t548;
t512 = qJD(2) * mrSges(3,1) - mrSges(3,3) * t549;
t509 = (-mrSges(3,1) * t528 + mrSges(3,2) * t525) * qJD(1);
t497 = Ifges(3,5) * qJD(2) + (Ifges(3,1) * t525 + Ifges(3,4) * t528) * qJD(1);
t496 = Ifges(3,6) * qJD(2) + (t525 * Ifges(3,4) + Ifges(3,2) * t528) * qJD(1);
t472 = Ifges(4,1) * t506 + Ifges(4,4) * t505 - Ifges(4,5) * t548;
t471 = Ifges(4,4) * t506 + Ifges(4,2) * t505 - Ifges(4,6) * t548;
t470 = Ifges(4,5) * t506 + Ifges(4,6) * t505 - Ifges(4,3) * t548;
t425 = Ifges(7,5) * t453 + Ifges(7,6) * t452 + Ifges(7,3) * t515;
t393 = mrSges(7,2) * t405 - mrSges(7,3) * t400 + Ifges(7,1) * t418 + Ifges(7,4) * t417 + Ifges(7,5) * t503 + t425 * t452 - t426 * t515;
t392 = -mrSges(7,1) * t405 + mrSges(7,3) * t401 + Ifges(7,4) * t418 + Ifges(7,2) * t417 + Ifges(7,6) * t503 - t425 * t453 + t427 * t515;
t379 = -mrSges(5,2) * t534 + mrSges(6,2) * t408 - mrSges(5,3) * t411 - mrSges(6,3) * t410 - pkin(9) * t391 - qJ(5) * t395 - t392 * t523 + t393 * t527 - t557 * t441 + t567 * t442 + t553 * t476 + t556 * t507 - t552 * t517;
t378 = mrSges(5,1) * t534 - mrSges(6,1) * t410 + mrSges(6,2) * t407 + mrSges(5,3) * t412 - pkin(4) * t395 - pkin(5) * t539 - pkin(9) * t542 - t527 * t392 - t523 * t393 + t566 * t441 + t557 * t442 + t553 * t477 + t555 * t507 + t551 * t517;
t376 = -mrSges(4,2) * t536 - mrSges(4,3) * t431 + Ifges(4,1) * t486 + Ifges(4,4) * t485 - Ifges(4,5) * t511 - pkin(8) * t383 - t524 * t378 + t379 * t559 + t505 * t470 + t471 * t548;
t375 = mrSges(4,1) * t536 + mrSges(4,3) * t432 + Ifges(4,4) * t486 + Ifges(4,2) * t485 - Ifges(4,6) * t511 - pkin(3) * t533 + pkin(8) * t543 + t378 * t559 + t524 * t379 - t506 * t470 - t472 * t548;
t1 = [Ifges(2,3) * qJDD(1) + mrSges(2,1) * t545 - mrSges(2,2) * t540 + t525 * (mrSges(3,2) * t498 - mrSges(3,3) * t480 + Ifges(3,1) * t510 + Ifges(3,4) * t511 + Ifges(3,5) * qJDD(2) - qJ(3) * t377 - qJD(2) * t496 - t521 * t375 + t522 * t376) + t528 * (t565 + Ifges(3,6) * qJDD(2) + (Ifges(3,2) + Ifges(4,3)) * t511 - pkin(3) * t383 + Ifges(3,4) * t510 + t505 * t472 - t506 * t471 - Ifges(4,5) * t486 + qJD(2) * t497 - mrSges(3,1) * t498 + mrSges(3,3) * t481 - Ifges(4,6) * t485 - mrSges(4,1) * t431 + mrSges(4,2) * t432 - pkin(2) * t377) + pkin(1) * (-m(3) * t498 + mrSges(3,1) * t511 - mrSges(3,2) * t510 + (-t512 * t525 + t513 * t528) * qJD(1) - t377) + pkin(7) * (t528 * (m(3) * t481 - qJDD(2) * mrSges(3,2) + mrSges(3,3) * t511 - qJD(2) * t512 + t509 * t548 + t544) - t525 * (m(3) * t480 + qJDD(2) * mrSges(3,1) - t510 * mrSges(3,3) + qJD(2) * t513 - t509 * t549 - t394)); Ifges(3,5) * t510 + Ifges(3,6) * t511 + Ifges(3,3) * qJDD(2) + mrSges(3,1) * t480 - mrSges(3,2) * t481 + t521 * t376 + t522 * t375 - pkin(2) * t394 + qJ(3) * t544 + (t525 * t496 - t528 * t497) * qJD(1); t394; -t565; t390; -t535;];
tauJ  = t1;
