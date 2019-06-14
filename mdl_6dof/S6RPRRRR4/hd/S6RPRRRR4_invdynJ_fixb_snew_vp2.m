% Calculate vector of inverse dynamics joint torques for with Newton-Euler
% S6RPRRRR4
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
%   pkin=[a2,a3,a4,a5,a6,d1,d3,d4,d5,d6,theta2]';
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
% Datum: 2019-05-06 03:17
% Revision: 89c353f7eff3bd693eda4e29f35b2761dbc3ada0 (2019-05-03)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauJ = S6RPRRRR4_invdynJ_fixb_snew_vp2(qJ, qJD, qJDD, g, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(6,1),zeros(3,1),zeros(11,1),zeros(7,1),zeros(7,3),zeros(7,6)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPRRRR4_invdynJ_fixb_snew_vp2: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RPRRRR4_invdynJ_fixb_snew_vp2: qJD has to be [6x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [6 1]), ...
  'S6RPRRRR4_invdynJ_fixb_snew_vp2: qJDD has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6RPRRRR4_invdynJ_fixb_snew_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [11 1]), ...
  'S6RPRRRR4_invdynJ_fixb_snew_vp2: pkin has to be [11x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6RPRRRR4_invdynJ_fixb_snew_vp2: m has to be [7x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [7,3]), ...
  'S6RPRRRR4_invdynJ_fixb_snew_vp2: mrSges has to be [7x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [7 6]), ...
  'S6RPRRRR4_invdynJ_fixb_snew_vp2: Ifges has to be [7x6] (double)'); 

%% Symbolic Calculation
% From invdyn_fixb_NewtonEuler_linkframe_tauJ_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-05-06 03:13:42
% EndTime: 2019-05-06 03:13:50
% DurationCPUTime: 8.28s
% Computational Cost: add. (120427->319), mult. (301285->405), div. (0->0), fcn. (244433->12), ass. (0->136)
t533 = qJD(1) ^ 2;
t559 = pkin(2) * t533;
t558 = pkin(7) * qJDD(1);
t527 = sin(qJ(1));
t532 = cos(qJ(1));
t546 = -g(1) * t532 - g(2) * t527;
t506 = -pkin(1) * t533 + qJDD(1) * qJ(2) + t546;
t521 = sin(pkin(11));
t522 = cos(pkin(11));
t553 = qJD(1) * qJD(2);
t551 = -g(3) * t522 - 0.2e1 * t521 * t553;
t483 = (t522 * t559 - t506 - t558) * t521 + t551;
t497 = -g(3) * t521 + (t506 + 0.2e1 * t553) * t522;
t519 = t522 ^ 2;
t484 = -t519 * t559 + t522 * t558 + t497;
t526 = sin(qJ(3));
t531 = cos(qJ(3));
t464 = t531 * t483 - t484 * t526;
t543 = t521 * t531 + t522 * t526;
t542 = -t521 * t526 + t522 * t531;
t504 = t542 * qJD(1);
t554 = qJD(3) * t504;
t495 = qJDD(1) * t543 + t554;
t505 = t543 * qJD(1);
t445 = (-t495 + t554) * pkin(8) + (t504 * t505 + qJDD(3)) * pkin(3) + t464;
t465 = t526 * t483 + t531 * t484;
t494 = -qJD(3) * t505 + qJDD(1) * t542;
t500 = qJD(3) * pkin(3) - pkin(8) * t505;
t503 = t504 ^ 2;
t453 = -pkin(3) * t503 + pkin(8) * t494 - qJD(3) * t500 + t465;
t525 = sin(qJ(4));
t530 = cos(qJ(4));
t427 = t530 * t445 - t453 * t525;
t489 = t504 * t530 - t505 * t525;
t461 = qJD(4) * t489 + t494 * t525 + t495 * t530;
t490 = t504 * t525 + t505 * t530;
t517 = qJDD(3) + qJDD(4);
t520 = qJD(3) + qJD(4);
t418 = (t489 * t520 - t461) * pkin(9) + (t489 * t490 + t517) * pkin(4) + t427;
t428 = t525 * t445 + t530 * t453;
t460 = -qJD(4) * t490 + t494 * t530 - t495 * t525;
t481 = pkin(4) * t520 - pkin(9) * t490;
t485 = t489 ^ 2;
t420 = -pkin(4) * t485 + pkin(9) * t460 - t481 * t520 + t428;
t524 = sin(qJ(5));
t529 = cos(qJ(5));
t416 = t524 * t418 + t529 * t420;
t475 = t489 * t524 + t490 * t529;
t433 = -qJD(5) * t475 + t460 * t529 - t461 * t524;
t474 = t489 * t529 - t490 * t524;
t450 = -mrSges(6,1) * t474 + mrSges(6,2) * t475;
t515 = qJD(5) + t520;
t467 = mrSges(6,1) * t515 - mrSges(6,3) * t475;
t514 = qJDD(5) + t517;
t452 = -pkin(5) * t474 - pkin(10) * t475;
t513 = t515 ^ 2;
t412 = -pkin(5) * t513 + pkin(10) * t514 + t452 * t474 + t416;
t552 = g(1) * t527 - t532 * g(2);
t545 = qJDD(2) - t552;
t556 = -t521 ^ 2 - t519;
t493 = (-pkin(2) * t522 - pkin(1)) * qJDD(1) + (pkin(7) * t556 - qJ(2)) * t533 + t545;
t456 = -pkin(3) * t494 - pkin(8) * t503 + t505 * t500 + t493;
t425 = -pkin(4) * t460 - pkin(9) * t485 + t490 * t481 + t456;
t434 = qJD(5) * t474 + t460 * t524 + t461 * t529;
t413 = t425 + (t475 * t515 - t433) * pkin(5) + (-t474 * t515 - t434) * pkin(10);
t523 = sin(qJ(6));
t528 = cos(qJ(6));
t409 = -t412 * t523 + t413 * t528;
t462 = -t475 * t523 + t515 * t528;
t423 = qJD(6) * t462 + t434 * t528 + t514 * t523;
t432 = qJDD(6) - t433;
t463 = t475 * t528 + t515 * t523;
t440 = -mrSges(7,1) * t462 + mrSges(7,2) * t463;
t471 = qJD(6) - t474;
t443 = -mrSges(7,2) * t471 + mrSges(7,3) * t462;
t406 = m(7) * t409 + mrSges(7,1) * t432 - mrSges(7,3) * t423 - t440 * t463 + t443 * t471;
t410 = t412 * t528 + t413 * t523;
t422 = -qJD(6) * t463 - t434 * t523 + t514 * t528;
t444 = mrSges(7,1) * t471 - mrSges(7,3) * t463;
t407 = m(7) * t410 - mrSges(7,2) * t432 + mrSges(7,3) * t422 + t440 * t462 - t444 * t471;
t547 = -t406 * t523 + t528 * t407;
t394 = m(6) * t416 - mrSges(6,2) * t514 + mrSges(6,3) * t433 + t450 * t474 - t467 * t515 + t547;
t415 = t418 * t529 - t420 * t524;
t466 = -mrSges(6,2) * t515 + mrSges(6,3) * t474;
t411 = -pkin(5) * t514 - pkin(10) * t513 + t452 * t475 - t415;
t540 = -m(7) * t411 + t422 * mrSges(7,1) - mrSges(7,2) * t423 + t462 * t443 - t444 * t463;
t402 = m(6) * t415 + mrSges(6,1) * t514 - mrSges(6,3) * t434 - t450 * t475 + t466 * t515 + t540;
t390 = t524 * t394 + t529 * t402;
t476 = -mrSges(5,1) * t489 + mrSges(5,2) * t490;
t479 = -mrSges(5,2) * t520 + mrSges(5,3) * t489;
t387 = m(5) * t427 + mrSges(5,1) * t517 - mrSges(5,3) * t461 - t476 * t490 + t479 * t520 + t390;
t480 = mrSges(5,1) * t520 - mrSges(5,3) * t490;
t548 = t529 * t394 - t402 * t524;
t388 = m(5) * t428 - mrSges(5,2) * t517 + mrSges(5,3) * t460 + t476 * t489 - t480 * t520 + t548;
t381 = t530 * t387 + t525 * t388;
t492 = -mrSges(4,1) * t504 + mrSges(4,2) * t505;
t498 = -qJD(3) * mrSges(4,2) + mrSges(4,3) * t504;
t379 = m(4) * t464 + qJDD(3) * mrSges(4,1) - mrSges(4,3) * t495 + qJD(3) * t498 - t492 * t505 + t381;
t499 = qJD(3) * mrSges(4,1) - mrSges(4,3) * t505;
t549 = -t387 * t525 + t530 * t388;
t380 = m(4) * t465 - qJDD(3) * mrSges(4,2) + mrSges(4,3) * t494 - qJD(3) * t499 + t492 * t504 + t549;
t557 = t531 * t379 + t526 * t380;
t396 = t528 * t406 + t523 * t407;
t550 = -t526 * t379 + t531 * t380;
t544 = -t522 * mrSges(3,1) + t521 * mrSges(3,2);
t541 = mrSges(3,3) * qJDD(1) + t533 * t544;
t539 = m(6) * t425 - t433 * mrSges(6,1) + t434 * mrSges(6,2) - t474 * t466 + t475 * t467 + t396;
t435 = Ifges(7,5) * t463 + Ifges(7,6) * t462 + Ifges(7,3) * t471;
t437 = Ifges(7,1) * t463 + Ifges(7,4) * t462 + Ifges(7,5) * t471;
t399 = -mrSges(7,1) * t411 + mrSges(7,3) * t410 + Ifges(7,4) * t423 + Ifges(7,2) * t422 + Ifges(7,6) * t432 - t435 * t463 + t437 * t471;
t436 = Ifges(7,4) * t463 + Ifges(7,2) * t462 + Ifges(7,6) * t471;
t400 = mrSges(7,2) * t411 - mrSges(7,3) * t409 + Ifges(7,1) * t423 + Ifges(7,4) * t422 + Ifges(7,5) * t432 + t435 * t462 - t436 * t471;
t447 = Ifges(6,4) * t475 + Ifges(6,2) * t474 + Ifges(6,6) * t515;
t448 = Ifges(6,1) * t475 + Ifges(6,4) * t474 + Ifges(6,5) * t515;
t538 = mrSges(6,1) * t415 - mrSges(6,2) * t416 + Ifges(6,5) * t434 + Ifges(6,6) * t433 + Ifges(6,3) * t514 + pkin(5) * t540 + pkin(10) * t547 + t528 * t399 + t523 * t400 + t475 * t447 - t474 * t448;
t537 = mrSges(7,1) * t409 - mrSges(7,2) * t410 + Ifges(7,5) * t423 + Ifges(7,6) * t422 + Ifges(7,3) * t432 + t436 * t463 - t437 * t462;
t536 = m(5) * t456 - t460 * mrSges(5,1) + t461 * mrSges(5,2) - t489 * t479 + t490 * t480 + t539;
t469 = Ifges(5,4) * t490 + Ifges(5,2) * t489 + Ifges(5,6) * t520;
t470 = Ifges(5,1) * t490 + Ifges(5,4) * t489 + Ifges(5,5) * t520;
t535 = mrSges(5,1) * t427 - mrSges(5,2) * t428 + Ifges(5,5) * t461 + Ifges(5,6) * t460 + Ifges(5,3) * t517 + pkin(4) * t390 + t490 * t469 - t489 * t470 + t538;
t534 = m(4) * t493 - t494 * mrSges(4,1) + t495 * mrSges(4,2) - t504 * t498 + t505 * t499 + t536;
t502 = -qJDD(1) * pkin(1) - qJ(2) * t533 + t545;
t496 = -t506 * t521 + t551;
t488 = Ifges(4,1) * t505 + Ifges(4,4) * t504 + Ifges(4,5) * qJD(3);
t487 = Ifges(4,4) * t505 + Ifges(4,2) * t504 + Ifges(4,6) * qJD(3);
t486 = Ifges(4,5) * t505 + Ifges(4,6) * t504 + Ifges(4,3) * qJD(3);
t468 = Ifges(5,5) * t490 + Ifges(5,6) * t489 + Ifges(5,3) * t520;
t446 = Ifges(6,5) * t475 + Ifges(6,6) * t474 + Ifges(6,3) * t515;
t391 = mrSges(3,3) * t533 * t556 + m(3) * t502 + qJDD(1) * t544 + t534;
t383 = -mrSges(6,1) * t425 + mrSges(6,3) * t416 + Ifges(6,4) * t434 + Ifges(6,2) * t433 + Ifges(6,6) * t514 - pkin(5) * t396 - t446 * t475 + t448 * t515 - t537;
t382 = mrSges(6,2) * t425 - mrSges(6,3) * t415 + Ifges(6,1) * t434 + Ifges(6,4) * t433 + Ifges(6,5) * t514 - pkin(10) * t396 - t399 * t523 + t400 * t528 + t446 * t474 - t447 * t515;
t375 = mrSges(5,2) * t456 - mrSges(5,3) * t427 + Ifges(5,1) * t461 + Ifges(5,4) * t460 + Ifges(5,5) * t517 - pkin(9) * t390 + t382 * t529 - t383 * t524 + t468 * t489 - t469 * t520;
t374 = -mrSges(5,1) * t456 + mrSges(5,3) * t428 + Ifges(5,4) * t461 + Ifges(5,2) * t460 + Ifges(5,6) * t517 - pkin(4) * t539 + pkin(9) * t548 + t524 * t382 + t529 * t383 - t490 * t468 + t520 * t470;
t373 = mrSges(4,2) * t493 - mrSges(4,3) * t464 + Ifges(4,1) * t495 + Ifges(4,4) * t494 + Ifges(4,5) * qJDD(3) - pkin(8) * t381 - qJD(3) * t487 - t374 * t525 + t375 * t530 + t486 * t504;
t372 = -mrSges(4,1) * t493 + mrSges(4,3) * t465 + Ifges(4,4) * t495 + Ifges(4,2) * t494 + Ifges(4,6) * qJDD(3) - pkin(3) * t536 + pkin(8) * t549 + qJD(3) * t488 + t530 * t374 + t525 * t375 - t505 * t486;
t1 = [Ifges(2,3) * qJDD(1) + mrSges(2,1) * t552 - mrSges(2,2) * t546 + t521 * (mrSges(3,2) * t502 - mrSges(3,3) * t496 + t531 * t373 - t526 * t372 - pkin(7) * t557 + (Ifges(3,1) * t521 + Ifges(3,4) * t522) * qJDD(1)) + t522 * (-mrSges(3,1) * t502 + mrSges(3,3) * t497 + t526 * t373 + t531 * t372 - pkin(2) * t534 + pkin(7) * t550 + (Ifges(3,4) * t521 + Ifges(3,2) * t522) * qJDD(1)) - pkin(1) * t391 + qJ(2) * ((m(3) * t497 + t522 * t541 + t550) * t522 + (-m(3) * t496 + t521 * t541 - t557) * t521); t391; mrSges(4,1) * t464 - mrSges(4,2) * t465 + Ifges(4,5) * t495 + Ifges(4,6) * t494 + Ifges(4,3) * qJDD(3) + pkin(3) * t381 + t505 * t487 - t504 * t488 + t535; t535; t538; t537;];
tauJ  = t1;
