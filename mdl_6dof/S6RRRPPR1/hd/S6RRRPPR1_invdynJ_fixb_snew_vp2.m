% Calculate vector of inverse dynamics joint torques for with Newton-Euler
% S6RRRPPR1
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
%   pkin=[a2,a3,a4,a5,a6,d1,d2,d3,d6,theta4,theta5]';
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
% Datum: 2019-05-07 04:11
% Revision: 89c353f7eff3bd693eda4e29f35b2761dbc3ada0 (2019-05-03)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauJ = S6RRRPPR1_invdynJ_fixb_snew_vp2(qJ, qJD, qJDD, g, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(6,1),zeros(3,1),zeros(11,1),zeros(7,1),zeros(7,3),zeros(7,6)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRRPPR1_invdynJ_fixb_snew_vp2: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RRRPPR1_invdynJ_fixb_snew_vp2: qJD has to be [6x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [6 1]), ...
  'S6RRRPPR1_invdynJ_fixb_snew_vp2: qJDD has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6RRRPPR1_invdynJ_fixb_snew_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [11 1]), ...
  'S6RRRPPR1_invdynJ_fixb_snew_vp2: pkin has to be [11x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6RRRPPR1_invdynJ_fixb_snew_vp2: m has to be [7x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [7,3]), ...
  'S6RRRPPR1_invdynJ_fixb_snew_vp2: mrSges has to be [7x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [7 6]), ...
  'S6RRRPPR1_invdynJ_fixb_snew_vp2: Ifges has to be [7x6] (double)'); 

%% Symbolic Calculation
% From invdyn_fixb_NewtonEuler_linkframe_tauJ_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-05-07 04:07:52
% EndTime: 2019-05-07 04:08:01
% DurationCPUTime: 9.11s
% Computational Cost: add. (129701->348), mult. (295935->446), div. (0->0), fcn. (220881->12), ass. (0->136)
t566 = -2 * qJD(4);
t540 = sin(qJ(2));
t544 = cos(qJ(2));
t559 = qJD(1) * qJD(2);
t523 = qJDD(1) * t540 + t544 * t559;
t546 = qJD(1) ^ 2;
t541 = sin(qJ(1));
t545 = cos(qJ(1));
t553 = -g(1) * t545 - g(2) * t541;
t520 = -pkin(1) * t546 + qJDD(1) * pkin(7) + t553;
t563 = t520 * t540;
t565 = pkin(2) * t546;
t481 = qJDD(2) * pkin(2) - pkin(8) * t523 - t563 + (pkin(8) * t559 + t540 * t565 - g(3)) * t544;
t507 = -g(3) * t540 + t544 * t520;
t524 = qJDD(1) * t544 - t540 * t559;
t562 = qJD(1) * t540;
t527 = qJD(2) * pkin(2) - pkin(8) * t562;
t534 = t544 ^ 2;
t482 = pkin(8) * t524 - qJD(2) * t527 - t534 * t565 + t507;
t539 = sin(qJ(3));
t543 = cos(qJ(3));
t457 = t543 * t481 - t482 * t539;
t517 = (-t539 * t540 + t543 * t544) * qJD(1);
t492 = qJD(3) * t517 + t523 * t543 + t524 * t539;
t518 = (t539 * t544 + t540 * t543) * qJD(1);
t532 = qJDD(2) + qJDD(3);
t533 = qJD(2) + qJD(3);
t438 = (t517 * t533 - t492) * qJ(4) + (t517 * t518 + t532) * pkin(3) + t457;
t458 = t539 * t481 + t543 * t482;
t491 = -qJD(3) * t518 - t523 * t539 + t524 * t543;
t509 = pkin(3) * t533 - qJ(4) * t518;
t513 = t517 ^ 2;
t442 = -pkin(3) * t513 + qJ(4) * t491 - t509 * t533 + t458;
t536 = sin(pkin(10));
t564 = cos(pkin(10));
t504 = t536 * t517 + t518 * t564;
t425 = t438 * t564 - t536 * t442 + t504 * t566;
t503 = -t564 * t517 + t518 * t536;
t426 = t536 * t438 + t564 * t442 + t503 * t566;
t466 = -t564 * t491 + t492 * t536;
t476 = mrSges(5,1) * t503 + mrSges(5,2) * t504;
t495 = mrSges(5,1) * t533 - mrSges(5,3) * t504;
t475 = pkin(4) * t503 - qJ(5) * t504;
t531 = t533 ^ 2;
t423 = -pkin(4) * t531 + qJ(5) * t532 - t475 * t503 + t426;
t558 = g(1) * t541 - t545 * g(2);
t551 = -qJDD(1) * pkin(1) - t558;
t493 = -pkin(2) * t524 + t527 * t562 + (-pkin(8) * t534 - pkin(7)) * t546 + t551;
t444 = -pkin(3) * t491 - qJ(4) * t513 + t518 * t509 + qJDD(4) + t493;
t467 = t536 * t491 + t492 * t564;
t429 = (t503 * t533 - t467) * qJ(5) + (t504 * t533 + t466) * pkin(4) + t444;
t535 = sin(pkin(11));
t537 = cos(pkin(11));
t488 = t504 * t537 + t533 * t535;
t418 = -0.2e1 * qJD(5) * t488 - t423 * t535 + t537 * t429;
t456 = t467 * t537 + t532 * t535;
t487 = -t504 * t535 + t533 * t537;
t416 = (t487 * t503 - t456) * pkin(9) + (t487 * t488 + t466) * pkin(5) + t418;
t419 = 0.2e1 * qJD(5) * t487 + t537 * t423 + t535 * t429;
t455 = -t467 * t535 + t532 * t537;
t470 = pkin(5) * t503 - pkin(9) * t488;
t486 = t487 ^ 2;
t417 = -pkin(5) * t486 + pkin(9) * t455 - t470 * t503 + t419;
t538 = sin(qJ(6));
t542 = cos(qJ(6));
t414 = t416 * t542 - t417 * t538;
t459 = t487 * t542 - t488 * t538;
t432 = qJD(6) * t459 + t455 * t538 + t456 * t542;
t460 = t487 * t538 + t488 * t542;
t439 = -mrSges(7,1) * t459 + mrSges(7,2) * t460;
t498 = qJD(6) + t503;
t445 = -mrSges(7,2) * t498 + mrSges(7,3) * t459;
t464 = qJDD(6) + t466;
t410 = m(7) * t414 + mrSges(7,1) * t464 - mrSges(7,3) * t432 - t439 * t460 + t445 * t498;
t415 = t416 * t538 + t417 * t542;
t431 = -qJD(6) * t460 + t455 * t542 - t456 * t538;
t446 = mrSges(7,1) * t498 - mrSges(7,3) * t460;
t411 = m(7) * t415 - mrSges(7,2) * t464 + mrSges(7,3) * t431 + t439 * t459 - t446 * t498;
t402 = t542 * t410 + t538 * t411;
t465 = -mrSges(6,1) * t487 + mrSges(6,2) * t488;
t468 = -mrSges(6,2) * t503 + mrSges(6,3) * t487;
t400 = m(6) * t418 + mrSges(6,1) * t466 - mrSges(6,3) * t456 - t465 * t488 + t468 * t503 + t402;
t469 = mrSges(6,1) * t503 - mrSges(6,3) * t488;
t554 = -t410 * t538 + t542 * t411;
t401 = m(6) * t419 - mrSges(6,2) * t466 + mrSges(6,3) * t455 + t465 * t487 - t469 * t503 + t554;
t555 = -t400 * t535 + t537 * t401;
t393 = m(5) * t426 - mrSges(5,2) * t532 - mrSges(5,3) * t466 - t476 * t503 - t495 * t533 + t555;
t422 = -t532 * pkin(4) - t531 * qJ(5) + t504 * t475 + qJDD(5) - t425;
t420 = -t455 * pkin(5) - t486 * pkin(9) + t488 * t470 + t422;
t550 = m(7) * t420 - t431 * mrSges(7,1) + mrSges(7,2) * t432 - t459 * t445 + t446 * t460;
t413 = m(6) * t422 - t455 * mrSges(6,1) + mrSges(6,2) * t456 - t487 * t468 + t469 * t488 + t550;
t494 = -mrSges(5,2) * t533 - mrSges(5,3) * t503;
t406 = m(5) * t425 + mrSges(5,1) * t532 - mrSges(5,3) * t467 - t476 * t504 + t494 * t533 - t413;
t386 = t536 * t393 + t564 * t406;
t505 = -mrSges(4,1) * t517 + mrSges(4,2) * t518;
t508 = -mrSges(4,2) * t533 + mrSges(4,3) * t517;
t383 = m(4) * t457 + mrSges(4,1) * t532 - mrSges(4,3) * t492 - t505 * t518 + t508 * t533 + t386;
t510 = mrSges(4,1) * t533 - mrSges(4,3) * t518;
t556 = t564 * t393 - t406 * t536;
t384 = m(4) * t458 - mrSges(4,2) * t532 + mrSges(4,3) * t491 + t505 * t517 - t510 * t533 + t556;
t378 = t543 * t383 + t539 * t384;
t396 = t537 * t400 + t535 * t401;
t561 = qJD(1) * t544;
t557 = -t383 * t539 + t543 * t384;
t394 = m(5) * t444 + mrSges(5,1) * t466 + t467 * mrSges(5,2) + t494 * t503 + t504 * t495 + t396;
t434 = Ifges(7,4) * t460 + Ifges(7,2) * t459 + Ifges(7,6) * t498;
t435 = Ifges(7,1) * t460 + Ifges(7,4) * t459 + Ifges(7,5) * t498;
t549 = mrSges(7,1) * t414 - mrSges(7,2) * t415 + Ifges(7,5) * t432 + Ifges(7,6) * t431 + Ifges(7,3) * t464 + t460 * t434 - t459 * t435;
t548 = m(4) * t493 - mrSges(4,1) * t491 + mrSges(4,2) * t492 - t508 * t517 + t510 * t518 + t394;
t433 = Ifges(7,5) * t460 + Ifges(7,6) * t459 + Ifges(7,3) * t498;
t403 = -mrSges(7,1) * t420 + mrSges(7,3) * t415 + Ifges(7,4) * t432 + Ifges(7,2) * t431 + Ifges(7,6) * t464 - t433 * t460 + t435 * t498;
t404 = mrSges(7,2) * t420 - mrSges(7,3) * t414 + Ifges(7,1) * t432 + Ifges(7,4) * t431 + Ifges(7,5) * t464 + t433 * t459 - t434 * t498;
t447 = Ifges(6,5) * t488 + Ifges(6,6) * t487 + Ifges(6,3) * t503;
t449 = Ifges(6,1) * t488 + Ifges(6,4) * t487 + Ifges(6,5) * t503;
t388 = -mrSges(6,1) * t422 + mrSges(6,3) * t419 + Ifges(6,4) * t456 + Ifges(6,2) * t455 + Ifges(6,6) * t466 - pkin(5) * t550 + pkin(9) * t554 + t542 * t403 + t538 * t404 - t488 * t447 + t503 * t449;
t448 = Ifges(6,4) * t488 + Ifges(6,2) * t487 + Ifges(6,6) * t503;
t390 = mrSges(6,2) * t422 - mrSges(6,3) * t418 + Ifges(6,1) * t456 + Ifges(6,4) * t455 + Ifges(6,5) * t466 - pkin(9) * t402 - t403 * t538 + t404 * t542 + t447 * t487 - t448 * t503;
t472 = Ifges(5,4) * t504 - Ifges(5,2) * t503 + Ifges(5,6) * t533;
t473 = Ifges(5,1) * t504 - Ifges(5,4) * t503 + Ifges(5,5) * t533;
t500 = Ifges(4,4) * t518 + Ifges(4,2) * t517 + Ifges(4,6) * t533;
t501 = Ifges(4,1) * t518 + Ifges(4,4) * t517 + Ifges(4,5) * t533;
t547 = mrSges(4,1) * t457 + mrSges(5,1) * t425 - mrSges(4,2) * t458 - mrSges(5,2) * t426 + pkin(3) * t386 - pkin(4) * t413 + qJ(5) * t555 + t537 * t388 + t535 * t390 + t504 * t472 + t503 * t473 - t517 * t501 - Ifges(5,6) * t466 + Ifges(5,5) * t467 + t518 * t500 + Ifges(4,6) * t491 + Ifges(4,5) * t492 + (Ifges(5,3) + Ifges(4,3)) * t532;
t526 = -qJD(2) * mrSges(3,2) + mrSges(3,3) * t561;
t525 = qJD(2) * mrSges(3,1) - mrSges(3,3) * t562;
t522 = (-mrSges(3,1) * t544 + mrSges(3,2) * t540) * qJD(1);
t519 = -pkin(7) * t546 + t551;
t516 = Ifges(3,5) * qJD(2) + (Ifges(3,1) * t540 + Ifges(3,4) * t544) * qJD(1);
t515 = Ifges(3,6) * qJD(2) + (Ifges(3,4) * t540 + Ifges(3,2) * t544) * qJD(1);
t506 = -g(3) * t544 - t563;
t499 = Ifges(4,5) * t518 + Ifges(4,6) * t517 + Ifges(4,3) * t533;
t471 = Ifges(5,5) * t504 - Ifges(5,6) * t503 + Ifges(5,3) * t533;
t379 = Ifges(5,6) * t532 + t533 * t473 - t504 * t471 + t487 * t449 - t488 * t448 + Ifges(5,4) * t467 - Ifges(6,5) * t456 - Ifges(6,6) * t455 - mrSges(5,1) * t444 + mrSges(5,3) * t426 - mrSges(6,1) * t418 + mrSges(6,2) * t419 + (-Ifges(5,2) - Ifges(6,3)) * t466 - pkin(5) * t402 - pkin(4) * t396 - t549;
t377 = mrSges(5,2) * t444 - mrSges(5,3) * t425 + Ifges(5,1) * t467 - Ifges(5,4) * t466 + Ifges(5,5) * t532 - qJ(5) * t396 - t388 * t535 + t390 * t537 - t471 * t503 - t472 * t533;
t376 = mrSges(4,2) * t493 - mrSges(4,3) * t457 + Ifges(4,1) * t492 + Ifges(4,4) * t491 + Ifges(4,5) * t532 - qJ(4) * t386 + t377 * t564 - t536 * t379 + t517 * t499 - t533 * t500;
t375 = -mrSges(4,1) * t493 + mrSges(4,3) * t458 + Ifges(4,4) * t492 + Ifges(4,2) * t491 + Ifges(4,6) * t532 - pkin(3) * t394 + qJ(4) * t556 + t536 * t377 + t379 * t564 - t518 * t499 + t533 * t501;
t1 = [Ifges(2,3) * qJDD(1) + mrSges(2,1) * t558 - mrSges(2,2) * t553 + t540 * (mrSges(3,2) * t519 - mrSges(3,3) * t506 + Ifges(3,1) * t523 + Ifges(3,4) * t524 + Ifges(3,5) * qJDD(2) - pkin(8) * t378 - qJD(2) * t515 - t539 * t375 + t543 * t376) + t544 * (-mrSges(3,1) * t519 + mrSges(3,3) * t507 + Ifges(3,4) * t523 + Ifges(3,2) * t524 + Ifges(3,6) * qJDD(2) - pkin(2) * t548 + pkin(8) * t557 + qJD(2) * t516 + t543 * t375 + t539 * t376) + pkin(1) * ((-t525 * t540 + t526 * t544) * qJD(1) - t548 - m(3) * t519 + mrSges(3,1) * t524 - mrSges(3,2) * t523) + pkin(7) * (t544 * (m(3) * t507 - qJDD(2) * mrSges(3,2) + mrSges(3,3) * t524 - qJD(2) * t525 + t522 * t561 + t557) - t540 * (m(3) * t506 + qJDD(2) * mrSges(3,1) - mrSges(3,3) * t523 + qJD(2) * t526 - t522 * t562 + t378)); pkin(2) * t378 + Ifges(3,5) * t523 + Ifges(3,6) * t524 + mrSges(3,1) * t506 - mrSges(3,2) * t507 + Ifges(3,3) * qJDD(2) + t547 + (t540 * t515 - t544 * t516) * qJD(1); t547; t394; t413; t549;];
tauJ  = t1;
