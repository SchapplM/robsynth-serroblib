% Calculate vector of inverse dynamics joint torques for with Newton-Euler
% S6RRRPRR2
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
%   pkin=[a2,a3,a4,a5,a6,d1,d2,d3,d5,d6,theta4]';
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
% Datum: 2019-05-07 10:01
% Revision: 89c353f7eff3bd693eda4e29f35b2761dbc3ada0 (2019-05-03)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauJ = S6RRRPRR2_invdynJ_fixb_snew_vp2(qJ, qJD, qJDD, g, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(6,1),zeros(3,1),zeros(11,1),zeros(7,1),zeros(7,3),zeros(7,6)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRRPRR2_invdynJ_fixb_snew_vp2: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RRRPRR2_invdynJ_fixb_snew_vp2: qJD has to be [6x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [6 1]), ...
  'S6RRRPRR2_invdynJ_fixb_snew_vp2: qJDD has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6RRRPRR2_invdynJ_fixb_snew_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [11 1]), ...
  'S6RRRPRR2_invdynJ_fixb_snew_vp2: pkin has to be [11x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6RRRPRR2_invdynJ_fixb_snew_vp2: m has to be [7x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [7,3]), ...
  'S6RRRPRR2_invdynJ_fixb_snew_vp2: mrSges has to be [7x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [7 6]), ...
  'S6RRRPRR2_invdynJ_fixb_snew_vp2: Ifges has to be [7x6] (double)'); 

%% Symbolic Calculation
% From invdyn_fixb_NewtonEuler_linkframe_tauJ_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-05-07 09:57:37
% EndTime: 2019-05-07 09:57:46
% DurationCPUTime: 9.34s
% Computational Cost: add. (137012->351), mult. (307930->446), div. (0->0), fcn. (231591->12), ass. (0->138)
t542 = sin(qJ(2));
t547 = cos(qJ(2));
t565 = qJD(1) * qJD(2);
t525 = qJDD(1) * t542 + t547 * t565;
t549 = qJD(1) ^ 2;
t543 = sin(qJ(1));
t548 = cos(qJ(1));
t558 = -g(1) * t548 - g(2) * t543;
t522 = -pkin(1) * t549 + qJDD(1) * pkin(7) + t558;
t569 = t542 * t522;
t570 = pkin(2) * t549;
t484 = qJDD(2) * pkin(2) - t525 * pkin(8) - t569 + (pkin(8) * t565 + t542 * t570 - g(3)) * t547;
t509 = -g(3) * t542 + t547 * t522;
t526 = qJDD(1) * t547 - t542 * t565;
t568 = qJD(1) * t542;
t529 = qJD(2) * pkin(2) - pkin(8) * t568;
t536 = t547 ^ 2;
t485 = pkin(8) * t526 - qJD(2) * t529 - t536 * t570 + t509;
t541 = sin(qJ(3));
t546 = cos(qJ(3));
t458 = t546 * t484 - t541 * t485;
t519 = (-t541 * t542 + t546 * t547) * qJD(1);
t493 = qJD(3) * t519 + t525 * t546 + t526 * t541;
t520 = (t541 * t547 + t542 * t546) * qJD(1);
t534 = qJDD(2) + qJDD(3);
t535 = qJD(2) + qJD(3);
t439 = (t519 * t535 - t493) * qJ(4) + (t519 * t520 + t534) * pkin(3) + t458;
t459 = t541 * t484 + t546 * t485;
t492 = -qJD(3) * t520 - t525 * t541 + t526 * t546;
t511 = pkin(3) * t535 - qJ(4) * t520;
t515 = t519 ^ 2;
t443 = -pkin(3) * t515 + qJ(4) * t492 - t511 * t535 + t459;
t537 = sin(pkin(11));
t538 = cos(pkin(11));
t506 = t519 * t537 + t520 * t538;
t423 = -0.2e1 * qJD(4) * t506 + t439 * t538 - t537 * t443;
t505 = t519 * t538 - t537 * t520;
t424 = 0.2e1 * qJD(4) * t505 + t537 * t439 + t538 * t443;
t468 = t492 * t538 - t537 * t493;
t478 = -mrSges(5,1) * t505 + mrSges(5,2) * t506;
t496 = mrSges(5,1) * t535 - mrSges(5,3) * t506;
t479 = -pkin(4) * t505 - pkin(9) * t506;
t533 = t535 ^ 2;
t421 = -pkin(4) * t533 + pkin(9) * t534 + t479 * t505 + t424;
t564 = t543 * g(1) - t548 * g(2);
t557 = -qJDD(1) * pkin(1) - t564;
t494 = -t526 * pkin(2) + t529 * t568 + (-pkin(8) * t536 - pkin(7)) * t549 + t557;
t448 = -t492 * pkin(3) - t515 * qJ(4) + t520 * t511 + qJDD(4) + t494;
t469 = t492 * t537 + t493 * t538;
t427 = (-t505 * t535 - t469) * pkin(9) + (t506 * t535 - t468) * pkin(4) + t448;
t540 = sin(qJ(5));
t545 = cos(qJ(5));
t416 = -t540 * t421 + t545 * t427;
t490 = -t506 * t540 + t535 * t545;
t446 = qJD(5) * t490 + t469 * t545 + t534 * t540;
t467 = qJDD(5) - t468;
t491 = t506 * t545 + t535 * t540;
t500 = qJD(5) - t505;
t414 = (t490 * t500 - t446) * pkin(10) + (t490 * t491 + t467) * pkin(5) + t416;
t417 = t545 * t421 + t540 * t427;
t445 = -qJD(5) * t491 - t469 * t540 + t534 * t545;
t473 = pkin(5) * t500 - pkin(10) * t491;
t487 = t490 ^ 2;
t415 = -pkin(5) * t487 + pkin(10) * t445 - t473 * t500 + t417;
t539 = sin(qJ(6));
t544 = cos(qJ(6));
t412 = t414 * t544 - t415 * t539;
t462 = t490 * t544 - t491 * t539;
t432 = qJD(6) * t462 + t445 * t539 + t446 * t544;
t463 = t490 * t539 + t491 * t544;
t440 = -mrSges(7,1) * t462 + mrSges(7,2) * t463;
t497 = qJD(6) + t500;
t449 = -mrSges(7,2) * t497 + mrSges(7,3) * t462;
t461 = qJDD(6) + t467;
t408 = m(7) * t412 + mrSges(7,1) * t461 - mrSges(7,3) * t432 - t440 * t463 + t449 * t497;
t413 = t414 * t539 + t415 * t544;
t431 = -qJD(6) * t463 + t445 * t544 - t446 * t539;
t450 = mrSges(7,1) * t497 - mrSges(7,3) * t463;
t409 = m(7) * t413 - mrSges(7,2) * t461 + mrSges(7,3) * t431 + t440 * t462 - t450 * t497;
t400 = t544 * t408 + t539 * t409;
t470 = -mrSges(6,1) * t490 + mrSges(6,2) * t491;
t471 = -mrSges(6,2) * t500 + mrSges(6,3) * t490;
t398 = m(6) * t416 + mrSges(6,1) * t467 - mrSges(6,3) * t446 - t470 * t491 + t471 * t500 + t400;
t472 = mrSges(6,1) * t500 - mrSges(6,3) * t491;
t560 = -t408 * t539 + t544 * t409;
t399 = m(6) * t417 - mrSges(6,2) * t467 + mrSges(6,3) * t445 + t470 * t490 - t472 * t500 + t560;
t561 = -t398 * t540 + t545 * t399;
t391 = m(5) * t424 - mrSges(5,2) * t534 + mrSges(5,3) * t468 + t478 * t505 - t496 * t535 + t561;
t495 = -mrSges(5,2) * t535 + mrSges(5,3) * t505;
t420 = -pkin(4) * t534 - pkin(9) * t533 + t506 * t479 - t423;
t418 = -pkin(5) * t445 - pkin(10) * t487 + t473 * t491 + t420;
t556 = m(7) * t418 - t431 * mrSges(7,1) + mrSges(7,2) * t432 - t462 * t449 + t450 * t463;
t553 = -m(6) * t420 + t445 * mrSges(6,1) - mrSges(6,2) * t446 + t490 * t471 - t472 * t491 - t556;
t404 = m(5) * t423 + mrSges(5,1) * t534 - mrSges(5,3) * t469 - t478 * t506 + t495 * t535 + t553;
t384 = t537 * t391 + t538 * t404;
t507 = -mrSges(4,1) * t519 + mrSges(4,2) * t520;
t510 = -mrSges(4,2) * t535 + mrSges(4,3) * t519;
t381 = m(4) * t458 + mrSges(4,1) * t534 - mrSges(4,3) * t493 - t507 * t520 + t510 * t535 + t384;
t512 = mrSges(4,1) * t535 - mrSges(4,3) * t520;
t562 = t538 * t391 - t404 * t537;
t382 = m(4) * t459 - mrSges(4,2) * t534 + mrSges(4,3) * t492 + t507 * t519 - t512 * t535 + t562;
t376 = t546 * t381 + t541 * t382;
t394 = t545 * t398 + t540 * t399;
t567 = qJD(1) * t547;
t563 = -t381 * t541 + t546 * t382;
t555 = -m(5) * t448 + mrSges(5,1) * t468 - t469 * mrSges(5,2) + t495 * t505 - t506 * t496 - t394;
t435 = Ifges(7,4) * t463 + Ifges(7,2) * t462 + Ifges(7,6) * t497;
t436 = Ifges(7,1) * t463 + Ifges(7,4) * t462 + Ifges(7,5) * t497;
t554 = -mrSges(7,1) * t412 + mrSges(7,2) * t413 - Ifges(7,5) * t432 - Ifges(7,6) * t431 - Ifges(7,3) * t461 - t463 * t435 + t462 * t436;
t552 = m(4) * t494 - mrSges(4,1) * t492 + mrSges(4,2) * t493 - t510 * t519 + t512 * t520 - t555;
t434 = Ifges(7,5) * t463 + Ifges(7,6) * t462 + Ifges(7,3) * t497;
t401 = -mrSges(7,1) * t418 + mrSges(7,3) * t413 + Ifges(7,4) * t432 + Ifges(7,2) * t431 + Ifges(7,6) * t461 - t434 * t463 + t436 * t497;
t402 = mrSges(7,2) * t418 - mrSges(7,3) * t412 + Ifges(7,1) * t432 + Ifges(7,4) * t431 + Ifges(7,5) * t461 + t434 * t462 - t435 * t497;
t451 = Ifges(6,5) * t491 + Ifges(6,6) * t490 + Ifges(6,3) * t500;
t453 = Ifges(6,1) * t491 + Ifges(6,4) * t490 + Ifges(6,5) * t500;
t386 = -mrSges(6,1) * t420 + mrSges(6,3) * t417 + Ifges(6,4) * t446 + Ifges(6,2) * t445 + Ifges(6,6) * t467 - pkin(5) * t556 + pkin(10) * t560 + t544 * t401 + t539 * t402 - t491 * t451 + t500 * t453;
t452 = Ifges(6,4) * t491 + Ifges(6,2) * t490 + Ifges(6,6) * t500;
t388 = mrSges(6,2) * t420 - mrSges(6,3) * t416 + Ifges(6,1) * t446 + Ifges(6,4) * t445 + Ifges(6,5) * t467 - pkin(10) * t400 - t401 * t539 + t402 * t544 + t451 * t490 - t452 * t500;
t475 = Ifges(5,4) * t506 + Ifges(5,2) * t505 + Ifges(5,6) * t535;
t476 = Ifges(5,1) * t506 + Ifges(5,4) * t505 + Ifges(5,5) * t535;
t502 = Ifges(4,4) * t520 + Ifges(4,2) * t519 + Ifges(4,6) * t535;
t503 = Ifges(4,1) * t520 + Ifges(4,4) * t519 + Ifges(4,5) * t535;
t551 = mrSges(4,1) * t458 + mrSges(5,1) * t423 - mrSges(4,2) * t459 - mrSges(5,2) * t424 + pkin(3) * t384 + pkin(4) * t553 + pkin(9) * t561 + t545 * t386 + t540 * t388 + t506 * t475 - t505 * t476 - t519 * t503 + Ifges(5,6) * t468 + Ifges(5,5) * t469 + t520 * t502 + Ifges(4,6) * t492 + Ifges(4,5) * t493 + (Ifges(5,3) + Ifges(4,3)) * t534;
t550 = mrSges(6,1) * t416 - mrSges(6,2) * t417 + Ifges(6,5) * t446 + Ifges(6,6) * t445 + Ifges(6,3) * t467 + pkin(5) * t400 + t491 * t452 - t490 * t453 - t554;
t528 = -qJD(2) * mrSges(3,2) + mrSges(3,3) * t567;
t527 = qJD(2) * mrSges(3,1) - mrSges(3,3) * t568;
t524 = (-mrSges(3,1) * t547 + mrSges(3,2) * t542) * qJD(1);
t521 = -t549 * pkin(7) + t557;
t518 = Ifges(3,5) * qJD(2) + (Ifges(3,1) * t542 + Ifges(3,4) * t547) * qJD(1);
t517 = Ifges(3,6) * qJD(2) + (Ifges(3,4) * t542 + Ifges(3,2) * t547) * qJD(1);
t508 = -t547 * g(3) - t569;
t501 = Ifges(4,5) * t520 + Ifges(4,6) * t519 + Ifges(4,3) * t535;
t474 = Ifges(5,5) * t506 + Ifges(5,6) * t505 + Ifges(5,3) * t535;
t377 = -mrSges(5,1) * t448 + mrSges(5,3) * t424 + Ifges(5,4) * t469 + Ifges(5,2) * t468 + Ifges(5,6) * t534 - pkin(4) * t394 - t506 * t474 + t535 * t476 - t550;
t375 = mrSges(5,2) * t448 - mrSges(5,3) * t423 + Ifges(5,1) * t469 + Ifges(5,4) * t468 + Ifges(5,5) * t534 - pkin(9) * t394 - t386 * t540 + t388 * t545 + t474 * t505 - t475 * t535;
t374 = mrSges(4,2) * t494 - mrSges(4,3) * t458 + Ifges(4,1) * t493 + Ifges(4,4) * t492 + Ifges(4,5) * t534 - qJ(4) * t384 + t375 * t538 - t377 * t537 + t501 * t519 - t502 * t535;
t373 = -mrSges(4,1) * t494 + mrSges(4,3) * t459 + Ifges(4,4) * t493 + Ifges(4,2) * t492 + Ifges(4,6) * t534 + pkin(3) * t555 + qJ(4) * t562 + t537 * t375 + t538 * t377 - t520 * t501 + t535 * t503;
t1 = [Ifges(2,3) * qJDD(1) + mrSges(2,1) * t564 - mrSges(2,2) * t558 + t542 * (mrSges(3,2) * t521 - mrSges(3,3) * t508 + Ifges(3,1) * t525 + Ifges(3,4) * t526 + Ifges(3,5) * qJDD(2) - pkin(8) * t376 - qJD(2) * t517 - t541 * t373 + t546 * t374) + t547 * (-mrSges(3,1) * t521 + mrSges(3,3) * t509 + Ifges(3,4) * t525 + Ifges(3,2) * t526 + Ifges(3,6) * qJDD(2) - pkin(2) * t552 + pkin(8) * t563 + qJD(2) * t518 + t546 * t373 + t541 * t374) + pkin(1) * (-t552 - m(3) * t521 + mrSges(3,1) * t526 - mrSges(3,2) * t525 + (-t527 * t542 + t528 * t547) * qJD(1)) + pkin(7) * (t547 * (m(3) * t509 - qJDD(2) * mrSges(3,2) + mrSges(3,3) * t526 - qJD(2) * t527 + t524 * t567 + t563) - t542 * (m(3) * t508 + qJDD(2) * mrSges(3,1) - mrSges(3,3) * t525 + qJD(2) * t528 - t524 * t568 + t376)); t551 + pkin(2) * t376 + Ifges(3,5) * t525 + Ifges(3,6) * t526 + mrSges(3,1) * t508 - mrSges(3,2) * t509 + Ifges(3,3) * qJDD(2) + (t542 * t517 - t547 * t518) * qJD(1); t551; -t555; t550; -t554;];
tauJ  = t1;
