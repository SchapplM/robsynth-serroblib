% Calculate vector of inverse dynamics joint torques for with Newton-Euler
% S6RPRRRR5
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
% Datum: 2019-05-06 03:34
% Revision: 89c353f7eff3bd693eda4e29f35b2761dbc3ada0 (2019-05-03)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauJ = S6RPRRRR5_invdynJ_fixb_snew_vp2(qJ, qJD, qJDD, g, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(6,1),zeros(3,1),zeros(11,1),zeros(7,1),zeros(7,3),zeros(7,6)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPRRRR5_invdynJ_fixb_snew_vp2: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RPRRRR5_invdynJ_fixb_snew_vp2: qJD has to be [6x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [6 1]), ...
  'S6RPRRRR5_invdynJ_fixb_snew_vp2: qJDD has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6RPRRRR5_invdynJ_fixb_snew_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [11 1]), ...
  'S6RPRRRR5_invdynJ_fixb_snew_vp2: pkin has to be [11x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6RPRRRR5_invdynJ_fixb_snew_vp2: m has to be [7x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [7,3]), ...
  'S6RPRRRR5_invdynJ_fixb_snew_vp2: mrSges has to be [7x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [7 6]), ...
  'S6RPRRRR5_invdynJ_fixb_snew_vp2: Ifges has to be [7x6] (double)'); 

%% Symbolic Calculation
% From invdyn_fixb_NewtonEuler_linkframe_tauJ_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-05-06 03:30:55
% EndTime: 2019-05-06 03:31:03
% DurationCPUTime: 7.97s
% Computational Cost: add. (114123->320), mult. (272946->405), div. (0->0), fcn. (218282->12), ass. (0->136)
t538 = qJD(1) ^ 2;
t564 = pkin(2) * t538;
t563 = pkin(7) * qJDD(1);
t532 = sin(qJ(1));
t537 = cos(qJ(1));
t551 = -g(1) * t537 - g(2) * t532;
t514 = -pkin(1) * t538 + qJDD(1) * qJ(2) + t551;
t526 = sin(pkin(11));
t527 = cos(pkin(11));
t558 = qJD(1) * qJD(2);
t556 = -t527 * g(3) - 0.2e1 * t526 * t558;
t487 = (t527 * t564 - t514 - t563) * t526 + t556;
t504 = -g(3) * t526 + (t514 + 0.2e1 * t558) * t527;
t524 = t527 ^ 2;
t488 = -t524 * t564 + t527 * t563 + t504;
t531 = sin(qJ(3));
t536 = cos(qJ(3));
t468 = t536 * t487 - t531 * t488;
t548 = t526 * t536 + t527 * t531;
t547 = -t526 * t531 + t527 * t536;
t512 = t547 * qJD(1);
t559 = t512 * qJD(3);
t502 = t548 * qJDD(1) + t559;
t513 = t548 * qJD(1);
t440 = (-t502 + t559) * pkin(8) + (t512 * t513 + qJDD(3)) * pkin(3) + t468;
t469 = t531 * t487 + t536 * t488;
t501 = -t513 * qJD(3) + t547 * qJDD(1);
t507 = qJD(3) * pkin(3) - pkin(8) * t513;
t511 = t512 ^ 2;
t446 = -pkin(3) * t511 + pkin(8) * t501 - qJD(3) * t507 + t469;
t530 = sin(qJ(4));
t535 = cos(qJ(4));
t433 = t530 * t440 + t535 * t446;
t496 = t512 * t530 + t513 * t535;
t463 = -t496 * qJD(4) + t501 * t535 - t530 * t502;
t495 = t512 * t535 - t530 * t513;
t477 = -mrSges(5,1) * t495 + mrSges(5,2) * t496;
t525 = qJD(3) + qJD(4);
t485 = mrSges(5,1) * t525 - mrSges(5,3) * t496;
t522 = qJDD(3) + qJDD(4);
t478 = -pkin(4) * t495 - pkin(9) * t496;
t521 = t525 ^ 2;
t421 = -pkin(4) * t521 + pkin(9) * t522 + t478 * t495 + t433;
t557 = t532 * g(1) - t537 * g(2);
t550 = qJDD(2) - t557;
t561 = -t526 ^ 2 - t524;
t500 = (-pkin(2) * t527 - pkin(1)) * qJDD(1) + (t561 * pkin(7) - qJ(2)) * t538 + t550;
t455 = -t501 * pkin(3) - t511 * pkin(8) + t513 * t507 + t500;
t464 = qJD(4) * t495 + t501 * t530 + t502 * t535;
t429 = (-t495 * t525 - t464) * pkin(9) + (t496 * t525 - t463) * pkin(4) + t455;
t529 = sin(qJ(5));
t534 = cos(qJ(5));
t416 = -t529 * t421 + t534 * t429;
t480 = -t496 * t529 + t525 * t534;
t443 = qJD(5) * t480 + t464 * t534 + t522 * t529;
t462 = qJDD(5) - t463;
t481 = t496 * t534 + t525 * t529;
t491 = qJD(5) - t495;
t414 = (t480 * t491 - t443) * pkin(10) + (t480 * t481 + t462) * pkin(5) + t416;
t417 = t534 * t421 + t529 * t429;
t442 = -qJD(5) * t481 - t464 * t529 + t522 * t534;
t472 = pkin(5) * t491 - pkin(10) * t481;
t479 = t480 ^ 2;
t415 = -pkin(5) * t479 + pkin(10) * t442 - t472 * t491 + t417;
t528 = sin(qJ(6));
t533 = cos(qJ(6));
t412 = t414 * t533 - t415 * t528;
t465 = t480 * t533 - t481 * t528;
t426 = qJD(6) * t465 + t442 * t528 + t443 * t533;
t466 = t480 * t528 + t481 * t533;
t438 = -mrSges(7,1) * t465 + mrSges(7,2) * t466;
t489 = qJD(6) + t491;
t447 = -mrSges(7,2) * t489 + mrSges(7,3) * t465;
t458 = qJDD(6) + t462;
t408 = m(7) * t412 + mrSges(7,1) * t458 - mrSges(7,3) * t426 - t438 * t466 + t447 * t489;
t413 = t414 * t528 + t415 * t533;
t425 = -qJD(6) * t466 + t442 * t533 - t443 * t528;
t448 = mrSges(7,1) * t489 - mrSges(7,3) * t466;
t409 = m(7) * t413 - mrSges(7,2) * t458 + mrSges(7,3) * t425 + t438 * t465 - t448 * t489;
t400 = t533 * t408 + t528 * t409;
t467 = -mrSges(6,1) * t480 + mrSges(6,2) * t481;
t470 = -mrSges(6,2) * t491 + mrSges(6,3) * t480;
t398 = m(6) * t416 + mrSges(6,1) * t462 - mrSges(6,3) * t443 - t467 * t481 + t470 * t491 + t400;
t471 = mrSges(6,1) * t491 - mrSges(6,3) * t481;
t552 = -t408 * t528 + t533 * t409;
t399 = m(6) * t417 - mrSges(6,2) * t462 + mrSges(6,3) * t442 + t467 * t480 - t471 * t491 + t552;
t553 = -t398 * t529 + t534 * t399;
t392 = m(5) * t433 - mrSges(5,2) * t522 + mrSges(5,3) * t463 + t477 * t495 - t485 * t525 + t553;
t432 = t440 * t535 - t530 * t446;
t484 = -mrSges(5,2) * t525 + mrSges(5,3) * t495;
t420 = -pkin(4) * t522 - pkin(9) * t521 + t496 * t478 - t432;
t418 = -pkin(5) * t442 - pkin(10) * t479 + t472 * t481 + t420;
t545 = m(7) * t418 - t425 * mrSges(7,1) + mrSges(7,2) * t426 - t465 * t447 + t448 * t466;
t541 = -m(6) * t420 + t442 * mrSges(6,1) - mrSges(6,2) * t443 + t480 * t470 - t471 * t481 - t545;
t404 = m(5) * t432 + mrSges(5,1) * t522 - mrSges(5,3) * t464 - t477 * t496 + t484 * t525 + t541;
t384 = t530 * t392 + t535 * t404;
t499 = -mrSges(4,1) * t512 + mrSges(4,2) * t513;
t505 = -qJD(3) * mrSges(4,2) + mrSges(4,3) * t512;
t382 = m(4) * t468 + qJDD(3) * mrSges(4,1) - mrSges(4,3) * t502 + qJD(3) * t505 - t499 * t513 + t384;
t506 = qJD(3) * mrSges(4,1) - mrSges(4,3) * t513;
t554 = t535 * t392 - t404 * t530;
t383 = m(4) * t469 - qJDD(3) * mrSges(4,2) + mrSges(4,3) * t501 - qJD(3) * t506 + t499 * t512 + t554;
t562 = t536 * t382 + t531 * t383;
t394 = t534 * t398 + t529 * t399;
t555 = -t531 * t382 + t536 * t383;
t549 = -mrSges(3,1) * t527 + mrSges(3,2) * t526;
t546 = mrSges(3,3) * qJDD(1) + t538 * t549;
t544 = m(5) * t455 - t463 * mrSges(5,1) + t464 * mrSges(5,2) - t495 * t484 + t496 * t485 + t394;
t435 = Ifges(7,4) * t466 + Ifges(7,2) * t465 + Ifges(7,6) * t489;
t436 = Ifges(7,1) * t466 + Ifges(7,4) * t465 + Ifges(7,5) * t489;
t543 = -mrSges(7,1) * t412 + mrSges(7,2) * t413 - Ifges(7,5) * t426 - Ifges(7,6) * t425 - Ifges(7,3) * t458 - t466 * t435 + t465 * t436;
t434 = Ifges(7,5) * t466 + Ifges(7,6) * t465 + Ifges(7,3) * t489;
t401 = -mrSges(7,1) * t418 + mrSges(7,3) * t413 + Ifges(7,4) * t426 + Ifges(7,2) * t425 + Ifges(7,6) * t458 - t434 * t466 + t436 * t489;
t402 = mrSges(7,2) * t418 - mrSges(7,3) * t412 + Ifges(7,1) * t426 + Ifges(7,4) * t425 + Ifges(7,5) * t458 + t434 * t465 - t435 * t489;
t449 = Ifges(6,5) * t481 + Ifges(6,6) * t480 + Ifges(6,3) * t491;
t451 = Ifges(6,1) * t481 + Ifges(6,4) * t480 + Ifges(6,5) * t491;
t386 = -mrSges(6,1) * t420 + mrSges(6,3) * t417 + Ifges(6,4) * t443 + Ifges(6,2) * t442 + Ifges(6,6) * t462 - pkin(5) * t545 + pkin(10) * t552 + t533 * t401 + t528 * t402 - t481 * t449 + t491 * t451;
t450 = Ifges(6,4) * t481 + Ifges(6,2) * t480 + Ifges(6,6) * t491;
t388 = mrSges(6,2) * t420 - mrSges(6,3) * t416 + Ifges(6,1) * t443 + Ifges(6,4) * t442 + Ifges(6,5) * t462 - pkin(10) * t400 - t401 * t528 + t402 * t533 + t449 * t480 - t450 * t491;
t474 = Ifges(5,4) * t496 + Ifges(5,2) * t495 + Ifges(5,6) * t525;
t475 = Ifges(5,1) * t496 + Ifges(5,4) * t495 + Ifges(5,5) * t525;
t542 = mrSges(5,1) * t432 - mrSges(5,2) * t433 + Ifges(5,5) * t464 + Ifges(5,6) * t463 + Ifges(5,3) * t522 + pkin(4) * t541 + pkin(9) * t553 + t534 * t386 + t529 * t388 + t496 * t474 - t475 * t495;
t540 = m(4) * t500 - t501 * mrSges(4,1) + t502 * mrSges(4,2) - t512 * t505 + t513 * t506 + t544;
t539 = mrSges(6,1) * t416 - mrSges(6,2) * t417 + Ifges(6,5) * t443 + Ifges(6,6) * t442 + Ifges(6,3) * t462 + pkin(5) * t400 + t481 * t450 - t480 * t451 - t543;
t510 = -qJDD(1) * pkin(1) - t538 * qJ(2) + t550;
t503 = -t526 * t514 + t556;
t494 = Ifges(4,1) * t513 + Ifges(4,4) * t512 + Ifges(4,5) * qJD(3);
t493 = Ifges(4,4) * t513 + Ifges(4,2) * t512 + Ifges(4,6) * qJD(3);
t492 = Ifges(4,5) * t513 + Ifges(4,6) * t512 + Ifges(4,3) * qJD(3);
t473 = Ifges(5,5) * t496 + Ifges(5,6) * t495 + Ifges(5,3) * t525;
t389 = t561 * t538 * mrSges(3,3) + m(3) * t510 + t549 * qJDD(1) + t540;
t378 = -mrSges(5,1) * t455 + mrSges(5,3) * t433 + Ifges(5,4) * t464 + Ifges(5,2) * t463 + Ifges(5,6) * t522 - pkin(4) * t394 - t496 * t473 + t525 * t475 - t539;
t377 = mrSges(5,2) * t455 - mrSges(5,3) * t432 + Ifges(5,1) * t464 + Ifges(5,4) * t463 + Ifges(5,5) * t522 - pkin(9) * t394 - t386 * t529 + t388 * t534 + t473 * t495 - t474 * t525;
t376 = mrSges(4,2) * t500 - mrSges(4,3) * t468 + Ifges(4,1) * t502 + Ifges(4,4) * t501 + Ifges(4,5) * qJDD(3) - pkin(8) * t384 - qJD(3) * t493 + t377 * t535 - t378 * t530 + t492 * t512;
t375 = -mrSges(4,1) * t500 + mrSges(4,3) * t469 + Ifges(4,4) * t502 + Ifges(4,2) * t501 + Ifges(4,6) * qJDD(3) - pkin(3) * t544 + pkin(8) * t554 + qJD(3) * t494 + t530 * t377 + t535 * t378 - t513 * t492;
t1 = [Ifges(2,3) * qJDD(1) + mrSges(2,1) * t557 - mrSges(2,2) * t551 + t526 * (mrSges(3,2) * t510 - mrSges(3,3) * t503 + t536 * t376 - t531 * t375 - pkin(7) * t562 + (Ifges(3,1) * t526 + Ifges(3,4) * t527) * qJDD(1)) + t527 * (-mrSges(3,1) * t510 + mrSges(3,3) * t504 + t531 * t376 + t536 * t375 - pkin(2) * t540 + pkin(7) * t555 + (Ifges(3,4) * t526 + Ifges(3,2) * t527) * qJDD(1)) - pkin(1) * t389 + qJ(2) * ((m(3) * t504 + t546 * t527 + t555) * t527 + (-m(3) * t503 + t546 * t526 - t562) * t526); t389; mrSges(4,1) * t468 - mrSges(4,2) * t469 + Ifges(4,5) * t502 + Ifges(4,6) * t501 + Ifges(4,3) * qJDD(3) + pkin(3) * t384 + t493 * t513 - t494 * t512 + t542; t542; t539; -t543;];
tauJ  = t1;
