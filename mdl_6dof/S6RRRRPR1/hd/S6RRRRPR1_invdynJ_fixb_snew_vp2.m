% Calculate vector of inverse dynamics joint torques for with Newton-Euler
% S6RRRRPR1
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
%   pkin=[a2,a3,a4,a5,a6,d1,d2,d3,d4,d6,theta5]';
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
% Datum: 2019-05-07 19:32
% Revision: 89c353f7eff3bd693eda4e29f35b2761dbc3ada0 (2019-05-03)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauJ = S6RRRRPR1_invdynJ_fixb_snew_vp2(qJ, qJD, qJDD, g, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(6,1),zeros(3,1),zeros(11,1),zeros(7,1),zeros(7,3),zeros(7,6)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRRRPR1_invdynJ_fixb_snew_vp2: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RRRRPR1_invdynJ_fixb_snew_vp2: qJD has to be [6x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [6 1]), ...
  'S6RRRRPR1_invdynJ_fixb_snew_vp2: qJDD has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6RRRRPR1_invdynJ_fixb_snew_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [11 1]), ...
  'S6RRRRPR1_invdynJ_fixb_snew_vp2: pkin has to be [11x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6RRRRPR1_invdynJ_fixb_snew_vp2: m has to be [7x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [7,3]), ...
  'S6RRRRPR1_invdynJ_fixb_snew_vp2: mrSges has to be [7x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [7 6]), ...
  'S6RRRRPR1_invdynJ_fixb_snew_vp2: Ifges has to be [7x6] (double)'); 

%% Symbolic Calculation
% From invdyn_fixb_NewtonEuler_linkframe_tauJ_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-05-07 19:27:39
% EndTime: 2019-05-07 19:27:50
% DurationCPUTime: 10.57s
% Computational Cost: add. (154862->351), mult. (358777->446), div. (0->0), fcn. (274461->12), ass. (0->140)
t564 = 2 * qJD(5);
t543 = qJD(1) ^ 2;
t563 = pkin(2) * t543;
t537 = sin(qJ(1));
t542 = cos(qJ(1));
t553 = -g(1) * t542 - g(2) * t537;
t513 = -pkin(1) * t543 + qJDD(1) * pkin(7) + t553;
t536 = sin(qJ(2));
t562 = t536 * t513;
t541 = cos(qJ(2));
t559 = qJD(1) * qJD(2);
t516 = qJDD(1) * t536 + t541 * t559;
t482 = qJDD(2) * pkin(2) - t516 * pkin(8) - t562 + (pkin(8) * t559 + t536 * t563 - g(3)) * t541;
t501 = -g(3) * t536 + t541 * t513;
t517 = qJDD(1) * t541 - t536 * t559;
t561 = qJD(1) * t536;
t520 = qJD(2) * pkin(2) - pkin(8) * t561;
t530 = t541 ^ 2;
t483 = pkin(8) * t517 - qJD(2) * t520 - t530 * t563 + t501;
t535 = sin(qJ(3));
t540 = cos(qJ(3));
t462 = t540 * t482 - t535 * t483;
t510 = (-t535 * t536 + t540 * t541) * qJD(1);
t487 = qJD(3) * t510 + t516 * t540 + t517 * t535;
t511 = (t535 * t541 + t536 * t540) * qJD(1);
t528 = qJDD(2) + qJDD(3);
t529 = qJD(2) + qJD(3);
t441 = (t510 * t529 - t487) * pkin(9) + (t510 * t511 + t528) * pkin(3) + t462;
t463 = t535 * t482 + t540 * t483;
t486 = -qJD(3) * t511 - t516 * t535 + t517 * t540;
t504 = pkin(3) * t529 - pkin(9) * t511;
t506 = t510 ^ 2;
t443 = -pkin(3) * t506 + pkin(9) * t486 - t504 * t529 + t463;
t534 = sin(qJ(4));
t539 = cos(qJ(4));
t421 = t539 * t441 - t534 * t443;
t497 = t510 * t539 - t511 * t534;
t459 = qJD(4) * t497 + t486 * t534 + t487 * t539;
t498 = t510 * t534 + t511 * t539;
t525 = qJDD(4) + t528;
t526 = qJD(4) + t529;
t417 = (t497 * t526 - t459) * qJ(5) + (t497 * t498 + t525) * pkin(4) + t421;
t422 = t534 * t441 + t539 * t443;
t458 = -qJD(4) * t498 + t486 * t539 - t487 * t534;
t490 = pkin(4) * t526 - qJ(5) * t498;
t496 = t497 ^ 2;
t419 = -pkin(4) * t496 + qJ(5) * t458 - t490 * t526 + t422;
t531 = sin(pkin(11));
t532 = cos(pkin(11));
t475 = t497 * t532 - t498 * t531;
t414 = t531 * t417 + t532 * t419 + t475 * t564;
t432 = t458 * t532 - t459 * t531;
t476 = t497 * t531 + t498 * t532;
t451 = -mrSges(6,1) * t475 + mrSges(6,2) * t476;
t467 = mrSges(6,1) * t526 - mrSges(6,3) * t476;
t452 = -pkin(5) * t475 - pkin(10) * t476;
t524 = t526 ^ 2;
t411 = -pkin(5) * t524 + pkin(10) * t525 + t452 * t475 + t414;
t558 = t537 * g(1) - t542 * g(2);
t551 = -qJDD(1) * pkin(1) - t558;
t488 = -t517 * pkin(2) + t520 * t561 + (-pkin(8) * t530 - pkin(7)) * t543 + t551;
t454 = -t486 * pkin(3) - t506 * pkin(9) + t511 * t504 + t488;
t424 = -t458 * pkin(4) - t496 * qJ(5) + t498 * t490 + qJDD(5) + t454;
t433 = t458 * t531 + t459 * t532;
t415 = t424 + (-t475 * t526 - t433) * pkin(10) + (t476 * t526 - t432) * pkin(5);
t533 = sin(qJ(6));
t538 = cos(qJ(6));
t408 = -t411 * t533 + t415 * t538;
t464 = -t476 * t533 + t526 * t538;
t427 = qJD(6) * t464 + t433 * t538 + t525 * t533;
t431 = qJDD(6) - t432;
t465 = t476 * t538 + t526 * t533;
t444 = -mrSges(7,1) * t464 + mrSges(7,2) * t465;
t474 = qJD(6) - t475;
t445 = -mrSges(7,2) * t474 + mrSges(7,3) * t464;
t405 = m(7) * t408 + mrSges(7,1) * t431 - mrSges(7,3) * t427 - t444 * t465 + t445 * t474;
t409 = t411 * t538 + t415 * t533;
t426 = -qJD(6) * t465 - t433 * t533 + t525 * t538;
t446 = mrSges(7,1) * t474 - mrSges(7,3) * t465;
t406 = m(7) * t409 - mrSges(7,2) * t431 + mrSges(7,3) * t426 + t444 * t464 - t446 * t474;
t554 = -t405 * t533 + t538 * t406;
t392 = m(6) * t414 - mrSges(6,2) * t525 + mrSges(6,3) * t432 + t451 * t475 - t467 * t526 + t554;
t552 = -t532 * t417 + t531 * t419;
t413 = -0.2e1 * qJD(5) * t476 - t552;
t466 = -mrSges(6,2) * t526 + mrSges(6,3) * t475;
t410 = -t525 * pkin(5) - t524 * pkin(10) + (t564 + t452) * t476 + t552;
t550 = -m(7) * t410 + t426 * mrSges(7,1) - mrSges(7,2) * t427 + t464 * t445 - t446 * t465;
t401 = m(6) * t413 + mrSges(6,1) * t525 - mrSges(6,3) * t433 - t451 * t476 + t466 * t526 + t550;
t389 = t531 * t392 + t532 * t401;
t477 = -mrSges(5,1) * t497 + mrSges(5,2) * t498;
t489 = -mrSges(5,2) * t526 + mrSges(5,3) * t497;
t386 = m(5) * t421 + mrSges(5,1) * t525 - mrSges(5,3) * t459 - t477 * t498 + t489 * t526 + t389;
t491 = mrSges(5,1) * t526 - mrSges(5,3) * t498;
t555 = t532 * t392 - t401 * t531;
t387 = m(5) * t422 - mrSges(5,2) * t525 + mrSges(5,3) * t458 + t477 * t497 - t491 * t526 + t555;
t380 = t539 * t386 + t534 * t387;
t499 = -mrSges(4,1) * t510 + mrSges(4,2) * t511;
t502 = -mrSges(4,2) * t529 + mrSges(4,3) * t510;
t377 = m(4) * t462 + mrSges(4,1) * t528 - mrSges(4,3) * t487 - t499 * t511 + t502 * t529 + t380;
t503 = mrSges(4,1) * t529 - mrSges(4,3) * t511;
t556 = -t386 * t534 + t539 * t387;
t378 = m(4) * t463 - mrSges(4,2) * t528 + mrSges(4,3) * t486 + t499 * t510 - t503 * t529 + t556;
t372 = t540 * t377 + t535 * t378;
t395 = t538 * t405 + t533 * t406;
t560 = qJD(1) * t541;
t557 = -t377 * t535 + t540 * t378;
t549 = -m(6) * t424 + t432 * mrSges(6,1) - t433 * mrSges(6,2) + t475 * t466 - t476 * t467 - t395;
t435 = Ifges(7,4) * t465 + Ifges(7,2) * t464 + Ifges(7,6) * t474;
t436 = Ifges(7,1) * t465 + Ifges(7,4) * t464 + Ifges(7,5) * t474;
t548 = mrSges(7,1) * t408 - mrSges(7,2) * t409 + Ifges(7,5) * t427 + Ifges(7,6) * t426 + Ifges(7,3) * t431 + t435 * t465 - t436 * t464;
t547 = -m(5) * t454 + t458 * mrSges(5,1) - t459 * mrSges(5,2) + t497 * t489 - t498 * t491 + t549;
t434 = Ifges(7,5) * t465 + Ifges(7,6) * t464 + Ifges(7,3) * t474;
t398 = -mrSges(7,1) * t410 + mrSges(7,3) * t409 + Ifges(7,4) * t427 + Ifges(7,2) * t426 + Ifges(7,6) * t431 - t434 * t465 + t436 * t474;
t399 = mrSges(7,2) * t410 - mrSges(7,3) * t408 + Ifges(7,1) * t427 + Ifges(7,4) * t426 + Ifges(7,5) * t431 + t434 * t464 - t435 * t474;
t448 = Ifges(6,4) * t476 + Ifges(6,2) * t475 + Ifges(6,6) * t526;
t449 = Ifges(6,1) * t476 + Ifges(6,4) * t475 + Ifges(6,5) * t526;
t469 = Ifges(5,4) * t498 + Ifges(5,2) * t497 + Ifges(5,6) * t526;
t470 = Ifges(5,1) * t498 + Ifges(5,4) * t497 + Ifges(5,5) * t526;
t546 = mrSges(5,1) * t421 + mrSges(6,1) * t413 - mrSges(5,2) * t422 - mrSges(6,2) * t414 + Ifges(6,6) * t432 + pkin(4) * t389 + pkin(5) * t550 + pkin(10) * t554 + t538 * t398 + t533 * t399 - t475 * t449 - t497 * t470 + Ifges(6,5) * t433 + t476 * t448 + Ifges(5,6) * t458 + Ifges(5,5) * t459 + t498 * t469 + (Ifges(6,3) + Ifges(5,3)) * t525;
t545 = m(4) * t488 - t486 * mrSges(4,1) + t487 * mrSges(4,2) - t510 * t502 + t511 * t503 - t547;
t493 = Ifges(4,4) * t511 + Ifges(4,2) * t510 + Ifges(4,6) * t529;
t494 = Ifges(4,1) * t511 + Ifges(4,4) * t510 + Ifges(4,5) * t529;
t544 = mrSges(4,1) * t462 - mrSges(4,2) * t463 + Ifges(4,5) * t487 + Ifges(4,6) * t486 + Ifges(4,3) * t528 + pkin(3) * t380 + t511 * t493 - t510 * t494 + t546;
t519 = -qJD(2) * mrSges(3,2) + mrSges(3,3) * t560;
t518 = qJD(2) * mrSges(3,1) - mrSges(3,3) * t561;
t515 = (-mrSges(3,1) * t541 + mrSges(3,2) * t536) * qJD(1);
t512 = -t543 * pkin(7) + t551;
t509 = Ifges(3,5) * qJD(2) + (Ifges(3,1) * t536 + Ifges(3,4) * t541) * qJD(1);
t508 = Ifges(3,6) * qJD(2) + (Ifges(3,4) * t536 + Ifges(3,2) * t541) * qJD(1);
t500 = -t541 * g(3) - t562;
t492 = Ifges(4,5) * t511 + Ifges(4,6) * t510 + Ifges(4,3) * t529;
t468 = Ifges(5,5) * t498 + Ifges(5,6) * t497 + Ifges(5,3) * t526;
t447 = Ifges(6,5) * t476 + Ifges(6,6) * t475 + Ifges(6,3) * t526;
t382 = -mrSges(6,1) * t424 + mrSges(6,3) * t414 + Ifges(6,4) * t433 + Ifges(6,2) * t432 + Ifges(6,6) * t525 - pkin(5) * t395 - t447 * t476 + t449 * t526 - t548;
t381 = mrSges(6,2) * t424 - mrSges(6,3) * t413 + Ifges(6,1) * t433 + Ifges(6,4) * t432 + Ifges(6,5) * t525 - pkin(10) * t395 - t398 * t533 + t399 * t538 + t447 * t475 - t448 * t526;
t373 = mrSges(5,2) * t454 - mrSges(5,3) * t421 + Ifges(5,1) * t459 + Ifges(5,4) * t458 + Ifges(5,5) * t525 - qJ(5) * t389 + t381 * t532 - t382 * t531 + t468 * t497 - t469 * t526;
t371 = -mrSges(5,1) * t454 + mrSges(5,3) * t422 + Ifges(5,4) * t459 + Ifges(5,2) * t458 + Ifges(5,6) * t525 + pkin(4) * t549 + qJ(5) * t555 + t531 * t381 + t532 * t382 - t498 * t468 + t526 * t470;
t370 = mrSges(4,2) * t488 - mrSges(4,3) * t462 + Ifges(4,1) * t487 + Ifges(4,4) * t486 + Ifges(4,5) * t528 - pkin(9) * t380 - t371 * t534 + t373 * t539 + t492 * t510 - t493 * t529;
t369 = -mrSges(4,1) * t488 + mrSges(4,3) * t463 + Ifges(4,4) * t487 + Ifges(4,2) * t486 + Ifges(4,6) * t528 + pkin(3) * t547 + pkin(9) * t556 + t539 * t371 + t534 * t373 - t511 * t492 + t529 * t494;
t1 = [Ifges(2,3) * qJDD(1) + mrSges(2,1) * t558 - mrSges(2,2) * t553 + t536 * (mrSges(3,2) * t512 - mrSges(3,3) * t500 + Ifges(3,1) * t516 + Ifges(3,4) * t517 + Ifges(3,5) * qJDD(2) - pkin(8) * t372 - qJD(2) * t508 - t535 * t369 + t540 * t370) + t541 * (-mrSges(3,1) * t512 + mrSges(3,3) * t501 + Ifges(3,4) * t516 + Ifges(3,2) * t517 + Ifges(3,6) * qJDD(2) - pkin(2) * t545 + pkin(8) * t557 + qJD(2) * t509 + t540 * t369 + t535 * t370) + pkin(1) * ((-t518 * t536 + t519 * t541) * qJD(1) - t545 + t517 * mrSges(3,1) - m(3) * t512 - t516 * mrSges(3,2)) + pkin(7) * (t541 * (m(3) * t501 - qJDD(2) * mrSges(3,2) + mrSges(3,3) * t517 - qJD(2) * t518 + t515 * t560 + t557) - t536 * (m(3) * t500 + qJDD(2) * mrSges(3,1) - mrSges(3,3) * t516 + qJD(2) * t519 - t515 * t561 + t372)); Ifges(3,3) * qJDD(2) + t544 + (t536 * t508 - t541 * t509) * qJD(1) + Ifges(3,6) * t517 + Ifges(3,5) * t516 + mrSges(3,1) * t500 - mrSges(3,2) * t501 + pkin(2) * t372; t544; t546; -t549; t548;];
tauJ  = t1;
