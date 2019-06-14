% Calculate vector of inverse dynamics joint torques for with Newton-Euler
% S6PRPRRR7
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
%   pkin=[a2,a3,a4,a5,a6,alpha2,alpha3,alpha4,d2,d4,d5,d6,theta1,theta3]';
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
% Datum: 2019-05-05 02:02
% Revision: 89c353f7eff3bd693eda4e29f35b2761dbc3ada0 (2019-05-03)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauJ = S6PRPRRR7_invdynJ_fixb_snew_vp2(qJ, qJD, qJDD, g, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(6,1),zeros(3,1),zeros(14,1),zeros(7,1),zeros(7,3),zeros(7,6)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6PRPRRR7_invdynJ_fixb_snew_vp2: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6PRPRRR7_invdynJ_fixb_snew_vp2: qJD has to be [6x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [6 1]), ...
  'S6PRPRRR7_invdynJ_fixb_snew_vp2: qJDD has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6PRPRRR7_invdynJ_fixb_snew_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [14 1]), ...
  'S6PRPRRR7_invdynJ_fixb_snew_vp2: pkin has to be [14x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6PRPRRR7_invdynJ_fixb_snew_vp2: m has to be [7x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [7,3]), ...
  'S6PRPRRR7_invdynJ_fixb_snew_vp2: mrSges has to be [7x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [7 6]), ...
  'S6PRPRRR7_invdynJ_fixb_snew_vp2: Ifges has to be [7x6] (double)'); 

%% Symbolic Calculation
% From invdyn_fixb_NewtonEuler_linkframe_tauJ_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-05-05 01:52:21
% EndTime: 2019-05-05 01:52:31
% DurationCPUTime: 9.15s
% Computational Cost: add. (133018->299), mult. (375646->425), div. (0->0), fcn. (320586->18), ass. (0->151)
t524 = sin(pkin(14));
t527 = sin(pkin(7));
t529 = cos(pkin(14));
t532 = cos(pkin(7));
t536 = sin(qJ(4));
t531 = cos(pkin(8));
t540 = cos(qJ(4));
t567 = t531 * t540;
t526 = sin(pkin(8));
t572 = t526 * t540;
t546 = t527 * (-t524 * t536 + t529 * t567) + t532 * t572;
t496 = t546 * qJD(2);
t568 = t531 * t536;
t573 = t526 * t536;
t548 = t532 * t573 + (t524 * t540 + t529 * t568) * t527;
t497 = t548 * qJD(2);
t480 = -t497 * qJD(4) + qJDD(2) * t546;
t525 = sin(pkin(13));
t530 = cos(pkin(13));
t519 = g(1) * t525 - g(2) * t530;
t523 = -g(3) + qJDD(1);
t528 = sin(pkin(6));
t533 = cos(pkin(6));
t583 = t519 * t533 + t523 * t528;
t520 = -g(1) * t530 - g(2) * t525;
t537 = sin(qJ(2));
t541 = cos(qJ(2));
t493 = t541 * t520 + t583 * t537;
t542 = qJD(2) ^ 2;
t579 = qJ(3) * t527;
t491 = -pkin(2) * t542 + qJDD(2) * t579 + t493;
t580 = pkin(10) * t526;
t555 = -pkin(3) * t529 - t524 * t580;
t565 = qJD(2) * t527;
t506 = t555 * t565;
t570 = t527 * t531;
t574 = t526 * t532;
t550 = pkin(10) * (t529 * t570 + t574);
t507 = qJD(2) * t550;
t492 = -t520 * t537 + t583 * t541;
t490 = qJDD(2) * pkin(2) + t542 * t579 + t492;
t509 = -t519 * t528 + t523 * t533;
t560 = qJD(3) * t565;
t569 = t529 * t532;
t571 = t527 * t529;
t561 = t490 * t569 + t509 * t571 - 0.2e1 * t524 * t560;
t449 = (pkin(3) * qJDD(2) + qJD(2) * t507) * t532 + (-t491 + (-pkin(10) * qJDD(2) * t531 - qJD(2) * t506) * t527) * t524 + t561;
t575 = t524 * t532;
t576 = t524 * t527;
t582 = 0.2e1 * t529;
t465 = t490 * t575 + t529 * t491 + t509 * t576 + t560 * t582;
t581 = pkin(10) * t524;
t511 = (pkin(3) * t532 - t570 * t581) * qJD(2);
t450 = (t506 * t571 - t511 * t532) * qJD(2) + qJDD(2) * t550 + t465;
t563 = t532 * t509 + qJDD(3);
t463 = (-t490 + t555 * qJDD(2) + (-t507 * t529 + t511 * t524) * qJD(2)) * t527 + t563;
t435 = -t536 * t450 + (t449 * t531 + t463 * t526) * t540;
t436 = t449 * t568 + t540 * t450 + t463 * t573;
t478 = -mrSges(5,1) * t496 + mrSges(5,2) * t497;
t551 = -t526 * t571 + t531 * t532;
t508 = qJD(2) * t551 + qJD(4);
t486 = mrSges(5,1) * t508 - mrSges(5,3) * t497;
t505 = qJDD(2) * t551 + qJDD(4);
t479 = -pkin(4) * t496 - pkin(11) * t497;
t504 = t508 ^ 2;
t432 = -pkin(4) * t504 + pkin(11) * t505 + t479 * t496 + t436;
t437 = -t526 * t449 + t531 * t463;
t481 = t496 * qJD(4) + qJDD(2) * t548;
t434 = (-t496 * t508 - t481) * pkin(11) + (t497 * t508 - t480) * pkin(4) + t437;
t535 = sin(qJ(5));
t539 = cos(qJ(5));
t428 = t539 * t432 + t535 * t434;
t483 = -t497 * t535 + t508 * t539;
t484 = t497 * t539 + t508 * t535;
t467 = -pkin(5) * t483 - pkin(12) * t484;
t477 = qJDD(5) - t480;
t495 = qJD(5) - t496;
t494 = t495 ^ 2;
t426 = -pkin(5) * t494 + pkin(12) * t477 + t467 * t483 + t428;
t431 = -t505 * pkin(4) - t504 * pkin(11) + t497 * t479 - t435;
t459 = -qJD(5) * t484 - t481 * t535 + t505 * t539;
t460 = qJD(5) * t483 + t481 * t539 + t505 * t535;
t429 = (-t483 * t495 - t460) * pkin(12) + (t484 * t495 - t459) * pkin(5) + t431;
t534 = sin(qJ(6));
t538 = cos(qJ(6));
t422 = -t426 * t534 + t429 * t538;
t469 = -t484 * t534 + t495 * t538;
t440 = qJD(6) * t469 + t460 * t538 + t477 * t534;
t470 = t484 * t538 + t495 * t534;
t448 = -mrSges(7,1) * t469 + mrSges(7,2) * t470;
t482 = qJD(6) - t483;
t451 = -mrSges(7,2) * t482 + mrSges(7,3) * t469;
t457 = qJDD(6) - t459;
t420 = m(7) * t422 + mrSges(7,1) * t457 - mrSges(7,3) * t440 - t448 * t470 + t451 * t482;
t423 = t426 * t538 + t429 * t534;
t439 = -qJD(6) * t470 - t460 * t534 + t477 * t538;
t452 = mrSges(7,1) * t482 - mrSges(7,3) * t470;
t421 = m(7) * t423 - mrSges(7,2) * t457 + mrSges(7,3) * t439 + t448 * t469 - t452 * t482;
t414 = -t420 * t534 + t538 * t421;
t466 = -mrSges(6,1) * t483 + mrSges(6,2) * t484;
t472 = mrSges(6,1) * t495 - mrSges(6,3) * t484;
t412 = m(6) * t428 - mrSges(6,2) * t477 + mrSges(6,3) * t459 + t466 * t483 - t472 * t495 + t414;
t427 = -t432 * t535 + t434 * t539;
t425 = -pkin(5) * t477 - pkin(12) * t494 + t467 * t484 - t427;
t424 = -m(7) * t425 + t439 * mrSges(7,1) - mrSges(7,2) * t440 + t469 * t451 - t452 * t470;
t471 = -mrSges(6,2) * t495 + mrSges(6,3) * t483;
t418 = m(6) * t427 + mrSges(6,1) * t477 - mrSges(6,3) * t460 - t466 * t484 + t471 * t495 + t424;
t558 = t539 * t412 - t418 * t535;
t403 = m(5) * t436 - mrSges(5,2) * t505 + mrSges(5,3) * t480 + t478 * t496 - t486 * t508 + t558;
t406 = t535 * t412 + t539 * t418;
t485 = -mrSges(5,2) * t508 + mrSges(5,3) * t496;
t405 = m(5) * t437 - mrSges(5,1) * t480 + mrSges(5,2) * t481 - t485 * t496 + t486 * t497 + t406;
t413 = t420 * t538 + t421 * t534;
t545 = -m(6) * t431 + t459 * mrSges(6,1) - mrSges(6,2) * t460 + t483 * t471 - t472 * t484 - t413;
t409 = m(5) * t435 + mrSges(5,1) * t505 - mrSges(5,3) * t481 - t478 * t497 + t485 * t508 + t545;
t393 = t403 * t568 - t526 * t405 + t409 * t567;
t464 = -t524 * t491 + t561;
t557 = -mrSges(4,1) * t529 + mrSges(4,2) * t524;
t510 = t557 * t565;
t553 = -mrSges(4,2) * t532 + mrSges(4,3) * t571;
t515 = t553 * qJD(2);
t554 = mrSges(4,1) * t532 - mrSges(4,3) * t576;
t390 = m(4) * t464 + t554 * qJDD(2) + (-t510 * t576 + t515 * t532) * qJD(2) + t393;
t397 = t540 * t403 - t536 * t409;
t514 = t554 * qJD(2);
t396 = m(4) * t465 + t553 * qJDD(2) + (t510 * t571 - t514 * t532) * qJD(2) + t397;
t566 = t390 * t569 + t396 * t575;
t392 = t403 * t573 + t531 * t405 + t409 * t572;
t559 = -t390 * t524 + t529 * t396;
t441 = Ifges(7,5) * t470 + Ifges(7,6) * t469 + Ifges(7,3) * t482;
t443 = Ifges(7,1) * t470 + Ifges(7,4) * t469 + Ifges(7,5) * t482;
t415 = -mrSges(7,1) * t425 + mrSges(7,3) * t423 + Ifges(7,4) * t440 + Ifges(7,2) * t439 + Ifges(7,6) * t457 - t441 * t470 + t443 * t482;
t442 = Ifges(7,4) * t470 + Ifges(7,2) * t469 + Ifges(7,6) * t482;
t416 = mrSges(7,2) * t425 - mrSges(7,3) * t422 + Ifges(7,1) * t440 + Ifges(7,4) * t439 + Ifges(7,5) * t457 + t441 * t469 - t442 * t482;
t453 = Ifges(6,5) * t484 + Ifges(6,6) * t483 + Ifges(6,3) * t495;
t454 = Ifges(6,4) * t484 + Ifges(6,2) * t483 + Ifges(6,6) * t495;
t398 = mrSges(6,2) * t431 - mrSges(6,3) * t427 + Ifges(6,1) * t460 + Ifges(6,4) * t459 + Ifges(6,5) * t477 - pkin(12) * t413 - t415 * t534 + t416 * t538 + t453 * t483 - t454 * t495;
t455 = Ifges(6,1) * t484 + Ifges(6,4) * t483 + Ifges(6,5) * t495;
t544 = mrSges(7,1) * t422 - mrSges(7,2) * t423 + Ifges(7,5) * t440 + Ifges(7,6) * t439 + Ifges(7,3) * t457 + t442 * t470 - t443 * t469;
t399 = -mrSges(6,1) * t431 + mrSges(6,3) * t428 + Ifges(6,4) * t460 + Ifges(6,2) * t459 + Ifges(6,6) * t477 - pkin(5) * t413 - t453 * t484 + t455 * t495 - t544;
t473 = Ifges(5,5) * t497 + Ifges(5,6) * t496 + Ifges(5,3) * t508;
t474 = Ifges(5,4) * t497 + Ifges(5,2) * t496 + Ifges(5,6) * t508;
t387 = mrSges(5,2) * t437 - mrSges(5,3) * t435 + Ifges(5,1) * t481 + Ifges(5,4) * t480 + Ifges(5,5) * t505 - pkin(11) * t406 + t398 * t539 - t399 * t535 + t473 * t496 - t474 * t508;
t475 = Ifges(5,1) * t497 + Ifges(5,4) * t496 + Ifges(5,5) * t508;
t543 = mrSges(6,1) * t427 - mrSges(6,2) * t428 + Ifges(6,5) * t460 + Ifges(6,6) * t459 + Ifges(6,3) * t477 + pkin(5) * t424 + pkin(12) * t414 + t538 * t415 + t534 * t416 + t484 * t454 - t483 * t455;
t388 = -mrSges(5,1) * t437 + mrSges(5,3) * t436 + Ifges(5,4) * t481 + Ifges(5,2) * t480 + Ifges(5,6) * t505 - pkin(4) * t406 - t497 * t473 + t508 * t475 - t543;
t549 = pkin(10) * t397 + t387 * t536 + t388 * t540;
t476 = -t527 * t490 + t563;
t391 = m(4) * t476 + (t557 * qJDD(2) + (t514 * t524 - t515 * t529) * qJD(2)) * t527 + t392;
t386 = mrSges(5,1) * t435 - mrSges(5,2) * t436 + Ifges(5,5) * t481 + Ifges(5,6) * t480 + Ifges(5,3) * t505 + pkin(4) * t545 + pkin(11) * t558 + t535 * t398 + t539 * t399 + t497 * t474 - t496 * t475;
t1 = [m(2) * t523 + t533 * (m(3) * t509 + t532 * t391 + (t390 * t529 + t396 * t524) * t527) + (t537 * (m(3) * t493 - mrSges(3,1) * t542 - qJDD(2) * mrSges(3,2) + t559) + t541 * (m(3) * t492 + qJDD(2) * mrSges(3,1) - mrSges(3,2) * t542 - t391 * t527 + t566)) * t528; Ifges(3,3) * qJDD(2) + mrSges(3,1) * t492 - mrSges(3,2) * t493 + t532 * (Ifges(4,3) * t532 * qJDD(2) + mrSges(4,1) * t464 - mrSges(4,2) * t465 + pkin(3) * t393 + t531 * t386) + pkin(2) * t566 + t549 * t574 + (t524 * (mrSges(4,2) * t476 - mrSges(4,3) * t464 + t540 * t387 - t536 * t388 - t392 * t580) + t529 * (-mrSges(4,1) * t476 + mrSges(4,3) * t465 - pkin(3) * t392 - t526 * t386) - pkin(2) * t391 + qJ(3) * t559 + (-t393 * t581 + t529 * t549) * t531 + ((Ifges(4,2) * t529 ^ 2 + (Ifges(4,1) * t524 + Ifges(4,4) * t582) * t524) * t527 + 0.2e1 * t532 * (Ifges(4,5) * t524 + Ifges(4,6) * t529)) * qJDD(2)) * t527; t391; t386; t543; t544;];
tauJ  = t1;
