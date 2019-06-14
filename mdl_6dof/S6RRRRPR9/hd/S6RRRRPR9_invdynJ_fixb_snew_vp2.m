% Calculate vector of inverse dynamics joint torques for with Newton-Euler
% S6RRRRPR9
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
% pkin [12x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,alpha2,d1,d2,d3,d4,d6,theta5]';
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
% Datum: 2019-05-07 22:28
% Revision: 89c353f7eff3bd693eda4e29f35b2761dbc3ada0 (2019-05-03)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauJ = S6RRRRPR9_invdynJ_fixb_snew_vp2(qJ, qJD, qJDD, g, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(6,1),zeros(3,1),zeros(12,1),zeros(7,1),zeros(7,3),zeros(7,6)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRRRPR9_invdynJ_fixb_snew_vp2: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RRRRPR9_invdynJ_fixb_snew_vp2: qJD has to be [6x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [6 1]), ...
  'S6RRRRPR9_invdynJ_fixb_snew_vp2: qJDD has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6RRRRPR9_invdynJ_fixb_snew_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [12 1]), ...
  'S6RRRRPR9_invdynJ_fixb_snew_vp2: pkin has to be [12x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6RRRRPR9_invdynJ_fixb_snew_vp2: m has to be [7x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [7,3]), ...
  'S6RRRRPR9_invdynJ_fixb_snew_vp2: mrSges has to be [7x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [7 6]), ...
  'S6RRRRPR9_invdynJ_fixb_snew_vp2: Ifges has to be [7x6] (double)'); 

%% Symbolic Calculation
% From invdyn_fixb_NewtonEuler_linkframe_tauJ_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-05-07 22:14:53
% EndTime: 2019-05-07 22:15:15
% DurationCPUTime: 17.08s
% Computational Cost: add. (276275->362), mult. (591202->471), div. (0->0), fcn. (483507->14), ass. (0->152)
t584 = cos(qJ(4));
t547 = sin(pkin(6));
t583 = pkin(8) * t547;
t549 = cos(pkin(6));
t582 = g(3) * t549;
t559 = qJD(1) ^ 2;
t554 = sin(qJ(1));
t558 = cos(qJ(1));
t572 = t554 * g(1) - g(2) * t558;
t531 = qJDD(1) * pkin(1) + t559 * t583 + t572;
t581 = t531 * t549;
t553 = sin(qJ(2));
t580 = t547 * t553;
t557 = cos(qJ(2));
t579 = t547 * t557;
t577 = qJD(1) * t547;
t534 = (-pkin(2) * t557 - pkin(9) * t553) * t577;
t543 = qJD(1) * t549 + qJD(2);
t541 = t543 ^ 2;
t542 = qJDD(1) * t549 + qJDD(2);
t576 = qJD(1) * t557;
t567 = -g(1) * t558 - g(2) * t554;
t575 = qJDD(1) * t547;
t532 = -pkin(1) * t559 + pkin(8) * t575 + t567;
t578 = t557 * t532 + t553 * t581;
t487 = -pkin(2) * t541 + pkin(9) * t542 + (-g(3) * t553 + t534 * t576) * t547 + t578;
t535 = (qJD(2) * t576 + qJDD(1) * t553) * t547;
t574 = t553 * t577;
t536 = -qJD(2) * t574 + t557 * t575;
t488 = -pkin(2) * t536 - pkin(9) * t535 - t582 + (-t531 + (pkin(2) * t553 - pkin(9) * t557) * t543 * qJD(1)) * t547;
t552 = sin(qJ(3));
t556 = cos(qJ(3));
t456 = -t487 * t552 + t556 * t488;
t523 = t543 * t556 - t552 * t574;
t504 = qJD(3) * t523 + t535 * t556 + t542 * t552;
t524 = t543 * t552 + t556 * t574;
t528 = qJDD(3) - t536;
t573 = t547 * t576;
t539 = qJD(3) - t573;
t443 = (t523 * t539 - t504) * pkin(10) + (t523 * t524 + t528) * pkin(3) + t456;
t457 = t556 * t487 + t552 * t488;
t503 = -qJD(3) * t524 - t535 * t552 + t542 * t556;
t513 = pkin(3) * t539 - pkin(10) * t524;
t522 = t523 ^ 2;
t450 = -pkin(3) * t522 + pkin(10) * t503 - t513 * t539 + t457;
t551 = sin(qJ(4));
t438 = t551 * t443 + t584 * t450;
t509 = t551 * t523 + t524 * t584;
t468 = qJD(4) * t509 - t584 * t503 + t504 * t551;
t508 = -t584 * t523 + t524 * t551;
t482 = mrSges(5,1) * t508 + mrSges(5,2) * t509;
t538 = qJD(4) + t539;
t495 = mrSges(5,1) * t538 - mrSges(5,3) * t509;
t527 = qJDD(4) + t528;
t481 = pkin(4) * t508 - qJ(5) * t509;
t537 = t538 ^ 2;
t432 = -pkin(4) * t537 + qJ(5) * t527 - t481 * t508 + t438;
t505 = -g(3) * t579 - t553 * t532 + t557 * t581;
t486 = -pkin(2) * t542 - pkin(9) * t541 + t534 * t574 - t505;
t452 = -pkin(3) * t503 - pkin(10) * t522 + t524 * t513 + t486;
t469 = -t508 * qJD(4) + t551 * t503 + t504 * t584;
t435 = (t508 * t538 - t469) * qJ(5) + (t509 * t538 + t468) * pkin(4) + t452;
t546 = sin(pkin(12));
t548 = cos(pkin(12));
t493 = t509 * t548 + t538 * t546;
t427 = -0.2e1 * qJD(5) * t493 - t432 * t546 + t548 * t435;
t459 = t469 * t548 + t527 * t546;
t492 = -t509 * t546 + t538 * t548;
t425 = (t492 * t508 - t459) * pkin(11) + (t492 * t493 + t468) * pkin(5) + t427;
t428 = 0.2e1 * qJD(5) * t492 + t548 * t432 + t546 * t435;
t458 = -t469 * t546 + t527 * t548;
t476 = pkin(5) * t508 - pkin(11) * t493;
t491 = t492 ^ 2;
t426 = -pkin(5) * t491 + pkin(11) * t458 - t476 * t508 + t428;
t550 = sin(qJ(6));
t555 = cos(qJ(6));
t423 = t425 * t555 - t426 * t550;
t472 = t492 * t555 - t493 * t550;
t441 = qJD(6) * t472 + t458 * t550 + t459 * t555;
t473 = t492 * t550 + t493 * t555;
t451 = -mrSges(7,1) * t472 + mrSges(7,2) * t473;
t507 = qJD(6) + t508;
t454 = -mrSges(7,2) * t507 + mrSges(7,3) * t472;
t467 = qJDD(6) + t468;
t419 = m(7) * t423 + mrSges(7,1) * t467 - mrSges(7,3) * t441 - t451 * t473 + t454 * t507;
t424 = t425 * t550 + t426 * t555;
t440 = -qJD(6) * t473 + t458 * t555 - t459 * t550;
t455 = mrSges(7,1) * t507 - mrSges(7,3) * t473;
t420 = m(7) * t424 - mrSges(7,2) * t467 + mrSges(7,3) * t440 + t451 * t472 - t455 * t507;
t411 = t555 * t419 + t550 * t420;
t474 = -mrSges(6,1) * t492 + mrSges(6,2) * t493;
t566 = -mrSges(6,2) * t508 + mrSges(6,3) * t492;
t409 = m(6) * t427 + t468 * mrSges(6,1) - t459 * mrSges(6,3) - t493 * t474 + t508 * t566 + t411;
t475 = mrSges(6,1) * t508 - mrSges(6,3) * t493;
t568 = -t419 * t550 + t555 * t420;
t410 = m(6) * t428 - mrSges(6,2) * t468 + mrSges(6,3) * t458 + t474 * t492 - t475 * t508 + t568;
t569 = -t409 * t546 + t548 * t410;
t402 = m(5) * t438 - mrSges(5,2) * t527 - mrSges(5,3) * t468 - t482 * t508 - t495 * t538 + t569;
t437 = t443 * t584 - t551 * t450;
t431 = -t527 * pkin(4) - t537 * qJ(5) + t509 * t481 + qJDD(5) - t437;
t429 = -t458 * pkin(5) - t491 * pkin(11) + t493 * t476 + t431;
t565 = m(7) * t429 - t440 * mrSges(7,1) + mrSges(7,2) * t441 - t472 * t454 + t455 * t473;
t422 = m(6) * t431 - t458 * mrSges(6,1) + mrSges(6,2) * t459 + t475 * t493 - t492 * t566 + t565;
t494 = -mrSges(5,2) * t538 - mrSges(5,3) * t508;
t415 = m(5) * t437 + mrSges(5,1) * t527 - mrSges(5,3) * t469 - t482 * t509 + t494 * t538 - t422;
t394 = t551 * t402 + t584 * t415;
t510 = -mrSges(4,1) * t523 + mrSges(4,2) * t524;
t511 = -mrSges(4,2) * t539 + mrSges(4,3) * t523;
t392 = m(4) * t456 + mrSges(4,1) * t528 - mrSges(4,3) * t504 - t510 * t524 + t511 * t539 + t394;
t512 = mrSges(4,1) * t539 - mrSges(4,3) * t524;
t570 = t584 * t402 - t415 * t551;
t393 = m(4) * t457 - mrSges(4,2) * t528 + mrSges(4,3) * t503 + t510 * t523 - t512 * t539 + t570;
t387 = t556 * t392 + t552 * t393;
t404 = t548 * t409 + t546 * t410;
t571 = -t392 * t552 + t556 * t393;
t564 = m(5) * t452 + t468 * mrSges(5,1) + mrSges(5,2) * t469 + t508 * t494 + t509 * t495 + t404;
t444 = Ifges(7,5) * t473 + Ifges(7,6) * t472 + Ifges(7,3) * t507;
t446 = Ifges(7,1) * t473 + Ifges(7,4) * t472 + Ifges(7,5) * t507;
t412 = -mrSges(7,1) * t429 + mrSges(7,3) * t424 + Ifges(7,4) * t441 + Ifges(7,2) * t440 + Ifges(7,6) * t467 - t444 * t473 + t446 * t507;
t445 = Ifges(7,4) * t473 + Ifges(7,2) * t472 + Ifges(7,6) * t507;
t413 = mrSges(7,2) * t429 - mrSges(7,3) * t423 + Ifges(7,1) * t441 + Ifges(7,4) * t440 + Ifges(7,5) * t467 + t444 * t472 - t445 * t507;
t460 = Ifges(6,5) * t493 + Ifges(6,6) * t492 + Ifges(6,3) * t508;
t462 = Ifges(6,1) * t493 + Ifges(6,4) * t492 + Ifges(6,5) * t508;
t396 = -mrSges(6,1) * t431 + mrSges(6,3) * t428 + Ifges(6,4) * t459 + Ifges(6,2) * t458 + Ifges(6,6) * t468 - pkin(5) * t565 + pkin(11) * t568 + t555 * t412 + t550 * t413 - t493 * t460 + t508 * t462;
t461 = Ifges(6,4) * t493 + Ifges(6,2) * t492 + Ifges(6,6) * t508;
t398 = mrSges(6,2) * t431 - mrSges(6,3) * t427 + Ifges(6,1) * t459 + Ifges(6,4) * t458 + Ifges(6,5) * t468 - pkin(11) * t411 - t412 * t550 + t413 * t555 + t460 * t492 - t461 * t508;
t478 = Ifges(5,4) * t509 - Ifges(5,2) * t508 + Ifges(5,6) * t538;
t479 = Ifges(5,1) * t509 - Ifges(5,4) * t508 + Ifges(5,5) * t538;
t563 = -mrSges(5,1) * t437 + mrSges(5,2) * t438 - Ifges(5,5) * t469 + Ifges(5,6) * t468 - Ifges(5,3) * t527 + pkin(4) * t422 - qJ(5) * t569 - t548 * t396 - t546 * t398 - t509 * t478 - t508 * t479;
t562 = mrSges(7,1) * t423 - mrSges(7,2) * t424 + Ifges(7,5) * t441 + Ifges(7,6) * t440 + Ifges(7,3) * t467 + t473 * t445 - t472 * t446;
t561 = -m(4) * t486 + t503 * mrSges(4,1) - mrSges(4,2) * t504 + t523 * t511 - t512 * t524 - t564;
t498 = Ifges(4,4) * t524 + Ifges(4,2) * t523 + Ifges(4,6) * t539;
t499 = Ifges(4,1) * t524 + Ifges(4,4) * t523 + Ifges(4,5) * t539;
t560 = mrSges(4,1) * t456 - mrSges(4,2) * t457 + Ifges(4,5) * t504 + Ifges(4,6) * t503 + Ifges(4,3) * t528 + pkin(3) * t394 + t524 * t498 - t523 * t499 - t563;
t533 = (-mrSges(3,1) * t557 + mrSges(3,2) * t553) * t577;
t530 = -mrSges(3,2) * t543 + mrSges(3,3) * t573;
t529 = mrSges(3,1) * t543 - mrSges(3,3) * t574;
t517 = -t531 * t547 - t582;
t516 = Ifges(3,5) * t543 + (Ifges(3,1) * t553 + Ifges(3,4) * t557) * t577;
t515 = Ifges(3,6) * t543 + (Ifges(3,4) * t553 + Ifges(3,2) * t557) * t577;
t514 = Ifges(3,3) * t543 + (Ifges(3,5) * t553 + Ifges(3,6) * t557) * t577;
t506 = -g(3) * t580 + t578;
t497 = Ifges(4,5) * t524 + Ifges(4,6) * t523 + Ifges(4,3) * t539;
t477 = Ifges(5,5) * t509 - Ifges(5,6) * t508 + Ifges(5,3) * t538;
t399 = m(3) * t505 + mrSges(3,1) * t542 - mrSges(3,3) * t535 + t530 * t543 - t533 * t574 + t561;
t388 = t538 * t479 + Ifges(5,6) * t527 - t509 * t477 + t492 * t462 - t493 * t461 + Ifges(5,4) * t469 - Ifges(6,6) * t458 - Ifges(6,5) * t459 - t562 - mrSges(5,1) * t452 + mrSges(5,3) * t438 + mrSges(6,2) * t428 - mrSges(6,1) * t427 - pkin(5) * t411 - pkin(4) * t404 + (-Ifges(5,2) - Ifges(6,3)) * t468;
t386 = m(3) * t506 - mrSges(3,2) * t542 + mrSges(3,3) * t536 - t529 * t543 + t533 * t573 + t571;
t385 = mrSges(5,2) * t452 - mrSges(5,3) * t437 + Ifges(5,1) * t469 - Ifges(5,4) * t468 + Ifges(5,5) * t527 - qJ(5) * t404 - t396 * t546 + t398 * t548 - t477 * t508 - t478 * t538;
t384 = mrSges(4,2) * t486 - mrSges(4,3) * t456 + Ifges(4,1) * t504 + Ifges(4,4) * t503 + Ifges(4,5) * t528 - pkin(10) * t394 + t385 * t584 - t551 * t388 + t523 * t497 - t539 * t498;
t383 = -mrSges(4,1) * t486 + mrSges(4,3) * t457 + Ifges(4,4) * t504 + Ifges(4,2) * t503 + Ifges(4,6) * t528 - pkin(3) * t564 + pkin(10) * t570 + t551 * t385 + t388 * t584 - t524 * t497 + t539 * t499;
t382 = Ifges(3,5) * t535 + Ifges(3,6) * t536 + Ifges(3,3) * t542 + mrSges(3,1) * t505 - mrSges(3,2) * t506 + t552 * t384 + t556 * t383 + pkin(2) * t561 + pkin(9) * t571 + (t515 * t553 - t516 * t557) * t577;
t1 = [Ifges(2,3) * qJDD(1) + mrSges(2,1) * t572 - mrSges(2,2) * t567 + (mrSges(3,2) * t517 - mrSges(3,3) * t505 + Ifges(3,1) * t535 + Ifges(3,4) * t536 + Ifges(3,5) * t542 - pkin(9) * t387 - t383 * t552 + t384 * t556 + t514 * t573 - t515 * t543) * t580 + (-mrSges(3,1) * t517 + mrSges(3,3) * t506 + Ifges(3,4) * t535 + Ifges(3,2) * t536 + Ifges(3,6) * t542 - pkin(2) * t387 - t514 * t574 + t543 * t516 - t560) * t579 + t549 * t382 + pkin(1) * ((t386 * t553 + t399 * t557) * t549 + (-m(3) * t517 + t536 * mrSges(3,1) - t535 * mrSges(3,2) + (-t529 * t553 + t530 * t557) * t577 - t387) * t547) + (t386 * t557 - t399 * t553) * t583; t382; t560; -t563; t422; t562;];
tauJ  = t1;
