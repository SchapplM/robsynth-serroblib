% Calculate vector of inverse dynamics joint torques for with Newton-Euler
% S6RRPRRR10
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
%   pkin=[a2,a3,a4,a5,a6,alpha2,d1,d2,d4,d5,d6,theta3]';
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
% Datum: 2019-05-06 23:51
% Revision: 89c353f7eff3bd693eda4e29f35b2761dbc3ada0 (2019-05-03)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauJ = S6RRPRRR10_invdynJ_fixb_snew_vp2(qJ, qJD, qJDD, g, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(6,1),zeros(3,1),zeros(12,1),zeros(7,1),zeros(7,3),zeros(7,6)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRPRRR10_invdynJ_fixb_snew_vp2: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RRPRRR10_invdynJ_fixb_snew_vp2: qJD has to be [6x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [6 1]), ...
  'S6RRPRRR10_invdynJ_fixb_snew_vp2: qJDD has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6RRPRRR10_invdynJ_fixb_snew_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [12 1]), ...
  'S6RRPRRR10_invdynJ_fixb_snew_vp2: pkin has to be [12x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6RRPRRR10_invdynJ_fixb_snew_vp2: m has to be [7x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [7,3]), ...
  'S6RRPRRR10_invdynJ_fixb_snew_vp2: mrSges has to be [7x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [7 6]), ...
  'S6RRPRRR10_invdynJ_fixb_snew_vp2: Ifges has to be [7x6] (double)'); 

%% Symbolic Calculation
% From invdyn_fixb_NewtonEuler_linkframe_tauJ_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-05-06 23:39:34
% EndTime: 2019-05-06 23:39:50
% DurationCPUTime: 15.46s
% Computational Cost: add. (238028->362), mult. (537833->471), div. (0->0), fcn. (440063->14), ass. (0->152)
t548 = sin(pkin(6));
t582 = pkin(8) * t548;
t550 = cos(pkin(6));
t581 = g(3) * t550;
t561 = qJD(1) ^ 2;
t555 = sin(qJ(1));
t560 = cos(qJ(1));
t571 = t555 * g(1) - g(2) * t560;
t533 = qJDD(1) * pkin(1) + t561 * t582 + t571;
t580 = t533 * t550;
t554 = sin(qJ(2));
t579 = t548 * t554;
t559 = cos(qJ(2));
t578 = t548 * t559;
t576 = qJD(1) * t548;
t535 = (-pkin(2) * t559 - qJ(3) * t554) * t576;
t544 = qJD(1) * t550 + qJD(2);
t542 = t544 ^ 2;
t543 = qJDD(1) * t550 + qJDD(2);
t575 = qJD(1) * t559;
t567 = -g(1) * t560 - g(2) * t555;
t574 = qJDD(1) * t548;
t534 = -pkin(1) * t561 + pkin(8) * t574 + t567;
t577 = t559 * t534 + t554 * t580;
t489 = -pkin(2) * t542 + qJ(3) * t543 + (-g(3) * t554 + t535 * t575) * t548 + t577;
t537 = (qJD(2) * t575 + qJDD(1) * t554) * t548;
t573 = t554 * t576;
t538 = -qJD(2) * t573 + t559 * t574;
t490 = -pkin(2) * t538 - t581 - qJ(3) * t537 + (-t533 + (pkin(2) * t554 - qJ(3) * t559) * t544 * qJD(1)) * t548;
t547 = sin(pkin(12));
t549 = cos(pkin(12));
t527 = t544 * t547 + t549 * t573;
t454 = -0.2e1 * qJD(3) * t527 - t489 * t547 + t549 * t490;
t514 = t537 * t549 + t543 * t547;
t526 = t544 * t549 - t547 * t573;
t572 = t548 * t575;
t445 = (-t526 * t572 - t514) * pkin(9) + (t526 * t527 - t538) * pkin(3) + t454;
t455 = 0.2e1 * qJD(3) * t526 + t549 * t489 + t547 * t490;
t513 = -t537 * t547 + t543 * t549;
t515 = -pkin(3) * t572 - pkin(9) * t527;
t525 = t526 ^ 2;
t452 = -pkin(3) * t525 + pkin(9) * t513 + t515 * t572 + t455;
t553 = sin(qJ(4));
t558 = cos(qJ(4));
t434 = t553 * t445 + t558 * t452;
t506 = t526 * t558 - t553 * t527;
t507 = t526 * t553 + t527 * t558;
t484 = -pkin(4) * t506 - pkin(10) * t507;
t530 = qJDD(4) - t538;
t540 = qJD(4) - t572;
t539 = t540 ^ 2;
t432 = -pkin(4) * t539 + pkin(10) * t530 + t484 * t506 + t434;
t502 = -g(3) * t578 - t554 * t534 + t559 * t580;
t488 = -pkin(2) * t543 - qJ(3) * t542 + t535 * t573 + qJDD(3) - t502;
t459 = -pkin(3) * t513 - pkin(9) * t525 + t527 * t515 + t488;
t474 = -t507 * qJD(4) + t513 * t558 - t553 * t514;
t475 = qJD(4) * t506 + t513 * t553 + t514 * t558;
t442 = (-t506 * t540 - t475) * pkin(10) + (t507 * t540 - t474) * pkin(4) + t459;
t552 = sin(qJ(5));
t557 = cos(qJ(5));
t427 = -t432 * t552 + t557 * t442;
t492 = -t507 * t552 + t540 * t557;
t458 = qJD(5) * t492 + t475 * t557 + t530 * t552;
t473 = qJDD(5) - t474;
t493 = t507 * t557 + t540 * t552;
t505 = qJD(5) - t506;
t425 = (t492 * t505 - t458) * pkin(11) + (t492 * t493 + t473) * pkin(5) + t427;
t428 = t557 * t432 + t552 * t442;
t457 = -qJD(5) * t493 - t475 * t552 + t530 * t557;
t478 = pkin(5) * t505 - pkin(11) * t493;
t491 = t492 ^ 2;
t426 = -pkin(5) * t491 + pkin(11) * t457 - t478 * t505 + t428;
t551 = sin(qJ(6));
t556 = cos(qJ(6));
t423 = t425 * t556 - t426 * t551;
t467 = t492 * t556 - t493 * t551;
t439 = qJD(6) * t467 + t457 * t551 + t458 * t556;
t468 = t492 * t551 + t493 * t556;
t453 = -mrSges(7,1) * t467 + mrSges(7,2) * t468;
t501 = qJD(6) + t505;
t460 = -mrSges(7,2) * t501 + mrSges(7,3) * t467;
t471 = qJDD(6) + t473;
t419 = m(7) * t423 + mrSges(7,1) * t471 - mrSges(7,3) * t439 - t453 * t468 + t460 * t501;
t424 = t425 * t551 + t426 * t556;
t438 = -qJD(6) * t468 + t457 * t556 - t458 * t551;
t461 = mrSges(7,1) * t501 - mrSges(7,3) * t468;
t420 = m(7) * t424 - mrSges(7,2) * t471 + mrSges(7,3) * t438 + t453 * t467 - t461 * t501;
t411 = t556 * t419 + t551 * t420;
t469 = -mrSges(6,1) * t492 + mrSges(6,2) * t493;
t476 = -mrSges(6,2) * t505 + mrSges(6,3) * t492;
t409 = m(6) * t427 + mrSges(6,1) * t473 - mrSges(6,3) * t458 - t469 * t493 + t476 * t505 + t411;
t477 = mrSges(6,1) * t505 - mrSges(6,3) * t493;
t568 = -t419 * t551 + t556 * t420;
t410 = m(6) * t428 - mrSges(6,2) * t473 + mrSges(6,3) * t457 + t469 * t492 - t477 * t505 + t568;
t405 = -t409 * t552 + t557 * t410;
t483 = -mrSges(5,1) * t506 + mrSges(5,2) * t507;
t495 = mrSges(5,1) * t540 - mrSges(5,3) * t507;
t402 = m(5) * t434 - mrSges(5,2) * t530 + mrSges(5,3) * t474 + t483 * t506 - t495 * t540 + t405;
t433 = t445 * t558 - t553 * t452;
t431 = -pkin(4) * t530 - pkin(10) * t539 + t507 * t484 - t433;
t429 = -pkin(5) * t457 - pkin(11) * t491 + t478 * t493 + t431;
t566 = m(7) * t429 - t438 * mrSges(7,1) + mrSges(7,2) * t439 - t467 * t460 + t461 * t468;
t421 = -m(6) * t431 + t457 * mrSges(6,1) - mrSges(6,2) * t458 + t492 * t476 - t477 * t493 - t566;
t494 = -mrSges(5,2) * t540 + mrSges(5,3) * t506;
t415 = m(5) * t433 + mrSges(5,1) * t530 - mrSges(5,3) * t475 - t483 * t507 + t494 * t540 + t421;
t396 = t553 * t402 + t558 * t415;
t508 = -mrSges(4,1) * t526 + mrSges(4,2) * t527;
t511 = mrSges(4,2) * t572 + mrSges(4,3) * t526;
t394 = m(4) * t454 - mrSges(4,1) * t538 - mrSges(4,3) * t514 - t508 * t527 - t511 * t572 + t396;
t512 = -mrSges(4,1) * t572 - mrSges(4,3) * t527;
t569 = t558 * t402 - t415 * t553;
t395 = m(4) * t455 + mrSges(4,2) * t538 + mrSges(4,3) * t513 + t508 * t526 + t512 * t572 + t569;
t389 = t549 * t394 + t547 * t395;
t404 = t557 * t409 + t552 * t410;
t570 = -t394 * t547 + t549 * t395;
t447 = Ifges(7,4) * t468 + Ifges(7,2) * t467 + Ifges(7,6) * t501;
t448 = Ifges(7,1) * t468 + Ifges(7,4) * t467 + Ifges(7,5) * t501;
t565 = -mrSges(7,1) * t423 + mrSges(7,2) * t424 - Ifges(7,5) * t439 - Ifges(7,6) * t438 - Ifges(7,3) * t471 - t468 * t447 + t467 * t448;
t564 = m(5) * t459 - t474 * mrSges(5,1) + mrSges(5,2) * t475 - t506 * t494 + t495 * t507 + t404;
t403 = m(4) * t488 - t513 * mrSges(4,1) + mrSges(4,2) * t514 - t526 * t511 + t512 * t527 + t564;
t446 = Ifges(7,5) * t468 + Ifges(7,6) * t467 + Ifges(7,3) * t501;
t412 = -mrSges(7,1) * t429 + mrSges(7,3) * t424 + Ifges(7,4) * t439 + Ifges(7,2) * t438 + Ifges(7,6) * t471 - t446 * t468 + t448 * t501;
t413 = mrSges(7,2) * t429 - mrSges(7,3) * t423 + Ifges(7,1) * t439 + Ifges(7,4) * t438 + Ifges(7,5) * t471 + t446 * t467 - t447 * t501;
t462 = Ifges(6,5) * t493 + Ifges(6,6) * t492 + Ifges(6,3) * t505;
t464 = Ifges(6,1) * t493 + Ifges(6,4) * t492 + Ifges(6,5) * t505;
t397 = -mrSges(6,1) * t431 + mrSges(6,3) * t428 + Ifges(6,4) * t458 + Ifges(6,2) * t457 + Ifges(6,6) * t473 - pkin(5) * t566 + pkin(11) * t568 + t556 * t412 + t551 * t413 - t493 * t462 + t505 * t464;
t463 = Ifges(6,4) * t493 + Ifges(6,2) * t492 + Ifges(6,6) * t505;
t398 = mrSges(6,2) * t431 - mrSges(6,3) * t427 + Ifges(6,1) * t458 + Ifges(6,4) * t457 + Ifges(6,5) * t473 - pkin(11) * t411 - t412 * t551 + t413 * t556 + t462 * t492 - t463 * t505;
t480 = Ifges(5,4) * t507 + Ifges(5,2) * t506 + Ifges(5,6) * t540;
t481 = Ifges(5,1) * t507 + Ifges(5,4) * t506 + Ifges(5,5) * t540;
t563 = mrSges(5,1) * t433 - mrSges(5,2) * t434 + Ifges(5,5) * t475 + Ifges(5,6) * t474 + Ifges(5,3) * t530 + pkin(4) * t421 + pkin(10) * t405 + t557 * t397 + t552 * t398 + t507 * t480 - t506 * t481;
t562 = mrSges(6,1) * t427 - mrSges(6,2) * t428 + Ifges(6,5) * t458 + Ifges(6,6) * t457 + Ifges(6,3) * t473 + pkin(5) * t411 + t493 * t463 - t492 * t464 - t565;
t536 = (-mrSges(3,1) * t559 + mrSges(3,2) * t554) * t576;
t532 = -mrSges(3,2) * t544 + mrSges(3,3) * t572;
t531 = mrSges(3,1) * t544 - mrSges(3,3) * t573;
t519 = -t533 * t548 - t581;
t518 = Ifges(3,5) * t544 + (Ifges(3,1) * t554 + Ifges(3,4) * t559) * t576;
t517 = Ifges(3,6) * t544 + (t554 * Ifges(3,4) + Ifges(3,2) * t559) * t576;
t516 = Ifges(3,3) * t544 + (Ifges(3,5) * t554 + Ifges(3,6) * t559) * t576;
t503 = -g(3) * t579 + t577;
t498 = Ifges(4,1) * t527 + Ifges(4,4) * t526 - Ifges(4,5) * t572;
t497 = Ifges(4,4) * t527 + Ifges(4,2) * t526 - Ifges(4,6) * t572;
t496 = Ifges(4,5) * t527 + Ifges(4,6) * t526 - Ifges(4,3) * t572;
t479 = Ifges(5,5) * t507 + Ifges(5,6) * t506 + Ifges(5,3) * t540;
t399 = m(3) * t502 + mrSges(3,1) * t543 - mrSges(3,3) * t537 + t532 * t544 - t536 * t573 - t403;
t390 = -mrSges(5,1) * t459 + mrSges(5,3) * t434 + Ifges(5,4) * t475 + Ifges(5,2) * t474 + Ifges(5,6) * t530 - pkin(4) * t404 - t507 * t479 + t540 * t481 - t562;
t388 = m(3) * t503 - mrSges(3,2) * t543 + mrSges(3,3) * t538 - t531 * t544 + t536 * t572 + t570;
t387 = mrSges(5,2) * t459 - mrSges(5,3) * t433 + Ifges(5,1) * t475 + Ifges(5,4) * t474 + Ifges(5,5) * t530 - pkin(10) * t404 - t397 * t552 + t398 * t557 + t479 * t506 - t480 * t540;
t386 = mrSges(4,2) * t488 - mrSges(4,3) * t454 + Ifges(4,1) * t514 + Ifges(4,4) * t513 - Ifges(4,5) * t538 - pkin(9) * t396 + t387 * t558 - t390 * t553 + t496 * t526 + t497 * t572;
t385 = -mrSges(4,1) * t488 + mrSges(4,3) * t455 + Ifges(4,4) * t514 + Ifges(4,2) * t513 - Ifges(4,6) * t538 - pkin(3) * t564 + pkin(9) * t569 + t553 * t387 + t558 * t390 - t527 * t496 - t498 * t572;
t384 = Ifges(3,5) * t537 + Ifges(3,6) * t538 + Ifges(3,3) * t543 + mrSges(3,1) * t502 - mrSges(3,2) * t503 + t547 * t386 + t549 * t385 - pkin(2) * t403 + qJ(3) * t570 + (t517 * t554 - t518 * t559) * t576;
t1 = [Ifges(2,3) * qJDD(1) + mrSges(2,1) * t571 - mrSges(2,2) * t567 + (mrSges(3,2) * t519 - mrSges(3,3) * t502 + Ifges(3,1) * t537 + Ifges(3,4) * t538 + Ifges(3,5) * t543 - qJ(3) * t389 - t385 * t547 + t549 * t386 + t516 * t572 - t544 * t517) * t579 + (-pkin(3) * t396 - t563 + (Ifges(3,2) + Ifges(4,3)) * t538 - pkin(2) * t389 + Ifges(3,6) * t543 + t544 * t518 + Ifges(3,4) * t537 + t526 * t498 - t527 * t497 - Ifges(4,5) * t514 - mrSges(3,1) * t519 - Ifges(4,6) * t513 - t516 * t573 - mrSges(4,1) * t454 + mrSges(4,2) * t455 + mrSges(3,3) * t503) * t578 + t550 * t384 + pkin(1) * ((t388 * t554 + t399 * t559) * t550 + (-m(3) * t519 + t538 * mrSges(3,1) - t537 * mrSges(3,2) + (-t531 * t554 + t532 * t559) * t576 - t389) * t548) + (t388 * t559 - t399 * t554) * t582; t384; t403; t563; t562; -t565;];
tauJ  = t1;
