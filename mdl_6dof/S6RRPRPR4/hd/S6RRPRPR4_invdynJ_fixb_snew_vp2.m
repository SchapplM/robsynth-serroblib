% Calculate vector of inverse dynamics joint torques for with Newton-Euler
% S6RRPRPR4
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
%   pkin=[a2,a3,a4,a5,a6,alpha2,d1,d2,d4,d6,theta3,theta5]';
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
% Datum: 2019-05-06 13:41
% Revision: 89c353f7eff3bd693eda4e29f35b2761dbc3ada0 (2019-05-03)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauJ = S6RRPRPR4_invdynJ_fixb_snew_vp2(qJ, qJD, qJDD, g, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(6,1),zeros(3,1),zeros(12,1),zeros(7,1),zeros(7,3),zeros(7,6)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRPRPR4_invdynJ_fixb_snew_vp2: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RRPRPR4_invdynJ_fixb_snew_vp2: qJD has to be [6x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [6 1]), ...
  'S6RRPRPR4_invdynJ_fixb_snew_vp2: qJDD has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6RRPRPR4_invdynJ_fixb_snew_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [12 1]), ...
  'S6RRPRPR4_invdynJ_fixb_snew_vp2: pkin has to be [12x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6RRPRPR4_invdynJ_fixb_snew_vp2: m has to be [7x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [7,3]), ...
  'S6RRPRPR4_invdynJ_fixb_snew_vp2: mrSges has to be [7x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [7 6]), ...
  'S6RRPRPR4_invdynJ_fixb_snew_vp2: Ifges has to be [7x6] (double)'); 

%% Symbolic Calculation
% From invdyn_fixb_NewtonEuler_linkframe_tauJ_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-05-06 13:32:44
% EndTime: 2019-05-06 13:33:01
% DurationCPUTime: 12.24s
% Computational Cost: add. (157033->361), mult. (417132->472), div. (0->0), fcn. (333362->14), ass. (0->151)
t588 = -2 * qJD(3);
t547 = sin(pkin(11));
t550 = cos(pkin(11));
t554 = sin(qJ(2));
t558 = cos(qJ(2));
t548 = sin(pkin(6));
t579 = qJD(1) * t548;
t529 = (t547 * t554 - t550 * t558) * t579;
t577 = qJD(1) * qJD(2);
t538 = (qJDD(1) * t554 + t558 * t577) * t548;
t551 = cos(pkin(6));
t541 = qJDD(1) * t551 + qJDD(2);
t542 = qJD(1) * t551 + qJD(2);
t560 = qJD(1) ^ 2;
t555 = sin(qJ(1));
t559 = cos(qJ(1));
t567 = -g(1) * t559 - g(2) * t555;
t585 = pkin(8) * t548;
t536 = -pkin(1) * t560 + qJDD(1) * t585 + t567;
t573 = t555 * g(1) - g(2) * t559;
t535 = qJDD(1) * pkin(1) + t560 * t585 + t573;
t583 = t535 * t551;
t568 = -t554 * t536 + t558 * t583;
t582 = t548 ^ 2 * t560;
t478 = t541 * pkin(2) - t538 * qJ(3) + (pkin(2) * t554 * t582 + (qJ(3) * qJD(1) * t542 - g(3)) * t548) * t558 + t568;
t581 = t548 * t554;
t506 = -g(3) * t581 + t558 * t536 + t554 * t583;
t575 = t554 * t579;
t532 = pkin(2) * t542 - qJ(3) * t575;
t539 = (qJDD(1) * t558 - t554 * t577) * t548;
t576 = t558 ^ 2 * t582;
t481 = -pkin(2) * t576 + qJ(3) * t539 - t532 * t542 + t506;
t530 = (t547 * t558 + t550 * t554) * t579;
t455 = t550 * t478 - t547 * t481 + t530 * t588;
t456 = t547 * t478 + t550 * t481 + t529 * t588;
t508 = pkin(3) * t529 - pkin(9) * t530;
t540 = t542 ^ 2;
t450 = -pkin(3) * t540 + pkin(9) * t541 - t508 * t529 + t456;
t521 = -t551 * g(3) - t548 * t535;
t492 = -t539 * pkin(2) - qJ(3) * t576 + t532 * t575 + qJDD(3) + t521;
t511 = -t538 * t547 + t539 * t550;
t512 = t538 * t550 + t539 * t547;
t459 = (t529 * t542 - t512) * pkin(9) + (t530 * t542 - t511) * pkin(3) + t492;
t553 = sin(qJ(4));
t557 = cos(qJ(4));
t442 = -t553 * t450 + t557 * t459;
t514 = -t530 * t553 + t542 * t557;
t489 = qJD(4) * t514 + t512 * t557 + t541 * t553;
t510 = qJDD(4) - t511;
t515 = t530 * t557 + t542 * t553;
t528 = qJD(4) + t529;
t439 = (t514 * t528 - t489) * qJ(5) + (t514 * t515 + t510) * pkin(4) + t442;
t443 = t557 * t450 + t553 * t459;
t488 = -qJD(4) * t515 - t512 * t553 + t541 * t557;
t499 = pkin(4) * t528 - qJ(5) * t515;
t513 = t514 ^ 2;
t441 = -pkin(4) * t513 + qJ(5) * t488 - t499 * t528 + t443;
t546 = sin(pkin(12));
t549 = cos(pkin(12));
t494 = t514 * t549 - t515 * t546;
t586 = 2 * qJD(5);
t436 = t546 * t439 + t549 * t441 + t494 * t586;
t495 = t514 * t546 + t515 * t549;
t472 = -pkin(5) * t494 - pkin(10) * t495;
t527 = t528 ^ 2;
t434 = -pkin(5) * t527 + pkin(10) * t510 + t472 * t494 + t436;
t449 = -t541 * pkin(3) - t540 * pkin(9) + t530 * t508 - t455;
t444 = -t488 * pkin(4) - t513 * qJ(5) + t515 * t499 + qJDD(5) + t449;
t463 = t488 * t549 - t489 * t546;
t464 = t488 * t546 + t489 * t549;
t437 = (-t494 * t528 - t464) * pkin(10) + (t495 * t528 - t463) * pkin(5) + t444;
t552 = sin(qJ(6));
t556 = cos(qJ(6));
t431 = -t434 * t552 + t437 * t556;
t474 = -t495 * t552 + t528 * t556;
t447 = qJD(6) * t474 + t464 * t556 + t510 * t552;
t475 = t495 * t556 + t528 * t552;
t460 = -mrSges(7,1) * t474 + mrSges(7,2) * t475;
t462 = qJDD(6) - t463;
t493 = qJD(6) - t494;
t465 = -mrSges(7,2) * t493 + mrSges(7,3) * t474;
t428 = m(7) * t431 + mrSges(7,1) * t462 - mrSges(7,3) * t447 - t460 * t475 + t465 * t493;
t432 = t434 * t556 + t437 * t552;
t446 = -qJD(6) * t475 - t464 * t552 + t510 * t556;
t466 = mrSges(7,1) * t493 - mrSges(7,3) * t475;
t429 = m(7) * t432 - mrSges(7,2) * t462 + mrSges(7,3) * t446 + t460 * t474 - t466 * t493;
t420 = -t428 * t552 + t556 * t429;
t471 = -mrSges(6,1) * t494 + mrSges(6,2) * t495;
t480 = mrSges(6,1) * t528 - mrSges(6,3) * t495;
t417 = m(6) * t436 - mrSges(6,2) * t510 + mrSges(6,3) * t463 + t471 * t494 - t480 * t528 + t420;
t566 = -t549 * t439 + t546 * t441;
t433 = -t510 * pkin(5) - t527 * pkin(10) + (t586 + t472) * t495 + t566;
t430 = -m(7) * t433 + t446 * mrSges(7,1) - mrSges(7,2) * t447 + t474 * t465 - t466 * t475;
t435 = -0.2e1 * qJD(5) * t495 - t566;
t479 = -mrSges(6,2) * t528 + mrSges(6,3) * t494;
t424 = m(6) * t435 + mrSges(6,1) * t510 - mrSges(6,3) * t464 - t471 * t495 + t479 * t528 + t430;
t412 = t546 * t417 + t549 * t424;
t451 = Ifges(7,5) * t475 + Ifges(7,6) * t474 + Ifges(7,3) * t493;
t453 = Ifges(7,1) * t475 + Ifges(7,4) * t474 + Ifges(7,5) * t493;
t421 = -mrSges(7,1) * t433 + mrSges(7,3) * t432 + Ifges(7,4) * t447 + Ifges(7,2) * t446 + Ifges(7,6) * t462 - t451 * t475 + t453 * t493;
t452 = Ifges(7,4) * t475 + Ifges(7,2) * t474 + Ifges(7,6) * t493;
t422 = mrSges(7,2) * t433 - mrSges(7,3) * t431 + Ifges(7,1) * t447 + Ifges(7,4) * t446 + Ifges(7,5) * t462 + t451 * t474 - t452 * t493;
t468 = Ifges(6,4) * t495 + Ifges(6,2) * t494 + Ifges(6,6) * t528;
t469 = Ifges(6,1) * t495 + Ifges(6,4) * t494 + Ifges(6,5) * t528;
t483 = Ifges(5,4) * t515 + Ifges(5,2) * t514 + Ifges(5,6) * t528;
t484 = Ifges(5,1) * t515 + Ifges(5,4) * t514 + Ifges(5,5) * t528;
t587 = Ifges(5,5) * t489 + Ifges(5,6) * t488 + t515 * t483 - t514 * t484 + mrSges(5,1) * t442 - mrSges(5,2) * t443 + Ifges(6,5) * t464 + Ifges(6,6) * t463 + t495 * t468 - t494 * t469 + mrSges(6,1) * t435 - mrSges(6,2) * t436 + t552 * t422 + t556 * t421 + pkin(5) * t430 + pkin(10) * t420 + pkin(4) * t412 + (Ifges(5,3) + Ifges(6,3)) * t510;
t580 = t548 * t558;
t507 = mrSges(4,1) * t529 + mrSges(4,2) * t530;
t517 = mrSges(4,1) * t542 - mrSges(4,3) * t530;
t496 = -mrSges(5,1) * t514 + mrSges(5,2) * t515;
t498 = -mrSges(5,2) * t528 + mrSges(5,3) * t514;
t410 = m(5) * t442 + mrSges(5,1) * t510 - mrSges(5,3) * t489 - t496 * t515 + t498 * t528 + t412;
t500 = mrSges(5,1) * t528 - mrSges(5,3) * t515;
t570 = t549 * t417 - t424 * t546;
t411 = m(5) * t443 - mrSges(5,2) * t510 + mrSges(5,3) * t488 + t496 * t514 - t500 * t528 + t570;
t571 = -t410 * t553 + t557 * t411;
t402 = m(4) * t456 - mrSges(4,2) * t541 + mrSges(4,3) * t511 - t507 * t529 - t517 * t542 + t571;
t516 = -mrSges(4,2) * t542 - mrSges(4,3) * t529;
t419 = t556 * t428 + t552 * t429;
t418 = m(6) * t444 - t463 * mrSges(6,1) + mrSges(6,2) * t464 - t494 * t479 + t480 * t495 + t419;
t562 = -m(5) * t449 + t488 * mrSges(5,1) - mrSges(5,2) * t489 + t514 * t498 - t500 * t515 - t418;
t414 = m(4) * t455 + mrSges(4,1) * t541 - mrSges(4,3) * t512 - t507 * t530 + t516 * t542 + t562;
t399 = t547 * t402 + t550 * t414;
t404 = t557 * t410 + t553 * t411;
t574 = t558 * t579;
t572 = t550 * t402 - t414 * t547;
t564 = -m(4) * t492 + t511 * mrSges(4,1) - t512 * mrSges(4,2) - t529 * t516 - t530 * t517 - t404;
t563 = mrSges(7,1) * t431 - mrSges(7,2) * t432 + Ifges(7,5) * t447 + Ifges(7,6) * t446 + Ifges(7,3) * t462 + t452 * t475 - t453 * t474;
t537 = (-mrSges(3,1) * t558 + mrSges(3,2) * t554) * t579;
t534 = -mrSges(3,2) * t542 + mrSges(3,3) * t574;
t533 = mrSges(3,1) * t542 - mrSges(3,3) * t575;
t520 = Ifges(3,5) * t542 + (Ifges(3,1) * t554 + Ifges(3,4) * t558) * t579;
t519 = Ifges(3,6) * t542 + (Ifges(3,4) * t554 + Ifges(3,2) * t558) * t579;
t518 = Ifges(3,3) * t542 + (Ifges(3,5) * t554 + Ifges(3,6) * t558) * t579;
t505 = -g(3) * t580 + t568;
t503 = Ifges(4,1) * t530 - Ifges(4,4) * t529 + Ifges(4,5) * t542;
t502 = Ifges(4,4) * t530 - Ifges(4,2) * t529 + Ifges(4,6) * t542;
t501 = Ifges(4,5) * t530 - Ifges(4,6) * t529 + Ifges(4,3) * t542;
t482 = Ifges(5,5) * t515 + Ifges(5,6) * t514 + Ifges(5,3) * t528;
t467 = Ifges(6,5) * t495 + Ifges(6,6) * t494 + Ifges(6,3) * t528;
t406 = -mrSges(6,1) * t444 + mrSges(6,3) * t436 + Ifges(6,4) * t464 + Ifges(6,2) * t463 + Ifges(6,6) * t510 - pkin(5) * t419 - t467 * t495 + t469 * t528 - t563;
t405 = mrSges(6,2) * t444 - mrSges(6,3) * t435 + Ifges(6,1) * t464 + Ifges(6,4) * t463 + Ifges(6,5) * t510 - pkin(10) * t419 - t421 * t552 + t422 * t556 + t467 * t494 - t468 * t528;
t398 = m(3) * t506 - mrSges(3,2) * t541 + mrSges(3,3) * t539 - t533 * t542 + t537 * t574 + t572;
t397 = m(3) * t505 + mrSges(3,1) * t541 - mrSges(3,3) * t538 + t534 * t542 - t537 * t575 + t399;
t396 = mrSges(5,2) * t449 - mrSges(5,3) * t442 + Ifges(5,1) * t489 + Ifges(5,4) * t488 + Ifges(5,5) * t510 - qJ(5) * t412 + t405 * t549 - t406 * t546 + t482 * t514 - t483 * t528;
t395 = -mrSges(5,1) * t449 + mrSges(5,3) * t443 + Ifges(5,4) * t489 + Ifges(5,2) * t488 + Ifges(5,6) * t510 - pkin(4) * t418 + qJ(5) * t570 + t546 * t405 + t549 * t406 - t515 * t482 + t528 * t484;
t394 = -mrSges(4,1) * t492 + mrSges(4,3) * t456 + Ifges(4,4) * t512 + Ifges(4,2) * t511 + Ifges(4,6) * t541 - pkin(3) * t404 - t530 * t501 + t542 * t503 - t587;
t393 = mrSges(4,2) * t492 - mrSges(4,3) * t455 + Ifges(4,1) * t512 + Ifges(4,4) * t511 + Ifges(4,5) * t541 - pkin(9) * t404 - t395 * t553 + t396 * t557 - t501 * t529 - t502 * t542;
t392 = Ifges(3,5) * t538 + Ifges(3,6) * t539 + mrSges(3,1) * t505 - mrSges(3,2) * t506 + Ifges(4,5) * t512 + Ifges(4,6) * t511 + t530 * t502 + t529 * t503 + mrSges(4,1) * t455 - mrSges(4,2) * t456 + t553 * t396 + t557 * t395 + pkin(3) * t562 + pkin(9) * t571 + pkin(2) * t399 + (Ifges(3,3) + Ifges(4,3)) * t541 + (t519 * t554 - t520 * t558) * t579;
t1 = [Ifges(2,3) * qJDD(1) + mrSges(2,1) * t573 - mrSges(2,2) * t567 + (mrSges(3,2) * t521 - mrSges(3,3) * t505 + Ifges(3,1) * t538 + Ifges(3,4) * t539 + Ifges(3,5) * t541 - qJ(3) * t399 + t393 * t550 - t394 * t547 + t518 * t574 - t519 * t542) * t581 + (-mrSges(3,1) * t521 + mrSges(3,3) * t506 + Ifges(3,4) * t538 + Ifges(3,2) * t539 + Ifges(3,6) * t541 + pkin(2) * t564 + qJ(3) * t572 + t547 * t393 + t550 * t394 - t518 * t575 + t542 * t520) * t580 + t551 * t392 + pkin(1) * ((t397 * t558 + t398 * t554) * t551 + (-m(3) * t521 + t539 * mrSges(3,1) - t538 * mrSges(3,2) + (-t533 * t554 + t534 * t558) * t579 + t564) * t548) + (-t397 * t554 + t398 * t558) * t585; t392; -t564; t587; t418; t563;];
tauJ  = t1;
