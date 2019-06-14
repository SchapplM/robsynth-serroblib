% Calculate vector of inverse dynamics joint torques for with Newton-Euler
% S6RRRPRR12
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
%   pkin=[a2,a3,a4,a5,a6,alpha2,d1,d2,d3,d5,d6,theta4]';
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
% Datum: 2019-05-07 15:13
% Revision: 89c353f7eff3bd693eda4e29f35b2761dbc3ada0 (2019-05-03)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauJ = S6RRRPRR12_invdynJ_fixb_snew_vp2(qJ, qJD, qJDD, g, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(6,1),zeros(3,1),zeros(12,1),zeros(7,1),zeros(7,3),zeros(7,6)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRRPRR12_invdynJ_fixb_snew_vp2: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RRRPRR12_invdynJ_fixb_snew_vp2: qJD has to be [6x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [6 1]), ...
  'S6RRRPRR12_invdynJ_fixb_snew_vp2: qJDD has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6RRRPRR12_invdynJ_fixb_snew_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [12 1]), ...
  'S6RRRPRR12_invdynJ_fixb_snew_vp2: pkin has to be [12x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6RRRPRR12_invdynJ_fixb_snew_vp2: m has to be [7x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [7,3]), ...
  'S6RRRPRR12_invdynJ_fixb_snew_vp2: mrSges has to be [7x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [7 6]), ...
  'S6RRRPRR12_invdynJ_fixb_snew_vp2: Ifges has to be [7x6] (double)'); 

%% Symbolic Calculation
% From invdyn_fixb_NewtonEuler_linkframe_tauJ_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-05-07 14:59:41
% EndTime: 2019-05-07 15:00:06
% DurationCPUTime: 18.98s
% Computational Cost: add. (300453->362), mult. (656223->472), div. (0->0), fcn. (530881->14), ass. (0->152)
t535 = sin(pkin(6));
t541 = sin(qJ(2));
t545 = cos(qJ(2));
t563 = qJD(1) * qJD(2);
t524 = (-qJDD(1) * t545 + t541 * t563) * t535;
t572 = cos(qJ(3));
t571 = pkin(8) * t535;
t537 = cos(pkin(6));
t570 = t537 * g(3);
t547 = qJD(1) ^ 2;
t542 = sin(qJ(1));
t546 = cos(qJ(1));
t560 = g(1) * t542 - g(2) * t546;
t519 = qJDD(1) * pkin(1) + t547 * t571 + t560;
t569 = t519 * t537;
t568 = t535 * t541;
t567 = t535 * t545;
t565 = qJD(1) * t535;
t522 = (-pkin(2) * t545 - pkin(9) * t541) * t565;
t531 = qJD(1) * t537 + qJD(2);
t529 = t531 ^ 2;
t530 = qJDD(1) * t537 + qJDD(2);
t564 = qJD(1) * t545;
t556 = -g(1) * t546 - g(2) * t542;
t520 = -pkin(1) * t547 + qJDD(1) * t571 + t556;
t566 = t520 * t545 + t541 * t569;
t471 = -t529 * pkin(2) + t530 * pkin(9) + (-g(3) * t541 + t522 * t564) * t535 + t566;
t523 = (qJDD(1) * t541 + t545 * t563) * t535;
t472 = t524 * pkin(2) - t523 * pkin(9) - t570 + (-t519 + (pkin(2) * t541 - pkin(9) * t545) * t531 * qJD(1)) * t535;
t540 = sin(qJ(3));
t451 = t471 * t572 + t472 * t540;
t562 = t541 * t565;
t512 = -t531 * t572 + t540 * t562;
t513 = t531 * t540 + t562 * t572;
t495 = pkin(3) * t512 - qJ(4) * t513;
t516 = qJDD(3) + t524;
t561 = t535 * t564;
t528 = qJD(3) - t561;
t527 = t528 ^ 2;
t441 = -pkin(3) * t527 + qJ(4) * t516 - t495 * t512 + t451;
t493 = -g(3) * t567 - t520 * t541 + t545 * t569;
t470 = -t530 * pkin(2) - t529 * pkin(9) + t522 * t562 - t493;
t491 = qJD(3) * t513 + t523 * t540 - t530 * t572;
t492 = -qJD(3) * t512 + t523 * t572 + t530 * t540;
t444 = (t512 * t528 - t492) * qJ(4) + (t513 * t528 + t491) * pkin(3) + t470;
t534 = sin(pkin(12));
t536 = cos(pkin(12));
t501 = t513 * t536 + t528 * t534;
t430 = -0.2e1 * qJD(4) * t501 - t534 * t441 + t444 * t536;
t478 = t492 * t536 + t516 * t534;
t500 = -t513 * t534 + t528 * t536;
t422 = (t500 * t512 - t478) * pkin(10) + (t500 * t501 + t491) * pkin(4) + t430;
t431 = 0.2e1 * qJD(4) * t500 + t441 * t536 + t444 * t534;
t477 = -t492 * t534 + t516 * t536;
t482 = pkin(4) * t512 - pkin(10) * t501;
t499 = t500 ^ 2;
t424 = -pkin(4) * t499 + pkin(10) * t477 - t482 * t512 + t431;
t539 = sin(qJ(5));
t544 = cos(qJ(5));
t417 = t422 * t544 - t539 * t424;
t474 = t500 * t544 - t501 * t539;
t447 = qJD(5) * t474 + t477 * t539 + t478 * t544;
t475 = t500 * t539 + t501 * t544;
t489 = qJDD(5) + t491;
t511 = qJD(5) + t512;
t415 = (t474 * t511 - t447) * pkin(11) + (t474 * t475 + t489) * pkin(5) + t417;
t418 = t422 * t539 + t424 * t544;
t446 = -qJD(5) * t475 + t477 * t544 - t478 * t539;
t461 = pkin(5) * t511 - pkin(11) * t475;
t473 = t474 ^ 2;
t416 = -pkin(5) * t473 + pkin(11) * t446 - t461 * t511 + t418;
t538 = sin(qJ(6));
t543 = cos(qJ(6));
t413 = t415 * t543 - t416 * t538;
t456 = t474 * t543 - t475 * t538;
t429 = qJD(6) * t456 + t446 * t538 + t447 * t543;
t457 = t474 * t538 + t475 * t543;
t438 = -mrSges(7,1) * t456 + mrSges(7,2) * t457;
t510 = qJD(6) + t511;
t448 = -mrSges(7,2) * t510 + mrSges(7,3) * t456;
t484 = qJDD(6) + t489;
t407 = m(7) * t413 + mrSges(7,1) * t484 - mrSges(7,3) * t429 - t438 * t457 + t448 * t510;
t414 = t415 * t538 + t416 * t543;
t428 = -qJD(6) * t457 + t446 * t543 - t447 * t538;
t449 = mrSges(7,1) * t510 - mrSges(7,3) * t457;
t408 = m(7) * t414 - mrSges(7,2) * t484 + mrSges(7,3) * t428 + t438 * t456 - t449 * t510;
t401 = t407 * t543 + t408 * t538;
t458 = -mrSges(6,1) * t474 + mrSges(6,2) * t475;
t459 = -mrSges(6,2) * t511 + mrSges(6,3) * t474;
t399 = m(6) * t417 + mrSges(6,1) * t489 - mrSges(6,3) * t447 - t458 * t475 + t459 * t511 + t401;
t460 = mrSges(6,1) * t511 - mrSges(6,3) * t475;
t557 = -t407 * t538 + t408 * t543;
t400 = m(6) * t418 - mrSges(6,2) * t489 + mrSges(6,3) * t446 + t458 * t474 - t460 * t511 + t557;
t395 = t399 * t544 + t400 * t539;
t479 = -mrSges(5,1) * t500 + mrSges(5,2) * t501;
t555 = -mrSges(5,2) * t512 + mrSges(5,3) * t500;
t393 = m(5) * t430 + t491 * mrSges(5,1) - t478 * mrSges(5,3) - t501 * t479 + t512 * t555 + t395;
t481 = mrSges(5,1) * t512 - mrSges(5,3) * t501;
t558 = -t399 * t539 + t400 * t544;
t394 = m(5) * t431 - mrSges(5,2) * t491 + mrSges(5,3) * t477 + t479 * t500 - t481 * t512 + t558;
t389 = -t393 * t534 + t394 * t536;
t496 = mrSges(4,1) * t512 + mrSges(4,2) * t513;
t503 = mrSges(4,1) * t528 - mrSges(4,3) * t513;
t387 = m(4) * t451 - mrSges(4,2) * t516 - mrSges(4,3) * t491 - t496 * t512 - t503 * t528 + t389;
t450 = -t471 * t540 + t472 * t572;
t440 = -pkin(3) * t516 - qJ(4) * t527 + t495 * t513 + qJDD(4) - t450;
t433 = -pkin(4) * t477 - pkin(10) * t499 + t482 * t501 + t440;
t419 = -pkin(5) * t446 - pkin(11) * t473 + t461 * t475 + t433;
t553 = m(7) * t419 - mrSges(7,1) * t428 + mrSges(7,2) * t429 - t448 * t456 + t449 * t457;
t550 = m(6) * t433 - mrSges(6,1) * t446 + mrSges(6,2) * t447 - t459 * t474 + t460 * t475 + t553;
t411 = m(5) * t440 - mrSges(5,1) * t477 + mrSges(5,2) * t478 + t481 * t501 - t500 * t555 + t550;
t502 = -mrSges(4,2) * t528 - mrSges(4,3) * t512;
t410 = m(4) * t450 + mrSges(4,1) * t516 - mrSges(4,3) * t492 - t496 * t513 + t502 * t528 - t411;
t383 = t387 * t540 + t410 * t572;
t559 = t387 * t572 - t410 * t540;
t388 = t393 * t536 + t394 * t534;
t435 = Ifges(7,4) * t457 + Ifges(7,2) * t456 + Ifges(7,6) * t510;
t436 = Ifges(7,1) * t457 + Ifges(7,4) * t456 + Ifges(7,5) * t510;
t552 = -mrSges(7,1) * t413 + mrSges(7,2) * t414 - Ifges(7,5) * t429 - Ifges(7,6) * t428 - Ifges(7,3) * t484 - t435 * t457 + t456 * t436;
t551 = -m(4) * t470 - mrSges(4,1) * t491 - mrSges(4,2) * t492 - t502 * t512 - t503 * t513 - t388;
t434 = Ifges(7,5) * t457 + Ifges(7,6) * t456 + Ifges(7,3) * t510;
t402 = -mrSges(7,1) * t419 + mrSges(7,3) * t414 + Ifges(7,4) * t429 + Ifges(7,2) * t428 + Ifges(7,6) * t484 - t434 * t457 + t436 * t510;
t403 = mrSges(7,2) * t419 - mrSges(7,3) * t413 + Ifges(7,1) * t429 + Ifges(7,4) * t428 + Ifges(7,5) * t484 + t434 * t456 - t435 * t510;
t452 = Ifges(6,5) * t475 + Ifges(6,6) * t474 + Ifges(6,3) * t511;
t454 = Ifges(6,1) * t475 + Ifges(6,4) * t474 + Ifges(6,5) * t511;
t390 = -mrSges(6,1) * t433 + mrSges(6,3) * t418 + Ifges(6,4) * t447 + Ifges(6,2) * t446 + Ifges(6,6) * t489 - pkin(5) * t553 + pkin(11) * t557 + t543 * t402 + t538 * t403 - t475 * t452 + t511 * t454;
t453 = Ifges(6,4) * t475 + Ifges(6,2) * t474 + Ifges(6,6) * t511;
t391 = mrSges(6,2) * t433 - mrSges(6,3) * t417 + Ifges(6,1) * t447 + Ifges(6,4) * t446 + Ifges(6,5) * t489 - pkin(11) * t401 - t402 * t538 + t403 * t543 + t452 * t474 - t453 * t511;
t462 = Ifges(5,5) * t501 + Ifges(5,6) * t500 + Ifges(5,3) * t512;
t464 = Ifges(5,1) * t501 + Ifges(5,4) * t500 + Ifges(5,5) * t512;
t380 = -mrSges(5,1) * t440 + mrSges(5,3) * t431 + Ifges(5,4) * t478 + Ifges(5,2) * t477 + Ifges(5,6) * t491 - pkin(4) * t550 + pkin(10) * t558 + t544 * t390 + t539 * t391 - t501 * t462 + t512 * t464;
t463 = Ifges(5,4) * t501 + Ifges(5,2) * t500 + Ifges(5,6) * t512;
t381 = mrSges(5,2) * t440 - mrSges(5,3) * t430 + Ifges(5,1) * t478 + Ifges(5,4) * t477 + Ifges(5,5) * t491 - pkin(10) * t395 - t390 * t539 + t391 * t544 + t462 * t500 - t463 * t512;
t486 = Ifges(4,4) * t513 - Ifges(4,2) * t512 + Ifges(4,6) * t528;
t487 = Ifges(4,1) * t513 - Ifges(4,4) * t512 + Ifges(4,5) * t528;
t549 = mrSges(4,1) * t450 - mrSges(4,2) * t451 + Ifges(4,5) * t492 - Ifges(4,6) * t491 + Ifges(4,3) * t516 - pkin(3) * t411 + qJ(4) * t389 + t536 * t380 + t534 * t381 + t513 * t486 + t512 * t487;
t548 = mrSges(6,1) * t417 - mrSges(6,2) * t418 + Ifges(6,5) * t447 + Ifges(6,6) * t446 + Ifges(6,3) * t489 + pkin(5) * t401 + t475 * t453 - t474 * t454 - t552;
t521 = (-mrSges(3,1) * t545 + mrSges(3,2) * t541) * t565;
t518 = -mrSges(3,2) * t531 + mrSges(3,3) * t561;
t517 = mrSges(3,1) * t531 - mrSges(3,3) * t562;
t507 = -t535 * t519 - t570;
t506 = Ifges(3,5) * t531 + (Ifges(3,1) * t541 + Ifges(3,4) * t545) * t565;
t505 = Ifges(3,6) * t531 + (Ifges(3,4) * t541 + Ifges(3,2) * t545) * t565;
t504 = Ifges(3,3) * t531 + (Ifges(3,5) * t541 + Ifges(3,6) * t545) * t565;
t494 = -g(3) * t568 + t566;
t485 = Ifges(4,5) * t513 - Ifges(4,6) * t512 + Ifges(4,3) * t528;
t384 = m(3) * t493 + mrSges(3,1) * t530 - mrSges(3,3) * t523 + t518 * t531 - t521 * t562 + t551;
t382 = m(3) * t494 - mrSges(3,2) * t530 - mrSges(3,3) * t524 - t517 * t531 + t521 * t561 + t559;
t379 = -pkin(4) * t395 - t548 - pkin(3) * t388 + (-Ifges(5,3) - Ifges(4,2)) * t491 + t528 * t487 + Ifges(4,6) * t516 - t513 * t485 - mrSges(5,1) * t430 + mrSges(5,2) * t431 + mrSges(4,3) * t451 - mrSges(4,1) * t470 - Ifges(5,6) * t477 - Ifges(5,5) * t478 + Ifges(4,4) * t492 + t500 * t464 - t501 * t463;
t378 = mrSges(4,2) * t470 - mrSges(4,3) * t450 + Ifges(4,1) * t492 - Ifges(4,4) * t491 + Ifges(4,5) * t516 - qJ(4) * t388 - t380 * t534 + t381 * t536 - t485 * t512 - t486 * t528;
t377 = Ifges(3,5) * t523 - Ifges(3,6) * t524 + Ifges(3,3) * t530 + mrSges(3,1) * t493 - mrSges(3,2) * t494 + t540 * t378 + t572 * t379 + pkin(2) * t551 + pkin(9) * t559 + (t505 * t541 - t506 * t545) * t565;
t1 = [Ifges(2,3) * qJDD(1) + mrSges(2,1) * t560 - mrSges(2,2) * t556 + (mrSges(3,2) * t507 - mrSges(3,3) * t493 + Ifges(3,1) * t523 - Ifges(3,4) * t524 + Ifges(3,5) * t530 - pkin(9) * t383 + t378 * t572 - t379 * t540 + t504 * t561 - t505 * t531) * t568 + (-mrSges(3,1) * t507 + mrSges(3,3) * t494 + Ifges(3,4) * t523 - Ifges(3,2) * t524 + Ifges(3,6) * t530 - pkin(2) * t383 - t504 * t562 + t531 * t506 - t549) * t567 + t537 * t377 + pkin(1) * ((t382 * t541 + t384 * t545) * t537 + (-m(3) * t507 - t524 * mrSges(3,1) - t523 * mrSges(3,2) + (-t517 * t541 + t518 * t545) * t565 - t383) * t535) + (t382 * t545 - t384 * t541) * t571; t377; t549; t411; t548; -t552;];
tauJ  = t1;
