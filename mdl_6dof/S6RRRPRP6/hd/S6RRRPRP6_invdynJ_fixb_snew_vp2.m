% Calculate vector of inverse dynamics joint torques for with Newton-Euler
% S6RRRPRP6
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
%   pkin=[a2,a3,a4,a5,a6,alpha2,d1,d2,d3,d5,theta4]';
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
% Datum: 2019-05-07 08:11
% Revision: 89c353f7eff3bd693eda4e29f35b2761dbc3ada0 (2019-05-03)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauJ = S6RRRPRP6_invdynJ_fixb_snew_vp2(qJ, qJD, qJDD, g, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(6,1),zeros(3,1),zeros(11,1),zeros(7,1),zeros(7,3),zeros(7,6)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRRPRP6_invdynJ_fixb_snew_vp2: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RRRPRP6_invdynJ_fixb_snew_vp2: qJD has to be [6x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [6 1]), ...
  'S6RRRPRP6_invdynJ_fixb_snew_vp2: qJDD has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6RRRPRP6_invdynJ_fixb_snew_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [11 1]), ...
  'S6RRRPRP6_invdynJ_fixb_snew_vp2: pkin has to be [11x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6RRRPRP6_invdynJ_fixb_snew_vp2: m has to be [7x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [7,3]), ...
  'S6RRRPRP6_invdynJ_fixb_snew_vp2: mrSges has to be [7x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [7 6]), ...
  'S6RRRPRP6_invdynJ_fixb_snew_vp2: Ifges has to be [7x6] (double)'); 

%% Symbolic Calculation
% From invdyn_fixb_NewtonEuler_linkframe_tauJ_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-05-07 08:01:18
% EndTime: 2019-05-07 08:01:30
% DurationCPUTime: 8.22s
% Computational Cost: add. (112719->339), mult. (248547->430), div. (0->0), fcn. (198870->12), ass. (0->143)
t585 = Ifges(6,1) + Ifges(7,1);
t577 = Ifges(6,4) + Ifges(7,4);
t576 = Ifges(6,5) + Ifges(7,5);
t584 = Ifges(6,2) + Ifges(7,2);
t575 = Ifges(6,6) + Ifges(7,6);
t583 = Ifges(6,3) + Ifges(7,3);
t535 = sin(pkin(6));
t540 = sin(qJ(2));
t544 = cos(qJ(2));
t561 = qJD(1) * qJD(2);
t526 = (-qJDD(1) * t544 + t540 * t561) * t535;
t564 = qJD(1) * t535;
t524 = (-pkin(2) * t544 - pkin(9) * t540) * t564;
t537 = cos(pkin(6));
t531 = qJD(1) * t537 + qJD(2);
t529 = t531 ^ 2;
t530 = qJDD(1) * t537 + qJDD(2);
t563 = qJD(1) * t544;
t546 = qJD(1) ^ 2;
t541 = sin(qJ(1));
t545 = cos(qJ(1));
t551 = -g(1) * t545 - g(2) * t541;
t580 = pkin(8) * t535;
t522 = -pkin(1) * t546 + qJDD(1) * t580 + t551;
t556 = t541 * g(1) - g(2) * t545;
t521 = qJDD(1) * pkin(1) + t546 * t580 + t556;
t572 = t521 * t537;
t565 = t544 * t522 + t540 * t572;
t481 = -pkin(2) * t529 + pkin(9) * t530 + (-g(3) * t540 + t524 * t563) * t535 + t565;
t525 = (qJDD(1) * t540 + t544 * t561) * t535;
t579 = g(3) * t537;
t482 = pkin(2) * t526 - pkin(9) * t525 - t579 + (-t521 + (pkin(2) * t540 - pkin(9) * t544) * t531 * qJD(1)) * t535;
t539 = sin(qJ(3));
t543 = cos(qJ(3));
t447 = -t481 * t539 + t543 * t482;
t558 = t540 * t564;
t514 = t531 * t543 - t539 * t558;
t495 = qJD(3) * t514 + t525 * t543 + t530 * t539;
t515 = t531 * t539 + t543 * t558;
t518 = qJDD(3) + t526;
t557 = t535 * t563;
t528 = qJD(3) - t557;
t436 = (t514 * t528 - t495) * qJ(4) + (t514 * t515 + t518) * pkin(3) + t447;
t448 = t543 * t481 + t539 * t482;
t494 = -qJD(3) * t515 - t525 * t539 + t530 * t543;
t505 = pkin(3) * t528 - qJ(4) * t515;
t513 = t514 ^ 2;
t439 = -pkin(3) * t513 + qJ(4) * t494 - t505 * t528 + t448;
t534 = sin(pkin(11));
t536 = cos(pkin(11));
t502 = t514 * t534 + t515 * t536;
t430 = -0.2e1 * qJD(4) * t502 + t436 * t536 - t534 * t439;
t501 = t514 * t536 - t515 * t534;
t431 = 0.2e1 * qJD(4) * t501 + t534 * t436 + t536 * t439;
t476 = -pkin(4) * t501 - pkin(10) * t502;
t527 = t528 ^ 2;
t429 = -pkin(4) * t527 + pkin(10) * t518 + t476 * t501 + t431;
t570 = t535 * t544;
t496 = -g(3) * t570 - t540 * t522 + t544 * t572;
t480 = -pkin(2) * t530 - pkin(9) * t529 + t524 * t558 - t496;
t440 = -pkin(3) * t494 - qJ(4) * t513 + t515 * t505 + qJDD(4) + t480;
t469 = t494 * t536 - t495 * t534;
t470 = t494 * t534 + t495 * t536;
t434 = (-t501 * t528 - t470) * pkin(10) + (t502 * t528 - t469) * pkin(4) + t440;
t538 = sin(qJ(5));
t542 = cos(qJ(5));
t424 = -t429 * t538 + t542 * t434;
t484 = -t502 * t538 + t528 * t542;
t445 = qJD(5) * t484 + t470 * t542 + t518 * t538;
t485 = t502 * t542 + t528 * t538;
t459 = -mrSges(7,1) * t484 + mrSges(7,2) * t485;
t460 = -mrSges(6,1) * t484 + mrSges(6,2) * t485;
t500 = qJD(5) - t501;
t462 = -mrSges(6,2) * t500 + mrSges(6,3) * t484;
t468 = qJDD(5) - t469;
t421 = -0.2e1 * qJD(6) * t485 + (t484 * t500 - t445) * qJ(6) + (t484 * t485 + t468) * pkin(5) + t424;
t461 = -mrSges(7,2) * t500 + mrSges(7,3) * t484;
t560 = m(7) * t421 + t468 * mrSges(7,1) + t500 * t461;
t411 = m(6) * t424 + mrSges(6,1) * t468 + t462 * t500 + (-t459 - t460) * t485 + (-mrSges(6,3) - mrSges(7,3)) * t445 + t560;
t425 = t542 * t429 + t538 * t434;
t444 = -qJD(5) * t485 - t470 * t538 + t518 * t542;
t463 = pkin(5) * t500 - qJ(6) * t485;
t483 = t484 ^ 2;
t423 = -pkin(5) * t483 + qJ(6) * t444 + 0.2e1 * qJD(6) * t484 - t463 * t500 + t425;
t559 = m(7) * t423 + t444 * mrSges(7,3) + t484 * t459;
t464 = mrSges(7,1) * t500 - mrSges(7,3) * t485;
t566 = -mrSges(6,1) * t500 + mrSges(6,3) * t485 - t464;
t578 = -mrSges(6,2) - mrSges(7,2);
t416 = m(6) * t425 + mrSges(6,3) * t444 + t460 * t484 + t468 * t578 + t500 * t566 + t559;
t409 = -t411 * t538 + t542 * t416;
t475 = -mrSges(5,1) * t501 + mrSges(5,2) * t502;
t487 = mrSges(5,1) * t528 - mrSges(5,3) * t502;
t405 = m(5) * t431 - mrSges(5,2) * t518 + mrSges(5,3) * t469 + t475 * t501 - t487 * t528 + t409;
t428 = -pkin(4) * t518 - pkin(10) * t527 + t502 * t476 - t430;
t426 = -pkin(5) * t444 - qJ(6) * t483 + t463 * t485 + qJDD(6) + t428;
t552 = -m(7) * t426 + t444 * mrSges(7,1) + t484 * t461;
t417 = -m(6) * t428 + t444 * mrSges(6,1) + t445 * t578 + t484 * t462 + t485 * t566 + t552;
t486 = -mrSges(5,2) * t528 + mrSges(5,3) * t501;
t413 = m(5) * t430 + mrSges(5,1) * t518 - mrSges(5,3) * t470 - t475 * t502 + t486 * t528 + t417;
t400 = t534 * t405 + t536 * t413;
t419 = mrSges(7,2) * t445 + t464 * t485 - t552;
t567 = t577 * t484 + t485 * t585 + t576 * t500;
t569 = -t484 * t575 - t485 * t576 - t500 * t583;
t401 = -mrSges(6,1) * t428 + mrSges(6,3) * t425 - mrSges(7,1) * t426 + mrSges(7,3) * t423 - pkin(5) * t419 + qJ(6) * t559 + (-qJ(6) * t464 + t567) * t500 + t569 * t485 + (-qJ(6) * mrSges(7,2) + t575) * t468 + t577 * t445 + t584 * t444;
t418 = -mrSges(7,3) * t445 - t459 * t485 + t560;
t568 = -t484 * t584 - t485 * t577 - t500 * t575;
t406 = mrSges(6,2) * t428 + mrSges(7,2) * t426 - mrSges(6,3) * t424 - mrSges(7,3) * t421 - qJ(6) * t418 + t577 * t444 + t445 * t585 + t576 * t468 - t569 * t484 + t568 * t500;
t472 = Ifges(5,4) * t502 + Ifges(5,2) * t501 + Ifges(5,6) * t528;
t473 = Ifges(5,1) * t502 + Ifges(5,4) * t501 + Ifges(5,5) * t528;
t489 = Ifges(4,4) * t515 + Ifges(4,2) * t514 + Ifges(4,6) * t528;
t490 = Ifges(4,1) * t515 + Ifges(4,4) * t514 + Ifges(4,5) * t528;
t582 = Ifges(4,5) * t495 + Ifges(4,6) * t494 + t515 * t489 - t514 * t490 + mrSges(4,1) * t447 - mrSges(4,2) * t448 + Ifges(5,5) * t470 + Ifges(5,6) * t469 + t502 * t472 - t501 * t473 + mrSges(5,1) * t430 - mrSges(5,2) * t431 + t538 * t406 + t542 * t401 + pkin(4) * t417 + pkin(10) * t409 + pkin(3) * t400 + (Ifges(4,3) + Ifges(5,3)) * t518;
t581 = mrSges(6,1) * t424 + mrSges(7,1) * t421 - mrSges(6,2) * t425 - mrSges(7,2) * t423 + pkin(5) * t418 + t444 * t575 + t445 * t576 + t468 * t583 - t484 * t567 - t485 * t568;
t571 = t535 * t540;
t503 = -mrSges(4,1) * t514 + mrSges(4,2) * t515;
t504 = -mrSges(4,2) * t528 + mrSges(4,3) * t514;
t398 = m(4) * t447 + mrSges(4,1) * t518 - mrSges(4,3) * t495 - t503 * t515 + t504 * t528 + t400;
t506 = mrSges(4,1) * t528 - mrSges(4,3) * t515;
t554 = t536 * t405 - t413 * t534;
t399 = m(4) * t448 - mrSges(4,2) * t518 + mrSges(4,3) * t494 + t503 * t514 - t506 * t528 + t554;
t393 = t543 * t398 + t539 * t399;
t408 = t542 * t411 + t538 * t416;
t555 = -t398 * t539 + t543 * t399;
t407 = m(5) * t440 - t469 * mrSges(5,1) + mrSges(5,2) * t470 - t501 * t486 + t487 * t502 + t408;
t548 = -m(4) * t480 + t494 * mrSges(4,1) - mrSges(4,2) * t495 + t514 * t504 - t506 * t515 - t407;
t523 = (-mrSges(3,1) * t544 + mrSges(3,2) * t540) * t564;
t520 = -mrSges(3,2) * t531 + mrSges(3,3) * t557;
t519 = mrSges(3,1) * t531 - mrSges(3,3) * t558;
t510 = -t521 * t535 - t579;
t509 = Ifges(3,5) * t531 + (Ifges(3,1) * t540 + Ifges(3,4) * t544) * t564;
t508 = Ifges(3,6) * t531 + (Ifges(3,4) * t540 + Ifges(3,2) * t544) * t564;
t507 = Ifges(3,3) * t531 + (Ifges(3,5) * t540 + Ifges(3,6) * t544) * t564;
t497 = -g(3) * t571 + t565;
t488 = Ifges(4,5) * t515 + Ifges(4,6) * t514 + Ifges(4,3) * t528;
t471 = Ifges(5,5) * t502 + Ifges(5,6) * t501 + Ifges(5,3) * t528;
t402 = m(3) * t496 + mrSges(3,1) * t530 - mrSges(3,3) * t525 + t520 * t531 - t523 * t558 + t548;
t394 = -mrSges(5,1) * t440 + mrSges(5,3) * t431 + Ifges(5,4) * t470 + Ifges(5,2) * t469 + Ifges(5,6) * t518 - pkin(4) * t408 - t502 * t471 + t528 * t473 - t581;
t392 = m(3) * t497 - mrSges(3,2) * t530 - mrSges(3,3) * t526 - t519 * t531 + t523 * t557 + t555;
t391 = mrSges(5,2) * t440 - mrSges(5,3) * t430 + Ifges(5,1) * t470 + Ifges(5,4) * t469 + Ifges(5,5) * t518 - pkin(10) * t408 - t401 * t538 + t406 * t542 + t471 * t501 - t472 * t528;
t390 = mrSges(4,2) * t480 - mrSges(4,3) * t447 + Ifges(4,1) * t495 + Ifges(4,4) * t494 + Ifges(4,5) * t518 - qJ(4) * t400 + t391 * t536 - t394 * t534 + t488 * t514 - t489 * t528;
t389 = -mrSges(4,1) * t480 + mrSges(4,3) * t448 + Ifges(4,4) * t495 + Ifges(4,2) * t494 + Ifges(4,6) * t518 - pkin(3) * t407 + qJ(4) * t554 + t534 * t391 + t536 * t394 - t515 * t488 + t528 * t490;
t388 = Ifges(3,5) * t525 - Ifges(3,6) * t526 + Ifges(3,3) * t530 + mrSges(3,1) * t496 - mrSges(3,2) * t497 + t539 * t390 + t543 * t389 + pkin(2) * t548 + pkin(9) * t555 + (t508 * t540 - t509 * t544) * t564;
t1 = [Ifges(2,3) * qJDD(1) + mrSges(2,1) * t556 - mrSges(2,2) * t551 + (mrSges(3,2) * t510 - mrSges(3,3) * t496 + Ifges(3,1) * t525 - Ifges(3,4) * t526 + Ifges(3,5) * t530 - pkin(9) * t393 - t389 * t539 + t390 * t543 + t507 * t557 - t508 * t531) * t571 + (-mrSges(3,1) * t510 + mrSges(3,3) * t497 + Ifges(3,4) * t525 - Ifges(3,2) * t526 + Ifges(3,6) * t530 - pkin(2) * t393 - t507 * t558 + t531 * t509 - t582) * t570 + t537 * t388 + pkin(1) * ((t392 * t540 + t402 * t544) * t537 + (-m(3) * t510 - t526 * mrSges(3,1) - t525 * mrSges(3,2) + (-t519 * t540 + t520 * t544) * t564 - t393) * t535) + (t392 * t544 - t402 * t540) * t580; t388; t582; t407; t581; t419;];
tauJ  = t1;
