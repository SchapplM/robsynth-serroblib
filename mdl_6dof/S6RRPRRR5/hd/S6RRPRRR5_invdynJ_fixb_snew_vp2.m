% Calculate vector of inverse dynamics joint torques for with Newton-Euler
% S6RRPRRR5
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
% Datum: 2019-05-06 21:30
% Revision: 89c353f7eff3bd693eda4e29f35b2761dbc3ada0 (2019-05-03)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauJ = S6RRPRRR5_invdynJ_fixb_snew_vp2(qJ, qJD, qJDD, g, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(6,1),zeros(3,1),zeros(12,1),zeros(7,1),zeros(7,3),zeros(7,6)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRPRRR5_invdynJ_fixb_snew_vp2: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RRPRRR5_invdynJ_fixb_snew_vp2: qJD has to be [6x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [6 1]), ...
  'S6RRPRRR5_invdynJ_fixb_snew_vp2: qJDD has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6RRPRRR5_invdynJ_fixb_snew_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [12 1]), ...
  'S6RRPRRR5_invdynJ_fixb_snew_vp2: pkin has to be [12x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6RRPRRR5_invdynJ_fixb_snew_vp2: m has to be [7x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [7,3]), ...
  'S6RRPRRR5_invdynJ_fixb_snew_vp2: mrSges has to be [7x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [7 6]), ...
  'S6RRPRRR5_invdynJ_fixb_snew_vp2: Ifges has to be [7x6] (double)'); 

%% Symbolic Calculation
% From invdyn_fixb_NewtonEuler_linkframe_tauJ_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-05-06 21:18:51
% EndTime: 2019-05-06 21:19:04
% DurationCPUTime: 12.70s
% Computational Cost: add. (179863->362), mult. (470401->471), div. (0->0), fcn. (377427->14), ass. (0->152)
t589 = -2 * qJD(3);
t549 = sin(pkin(12));
t551 = cos(pkin(12));
t556 = sin(qJ(2));
t561 = cos(qJ(2));
t550 = sin(pkin(6));
t583 = qJD(1) * t550;
t532 = (t549 * t556 - t551 * t561) * t583;
t581 = qJD(1) * qJD(2);
t541 = (qJDD(1) * t556 + t561 * t581) * t550;
t552 = cos(pkin(6));
t544 = qJDD(1) * t552 + qJDD(2);
t545 = qJD(1) * t552 + qJD(2);
t563 = qJD(1) ^ 2;
t557 = sin(qJ(1));
t562 = cos(qJ(1));
t571 = -g(1) * t562 - g(2) * t557;
t588 = pkin(8) * t550;
t539 = -pkin(1) * t563 + qJDD(1) * t588 + t571;
t577 = t557 * g(1) - g(2) * t562;
t538 = qJDD(1) * pkin(1) + t563 * t588 + t577;
t587 = t538 * t552;
t572 = -t556 * t539 + t561 * t587;
t586 = t550 ^ 2 * t563;
t475 = t544 * pkin(2) - t541 * qJ(3) + (pkin(2) * t556 * t586 + (qJ(3) * qJD(1) * t545 - g(3)) * t550) * t561 + t572;
t585 = t550 * t556;
t505 = -g(3) * t585 + t561 * t539 + t556 * t587;
t579 = t556 * t583;
t535 = pkin(2) * t545 - qJ(3) * t579;
t542 = (qJDD(1) * t561 - t556 * t581) * t550;
t580 = t561 ^ 2 * t586;
t479 = -pkin(2) * t580 + qJ(3) * t542 - t535 * t545 + t505;
t533 = (t549 * t561 + t551 * t556) * t583;
t455 = t551 * t475 - t549 * t479 + t533 * t589;
t584 = t550 * t561;
t456 = t549 * t475 + t551 * t479 + t532 * t589;
t506 = mrSges(4,1) * t532 + mrSges(4,2) * t533;
t511 = -t541 * t549 + t542 * t551;
t519 = mrSges(4,1) * t545 - mrSges(4,3) * t533;
t507 = pkin(3) * t532 - pkin(9) * t533;
t543 = t545 ^ 2;
t449 = -pkin(3) * t543 + pkin(9) * t544 - t507 * t532 + t456;
t523 = -t552 * g(3) - t550 * t538;
t491 = -t542 * pkin(2) - qJ(3) * t580 + t535 * t579 + qJDD(3) + t523;
t512 = t541 * t551 + t542 * t549;
t458 = (t532 * t545 - t512) * pkin(9) + (t533 * t545 - t511) * pkin(3) + t491;
t555 = sin(qJ(4));
t560 = cos(qJ(4));
t444 = t560 * t449 + t555 * t458;
t516 = -t555 * t533 + t545 * t560;
t517 = t533 * t560 + t545 * t555;
t493 = -pkin(4) * t516 - pkin(10) * t517;
t510 = qJDD(4) - t511;
t531 = qJD(4) + t532;
t530 = t531 ^ 2;
t434 = -pkin(4) * t530 + pkin(10) * t510 + t493 * t516 + t444;
t448 = -t544 * pkin(3) - t543 * pkin(9) + t533 * t507 - t455;
t488 = -t517 * qJD(4) - t555 * t512 + t544 * t560;
t489 = qJD(4) * t516 + t512 * t560 + t544 * t555;
t437 = (-t516 * t531 - t489) * pkin(10) + (t517 * t531 - t488) * pkin(4) + t448;
t554 = sin(qJ(5));
t559 = cos(qJ(5));
t429 = -t554 * t434 + t559 * t437;
t496 = -t517 * t554 + t531 * t559;
t461 = qJD(5) * t496 + t489 * t559 + t510 * t554;
t487 = qJDD(5) - t488;
t497 = t517 * t559 + t531 * t554;
t515 = qJD(5) - t516;
t427 = (t496 * t515 - t461) * pkin(11) + (t496 * t497 + t487) * pkin(5) + t429;
t430 = t559 * t434 + t554 * t437;
t460 = -qJD(5) * t497 - t489 * t554 + t510 * t559;
t478 = pkin(5) * t515 - pkin(11) * t497;
t495 = t496 ^ 2;
t428 = -pkin(5) * t495 + pkin(11) * t460 - t478 * t515 + t430;
t553 = sin(qJ(6));
t558 = cos(qJ(6));
t425 = t427 * t558 - t428 * t553;
t468 = t496 * t558 - t497 * t553;
t442 = qJD(6) * t468 + t460 * t553 + t461 * t558;
t469 = t496 * t553 + t497 * t558;
t454 = -mrSges(7,1) * t468 + mrSges(7,2) * t469;
t513 = qJD(6) + t515;
t462 = -mrSges(7,2) * t513 + mrSges(7,3) * t468;
t485 = qJDD(6) + t487;
t421 = m(7) * t425 + mrSges(7,1) * t485 - mrSges(7,3) * t442 - t454 * t469 + t462 * t513;
t426 = t427 * t553 + t428 * t558;
t441 = -qJD(6) * t469 + t460 * t558 - t461 * t553;
t463 = mrSges(7,1) * t513 - mrSges(7,3) * t469;
t422 = m(7) * t426 - mrSges(7,2) * t485 + mrSges(7,3) * t441 + t454 * t468 - t463 * t513;
t413 = t558 * t421 + t553 * t422;
t470 = -mrSges(6,1) * t496 + mrSges(6,2) * t497;
t476 = -mrSges(6,2) * t515 + mrSges(6,3) * t496;
t411 = m(6) * t429 + mrSges(6,1) * t487 - mrSges(6,3) * t461 - t470 * t497 + t476 * t515 + t413;
t477 = mrSges(6,1) * t515 - mrSges(6,3) * t497;
t574 = -t421 * t553 + t558 * t422;
t412 = m(6) * t430 - mrSges(6,2) * t487 + mrSges(6,3) * t460 + t470 * t496 - t477 * t515 + t574;
t409 = -t411 * t554 + t559 * t412;
t492 = -mrSges(5,1) * t516 + mrSges(5,2) * t517;
t499 = mrSges(5,1) * t531 - mrSges(5,3) * t517;
t407 = m(5) * t444 - mrSges(5,2) * t510 + mrSges(5,3) * t488 + t492 * t516 - t499 * t531 + t409;
t443 = -t555 * t449 + t458 * t560;
t433 = -pkin(4) * t510 - pkin(10) * t530 + t517 * t493 - t443;
t431 = -pkin(5) * t460 - pkin(11) * t495 + t478 * t497 + t433;
t569 = m(7) * t431 - t441 * mrSges(7,1) + mrSges(7,2) * t442 - t468 * t462 + t463 * t469;
t423 = -m(6) * t433 + t460 * mrSges(6,1) - mrSges(6,2) * t461 + t496 * t476 - t477 * t497 - t569;
t498 = -mrSges(5,2) * t531 + mrSges(5,3) * t516;
t417 = m(5) * t443 + mrSges(5,1) * t510 - mrSges(5,3) * t489 - t492 * t517 + t498 * t531 + t423;
t575 = t560 * t407 - t417 * t555;
t398 = m(4) * t456 - mrSges(4,2) * t544 + mrSges(4,3) * t511 - t506 * t532 - t519 * t545 + t575;
t518 = -mrSges(4,2) * t545 - mrSges(4,3) * t532;
t408 = t411 * t559 + t412 * t554;
t566 = -m(5) * t448 + t488 * mrSges(5,1) - mrSges(5,2) * t489 + t516 * t498 - t499 * t517 - t408;
t404 = m(4) * t455 + mrSges(4,1) * t544 - mrSges(4,3) * t512 - t506 * t533 + t518 * t545 + t566;
t394 = t549 * t398 + t551 * t404;
t400 = t555 * t407 + t560 * t417;
t578 = t561 * t583;
t576 = t551 * t398 - t404 * t549;
t568 = -m(4) * t491 + t511 * mrSges(4,1) - t512 * mrSges(4,2) - t532 * t518 - t533 * t519 - t400;
t451 = Ifges(7,4) * t469 + Ifges(7,2) * t468 + Ifges(7,6) * t513;
t452 = Ifges(7,1) * t469 + Ifges(7,4) * t468 + Ifges(7,5) * t513;
t567 = -mrSges(7,1) * t425 + mrSges(7,2) * t426 - Ifges(7,5) * t442 - Ifges(7,6) * t441 - Ifges(7,3) * t485 - t469 * t451 + t468 * t452;
t450 = Ifges(7,5) * t469 + Ifges(7,6) * t468 + Ifges(7,3) * t513;
t414 = -mrSges(7,1) * t431 + mrSges(7,3) * t426 + Ifges(7,4) * t442 + Ifges(7,2) * t441 + Ifges(7,6) * t485 - t450 * t469 + t452 * t513;
t415 = mrSges(7,2) * t431 - mrSges(7,3) * t425 + Ifges(7,1) * t442 + Ifges(7,4) * t441 + Ifges(7,5) * t485 + t450 * t468 - t451 * t513;
t464 = Ifges(6,5) * t497 + Ifges(6,6) * t496 + Ifges(6,3) * t515;
t466 = Ifges(6,1) * t497 + Ifges(6,4) * t496 + Ifges(6,5) * t515;
t401 = -mrSges(6,1) * t433 + mrSges(6,3) * t430 + Ifges(6,4) * t461 + Ifges(6,2) * t460 + Ifges(6,6) * t487 - pkin(5) * t569 + pkin(11) * t574 + t558 * t414 + t553 * t415 - t497 * t464 + t515 * t466;
t465 = Ifges(6,4) * t497 + Ifges(6,2) * t496 + Ifges(6,6) * t515;
t402 = mrSges(6,2) * t433 - mrSges(6,3) * t429 + Ifges(6,1) * t461 + Ifges(6,4) * t460 + Ifges(6,5) * t487 - pkin(11) * t413 - t414 * t553 + t415 * t558 + t464 * t496 - t465 * t515;
t481 = Ifges(5,4) * t517 + Ifges(5,2) * t516 + Ifges(5,6) * t531;
t482 = Ifges(5,1) * t517 + Ifges(5,4) * t516 + Ifges(5,5) * t531;
t565 = mrSges(5,1) * t443 - mrSges(5,2) * t444 + Ifges(5,5) * t489 + Ifges(5,6) * t488 + Ifges(5,3) * t510 + pkin(4) * t423 + pkin(10) * t409 + t559 * t401 + t554 * t402 + t517 * t481 - t516 * t482;
t564 = mrSges(6,1) * t429 - mrSges(6,2) * t430 + Ifges(6,5) * t461 + Ifges(6,6) * t460 + Ifges(6,3) * t487 + pkin(5) * t413 + t497 * t465 - t496 * t466 - t567;
t540 = (-mrSges(3,1) * t561 + mrSges(3,2) * t556) * t583;
t537 = -mrSges(3,2) * t545 + mrSges(3,3) * t578;
t536 = mrSges(3,1) * t545 - mrSges(3,3) * t579;
t522 = Ifges(3,5) * t545 + (Ifges(3,1) * t556 + Ifges(3,4) * t561) * t583;
t521 = Ifges(3,6) * t545 + (Ifges(3,4) * t556 + Ifges(3,2) * t561) * t583;
t520 = Ifges(3,3) * t545 + (Ifges(3,5) * t556 + Ifges(3,6) * t561) * t583;
t504 = -g(3) * t584 + t572;
t502 = Ifges(4,1) * t533 - Ifges(4,4) * t532 + Ifges(4,5) * t545;
t501 = Ifges(4,4) * t533 - Ifges(4,2) * t532 + Ifges(4,6) * t545;
t500 = Ifges(4,5) * t533 - Ifges(4,6) * t532 + Ifges(4,3) * t545;
t480 = Ifges(5,5) * t517 + Ifges(5,6) * t516 + Ifges(5,3) * t531;
t395 = -mrSges(5,1) * t448 + mrSges(5,3) * t444 + Ifges(5,4) * t489 + Ifges(5,2) * t488 + Ifges(5,6) * t510 - pkin(4) * t408 - t517 * t480 + t531 * t482 - t564;
t393 = m(3) * t505 - mrSges(3,2) * t544 + mrSges(3,3) * t542 - t536 * t545 + t540 * t578 + t576;
t392 = m(3) * t504 + mrSges(3,1) * t544 - mrSges(3,3) * t541 + t537 * t545 - t540 * t579 + t394;
t391 = mrSges(5,2) * t448 - mrSges(5,3) * t443 + Ifges(5,1) * t489 + Ifges(5,4) * t488 + Ifges(5,5) * t510 - pkin(10) * t408 - t401 * t554 + t402 * t559 + t480 * t516 - t481 * t531;
t390 = -mrSges(4,1) * t491 + mrSges(4,3) * t456 + Ifges(4,4) * t512 + Ifges(4,2) * t511 + Ifges(4,6) * t544 - pkin(3) * t400 - t533 * t500 + t545 * t502 - t565;
t389 = mrSges(4,2) * t491 - mrSges(4,3) * t455 + Ifges(4,1) * t512 + Ifges(4,4) * t511 + Ifges(4,5) * t544 - pkin(9) * t400 + t391 * t560 - t395 * t555 - t500 * t532 - t501 * t545;
t388 = Ifges(3,5) * t541 + Ifges(3,6) * t542 + mrSges(3,1) * t504 - mrSges(3,2) * t505 + Ifges(4,5) * t512 + Ifges(4,6) * t511 + t533 * t501 + t532 * t502 + mrSges(4,1) * t455 - mrSges(4,2) * t456 + t555 * t391 + t560 * t395 + pkin(3) * t566 + pkin(9) * t575 + pkin(2) * t394 + (Ifges(3,3) + Ifges(4,3)) * t544 + (t521 * t556 - t522 * t561) * t583;
t1 = [Ifges(2,3) * qJDD(1) + mrSges(2,1) * t577 - mrSges(2,2) * t571 + (mrSges(3,2) * t523 - mrSges(3,3) * t504 + Ifges(3,1) * t541 + Ifges(3,4) * t542 + Ifges(3,5) * t544 - qJ(3) * t394 + t389 * t551 - t390 * t549 + t520 * t578 - t521 * t545) * t585 + (-mrSges(3,1) * t523 + mrSges(3,3) * t505 + Ifges(3,4) * t541 + Ifges(3,2) * t542 + Ifges(3,6) * t544 + pkin(2) * t568 + qJ(3) * t576 + t549 * t389 + t551 * t390 - t520 * t579 + t545 * t522) * t584 + t552 * t388 + pkin(1) * ((t392 * t561 + t393 * t556) * t552 + (-m(3) * t523 + t542 * mrSges(3,1) - t541 * mrSges(3,2) + (-t536 * t556 + t537 * t561) * t583 + t568) * t550) + (-t392 * t556 + t393 * t561) * t588; t388; -t568; t565; t564; -t567;];
tauJ  = t1;
