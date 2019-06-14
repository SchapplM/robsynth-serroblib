% Calculate vector of inverse dynamics joint torques for with Newton-Euler
% S6RRPPPR1
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
% pkin [10x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d2,d6,theta3,theta4]';
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
% Datum: 2019-05-06 08:15
% Revision: 89c353f7eff3bd693eda4e29f35b2761dbc3ada0 (2019-05-03)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauJ = S6RRPPPR1_invdynJ_fixb_snew_vp2(qJ, qJD, qJDD, g, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(6,1),zeros(3,1),zeros(10,1),zeros(7,1),zeros(7,3),zeros(7,6)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRPPPR1_invdynJ_fixb_snew_vp2: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RRPPPR1_invdynJ_fixb_snew_vp2: qJD has to be [6x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [6 1]), ...
  'S6RRPPPR1_invdynJ_fixb_snew_vp2: qJDD has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6RRPPPR1_invdynJ_fixb_snew_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S6RRPPPR1_invdynJ_fixb_snew_vp2: pkin has to be [10x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6RRPPPR1_invdynJ_fixb_snew_vp2: m has to be [7x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [7,3]), ...
  'S6RRPPPR1_invdynJ_fixb_snew_vp2: mrSges has to be [7x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [7 6]), ...
  'S6RRPPPR1_invdynJ_fixb_snew_vp2: Ifges has to be [7x6] (double)'); 

%% Symbolic Calculation
% From invdyn_fixb_NewtonEuler_linkframe_tauJ_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-05-06 08:10:58
% EndTime: 2019-05-06 08:11:05
% DurationCPUTime: 3.95s
% Computational Cost: add. (33437->330), mult. (79493->404), div. (0->0), fcn. (55103->10), ass. (0->133)
t587 = -2 * qJD(4);
t586 = Ifges(5,1) + Ifges(6,1);
t579 = Ifges(5,4) - Ifges(6,5);
t578 = -Ifges(5,5) - Ifges(6,4);
t585 = Ifges(5,2) + Ifges(6,3);
t584 = -Ifges(6,2) - Ifges(5,3);
t577 = Ifges(5,6) - Ifges(6,6);
t543 = sin(qJ(2));
t546 = cos(qJ(2));
t563 = qJD(1) * qJD(2);
t530 = qJDD(1) * t543 + t546 * t563;
t531 = qJDD(1) * t546 - t543 * t563;
t541 = sin(pkin(9));
t576 = cos(pkin(9));
t503 = t530 * t576 + t531 * t541;
t540 = sin(pkin(10));
t575 = cos(pkin(10));
t491 = qJDD(2) * t540 + t503 * t575;
t521 = (t541 * t546 + t543 * t576) * qJD(1);
t507 = -qJD(2) * t575 + t521 * t540;
t566 = qJD(1) * t546;
t567 = qJD(1) * t543;
t520 = t541 * t567 - t566 * t576;
t574 = t507 * t520;
t583 = (-t491 + t574) * qJ(5);
t549 = qJD(1) ^ 2;
t544 = sin(qJ(1));
t547 = cos(qJ(1));
t557 = -g(1) * t547 - g(2) * t544;
t527 = -pkin(1) * t549 + qJDD(1) * pkin(7) + t557;
t573 = t543 * t527;
t581 = pkin(2) * t549;
t476 = qJDD(2) * pkin(2) - t530 * qJ(3) - t573 + (qJ(3) * t563 + t543 * t581 - g(3)) * t546;
t510 = -g(3) * t543 + t527 * t546;
t532 = qJD(2) * pkin(2) - qJ(3) * t567;
t539 = t546 ^ 2;
t480 = qJ(3) * t531 - qJD(2) * t532 - t539 * t581 + t510;
t454 = -0.2e1 * qJD(3) * t520 + t476 * t541 + t480 * t576;
t495 = pkin(3) * t520 - qJ(4) * t521;
t548 = qJD(2) ^ 2;
t439 = -pkin(3) * t548 + qJDD(2) * qJ(4) - t495 * t520 + t454;
t562 = t544 * g(1) - g(2) * t547;
t554 = -qJDD(1) * pkin(1) - t562;
t482 = -t531 * pkin(2) + qJDD(3) + t532 * t567 + (-qJ(3) * t539 - pkin(7)) * t549 + t554;
t502 = t530 * t541 - t531 * t576;
t441 = (qJD(2) * t520 - t503) * qJ(4) + (qJD(2) * t521 + t502) * pkin(3) + t482;
t508 = qJD(2) * t540 + t521 * t575;
t433 = -t439 * t540 + t441 * t575 + t508 * t587;
t582 = 2 * qJD(5);
t580 = -mrSges(5,3) - mrSges(6,2);
t496 = mrSges(4,1) * t520 + mrSges(4,2) * t521;
t512 = qJD(2) * mrSges(4,1) - mrSges(4,3) * t521;
t434 = t439 * t575 + t441 * t540 + t507 * t587;
t485 = mrSges(5,1) * t520 - mrSges(5,3) * t508;
t490 = -qJDD(2) * t575 + t503 * t540;
t477 = pkin(4) * t507 - qJ(5) * t508;
t519 = t520 ^ 2;
t430 = -pkin(4) * t519 + qJ(5) * t502 - t477 * t507 + t520 * t582 + t434;
t486 = -mrSges(6,1) * t520 + mrSges(6,2) * t508;
t431 = -t502 * pkin(4) - t519 * qJ(5) + t477 * t508 + qJDD(5) - t433;
t425 = (-t491 - t574) * pkin(8) + (t507 * t508 - t502) * pkin(5) + t431;
t487 = -pkin(5) * t520 - pkin(8) * t508;
t506 = t507 ^ 2;
t426 = -pkin(5) * t506 + pkin(8) * t490 + t487 * t520 + t430;
t542 = sin(qJ(6));
t545 = cos(qJ(6));
t423 = t425 * t545 - t426 * t542;
t472 = t507 * t545 - t508 * t542;
t447 = qJD(6) * t472 + t490 * t542 + t491 * t545;
t473 = t507 * t542 + t508 * t545;
t455 = -mrSges(7,1) * t472 + mrSges(7,2) * t473;
t518 = qJD(6) - t520;
t458 = -mrSges(7,2) * t518 + mrSges(7,3) * t472;
t501 = qJDD(6) - t502;
t419 = m(7) * t423 + mrSges(7,1) * t501 - mrSges(7,3) * t447 - t455 * t473 + t458 * t518;
t424 = t425 * t542 + t426 * t545;
t446 = -qJD(6) * t473 + t490 * t545 - t491 * t542;
t459 = mrSges(7,1) * t518 - mrSges(7,3) * t473;
t420 = m(7) * t424 - mrSges(7,2) * t501 + mrSges(7,3) * t446 + t455 * t472 - t459 * t518;
t559 = -t542 * t419 + t420 * t545;
t553 = m(6) * t430 + mrSges(6,3) * t502 + t486 * t520 + t559;
t478 = mrSges(6,1) * t507 - mrSges(6,3) * t508;
t568 = -mrSges(5,1) * t507 - mrSges(5,2) * t508 - t478;
t408 = m(5) * t434 - t502 * mrSges(5,2) - t520 * t485 + t490 * t580 + t507 * t568 + t553;
t484 = -mrSges(5,2) * t520 - mrSges(5,3) * t507;
t412 = t545 * t419 + t542 * t420;
t483 = -mrSges(6,2) * t507 + mrSges(6,3) * t520;
t552 = -m(6) * t431 + mrSges(6,1) * t502 + t483 * t520 - t412;
t410 = m(5) * t433 + t502 * mrSges(5,1) + t520 * t484 + t491 * t580 + t508 * t568 + t552;
t560 = t408 * t575 - t410 * t540;
t403 = m(4) * t454 - qJDD(2) * mrSges(4,2) - mrSges(4,3) * t502 - qJD(2) * t512 - t496 * t520 + t560;
t569 = t476 * t576 - t480 * t541;
t551 = qJDD(2) * pkin(3) + t548 * qJ(4) - t495 * t521 - qJDD(4) + t569;
t565 = qJD(3) * t521;
t438 = 0.2e1 * t565 - t551;
t432 = -0.2e1 * qJD(5) * t508 + t583 + (t508 * t520 + t490) * pkin(4) + t438;
t515 = -0.2e1 * t565;
t429 = -t506 * pkin(8) + t515 + (-pkin(4) - pkin(5)) * t490 - t583 + (-pkin(4) * t520 + t487 + t582) * t508 + t551;
t555 = -m(7) * t429 + mrSges(7,1) * t446 - mrSges(7,2) * t447 + t458 * t472 - t459 * t473;
t421 = m(6) * t432 + mrSges(6,1) * t490 - mrSges(6,3) * t491 + t483 * t507 - t486 * t508 + t555;
t417 = m(5) * t438 + mrSges(5,1) * t490 + mrSges(5,2) * t491 + t484 * t507 + t485 * t508 + t421;
t453 = t515 + t569;
t511 = -qJD(2) * mrSges(4,2) - mrSges(4,3) * t520;
t416 = m(4) * t453 + qJDD(2) * mrSges(4,1) - mrSges(4,3) * t503 + qJD(2) * t511 - t496 * t521 - t417;
t400 = t403 * t541 + t416 * t576;
t405 = t408 * t540 + t410 * t575;
t572 = t507 * t585 - t508 * t579 - t520 * t577;
t571 = t507 * t577 + t508 * t578 + t520 * t584;
t570 = -t507 * t579 + t508 * t586 - t520 * t578;
t561 = t403 * t576 - t416 * t541;
t449 = Ifges(7,4) * t473 + Ifges(7,2) * t472 + Ifges(7,6) * t518;
t450 = Ifges(7,1) * t473 + Ifges(7,4) * t472 + Ifges(7,5) * t518;
t550 = mrSges(7,1) * t423 - mrSges(7,2) * t424 + Ifges(7,5) * t447 + Ifges(7,6) * t446 + Ifges(7,3) * t501 + t449 * t473 - t472 * t450;
t404 = m(4) * t482 + t502 * mrSges(4,1) + t503 * mrSges(4,2) + t520 * t511 + t521 * t512 + t405;
t534 = -qJD(2) * mrSges(3,2) + mrSges(3,3) * t566;
t533 = qJD(2) * mrSges(3,1) - mrSges(3,3) * t567;
t529 = (-mrSges(3,1) * t546 + mrSges(3,2) * t543) * qJD(1);
t526 = -t549 * pkin(7) + t554;
t524 = Ifges(3,5) * qJD(2) + (Ifges(3,1) * t543 + Ifges(3,4) * t546) * qJD(1);
t523 = Ifges(3,6) * qJD(2) + (Ifges(3,4) * t543 + Ifges(3,2) * t546) * qJD(1);
t509 = -t546 * g(3) - t573;
t494 = Ifges(4,1) * t521 - Ifges(4,4) * t520 + Ifges(4,5) * qJD(2);
t493 = Ifges(4,4) * t521 - Ifges(4,2) * t520 + Ifges(4,6) * qJD(2);
t492 = Ifges(4,5) * t521 - Ifges(4,6) * t520 + Ifges(4,3) * qJD(2);
t448 = Ifges(7,5) * t473 + Ifges(7,6) * t472 + Ifges(7,3) * t518;
t414 = mrSges(7,2) * t429 - mrSges(7,3) * t423 + Ifges(7,1) * t447 + Ifges(7,4) * t446 + Ifges(7,5) * t501 + t448 * t472 - t449 * t518;
t413 = -mrSges(7,1) * t429 + mrSges(7,3) * t424 + Ifges(7,4) * t447 + Ifges(7,2) * t446 + Ifges(7,6) * t501 - t448 * t473 + t450 * t518;
t411 = t491 * mrSges(6,2) + t508 * t478 - t552;
t399 = mrSges(5,2) * t438 + mrSges(6,2) * t431 - mrSges(5,3) * t433 - mrSges(6,3) * t432 - pkin(8) * t412 - qJ(5) * t421 - t542 * t413 + t545 * t414 - t490 * t579 + t491 * t586 - t502 * t578 + t507 * t571 + t520 * t572;
t398 = -mrSges(5,1) * t438 - mrSges(6,1) * t432 + mrSges(6,2) * t430 + mrSges(5,3) * t434 - pkin(4) * t421 - pkin(5) * t555 - pkin(8) * t559 - t545 * t413 - t542 * t414 - t490 * t585 + t491 * t579 + t502 * t577 + t508 * t571 + t520 * t570;
t397 = -t521 * t492 + Ifges(4,4) * t503 + qJD(2) * t494 - mrSges(4,1) * t482 + t572 * t508 + mrSges(4,3) * t454 - mrSges(5,1) * t433 + mrSges(5,2) * t434 - mrSges(6,3) * t430 + mrSges(6,1) * t431 + (qJ(5) * t478 - t570) * t507 + (-Ifges(4,2) + t584) * t502 + t578 * t491 + (mrSges(6,2) * qJ(5) + t577) * t490 + Ifges(4,6) * qJDD(2) + pkin(5) * t412 + pkin(4) * t411 - pkin(3) * t405 + t550 - qJ(5) * t553;
t396 = mrSges(4,2) * t482 - mrSges(4,3) * t453 + Ifges(4,1) * t503 - Ifges(4,4) * t502 + Ifges(4,5) * qJDD(2) - qJ(4) * t405 - qJD(2) * t493 - t398 * t540 + t399 * t575 - t492 * t520;
t1 = [Ifges(2,3) * qJDD(1) + mrSges(2,1) * t562 - mrSges(2,2) * t557 + t543 * (mrSges(3,2) * t526 - mrSges(3,3) * t509 + Ifges(3,1) * t530 + Ifges(3,4) * t531 + Ifges(3,5) * qJDD(2) - qJ(3) * t400 - qJD(2) * t523 + t576 * t396 - t541 * t397) + t546 * (-mrSges(3,1) * t526 + mrSges(3,3) * t510 + Ifges(3,4) * t530 + Ifges(3,2) * t531 + Ifges(3,6) * qJDD(2) - pkin(2) * t404 + qJ(3) * t561 + qJD(2) * t524 + t541 * t396 + t576 * t397) + pkin(1) * (-m(3) * t526 + t531 * mrSges(3,1) - t530 * mrSges(3,2) + (-t533 * t543 + t534 * t546) * qJD(1) - t404) + pkin(7) * (t546 * (m(3) * t510 - qJDD(2) * mrSges(3,2) + mrSges(3,3) * t531 - qJD(2) * t533 + t529 * t566 + t561) - t543 * (m(3) * t509 + qJDD(2) * mrSges(3,1) - mrSges(3,3) * t530 + qJD(2) * t534 - t529 * t567 + t400)); Ifges(3,5) * t530 + Ifges(3,6) * t531 + mrSges(3,1) * t509 - mrSges(3,2) * t510 + Ifges(4,5) * t503 - Ifges(4,6) * t502 + t521 * t493 + t520 * t494 + mrSges(4,1) * t453 - mrSges(4,2) * t454 + t540 * t399 + t575 * t398 - pkin(3) * t417 + qJ(4) * t560 + pkin(2) * t400 + (Ifges(3,3) + Ifges(4,3)) * qJDD(2) + (t543 * t523 - t546 * t524) * qJD(1); t404; t417; t411; t550;];
tauJ  = t1;
