% Calculate vector of inverse dynamics joint torques for with Newton-Euler
% S6RRRRRP7
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
%   pkin=[a2,a3,a4,a5,a6,alpha2,d1,d2,d3,d4,d5]';
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
% Datum: 2019-05-08 05:36
% Revision: 89c353f7eff3bd693eda4e29f35b2761dbc3ada0 (2019-05-03)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauJ = S6RRRRRP7_invdynJ_fixb_snew_vp2(qJ, qJD, qJDD, g, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(6,1),zeros(3,1),zeros(11,1),zeros(7,1),zeros(7,3),zeros(7,6)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRRRRP7_invdynJ_fixb_snew_vp2: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RRRRRP7_invdynJ_fixb_snew_vp2: qJD has to be [6x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [6 1]), ...
  'S6RRRRRP7_invdynJ_fixb_snew_vp2: qJDD has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6RRRRRP7_invdynJ_fixb_snew_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [11 1]), ...
  'S6RRRRRP7_invdynJ_fixb_snew_vp2: pkin has to be [11x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6RRRRRP7_invdynJ_fixb_snew_vp2: m has to be [7x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [7,3]), ...
  'S6RRRRRP7_invdynJ_fixb_snew_vp2: mrSges has to be [7x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [7 6]), ...
  'S6RRRRRP7_invdynJ_fixb_snew_vp2: Ifges has to be [7x6] (double)'); 

%% Symbolic Calculation
% From invdyn_fixb_NewtonEuler_linkframe_tauJ_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-05-08 05:22:30
% EndTime: 2019-05-08 05:22:49
% DurationCPUTime: 8.70s
% Computational Cost: add. (127573->341), mult. (270949->428), div. (0->0), fcn. (218965->12), ass. (0->146)
t588 = Ifges(6,1) + Ifges(7,1);
t581 = Ifges(6,4) + Ifges(7,4);
t580 = Ifges(6,5) + Ifges(7,5);
t587 = Ifges(6,2) + Ifges(7,2);
t579 = Ifges(6,6) + Ifges(7,6);
t586 = Ifges(6,3) + Ifges(7,3);
t539 = cos(pkin(6));
t535 = qJD(1) * t539 + qJD(2);
t542 = sin(qJ(3));
t547 = cos(qJ(3));
t543 = sin(qJ(2));
t538 = sin(pkin(6));
t569 = qJD(1) * t538;
t564 = t543 * t569;
t516 = t535 * t542 + t547 * t564;
t548 = cos(qJ(2));
t568 = qJD(1) * t548;
t527 = (qJD(2) * t568 + qJDD(1) * t543) * t538;
t534 = qJDD(1) * t539 + qJDD(2);
t496 = -qJD(3) * t516 - t527 * t542 + t534 * t547;
t515 = t535 * t547 - t542 * t564;
t497 = qJD(3) * t515 + t527 * t547 + t534 * t542;
t541 = sin(qJ(4));
t546 = cos(qJ(4));
t501 = t515 * t546 - t516 * t541;
t463 = qJD(4) * t501 + t496 * t541 + t497 * t546;
t502 = t515 * t541 + t516 * t546;
t563 = t538 * t568;
t531 = qJD(3) - t563;
t530 = qJD(4) + t531;
t540 = sin(qJ(5));
t545 = cos(qJ(5));
t486 = -t502 * t540 + t530 * t545;
t567 = qJDD(1) * t538;
t528 = -qJD(2) * t564 + t548 * t567;
t520 = qJDD(3) - t528;
t519 = qJDD(4) + t520;
t443 = qJD(5) * t486 + t463 * t545 + t519 * t540;
t487 = t502 * t545 + t530 * t540;
t466 = -mrSges(7,1) * t486 + mrSges(7,2) * t487;
t526 = (-pkin(2) * t548 - pkin(9) * t543) * t569;
t533 = t535 ^ 2;
t550 = qJD(1) ^ 2;
t544 = sin(qJ(1));
t549 = cos(qJ(1));
t557 = -g(1) * t549 - g(2) * t544;
t524 = -pkin(1) * t550 + pkin(8) * t567 + t557;
t562 = t544 * g(1) - g(2) * t549;
t584 = pkin(8) * t538;
t523 = qJDD(1) * pkin(1) + t550 * t584 + t562;
t577 = t523 * t539;
t570 = t548 * t524 + t543 * t577;
t483 = -pkin(2) * t533 + pkin(9) * t534 + (-g(3) * t543 + t526 * t568) * t538 + t570;
t583 = g(3) * t539;
t484 = -pkin(2) * t528 - pkin(9) * t527 - t583 + (-t523 + (pkin(2) * t543 - pkin(9) * t548) * t535 * qJD(1)) * t538;
t446 = -t483 * t542 + t547 * t484;
t435 = (t515 * t531 - t497) * pkin(10) + (t515 * t516 + t520) * pkin(3) + t446;
t447 = t547 * t483 + t542 * t484;
t506 = pkin(3) * t531 - pkin(10) * t516;
t514 = t515 ^ 2;
t438 = -pkin(3) * t514 + pkin(10) * t496 - t506 * t531 + t447;
t433 = t541 * t435 + t546 * t438;
t478 = -pkin(4) * t501 - pkin(11) * t502;
t529 = t530 ^ 2;
t427 = -pkin(4) * t529 + pkin(11) * t519 + t478 * t501 + t433;
t575 = t538 * t548;
t498 = -g(3) * t575 - t543 * t524 + t548 * t577;
t482 = -pkin(2) * t534 - pkin(9) * t533 + t526 * t564 - t498;
t444 = -pkin(3) * t496 - pkin(10) * t514 + t516 * t506 + t482;
t462 = -qJD(4) * t502 + t496 * t546 - t497 * t541;
t430 = (-t501 * t530 - t463) * pkin(11) + (t502 * t530 - t462) * pkin(4) + t444;
t422 = -t427 * t540 + t545 * t430;
t461 = qJDD(5) - t462;
t500 = qJD(5) - t501;
t418 = -0.2e1 * qJD(6) * t487 + (t486 * t500 - t443) * qJ(6) + (t486 * t487 + t461) * pkin(5) + t422;
t468 = -mrSges(7,2) * t500 + mrSges(7,3) * t486;
t566 = m(7) * t418 + t461 * mrSges(7,1) + t500 * t468;
t416 = -mrSges(7,3) * t443 - t466 * t487 + t566;
t423 = t545 * t427 + t540 * t430;
t442 = -qJD(5) * t487 - t463 * t540 + t519 * t545;
t470 = pkin(5) * t500 - qJ(6) * t487;
t485 = t486 ^ 2;
t421 = -pkin(5) * t485 + qJ(6) * t442 + 0.2e1 * qJD(6) * t486 - t470 * t500 + t423;
t572 = t581 * t486 + t487 * t588 + t580 * t500;
t573 = -t486 * t587 - t487 * t581 - t500 * t579;
t585 = mrSges(6,1) * t422 + mrSges(7,1) * t418 - mrSges(6,2) * t423 - mrSges(7,2) * t421 + pkin(5) * t416 + t442 * t579 + t443 * t580 + t461 * t586 - t486 * t572 - t487 * t573;
t582 = -mrSges(6,2) - mrSges(7,2);
t576 = t538 * t543;
t477 = -mrSges(5,1) * t501 + mrSges(5,2) * t502;
t489 = mrSges(5,1) * t530 - mrSges(5,3) * t502;
t467 = -mrSges(6,1) * t486 + mrSges(6,2) * t487;
t469 = -mrSges(6,2) * t500 + mrSges(6,3) * t486;
t408 = m(6) * t422 + mrSges(6,1) * t461 + t469 * t500 + (-t466 - t467) * t487 + (-mrSges(6,3) - mrSges(7,3)) * t443 + t566;
t565 = m(7) * t421 + t442 * mrSges(7,3) + t486 * t466;
t471 = mrSges(7,1) * t500 - mrSges(7,3) * t487;
t571 = -mrSges(6,1) * t500 + mrSges(6,3) * t487 - t471;
t411 = m(6) * t423 + mrSges(6,3) * t442 + t461 * t582 + t467 * t486 + t500 * t571 + t565;
t559 = -t408 * t540 + t545 * t411;
t401 = m(5) * t433 - mrSges(5,2) * t519 + mrSges(5,3) * t462 + t477 * t501 - t489 * t530 + t559;
t432 = t435 * t546 - t541 * t438;
t488 = -mrSges(5,2) * t530 + mrSges(5,3) * t501;
t426 = -pkin(4) * t519 - pkin(11) * t529 + t502 * t478 - t432;
t424 = -pkin(5) * t442 - qJ(6) * t485 + t470 * t487 + qJDD(6) + t426;
t558 = -m(7) * t424 + t442 * mrSges(7,1) + t486 * t468;
t553 = -m(6) * t426 + t442 * mrSges(6,1) + t443 * t582 + t486 * t469 + t487 * t571 + t558;
t413 = m(5) * t432 + mrSges(5,1) * t519 - mrSges(5,3) * t463 - t477 * t502 + t488 * t530 + t553;
t395 = t541 * t401 + t546 * t413;
t503 = -mrSges(4,1) * t515 + mrSges(4,2) * t516;
t504 = -mrSges(4,2) * t531 + mrSges(4,3) * t515;
t393 = m(4) * t446 + mrSges(4,1) * t520 - mrSges(4,3) * t497 - t503 * t516 + t504 * t531 + t395;
t505 = mrSges(4,1) * t531 - mrSges(4,3) * t516;
t560 = t546 * t401 - t413 * t541;
t394 = m(4) * t447 - mrSges(4,2) * t520 + mrSges(4,3) * t496 + t503 * t515 - t505 * t531 + t560;
t388 = t547 * t393 + t542 * t394;
t405 = t545 * t408 + t540 * t411;
t574 = -t486 * t579 - t487 * t580 - t500 * t586;
t561 = -t393 * t542 + t547 * t394;
t556 = m(5) * t444 - t462 * mrSges(5,1) + mrSges(5,2) * t463 - t501 * t488 + t489 * t502 + t405;
t419 = mrSges(7,2) * t443 + t471 * t487 - t558;
t397 = -mrSges(6,1) * t426 + mrSges(6,3) * t423 - mrSges(7,1) * t424 + mrSges(7,3) * t421 - pkin(5) * t419 + qJ(6) * t565 + (-qJ(6) * t471 + t572) * t500 + t574 * t487 + (-mrSges(7,2) * qJ(6) + t579) * t461 + t581 * t443 + t587 * t442;
t403 = mrSges(6,2) * t426 + mrSges(7,2) * t424 - mrSges(6,3) * t422 - mrSges(7,3) * t418 - qJ(6) * t416 + t581 * t442 + t443 * t588 + t580 * t461 - t574 * t486 + t573 * t500;
t474 = Ifges(5,4) * t502 + Ifges(5,2) * t501 + Ifges(5,6) * t530;
t475 = Ifges(5,1) * t502 + Ifges(5,4) * t501 + Ifges(5,5) * t530;
t554 = -mrSges(5,1) * t432 + mrSges(5,2) * t433 - Ifges(5,5) * t463 - Ifges(5,6) * t462 - Ifges(5,3) * t519 - pkin(4) * t553 - pkin(11) * t559 - t545 * t397 - t540 * t403 - t502 * t474 + t501 * t475;
t552 = -m(4) * t482 + t496 * mrSges(4,1) - mrSges(4,2) * t497 + t515 * t504 - t505 * t516 - t556;
t491 = Ifges(4,4) * t516 + Ifges(4,2) * t515 + Ifges(4,6) * t531;
t492 = Ifges(4,1) * t516 + Ifges(4,4) * t515 + Ifges(4,5) * t531;
t551 = mrSges(4,1) * t446 - mrSges(4,2) * t447 + Ifges(4,5) * t497 + Ifges(4,6) * t496 + Ifges(4,3) * t520 + pkin(3) * t395 + t516 * t491 - t515 * t492 - t554;
t525 = (-mrSges(3,1) * t548 + mrSges(3,2) * t543) * t569;
t522 = -mrSges(3,2) * t535 + mrSges(3,3) * t563;
t521 = mrSges(3,1) * t535 - mrSges(3,3) * t564;
t510 = -t523 * t538 - t583;
t509 = Ifges(3,5) * t535 + (Ifges(3,1) * t543 + Ifges(3,4) * t548) * t569;
t508 = Ifges(3,6) * t535 + (Ifges(3,4) * t543 + Ifges(3,2) * t548) * t569;
t507 = Ifges(3,3) * t535 + (Ifges(3,5) * t543 + Ifges(3,6) * t548) * t569;
t499 = -g(3) * t576 + t570;
t490 = Ifges(4,5) * t516 + Ifges(4,6) * t515 + Ifges(4,3) * t531;
t473 = Ifges(5,5) * t502 + Ifges(5,6) * t501 + Ifges(5,3) * t530;
t398 = m(3) * t498 + mrSges(3,1) * t534 - mrSges(3,3) * t527 + t522 * t535 - t525 * t564 + t552;
t389 = -mrSges(5,1) * t444 + mrSges(5,3) * t433 + Ifges(5,4) * t463 + Ifges(5,2) * t462 + Ifges(5,6) * t519 - pkin(4) * t405 - t502 * t473 + t530 * t475 - t585;
t387 = m(3) * t499 - mrSges(3,2) * t534 + mrSges(3,3) * t528 - t521 * t535 + t525 * t563 + t561;
t386 = mrSges(5,2) * t444 - mrSges(5,3) * t432 + Ifges(5,1) * t463 + Ifges(5,4) * t462 + Ifges(5,5) * t519 - pkin(11) * t405 - t397 * t540 + t403 * t545 + t473 * t501 - t474 * t530;
t385 = mrSges(4,2) * t482 - mrSges(4,3) * t446 + Ifges(4,1) * t497 + Ifges(4,4) * t496 + Ifges(4,5) * t520 - pkin(10) * t395 + t386 * t546 - t389 * t541 + t490 * t515 - t491 * t531;
t384 = -mrSges(4,1) * t482 + mrSges(4,3) * t447 + Ifges(4,4) * t497 + Ifges(4,2) * t496 + Ifges(4,6) * t520 - pkin(3) * t556 + pkin(10) * t560 + t541 * t386 + t546 * t389 - t516 * t490 + t531 * t492;
t383 = Ifges(3,5) * t527 + Ifges(3,6) * t528 + Ifges(3,3) * t534 + mrSges(3,1) * t498 - mrSges(3,2) * t499 + t542 * t385 + t547 * t384 + pkin(2) * t552 + pkin(9) * t561 + (t508 * t543 - t509 * t548) * t569;
t1 = [Ifges(2,3) * qJDD(1) + mrSges(2,1) * t562 - mrSges(2,2) * t557 + (mrSges(3,2) * t510 - mrSges(3,3) * t498 + Ifges(3,1) * t527 + Ifges(3,4) * t528 + Ifges(3,5) * t534 - pkin(9) * t388 - t384 * t542 + t385 * t547 + t507 * t563 - t508 * t535) * t576 + (-mrSges(3,1) * t510 + mrSges(3,3) * t499 + Ifges(3,4) * t527 + Ifges(3,2) * t528 + Ifges(3,6) * t534 - pkin(2) * t388 - t507 * t564 + t535 * t509 - t551) * t575 + t539 * t383 + pkin(1) * ((t387 * t543 + t398 * t548) * t539 + (-m(3) * t510 + t528 * mrSges(3,1) - t527 * mrSges(3,2) + (-t521 * t543 + t522 * t548) * t569 - t388) * t538) + (t387 * t548 - t398 * t543) * t584; t383; t551; -t554; t585; t419;];
tauJ  = t1;
