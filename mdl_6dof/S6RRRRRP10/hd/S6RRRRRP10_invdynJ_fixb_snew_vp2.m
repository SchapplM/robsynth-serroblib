% Calculate vector of inverse dynamics joint torques for with Newton-Euler
% S6RRRRRP10
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
% Datum: 2019-05-08 06:39
% Revision: 89c353f7eff3bd693eda4e29f35b2761dbc3ada0 (2019-05-03)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauJ = S6RRRRRP10_invdynJ_fixb_snew_vp2(qJ, qJD, qJDD, g, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(6,1),zeros(3,1),zeros(11,1),zeros(7,1),zeros(7,3),zeros(7,6)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRRRRP10_invdynJ_fixb_snew_vp2: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RRRRRP10_invdynJ_fixb_snew_vp2: qJD has to be [6x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [6 1]), ...
  'S6RRRRRP10_invdynJ_fixb_snew_vp2: qJDD has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6RRRRRP10_invdynJ_fixb_snew_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [11 1]), ...
  'S6RRRRRP10_invdynJ_fixb_snew_vp2: pkin has to be [11x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6RRRRRP10_invdynJ_fixb_snew_vp2: m has to be [7x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [7,3]), ...
  'S6RRRRRP10_invdynJ_fixb_snew_vp2: mrSges has to be [7x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [7 6]), ...
  'S6RRRRRP10_invdynJ_fixb_snew_vp2: Ifges has to be [7x6] (double)'); 

%% Symbolic Calculation
% From invdyn_fixb_NewtonEuler_linkframe_tauJ_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-05-08 06:20:20
% EndTime: 2019-05-08 06:20:43
% DurationCPUTime: 9.10s
% Computational Cost: add. (128512->339), mult. (272635->429), div. (0->0), fcn. (218114->12), ass. (0->145)
t586 = Ifges(6,1) + Ifges(7,1);
t575 = Ifges(6,4) - Ifges(7,5);
t584 = Ifges(7,4) + Ifges(6,5);
t585 = Ifges(6,2) + Ifges(7,3);
t583 = Ifges(6,6) - Ifges(7,6);
t582 = -Ifges(6,3) - Ifges(7,2);
t541 = cos(pkin(6));
t537 = qJD(1) * t541 + qJD(2);
t544 = sin(qJ(3));
t548 = cos(qJ(3));
t545 = sin(qJ(2));
t540 = sin(pkin(6));
t568 = qJD(1) * t540;
t564 = t545 * t568;
t519 = t537 * t544 + t548 * t564;
t549 = cos(qJ(2));
t567 = qJD(1) * t549;
t563 = t540 * t567;
t534 = qJD(3) - t563;
t543 = sin(qJ(4));
t547 = cos(qJ(4));
t504 = -t519 * t543 + t534 * t547;
t505 = t519 * t547 + t534 * t543;
t542 = sin(qJ(5));
t579 = cos(qJ(5));
t478 = -t504 * t579 + t505 * t542;
t479 = t542 * t504 + t505 * t579;
t518 = t537 * t548 - t544 * t564;
t517 = qJD(4) - t518;
t515 = qJD(5) + t517;
t581 = t585 * t478 - t575 * t479 - t583 * t515;
t580 = -t575 * t478 + t586 * t479 + t584 * t515;
t566 = qJD(1) * qJD(2);
t531 = (-qJDD(1) * t549 + t545 * t566) * t540;
t578 = pkin(8) * t540;
t577 = t541 * g(3);
t576 = -mrSges(6,3) - mrSges(7,2);
t551 = qJD(1) ^ 2;
t546 = sin(qJ(1));
t550 = cos(qJ(1));
t562 = t546 * g(1) - g(2) * t550;
t526 = qJDD(1) * pkin(1) + t551 * t578 + t562;
t574 = t526 * t541;
t573 = t540 * t545;
t572 = t540 * t549;
t529 = (-pkin(2) * t549 - pkin(9) * t545) * t568;
t535 = t537 ^ 2;
t536 = qJDD(1) * t541 + qJDD(2);
t558 = -g(1) * t550 - g(2) * t546;
t527 = -pkin(1) * t551 + qJDD(1) * t578 + t558;
t569 = t549 * t527 + t545 * t574;
t476 = -t535 * pkin(2) + t536 * pkin(9) + (-g(3) * t545 + t529 * t567) * t540 + t569;
t530 = (qJDD(1) * t545 + t549 * t566) * t540;
t477 = t531 * pkin(2) - t530 * pkin(9) - t577 + (-t526 + (pkin(2) * t545 - pkin(9) * t549) * t537 * qJD(1)) * t540;
t444 = t548 * t476 + t544 * t477;
t501 = -pkin(3) * t518 - pkin(10) * t519;
t523 = qJDD(3) + t531;
t532 = t534 ^ 2;
t437 = -pkin(3) * t532 + pkin(10) * t523 + t501 * t518 + t444;
t498 = -g(3) * t572 - t545 * t527 + t549 * t574;
t475 = -t536 * pkin(2) - t535 * pkin(9) + t529 * t564 - t498;
t496 = -t519 * qJD(3) - t544 * t530 + t536 * t548;
t497 = qJD(3) * t518 + t530 * t548 + t536 * t544;
t442 = (-t518 * t534 - t497) * pkin(10) + (t519 * t534 - t496) * pkin(3) + t475;
t423 = -t543 * t437 + t547 * t442;
t462 = qJD(4) * t504 + t497 * t547 + t523 * t543;
t494 = qJDD(4) - t496;
t420 = (t504 * t517 - t462) * pkin(11) + (t504 * t505 + t494) * pkin(4) + t423;
t424 = t547 * t437 + t543 * t442;
t461 = -qJD(4) * t505 - t497 * t543 + t523 * t547;
t484 = pkin(4) * t517 - pkin(11) * t505;
t503 = t504 ^ 2;
t422 = -pkin(4) * t503 + pkin(11) * t461 - t484 * t517 + t424;
t416 = t542 * t420 + t422 * t579;
t433 = qJD(5) * t479 - t461 * t579 + t462 * t542;
t465 = mrSges(6,1) * t515 - mrSges(6,3) * t479;
t489 = qJDD(5) + t494;
t454 = pkin(5) * t478 - qJ(6) * t479;
t513 = t515 ^ 2;
t412 = -pkin(5) * t513 + qJ(6) * t489 + 0.2e1 * qJD(6) * t515 - t454 * t478 + t416;
t466 = -mrSges(7,1) * t515 + mrSges(7,2) * t479;
t565 = m(7) * t412 + t489 * mrSges(7,3) + t515 * t466;
t455 = mrSges(7,1) * t478 - mrSges(7,3) * t479;
t570 = -mrSges(6,1) * t478 - mrSges(6,2) * t479 - t455;
t399 = m(6) * t416 - t489 * mrSges(6,2) + t433 * t576 - t515 * t465 + t478 * t570 + t565;
t415 = t420 * t579 - t542 * t422;
t434 = -t478 * qJD(5) + t542 * t461 + t462 * t579;
t464 = -mrSges(6,2) * t515 - mrSges(6,3) * t478;
t413 = -t489 * pkin(5) - t513 * qJ(6) + t479 * t454 + qJDD(6) - t415;
t463 = -mrSges(7,2) * t478 + mrSges(7,3) * t515;
t559 = -m(7) * t413 + t489 * mrSges(7,1) + t515 * t463;
t401 = m(6) * t415 + t489 * mrSges(6,1) + t434 * t576 + t515 * t464 + t479 * t570 + t559;
t396 = t542 * t399 + t401 * t579;
t480 = -mrSges(5,1) * t504 + mrSges(5,2) * t505;
t482 = -mrSges(5,2) * t517 + mrSges(5,3) * t504;
t392 = m(5) * t423 + mrSges(5,1) * t494 - mrSges(5,3) * t462 - t480 * t505 + t482 * t517 + t396;
t483 = mrSges(5,1) * t517 - mrSges(5,3) * t505;
t560 = t399 * t579 - t401 * t542;
t393 = m(5) * t424 - mrSges(5,2) * t494 + mrSges(5,3) * t461 + t480 * t504 - t483 * t517 + t560;
t390 = -t392 * t543 + t547 * t393;
t500 = -mrSges(4,1) * t518 + mrSges(4,2) * t519;
t507 = mrSges(4,1) * t534 - mrSges(4,3) * t519;
t388 = m(4) * t444 - mrSges(4,2) * t523 + mrSges(4,3) * t496 + t500 * t518 - t507 * t534 + t390;
t443 = -t544 * t476 + t548 * t477;
t436 = -t523 * pkin(3) - t532 * pkin(10) + t519 * t501 - t443;
t425 = -t461 * pkin(4) - t503 * pkin(11) + t505 * t484 + t436;
t418 = -0.2e1 * qJD(6) * t479 + (t478 * t515 - t434) * qJ(6) + (t479 * t515 + t433) * pkin(5) + t425;
t409 = m(7) * t418 + t433 * mrSges(7,1) - t434 * mrSges(7,3) + t478 * t463 - t479 * t466;
t556 = m(6) * t425 + t433 * mrSges(6,1) + mrSges(6,2) * t434 + t478 * t464 + t465 * t479 + t409;
t404 = -m(5) * t436 + t461 * mrSges(5,1) - mrSges(5,2) * t462 + t504 * t482 - t483 * t505 - t556;
t506 = -mrSges(4,2) * t534 + mrSges(4,3) * t518;
t403 = m(4) * t443 + mrSges(4,1) * t523 - mrSges(4,3) * t497 - t500 * t519 + t506 * t534 + t404;
t384 = t544 * t388 + t548 * t403;
t571 = t583 * t478 - t584 * t479 + t582 * t515;
t561 = t548 * t388 - t403 * t544;
t389 = t392 * t547 + t393 * t543;
t555 = -m(4) * t475 + t496 * mrSges(4,1) - mrSges(4,2) * t497 + t518 * t506 - t507 * t519 - t389;
t408 = t434 * mrSges(7,2) + t479 * t455 - t559;
t554 = -mrSges(6,1) * t415 + mrSges(7,1) * t413 + mrSges(6,2) * t416 - mrSges(7,3) * t412 + pkin(5) * t408 - qJ(6) * t565 + t582 * t489 + t581 * t479 + (qJ(6) * t455 - t580) * t478 - t584 * t434 + (qJ(6) * mrSges(7,2) + t583) * t433;
t394 = -mrSges(6,1) * t425 - mrSges(7,1) * t418 + mrSges(7,2) * t412 + mrSges(6,3) * t416 - pkin(5) * t409 - t585 * t433 + t575 * t434 + t571 * t479 + t583 * t489 + t580 * t515;
t395 = mrSges(6,2) * t425 + mrSges(7,2) * t413 - mrSges(6,3) * t415 - mrSges(7,3) * t418 - qJ(6) * t409 - t575 * t433 + t586 * t434 + t571 * t478 + t584 * t489 + t581 * t515;
t467 = Ifges(5,5) * t505 + Ifges(5,6) * t504 + Ifges(5,3) * t517;
t469 = Ifges(5,1) * t505 + Ifges(5,4) * t504 + Ifges(5,5) * t517;
t381 = -mrSges(5,1) * t436 + mrSges(5,3) * t424 + Ifges(5,4) * t462 + Ifges(5,2) * t461 + Ifges(5,6) * t494 - pkin(4) * t556 + pkin(11) * t560 + t394 * t579 + t542 * t395 - t505 * t467 + t517 * t469;
t468 = Ifges(5,4) * t505 + Ifges(5,2) * t504 + Ifges(5,6) * t517;
t382 = mrSges(5,2) * t436 - mrSges(5,3) * t423 + Ifges(5,1) * t462 + Ifges(5,4) * t461 + Ifges(5,5) * t494 - pkin(11) * t396 - t542 * t394 + t395 * t579 + t504 * t467 - t517 * t468;
t491 = Ifges(4,4) * t519 + Ifges(4,2) * t518 + Ifges(4,6) * t534;
t492 = Ifges(4,1) * t519 + Ifges(4,4) * t518 + Ifges(4,5) * t534;
t553 = mrSges(4,1) * t443 - mrSges(4,2) * t444 + Ifges(4,5) * t497 + Ifges(4,6) * t496 + Ifges(4,3) * t523 + pkin(3) * t404 + pkin(10) * t390 + t547 * t381 + t543 * t382 + t519 * t491 - t518 * t492;
t552 = mrSges(5,1) * t423 - mrSges(5,2) * t424 + Ifges(5,5) * t462 + Ifges(5,6) * t461 + Ifges(5,3) * t494 + pkin(4) * t396 + t505 * t468 - t504 * t469 - t554;
t528 = (-mrSges(3,1) * t549 + mrSges(3,2) * t545) * t568;
t525 = -mrSges(3,2) * t537 + mrSges(3,3) * t563;
t524 = mrSges(3,1) * t537 - mrSges(3,3) * t564;
t511 = -t540 * t526 - t577;
t510 = Ifges(3,5) * t537 + (Ifges(3,1) * t545 + Ifges(3,4) * t549) * t568;
t509 = Ifges(3,6) * t537 + (Ifges(3,4) * t545 + Ifges(3,2) * t549) * t568;
t508 = Ifges(3,3) * t537 + (Ifges(3,5) * t545 + Ifges(3,6) * t549) * t568;
t499 = -g(3) * t573 + t569;
t490 = Ifges(4,5) * t519 + Ifges(4,6) * t518 + Ifges(4,3) * t534;
t385 = m(3) * t498 + mrSges(3,1) * t536 - mrSges(3,3) * t530 + t525 * t537 - t528 * t564 + t555;
t383 = m(3) * t499 - mrSges(3,2) * t536 - mrSges(3,3) * t531 - t524 * t537 + t528 * t563 + t561;
t380 = -mrSges(4,1) * t475 + mrSges(4,3) * t444 + Ifges(4,4) * t497 + Ifges(4,2) * t496 + Ifges(4,6) * t523 - pkin(3) * t389 - t519 * t490 + t534 * t492 - t552;
t379 = mrSges(4,2) * t475 - mrSges(4,3) * t443 + Ifges(4,1) * t497 + Ifges(4,4) * t496 + Ifges(4,5) * t523 - pkin(10) * t389 - t381 * t543 + t382 * t547 + t490 * t518 - t491 * t534;
t378 = Ifges(3,5) * t530 - Ifges(3,6) * t531 + Ifges(3,3) * t536 + mrSges(3,1) * t498 - mrSges(3,2) * t499 + t544 * t379 + t548 * t380 + pkin(2) * t555 + pkin(9) * t561 + (t509 * t545 - t510 * t549) * t568;
t1 = [Ifges(2,3) * qJDD(1) + mrSges(2,1) * t562 - mrSges(2,2) * t558 + (mrSges(3,2) * t511 - mrSges(3,3) * t498 + Ifges(3,1) * t530 - Ifges(3,4) * t531 + Ifges(3,5) * t536 - pkin(9) * t384 + t379 * t548 - t380 * t544 + t508 * t563 - t509 * t537) * t573 + (-mrSges(3,1) * t511 + mrSges(3,3) * t499 + Ifges(3,4) * t530 - Ifges(3,2) * t531 + Ifges(3,6) * t536 - pkin(2) * t384 - t508 * t564 + t537 * t510 - t553) * t572 + t541 * t378 + pkin(1) * ((t383 * t545 + t385 * t549) * t541 + (-m(3) * t511 - t531 * mrSges(3,1) - t530 * mrSges(3,2) + (-t524 * t545 + t525 * t549) * t568 - t384) * t540) + (t383 * t549 - t385 * t545) * t578; t378; t553; t552; -t554; t408;];
tauJ  = t1;
