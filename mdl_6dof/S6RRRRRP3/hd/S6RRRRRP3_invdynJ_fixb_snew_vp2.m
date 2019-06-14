% Calculate vector of inverse dynamics joint torques for with Newton-Euler
% S6RRRRRP3
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
%   pkin=[a2,a3,a4,a5,a6,d1,d2,d3,d4,d5]';
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
% Datum: 2019-05-08 04:45
% Revision: 89c353f7eff3bd693eda4e29f35b2761dbc3ada0 (2019-05-03)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauJ = S6RRRRRP3_invdynJ_fixb_snew_vp2(qJ, qJD, qJDD, g, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(6,1),zeros(3,1),zeros(10,1),zeros(7,1),zeros(7,3),zeros(7,6)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRRRRP3_invdynJ_fixb_snew_vp2: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RRRRRP3_invdynJ_fixb_snew_vp2: qJD has to be [6x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [6 1]), ...
  'S6RRRRRP3_invdynJ_fixb_snew_vp2: qJDD has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6RRRRRP3_invdynJ_fixb_snew_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S6RRRRRP3_invdynJ_fixb_snew_vp2: pkin has to be [10x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6RRRRRP3_invdynJ_fixb_snew_vp2: m has to be [7x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [7,3]), ...
  'S6RRRRRP3_invdynJ_fixb_snew_vp2: mrSges has to be [7x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [7 6]), ...
  'S6RRRRRP3_invdynJ_fixb_snew_vp2: Ifges has to be [7x6] (double)'); 

%% Symbolic Calculation
% From invdyn_fixb_NewtonEuler_linkframe_tauJ_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-05-08 04:38:54
% EndTime: 2019-05-08 04:39:08
% DurationCPUTime: 7.02s
% Computational Cost: add. (76952->328), mult. (153180->402), div. (0->0), fcn. (111752->10), ass. (0->130)
t570 = Ifges(6,4) + Ifges(7,4);
t579 = Ifges(6,2) + Ifges(7,2);
t575 = Ifges(6,6) + Ifges(7,6);
t540 = sin(qJ(3));
t545 = cos(qJ(3));
t546 = cos(qJ(2));
t565 = qJD(1) * t546;
t541 = sin(qJ(2));
t566 = qJD(1) * t541;
t519 = -t540 * t566 + t545 * t565;
t564 = qJD(1) * qJD(2);
t526 = t541 * qJDD(1) + t546 * t564;
t527 = t546 * qJDD(1) - t541 * t564;
t495 = t519 * qJD(3) + t545 * t526 + t540 * t527;
t520 = (t546 * t540 + t541 * t545) * qJD(1);
t536 = qJD(2) + qJD(3);
t539 = sin(qJ(4));
t544 = cos(qJ(4));
t507 = t544 * t520 + t539 * t536;
t535 = qJDD(2) + qJDD(3);
t464 = -t507 * qJD(4) - t539 * t495 + t544 * t535;
t506 = -t539 * t520 + t544 * t536;
t465 = t506 * qJD(4) + t544 * t495 + t539 * t535;
t538 = sin(qJ(5));
t543 = cos(qJ(5));
t477 = t543 * t506 - t538 * t507;
t438 = t477 * qJD(5) + t538 * t464 + t543 * t465;
t478 = t538 * t506 + t543 * t507;
t457 = -t477 * mrSges(7,1) + t478 * mrSges(7,2);
t494 = -t520 * qJD(3) - t540 * t526 + t545 * t527;
t530 = qJD(2) * pkin(2) - pkin(8) * t566;
t537 = t546 ^ 2;
t548 = qJD(1) ^ 2;
t542 = sin(qJ(1));
t547 = cos(qJ(1));
t561 = t542 * g(1) - t547 * g(2);
t555 = -qJDD(1) * pkin(1) - t561;
t496 = -t527 * pkin(2) + t530 * t566 + (-pkin(8) * t537 - pkin(7)) * t548 + t555;
t443 = (-t519 * t536 - t495) * pkin(9) + (t520 * t536 - t494) * pkin(3) + t496;
t556 = -t547 * g(1) - t542 * g(2);
t522 = -t548 * pkin(1) + qJDD(1) * pkin(7) + t556;
t569 = t541 * t522;
t571 = pkin(2) * t548;
t484 = qJDD(2) * pkin(2) - t526 * pkin(8) - t569 + (pkin(8) * t564 + t541 * t571 - g(3)) * t546;
t509 = -t541 * g(3) + t546 * t522;
t485 = t527 * pkin(8) - qJD(2) * t530 - t537 * t571 + t509;
t461 = t540 * t484 + t545 * t485;
t504 = -t519 * pkin(3) - t520 * pkin(9);
t534 = t536 ^ 2;
t447 = -t534 * pkin(3) + t535 * pkin(9) + t519 * t504 + t461;
t426 = t544 * t443 - t539 * t447;
t493 = qJDD(4) - t494;
t515 = qJD(4) - t519;
t423 = (t506 * t515 - t465) * pkin(10) + (t506 * t507 + t493) * pkin(4) + t426;
t427 = t539 * t443 + t544 * t447;
t499 = t515 * pkin(4) - t507 * pkin(10);
t505 = t506 ^ 2;
t425 = -t505 * pkin(4) + t464 * pkin(10) - t515 * t499 + t427;
t417 = t543 * t423 - t538 * t425;
t489 = qJDD(5) + t493;
t513 = qJD(5) + t515;
t412 = -0.2e1 * qJD(6) * t478 + (t477 * t513 - t438) * qJ(6) + (t477 * t478 + t489) * pkin(5) + t417;
t466 = -t513 * mrSges(7,2) + t477 * mrSges(7,3);
t563 = m(7) * t412 + t489 * mrSges(7,1) + t513 * t466;
t409 = -t438 * mrSges(7,3) - t478 * t457 + t563;
t418 = t538 * t423 + t543 * t425;
t437 = -t478 * qJD(5) + t543 * t464 - t538 * t465;
t468 = t513 * pkin(5) - t478 * qJ(6);
t476 = t477 ^ 2;
t414 = -t476 * pkin(5) + t437 * qJ(6) + 0.2e1 * qJD(6) * t477 - t513 * t468 + t418;
t576 = Ifges(6,5) + Ifges(7,5);
t577 = Ifges(6,1) + Ifges(7,1);
t567 = t570 * t477 + t577 * t478 + t576 * t513;
t573 = t579 * t477 + t570 * t478 + t575 * t513;
t574 = Ifges(6,3) + Ifges(7,3);
t578 = mrSges(6,1) * t417 + mrSges(7,1) * t412 - mrSges(6,2) * t418 - mrSges(7,2) * t414 + pkin(5) * t409 + t575 * t437 + t576 * t438 - t567 * t477 + t573 * t478 + t574 * t489;
t458 = -t477 * mrSges(6,1) + t478 * mrSges(6,2);
t467 = -t513 * mrSges(6,2) + t477 * mrSges(6,3);
t401 = m(6) * t417 + t489 * mrSges(6,1) + t513 * t467 + (-t457 - t458) * t478 + (-mrSges(6,3) - mrSges(7,3)) * t438 + t563;
t469 = t513 * mrSges(7,1) - t478 * mrSges(7,3);
t470 = t513 * mrSges(6,1) - t478 * mrSges(6,3);
t562 = m(7) * t414 + t437 * mrSges(7,3) + t477 * t457;
t404 = m(6) * t418 + t437 * mrSges(6,3) + t477 * t458 + (-t469 - t470) * t513 + (-mrSges(6,2) - mrSges(7,2)) * t489 + t562;
t399 = t543 * t401 + t538 * t404;
t472 = Ifges(5,4) * t507 + Ifges(5,2) * t506 + Ifges(5,6) * t515;
t473 = Ifges(5,1) * t507 + Ifges(5,4) * t506 + Ifges(5,5) * t515;
t572 = mrSges(5,1) * t426 - mrSges(5,2) * t427 + Ifges(5,5) * t465 + Ifges(5,6) * t464 + Ifges(5,3) * t493 + pkin(4) * t399 + t507 * t472 - t506 * t473 + t578;
t503 = -t519 * mrSges(4,1) + t520 * mrSges(4,2);
t511 = t536 * mrSges(4,1) - t520 * mrSges(4,3);
t482 = -t506 * mrSges(5,1) + t507 * mrSges(5,2);
t497 = -t515 * mrSges(5,2) + t506 * mrSges(5,3);
t396 = m(5) * t426 + t493 * mrSges(5,1) - t465 * mrSges(5,3) - t507 * t482 + t515 * t497 + t399;
t498 = t515 * mrSges(5,1) - t507 * mrSges(5,3);
t557 = -t538 * t401 + t543 * t404;
t397 = m(5) * t427 - t493 * mrSges(5,2) + t464 * mrSges(5,3) + t506 * t482 - t515 * t498 + t557;
t558 = -t539 * t396 + t544 * t397;
t389 = m(4) * t461 - t535 * mrSges(4,2) + t494 * mrSges(4,3) + t519 * t503 - t536 * t511 + t558;
t460 = t545 * t484 - t540 * t485;
t510 = -t536 * mrSges(4,2) + t519 * mrSges(4,3);
t446 = -t535 * pkin(3) - t534 * pkin(9) + t520 * t504 - t460;
t428 = -t464 * pkin(4) - t505 * pkin(10) + t507 * t499 + t446;
t420 = -t437 * pkin(5) - t476 * qJ(6) + t478 * t468 + qJDD(6) + t428;
t415 = m(7) * t420 - t437 * mrSges(7,1) + t438 * mrSges(7,2) - t477 * t466 + t478 * t469;
t552 = m(6) * t428 - t437 * mrSges(6,1) + t438 * mrSges(6,2) - t477 * t467 + t478 * t470 + t415;
t550 = -m(5) * t446 + t464 * mrSges(5,1) - t465 * mrSges(5,2) + t506 * t497 - t507 * t498 - t552;
t406 = m(4) * t460 + t535 * mrSges(4,1) - t495 * mrSges(4,3) - t520 * t503 + t536 * t510 + t550;
t386 = t540 * t389 + t545 * t406;
t391 = t544 * t396 + t539 * t397;
t568 = -t575 * t477 - t576 * t478 - t574 * t513;
t559 = t545 * t389 - t540 * t406;
t392 = -mrSges(6,1) * t428 + mrSges(6,3) * t418 - mrSges(7,1) * t420 + mrSges(7,3) * t414 - pkin(5) * t415 + qJ(6) * t562 + (-qJ(6) * t469 + t567) * t513 + (-qJ(6) * mrSges(7,2) + t575) * t489 + t568 * t478 + t570 * t438 + t579 * t437;
t398 = mrSges(6,2) * t428 + mrSges(7,2) * t420 - mrSges(6,3) * t417 - mrSges(7,3) * t412 - qJ(6) * t409 + t570 * t437 + t577 * t438 - t568 * t477 + t576 * t489 - t573 * t513;
t471 = Ifges(5,5) * t507 + Ifges(5,6) * t506 + Ifges(5,3) * t515;
t383 = -mrSges(5,1) * t446 + mrSges(5,3) * t427 + Ifges(5,4) * t465 + Ifges(5,2) * t464 + Ifges(5,6) * t493 - pkin(4) * t552 + pkin(10) * t557 + t543 * t392 + t538 * t398 - t507 * t471 + t515 * t473;
t385 = mrSges(5,2) * t446 - mrSges(5,3) * t426 + Ifges(5,1) * t465 + Ifges(5,4) * t464 + Ifges(5,5) * t493 - pkin(10) * t399 - t538 * t392 + t543 * t398 + t506 * t471 - t515 * t472;
t501 = Ifges(4,4) * t520 + Ifges(4,2) * t519 + Ifges(4,6) * t536;
t502 = Ifges(4,1) * t520 + Ifges(4,4) * t519 + Ifges(4,5) * t536;
t554 = mrSges(4,1) * t460 - mrSges(4,2) * t461 + Ifges(4,5) * t495 + Ifges(4,6) * t494 + Ifges(4,3) * t535 + pkin(3) * t550 + pkin(9) * t558 + t544 * t383 + t539 * t385 + t520 * t501 - t519 * t502;
t553 = m(4) * t496 - t494 * mrSges(4,1) + t495 * mrSges(4,2) - t519 * t510 + t520 * t511 + t391;
t529 = -qJD(2) * mrSges(3,2) + mrSges(3,3) * t565;
t528 = qJD(2) * mrSges(3,1) - mrSges(3,3) * t566;
t525 = (-t546 * mrSges(3,1) + t541 * mrSges(3,2)) * qJD(1);
t521 = -t548 * pkin(7) + t555;
t518 = Ifges(3,5) * qJD(2) + (t541 * Ifges(3,1) + t546 * Ifges(3,4)) * qJD(1);
t517 = Ifges(3,6) * qJD(2) + (t541 * Ifges(3,4) + t546 * Ifges(3,2)) * qJD(1);
t508 = -t546 * g(3) - t569;
t500 = Ifges(4,5) * t520 + Ifges(4,6) * t519 + Ifges(4,3) * t536;
t381 = -mrSges(4,1) * t496 + mrSges(4,3) * t461 + Ifges(4,4) * t495 + Ifges(4,2) * t494 + Ifges(4,6) * t535 - pkin(3) * t391 - t520 * t500 + t536 * t502 - t572;
t380 = mrSges(4,2) * t496 - mrSges(4,3) * t460 + Ifges(4,1) * t495 + Ifges(4,4) * t494 + Ifges(4,5) * t535 - pkin(9) * t391 - t539 * t383 + t544 * t385 + t519 * t500 - t536 * t501;
t1 = [Ifges(2,3) * qJDD(1) + mrSges(2,1) * t561 - mrSges(2,2) * t556 + t541 * (mrSges(3,2) * t521 - mrSges(3,3) * t508 + Ifges(3,1) * t526 + Ifges(3,4) * t527 + Ifges(3,5) * qJDD(2) - pkin(8) * t386 - qJD(2) * t517 + t545 * t380 - t540 * t381) + t546 * (-mrSges(3,1) * t521 + mrSges(3,3) * t509 + Ifges(3,4) * t526 + Ifges(3,2) * t527 + Ifges(3,6) * qJDD(2) - pkin(2) * t553 + pkin(8) * t559 + qJD(2) * t518 + t540 * t380 + t545 * t381) + pkin(1) * (-m(3) * t521 + t527 * mrSges(3,1) - t526 * mrSges(3,2) + (-t528 * t541 + t529 * t546) * qJD(1) - t553) + pkin(7) * (t546 * (m(3) * t509 - qJDD(2) * mrSges(3,2) + t527 * mrSges(3,3) - qJD(2) * t528 + t525 * t565 + t559) - t541 * (m(3) * t508 + qJDD(2) * mrSges(3,1) - t526 * mrSges(3,3) + qJD(2) * t529 - t525 * t566 + t386)); t554 + Ifges(3,3) * qJDD(2) + (t541 * t517 - t546 * t518) * qJD(1) + Ifges(3,5) * t526 + Ifges(3,6) * t527 + mrSges(3,1) * t508 - mrSges(3,2) * t509 + pkin(2) * t386; t554; t572; t578; t415;];
tauJ  = t1;
