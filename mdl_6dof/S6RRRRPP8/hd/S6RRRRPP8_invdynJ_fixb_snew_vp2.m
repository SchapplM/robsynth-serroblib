% Calculate vector of inverse dynamics joint torques for with Newton-Euler
% S6RRRRPP8
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
%   pkin=[a2,a3,a4,a5,a6,alpha2,d1,d2,d3,d4]';
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
% Datum: 2019-05-07 19:10
% Revision: 89c353f7eff3bd693eda4e29f35b2761dbc3ada0 (2019-05-03)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauJ = S6RRRRPP8_invdynJ_fixb_snew_vp2(qJ, qJD, qJDD, g, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(6,1),zeros(3,1),zeros(10,1),zeros(7,1),zeros(7,3),zeros(7,6)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRRRPP8_invdynJ_fixb_snew_vp2: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RRRRPP8_invdynJ_fixb_snew_vp2: qJD has to be [6x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [6 1]), ...
  'S6RRRRPP8_invdynJ_fixb_snew_vp2: qJDD has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6RRRRPP8_invdynJ_fixb_snew_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S6RRRRPP8_invdynJ_fixb_snew_vp2: pkin has to be [10x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6RRRRPP8_invdynJ_fixb_snew_vp2: m has to be [7x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [7,3]), ...
  'S6RRRRPP8_invdynJ_fixb_snew_vp2: mrSges has to be [7x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [7 6]), ...
  'S6RRRRPP8_invdynJ_fixb_snew_vp2: Ifges has to be [7x6] (double)'); 

%% Symbolic Calculation
% From invdyn_fixb_NewtonEuler_linkframe_tauJ_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-05-07 18:58:02
% EndTime: 2019-05-07 18:58:11
% DurationCPUTime: 4.63s
% Computational Cost: add. (45783->316), mult. (97086->388), div. (0->0), fcn. (75169->10), ass. (0->132)
t593 = Ifges(5,1) + Ifges(6,1) + Ifges(7,1);
t570 = Ifges(5,4) - Ifges(6,5) - Ifges(7,4);
t569 = Ifges(5,5) + Ifges(6,4) - Ifges(7,5);
t592 = -Ifges(5,2) - Ifges(6,3) - Ifges(7,2);
t568 = -Ifges(6,6) + Ifges(7,6) + Ifges(5,6);
t591 = Ifges(5,3) + Ifges(6,2) + Ifges(7,3);
t540 = sin(pkin(6));
t544 = sin(qJ(2));
t546 = cos(qJ(2));
t571 = qJD(1) * qJD(2);
t530 = (-qJDD(1) * t546 + t544 * t571) * t540;
t573 = qJD(1) * t540;
t528 = (-pkin(2) * t546 - pkin(9) * t544) * t573;
t541 = cos(pkin(6));
t537 = qJD(1) * t541 + qJD(2);
t536 = t537 ^ 2;
t557 = qJDD(1) * t541 + qJDD(2);
t572 = qJD(1) * t546;
t548 = qJD(1) ^ 2;
t545 = sin(qJ(1));
t547 = cos(qJ(1));
t555 = -g(1) * t547 - g(2) * t545;
t582 = pkin(8) * t540;
t526 = -pkin(1) * t548 + qJDD(1) * t582 + t555;
t560 = t545 * g(1) - g(2) * t547;
t525 = qJDD(1) * pkin(1) + t548 * t582 + t560;
t578 = t525 * t541;
t574 = t546 * t526 + t544 * t578;
t472 = t557 * pkin(9) - t536 * pkin(2) + (-g(3) * t544 + t528 * t572) * t540 + t574;
t529 = (qJDD(1) * t544 + t546 * t571) * t540;
t581 = t541 * g(3);
t473 = t530 * pkin(2) - t529 * pkin(9) - t581 + (-t525 + (pkin(2) * t544 - pkin(9) * t546) * t537 * qJD(1)) * t540;
t543 = sin(qJ(3));
t585 = cos(qJ(3));
t441 = t585 * t472 + t543 * t473;
t562 = t544 * t573;
t518 = t585 * t537 - t543 * t562;
t519 = t543 * t537 + t585 * t562;
t502 = -pkin(3) * t518 - pkin(10) * t519;
t522 = qJDD(3) + t530;
t561 = t540 * t572;
t535 = qJD(3) - t561;
t534 = t535 ^ 2;
t437 = -pkin(3) * t534 + pkin(10) * t522 + t502 * t518 + t441;
t576 = t540 * t546;
t499 = -g(3) * t576 - t544 * t526 + t546 * t578;
t471 = -t557 * pkin(2) - t536 * pkin(9) + t528 * t562 - t499;
t497 = -qJD(3) * t519 - t529 * t543 + t585 * t557;
t498 = t518 * qJD(3) + t585 * t529 + t543 * t557;
t439 = (-t518 * t535 - t498) * pkin(10) + (t519 * t535 - t497) * pkin(3) + t471;
t542 = sin(qJ(4));
t584 = cos(qJ(4));
t432 = -t542 * t437 + t584 * t439;
t504 = t542 * t519 - t584 * t535;
t505 = t584 * t519 + t542 * t535;
t476 = pkin(4) * t504 - qJ(5) * t505;
t495 = qJDD(4) - t497;
t516 = qJD(4) - t518;
t515 = t516 ^ 2;
t430 = -t495 * pkin(4) - t515 * qJ(5) + t505 * t476 + qJDD(5) - t432;
t481 = -mrSges(6,2) * t504 + mrSges(6,3) * t516;
t590 = -m(6) * t430 + t495 * mrSges(6,1) + t516 * t481;
t449 = -t504 * qJD(4) + t584 * t498 + t542 * t522;
t440 = -t543 * t472 + t585 * t473;
t552 = t522 * pkin(3) + t534 * pkin(10) - t519 * t502 + t440;
t579 = t504 * t516;
t589 = (-t449 + t579) * qJ(5) - t552;
t482 = mrSges(7,2) * t516 + mrSges(7,3) * t504;
t587 = -0.2e1 * t505;
t423 = qJD(6) * t587 + (-t449 - t579) * qJ(6) + (t504 * t505 - t495) * pkin(5) + t430;
t478 = -mrSges(7,1) * t504 + mrSges(7,2) * t505;
t556 = -m(7) * t423 + t449 * mrSges(7,3) + t505 * t478;
t421 = -t495 * mrSges(7,1) - t516 * t482 - t556;
t477 = mrSges(6,1) * t504 - mrSges(6,3) * t505;
t418 = t449 * mrSges(6,2) + t505 * t477 + t421 - t590;
t433 = t584 * t437 + t542 * t439;
t586 = 2 * qJD(5);
t429 = -pkin(4) * t515 + t495 * qJ(5) - t504 * t476 + t516 * t586 + t433;
t448 = t505 * qJD(4) + t542 * t498 - t584 * t522;
t484 = -pkin(5) * t516 - qJ(6) * t505;
t503 = t504 ^ 2;
t425 = -pkin(5) * t503 + qJ(6) * t448 + 0.2e1 * qJD(6) * t504 + t484 * t516 + t429;
t485 = -mrSges(7,1) * t516 - mrSges(7,3) * t505;
t487 = -mrSges(6,1) * t516 + mrSges(6,2) * t505;
t566 = m(7) * t425 + t448 * mrSges(7,3) + t504 * t478;
t554 = m(6) * t429 + t495 * mrSges(6,3) + t516 * t487 + t566;
t563 = -t570 * t504 + t505 * t593 + t569 * t516;
t564 = t504 * t592 + t505 * t570 + t516 * t568;
t588 = -t568 * t448 + t569 * t449 + t591 * t495 + t563 * t504 + t564 * t505 + mrSges(5,1) * t432 - mrSges(6,1) * t430 - mrSges(7,1) * t423 - mrSges(5,2) * t433 + mrSges(7,2) * t425 + mrSges(6,3) * t429 - pkin(4) * t418 - pkin(5) * t421 + qJ(5) * (-t448 * mrSges(6,2) + t495 * mrSges(7,2) - t504 * t477 + t516 * t485 + t554);
t580 = -mrSges(5,3) - mrSges(6,2);
t577 = t540 * t544;
t483 = -mrSges(5,2) * t516 - mrSges(5,3) * t504;
t575 = -mrSges(5,1) * t504 - mrSges(5,2) * t505 - t477;
t414 = m(5) * t432 + (t482 + t483) * t516 + t575 * t505 + (mrSges(5,1) + mrSges(7,1)) * t495 + t580 * t449 + t556 + t590;
t486 = mrSges(5,1) * t516 - mrSges(5,3) * t505;
t416 = m(5) * t433 + (t485 - t486) * t516 + t575 * t504 + (-mrSges(5,2) + mrSges(7,2)) * t495 + t580 * t448 + t554;
t411 = -t414 * t542 + t584 * t416;
t501 = -mrSges(4,1) * t518 + mrSges(4,2) * t519;
t507 = mrSges(4,1) * t535 - mrSges(4,3) * t519;
t409 = m(4) * t441 - mrSges(4,2) * t522 + mrSges(4,3) * t497 + t501 * t518 - t507 * t535 + t411;
t427 = -t503 * qJ(6) + qJDD(6) + (-pkin(4) - pkin(5)) * t448 + (-pkin(4) * t516 + t484 + t586) * t505 - t589;
t422 = m(7) * t427 - t448 * mrSges(7,1) + t449 * mrSges(7,2) - t504 * t482 + t505 * t485;
t431 = qJD(5) * t587 + (t505 * t516 + t448) * pkin(4) + t589;
t420 = m(6) * t431 + mrSges(6,1) * t448 - t449 * mrSges(6,3) + t481 * t504 - t505 * t487 - t422;
t417 = m(5) * t552 - t448 * mrSges(5,1) - mrSges(5,2) * t449 - t504 * t483 - t486 * t505 - t420;
t506 = -mrSges(4,2) * t535 + mrSges(4,3) * t518;
t413 = m(4) * t440 + mrSges(4,1) * t522 - mrSges(4,3) * t498 - t501 * t519 + t506 * t535 + t417;
t404 = t543 * t409 + t585 * t413;
t565 = t504 * t568 - t505 * t569 - t516 * t591;
t559 = t585 * t409 - t543 * t413;
t410 = t584 * t414 + t542 * t416;
t551 = -m(4) * t471 + t497 * mrSges(4,1) - t498 * mrSges(4,2) + t518 * t506 - t519 * t507 - t410;
t402 = mrSges(5,1) * t552 + mrSges(5,3) * t433 - mrSges(6,1) * t431 + mrSges(6,2) * t429 + mrSges(7,1) * t427 - mrSges(7,3) * t425 + pkin(5) * t422 - qJ(6) * t566 - pkin(4) * t420 + (-qJ(6) * t485 + t563) * t516 + t565 * t505 + (-mrSges(7,2) * qJ(6) + t568) * t495 + t570 * t449 + t592 * t448;
t405 = -mrSges(5,2) * t552 + mrSges(6,2) * t430 + mrSges(7,2) * t427 - mrSges(5,3) * t432 - mrSges(6,3) * t431 - mrSges(7,3) * t423 - qJ(5) * t420 - qJ(6) * t421 - t570 * t448 + t449 * t593 + t569 * t495 + t565 * t504 - t564 * t516;
t492 = Ifges(4,4) * t519 + Ifges(4,2) * t518 + Ifges(4,6) * t535;
t493 = Ifges(4,1) * t519 + Ifges(4,4) * t518 + Ifges(4,5) * t535;
t549 = mrSges(4,1) * t440 - mrSges(4,2) * t441 + Ifges(4,5) * t498 + Ifges(4,6) * t497 + Ifges(4,3) * t522 + pkin(3) * t417 + pkin(10) * t411 + t584 * t402 + t542 * t405 + t519 * t492 - t518 * t493;
t527 = (-mrSges(3,1) * t546 + mrSges(3,2) * t544) * t573;
t524 = -mrSges(3,2) * t537 + mrSges(3,3) * t561;
t523 = mrSges(3,1) * t537 - mrSges(3,3) * t562;
t511 = -t540 * t525 - t581;
t510 = Ifges(3,5) * t537 + (Ifges(3,1) * t544 + Ifges(3,4) * t546) * t573;
t509 = Ifges(3,6) * t537 + (Ifges(3,4) * t544 + Ifges(3,2) * t546) * t573;
t508 = Ifges(3,3) * t537 + (Ifges(3,5) * t544 + Ifges(3,6) * t546) * t573;
t500 = -g(3) * t577 + t574;
t491 = Ifges(4,5) * t519 + Ifges(4,6) * t518 + Ifges(4,3) * t535;
t406 = m(3) * t499 + t557 * mrSges(3,1) - t529 * mrSges(3,3) + t537 * t524 - t527 * t562 + t551;
t403 = m(3) * t500 - t557 * mrSges(3,2) - t530 * mrSges(3,3) - t537 * t523 + t527 * t561 + t559;
t401 = -mrSges(4,1) * t471 + mrSges(4,3) * t441 + Ifges(4,4) * t498 + Ifges(4,2) * t497 + Ifges(4,6) * t522 - pkin(3) * t410 - t519 * t491 + t535 * t493 - t588;
t400 = mrSges(4,2) * t471 - mrSges(4,3) * t440 + Ifges(4,1) * t498 + Ifges(4,4) * t497 + Ifges(4,5) * t522 - pkin(10) * t410 - t542 * t402 + t584 * t405 + t518 * t491 - t535 * t492;
t399 = Ifges(3,5) * t529 - Ifges(3,6) * t530 + Ifges(3,3) * t557 + mrSges(3,1) * t499 - mrSges(3,2) * t500 + t543 * t400 + t585 * t401 + pkin(2) * t551 + pkin(9) * t559 + (t509 * t544 - t510 * t546) * t573;
t1 = [Ifges(2,3) * qJDD(1) + mrSges(2,1) * t560 - mrSges(2,2) * t555 + (mrSges(3,2) * t511 - mrSges(3,3) * t499 + Ifges(3,1) * t529 - Ifges(3,4) * t530 + Ifges(3,5) * t557 - pkin(9) * t404 + t585 * t400 - t543 * t401 + t508 * t561 - t537 * t509) * t577 + (-mrSges(3,1) * t511 + mrSges(3,3) * t500 + Ifges(3,4) * t529 - Ifges(3,2) * t530 + Ifges(3,6) * t557 - pkin(2) * t404 - t508 * t562 + t537 * t510 - t549) * t576 + t541 * t399 + pkin(1) * ((t403 * t544 + t406 * t546) * t541 + (-m(3) * t511 - t530 * mrSges(3,1) - t529 * mrSges(3,2) + (-t523 * t544 + t524 * t546) * t573 - t404) * t540) + (t403 * t546 - t406 * t544) * t582; t399; t549; t588; t418; t422;];
tauJ  = t1;
