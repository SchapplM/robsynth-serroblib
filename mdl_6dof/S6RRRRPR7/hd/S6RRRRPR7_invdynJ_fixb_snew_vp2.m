% Calculate vector of inverse dynamics joint torques for with Newton-Euler
% S6RRRRPR7
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
%   pkin=[a2,a3,a4,a5,a6,alpha2,d1,d2,d3,d4,d6,theta5]';
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
% Datum: 2019-05-07 21:26
% Revision: 89c353f7eff3bd693eda4e29f35b2761dbc3ada0 (2019-05-03)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauJ = S6RRRRPR7_invdynJ_fixb_snew_vp2(qJ, qJD, qJDD, g, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(6,1),zeros(3,1),zeros(12,1),zeros(7,1),zeros(7,3),zeros(7,6)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRRRPR7_invdynJ_fixb_snew_vp2: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RRRRPR7_invdynJ_fixb_snew_vp2: qJD has to be [6x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [6 1]), ...
  'S6RRRRPR7_invdynJ_fixb_snew_vp2: qJDD has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6RRRRPR7_invdynJ_fixb_snew_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [12 1]), ...
  'S6RRRRPR7_invdynJ_fixb_snew_vp2: pkin has to be [12x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6RRRRPR7_invdynJ_fixb_snew_vp2: m has to be [7x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [7,3]), ...
  'S6RRRRPR7_invdynJ_fixb_snew_vp2: mrSges has to be [7x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [7 6]), ...
  'S6RRRRPR7_invdynJ_fixb_snew_vp2: Ifges has to be [7x6] (double)'); 

%% Symbolic Calculation
% From invdyn_fixb_NewtonEuler_linkframe_tauJ_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-05-07 21:13:16
% EndTime: 2019-05-07 21:13:41
% DurationCPUTime: 18.61s
% Computational Cost: add. (289923->363), mult. (635267->471), div. (0->0), fcn. (521486->14), ass. (0->154)
t592 = 2 * qJD(5);
t554 = sin(pkin(6));
t591 = pkin(8) * t554;
t556 = cos(pkin(6));
t590 = g(3) * t556;
t567 = qJD(1) ^ 2;
t561 = sin(qJ(1));
t566 = cos(qJ(1));
t580 = t561 * g(1) - g(2) * t566;
t538 = qJDD(1) * pkin(1) + t567 * t591 + t580;
t589 = t538 * t556;
t560 = sin(qJ(2));
t588 = t554 * t560;
t565 = cos(qJ(2));
t587 = t554 * t565;
t585 = qJD(1) * t554;
t541 = (-pkin(2) * t565 - pkin(9) * t560) * t585;
t550 = qJD(1) * t556 + qJD(2);
t548 = t550 ^ 2;
t549 = qJDD(1) * t556 + qJDD(2);
t584 = qJD(1) * t565;
t575 = -g(1) * t566 - g(2) * t561;
t583 = qJDD(1) * t554;
t539 = -pkin(1) * t567 + pkin(8) * t583 + t575;
t586 = t565 * t539 + t560 * t589;
t499 = -pkin(2) * t548 + pkin(9) * t549 + (-g(3) * t560 + t541 * t584) * t554 + t586;
t542 = (qJD(2) * t584 + qJDD(1) * t560) * t554;
t582 = t560 * t585;
t543 = -qJD(2) * t582 + t565 * t583;
t500 = -pkin(2) * t543 - pkin(9) * t542 - t590 + (-t538 + (pkin(2) * t560 - pkin(9) * t565) * t550 * qJD(1)) * t554;
t559 = sin(qJ(3));
t564 = cos(qJ(3));
t473 = -t499 * t559 + t564 * t500;
t530 = t550 * t564 - t559 * t582;
t511 = qJD(3) * t530 + t542 * t564 + t549 * t559;
t531 = t550 * t559 + t564 * t582;
t535 = qJDD(3) - t543;
t581 = t554 * t584;
t546 = qJD(3) - t581;
t451 = (t530 * t546 - t511) * pkin(10) + (t530 * t531 + t535) * pkin(3) + t473;
t474 = t564 * t499 + t559 * t500;
t510 = -qJD(3) * t531 - t542 * t559 + t549 * t564;
t520 = pkin(3) * t546 - pkin(10) * t531;
t529 = t530 ^ 2;
t461 = -pkin(3) * t529 + pkin(10) * t510 - t520 * t546 + t474;
t558 = sin(qJ(4));
t563 = cos(qJ(4));
t438 = t563 * t451 - t461 * t558;
t515 = t530 * t563 - t531 * t558;
t479 = qJD(4) * t515 + t510 * t558 + t511 * t563;
t516 = t530 * t558 + t531 * t563;
t534 = qJDD(4) + t535;
t545 = qJD(4) + t546;
t434 = (t515 * t545 - t479) * qJ(5) + (t515 * t516 + t534) * pkin(4) + t438;
t439 = t558 * t451 + t563 * t461;
t478 = -qJD(4) * t516 + t510 * t563 - t511 * t558;
t502 = pkin(4) * t545 - qJ(5) * t516;
t514 = t515 ^ 2;
t436 = -pkin(4) * t514 + qJ(5) * t478 - t502 * t545 + t439;
t553 = sin(pkin(12));
t555 = cos(pkin(12));
t492 = t515 * t555 - t516 * t553;
t431 = t553 * t434 + t555 * t436 + t492 * t592;
t457 = t478 * t555 - t479 * t553;
t493 = t515 * t553 + t516 * t555;
t471 = -mrSges(6,1) * t492 + mrSges(6,2) * t493;
t484 = mrSges(6,1) * t545 - mrSges(6,3) * t493;
t472 = -pkin(5) * t492 - pkin(11) * t493;
t544 = t545 ^ 2;
t428 = -pkin(5) * t544 + pkin(11) * t534 + t472 * t492 + t431;
t512 = -g(3) * t587 - t560 * t539 + t565 * t589;
t498 = -pkin(2) * t549 - pkin(9) * t548 + t541 * t582 - t512;
t465 = -pkin(3) * t510 - pkin(10) * t529 + t531 * t520 + t498;
t441 = -pkin(4) * t478 - qJ(5) * t514 + t516 * t502 + qJDD(5) + t465;
t458 = t478 * t553 + t479 * t555;
t432 = (-t492 * t545 - t458) * pkin(11) + (t493 * t545 - t457) * pkin(5) + t441;
t557 = sin(qJ(6));
t562 = cos(qJ(6));
t425 = -t428 * t557 + t432 * t562;
t481 = -t493 * t557 + t545 * t562;
t444 = qJD(6) * t481 + t458 * t562 + t534 * t557;
t456 = qJDD(6) - t457;
t482 = t493 * t562 + t545 * t557;
t462 = -mrSges(7,1) * t481 + mrSges(7,2) * t482;
t491 = qJD(6) - t492;
t463 = -mrSges(7,2) * t491 + mrSges(7,3) * t481;
t421 = m(7) * t425 + mrSges(7,1) * t456 - mrSges(7,3) * t444 - t462 * t482 + t463 * t491;
t426 = t428 * t562 + t432 * t557;
t443 = -qJD(6) * t482 - t458 * t557 + t534 * t562;
t464 = mrSges(7,1) * t491 - mrSges(7,3) * t482;
t422 = m(7) * t426 - mrSges(7,2) * t456 + mrSges(7,3) * t443 + t462 * t481 - t464 * t491;
t576 = -t421 * t557 + t562 * t422;
t407 = m(6) * t431 - mrSges(6,2) * t534 + mrSges(6,3) * t457 + t471 * t492 - t484 * t545 + t576;
t574 = -t434 * t555 + t436 * t553;
t430 = -0.2e1 * qJD(5) * t493 - t574;
t483 = -mrSges(6,2) * t545 + mrSges(6,3) * t492;
t427 = -pkin(5) * t534 - pkin(11) * t544 + (t592 + t472) * t493 + t574;
t573 = -m(7) * t427 + t443 * mrSges(7,1) - mrSges(7,2) * t444 + t481 * t463 - t464 * t482;
t417 = m(6) * t430 + mrSges(6,1) * t534 - mrSges(6,3) * t458 - t471 * t493 + t483 * t545 + t573;
t403 = t553 * t407 + t555 * t417;
t494 = -mrSges(5,1) * t515 + mrSges(5,2) * t516;
t501 = -mrSges(5,2) * t545 + mrSges(5,3) * t515;
t400 = m(5) * t438 + mrSges(5,1) * t534 - mrSges(5,3) * t479 - t494 * t516 + t501 * t545 + t403;
t503 = mrSges(5,1) * t545 - mrSges(5,3) * t516;
t577 = t555 * t407 - t417 * t553;
t401 = m(5) * t439 - mrSges(5,2) * t534 + mrSges(5,3) * t478 + t494 * t515 - t503 * t545 + t577;
t394 = t563 * t400 + t558 * t401;
t517 = -mrSges(4,1) * t530 + mrSges(4,2) * t531;
t518 = -mrSges(4,2) * t546 + mrSges(4,3) * t530;
t392 = m(4) * t473 + mrSges(4,1) * t535 - mrSges(4,3) * t511 - t517 * t531 + t518 * t546 + t394;
t519 = mrSges(4,1) * t546 - mrSges(4,3) * t531;
t578 = -t400 * t558 + t563 * t401;
t393 = m(4) * t474 - mrSges(4,2) * t535 + mrSges(4,3) * t510 + t517 * t530 - t519 * t546 + t578;
t387 = t564 * t392 + t559 * t393;
t410 = t562 * t421 + t557 * t422;
t579 = -t392 * t559 + t564 * t393;
t408 = m(6) * t441 - t457 * mrSges(6,1) + t458 * mrSges(6,2) - t492 * t483 + t493 * t484 + t410;
t572 = m(5) * t465 - t478 * mrSges(5,1) + t479 * mrSges(5,2) - t515 * t501 + t516 * t503 + t408;
t446 = Ifges(7,4) * t482 + Ifges(7,2) * t481 + Ifges(7,6) * t491;
t447 = Ifges(7,1) * t482 + Ifges(7,4) * t481 + Ifges(7,5) * t491;
t571 = mrSges(7,1) * t425 - mrSges(7,2) * t426 + Ifges(7,5) * t444 + Ifges(7,6) * t443 + Ifges(7,3) * t456 + t446 * t482 - t447 * t481;
t445 = Ifges(7,5) * t482 + Ifges(7,6) * t481 + Ifges(7,3) * t491;
t414 = -mrSges(7,1) * t427 + mrSges(7,3) * t426 + Ifges(7,4) * t444 + Ifges(7,2) * t443 + Ifges(7,6) * t456 - t445 * t482 + t447 * t491;
t415 = mrSges(7,2) * t427 - mrSges(7,3) * t425 + Ifges(7,1) * t444 + Ifges(7,4) * t443 + Ifges(7,5) * t456 + t445 * t481 - t446 * t491;
t467 = Ifges(6,4) * t493 + Ifges(6,2) * t492 + Ifges(6,6) * t545;
t468 = Ifges(6,1) * t493 + Ifges(6,4) * t492 + Ifges(6,5) * t545;
t486 = Ifges(5,4) * t516 + Ifges(5,2) * t515 + Ifges(5,6) * t545;
t487 = Ifges(5,1) * t516 + Ifges(5,4) * t515 + Ifges(5,5) * t545;
t570 = -mrSges(5,1) * t438 - mrSges(6,1) * t430 + mrSges(5,2) * t439 + mrSges(6,2) * t431 - Ifges(6,6) * t457 - pkin(4) * t403 - pkin(5) * t573 - pkin(11) * t576 - t562 * t414 - t557 * t415 + t492 * t468 + t515 * t487 - Ifges(6,5) * t458 - t493 * t467 - Ifges(5,6) * t478 - Ifges(5,5) * t479 - t516 * t486 + (-Ifges(6,3) - Ifges(5,3)) * t534;
t569 = -m(4) * t498 + t510 * mrSges(4,1) - t511 * mrSges(4,2) + t530 * t518 - t531 * t519 - t572;
t505 = Ifges(4,4) * t531 + Ifges(4,2) * t530 + Ifges(4,6) * t546;
t506 = Ifges(4,1) * t531 + Ifges(4,4) * t530 + Ifges(4,5) * t546;
t568 = mrSges(4,1) * t473 - mrSges(4,2) * t474 + Ifges(4,5) * t511 + Ifges(4,6) * t510 + Ifges(4,3) * t535 + pkin(3) * t394 + t531 * t505 - t530 * t506 - t570;
t540 = (-mrSges(3,1) * t565 + mrSges(3,2) * t560) * t585;
t537 = -mrSges(3,2) * t550 + mrSges(3,3) * t581;
t536 = mrSges(3,1) * t550 - mrSges(3,3) * t582;
t524 = -t538 * t554 - t590;
t523 = Ifges(3,5) * t550 + (Ifges(3,1) * t560 + Ifges(3,4) * t565) * t585;
t522 = Ifges(3,6) * t550 + (Ifges(3,4) * t560 + Ifges(3,2) * t565) * t585;
t521 = Ifges(3,3) * t550 + (Ifges(3,5) * t560 + Ifges(3,6) * t565) * t585;
t513 = -g(3) * t588 + t586;
t504 = Ifges(4,5) * t531 + Ifges(4,6) * t530 + Ifges(4,3) * t546;
t485 = Ifges(5,5) * t516 + Ifges(5,6) * t515 + Ifges(5,3) * t545;
t466 = Ifges(6,5) * t493 + Ifges(6,6) * t492 + Ifges(6,3) * t545;
t404 = m(3) * t512 + t549 * mrSges(3,1) - t542 * mrSges(3,3) + t550 * t537 - t540 * t582 + t569;
t396 = -mrSges(6,1) * t441 + mrSges(6,3) * t431 + Ifges(6,4) * t458 + Ifges(6,2) * t457 + Ifges(6,6) * t534 - pkin(5) * t410 - t466 * t493 + t468 * t545 - t571;
t395 = mrSges(6,2) * t441 - mrSges(6,3) * t430 + Ifges(6,1) * t458 + Ifges(6,4) * t457 + Ifges(6,5) * t534 - pkin(11) * t410 - t414 * t557 + t415 * t562 + t466 * t492 - t467 * t545;
t388 = mrSges(5,2) * t465 - mrSges(5,3) * t438 + Ifges(5,1) * t479 + Ifges(5,4) * t478 + Ifges(5,5) * t534 - qJ(5) * t403 + t395 * t555 - t396 * t553 + t485 * t515 - t486 * t545;
t386 = m(3) * t513 - mrSges(3,2) * t549 + mrSges(3,3) * t543 - t536 * t550 + t540 * t581 + t579;
t385 = -mrSges(5,1) * t465 + mrSges(5,3) * t439 + Ifges(5,4) * t479 + Ifges(5,2) * t478 + Ifges(5,6) * t534 - pkin(4) * t408 + qJ(5) * t577 + t553 * t395 + t555 * t396 - t516 * t485 + t545 * t487;
t384 = mrSges(4,2) * t498 - mrSges(4,3) * t473 + Ifges(4,1) * t511 + Ifges(4,4) * t510 + Ifges(4,5) * t535 - pkin(10) * t394 - t385 * t558 + t388 * t563 + t504 * t530 - t505 * t546;
t383 = -mrSges(4,1) * t498 + mrSges(4,3) * t474 + Ifges(4,4) * t511 + Ifges(4,2) * t510 + Ifges(4,6) * t535 - pkin(3) * t572 + pkin(10) * t578 + t563 * t385 + t558 * t388 - t531 * t504 + t546 * t506;
t382 = Ifges(3,5) * t542 + Ifges(3,6) * t543 + Ifges(3,3) * t549 + mrSges(3,1) * t512 - mrSges(3,2) * t513 + t559 * t384 + t564 * t383 + pkin(2) * t569 + pkin(9) * t579 + (t522 * t560 - t523 * t565) * t585;
t1 = [Ifges(2,3) * qJDD(1) + mrSges(2,1) * t580 - mrSges(2,2) * t575 + (mrSges(3,2) * t524 - mrSges(3,3) * t512 + Ifges(3,1) * t542 + Ifges(3,4) * t543 + Ifges(3,5) * t549 - pkin(9) * t387 - t383 * t559 + t384 * t564 + t521 * t581 - t522 * t550) * t588 + (-mrSges(3,1) * t524 + mrSges(3,3) * t513 + Ifges(3,4) * t542 + Ifges(3,2) * t543 + Ifges(3,6) * t549 - pkin(2) * t387 - t521 * t582 + t550 * t523 - t568) * t587 + t556 * t382 + pkin(1) * ((t386 * t560 + t404 * t565) * t556 + (-m(3) * t524 + t543 * mrSges(3,1) - t542 * mrSges(3,2) + (-t536 * t560 + t537 * t565) * t585 - t387) * t554) + (t386 * t565 - t404 * t560) * t591; t382; t568; -t570; t408; t571;];
tauJ  = t1;
