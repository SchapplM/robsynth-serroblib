% Calculate vector of inverse dynamics joint torques for with Newton-Euler
% S6RRRPRR13
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
% pkin [13x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,alpha2,alpha3,d1,d2,d3,d5,d6,theta4]';
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
% Datum: 2019-05-07 15:59
% Revision: 89c353f7eff3bd693eda4e29f35b2761dbc3ada0 (2019-05-03)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauJ = S6RRRPRR13_invdynJ_fixb_snew_vp2(qJ, qJD, qJDD, g, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(6,1),zeros(3,1),zeros(13,1),zeros(7,1),zeros(7,3),zeros(7,6)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRRPRR13_invdynJ_fixb_snew_vp2: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RRRPRR13_invdynJ_fixb_snew_vp2: qJD has to be [6x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [6 1]), ...
  'S6RRRPRR13_invdynJ_fixb_snew_vp2: qJDD has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6RRRPRR13_invdynJ_fixb_snew_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [13 1]), ...
  'S6RRRPRR13_invdynJ_fixb_snew_vp2: pkin has to be [13x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6RRRPRR13_invdynJ_fixb_snew_vp2: m has to be [7x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [7,3]), ...
  'S6RRRPRR13_invdynJ_fixb_snew_vp2: mrSges has to be [7x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [7 6]), ...
  'S6RRRPRR13_invdynJ_fixb_snew_vp2: Ifges has to be [7x6] (double)'); 

%% Symbolic Calculation
% From invdyn_fixb_NewtonEuler_linkframe_tauJ_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-05-07 15:37:33
% EndTime: 2019-05-07 15:38:15
% DurationCPUTime: 37.24s
% Computational Cost: add. (589972->380), mult. (1476127->509), div. (0->0), fcn. (1251541->16), ass. (0->168)
t617 = cos(qJ(3));
t574 = sin(pkin(6));
t616 = pkin(9) * t574;
t573 = sin(pkin(7));
t615 = pkin(10) * t573;
t576 = cos(pkin(7));
t614 = pkin(10) * t576;
t577 = cos(pkin(6));
t613 = g(3) * t577;
t587 = qJD(1) ^ 2;
t582 = sin(qJ(1));
t586 = cos(qJ(1));
t598 = g(1) * t582 - g(2) * t586;
t561 = qJDD(1) * pkin(1) + t587 * t616 + t598;
t612 = t561 * t577;
t580 = sin(qJ(3));
t611 = t573 * t580;
t581 = sin(qJ(2));
t610 = t574 * t581;
t585 = cos(qJ(2));
t609 = t574 * t585;
t608 = t576 * t580;
t570 = qJD(1) * t577 + qJD(2);
t604 = qJD(1) * t585;
t599 = t574 * t604;
t594 = t576 * t599;
t553 = (t570 * t573 + t594) * pkin(10);
t606 = qJD(1) * t574;
t558 = (-pkin(2) * t585 - t581 * t615) * t606;
t603 = qJD(1) * qJD(2);
t564 = (qJDD(1) * t581 + t585 * t603) * t574;
t569 = qJDD(1) * t577 + qJDD(2);
t593 = -g(1) * t586 - g(2) * t582;
t562 = -pkin(1) * t587 + qJDD(1) * t616 + t593;
t595 = -t562 * t581 + t585 * t612;
t605 = qJD(1) * t581;
t511 = -t564 * t614 + pkin(2) * t569 + t553 * t570 + (-g(3) * t585 - t558 * t605) * t574 + t595;
t600 = t574 * t605;
t557 = pkin(2) * t570 - t600 * t614;
t565 = (qJDD(1) * t585 - t581 * t603) * t574;
t592 = t565 * t576 + t569 * t573;
t607 = t562 * t585 + t581 * t612;
t512 = -t557 * t570 + (-g(3) * t581 + t558 * t604) * t574 + t592 * pkin(10) + t607;
t519 = -t564 * t615 - pkin(2) * t565 - t613 + (-t561 + (-t553 * t585 + t557 * t581) * qJD(1)) * t574;
t485 = t511 * t608 + t512 * t617 + t519 * t611;
t602 = t573 * t617;
t543 = -t570 * t602 + t580 * t600 - t594 * t617;
t544 = t570 * t611 + (t581 * t617 + t585 * t608) * t606;
t530 = pkin(3) * t543 - qJ(4) * t544;
t545 = -t565 * t573 + t569 * t576 + qJDD(3);
t555 = t570 * t576 - t573 * t599 + qJD(3);
t552 = t555 ^ 2;
t472 = -pkin(3) * t552 + qJ(4) * t545 - t530 * t543 + t485;
t495 = -t511 * t573 + t519 * t576;
t601 = t576 * t617;
t528 = qJD(3) * t544 + t564 * t580 - t565 * t601 - t569 * t602;
t529 = -t543 * qJD(3) + t564 * t617 + t580 * t592;
t475 = (t543 * t555 - t529) * qJ(4) + (t544 * t555 + t528) * pkin(3) + t495;
t572 = sin(pkin(13));
t575 = cos(pkin(13));
t536 = t544 * t575 + t555 * t572;
t464 = -0.2e1 * qJD(4) * t536 - t472 * t572 + t475 * t575;
t516 = t529 * t575 + t545 * t572;
t535 = -t544 * t572 + t555 * t575;
t461 = (t535 * t543 - t516) * pkin(11) + (t535 * t536 + t528) * pkin(4) + t464;
t465 = 0.2e1 * qJD(4) * t535 + t472 * t575 + t475 * t572;
t515 = -t529 * t572 + t545 * t575;
t523 = pkin(4) * t543 - pkin(11) * t536;
t534 = t535 ^ 2;
t463 = -pkin(4) * t534 + pkin(11) * t515 - t523 * t543 + t465;
t579 = sin(qJ(5));
t584 = cos(qJ(5));
t458 = t461 * t579 + t463 * t584;
t509 = t535 * t584 - t536 * t579;
t510 = t535 * t579 + t536 * t584;
t494 = -pkin(5) * t509 - pkin(12) * t510;
t527 = qJDD(5) + t528;
t542 = qJD(5) + t543;
t541 = t542 ^ 2;
t456 = -pkin(5) * t541 + pkin(12) * t527 + t494 * t509 + t458;
t484 = t511 * t601 - t512 * t580 + t519 * t602;
t471 = -t545 * pkin(3) - t552 * qJ(4) + t530 * t544 + qJDD(4) - t484;
t466 = -t515 * pkin(4) - t534 * pkin(11) + t523 * t536 + t471;
t482 = -qJD(5) * t510 + t515 * t584 - t516 * t579;
t483 = qJD(5) * t509 + t515 * t579 + t516 * t584;
t459 = (-t509 * t542 - t483) * pkin(12) + (t510 * t542 - t482) * pkin(5) + t466;
t578 = sin(qJ(6));
t583 = cos(qJ(6));
t453 = -t456 * t578 + t459 * t583;
t496 = -t510 * t578 + t542 * t583;
t469 = qJD(6) * t496 + t483 * t583 + t527 * t578;
t481 = qJDD(6) - t482;
t497 = t510 * t583 + t542 * t578;
t486 = -mrSges(7,1) * t496 + mrSges(7,2) * t497;
t508 = qJD(6) - t509;
t487 = -mrSges(7,2) * t508 + mrSges(7,3) * t496;
t450 = m(7) * t453 + mrSges(7,1) * t481 - mrSges(7,3) * t469 - t486 * t497 + t487 * t508;
t454 = t456 * t583 + t459 * t578;
t468 = -qJD(6) * t497 - t483 * t578 + t527 * t583;
t488 = mrSges(7,1) * t508 - mrSges(7,3) * t497;
t451 = m(7) * t454 - mrSges(7,2) * t481 + mrSges(7,3) * t468 + t486 * t496 - t488 * t508;
t442 = -t450 * t578 + t451 * t583;
t493 = -mrSges(6,1) * t509 + mrSges(6,2) * t510;
t499 = mrSges(6,1) * t542 - mrSges(6,3) * t510;
t436 = m(6) * t458 - mrSges(6,2) * t527 + mrSges(6,3) * t482 + t493 * t509 - t499 * t542 + t442;
t457 = t461 * t584 - t463 * t579;
t455 = -pkin(5) * t527 - pkin(12) * t541 + t494 * t510 - t457;
t452 = -m(7) * t455 + mrSges(7,1) * t468 - mrSges(7,2) * t469 + t487 * t496 - t488 * t497;
t498 = -mrSges(6,2) * t542 + mrSges(6,3) * t509;
t446 = m(6) * t457 + mrSges(6,1) * t527 - mrSges(6,3) * t483 - t493 * t510 + t498 * t542 + t452;
t433 = t436 * t579 + t446 * t584;
t514 = -mrSges(5,1) * t535 + mrSges(5,2) * t536;
t521 = -mrSges(5,2) * t543 + mrSges(5,3) * t535;
t431 = m(5) * t464 + mrSges(5,1) * t528 - mrSges(5,3) * t516 - t514 * t536 + t521 * t543 + t433;
t522 = mrSges(5,1) * t543 - mrSges(5,3) * t536;
t596 = t436 * t584 - t446 * t579;
t432 = m(5) * t465 - mrSges(5,2) * t528 + mrSges(5,3) * t515 + t514 * t535 - t522 * t543 + t596;
t425 = t431 * t575 + t432 * t572;
t441 = t450 * t583 + t451 * t578;
t531 = mrSges(4,1) * t543 + mrSges(4,2) * t544;
t538 = mrSges(4,1) * t555 - mrSges(4,3) * t544;
t597 = -t431 * t572 + t432 * t575;
t422 = m(4) * t485 - mrSges(4,2) * t545 - mrSges(4,3) * t528 - t531 * t543 - t538 * t555 + t597;
t537 = -mrSges(4,2) * t555 - mrSges(4,3) * t543;
t424 = m(4) * t495 + mrSges(4,1) * t528 + mrSges(4,2) * t529 + t537 * t543 + t538 * t544 + t425;
t590 = m(6) * t466 - mrSges(6,1) * t482 + mrSges(6,2) * t483 - t498 * t509 + t499 * t510 + t441;
t440 = m(5) * t471 - mrSges(5,1) * t515 + mrSges(5,2) * t516 - t521 * t535 + t522 * t536 + t590;
t439 = m(4) * t484 + mrSges(4,1) * t545 - mrSges(4,3) * t529 - t531 * t544 + t537 * t555 - t440;
t413 = t422 * t611 + t424 * t576 + t439 * t602;
t418 = t422 * t617 - t439 * t580;
t414 = t422 * t608 - t424 * t573 + t439 * t601;
t476 = Ifges(7,5) * t497 + Ifges(7,6) * t496 + Ifges(7,3) * t508;
t478 = Ifges(7,1) * t497 + Ifges(7,4) * t496 + Ifges(7,5) * t508;
t443 = -mrSges(7,1) * t455 + mrSges(7,3) * t454 + Ifges(7,4) * t469 + Ifges(7,2) * t468 + Ifges(7,6) * t481 - t476 * t497 + t478 * t508;
t477 = Ifges(7,4) * t497 + Ifges(7,2) * t496 + Ifges(7,6) * t508;
t444 = mrSges(7,2) * t455 - mrSges(7,3) * t453 + Ifges(7,1) * t469 + Ifges(7,4) * t468 + Ifges(7,5) * t481 + t476 * t496 - t477 * t508;
t489 = Ifges(6,5) * t510 + Ifges(6,6) * t509 + Ifges(6,3) * t542;
t490 = Ifges(6,4) * t510 + Ifges(6,2) * t509 + Ifges(6,6) * t542;
t426 = mrSges(6,2) * t466 - mrSges(6,3) * t457 + Ifges(6,1) * t483 + Ifges(6,4) * t482 + Ifges(6,5) * t527 - pkin(12) * t441 - t443 * t578 + t444 * t583 + t489 * t509 - t490 * t542;
t491 = Ifges(6,1) * t510 + Ifges(6,4) * t509 + Ifges(6,5) * t542;
t589 = mrSges(7,1) * t453 - mrSges(7,2) * t454 + Ifges(7,5) * t469 + Ifges(7,6) * t468 + Ifges(7,3) * t481 + t477 * t497 - t478 * t496;
t427 = -mrSges(6,1) * t466 + mrSges(6,3) * t458 + Ifges(6,4) * t483 + Ifges(6,2) * t482 + Ifges(6,6) * t527 - pkin(5) * t441 - t489 * t510 + t491 * t542 - t589;
t500 = Ifges(5,5) * t536 + Ifges(5,6) * t535 + Ifges(5,3) * t543;
t502 = Ifges(5,1) * t536 + Ifges(5,4) * t535 + Ifges(5,5) * t543;
t415 = -mrSges(5,1) * t471 + mrSges(5,3) * t465 + Ifges(5,4) * t516 + Ifges(5,2) * t515 + Ifges(5,6) * t528 - pkin(4) * t590 + pkin(11) * t596 + t579 * t426 + t584 * t427 - t536 * t500 + t543 * t502;
t501 = Ifges(5,4) * t536 + Ifges(5,2) * t535 + Ifges(5,6) * t543;
t416 = mrSges(5,2) * t471 - mrSges(5,3) * t464 + Ifges(5,1) * t516 + Ifges(5,4) * t515 + Ifges(5,5) * t528 - pkin(11) * t433 + t426 * t584 - t427 * t579 + t500 * t535 - t501 * t543;
t524 = Ifges(4,5) * t544 - Ifges(4,6) * t543 + Ifges(4,3) * t555;
t525 = Ifges(4,4) * t544 - Ifges(4,2) * t543 + Ifges(4,6) * t555;
t410 = mrSges(4,2) * t495 - mrSges(4,3) * t484 + Ifges(4,1) * t529 - Ifges(4,4) * t528 + Ifges(4,5) * t545 - qJ(4) * t425 - t415 * t572 + t416 * t575 - t524 * t543 - t525 * t555;
t526 = Ifges(4,1) * t544 - Ifges(4,4) * t543 + Ifges(4,5) * t555;
t588 = mrSges(6,1) * t457 - mrSges(6,2) * t458 + Ifges(6,5) * t483 + Ifges(6,6) * t482 + Ifges(6,3) * t527 + pkin(5) * t452 + pkin(12) * t442 + t583 * t443 + t578 * t444 + t510 * t490 - t509 * t491;
t411 = -pkin(3) * t425 - pkin(4) * t433 - t588 + t555 * t526 - t544 * t524 + Ifges(4,6) * t545 + t535 * t502 - t536 * t501 + Ifges(4,4) * t529 - Ifges(5,6) * t515 - Ifges(5,5) * t516 - mrSges(4,1) * t495 + mrSges(4,3) * t485 - mrSges(5,1) * t464 + mrSges(5,2) * t465 + (-Ifges(5,3) - Ifges(4,2)) * t528;
t591 = pkin(10) * t418 + t410 * t580 + t411 * t617;
t563 = (-mrSges(3,1) * t585 + mrSges(3,2) * t581) * t606;
t560 = -mrSges(3,2) * t570 + mrSges(3,3) * t599;
t559 = mrSges(3,1) * t570 - mrSges(3,3) * t600;
t549 = -t561 * t574 - t613;
t548 = Ifges(3,5) * t570 + (Ifges(3,1) * t581 + Ifges(3,4) * t585) * t606;
t547 = Ifges(3,6) * t570 + (Ifges(3,4) * t581 + Ifges(3,2) * t585) * t606;
t546 = Ifges(3,3) * t570 + (Ifges(3,5) * t581 + Ifges(3,6) * t585) * t606;
t540 = -g(3) * t610 + t607;
t539 = -g(3) * t609 + t595;
t417 = m(3) * t540 - mrSges(3,2) * t569 + mrSges(3,3) * t565 - t559 * t570 + t563 * t599 + t418;
t412 = m(3) * t539 + mrSges(3,1) * t569 - mrSges(3,3) * t564 + t560 * t570 - t563 * t600 + t414;
t409 = mrSges(4,1) * t484 - mrSges(4,2) * t485 + Ifges(4,5) * t529 - Ifges(4,6) * t528 + Ifges(4,3) * t545 - pkin(3) * t440 + qJ(4) * t597 + t575 * t415 + t572 * t416 + t544 * t525 + t543 * t526;
t408 = mrSges(3,1) * t539 - mrSges(3,2) * t540 + Ifges(3,5) * t564 + Ifges(3,6) * t565 + Ifges(3,3) * t569 + pkin(2) * t414 + t576 * t409 + (t547 * t581 - t548 * t585) * t606 + t591 * t573;
t1 = [Ifges(2,3) * qJDD(1) + mrSges(2,1) * t598 - mrSges(2,2) * t593 + (t546 * t599 + mrSges(3,2) * t549 - mrSges(3,3) * t539 + t617 * t410 + Ifges(3,1) * t564 + Ifges(3,4) * t565 + Ifges(3,5) * t569 - t580 * t411 - t570 * t547 + (-t413 * t573 - t414 * t576) * pkin(10)) * t610 + (-mrSges(3,1) * t549 + mrSges(3,3) * t540 + Ifges(3,4) * t564 + Ifges(3,2) * t565 + Ifges(3,6) * t569 - pkin(2) * t413 - t573 * t409 - t546 * t600 + t570 * t548 + t576 * t591) * t609 + t577 * t408 + pkin(1) * ((t412 * t585 + t417 * t581) * t577 + (-m(3) * t549 + t565 * mrSges(3,1) - t564 * mrSges(3,2) + (-t559 * t581 + t560 * t585) * t606 - t413) * t574) + (-t412 * t581 + t417 * t585) * t616; t408; t409; t440; t588; t589;];
tauJ  = t1;
