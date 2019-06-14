% Calculate vector of inverse dynamics base forces with Newton-Euler for
% S6PPRRPR1
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
%   pkin=[a2,a3,a4,a5,a6,alpha2,alpha3,d3,d4,d6,theta1,theta2,theta5]';
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
% tauB [6x1]
%   base forces of inverse dynamics (contains inertial, gravitational coriolis and centrifugal forces)

% Quelle: HybrDyn-Toolbox
% Datum: 2019-05-04 20:07
% Revision: 89c353f7eff3bd693eda4e29f35b2761dbc3ada0 (2019-05-03)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauB = S6PPRRPR1_invdynB_fixb_snew_vp2(qJ, qJD, qJDD, g, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(6,1),zeros(3,1),zeros(13,1),zeros(7,1),zeros(7,3),zeros(7,6)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6PPRRPR1_invdynB_fixb_snew_vp2: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6PPRRPR1_invdynB_fixb_snew_vp2: qJD has to be [6x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [6 1]), ...
  'S6PPRRPR1_invdynB_fixb_snew_vp2: qJDD has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6PPRRPR1_invdynB_fixb_snew_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [13 1]), ...
  'S6PPRRPR1_invdynB_fixb_snew_vp2: pkin has to be [13x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6PPRRPR1_invdynB_fixb_snew_vp2: m has to be [7x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [7,3]), ...
  'S6PPRRPR1_invdynB_fixb_snew_vp2: mrSges has to be [7x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [7 6]), ...
  'S6PPRRPR1_invdynB_fixb_snew_vp2: Ifges has to be [7x6] (double)'); 

%% Symbolic Calculation
% From invdyn_fixb_NewtonEuler_linkframe_tauB_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-05-04 20:04:09
% EndTime: 2019-05-04 20:04:31
% DurationCPUTime: 21.77s
% Computational Cost: add. (375260->293), mult. (710761->387), div. (0->0), fcn. (554702->16), ass. (0->132)
t574 = sin(pkin(11));
t579 = cos(pkin(11));
t566 = -g(1) * t579 - g(2) * t574;
t573 = sin(pkin(12));
t578 = cos(pkin(12));
t565 = g(1) * t574 - g(2) * t579;
t571 = -g(3) + qJDD(1);
t576 = sin(pkin(6));
t581 = cos(pkin(6));
t595 = t565 * t581 + t571 * t576;
t532 = -t573 * t566 + t578 * t595;
t533 = t578 * t566 + t573 * t595;
t543 = -t565 * t576 + t571 * t581 + qJDD(2);
t587 = cos(qJ(3));
t580 = cos(pkin(7));
t584 = sin(qJ(3));
t606 = t580 * t584;
t575 = sin(pkin(7));
t607 = t575 * t584;
t516 = t532 * t606 + t587 * t533 + t543 * t607;
t589 = qJD(3) ^ 2;
t514 = -pkin(3) * t589 + qJDD(3) * pkin(9) + t516;
t526 = -t532 * t575 + t543 * t580;
t583 = sin(qJ(4));
t586 = cos(qJ(4));
t507 = t586 * t514 + t583 * t526;
t562 = (-mrSges(5,1) * t586 + mrSges(5,2) * t583) * qJD(3);
t602 = qJD(3) * qJD(4);
t570 = t583 * t602;
t564 = qJDD(3) * t586 - t570;
t604 = qJD(3) * t583;
t567 = qJD(4) * mrSges(5,1) - mrSges(5,3) * t604;
t561 = (-pkin(4) * t586 - qJ(5) * t583) * qJD(3);
t588 = qJD(4) ^ 2;
t603 = qJD(3) * t586;
t505 = -pkin(4) * t588 + qJDD(4) * qJ(5) + t561 * t603 + t507;
t515 = -t584 * t533 + (t532 * t580 + t543 * t575) * t587;
t513 = -qJDD(3) * pkin(3) - t589 * pkin(9) - t515;
t601 = t586 * t602;
t563 = qJDD(3) * t583 + t601;
t510 = (-t563 - t601) * qJ(5) + (-t564 + t570) * pkin(4) + t513;
t572 = sin(pkin(13));
t577 = cos(pkin(13));
t558 = qJD(4) * t572 + t577 * t604;
t500 = -0.2e1 * qJD(5) * t558 - t572 * t505 + t577 * t510;
t547 = qJDD(4) * t572 + t563 * t577;
t557 = qJD(4) * t577 - t572 * t604;
t498 = (-t557 * t603 - t547) * pkin(10) + (t557 * t558 - t564) * pkin(5) + t500;
t501 = 0.2e1 * qJD(5) * t557 + t577 * t505 + t572 * t510;
t546 = qJDD(4) * t577 - t563 * t572;
t548 = -pkin(5) * t603 - pkin(10) * t558;
t556 = t557 ^ 2;
t499 = -pkin(5) * t556 + pkin(10) * t546 + t548 * t603 + t501;
t582 = sin(qJ(6));
t585 = cos(qJ(6));
t496 = t498 * t585 - t499 * t582;
t538 = t557 * t585 - t558 * t582;
t519 = qJD(6) * t538 + t546 * t582 + t547 * t585;
t539 = t557 * t582 + t558 * t585;
t525 = -mrSges(7,1) * t538 + mrSges(7,2) * t539;
t569 = qJD(6) - t603;
t528 = -mrSges(7,2) * t569 + mrSges(7,3) * t538;
t560 = qJDD(6) - t564;
t492 = m(7) * t496 + mrSges(7,1) * t560 - mrSges(7,3) * t519 - t525 * t539 + t528 * t569;
t497 = t498 * t582 + t499 * t585;
t518 = -qJD(6) * t539 + t546 * t585 - t547 * t582;
t529 = mrSges(7,1) * t569 - mrSges(7,3) * t539;
t493 = m(7) * t497 - mrSges(7,2) * t560 + mrSges(7,3) * t518 + t525 * t538 - t529 * t569;
t486 = t585 * t492 + t582 * t493;
t540 = -mrSges(6,1) * t557 + mrSges(6,2) * t558;
t544 = mrSges(6,2) * t603 + mrSges(6,3) * t557;
t484 = m(6) * t500 - mrSges(6,1) * t564 - mrSges(6,3) * t547 - t540 * t558 - t544 * t603 + t486;
t545 = -mrSges(6,1) * t603 - mrSges(6,3) * t558;
t597 = -t492 * t582 + t585 * t493;
t485 = m(6) * t501 + mrSges(6,2) * t564 + mrSges(6,3) * t546 + t540 * t557 + t545 * t603 + t597;
t598 = -t484 * t572 + t577 * t485;
t481 = m(5) * t507 - qJDD(4) * mrSges(5,2) + mrSges(5,3) * t564 - qJD(4) * t567 + t562 * t603 + t598;
t506 = -t583 * t514 + t526 * t586;
t568 = -qJD(4) * mrSges(5,2) + mrSges(5,3) * t603;
t504 = -qJDD(4) * pkin(4) - qJ(5) * t588 + t561 * t604 + qJDD(5) - t506;
t502 = -pkin(5) * t546 - pkin(10) * t556 + t548 * t558 + t504;
t592 = m(7) * t502 - t518 * mrSges(7,1) + mrSges(7,2) * t519 - t538 * t528 + t529 * t539;
t590 = -m(6) * t504 + t546 * mrSges(6,1) - mrSges(6,2) * t547 + t557 * t544 - t545 * t558 - t592;
t495 = m(5) * t506 + qJDD(4) * mrSges(5,1) - mrSges(5,3) * t563 + qJD(4) * t568 - t562 * t604 + t590;
t599 = t586 * t481 - t495 * t583;
t470 = m(4) * t516 - mrSges(4,1) * t589 - qJDD(3) * mrSges(4,2) + t599;
t473 = t583 * t481 + t586 * t495;
t472 = m(4) * t526 + t473;
t482 = t484 * t577 + t485 * t572;
t591 = -m(5) * t513 + t564 * mrSges(5,1) - mrSges(5,2) * t563 - t567 * t604 + t568 * t603 - t482;
t478 = m(4) * t515 + qJDD(3) * mrSges(4,1) - mrSges(4,2) * t589 + t591;
t608 = t478 * t587;
t459 = t470 * t606 - t472 * t575 + t580 * t608;
t455 = m(3) * t532 + t459;
t465 = t587 * t470 - t478 * t584;
t464 = m(3) * t533 + t465;
t611 = t455 * t578 + t464 * t573;
t458 = t470 * t607 + t580 * t472 + t575 * t608;
t457 = m(3) * t543 + t458;
t445 = -t457 * t576 + t611 * t581;
t443 = m(2) * t565 + t445;
t451 = -t455 * t573 + t578 * t464;
t450 = m(2) * t566 + t451;
t605 = t579 * t443 + t574 * t450;
t444 = t581 * t457 + t611 * t576;
t600 = -t443 * t574 + t579 * t450;
t520 = Ifges(7,5) * t539 + Ifges(7,6) * t538 + Ifges(7,3) * t569;
t522 = Ifges(7,1) * t539 + Ifges(7,4) * t538 + Ifges(7,5) * t569;
t487 = -mrSges(7,1) * t502 + mrSges(7,3) * t497 + Ifges(7,4) * t519 + Ifges(7,2) * t518 + Ifges(7,6) * t560 - t520 * t539 + t522 * t569;
t521 = Ifges(7,4) * t539 + Ifges(7,2) * t538 + Ifges(7,6) * t569;
t488 = mrSges(7,2) * t502 - mrSges(7,3) * t496 + Ifges(7,1) * t519 + Ifges(7,4) * t518 + Ifges(7,5) * t560 + t520 * t538 - t521 * t569;
t534 = Ifges(6,5) * t558 + Ifges(6,6) * t557 - Ifges(6,3) * t603;
t536 = Ifges(6,1) * t558 + Ifges(6,4) * t557 - Ifges(6,5) * t603;
t474 = -mrSges(6,1) * t504 + mrSges(6,3) * t501 + Ifges(6,4) * t547 + Ifges(6,2) * t546 - Ifges(6,6) * t564 - pkin(5) * t592 + pkin(10) * t597 + t585 * t487 + t582 * t488 - t558 * t534 - t536 * t603;
t535 = Ifges(6,4) * t558 + Ifges(6,2) * t557 - Ifges(6,6) * t603;
t475 = mrSges(6,2) * t504 - mrSges(6,3) * t500 + Ifges(6,1) * t547 + Ifges(6,4) * t546 - Ifges(6,5) * t564 - pkin(10) * t486 - t487 * t582 + t488 * t585 + t534 * t557 + t535 * t603;
t552 = Ifges(5,3) * qJD(4) + (Ifges(5,5) * t583 + Ifges(5,6) * t586) * qJD(3);
t553 = Ifges(5,6) * qJD(4) + (Ifges(5,4) * t583 + Ifges(5,2) * t586) * qJD(3);
t460 = mrSges(5,2) * t513 - mrSges(5,3) * t506 + Ifges(5,1) * t563 + Ifges(5,4) * t564 + Ifges(5,5) * qJDD(4) - qJ(5) * t482 - qJD(4) * t553 - t474 * t572 + t475 * t577 + t552 * t603;
t554 = Ifges(5,5) * qJD(4) + (Ifges(5,1) * t583 + Ifges(5,4) * t586) * qJD(3);
t466 = Ifges(5,4) * t563 + Ifges(5,6) * qJDD(4) - t552 * t604 + qJD(4) * t554 - mrSges(5,1) * t513 + mrSges(5,3) * t507 - Ifges(6,5) * t547 - Ifges(6,6) * t546 - t558 * t535 + t557 * t536 - mrSges(6,1) * t500 + mrSges(6,2) * t501 - Ifges(7,5) * t519 - Ifges(7,6) * t518 - Ifges(7,3) * t560 - t539 * t521 + t538 * t522 - mrSges(7,1) * t496 + mrSges(7,2) * t497 - pkin(5) * t486 - pkin(4) * t482 + (Ifges(5,2) + Ifges(6,3)) * t564;
t447 = mrSges(4,2) * t526 - mrSges(4,3) * t515 + Ifges(4,5) * qJDD(3) - Ifges(4,6) * t589 - pkin(9) * t473 + t460 * t586 - t466 * t583;
t452 = Ifges(4,6) * qJDD(3) + t589 * Ifges(4,5) - mrSges(4,1) * t526 + mrSges(4,3) * t516 - Ifges(5,5) * t563 - Ifges(5,6) * t564 - Ifges(5,3) * qJDD(4) - mrSges(5,1) * t506 + mrSges(5,2) * t507 - t572 * t475 - t577 * t474 - pkin(4) * t590 - qJ(5) * t598 - pkin(3) * t473 + (-t553 * t583 + t554 * t586) * qJD(3);
t594 = pkin(8) * t465 + t447 * t584 + t452 * t587;
t446 = mrSges(4,1) * t515 - mrSges(4,2) * t516 + Ifges(4,3) * qJDD(3) + pkin(3) * t591 + pkin(9) * t599 + t583 * t460 + t586 * t466;
t440 = -mrSges(3,1) * t543 + mrSges(3,3) * t533 - pkin(2) * t458 - t575 * t446 + t580 * t594;
t441 = mrSges(3,2) * t543 - mrSges(3,3) * t532 + t587 * t447 - t584 * t452 + (-t458 * t575 - t459 * t580) * pkin(8);
t593 = qJ(2) * t451 + t440 * t578 + t441 * t573;
t439 = mrSges(3,1) * t532 - mrSges(3,2) * t533 + pkin(2) * t459 + t580 * t446 + t575 * t594;
t438 = mrSges(2,2) * t571 - mrSges(2,3) * t565 - t573 * t440 + t578 * t441 + (-t444 * t576 - t445 * t581) * qJ(2);
t437 = -mrSges(2,1) * t571 + mrSges(2,3) * t566 - pkin(1) * t444 - t576 * t439 + t581 * t593;
t1 = [-m(1) * g(1) + t600; -m(1) * g(2) + t605; -m(1) * g(3) + m(2) * t571 + t444; -mrSges(1,2) * g(3) + mrSges(1,3) * g(2) - qJ(1) * t605 - t574 * t437 + t579 * t438; mrSges(1,1) * g(3) - mrSges(1,3) * g(1) + qJ(1) * t600 + t579 * t437 + t574 * t438; -mrSges(1,1) * g(2) + mrSges(2,1) * t565 + mrSges(1,2) * g(1) - mrSges(2,2) * t566 + pkin(1) * t445 + t581 * t439 + t576 * t593;];
tauB  = t1;
