% Calculate vector of inverse dynamics base forces with Newton-Euler for
% S6RRRRRR10V2
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
% pkin [6x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,d1,d4,d6]';
% m [7x1]
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
% Datum: 2020-06-19 10:42
% Revision: caa0dbda1e8a16d11faaa29ba3bbef6afcd619f7 (2020-05-25)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauB = S6RRRRRR10V2_invdynB_fixb_snew_vp2(qJ, qJD, qJDD, g, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(6,1),zeros(3,1),zeros(6,1),zeros(7,1),zeros(7,3),zeros(7,6)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRRRRR10V2_invdynB_fixb_snew_vp2: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RRRRRR10V2_invdynB_fixb_snew_vp2: qJD has to be [6x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [6 1]), ...
  'S6RRRRRR10V2_invdynB_fixb_snew_vp2: qJDD has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6RRRRRR10V2_invdynB_fixb_snew_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [6 1]), ...
  'S6RRRRRR10V2_invdynB_fixb_snew_vp2: pkin has to be [6x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6RRRRRR10V2_invdynB_fixb_snew_vp2: m has to be [7x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [7,3]), ...
  'S6RRRRRR10V2_invdynB_fixb_snew_vp2: mrSges has to be [7x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [7 6]), ...
  'S6RRRRRR10V2_invdynB_fixb_snew_vp2: Ifges has to be [7x6] (double)'); 

%% Symbolic Calculation
% From invdyn_fixb_NewtonEuler_linkframe_tauB_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2020-06-19 10:17:39
% EndTime: 2020-06-19 10:18:07
% DurationCPUTime: 22.74s
% Computational Cost: add. (156436->353), mult. (308663->449), div. (0->0), fcn. (235790->12), ass. (0->134)
t567 = sin(qJ(3));
t568 = sin(qJ(2));
t573 = cos(qJ(3));
t574 = cos(qJ(2));
t545 = (t567 * t568 - t573 * t574) * qJD(1);
t569 = sin(qJ(1));
t575 = cos(qJ(1));
t558 = -g(1) * t575 - g(2) * t569;
t576 = qJD(1) ^ 2;
t552 = -pkin(1) * t576 + t558;
t539 = -t574 * g(3) - t568 * t552;
t551 = (-mrSges(3,1) * t574 + mrSges(3,2) * t568) * qJD(1);
t587 = qJD(1) * qJD(2);
t553 = qJDD(1) * t568 + t574 * t587;
t588 = qJD(1) * t574;
t556 = -qJD(2) * mrSges(3,2) + mrSges(3,3) * t588;
t589 = qJD(1) * t568;
t532 = (t568 * t574 * t576 + qJDD(2)) * pkin(2) + t539;
t540 = -t568 * g(3) + t574 * t552;
t533 = (-t574 ^ 2 * t576 - qJD(2) ^ 2) * pkin(2) + t540;
t512 = t567 * t532 + t573 * t533;
t546 = (t567 * t574 + t568 * t573) * qJD(1);
t586 = t568 * t587;
t554 = qJDD(1) * t574 - t586;
t519 = -qJD(3) * t546 - t553 * t567 + t554 * t573;
t527 = mrSges(4,1) * t545 + mrSges(4,2) * t546;
t563 = qJD(2) + qJD(3);
t538 = mrSges(4,1) * t563 - mrSges(4,3) * t546;
t562 = qJDD(2) + qJDD(3);
t529 = pkin(3) * t545 - pkin(5) * t546;
t561 = t563 ^ 2;
t497 = -pkin(3) * t561 + pkin(5) * t562 - t529 * t545 + t512;
t566 = sin(qJ(4));
t572 = cos(qJ(4));
t520 = -qJD(3) * t545 + t553 * t573 + t554 * t567;
t557 = g(1) * t569 - t575 * g(2);
t550 = -qJDD(1) * pkin(1) - t557;
t530 = t550 + (-t554 + t586) * pkin(2);
t578 = (t545 * t563 - t520) * pkin(5) + (t546 * t563 - t519) * pkin(3) + t530;
t481 = t572 * t497 + t566 * t578;
t511 = t532 * t573 - t533 * t567;
t496 = -pkin(3) * t562 - pkin(5) * t561 + t529 * t546 - t511;
t565 = sin(qJ(5));
t571 = cos(qJ(5));
t476 = t571 * t481 + t565 * t496;
t535 = -t546 * t566 + t563 * t572;
t502 = qJD(4) * t535 + t520 * t572 + t562 * t566;
t536 = t546 * t572 + t563 * t566;
t541 = qJD(4) + t545;
t516 = t536 * t571 + t541 * t565;
t518 = qJDD(4) - t519;
t486 = -qJD(5) * t516 - t502 * t565 + t518 * t571;
t515 = -t536 * t565 + t541 * t571;
t498 = -mrSges(6,1) * t515 + mrSges(6,2) * t516;
t501 = -qJD(4) * t536 - t520 * t566 + t562 * t572;
t500 = qJDD(5) - t501;
t534 = qJD(5) - t535;
t506 = mrSges(6,1) * t534 - mrSges(6,3) * t516;
t472 = (-t515 * t516 + t500) * pkin(6) + t476;
t480 = t566 * t497 - t572 * t578;
t487 = qJD(5) * t515 + t502 * t571 + t518 * t565;
t473 = (-t515 * t534 - t487) * pkin(6) + t480;
t564 = sin(qJ(6));
t570 = cos(qJ(6));
t470 = -t472 * t564 + t473 * t570;
t503 = -t516 * t564 + t534 * t570;
t478 = qJD(6) * t503 + t487 * t570 + t500 * t564;
t485 = qJDD(6) - t486;
t504 = t516 * t570 + t534 * t564;
t489 = -mrSges(7,1) * t503 + mrSges(7,2) * t504;
t514 = qJD(6) - t515;
t490 = -mrSges(7,2) * t514 + mrSges(7,3) * t503;
t468 = m(7) * t470 + mrSges(7,1) * t485 - mrSges(7,3) * t478 - t489 * t504 + t490 * t514;
t471 = t472 * t570 + t473 * t564;
t477 = -qJD(6) * t504 - t487 * t564 + t500 * t570;
t491 = mrSges(7,1) * t514 - mrSges(7,3) * t504;
t469 = m(7) * t471 - mrSges(7,2) * t485 + mrSges(7,3) * t477 + t489 * t503 - t491 * t514;
t583 = -t468 * t564 + t570 * t469;
t463 = m(6) * t476 - mrSges(6,2) * t500 + mrSges(6,3) * t486 + t498 * t515 - t506 * t534 + t583;
t475 = -t565 * t481 + t571 * t496;
t474 = (-t516 ^ 2 - t534 ^ 2) * pkin(6) - t475;
t505 = -mrSges(6,2) * t534 + mrSges(6,3) * t515;
t466 = m(6) * t475 - m(7) * t474 + mrSges(6,1) * t500 + mrSges(7,1) * t477 - mrSges(7,2) * t478 - mrSges(6,3) * t487 + t490 * t503 - t491 * t504 - t498 * t516 + t505 * t534;
t513 = -mrSges(5,1) * t535 + mrSges(5,2) * t536;
t522 = mrSges(5,1) * t541 - mrSges(5,3) * t536;
t459 = m(5) * t481 - mrSges(5,2) * t518 + mrSges(5,3) * t501 + t463 * t571 - t466 * t565 + t513 * t535 - t522 * t541;
t521 = -mrSges(5,2) * t541 + mrSges(5,3) * t535;
t582 = t468 * t570 + t469 * t564;
t461 = mrSges(5,1) * t518 + mrSges(6,1) * t486 - mrSges(6,2) * t487 - mrSges(5,3) * t502 + t505 * t515 - t506 * t516 - t513 * t536 + t521 * t541 + (-m(5) - m(6)) * t480 - t582;
t584 = t572 * t459 - t461 * t566;
t451 = m(4) * t512 - mrSges(4,2) * t562 + mrSges(4,3) * t519 - t527 * t545 - t538 * t563 + t584;
t537 = -mrSges(4,2) * t563 - mrSges(4,3) * t545;
t579 = -m(5) * t496 + t501 * mrSges(5,1) - mrSges(5,2) * t502 - t463 * t565 - t466 * t571 + t535 * t521 - t522 * t536;
t456 = m(4) * t511 + mrSges(4,1) * t562 - mrSges(4,3) * t520 - t527 * t546 + t537 * t563 + t579;
t590 = t567 * t451 + t573 * t456;
t445 = m(3) * t539 + qJDD(2) * mrSges(3,1) - mrSges(3,3) * t553 + qJD(2) * t556 - t551 * t589 + t590;
t555 = qJD(2) * mrSges(3,1) - mrSges(3,3) * t589;
t446 = m(3) * t540 - qJDD(2) * mrSges(3,2) + mrSges(3,3) * t554 - qJD(2) * t555 + t451 * t573 - t456 * t567 + t551 * t588;
t440 = m(2) * t558 - mrSges(2,1) * t576 - qJDD(1) * mrSges(2,2) - t445 * t568 + t446 * t574;
t452 = t566 * t459 + t572 * t461;
t580 = m(4) * t530 - t519 * mrSges(4,1) + mrSges(4,2) * t520 + t545 * t537 + t538 * t546 + t452;
t577 = -m(3) * t550 + t554 * mrSges(3,1) - mrSges(3,2) * t553 - t555 * t589 + t556 * t588 - t580;
t449 = m(2) * t557 + qJDD(1) * mrSges(2,1) - mrSges(2,2) * t576 + t577;
t592 = t569 * t440 + t575 * t449;
t591 = t574 * t445 + t568 * t446;
t585 = t575 * t440 - t449 * t569;
t544 = Ifges(3,5) * qJD(2) + (Ifges(3,1) * t568 + Ifges(3,4) * t574) * qJD(1);
t543 = Ifges(3,6) * qJD(2) + (Ifges(3,4) * t568 + Ifges(3,2) * t574) * qJD(1);
t542 = Ifges(3,3) * qJD(2) + (Ifges(3,5) * t568 + Ifges(3,6) * t574) * qJD(1);
t525 = Ifges(4,1) * t546 - Ifges(4,4) * t545 + Ifges(4,5) * t563;
t524 = Ifges(4,4) * t546 - Ifges(4,2) * t545 + Ifges(4,6) * t563;
t523 = Ifges(4,5) * t546 - Ifges(4,6) * t545 + Ifges(4,3) * t563;
t509 = Ifges(5,1) * t536 + Ifges(5,4) * t535 + Ifges(5,5) * t541;
t508 = Ifges(5,4) * t536 + Ifges(5,2) * t535 + Ifges(5,6) * t541;
t507 = Ifges(5,5) * t536 + Ifges(5,6) * t535 + Ifges(5,3) * t541;
t494 = Ifges(6,1) * t516 + Ifges(6,4) * t515 + Ifges(6,5) * t534;
t493 = Ifges(6,4) * t516 + Ifges(6,2) * t515 + Ifges(6,6) * t534;
t492 = Ifges(6,5) * t516 + Ifges(6,6) * t515 + Ifges(6,3) * t534;
t484 = Ifges(7,1) * t504 + Ifges(7,4) * t503 + Ifges(7,5) * t514;
t483 = Ifges(7,4) * t504 + Ifges(7,2) * t503 + Ifges(7,6) * t514;
t482 = Ifges(7,5) * t504 + Ifges(7,6) * t503 + Ifges(7,3) * t514;
t465 = mrSges(7,2) * t474 - mrSges(7,3) * t470 + Ifges(7,1) * t478 + Ifges(7,4) * t477 + Ifges(7,5) * t485 + t482 * t503 - t483 * t514;
t464 = -mrSges(7,1) * t474 + mrSges(7,3) * t471 + Ifges(7,4) * t478 + Ifges(7,2) * t477 + Ifges(7,6) * t485 - t482 * t504 + t484 * t514;
t462 = -mrSges(6,1) * t480 - mrSges(7,1) * t470 + mrSges(7,2) * t471 + mrSges(6,3) * t476 + Ifges(6,4) * t487 - Ifges(7,5) * t478 + Ifges(6,2) * t486 + Ifges(6,6) * t500 - Ifges(7,6) * t477 - Ifges(7,3) * t485 - t483 * t504 + t484 * t503 - t492 * t516 + t494 * t534;
t454 = mrSges(6,2) * t480 - mrSges(6,3) * t475 + Ifges(6,1) * t487 + Ifges(6,4) * t486 + Ifges(6,5) * t500 - pkin(6) * t582 - t564 * t464 + t570 * t465 + t515 * t492 - t534 * t493;
t453 = Ifges(5,4) * t502 + Ifges(5,2) * t501 + Ifges(5,6) * t518 - t536 * t507 + t541 * t509 - mrSges(5,1) * t496 + mrSges(5,3) * t481 - Ifges(6,5) * t487 - Ifges(6,6) * t486 - Ifges(6,3) * t500 - t516 * t493 + t515 * t494 - mrSges(6,1) * t475 + mrSges(6,2) * t476 - t564 * t465 - t570 * t464 - pkin(6) * t583;
t447 = mrSges(5,2) * t496 + mrSges(5,3) * t480 + Ifges(5,1) * t502 + Ifges(5,4) * t501 + Ifges(5,5) * t518 + t454 * t571 - t462 * t565 + t507 * t535 - t508 * t541;
t442 = Ifges(4,4) * t520 + Ifges(4,2) * t519 + Ifges(4,6) * t562 - t546 * t523 + t563 * t525 - mrSges(4,1) * t530 + mrSges(4,3) * t512 - Ifges(5,5) * t502 - Ifges(5,6) * t501 - Ifges(5,3) * t518 - t536 * t508 + t535 * t509 + mrSges(5,1) * t480 + mrSges(5,2) * t481 - t565 * t454 - t571 * t462 - pkin(3) * t452;
t441 = mrSges(4,2) * t530 - mrSges(4,3) * t511 + Ifges(4,1) * t520 + Ifges(4,4) * t519 + Ifges(4,5) * t562 - pkin(5) * t452 + t447 * t572 - t453 * t566 - t523 * t545 - t524 * t563;
t437 = mrSges(3,2) * t550 - mrSges(3,3) * t539 + Ifges(3,1) * t553 + Ifges(3,4) * t554 + Ifges(3,5) * qJDD(2) - qJD(2) * t543 + t441 * t573 - t442 * t567 + t542 * t588;
t436 = -mrSges(3,1) * t550 + mrSges(3,3) * t540 + Ifges(3,4) * t553 + Ifges(3,2) * t554 + Ifges(3,6) * qJDD(2) - pkin(2) * t580 + qJD(2) * t544 + t567 * t441 + t573 * t442 - t542 * t589;
t435 = -pkin(1) * t591 + mrSges(2,3) * t558 - pkin(2) * t590 + mrSges(3,2) * t540 - Ifges(3,5) * t553 - Ifges(3,6) * t554 - Ifges(3,3) * qJDD(2) - mrSges(3,1) * t539 - t566 * t447 - t572 * t453 - pkin(3) * t579 - pkin(5) * t584 - Ifges(4,5) * t520 - Ifges(4,6) * t519 - Ifges(4,3) * t562 - mrSges(4,1) * t511 + mrSges(4,2) * t512 + mrSges(2,1) * g(3) - t546 * t524 - t545 * t525 + t576 * Ifges(2,5) + Ifges(2,6) * qJDD(1) + (-t543 * t568 + t544 * t574) * qJD(1);
t434 = -mrSges(2,2) * g(3) - mrSges(2,3) * t557 + Ifges(2,5) * qJDD(1) - Ifges(2,6) * t576 - t436 * t568 + t437 * t574;
t1 = [-m(1) * g(1) + t585; -m(1) * g(2) + t592; (-m(1) - m(2)) * g(3) + t591; -mrSges(1,2) * g(3) + mrSges(1,3) * g(2) - pkin(4) * t592 + t575 * t434 - t569 * t435; mrSges(1,1) * g(3) - mrSges(1,3) * g(1) + pkin(4) * t585 + t569 * t434 + t575 * t435; -mrSges(1,1) * g(2) + mrSges(2,1) * t557 + mrSges(1,2) * g(1) - mrSges(2,2) * t558 + Ifges(2,3) * qJDD(1) + pkin(1) * t577 + t574 * t436 + t568 * t437;];
tauB = t1;
