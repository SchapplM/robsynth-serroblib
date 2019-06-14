% Calculate vector of inverse dynamics joint torques for with Newton-Euler
% S6RRRRRR6
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
%   pkin=[a2,a3,a4,a5,a6,alpha2,d1,d2,d3,d4,d5,d6]';
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
% Datum: 2019-05-08 11:29
% Revision: 89c353f7eff3bd693eda4e29f35b2761dbc3ada0 (2019-05-03)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauJ = S6RRRRRR6_invdynJ_fixb_snew_vp2(qJ, qJD, qJDD, g, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(6,1),zeros(3,1),zeros(12,1),zeros(7,1),zeros(7,3),zeros(7,6)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRRRRR6_invdynJ_fixb_snew_vp2: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RRRRRR6_invdynJ_fixb_snew_vp2: qJD has to be [6x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [6 1]), ...
  'S6RRRRRR6_invdynJ_fixb_snew_vp2: qJDD has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6RRRRRR6_invdynJ_fixb_snew_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [12 1]), ...
  'S6RRRRRR6_invdynJ_fixb_snew_vp2: pkin has to be [12x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6RRRRRR6_invdynJ_fixb_snew_vp2: m has to be [7x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [7,3]), ...
  'S6RRRRRR6_invdynJ_fixb_snew_vp2: mrSges has to be [7x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [7 6]), ...
  'S6RRRRRR6_invdynJ_fixb_snew_vp2: Ifges has to be [7x6] (double)'); 

%% Symbolic Calculation
% From invdyn_fixb_NewtonEuler_linkframe_tauJ_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-05-08 11:10:44
% EndTime: 2019-05-08 11:11:08
% DurationCPUTime: 17.56s
% Computational Cost: add. (289564->364), mult. (611835->470), div. (0->0), fcn. (502797->14), ass. (0->155)
t555 = sin(pkin(6));
t594 = pkin(8) * t555;
t556 = cos(pkin(6));
t593 = t556 * g(3);
t569 = qJD(1) ^ 2;
t562 = sin(qJ(1));
t568 = cos(qJ(1));
t583 = t562 * g(1) - g(2) * t568;
t540 = qJDD(1) * pkin(1) + t569 * t594 + t583;
t592 = t540 * t556;
t561 = sin(qJ(2));
t591 = t555 * t561;
t567 = cos(qJ(2));
t590 = t555 * t567;
t588 = qJD(1) * t555;
t543 = (-pkin(2) * t567 - pkin(9) * t561) * t588;
t552 = qJD(1) * t556 + qJD(2);
t550 = t552 ^ 2;
t551 = qJDD(1) * t556 + qJDD(2);
t587 = qJD(1) * t567;
t578 = -g(1) * t568 - g(2) * t562;
t586 = qJDD(1) * t555;
t541 = -pkin(1) * t569 + pkin(8) * t586 + t578;
t589 = t567 * t541 + t561 * t592;
t496 = -t550 * pkin(2) + t551 * pkin(9) + (-g(3) * t561 + t543 * t587) * t555 + t589;
t544 = (qJD(2) * t587 + qJDD(1) * t561) * t555;
t585 = t561 * t588;
t545 = -qJD(2) * t585 + t567 * t586;
t497 = -t545 * pkin(2) - t544 * pkin(9) - t593 + (-t540 + (pkin(2) * t561 - pkin(9) * t567) * t552 * qJD(1)) * t555;
t560 = sin(qJ(3));
t566 = cos(qJ(3));
t464 = -t560 * t496 + t566 * t497;
t532 = t552 * t566 - t560 * t585;
t511 = qJD(3) * t532 + t544 * t566 + t551 * t560;
t533 = t552 * t560 + t566 * t585;
t537 = qJDD(3) - t545;
t584 = t555 * t587;
t548 = qJD(3) - t584;
t449 = (t532 * t548 - t511) * pkin(10) + (t532 * t533 + t537) * pkin(3) + t464;
t465 = t566 * t496 + t560 * t497;
t510 = -qJD(3) * t533 - t544 * t560 + t551 * t566;
t522 = pkin(3) * t548 - pkin(10) * t533;
t531 = t532 ^ 2;
t456 = -pkin(3) * t531 + pkin(10) * t510 - t522 * t548 + t465;
t559 = sin(qJ(4));
t565 = cos(qJ(4));
t446 = t559 * t449 + t565 * t456;
t518 = t532 * t559 + t533 * t565;
t476 = -t518 * qJD(4) + t510 * t565 - t559 * t511;
t517 = t532 * t565 - t559 * t533;
t490 = -mrSges(5,1) * t517 + mrSges(5,2) * t518;
t547 = qJD(4) + t548;
t502 = mrSges(5,1) * t547 - mrSges(5,3) * t518;
t536 = qJDD(4) + t537;
t491 = -pkin(4) * t517 - pkin(11) * t518;
t546 = t547 ^ 2;
t435 = -pkin(4) * t546 + pkin(11) * t536 + t491 * t517 + t446;
t512 = -g(3) * t590 - t561 * t541 + t567 * t592;
t495 = -t551 * pkin(2) - t550 * pkin(9) + t543 * t585 - t512;
t461 = -t510 * pkin(3) - t531 * pkin(10) + t533 * t522 + t495;
t477 = qJD(4) * t517 + t510 * t559 + t511 * t565;
t438 = (-t517 * t547 - t477) * pkin(11) + (t518 * t547 - t476) * pkin(4) + t461;
t558 = sin(qJ(5));
t564 = cos(qJ(5));
t430 = -t558 * t435 + t564 * t438;
t499 = -t518 * t558 + t547 * t564;
t460 = qJD(5) * t499 + t477 * t564 + t536 * t558;
t475 = qJDD(5) - t476;
t500 = t518 * t564 + t547 * t558;
t516 = qJD(5) - t517;
t428 = (t499 * t516 - t460) * pkin(12) + (t499 * t500 + t475) * pkin(5) + t430;
t431 = t564 * t435 + t558 * t438;
t459 = -qJD(5) * t500 - t477 * t558 + t536 * t564;
t485 = pkin(5) * t516 - pkin(12) * t500;
t498 = t499 ^ 2;
t429 = -pkin(5) * t498 + pkin(12) * t459 - t485 * t516 + t431;
t557 = sin(qJ(6));
t563 = cos(qJ(6));
t426 = t428 * t563 - t429 * t557;
t480 = t499 * t563 - t500 * t557;
t443 = qJD(6) * t480 + t459 * t557 + t460 * t563;
t481 = t499 * t557 + t500 * t563;
t457 = -mrSges(7,1) * t480 + mrSges(7,2) * t481;
t514 = qJD(6) + t516;
t462 = -mrSges(7,2) * t514 + mrSges(7,3) * t480;
t471 = qJDD(6) + t475;
t421 = m(7) * t426 + mrSges(7,1) * t471 - mrSges(7,3) * t443 - t457 * t481 + t462 * t514;
t427 = t428 * t557 + t429 * t563;
t442 = -qJD(6) * t481 + t459 * t563 - t460 * t557;
t463 = mrSges(7,1) * t514 - mrSges(7,3) * t481;
t422 = m(7) * t427 - mrSges(7,2) * t471 + mrSges(7,3) * t442 + t457 * t480 - t463 * t514;
t413 = t563 * t421 + t557 * t422;
t482 = -mrSges(6,1) * t499 + mrSges(6,2) * t500;
t483 = -mrSges(6,2) * t516 + mrSges(6,3) * t499;
t411 = m(6) * t430 + mrSges(6,1) * t475 - mrSges(6,3) * t460 - t482 * t500 + t483 * t516 + t413;
t484 = mrSges(6,1) * t516 - mrSges(6,3) * t500;
t579 = -t421 * t557 + t563 * t422;
t412 = m(6) * t431 - mrSges(6,2) * t475 + mrSges(6,3) * t459 + t482 * t499 - t484 * t516 + t579;
t580 = -t411 * t558 + t564 * t412;
t404 = m(5) * t446 - mrSges(5,2) * t536 + mrSges(5,3) * t476 + t490 * t517 - t502 * t547 + t580;
t445 = t449 * t565 - t559 * t456;
t501 = -mrSges(5,2) * t547 + mrSges(5,3) * t517;
t434 = -pkin(4) * t536 - pkin(11) * t546 + t518 * t491 - t445;
t432 = -pkin(5) * t459 - pkin(12) * t498 + t485 * t500 + t434;
t577 = m(7) * t432 - t442 * mrSges(7,1) + mrSges(7,2) * t443 - t480 * t462 + t463 * t481;
t573 = -m(6) * t434 + t459 * mrSges(6,1) - mrSges(6,2) * t460 + t499 * t483 - t484 * t500 - t577;
t417 = m(5) * t445 + mrSges(5,1) * t536 - mrSges(5,3) * t477 - t490 * t518 + t501 * t547 + t573;
t396 = t559 * t404 + t565 * t417;
t519 = -mrSges(4,1) * t532 + mrSges(4,2) * t533;
t520 = -mrSges(4,2) * t548 + mrSges(4,3) * t532;
t394 = m(4) * t464 + mrSges(4,1) * t537 - mrSges(4,3) * t511 - t519 * t533 + t520 * t548 + t396;
t521 = mrSges(4,1) * t548 - mrSges(4,3) * t533;
t581 = t565 * t404 - t417 * t559;
t395 = m(4) * t465 - mrSges(4,2) * t537 + mrSges(4,3) * t510 + t519 * t532 - t521 * t548 + t581;
t389 = t566 * t394 + t560 * t395;
t406 = t564 * t411 + t558 * t412;
t582 = -t394 * t560 + t566 * t395;
t451 = Ifges(7,4) * t481 + Ifges(7,2) * t480 + Ifges(7,6) * t514;
t452 = Ifges(7,1) * t481 + Ifges(7,4) * t480 + Ifges(7,5) * t514;
t576 = -mrSges(7,1) * t426 + mrSges(7,2) * t427 - Ifges(7,5) * t443 - Ifges(7,6) * t442 - Ifges(7,3) * t471 - t481 * t451 + t480 * t452;
t575 = m(5) * t461 - t476 * mrSges(5,1) + mrSges(5,2) * t477 - t517 * t501 + t502 * t518 + t406;
t450 = Ifges(7,5) * t481 + Ifges(7,6) * t480 + Ifges(7,3) * t514;
t414 = -mrSges(7,1) * t432 + mrSges(7,3) * t427 + Ifges(7,4) * t443 + Ifges(7,2) * t442 + Ifges(7,6) * t471 - t450 * t481 + t452 * t514;
t415 = mrSges(7,2) * t432 - mrSges(7,3) * t426 + Ifges(7,1) * t443 + Ifges(7,4) * t442 + Ifges(7,5) * t471 + t450 * t480 - t451 * t514;
t466 = Ifges(6,5) * t500 + Ifges(6,6) * t499 + Ifges(6,3) * t516;
t468 = Ifges(6,1) * t500 + Ifges(6,4) * t499 + Ifges(6,5) * t516;
t398 = -mrSges(6,1) * t434 + mrSges(6,3) * t431 + Ifges(6,4) * t460 + Ifges(6,2) * t459 + Ifges(6,6) * t475 - pkin(5) * t577 + pkin(12) * t579 + t563 * t414 + t557 * t415 - t500 * t466 + t516 * t468;
t467 = Ifges(6,4) * t500 + Ifges(6,2) * t499 + Ifges(6,6) * t516;
t400 = mrSges(6,2) * t434 - mrSges(6,3) * t430 + Ifges(6,1) * t460 + Ifges(6,4) * t459 + Ifges(6,5) * t475 - pkin(12) * t413 - t414 * t557 + t415 * t563 + t466 * t499 - t467 * t516;
t487 = Ifges(5,4) * t518 + Ifges(5,2) * t517 + Ifges(5,6) * t547;
t488 = Ifges(5,1) * t518 + Ifges(5,4) * t517 + Ifges(5,5) * t547;
t574 = -mrSges(5,1) * t445 + mrSges(5,2) * t446 - Ifges(5,5) * t477 - Ifges(5,6) * t476 - Ifges(5,3) * t536 - pkin(4) * t573 - pkin(11) * t580 - t564 * t398 - t558 * t400 - t518 * t487 + t517 * t488;
t572 = -m(4) * t495 + t510 * mrSges(4,1) - mrSges(4,2) * t511 + t532 * t520 - t521 * t533 - t575;
t571 = mrSges(6,1) * t430 - mrSges(6,2) * t431 + Ifges(6,5) * t460 + Ifges(6,6) * t459 + Ifges(6,3) * t475 + pkin(5) * t413 + t500 * t467 - t499 * t468 - t576;
t505 = Ifges(4,4) * t533 + Ifges(4,2) * t532 + Ifges(4,6) * t548;
t506 = Ifges(4,1) * t533 + Ifges(4,4) * t532 + Ifges(4,5) * t548;
t570 = mrSges(4,1) * t464 - mrSges(4,2) * t465 + Ifges(4,5) * t511 + Ifges(4,6) * t510 + Ifges(4,3) * t537 + pkin(3) * t396 + t533 * t505 - t532 * t506 - t574;
t542 = (-mrSges(3,1) * t567 + mrSges(3,2) * t561) * t588;
t539 = -mrSges(3,2) * t552 + mrSges(3,3) * t584;
t538 = mrSges(3,1) * t552 - mrSges(3,3) * t585;
t526 = -t555 * t540 - t593;
t525 = Ifges(3,5) * t552 + (Ifges(3,1) * t561 + Ifges(3,4) * t567) * t588;
t524 = Ifges(3,6) * t552 + (Ifges(3,4) * t561 + Ifges(3,2) * t567) * t588;
t523 = Ifges(3,3) * t552 + (Ifges(3,5) * t561 + Ifges(3,6) * t567) * t588;
t513 = -g(3) * t591 + t589;
t504 = Ifges(4,5) * t533 + Ifges(4,6) * t532 + Ifges(4,3) * t548;
t486 = Ifges(5,5) * t518 + Ifges(5,6) * t517 + Ifges(5,3) * t547;
t401 = m(3) * t512 + mrSges(3,1) * t551 - mrSges(3,3) * t544 + t539 * t552 - t542 * t585 + t572;
t390 = -mrSges(5,1) * t461 + mrSges(5,3) * t446 + Ifges(5,4) * t477 + Ifges(5,2) * t476 + Ifges(5,6) * t536 - pkin(4) * t406 - t518 * t486 + t547 * t488 - t571;
t388 = m(3) * t513 - mrSges(3,2) * t551 + mrSges(3,3) * t545 - t538 * t552 + t542 * t584 + t582;
t387 = mrSges(5,2) * t461 - mrSges(5,3) * t445 + Ifges(5,1) * t477 + Ifges(5,4) * t476 + Ifges(5,5) * t536 - pkin(11) * t406 - t398 * t558 + t400 * t564 + t486 * t517 - t487 * t547;
t386 = mrSges(4,2) * t495 - mrSges(4,3) * t464 + Ifges(4,1) * t511 + Ifges(4,4) * t510 + Ifges(4,5) * t537 - pkin(10) * t396 + t387 * t565 - t390 * t559 + t504 * t532 - t505 * t548;
t385 = -mrSges(4,1) * t495 + mrSges(4,3) * t465 + Ifges(4,4) * t511 + Ifges(4,2) * t510 + Ifges(4,6) * t537 - pkin(3) * t575 + pkin(10) * t581 + t559 * t387 + t565 * t390 - t533 * t504 + t548 * t506;
t384 = Ifges(3,5) * t544 + Ifges(3,6) * t545 + Ifges(3,3) * t551 + mrSges(3,1) * t512 - mrSges(3,2) * t513 + t560 * t386 + t566 * t385 + pkin(2) * t572 + pkin(9) * t582 + (t524 * t561 - t525 * t567) * t588;
t1 = [Ifges(2,3) * qJDD(1) + mrSges(2,1) * t583 - mrSges(2,2) * t578 + (mrSges(3,2) * t526 - mrSges(3,3) * t512 + Ifges(3,1) * t544 + Ifges(3,4) * t545 + Ifges(3,5) * t551 - pkin(9) * t389 - t385 * t560 + t386 * t566 + t523 * t584 - t524 * t552) * t591 + (-mrSges(3,1) * t526 + mrSges(3,3) * t513 + Ifges(3,4) * t544 + Ifges(3,2) * t545 + Ifges(3,6) * t551 - pkin(2) * t389 - t523 * t585 + t552 * t525 - t570) * t590 + t556 * t384 + pkin(1) * ((t388 * t561 + t401 * t567) * t556 + (-m(3) * t526 + t545 * mrSges(3,1) - t544 * mrSges(3,2) + (-t538 * t561 + t539 * t567) * t588 - t389) * t555) + (t388 * t567 - t401 * t561) * t594; t384; t570; -t574; t571; -t576;];
tauJ  = t1;
