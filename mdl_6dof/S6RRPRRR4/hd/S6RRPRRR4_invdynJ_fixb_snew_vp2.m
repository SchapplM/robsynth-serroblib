% Calculate vector of inverse dynamics joint torques for with Newton-Euler
% S6RRPRRR4
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
%   pkin=[a2,a3,a4,a5,a6,alpha2,d1,d2,d4,d5,d6,theta3]';
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
% Datum: 2019-05-06 20:52
% Revision: 89c353f7eff3bd693eda4e29f35b2761dbc3ada0 (2019-05-03)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauJ = S6RRPRRR4_invdynJ_fixb_snew_vp2(qJ, qJD, qJDD, g, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(6,1),zeros(3,1),zeros(12,1),zeros(7,1),zeros(7,3),zeros(7,6)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRPRRR4_invdynJ_fixb_snew_vp2: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RRPRRR4_invdynJ_fixb_snew_vp2: qJD has to be [6x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [6 1]), ...
  'S6RRPRRR4_invdynJ_fixb_snew_vp2: qJDD has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6RRPRRR4_invdynJ_fixb_snew_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [12 1]), ...
  'S6RRPRRR4_invdynJ_fixb_snew_vp2: pkin has to be [12x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6RRPRRR4_invdynJ_fixb_snew_vp2: m has to be [7x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [7,3]), ...
  'S6RRPRRR4_invdynJ_fixb_snew_vp2: mrSges has to be [7x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [7 6]), ...
  'S6RRPRRR4_invdynJ_fixb_snew_vp2: Ifges has to be [7x6] (double)'); 

%% Symbolic Calculation
% From invdyn_fixb_NewtonEuler_linkframe_tauJ_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-05-06 20:41:35
% EndTime: 2019-05-06 20:41:52
% DurationCPUTime: 12.18s
% Computational Cost: add. (167462->362), mult. (437025->471), div. (0->0), fcn. (351747->14), ass. (0->151)
t557 = sin(pkin(6));
t563 = sin(qJ(2));
t568 = cos(qJ(2));
t589 = qJD(1) * qJD(2);
t547 = (qJDD(1) * t563 + t568 * t589) * t557;
t559 = cos(pkin(6));
t551 = qJDD(1) * t559 + qJDD(2);
t552 = qJD(1) * t559 + qJD(2);
t570 = qJD(1) ^ 2;
t564 = sin(qJ(1));
t569 = cos(qJ(1));
t578 = -g(1) * t569 - g(2) * t564;
t596 = pkin(8) * t557;
t545 = -pkin(1) * t570 + qJDD(1) * t596 + t578;
t585 = t564 * g(1) - g(2) * t569;
t544 = qJDD(1) * pkin(1) + t570 * t596 + t585;
t595 = t544 * t559;
t579 = -t563 * t545 + t568 * t595;
t594 = t557 ^ 2 * t570;
t484 = t551 * pkin(2) - t547 * qJ(3) + (pkin(2) * t563 * t594 + (qJ(3) * qJD(1) * t552 - g(3)) * t557) * t568 + t579;
t593 = t557 * t563;
t511 = -g(3) * t593 + t568 * t545 + t563 * t595;
t591 = qJD(1) * t557;
t587 = t563 * t591;
t541 = pkin(2) * t552 - qJ(3) * t587;
t548 = (qJDD(1) * t568 - t563 * t589) * t557;
t588 = t568 ^ 2 * t594;
t487 = -pkin(2) * t588 + qJ(3) * t548 - t541 * t552 + t511;
t556 = sin(pkin(12));
t558 = cos(pkin(12));
t538 = (t556 * t568 + t558 * t563) * t591;
t465 = -0.2e1 * qJD(3) * t538 + t558 * t484 - t556 * t487;
t592 = t557 * t568;
t586 = t568 * t591;
t537 = -t556 * t587 + t558 * t586;
t466 = 0.2e1 * qJD(3) * t537 + t556 * t484 + t558 * t487;
t512 = -mrSges(4,1) * t537 + mrSges(4,2) * t538;
t518 = -t556 * t547 + t548 * t558;
t524 = mrSges(4,1) * t552 - mrSges(4,3) * t538;
t514 = -pkin(3) * t537 - pkin(9) * t538;
t550 = t552 ^ 2;
t453 = -pkin(3) * t550 + pkin(9) * t551 + t514 * t537 + t466;
t528 = -t559 * g(3) - t557 * t544;
t497 = -t548 * pkin(2) - qJ(3) * t588 + t541 * t587 + qJDD(3) + t528;
t519 = t547 * t558 + t548 * t556;
t469 = (-t537 * t552 - t519) * pkin(9) + (t538 * t552 - t518) * pkin(3) + t497;
t562 = sin(qJ(4));
t567 = cos(qJ(4));
t445 = -t562 * t453 + t567 * t469;
t521 = -t538 * t562 + t552 * t567;
t495 = qJD(4) * t521 + t519 * t567 + t551 * t562;
t517 = qJDD(4) - t518;
t522 = t538 * t567 + t552 * t562;
t536 = qJD(4) - t537;
t442 = (t521 * t536 - t495) * pkin(10) + (t521 * t522 + t517) * pkin(4) + t445;
t446 = t567 * t453 + t562 * t469;
t494 = -qJD(4) * t522 - t519 * t562 + t551 * t567;
t505 = pkin(4) * t536 - pkin(10) * t522;
t520 = t521 ^ 2;
t444 = -pkin(4) * t520 + pkin(10) * t494 - t505 * t536 + t446;
t561 = sin(qJ(5));
t566 = cos(qJ(5));
t439 = t561 * t442 + t566 * t444;
t500 = t521 * t561 + t522 * t566;
t462 = -qJD(5) * t500 + t494 * t566 - t495 * t561;
t499 = t521 * t566 - t522 * t561;
t477 = -mrSges(6,1) * t499 + mrSges(6,2) * t500;
t531 = qJD(5) + t536;
t486 = mrSges(6,1) * t531 - mrSges(6,3) * t500;
t515 = qJDD(5) + t517;
t478 = -pkin(5) * t499 - pkin(11) * t500;
t530 = t531 ^ 2;
t436 = -pkin(5) * t530 + pkin(11) * t515 + t478 * t499 + t439;
t452 = -t551 * pkin(3) - t550 * pkin(9) + t538 * t514 - t465;
t447 = -t494 * pkin(4) - t520 * pkin(10) + t522 * t505 + t452;
t463 = qJD(5) * t499 + t494 * t561 + t495 * t566;
t440 = (-t499 * t531 - t463) * pkin(11) + (t500 * t531 - t462) * pkin(5) + t447;
t560 = sin(qJ(6));
t565 = cos(qJ(6));
t433 = -t436 * t560 + t440 * t565;
t480 = -t500 * t560 + t531 * t565;
t450 = qJD(6) * t480 + t463 * t565 + t515 * t560;
t461 = qJDD(6) - t462;
t481 = t500 * t565 + t531 * t560;
t470 = -mrSges(7,1) * t480 + mrSges(7,2) * t481;
t498 = qJD(6) - t499;
t471 = -mrSges(7,2) * t498 + mrSges(7,3) * t480;
t429 = m(7) * t433 + mrSges(7,1) * t461 - mrSges(7,3) * t450 - t470 * t481 + t471 * t498;
t434 = t436 * t565 + t440 * t560;
t449 = -qJD(6) * t481 - t463 * t560 + t515 * t565;
t472 = mrSges(7,1) * t498 - mrSges(7,3) * t481;
t430 = m(7) * t434 - mrSges(7,2) * t461 + mrSges(7,3) * t449 + t470 * t480 - t472 * t498;
t581 = -t429 * t560 + t565 * t430;
t416 = m(6) * t439 - mrSges(6,2) * t515 + mrSges(6,3) * t462 + t477 * t499 - t486 * t531 + t581;
t438 = t442 * t566 - t444 * t561;
t485 = -mrSges(6,2) * t531 + mrSges(6,3) * t499;
t435 = -pkin(5) * t515 - pkin(11) * t530 + t478 * t500 - t438;
t577 = -m(7) * t435 + t449 * mrSges(7,1) - mrSges(7,2) * t450 + t480 * t471 - t472 * t481;
t425 = m(6) * t438 + mrSges(6,1) * t515 - mrSges(6,3) * t463 - t477 * t500 + t485 * t531 + t577;
t411 = t561 * t416 + t566 * t425;
t501 = -mrSges(5,1) * t521 + mrSges(5,2) * t522;
t503 = -mrSges(5,2) * t536 + mrSges(5,3) * t521;
t409 = m(5) * t445 + mrSges(5,1) * t517 - mrSges(5,3) * t495 - t501 * t522 + t503 * t536 + t411;
t504 = mrSges(5,1) * t536 - mrSges(5,3) * t522;
t582 = t566 * t416 - t425 * t561;
t410 = m(5) * t446 - mrSges(5,2) * t517 + mrSges(5,3) * t494 + t501 * t521 - t504 * t536 + t582;
t583 = -t409 * t562 + t567 * t410;
t401 = m(4) * t466 - mrSges(4,2) * t551 + mrSges(4,3) * t518 + t512 * t537 - t524 * t552 + t583;
t523 = -mrSges(4,2) * t552 + mrSges(4,3) * t537;
t418 = t565 * t429 + t560 * t430;
t575 = m(6) * t447 - t462 * mrSges(6,1) + mrSges(6,2) * t463 - t499 * t485 + t486 * t500 + t418;
t572 = -m(5) * t452 + t494 * mrSges(5,1) - mrSges(5,2) * t495 + t521 * t503 - t504 * t522 - t575;
t413 = m(4) * t465 + mrSges(4,1) * t551 - mrSges(4,3) * t519 - t512 * t538 + t523 * t552 + t572;
t398 = t556 * t401 + t558 * t413;
t403 = t567 * t409 + t562 * t410;
t584 = t558 * t401 - t413 * t556;
t576 = -m(4) * t497 + t518 * mrSges(4,1) - t519 * mrSges(4,2) + t537 * t523 - t538 * t524 - t403;
t454 = Ifges(7,5) * t481 + Ifges(7,6) * t480 + Ifges(7,3) * t498;
t456 = Ifges(7,1) * t481 + Ifges(7,4) * t480 + Ifges(7,5) * t498;
t422 = -mrSges(7,1) * t435 + mrSges(7,3) * t434 + Ifges(7,4) * t450 + Ifges(7,2) * t449 + Ifges(7,6) * t461 - t454 * t481 + t456 * t498;
t455 = Ifges(7,4) * t481 + Ifges(7,2) * t480 + Ifges(7,6) * t498;
t423 = mrSges(7,2) * t435 - mrSges(7,3) * t433 + Ifges(7,1) * t450 + Ifges(7,4) * t449 + Ifges(7,5) * t461 + t454 * t480 - t455 * t498;
t474 = Ifges(6,4) * t500 + Ifges(6,2) * t499 + Ifges(6,6) * t531;
t475 = Ifges(6,1) * t500 + Ifges(6,4) * t499 + Ifges(6,5) * t531;
t574 = -mrSges(6,1) * t438 + mrSges(6,2) * t439 - Ifges(6,5) * t463 - Ifges(6,6) * t462 - Ifges(6,3) * t515 - pkin(5) * t577 - pkin(11) * t581 - t565 * t422 - t560 * t423 - t500 * t474 + t499 * t475;
t573 = mrSges(7,1) * t433 - mrSges(7,2) * t434 + Ifges(7,5) * t450 + Ifges(7,6) * t449 + Ifges(7,3) * t461 + t455 * t481 - t456 * t480;
t489 = Ifges(5,4) * t522 + Ifges(5,2) * t521 + Ifges(5,6) * t536;
t490 = Ifges(5,1) * t522 + Ifges(5,4) * t521 + Ifges(5,5) * t536;
t571 = mrSges(5,1) * t445 - mrSges(5,2) * t446 + Ifges(5,5) * t495 + Ifges(5,6) * t494 + Ifges(5,3) * t517 + pkin(4) * t411 + t522 * t489 - t521 * t490 - t574;
t546 = (-mrSges(3,1) * t568 + mrSges(3,2) * t563) * t591;
t543 = -mrSges(3,2) * t552 + mrSges(3,3) * t586;
t542 = mrSges(3,1) * t552 - mrSges(3,3) * t587;
t527 = Ifges(3,5) * t552 + (Ifges(3,1) * t563 + Ifges(3,4) * t568) * t591;
t526 = Ifges(3,6) * t552 + (Ifges(3,4) * t563 + Ifges(3,2) * t568) * t591;
t525 = Ifges(3,3) * t552 + (Ifges(3,5) * t563 + Ifges(3,6) * t568) * t591;
t510 = -g(3) * t592 + t579;
t508 = Ifges(4,1) * t538 + Ifges(4,4) * t537 + Ifges(4,5) * t552;
t507 = Ifges(4,4) * t538 + Ifges(4,2) * t537 + Ifges(4,6) * t552;
t506 = Ifges(4,5) * t538 + Ifges(4,6) * t537 + Ifges(4,3) * t552;
t488 = Ifges(5,5) * t522 + Ifges(5,6) * t521 + Ifges(5,3) * t536;
t473 = Ifges(6,5) * t500 + Ifges(6,6) * t499 + Ifges(6,3) * t531;
t405 = -mrSges(6,1) * t447 + mrSges(6,3) * t439 + Ifges(6,4) * t463 + Ifges(6,2) * t462 + Ifges(6,6) * t515 - pkin(5) * t418 - t473 * t500 + t475 * t531 - t573;
t404 = mrSges(6,2) * t447 - mrSges(6,3) * t438 + Ifges(6,1) * t463 + Ifges(6,4) * t462 + Ifges(6,5) * t515 - pkin(11) * t418 - t422 * t560 + t423 * t565 + t473 * t499 - t474 * t531;
t397 = m(3) * t511 - mrSges(3,2) * t551 + mrSges(3,3) * t548 - t542 * t552 + t546 * t586 + t584;
t396 = m(3) * t510 + mrSges(3,1) * t551 - mrSges(3,3) * t547 + t543 * t552 - t546 * t587 + t398;
t395 = mrSges(5,2) * t452 - mrSges(5,3) * t445 + Ifges(5,1) * t495 + Ifges(5,4) * t494 + Ifges(5,5) * t517 - pkin(10) * t411 + t404 * t566 - t405 * t561 + t488 * t521 - t489 * t536;
t394 = -mrSges(5,1) * t452 + mrSges(5,3) * t446 + Ifges(5,4) * t495 + Ifges(5,2) * t494 + Ifges(5,6) * t517 - pkin(4) * t575 + pkin(10) * t582 + t561 * t404 + t566 * t405 - t522 * t488 + t536 * t490;
t393 = -mrSges(4,1) * t497 + mrSges(4,3) * t466 + Ifges(4,4) * t519 + Ifges(4,2) * t518 + Ifges(4,6) * t551 - pkin(3) * t403 - t538 * t506 + t552 * t508 - t571;
t392 = mrSges(4,2) * t497 - mrSges(4,3) * t465 + Ifges(4,1) * t519 + Ifges(4,4) * t518 + Ifges(4,5) * t551 - pkin(9) * t403 - t394 * t562 + t395 * t567 + t506 * t537 - t507 * t552;
t391 = Ifges(3,5) * t547 + Ifges(3,6) * t548 + mrSges(3,1) * t510 - mrSges(3,2) * t511 + Ifges(4,5) * t519 + Ifges(4,6) * t518 + t538 * t507 - t537 * t508 + mrSges(4,1) * t465 - mrSges(4,2) * t466 + t562 * t395 + t567 * t394 + pkin(3) * t572 + pkin(9) * t583 + pkin(2) * t398 + (Ifges(3,3) + Ifges(4,3)) * t551 + (t526 * t563 - t527 * t568) * t591;
t1 = [Ifges(2,3) * qJDD(1) + mrSges(2,1) * t585 - mrSges(2,2) * t578 + (mrSges(3,2) * t528 - mrSges(3,3) * t510 + Ifges(3,1) * t547 + Ifges(3,4) * t548 + Ifges(3,5) * t551 - qJ(3) * t398 + t392 * t558 - t393 * t556 + t525 * t586 - t526 * t552) * t593 + (-mrSges(3,1) * t528 + mrSges(3,3) * t511 + Ifges(3,4) * t547 + Ifges(3,2) * t548 + Ifges(3,6) * t551 + pkin(2) * t576 + qJ(3) * t584 + t556 * t392 + t558 * t393 - t525 * t587 + t552 * t527) * t592 + t559 * t391 + pkin(1) * ((t396 * t568 + t397 * t563) * t559 + (-m(3) * t528 + t548 * mrSges(3,1) - t547 * mrSges(3,2) + (-t542 * t563 + t543 * t568) * t591 + t576) * t557) + (-t396 * t563 + t397 * t568) * t596; t391; -t576; t571; -t574; t573;];
tauJ  = t1;
