% Calculate vector of inverse dynamics joint torques for with Newton-Euler
% S6RRRPRR9
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
% Datum: 2019-05-07 13:24
% Revision: 89c353f7eff3bd693eda4e29f35b2761dbc3ada0 (2019-05-03)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauJ = S6RRRPRR9_invdynJ_fixb_snew_vp2(qJ, qJD, qJDD, g, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(6,1),zeros(3,1),zeros(13,1),zeros(7,1),zeros(7,3),zeros(7,6)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRRPRR9_invdynJ_fixb_snew_vp2: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RRRPRR9_invdynJ_fixb_snew_vp2: qJD has to be [6x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [6 1]), ...
  'S6RRRPRR9_invdynJ_fixb_snew_vp2: qJDD has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6RRRPRR9_invdynJ_fixb_snew_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [13 1]), ...
  'S6RRRPRR9_invdynJ_fixb_snew_vp2: pkin has to be [13x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6RRRPRR9_invdynJ_fixb_snew_vp2: m has to be [7x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [7,3]), ...
  'S6RRRPRR9_invdynJ_fixb_snew_vp2: mrSges has to be [7x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [7 6]), ...
  'S6RRRPRR9_invdynJ_fixb_snew_vp2: Ifges has to be [7x6] (double)'); 

%% Symbolic Calculation
% From invdyn_fixb_NewtonEuler_linkframe_tauJ_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-05-07 12:59:06
% EndTime: 2019-05-07 12:59:53
% DurationCPUTime: 35.18s
% Computational Cost: add. (537981->379), mult. (1411287->510), div. (0->0), fcn. (1198716->16), ass. (0->168)
t564 = cos(pkin(6));
t557 = qJD(1) * t564 + qJD(2);
t560 = sin(pkin(7));
t563 = cos(pkin(7));
t561 = sin(pkin(6));
t573 = cos(qJ(2));
t591 = qJD(1) * t573;
t587 = t561 * t591;
t544 = (t557 * t560 + t563 * t587) * pkin(10);
t568 = sin(qJ(2));
t593 = qJD(1) * t561;
t605 = pkin(10) * t560;
t548 = (-pkin(2) * t573 - t568 * t605) * t593;
t589 = qJD(1) * qJD(2);
t554 = (qJDD(1) * t568 + t573 * t589) * t561;
t556 = qJDD(1) * t564 + qJDD(2);
t575 = qJD(1) ^ 2;
t569 = sin(qJ(1));
t574 = cos(qJ(1));
t581 = -g(1) * t574 - g(2) * t569;
t606 = pkin(9) * t561;
t552 = -pkin(1) * t575 + qJDD(1) * t606 + t581;
t586 = t569 * g(1) - g(2) * t574;
t551 = qJDD(1) * pkin(1) + t575 * t606 + t586;
t602 = t551 * t564;
t582 = -t568 * t552 + t573 * t602;
t592 = qJD(1) * t568;
t604 = pkin(10) * t563;
t504 = -t554 * t604 + t556 * pkin(2) + t557 * t544 + (-g(3) * t573 - t548 * t592) * t561 + t582;
t588 = t561 * t592;
t547 = pkin(2) * t557 - t588 * t604;
t555 = (qJDD(1) * t573 - t568 * t589) * t561;
t580 = t555 * t563 + t556 * t560;
t594 = t573 * t552 + t568 * t602;
t505 = -t557 * t547 + (-g(3) * t568 + t548 * t591) * t561 + t580 * pkin(10) + t594;
t603 = t564 * g(3);
t509 = -t554 * t605 - t555 * pkin(2) - t603 + (-t551 + (-t544 * t573 + t547 * t568) * qJD(1)) * t561;
t567 = sin(qJ(3));
t572 = cos(qJ(3));
t596 = t563 * t572;
t600 = t560 * t572;
t470 = t504 * t596 - t567 * t505 + t509 * t600;
t595 = t563 * t573;
t534 = t557 * t600 + (-t567 * t568 + t572 * t595) * t593;
t520 = t534 * qJD(3) + t572 * t554 + t567 * t580;
t601 = t560 * t567;
t535 = t557 * t601 + (t567 * t595 + t568 * t572) * t593;
t536 = -t555 * t560 + t556 * t563 + qJDD(3);
t545 = t557 * t563 - t560 * t587 + qJD(3);
t460 = (t534 * t545 - t520) * qJ(4) + (t534 * t535 + t536) * pkin(3) + t470;
t597 = t563 * t567;
t471 = t504 * t597 + t572 * t505 + t509 * t601;
t519 = -t535 * qJD(3) - t567 * t554 + t572 * t580;
t529 = pkin(3) * t545 - qJ(4) * t535;
t533 = t534 ^ 2;
t463 = -pkin(3) * t533 + qJ(4) * t519 - t529 * t545 + t471;
t559 = sin(pkin(13));
t562 = cos(pkin(13));
t526 = t534 * t559 + t535 * t562;
t452 = -0.2e1 * qJD(4) * t526 + t562 * t460 - t559 * t463;
t599 = t561 * t568;
t598 = t561 * t573;
t525 = t534 * t562 - t535 * t559;
t453 = 0.2e1 * qJD(4) * t525 + t559 * t460 + t562 * t463;
t491 = t519 * t562 - t520 * t559;
t499 = -mrSges(5,1) * t525 + mrSges(5,2) * t526;
t514 = mrSges(5,1) * t545 - mrSges(5,3) * t526;
t500 = -pkin(4) * t525 - pkin(11) * t526;
t543 = t545 ^ 2;
t451 = -pkin(4) * t543 + pkin(11) * t536 + t500 * t525 + t453;
t483 = -t560 * t504 + t563 * t509;
t469 = -t519 * pkin(3) - t533 * qJ(4) + t535 * t529 + qJDD(4) + t483;
t492 = t519 * t559 + t520 * t562;
t455 = (-t525 * t545 - t492) * pkin(11) + (t526 * t545 - t491) * pkin(4) + t469;
t566 = sin(qJ(5));
t571 = cos(qJ(5));
t447 = t571 * t451 + t566 * t455;
t511 = -t526 * t566 + t545 * t571;
t512 = t526 * t571 + t545 * t566;
t486 = -pkin(5) * t511 - pkin(12) * t512;
t488 = qJDD(5) - t491;
t524 = qJD(5) - t525;
t523 = t524 ^ 2;
t445 = -pkin(5) * t523 + pkin(12) * t488 + t486 * t511 + t447;
t450 = -t536 * pkin(4) - t543 * pkin(11) + t526 * t500 - t452;
t475 = -qJD(5) * t512 - t492 * t566 + t536 * t571;
t476 = qJD(5) * t511 + t492 * t571 + t536 * t566;
t448 = (-t511 * t524 - t476) * pkin(12) + (t512 * t524 - t475) * pkin(5) + t450;
t565 = sin(qJ(6));
t570 = cos(qJ(6));
t441 = -t445 * t565 + t448 * t570;
t489 = -t512 * t565 + t524 * t570;
t458 = qJD(6) * t489 + t476 * t570 + t488 * t565;
t490 = t512 * t570 + t524 * t565;
t472 = -mrSges(7,1) * t489 + mrSges(7,2) * t490;
t474 = qJDD(6) - t475;
t510 = qJD(6) - t511;
t477 = -mrSges(7,2) * t510 + mrSges(7,3) * t489;
t439 = m(7) * t441 + mrSges(7,1) * t474 - mrSges(7,3) * t458 - t472 * t490 + t477 * t510;
t442 = t445 * t570 + t448 * t565;
t457 = -qJD(6) * t490 - t476 * t565 + t488 * t570;
t478 = mrSges(7,1) * t510 - mrSges(7,3) * t490;
t440 = m(7) * t442 - mrSges(7,2) * t474 + mrSges(7,3) * t457 + t472 * t489 - t478 * t510;
t433 = -t439 * t565 + t570 * t440;
t485 = -mrSges(6,1) * t511 + mrSges(6,2) * t512;
t494 = mrSges(6,1) * t524 - mrSges(6,3) * t512;
t431 = m(6) * t447 - mrSges(6,2) * t488 + mrSges(6,3) * t475 + t485 * t511 - t494 * t524 + t433;
t446 = -t451 * t566 + t455 * t571;
t444 = -pkin(5) * t488 - pkin(12) * t523 + t486 * t512 - t446;
t443 = -m(7) * t444 + t457 * mrSges(7,1) - mrSges(7,2) * t458 + t489 * t477 - t478 * t490;
t493 = -mrSges(6,2) * t524 + mrSges(6,3) * t511;
t437 = m(6) * t446 + mrSges(6,1) * t488 - mrSges(6,3) * t476 - t485 * t512 + t493 * t524 + t443;
t584 = t571 * t431 - t437 * t566;
t422 = m(5) * t453 - mrSges(5,2) * t536 + mrSges(5,3) * t491 + t499 * t525 - t514 * t545 + t584;
t513 = -mrSges(5,2) * t545 + mrSges(5,3) * t525;
t432 = t439 * t570 + t440 * t565;
t578 = -m(6) * t450 + t475 * mrSges(6,1) - mrSges(6,2) * t476 + t511 * t493 - t494 * t512 - t432;
t428 = m(5) * t452 + mrSges(5,1) * t536 - mrSges(5,3) * t492 - t499 * t526 + t513 * t545 + t578;
t417 = t559 * t422 + t562 * t428;
t426 = t566 * t431 + t571 * t437;
t527 = -mrSges(4,1) * t534 + mrSges(4,2) * t535;
t528 = -mrSges(4,2) * t545 + mrSges(4,3) * t534;
t415 = m(4) * t470 + mrSges(4,1) * t536 - mrSges(4,3) * t520 - t527 * t535 + t528 * t545 + t417;
t530 = mrSges(4,1) * t545 - mrSges(4,3) * t535;
t585 = t562 * t422 - t428 * t559;
t416 = m(4) * t471 - mrSges(4,2) * t536 + mrSges(4,3) * t519 + t527 * t534 - t530 * t545 + t585;
t425 = m(5) * t469 - mrSges(5,1) * t491 + t492 * mrSges(5,2) - t513 * t525 + t526 * t514 + t426;
t424 = m(4) * t483 - mrSges(4,1) * t519 + mrSges(4,2) * t520 - t528 * t534 + t530 * t535 + t425;
t404 = t415 * t600 + t416 * t601 + t563 * t424;
t408 = -t415 * t567 + t572 * t416;
t405 = t415 * t596 + t416 * t597 - t424 * t560;
t464 = Ifges(7,5) * t490 + Ifges(7,6) * t489 + Ifges(7,3) * t510;
t466 = Ifges(7,1) * t490 + Ifges(7,4) * t489 + Ifges(7,5) * t510;
t434 = -mrSges(7,1) * t444 + mrSges(7,3) * t442 + Ifges(7,4) * t458 + Ifges(7,2) * t457 + Ifges(7,6) * t474 - t464 * t490 + t466 * t510;
t465 = Ifges(7,4) * t490 + Ifges(7,2) * t489 + Ifges(7,6) * t510;
t435 = mrSges(7,2) * t444 - mrSges(7,3) * t441 + Ifges(7,1) * t458 + Ifges(7,4) * t457 + Ifges(7,5) * t474 + t464 * t489 - t465 * t510;
t479 = Ifges(6,5) * t512 + Ifges(6,6) * t511 + Ifges(6,3) * t524;
t480 = Ifges(6,4) * t512 + Ifges(6,2) * t511 + Ifges(6,6) * t524;
t418 = mrSges(6,2) * t450 - mrSges(6,3) * t446 + Ifges(6,1) * t476 + Ifges(6,4) * t475 + Ifges(6,5) * t488 - pkin(12) * t432 - t434 * t565 + t435 * t570 + t479 * t511 - t480 * t524;
t481 = Ifges(6,1) * t512 + Ifges(6,4) * t511 + Ifges(6,5) * t524;
t577 = mrSges(7,1) * t441 - mrSges(7,2) * t442 + Ifges(7,5) * t458 + Ifges(7,6) * t457 + Ifges(7,3) * t474 + t465 * t490 - t466 * t489;
t419 = -mrSges(6,1) * t450 + mrSges(6,3) * t447 + Ifges(6,4) * t476 + Ifges(6,2) * t475 + Ifges(6,6) * t488 - pkin(5) * t432 - t479 * t512 + t481 * t524 - t577;
t495 = Ifges(5,5) * t526 + Ifges(5,6) * t525 + Ifges(5,3) * t545;
t496 = Ifges(5,4) * t526 + Ifges(5,2) * t525 + Ifges(5,6) * t545;
t406 = mrSges(5,2) * t469 - mrSges(5,3) * t452 + Ifges(5,1) * t492 + Ifges(5,4) * t491 + Ifges(5,5) * t536 - pkin(11) * t426 + t418 * t571 - t419 * t566 + t495 * t525 - t496 * t545;
t497 = Ifges(5,1) * t526 + Ifges(5,4) * t525 + Ifges(5,5) * t545;
t576 = mrSges(6,1) * t446 - mrSges(6,2) * t447 + Ifges(6,5) * t476 + Ifges(6,6) * t475 + Ifges(6,3) * t488 + pkin(5) * t443 + pkin(12) * t433 + t570 * t434 + t565 * t435 + t512 * t480 - t511 * t481;
t409 = -mrSges(5,1) * t469 + mrSges(5,3) * t453 + Ifges(5,4) * t492 + Ifges(5,2) * t491 + Ifges(5,6) * t536 - pkin(4) * t426 - t526 * t495 + t545 * t497 - t576;
t515 = Ifges(4,5) * t535 + Ifges(4,6) * t534 + Ifges(4,3) * t545;
t517 = Ifges(4,1) * t535 + Ifges(4,4) * t534 + Ifges(4,5) * t545;
t400 = -mrSges(4,1) * t483 + mrSges(4,3) * t471 + Ifges(4,4) * t520 + Ifges(4,2) * t519 + Ifges(4,6) * t536 - pkin(3) * t425 + qJ(4) * t585 + t559 * t406 + t562 * t409 - t535 * t515 + t545 * t517;
t516 = Ifges(4,4) * t535 + Ifges(4,2) * t534 + Ifges(4,6) * t545;
t401 = mrSges(4,2) * t483 - mrSges(4,3) * t470 + Ifges(4,1) * t520 + Ifges(4,4) * t519 + Ifges(4,5) * t536 - qJ(4) * t417 + t406 * t562 - t409 * t559 + t515 * t534 - t516 * t545;
t579 = pkin(10) * t408 + t400 * t572 + t401 * t567;
t553 = (-mrSges(3,1) * t573 + mrSges(3,2) * t568) * t593;
t550 = -mrSges(3,2) * t557 + mrSges(3,3) * t587;
t549 = mrSges(3,1) * t557 - mrSges(3,3) * t588;
t540 = -t561 * t551 - t603;
t539 = Ifges(3,5) * t557 + (Ifges(3,1) * t568 + Ifges(3,4) * t573) * t593;
t538 = Ifges(3,6) * t557 + (Ifges(3,4) * t568 + Ifges(3,2) * t573) * t593;
t537 = Ifges(3,3) * t557 + (Ifges(3,5) * t568 + Ifges(3,6) * t573) * t593;
t532 = -g(3) * t599 + t594;
t531 = -g(3) * t598 + t582;
t407 = m(3) * t532 - mrSges(3,2) * t556 + mrSges(3,3) * t555 - t549 * t557 + t553 * t587 + t408;
t403 = m(3) * t531 + mrSges(3,1) * t556 - mrSges(3,3) * t554 + t550 * t557 - t553 * t588 + t405;
t402 = Ifges(4,5) * t520 + Ifges(4,6) * t519 + t535 * t516 - t534 * t517 + mrSges(4,1) * t470 - mrSges(4,2) * t471 + Ifges(5,5) * t492 + Ifges(5,6) * t491 + t526 * t496 - t525 * t497 + mrSges(5,1) * t452 - mrSges(5,2) * t453 + t566 * t418 + t571 * t419 + pkin(4) * t578 + pkin(11) * t584 + pkin(3) * t417 + (Ifges(4,3) + Ifges(5,3)) * t536;
t399 = mrSges(3,1) * t531 - mrSges(3,2) * t532 + Ifges(3,5) * t554 + Ifges(3,6) * t555 + Ifges(3,3) * t556 + pkin(2) * t405 + t563 * t402 + (t538 * t568 - t539 * t573) * t593 + t579 * t560;
t1 = [Ifges(2,3) * qJDD(1) + mrSges(2,1) * t586 - mrSges(2,2) * t581 + (t537 * t587 + mrSges(3,2) * t540 - mrSges(3,3) * t531 + Ifges(3,1) * t554 + Ifges(3,4) * t555 + Ifges(3,5) * t556 - t567 * t400 + t572 * t401 - t557 * t538 + (-t404 * t560 - t405 * t563) * pkin(10)) * t599 + (-mrSges(3,1) * t540 + mrSges(3,3) * t532 + Ifges(3,4) * t554 + Ifges(3,2) * t555 + Ifges(3,6) * t556 - pkin(2) * t404 - t560 * t402 - t537 * t588 + t557 * t539 + t563 * t579) * t598 + t564 * t399 + pkin(1) * ((t403 * t573 + t407 * t568) * t564 + (-m(3) * t540 + t555 * mrSges(3,1) - t554 * mrSges(3,2) + (-t549 * t568 + t550 * t573) * t593 - t404) * t561) + (-t403 * t568 + t407 * t573) * t606; t399; t402; t425; t576; t577;];
tauJ  = t1;
