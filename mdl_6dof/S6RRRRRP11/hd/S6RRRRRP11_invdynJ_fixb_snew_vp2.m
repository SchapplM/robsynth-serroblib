% Calculate vector of inverse dynamics joint torques for with Newton-Euler
% S6RRRRRP11
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
%   pkin=[a2,a3,a4,a5,a6,alpha2,alpha3,d1,d2,d3,d4,d5]';
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
% Datum: 2019-05-08 07:15
% Revision: 89c353f7eff3bd693eda4e29f35b2761dbc3ada0 (2019-05-03)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauJ = S6RRRRRP11_invdynJ_fixb_snew_vp2(qJ, qJD, qJDD, g, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(6,1),zeros(3,1),zeros(12,1),zeros(7,1),zeros(7,3),zeros(7,6)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRRRRP11_invdynJ_fixb_snew_vp2: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RRRRRP11_invdynJ_fixb_snew_vp2: qJD has to be [6x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [6 1]), ...
  'S6RRRRRP11_invdynJ_fixb_snew_vp2: qJDD has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6RRRRRP11_invdynJ_fixb_snew_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [12 1]), ...
  'S6RRRRRP11_invdynJ_fixb_snew_vp2: pkin has to be [12x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6RRRRRP11_invdynJ_fixb_snew_vp2: m has to be [7x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [7,3]), ...
  'S6RRRRRP11_invdynJ_fixb_snew_vp2: mrSges has to be [7x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [7 6]), ...
  'S6RRRRRP11_invdynJ_fixb_snew_vp2: Ifges has to be [7x6] (double)'); 

%% Symbolic Calculation
% From invdyn_fixb_NewtonEuler_linkframe_tauJ_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-05-08 06:45:25
% EndTime: 2019-05-08 06:45:50
% DurationCPUTime: 18.39s
% Computational Cost: add. (292844->358), mult. (722209->468), div. (0->0), fcn. (608051->14), ass. (0->161)
t611 = Ifges(6,1) + Ifges(7,1);
t602 = Ifges(6,4) + Ifges(7,4);
t601 = Ifges(6,5) + Ifges(7,5);
t610 = Ifges(6,2) + Ifges(7,2);
t600 = Ifges(6,6) + Ifges(7,6);
t609 = Ifges(6,3) + Ifges(7,3);
t553 = cos(pkin(6));
t548 = qJD(1) * t553 + qJD(2);
t550 = sin(pkin(7));
t552 = cos(pkin(7));
t551 = sin(pkin(6));
t562 = cos(qJ(2));
t584 = qJD(1) * t562;
t578 = t551 * t584;
t535 = (t548 * t550 + t552 * t578) * pkin(10);
t557 = sin(qJ(2));
t586 = qJD(1) * t551;
t606 = pkin(10) * t550;
t539 = (-pkin(2) * t562 - t557 * t606) * t586;
t583 = qJD(1) * qJD(2);
t545 = (qJDD(1) * t557 + t562 * t583) * t551;
t547 = qJDD(1) * t553 + qJDD(2);
t564 = qJD(1) ^ 2;
t558 = sin(qJ(1));
t563 = cos(qJ(1));
t573 = -g(1) * t563 - g(2) * t558;
t607 = pkin(9) * t551;
t543 = -pkin(1) * t564 + qJDD(1) * t607 + t573;
t577 = t558 * g(1) - g(2) * t563;
t542 = qJDD(1) * pkin(1) + t564 * t607 + t577;
t598 = t542 * t553;
t575 = -t543 * t557 + t562 * t598;
t585 = qJD(1) * t557;
t605 = pkin(10) * t552;
t495 = -t545 * t605 + pkin(2) * t547 + t535 * t548 + (-g(3) * t562 - t539 * t585) * t551 + t575;
t579 = t551 * t585;
t538 = pkin(2) * t548 - t579 * t605;
t546 = (qJDD(1) * t562 - t557 * t583) * t551;
t571 = t546 * t552 + t547 * t550;
t587 = t562 * t543 + t557 * t598;
t496 = -t538 * t548 + (-g(3) * t557 + t539 * t584) * t551 + t571 * pkin(10) + t587;
t604 = g(3) * t553;
t501 = -t545 * t606 - pkin(2) * t546 - t604 + (-t542 + (-t535 * t562 + t538 * t557) * qJD(1)) * t551;
t556 = sin(qJ(3));
t561 = cos(qJ(3));
t462 = -t556 * t496 + (t495 * t552 + t501 * t550) * t561;
t592 = t552 * t562;
t597 = t550 * t556;
t526 = t548 * t597 + (t556 * t592 + t557 * t561) * t586;
t512 = -qJD(3) * t526 - t545 * t556 + t571 * t561;
t596 = t550 * t561;
t525 = (-t556 * t557 + t561 * t592) * t586 + t548 * t596;
t513 = qJD(3) * t525 + t545 * t561 + t571 * t556;
t536 = t548 * t552 - t550 * t578 + qJD(3);
t555 = sin(qJ(4));
t560 = cos(qJ(4));
t517 = -t526 * t555 + t536 * t560;
t527 = -t546 * t550 + t547 * t552 + qJDD(3);
t481 = qJD(4) * t517 + t513 * t560 + t527 * t555;
t518 = t526 * t560 + t536 * t555;
t524 = qJD(4) - t525;
t554 = sin(qJ(5));
t559 = cos(qJ(5));
t504 = -t518 * t554 + t524 * t559;
t511 = qJDD(4) - t512;
t461 = qJD(5) * t504 + t481 * t559 + t511 * t554;
t505 = t518 * t559 + t524 * t554;
t475 = -mrSges(7,1) * t504 + mrSges(7,2) * t505;
t593 = t552 * t556;
t463 = t495 * t593 + t561 * t496 + t501 * t597;
t515 = -pkin(3) * t525 - pkin(11) * t526;
t534 = t536 ^ 2;
t454 = -pkin(3) * t534 + pkin(11) * t527 + t515 * t525 + t463;
t473 = -t495 * t550 + t552 * t501;
t456 = (-t525 * t536 - t513) * pkin(11) + (t526 * t536 - t512) * pkin(3) + t473;
t450 = t560 * t454 + t555 * t456;
t498 = -pkin(4) * t517 - pkin(12) * t518;
t523 = t524 ^ 2;
t445 = -pkin(4) * t523 + pkin(12) * t511 + t498 * t517 + t450;
t453 = -pkin(3) * t527 - pkin(11) * t534 + t526 * t515 - t462;
t480 = -qJD(4) * t518 - t513 * t555 + t527 * t560;
t448 = (-t517 * t524 - t481) * pkin(12) + (t518 * t524 - t480) * pkin(4) + t453;
t440 = -t445 * t554 + t559 * t448;
t479 = qJDD(5) - t480;
t516 = qJD(5) - t517;
t436 = -0.2e1 * qJD(6) * t505 + (t504 * t516 - t461) * qJ(6) + (t504 * t505 + t479) * pkin(5) + t440;
t483 = -mrSges(7,2) * t516 + mrSges(7,3) * t504;
t581 = m(7) * t436 + t479 * mrSges(7,1) + t516 * t483;
t434 = -mrSges(7,3) * t461 - t475 * t505 + t581;
t441 = t559 * t445 + t554 * t448;
t460 = -qJD(5) * t505 - t481 * t554 + t511 * t559;
t485 = pkin(5) * t516 - qJ(6) * t505;
t503 = t504 ^ 2;
t439 = -pkin(5) * t503 + qJ(6) * t460 + 0.2e1 * qJD(6) * t504 - t485 * t516 + t441;
t589 = t602 * t504 + t611 * t505 + t601 * t516;
t590 = -t610 * t504 - t602 * t505 - t600 * t516;
t608 = mrSges(6,1) * t440 + mrSges(7,1) * t436 - mrSges(6,2) * t441 - mrSges(7,2) * t439 + pkin(5) * t434 + t600 * t460 + t601 * t461 + t609 * t479 - t589 * t504 - t590 * t505;
t603 = -mrSges(6,2) - mrSges(7,2);
t595 = t551 * t557;
t594 = t551 * t562;
t476 = -mrSges(6,1) * t504 + mrSges(6,2) * t505;
t484 = -mrSges(6,2) * t516 + mrSges(6,3) * t504;
t428 = m(6) * t440 + mrSges(6,1) * t479 + t484 * t516 + (-t475 - t476) * t505 + (-mrSges(6,3) - mrSges(7,3)) * t461 + t581;
t580 = m(7) * t439 + t460 * mrSges(7,3) + t504 * t475;
t486 = mrSges(7,1) * t516 - mrSges(7,3) * t505;
t588 = -mrSges(6,1) * t516 + mrSges(6,3) * t505 - t486;
t430 = m(6) * t441 + mrSges(6,3) * t460 + t476 * t504 + t603 * t479 + t588 * t516 + t580;
t427 = -t428 * t554 + t559 * t430;
t497 = -mrSges(5,1) * t517 + mrSges(5,2) * t518;
t507 = mrSges(5,1) * t524 - mrSges(5,3) * t518;
t424 = m(5) * t450 - mrSges(5,2) * t511 + mrSges(5,3) * t480 + t497 * t517 - t507 * t524 + t427;
t449 = -t555 * t454 + t456 * t560;
t444 = -pkin(4) * t511 - pkin(12) * t523 + t518 * t498 - t449;
t442 = -pkin(5) * t460 - qJ(6) * t503 + t485 * t505 + qJDD(6) + t444;
t574 = -m(7) * t442 + t460 * mrSges(7,1) + t504 * t483;
t433 = -m(6) * t444 + t460 * mrSges(6,1) + t603 * t461 + t504 * t484 + t588 * t505 + t574;
t506 = -mrSges(5,2) * t524 + mrSges(5,3) * t517;
t432 = m(5) * t449 + mrSges(5,1) * t511 - mrSges(5,3) * t481 - t497 * t518 + t506 * t524 + t433;
t417 = t555 * t424 + t560 * t432;
t591 = -t600 * t504 - t601 * t505 - t609 * t516;
t514 = -mrSges(4,1) * t525 + mrSges(4,2) * t526;
t520 = mrSges(4,1) * t536 - mrSges(4,3) * t526;
t576 = t560 * t424 - t432 * t555;
t414 = m(4) * t463 - mrSges(4,2) * t527 + mrSges(4,3) * t512 + t514 * t525 - t520 * t536 + t576;
t519 = -mrSges(4,2) * t536 + mrSges(4,3) * t525;
t416 = m(4) * t473 - mrSges(4,1) * t512 + mrSges(4,2) * t513 - t519 * t525 + t520 * t526 + t417;
t426 = t559 * t428 + t430 * t554;
t566 = -m(5) * t453 + t480 * mrSges(5,1) - mrSges(5,2) * t481 + t517 * t506 - t507 * t518 - t426;
t421 = m(4) * t462 + mrSges(4,1) * t527 - mrSges(4,3) * t513 - t514 * t526 + t519 * t536 + t566;
t405 = t414 * t597 + t552 * t416 + t421 * t596;
t409 = t561 * t414 - t421 * t556;
t406 = t552 * t561 * t421 + t414 * t593 - t416 * t550;
t437 = mrSges(7,2) * t461 + t486 * t505 - t574;
t418 = -mrSges(6,1) * t444 + mrSges(6,3) * t441 - mrSges(7,1) * t442 + mrSges(7,3) * t439 - pkin(5) * t437 + qJ(6) * t580 + (-qJ(6) * t486 + t589) * t516 + t591 * t505 + (-qJ(6) * mrSges(7,2) + t600) * t479 + t602 * t461 + t610 * t460;
t425 = mrSges(6,2) * t444 + mrSges(7,2) * t442 - mrSges(6,3) * t440 - mrSges(7,3) * t436 - qJ(6) * t434 + t602 * t460 + t611 * t461 + t601 * t479 - t591 * t504 + t590 * t516;
t488 = Ifges(5,5) * t518 + Ifges(5,6) * t517 + Ifges(5,3) * t524;
t489 = Ifges(5,4) * t518 + Ifges(5,2) * t517 + Ifges(5,6) * t524;
t407 = mrSges(5,2) * t453 - mrSges(5,3) * t449 + Ifges(5,1) * t481 + Ifges(5,4) * t480 + Ifges(5,5) * t511 - pkin(12) * t426 - t418 * t554 + t425 * t559 + t488 * t517 - t489 * t524;
t490 = Ifges(5,1) * t518 + Ifges(5,4) * t517 + Ifges(5,5) * t524;
t410 = -mrSges(5,1) * t453 + mrSges(5,3) * t450 + Ifges(5,4) * t481 + Ifges(5,2) * t480 + Ifges(5,6) * t511 - pkin(4) * t426 - t518 * t488 + t524 * t490 - t608;
t508 = Ifges(4,5) * t526 + Ifges(4,6) * t525 + Ifges(4,3) * t536;
t509 = Ifges(4,4) * t526 + Ifges(4,2) * t525 + Ifges(4,6) * t536;
t402 = mrSges(4,2) * t473 - mrSges(4,3) * t462 + Ifges(4,1) * t513 + Ifges(4,4) * t512 + Ifges(4,5) * t527 - pkin(11) * t417 + t407 * t560 - t410 * t555 + t508 * t525 - t509 * t536;
t510 = Ifges(4,1) * t526 + Ifges(4,4) * t525 + Ifges(4,5) * t536;
t565 = mrSges(5,1) * t449 - mrSges(5,2) * t450 + Ifges(5,5) * t481 + Ifges(5,6) * t480 + Ifges(5,3) * t511 + pkin(4) * t433 + pkin(12) * t427 + t559 * t418 + t554 * t425 + t518 * t489 - t517 * t490;
t403 = -mrSges(4,1) * t473 + mrSges(4,3) * t463 + Ifges(4,4) * t513 + Ifges(4,2) * t512 + Ifges(4,6) * t527 - pkin(3) * t417 - t526 * t508 + t536 * t510 - t565;
t568 = pkin(10) * t409 + t402 * t556 + t403 * t561;
t544 = (-mrSges(3,1) * t562 + mrSges(3,2) * t557) * t586;
t541 = -mrSges(3,2) * t548 + mrSges(3,3) * t578;
t540 = mrSges(3,1) * t548 - mrSges(3,3) * t579;
t531 = -t542 * t551 - t604;
t530 = Ifges(3,5) * t548 + (Ifges(3,1) * t557 + Ifges(3,4) * t562) * t586;
t529 = Ifges(3,6) * t548 + (Ifges(3,4) * t557 + Ifges(3,2) * t562) * t586;
t528 = Ifges(3,3) * t548 + (Ifges(3,5) * t557 + Ifges(3,6) * t562) * t586;
t522 = -g(3) * t595 + t587;
t521 = -g(3) * t594 + t575;
t408 = m(3) * t522 - mrSges(3,2) * t547 + mrSges(3,3) * t546 - t540 * t548 + t544 * t578 + t409;
t404 = m(3) * t521 + mrSges(3,1) * t547 - mrSges(3,3) * t545 + t541 * t548 - t544 * t579 + t406;
t401 = mrSges(4,1) * t462 - mrSges(4,2) * t463 + Ifges(4,5) * t513 + Ifges(4,6) * t512 + Ifges(4,3) * t527 + pkin(3) * t566 + pkin(11) * t576 + t555 * t407 + t560 * t410 + t526 * t509 - t525 * t510;
t400 = mrSges(3,1) * t521 - mrSges(3,2) * t522 + Ifges(3,5) * t545 + Ifges(3,6) * t546 + Ifges(3,3) * t547 + pkin(2) * t406 + t401 * t552 + (t529 * t557 - t530 * t562) * t586 + t568 * t550;
t1 = [Ifges(2,3) * qJDD(1) + mrSges(2,1) * t577 - mrSges(2,2) * t573 + (t528 * t578 + mrSges(3,2) * t531 - mrSges(3,3) * t521 + Ifges(3,1) * t545 + Ifges(3,4) * t546 + Ifges(3,5) * t547 + t402 * t561 - t403 * t556 - t529 * t548 + (-t405 * t550 - t406 * t552) * pkin(10)) * t595 + (-mrSges(3,1) * t531 + mrSges(3,3) * t522 + Ifges(3,4) * t545 + Ifges(3,2) * t546 + Ifges(3,6) * t547 - pkin(2) * t405 - t401 * t550 - t528 * t579 + t530 * t548 + t568 * t552) * t594 + t553 * t400 + pkin(1) * ((t404 * t562 + t408 * t557) * t553 + (-m(3) * t531 + t546 * mrSges(3,1) - t545 * mrSges(3,2) + (-t540 * t557 + t541 * t562) * t586 - t405) * t551) + (-t404 * t557 + t408 * t562) * t607; t400; t401; t565; t608; t437;];
tauJ  = t1;
