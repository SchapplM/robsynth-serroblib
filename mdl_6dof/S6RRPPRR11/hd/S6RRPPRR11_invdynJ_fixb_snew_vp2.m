% Calculate vector of inverse dynamics joint torques for with Newton-Euler
% S6RRPPRR11
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
% pkin [11x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,alpha2,d1,d2,d5,d6,theta4]';
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
% Datum: 2019-05-06 12:05
% Revision: 89c353f7eff3bd693eda4e29f35b2761dbc3ada0 (2019-05-03)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauJ = S6RRPPRR11_invdynJ_fixb_snew_vp2(qJ, qJD, qJDD, g, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(6,1),zeros(3,1),zeros(11,1),zeros(7,1),zeros(7,3),zeros(7,6)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRPPRR11_invdynJ_fixb_snew_vp2: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RRPPRR11_invdynJ_fixb_snew_vp2: qJD has to be [6x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [6 1]), ...
  'S6RRPPRR11_invdynJ_fixb_snew_vp2: qJDD has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6RRPPRR11_invdynJ_fixb_snew_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [11 1]), ...
  'S6RRPPRR11_invdynJ_fixb_snew_vp2: pkin has to be [11x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6RRPPRR11_invdynJ_fixb_snew_vp2: m has to be [7x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [7,3]), ...
  'S6RRPPRR11_invdynJ_fixb_snew_vp2: mrSges has to be [7x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [7 6]), ...
  'S6RRPPRR11_invdynJ_fixb_snew_vp2: Ifges has to be [7x6] (double)'); 

%% Symbolic Calculation
% From invdyn_fixb_NewtonEuler_linkframe_tauJ_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-05-06 11:59:53
% EndTime: 2019-05-06 12:00:00
% DurationCPUTime: 7.06s
% Computational Cost: add. (69275->345), mult. (163427->433), div. (0->0), fcn. (121325->12), ass. (0->149)
t611 = -2 * qJD(3);
t610 = Ifges(3,1) + Ifges(4,2);
t603 = Ifges(3,4) + Ifges(4,6);
t602 = Ifges(3,5) - Ifges(4,4);
t609 = Ifges(3,2) + Ifges(4,3);
t601 = Ifges(3,6) - Ifges(4,5);
t608 = Ifges(3,3) + Ifges(4,1);
t560 = cos(pkin(6));
t552 = qJD(1) * t560 + qJD(2);
t563 = sin(qJ(2));
t558 = sin(pkin(6));
t591 = qJD(1) * t558;
t585 = t563 * t591;
t607 = (pkin(2) * t552 + t611) * t585;
t569 = qJD(1) ^ 2;
t564 = sin(qJ(1));
t568 = cos(qJ(1));
t579 = -g(1) * t568 - g(2) * t564;
t588 = qJDD(1) * t558;
t536 = -pkin(1) * t569 + pkin(8) * t588 + t579;
t567 = cos(qJ(2));
t597 = t558 * t563;
t583 = t564 * g(1) - g(2) * t568;
t606 = pkin(8) * t558;
t535 = qJDD(1) * pkin(1) + t569 * t606 + t583;
t599 = t535 * t560;
t494 = -g(3) * t597 + t567 * t536 + t563 * t599;
t537 = (-t567 * pkin(2) - t563 * qJ(3)) * t591;
t550 = t552 ^ 2;
t551 = qJDD(1) * t560 + qJDD(2);
t590 = qJD(1) * t567;
t584 = t558 * t590;
t481 = pkin(2) * t550 - t551 * qJ(3) - t537 * t584 + t552 * t611 - t494;
t605 = g(3) * t560;
t604 = mrSges(3,1) - mrSges(4,2);
t600 = -pkin(2) - qJ(4);
t598 = t558 ^ 2 * t569;
t596 = t558 * t567;
t532 = pkin(3) * t585 - qJ(4) * t552;
t540 = (qJD(2) * t590 + qJDD(1) * t563) * t558;
t541 = -qJD(2) * t585 + t567 * t588;
t586 = t567 ^ 2 * t598;
t461 = -pkin(3) * t586 - t605 - qJ(3) * t540 + t600 * t541 + (-t535 + (-qJ(3) * t552 * t567 - t532 * t563) * qJD(1)) * t558 + t607;
t592 = g(3) * t596 + t563 * t536;
t577 = -qJ(3) * t550 + t537 * t585 + qJDD(3) + t592;
t464 = pkin(3) * t540 + t600 * t551 + (-pkin(3) * t552 * t591 - qJ(4) * t563 * t598 - t599) * t567 + t577;
t557 = sin(pkin(11));
t559 = cos(pkin(11));
t523 = t552 * t559 - t557 * t584;
t446 = -0.2e1 * qJD(4) * t523 - t461 * t557 + t559 * t464;
t503 = -t541 * t557 + t551 * t559;
t522 = -t552 * t557 - t559 * t584;
t443 = (t522 * t585 - t503) * pkin(9) + (t522 * t523 + t540) * pkin(4) + t446;
t447 = 0.2e1 * qJD(4) * t522 + t559 * t461 + t557 * t464;
t502 = -t541 * t559 - t551 * t557;
t504 = pkin(4) * t585 - pkin(9) * t523;
t521 = t522 ^ 2;
t445 = -pkin(4) * t521 + pkin(9) * t502 - t504 * t585 + t447;
t562 = sin(qJ(5));
t566 = cos(qJ(5));
t440 = t562 * t443 + t566 * t445;
t497 = t522 * t562 + t523 * t566;
t470 = -qJD(5) * t497 + t502 * t566 - t503 * t562;
t496 = t522 * t566 - t523 * t562;
t479 = -mrSges(6,1) * t496 + mrSges(6,2) * t497;
t545 = qJD(5) + t585;
t487 = mrSges(6,1) * t545 - mrSges(6,3) * t497;
t529 = qJDD(5) + t540;
t480 = -pkin(5) * t496 - pkin(10) * t497;
t543 = t545 ^ 2;
t437 = -pkin(5) * t543 + pkin(10) * t529 + t480 * t496 + t440;
t460 = pkin(3) * t541 - qJ(4) * t586 + t552 * t532 + qJDD(4) - t481;
t449 = -pkin(4) * t502 - pkin(9) * t521 + t523 * t504 + t460;
t471 = qJD(5) * t496 + t502 * t562 + t503 * t566;
t441 = (-t496 * t545 - t471) * pkin(10) + (t497 * t545 - t470) * pkin(5) + t449;
t561 = sin(qJ(6));
t565 = cos(qJ(6));
t434 = -t437 * t561 + t441 * t565;
t484 = -t497 * t561 + t545 * t565;
t452 = qJD(6) * t484 + t471 * t565 + t529 * t561;
t485 = t497 * t565 + t545 * t561;
t465 = -mrSges(7,1) * t484 + mrSges(7,2) * t485;
t469 = qJDD(6) - t470;
t495 = qJD(6) - t496;
t472 = -mrSges(7,2) * t495 + mrSges(7,3) * t484;
t431 = m(7) * t434 + mrSges(7,1) * t469 - mrSges(7,3) * t452 - t465 * t485 + t472 * t495;
t435 = t437 * t565 + t441 * t561;
t451 = -qJD(6) * t485 - t471 * t561 + t529 * t565;
t473 = mrSges(7,1) * t495 - mrSges(7,3) * t485;
t432 = m(7) * t435 - mrSges(7,2) * t469 + mrSges(7,3) * t451 + t465 * t484 - t473 * t495;
t580 = -t431 * t561 + t565 * t432;
t418 = m(6) * t440 - mrSges(6,2) * t529 + mrSges(6,3) * t470 + t479 * t496 - t487 * t545 + t580;
t439 = t443 * t566 - t445 * t562;
t486 = -mrSges(6,2) * t545 + mrSges(6,3) * t496;
t436 = -pkin(5) * t529 - pkin(10) * t543 + t480 * t497 - t439;
t575 = -m(7) * t436 + t451 * mrSges(7,1) - mrSges(7,2) * t452 + t484 * t472 - t473 * t485;
t427 = m(6) * t439 + mrSges(6,1) * t529 - mrSges(6,3) * t471 - t479 * t497 + t486 * t545 + t575;
t414 = t562 * t418 + t566 * t427;
t421 = t565 * t431 + t561 * t432;
t595 = (t563 * t602 + t567 * t601) * t591 + t608 * t552;
t594 = (t563 * t603 + t567 * t609) * t591 + t601 * t552;
t593 = (t610 * t563 + t567 * t603) * t591 + t602 * t552;
t587 = t567 * t599;
t498 = -mrSges(5,1) * t522 + mrSges(5,2) * t523;
t500 = -mrSges(5,2) * t585 + mrSges(5,3) * t522;
t412 = m(5) * t446 + mrSges(5,1) * t540 - mrSges(5,3) * t503 - t498 * t523 + t500 * t585 + t414;
t501 = mrSges(5,1) * t585 - mrSges(5,3) * t523;
t581 = t566 * t418 - t427 * t562;
t413 = m(5) * t447 - mrSges(5,2) * t540 + mrSges(5,3) * t502 + t498 * t522 - t501 * t585 + t581;
t582 = -t557 * t412 + t559 * t413;
t511 = -t535 * t558 - t605;
t408 = t412 * t559 + t413 * t557;
t482 = -pkin(2) * t541 + (-t552 * t584 - t540) * qJ(3) + t511 + t607;
t533 = -mrSges(4,1) * t584 - mrSges(4,3) * t552;
t578 = -m(4) * t482 + t540 * mrSges(4,3) - t533 * t584 - t582;
t483 = -pkin(2) * t551 + t577 - t587;
t576 = -m(4) * t483 - t540 * mrSges(4,1) - t408;
t573 = m(6) * t449 - t470 * mrSges(6,1) + t471 * mrSges(6,2) - t496 * t486 + t497 * t487 + t421;
t453 = Ifges(7,5) * t485 + Ifges(7,6) * t484 + Ifges(7,3) * t495;
t455 = Ifges(7,1) * t485 + Ifges(7,4) * t484 + Ifges(7,5) * t495;
t424 = -mrSges(7,1) * t436 + mrSges(7,3) * t435 + Ifges(7,4) * t452 + Ifges(7,2) * t451 + Ifges(7,6) * t469 - t453 * t485 + t455 * t495;
t454 = Ifges(7,4) * t485 + Ifges(7,2) * t484 + Ifges(7,6) * t495;
t425 = mrSges(7,2) * t436 - mrSges(7,3) * t434 + Ifges(7,1) * t452 + Ifges(7,4) * t451 + Ifges(7,5) * t469 + t453 * t484 - t454 * t495;
t475 = Ifges(6,4) * t497 + Ifges(6,2) * t496 + Ifges(6,6) * t545;
t476 = Ifges(6,1) * t497 + Ifges(6,4) * t496 + Ifges(6,5) * t545;
t572 = mrSges(6,1) * t439 - mrSges(6,2) * t440 + Ifges(6,5) * t471 + Ifges(6,6) * t470 + Ifges(6,3) * t529 + pkin(5) * t575 + pkin(10) * t580 + t565 * t424 + t561 * t425 + t497 * t475 - t496 * t476;
t571 = mrSges(7,1) * t434 - mrSges(7,2) * t435 + Ifges(7,5) * t452 + Ifges(7,6) * t451 + Ifges(7,3) * t469 + t454 * t485 - t455 * t484;
t419 = m(5) * t460 - t502 * mrSges(5,1) + t503 * mrSges(5,2) - t522 * t500 + t523 * t501 + t573;
t534 = mrSges(4,1) * t585 + mrSges(4,2) * t552;
t538 = (t567 * mrSges(4,2) - t563 * mrSges(4,3)) * t591;
t570 = -m(4) * t481 + t551 * mrSges(4,3) + t552 * t534 + t538 * t584 + t419;
t539 = (-t567 * mrSges(3,1) + t563 * mrSges(3,2)) * t591;
t531 = -mrSges(3,2) * t552 + mrSges(3,3) * t584;
t530 = mrSges(3,1) * t552 - mrSges(3,3) * t585;
t493 = t587 - t592;
t490 = Ifges(5,1) * t523 + Ifges(5,4) * t522 + Ifges(5,5) * t585;
t489 = Ifges(5,4) * t523 + Ifges(5,2) * t522 + Ifges(5,6) * t585;
t488 = Ifges(5,5) * t523 + Ifges(5,6) * t522 + Ifges(5,3) * t585;
t474 = Ifges(6,5) * t497 + Ifges(6,6) * t496 + Ifges(6,3) * t545;
t415 = (mrSges(3,3) + mrSges(4,1)) * t541 - t552 * t530 - t551 * mrSges(3,2) + m(3) * t494 + t570 + t539 * t584;
t410 = -mrSges(6,1) * t449 + mrSges(6,3) * t440 + Ifges(6,4) * t471 + Ifges(6,2) * t470 + Ifges(6,6) * t529 - pkin(5) * t421 - t474 * t497 + t476 * t545 - t571;
t409 = mrSges(6,2) * t449 - mrSges(6,3) * t439 + Ifges(6,1) * t471 + Ifges(6,4) * t470 + Ifges(6,5) * t529 - pkin(10) * t421 - t424 * t561 + t425 * t565 + t474 * t496 - t475 * t545;
t407 = mrSges(4,2) * t551 + t533 * t552 + t538 * t585 - t576;
t406 = t541 * mrSges(4,2) - t534 * t585 - t578;
t405 = m(3) * t493 - mrSges(3,3) * t540 + (t531 - t533) * t552 + t604 * t551 + (-t538 - t539) * t585 + t576;
t404 = mrSges(5,2) * t460 - mrSges(5,3) * t446 + Ifges(5,1) * t503 + Ifges(5,4) * t502 + Ifges(5,5) * t540 - pkin(9) * t414 + t409 * t566 - t410 * t562 + t488 * t522 - t489 * t585;
t403 = -mrSges(5,1) * t460 + mrSges(5,3) * t447 + Ifges(5,4) * t503 + Ifges(5,2) * t502 + Ifges(5,6) * t540 - pkin(4) * t573 + pkin(9) * t581 + t562 * t409 + t566 * t410 - t523 * t488 + t490 * t585;
t402 = mrSges(3,1) * t493 - mrSges(3,2) * t494 + mrSges(4,2) * t483 - mrSges(4,3) * t481 + t559 * t404 - t557 * t403 - qJ(4) * t408 - pkin(2) * t407 + qJ(3) * t570 + t608 * t551 + (qJ(3) * mrSges(4,1) + t601) * t541 + t602 * t540 + (t594 * t563 - t593 * t567) * t591;
t1 = [Ifges(2,3) * qJDD(1) + mrSges(2,1) * t583 - mrSges(2,2) * t579 + (-t594 * t552 + t602 * t551 + t603 * t541 + (Ifges(5,3) + t610) * t540 + t595 * t584 - t522 * t490 + t523 * t489 + Ifges(5,6) * t502 + Ifges(5,5) * t503 + mrSges(3,2) * t511 - mrSges(3,3) * t493 - mrSges(4,3) * t482 + mrSges(4,1) * t483 + mrSges(5,1) * t446 - mrSges(5,2) * t447 + t572 + pkin(4) * t414 + pkin(3) * t408 - qJ(3) * t406) * t597 + (-mrSges(3,1) * t511 - mrSges(4,1) * t481 + mrSges(4,2) * t482 + mrSges(3,3) * t494 - pkin(2) * t406 + pkin(3) * t419 - qJ(4) * t582 - t559 * t403 - t557 * t404 + t603 * t540 + t609 * t541 + t601 * t551 + t593 * t552 - t595 * t585) * t596 + t560 * t402 + pkin(1) * ((t405 * t567 + t415 * t563) * t560 + (-m(3) * t511 - t540 * mrSges(3,2) + t604 * t541 + (t531 * t567 + (-t530 + t534) * t563) * t591 + t578) * t558) + (-t405 * t563 + t415 * t567) * t606; t402; t407; t419; t572; t571;];
tauJ  = t1;
