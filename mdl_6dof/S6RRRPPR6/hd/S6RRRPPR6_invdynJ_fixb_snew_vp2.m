% Calculate vector of inverse dynamics joint torques for with Newton-Euler
% S6RRRPPR6
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
%   pkin=[a2,a3,a4,a5,a6,alpha2,d1,d2,d3,d6,theta4]';
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
% Datum: 2019-05-07 05:37
% Revision: 89c353f7eff3bd693eda4e29f35b2761dbc3ada0 (2019-05-03)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauJ = S6RRRPPR6_invdynJ_fixb_snew_vp2(qJ, qJD, qJDD, g, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(6,1),zeros(3,1),zeros(11,1),zeros(7,1),zeros(7,3),zeros(7,6)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRRPPR6_invdynJ_fixb_snew_vp2: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RRRPPR6_invdynJ_fixb_snew_vp2: qJD has to be [6x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [6 1]), ...
  'S6RRRPPR6_invdynJ_fixb_snew_vp2: qJDD has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6RRRPPR6_invdynJ_fixb_snew_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [11 1]), ...
  'S6RRRPPR6_invdynJ_fixb_snew_vp2: pkin has to be [11x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6RRRPPR6_invdynJ_fixb_snew_vp2: m has to be [7x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [7,3]), ...
  'S6RRRPPR6_invdynJ_fixb_snew_vp2: mrSges has to be [7x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [7 6]), ...
  'S6RRRPPR6_invdynJ_fixb_snew_vp2: Ifges has to be [7x6] (double)'); 

%% Symbolic Calculation
% From invdyn_fixb_NewtonEuler_linkframe_tauJ_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-05-07 05:27:52
% EndTime: 2019-05-07 05:28:03
% DurationCPUTime: 6.88s
% Computational Cost: add. (88326->344), mult. (195708->430), div. (0->0), fcn. (155349->12), ass. (0->148)
t608 = -2 * qJD(4);
t607 = Ifges(5,1) + Ifges(6,2);
t606 = Ifges(6,1) + Ifges(5,3);
t600 = Ifges(5,4) + Ifges(6,6);
t599 = Ifges(5,5) - Ifges(6,4);
t605 = -Ifges(5,2) - Ifges(6,3);
t598 = Ifges(5,6) - Ifges(6,5);
t555 = sin(pkin(6));
t559 = sin(qJ(2));
t563 = cos(qJ(2));
t583 = qJD(1) * qJD(2);
t545 = (-qJDD(1) * t563 + t559 * t583) * t555;
t586 = qJD(1) * t555;
t543 = (-pkin(2) * t563 - pkin(9) * t559) * t586;
t556 = cos(pkin(6));
t551 = qJD(1) * t556 + qJD(2);
t549 = t551 ^ 2;
t550 = qJDD(1) * t556 + qJDD(2);
t585 = qJD(1) * t563;
t565 = qJD(1) ^ 2;
t560 = sin(qJ(1));
t564 = cos(qJ(1));
t576 = -g(1) * t564 - g(2) * t560;
t602 = pkin(8) * t555;
t541 = -pkin(1) * t565 + qJDD(1) * t602 + t576;
t579 = t560 * g(1) - g(2) * t564;
t540 = qJDD(1) * pkin(1) + t565 * t602 + t579;
t595 = t540 * t556;
t587 = t563 * t541 + t559 * t595;
t491 = -t549 * pkin(2) + t550 * pkin(9) + (-g(3) * t559 + t543 * t585) * t555 + t587;
t544 = (qJDD(1) * t559 + t563 * t583) * t555;
t601 = t556 * g(3);
t492 = t545 * pkin(2) - t544 * pkin(9) - t601 + (-t540 + (pkin(2) * t559 - pkin(9) * t563) * t551 * qJD(1)) * t555;
t558 = sin(qJ(3));
t562 = cos(qJ(3));
t456 = -t558 * t491 + t562 * t492;
t581 = t559 * t586;
t532 = t551 * t562 - t558 * t581;
t510 = qJD(3) * t532 + t544 * t562 + t550 * t558;
t533 = t551 * t558 + t562 * t581;
t537 = qJDD(3) + t545;
t580 = t555 * t585;
t547 = -qJD(3) + t580;
t445 = (-t532 * t547 - t510) * qJ(4) + (t532 * t533 + t537) * pkin(3) + t456;
t457 = t562 * t491 + t558 * t492;
t509 = -qJD(3) * t533 - t544 * t558 + t550 * t562;
t522 = -pkin(3) * t547 - qJ(4) * t533;
t531 = t532 ^ 2;
t448 = -pkin(3) * t531 + qJ(4) * t509 + t522 * t547 + t457;
t554 = sin(pkin(11));
t597 = cos(pkin(11));
t519 = t554 * t532 + t533 * t597;
t440 = t445 * t597 - t554 * t448 + t519 * t608;
t474 = t554 * t509 + t510 * t597;
t518 = -t532 * t597 + t533 * t554;
t485 = mrSges(5,1) * t518 + mrSges(5,2) * t519;
t496 = mrSges(6,1) * t518 + mrSges(6,3) * t547;
t498 = mrSges(5,2) * t547 - mrSges(5,3) * t518;
t484 = pkin(4) * t518 - qJ(5) * t519;
t546 = t547 ^ 2;
t439 = -t537 * pkin(4) - t546 * qJ(5) + t519 * t484 + qJDD(5) - t440;
t596 = t518 * t547;
t434 = (t518 * t519 - t537) * pkin(10) + (t474 - t596) * pkin(5) + t439;
t473 = -t509 * t597 + t510 * t554;
t500 = pkin(5) * t519 + pkin(10) * t547;
t517 = t518 ^ 2;
t593 = t555 * t563;
t511 = -g(3) * t593 - t559 * t541 + t563 * t595;
t490 = -t550 * pkin(2) - t549 * pkin(9) + t543 * t581 - t511;
t450 = -t509 * pkin(3) - t531 * qJ(4) + t533 * t522 + qJDD(4) + t490;
t603 = -2 * qJD(5);
t567 = (-t474 - t596) * qJ(5) + t450 + (-t547 * pkin(4) + t603) * t519;
t437 = t567 - t519 * t500 - t517 * pkin(5) + (pkin(4) + pkin(10)) * t473;
t557 = sin(qJ(6));
t561 = cos(qJ(6));
t432 = t434 * t561 - t437 * t557;
t494 = t518 * t561 + t547 * t557;
t455 = qJD(6) * t494 + t473 * t557 + t537 * t561;
t495 = t518 * t557 - t547 * t561;
t464 = -mrSges(7,1) * t494 + mrSges(7,2) * t495;
t516 = qJD(6) + t519;
t465 = -mrSges(7,2) * t516 + mrSges(7,3) * t494;
t472 = qJDD(6) + t474;
t429 = m(7) * t432 + mrSges(7,1) * t472 - mrSges(7,3) * t455 - t464 * t495 + t465 * t516;
t433 = t434 * t557 + t437 * t561;
t454 = -qJD(6) * t495 + t473 * t561 - t537 * t557;
t466 = mrSges(7,1) * t516 - mrSges(7,3) * t495;
t430 = m(7) * t433 - mrSges(7,2) * t472 + mrSges(7,3) * t454 + t464 * t494 - t466 * t516;
t420 = t561 * t429 + t557 * t430;
t486 = -mrSges(6,2) * t518 - mrSges(6,3) * t519;
t571 = -m(6) * t439 - t474 * mrSges(6,1) - t519 * t486 - t420;
t415 = m(5) * t440 - t474 * mrSges(5,3) - t519 * t485 + (t496 - t498) * t547 + (mrSges(5,1) - mrSges(6,2)) * t537 + t571;
t514 = t518 * t608;
t591 = t554 * t445 + t597 * t448;
t441 = t514 + t591;
t499 = -mrSges(5,1) * t547 - mrSges(5,3) * t519;
t573 = t546 * pkin(4) - t537 * qJ(5) - t591;
t438 = 0.2e1 * qJD(5) * t547 + ((2 * qJD(4)) + t484) * t518 + t573;
t497 = mrSges(6,1) * t519 - mrSges(6,2) * t547;
t436 = -t473 * pkin(5) - t517 * pkin(10) - t518 * t484 + t514 + (t603 - t500) * t547 - t573;
t572 = -m(7) * t436 + t454 * mrSges(7,1) - t455 * mrSges(7,2) + t494 * t465 - t495 * t466;
t569 = -m(6) * t438 + t537 * mrSges(6,3) - t547 * t497 - t572;
t425 = m(5) * t441 - t537 * mrSges(5,2) + t547 * t499 + (-t485 - t486) * t518 + (-mrSges(5,3) - mrSges(6,1)) * t473 + t569;
t413 = t597 * t415 + t554 * t425;
t418 = t537 * mrSges(6,2) - t547 * t496 - t571;
t458 = Ifges(7,5) * t495 + Ifges(7,6) * t494 + Ifges(7,3) * t516;
t460 = Ifges(7,1) * t495 + Ifges(7,4) * t494 + Ifges(7,5) * t516;
t423 = -mrSges(7,1) * t436 + mrSges(7,3) * t433 + Ifges(7,4) * t455 + Ifges(7,2) * t454 + Ifges(7,6) * t472 - t458 * t495 + t460 * t516;
t459 = Ifges(7,4) * t495 + Ifges(7,2) * t494 + Ifges(7,6) * t516;
t424 = mrSges(7,2) * t436 - mrSges(7,3) * t432 + Ifges(7,1) * t455 + Ifges(7,4) * t454 + Ifges(7,5) * t472 + t458 * t494 - t459 * t516;
t504 = Ifges(4,4) * t533 + Ifges(4,2) * t532 - Ifges(4,6) * t547;
t505 = Ifges(4,1) * t533 + Ifges(4,4) * t532 - Ifges(4,5) * t547;
t588 = -t600 * t518 + t607 * t519 - t599 * t547;
t589 = t605 * t518 + t600 * t519 - t598 * t547;
t604 = -t473 * t598 + t474 * t599 + t518 * t588 + t519 * t589 + (Ifges(4,3) + t606) * t537 + mrSges(4,1) * t456 + mrSges(5,1) * t440 + mrSges(6,2) * t439 + Ifges(4,5) * t510 + Ifges(4,6) * t509 + pkin(3) * t413 + qJ(5) * (-t473 * mrSges(6,1) - t518 * t486 + t569) + t533 * t504 + t561 * t424 - mrSges(4,2) * t457 - mrSges(5,2) * t441 - mrSges(6,3) * t438 - pkin(4) * t418 - pkin(10) * t420 - t532 * t505 - t557 * t423;
t594 = t555 * t559;
t520 = -mrSges(4,1) * t532 + mrSges(4,2) * t533;
t521 = mrSges(4,2) * t547 + mrSges(4,3) * t532;
t411 = m(4) * t456 + mrSges(4,1) * t537 - mrSges(4,3) * t510 - t520 * t533 - t521 * t547 + t413;
t523 = -mrSges(4,1) * t547 - mrSges(4,3) * t533;
t577 = -t415 * t554 + t597 * t425;
t412 = m(4) * t457 - mrSges(4,2) * t537 + mrSges(4,3) * t509 + t520 * t532 + t523 * t547 + t577;
t406 = t562 * t411 + t558 * t412;
t592 = -t557 * t429 + t561 * t430;
t590 = t598 * t518 - t599 * t519 + t606 * t547;
t578 = -t411 * t558 + t562 * t412;
t443 = t473 * pkin(4) + t567;
t419 = m(6) * t443 - t473 * mrSges(6,2) - t474 * mrSges(6,3) - t518 * t496 - t519 * t497 + t592;
t570 = mrSges(7,1) * t432 - mrSges(7,2) * t433 + Ifges(7,5) * t455 + Ifges(7,6) * t454 + Ifges(7,3) * t472 + t495 * t459 - t494 * t460;
t417 = m(5) * t450 + t473 * mrSges(5,1) + t474 * mrSges(5,2) + t518 * t498 + t519 * t499 + t419;
t568 = -m(4) * t490 + t509 * mrSges(4,1) - t510 * mrSges(4,2) + t532 * t521 - t533 * t523 - t417;
t542 = (-mrSges(3,1) * t563 + mrSges(3,2) * t559) * t586;
t539 = -mrSges(3,2) * t551 + mrSges(3,3) * t580;
t538 = mrSges(3,1) * t551 - mrSges(3,3) * t581;
t527 = -t555 * t540 - t601;
t526 = Ifges(3,5) * t551 + (Ifges(3,1) * t559 + Ifges(3,4) * t563) * t586;
t525 = Ifges(3,6) * t551 + (Ifges(3,4) * t559 + Ifges(3,2) * t563) * t586;
t524 = Ifges(3,3) * t551 + (Ifges(3,5) * t559 + Ifges(3,6) * t563) * t586;
t512 = -g(3) * t594 + t587;
t503 = Ifges(4,5) * t533 + Ifges(4,6) * t532 - Ifges(4,3) * t547;
t416 = m(3) * t511 + t550 * mrSges(3,1) - t544 * mrSges(3,3) + t551 * t539 - t542 * t581 + t568;
t407 = mrSges(6,1) * t439 + mrSges(5,2) * t450 - mrSges(5,3) * t440 - mrSges(6,3) * t443 + pkin(5) * t420 - qJ(5) * t419 - t600 * t473 + t607 * t474 + t590 * t518 + t599 * t537 + t589 * t547 + t570;
t405 = m(3) * t512 - mrSges(3,2) * t550 - mrSges(3,3) * t545 - t538 * t551 + t542 * t580 + t578;
t404 = -mrSges(5,1) * t450 - mrSges(6,1) * t438 + mrSges(6,2) * t443 + mrSges(5,3) * t441 - pkin(4) * t419 - pkin(5) * t572 - pkin(10) * t592 - t561 * t423 - t557 * t424 + t605 * t473 + t600 * t474 + t590 * t519 + t598 * t537 - t588 * t547;
t403 = mrSges(4,2) * t490 - mrSges(4,3) * t456 + Ifges(4,1) * t510 + Ifges(4,4) * t509 + Ifges(4,5) * t537 - qJ(4) * t413 - t554 * t404 + t407 * t597 + t532 * t503 + t547 * t504;
t402 = -mrSges(4,1) * t490 + mrSges(4,3) * t457 + Ifges(4,4) * t510 + Ifges(4,2) * t509 + Ifges(4,6) * t537 - pkin(3) * t417 + qJ(4) * t577 + t404 * t597 + t554 * t407 - t533 * t503 - t547 * t505;
t401 = Ifges(3,5) * t544 - Ifges(3,6) * t545 + Ifges(3,3) * t550 + mrSges(3,1) * t511 - mrSges(3,2) * t512 + t558 * t403 + t562 * t402 + pkin(2) * t568 + pkin(9) * t578 + (t525 * t559 - t526 * t563) * t586;
t1 = [Ifges(2,3) * qJDD(1) + mrSges(2,1) * t579 - mrSges(2,2) * t576 + (mrSges(3,2) * t527 - mrSges(3,3) * t511 + Ifges(3,1) * t544 - Ifges(3,4) * t545 + Ifges(3,5) * t550 - pkin(9) * t406 - t402 * t558 + t403 * t562 + t524 * t580 - t525 * t551) * t594 + (-mrSges(3,1) * t527 + mrSges(3,3) * t512 + Ifges(3,4) * t544 - Ifges(3,2) * t545 + Ifges(3,6) * t550 - pkin(2) * t406 - t524 * t581 + t551 * t526 - t604) * t593 + t556 * t401 + pkin(1) * ((t405 * t559 + t416 * t563) * t556 + (-m(3) * t527 - t545 * mrSges(3,1) - t544 * mrSges(3,2) + (-t538 * t559 + t539 * t563) * t586 - t406) * t555) + (t405 * t563 - t416 * t559) * t602; t401; t604; t417; t418; t570;];
tauJ  = t1;
