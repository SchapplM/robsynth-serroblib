% Calculate vector of inverse dynamics joint torques for with Newton-Euler
% S6RRRPPP1
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
% pkin [10x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,alpha4,d1,d2,d3,theta4]';
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
% Datum: 2019-05-07 04:05
% Revision: 89c353f7eff3bd693eda4e29f35b2761dbc3ada0 (2019-05-03)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauJ = S6RRRPPP1_invdynJ_fixb_snew_vp2(qJ, qJD, qJDD, g, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(6,1),zeros(3,1),zeros(10,1),zeros(7,1),zeros(7,3),zeros(7,6)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRRPPP1_invdynJ_fixb_snew_vp2: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RRRPPP1_invdynJ_fixb_snew_vp2: qJD has to be [6x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [6 1]), ...
  'S6RRRPPP1_invdynJ_fixb_snew_vp2: qJDD has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6RRRPPP1_invdynJ_fixb_snew_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S6RRRPPP1_invdynJ_fixb_snew_vp2: pkin has to be [10x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6RRRPPP1_invdynJ_fixb_snew_vp2: m has to be [7x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [7,3]), ...
  'S6RRRPPP1_invdynJ_fixb_snew_vp2: mrSges has to be [7x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [7 6]), ...
  'S6RRRPPP1_invdynJ_fixb_snew_vp2: Ifges has to be [7x6] (double)'); 

%% Symbolic Calculation
% From invdyn_fixb_NewtonEuler_linkframe_tauJ_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-05-07 03:58:08
% EndTime: 2019-05-07 03:58:17
% DurationCPUTime: 6.36s
% Computational Cost: add. (66845->330), mult. (141173->398), div. (0->0), fcn. (103193->10), ass. (0->138)
t597 = -2 * qJD(4);
t596 = Ifges(5,1) + Ifges(6,2) + Ifges(7,3);
t595 = Ifges(5,4) + Ifges(6,6) - Ifges(7,6);
t575 = Ifges(5,5) - Ifges(6,4) + Ifges(7,5);
t594 = -Ifges(5,2) - Ifges(6,3) - Ifges(7,2);
t574 = Ifges(5,6) - Ifges(6,5) - Ifges(7,4);
t593 = Ifges(5,3) + Ifges(6,1) + Ifges(7,1);
t542 = qJD(1) ^ 2;
t538 = sin(qJ(1));
t540 = cos(qJ(1));
t565 = t538 * g(1) - t540 * g(2);
t519 = -qJDD(1) * pkin(1) - t542 * pkin(8) - t565;
t537 = sin(qJ(2));
t588 = cos(qJ(2));
t567 = t588 * qJD(1);
t558 = qJD(2) * t567;
t527 = t537 * qJDD(1) + t558;
t578 = qJD(1) * t537;
t566 = qJD(2) * t578;
t528 = t588 * qJDD(1) - t566;
t482 = (-t558 - t527) * pkin(9) + (-t528 + t566) * pkin(2) + t519;
t557 = -g(1) * t540 - g(2) * t538;
t520 = -pkin(1) * t542 + qJDD(1) * pkin(8) + t557;
t510 = -t537 * g(3) + t588 * t520;
t526 = (-t588 * pkin(2) - pkin(9) * t537) * qJD(1);
t541 = qJD(2) ^ 2;
t491 = -t541 * pkin(2) + qJDD(2) * pkin(9) + t526 * t567 + t510;
t536 = sin(qJ(3));
t539 = cos(qJ(3));
t457 = t539 * t482 - t536 * t491;
t523 = qJD(2) * t539 - t536 * t578;
t499 = qJD(3) * t523 + qJDD(2) * t536 + t527 * t539;
t524 = qJD(2) * t536 + t539 * t578;
t535 = sin(pkin(6));
t582 = qJ(4) * t535;
t500 = -pkin(3) * t523 - t524 * t582;
t555 = t567 - qJD(3);
t553 = t555 * t535;
t584 = cos(pkin(6));
t547 = t523 * t584 - t553;
t501 = t547 * qJ(4);
t522 = qJDD(3) - t528;
t564 = qJ(4) * t584;
t436 = t522 * pkin(3) - t499 * t564 - t524 * t500 - t555 * t501 + t457;
t458 = t536 * t482 + t539 * t491;
t506 = -t555 * pkin(3) - t524 * t564;
t498 = -qJD(3) * t524 + qJDD(2) * t539 - t527 * t536;
t552 = t584 * t498 + t522 * t535;
t437 = t552 * qJ(4) + t523 * t500 + t555 * t506 + t458;
t509 = -t588 * g(3) - t537 * t520;
t490 = -qJDD(2) * pkin(2) - pkin(9) * t541 + t526 * t578 - t509;
t440 = -pkin(3) * t498 - t499 * t582 - t501 * t523 + t506 * t524 + t490;
t534 = sin(pkin(10));
t583 = cos(pkin(10));
t488 = t583 * t524 + t547 * t534;
t556 = t584 * t583;
t562 = t535 * t583;
t430 = t436 * t556 - t534 * t437 + t440 * t562 + t488 * t597;
t487 = -t523 * t556 + t524 * t534 + t583 * t553;
t460 = pkin(4) * t487 - qJ(5) * t488;
t481 = -t535 * t498 + t584 * t522;
t504 = t523 * t535 + t584 * t555;
t503 = t504 ^ 2;
t427 = -t481 * pkin(4) - t503 * qJ(5) + t488 * t460 + qJDD(5) - t430;
t462 = -mrSges(6,2) * t487 - mrSges(6,3) * t488;
t467 = t583 * t499 + t552 * t534;
t592 = -m(6) * t427 - t467 * mrSges(6,1) - t488 * t462;
t459 = -mrSges(7,2) * t488 + mrSges(7,3) * t487;
t581 = t487 * t504;
t589 = 2 * qJD(6);
t421 = t504 * t589 + (t487 * t488 - t481) * qJ(6) + (t467 - t581) * pkin(5) + t427;
t473 = -mrSges(7,1) * t487 - mrSges(7,2) * t504;
t559 = -m(7) * t421 + t481 * mrSges(7,3) - t504 * t473;
t419 = mrSges(7,1) * t467 + t459 * t488 - t559;
t472 = mrSges(6,1) * t487 + mrSges(6,3) * t504;
t417 = mrSges(6,2) * t481 - t472 * t504 + t419 - t592;
t466 = -t498 * t556 + t499 * t534 - t522 * t562;
t470 = pkin(5) * t488 + qJ(6) * t504;
t484 = t487 * t597;
t486 = t487 ^ 2;
t563 = t534 * t584;
t571 = t535 * t534 * t440 + t436 * t563 + t583 * t437;
t549 = pkin(4) * t503 - qJ(5) * t481 - t571;
t590 = -2 * qJD(5);
t423 = -pkin(5) * t466 - qJ(6) * t486 - t460 * t487 + qJDD(6) + t484 + (t590 - t470) * t504 - t549;
t426 = 0.2e1 * qJD(5) * t504 + ((2 * qJD(4)) + t460) * t487 + t549;
t431 = t484 + t571;
t474 = mrSges(6,1) * t488 - mrSges(6,2) * t504;
t471 = mrSges(7,1) * t488 + mrSges(7,3) * t504;
t573 = m(7) * t423 + t481 * mrSges(7,2) - t504 * t471;
t550 = -m(6) * t426 + t481 * mrSges(6,3) - t504 * t474 + t573;
t568 = t487 * t595 - t488 * t596 + t504 * t575;
t569 = t487 * t594 + t488 * t595 - t504 * t574;
t580 = -t459 - t462;
t402 = mrSges(5,1) * t430 - mrSges(5,2) * t431 + mrSges(6,2) * t427 - mrSges(6,3) * t426 + mrSges(7,2) * t423 - mrSges(7,3) * t421 - qJ(6) * t419 - pkin(4) * t417 + qJ(5) * t550 + t569 * t488 + t593 * t481 + t575 * t467 + (qJ(5) * t580 - t568) * t487 + (qJ(5) * (-mrSges(6,1) - mrSges(7,1)) - t574) * t466;
t432 = -t436 * t535 + t584 * t440 + qJDD(4);
t545 = (-t467 - t581) * qJ(5) + t432 + (-pkin(4) * t504 + t590) * t488;
t429 = pkin(4) * t466 + t545;
t425 = -pkin(5) * t486 + t487 * t589 - t470 * t488 + (pkin(4) + qJ(6)) * t466 + t545;
t572 = m(7) * t425 + t466 * mrSges(7,3) + t487 * t473;
t554 = m(6) * t429 - t467 * mrSges(6,3) - t488 * t474 + t572;
t418 = -mrSges(6,2) * t466 - mrSges(7,2) * t467 - t471 * t488 - t472 * t487 + t554;
t420 = -mrSges(7,1) * t466 - t459 * t487 + t573;
t570 = t487 * t574 - t488 * t575 + t504 * t593;
t403 = -mrSges(5,1) * t432 + mrSges(5,3) * t431 - mrSges(6,1) * t426 + mrSges(6,2) * t429 + mrSges(7,1) * t423 - mrSges(7,3) * t425 + pkin(5) * t420 - qJ(6) * t572 - pkin(4) * t418 + t568 * t504 + (qJ(6) * t471 + t570) * t488 + t574 * t481 + (qJ(6) * mrSges(7,2) + t595) * t467 + t594 * t466;
t461 = mrSges(5,1) * t487 + mrSges(5,2) * t488;
t579 = mrSges(5,2) * t504 - mrSges(5,3) * t487 - t472;
t585 = -mrSges(7,1) - mrSges(5,3);
t586 = mrSges(5,1) - mrSges(6,2);
t412 = m(5) * t430 - t579 * t504 + (-t459 - t461) * t488 + t586 * t481 + t585 * t467 + t559 + t592;
t469 = -mrSges(5,1) * t504 - mrSges(5,3) * t488;
t415 = m(5) * t431 - mrSges(5,2) * t481 + t469 * t504 + (-t461 + t580) * t487 + (-mrSges(6,1) + t585) * t466 + t550;
t416 = m(5) * t432 + (t469 - t471) * t488 + t579 * t487 + (mrSges(5,2) - mrSges(7,2)) * t467 + t586 * t466 + t554;
t406 = t412 * t556 + t415 * t563 - t535 * t416;
t407 = mrSges(6,1) * t427 + mrSges(7,1) * t421 + mrSges(5,2) * t432 - mrSges(7,2) * t425 - mrSges(5,3) * t430 - mrSges(6,3) * t429 + pkin(5) * t419 - qJ(5) * t418 - t466 * t595 + t467 * t596 + t575 * t481 + t570 * t487 + t569 * t504;
t410 = -t534 * t412 + t583 * t415;
t494 = Ifges(4,4) * t524 + Ifges(4,2) * t523 - Ifges(4,6) * t555;
t495 = Ifges(4,1) * t524 + Ifges(4,4) * t523 - Ifges(4,5) * t555;
t591 = mrSges(4,1) * t457 - mrSges(4,2) * t458 + Ifges(4,5) * t499 + Ifges(4,6) * t498 + Ifges(4,3) * t522 + pkin(3) * t406 + t584 * t402 + t524 * t494 - t523 * t495 + (qJ(4) * t410 + t583 * t403 + t534 * t407) * t535;
t405 = (t583 * t412 + t415 * t534) * t535 + t584 * t416;
t502 = -mrSges(4,1) * t523 + mrSges(4,2) * t524;
t507 = t555 * mrSges(4,2) + t523 * mrSges(4,3);
t404 = m(4) * t457 + t522 * mrSges(4,1) - t499 * mrSges(4,3) - t524 * t502 - t555 * t507 + t406;
t508 = -t555 * mrSges(4,1) - t524 * mrSges(4,3);
t409 = m(4) * t458 - t522 * mrSges(4,2) + t498 * mrSges(4,3) + t523 * t502 + t555 * t508 + t410;
t560 = -t536 * t404 + t539 * t409;
t401 = t539 * t404 + t536 * t409;
t544 = -m(4) * t490 + t498 * mrSges(4,1) - t499 * mrSges(4,2) + t523 * t507 - t524 * t508 - t405;
t531 = -qJD(2) * mrSges(3,2) + mrSges(3,3) * t567;
t530 = qJD(2) * mrSges(3,1) - mrSges(3,3) * t578;
t525 = (-mrSges(3,1) * t588 + mrSges(3,2) * t537) * qJD(1);
t518 = Ifges(3,5) * qJD(2) + (Ifges(3,1) * t537 + t588 * Ifges(3,4)) * qJD(1);
t517 = Ifges(3,6) * qJD(2) + (Ifges(3,4) * t537 + t588 * Ifges(3,2)) * qJD(1);
t516 = Ifges(3,3) * qJD(2) + (Ifges(3,5) * t537 + t588 * Ifges(3,6)) * qJD(1);
t493 = Ifges(4,5) * t524 + Ifges(4,6) * t523 - Ifges(4,3) * t555;
t400 = Ifges(4,1) * t499 + Ifges(4,4) * t498 + Ifges(4,5) * t522 + t523 * t493 + t555 * t494 + mrSges(4,2) * t490 - mrSges(4,3) * t457 + t583 * t407 - t534 * t403 + (-t535 * t405 - t584 * t406) * qJ(4);
t399 = -mrSges(4,1) * t490 + mrSges(4,3) * t458 + Ifges(4,4) * t499 + Ifges(4,2) * t498 + Ifges(4,6) * t522 - pkin(3) * t405 - t535 * t402 + t403 * t556 + t407 * t563 + t410 * t564 - t524 * t493 - t555 * t495;
t1 = [Ifges(2,3) * qJDD(1) + mrSges(2,1) * t565 - mrSges(2,2) * t557 + t537 * (mrSges(3,2) * t519 - mrSges(3,3) * t509 + Ifges(3,1) * t527 + Ifges(3,4) * t528 + Ifges(3,5) * qJDD(2) - pkin(9) * t401 - qJD(2) * t517 - t536 * t399 + t539 * t400 + t516 * t567) + t588 * (-mrSges(3,1) * t519 + mrSges(3,3) * t510 + Ifges(3,4) * t527 + Ifges(3,2) * t528 + Ifges(3,6) * qJDD(2) - pkin(2) * t401 + qJD(2) * t518 - t516 * t578 - t591) + pkin(1) * (-m(3) * t519 + t528 * mrSges(3,1) - t527 * mrSges(3,2) + (-t530 * t537 + t588 * t531) * qJD(1) - t401) + pkin(8) * (t588 * (m(3) * t510 - qJDD(2) * mrSges(3,2) + t528 * mrSges(3,3) - qJD(2) * t530 + t525 * t567 + t560) - t537 * (m(3) * t509 + qJDD(2) * mrSges(3,1) - t527 * mrSges(3,3) + qJD(2) * t531 - t525 * t578 + t544)); mrSges(3,1) * t509 - mrSges(3,2) * t510 + Ifges(3,5) * t527 + Ifges(3,6) * t528 + Ifges(3,3) * qJDD(2) + pkin(2) * t544 + pkin(9) * t560 + t539 * t399 + t536 * t400 + t517 * t578 - t518 * t567; t591; t416; t417; t420;];
tauJ  = t1;
