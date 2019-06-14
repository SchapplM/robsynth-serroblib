% Calculate vector of inverse dynamics joint torques for with Newton-Euler
% S6RRRPRP11
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
%   pkin=[a2,a3,a4,a5,a6,alpha2,d1,d2,d3,d5]';
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
% Datum: 2019-05-07 09:21
% Revision: 89c353f7eff3bd693eda4e29f35b2761dbc3ada0 (2019-05-03)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauJ = S6RRRPRP11_invdynJ_fixb_snew_vp2(qJ, qJD, qJDD, g, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(6,1),zeros(3,1),zeros(10,1),zeros(7,1),zeros(7,3),zeros(7,6)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRRPRP11_invdynJ_fixb_snew_vp2: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RRRPRP11_invdynJ_fixb_snew_vp2: qJD has to be [6x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [6 1]), ...
  'S6RRRPRP11_invdynJ_fixb_snew_vp2: qJDD has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6RRRPRP11_invdynJ_fixb_snew_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S6RRRPRP11_invdynJ_fixb_snew_vp2: pkin has to be [10x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6RRRPRP11_invdynJ_fixb_snew_vp2: m has to be [7x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [7,3]), ...
  'S6RRRPRP11_invdynJ_fixb_snew_vp2: mrSges has to be [7x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [7 6]), ...
  'S6RRRPRP11_invdynJ_fixb_snew_vp2: Ifges has to be [7x6] (double)'); 

%% Symbolic Calculation
% From invdyn_fixb_NewtonEuler_linkframe_tauJ_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-05-07 09:08:23
% EndTime: 2019-05-07 09:08:36
% DurationCPUTime: 4.97s
% Computational Cost: add. (41218->323), mult. (88342->386), div. (0->0), fcn. (67317->10), ass. (0->137)
t602 = Ifges(6,4) + Ifges(7,4);
t617 = Ifges(6,2) + Ifges(7,2);
t612 = Ifges(6,6) + Ifges(7,6);
t616 = Ifges(4,1) + Ifges(5,2);
t615 = Ifges(6,1) + Ifges(7,1);
t601 = Ifges(4,5) - Ifges(5,4);
t614 = Ifges(6,5) + Ifges(7,5);
t613 = -Ifges(4,2) - Ifges(5,3);
t600 = Ifges(4,6) - Ifges(5,5);
t599 = -Ifges(5,6) - Ifges(4,4);
t611 = Ifges(4,3) + Ifges(5,1);
t610 = Ifges(6,3) + Ifges(7,3);
t555 = cos(pkin(6));
t551 = qJD(1) * t555 + qJD(2);
t557 = sin(qJ(3));
t558 = sin(qJ(2));
t554 = sin(pkin(6));
t585 = qJD(1) * t554;
t578 = t558 * t585;
t606 = cos(qJ(3));
t528 = -t551 * t606 + t557 * t578;
t561 = cos(qJ(2));
t584 = qJD(1) * t561;
t577 = t554 * t584;
t547 = -qJD(3) + t577;
t556 = sin(qJ(5));
t560 = cos(qJ(5));
t511 = t528 * t560 + t547 * t556;
t512 = t528 * t556 - t547 * t560;
t529 = t551 * t557 + t578 * t606;
t526 = qJD(5) + t529;
t609 = t511 * t617 + t512 * t602 + t526 * t612;
t582 = qJD(1) * qJD(2);
t542 = (-qJDD(1) * t561 + t558 * t582) * t554;
t540 = (-pkin(2) * t561 - pkin(9) * t558) * t585;
t549 = t551 ^ 2;
t550 = qJDD(1) * t555 + qJDD(2);
t563 = qJD(1) ^ 2;
t559 = sin(qJ(1));
t562 = cos(qJ(1));
t574 = -g(1) * t562 - g(2) * t559;
t605 = pkin(8) * t554;
t538 = -pkin(1) * t563 + qJDD(1) * t605 + t574;
t576 = g(1) * t559 - g(2) * t562;
t537 = qJDD(1) * pkin(1) + t563 * t605 + t576;
t596 = t537 * t555;
t586 = t538 * t561 + t558 * t596;
t473 = -t549 * pkin(2) + t550 * pkin(9) + (-g(3) * t558 + t540 * t584) * t554 + t586;
t541 = (qJDD(1) * t558 + t561 * t582) * t554;
t604 = t555 * g(3);
t474 = t542 * pkin(2) - t541 * pkin(9) - t604 + (-t537 + (pkin(2) * t558 - pkin(9) * t561) * t551 * qJD(1)) * t554;
t448 = t473 * t606 + t474 * t557;
t505 = pkin(3) * t528 - qJ(4) * t529;
t534 = qJDD(3) + t542;
t546 = t547 ^ 2;
t444 = pkin(3) * t546 - qJ(4) * t534 + 0.2e1 * qJD(4) * t547 + t505 * t528 - t448;
t501 = qJD(3) * t529 + t541 * t557 - t550 * t606;
t517 = pkin(4) * t529 + pkin(10) * t547;
t527 = t528 ^ 2;
t442 = -pkin(4) * t501 - pkin(10) * t527 - t517 * t547 - t444;
t458 = -qJD(5) * t512 + t501 * t560 - t534 * t556;
t459 = qJD(5) * t511 + t501 * t556 + t534 * t560;
t481 = -mrSges(7,2) * t526 + mrSges(7,3) * t511;
t482 = -mrSges(6,2) * t526 + mrSges(6,3) * t511;
t485 = mrSges(6,1) * t526 - mrSges(6,3) * t512;
t483 = pkin(5) * t526 - qJ(6) * t512;
t510 = t511 ^ 2;
t436 = -pkin(5) * t458 - qJ(6) * t510 + t483 * t512 + qJDD(6) + t442;
t484 = mrSges(7,1) * t526 - mrSges(7,3) * t512;
t579 = m(7) * t436 + mrSges(7,2) * t459 + t484 * t512;
t608 = -m(6) * t442 - t459 * mrSges(6,2) + (mrSges(6,1) + mrSges(7,1)) * t458 - t512 * t485 + (t481 + t482) * t511 - t579;
t447 = -t473 * t557 + t474 * t606;
t445 = -t534 * pkin(3) - t546 * qJ(4) + t505 * t529 + qJDD(4) - t447;
t502 = -qJD(3) * t528 + t541 * t606 + t550 * t557;
t597 = t528 * t547;
t439 = (t528 * t529 - t534) * pkin(10) + (t502 - t597) * pkin(4) + t445;
t594 = t554 * t561;
t503 = -g(3) * t594 - t538 * t558 + t561 * t596;
t472 = -t550 * pkin(2) - t549 * pkin(9) + t540 * t578 - t503;
t566 = (-t502 - t597) * qJ(4) + t472 + (-pkin(3) * t547 - 0.2e1 * qJD(4)) * t529;
t443 = -t527 * pkin(4) - t529 * t517 + (pkin(3) + pkin(10)) * t501 + t566;
t434 = t439 * t556 + t443 * t560;
t430 = -pkin(5) * t510 + qJ(6) * t458 + 0.2e1 * qJD(6) * t511 - t483 * t526 + t434;
t431 = -t458 * mrSges(7,1) - t511 * t481 + t579;
t498 = qJDD(5) + t502;
t476 = -mrSges(7,1) * t511 + mrSges(7,2) * t512;
t580 = m(7) * t430 + mrSges(7,3) * t458 + t476 * t511;
t591 = -t511 * t602 - t512 * t615 - t526 * t614;
t592 = -t511 * t612 - t512 * t614 - t526 * t610;
t408 = -mrSges(6,1) * t442 + mrSges(6,3) * t434 - mrSges(7,1) * t436 + mrSges(7,3) * t430 - pkin(5) * t431 + qJ(6) * t580 + (-qJ(6) * t484 - t591) * t526 + t592 * t512 + (-mrSges(7,2) * qJ(6) + t612) * t498 + t602 * t459 + t617 * t458;
t515 = mrSges(5,1) * t528 + mrSges(5,3) * t547;
t433 = t439 * t560 - t556 * t443;
t477 = -mrSges(6,1) * t511 + mrSges(6,2) * t512;
t428 = -0.2e1 * qJD(6) * t512 + (t511 * t526 - t459) * qJ(6) + (t511 * t512 + t498) * pkin(5) + t433;
t581 = m(7) * t428 + mrSges(7,1) * t498 + t481 * t526;
t420 = m(6) * t433 + t498 * mrSges(6,1) + t526 * t482 + (-t476 - t477) * t512 + (-mrSges(6,3) - mrSges(7,3)) * t459 + t581;
t422 = m(6) * t434 + t458 * mrSges(6,3) + t511 * t477 + (-t484 - t485) * t526 + (-mrSges(6,2) - mrSges(7,2)) * t498 + t580;
t415 = t560 * t420 + t556 * t422;
t507 = -mrSges(5,2) * t528 - mrSges(5,3) * t529;
t569 = -m(5) * t445 - mrSges(5,1) * t502 - t507 * t529 - t415;
t413 = t534 * mrSges(5,2) - t515 * t547 - t569;
t425 = -t459 * mrSges(7,3) - t512 * t476 + t581;
t414 = mrSges(6,2) * t442 + mrSges(7,2) * t436 - mrSges(6,3) * t433 - mrSges(7,3) * t428 - qJ(6) * t425 + t458 * t602 + t459 * t615 + t498 * t614 - t511 * t592 - t526 * t609;
t516 = mrSges(5,1) * t529 - mrSges(5,2) * t547;
t568 = -m(5) * t444 + mrSges(5,3) * t534 - t516 * t547 - t608;
t587 = t528 * t599 + t529 * t616 - t547 * t601;
t588 = t528 * t613 - t529 * t599 - t547 * t600;
t607 = -t600 * t501 + t601 * t502 + t587 * t528 + t588 * t529 + t611 * t534 + mrSges(4,1) * t447 - mrSges(4,2) * t448 + mrSges(5,2) * t445 - mrSges(5,3) * t444 - pkin(3) * t413 - pkin(10) * t415 + qJ(4) * (-t501 * mrSges(5,1) - t528 * t507 + t568) - t556 * t408 + t560 * t414;
t595 = t554 * t558;
t506 = mrSges(4,1) * t528 + mrSges(4,2) * t529;
t513 = mrSges(4,2) * t547 - mrSges(4,3) * t528;
t411 = m(4) * t447 - t502 * mrSges(4,3) - t529 * t506 + (-t513 + t515) * t547 + (mrSges(4,1) - mrSges(5,2)) * t534 + t569;
t514 = -mrSges(4,1) * t547 - mrSges(4,3) * t529;
t418 = t568 + (-t506 - t507) * t528 + t547 * t514 - t534 * mrSges(4,2) + m(4) * t448 + (-mrSges(4,3) - mrSges(5,1)) * t501;
t407 = t411 * t606 + t418 * t557;
t593 = -t420 * t556 + t422 * t560;
t589 = t528 * t600 - t529 * t601 + t547 * t611;
t575 = -t411 * t557 + t418 * t606;
t446 = t501 * pkin(3) + t566;
t572 = -m(5) * t446 + mrSges(5,2) * t501 + t515 * t528 - t593;
t567 = -m(4) * t472 - t501 * mrSges(4,1) - t528 * t513 + (-t514 + t516) * t529 + (-mrSges(4,2) + mrSges(5,3)) * t502 + t572;
t565 = mrSges(6,1) * t433 + mrSges(7,1) * t428 - mrSges(6,2) * t434 - mrSges(7,2) * t430 + pkin(5) * t425 + t458 * t612 + t459 * t614 + t498 * t610 + t591 * t511 + t512 * t609;
t539 = (-mrSges(3,1) * t561 + mrSges(3,2) * t558) * t585;
t536 = -mrSges(3,2) * t551 + mrSges(3,3) * t577;
t535 = mrSges(3,1) * t551 - mrSges(3,3) * t578;
t522 = -t554 * t537 - t604;
t521 = Ifges(3,5) * t551 + (Ifges(3,1) * t558 + Ifges(3,4) * t561) * t585;
t520 = Ifges(3,6) * t551 + (Ifges(3,4) * t558 + Ifges(3,2) * t561) * t585;
t519 = Ifges(3,3) * t551 + (Ifges(3,5) * t558 + Ifges(3,6) * t561) * t585;
t504 = -g(3) * t595 + t586;
t412 = -t502 * mrSges(5,3) - t529 * t516 - t572;
t409 = m(3) * t503 + t550 * mrSges(3,1) - t541 * mrSges(3,3) + t551 * t536 - t539 * t578 + t567;
t406 = m(3) * t504 - mrSges(3,2) * t550 - mrSges(3,3) * t542 - t535 * t551 + t539 * t577 + t575;
t405 = mrSges(5,1) * t445 + mrSges(4,2) * t472 - mrSges(4,3) * t447 - mrSges(5,3) * t446 + pkin(4) * t415 - qJ(4) * t412 + t501 * t599 + t502 * t616 + t528 * t589 + t534 * t601 + t547 * t588 + t565;
t404 = -mrSges(4,1) * t472 - mrSges(5,1) * t444 + mrSges(5,2) * t446 + mrSges(4,3) * t448 - pkin(3) * t412 - pkin(4) * t608 - pkin(10) * t593 - t560 * t408 - t556 * t414 + t501 * t613 - t502 * t599 + t529 * t589 + t534 * t600 - t547 * t587;
t403 = Ifges(3,5) * t541 - Ifges(3,6) * t542 + Ifges(3,3) * t550 + mrSges(3,1) * t503 - mrSges(3,2) * t504 + t557 * t405 + t606 * t404 + pkin(2) * t567 + pkin(9) * t575 + (t520 * t558 - t521 * t561) * t585;
t1 = [Ifges(2,3) * qJDD(1) + mrSges(2,1) * t576 - mrSges(2,2) * t574 + (mrSges(3,2) * t522 - mrSges(3,3) * t503 + Ifges(3,1) * t541 - Ifges(3,4) * t542 + Ifges(3,5) * t550 - pkin(9) * t407 - t404 * t557 + t405 * t606 + t519 * t577 - t520 * t551) * t595 + (-mrSges(3,1) * t522 + mrSges(3,3) * t504 + Ifges(3,4) * t541 - Ifges(3,2) * t542 + Ifges(3,6) * t550 - pkin(2) * t407 - t519 * t578 + t551 * t521 - t607) * t594 + t555 * t403 + pkin(1) * ((t406 * t558 + t409 * t561) * t555 + (-m(3) * t522 - t542 * mrSges(3,1) - t541 * mrSges(3,2) + (-t535 * t558 + t536 * t561) * t585 - t407) * t554) + (t406 * t561 - t409 * t558) * t605; t403; t607; t413; t565; t431;];
tauJ  = t1;
