% Calculate vector of inverse dynamics joint torques for with Newton-Euler
% S6RPRRRP11
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
%   pkin=[a2,a3,a4,a5,a6,alpha2,alpha3,d1,d3,d4,d5,theta2]';
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
% Datum: 2019-05-06 02:13
% Revision: 89c353f7eff3bd693eda4e29f35b2761dbc3ada0 (2019-05-03)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauJ = S6RPRRRP11_invdynJ_fixb_snew_vp2(qJ, qJD, qJDD, g, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(6,1),zeros(3,1),zeros(12,1),zeros(7,1),zeros(7,3),zeros(7,6)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPRRRP11_invdynJ_fixb_snew_vp2: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RPRRRP11_invdynJ_fixb_snew_vp2: qJD has to be [6x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [6 1]), ...
  'S6RPRRRP11_invdynJ_fixb_snew_vp2: qJDD has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6RPRRRP11_invdynJ_fixb_snew_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [12 1]), ...
  'S6RPRRRP11_invdynJ_fixb_snew_vp2: pkin has to be [12x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6RPRRRP11_invdynJ_fixb_snew_vp2: m has to be [7x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [7,3]), ...
  'S6RPRRRP11_invdynJ_fixb_snew_vp2: mrSges has to be [7x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [7 6]), ...
  'S6RPRRRP11_invdynJ_fixb_snew_vp2: Ifges has to be [7x6] (double)'); 

%% Symbolic Calculation
% From invdyn_fixb_NewtonEuler_linkframe_tauJ_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-05-06 02:03:29
% EndTime: 2019-05-06 02:03:43
% DurationCPUTime: 13.29s
% Computational Cost: add. (178830->334), mult. (555375->449), div. (0->0), fcn. (473230->14), ass. (0->159)
t557 = sin(pkin(12));
t559 = sin(pkin(6));
t560 = cos(pkin(12));
t562 = cos(pkin(6));
t565 = sin(qJ(3));
t561 = cos(pkin(7));
t569 = cos(qJ(3));
t604 = t561 * t569;
t558 = sin(pkin(7));
t608 = t558 * t569;
t575 = t559 * (-t557 * t565 + t560 * t604) + t562 * t608;
t534 = t575 * qJD(1);
t605 = t561 * t565;
t609 = t558 * t565;
t577 = t562 * t609 + (t557 * t569 + t560 * t605) * t559;
t535 = t577 * qJD(1);
t523 = -t535 * qJD(3) + qJDD(1) * t575;
t625 = Ifges(6,1) + Ifges(7,1);
t617 = Ifges(6,4) + Ifges(7,4);
t616 = Ifges(6,5) + Ifges(7,5);
t624 = Ifges(6,2) + Ifges(7,2);
t615 = Ifges(6,6) + Ifges(7,6);
t623 = Ifges(6,3) + Ifges(7,3);
t606 = t559 * t561;
t546 = (t558 * t562 + t560 * t606) * qJD(1) * pkin(9);
t571 = qJD(1) ^ 2;
t566 = sin(qJ(1));
t570 = cos(qJ(1));
t588 = -g(1) * t570 - g(2) * t566;
t612 = qJ(2) * t559;
t550 = -pkin(1) * t571 + qJDD(1) * t612 + t588;
t619 = pkin(9) * t558;
t584 = -pkin(2) * t560 - t557 * t619;
t599 = qJD(1) * t559;
t613 = pkin(9) * qJDD(1);
t579 = qJD(1) * t584 * t599 + t561 * t613;
t593 = qJD(2) * t599;
t607 = t559 * t560;
t592 = t566 * g(1) - g(2) * t570;
t549 = qJDD(1) * pkin(1) + t571 * t612 + t592;
t611 = t549 * t562;
t585 = -g(3) * t607 - 0.2e1 * t557 * t593 + t560 * t611;
t504 = (pkin(2) * qJDD(1) + qJD(1) * t546) * t562 + (-t559 * t579 - t550) * t557 + t585;
t620 = pkin(9) * t557;
t551 = (pkin(2) * t562 - t606 * t620) * qJD(1);
t621 = 0.2e1 * t560;
t594 = t560 * t550 + t557 * t611 + t593 * t621;
t505 = (-qJD(1) * t551 + t558 * t613) * t562 + (-g(3) * t557 + t560 * t579) * t559 + t594;
t591 = -g(3) * t562 + qJDD(2);
t514 = (-t549 + t584 * qJDD(1) + (-t546 * t560 + t551 * t557) * qJD(1)) * t559 + t591;
t466 = -t565 * t505 + (t504 * t561 + t514 * t558) * t569;
t524 = t534 * qJD(3) + qJDD(1) * t577;
t580 = -t558 * t607 + t561 * t562;
t547 = qJD(1) * t580 + qJD(3);
t564 = sin(qJ(4));
t568 = cos(qJ(4));
t528 = -t535 * t564 + t547 * t568;
t544 = qJDD(1) * t580 + qJDD(3);
t500 = qJD(4) * t528 + t524 * t568 + t544 * t564;
t529 = t535 * t568 + t547 * t564;
t533 = qJD(4) - t534;
t563 = sin(qJ(5));
t567 = cos(qJ(5));
t512 = -t529 * t563 + t533 * t567;
t520 = qJDD(4) - t523;
t472 = qJD(5) * t512 + t500 * t567 + t520 * t563;
t513 = t529 * t567 + t533 * t563;
t484 = -mrSges(7,1) * t512 + mrSges(7,2) * t513;
t467 = t504 * t605 + t569 * t505 + t514 * t609;
t522 = -pkin(3) * t534 - pkin(10) * t535;
t543 = t547 ^ 2;
t463 = -pkin(3) * t543 + pkin(10) * t544 + t522 * t534 + t467;
t482 = -t558 * t504 + t561 * t514;
t465 = (-t534 * t547 - t524) * pkin(10) + (t535 * t547 - t523) * pkin(3) + t482;
t459 = t568 * t463 + t564 * t465;
t507 = -pkin(4) * t528 - pkin(11) * t529;
t532 = t533 ^ 2;
t454 = -pkin(4) * t532 + pkin(11) * t520 + t507 * t528 + t459;
t462 = -pkin(3) * t544 - pkin(10) * t543 + t535 * t522 - t466;
t499 = -qJD(4) * t529 - t524 * t564 + t544 * t568;
t457 = (-t528 * t533 - t500) * pkin(11) + (t529 * t533 - t499) * pkin(4) + t462;
t449 = -t454 * t563 + t567 * t457;
t497 = qJDD(5) - t499;
t525 = qJD(5) - t528;
t445 = -0.2e1 * qJD(6) * t513 + (t512 * t525 - t472) * qJ(6) + (t512 * t513 + t497) * pkin(5) + t449;
t487 = -mrSges(7,2) * t525 + mrSges(7,3) * t512;
t596 = m(7) * t445 + t497 * mrSges(7,1) + t525 * t487;
t443 = -mrSges(7,3) * t472 - t484 * t513 + t596;
t450 = t567 * t454 + t563 * t457;
t471 = -qJD(5) * t513 - t500 * t563 + t520 * t567;
t489 = pkin(5) * t525 - qJ(6) * t513;
t511 = t512 ^ 2;
t448 = -pkin(5) * t511 + qJ(6) * t471 + 0.2e1 * qJD(6) * t512 - t489 * t525 + t450;
t601 = t617 * t512 + t625 * t513 + t616 * t525;
t602 = -t624 * t512 - t617 * t513 - t615 * t525;
t622 = mrSges(6,1) * t449 + mrSges(7,1) * t445 - mrSges(6,2) * t450 - mrSges(7,2) * t448 + pkin(5) * t443 + t471 * t615 + t472 * t616 + t623 * t497 - t512 * t601 - t513 * t602;
t618 = -mrSges(6,2) - mrSges(7,2);
t610 = t557 * t559;
t485 = -mrSges(6,1) * t512 + mrSges(6,2) * t513;
t488 = -mrSges(6,2) * t525 + mrSges(6,3) * t512;
t437 = m(6) * t449 + mrSges(6,1) * t497 + t488 * t525 + (-t484 - t485) * t513 + (-mrSges(6,3) - mrSges(7,3)) * t472 + t596;
t595 = m(7) * t448 + t471 * mrSges(7,3) + t512 * t484;
t490 = mrSges(7,1) * t525 - mrSges(7,3) * t513;
t600 = -mrSges(6,1) * t525 + mrSges(6,3) * t513 - t490;
t439 = m(6) * t450 + mrSges(6,3) * t471 + t485 * t512 + t497 * t618 + t525 * t600 + t595;
t436 = -t437 * t563 + t567 * t439;
t506 = -mrSges(5,1) * t528 + mrSges(5,2) * t529;
t516 = mrSges(5,1) * t533 - mrSges(5,3) * t529;
t433 = m(5) * t459 - mrSges(5,2) * t520 + mrSges(5,3) * t499 + t506 * t528 - t516 * t533 + t436;
t458 = -t564 * t463 + t465 * t568;
t453 = -pkin(4) * t520 - pkin(11) * t532 + t529 * t507 - t458;
t451 = -pkin(5) * t471 - qJ(6) * t511 + t489 * t513 + qJDD(6) + t453;
t589 = -m(7) * t451 + t471 * mrSges(7,1) + t512 * t487;
t442 = -m(6) * t453 + t471 * mrSges(6,1) + t472 * t618 + t512 * t488 + t513 * t600 + t589;
t515 = -mrSges(5,2) * t533 + mrSges(5,3) * t528;
t441 = m(5) * t458 + mrSges(5,1) * t520 - mrSges(5,3) * t500 - t506 * t529 + t515 * t533 + t442;
t426 = t564 * t433 + t568 * t441;
t603 = -t615 * t512 - t616 * t513 - t623 * t525;
t521 = -mrSges(4,1) * t534 + mrSges(4,2) * t535;
t531 = mrSges(4,1) * t547 - mrSges(4,3) * t535;
t590 = t568 * t433 - t441 * t564;
t423 = m(4) * t467 - mrSges(4,2) * t544 + mrSges(4,3) * t523 + t521 * t534 - t531 * t547 + t590;
t530 = -mrSges(4,2) * t547 + mrSges(4,3) * t534;
t425 = m(4) * t482 - mrSges(4,1) * t523 + mrSges(4,2) * t524 - t530 * t534 + t531 * t535 + t426;
t435 = t437 * t567 + t439 * t563;
t573 = -m(5) * t462 + t499 * mrSges(5,1) - mrSges(5,2) * t500 + t528 * t515 - t516 * t529 - t435;
t430 = m(4) * t466 + mrSges(4,1) * t544 - mrSges(4,3) * t524 - t521 * t535 + t530 * t547 + t573;
t414 = t423 * t609 + t561 * t425 + t430 * t608;
t418 = t569 * t423 - t430 * t565;
t415 = t423 * t605 - t425 * t558 + t430 * t604;
t587 = -mrSges(3,1) * t560 + mrSges(3,2) * t557;
t583 = mrSges(3,1) * t562 - mrSges(3,3) * t610;
t582 = -mrSges(3,2) * t562 + mrSges(3,3) * t607;
t446 = mrSges(7,2) * t472 + t490 * t513 - t589;
t427 = -mrSges(6,1) * t453 + mrSges(6,3) * t450 - mrSges(7,1) * t451 + mrSges(7,3) * t448 - pkin(5) * t446 + qJ(6) * t595 + (-qJ(6) * t490 + t601) * t525 + t603 * t513 + (-mrSges(7,2) * qJ(6) + t615) * t497 + t617 * t472 + t624 * t471;
t434 = mrSges(6,2) * t453 + mrSges(7,2) * t451 - mrSges(6,3) * t449 - mrSges(7,3) * t445 - qJ(6) * t443 + t617 * t471 + t625 * t472 + t616 * t497 - t603 * t512 + t602 * t525;
t493 = Ifges(5,5) * t529 + Ifges(5,6) * t528 + Ifges(5,3) * t533;
t494 = Ifges(5,4) * t529 + Ifges(5,2) * t528 + Ifges(5,6) * t533;
t416 = mrSges(5,2) * t462 - mrSges(5,3) * t458 + Ifges(5,1) * t500 + Ifges(5,4) * t499 + Ifges(5,5) * t520 - pkin(11) * t435 - t427 * t563 + t434 * t567 + t493 * t528 - t494 * t533;
t495 = Ifges(5,1) * t529 + Ifges(5,4) * t528 + Ifges(5,5) * t533;
t419 = -mrSges(5,1) * t462 + mrSges(5,3) * t459 + Ifges(5,4) * t500 + Ifges(5,2) * t499 + Ifges(5,6) * t520 - pkin(4) * t435 - t529 * t493 + t533 * t495 - t622;
t517 = Ifges(4,5) * t535 + Ifges(4,6) * t534 + Ifges(4,3) * t547;
t518 = Ifges(4,4) * t535 + Ifges(4,2) * t534 + Ifges(4,6) * t547;
t410 = mrSges(4,2) * t482 - mrSges(4,3) * t466 + Ifges(4,1) * t524 + Ifges(4,4) * t523 + Ifges(4,5) * t544 - pkin(10) * t426 + t416 * t568 - t419 * t564 + t517 * t534 - t518 * t547;
t519 = Ifges(4,1) * t535 + Ifges(4,4) * t534 + Ifges(4,5) * t547;
t572 = mrSges(5,1) * t458 - mrSges(5,2) * t459 + Ifges(5,5) * t500 + Ifges(5,6) * t499 + Ifges(5,3) * t520 + pkin(4) * t442 + pkin(11) * t436 + t567 * t427 + t563 * t434 + t529 * t494 - t528 * t495;
t411 = -mrSges(4,1) * t482 + mrSges(4,3) * t467 + Ifges(4,4) * t524 + Ifges(4,2) * t523 + Ifges(4,6) * t544 - pkin(3) * t426 - t535 * t517 + t547 * t519 - t572;
t578 = pkin(9) * t418 + t410 * t565 + t411 * t569;
t553 = t582 * qJD(1);
t552 = t583 * qJD(1);
t548 = t587 * t599;
t536 = -t549 * t559 + t591;
t527 = -g(3) * t610 + t594;
t526 = -t550 * t557 + t585;
t417 = m(3) * t527 + t582 * qJDD(1) + (t548 * t607 - t552 * t562) * qJD(1) + t418;
t413 = m(3) * t536 + (t587 * qJDD(1) + (t552 * t557 - t553 * t560) * qJD(1)) * t559 + t414;
t412 = m(3) * t526 + t583 * qJDD(1) + (-t548 * t610 + t553 * t562) * qJD(1) + t415;
t409 = mrSges(4,1) * t466 - mrSges(4,2) * t467 + Ifges(4,5) * t524 + Ifges(4,6) * t523 + Ifges(4,3) * t544 + pkin(3) * t573 + pkin(10) * t590 + t564 * t416 + t568 * t419 + t535 * t518 - t534 * t519;
t1 = [Ifges(2,3) * qJDD(1) + mrSges(2,1) * t592 - mrSges(2,2) * t588 + (mrSges(3,1) * t526 - mrSges(3,2) * t527 + pkin(2) * t415 + t561 * t409 + pkin(1) * (t412 * t560 + t417 * t557) + Ifges(3,3) * t562 * qJDD(1) + t578 * t558) * t562 + (t557 * (mrSges(3,2) * t536 - mrSges(3,3) * t526 + t569 * t410 - t565 * t411 - t414 * t619) + t560 * (-mrSges(3,1) * t536 + mrSges(3,3) * t527 - pkin(2) * t414 - t558 * t409) - pkin(1) * t413 + qJ(2) * (-t412 * t557 + t417 * t560) + (-t415 * t620 + t560 * t578) * t561 + ((Ifges(3,2) * t560 ^ 2 + (Ifges(3,1) * t557 + Ifges(3,4) * t621) * t557) * t559 + 0.2e1 * t562 * (Ifges(3,5) * t557 + Ifges(3,6) * t560)) * qJDD(1)) * t559; t413; t409; t572; t622; t446;];
tauJ  = t1;
