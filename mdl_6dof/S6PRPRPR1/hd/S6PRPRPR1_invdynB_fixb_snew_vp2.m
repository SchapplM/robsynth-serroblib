% Calculate vector of inverse dynamics base forces with Newton-Euler for
% S6PRPRPR1
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
%   pkin=[a2,a3,a4,a5,a6,alpha2,d2,d4,d6,theta1,theta3,theta5]';
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
% tauB [6x1]
%   base forces of inverse dynamics (contains inertial, gravitational coriolis and centrifugal forces)

% Quelle: HybrDyn-Toolbox
% Datum: 2019-05-04 22:12
% Revision: 89c353f7eff3bd693eda4e29f35b2761dbc3ada0 (2019-05-03)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauB = S6PRPRPR1_invdynB_fixb_snew_vp2(qJ, qJD, qJDD, g, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(6,1),zeros(3,1),zeros(12,1),zeros(7,1),zeros(7,3),zeros(7,6)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6PRPRPR1_invdynB_fixb_snew_vp2: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6PRPRPR1_invdynB_fixb_snew_vp2: qJD has to be [6x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [6 1]), ...
  'S6PRPRPR1_invdynB_fixb_snew_vp2: qJDD has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6PRPRPR1_invdynB_fixb_snew_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [12 1]), ...
  'S6PRPRPR1_invdynB_fixb_snew_vp2: pkin has to be [12x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6PRPRPR1_invdynB_fixb_snew_vp2: m has to be [7x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [7,3]), ...
  'S6PRPRPR1_invdynB_fixb_snew_vp2: mrSges has to be [7x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [7 6]), ...
  'S6PRPRPR1_invdynB_fixb_snew_vp2: Ifges has to be [7x6] (double)'); 

%% Symbolic Calculation
% From invdyn_fixb_NewtonEuler_linkframe_tauB_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-05-04 22:08:12
% EndTime: 2019-05-04 22:08:25
% DurationCPUTime: 10.38s
% Computational Cost: add. (176632->299), mult. (352154->388), div. (0->0), fcn. (248931->14), ass. (0->130)
t596 = sin(pkin(12));
t600 = cos(pkin(12));
t605 = sin(qJ(4));
t608 = cos(qJ(4));
t572 = (t596 * t605 - t600 * t608) * qJD(2);
t634 = 2 * qJD(5);
t599 = sin(pkin(6));
t606 = sin(qJ(2));
t633 = t599 * t606;
t609 = cos(qJ(2));
t632 = t599 * t609;
t603 = cos(pkin(6));
t631 = t603 * t606;
t630 = t603 * t609;
t598 = sin(pkin(10));
t602 = cos(pkin(10));
t586 = g(1) * t598 - g(2) * t602;
t587 = -g(1) * t602 - g(2) * t598;
t595 = -g(3) + qJDD(1);
t552 = t586 * t630 - t587 * t606 + t595 * t632;
t547 = qJDD(2) * pkin(2) + t552;
t553 = t586 * t631 + t609 * t587 + t595 * t633;
t611 = qJD(2) ^ 2;
t548 = -pkin(2) * t611 + t553;
t597 = sin(pkin(11));
t601 = cos(pkin(11));
t533 = t597 * t547 + t601 * t548;
t531 = -pkin(3) * t611 + qJDD(2) * pkin(8) + t533;
t569 = -t586 * t599 + t603 * t595;
t566 = qJDD(3) + t569;
t527 = -t605 * t531 + t608 * t566;
t626 = qJD(2) * qJD(4);
t624 = t608 * t626;
t584 = qJDD(2) * t605 + t624;
t524 = (-t584 + t624) * qJ(5) + (t605 * t608 * t611 + qJDD(4)) * pkin(4) + t527;
t528 = t608 * t531 + t605 * t566;
t585 = qJDD(2) * t608 - t605 * t626;
t628 = qJD(2) * t605;
t588 = qJD(4) * pkin(4) - qJ(5) * t628;
t594 = t608 ^ 2;
t525 = -pkin(4) * t594 * t611 + qJ(5) * t585 - qJD(4) * t588 + t528;
t520 = t596 * t524 + t600 * t525 - t572 * t634;
t573 = (t596 * t608 + t600 * t605) * qJD(2);
t555 = mrSges(6,1) * t572 + mrSges(6,2) * t573;
t559 = -t584 * t596 + t585 * t600;
t568 = qJD(4) * mrSges(6,1) - mrSges(6,3) * t573;
t556 = pkin(5) * t572 - pkin(9) * t573;
t610 = qJD(4) ^ 2;
t518 = -pkin(5) * t610 + qJDD(4) * pkin(9) - t556 * t572 + t520;
t532 = t601 * t547 - t597 * t548;
t616 = -qJDD(2) * pkin(3) - t532;
t526 = -t585 * pkin(4) + qJDD(5) + t588 * t628 + (-qJ(5) * t594 - pkin(8)) * t611 + t616;
t560 = t584 * t600 + t585 * t596;
t521 = (qJD(4) * t572 - t560) * pkin(9) + (qJD(4) * t573 - t559) * pkin(5) + t526;
t604 = sin(qJ(6));
t607 = cos(qJ(6));
t515 = -t518 * t604 + t521 * t607;
t563 = qJD(4) * t607 - t573 * t604;
t540 = qJD(6) * t563 + qJDD(4) * t604 + t560 * t607;
t564 = qJD(4) * t604 + t573 * t607;
t541 = -mrSges(7,1) * t563 + mrSges(7,2) * t564;
t571 = qJD(6) + t572;
t542 = -mrSges(7,2) * t571 + mrSges(7,3) * t563;
t558 = qJDD(6) - t559;
t513 = m(7) * t515 + mrSges(7,1) * t558 - mrSges(7,3) * t540 - t541 * t564 + t542 * t571;
t516 = t518 * t607 + t521 * t604;
t539 = -qJD(6) * t564 + qJDD(4) * t607 - t560 * t604;
t543 = mrSges(7,1) * t571 - mrSges(7,3) * t564;
t514 = m(7) * t516 - mrSges(7,2) * t558 + mrSges(7,3) * t539 + t541 * t563 - t543 * t571;
t619 = -t513 * t604 + t607 * t514;
t504 = m(6) * t520 - qJDD(4) * mrSges(6,2) + mrSges(6,3) * t559 - qJD(4) * t568 - t555 * t572 + t619;
t618 = -t600 * t524 + t596 * t525;
t519 = -0.2e1 * qJD(5) * t573 - t618;
t567 = -qJD(4) * mrSges(6,2) - mrSges(6,3) * t572;
t517 = -qJDD(4) * pkin(5) - t610 * pkin(9) + (t634 + t556) * t573 + t618;
t614 = -m(7) * t517 + t539 * mrSges(7,1) - mrSges(7,2) * t540 + t563 * t542 - t543 * t564;
t509 = m(6) * t519 + qJDD(4) * mrSges(6,1) - mrSges(6,3) * t560 + qJD(4) * t567 - t555 * t573 + t614;
t499 = t596 * t504 + t600 * t509;
t583 = (-mrSges(5,1) * t608 + mrSges(5,2) * t605) * qJD(2);
t627 = qJD(2) * t608;
t590 = -qJD(4) * mrSges(5,2) + mrSges(5,3) * t627;
t497 = m(5) * t527 + qJDD(4) * mrSges(5,1) - mrSges(5,3) * t584 + qJD(4) * t590 - t583 * t628 + t499;
t589 = qJD(4) * mrSges(5,1) - mrSges(5,3) * t628;
t620 = t600 * t504 - t509 * t596;
t498 = m(5) * t528 - qJDD(4) * mrSges(5,2) + mrSges(5,3) * t585 - qJD(4) * t589 + t583 * t627 + t620;
t621 = -t497 * t605 + t608 * t498;
t488 = m(4) * t533 - mrSges(4,1) * t611 - qJDD(2) * mrSges(4,2) + t621;
t530 = -t611 * pkin(8) + t616;
t505 = t607 * t513 + t604 * t514;
t613 = m(6) * t526 - t559 * mrSges(6,1) + mrSges(6,2) * t560 + t572 * t567 + t568 * t573 + t505;
t612 = -m(5) * t530 + t585 * mrSges(5,1) - mrSges(5,2) * t584 - t589 * t628 + t590 * t627 - t613;
t501 = m(4) * t532 + qJDD(2) * mrSges(4,1) - mrSges(4,2) * t611 + t612;
t485 = t597 * t488 + t601 * t501;
t483 = m(3) * t552 + qJDD(2) * mrSges(3,1) - mrSges(3,2) * t611 + t485;
t622 = t601 * t488 - t501 * t597;
t484 = m(3) * t553 - mrSges(3,1) * t611 - qJDD(2) * mrSges(3,2) + t622;
t491 = t608 * t497 + t605 * t498;
t625 = m(4) * t566 + t491;
t490 = m(3) * t569 + t625;
t470 = t483 * t630 + t484 * t631 - t490 * t599;
t468 = m(2) * t586 + t470;
t474 = -t483 * t606 + t609 * t484;
t473 = m(2) * t587 + t474;
t629 = t602 * t468 + t598 * t473;
t469 = t483 * t632 + t484 * t633 + t603 * t490;
t623 = -t468 * t598 + t602 * t473;
t534 = Ifges(7,5) * t564 + Ifges(7,6) * t563 + Ifges(7,3) * t571;
t536 = Ifges(7,1) * t564 + Ifges(7,4) * t563 + Ifges(7,5) * t571;
t506 = -mrSges(7,1) * t517 + mrSges(7,3) * t516 + Ifges(7,4) * t540 + Ifges(7,2) * t539 + Ifges(7,6) * t558 - t534 * t564 + t536 * t571;
t535 = Ifges(7,4) * t564 + Ifges(7,2) * t563 + Ifges(7,6) * t571;
t507 = mrSges(7,2) * t517 - mrSges(7,3) * t515 + Ifges(7,1) * t540 + Ifges(7,4) * t539 + Ifges(7,5) * t558 + t534 * t563 - t535 * t571;
t549 = Ifges(6,5) * t573 - Ifges(6,6) * t572 + Ifges(6,3) * qJD(4);
t550 = Ifges(6,4) * t573 - Ifges(6,2) * t572 + Ifges(6,6) * qJD(4);
t492 = mrSges(6,2) * t526 - mrSges(6,3) * t519 + Ifges(6,1) * t560 + Ifges(6,4) * t559 + Ifges(6,5) * qJDD(4) - pkin(9) * t505 - qJD(4) * t550 - t506 * t604 + t507 * t607 - t549 * t572;
t551 = Ifges(6,1) * t573 - Ifges(6,4) * t572 + Ifges(6,5) * qJD(4);
t493 = -mrSges(6,1) * t526 - mrSges(7,1) * t515 + mrSges(7,2) * t516 + mrSges(6,3) * t520 + Ifges(6,4) * t560 - Ifges(7,5) * t540 + Ifges(6,2) * t559 + Ifges(6,6) * qJDD(4) - Ifges(7,6) * t539 - Ifges(7,3) * t558 - pkin(5) * t505 + qJD(4) * t551 - t535 * t564 + t536 * t563 - t549 * t573;
t576 = Ifges(5,3) * qJD(4) + (Ifges(5,5) * t605 + Ifges(5,6) * t608) * qJD(2);
t578 = Ifges(5,5) * qJD(4) + (Ifges(5,1) * t605 + Ifges(5,4) * t608) * qJD(2);
t476 = -mrSges(5,1) * t530 + mrSges(5,3) * t528 + Ifges(5,4) * t584 + Ifges(5,2) * t585 + Ifges(5,6) * qJDD(4) - pkin(4) * t613 + qJ(5) * t620 + qJD(4) * t578 + t596 * t492 + t600 * t493 - t576 * t628;
t577 = Ifges(5,6) * qJD(4) + (Ifges(5,4) * t605 + Ifges(5,2) * t608) * qJD(2);
t477 = mrSges(5,2) * t530 - mrSges(5,3) * t527 + Ifges(5,1) * t584 + Ifges(5,4) * t585 + Ifges(5,5) * qJDD(4) - qJ(5) * t499 - qJD(4) * t577 + t492 * t600 - t493 * t596 + t576 * t627;
t466 = mrSges(4,2) * t566 - mrSges(4,3) * t532 + Ifges(4,5) * qJDD(2) - Ifges(4,6) * t611 - pkin(8) * t491 - t476 * t605 + t477 * t608;
t475 = -pkin(3) * t491 + mrSges(4,3) * t533 - mrSges(4,1) * t566 + Ifges(4,6) * qJDD(2) - pkin(4) * t499 - Ifges(5,5) * t584 - Ifges(5,6) * t585 - mrSges(5,1) * t527 + mrSges(5,2) * t528 - Ifges(6,5) * t560 - Ifges(6,6) * t559 - mrSges(6,1) * t519 + mrSges(6,2) * t520 - t604 * t507 - t607 * t506 - pkin(5) * t614 - pkin(9) * t619 - t573 * t550 - t572 * t551 + t611 * Ifges(4,5) + (-Ifges(5,3) - Ifges(6,3)) * qJDD(4) + (-t577 * t605 + t578 * t608) * qJD(2);
t463 = -mrSges(3,1) * t569 + mrSges(3,3) * t553 + t611 * Ifges(3,5) + Ifges(3,6) * qJDD(2) - pkin(2) * t625 + qJ(3) * t622 + t597 * t466 + t601 * t475;
t464 = mrSges(3,2) * t569 - mrSges(3,3) * t552 + Ifges(3,5) * qJDD(2) - Ifges(3,6) * t611 - qJ(3) * t485 + t466 * t601 - t475 * t597;
t615 = pkin(7) * t474 + t463 * t609 + t464 * t606;
t465 = mrSges(3,1) * t552 - mrSges(3,2) * t553 + mrSges(4,1) * t532 - mrSges(4,2) * t533 + t605 * t477 + t608 * t476 + pkin(3) * t612 + pkin(8) * t621 + pkin(2) * t485 + (Ifges(3,3) + Ifges(4,3)) * qJDD(2);
t462 = mrSges(2,2) * t595 - mrSges(2,3) * t586 - t606 * t463 + t609 * t464 + (-t469 * t599 - t470 * t603) * pkin(7);
t461 = -mrSges(2,1) * t595 + mrSges(2,3) * t587 - pkin(1) * t469 - t599 * t465 + t615 * t603;
t1 = [-m(1) * g(1) + t623; -m(1) * g(2) + t629; -m(1) * g(3) + m(2) * t595 + t469; -mrSges(1,2) * g(3) + mrSges(1,3) * g(2) - qJ(1) * t629 - t598 * t461 + t602 * t462; mrSges(1,1) * g(3) - mrSges(1,3) * g(1) + qJ(1) * t623 + t602 * t461 + t598 * t462; -mrSges(1,1) * g(2) + mrSges(2,1) * t586 + mrSges(1,2) * g(1) - mrSges(2,2) * t587 + pkin(1) * t470 + t603 * t465 + t615 * t599;];
tauB  = t1;
