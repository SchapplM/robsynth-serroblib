% Calculate vector of inverse dynamics base forces with Newton-Euler for
% S6PRPPRR1
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
%   pkin=[a2,a3,a4,a5,a6,alpha2,d2,d5,d6,theta1,theta3,theta4]';
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
% Datum: 2019-05-04 21:43
% Revision: 89c353f7eff3bd693eda4e29f35b2761dbc3ada0 (2019-05-03)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauB = S6PRPPRR1_invdynB_fixb_snew_vp2(qJ, qJD, qJDD, g, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(6,1),zeros(3,1),zeros(12,1),zeros(7,1),zeros(7,3),zeros(7,6)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6PRPPRR1_invdynB_fixb_snew_vp2: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6PRPPRR1_invdynB_fixb_snew_vp2: qJD has to be [6x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [6 1]), ...
  'S6PRPPRR1_invdynB_fixb_snew_vp2: qJDD has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6PRPPRR1_invdynB_fixb_snew_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [12 1]), ...
  'S6PRPPRR1_invdynB_fixb_snew_vp2: pkin has to be [12x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6PRPPRR1_invdynB_fixb_snew_vp2: m has to be [7x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [7,3]), ...
  'S6PRPPRR1_invdynB_fixb_snew_vp2: mrSges has to be [7x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [7 6]), ...
  'S6PRPPRR1_invdynB_fixb_snew_vp2: Ifges has to be [7x6] (double)'); 

%% Symbolic Calculation
% From invdyn_fixb_NewtonEuler_linkframe_tauB_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-05-04 21:40:16
% EndTime: 2019-05-04 21:40:26
% DurationCPUTime: 9.66s
% Computational Cost: add. (164134->277), mult. (333003->354), div. (0->0), fcn. (245761->14), ass. (0->131)
t612 = qJD(2) ^ 2;
t597 = sin(pkin(12));
t601 = cos(pkin(12));
t606 = sin(qJ(5));
t609 = cos(qJ(5));
t618 = t597 * t606 - t601 * t609;
t577 = t618 * qJD(2);
t619 = t597 * t609 + t601 * t606;
t578 = t619 * qJD(2);
t632 = t578 * qJD(5);
t565 = -qJDD(2) * t618 - t632;
t643 = pkin(4) * t601;
t642 = mrSges(5,2) * t597;
t595 = t601 ^ 2;
t641 = t595 * t612;
t600 = sin(pkin(6));
t607 = sin(qJ(2));
t640 = t600 * t607;
t610 = cos(qJ(2));
t639 = t600 * t610;
t604 = cos(pkin(6));
t638 = t604 * t607;
t637 = t604 * t610;
t599 = sin(pkin(10));
t603 = cos(pkin(10));
t584 = g(1) * t599 - g(2) * t603;
t585 = -g(1) * t603 - g(2) * t599;
t596 = -g(3) + qJDD(1);
t558 = t584 * t637 - t585 * t607 + t596 * t639;
t553 = qJDD(2) * pkin(2) + t558;
t559 = t584 * t638 + t585 * t610 + t596 * t640;
t554 = -pkin(2) * t612 + t559;
t598 = sin(pkin(11));
t602 = cos(pkin(11));
t539 = t553 * t598 + t554 * t602;
t537 = -pkin(3) * t612 + qJDD(2) * qJ(4) + t539;
t575 = -t584 * t600 + t596 * t604;
t572 = qJDD(3) + t575;
t631 = qJD(2) * qJD(4);
t635 = t572 * t601 - 0.2e1 * t597 * t631;
t530 = (-pkin(8) * qJDD(2) + t612 * t643 - t537) * t597 + t635;
t533 = t597 * t572 + (t537 + 0.2e1 * t631) * t601;
t630 = qJDD(2) * t601;
t531 = -pkin(4) * t641 + pkin(8) * t630 + t533;
t526 = t530 * t606 + t531 * t609;
t561 = mrSges(6,1) * t577 + mrSges(6,2) * t578;
t574 = qJD(5) * mrSges(6,1) - mrSges(6,3) * t578;
t564 = pkin(5) * t577 - pkin(9) * t578;
t611 = qJD(5) ^ 2;
t524 = -pkin(5) * t611 + qJDD(5) * pkin(9) - t564 * t577 + t526;
t594 = t597 ^ 2;
t538 = t602 * t553 - t554 * t598;
t620 = qJDD(4) - t538;
t534 = (-pkin(3) - t643) * qJDD(2) + (-qJ(4) + (-t594 - t595) * pkin(8)) * t612 + t620;
t633 = t577 * qJD(5);
t566 = qJDD(2) * t619 - t633;
t527 = (-t566 + t633) * pkin(9) + (-t565 + t632) * pkin(5) + t534;
t605 = sin(qJ(6));
t608 = cos(qJ(6));
t521 = -t524 * t605 + t527 * t608;
t569 = qJD(5) * t608 - t578 * t605;
t546 = qJD(6) * t569 + qJDD(5) * t605 + t566 * t608;
t570 = qJD(5) * t605 + t578 * t608;
t547 = -mrSges(7,1) * t569 + mrSges(7,2) * t570;
t576 = qJD(6) + t577;
t548 = -mrSges(7,2) * t576 + mrSges(7,3) * t569;
t563 = qJDD(6) - t565;
t519 = m(7) * t521 + mrSges(7,1) * t563 - mrSges(7,3) * t546 - t547 * t570 + t548 * t576;
t522 = t524 * t608 + t527 * t605;
t545 = -qJD(6) * t570 + qJDD(5) * t608 - t566 * t605;
t549 = mrSges(7,1) * t576 - mrSges(7,3) * t570;
t520 = m(7) * t522 - mrSges(7,2) * t563 + mrSges(7,3) * t545 + t547 * t569 - t549 * t576;
t624 = -t519 * t605 + t520 * t608;
t510 = m(6) * t526 - qJDD(5) * mrSges(6,2) + mrSges(6,3) * t565 - qJD(5) * t574 - t561 * t577 + t624;
t525 = t530 * t609 - t531 * t606;
t573 = -qJD(5) * mrSges(6,2) - mrSges(6,3) * t577;
t523 = -qJDD(5) * pkin(5) - pkin(9) * t611 + t564 * t578 - t525;
t615 = -m(7) * t523 + mrSges(7,1) * t545 - mrSges(7,2) * t546 + t548 * t569 - t549 * t570;
t515 = m(6) * t525 + qJDD(5) * mrSges(6,1) - mrSges(6,3) * t566 + qJD(5) * t573 - t561 * t578 + t615;
t505 = t510 * t606 + t515 * t609;
t532 = -t537 * t597 + t635;
t617 = mrSges(5,3) * qJDD(2) + t612 * (-mrSges(5,1) * t601 + t642);
t503 = m(5) * t532 - t597 * t617 + t505;
t625 = t510 * t609 - t606 * t515;
t504 = m(5) * t533 + t601 * t617 + t625;
t626 = -t503 * t597 + t504 * t601;
t494 = m(4) * t539 - mrSges(4,1) * t612 - qJDD(2) * mrSges(4,2) + t626;
t536 = -qJDD(2) * pkin(3) - t612 * qJ(4) + t620;
t511 = t519 * t608 + t520 * t605;
t614 = m(6) * t534 - mrSges(6,1) * t565 + t566 * mrSges(6,2) + t573 * t577 + t578 * t574 + t511;
t613 = -m(5) * t536 + mrSges(5,1) * t630 - t614 + (t594 * t612 + t641) * mrSges(5,3);
t507 = t613 + (mrSges(4,1) - t642) * qJDD(2) + m(4) * t538 - t612 * mrSges(4,2);
t491 = t494 * t598 + t507 * t602;
t489 = m(3) * t558 + qJDD(2) * mrSges(3,1) - mrSges(3,2) * t612 + t491;
t627 = t494 * t602 - t507 * t598;
t490 = m(3) * t559 - mrSges(3,1) * t612 - qJDD(2) * mrSges(3,2) + t627;
t497 = t503 * t601 + t504 * t597;
t629 = m(4) * t572 + t497;
t496 = m(3) * t575 + t629;
t476 = t489 * t637 + t490 * t638 - t496 * t600;
t474 = m(2) * t584 + t476;
t480 = -t489 * t607 + t490 * t610;
t479 = m(2) * t585 + t480;
t636 = t474 * t603 + t479 * t599;
t621 = Ifges(5,5) * t597 + Ifges(5,6) * t601;
t634 = t612 * t621;
t475 = t489 * t639 + t490 * t640 + t496 * t604;
t628 = -t474 * t599 + t479 * t603;
t623 = Ifges(5,1) * t597 + Ifges(5,4) * t601;
t622 = Ifges(5,4) * t597 + Ifges(5,2) * t601;
t540 = Ifges(7,5) * t570 + Ifges(7,6) * t569 + Ifges(7,3) * t576;
t542 = Ifges(7,1) * t570 + Ifges(7,4) * t569 + Ifges(7,5) * t576;
t512 = -mrSges(7,1) * t523 + mrSges(7,3) * t522 + Ifges(7,4) * t546 + Ifges(7,2) * t545 + Ifges(7,6) * t563 - t540 * t570 + t542 * t576;
t541 = Ifges(7,4) * t570 + Ifges(7,2) * t569 + Ifges(7,6) * t576;
t513 = mrSges(7,2) * t523 - mrSges(7,3) * t521 + Ifges(7,1) * t546 + Ifges(7,4) * t545 + Ifges(7,5) * t563 + t540 * t569 - t541 * t576;
t555 = Ifges(6,5) * t578 - Ifges(6,6) * t577 + Ifges(6,3) * qJD(5);
t556 = Ifges(6,4) * t578 - Ifges(6,2) * t577 + Ifges(6,6) * qJD(5);
t498 = mrSges(6,2) * t534 - mrSges(6,3) * t525 + Ifges(6,1) * t566 + Ifges(6,4) * t565 + Ifges(6,5) * qJDD(5) - pkin(9) * t511 - qJD(5) * t556 - t512 * t605 + t513 * t608 - t555 * t577;
t557 = Ifges(6,1) * t578 - Ifges(6,4) * t577 + Ifges(6,5) * qJD(5);
t499 = -mrSges(6,1) * t534 - mrSges(7,1) * t521 + mrSges(7,2) * t522 + mrSges(6,3) * t526 + Ifges(6,4) * t566 - Ifges(7,5) * t546 + Ifges(6,2) * t565 + Ifges(6,6) * qJDD(5) - Ifges(7,6) * t545 - Ifges(7,3) * t563 - pkin(5) * t511 + qJD(5) * t557 - t541 * t570 + t542 * t569 - t555 * t578;
t482 = -mrSges(5,1) * t536 + mrSges(5,3) * t533 - pkin(4) * t614 + pkin(8) * t625 + qJDD(2) * t622 + t606 * t498 + t609 * t499 - t597 * t634;
t483 = mrSges(5,2) * t536 - mrSges(5,3) * t532 - pkin(8) * t505 + qJDD(2) * t623 + t609 * t498 - t606 * t499 + t601 * t634;
t472 = mrSges(4,2) * t572 - mrSges(4,3) * t538 + Ifges(4,5) * qJDD(2) - Ifges(4,6) * t612 - qJ(4) * t497 - t482 * t597 + t483 * t601;
t481 = -pkin(3) * t497 + mrSges(4,3) * t539 - mrSges(4,1) * t572 - pkin(4) * t505 + mrSges(5,2) * t533 - mrSges(5,1) * t532 - t605 * t513 - t608 * t512 - pkin(5) * t615 - pkin(9) * t624 - Ifges(6,3) * qJDD(5) - t578 * t556 - t577 * t557 - mrSges(6,1) * t525 + mrSges(6,2) * t526 - Ifges(6,5) * t566 - Ifges(6,6) * t565 + (Ifges(4,6) - t621) * qJDD(2) + (-t597 * t622 + t601 * t623 + Ifges(4,5)) * t612;
t469 = -mrSges(3,1) * t575 + mrSges(3,3) * t559 + t612 * Ifges(3,5) + Ifges(3,6) * qJDD(2) - pkin(2) * t629 + qJ(3) * t627 + t598 * t472 + t602 * t481;
t470 = mrSges(3,2) * t575 - mrSges(3,3) * t558 + Ifges(3,5) * qJDD(2) - Ifges(3,6) * t612 - qJ(3) * t491 + t472 * t602 - t481 * t598;
t616 = pkin(7) * t480 + t469 * t610 + t470 * t607;
t471 = mrSges(3,1) * t558 - mrSges(3,2) * t559 + mrSges(4,1) * t538 - mrSges(4,2) * t539 + t597 * t483 + t601 * t482 + pkin(3) * t613 + qJ(4) * t626 + pkin(2) * t491 + (-pkin(3) * t642 + Ifges(3,3) + Ifges(4,3)) * qJDD(2);
t468 = mrSges(2,2) * t596 - mrSges(2,3) * t584 - t607 * t469 + t610 * t470 + (-t475 * t600 - t476 * t604) * pkin(7);
t467 = -mrSges(2,1) * t596 + mrSges(2,3) * t585 - pkin(1) * t475 - t600 * t471 + t604 * t616;
t1 = [-m(1) * g(1) + t628; -m(1) * g(2) + t636; -m(1) * g(3) + m(2) * t596 + t475; -mrSges(1,2) * g(3) + mrSges(1,3) * g(2) - qJ(1) * t636 - t467 * t599 + t468 * t603; mrSges(1,1) * g(3) - mrSges(1,3) * g(1) + qJ(1) * t628 + t603 * t467 + t599 * t468; -mrSges(1,1) * g(2) + mrSges(2,1) * t584 + mrSges(1,2) * g(1) - mrSges(2,2) * t585 + pkin(1) * t476 + t604 * t471 + t600 * t616;];
tauB  = t1;
