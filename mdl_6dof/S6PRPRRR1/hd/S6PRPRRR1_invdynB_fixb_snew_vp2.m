% Calculate vector of inverse dynamics base forces with Newton-Euler for
% S6PRPRRR1
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
%   pkin=[a2,a3,a4,a5,a6,alpha2,d2,d4,d5,d6,theta1,theta3]';
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
% Datum: 2019-05-05 00:09
% Revision: 89c353f7eff3bd693eda4e29f35b2761dbc3ada0 (2019-05-03)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauB = S6PRPRRR1_invdynB_fixb_snew_vp2(qJ, qJD, qJDD, g, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(6,1),zeros(3,1),zeros(12,1),zeros(7,1),zeros(7,3),zeros(7,6)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6PRPRRR1_invdynB_fixb_snew_vp2: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6PRPRRR1_invdynB_fixb_snew_vp2: qJD has to be [6x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [6 1]), ...
  'S6PRPRRR1_invdynB_fixb_snew_vp2: qJDD has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6PRPRRR1_invdynB_fixb_snew_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [12 1]), ...
  'S6PRPRRR1_invdynB_fixb_snew_vp2: pkin has to be [12x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6PRPRRR1_invdynB_fixb_snew_vp2: m has to be [7x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [7,3]), ...
  'S6PRPRRR1_invdynB_fixb_snew_vp2: mrSges has to be [7x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [7 6]), ...
  'S6PRPRRR1_invdynB_fixb_snew_vp2: Ifges has to be [7x6] (double)'); 

%% Symbolic Calculation
% From invdyn_fixb_NewtonEuler_linkframe_tauB_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-05-05 00:05:38
% EndTime: 2019-05-05 00:05:49
% DurationCPUTime: 10.65s
% Computational Cost: add. (193827->299), mult. (364738->387), div. (0->0), fcn. (260299->14), ass. (0->130)
t613 = sin(qJ(5));
t614 = sin(qJ(4));
t617 = cos(qJ(5));
t618 = cos(qJ(4));
t584 = (t613 * t614 - t617 * t618) * qJD(2);
t608 = sin(pkin(6));
t615 = sin(qJ(2));
t641 = t608 * t615;
t619 = cos(qJ(2));
t640 = t608 * t619;
t611 = cos(pkin(6));
t639 = t611 * t615;
t638 = t611 * t619;
t607 = sin(pkin(11));
t610 = cos(pkin(11));
t593 = g(1) * t607 - g(2) * t610;
t594 = -g(1) * t610 - g(2) * t607;
t605 = -g(3) + qJDD(1);
t564 = t593 * t638 - t594 * t615 + t605 * t640;
t559 = qJDD(2) * pkin(2) + t564;
t565 = t593 * t639 + t619 * t594 + t605 * t641;
t620 = qJD(2) ^ 2;
t560 = -pkin(2) * t620 + t565;
t606 = sin(pkin(12));
t609 = cos(pkin(12));
t541 = t606 * t559 + t609 * t560;
t539 = -pkin(3) * t620 + qJDD(2) * pkin(8) + t541;
t577 = -t593 * t608 + t611 * t605;
t576 = qJDD(3) + t577;
t535 = -t614 * t539 + t618 * t576;
t634 = qJD(2) * qJD(4);
t632 = t618 * t634;
t591 = qJDD(2) * t614 + t632;
t532 = (-t591 + t632) * pkin(9) + (t614 * t618 * t620 + qJDD(4)) * pkin(4) + t535;
t536 = t618 * t539 + t614 * t576;
t592 = qJDD(2) * t618 - t614 * t634;
t636 = qJD(2) * t614;
t599 = qJD(4) * pkin(4) - pkin(9) * t636;
t604 = t618 ^ 2;
t533 = -pkin(4) * t604 * t620 + pkin(9) * t592 - qJD(4) * t599 + t536;
t528 = t613 * t532 + t617 * t533;
t585 = (t613 * t618 + t614 * t617) * qJD(2);
t552 = -qJD(5) * t585 - t591 * t613 + t592 * t617;
t567 = mrSges(6,1) * t584 + mrSges(6,2) * t585;
t603 = qJD(4) + qJD(5);
t575 = mrSges(6,1) * t603 - mrSges(6,3) * t585;
t602 = qJDD(4) + qJDD(5);
t568 = pkin(5) * t584 - pkin(10) * t585;
t601 = t603 ^ 2;
t526 = -pkin(5) * t601 + pkin(10) * t602 - t568 * t584 + t528;
t540 = t609 * t559 - t606 * t560;
t625 = -qJDD(2) * pkin(3) - t540;
t534 = -t592 * pkin(4) + t599 * t636 + (-pkin(9) * t604 - pkin(8)) * t620 + t625;
t553 = -qJD(5) * t584 + t591 * t617 + t592 * t613;
t529 = (t584 * t603 - t553) * pkin(10) + (t585 * t603 - t552) * pkin(5) + t534;
t612 = sin(qJ(6));
t616 = cos(qJ(6));
t523 = -t526 * t612 + t529 * t616;
t569 = -t585 * t612 + t603 * t616;
t544 = qJD(6) * t569 + t553 * t616 + t602 * t612;
t570 = t585 * t616 + t603 * t612;
t549 = -mrSges(7,1) * t569 + mrSges(7,2) * t570;
t551 = qJDD(6) - t552;
t578 = qJD(6) + t584;
t554 = -mrSges(7,2) * t578 + mrSges(7,3) * t569;
t521 = m(7) * t523 + mrSges(7,1) * t551 - mrSges(7,3) * t544 - t549 * t570 + t554 * t578;
t524 = t526 * t616 + t529 * t612;
t543 = -qJD(6) * t570 - t553 * t612 + t602 * t616;
t555 = mrSges(7,1) * t578 - mrSges(7,3) * t570;
t522 = m(7) * t524 - mrSges(7,2) * t551 + mrSges(7,3) * t543 + t549 * t569 - t555 * t578;
t627 = -t521 * t612 + t616 * t522;
t512 = m(6) * t528 - mrSges(6,2) * t602 + mrSges(6,3) * t552 - t567 * t584 - t575 * t603 + t627;
t527 = t532 * t617 - t533 * t613;
t574 = -mrSges(6,2) * t603 - mrSges(6,3) * t584;
t525 = -pkin(5) * t602 - pkin(10) * t601 + t568 * t585 - t527;
t623 = -m(7) * t525 + t543 * mrSges(7,1) - mrSges(7,2) * t544 + t569 * t554 - t555 * t570;
t517 = m(6) * t527 + mrSges(6,1) * t602 - mrSges(6,3) * t553 - t567 * t585 + t574 * t603 + t623;
t507 = t613 * t512 + t617 * t517;
t590 = (-mrSges(5,1) * t618 + mrSges(5,2) * t614) * qJD(2);
t635 = qJD(2) * t618;
t596 = -qJD(4) * mrSges(5,2) + mrSges(5,3) * t635;
t505 = m(5) * t535 + qJDD(4) * mrSges(5,1) - mrSges(5,3) * t591 + qJD(4) * t596 - t590 * t636 + t507;
t595 = qJD(4) * mrSges(5,1) - mrSges(5,3) * t636;
t628 = t617 * t512 - t517 * t613;
t506 = m(5) * t536 - qJDD(4) * mrSges(5,2) + mrSges(5,3) * t592 - qJD(4) * t595 + t590 * t635 + t628;
t629 = -t505 * t614 + t618 * t506;
t496 = m(4) * t541 - mrSges(4,1) * t620 - qJDD(2) * mrSges(4,2) + t629;
t538 = -t620 * pkin(8) + t625;
t513 = t616 * t521 + t612 * t522;
t622 = m(6) * t534 - t552 * mrSges(6,1) + mrSges(6,2) * t553 + t584 * t574 + t575 * t585 + t513;
t621 = -m(5) * t538 + t592 * mrSges(5,1) - mrSges(5,2) * t591 - t595 * t636 + t596 * t635 - t622;
t509 = m(4) * t540 + qJDD(2) * mrSges(4,1) - mrSges(4,2) * t620 + t621;
t493 = t606 * t496 + t609 * t509;
t491 = m(3) * t564 + qJDD(2) * mrSges(3,1) - mrSges(3,2) * t620 + t493;
t630 = t609 * t496 - t509 * t606;
t492 = m(3) * t565 - mrSges(3,1) * t620 - qJDD(2) * mrSges(3,2) + t630;
t499 = t618 * t505 + t614 * t506;
t633 = m(4) * t576 + t499;
t498 = m(3) * t577 + t633;
t478 = t491 * t638 + t492 * t639 - t498 * t608;
t476 = m(2) * t593 + t478;
t482 = -t491 * t615 + t619 * t492;
t481 = m(2) * t594 + t482;
t637 = t610 * t476 + t607 * t481;
t477 = t491 * t640 + t492 * t641 + t611 * t498;
t631 = -t476 * t607 + t610 * t481;
t545 = Ifges(7,5) * t570 + Ifges(7,6) * t569 + Ifges(7,3) * t578;
t547 = Ifges(7,1) * t570 + Ifges(7,4) * t569 + Ifges(7,5) * t578;
t514 = -mrSges(7,1) * t525 + mrSges(7,3) * t524 + Ifges(7,4) * t544 + Ifges(7,2) * t543 + Ifges(7,6) * t551 - t545 * t570 + t547 * t578;
t546 = Ifges(7,4) * t570 + Ifges(7,2) * t569 + Ifges(7,6) * t578;
t515 = mrSges(7,2) * t525 - mrSges(7,3) * t523 + Ifges(7,1) * t544 + Ifges(7,4) * t543 + Ifges(7,5) * t551 + t545 * t569 - t546 * t578;
t561 = Ifges(6,5) * t585 - Ifges(6,6) * t584 + Ifges(6,3) * t603;
t562 = Ifges(6,4) * t585 - Ifges(6,2) * t584 + Ifges(6,6) * t603;
t500 = mrSges(6,2) * t534 - mrSges(6,3) * t527 + Ifges(6,1) * t553 + Ifges(6,4) * t552 + Ifges(6,5) * t602 - pkin(10) * t513 - t514 * t612 + t515 * t616 - t561 * t584 - t562 * t603;
t563 = Ifges(6,1) * t585 - Ifges(6,4) * t584 + Ifges(6,5) * t603;
t501 = -mrSges(6,1) * t534 - mrSges(7,1) * t523 + mrSges(7,2) * t524 + mrSges(6,3) * t528 + Ifges(6,4) * t553 - Ifges(7,5) * t544 + Ifges(6,2) * t552 + Ifges(6,6) * t602 - Ifges(7,6) * t543 - Ifges(7,3) * t551 - pkin(5) * t513 - t546 * t570 + t547 * t569 - t561 * t585 + t563 * t603;
t581 = Ifges(5,3) * qJD(4) + (Ifges(5,5) * t614 + Ifges(5,6) * t618) * qJD(2);
t583 = Ifges(5,5) * qJD(4) + (Ifges(5,1) * t614 + Ifges(5,4) * t618) * qJD(2);
t484 = -mrSges(5,1) * t538 + mrSges(5,3) * t536 + Ifges(5,4) * t591 + Ifges(5,2) * t592 + Ifges(5,6) * qJDD(4) - pkin(4) * t622 + pkin(9) * t628 + qJD(4) * t583 + t613 * t500 + t617 * t501 - t581 * t636;
t582 = Ifges(5,6) * qJD(4) + (Ifges(5,4) * t614 + Ifges(5,2) * t618) * qJD(2);
t485 = mrSges(5,2) * t538 - mrSges(5,3) * t535 + Ifges(5,1) * t591 + Ifges(5,4) * t592 + Ifges(5,5) * qJDD(4) - pkin(9) * t507 - qJD(4) * t582 + t500 * t617 - t501 * t613 + t581 * t635;
t474 = mrSges(4,2) * t576 - mrSges(4,3) * t540 + Ifges(4,5) * qJDD(2) - Ifges(4,6) * t620 - pkin(8) * t499 - t484 * t614 + t485 * t618;
t483 = Ifges(4,6) * qJDD(2) - pkin(3) * t499 - mrSges(4,1) * t576 + mrSges(4,3) * t541 - pkin(4) * t507 - Ifges(5,5) * t591 - Ifges(5,6) * t592 - Ifges(5,3) * qJDD(4) - mrSges(5,1) * t535 + mrSges(5,2) * t536 - mrSges(6,1) * t527 + mrSges(6,2) * t528 - t612 * t515 - t616 * t514 - pkin(5) * t623 - pkin(10) * t627 - Ifges(6,5) * t553 - Ifges(6,6) * t552 - Ifges(6,3) * t602 - t585 * t562 - t584 * t563 + t620 * Ifges(4,5) + (-t582 * t614 + t583 * t618) * qJD(2);
t471 = -mrSges(3,1) * t577 + mrSges(3,3) * t565 + t620 * Ifges(3,5) + Ifges(3,6) * qJDD(2) - pkin(2) * t633 + qJ(3) * t630 + t606 * t474 + t609 * t483;
t472 = mrSges(3,2) * t577 - mrSges(3,3) * t564 + Ifges(3,5) * qJDD(2) - Ifges(3,6) * t620 - qJ(3) * t493 + t474 * t609 - t483 * t606;
t624 = pkin(7) * t482 + t471 * t619 + t472 * t615;
t473 = mrSges(3,1) * t564 - mrSges(3,2) * t565 + mrSges(4,1) * t540 - mrSges(4,2) * t541 + t614 * t485 + t618 * t484 + pkin(3) * t621 + pkin(8) * t629 + pkin(2) * t493 + (Ifges(3,3) + Ifges(4,3)) * qJDD(2);
t470 = mrSges(2,2) * t605 - mrSges(2,3) * t593 - t615 * t471 + t619 * t472 + (-t477 * t608 - t478 * t611) * pkin(7);
t469 = -mrSges(2,1) * t605 + mrSges(2,3) * t594 - pkin(1) * t477 - t608 * t473 + t624 * t611;
t1 = [-m(1) * g(1) + t631; -m(1) * g(2) + t637; -m(1) * g(3) + m(2) * t605 + t477; -mrSges(1,2) * g(3) + mrSges(1,3) * g(2) - qJ(1) * t637 - t607 * t469 + t610 * t470; mrSges(1,1) * g(3) - mrSges(1,3) * g(1) + qJ(1) * t631 + t610 * t469 + t607 * t470; -mrSges(1,1) * g(2) + mrSges(2,1) * t593 + mrSges(1,2) * g(1) - mrSges(2,2) * t594 + pkin(1) * t478 + t611 * t473 + t624 * t608;];
tauB  = t1;
