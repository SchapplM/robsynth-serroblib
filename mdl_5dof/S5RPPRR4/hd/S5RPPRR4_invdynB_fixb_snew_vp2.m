% Calculate vector of inverse dynamics base forces with Newton-Euler for
% S5RPPRR4
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% qJD [5x1]
%   Generalized joint velocities
% qJDD [5x1]
%   Generalized joint accelerations
% g [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [9x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d1,d4,d5,theta2,theta3]';
% m [6x1]
%   mass of all robot links (including the base)
% mrSges [6x3]
%  first moment of all robot links (mass times center of mass in body frames)
%  rows: links of the robot (starting with base)
%  columns: x-, y-, z-coordinates
% Ifges [6x6]
%   inertia of all robot links about their respective body frame origins, in body frames
%   rows: links of the robot (starting with base)
%   columns: xx, yy, zz, xy, xz, yz (see inertial_parameters_convert_par1_par2.m)
% 
% Output:
% tauB [6x1]
%   base forces of inverse dynamics (contains inertial, gravitational coriolis and centrifugal forces)

% Quelle: HybrDyn-Toolbox
% Datum: 2022-01-23 09:17
% Revision: fd3771346c4aea32fdeb66112c511235427c26a7 (2022-01-20)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauB = S5RPPRR4_invdynB_fixb_snew_vp2(qJ, qJD, qJDD, g, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(5,1),zeros(3,1),zeros(9,1),zeros(6,1),zeros(6,3),zeros(6,6)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RPPRR4_invdynB_fixb_snew_vp2: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RPPRR4_invdynB_fixb_snew_vp2: qJD has to be [5x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [5 1]), ...
  'S5RPPRR4_invdynB_fixb_snew_vp2: qJDD has to be [5x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S5RPPRR4_invdynB_fixb_snew_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S5RPPRR4_invdynB_fixb_snew_vp2: pkin has to be [9x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5RPPRR4_invdynB_fixb_snew_vp2: m has to be [6x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [6,3]), ...
  'S5RPPRR4_invdynB_fixb_snew_vp2: mrSges has to be [6x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [6 6]), ...
  'S5RPPRR4_invdynB_fixb_snew_vp2: Ifges has to be [6x6] (double)'); 

%% Symbolic Calculation
% From invdyn_fixb_NewtonEuler_linkframe_tauB_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2022-01-23 09:16:37
% EndTime: 2022-01-23 09:16:43
% DurationCPUTime: 6.23s
% Computational Cost: add. (59829->279), mult. (163866->381), div. (0->0), fcn. (112502->10), ass. (0->128)
t600 = sin(qJ(1));
t603 = cos(qJ(1));
t580 = -t603 * g(1) - t600 * g(2);
t604 = qJD(1) ^ 2;
t644 = -t604 * pkin(1) + qJDD(1) * qJ(2) + (2 * qJD(1) * qJD(2)) + t580;
t579 = t600 * g(1) - t603 * g(2);
t610 = -t604 * qJ(2) + qJDD(2) - t579;
t595 = sin(pkin(8));
t597 = cos(pkin(8));
t617 = -pkin(2) * t597 - qJ(3) * t595;
t634 = t595 * qJD(1);
t643 = (-pkin(1) + t617) * qJDD(1) + t610 - 0.2e1 * qJD(3) * t634;
t554 = -t597 * g(3) - t644 * t595;
t642 = mrSges(3,2) * t595;
t641 = Ifges(3,6) * t597;
t593 = t595 ^ 2;
t640 = t593 * t604;
t594 = sin(pkin(9));
t639 = t594 * t595;
t596 = cos(pkin(9));
t638 = t595 * t596;
t555 = -t595 * g(3) + t644 * t597;
t574 = (-mrSges(3,1) * t597 + t642) * qJD(1);
t573 = t617 * qJD(1);
t633 = t597 * qJD(1);
t544 = t573 * t633 + t555;
t615 = -pkin(3) * t597 - pkin(6) * t638;
t636 = t643 * t596;
t523 = t615 * qJDD(1) + (-t544 + (-pkin(3) * t593 * t596 + pkin(6) * t595 * t597) * t604) * t594 + t636;
t532 = t596 * t544 + t643 * t594;
t572 = t615 * qJD(1);
t631 = qJDD(1) * t595;
t627 = t594 * t631;
t629 = t594 ^ 2 * t640;
t524 = -pkin(3) * t629 - pkin(6) * t627 + t572 * t633 + t532;
t599 = sin(qJ(4));
t602 = cos(qJ(4));
t512 = t602 * t523 - t599 * t524;
t612 = (-t594 * t602 - t596 * t599) * t595;
t562 = qJD(1) * t612;
t611 = (-t594 * t599 + t596 * t602) * t595;
t548 = t562 * qJD(4) + qJDD(1) * t611;
t563 = qJD(1) * t611;
t630 = t597 * qJDD(1);
t583 = qJDD(4) - t630;
t584 = qJD(4) - t633;
t510 = (t562 * t584 - t548) * pkin(7) + (t562 * t563 + t583) * pkin(4) + t512;
t513 = t599 * t523 + t602 * t524;
t547 = -t563 * qJD(4) + qJDD(1) * t612;
t553 = t584 * pkin(4) - t563 * pkin(7);
t561 = t562 ^ 2;
t511 = -t561 * pkin(4) + t547 * pkin(7) - t584 * t553 + t513;
t598 = sin(qJ(5));
t601 = cos(qJ(5));
t508 = t601 * t510 - t598 * t511;
t541 = t601 * t562 - t598 * t563;
t519 = t541 * qJD(5) + t598 * t547 + t601 * t548;
t542 = t598 * t562 + t601 * t563;
t530 = -t541 * mrSges(6,1) + t542 * mrSges(6,2);
t582 = qJD(5) + t584;
t534 = -t582 * mrSges(6,2) + t541 * mrSges(6,3);
t578 = qJDD(5) + t583;
t506 = m(6) * t508 + t578 * mrSges(6,1) - t519 * mrSges(6,3) - t542 * t530 + t582 * t534;
t509 = t598 * t510 + t601 * t511;
t518 = -t542 * qJD(5) + t601 * t547 - t598 * t548;
t535 = t582 * mrSges(6,1) - t542 * mrSges(6,3);
t507 = m(6) * t509 - t578 * mrSges(6,2) + t518 * mrSges(6,3) + t541 * t530 - t582 * t535;
t498 = t601 * t506 + t598 * t507;
t545 = -t562 * mrSges(5,1) + t563 * mrSges(5,2);
t549 = -t584 * mrSges(5,2) + t562 * mrSges(5,3);
t496 = m(5) * t512 + t583 * mrSges(5,1) - t548 * mrSges(5,3) - t563 * t545 + t584 * t549 + t498;
t550 = t584 * mrSges(5,1) - t563 * mrSges(5,3);
t622 = -t598 * t506 + t601 * t507;
t497 = m(5) * t513 - t583 * mrSges(5,2) + t547 * mrSges(5,3) + t562 * t545 - t584 * t550 + t622;
t492 = t602 * t496 + t599 * t497;
t531 = -t594 * t544 + t636;
t620 = mrSges(4,1) * t594 + mrSges(4,2) * t596;
t566 = t620 * t634;
t613 = mrSges(4,2) * t597 - mrSges(4,3) * t639;
t569 = t613 * qJD(1);
t614 = -mrSges(4,1) * t597 - mrSges(4,3) * t638;
t490 = m(4) * t531 + t614 * qJDD(1) + (-t566 * t638 - t569 * t597) * qJD(1) + t492;
t570 = t614 * qJD(1);
t623 = -t599 * t496 + t602 * t497;
t491 = m(4) * t532 + t613 * qJDD(1) + (-t566 * t639 + t570 * t597) * qJD(1) + t623;
t624 = -t594 * t490 + t596 * t491;
t485 = m(3) * t555 + (qJDD(1) * mrSges(3,3) + qJD(1) * t574) * t597 + t624;
t543 = t573 * t634 + qJDD(3) - t554;
t533 = t596 * t572 * t634 + pkin(3) * t627 - pkin(6) * t629 + t543;
t515 = -t547 * pkin(4) - t561 * pkin(7) + t563 * t553 + t533;
t616 = m(6) * t515 - t518 * mrSges(6,1) + t519 * mrSges(6,2) - t541 * t534 + t542 * t535;
t606 = m(5) * t533 - t547 * mrSges(5,1) + t548 * mrSges(5,2) - t562 * t549 + t563 * t550 + t616;
t605 = -m(4) * t543 - t606;
t502 = t605 + m(3) * t554 + ((-mrSges(3,3) - t620) * qJDD(1) + (-t569 * t594 - t570 * t596 - t574) * qJD(1)) * t595;
t625 = t597 * t485 - t595 * t502;
t479 = m(2) * t580 - t604 * mrSges(2,1) - qJDD(1) * mrSges(2,2) + t625;
t486 = t596 * t490 + t594 * t491;
t568 = -qJDD(1) * pkin(1) + t610;
t607 = -m(3) * t568 + mrSges(3,1) * t630 - t486 + (t597 ^ 2 * t604 + t640) * mrSges(3,3);
t482 = m(2) * t579 - t604 * mrSges(2,2) + (mrSges(2,1) - t642) * qJDD(1) + t607;
t637 = t600 * t479 + t603 * t482;
t480 = t595 * t485 + t597 * t502;
t626 = t603 * t479 - t600 * t482;
t619 = Ifges(3,1) * t595 + Ifges(3,4) * t597;
t618 = Ifges(4,5) * t596 - Ifges(4,6) * t594;
t609 = -Ifges(4,5) * t597 + (Ifges(4,1) * t596 - Ifges(4,4) * t594) * t595;
t608 = -Ifges(4,6) * t597 + (Ifges(4,4) * t596 - Ifges(4,2) * t594) * t595;
t575 = (Ifges(3,5) * t595 + t641) * qJD(1);
t559 = t609 * qJD(1);
t558 = t608 * qJD(1);
t557 = (-Ifges(4,3) * t597 + t618 * t595) * qJD(1);
t538 = Ifges(5,1) * t563 + Ifges(5,4) * t562 + Ifges(5,5) * t584;
t537 = Ifges(5,4) * t563 + Ifges(5,2) * t562 + Ifges(5,6) * t584;
t536 = Ifges(5,5) * t563 + Ifges(5,6) * t562 + Ifges(5,3) * t584;
t527 = Ifges(6,1) * t542 + Ifges(6,4) * t541 + Ifges(6,5) * t582;
t526 = Ifges(6,4) * t542 + Ifges(6,2) * t541 + Ifges(6,6) * t582;
t525 = Ifges(6,5) * t542 + Ifges(6,6) * t541 + Ifges(6,3) * t582;
t500 = mrSges(6,2) * t515 - mrSges(6,3) * t508 + Ifges(6,1) * t519 + Ifges(6,4) * t518 + Ifges(6,5) * t578 + t541 * t525 - t582 * t526;
t499 = -mrSges(6,1) * t515 + mrSges(6,3) * t509 + Ifges(6,4) * t519 + Ifges(6,2) * t518 + Ifges(6,6) * t578 - t542 * t525 + t582 * t527;
t488 = mrSges(5,2) * t533 - mrSges(5,3) * t512 + Ifges(5,1) * t548 + Ifges(5,4) * t547 + Ifges(5,5) * t583 - pkin(7) * t498 - t598 * t499 + t601 * t500 + t562 * t536 - t584 * t537;
t487 = -mrSges(5,1) * t533 + mrSges(5,3) * t513 + Ifges(5,4) * t548 + Ifges(5,2) * t547 + Ifges(5,6) * t583 - pkin(4) * t616 + pkin(7) * t622 + t601 * t499 + t598 * t500 - t563 * t536 + t584 * t538;
t476 = mrSges(4,2) * t543 - mrSges(4,3) * t531 - pkin(6) * t492 - t599 * t487 + t602 * t488 + (-t557 * t639 + t558 * t597) * qJD(1) + t609 * qJDD(1);
t475 = -mrSges(4,1) * t543 + mrSges(4,3) * t532 + t599 * t488 + t602 * t487 - pkin(3) * t606 + pkin(6) * t623 + (-t557 * t638 - t597 * t559) * qJD(1) + t608 * qJDD(1);
t474 = -Ifges(5,3) * t583 - Ifges(6,3) * t578 + t562 * t538 - t563 * t537 - mrSges(3,1) * t568 - Ifges(5,6) * t547 - Ifges(5,5) * t548 + mrSges(3,3) * t555 + t541 * t527 - t542 * t526 - mrSges(4,1) * t531 + mrSges(4,2) * t532 - Ifges(6,6) * t518 - Ifges(6,5) * t519 - mrSges(5,1) * t512 + mrSges(5,2) * t513 - mrSges(6,1) * t508 + mrSges(6,2) * t509 - pkin(4) * t498 - pkin(3) * t492 - pkin(2) * t486 + (Ifges(3,2) + Ifges(4,3)) * t630 + ((Ifges(3,4) - t618) * qJDD(1) + (-t558 * t596 - t559 * t594 - t575) * qJD(1)) * t595;
t473 = mrSges(3,2) * t568 - mrSges(3,3) * t554 - qJ(3) * t486 + t619 * qJDD(1) - t594 * t475 + t596 * t476 + t575 * t633;
t472 = t604 * Ifges(2,5) + mrSges(2,1) * g(3) + mrSges(2,3) * t580 - mrSges(3,1) * t554 + mrSges(3,2) * t555 - t594 * t476 - t596 * t475 - pkin(2) * t605 - qJ(3) * t624 - pkin(1) * t480 + (-t641 + Ifges(2,6) + (pkin(2) * t620 - Ifges(3,5)) * t595) * qJDD(1) + (-pkin(2) * (-t569 * t639 - t570 * t638) + (-t595 * (Ifges(3,4) * t595 + Ifges(3,2) * t597) + t597 * t619) * qJD(1)) * qJD(1);
t471 = -mrSges(2,2) * g(3) - mrSges(2,3) * t579 + Ifges(2,5) * qJDD(1) - t604 * Ifges(2,6) - qJ(2) * t480 + t597 * t473 - t595 * t474;
t1 = [-m(1) * g(1) + t626; -m(1) * g(2) + t637; (-m(1) - m(2)) * g(3) + t480; -mrSges(1,2) * g(3) + mrSges(1,3) * g(2) - pkin(5) * t637 + t603 * t471 - t600 * t472; mrSges(1,1) * g(3) - mrSges(1,3) * g(1) + pkin(5) * t626 + t600 * t471 + t603 * t472; -mrSges(1,1) * g(2) + mrSges(1,2) * g(1) + Ifges(2,3) * qJDD(1) + mrSges(2,1) * t579 - mrSges(2,2) * t580 + t595 * t473 + t597 * t474 + pkin(1) * (-mrSges(3,2) * t631 + t607) + qJ(2) * t625;];
tauB = t1;
