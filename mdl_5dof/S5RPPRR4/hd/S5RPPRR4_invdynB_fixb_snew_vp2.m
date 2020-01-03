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
% m_mdh [6x1]
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
% Datum: 2020-01-03 11:32
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
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
% StartTime: 2020-01-03 11:31:04
% EndTime: 2020-01-03 11:31:13
% DurationCPUTime: 6.66s
% Computational Cost: add. (59829->279), mult. (163866->381), div. (0->0), fcn. (112502->10), ass. (0->128)
t608 = sin(qJ(1));
t611 = cos(qJ(1));
t585 = -t608 * g(2) + t611 * g(3);
t612 = qJD(1) ^ 2;
t652 = -t612 * pkin(1) + qJDD(1) * qJ(2) + (2 * qJD(1) * qJD(2)) + t585;
t586 = -t611 * g(2) - t608 * g(3);
t620 = -t612 * qJ(2) + qJDD(2) - t586;
t603 = sin(pkin(8));
t605 = cos(pkin(8));
t625 = -pkin(2) * t605 - qJ(3) * t603;
t642 = t603 * qJD(1);
t651 = (-pkin(1) + t625) * qJDD(1) + t620 - 0.2e1 * qJD(3) * t642;
t560 = -t605 * g(1) - t652 * t603;
t650 = mrSges(3,2) * t603;
t649 = Ifges(3,6) * t605;
t601 = t603 ^ 2;
t648 = t601 * t612;
t602 = sin(pkin(9));
t647 = t602 * t603;
t604 = cos(pkin(9));
t646 = t603 * t604;
t561 = -t603 * g(1) + t652 * t605;
t580 = (-mrSges(3,1) * t605 + t650) * qJD(1);
t579 = t625 * qJD(1);
t641 = t605 * qJD(1);
t550 = t579 * t641 + t561;
t623 = -pkin(3) * t605 - pkin(6) * t646;
t644 = t651 * t604;
t529 = t623 * qJDD(1) + (-t550 + (-pkin(3) * t601 * t604 + pkin(6) * t603 * t605) * t612) * t602 + t644;
t538 = t604 * t550 + t651 * t602;
t578 = t623 * qJD(1);
t639 = qJDD(1) * t603;
t635 = t602 * t639;
t637 = t602 ^ 2 * t648;
t530 = -pkin(3) * t637 - pkin(6) * t635 + t578 * t641 + t538;
t607 = sin(qJ(4));
t610 = cos(qJ(4));
t518 = t610 * t529 - t607 * t530;
t619 = (-t602 * t610 - t604 * t607) * t603;
t568 = qJD(1) * t619;
t618 = (-t602 * t607 + t604 * t610) * t603;
t554 = t568 * qJD(4) + qJDD(1) * t618;
t569 = qJD(1) * t618;
t638 = t605 * qJDD(1);
t589 = qJDD(4) - t638;
t590 = qJD(4) - t641;
t516 = (t568 * t590 - t554) * pkin(7) + (t568 * t569 + t589) * pkin(4) + t518;
t519 = t607 * t529 + t610 * t530;
t553 = -t569 * qJD(4) + qJDD(1) * t619;
t559 = t590 * pkin(4) - t569 * pkin(7);
t567 = t568 ^ 2;
t517 = -t567 * pkin(4) + t553 * pkin(7) - t590 * t559 + t519;
t606 = sin(qJ(5));
t609 = cos(qJ(5));
t514 = t609 * t516 - t606 * t517;
t547 = t609 * t568 - t606 * t569;
t525 = t547 * qJD(5) + t606 * t553 + t609 * t554;
t548 = t606 * t568 + t609 * t569;
t536 = -t547 * mrSges(6,1) + t548 * mrSges(6,2);
t588 = qJD(5) + t590;
t540 = -t588 * mrSges(6,2) + t547 * mrSges(6,3);
t584 = qJDD(5) + t589;
t512 = m(6) * t514 + t584 * mrSges(6,1) - t525 * mrSges(6,3) - t548 * t536 + t588 * t540;
t515 = t606 * t516 + t609 * t517;
t524 = -t548 * qJD(5) + t609 * t553 - t606 * t554;
t541 = t588 * mrSges(6,1) - t548 * mrSges(6,3);
t513 = m(6) * t515 - t584 * mrSges(6,2) + t524 * mrSges(6,3) + t547 * t536 - t588 * t541;
t504 = t609 * t512 + t606 * t513;
t551 = -t568 * mrSges(5,1) + t569 * mrSges(5,2);
t555 = -t590 * mrSges(5,2) + t568 * mrSges(5,3);
t502 = m(5) * t518 + t589 * mrSges(5,1) - t554 * mrSges(5,3) - t569 * t551 + t590 * t555 + t504;
t556 = t590 * mrSges(5,1) - t569 * mrSges(5,3);
t630 = -t606 * t512 + t609 * t513;
t503 = m(5) * t519 - t589 * mrSges(5,2) + t553 * mrSges(5,3) + t568 * t551 - t590 * t556 + t630;
t498 = t610 * t502 + t607 * t503;
t537 = -t602 * t550 + t644;
t628 = mrSges(4,1) * t602 + mrSges(4,2) * t604;
t572 = t628 * t642;
t621 = mrSges(4,2) * t605 - mrSges(4,3) * t647;
t575 = t621 * qJD(1);
t622 = -mrSges(4,1) * t605 - mrSges(4,3) * t646;
t496 = m(4) * t537 + t622 * qJDD(1) + (-t572 * t646 - t575 * t605) * qJD(1) + t498;
t576 = t622 * qJD(1);
t631 = -t607 * t502 + t610 * t503;
t497 = m(4) * t538 + t621 * qJDD(1) + (-t572 * t647 + t576 * t605) * qJD(1) + t631;
t632 = -t602 * t496 + t604 * t497;
t491 = m(3) * t561 + (qJDD(1) * mrSges(3,3) + qJD(1) * t580) * t605 + t632;
t549 = t579 * t642 + qJDD(3) - t560;
t539 = t604 * t578 * t642 + pkin(3) * t635 - pkin(6) * t637 + t549;
t521 = -t553 * pkin(4) - t567 * pkin(7) + t569 * t559 + t539;
t624 = m(6) * t521 - t524 * mrSges(6,1) + t525 * mrSges(6,2) - t547 * t540 + t548 * t541;
t614 = m(5) * t539 - t553 * mrSges(5,1) + t554 * mrSges(5,2) - t568 * t555 + t569 * t556 + t624;
t613 = -m(4) * t549 - t614;
t508 = t613 + ((-mrSges(3,3) - t628) * qJDD(1) + (-t575 * t602 - t576 * t604 - t580) * qJD(1)) * t603 + m(3) * t560;
t633 = t605 * t491 - t603 * t508;
t484 = m(2) * t585 - t612 * mrSges(2,1) - qJDD(1) * mrSges(2,2) + t633;
t492 = t604 * t496 + t602 * t497;
t574 = -qJDD(1) * pkin(1) + t620;
t615 = -m(3) * t574 + mrSges(3,1) * t638 - t492 + (t605 ^ 2 * t612 + t648) * mrSges(3,3);
t488 = m(2) * t586 - t612 * mrSges(2,2) + (mrSges(2,1) - t650) * qJDD(1) + t615;
t645 = t608 * t484 + t611 * t488;
t485 = t603 * t491 + t605 * t508;
t634 = -t611 * t484 + t608 * t488;
t627 = Ifges(3,1) * t603 + Ifges(3,4) * t605;
t626 = Ifges(4,5) * t604 - Ifges(4,6) * t602;
t617 = -Ifges(4,5) * t605 + (Ifges(4,1) * t604 - Ifges(4,4) * t602) * t603;
t616 = -Ifges(4,6) * t605 + (Ifges(4,4) * t604 - Ifges(4,2) * t602) * t603;
t581 = (Ifges(3,5) * t603 + t649) * qJD(1);
t565 = t617 * qJD(1);
t564 = t616 * qJD(1);
t563 = (-Ifges(4,3) * t605 + t626 * t603) * qJD(1);
t544 = Ifges(5,1) * t569 + Ifges(5,4) * t568 + Ifges(5,5) * t590;
t543 = Ifges(5,4) * t569 + Ifges(5,2) * t568 + Ifges(5,6) * t590;
t542 = Ifges(5,5) * t569 + Ifges(5,6) * t568 + Ifges(5,3) * t590;
t533 = Ifges(6,1) * t548 + Ifges(6,4) * t547 + Ifges(6,5) * t588;
t532 = Ifges(6,4) * t548 + Ifges(6,2) * t547 + Ifges(6,6) * t588;
t531 = Ifges(6,5) * t548 + Ifges(6,6) * t547 + Ifges(6,3) * t588;
t506 = mrSges(6,2) * t521 - mrSges(6,3) * t514 + Ifges(6,1) * t525 + Ifges(6,4) * t524 + Ifges(6,5) * t584 + t547 * t531 - t588 * t532;
t505 = -mrSges(6,1) * t521 + mrSges(6,3) * t515 + Ifges(6,4) * t525 + Ifges(6,2) * t524 + Ifges(6,6) * t584 - t548 * t531 + t588 * t533;
t494 = mrSges(5,2) * t539 - mrSges(5,3) * t518 + Ifges(5,1) * t554 + Ifges(5,4) * t553 + Ifges(5,5) * t589 - pkin(7) * t504 - t606 * t505 + t609 * t506 + t568 * t542 - t590 * t543;
t493 = -mrSges(5,1) * t539 + mrSges(5,3) * t519 + Ifges(5,4) * t554 + Ifges(5,2) * t553 + Ifges(5,6) * t589 - pkin(4) * t624 + pkin(7) * t630 + t609 * t505 + t606 * t506 - t569 * t542 + t590 * t544;
t482 = mrSges(4,2) * t549 - mrSges(4,3) * t537 - pkin(6) * t498 - t607 * t493 + t610 * t494 + (-t563 * t647 + t564 * t605) * qJD(1) + t617 * qJDD(1);
t481 = -mrSges(4,1) * t549 + mrSges(4,3) * t538 + t607 * t494 + t610 * t493 - pkin(3) * t614 + pkin(6) * t631 + (-t563 * t646 - t605 * t565) * qJD(1) + t616 * qJDD(1);
t480 = -Ifges(5,3) * t589 - mrSges(3,1) * t574 - Ifges(6,3) * t584 + mrSges(3,3) * t561 + t568 * t544 - t569 * t543 - t548 * t532 - Ifges(5,6) * t553 - Ifges(5,5) * t554 + mrSges(4,2) * t538 + t547 * t533 - mrSges(4,1) * t537 - Ifges(6,6) * t524 - Ifges(6,5) * t525 - mrSges(5,1) * t518 + mrSges(5,2) * t519 - mrSges(6,1) * t514 + mrSges(6,2) * t515 - pkin(4) * t504 - pkin(3) * t498 - pkin(2) * t492 + (Ifges(3,2) + Ifges(4,3)) * t638 + ((Ifges(3,4) - t626) * qJDD(1) + (-t564 * t604 - t565 * t602 - t581) * qJD(1)) * t603;
t479 = mrSges(3,2) * t574 - mrSges(3,3) * t560 - qJ(3) * t492 + t627 * qJDD(1) - t602 * t481 + t604 * t482 + t581 * t641;
t478 = t612 * Ifges(2,5) + mrSges(2,1) * g(1) + mrSges(2,3) * t585 - mrSges(3,1) * t560 + mrSges(3,2) * t561 - t602 * t482 - t604 * t481 - pkin(2) * t613 - qJ(3) * t632 - pkin(1) * t485 + (-t649 + Ifges(2,6) + (pkin(2) * t628 - Ifges(3,5)) * t603) * qJDD(1) + (-pkin(2) * (-t575 * t647 - t576 * t646) + (-t603 * (Ifges(3,4) * t603 + Ifges(3,2) * t605) + t605 * t627) * qJD(1)) * qJD(1);
t477 = -mrSges(2,2) * g(1) - mrSges(2,3) * t586 + Ifges(2,5) * qJDD(1) - t612 * Ifges(2,6) - qJ(2) * t485 + t605 * t479 - t603 * t480;
t1 = [(-m(1) - m(2)) * g(1) + t485; -m(1) * g(2) + t645; -m(1) * g(3) + t634; -mrSges(1,2) * g(3) + mrSges(1,3) * g(2) + Ifges(2,3) * qJDD(1) + mrSges(2,1) * t586 - mrSges(2,2) * t585 + t603 * t479 + t605 * t480 + pkin(1) * (-mrSges(3,2) * t639 + t615) + qJ(2) * t633; mrSges(1,1) * g(3) - mrSges(1,3) * g(1) - pkin(5) * t634 + t608 * t477 + t611 * t478; -mrSges(1,1) * g(2) + mrSges(1,2) * g(1) + pkin(5) * t645 - t611 * t477 + t608 * t478;];
tauB = t1;
