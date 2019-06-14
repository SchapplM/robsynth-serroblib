% Calculate vector of inverse dynamics base forces with Newton-Euler for
% S6RPRPRR10
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
%   pkin=[a2,a3,a4,a5,a6,d1,d3,d5,d6,theta4]';
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
% Datum: 2019-05-05 20:04
% Revision: 89c353f7eff3bd693eda4e29f35b2761dbc3ada0 (2019-05-03)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauB = S6RPRPRR10_invdynB_fixb_snew_vp2(qJ, qJD, qJDD, g, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(6,1),zeros(3,1),zeros(10,1),zeros(7,1),zeros(7,3),zeros(7,6)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPRPRR10_invdynB_fixb_snew_vp2: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RPRPRR10_invdynB_fixb_snew_vp2: qJD has to be [6x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [6 1]), ...
  'S6RPRPRR10_invdynB_fixb_snew_vp2: qJDD has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6RPRPRR10_invdynB_fixb_snew_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S6RPRPRR10_invdynB_fixb_snew_vp2: pkin has to be [10x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6RPRPRR10_invdynB_fixb_snew_vp2: m has to be [7x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [7,3]), ...
  'S6RPRPRR10_invdynB_fixb_snew_vp2: mrSges has to be [7x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [7 6]), ...
  'S6RPRPRR10_invdynB_fixb_snew_vp2: Ifges has to be [7x6] (double)'); 

%% Symbolic Calculation
% From invdyn_fixb_NewtonEuler_linkframe_tauB_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-05-05 20:02:04
% EndTime: 2019-05-05 20:02:14
% DurationCPUTime: 7.93s
% Computational Cost: add. (127010->337), mult. (268375->413), div. (0->0), fcn. (178861->10), ass. (0->131)
t636 = sin(qJ(1));
t640 = cos(qJ(1));
t621 = -t640 * g(1) - t636 * g(2);
t666 = -qJDD(1) * qJ(2) - (2 * qJD(2) * qJD(1)) - t621;
t665 = -pkin(1) - pkin(7);
t664 = mrSges(2,1) - mrSges(3,2);
t663 = -Ifges(3,4) + Ifges(2,5);
t662 = (Ifges(3,5) - Ifges(2,6));
t620 = t636 * g(1) - t640 * g(2);
t642 = qJD(1) ^ 2;
t648 = -t642 * qJ(2) + qJDD(2) - t620;
t597 = t665 * qJDD(1) + t648;
t635 = sin(qJ(3));
t639 = cos(qJ(3));
t587 = -t639 * g(3) + t635 * t597;
t615 = (mrSges(4,1) * t635 + mrSges(4,2) * t639) * qJD(1);
t659 = qJD(1) * qJD(3);
t656 = t639 * t659;
t616 = t635 * qJDD(1) + t656;
t660 = qJD(1) * t639;
t619 = qJD(3) * mrSges(4,1) - mrSges(4,3) * t660;
t624 = t635 * qJD(1);
t591 = t665 * t642 - t666;
t657 = t635 * t659;
t617 = t639 * qJDD(1) - t657;
t570 = (-t617 + t657) * qJ(4) + (t616 + t656) * pkin(3) + t591;
t614 = (pkin(3) * t635 - qJ(4) * t639) * qJD(1);
t641 = qJD(3) ^ 2;
t573 = -t641 * pkin(3) + qJDD(3) * qJ(4) - t614 * t624 + t587;
t631 = sin(pkin(10));
t632 = cos(pkin(10));
t611 = t631 * qJD(3) + t632 * t660;
t553 = -0.2e1 * qJD(4) * t611 + t632 * t570 - t631 * t573;
t595 = t631 * qJDD(3) + t632 * t617;
t610 = t632 * qJD(3) - t631 * t660;
t544 = (t610 * t624 - t595) * pkin(8) + (t610 * t611 + t616) * pkin(4) + t553;
t554 = 0.2e1 * qJD(4) * t610 + t631 * t570 + t632 * t573;
t594 = t632 * qJDD(3) - t631 * t617;
t596 = pkin(4) * t624 - t611 * pkin(8);
t609 = t610 ^ 2;
t546 = -t609 * pkin(4) + t594 * pkin(8) - t596 * t624 + t554;
t634 = sin(qJ(5));
t638 = cos(qJ(5));
t534 = t638 * t544 - t634 * t546;
t583 = t638 * t610 - t634 * t611;
t560 = t583 * qJD(5) + t634 * t594 + t638 * t595;
t584 = t634 * t610 + t638 * t611;
t613 = qJDD(5) + t616;
t623 = t624 + qJD(5);
t532 = (t583 * t623 - t560) * pkin(9) + (t583 * t584 + t613) * pkin(5) + t534;
t535 = t634 * t544 + t638 * t546;
t559 = -t584 * qJD(5) + t638 * t594 - t634 * t595;
t576 = t623 * pkin(5) - t584 * pkin(9);
t582 = t583 ^ 2;
t533 = -t582 * pkin(5) + t559 * pkin(9) - t623 * t576 + t535;
t633 = sin(qJ(6));
t637 = cos(qJ(6));
t530 = t637 * t532 - t633 * t533;
t565 = t637 * t583 - t633 * t584;
t541 = t565 * qJD(6) + t633 * t559 + t637 * t560;
t566 = t633 * t583 + t637 * t584;
t552 = -t565 * mrSges(7,1) + t566 * mrSges(7,2);
t622 = qJD(6) + t623;
t556 = -t622 * mrSges(7,2) + t565 * mrSges(7,3);
t606 = qJDD(6) + t613;
t528 = m(7) * t530 + t606 * mrSges(7,1) - t541 * mrSges(7,3) - t566 * t552 + t622 * t556;
t531 = t633 * t532 + t637 * t533;
t540 = -t566 * qJD(6) + t637 * t559 - t633 * t560;
t557 = t622 * mrSges(7,1) - t566 * mrSges(7,3);
t529 = m(7) * t531 - t606 * mrSges(7,2) + t540 * mrSges(7,3) + t565 * t552 - t622 * t557;
t521 = t637 * t528 + t633 * t529;
t567 = -t583 * mrSges(6,1) + t584 * mrSges(6,2);
t574 = -t623 * mrSges(6,2) + t583 * mrSges(6,3);
t519 = m(6) * t534 + t613 * mrSges(6,1) - t560 * mrSges(6,3) - t584 * t567 + t623 * t574 + t521;
t575 = t623 * mrSges(6,1) - t584 * mrSges(6,3);
t651 = -t633 * t528 + t637 * t529;
t520 = m(6) * t535 - t613 * mrSges(6,2) + t559 * mrSges(6,3) + t583 * t567 - t623 * t575 + t651;
t515 = t638 * t519 + t634 * t520;
t585 = -t610 * mrSges(5,1) + t611 * mrSges(5,2);
t592 = -mrSges(5,2) * t624 + t610 * mrSges(5,3);
t513 = m(5) * t553 + t616 * mrSges(5,1) - t595 * mrSges(5,3) - t611 * t585 + t592 * t624 + t515;
t593 = mrSges(5,1) * t624 - t611 * mrSges(5,3);
t652 = -t634 * t519 + t638 * t520;
t514 = m(5) * t554 - t616 * mrSges(5,2) + t594 * mrSges(5,3) + t610 * t585 - t593 * t624 + t652;
t653 = -t631 * t513 + t632 * t514;
t506 = m(4) * t587 - qJDD(3) * mrSges(4,2) - t616 * mrSges(4,3) - qJD(3) * t619 - t615 * t624 + t653;
t586 = t635 * g(3) + t639 * t597;
t618 = -qJD(3) * mrSges(4,2) - mrSges(4,3) * t624;
t572 = -qJDD(3) * pkin(3) - t641 * qJ(4) + t614 * t660 + qJDD(4) - t586;
t555 = -t594 * pkin(4) - t609 * pkin(8) + t611 * t596 + t572;
t537 = -t559 * pkin(5) - t582 * pkin(9) + t584 * t576 + t555;
t650 = m(7) * t537 - t540 * mrSges(7,1) + t541 * mrSges(7,2) - t565 * t556 + t566 * t557;
t644 = m(6) * t555 - t559 * mrSges(6,1) + t560 * mrSges(6,2) - t583 * t574 + t584 * t575 + t650;
t643 = -m(5) * t572 + t594 * mrSges(5,1) - t595 * mrSges(5,2) + t610 * t592 - t611 * t593 - t644;
t524 = m(4) * t586 + qJDD(3) * mrSges(4,1) - t617 * mrSges(4,3) + qJD(3) * t618 - t615 * t660 + t643;
t501 = t635 * t506 + t639 * t524;
t599 = -qJDD(1) * pkin(1) + t648;
t647 = -m(3) * t599 + (t642 * mrSges(3,3)) - t501;
t499 = m(2) * t620 - (t642 * mrSges(2,2)) + t664 * qJDD(1) + t647;
t598 = t642 * pkin(1) + t666;
t507 = t632 * t513 + t631 * t514;
t646 = -m(4) * t591 - t616 * mrSges(4,1) - t617 * mrSges(4,2) - t618 * t624 - t619 * t660 - t507;
t645 = -m(3) * t598 + (t642 * mrSges(3,2)) + qJDD(1) * mrSges(3,3) - t646;
t504 = m(2) * t621 - (t642 * mrSges(2,1)) - qJDD(1) * mrSges(2,2) + t645;
t661 = t640 * t499 + t636 * t504;
t655 = -t636 * t499 + t640 * t504;
t654 = t639 * t506 - t635 * t524;
t605 = (Ifges(4,5) * qJD(3)) + (Ifges(4,1) * t639 - Ifges(4,4) * t635) * qJD(1);
t604 = (Ifges(4,6) * qJD(3)) + (Ifges(4,4) * t639 - Ifges(4,2) * t635) * qJD(1);
t603 = (Ifges(4,3) * qJD(3)) + (Ifges(4,5) * t639 - Ifges(4,6) * t635) * qJD(1);
t579 = Ifges(5,1) * t611 + Ifges(5,4) * t610 + Ifges(5,5) * t624;
t578 = Ifges(5,4) * t611 + Ifges(5,2) * t610 + Ifges(5,6) * t624;
t577 = Ifges(5,5) * t611 + Ifges(5,6) * t610 + Ifges(5,3) * t624;
t563 = Ifges(6,1) * t584 + Ifges(6,4) * t583 + Ifges(6,5) * t623;
t562 = Ifges(6,4) * t584 + Ifges(6,2) * t583 + Ifges(6,6) * t623;
t561 = Ifges(6,5) * t584 + Ifges(6,6) * t583 + Ifges(6,3) * t623;
t549 = Ifges(7,1) * t566 + Ifges(7,4) * t565 + Ifges(7,5) * t622;
t548 = Ifges(7,4) * t566 + Ifges(7,2) * t565 + Ifges(7,6) * t622;
t547 = Ifges(7,5) * t566 + Ifges(7,6) * t565 + Ifges(7,3) * t622;
t523 = mrSges(7,2) * t537 - mrSges(7,3) * t530 + Ifges(7,1) * t541 + Ifges(7,4) * t540 + Ifges(7,5) * t606 + t565 * t547 - t622 * t548;
t522 = -mrSges(7,1) * t537 + mrSges(7,3) * t531 + Ifges(7,4) * t541 + Ifges(7,2) * t540 + Ifges(7,6) * t606 - t566 * t547 + t622 * t549;
t509 = mrSges(6,2) * t555 - mrSges(6,3) * t534 + Ifges(6,1) * t560 + Ifges(6,4) * t559 + Ifges(6,5) * t613 - pkin(9) * t521 - t633 * t522 + t637 * t523 + t583 * t561 - t623 * t562;
t508 = -mrSges(6,1) * t555 + mrSges(6,3) * t535 + Ifges(6,4) * t560 + Ifges(6,2) * t559 + Ifges(6,6) * t613 - pkin(5) * t650 + pkin(9) * t651 + t637 * t522 + t633 * t523 - t584 * t561 + t623 * t563;
t500 = -m(3) * g(3) + t654;
t497 = mrSges(5,2) * t572 - mrSges(5,3) * t553 + Ifges(5,1) * t595 + Ifges(5,4) * t594 + Ifges(5,5) * t616 - pkin(8) * t515 - t634 * t508 + t638 * t509 + t610 * t577 - t578 * t624;
t496 = -mrSges(5,1) * t572 + mrSges(5,3) * t554 + Ifges(5,4) * t595 + Ifges(5,2) * t594 + Ifges(5,6) * t616 - pkin(4) * t644 + pkin(8) * t652 + t638 * t508 + t634 * t509 - t611 * t577 + t579 * t624;
t495 = Ifges(4,6) * qJDD(3) + (-Ifges(5,3) - Ifges(4,2)) * t616 + Ifges(4,4) * t617 - Ifges(7,3) * t606 + t610 * t579 - t611 * t578 - Ifges(6,3) * t613 - Ifges(5,6) * t594 - Ifges(5,5) * t595 + qJD(3) * t605 + t583 * t563 - t584 * t562 + mrSges(4,3) * t587 - mrSges(4,1) * t591 - Ifges(6,5) * t560 + t565 * t549 - t566 * t548 - mrSges(5,1) * t553 + mrSges(5,2) * t554 - Ifges(6,6) * t559 - Ifges(7,6) * t540 - Ifges(7,5) * t541 - mrSges(6,1) * t534 + mrSges(6,2) * t535 + mrSges(7,2) * t531 - mrSges(7,1) * t530 - pkin(5) * t521 - pkin(4) * t515 - pkin(3) * t507 - t603 * t660;
t494 = mrSges(4,2) * t591 - mrSges(4,3) * t586 + Ifges(4,1) * t617 - Ifges(4,4) * t616 + Ifges(4,5) * qJDD(3) - qJ(4) * t507 - qJD(3) * t604 - t631 * t496 + t632 * t497 - t603 * t624;
t493 = -qJ(2) * t500 - mrSges(2,3) * t620 + pkin(2) * t501 + mrSges(3,1) * t599 + qJ(4) * t653 + t631 * t497 + t632 * t496 + pkin(3) * t643 + Ifges(4,5) * t617 - Ifges(4,6) * t616 + Ifges(4,3) * qJDD(3) + mrSges(4,1) * t586 - mrSges(4,2) * t587 + (t662 * t642) + t663 * qJDD(1) + (t639 * t604 + t635 * t605) * qJD(1) + (-mrSges(2,2) + mrSges(3,3)) * g(3);
t492 = -mrSges(3,1) * t598 + mrSges(2,3) * t621 - pkin(1) * t500 - pkin(2) * t646 - pkin(7) * t654 + t664 * g(3) - t662 * qJDD(1) - t635 * t494 - t639 * t495 + t663 * t642;
t1 = [-m(1) * g(1) + t655; -m(1) * g(2) + t661; (-m(1) - m(2) - m(3)) * g(3) + t654; -mrSges(1,2) * g(3) + mrSges(1,3) * g(2) - pkin(6) * t661 - t636 * t492 + t640 * t493; mrSges(1,1) * g(3) - mrSges(1,3) * g(1) + pkin(6) * t655 + t640 * t492 + t636 * t493; pkin(1) * t647 + qJ(2) * t645 + t639 * t494 - t635 * t495 - pkin(7) * t501 + mrSges(2,1) * t620 - mrSges(2,2) * t621 + mrSges(3,2) * t599 - mrSges(3,3) * t598 - mrSges(1,1) * g(2) + mrSges(1,2) * g(1) + (-pkin(1) * mrSges(3,2) + Ifges(3,1) + Ifges(2,3)) * qJDD(1);];
tauB  = t1;
