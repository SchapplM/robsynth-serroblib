% Calculate vector of inverse dynamics base forces with Newton-Euler for
% S6RPRPRP3
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
%   pkin=[a2,a3,a4,a5,a6,d1,d3,d5,theta2,theta4]';
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
% Datum: 2019-05-05 17:41
% Revision: 89c353f7eff3bd693eda4e29f35b2761dbc3ada0 (2019-05-03)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauB = S6RPRPRP3_invdynB_fixb_snew_vp2(qJ, qJD, qJDD, g, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(6,1),zeros(3,1),zeros(10,1),zeros(7,1),zeros(7,3),zeros(7,6)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPRPRP3_invdynB_fixb_snew_vp2: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RPRPRP3_invdynB_fixb_snew_vp2: qJD has to be [6x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [6 1]), ...
  'S6RPRPRP3_invdynB_fixb_snew_vp2: qJDD has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6RPRPRP3_invdynB_fixb_snew_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S6RPRPRP3_invdynB_fixb_snew_vp2: pkin has to be [10x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6RPRPRP3_invdynB_fixb_snew_vp2: m has to be [7x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [7,3]), ...
  'S6RPRPRP3_invdynB_fixb_snew_vp2: mrSges has to be [7x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [7 6]), ...
  'S6RPRPRP3_invdynB_fixb_snew_vp2: Ifges has to be [7x6] (double)'); 

%% Symbolic Calculation
% From invdyn_fixb_NewtonEuler_linkframe_tauB_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-05-05 17:37:52
% EndTime: 2019-05-05 17:38:02
% DurationCPUTime: 8.01s
% Computational Cost: add. (95098->318), mult. (195445->387), div. (0->0), fcn. (125838->10), ass. (0->124)
t662 = Ifges(6,1) + Ifges(7,1);
t656 = Ifges(6,4) - Ifges(7,5);
t661 = -Ifges(6,5) - Ifges(7,4);
t660 = Ifges(6,2) + Ifges(7,3);
t654 = Ifges(6,6) - Ifges(7,6);
t659 = -Ifges(6,3) - Ifges(7,2);
t658 = cos(qJ(5));
t657 = -mrSges(6,3) - mrSges(7,2);
t629 = sin(qJ(1));
t631 = cos(qJ(1));
t614 = t629 * g(1) - t631 * g(2);
t606 = qJDD(1) * pkin(1) + t614;
t615 = -t631 * g(1) - t629 * g(2);
t633 = qJD(1) ^ 2;
t609 = -t633 * pkin(1) + t615;
t624 = sin(pkin(9));
t626 = cos(pkin(9));
t579 = t624 * t606 + t626 * t609;
t570 = -t633 * pkin(2) + qJDD(1) * pkin(7) + t579;
t622 = -g(3) + qJDD(2);
t628 = sin(qJ(3));
t630 = cos(qJ(3));
t562 = t630 * t570 + t628 * t622;
t608 = (-mrSges(4,1) * t630 + mrSges(4,2) * t628) * qJD(1);
t646 = qJD(1) * qJD(3);
t619 = t628 * t646;
t611 = t630 * qJDD(1) - t619;
t648 = qJD(1) * t628;
t612 = qJD(3) * mrSges(4,1) - mrSges(4,3) * t648;
t578 = t626 * t606 - t624 * t609;
t569 = -qJDD(1) * pkin(2) - t633 * pkin(7) - t578;
t643 = t630 * t646;
t610 = t628 * qJDD(1) + t643;
t552 = (-t610 - t643) * qJ(4) + (-t611 + t619) * pkin(3) + t569;
t607 = (-pkin(3) * t630 - qJ(4) * t628) * qJD(1);
t632 = qJD(3) ^ 2;
t647 = t630 * qJD(1);
t558 = -t632 * pkin(3) + qJDD(3) * qJ(4) + t607 * t647 + t562;
t623 = sin(pkin(10));
t625 = cos(pkin(10));
t603 = t623 * qJD(3) + t625 * t648;
t533 = -0.2e1 * qJD(4) * t603 + t625 * t552 - t623 * t558;
t586 = t623 * qJDD(3) + t625 * t610;
t602 = t625 * qJD(3) - t623 * t648;
t530 = (-t602 * t647 - t586) * pkin(8) + (t602 * t603 - t611) * pkin(4) + t533;
t534 = 0.2e1 * qJD(4) * t602 + t623 * t552 + t625 * t558;
t585 = t625 * qJDD(3) - t623 * t610;
t587 = -pkin(4) * t647 - t603 * pkin(8);
t601 = t602 ^ 2;
t532 = -t601 * pkin(4) + t585 * pkin(8) + t587 * t647 + t534;
t627 = sin(qJ(5));
t526 = t627 * t530 + t658 * t532;
t577 = t627 * t602 + t658 * t603;
t539 = t577 * qJD(5) - t658 * t585 + t627 * t586;
t617 = qJD(5) - t647;
t564 = t617 * mrSges(6,1) - t577 * mrSges(6,3);
t576 = -t658 * t602 + t627 * t603;
t605 = qJDD(5) - t611;
t555 = t576 * pkin(5) - t577 * qJ(6);
t616 = t617 ^ 2;
t523 = -t616 * pkin(5) + t605 * qJ(6) + 0.2e1 * qJD(6) * t617 - t576 * t555 + t526;
t565 = -t617 * mrSges(7,1) + t577 * mrSges(7,2);
t645 = m(7) * t523 + t605 * mrSges(7,3) + t617 * t565;
t556 = t576 * mrSges(7,1) - t577 * mrSges(7,3);
t649 = -t576 * mrSges(6,1) - t577 * mrSges(6,2) - t556;
t516 = m(6) * t526 - t605 * mrSges(6,2) + t657 * t539 - t617 * t564 + t649 * t576 + t645;
t525 = t658 * t530 - t627 * t532;
t540 = -t576 * qJD(5) + t627 * t585 + t658 * t586;
t563 = -t617 * mrSges(6,2) - t576 * mrSges(6,3);
t524 = -t605 * pkin(5) - t616 * qJ(6) + t577 * t555 + qJDD(6) - t525;
t566 = -t576 * mrSges(7,2) + t617 * mrSges(7,3);
t637 = -m(7) * t524 + t605 * mrSges(7,1) + t617 * t566;
t518 = m(6) * t525 + t605 * mrSges(6,1) + t657 * t540 + t617 * t563 + t649 * t577 + t637;
t511 = t627 * t516 + t658 * t518;
t580 = -t602 * mrSges(5,1) + t603 * mrSges(5,2);
t583 = mrSges(5,2) * t647 + t602 * mrSges(5,3);
t509 = m(5) * t533 - t611 * mrSges(5,1) - t586 * mrSges(5,3) - t603 * t580 - t583 * t647 + t511;
t584 = -mrSges(5,1) * t647 - t603 * mrSges(5,3);
t638 = t658 * t516 - t627 * t518;
t510 = m(5) * t534 + t611 * mrSges(5,2) + t585 * mrSges(5,3) + t602 * t580 + t584 * t647 + t638;
t639 = -t623 * t509 + t625 * t510;
t506 = m(4) * t562 - qJDD(3) * mrSges(4,2) + t611 * mrSges(4,3) - qJD(3) * t612 + t608 * t647 + t639;
t561 = -t628 * t570 + t630 * t622;
t613 = -qJD(3) * mrSges(4,2) + mrSges(4,3) * t647;
t554 = -qJDD(3) * pkin(3) - t632 * qJ(4) + t607 * t648 + qJDD(4) - t561;
t535 = -t585 * pkin(4) - t601 * pkin(8) + t603 * t587 + t554;
t528 = -0.2e1 * qJD(6) * t577 + (t576 * t617 - t540) * qJ(6) + (t577 * t617 + t539) * pkin(5) + t535;
t521 = m(7) * t528 + t539 * mrSges(7,1) - t540 * mrSges(7,3) - t577 * t565 + t576 * t566;
t636 = m(6) * t535 + t539 * mrSges(6,1) + t540 * mrSges(6,2) + t576 * t563 + t577 * t564 + t521;
t634 = -m(5) * t554 + t585 * mrSges(5,1) - t586 * mrSges(5,2) + t602 * t583 - t603 * t584 - t636;
t520 = m(4) * t561 + qJDD(3) * mrSges(4,1) - t610 * mrSges(4,3) + qJD(3) * t613 - t608 * t648 + t634;
t640 = t630 * t506 - t628 * t520;
t500 = m(3) * t579 - t633 * mrSges(3,1) - qJDD(1) * mrSges(3,2) + t640;
t507 = t625 * t509 + t623 * t510;
t635 = -m(4) * t569 + t611 * mrSges(4,1) - t610 * mrSges(4,2) - t612 * t648 + t613 * t647 - t507;
t503 = m(3) * t578 + qJDD(1) * mrSges(3,1) - t633 * mrSges(3,2) + t635;
t494 = t624 * t500 + t626 * t503;
t492 = m(2) * t614 + qJDD(1) * mrSges(2,1) - t633 * mrSges(2,2) + t494;
t641 = t626 * t500 - t624 * t503;
t493 = m(2) * t615 - t633 * mrSges(2,1) - qJDD(1) * mrSges(2,2) + t641;
t653 = t631 * t492 + t629 * t493;
t501 = t628 * t506 + t630 * t520;
t652 = t660 * t576 - t656 * t577 - t654 * t617;
t651 = t654 * t576 + t661 * t577 + t659 * t617;
t650 = -t656 * t576 + t662 * t577 - t661 * t617;
t644 = m(3) * t622 + t501;
t642 = -t629 * t492 + t631 * t493;
t597 = Ifges(4,5) * qJD(3) + (Ifges(4,1) * t628 + Ifges(4,4) * t630) * qJD(1);
t596 = Ifges(4,6) * qJD(3) + (Ifges(4,4) * t628 + Ifges(4,2) * t630) * qJD(1);
t595 = Ifges(4,3) * qJD(3) + (Ifges(4,5) * t628 + Ifges(4,6) * t630) * qJD(1);
t573 = Ifges(5,1) * t603 + Ifges(5,4) * t602 - Ifges(5,5) * t647;
t572 = Ifges(5,4) * t603 + Ifges(5,2) * t602 - Ifges(5,6) * t647;
t571 = Ifges(5,5) * t603 + Ifges(5,6) * t602 - Ifges(5,3) * t647;
t513 = mrSges(6,2) * t535 + mrSges(7,2) * t524 - mrSges(6,3) * t525 - mrSges(7,3) * t528 - qJ(6) * t521 - t656 * t539 + t662 * t540 + t651 * t576 - t605 * t661 + t652 * t617;
t512 = -mrSges(6,1) * t535 - mrSges(7,1) * t528 + mrSges(7,2) * t523 + mrSges(6,3) * t526 - pkin(5) * t521 - t660 * t539 + t656 * t540 + t651 * t577 + t654 * t605 + t650 * t617;
t497 = mrSges(5,2) * t554 - mrSges(5,3) * t533 + Ifges(5,1) * t586 + Ifges(5,4) * t585 - Ifges(5,5) * t611 - pkin(8) * t511 - t627 * t512 + t658 * t513 + t602 * t571 + t572 * t647;
t496 = -mrSges(5,1) * t554 + mrSges(5,3) * t534 + Ifges(5,4) * t586 + Ifges(5,2) * t585 - Ifges(5,6) * t611 - pkin(4) * t636 + pkin(8) * t638 + t658 * t512 + t627 * t513 - t603 * t571 - t573 * t647;
t495 = (pkin(5) * mrSges(7,2) + t661) * t540 + t659 * t605 + t602 * t573 - t603 * t572 + Ifges(4,4) * t610 - Ifges(5,6) * t585 - Ifges(5,5) * t586 + qJD(3) * t597 - mrSges(4,1) * t569 + mrSges(4,3) * t562 - mrSges(5,1) * t533 + mrSges(5,2) * t534 - mrSges(6,1) * t525 + mrSges(6,2) * t526 - mrSges(7,3) * t523 + mrSges(7,1) * t524 - pkin(4) * t511 - pkin(3) * t507 + Ifges(4,6) * qJDD(3) - pkin(5) * t637 - qJ(6) * t645 - t595 * t648 + (qJ(6) * t556 - t650) * t576 + (pkin(5) * t556 + t652) * t577 + (qJ(6) * mrSges(7,2) + t654) * t539 + (Ifges(4,2) + Ifges(5,3)) * t611;
t488 = mrSges(4,2) * t569 - mrSges(4,3) * t561 + Ifges(4,1) * t610 + Ifges(4,4) * t611 + Ifges(4,5) * qJDD(3) - qJ(4) * t507 - qJD(3) * t596 - t623 * t496 + t625 * t497 + t595 * t647;
t487 = Ifges(3,6) * qJDD(1) + t633 * Ifges(3,5) - mrSges(3,1) * t622 + mrSges(3,3) * t579 - Ifges(4,5) * t610 - Ifges(4,6) * t611 - Ifges(4,3) * qJDD(3) - mrSges(4,1) * t561 + mrSges(4,2) * t562 - t623 * t497 - t625 * t496 - pkin(3) * t634 - qJ(4) * t639 - pkin(2) * t501 + (-t628 * t596 + t630 * t597) * qJD(1);
t486 = mrSges(3,2) * t622 - mrSges(3,3) * t578 + Ifges(3,5) * qJDD(1) - t633 * Ifges(3,6) - pkin(7) * t501 + t630 * t488 - t628 * t495;
t485 = -mrSges(2,2) * g(3) - mrSges(2,3) * t614 + Ifges(2,5) * qJDD(1) - t633 * Ifges(2,6) - qJ(2) * t494 + t626 * t486 - t624 * t487;
t484 = mrSges(2,1) * g(3) + mrSges(2,3) * t615 + t633 * Ifges(2,5) + Ifges(2,6) * qJDD(1) - pkin(1) * t644 + qJ(2) * t641 + t624 * t486 + t626 * t487;
t1 = [-m(1) * g(1) + t642; -m(1) * g(2) + t653; (-m(1) - m(2)) * g(3) + t644; -mrSges(1,2) * g(3) + mrSges(1,3) * g(2) - pkin(6) * t653 - t629 * t484 + t631 * t485; mrSges(1,1) * g(3) - mrSges(1,3) * g(1) + pkin(6) * t642 + t631 * t484 + t629 * t485; pkin(1) * t494 + mrSges(2,1) * t614 - mrSges(2,2) * t615 + t628 * t488 + t630 * t495 + pkin(2) * t635 + pkin(7) * t640 + mrSges(3,1) * t578 - mrSges(3,2) * t579 - mrSges(1,1) * g(2) + mrSges(1,2) * g(1) + (Ifges(2,3) + Ifges(3,3)) * qJDD(1);];
tauB  = t1;
