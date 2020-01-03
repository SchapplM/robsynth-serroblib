% Calculate vector of inverse dynamics base forces with Newton-Euler for
% S5RRPPP1
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
% pkin [8x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,alpha3,d1,d2,theta3]';
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
% Datum: 2019-12-31 19:25
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauB = S5RRPPP1_invdynB_fixb_snew_vp2(qJ, qJD, qJDD, g, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(5,1),zeros(3,1),zeros(8,1),zeros(6,1),zeros(6,3),zeros(6,6)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RRPPP1_invdynB_fixb_snew_vp2: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RRPPP1_invdynB_fixb_snew_vp2: qJD has to be [5x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [5 1]), ...
  'S5RRPPP1_invdynB_fixb_snew_vp2: qJDD has to be [5x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S5RRPPP1_invdynB_fixb_snew_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S5RRPPP1_invdynB_fixb_snew_vp2: pkin has to be [8x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5RRPPP1_invdynB_fixb_snew_vp2: m has to be [6x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [6,3]), ...
  'S5RRPPP1_invdynB_fixb_snew_vp2: mrSges has to be [6x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [6 6]), ...
  'S5RRPPP1_invdynB_fixb_snew_vp2: Ifges has to be [6x6] (double)'); 

%% Symbolic Calculation
% From invdyn_fixb_NewtonEuler_linkframe_tauB_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 19:23:43
% EndTime: 2019-12-31 19:23:51
% DurationCPUTime: 4.55s
% Computational Cost: add. (39974->293), mult. (100460->358), div. (0->0), fcn. (69794->8), ass. (0->120)
t667 = -2 * qJD(3);
t666 = Ifges(4,1) + Ifges(5,2) + Ifges(6,3);
t665 = Ifges(4,4) + Ifges(5,6) - Ifges(6,6);
t644 = Ifges(4,5) - Ifges(5,4) + Ifges(6,5);
t664 = -Ifges(4,2) - Ifges(5,3) - Ifges(6,2);
t643 = Ifges(4,6) - Ifges(5,5) - Ifges(6,4);
t663 = Ifges(4,3) + Ifges(5,1) + Ifges(6,1);
t616 = cos(qJ(2));
t657 = cos(pkin(5));
t632 = qJD(1) * t657;
t613 = sin(pkin(5));
t648 = qJD(2) * t613;
t587 = (t616 * t632 + t648) * qJ(3);
t615 = sin(qJ(1));
t617 = cos(qJ(1));
t606 = t615 * g(1) - t617 * g(2);
t618 = qJD(1) ^ 2;
t591 = -qJDD(1) * pkin(1) - t618 * pkin(7) - t606;
t614 = sin(qJ(2));
t655 = qJ(3) * t614;
t599 = qJD(2) * pkin(2) - t632 * t655;
t646 = qJD(1) * qJD(2);
t601 = t614 * qJDD(1) + t616 * t646;
t602 = t616 * qJDD(1) - t614 * t646;
t539 = -t601 * t613 * qJ(3) - t602 * pkin(2) + (-t587 * t616 + t599 * t614) * qJD(1) + t591;
t607 = -t617 * g(1) - t615 * g(2);
t592 = -t618 * pkin(1) + qJDD(1) * pkin(7) + t607;
t593 = (-pkin(2) * t616 - t613 * t655) * qJD(1);
t636 = qJ(3) * t657;
t660 = t616 * g(3);
t540 = -t601 * t636 + qJDD(2) * pkin(2) - t660 + qJD(2) * t587 + (-qJD(1) * t593 - t592) * t614;
t580 = -t614 * g(3) + t616 * t592;
t625 = qJDD(2) * t613 + t657 * t602;
t649 = qJD(1) * t616;
t541 = t625 * qJ(3) - qJD(2) * t599 + t593 * t649 + t580;
t612 = sin(pkin(8));
t635 = t612 * t657;
t656 = cos(pkin(8));
t578 = t612 * t648 + (t656 * t614 + t616 * t635) * qJD(1);
t628 = t657 * t656;
t634 = t613 * t656;
t529 = t539 * t634 + t540 * t628 - t612 * t541 + t578 * t667;
t650 = qJD(1) * t614;
t577 = -qJD(2) * t634 + t612 * t650 - t628 * t649;
t556 = -t578 * mrSges(6,2) + t577 * mrSges(6,3);
t558 = t577 * mrSges(4,1) + t578 * mrSges(4,2);
t571 = t656 * t601 + t625 * t612;
t584 = t657 * qJDD(2) - t613 * t602;
t596 = -t657 * qJD(2) + t613 * t649;
t557 = t577 * pkin(3) - t578 * qJ(4);
t595 = t596 ^ 2;
t526 = -t584 * pkin(3) - t595 * qJ(4) + t578 * t557 + qJDD(4) - t529;
t559 = -t577 * mrSges(5,2) - t578 * mrSges(5,3);
t654 = t577 * t596;
t661 = 2 * qJD(5);
t520 = t596 * t661 + (t577 * t578 - t584) * qJ(5) + (t571 - t654) * pkin(4) + t526;
t565 = -t577 * mrSges(6,1) - t596 * mrSges(6,2);
t629 = m(6) * t520 - t584 * mrSges(6,3) + t596 * t565;
t622 = -m(5) * t526 - t571 * mrSges(5,1) - t578 * t559 - t629;
t564 = t577 * mrSges(5,1) + t596 * mrSges(5,3);
t651 = t596 * mrSges(4,2) - t577 * mrSges(4,3) - t564;
t658 = -mrSges(6,1) - mrSges(4,3);
t659 = mrSges(4,1) - mrSges(5,2);
t513 = m(4) * t529 - t651 * t596 + t659 * t584 + (-t556 - t558) * t578 + t658 * t571 + t622;
t574 = t577 * t667;
t641 = t613 * t612 * t539 + t540 * t635 + t656 * t541;
t530 = t574 + t641;
t561 = -t596 * mrSges(4,1) - t578 * mrSges(4,3);
t570 = -qJDD(2) * t634 + t612 * t601 - t602 * t628;
t623 = t595 * pkin(3) - t584 * qJ(4) - t641;
t525 = 0.2e1 * qJD(4) * t596 + ((2 * qJD(3)) + t557) * t577 + t623;
t566 = t578 * mrSges(5,1) - t596 * mrSges(5,2);
t562 = t578 * pkin(4) + t596 * qJ(5);
t576 = t577 ^ 2;
t662 = -0.2e1 * qJD(4);
t522 = -t570 * pkin(4) - t576 * qJ(5) - t577 * t557 + qJDD(5) + t574 + (t662 - t562) * t596 - t623;
t563 = t578 * mrSges(6,1) + t596 * mrSges(6,3);
t637 = -m(6) * t522 - t584 * mrSges(6,2) + t596 * t563;
t624 = -m(5) * t525 + t584 * mrSges(5,3) - t596 * t566 - t637;
t652 = -t556 - t559;
t516 = m(4) * t530 - t584 * mrSges(4,2) + t596 * t561 + (-t558 + t652) * t577 + (-mrSges(5,1) + t658) * t570 + t624;
t531 = t657 * t539 - t613 * t540 + qJDD(3);
t620 = (-t571 - t654) * qJ(4) + t531 + (-t596 * pkin(3) + t662) * t578;
t528 = t570 * pkin(3) + t620;
t524 = -t576 * pkin(4) + t577 * t661 - t578 * t562 + (pkin(3) + qJ(5)) * t570 + t620;
t642 = m(6) * t524 + t570 * mrSges(6,3) + t577 * t565;
t627 = m(5) * t528 - t571 * mrSges(5,3) - t578 * t566 + t642;
t517 = m(4) * t531 + (t561 - t563) * t578 + t651 * t577 + (mrSges(4,2) - mrSges(6,2)) * t571 + t659 * t570 + t627;
t505 = (t656 * t513 + t516 * t612) * t613 + t657 * t517;
t506 = t513 * t628 + t516 * t635 - t613 * t517;
t579 = -t614 * t592 - t660;
t600 = (-mrSges(3,1) * t616 + mrSges(3,2) * t614) * qJD(1);
t605 = -qJD(2) * mrSges(3,2) + mrSges(3,3) * t649;
t504 = m(3) * t579 + qJDD(2) * mrSges(3,1) - t601 * mrSges(3,3) + qJD(2) * t605 - t600 * t650 + t506;
t511 = -t612 * t513 + t656 * t516;
t604 = qJD(2) * mrSges(3,1) - mrSges(3,3) * t650;
t510 = m(3) * t580 - qJDD(2) * mrSges(3,2) + t602 * mrSges(3,3) - qJD(2) * t604 + t600 * t649 + t511;
t630 = -t614 * t504 + t616 * t510;
t497 = m(2) * t607 - t618 * mrSges(2,1) - qJDD(1) * mrSges(2,2) + t630;
t619 = -m(3) * t591 + t602 * mrSges(3,1) - t601 * mrSges(3,2) - t604 * t650 + t605 * t649 - t505;
t502 = m(2) * t606 + qJDD(1) * mrSges(2,1) - t618 * mrSges(2,2) + t619;
t653 = t615 * t497 + t617 * t502;
t498 = t616 * t504 + t614 * t510;
t640 = t643 * t577 - t644 * t578 + t663 * t596;
t639 = t664 * t577 + t665 * t578 - t643 * t596;
t638 = t665 * t577 - t666 * t578 + t644 * t596;
t631 = t617 * t497 - t615 * t502;
t590 = Ifges(3,5) * qJD(2) + (Ifges(3,1) * t614 + Ifges(3,4) * t616) * qJD(1);
t589 = Ifges(3,6) * qJD(2) + (Ifges(3,4) * t614 + Ifges(3,2) * t616) * qJD(1);
t588 = Ifges(3,3) * qJD(2) + (Ifges(3,5) * t614 + Ifges(3,6) * t616) * qJD(1);
t519 = t571 * mrSges(6,1) + t578 * t556 + t629;
t518 = -t570 * mrSges(5,2) - t571 * mrSges(6,2) - t578 * t563 - t577 * t564 + t627;
t507 = mrSges(5,1) * t526 + mrSges(6,1) * t520 + mrSges(4,2) * t531 - mrSges(6,2) * t524 - mrSges(4,3) * t529 - mrSges(5,3) * t528 + pkin(4) * t519 - qJ(4) * t518 - t570 * t665 + t666 * t571 + t640 * t577 + t644 * t584 + t639 * t596;
t500 = -mrSges(4,1) * t531 + mrSges(4,3) * t530 - mrSges(5,1) * t525 + mrSges(5,2) * t528 + mrSges(6,1) * t522 - mrSges(6,3) * t524 - pkin(4) * (t577 * t556 + t637) - qJ(5) * t642 - pkin(3) * t518 + t638 * t596 + t643 * t584 + (qJ(5) * t563 + t640) * t578 + (qJ(5) * mrSges(6,2) + t665) * t571 + (-pkin(4) * mrSges(6,1) + t664) * t570;
t499 = mrSges(4,1) * t529 - mrSges(4,2) * t530 + mrSges(5,2) * t526 - mrSges(5,3) * t525 + mrSges(6,2) * t522 - mrSges(6,3) * t520 - qJ(5) * t519 + pkin(3) * (t596 * t564 + t622) + qJ(4) * t624 + (-pkin(3) * mrSges(5,2) + t663) * t584 + (-pkin(3) * t556 + t639) * t578 + (-pkin(3) * mrSges(6,1) + t644) * t571 + (qJ(4) * t652 - t638) * t577 + (qJ(4) * (-mrSges(5,1) - mrSges(6,1)) - t643) * t570;
t494 = t588 * t649 + mrSges(3,2) * t591 - mrSges(3,3) * t579 + t656 * t507 + Ifges(3,1) * t601 + Ifges(3,4) * t602 + Ifges(3,5) * qJDD(2) - qJD(2) * t589 - t612 * t500 + (-t505 * t613 - t657 * t506) * qJ(3);
t493 = -mrSges(3,1) * t591 + mrSges(3,3) * t580 + Ifges(3,4) * t601 + Ifges(3,2) * t602 + Ifges(3,6) * qJDD(2) - pkin(2) * t505 + qJD(2) * t590 - t613 * t499 + t500 * t628 + t507 * t635 + t511 * t636 - t588 * t650;
t492 = mrSges(2,1) * g(3) - mrSges(3,1) * t579 + mrSges(3,2) * t580 + mrSges(2,3) * t607 - t657 * t499 + t618 * Ifges(2,5) - Ifges(3,5) * t601 + Ifges(2,6) * qJDD(1) - Ifges(3,6) * t602 - Ifges(3,3) * qJDD(2) - pkin(1) * t498 - pkin(2) * t506 + (-t614 * t589 + t616 * t590) * qJD(1) + (-qJ(3) * t511 - t656 * t500 - t612 * t507) * t613;
t491 = -mrSges(2,2) * g(3) - mrSges(2,3) * t606 + Ifges(2,5) * qJDD(1) - t618 * Ifges(2,6) - pkin(7) * t498 - t614 * t493 + t616 * t494;
t1 = [-m(1) * g(1) + t631; -m(1) * g(2) + t653; (-m(1) - m(2)) * g(3) + t498; -mrSges(1,2) * g(3) + mrSges(1,3) * g(2) - pkin(6) * t653 + t617 * t491 - t615 * t492; mrSges(1,1) * g(3) - mrSges(1,3) * g(1) + pkin(6) * t631 + t615 * t491 + t617 * t492; -mrSges(1,1) * g(2) + mrSges(2,1) * t606 + mrSges(1,2) * g(1) - mrSges(2,2) * t607 + Ifges(2,3) * qJDD(1) + pkin(1) * t619 + pkin(7) * t630 + t616 * t493 + t614 * t494;];
tauB = t1;
