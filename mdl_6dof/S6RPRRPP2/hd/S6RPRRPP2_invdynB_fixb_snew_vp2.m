% Calculate vector of inverse dynamics base forces with Newton-Euler for
% S6RPRRPP2
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
% pkin [9x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d3,d4,theta2]';
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
% Datum: 2019-05-05 21:24
% Revision: 89c353f7eff3bd693eda4e29f35b2761dbc3ada0 (2019-05-03)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauB = S6RPRRPP2_invdynB_fixb_snew_vp2(qJ, qJD, qJDD, g, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(6,1),zeros(3,1),zeros(9,1),zeros(7,1),zeros(7,3),zeros(7,6)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPRRPP2_invdynB_fixb_snew_vp2: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RPRRPP2_invdynB_fixb_snew_vp2: qJD has to be [6x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [6 1]), ...
  'S6RPRRPP2_invdynB_fixb_snew_vp2: qJDD has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6RPRRPP2_invdynB_fixb_snew_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S6RPRRPP2_invdynB_fixb_snew_vp2: pkin has to be [9x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6RPRRPP2_invdynB_fixb_snew_vp2: m has to be [7x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [7,3]), ...
  'S6RPRRPP2_invdynB_fixb_snew_vp2: mrSges has to be [7x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [7 6]), ...
  'S6RPRRPP2_invdynB_fixb_snew_vp2: Ifges has to be [7x6] (double)'); 

%% Symbolic Calculation
% From invdyn_fixb_NewtonEuler_linkframe_tauB_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-05-05 21:20:57
% EndTime: 2019-05-05 21:21:03
% DurationCPUTime: 3.89s
% Computational Cost: add. (37469->298), mult. (70645->347), div. (0->0), fcn. (41218->8), ass. (0->114)
t678 = Ifges(5,1) + Ifges(6,1) + Ifges(7,1);
t663 = Ifges(5,4) - Ifges(6,5) - Ifges(7,4);
t677 = Ifges(5,5) + Ifges(6,4) - Ifges(7,5);
t676 = -Ifges(5,2) - Ifges(6,3) - Ifges(7,2);
t661 = Ifges(5,6) - Ifges(6,6) + Ifges(7,6);
t675 = Ifges(7,3) + Ifges(5,3) + Ifges(6,2);
t636 = sin(qJ(4));
t637 = sin(qJ(3));
t666 = qJD(1) * t637;
t672 = cos(qJ(4));
t610 = -qJD(3) * t672 + t636 * t666;
t639 = cos(qJ(3));
t664 = qJD(1) * qJD(3);
t655 = t639 * t664;
t616 = qJDD(1) * t637 + t655;
t580 = -t610 * qJD(4) + t636 * qJDD(3) + t616 * t672;
t638 = sin(qJ(1));
t640 = cos(qJ(1));
t621 = t638 * g(1) - g(2) * t640;
t612 = qJDD(1) * pkin(1) + t621;
t622 = -g(1) * t640 - g(2) * t638;
t642 = qJD(1) ^ 2;
t614 = -pkin(1) * t642 + t622;
t634 = sin(pkin(9));
t635 = cos(pkin(9));
t582 = t634 * t612 + t635 * t614;
t559 = -pkin(2) * t642 + qJDD(1) * pkin(7) + t582;
t633 = -g(3) + qJDD(2);
t550 = -t637 * t559 + t639 * t633;
t615 = (-pkin(3) * t639 - pkin(8) * t637) * qJD(1);
t641 = qJD(3) ^ 2;
t646 = qJDD(3) * pkin(3) + pkin(8) * t641 - t615 * t666 + t550;
t665 = qJD(1) * t639;
t624 = -qJD(4) + t665;
t669 = t610 * t624;
t674 = -(t580 + t669) * qJ(5) - t646;
t673 = -2 * qJD(5);
t671 = -mrSges(5,3) - mrSges(6,2);
t589 = -mrSges(7,2) * t624 + mrSges(7,3) * t610;
t670 = t589 * t624;
t551 = t639 * t559 + t637 * t633;
t613 = (-mrSges(4,1) * t639 + mrSges(4,2) * t637) * qJD(1);
t628 = t637 * t664;
t617 = qJDD(1) * t639 - t628;
t618 = qJD(3) * mrSges(4,1) - mrSges(4,3) * t666;
t581 = t612 * t635 - t634 * t614;
t558 = -qJDD(1) * pkin(2) - pkin(7) * t642 - t581;
t545 = (-t616 - t655) * pkin(8) + (-t617 + t628) * pkin(3) + t558;
t549 = -pkin(3) * t641 + qJDD(3) * pkin(8) + t615 * t665 + t551;
t542 = t545 * t672 - t636 * t549;
t590 = mrSges(5,2) * t624 - mrSges(5,3) * t610;
t609 = -qJDD(4) + t617;
t611 = t636 * qJD(3) + t666 * t672;
t585 = pkin(4) * t610 - qJ(5) * t611;
t623 = t624 ^ 2;
t540 = t609 * pkin(4) - t623 * qJ(5) + t611 * t585 + qJDD(5) - t542;
t595 = -mrSges(6,2) * t610 - mrSges(6,3) * t624;
t533 = -0.2e1 * qJD(6) * t611 + (-t580 + t669) * qJ(6) + (t610 * t611 + t609) * pkin(5) + t540;
t587 = -mrSges(7,1) * t610 + mrSges(7,2) * t611;
t649 = -m(7) * t533 + t580 * mrSges(7,3) + t611 * t587;
t645 = -m(6) * t540 - t609 * mrSges(6,1) - t624 * t595 + t649;
t586 = mrSges(6,1) * t610 - mrSges(6,3) * t611;
t667 = -mrSges(5,1) * t610 - mrSges(5,2) * t611 - t586;
t529 = m(5) * t542 + (-t589 - t590) * t624 + t667 * t611 + (-mrSges(5,1) - mrSges(7,1)) * t609 + t671 * t580 + t645;
t543 = t636 * t545 + t672 * t549;
t579 = qJD(4) * t611 - qJDD(3) * t672 + t616 * t636;
t592 = mrSges(7,1) * t624 - mrSges(7,3) * t611;
t593 = -mrSges(5,1) * t624 - mrSges(5,3) * t611;
t539 = -pkin(4) * t623 - t609 * qJ(5) - t610 * t585 + t624 * t673 + t543;
t594 = mrSges(6,1) * t624 + mrSges(6,2) * t611;
t591 = pkin(5) * t624 - qJ(6) * t611;
t608 = t610 ^ 2;
t535 = -pkin(5) * t608 + qJ(6) * t579 + 0.2e1 * qJD(6) * t610 - t591 * t624 + t539;
t660 = m(7) * t535 + t579 * mrSges(7,3) + t610 * t587;
t647 = m(6) * t539 - t609 * mrSges(6,3) - t624 * t594 + t660;
t530 = m(5) * t543 + (-t592 + t593) * t624 + t667 * t610 + (mrSges(5,2) - mrSges(7,2)) * t609 + t671 * t579 + t647;
t651 = -t529 * t636 + t672 * t530;
t524 = m(4) * t551 - qJDD(3) * mrSges(4,2) + mrSges(4,3) * t617 - qJD(3) * t618 + t613 * t665 + t651;
t619 = -qJD(3) * mrSges(4,2) + mrSges(4,3) * t665;
t541 = t611 * t673 + (-t611 * t624 + t579) * pkin(4) + t674;
t537 = -qJ(6) * t608 + qJDD(6) + (-pkin(4) - pkin(5)) * t579 + (pkin(4) * t624 + (2 * qJD(5)) + t591) * t611 - t674;
t648 = -m(7) * t537 + t579 * mrSges(7,1) - t580 * mrSges(7,2) + t610 * t589 - t611 * t592;
t531 = m(6) * t541 + t579 * mrSges(6,1) - t580 * mrSges(6,3) - t611 * t594 + t610 * t595 + t648;
t643 = m(5) * t646 - t579 * mrSges(5,1) - t580 * mrSges(5,2) - t610 * t590 - t611 * t593 - t531;
t527 = m(4) * t550 + qJDD(3) * mrSges(4,1) - t616 * mrSges(4,3) + qJD(3) * t619 - t613 * t666 + t643;
t652 = t639 * t524 - t527 * t637;
t517 = m(3) * t582 - mrSges(3,1) * t642 - qJDD(1) * mrSges(3,2) + t652;
t525 = t529 * t672 + t636 * t530;
t644 = -m(4) * t558 + t617 * mrSges(4,1) - t616 * mrSges(4,2) - t618 * t666 + t619 * t665 - t525;
t521 = m(3) * t581 + qJDD(1) * mrSges(3,1) - t642 * mrSges(3,2) + t644;
t513 = t634 * t517 + t635 * t521;
t511 = m(2) * t621 + qJDD(1) * mrSges(2,1) - mrSges(2,2) * t642 + t513;
t653 = t635 * t517 - t521 * t634;
t512 = m(2) * t622 - mrSges(2,1) * t642 - qJDD(1) * mrSges(2,2) + t653;
t668 = t640 * t511 + t638 * t512;
t518 = t637 * t524 + t639 * t527;
t659 = t610 * t661 - t611 * t677 + t624 * t675;
t658 = t610 * t676 + t611 * t663 - t624 * t661;
t657 = t663 * t610 - t611 * t678 + t677 * t624;
t656 = m(3) * t633 + t518;
t654 = -t511 * t638 + t640 * t512;
t602 = Ifges(4,5) * qJD(3) + (Ifges(4,1) * t637 + Ifges(4,4) * t639) * qJD(1);
t601 = Ifges(4,6) * qJD(3) + (Ifges(4,4) * t637 + Ifges(4,2) * t639) * qJD(1);
t600 = Ifges(4,3) * qJD(3) + (Ifges(4,5) * t637 + Ifges(4,6) * t639) * qJD(1);
t532 = mrSges(7,1) * t609 - t649 + t670;
t519 = -mrSges(5,2) * t646 + mrSges(6,2) * t540 + mrSges(7,2) * t537 - mrSges(5,3) * t542 - mrSges(6,3) * t541 - mrSges(7,3) * t533 - qJ(5) * t531 - qJ(6) * t532 - t663 * t579 + t580 * t678 - t609 * t677 + t659 * t610 + t658 * t624;
t514 = mrSges(5,1) * t646 + mrSges(5,3) * t543 - mrSges(6,1) * t541 + mrSges(6,2) * t539 + mrSges(7,1) * t537 - mrSges(7,3) * t535 - pkin(5) * t648 - qJ(6) * t660 - pkin(4) * t531 + (qJ(6) * t592 + t657) * t624 + t659 * t611 + (qJ(6) * mrSges(7,2) - t661) * t609 + t663 * t580 + t676 * t579;
t507 = Ifges(4,6) * qJDD(3) - pkin(4) * (t645 - t670) - qJ(5) * (-t592 * t624 + t647) + Ifges(4,4) * t616 + Ifges(4,2) * t617 + qJD(3) * t602 + mrSges(4,3) * t551 - mrSges(4,1) * t558 + mrSges(6,1) * t540 - mrSges(5,1) * t542 + mrSges(5,2) * t543 + pkin(5) * t532 + mrSges(7,1) * t533 - mrSges(7,2) * t535 - mrSges(6,3) * t539 - pkin(3) * t525 - t600 * t666 + (pkin(4) * t586 - t658) * t611 + (qJ(5) * t586 + t657) * t610 + (pkin(4) * mrSges(7,1) + qJ(5) * mrSges(7,2) + t675) * t609 + (pkin(4) * mrSges(6,2) - t677) * t580 + (qJ(5) * mrSges(6,2) + t661) * t579;
t506 = mrSges(4,2) * t558 - mrSges(4,3) * t550 + Ifges(4,1) * t616 + Ifges(4,4) * t617 + Ifges(4,5) * qJDD(3) - pkin(8) * t525 - qJD(3) * t601 - t636 * t514 + t519 * t672 + t600 * t665;
t505 = Ifges(3,6) * qJDD(1) + t642 * Ifges(3,5) - mrSges(3,1) * t633 + mrSges(3,3) * t582 - Ifges(4,5) * t616 - Ifges(4,6) * t617 - Ifges(4,3) * qJDD(3) - mrSges(4,1) * t550 + mrSges(4,2) * t551 - t636 * t519 - t672 * t514 - pkin(3) * t643 - pkin(8) * t651 - pkin(2) * t518 + (-t601 * t637 + t602 * t639) * qJD(1);
t504 = mrSges(3,2) * t633 - mrSges(3,3) * t581 + Ifges(3,5) * qJDD(1) - Ifges(3,6) * t642 - pkin(7) * t518 + t506 * t639 - t507 * t637;
t503 = -mrSges(2,2) * g(3) - mrSges(2,3) * t621 + Ifges(2,5) * qJDD(1) - Ifges(2,6) * t642 - qJ(2) * t513 + t504 * t635 - t505 * t634;
t502 = mrSges(2,1) * g(3) + mrSges(2,3) * t622 + t642 * Ifges(2,5) + Ifges(2,6) * qJDD(1) - pkin(1) * t656 + qJ(2) * t653 + t634 * t504 + t635 * t505;
t1 = [-m(1) * g(1) + t654; -m(1) * g(2) + t668; (-m(1) - m(2)) * g(3) + t656; -mrSges(1,2) * g(3) + mrSges(1,3) * g(2) - pkin(6) * t668 - t638 * t502 + t640 * t503; mrSges(1,1) * g(3) - mrSges(1,3) * g(1) + pkin(6) * t654 + t640 * t502 + t638 * t503; pkin(1) * t513 + mrSges(2,1) * t621 - mrSges(2,2) * t622 + t637 * t506 + t639 * t507 + pkin(2) * t644 + pkin(7) * t652 + mrSges(3,1) * t581 - mrSges(3,2) * t582 - mrSges(1,1) * g(2) + mrSges(1,2) * g(1) + (Ifges(2,3) + Ifges(3,3)) * qJDD(1);];
tauB  = t1;
