% Calculate vector of inverse dynamics base forces with Newton-Euler for
% S6RPPRRP1
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
%   pkin=[a2,a3,a4,a5,a6,d1,d4,d5,theta2,theta3]';
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
% Datum: 2019-05-05 14:45
% Revision: 89c353f7eff3bd693eda4e29f35b2761dbc3ada0 (2019-05-03)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauB = S6RPPRRP1_invdynB_fixb_snew_vp2(qJ, qJD, qJDD, g, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(6,1),zeros(3,1),zeros(10,1),zeros(7,1),zeros(7,3),zeros(7,6)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPPRRP1_invdynB_fixb_snew_vp2: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RPPRRP1_invdynB_fixb_snew_vp2: qJD has to be [6x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [6 1]), ...
  'S6RPPRRP1_invdynB_fixb_snew_vp2: qJDD has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6RPPRRP1_invdynB_fixb_snew_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S6RPPRRP1_invdynB_fixb_snew_vp2: pkin has to be [10x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6RPPRRP1_invdynB_fixb_snew_vp2: m has to be [7x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [7,3]), ...
  'S6RPPRRP1_invdynB_fixb_snew_vp2: mrSges has to be [7x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [7 6]), ...
  'S6RPPRRP1_invdynB_fixb_snew_vp2: Ifges has to be [7x6] (double)'); 

%% Symbolic Calculation
% From invdyn_fixb_NewtonEuler_linkframe_tauB_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-05-05 14:43:55
% EndTime: 2019-05-05 14:44:02
% DurationCPUTime: 6.90s
% Computational Cost: add. (78412->300), mult. (173607->358), div. (0->0), fcn. (119029->10), ass. (0->128)
t685 = Ifges(6,1) + Ifges(7,1);
t679 = Ifges(6,4) + Ifges(7,4);
t678 = Ifges(6,5) + Ifges(7,5);
t684 = Ifges(6,2) + Ifges(7,2);
t683 = Ifges(6,6) + Ifges(7,6);
t682 = Ifges(6,3) + Ifges(7,3);
t644 = qJD(1) ^ 2;
t633 = sin(pkin(10));
t635 = cos(pkin(10));
t638 = sin(qJ(4));
t641 = cos(qJ(4));
t649 = t633 * t638 - t635 * t641;
t610 = t649 * qJD(1);
t650 = t633 * t641 + t635 * t638;
t611 = t650 * qJD(1);
t666 = t611 * qJD(4);
t600 = -qJDD(1) * t649 - t666;
t681 = pkin(3) * t635;
t680 = -mrSges(6,2) - mrSges(7,2);
t676 = mrSges(4,2) * t633;
t631 = t635 ^ 2;
t675 = t631 * t644;
t639 = sin(qJ(1));
t642 = cos(qJ(1));
t619 = t639 * g(1) - g(2) * t642;
t617 = qJDD(1) * pkin(1) + t619;
t620 = -g(1) * t642 - g(2) * t639;
t618 = -pkin(1) * t644 + t620;
t634 = sin(pkin(9));
t636 = cos(pkin(9));
t603 = t634 * t617 + t636 * t618;
t593 = -pkin(2) * t644 + qJDD(1) * qJ(3) + t603;
t632 = -g(3) + qJDD(2);
t665 = qJD(1) * qJD(3);
t669 = t635 * t632 - 0.2e1 * t633 * t665;
t572 = (-pkin(7) * qJDD(1) + t644 * t681 - t593) * t633 + t669;
t579 = t633 * t632 + (t593 + 0.2e1 * t665) * t635;
t664 = qJDD(1) * t635;
t575 = -pkin(3) * t675 + pkin(7) * t664 + t579;
t556 = t638 * t572 + t641 * t575;
t596 = mrSges(5,1) * t610 + mrSges(5,2) * t611;
t608 = qJD(4) * mrSges(5,1) - mrSges(5,3) * t611;
t599 = pkin(4) * t610 - pkin(8) * t611;
t643 = qJD(4) ^ 2;
t551 = -pkin(4) * t643 + qJDD(4) * pkin(8) - t599 * t610 + t556;
t630 = t633 ^ 2;
t602 = t636 * t617 - t634 * t618;
t651 = qJDD(3) - t602;
t577 = (-pkin(2) - t681) * qJDD(1) + (-qJ(3) + (-t630 - t631) * pkin(7)) * t644 + t651;
t667 = t610 * qJD(4);
t601 = qJDD(1) * t650 - t667;
t554 = (-t601 + t667) * pkin(8) + (-t600 + t666) * pkin(4) + t577;
t637 = sin(qJ(5));
t640 = cos(qJ(5));
t546 = -t637 * t551 + t640 * t554;
t605 = qJD(4) * t640 - t611 * t637;
t574 = qJD(5) * t605 + qJDD(4) * t637 + t601 * t640;
t606 = qJD(4) * t637 + t611 * t640;
t580 = -mrSges(7,1) * t605 + mrSges(7,2) * t606;
t581 = -mrSges(6,1) * t605 + mrSges(6,2) * t606;
t609 = qJD(5) + t610;
t584 = -mrSges(6,2) * t609 + mrSges(6,3) * t605;
t598 = qJDD(5) - t600;
t543 = -0.2e1 * qJD(6) * t606 + (t605 * t609 - t574) * qJ(6) + (t605 * t606 + t598) * pkin(5) + t546;
t583 = -mrSges(7,2) * t609 + mrSges(7,3) * t605;
t663 = m(7) * t543 + t598 * mrSges(7,1) + t609 * t583;
t535 = m(6) * t546 + t598 * mrSges(6,1) + t609 * t584 + (-t580 - t581) * t606 + (-mrSges(6,3) - mrSges(7,3)) * t574 + t663;
t547 = t640 * t551 + t637 * t554;
t573 = -qJD(5) * t606 + qJDD(4) * t640 - t601 * t637;
t585 = pkin(5) * t609 - qJ(6) * t606;
t604 = t605 ^ 2;
t545 = -pkin(5) * t604 + qJ(6) * t573 + 0.2e1 * qJD(6) * t605 - t585 * t609 + t547;
t662 = m(7) * t545 + t573 * mrSges(7,3) + t605 * t580;
t586 = mrSges(7,1) * t609 - mrSges(7,3) * t606;
t670 = -mrSges(6,1) * t609 + mrSges(6,3) * t606 - t586;
t538 = m(6) * t547 + t573 * mrSges(6,3) + t605 * t581 + t598 * t680 + t609 * t670 + t662;
t656 = -t535 * t637 + t640 * t538;
t531 = m(5) * t556 - qJDD(4) * mrSges(5,2) + mrSges(5,3) * t600 - qJD(4) * t608 - t596 * t610 + t656;
t555 = t572 * t641 - t638 * t575;
t607 = -qJD(4) * mrSges(5,2) - mrSges(5,3) * t610;
t550 = -qJDD(4) * pkin(4) - pkin(8) * t643 + t611 * t599 - t555;
t548 = -pkin(5) * t573 - qJ(6) * t604 + t585 * t606 + qJDD(6) + t550;
t655 = m(7) * t548 - t573 * mrSges(7,1) - t605 * t583;
t646 = -m(6) * t550 + t573 * mrSges(6,1) + t574 * t680 + t605 * t584 + t606 * t670 - t655;
t540 = m(5) * t555 + qJDD(4) * mrSges(5,1) - t601 * mrSges(5,3) + qJD(4) * t607 - t611 * t596 + t646;
t525 = t638 * t531 + t641 * t540;
t578 = -t593 * t633 + t669;
t648 = mrSges(4,3) * qJDD(1) + t644 * (-mrSges(4,1) * t635 + t676);
t523 = m(4) * t578 - t633 * t648 + t525;
t657 = t641 * t531 - t638 * t540;
t524 = m(4) * t579 + t635 * t648 + t657;
t658 = -t523 * t633 + t635 * t524;
t517 = m(3) * t603 - mrSges(3,1) * t644 - qJDD(1) * mrSges(3,2) + t658;
t589 = -qJDD(1) * pkin(2) - t644 * qJ(3) + t651;
t533 = t640 * t535 + t637 * t538;
t647 = m(5) * t577 - t600 * mrSges(5,1) + t601 * mrSges(5,2) + t610 * t607 + t611 * t608 + t533;
t645 = -m(4) * t589 + mrSges(4,1) * t664 - t647 + (t630 * t644 + t675) * mrSges(4,3);
t528 = t645 - t644 * mrSges(3,2) + m(3) * t602 + (mrSges(3,1) - t676) * qJDD(1);
t513 = t634 * t517 + t636 * t528;
t511 = m(2) * t619 + qJDD(1) * mrSges(2,1) - mrSges(2,2) * t644 + t513;
t659 = t636 * t517 - t528 * t634;
t512 = m(2) * t620 - mrSges(2,1) * t644 - qJDD(1) * mrSges(2,2) + t659;
t674 = t642 * t511 + t639 * t512;
t518 = t635 * t523 + t633 * t524;
t673 = t683 * t605 + t678 * t606 + t682 * t609;
t672 = -t684 * t605 - t679 * t606 - t683 * t609;
t671 = t679 * t605 + t685 * t606 + t678 * t609;
t652 = Ifges(4,5) * t633 + Ifges(4,6) * t635;
t668 = t644 * t652;
t661 = m(3) * t632 + t518;
t660 = -t511 * t639 + t642 * t512;
t654 = Ifges(4,1) * t633 + Ifges(4,4) * t635;
t653 = Ifges(4,4) * t633 + Ifges(4,2) * t635;
t592 = Ifges(5,1) * t611 - Ifges(5,4) * t610 + Ifges(5,5) * qJD(4);
t591 = Ifges(5,4) * t611 - Ifges(5,2) * t610 + Ifges(5,6) * qJD(4);
t590 = Ifges(5,5) * t611 - Ifges(5,6) * t610 + Ifges(5,3) * qJD(4);
t541 = -t574 * mrSges(7,3) - t606 * t580 + t663;
t532 = mrSges(6,2) * t550 + mrSges(7,2) * t548 - mrSges(6,3) * t546 - mrSges(7,3) * t543 - qJ(6) * t541 + t679 * t573 + t685 * t574 + t678 * t598 + t673 * t605 + t672 * t609;
t526 = -mrSges(6,1) * t550 + mrSges(6,3) * t547 - mrSges(7,1) * t548 + mrSges(7,3) * t545 - pkin(5) * t655 + qJ(6) * t662 + (-qJ(6) * t586 + t671) * t609 + (-pkin(5) * t586 - t673) * t606 + (-mrSges(7,2) * qJ(6) + t683) * t598 + (-mrSges(7,2) * pkin(5) + t679) * t574 + t684 * t573;
t519 = -mrSges(5,1) * t577 - mrSges(6,1) * t546 - mrSges(7,1) * t543 + mrSges(6,2) * t547 + mrSges(7,2) * t545 + mrSges(5,3) * t556 + Ifges(5,4) * t601 + Ifges(5,2) * t600 + Ifges(5,6) * qJDD(4) - pkin(4) * t533 - pkin(5) * t541 + qJD(4) * t592 - t611 * t590 + t672 * t606 + t671 * t605 - t682 * t598 - t678 * t574 - t683 * t573;
t514 = mrSges(5,2) * t577 - mrSges(5,3) * t555 + Ifges(5,1) * t601 + Ifges(5,4) * t600 + Ifges(5,5) * qJDD(4) - pkin(8) * t533 - qJD(4) * t591 - t526 * t637 + t532 * t640 - t590 * t610;
t507 = mrSges(4,2) * t589 - mrSges(4,3) * t578 - pkin(7) * t525 + qJDD(1) * t654 + t641 * t514 - t638 * t519 + t635 * t668;
t506 = -mrSges(4,1) * t589 + mrSges(4,3) * t579 - pkin(3) * t647 + pkin(7) * t657 + qJDD(1) * t653 + t638 * t514 + t641 * t519 - t633 * t668;
t505 = -pkin(2) * t518 - mrSges(3,1) * t632 + mrSges(3,3) * t603 - pkin(3) * t525 - mrSges(4,1) * t578 + mrSges(4,2) * t579 - t637 * t532 - t640 * t526 - pkin(4) * t646 - pkin(8) * t656 - mrSges(5,1) * t555 + mrSges(5,2) * t556 - Ifges(5,5) * t601 - Ifges(5,6) * t600 - Ifges(5,3) * qJDD(4) - t611 * t591 - t610 * t592 + (Ifges(3,6) - t652) * qJDD(1) + (-t633 * t653 + t635 * t654 + Ifges(3,5)) * t644;
t504 = mrSges(3,2) * t632 - mrSges(3,3) * t602 + Ifges(3,5) * qJDD(1) - Ifges(3,6) * t644 - qJ(3) * t518 - t506 * t633 + t507 * t635;
t503 = -mrSges(2,2) * g(3) - mrSges(2,3) * t619 + Ifges(2,5) * qJDD(1) - Ifges(2,6) * t644 - qJ(2) * t513 + t504 * t636 - t505 * t634;
t502 = mrSges(2,1) * g(3) + mrSges(2,3) * t620 + t644 * Ifges(2,5) + Ifges(2,6) * qJDD(1) - pkin(1) * t661 + qJ(2) * t659 + t634 * t504 + t636 * t505;
t1 = [-m(1) * g(1) + t660; -m(1) * g(2) + t674; (-m(1) - m(2)) * g(3) + t661; -mrSges(1,2) * g(3) + mrSges(1,3) * g(2) - pkin(6) * t674 - t639 * t502 + t642 * t503; mrSges(1,1) * g(3) - mrSges(1,3) * g(1) + pkin(6) * t660 + t642 * t502 + t639 * t503; pkin(1) * t513 + mrSges(2,1) * t619 - mrSges(2,2) * t620 + t635 * t506 + pkin(2) * t645 + qJ(3) * t658 + t633 * t507 + mrSges(3,1) * t602 - mrSges(3,2) * t603 - mrSges(1,1) * g(2) + mrSges(1,2) * g(1) + (-pkin(2) * t676 + Ifges(2,3) + Ifges(3,3)) * qJDD(1);];
tauB  = t1;
