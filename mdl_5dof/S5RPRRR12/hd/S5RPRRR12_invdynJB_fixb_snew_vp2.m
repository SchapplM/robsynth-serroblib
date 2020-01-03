% Calculate vector of inverse dynamics joint torques and base forces with Newton-Euler
% S5RPRRR12
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
%   pkin=[a2,a3,a4,a5,d1,d3,d4,d5]';
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
% tauJB [(6+5)x1]
%   joint torques and base forces of inverse dynamics (contains inertial, gravitational coriolis and centrifugal forces)

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-31 19:13
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauJB = S5RPRRR12_invdynJB_fixb_snew_vp2(qJ, qJD, qJDD, g, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(5,1),zeros(3,1),zeros(8,1),zeros(6,1),zeros(6,3),zeros(6,6)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RPRRR12_invdynJB_fixb_snew_vp2: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RPRRR12_invdynJB_fixb_snew_vp2: qJD has to be [5x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [5 1]), ...
  'S5RPRRR12_invdynJB_fixb_snew_vp2: qJDD has to be [5x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S5RPRRR12_invdynJB_fixb_snew_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S5RPRRR12_invdynJB_fixb_snew_vp2: pkin has to be [8x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5RPRRR12_invdynJB_fixb_snew_vp2: m has to be [6x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [6,3]), ...
  'S5RPRRR12_invdynJB_fixb_snew_vp2: mrSges has to be [6x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [6 6]), ...
  'S5RPRRR12_invdynJB_fixb_snew_vp2: Ifges has to be [6x6] (double)'); 

%% Symbolic Calculation
% From invdyn_fixb_NewtonEuler_linkframe_tauJB_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 19:12:44
% EndTime: 2019-12-31 19:12:47
% DurationCPUTime: 2.91s
% Computational Cost: add. (31954->265), mult. (62561->326), div. (0->0), fcn. (38961->8), ass. (0->109)
t655 = sin(qJ(1));
t659 = cos(qJ(1));
t636 = -t659 * g(1) - t655 * g(2);
t671 = qJDD(1) * qJ(2) + (2 * qJD(2) * qJD(1)) + t636;
t653 = sin(qJ(4));
t654 = sin(qJ(3));
t657 = cos(qJ(4));
t658 = cos(qJ(3));
t620 = (t653 * t658 + t654 * t657) * qJD(1);
t687 = -pkin(1) - pkin(6);
t686 = mrSges(2,1) - mrSges(3,2);
t685 = Ifges(2,5) - Ifges(3,4);
t684 = (-Ifges(2,6) + Ifges(3,5));
t635 = t655 * g(1) - t659 * g(2);
t660 = qJD(1) ^ 2;
t670 = -t660 * qJ(2) + qJDD(2) - t635;
t610 = qJDD(1) * t687 + t670;
t600 = t654 * g(3) + t658 * t610;
t680 = qJD(1) * qJD(3);
t678 = t654 * t680;
t630 = qJDD(1) * t658 - t678;
t579 = (-t630 - t678) * pkin(7) + (-t654 * t658 * t660 + qJDD(3)) * pkin(3) + t600;
t601 = -g(3) * t658 + t654 * t610;
t629 = -qJDD(1) * t654 - t658 * t680;
t681 = qJD(1) * t658;
t634 = qJD(3) * pkin(3) - pkin(7) * t681;
t649 = t654 ^ 2;
t580 = -pkin(3) * t649 * t660 + pkin(7) * t629 - qJD(3) * t634 + t601;
t569 = t579 * t653 + t580 * t657;
t621 = (-t653 * t654 + t657 * t658) * qJD(1);
t589 = -qJD(4) * t621 + t629 * t657 - t630 * t653;
t597 = mrSges(5,1) * t620 + mrSges(5,2) * t621;
t643 = qJD(3) + qJD(4);
t608 = mrSges(5,1) * t643 - mrSges(5,3) * t621;
t642 = qJDD(3) + qJDD(4);
t584 = -t629 * pkin(3) + t634 * t681 + (-pkin(7) * t649 + t687) * t660 + t671;
t590 = -qJD(4) * t620 + t629 * t653 + t630 * t657;
t564 = (t620 * t643 - t590) * pkin(8) + (t621 * t643 - t589) * pkin(4) + t584;
t598 = pkin(4) * t620 - pkin(8) * t621;
t641 = t643 ^ 2;
t566 = -pkin(4) * t641 + pkin(8) * t642 - t598 * t620 + t569;
t652 = sin(qJ(5));
t656 = cos(qJ(5));
t562 = t564 * t656 - t566 * t652;
t602 = -t621 * t652 + t643 * t656;
t572 = qJD(5) * t602 + t590 * t656 + t642 * t652;
t603 = t621 * t656 + t643 * t652;
t581 = -mrSges(6,1) * t602 + mrSges(6,2) * t603;
t588 = qJDD(5) - t589;
t616 = qJD(5) + t620;
t591 = -mrSges(6,2) * t616 + mrSges(6,3) * t602;
t559 = m(6) * t562 + t588 * mrSges(6,1) - t572 * mrSges(6,3) - t581 * t603 + t591 * t616;
t563 = t564 * t652 + t566 * t656;
t571 = -qJD(5) * t603 - t590 * t652 + t642 * t656;
t592 = mrSges(6,1) * t616 - mrSges(6,3) * t603;
t560 = m(6) * t563 - t588 * mrSges(6,2) + t571 * mrSges(6,3) + t581 * t602 - t592 * t616;
t674 = -t559 * t652 + t560 * t656;
t547 = m(5) * t569 - mrSges(5,2) * t642 + t589 * mrSges(5,3) - t597 * t620 - t608 * t643 + t674;
t568 = t579 * t657 - t580 * t653;
t607 = -mrSges(5,2) * t643 - mrSges(5,3) * t620;
t565 = -pkin(4) * t642 - pkin(8) * t641 + t598 * t621 - t568;
t668 = -m(6) * t565 + t571 * mrSges(6,1) - t572 * mrSges(6,2) + t591 * t602 - t592 * t603;
t555 = m(5) * t568 + mrSges(5,1) * t642 - t590 * mrSges(5,3) - t597 * t621 + t607 * t643 + t668;
t540 = t547 * t653 + t555 * t657;
t628 = (mrSges(4,1) * t654 + mrSges(4,2) * t658) * qJD(1);
t682 = qJD(1) * t654;
t632 = -qJD(3) * mrSges(4,2) - mrSges(4,3) * t682;
t537 = m(4) * t600 + qJDD(3) * mrSges(4,1) - mrSges(4,3) * t630 + qJD(3) * t632 - t628 * t681 + t540;
t633 = qJD(3) * mrSges(4,1) - mrSges(4,3) * t681;
t675 = t547 * t657 - t555 * t653;
t538 = m(4) * t601 - qJDD(3) * mrSges(4,2) + mrSges(4,3) * t629 - qJD(3) * t633 - t628 * t682 + t675;
t533 = t658 * t537 + t654 * t538;
t615 = -qJDD(1) * pkin(1) + t670;
t669 = -m(3) * t615 + (t660 * mrSges(3,3)) - t533;
t529 = m(2) * t635 - (t660 * mrSges(2,2)) + qJDD(1) * t686 + t669;
t613 = t660 * pkin(1) - t671;
t609 = t660 * t687 + t671;
t549 = t559 * t656 + t560 * t652;
t667 = m(5) * t584 - t589 * mrSges(5,1) + t590 * mrSges(5,2) + t607 * t620 + t621 * t608 + t549;
t665 = -m(4) * t609 + mrSges(4,1) * t629 - t630 * mrSges(4,2) - t632 * t682 - t633 * t681 - t667;
t662 = -m(3) * t613 + (t660 * mrSges(3,2)) + qJDD(1) * mrSges(3,3) - t665;
t543 = m(2) * t636 - (mrSges(2,1) * t660) - qJDD(1) * mrSges(2,2) + t662;
t683 = t529 * t659 + t543 * t655;
t677 = -t529 * t655 + t543 * t659;
t676 = -t654 * t537 + t538 * t658;
t573 = Ifges(6,5) * t603 + Ifges(6,6) * t602 + Ifges(6,3) * t616;
t575 = Ifges(6,1) * t603 + Ifges(6,4) * t602 + Ifges(6,5) * t616;
t552 = -mrSges(6,1) * t565 + mrSges(6,3) * t563 + Ifges(6,4) * t572 + Ifges(6,2) * t571 + Ifges(6,6) * t588 - t573 * t603 + t575 * t616;
t574 = Ifges(6,4) * t603 + Ifges(6,2) * t602 + Ifges(6,6) * t616;
t553 = mrSges(6,2) * t565 - mrSges(6,3) * t562 + Ifges(6,1) * t572 + Ifges(6,4) * t571 + Ifges(6,5) * t588 + t573 * t602 - t574 * t616;
t594 = Ifges(5,4) * t621 - Ifges(5,2) * t620 + Ifges(5,6) * t643;
t595 = Ifges(5,1) * t621 - Ifges(5,4) * t620 + Ifges(5,5) * t643;
t666 = mrSges(5,1) * t568 - mrSges(5,2) * t569 + Ifges(5,5) * t590 + Ifges(5,6) * t589 + Ifges(5,3) * t642 + pkin(4) * t668 + pkin(8) * t674 + t552 * t656 + t553 * t652 + t594 * t621 + t620 * t595;
t593 = Ifges(5,5) * t621 - Ifges(5,6) * t620 + Ifges(5,3) * t643;
t534 = mrSges(5,2) * t584 - mrSges(5,3) * t568 + Ifges(5,1) * t590 + Ifges(5,4) * t589 + Ifges(5,5) * t642 - pkin(8) * t549 - t552 * t652 + t553 * t656 - t593 * t620 - t594 * t643;
t663 = mrSges(6,1) * t562 - mrSges(6,2) * t563 + Ifges(6,5) * t572 + Ifges(6,6) * t571 + Ifges(6,3) * t588 + t574 * t603 - t575 * t602;
t535 = -mrSges(5,1) * t584 + mrSges(5,3) * t569 + Ifges(5,4) * t590 + Ifges(5,2) * t589 + Ifges(5,6) * t642 - pkin(4) * t549 - t593 * t621 + t595 * t643 - t663;
t617 = (Ifges(4,3) * qJD(3)) + (Ifges(4,5) * t658 - Ifges(4,6) * t654) * qJD(1);
t619 = Ifges(4,5) * qJD(3) + (Ifges(4,1) * t658 - Ifges(4,4) * t654) * qJD(1);
t525 = -mrSges(4,1) * t609 + mrSges(4,3) * t601 + Ifges(4,4) * t630 + Ifges(4,2) * t629 + Ifges(4,6) * qJDD(3) - pkin(3) * t667 + pkin(7) * t675 + qJD(3) * t619 + t653 * t534 + t657 * t535 - t617 * t681;
t618 = Ifges(4,6) * qJD(3) + (Ifges(4,4) * t658 - Ifges(4,2) * t654) * qJD(1);
t527 = mrSges(4,2) * t609 - mrSges(4,3) * t600 + Ifges(4,1) * t630 + Ifges(4,4) * t629 + Ifges(4,5) * qJDD(3) - pkin(7) * t540 - qJD(3) * t618 + t534 * t657 - t535 * t653 - t617 * t682;
t531 = qJDD(1) * mrSges(3,2) - t669;
t664 = mrSges(2,1) * t635 - mrSges(2,2) * t636 + mrSges(3,2) * t615 - mrSges(3,3) * t613 - pkin(1) * t531 - pkin(6) * t533 + qJ(2) * t662 - t654 * t525 + t658 * t527 + (Ifges(3,1) + Ifges(2,3)) * qJDD(1);
t661 = mrSges(4,1) * t600 - mrSges(4,2) * t601 + Ifges(4,5) * t630 + Ifges(4,6) * t629 + Ifges(4,3) * qJDD(3) + pkin(3) * t540 + t618 * t681 + t619 * t682 + t666;
t532 = -m(3) * g(3) + t676;
t524 = t661 + t685 * qJDD(1) + (-mrSges(2,2) + mrSges(3,3)) * g(3) + (t684 * t660) - mrSges(2,3) * t635 + mrSges(3,1) * t615 + pkin(2) * t533 - qJ(2) * t532;
t523 = -mrSges(3,1) * t613 + mrSges(2,3) * t636 - pkin(1) * t532 - pkin(2) * t665 - pkin(6) * t676 + g(3) * t686 - qJDD(1) * t684 - t658 * t525 - t654 * t527 + t660 * t685;
t1 = [-m(1) * g(1) + t677; -m(1) * g(2) + t683; (-m(1) - m(2) - m(3)) * g(3) + t676; -mrSges(1,2) * g(3) + mrSges(1,3) * g(2) - pkin(5) * t683 - t523 * t655 + t524 * t659; mrSges(1,1) * g(3) - mrSges(1,3) * g(1) + pkin(5) * t677 + t659 * t523 + t655 * t524; -mrSges(1,1) * g(2) + mrSges(1,2) * g(1) + t664; t664; t531; t661; t666; t663;];
tauJB = t1;
