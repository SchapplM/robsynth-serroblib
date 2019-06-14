% Calculate vector of inverse dynamics joint torques and base forces with Newton-Euler
% S6RPPRPR4
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
%   pkin=[a2,a3,a4,a5,a6,d1,d4,d6,theta3,theta5]';
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
% tauJB [(6+6)x1]
%   joint torques and base forces of inverse dynamics (contains inertial, gravitational coriolis and centrifugal forces)

% Quelle: HybrDyn-Toolbox
% Datum: 2019-05-05 14:16
% Revision: 89c353f7eff3bd693eda4e29f35b2761dbc3ada0 (2019-05-03)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauJB = S6RPPRPR4_invdynJB_fixb_snew_vp2(qJ, qJD, qJDD, g, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(6,1),zeros(3,1),zeros(10,1),zeros(7,1),zeros(7,3),zeros(7,6)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPPRPR4_invdynJB_fixb_snew_vp2: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RPPRPR4_invdynJB_fixb_snew_vp2: qJD has to be [6x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [6 1]), ...
  'S6RPPRPR4_invdynJB_fixb_snew_vp2: qJDD has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6RPPRPR4_invdynJB_fixb_snew_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S6RPPRPR4_invdynJB_fixb_snew_vp2: pkin has to be [10x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6RPPRPR4_invdynJB_fixb_snew_vp2: m has to be [7x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [7,3]), ...
  'S6RPPRPR4_invdynJB_fixb_snew_vp2: mrSges has to be [7x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [7 6]), ...
  'S6RPPRPR4_invdynJB_fixb_snew_vp2: Ifges has to be [7x6] (double)'); 

%% Symbolic Calculation
% From invdyn_fixb_NewtonEuler_linkframe_tauJB_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-05-05 14:15:16
% EndTime: 2019-05-05 14:15:22
% DurationCPUTime: 5.20s
% Computational Cost: add. (67605->298), mult. (133771->365), div. (0->0), fcn. (73190->10), ass. (0->122)
t700 = sin(pkin(10));
t702 = cos(pkin(10));
t705 = sin(qJ(4));
t708 = cos(qJ(4));
t670 = (t700 * t705 - t702 * t708) * qJD(1);
t711 = qJD(1) ^ 2;
t706 = sin(qJ(1));
t709 = cos(qJ(1));
t685 = -t709 * g(1) - t706 * g(2);
t716 = qJDD(1) * qJ(2) + 0.2e1 * qJD(2) * qJD(1) + t685;
t734 = -pkin(1) - pkin(2);
t660 = t711 * t734 + t716;
t684 = t706 * g(1) - t709 * g(2);
t715 = -t711 * qJ(2) + qJDD(2) - t684;
t665 = qJDD(1) * t734 + t715;
t701 = sin(pkin(9));
t703 = cos(pkin(9));
t642 = t703 * t660 + t701 * t665;
t638 = -pkin(3) * t711 - qJDD(1) * pkin(7) + t642;
t697 = g(3) + qJDD(3);
t634 = -t705 * t638 + t708 * t697;
t726 = qJD(1) * qJD(4);
t725 = t708 * t726;
t678 = -qJDD(1) * t705 - t725;
t624 = (-t678 - t725) * qJ(5) + (t705 * t708 * t711 + qJDD(4)) * pkin(4) + t634;
t635 = t708 * t638 + t705 * t697;
t679 = -qJDD(1) * t708 + t705 * t726;
t728 = qJD(1) * t705;
t681 = qJD(4) * pkin(4) + qJ(5) * t728;
t696 = t708 ^ 2;
t625 = -pkin(4) * t696 * t711 + qJ(5) * t679 - qJD(4) * t681 + t635;
t735 = 2 * qJD(5);
t620 = t700 * t624 + t702 * t625 + t670 * t735;
t671 = (t700 * t708 + t702 * t705) * qJD(1);
t650 = -pkin(5) * t670 + pkin(8) * t671;
t710 = qJD(4) ^ 2;
t618 = -pkin(5) * t710 + qJDD(4) * pkin(8) + t650 * t670 + t620;
t641 = -t701 * t660 + t703 * t665;
t721 = qJDD(1) * pkin(3) - t641;
t626 = -t681 * t728 - t679 * pkin(4) + qJDD(5) + (-qJ(5) * t696 - pkin(7)) * t711 + t721;
t653 = -t678 * t700 + t679 * t702;
t654 = t678 * t702 + t679 * t700;
t621 = (-qJD(4) * t670 - t654) * pkin(8) + (-qJD(4) * t671 - t653) * pkin(5) + t626;
t704 = sin(qJ(6));
t707 = cos(qJ(6));
t615 = -t618 * t704 + t621 * t707;
t657 = qJD(4) * t707 + t671 * t704;
t633 = qJD(6) * t657 + qJDD(4) * t704 + t654 * t707;
t658 = qJD(4) * t704 - t671 * t707;
t639 = -mrSges(7,1) * t657 + mrSges(7,2) * t658;
t668 = qJD(6) - t670;
t643 = -mrSges(7,2) * t668 + mrSges(7,3) * t657;
t652 = qJDD(6) - t653;
t612 = m(7) * t615 + mrSges(7,1) * t652 - mrSges(7,3) * t633 - t639 * t658 + t643 * t668;
t616 = t618 * t707 + t621 * t704;
t632 = -qJD(6) * t658 + qJDD(4) * t707 - t654 * t704;
t644 = mrSges(7,1) * t668 - mrSges(7,3) * t658;
t613 = m(7) * t616 - mrSges(7,2) * t652 + mrSges(7,3) * t632 + t639 * t657 - t644 * t668;
t604 = -t612 * t704 + t707 * t613;
t649 = -mrSges(6,1) * t670 - mrSges(6,2) * t671;
t662 = qJD(4) * mrSges(6,1) + mrSges(6,3) * t671;
t601 = m(6) * t620 - qJDD(4) * mrSges(6,2) + mrSges(6,3) * t653 - qJD(4) * t662 + t649 * t670 + t604;
t720 = -t702 * t624 + t700 * t625;
t617 = -qJDD(4) * pkin(5) - t710 * pkin(8) + (-(2 * qJD(5)) - t650) * t671 + t720;
t614 = -m(7) * t617 + t632 * mrSges(7,1) - mrSges(7,2) * t633 + t657 * t643 - t644 * t658;
t619 = t671 * t735 - t720;
t661 = -qJD(4) * mrSges(6,2) + mrSges(6,3) * t670;
t608 = m(6) * t619 + qJDD(4) * mrSges(6,1) - mrSges(6,3) * t654 + qJD(4) * t661 + t649 * t671 + t614;
t596 = t700 * t601 + t702 * t608;
t627 = Ifges(7,5) * t658 + Ifges(7,6) * t657 + Ifges(7,3) * t668;
t629 = Ifges(7,1) * t658 + Ifges(7,4) * t657 + Ifges(7,5) * t668;
t605 = -mrSges(7,1) * t617 + mrSges(7,3) * t616 + Ifges(7,4) * t633 + Ifges(7,2) * t632 + Ifges(7,6) * t652 - t627 * t658 + t629 * t668;
t628 = Ifges(7,4) * t658 + Ifges(7,2) * t657 + Ifges(7,6) * t668;
t606 = mrSges(7,2) * t617 - mrSges(7,3) * t615 + Ifges(7,1) * t633 + Ifges(7,4) * t632 + Ifges(7,5) * t652 + t627 * t657 - t628 * t668;
t646 = -Ifges(6,4) * t671 + Ifges(6,2) * t670 + Ifges(6,6) * qJD(4);
t647 = -Ifges(6,1) * t671 + Ifges(6,4) * t670 + Ifges(6,5) * qJD(4);
t673 = Ifges(5,6) * qJD(4) + (-Ifges(5,4) * t705 - Ifges(5,2) * t708) * qJD(1);
t674 = Ifges(5,5) * qJD(4) + (-Ifges(5,1) * t705 - Ifges(5,4) * t708) * qJD(1);
t736 = -(t673 * t705 - t674 * t708) * qJD(1) + (Ifges(5,3) + Ifges(6,3)) * qJDD(4) + mrSges(5,1) * t634 + mrSges(6,1) * t619 - mrSges(5,2) * t635 - mrSges(6,2) * t620 + Ifges(5,5) * t678 + Ifges(6,5) * t654 + Ifges(5,6) * t679 + Ifges(6,6) * t653 + pkin(4) * t596 + pkin(5) * t614 + pkin(8) * t604 + t707 * t605 + t704 * t606 - t671 * t646 - t670 * t647;
t733 = mrSges(2,1) + mrSges(3,1);
t732 = Ifges(3,4) + Ifges(2,5);
t731 = Ifges(2,6) - Ifges(3,6);
t666 = -pkin(1) * t711 + t716;
t677 = (mrSges(5,1) * t708 - mrSges(5,2) * t705) * qJD(1);
t727 = qJD(1) * t708;
t683 = -qJD(4) * mrSges(5,2) - mrSges(5,3) * t727;
t594 = m(5) * t634 + qJDD(4) * mrSges(5,1) - mrSges(5,3) * t678 + qJD(4) * t683 + t677 * t728 + t596;
t682 = qJD(4) * mrSges(5,1) + mrSges(5,3) * t728;
t722 = t702 * t601 - t608 * t700;
t595 = m(5) * t635 - qJDD(4) * mrSges(5,2) + mrSges(5,3) * t679 - qJD(4) * t682 - t677 * t727 + t722;
t590 = -t594 * t705 + t708 * t595;
t586 = m(4) * t642 - mrSges(4,1) * t711 + qJDD(1) * mrSges(4,2) + t590;
t603 = t707 * t612 + t704 * t613;
t602 = m(6) * t626 - t653 * mrSges(6,1) + mrSges(6,2) * t654 - t670 * t661 - t662 * t671 + t603;
t637 = -t711 * pkin(7) + t721;
t598 = -m(5) * t637 + t679 * mrSges(5,1) - mrSges(5,2) * t678 + t682 * t728 - t683 * t727 - t602;
t597 = m(4) * t641 - qJDD(1) * mrSges(4,1) - mrSges(4,2) * t711 + t598;
t723 = t703 * t586 - t701 * t597;
t717 = m(3) * t666 + qJDD(1) * mrSges(3,3) + t723;
t579 = m(2) * t685 - qJDD(1) * mrSges(2,2) - t711 * t733 + t717;
t584 = t586 * t701 + t597 * t703;
t669 = -qJDD(1) * pkin(1) + t715;
t583 = m(3) * t669 - qJDD(1) * mrSges(3,1) - t711 * mrSges(3,3) + t584;
t580 = m(2) * t684 + qJDD(1) * mrSges(2,1) - mrSges(2,2) * t711 - t583;
t729 = t706 * t579 + t709 * t580;
t724 = t709 * t579 - t580 * t706;
t589 = t708 * t594 + t705 * t595;
t588 = m(4) * t697 + t589;
t714 = mrSges(7,1) * t615 - mrSges(7,2) * t616 + Ifges(7,5) * t633 + Ifges(7,6) * t632 + Ifges(7,3) * t652 + t628 * t658 - t629 * t657;
t645 = -Ifges(6,5) * t671 + Ifges(6,6) * t670 + Ifges(6,3) * qJD(4);
t591 = mrSges(6,2) * t626 - mrSges(6,3) * t619 + Ifges(6,1) * t654 + Ifges(6,4) * t653 + Ifges(6,5) * qJDD(4) - pkin(8) * t603 - qJD(4) * t646 - t605 * t704 + t606 * t707 + t645 * t670;
t592 = -mrSges(6,1) * t626 + mrSges(6,3) * t620 + Ifges(6,4) * t654 + Ifges(6,2) * t653 + Ifges(6,6) * qJDD(4) - pkin(5) * t603 + qJD(4) * t647 + t645 * t671 - t714;
t672 = Ifges(5,3) * qJD(4) + (-Ifges(5,5) * t705 - Ifges(5,6) * t708) * qJD(1);
t574 = -mrSges(5,1) * t637 + mrSges(5,3) * t635 + Ifges(5,4) * t678 + Ifges(5,2) * t679 + Ifges(5,6) * qJDD(4) - pkin(4) * t602 + qJ(5) * t722 + qJD(4) * t674 + t700 * t591 + t702 * t592 + t672 * t728;
t575 = mrSges(5,2) * t637 - mrSges(5,3) * t634 + Ifges(5,1) * t678 + Ifges(5,4) * t679 + Ifges(5,5) * qJDD(4) - qJ(5) * t596 - qJD(4) * t673 + t591 * t702 - t592 * t700 - t672 * t727;
t713 = -mrSges(3,1) * t669 - mrSges(4,1) * t641 - mrSges(2,2) * t685 - pkin(2) * t584 - pkin(3) * t598 - pkin(7) * t590 - t574 * t708 - t575 * t705 + qJ(2) * (-mrSges(3,1) * t711 + t717) - pkin(1) * t583 + mrSges(4,2) * t642 + mrSges(3,3) * t666 + mrSges(2,1) * t684 + (Ifges(3,2) + Ifges(2,3) + Ifges(4,3)) * qJDD(1);
t587 = -m(3) * g(3) - t588;
t573 = -mrSges(4,1) * t697 + mrSges(4,3) * t642 + t711 * Ifges(4,5) - Ifges(4,6) * qJDD(1) - pkin(3) * t589 - t736;
t572 = mrSges(4,2) * t697 - mrSges(4,3) * t641 - Ifges(4,5) * qJDD(1) - Ifges(4,6) * t711 - pkin(7) * t589 - t574 * t705 + t575 * t708;
t571 = mrSges(3,2) * t669 - mrSges(2,3) * t684 - qJ(2) * t587 - qJ(3) * t584 + t703 * t572 - t701 * t573 - t731 * t711 + t732 * qJDD(1) + (-mrSges(2,2) + mrSges(3,3)) * g(3);
t570 = mrSges(3,2) * t666 + mrSges(2,3) * t685 - pkin(1) * t587 + pkin(2) * t588 + g(3) * t733 - qJ(3) * t723 + qJDD(1) * t731 - t701 * t572 - t703 * t573 + t711 * t732;
t1 = [-m(1) * g(1) + t724; -m(1) * g(2) + t729; (-m(1) - m(2) - m(3)) * g(3) - t588; -mrSges(1,2) * g(3) + mrSges(1,3) * g(2) - pkin(6) * t729 - t706 * t570 + t709 * t571; mrSges(1,1) * g(3) - mrSges(1,3) * g(1) + pkin(6) * t724 + t709 * t570 + t706 * t571; -mrSges(1,1) * g(2) + mrSges(1,2) * g(1) + t713; t713; t583; t588; t736; t602; t714;];
tauJB  = t1;
