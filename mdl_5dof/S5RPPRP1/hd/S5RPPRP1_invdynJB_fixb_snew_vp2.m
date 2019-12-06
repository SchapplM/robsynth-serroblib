% Calculate vector of inverse dynamics joint torques and base forces with Newton-Euler
% S5RPPRP1
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
%   pkin=[a2,a3,a4,a5,d1,d4,theta2,theta3]';
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
% Datum: 2019-12-05 17:36
% Revision: 77da58f92bca3eff71542919beafa37024070d86 (2019-12-05)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauJB = S5RPPRP1_invdynJB_fixb_snew_vp2(qJ, qJD, qJDD, g, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(5,1),zeros(3,1),zeros(8,1),zeros(6,1),zeros(6,3),zeros(6,6)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RPPRP1_invdynJB_fixb_snew_vp2: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RPPRP1_invdynJB_fixb_snew_vp2: qJD has to be [5x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [5 1]), ...
  'S5RPPRP1_invdynJB_fixb_snew_vp2: qJDD has to be [5x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S5RPPRP1_invdynJB_fixb_snew_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S5RPPRP1_invdynJB_fixb_snew_vp2: pkin has to be [8x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5RPPRP1_invdynJB_fixb_snew_vp2: m has to be [6x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [6,3]), ...
  'S5RPPRP1_invdynJB_fixb_snew_vp2: mrSges has to be [6x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [6 6]), ...
  'S5RPPRP1_invdynJB_fixb_snew_vp2: Ifges has to be [6x6] (double)'); 

%% Symbolic Calculation
% From invdyn_fixb_NewtonEuler_linkframe_tauJB_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-05 17:35:51
% EndTime: 2019-12-05 17:35:53
% DurationCPUTime: 2.21s
% Computational Cost: add. (17596->229), mult. (36515->286), div. (0->0), fcn. (20211->8), ass. (0->107)
t658 = sin(qJ(1));
t660 = cos(qJ(1));
t635 = t660 * g(2) + t658 * g(3);
t630 = qJDD(1) * pkin(1) + t635;
t634 = t658 * g(2) - g(3) * t660;
t661 = qJD(1) ^ 2;
t631 = -pkin(1) * t661 + t634;
t654 = sin(pkin(7));
t656 = cos(pkin(7));
t603 = t654 * t630 + t656 * t631;
t709 = -pkin(2) * t661 + qJDD(1) * qJ(3) + (2 * qJD(1) * qJD(3)) + t603;
t708 = Ifges(5,1) + Ifges(6,1);
t701 = Ifges(5,4) + Ifges(6,4);
t700 = Ifges(5,5) + Ifges(6,5);
t707 = Ifges(5,2) + Ifges(6,2);
t699 = Ifges(5,6) + Ifges(6,6);
t706 = Ifges(5,3) + Ifges(6,3);
t652 = -g(1) + qJDD(2);
t653 = sin(pkin(8));
t655 = cos(pkin(8));
t592 = t652 * t655 - t709 * t653;
t657 = sin(qJ(4));
t659 = cos(qJ(4));
t688 = qJD(1) * t655;
t637 = qJD(4) - t688;
t689 = qJD(1) * t653;
t691 = (-t701 * t657 + t708 * t659) * t689 + t700 * t637;
t692 = (-t707 * t657 + t701 * t659) * t689 + t699 * t637;
t705 = t691 * t657 + t692 * t659;
t619 = (mrSges(6,1) * t657 + mrSges(6,2) * t659) * t689;
t685 = qJD(1) * qJD(4);
t622 = (qJDD(1) * t659 - t657 * t685) * t653;
t680 = t659 * t689;
t593 = t653 * t652 + t709 * t655;
t671 = -pkin(3) * t655 - pkin(6) * t653;
t629 = t671 * qJD(1);
t591 = t629 * t688 + t593;
t602 = t630 * t656 - t654 * t631;
t666 = -qJ(3) * t661 + qJDD(3) - t602;
t596 = (-pkin(2) + t671) * qJDD(1) + t666;
t595 = t659 * t596;
t684 = qJDD(1) * t655;
t636 = qJDD(4) - t684;
t672 = -0.2e1 * qJD(5) * t689;
t695 = t653 ^ 2 * t661;
t583 = t659 * t672 + pkin(4) * t636 - qJ(5) * t622 + t595 + (-pkin(4) * t659 * t695 - qJ(5) * t637 * t689 - t591) * t657;
t681 = t657 * t689;
t614 = -mrSges(6,2) * t637 - mrSges(6,3) * t681;
t682 = m(6) * t583 + t636 * mrSges(6,1) + t637 * t614;
t580 = -mrSges(6,3) * t622 - t619 * t680 + t682;
t587 = t659 * t591 + t657 * t596;
t616 = pkin(4) * t637 - qJ(5) * t680;
t621 = (-qJDD(1) * t657 - t659 * t685) * t653;
t683 = t657 ^ 2 * t695;
t585 = -pkin(4) * t683 + qJ(5) * t621 - t616 * t637 + t657 * t672 + t587;
t586 = -t591 * t657 + t595;
t704 = mrSges(5,1) * t586 + mrSges(6,1) * t583 - mrSges(5,2) * t587 - mrSges(6,2) * t585 + pkin(4) * t580 + t699 * t621 + t700 * t622 + t706 * t636;
t703 = t655 ^ 2;
t702 = -mrSges(5,2) - mrSges(6,2);
t697 = mrSges(4,2) * t653;
t696 = mrSges(6,2) * t622;
t627 = (-mrSges(4,1) * t655 + t697) * qJD(1);
t615 = -mrSges(5,2) * t637 - mrSges(5,3) * t681;
t670 = (-t619 - (mrSges(5,1) * t657 + mrSges(5,2) * t659) * t689) * t689;
t576 = m(5) * t586 + mrSges(5,1) * t636 + t615 * t637 + (-mrSges(5,3) - mrSges(6,3)) * t622 + t659 * t670 + t682;
t617 = mrSges(6,1) * t637 - mrSges(6,3) * t680;
t690 = -mrSges(5,1) * t637 + mrSges(5,3) * t680 - t617;
t694 = m(6) * t585 + t621 * mrSges(6,3);
t577 = m(5) * t587 + mrSges(5,3) * t621 + t702 * t636 + t690 * t637 + t657 * t670 + t694;
t674 = -t576 * t657 + t659 * t577;
t687 = qJDD(1) * mrSges(4,3);
t570 = m(4) * t593 + (qJD(1) * t627 + t687) * t655 + t674;
t665 = t690 * t659 + (-t614 - t615) * t657;
t590 = t629 * t689 - t592;
t588 = -pkin(4) * t621 - qJ(5) * t683 + t616 * t680 + qJDD(5) + t590;
t678 = m(6) * t588 - t621 * mrSges(6,1);
t667 = -m(5) * t590 + t621 * mrSges(5,1) - t678;
t579 = m(4) * t592 + t702 * t622 + (-t687 + (-t627 + t665) * qJD(1)) * t653 + t667;
t675 = t655 * t570 - t579 * t653;
t561 = m(3) * t603 - mrSges(3,1) * t661 - qJDD(1) * mrSges(3,2) + t675;
t574 = t576 * t659 + t577 * t657;
t599 = -qJDD(1) * pkin(2) + t666;
t664 = -m(4) * t599 + mrSges(4,1) * t684 - t574 + (t661 * t703 + t695) * mrSges(4,3);
t567 = m(3) * t602 - mrSges(3,2) * t661 + (mrSges(3,1) - t697) * qJDD(1) + t664;
t556 = t654 * t561 + t656 * t567;
t564 = t653 * t570 + t655 * t579;
t693 = (t699 * t657 - t700 * t659) * t689 - t706 * t637;
t562 = m(3) * t652 + t564;
t676 = t656 * t561 - t654 * t567;
t553 = m(2) * t634 - mrSges(2,1) * t661 - qJDD(1) * mrSges(2,2) + t676;
t554 = m(2) * t635 + qJDD(1) * mrSges(2,1) - mrSges(2,2) * t661 + t556;
t677 = t660 * t553 - t554 * t658;
t669 = Ifges(4,5) * t653 + Ifges(4,6) * t655;
t668 = -t553 * t658 - t554 * t660;
t581 = t696 + (t614 * t657 + t617 * t659) * t689 + t678;
t565 = -mrSges(5,1) * t590 + mrSges(5,3) * t587 - mrSges(6,1) * t588 + mrSges(6,3) * t585 - pkin(4) * t581 + qJ(5) * t694 + (-qJ(5) * t617 + t691) * t637 + (-mrSges(6,2) * qJ(5) + t699) * t636 + t701 * t622 + t707 * t621 + (-qJ(5) * t619 * t657 + t693 * t659) * t689;
t573 = mrSges(5,2) * t590 + mrSges(6,2) * t588 - mrSges(5,3) * t586 - mrSges(6,3) * t583 - qJ(5) * t580 + t701 * t621 + t708 * t622 + t700 * t636 - t692 * t637 + t693 * t681;
t628 = t669 * qJD(1);
t551 = t628 * t688 + mrSges(4,2) * t599 - mrSges(4,3) * t592 - pkin(6) * t574 - t565 * t657 + t573 * t659 + (Ifges(4,1) * t653 + Ifges(4,4) * t655) * qJDD(1);
t558 = Ifges(4,2) * t684 - mrSges(4,1) * t599 + mrSges(4,3) * t593 - pkin(3) * t574 + (Ifges(4,4) * qJDD(1) + (-t628 - t705) * qJD(1)) * t653 - t704;
t572 = qJDD(1) * t697 - t664;
t662 = mrSges(2,1) * t635 + mrSges(3,1) * t602 - mrSges(2,2) * t634 - mrSges(3,2) * t603 + pkin(1) * t556 - pkin(2) * t572 + qJ(3) * t675 + t653 * t551 + t655 * t558 + (Ifges(2,3) + Ifges(3,3)) * qJDD(1);
t549 = t661 * Ifges(3,5) - mrSges(3,1) * t652 + mrSges(3,3) * t603 - mrSges(4,1) * t592 + mrSges(4,2) * t593 - t657 * t573 - t659 * t565 - pkin(3) * (-mrSges(5,2) * t622 + t667 - t696) - pkin(6) * t674 - pkin(2) * t564 + (Ifges(3,6) - t669) * qJDD(1) + (Ifges(4,4) * t703 * qJD(1) + (-pkin(3) * t665 + (-Ifges(4,4) * t653 + (Ifges(4,1) - Ifges(4,2)) * t655) * qJD(1)) * t653) * qJD(1);
t548 = mrSges(3,2) * t652 - mrSges(3,3) * t602 + Ifges(3,5) * qJDD(1) - Ifges(3,6) * t661 - qJ(3) * t564 + t551 * t655 - t558 * t653;
t547 = -mrSges(2,2) * g(1) - mrSges(2,3) * t635 + Ifges(2,5) * qJDD(1) - Ifges(2,6) * t661 - qJ(2) * t556 + t548 * t656 - t549 * t654;
t546 = mrSges(2,1) * g(1) + mrSges(2,3) * t634 + t661 * Ifges(2,5) + Ifges(2,6) * qJDD(1) - pkin(1) * t562 + qJ(2) * t676 + t654 * t548 + t656 * t549;
t1 = [(-m(1) - m(2)) * g(1) + t562; -m(1) * g(2) + t668; -m(1) * g(3) + t677; -mrSges(1,2) * g(3) + mrSges(1,3) * g(2) + t662; mrSges(1,1) * g(3) - mrSges(1,3) * g(1) - pkin(5) * t677 - t660 * t546 - t658 * t547; -mrSges(1,1) * g(2) + mrSges(1,2) * g(1) + pkin(5) * t668 - t658 * t546 + t660 * t547; t662; t562; t572; t705 * t689 + t704; t581;];
tauJB = t1;
