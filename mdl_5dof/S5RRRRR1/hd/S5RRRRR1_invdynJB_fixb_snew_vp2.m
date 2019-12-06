% Calculate vector of inverse dynamics joint torques and base forces with Newton-Euler
% S5RRRRR1
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
% pkin [6x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d1,d5]';
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
% Datum: 2019-12-05 18:52
% Revision: 77da58f92bca3eff71542919beafa37024070d86 (2019-12-05)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauJB = S5RRRRR1_invdynJB_fixb_snew_vp2(qJ, qJD, qJDD, g, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(5,1),zeros(3,1),zeros(6,1),zeros(6,1),zeros(6,3),zeros(6,6)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RRRRR1_invdynJB_fixb_snew_vp2: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RRRRR1_invdynJB_fixb_snew_vp2: qJD has to be [5x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [5 1]), ...
  'S5RRRRR1_invdynJB_fixb_snew_vp2: qJDD has to be [5x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S5RRRRR1_invdynJB_fixb_snew_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [6 1]), ...
  'S5RRRRR1_invdynJB_fixb_snew_vp2: pkin has to be [6x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5RRRRR1_invdynJB_fixb_snew_vp2: m has to be [6x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [6,3]), ...
  'S5RRRRR1_invdynJB_fixb_snew_vp2: mrSges has to be [6x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [6 6]), ...
  'S5RRRRR1_invdynJB_fixb_snew_vp2: Ifges has to be [6x6] (double)'); 

%% Symbolic Calculation
% From invdyn_fixb_NewtonEuler_linkframe_tauJB_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-05 18:50:21
% EndTime: 2019-12-05 18:50:30
% DurationCPUTime: 6.65s
% Computational Cost: add. (71792->294), mult. (157349->376), div. (0->0), fcn. (119548->10), ass. (0->119)
t704 = sin(qJ(1));
t709 = cos(qJ(1));
t686 = -t709 * g(1) - t704 * g(2);
t710 = qJD(1) ^ 2;
t679 = -t710 * pkin(1) + t686;
t703 = sin(qJ(2));
t708 = cos(qJ(2));
t663 = t708 * g(3) - t703 * t679;
t678 = (mrSges(3,1) * t708 - mrSges(3,2) * t703) * qJD(1);
t724 = qJD(1) * qJD(2);
t680 = -t703 * qJDD(1) - t708 * t724;
t725 = qJD(1) * t708;
t684 = -qJD(2) * mrSges(3,2) - mrSges(3,3) * t725;
t726 = qJD(1) * t703;
t657 = (t703 * t708 * t710 + qJDD(2)) * pkin(2) + t663;
t664 = t703 * g(3) + t708 * t679;
t658 = (-t708 ^ 2 * t710 - qJD(2) ^ 2) * pkin(2) + t664;
t702 = sin(qJ(3));
t707 = cos(qJ(3));
t632 = t707 * t657 - t702 * t658;
t669 = (t702 * t703 - t707 * t708) * qJD(1);
t723 = t703 * t724;
t681 = -t708 * qJDD(1) + t723;
t641 = t669 * qJD(3) + t707 * t680 + t702 * t681;
t670 = (-t702 * t708 - t703 * t707) * qJD(1);
t651 = -t669 * mrSges(4,1) + t670 * mrSges(4,2);
t696 = qJD(2) + qJD(3);
t659 = -t696 * mrSges(4,2) + t669 * mrSges(4,3);
t695 = qJDD(2) + qJDD(3);
t617 = (t669 * t670 + t695) * pkin(3) + t632;
t633 = t702 * t657 + t707 * t658;
t624 = (-t669 ^ 2 - t696 ^ 2) * pkin(3) + t633;
t701 = sin(qJ(4));
t706 = cos(qJ(4));
t604 = t701 * t617 + t706 * t624;
t640 = -t670 * qJD(3) - t702 * t680 + t707 * t681;
t650 = t701 * t669 + t706 * t670;
t613 = -t650 * qJD(4) + t706 * t640 - t701 * t641;
t649 = t706 * t669 - t701 * t670;
t629 = -t649 * mrSges(5,1) + t650 * mrSges(5,2);
t691 = qJD(4) + t696;
t643 = t691 * mrSges(5,1) - t650 * mrSges(5,3);
t690 = qJDD(4) + t695;
t614 = t649 * qJD(4) + t701 * t640 + t706 * t641;
t685 = t704 * g(1) - t709 * g(2);
t677 = qJDD(1) * pkin(1) + t685;
t655 = (-t681 - t723) * pkin(2) + t677;
t620 = t655 + (t670 * t696 - t640) * pkin(3);
t596 = (-t649 * t691 - t614) * pkin(6) + (t650 * t691 - t613) * pkin(4) + t620;
t630 = -t649 * pkin(4) - t650 * pkin(6);
t689 = t691 ^ 2;
t598 = -t689 * pkin(4) + t690 * pkin(6) + t649 * t630 + t604;
t700 = sin(qJ(5));
t705 = cos(qJ(5));
t594 = t705 * t596 - t700 * t598;
t635 = -t700 * t650 + t705 * t691;
t601 = t635 * qJD(5) + t705 * t614 + t700 * t690;
t611 = qJDD(5) - t613;
t636 = t705 * t650 + t700 * t691;
t618 = -t635 * mrSges(6,1) + t636 * mrSges(6,2);
t647 = qJD(5) - t649;
t621 = -t647 * mrSges(6,2) + t635 * mrSges(6,3);
t591 = m(6) * t594 + t611 * mrSges(6,1) - t601 * mrSges(6,3) - t636 * t618 + t647 * t621;
t595 = t700 * t596 + t705 * t598;
t600 = -t636 * qJD(5) - t700 * t614 + t705 * t690;
t622 = t647 * mrSges(6,1) - t636 * mrSges(6,3);
t592 = m(6) * t595 - t611 * mrSges(6,2) + t600 * mrSges(6,3) + t635 * t618 - t647 * t622;
t721 = -t700 * t591 + t705 * t592;
t579 = m(5) * t604 - t690 * mrSges(5,2) + t613 * mrSges(5,3) + t649 * t629 - t691 * t643 + t721;
t603 = t706 * t617 - t701 * t624;
t642 = -t691 * mrSges(5,2) + t649 * mrSges(5,3);
t597 = -t690 * pkin(4) - t689 * pkin(6) + t650 * t630 - t603;
t719 = -m(6) * t597 + t600 * mrSges(6,1) - t601 * mrSges(6,2) + t635 * t621 - t636 * t622;
t587 = m(5) * t603 + t690 * mrSges(5,1) - t614 * mrSges(5,3) - t650 * t629 + t691 * t642 + t719;
t727 = t701 * t579 + t706 * t587;
t572 = m(4) * t632 + t695 * mrSges(4,1) - t641 * mrSges(4,3) - t670 * t651 + t696 * t659 + t727;
t660 = t696 * mrSges(4,1) - t670 * mrSges(4,3);
t573 = m(4) * t633 - t695 * mrSges(4,2) + t640 * mrSges(4,3) + t706 * t579 - t701 * t587 + t669 * t651 - t696 * t660;
t728 = t707 * t572 + t702 * t573;
t565 = m(3) * t663 + qJDD(2) * mrSges(3,1) - t680 * mrSges(3,3) + qJD(2) * t684 + t678 * t726 + t728;
t683 = qJD(2) * mrSges(3,1) + mrSges(3,3) * t726;
t566 = m(3) * t664 - qJDD(2) * mrSges(3,2) + t681 * mrSges(3,3) - qJD(2) * t683 - t702 * t572 + t707 * t573 - t678 * t725;
t562 = m(2) * t686 - t710 * mrSges(2,1) - qJDD(1) * mrSges(2,2) - t703 * t565 + t708 * t566;
t581 = t705 * t591 + t700 * t592;
t718 = m(5) * t620 - t613 * mrSges(5,1) + t614 * mrSges(5,2) - t649 * t642 + t650 * t643 + t581;
t714 = m(4) * t655 - t640 * mrSges(4,1) + t641 * mrSges(4,2) - t669 * t659 + t670 * t660 + t718;
t712 = m(3) * t677 - t681 * mrSges(3,1) + t680 * mrSges(3,2) - t683 * t726 + t684 * t725 + t714;
t576 = m(2) * t685 + qJDD(1) * mrSges(2,1) - t710 * mrSges(2,2) + t712;
t729 = t704 * t562 + t709 * t576;
t722 = t709 * t562 - t704 * t576;
t720 = -t708 * t565 - t703 * t566;
t605 = Ifges(6,5) * t636 + Ifges(6,6) * t635 + Ifges(6,3) * t647;
t607 = Ifges(6,1) * t636 + Ifges(6,4) * t635 + Ifges(6,5) * t647;
t584 = -mrSges(6,1) * t597 + mrSges(6,3) * t595 + Ifges(6,4) * t601 + Ifges(6,2) * t600 + Ifges(6,6) * t611 - t636 * t605 + t647 * t607;
t606 = Ifges(6,4) * t636 + Ifges(6,2) * t635 + Ifges(6,6) * t647;
t585 = mrSges(6,2) * t597 - mrSges(6,3) * t594 + Ifges(6,1) * t601 + Ifges(6,4) * t600 + Ifges(6,5) * t611 + t635 * t605 - t647 * t606;
t625 = Ifges(5,5) * t650 + Ifges(5,6) * t649 + Ifges(5,3) * t691;
t626 = Ifges(5,4) * t650 + Ifges(5,2) * t649 + Ifges(5,6) * t691;
t568 = mrSges(5,2) * t620 - mrSges(5,3) * t603 + Ifges(5,1) * t614 + Ifges(5,4) * t613 + Ifges(5,5) * t690 - pkin(6) * t581 - t700 * t584 + t705 * t585 + t649 * t625 - t691 * t626;
t627 = Ifges(5,1) * t650 + Ifges(5,4) * t649 + Ifges(5,5) * t691;
t715 = mrSges(6,1) * t594 - mrSges(6,2) * t595 + Ifges(6,5) * t601 + Ifges(6,6) * t600 + Ifges(6,3) * t611 + t636 * t606 - t635 * t607;
t569 = -mrSges(5,1) * t620 + mrSges(5,3) * t604 + Ifges(5,4) * t614 + Ifges(5,2) * t613 + Ifges(5,6) * t690 - pkin(4) * t581 - t650 * t625 + t691 * t627 - t715;
t644 = Ifges(4,5) * t670 + Ifges(4,6) * t669 + Ifges(4,3) * t696;
t646 = Ifges(4,1) * t670 + Ifges(4,4) * t669 + Ifges(4,5) * t696;
t563 = -mrSges(4,1) * t655 + mrSges(4,3) * t633 + Ifges(4,4) * t641 + Ifges(4,2) * t640 + Ifges(4,6) * t695 - pkin(3) * t718 + t701 * t568 + t706 * t569 - t670 * t644 + t696 * t646;
t645 = Ifges(4,4) * t670 + Ifges(4,2) * t669 + Ifges(4,6) * t696;
t564 = mrSges(4,2) * t655 - mrSges(4,3) * t632 + Ifges(4,1) * t641 + Ifges(4,4) * t640 + Ifges(4,5) * t695 + t706 * t568 - t701 * t569 + t669 * t644 - t696 * t645;
t666 = Ifges(3,3) * qJD(2) + (-Ifges(3,5) * t703 - Ifges(3,6) * t708) * qJD(1);
t668 = Ifges(3,5) * qJD(2) + (-Ifges(3,1) * t703 - Ifges(3,4) * t708) * qJD(1);
t558 = -mrSges(3,1) * t677 + mrSges(3,3) * t664 + Ifges(3,4) * t680 + Ifges(3,2) * t681 + Ifges(3,6) * qJDD(2) - pkin(2) * t714 + qJD(2) * t668 + t707 * t563 + t702 * t564 + t666 * t726;
t667 = Ifges(3,6) * qJD(2) + (-Ifges(3,4) * t703 - Ifges(3,2) * t708) * qJD(1);
t559 = mrSges(3,2) * t677 - mrSges(3,3) * t663 + Ifges(3,1) * t680 + Ifges(3,4) * t681 + Ifges(3,5) * qJDD(2) - qJD(2) * t667 - t702 * t563 + t707 * t564 - t666 * t725;
t717 = mrSges(2,1) * t685 - mrSges(2,2) * t686 + Ifges(2,3) * qJDD(1) + pkin(1) * t712 - t708 * t558 - t703 * t559;
t716 = mrSges(5,1) * t603 - mrSges(5,2) * t604 + Ifges(5,5) * t614 + Ifges(5,6) * t613 + Ifges(5,3) * t690 + pkin(4) * t719 + pkin(6) * t721 + t705 * t584 + t700 * t585 + t650 * t626 - t649 * t627;
t713 = mrSges(4,1) * t632 - mrSges(4,2) * t633 + Ifges(4,5) * t641 + Ifges(4,6) * t640 + Ifges(4,3) * t695 + pkin(3) * t727 + t670 * t645 - t669 * t646 + t716;
t711 = mrSges(3,1) * t663 - mrSges(3,2) * t664 + Ifges(3,5) * t680 + Ifges(3,6) * t681 + Ifges(3,3) * qJDD(2) + pkin(2) * t728 - t667 * t726 + t668 * t725 + t713;
t557 = mrSges(2,1) * g(3) + mrSges(2,3) * t686 + t710 * Ifges(2,5) + Ifges(2,6) * qJDD(1) - pkin(1) * t720 + t711;
t556 = -mrSges(2,2) * g(3) - mrSges(2,3) * t685 + Ifges(2,5) * qJDD(1) - t710 * Ifges(2,6) - t703 * t558 + t708 * t559;
t1 = [-m(1) * g(1) + t722; -m(1) * g(2) + t729; (-m(1) - m(2)) * g(3) + t720; -mrSges(1,2) * g(3) + mrSges(1,3) * g(2) - pkin(5) * t729 + t709 * t556 - t704 * t557; mrSges(1,1) * g(3) - mrSges(1,3) * g(1) + pkin(5) * t722 + t704 * t556 + t709 * t557; -mrSges(1,1) * g(2) + mrSges(1,2) * g(1) + t717; t717; t711; t713; t716; t715;];
tauJB = t1;
