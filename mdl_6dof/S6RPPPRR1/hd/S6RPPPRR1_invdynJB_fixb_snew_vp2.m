% Calculate vector of inverse dynamics joint torques and base forces with Newton-Euler
% S6RPPPRR1
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
%   pkin=[a2,a3,a4,a5,a6,d1,d5,d6,theta2]';
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
% Datum: 2019-05-05 13:32
% Revision: 89c353f7eff3bd693eda4e29f35b2761dbc3ada0 (2019-05-03)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauJB = S6RPPPRR1_invdynJB_fixb_snew_vp2(qJ, qJD, qJDD, g, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(6,1),zeros(3,1),zeros(9,1),zeros(7,1),zeros(7,3),zeros(7,6)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPPPRR1_invdynJB_fixb_snew_vp2: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RPPPRR1_invdynJB_fixb_snew_vp2: qJD has to be [6x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [6 1]), ...
  'S6RPPPRR1_invdynJB_fixb_snew_vp2: qJDD has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6RPPPRR1_invdynJB_fixb_snew_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S6RPPPRR1_invdynJB_fixb_snew_vp2: pkin has to be [9x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6RPPPRR1_invdynJB_fixb_snew_vp2: m has to be [7x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [7,3]), ...
  'S6RPPPRR1_invdynJB_fixb_snew_vp2: mrSges has to be [7x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [7 6]), ...
  'S6RPPPRR1_invdynJB_fixb_snew_vp2: Ifges has to be [7x6] (double)'); 

%% Symbolic Calculation
% From invdyn_fixb_NewtonEuler_linkframe_tauJB_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-05-05 13:32:04
% EndTime: 2019-05-05 13:32:06
% DurationCPUTime: 1.84s
% Computational Cost: add. (19129->245), mult. (33622->286), div. (0->0), fcn. (16242->8), ass. (0->101)
t670 = sin(qJ(1));
t673 = cos(qJ(1));
t640 = t670 * g(1) - t673 * g(2);
t631 = qJDD(1) * pkin(1) + t640;
t641 = -t673 * g(1) - t670 * g(2);
t675 = qJD(1) ^ 2;
t633 = -t675 * pkin(1) + t641;
t666 = sin(pkin(9));
t667 = cos(pkin(9));
t608 = t667 * t631 - t666 * t633;
t599 = -qJDD(1) * pkin(2) - t675 * qJ(3) + qJDD(3) - t608;
t595 = -qJDD(1) * qJ(4) - (2 * qJD(4) * qJD(1)) + t599;
t609 = t666 * t631 + t667 * t633;
t704 = qJDD(1) * qJ(3) + (2 * qJD(3) * qJD(1)) + t609;
t703 = mrSges(3,1) - mrSges(4,2);
t702 = t675 * mrSges(5,3);
t661 = -g(3) + qJDD(2);
t669 = sin(qJ(5));
t701 = t669 * t661;
t597 = t675 * pkin(2) - t704;
t596 = qJDD(4) + (-pkin(2) - qJ(4)) * t675 + t704;
t592 = -qJDD(1) * pkin(7) + t596;
t672 = cos(qJ(5));
t588 = t669 * t592 + t672 * t661;
t632 = (mrSges(6,1) * t669 + mrSges(6,2) * t672) * qJD(1);
t697 = qJD(1) * qJD(5);
t690 = t672 * t697;
t635 = -t669 * qJDD(1) - t690;
t699 = qJD(1) * t672;
t639 = qJD(5) * mrSges(6,1) - mrSges(6,3) * t699;
t591 = -t675 * pkin(7) - t595;
t691 = t669 * t697;
t636 = t672 * qJDD(1) - t691;
t583 = (-t636 + t691) * pkin(8) + (-t635 + t690) * pkin(5) + t591;
t634 = (pkin(5) * t669 - pkin(8) * t672) * qJD(1);
t674 = qJD(5) ^ 2;
t698 = t669 * qJD(1);
t585 = -t674 * pkin(5) + qJDD(5) * pkin(8) - t634 * t698 + t588;
t668 = sin(qJ(6));
t671 = cos(qJ(6));
t581 = t671 * t583 - t668 * t585;
t629 = t671 * qJD(5) - t668 * t699;
t606 = t629 * qJD(6) + t668 * qJDD(5) + t671 * t636;
t630 = t668 * qJD(5) + t671 * t699;
t610 = -t629 * mrSges(7,1) + t630 * mrSges(7,2);
t642 = qJD(6) + t698;
t611 = -t642 * mrSges(7,2) + t629 * mrSges(7,3);
t628 = qJDD(6) - t635;
t578 = m(7) * t581 + t628 * mrSges(7,1) - t606 * mrSges(7,3) - t630 * t610 + t642 * t611;
t582 = t668 * t583 + t671 * t585;
t605 = -t630 * qJD(6) + t671 * qJDD(5) - t668 * t636;
t612 = t642 * mrSges(7,1) - t630 * mrSges(7,3);
t579 = m(7) * t582 - t628 * mrSges(7,2) + t605 * mrSges(7,3) + t629 * t610 - t642 * t612;
t686 = -t668 * t578 + t671 * t579;
t566 = m(6) * t588 - qJDD(5) * mrSges(6,2) + t635 * mrSges(6,3) - qJD(5) * t639 - t632 * t698 + t686;
t587 = t672 * t592 - t701;
t638 = -qJD(5) * mrSges(6,2) - mrSges(6,3) * t698;
t584 = -qJDD(5) * pkin(5) - t674 * pkin(8) + t701 + (qJD(1) * t634 - t592) * t672;
t680 = -m(7) * t584 + t605 * mrSges(7,1) - t606 * mrSges(7,2) + t629 * t611 - t630 * t612;
t574 = m(6) * t587 + qJDD(5) * mrSges(6,1) - t636 * mrSges(6,3) + qJD(5) * t638 - t632 * t699 + t680;
t558 = t669 * t566 + t672 * t574;
t685 = m(5) * t596 + qJDD(1) * mrSges(5,2) + t558;
t681 = -m(4) * t597 + t675 * mrSges(4,2) + qJDD(1) * mrSges(4,3) + t685;
t553 = m(3) * t609 - qJDD(1) * mrSges(3,2) + (-mrSges(3,1) - mrSges(5,3)) * t675 + t681;
t568 = t671 * t578 + t668 * t579;
t682 = -m(6) * t591 + t635 * mrSges(6,1) - t636 * mrSges(6,2) - t638 * t698 - t639 * t699 - t568;
t563 = m(5) * t595 - t675 * mrSges(5,2) - qJDD(1) * mrSges(5,3) + t682;
t678 = -m(4) * t599 + t675 * mrSges(4,3) - t563;
t560 = m(3) * t608 - t675 * mrSges(3,2) + qJDD(1) * t703 + t678;
t547 = t666 * t553 + t667 * t560;
t544 = m(2) * t640 + qJDD(1) * mrSges(2,1) - t675 * mrSges(2,2) + t547;
t688 = t667 * t553 - t666 * t560;
t545 = m(2) * t641 - t675 * mrSges(2,1) - qJDD(1) * mrSges(2,2) + t688;
t700 = t673 * t544 + t670 * t545;
t694 = Ifges(3,5) - Ifges(4,4) + Ifges(5,5);
t693 = Ifges(3,6) - Ifges(4,5) - Ifges(5,4);
t689 = -t670 * t544 + t673 * t545;
t687 = t672 * t566 - t669 * t574;
t684 = m(5) * t661 + t687;
t557 = m(4) * t661 + t684;
t556 = m(3) * t661 + t557;
t600 = Ifges(7,5) * t630 + Ifges(7,6) * t629 + Ifges(7,3) * t642;
t602 = Ifges(7,1) * t630 + Ifges(7,4) * t629 + Ifges(7,5) * t642;
t571 = -mrSges(7,1) * t584 + mrSges(7,3) * t582 + Ifges(7,4) * t606 + Ifges(7,2) * t605 + Ifges(7,6) * t628 - t630 * t600 + t642 * t602;
t601 = Ifges(7,4) * t630 + Ifges(7,2) * t629 + Ifges(7,6) * t642;
t572 = mrSges(7,2) * t584 - mrSges(7,3) * t581 + Ifges(7,1) * t606 + Ifges(7,4) * t605 + Ifges(7,5) * t628 + t629 * t600 - t642 * t601;
t620 = (Ifges(6,6) * qJD(5)) + (Ifges(6,4) * t672 - Ifges(6,2) * t669) * qJD(1);
t621 = (Ifges(6,5) * qJD(5)) + (Ifges(6,1) * t672 - Ifges(6,4) * t669) * qJD(1);
t679 = mrSges(6,1) * t587 - mrSges(6,2) * t588 + Ifges(6,5) * t636 + Ifges(6,6) * t635 + Ifges(6,3) * qJDD(5) + pkin(5) * t680 + pkin(8) * t686 + t671 * t571 + t668 * t572 + t620 * t699 + t621 * t698;
t677 = mrSges(7,1) * t581 - mrSges(7,2) * t582 + Ifges(7,5) * t606 + Ifges(7,6) * t605 + Ifges(7,3) * t628 + t630 * t601 - t629 * t602;
t619 = Ifges(6,3) * qJD(5) + (Ifges(6,5) * t672 - Ifges(6,6) * t669) * qJD(1);
t549 = mrSges(6,2) * t591 - mrSges(6,3) * t587 + Ifges(6,1) * t636 + Ifges(6,4) * t635 + Ifges(6,5) * qJDD(5) - pkin(8) * t568 - qJD(5) * t620 - t668 * t571 + t671 * t572 - t619 * t698;
t550 = -mrSges(6,1) * t591 + mrSges(6,3) * t588 + Ifges(6,4) * t636 + Ifges(6,2) * t635 + Ifges(6,6) * qJDD(5) - pkin(5) * t568 + qJD(5) * t621 - t619 * t699 - t677;
t562 = qJDD(1) * mrSges(4,2) - t678;
t676 = -mrSges(2,2) * t641 - mrSges(3,2) * t609 - mrSges(4,3) * t597 - mrSges(5,3) * t595 + pkin(1) * t547 - pkin(7) * t558 - qJ(4) * t563 + t672 * t549 - t669 * t550 + qJ(3) * (t681 - t702) - pkin(2) * t562 + mrSges(5,2) * t596 + mrSges(4,2) * t599 + mrSges(3,1) * t608 + mrSges(2,1) * t640 + (Ifges(2,3) + Ifges(3,3) + Ifges(4,1) + Ifges(5,1)) * qJDD(1);
t555 = t685 - t702;
t540 = t679 + t693 * qJDD(1) + t694 * t675 + (-mrSges(5,3) - t703) * t661 - qJ(4) * t684 + mrSges(3,3) * t609 + mrSges(5,1) * t596 - mrSges(4,1) * t597 - pkin(2) * t557 + pkin(4) * t558 + pkin(3) * t555;
t539 = -qJ(3) * t557 - mrSges(3,3) * t608 + pkin(3) * t563 + mrSges(4,1) * t599 + t672 * t550 + pkin(4) * t682 + pkin(7) * t687 + t669 * t549 + mrSges(5,1) * t595 - t693 * t675 + (mrSges(3,2) - mrSges(5,2) - mrSges(4,3)) * t661 + t694 * qJDD(1);
t538 = -mrSges(2,2) * g(3) - mrSges(2,3) * t640 + Ifges(2,5) * qJDD(1) - t675 * Ifges(2,6) - qJ(2) * t547 + t667 * t539 - t666 * t540;
t537 = mrSges(2,1) * g(3) + mrSges(2,3) * t641 + t675 * Ifges(2,5) + Ifges(2,6) * qJDD(1) - pkin(1) * t556 + qJ(2) * t688 + t666 * t539 + t667 * t540;
t1 = [-m(1) * g(1) + t689; -m(1) * g(2) + t700; (-m(1) - m(2)) * g(3) + t556; -mrSges(1,2) * g(3) + mrSges(1,3) * g(2) - pkin(6) * t700 - t670 * t537 + t673 * t538; mrSges(1,1) * g(3) - mrSges(1,3) * g(1) + pkin(6) * t689 + t673 * t537 + t670 * t538; -mrSges(1,1) * g(2) + mrSges(1,2) * g(1) + t676; t676; t556; t562; t555; t679; t677;];
tauJB  = t1;
