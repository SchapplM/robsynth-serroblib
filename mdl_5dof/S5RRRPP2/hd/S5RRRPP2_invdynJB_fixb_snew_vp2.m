% Calculate vector of inverse dynamics joint torques and base forces with Newton-Euler
% S5RRRPP2
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
% pkin [7x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d1,d2,d3]';
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
% Datum: 2019-12-31 20:52
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauJB = S5RRRPP2_invdynJB_fixb_snew_vp2(qJ, qJD, qJDD, g, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(5,1),zeros(3,1),zeros(7,1),zeros(6,1),zeros(6,3),zeros(6,6)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RRRPP2_invdynJB_fixb_snew_vp2: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RRRPP2_invdynJB_fixb_snew_vp2: qJD has to be [5x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [5 1]), ...
  'S5RRRPP2_invdynJB_fixb_snew_vp2: qJDD has to be [5x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S5RRRPP2_invdynJB_fixb_snew_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [7 1]), ...
  'S5RRRPP2_invdynJB_fixb_snew_vp2: pkin has to be [7x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5RRRPP2_invdynJB_fixb_snew_vp2: m has to be [6x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [6,3]), ...
  'S5RRRPP2_invdynJB_fixb_snew_vp2: mrSges has to be [6x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [6 6]), ...
  'S5RRRPP2_invdynJB_fixb_snew_vp2: Ifges has to be [6x6] (double)'); 

%% Symbolic Calculation
% From invdyn_fixb_NewtonEuler_linkframe_tauJB_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 20:51:38
% EndTime: 2019-12-31 20:51:40
% DurationCPUTime: 1.49s
% Computational Cost: add. (14892->223), mult. (18704->269), div. (0->0), fcn. (8064->6), ass. (0->90)
t706 = Ifges(4,1) + Ifges(5,1) + Ifges(6,1);
t693 = Ifges(4,4) - Ifges(5,5) - Ifges(6,4);
t692 = Ifges(4,5) + Ifges(5,4) - Ifges(6,5);
t705 = Ifges(4,2) + Ifges(5,3) + Ifges(6,2);
t691 = Ifges(4,6) - Ifges(5,6) + Ifges(6,6);
t704 = Ifges(4,3) + Ifges(5,2) + Ifges(6,3);
t657 = qJD(1) + qJD(2);
t665 = sin(qJ(3));
t668 = cos(qJ(3));
t625 = (-mrSges(5,1) * t668 - mrSges(5,3) * t665) * t657;
t656 = qJDD(1) + qJDD(2);
t694 = qJD(3) * t668;
t685 = t657 * t694;
t628 = t665 * t656 + t685;
t667 = sin(qJ(1));
t670 = cos(qJ(1));
t648 = t667 * g(1) - t670 * g(2);
t638 = qJDD(1) * pkin(1) + t648;
t649 = -t670 * g(1) - t667 * g(2);
t672 = qJD(1) ^ 2;
t639 = -t672 * pkin(1) + t649;
t666 = sin(qJ(2));
t669 = cos(qJ(2));
t599 = t666 * t638 + t669 * t639;
t655 = t657 ^ 2;
t596 = -t655 * pkin(2) + t656 * pkin(7) + t599;
t591 = -t668 * g(3) - t665 * t596;
t624 = (-pkin(3) * t668 - qJ(4) * t665) * t657;
t671 = qJD(3) ^ 2;
t697 = t657 * t665;
t590 = -qJDD(3) * pkin(3) - t671 * qJ(4) + t624 * t697 + qJDD(4) - t591;
t689 = -0.2e1 * qJD(5) * t657;
t586 = t665 * t689 + (-t628 + t685) * qJ(5) + (-t655 * t665 * t668 - qJDD(3)) * pkin(4) + t590;
t626 = (mrSges(6,1) * t668 + mrSges(6,2) * t665) * t657;
t696 = t657 * t668;
t644 = qJD(3) * mrSges(6,2) - mrSges(6,3) * t696;
t581 = m(6) * t586 - qJDD(3) * mrSges(6,1) - t628 * mrSges(6,3) - qJD(3) * t644 - t626 * t697;
t646 = mrSges(5,2) * t696 + qJD(3) * mrSges(5,3);
t676 = -m(5) * t590 + qJDD(3) * mrSges(5,1) + qJD(3) * t646 - t581;
t579 = t628 * mrSges(5,2) + t625 * t697 - t676;
t592 = -t665 * g(3) + t668 * t596;
t700 = 2 * qJD(4);
t589 = -t671 * pkin(3) + qJDD(3) * qJ(4) + qJD(3) * t700 + t624 * t696 + t592;
t629 = -qJD(3) * t697 + t668 * t656;
t640 = -qJD(3) * pkin(4) - qJ(5) * t697;
t664 = t668 ^ 2;
t585 = -t664 * t655 * pkin(4) - t629 * qJ(5) + qJD(3) * t640 + t668 * t689 + t589;
t643 = -qJD(3) * mrSges(5,1) + mrSges(5,2) * t697;
t641 = -qJD(3) * mrSges(6,1) - mrSges(6,3) * t697;
t681 = m(6) * t585 + qJDD(3) * mrSges(6,2) - t629 * mrSges(6,3) + qJD(3) * t641;
t678 = m(5) * t589 + qJDD(3) * mrSges(5,3) + qJD(3) * t643 + t625 * t696 + t681;
t686 = (t706 * t665 + t693 * t668) * t657 + t692 * qJD(3);
t687 = (-t693 * t665 - t705 * t668) * t657 - t691 * qJD(3);
t703 = -(t687 * t665 + t686 * t668) * t657 + t704 * qJDD(3) + t692 * t628 + t691 * t629 + mrSges(4,1) * t591 - mrSges(5,1) * t590 - mrSges(6,1) * t586 - mrSges(4,2) * t592 + mrSges(6,2) * t585 + mrSges(5,3) * t589 - pkin(3) * t579 - pkin(4) * t581 + qJ(4) * (t629 * mrSges(5,2) - t626 * t696 + t678);
t699 = t655 * pkin(7);
t698 = mrSges(4,3) + mrSges(5,2);
t627 = (-mrSges(4,1) * t668 + mrSges(4,2) * t665) * t657;
t642 = qJD(3) * mrSges(4,1) - mrSges(4,3) * t697;
t575 = m(4) * t592 - qJDD(3) * mrSges(4,2) - qJD(3) * t642 + (-t626 + t627) * t696 + t698 * t629 + t678;
t645 = -qJD(3) * mrSges(4,2) + mrSges(4,3) * t696;
t576 = m(4) * t591 + qJDD(3) * mrSges(4,1) + qJD(3) * t645 + (-t625 - t627) * t697 - t698 * t628 + t676;
t682 = t668 * t575 - t665 * t576;
t566 = m(3) * t599 - t655 * mrSges(3,1) - t656 * mrSges(3,2) + t682;
t598 = t669 * t638 - t666 * t639;
t680 = t656 * pkin(2) + t598;
t679 = -t628 * qJ(4) - t680;
t583 = qJDD(5) + (-qJ(5) * t664 + pkin(7)) * t655 + (pkin(3) + pkin(4)) * t629 + (qJ(4) * t694 + (-pkin(3) * qJD(3) + t640 + t700) * t665) * t657 - t679;
t580 = m(6) * t583 + t629 * mrSges(6,1) + t628 * mrSges(6,2) + t641 * t697 + t644 * t696;
t587 = -t629 * pkin(3) - t699 + (-0.2e1 * qJD(4) * t665 + (pkin(3) * t665 - qJ(4) * t668) * qJD(3)) * t657 + t679;
t577 = m(5) * t587 - t629 * mrSges(5,1) - t628 * mrSges(5,3) - t643 * t697 - t646 * t696 - t580;
t595 = -t680 - t699;
t673 = -m(4) * t595 + t629 * mrSges(4,1) - t628 * mrSges(4,2) - t642 * t697 + t645 * t696 - t577;
t570 = m(3) * t598 + t656 * mrSges(3,1) - t655 * mrSges(3,2) + t673;
t561 = t666 * t566 + t669 * t570;
t558 = m(2) * t648 + qJDD(1) * mrSges(2,1) - t672 * mrSges(2,2) + t561;
t683 = t669 * t566 - t666 * t570;
t559 = m(2) * t649 - t672 * mrSges(2,1) - qJDD(1) * mrSges(2,2) + t683;
t695 = t670 * t558 + t667 * t559;
t568 = t665 * t575 + t668 * t576;
t688 = (-t692 * t665 - t691 * t668) * t657 - t704 * qJD(3);
t684 = -t667 * t558 + t670 * t559;
t554 = -mrSges(4,1) * t595 + mrSges(4,3) * t592 - mrSges(5,1) * t587 + mrSges(5,2) * t589 + mrSges(6,1) * t583 - mrSges(6,3) * t585 + pkin(4) * t580 - qJ(5) * t681 - pkin(3) * t577 + (qJ(5) * t626 * t668 + t688 * t665) * t657 + t705 * t629 + t693 * t628 + t691 * qJDD(3) + t686 * qJD(3);
t563 = mrSges(4,2) * t595 + mrSges(5,2) * t590 + mrSges(6,2) * t583 - mrSges(4,3) * t591 - mrSges(5,3) * t587 - mrSges(6,3) * t586 - qJ(4) * t577 - qJ(5) * t581 + t687 * qJD(3) + t692 * qJDD(3) + t706 * t628 + t693 * t629 - t688 * t696;
t677 = mrSges(3,1) * t598 - mrSges(3,2) * t599 + Ifges(3,3) * t656 + pkin(2) * t673 + pkin(7) * t682 + t668 * t554 + t665 * t563;
t675 = mrSges(2,1) * t648 - mrSges(2,2) * t649 + Ifges(2,3) * qJDD(1) + pkin(1) * t561 + t677;
t552 = mrSges(3,1) * g(3) + mrSges(3,3) * t599 + t655 * Ifges(3,5) + Ifges(3,6) * t656 - pkin(2) * t568 - t703;
t551 = -mrSges(3,2) * g(3) - mrSges(3,3) * t598 + Ifges(3,5) * t656 - t655 * Ifges(3,6) - pkin(7) * t568 - t665 * t554 + t668 * t563;
t550 = -mrSges(2,2) * g(3) - mrSges(2,3) * t648 + Ifges(2,5) * qJDD(1) - t672 * Ifges(2,6) - pkin(6) * t561 + t669 * t551 - t666 * t552;
t549 = Ifges(2,6) * qJDD(1) + t672 * Ifges(2,5) + mrSges(2,1) * g(3) + mrSges(2,3) * t649 + t666 * t551 + t669 * t552 - pkin(1) * (-m(3) * g(3) + t568) + pkin(6) * t683;
t1 = [-m(1) * g(1) + t684; -m(1) * g(2) + t695; (-m(1) - m(2) - m(3)) * g(3) + t568; -mrSges(1,2) * g(3) + mrSges(1,3) * g(2) - pkin(5) * t695 - t667 * t549 + t670 * t550; mrSges(1,1) * g(3) - mrSges(1,3) * g(1) + pkin(5) * t684 + t670 * t549 + t667 * t550; -mrSges(1,1) * g(2) + mrSges(1,2) * g(1) + t675; t675; t677; t703; t579; t580;];
tauJB = t1;
