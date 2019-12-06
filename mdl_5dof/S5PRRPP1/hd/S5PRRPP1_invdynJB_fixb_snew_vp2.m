% Calculate vector of inverse dynamics joint torques and base forces with Newton-Euler
% S5PRRPP1
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
%   pkin=[a2,a3,a4,a5,d2,d3,theta1,theta4]';
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
% Datum: 2019-12-05 16:07
% Revision: 77da58f92bca3eff71542919beafa37024070d86 (2019-12-05)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauJB = S5PRRPP1_invdynJB_fixb_snew_vp2(qJ, qJD, qJDD, g, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(5,1),zeros(3,1),zeros(8,1),zeros(6,1),zeros(6,3),zeros(6,6)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5PRRPP1_invdynJB_fixb_snew_vp2: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5PRRPP1_invdynJB_fixb_snew_vp2: qJD has to be [5x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [5 1]), ...
  'S5PRRPP1_invdynJB_fixb_snew_vp2: qJDD has to be [5x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S5PRRPP1_invdynJB_fixb_snew_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S5PRRPP1_invdynJB_fixb_snew_vp2: pkin has to be [8x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5PRRPP1_invdynJB_fixb_snew_vp2: m has to be [6x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [6,3]), ...
  'S5PRRPP1_invdynJB_fixb_snew_vp2: mrSges has to be [6x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [6 6]), ...
  'S5PRRPP1_invdynJB_fixb_snew_vp2: Ifges has to be [6x6] (double)'); 

%% Symbolic Calculation
% From invdyn_fixb_NewtonEuler_linkframe_tauJB_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-05 16:06:30
% EndTime: 2019-12-05 16:06:33
% DurationCPUTime: 2.40s
% Computational Cost: add. (21864->235), mult. (47067->292), div. (0->0), fcn. (28649->8), ass. (0->100)
t667 = Ifges(5,1) + Ifges(6,1);
t661 = Ifges(5,4) - Ifges(6,5);
t660 = Ifges(5,5) + Ifges(6,4);
t666 = -Ifges(5,2) - Ifges(6,3);
t665 = -Ifges(6,2) - Ifges(5,3);
t659 = Ifges(5,6) - Ifges(6,6);
t626 = sin(pkin(7));
t627 = cos(pkin(7));
t610 = t626 * g(1) - t627 * g(2);
t611 = -t627 * g(1) - t626 * g(2);
t629 = sin(qJ(2));
t631 = cos(qJ(2));
t587 = t629 * t610 + t631 * t611;
t633 = qJD(2) ^ 2;
t582 = -t633 * pkin(2) + qJDD(2) * pkin(6) + t587;
t624 = -g(3) + qJDD(1);
t628 = sin(qJ(3));
t630 = cos(qJ(3));
t563 = -t628 * t582 + t630 * t624;
t650 = qJD(2) * qJD(3);
t646 = t630 * t650;
t608 = t628 * qJDD(2) + t646;
t560 = (-t608 + t646) * qJ(4) + (t628 * t630 * t633 + qJDD(3)) * pkin(3) + t563;
t564 = t630 * t582 + t628 * t624;
t609 = t630 * qJDD(2) - t628 * t650;
t652 = qJD(2) * t628;
t612 = qJD(3) * pkin(3) - qJ(4) * t652;
t623 = t630 ^ 2;
t561 = -t623 * t633 * pkin(3) + t609 * qJ(4) - qJD(3) * t612 + t564;
t625 = sin(pkin(8));
t651 = qJD(2) * t630;
t658 = cos(pkin(8));
t595 = t625 * t652 - t658 * t651;
t663 = -2 * qJD(4);
t557 = t625 * t560 + t658 * t561 + t595 * t663;
t583 = t625 * t608 - t658 * t609;
t596 = (t625 * t630 + t658 * t628) * qJD(2);
t591 = qJD(3) * mrSges(5,1) - t596 * mrSges(5,3);
t574 = t595 * pkin(4) - t596 * qJ(5);
t632 = qJD(3) ^ 2;
t552 = -t632 * pkin(4) + qJDD(3) * qJ(5) + 0.2e1 * qJD(5) * qJD(3) - t595 * t574 + t557;
t592 = -qJD(3) * mrSges(6,1) + t596 * mrSges(6,2);
t648 = m(6) * t552 + qJDD(3) * mrSges(6,3) + qJD(3) * t592;
t575 = t595 * mrSges(6,1) - t596 * mrSges(6,3);
t653 = -t595 * mrSges(5,1) - t596 * mrSges(5,2) - t575;
t662 = -mrSges(5,3) - mrSges(6,2);
t544 = m(5) * t557 - qJDD(3) * mrSges(5,2) - qJD(3) * t591 + t662 * t583 + t653 * t595 + t648;
t637 = t658 * t560 - t625 * t561;
t556 = t596 * t663 + t637;
t584 = t658 * t608 + t625 * t609;
t590 = -qJD(3) * mrSges(5,2) - t595 * mrSges(5,3);
t553 = -qJDD(3) * pkin(4) - t632 * qJ(5) + qJDD(5) + ((2 * qJD(4)) + t574) * t596 - t637;
t593 = -t595 * mrSges(6,2) + qJD(3) * mrSges(6,3);
t640 = -m(6) * t553 + qJDD(3) * mrSges(6,1) + qJD(3) * t593;
t545 = m(5) * t556 + qJDD(3) * mrSges(5,1) + qJD(3) * t590 + t662 * t584 + t653 * t596 + t640;
t538 = t625 * t544 + t658 * t545;
t549 = t584 * mrSges(6,2) + t596 * t575 - t640;
t598 = Ifges(4,6) * qJD(3) + (Ifges(4,4) * t628 + Ifges(4,2) * t630) * qJD(2);
t599 = Ifges(4,5) * qJD(3) + (Ifges(4,1) * t628 + Ifges(4,4) * t630) * qJD(2);
t654 = t660 * qJD(3) - t661 * t595 + t667 * t596;
t655 = t659 * qJD(3) + t666 * t595 + t661 * t596;
t664 = (t628 * t598 - t630 * t599) * qJD(2) + (Ifges(4,3) - t665) * qJDD(3) - t659 * t583 + t660 * t584 + t654 * t595 + t655 * t596 + mrSges(4,1) * t563 + mrSges(5,1) * t556 - mrSges(6,1) * t553 - mrSges(4,2) * t564 - mrSges(5,2) * t557 + mrSges(6,3) * t552 + Ifges(4,5) * t608 + Ifges(4,6) * t609 + pkin(3) * t538 - pkin(4) * t549 + qJ(5) * (-t583 * mrSges(6,2) - t595 * t575 + t648);
t607 = (-mrSges(4,1) * t630 + mrSges(4,2) * t628) * qJD(2);
t614 = -qJD(3) * mrSges(4,2) + mrSges(4,3) * t651;
t534 = m(4) * t563 + qJDD(3) * mrSges(4,1) - t608 * mrSges(4,3) + qJD(3) * t614 - t607 * t652 + t538;
t613 = qJD(3) * mrSges(4,1) - mrSges(4,3) * t652;
t642 = t658 * t544 - t625 * t545;
t535 = m(4) * t564 - qJDD(3) * mrSges(4,2) + t609 * mrSges(4,3) - qJD(3) * t613 + t607 * t651 + t642;
t643 = -t628 * t534 + t630 * t535;
t528 = m(3) * t587 - t633 * mrSges(3,1) - qJDD(2) * mrSges(3,2) + t643;
t586 = t631 * t610 - t629 * t611;
t638 = -qJDD(2) * pkin(2) - t586;
t562 = -t609 * pkin(3) + qJDD(4) + t612 * t652 + (-qJ(4) * t623 - pkin(6)) * t633 + t638;
t555 = -0.2e1 * qJD(5) * t596 + (qJD(3) * t595 - t584) * qJ(5) + (qJD(3) * t596 + t583) * pkin(4) + t562;
t550 = m(6) * t555 + t583 * mrSges(6,1) - t584 * mrSges(6,3) - t596 * t592 + t595 * t593;
t547 = m(5) * t562 + t583 * mrSges(5,1) + t584 * mrSges(5,2) + t595 * t590 + t596 * t591 + t550;
t581 = -t633 * pkin(6) + t638;
t635 = -m(4) * t581 + t609 * mrSges(4,1) - t608 * mrSges(4,2) - t613 * t652 + t614 * t651 - t547;
t540 = m(3) * t586 + qJDD(2) * mrSges(3,1) - t633 * mrSges(3,2) + t635;
t525 = t629 * t528 + t631 * t540;
t523 = m(2) * t610 + t525;
t644 = t631 * t528 - t629 * t540;
t524 = m(2) * t611 + t644;
t657 = t627 * t523 + t626 * t524;
t530 = t630 * t534 + t628 * t535;
t656 = t665 * qJD(3) + t659 * t595 - t660 * t596;
t647 = m(3) * t624 + t530;
t645 = -t626 * t523 + t627 * t524;
t641 = m(2) * t624 + t647;
t536 = -mrSges(5,1) * t562 - mrSges(6,1) * t555 + mrSges(6,2) * t552 + mrSges(5,3) * t557 - pkin(4) * t550 + t654 * qJD(3) + t659 * qJDD(3) + t666 * t583 + t661 * t584 + t656 * t596;
t537 = mrSges(5,2) * t562 + mrSges(6,2) * t553 - mrSges(5,3) * t556 - mrSges(6,3) * t555 - qJ(5) * t550 - t655 * qJD(3) + t660 * qJDD(3) - t661 * t583 + t667 * t584 + t656 * t595;
t597 = Ifges(4,3) * qJD(3) + (Ifges(4,5) * t628 + Ifges(4,6) * t630) * qJD(2);
t517 = -mrSges(4,1) * t581 + mrSges(4,3) * t564 + Ifges(4,4) * t608 + Ifges(4,2) * t609 + Ifges(4,6) * qJDD(3) - pkin(3) * t547 + qJ(4) * t642 + qJD(3) * t599 + t658 * t536 + t625 * t537 - t597 * t652;
t519 = mrSges(4,2) * t581 - mrSges(4,3) * t563 + Ifges(4,1) * t608 + Ifges(4,4) * t609 + Ifges(4,5) * qJDD(3) - qJ(4) * t538 - qJD(3) * t598 - t625 * t536 + t658 * t537 + t597 * t651;
t636 = mrSges(3,1) * t586 - mrSges(3,2) * t587 + Ifges(3,3) * qJDD(2) + pkin(2) * t635 + pkin(6) * t643 + t630 * t517 + t628 * t519;
t515 = -mrSges(3,1) * t624 + mrSges(3,3) * t587 + t633 * Ifges(3,5) + Ifges(3,6) * qJDD(2) - pkin(2) * t530 - t664;
t514 = mrSges(3,2) * t624 - mrSges(3,3) * t586 + Ifges(3,5) * qJDD(2) - t633 * Ifges(3,6) - pkin(6) * t530 - t628 * t517 + t630 * t519;
t513 = mrSges(2,2) * t624 - mrSges(2,3) * t610 - pkin(5) * t525 + t631 * t514 - t629 * t515;
t512 = -mrSges(2,1) * t624 + mrSges(2,3) * t611 - pkin(1) * t647 + pkin(5) * t644 + t629 * t514 + t631 * t515;
t1 = [-m(1) * g(1) + t645; -m(1) * g(2) + t657; -m(1) * g(3) + t641; -mrSges(1,2) * g(3) + mrSges(1,3) * g(2) - qJ(1) * t657 - t626 * t512 + t627 * t513; mrSges(1,1) * g(3) - mrSges(1,3) * g(1) + qJ(1) * t645 + t627 * t512 + t626 * t513; -mrSges(1,1) * g(2) + mrSges(2,1) * t610 + mrSges(1,2) * g(1) - mrSges(2,2) * t611 + pkin(1) * t525 + t636; t641; t636; t664; t547; t549;];
tauJB = t1;
