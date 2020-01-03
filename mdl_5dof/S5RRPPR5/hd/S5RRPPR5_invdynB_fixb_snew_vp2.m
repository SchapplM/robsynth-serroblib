% Calculate vector of inverse dynamics base forces with Newton-Euler for
% S5RRPPR5
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
%   pkin=[a2,a3,a4,a5,d1,d2,d5,theta3]';
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
% Datum: 2019-12-31 19:30
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauB = S5RRPPR5_invdynB_fixb_snew_vp2(qJ, qJD, qJDD, g, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(5,1),zeros(3,1),zeros(8,1),zeros(6,1),zeros(6,3),zeros(6,6)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RRPPR5_invdynB_fixb_snew_vp2: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RRPPR5_invdynB_fixb_snew_vp2: qJD has to be [5x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [5 1]), ...
  'S5RRPPR5_invdynB_fixb_snew_vp2: qJDD has to be [5x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S5RRPPR5_invdynB_fixb_snew_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S5RRPPR5_invdynB_fixb_snew_vp2: pkin has to be [8x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5RRPPR5_invdynB_fixb_snew_vp2: m has to be [6x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [6,3]), ...
  'S5RRPPR5_invdynB_fixb_snew_vp2: mrSges has to be [6x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [6 6]), ...
  'S5RRPPR5_invdynB_fixb_snew_vp2: Ifges has to be [6x6] (double)'); 

%% Symbolic Calculation
% From invdyn_fixb_NewtonEuler_linkframe_tauB_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 19:29:11
% EndTime: 2019-12-31 19:29:16
% DurationCPUTime: 3.41s
% Computational Cost: add. (30611->292), mult. (71176->356), div. (0->0), fcn. (45460->8), ass. (0->115)
t639 = -2 * qJD(3);
t638 = Ifges(4,1) + Ifges(5,1);
t630 = Ifges(4,4) - Ifges(5,5);
t629 = Ifges(4,5) + Ifges(5,4);
t637 = Ifges(4,2) + Ifges(5,3);
t636 = -Ifges(5,2) - Ifges(4,3);
t628 = Ifges(4,6) - Ifges(5,6);
t595 = sin(qJ(2));
t598 = cos(qJ(2));
t616 = qJD(1) * qJD(2);
t577 = t598 * qJDD(1) - t595 * t616;
t620 = qJD(1) * t595;
t578 = qJD(2) * pkin(2) - qJ(3) * t620;
t592 = t598 ^ 2;
t601 = qJD(1) ^ 2;
t596 = sin(qJ(1));
t599 = cos(qJ(1));
t581 = t596 * g(1) - t599 * g(2);
t609 = qJDD(1) * pkin(1) + t581;
t529 = -t577 * pkin(2) + t578 * t620 - (qJ(3) * t592 + pkin(6)) * t601 + qJDD(3) - t609;
t576 = t595 * qJDD(1) + t598 * t616;
t593 = sin(pkin(8));
t627 = cos(pkin(8));
t550 = t627 * t576 + t593 * t577;
t619 = qJD(1) * t598;
t564 = t593 * t620 - t627 * t619;
t618 = qJD(2) * t564;
t635 = t529 + (-t550 + t618) * qJ(4);
t582 = -t599 * g(1) - t596 * g(2);
t571 = -t601 * pkin(1) + qJDD(1) * pkin(6) + t582;
t626 = t595 * t571;
t632 = pkin(2) * t601;
t525 = qJDD(2) * pkin(2) - t576 * qJ(3) - t626 + (qJ(3) * t616 + t595 * t632 - g(3)) * t598;
t554 = -t595 * g(3) + t598 * t571;
t526 = t577 * qJ(3) - qJD(2) * t578 - t592 * t632 + t554;
t565 = (t593 * t598 + t627 * t595) * qJD(1);
t510 = t627 * t525 - t593 * t526 + t565 * t639;
t633 = 2 * qJD(4);
t631 = -mrSges(4,3) - mrSges(5,2);
t511 = t593 * t525 + t627 * t526 + t564 * t639;
t549 = t593 * t576 - t627 * t577;
t556 = qJD(2) * mrSges(4,1) - t565 * mrSges(4,3);
t543 = t564 * pkin(3) - t565 * qJ(4);
t600 = qJD(2) ^ 2;
t506 = -t600 * pkin(3) + qJDD(2) * qJ(4) + qJD(2) * t633 - t564 * t543 + t511;
t557 = -qJD(2) * mrSges(5,1) + t565 * mrSges(5,2);
t507 = -qJDD(2) * pkin(3) - t600 * qJ(4) + t565 * t543 + qJDD(4) - t510;
t501 = (-t550 - t618) * pkin(7) + (t564 * t565 - qJDD(2)) * pkin(4) + t507;
t559 = -qJD(2) * pkin(4) - t565 * pkin(7);
t563 = t564 ^ 2;
t502 = -t563 * pkin(4) + t549 * pkin(7) + qJD(2) * t559 + t506;
t594 = sin(qJ(5));
t597 = cos(qJ(5));
t499 = t597 * t501 - t594 * t502;
t538 = t597 * t564 - t594 * t565;
t515 = t538 * qJD(5) + t594 * t549 + t597 * t550;
t539 = t594 * t564 + t597 * t565;
t521 = -t538 * mrSges(6,1) + t539 * mrSges(6,2);
t588 = -qJD(2) + qJD(5);
t530 = -t588 * mrSges(6,2) + t538 * mrSges(6,3);
t587 = -qJDD(2) + qJDD(5);
t497 = m(6) * t499 + t587 * mrSges(6,1) - t515 * mrSges(6,3) - t539 * t521 + t588 * t530;
t500 = t594 * t501 + t597 * t502;
t514 = -t539 * qJD(5) + t597 * t549 - t594 * t550;
t531 = t588 * mrSges(6,1) - t539 * mrSges(6,3);
t498 = m(6) * t500 - t587 * mrSges(6,2) + t514 * mrSges(6,3) + t538 * t521 - t588 * t531;
t611 = -t594 * t497 + t597 * t498;
t606 = m(5) * t506 + qJDD(2) * mrSges(5,3) + qJD(2) * t557 + t611;
t544 = t564 * mrSges(5,1) - t565 * mrSges(5,3);
t621 = -t564 * mrSges(4,1) - t565 * mrSges(4,2) - t544;
t488 = m(4) * t511 - qJDD(2) * mrSges(4,2) - qJD(2) * t556 + t631 * t549 + t621 * t564 + t606;
t555 = -qJD(2) * mrSges(4,2) - t564 * mrSges(4,3);
t490 = t597 * t497 + t594 * t498;
t558 = -t564 * mrSges(5,2) + qJD(2) * mrSges(5,3);
t604 = -m(5) * t507 + qJDD(2) * mrSges(5,1) + qJD(2) * t558 - t490;
t489 = m(4) * t510 + qJDD(2) * mrSges(4,1) + qJD(2) * t555 + t631 * t550 + t621 * t565 + t604;
t484 = t593 * t488 + t627 * t489;
t553 = -t598 * g(3) - t626;
t575 = (-mrSges(3,1) * t598 + mrSges(3,2) * t595) * qJD(1);
t580 = -qJD(2) * mrSges(3,2) + mrSges(3,3) * t619;
t482 = m(3) * t553 + qJDD(2) * mrSges(3,1) - t576 * mrSges(3,3) + qJD(2) * t580 - t575 * t620 + t484;
t579 = qJD(2) * mrSges(3,1) - mrSges(3,3) * t620;
t612 = t627 * t488 - t593 * t489;
t483 = m(3) * t554 - qJDD(2) * mrSges(3,2) + t577 * mrSges(3,3) - qJD(2) * t579 + t575 * t619 + t612;
t613 = -t595 * t482 + t598 * t483;
t475 = m(2) * t582 - t601 * mrSges(2,1) - qJDD(1) * mrSges(2,2) + t613;
t570 = -t601 * pkin(6) - t609;
t509 = -0.2e1 * qJD(4) * t565 + (qJD(2) * t565 + t549) * pkin(3) + t635;
t504 = -t563 * pkin(7) + (-pkin(3) - pkin(4)) * t549 + (-pkin(3) * qJD(2) + t559 + t633) * t565 - t635;
t607 = -m(6) * t504 + t514 * mrSges(6,1) - t515 * mrSges(6,2) + t538 * t530 - t539 * t531;
t495 = m(5) * t509 + t549 * mrSges(5,1) - t550 * mrSges(5,3) - t565 * t557 + t564 * t558 + t607;
t603 = m(4) * t529 + t549 * mrSges(4,1) + t550 * mrSges(4,2) + t564 * t555 + t565 * t556 + t495;
t602 = -m(3) * t570 + t577 * mrSges(3,1) - t576 * mrSges(3,2) - t579 * t620 + t580 * t619 - t603;
t494 = m(2) * t581 + qJDD(1) * mrSges(2,1) - t601 * mrSges(2,2) + t602;
t625 = t596 * t475 + t599 * t494;
t476 = t598 * t482 + t595 * t483;
t624 = -t628 * qJD(2) + t637 * t564 - t630 * t565;
t623 = t636 * qJD(2) + t628 * t564 - t629 * t565;
t622 = t629 * qJD(2) - t630 * t564 + t638 * t565;
t614 = t599 * t475 - t596 * t494;
t568 = Ifges(3,5) * qJD(2) + (Ifges(3,1) * t595 + Ifges(3,4) * t598) * qJD(1);
t567 = Ifges(3,6) * qJD(2) + (Ifges(3,4) * t595 + Ifges(3,2) * t598) * qJD(1);
t566 = Ifges(3,3) * qJD(2) + (Ifges(3,5) * t595 + Ifges(3,6) * t598) * qJD(1);
t518 = Ifges(6,1) * t539 + Ifges(6,4) * t538 + Ifges(6,5) * t588;
t517 = Ifges(6,4) * t539 + Ifges(6,2) * t538 + Ifges(6,6) * t588;
t516 = Ifges(6,5) * t539 + Ifges(6,6) * t538 + Ifges(6,3) * t588;
t492 = mrSges(6,2) * t504 - mrSges(6,3) * t499 + Ifges(6,1) * t515 + Ifges(6,4) * t514 + Ifges(6,5) * t587 + t538 * t516 - t588 * t517;
t491 = -mrSges(6,1) * t504 + mrSges(6,3) * t500 + Ifges(6,4) * t515 + Ifges(6,2) * t514 + Ifges(6,6) * t587 - t539 * t516 + t588 * t518;
t478 = mrSges(4,2) * t529 + mrSges(5,2) * t507 - mrSges(4,3) * t510 - mrSges(5,3) * t509 - pkin(7) * t490 - qJ(4) * t495 + t624 * qJD(2) + t629 * qJDD(2) - t594 * t491 + t597 * t492 - t630 * t549 + t638 * t550 + t623 * t564;
t477 = -mrSges(4,1) * t529 - mrSges(5,1) * t509 + mrSges(5,2) * t506 + mrSges(4,3) * t511 - pkin(3) * t495 - pkin(4) * t607 - pkin(7) * t611 + t622 * qJD(2) + t628 * qJDD(2) - t597 * t491 - t594 * t492 - t637 * t549 + t630 * t550 + t623 * t565;
t472 = mrSges(3,2) * t570 - mrSges(3,3) * t553 + Ifges(3,1) * t576 + Ifges(3,4) * t577 + Ifges(3,5) * qJDD(2) - qJ(3) * t484 - qJD(2) * t567 - t593 * t477 + t627 * t478 + t566 * t619;
t471 = -mrSges(3,1) * t570 + mrSges(3,3) * t554 + Ifges(3,4) * t576 + Ifges(3,2) * t577 + Ifges(3,6) * qJDD(2) - pkin(2) * t603 + qJ(3) * t612 + qJD(2) * t568 + t627 * t477 + t593 * t478 - t566 * t620;
t470 = -qJ(4) * t606 - pkin(3) * t604 + (qJ(4) * mrSges(5,2) + t628) * t549 + (pkin(3) * mrSges(5,2) - t629) * t550 + (qJ(4) * t544 - t622) * t564 + (pkin(3) * t544 + t624) * t565 + (-Ifges(3,3) + t636) * qJDD(2) + t601 * Ifges(2,5) + mrSges(5,1) * t507 - mrSges(4,1) * t510 + mrSges(4,2) * t511 + mrSges(6,1) * t499 - mrSges(6,2) * t500 + Ifges(6,6) * t514 + Ifges(6,5) * t515 + Ifges(2,6) * qJDD(1) + mrSges(2,1) * g(3) - mrSges(3,1) * t553 + mrSges(3,2) * t554 - pkin(2) * t484 - mrSges(5,3) * t506 + pkin(4) * t490 - Ifges(3,5) * t576 - Ifges(3,6) * t577 + mrSges(2,3) * t582 + Ifges(6,3) * t587 + (-t595 * t567 + t598 * t568) * qJD(1) - pkin(1) * t476 - t538 * t518 + t539 * t517;
t469 = -mrSges(2,2) * g(3) - mrSges(2,3) * t581 + Ifges(2,5) * qJDD(1) - t601 * Ifges(2,6) - pkin(6) * t476 - t595 * t471 + t598 * t472;
t1 = [-m(1) * g(1) + t614; -m(1) * g(2) + t625; (-m(1) - m(2)) * g(3) + t476; -mrSges(1,2) * g(3) + mrSges(1,3) * g(2) - pkin(5) * t625 + t599 * t469 - t596 * t470; mrSges(1,1) * g(3) - mrSges(1,3) * g(1) + pkin(5) * t614 + t596 * t469 + t599 * t470; -mrSges(1,1) * g(2) + mrSges(2,1) * t581 + mrSges(1,2) * g(1) - mrSges(2,2) * t582 + Ifges(2,3) * qJDD(1) + pkin(1) * t602 + pkin(6) * t613 + t598 * t471 + t595 * t472;];
tauB = t1;
