% Calculate vector of inverse dynamics base forces with Newton-Euler for
% S5RRPRP8
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
%   pkin=[a2,a3,a4,a5,d1,d2,d4]';
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
% Datum: 2019-12-31 20:04
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauB = S5RRPRP8_invdynB_fixb_snew_vp2(qJ, qJD, qJDD, g, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(5,1),zeros(3,1),zeros(7,1),zeros(6,1),zeros(6,3),zeros(6,6)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RRPRP8_invdynB_fixb_snew_vp2: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RRPRP8_invdynB_fixb_snew_vp2: qJD has to be [5x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [5 1]), ...
  'S5RRPRP8_invdynB_fixb_snew_vp2: qJDD has to be [5x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S5RRPRP8_invdynB_fixb_snew_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [7 1]), ...
  'S5RRPRP8_invdynB_fixb_snew_vp2: pkin has to be [7x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5RRPRP8_invdynB_fixb_snew_vp2: m has to be [6x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [6,3]), ...
  'S5RRPRP8_invdynB_fixb_snew_vp2: mrSges has to be [6x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [6 6]), ...
  'S5RRPRP8_invdynB_fixb_snew_vp2: Ifges has to be [6x6] (double)'); 

%% Symbolic Calculation
% From invdyn_fixb_NewtonEuler_linkframe_tauB_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 20:03:31
% EndTime: 2019-12-31 20:03:34
% DurationCPUTime: 1.90s
% Computational Cost: add. (13536->268), mult. (28423->314), div. (0->0), fcn. (15393->6), ass. (0->104)
t644 = Ifges(3,1) + Ifges(4,1);
t643 = Ifges(5,1) + Ifges(6,1);
t636 = Ifges(3,4) - Ifges(4,5);
t635 = Ifges(5,4) + Ifges(6,4);
t634 = Ifges(3,5) + Ifges(4,4);
t633 = Ifges(5,5) + Ifges(6,5);
t642 = Ifges(3,2) + Ifges(4,3);
t641 = Ifges(5,2) + Ifges(6,2);
t632 = Ifges(3,6) - Ifges(4,6);
t631 = Ifges(5,6) + Ifges(6,6);
t640 = Ifges(3,3) + Ifges(4,2);
t639 = Ifges(5,3) + Ifges(6,3);
t638 = 2 * qJD(3);
t637 = mrSges(3,3) + mrSges(4,2);
t603 = cos(qJ(2));
t606 = qJD(1) ^ 2;
t630 = t603 ^ 2 * t606;
t601 = sin(qJ(1));
t604 = cos(qJ(1));
t582 = -t604 * g(1) - t601 * g(2);
t560 = -t606 * pkin(1) + qJDD(1) * pkin(6) + t582;
t600 = sin(qJ(2));
t539 = -t600 * g(3) + t603 * t560;
t571 = (-mrSges(3,1) * t603 + mrSges(3,2) * t600) * qJD(1);
t620 = qJD(1) * qJD(2);
t617 = t600 * t620;
t573 = t603 * qJDD(1) - t617;
t622 = qJD(1) * t600;
t576 = qJD(2) * mrSges(3,1) - mrSges(3,3) * t622;
t569 = (-pkin(2) * t603 - qJ(3) * t600) * qJD(1);
t605 = qJD(2) ^ 2;
t621 = qJD(1) * t603;
t516 = -t605 * pkin(2) + qJDD(2) * qJ(3) + qJD(2) * t638 + t569 * t621 + t539;
t570 = (-mrSges(4,1) * t603 - mrSges(4,3) * t600) * qJD(1);
t577 = -qJD(2) * mrSges(4,1) + mrSges(4,2) * t622;
t580 = -qJD(2) * pkin(3) - pkin(7) * t622;
t511 = -pkin(3) * t630 - t573 * pkin(7) + qJD(2) * t580 + t516;
t538 = -t603 * g(3) - t600 * t560;
t522 = -qJDD(2) * pkin(2) - t605 * qJ(3) + t569 * t622 + qJDD(3) - t538;
t616 = t603 * t620;
t572 = t600 * qJDD(1) + t616;
t512 = (-t572 + t616) * pkin(7) + (-t600 * t603 * t606 - qJDD(2)) * pkin(3) + t522;
t599 = sin(qJ(4));
t602 = cos(qJ(4));
t504 = -t599 * t511 + t602 * t512;
t557 = (-t599 * t600 - t602 * t603) * qJD(1);
t524 = t557 * qJD(4) + t602 * t572 - t599 * t573;
t558 = (-t599 * t603 + t600 * t602) * qJD(1);
t535 = -t557 * mrSges(6,1) + t558 * mrSges(6,2);
t536 = -t557 * mrSges(5,1) + t558 * mrSges(5,2);
t592 = -qJD(2) + qJD(4);
t541 = -t592 * mrSges(5,2) + t557 * mrSges(5,3);
t591 = -qJDD(2) + qJDD(4);
t499 = -0.2e1 * qJD(5) * t558 + (t557 * t592 - t524) * qJ(5) + (t557 * t558 + t591) * pkin(4) + t504;
t540 = -t592 * mrSges(6,2) + t557 * mrSges(6,3);
t619 = m(6) * t499 + t591 * mrSges(6,1) + t592 * t540;
t494 = m(5) * t504 + t591 * mrSges(5,1) + t592 * t541 + (-t535 - t536) * t558 + (-mrSges(5,3) - mrSges(6,3)) * t524 + t619;
t505 = t602 * t511 + t599 * t512;
t523 = -t558 * qJD(4) - t599 * t572 - t602 * t573;
t543 = t592 * mrSges(6,1) - t558 * mrSges(6,3);
t544 = t592 * mrSges(5,1) - t558 * mrSges(5,3);
t542 = t592 * pkin(4) - t558 * qJ(5);
t550 = t557 ^ 2;
t501 = -t550 * pkin(4) + t523 * qJ(5) + 0.2e1 * qJD(5) * t557 - t592 * t542 + t505;
t618 = m(6) * t501 + t523 * mrSges(6,3) + t557 * t535;
t496 = m(5) * t505 + t523 * mrSges(5,3) + t557 * t536 + (-t543 - t544) * t592 + (-mrSges(5,2) - mrSges(6,2)) * t591 + t618;
t613 = -t599 * t494 + t602 * t496;
t610 = m(4) * t516 + qJDD(2) * mrSges(4,3) + qJD(2) * t577 + t570 * t621 + t613;
t487 = m(3) * t539 - qJDD(2) * mrSges(3,2) - qJD(2) * t576 + t571 * t621 + t637 * t573 + t610;
t578 = -qJD(2) * mrSges(3,2) + mrSges(3,3) * t621;
t490 = t602 * t494 + t599 * t496;
t579 = mrSges(4,2) * t621 + qJD(2) * mrSges(4,3);
t609 = -m(4) * t522 + qJDD(2) * mrSges(4,1) + qJD(2) * t579 - t490;
t488 = m(3) * t538 + qJDD(2) * mrSges(3,1) + qJD(2) * t578 - t637 * t572 + (-t570 - t571) * t622 + t609;
t614 = t603 * t487 - t600 * t488;
t481 = m(2) * t582 - t606 * mrSges(2,1) - qJDD(1) * mrSges(2,2) + t614;
t581 = t601 * g(1) - t604 * g(2);
t559 = -qJDD(1) * pkin(1) - t606 * pkin(6) - t581;
t611 = -t573 * pkin(2) + t559 + (-t572 - t616) * qJ(3);
t513 = (pkin(2) * qJD(2) - (2 * qJD(3))) * t622 + t611;
t507 = -pkin(2) * t617 + t573 * pkin(3) - pkin(7) * t630 - t611 + (t580 + t638) * t622;
t503 = -t523 * pkin(4) - t550 * qJ(5) + t558 * t542 + qJDD(5) + t507;
t612 = m(6) * t503 - t523 * mrSges(6,1) + t524 * mrSges(6,2) - t557 * t540 + t558 * t543;
t608 = -m(5) * t507 + t523 * mrSges(5,1) - t524 * mrSges(5,2) + t557 * t541 - t558 * t544 - t612;
t493 = m(4) * t513 - t573 * mrSges(4,1) - t572 * mrSges(4,3) - t577 * t622 - t579 * t621 + t608;
t607 = -m(3) * t559 + t573 * mrSges(3,1) - t572 * mrSges(3,2) - t576 * t622 + t578 * t621 - t493;
t492 = m(2) * t581 + qJDD(1) * mrSges(2,1) - t606 * mrSges(2,2) + t607;
t629 = t601 * t481 + t604 * t492;
t482 = t600 * t487 + t603 * t488;
t628 = -t631 * t557 - t633 * t558 - t639 * t592;
t627 = t641 * t557 + t635 * t558 + t631 * t592;
t626 = -t635 * t557 - t643 * t558 - t633 * t592;
t625 = t640 * qJD(2) + (t634 * t600 + t632 * t603) * qJD(1);
t624 = -t632 * qJD(2) + (-t636 * t600 - t642 * t603) * qJD(1);
t623 = t634 * qJD(2) + (t644 * t600 + t636 * t603) * qJD(1);
t615 = t604 * t481 - t601 * t492;
t497 = -t524 * mrSges(6,3) - t558 * t535 + t619;
t489 = mrSges(5,2) * t507 + mrSges(6,2) * t503 - mrSges(5,3) * t504 - mrSges(6,3) * t499 - qJ(5) * t497 + t635 * t523 + t643 * t524 - t628 * t557 + t633 * t591 - t627 * t592;
t483 = -mrSges(5,1) * t507 + mrSges(5,3) * t505 - mrSges(6,1) * t503 + mrSges(6,3) * t501 - pkin(4) * t612 + qJ(5) * t618 + (-qJ(5) * t543 - t626) * t592 + (-qJ(5) * mrSges(6,2) + t631) * t591 + t628 * t558 + t635 * t524 + t641 * t523;
t478 = mrSges(3,2) * t559 + mrSges(4,2) * t522 - mrSges(3,3) * t538 - mrSges(4,3) * t513 - pkin(7) * t490 - qJ(3) * t493 + t624 * qJD(2) + t634 * qJDD(2) - t599 * t483 + t602 * t489 + t644 * t572 + t636 * t573 + t625 * t621;
t477 = -mrSges(3,1) * t559 - mrSges(4,1) * t513 + mrSges(4,2) * t516 + mrSges(3,3) * t539 - pkin(2) * t493 - pkin(3) * t608 - pkin(7) * t613 + t623 * qJD(2) + t632 * qJDD(2) - t602 * t483 - t599 * t489 + t636 * t572 + t642 * t573 - t625 * t622;
t476 = t639 * t591 - t640 * qJDD(2) + t631 * t523 + (-qJ(3) * mrSges(4,2) - t632) * t573 + t633 * t524 + (pkin(2) * mrSges(4,2) - t634) * t572 + (t623 * t603 + (pkin(2) * t570 + t624) * t600) * qJD(1) + t626 * t557 + t627 * t558 - pkin(2) * t609 - qJ(3) * t610 - pkin(1) * t482 + Ifges(2,6) * qJDD(1) + mrSges(2,1) * g(3) + pkin(3) * t490 + pkin(4) * t497 + mrSges(6,1) * t499 - mrSges(6,2) * t501 + mrSges(5,1) * t504 - mrSges(5,2) * t505 - mrSges(4,3) * t516 + mrSges(4,1) * t522 - mrSges(3,1) * t538 + mrSges(3,2) * t539 + mrSges(2,3) * t582 + t606 * Ifges(2,5);
t475 = -mrSges(2,2) * g(3) - mrSges(2,3) * t581 + Ifges(2,5) * qJDD(1) - t606 * Ifges(2,6) - pkin(6) * t482 - t600 * t477 + t603 * t478;
t1 = [-m(1) * g(1) + t615; -m(1) * g(2) + t629; (-m(1) - m(2)) * g(3) + t482; -mrSges(1,2) * g(3) + mrSges(1,3) * g(2) - pkin(5) * t629 + t604 * t475 - t601 * t476; mrSges(1,1) * g(3) - mrSges(1,3) * g(1) + pkin(5) * t615 + t601 * t475 + t604 * t476; -mrSges(1,1) * g(2) + mrSges(2,1) * t581 + mrSges(1,2) * g(1) - mrSges(2,2) * t582 + Ifges(2,3) * qJDD(1) + pkin(1) * t607 + pkin(6) * t614 + t603 * t477 + t600 * t478;];
tauB = t1;
