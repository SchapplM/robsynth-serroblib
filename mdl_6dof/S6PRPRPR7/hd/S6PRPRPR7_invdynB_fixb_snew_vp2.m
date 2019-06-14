% Calculate vector of inverse dynamics base forces with Newton-Euler for
% S6PRPRPR7
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
%   pkin=[a2,a3,a4,a5,a6,alpha2,d2,d4,d6,theta1]';
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
% Datum: 2019-05-04 23:22
% Revision: 89c353f7eff3bd693eda4e29f35b2761dbc3ada0 (2019-05-03)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauB = S6PRPRPR7_invdynB_fixb_snew_vp2(qJ, qJD, qJDD, g, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(6,1),zeros(3,1),zeros(10,1),zeros(7,1),zeros(7,3),zeros(7,6)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6PRPRPR7_invdynB_fixb_snew_vp2: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6PRPRPR7_invdynB_fixb_snew_vp2: qJD has to be [6x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [6 1]), ...
  'S6PRPRPR7_invdynB_fixb_snew_vp2: qJDD has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6PRPRPR7_invdynB_fixb_snew_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S6PRPRPR7_invdynB_fixb_snew_vp2: pkin has to be [10x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6PRPRPR7_invdynB_fixb_snew_vp2: m has to be [7x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [7,3]), ...
  'S6PRPRPR7_invdynB_fixb_snew_vp2: mrSges has to be [7x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [7 6]), ...
  'S6PRPRPR7_invdynB_fixb_snew_vp2: Ifges has to be [7x6] (double)'); 

%% Symbolic Calculation
% From invdyn_fixb_NewtonEuler_linkframe_tauB_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-05-04 23:19:09
% EndTime: 2019-05-04 23:19:14
% DurationCPUTime: 3.10s
% Computational Cost: add. (32944->274), mult. (61818->331), div. (0->0), fcn. (34861->10), ass. (0->121)
t650 = Ifges(5,1) + Ifges(6,2);
t641 = Ifges(5,4) + Ifges(6,6);
t639 = (Ifges(5,5) - Ifges(6,4));
t649 = Ifges(5,2) + Ifges(6,3);
t637 = (Ifges(5,6) - Ifges(6,5));
t648 = (-Ifges(5,3) - Ifges(6,1));
t593 = sin(pkin(10));
t595 = cos(pkin(10));
t572 = g(1) * t593 - g(2) * t595;
t573 = -g(1) * t595 - g(2) * t593;
t590 = -g(3) + qJDD(1);
t602 = cos(qJ(2));
t596 = cos(pkin(6));
t599 = sin(qJ(2));
t634 = t596 * t599;
t594 = sin(pkin(6));
t635 = t594 * t599;
t533 = t572 * t634 + t602 * t573 + t590 * t635;
t618 = qJDD(2) * qJ(3) + (2 * qJD(3) * qJD(2)) + t533;
t604 = qJD(2) ^ 2;
t644 = (-pkin(2) - pkin(8));
t626 = t644 * t604;
t528 = t626 + t618;
t598 = sin(qJ(4));
t601 = cos(qJ(4));
t628 = qJD(2) * qJD(4);
t625 = t601 * t628;
t569 = qJDD(2) * t598 + t625;
t584 = t598 * t628;
t570 = qJDD(2) * t601 - t584;
t630 = qJD(2) * t598;
t574 = -(qJD(4) * mrSges(5,2)) - mrSges(5,3) * t630;
t629 = qJD(2) * t601;
t575 = (qJD(4) * mrSges(5,1)) - mrSges(5,3) * t629;
t645 = -2 * qJD(5);
t610 = pkin(4) * t625 + t629 * t645 + t618 + (-t570 + t584) * qJ(5);
t523 = t569 * pkin(4) + t610 + t626;
t576 = mrSges(6,1) * t630 - (qJD(4) * mrSges(6,3));
t577 = mrSges(6,1) * t629 + (qJD(4) * mrSges(6,2));
t532 = -t599 * t573 + (t572 * t596 + t590 * t594) * t602;
t608 = -t604 * qJ(3) + qJDD(3) - t532;
t529 = t644 * qJDD(2) + t608;
t547 = -t572 * t594 + t590 * t596;
t524 = t601 * t529 - t598 * t547;
t566 = (pkin(4) * t598 - qJ(5) * t601) * qJD(2);
t603 = qJD(4) ^ 2;
t521 = -qJDD(4) * pkin(4) - t603 * qJ(5) + t566 * t629 + qJDD(5) - t524;
t518 = (t598 * t601 * t604 - qJDD(4)) * pkin(9) + (t570 + t584) * pkin(5) + t521;
t579 = pkin(5) * t629 - qJD(4) * pkin(9);
t589 = t598 ^ 2;
t519 = -t579 * t629 + (pkin(4) + pkin(9)) * t569 + (-pkin(5) * t589 + t644) * t604 + t610;
t597 = sin(qJ(6));
t600 = cos(qJ(6));
t514 = t518 * t600 - t519 * t597;
t564 = -qJD(4) * t597 + t600 * t630;
t540 = qJD(6) * t564 + qJDD(4) * t600 + t569 * t597;
t565 = qJD(4) * t600 + t597 * t630;
t541 = -mrSges(7,1) * t564 + mrSges(7,2) * t565;
t582 = qJD(6) + t629;
t544 = -mrSges(7,2) * t582 + mrSges(7,3) * t564;
t561 = qJDD(6) + t570;
t512 = m(7) * t514 + mrSges(7,1) * t561 - mrSges(7,3) * t540 - t541 * t565 + t544 * t582;
t515 = t518 * t597 + t519 * t600;
t539 = -qJD(6) * t565 - qJDD(4) * t597 + t569 * t600;
t545 = mrSges(7,1) * t582 - mrSges(7,3) * t565;
t513 = m(7) * t515 - mrSges(7,2) * t561 + mrSges(7,3) * t539 + t541 * t564 - t545 * t582;
t622 = -t597 * t512 + t600 * t513;
t606 = m(6) * t523 - t570 * mrSges(6,3) - (t576 * t598 + t577 * t601) * qJD(2) + t622;
t642 = mrSges(5,1) - mrSges(6,2);
t647 = -m(5) * t528 - t570 * mrSges(5,2) - t642 * t569 - t574 * t630 - t575 * t629 - t606;
t643 = mrSges(3,1) - mrSges(4,2);
t640 = (Ifges(3,5) - Ifges(4,4));
t638 = Ifges(3,6) - Ifges(4,5);
t506 = t600 * t512 + t597 * t513;
t611 = -m(6) * t521 - t570 * mrSges(6,1) - t506;
t567 = (-mrSges(6,2) * t598 - mrSges(6,3) * t601) * qJD(2);
t620 = qJD(2) * (-t567 - (mrSges(5,1) * t598 + mrSges(5,2) * t601) * qJD(2));
t504 = m(5) * t524 - t570 * mrSges(5,3) + t642 * qJDD(4) + (t574 - t576) * qJD(4) + t601 * t620 + t611;
t525 = t598 * t529 + t601 * t547;
t609 = -t603 * pkin(4) + qJDD(4) * qJ(5) - t566 * t630 + t525;
t520 = (qJD(4) * t645) - t609;
t517 = -t589 * t604 * pkin(9) - t569 * pkin(5) + ((2 * qJD(5)) + t579) * qJD(4) + t609;
t613 = -m(7) * t517 + t539 * mrSges(7,1) - t540 * mrSges(7,2) + t564 * t544 - t565 * t545;
t607 = -m(6) * t520 + qJDD(4) * mrSges(6,3) + qJD(4) * t577 - t613;
t510 = m(5) * t525 - qJDD(4) * mrSges(5,2) - qJD(4) * t575 + (-mrSges(5,3) - mrSges(6,1)) * t569 + t598 * t620 + t607;
t499 = t601 * t504 + t598 * t510;
t531 = -qJDD(2) * pkin(2) + t608;
t612 = -m(4) * t531 + (t604 * mrSges(4,3)) - t499;
t495 = m(3) * t532 - (t604 * mrSges(3,2)) + t643 * qJDD(2) + t612;
t636 = t495 * t602;
t623 = -t504 * t598 + t601 * t510;
t498 = m(4) * t547 + t623;
t497 = m(3) * t547 + t498;
t530 = t604 * pkin(2) - t618;
t605 = -m(4) * t530 + (t604 * mrSges(4,2)) + qJDD(2) * mrSges(4,3) - t647;
t503 = m(3) * t533 - (t604 * mrSges(3,1)) - qJDD(2) * mrSges(3,2) + t605;
t486 = -t497 * t594 + t503 * t634 + t596 * t636;
t484 = m(2) * t572 + t486;
t492 = -t495 * t599 + t602 * t503;
t491 = m(2) * t573 + t492;
t633 = t595 * t484 + t593 * t491;
t632 = -(t637 * qJD(4)) + (t598 * t649 - t601 * t641) * qJD(2);
t631 = (t639 * qJD(4)) + (-t641 * t598 + t601 * t650) * qJD(2);
t485 = t596 * t497 + t503 * t635 + t594 * t636;
t624 = -t484 * t593 + t595 * t491;
t621 = qJD(2) * ((t648 * qJD(4)) + (t598 * t637 - t601 * t639) * qJD(2));
t505 = -t569 * mrSges(6,2) + t606;
t534 = Ifges(7,5) * t565 + Ifges(7,6) * t564 + Ifges(7,3) * t582;
t536 = Ifges(7,1) * t565 + Ifges(7,4) * t564 + Ifges(7,5) * t582;
t507 = -mrSges(7,1) * t517 + mrSges(7,3) * t515 + Ifges(7,4) * t540 + Ifges(7,2) * t539 + Ifges(7,6) * t561 - t534 * t565 + t536 * t582;
t535 = Ifges(7,4) * t565 + Ifges(7,2) * t564 + Ifges(7,6) * t582;
t508 = mrSges(7,2) * t517 - mrSges(7,3) * t514 + Ifges(7,1) * t540 + Ifges(7,4) * t539 + Ifges(7,5) * t561 + t534 * t564 - t535 * t582;
t487 = -mrSges(5,1) * t528 - mrSges(6,1) * t520 + mrSges(6,2) * t523 + mrSges(5,3) * t525 - pkin(4) * t505 - pkin(5) * t613 - pkin(9) * t622 + t631 * qJD(4) + t637 * qJDD(4) - t600 * t507 - t597 * t508 - t569 * t649 + t641 * t570 + t601 * t621;
t488 = mrSges(6,1) * t521 + mrSges(7,1) * t514 + mrSges(5,2) * t528 - mrSges(7,2) * t515 - mrSges(5,3) * t524 - mrSges(6,3) * t523 + Ifges(7,5) * t540 + Ifges(7,6) * t539 + Ifges(7,3) * t561 + pkin(5) * t506 - qJ(5) * t505 + t565 * t535 - t564 * t536 + t650 * t570 - t641 * t569 + t639 * qJDD(4) + t632 * qJD(4) + t598 * t621;
t481 = -mrSges(4,1) * t530 + mrSges(3,3) * t533 - pkin(2) * t498 - pkin(3) * t647 - pkin(8) * t623 + t638 * qJDD(2) - t601 * t487 - t598 * t488 - t643 * t547 + (t640 * t604);
t482 = pkin(4) * (-qJD(4) * t576 + t611) + qJ(5) * t607 - t597 * t507 + t600 * t508 - mrSges(3,3) * t532 + mrSges(5,1) * t524 - mrSges(5,2) * t525 + mrSges(4,1) * t531 - mrSges(6,3) * t520 + mrSges(6,2) * t521 - pkin(9) * t506 + pkin(3) * t499 - qJ(3) * t498 - t638 * t604 + t639 * t570 + (-mrSges(6,1) * qJ(5) - t637) * t569 + (mrSges(3,2) - mrSges(4,3)) * t547 + (-mrSges(6,2) * pkin(4) - t648) * qJDD(4) + t640 * qJDD(2) + ((-pkin(4) * t567 - t632) * t601 + (-qJ(5) * t567 + t631) * t598) * qJD(2);
t614 = pkin(7) * t492 + t481 * t602 + t482 * t599;
t480 = mrSges(3,1) * t532 - mrSges(3,2) * t533 + mrSges(4,2) * t531 - mrSges(4,3) * t530 + t601 * t488 - t598 * t487 - pkin(8) * t499 + pkin(2) * t612 + qJ(3) * t605 + (-mrSges(4,2) * pkin(2) + Ifges(4,1) + Ifges(3,3)) * qJDD(2);
t479 = mrSges(2,2) * t590 - mrSges(2,3) * t572 - t599 * t481 + t602 * t482 + (-t485 * t594 - t486 * t596) * pkin(7);
t478 = -mrSges(2,1) * t590 + mrSges(2,3) * t573 - pkin(1) * t485 - t594 * t480 + t614 * t596;
t1 = [-m(1) * g(1) + t624; -m(1) * g(2) + t633; -m(1) * g(3) + m(2) * t590 + t485; -mrSges(1,2) * g(3) + mrSges(1,3) * g(2) - qJ(1) * t633 - t593 * t478 + t595 * t479; mrSges(1,1) * g(3) - mrSges(1,3) * g(1) + qJ(1) * t624 + t595 * t478 + t593 * t479; -mrSges(1,1) * g(2) + mrSges(2,1) * t572 + mrSges(1,2) * g(1) - mrSges(2,2) * t573 + pkin(1) * t486 + t596 * t480 + t614 * t594;];
tauB  = t1;
