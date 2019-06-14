% Calculate vector of inverse dynamics joint torques for with Newton-Euler
% S6RRPPRR9
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
%   pkin=[a2,a3,a4,a5,a6,alpha2,d1,d2,d5,d6]';
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
% tauJ [6x1]
%   joint torques of inverse dynamics (contains inertial, gravitational coriolis and centrifugal forces)

% Quelle: HybrDyn-Toolbox
% Datum: 2019-05-06 11:37
% Revision: 89c353f7eff3bd693eda4e29f35b2761dbc3ada0 (2019-05-03)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauJ = S6RRPPRR9_invdynJ_fixb_snew_vp2(qJ, qJD, qJDD, g, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(6,1),zeros(3,1),zeros(10,1),zeros(7,1),zeros(7,3),zeros(7,6)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRPPRR9_invdynJ_fixb_snew_vp2: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RRPPRR9_invdynJ_fixb_snew_vp2: qJD has to be [6x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [6 1]), ...
  'S6RRPPRR9_invdynJ_fixb_snew_vp2: qJDD has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6RRPPRR9_invdynJ_fixb_snew_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S6RRPPRR9_invdynJ_fixb_snew_vp2: pkin has to be [10x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6RRPPRR9_invdynJ_fixb_snew_vp2: m has to be [7x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [7,3]), ...
  'S6RRPPRR9_invdynJ_fixb_snew_vp2: mrSges has to be [7x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [7 6]), ...
  'S6RRPPRR9_invdynJ_fixb_snew_vp2: Ifges has to be [7x6] (double)'); 

%% Symbolic Calculation
% From invdyn_fixb_NewtonEuler_linkframe_tauJ_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-05-06 11:32:09
% EndTime: 2019-05-06 11:32:14
% DurationCPUTime: 3.21s
% Computational Cost: add. (22419->329), mult. (52266->401), div. (0->0), fcn. (34952->10), ass. (0->150)
t681 = -2 * qJD(3);
t680 = -2 * qJD(4);
t679 = Ifges(3,1) + Ifges(5,3) + Ifges(4,2);
t654 = Ifges(3,4) + Ifges(4,6) - Ifges(5,6);
t653 = Ifges(3,5) - Ifges(4,4) + Ifges(5,5);
t678 = Ifges(3,2) + Ifges(4,3) + Ifges(5,2);
t652 = Ifges(3,6) - Ifges(4,5) - Ifges(5,4);
t677 = Ifges(3,3) + Ifges(4,1) + Ifges(5,1);
t622 = qJD(1) ^ 2;
t617 = sin(qJ(1));
t621 = cos(qJ(1));
t638 = -g(1) * t621 - g(2) * t617;
t612 = sin(pkin(6));
t655 = qJDD(1) * t612;
t579 = -pkin(1) * t622 + pkin(8) * t655 + t638;
t616 = sin(qJ(2));
t620 = cos(qJ(2));
t663 = t612 * t620;
t642 = t617 * g(1) - g(2) * t621;
t675 = pkin(8) * t612;
t578 = qJDD(1) * pkin(1) + t622 * t675 + t642;
t613 = cos(pkin(6));
t666 = t578 * t613;
t529 = -g(3) * t663 - t616 * t579 + t620 * t666;
t659 = qJD(1) * t612;
t580 = (-t620 * pkin(2) - t616 * qJ(3)) * t659;
t606 = qJD(1) * t613 + qJD(2);
t604 = t606 ^ 2;
t605 = qJDD(1) * t613 + qJDD(2);
t645 = t616 * t659;
t512 = -t605 * pkin(2) - t604 * qJ(3) + t580 * t645 + qJDD(3) - t529;
t658 = qJD(1) * t620;
t585 = (qJD(2) * t658 + qJDD(1) * t616) * t612;
t644 = t612 * t658;
t639 = t606 * t644;
t665 = t612 ^ 2 * t622;
t649 = t620 * t665;
t499 = t512 - (t616 * t649 + t605) * qJ(4) - (-t585 + t639) * pkin(3) + t606 * t680;
t662 = t620 * t579 + t616 * t666;
t676 = pkin(2) * t604 - t605 * qJ(3) - t580 * t644 + t606 * t681 - t662;
t674 = g(3) * t613;
t673 = mrSges(3,1) - mrSges(4,2);
t672 = -mrSges(3,3) - mrSges(5,1);
t671 = t585 * mrSges(5,1);
t670 = qJ(3) * t606;
t573 = pkin(3) * t645 - qJ(4) * t606;
t669 = t573 * t616;
t576 = mrSges(5,1) * t644 + mrSges(5,2) * t606;
t668 = t576 * t620;
t667 = t578 * t612;
t664 = t612 * t616;
t584 = pkin(4) * t644 - pkin(9) * t606;
t610 = t616 ^ 2;
t611 = t620 ^ 2;
t586 = -qJD(2) * t645 + t620 * t655;
t634 = t645 * t681 - t674 + (t606 * t645 - t586) * pkin(2);
t632 = -t586 * qJ(4) + t644 * t680 + t634;
t492 = (-pkin(3) * t611 - pkin(4) * t610) * t665 + (pkin(9) - qJ(3)) * t585 + (-t578 + (-t669 + (-t584 - t670) * t620) * qJD(1)) * t612 + t632;
t650 = t611 * t665;
t625 = -qJ(4) * t650 + t606 * t573 + qJDD(4) - t676;
t496 = -pkin(9) * t605 + (pkin(3) + pkin(4)) * t586 + (pkin(9) * t649 + (pkin(4) * qJD(1) * t606 - g(3)) * t612) * t616 + t625;
t615 = sin(qJ(5));
t619 = cos(qJ(5));
t489 = t619 * t492 + t615 * t496;
t561 = t606 * t619 + t615 * t645;
t527 = -qJD(5) * t561 + t585 * t619 - t605 * t615;
t560 = -t606 * t615 + t619 * t645;
t531 = -mrSges(6,1) * t560 + mrSges(6,2) * t561;
t593 = qJD(5) + t644;
t536 = mrSges(6,1) * t593 - mrSges(6,3) * t561;
t569 = qJDD(5) + t586;
t532 = -pkin(5) * t560 - pkin(10) * t561;
t590 = t593 ^ 2;
t486 = -pkin(5) * t590 + pkin(10) * t569 + t532 * t560 + t489;
t495 = -pkin(9) * t610 * t665 - pkin(4) * t585 + t606 * t584 - t499;
t528 = qJD(5) * t560 + t585 * t615 + t605 * t619;
t490 = (-t560 * t593 - t528) * pkin(10) + t495 + (t561 * t593 - t527) * pkin(5);
t614 = sin(qJ(6));
t618 = cos(qJ(6));
t483 = -t486 * t614 + t490 * t618;
t533 = -t561 * t614 + t593 * t618;
t504 = qJD(6) * t533 + t528 * t618 + t569 * t614;
t534 = t561 * t618 + t593 * t614;
t513 = -mrSges(7,1) * t533 + mrSges(7,2) * t534;
t559 = qJD(6) - t560;
t515 = -mrSges(7,2) * t559 + mrSges(7,3) * t533;
t524 = qJDD(6) - t527;
t480 = m(7) * t483 + mrSges(7,1) * t524 - mrSges(7,3) * t504 - t513 * t534 + t515 * t559;
t484 = t486 * t618 + t490 * t614;
t503 = -qJD(6) * t534 - t528 * t614 + t569 * t618;
t516 = mrSges(7,1) * t559 - mrSges(7,3) * t534;
t481 = m(7) * t484 - mrSges(7,2) * t524 + mrSges(7,3) * t503 + t513 * t533 - t516 * t559;
t640 = -t480 * t614 + t618 * t481;
t468 = m(6) * t489 - mrSges(6,2) * t569 + mrSges(6,3) * t527 + t531 * t560 - t536 * t593 + t640;
t488 = -t492 * t615 + t496 * t619;
t535 = -mrSges(6,2) * t593 + mrSges(6,3) * t560;
t485 = -pkin(5) * t569 - pkin(10) * t590 + t532 * t561 - t488;
t630 = -m(7) * t485 + t503 * mrSges(7,1) - mrSges(7,2) * t504 + t533 * t515 - t516 * t534;
t476 = m(6) * t488 + mrSges(6,1) * t569 - mrSges(6,3) * t528 - t531 * t561 + t535 * t593 + t630;
t462 = t615 * t468 + t619 * t476;
t470 = t618 * t480 + t614 * t481;
t574 = mrSges(5,1) * t645 - mrSges(5,3) * t606;
t577 = mrSges(4,1) * t645 + mrSges(4,2) * t606;
t661 = -t574 - t577;
t581 = (t620 * mrSges(4,2) - t616 * mrSges(4,3)) * t659;
t583 = (-t616 * mrSges(5,2) - t620 * mrSges(5,3)) * t659;
t660 = t581 + t583;
t651 = g(3) * t664;
t648 = (t616 * t653 + t620 * t652) * t659 + t677 * t606;
t647 = (t616 * t654 + t620 * t678) * t659 + t652 * t606;
t646 = (t616 * t679 + t654 * t620) * t659 + t653 * t606;
t643 = t583 * t659;
t641 = t619 * t468 - t615 * t476;
t498 = -pkin(3) * t650 - qJ(3) * t585 + (-t578 + (-t620 * t670 - t669) * qJD(1)) * t612 + t632;
t637 = -m(5) * t498 + t586 * mrSges(5,3) - t641;
t501 = pkin(3) * t586 + t625 - t651;
t636 = m(5) * t501 + t605 * mrSges(5,2) + t606 * t574 + t620 * t643 + t462;
t633 = -m(6) * t495 + t527 * mrSges(6,1) - t528 * mrSges(6,2) + t560 * t535 - t561 * t536 - t470;
t507 = -t667 + (-t585 - t639) * qJ(3) + t634;
t575 = -mrSges(4,1) * t644 - mrSges(4,3) * t606;
t631 = m(4) * t507 - t585 * mrSges(4,3) + t575 * t644 - t637;
t506 = t651 + t676;
t628 = -m(4) * t506 + t605 * mrSges(4,3) + t606 * t577 + t581 * t644 + t636;
t627 = -m(5) * t499 + t605 * mrSges(5,3) + t606 * t576 - t633;
t508 = Ifges(7,5) * t534 + Ifges(7,6) * t533 + Ifges(7,3) * t559;
t510 = Ifges(7,1) * t534 + Ifges(7,4) * t533 + Ifges(7,5) * t559;
t473 = -mrSges(7,1) * t485 + mrSges(7,3) * t484 + Ifges(7,4) * t504 + Ifges(7,2) * t503 + Ifges(7,6) * t524 - t508 * t534 + t510 * t559;
t509 = Ifges(7,4) * t534 + Ifges(7,2) * t533 + Ifges(7,6) * t559;
t474 = mrSges(7,2) * t485 - mrSges(7,3) * t483 + Ifges(7,1) * t504 + Ifges(7,4) * t503 + Ifges(7,5) * t524 + t508 * t533 - t509 * t559;
t518 = Ifges(6,4) * t561 + Ifges(6,2) * t560 + Ifges(6,6) * t593;
t519 = Ifges(6,1) * t561 + Ifges(6,4) * t560 + Ifges(6,5) * t593;
t626 = mrSges(6,1) * t488 - mrSges(6,2) * t489 + Ifges(6,5) * t528 + Ifges(6,6) * t527 + Ifges(6,3) * t569 + pkin(5) * t630 + pkin(10) * t640 + t618 * t473 + t614 * t474 + t561 * t518 - t560 * t519;
t624 = m(4) * t512 + t585 * mrSges(4,1) - t627;
t623 = mrSges(7,1) * t483 - mrSges(7,2) * t484 + Ifges(7,5) * t504 + Ifges(7,6) * t503 + Ifges(7,3) * t524 + t509 * t534 - t510 * t533;
t582 = (-t620 * mrSges(3,1) + t616 * mrSges(3,2)) * t659;
t572 = -mrSges(3,2) * t606 + mrSges(3,3) * t644;
t571 = mrSges(3,1) * t606 - mrSges(3,3) * t645;
t546 = -t667 - t674;
t530 = -t651 + t662;
t517 = Ifges(6,5) * t561 + Ifges(6,6) * t560 + Ifges(6,3) * t593;
t465 = t616 * t643 - t627 + t671;
t464 = t605 * mrSges(4,2) + t606 * t575 + t660 * t645 + t624 + t671;
t463 = (-t582 - t660) * t645 + m(3) * t529 + (t572 - t575) * t606 + t672 * t585 + t673 * t605 - t624;
t461 = mrSges(5,1) * t586 + t636;
t460 = t586 * mrSges(4,2) - t585 * mrSges(5,2) + (t661 * t616 - t668) * t659 + t631;
t459 = t582 * t644 + t628 + m(3) * t530 - mrSges(3,2) * t605 - t571 * t606 + (mrSges(4,1) - t672) * t586;
t458 = -mrSges(6,1) * t495 + mrSges(6,3) * t489 + Ifges(6,4) * t528 + Ifges(6,2) * t527 + Ifges(6,6) * t569 - pkin(5) * t470 - t517 * t561 + t519 * t593 - t623;
t457 = mrSges(6,2) * t495 - mrSges(6,3) * t488 + Ifges(6,1) * t528 + Ifges(6,4) * t527 + Ifges(6,5) * t569 - pkin(10) * t470 - t473 * t614 + t474 * t618 + t517 * t560 - t518 * t593;
t456 = t619 * t457 - t615 * t458 + qJ(3) * t628 + mrSges(3,1) * t529 - mrSges(3,2) * t530 - mrSges(4,3) * t506 + mrSges(4,2) * t512 - mrSges(5,3) * t499 + mrSges(5,2) * t501 - qJ(4) * t465 - pkin(2) * t464 - pkin(9) * t462 + t677 * t605 + t653 * t585 + (qJ(3) * (mrSges(4,1) + mrSges(5,1)) + t652) * t586 + (t647 * t616 - t646 * t620) * t659;
t1 = [Ifges(2,3) * qJDD(1) + mrSges(2,1) * t642 - mrSges(2,2) * t638 + (t619 * t458 + pkin(9) * t641 + t615 * t457 + mrSges(3,2) * t546 - mrSges(3,3) * t529 + pkin(4) * t633 - mrSges(4,3) * t507 + mrSges(4,1) * t512 - mrSges(5,2) * t498 + mrSges(5,1) * t499 + pkin(3) * t465 - qJ(3) * t460 - t647 * t606 + t653 * t605 + t654 * t586 + t679 * t585 + t648 * t644) * t664 + (qJ(4) * t637 + (qJ(4) * t668 + (qJ(4) * t574 - t648) * t616) * t659 + t646 * t606 + t652 * t605 + t678 * t586 + (mrSges(5,2) * qJ(4) + t654) * t585 - mrSges(3,1) * t546 + mrSges(3,3) * t530 - mrSges(4,1) * t506 + mrSges(4,2) * t507 - mrSges(5,3) * t498 + mrSges(5,1) * t501 + pkin(4) * t462 + t626 + pkin(3) * t461 - pkin(2) * t460) * t663 + t613 * t456 + pkin(1) * ((t459 * t616 + t463 * t620) * t613 + (-m(3) * t546 + t673 * t586 + (-mrSges(3,2) + mrSges(5,2)) * t585 + ((t572 + t576) * t620 + (-t571 - t661) * t616) * t659 - t631) * t612) + (t459 * t620 - t463 * t616) * t675; t456; t464; t461; t626; t623;];
tauJ  = t1;
