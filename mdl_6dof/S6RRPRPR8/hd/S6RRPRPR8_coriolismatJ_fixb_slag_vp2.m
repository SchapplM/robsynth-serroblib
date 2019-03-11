% Calculate matrix of centrifugal and coriolis load on the joints for
% S6RRPRPR8
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% qJD [6x1]
%   Generalized joint velocities
% pkin [10x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d2,d4,d6,theta3]';
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
% Cq [6x6]
%   matrix of coriolis and centrifugal joint torques

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-09 10:54
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Cq = S6RRPRPR8_coriolismatJ_fixb_slag_vp2(qJ, qJD, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(10,1),zeros(7,1),zeros(7,3),zeros(7,6)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRPRPR8_coriolismatJ_fixb_slag_vp2: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RRPRPR8_coriolismatJ_fixb_slag_vp2: qJD has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S6RRPRPR8_coriolismatJ_fixb_slag_vp2: pkin has to be [10x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6RRPRPR8_coriolismatJ_fixb_slag_vp2: m has to be [7x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [7,3]), ...
  'S6RRPRPR8_coriolismatJ_fixb_slag_vp2: mrSges has to be [7x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [7 6]), ...
  'S6RRPRPR8_coriolismatJ_fixb_slag_vp2: Ifges has to be [7x6] (double)'); 

%% Symbolic Calculation
% From coriolismat_joint_fixb_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 10:49:44
% EndTime: 2019-03-09 10:50:13
% DurationCPUTime: 17.24s
% Computational Cost: add. (25700->757), mult. (54105->1025), div. (0->0), fcn. (57251->8), ass. (0->362)
t470 = cos(pkin(10));
t624 = pkin(8) + qJ(3);
t443 = t624 * t470;
t469 = sin(pkin(10));
t526 = t624 * t469;
t635 = sin(qJ(4));
t636 = cos(qJ(4));
t359 = t443 * t635 + t636 * t526;
t496 = t469 * t636 + t470 * t635;
t284 = -pkin(9) * t496 + t359;
t471 = sin(qJ(6));
t473 = cos(qJ(6));
t360 = t443 * t636 - t526 * t635;
t427 = t635 * t469 - t636 * t470;
t487 = t427 * pkin(9) + t360;
t140 = -t284 * t473 + t471 * t487;
t334 = t427 * t473 - t471 * t496;
t504 = t427 * t471 + t473 * t496;
t551 = Ifges(7,5) * t334 - Ifges(7,6) * t504;
t695 = t284 * t471 + t473 * t487;
t746 = -t695 * mrSges(7,1) + t140 * mrSges(7,2) + t551;
t749 = qJD(6) * t746;
t472 = sin(qJ(2));
t401 = t496 * t472;
t403 = t427 * t472;
t278 = t401 * t473 + t403 * t471;
t266 = Ifges(7,4) * t278;
t474 = cos(qJ(2));
t505 = t401 * t471 - t403 * t473;
t119 = Ifges(7,1) * t505 + t474 * Ifges(7,5) + t266;
t700 = t505 * mrSges(7,1);
t128 = t278 * mrSges(7,2) + t700;
t133 = Ifges(7,2) * t505 - t266;
t465 = t472 * pkin(7);
t562 = t469 * t472;
t438 = pkin(3) * t562 + t465;
t571 = qJ(5) * t403;
t520 = -t438 - t571;
t680 = pkin(4) + pkin(5);
t170 = -t401 * t680 + t520;
t723 = t334 * mrSges(7,2);
t175 = -t504 * mrSges(7,1) - t723;
t325 = Ifges(7,4) * t334;
t178 = Ifges(7,2) * t504 - t325;
t182 = Ifges(7,1) * t504 + t325;
t457 = -pkin(3) * t470 - pkin(2);
t569 = qJ(5) * t496;
t503 = -t457 + t569;
t232 = -t427 * t680 + t503;
t574 = t474 * mrSges(7,2);
t724 = t278 * mrSges(7,3);
t228 = -t574 + t724;
t673 = t228 / 0.2e1;
t748 = -(t178 / 0.4e1 - mrSges(7,3) * t140 / 0.2e1 - t182 / 0.4e1) * t278 + (t119 / 0.4e1 - t133 / 0.4e1) * t334 + t232 * t128 / 0.2e1 - t170 * t175 / 0.2e1 + t551 * t474 / 0.4e1 - t140 * t673;
t747 = t232 * t175;
t703 = Ifges(7,6) * t505;
t727 = Ifges(7,5) * t278;
t552 = t727 - t703;
t743 = t170 * t128 + t552 * t474 / 0.2e1;
t442 = -pkin(2) * t474 - qJ(3) * t472 - pkin(1);
t425 = t470 * t442;
t560 = t470 * t472;
t337 = -pkin(8) * t560 + t425 + (-pkin(7) * t469 - pkin(3)) * t474;
t559 = t470 * t474;
t373 = pkin(7) * t559 + t469 * t442;
t355 = -pkin(8) * t562 + t373;
t550 = -t636 * t337 + t635 * t355;
t144 = -pkin(9) * t403 - t550;
t466 = t474 * pkin(4);
t115 = pkin(5) * t474 - t144 + t466;
t186 = t337 * t635 + t355 * t636;
t554 = t474 * qJ(5);
t166 = t186 - t554;
t632 = t401 * pkin(9);
t127 = t166 + t632;
t57 = t115 * t473 - t127 * t471;
t145 = t186 + t632;
t73 = t144 * t473 + t145 * t471;
t741 = t57 + t73;
t176 = -mrSges(7,1) * t334 + mrSges(7,2) * t504;
t740 = m(7) * t232 + t176;
t444 = pkin(2) * t472 - qJ(3) * t474;
t385 = pkin(7) * t562 + t470 * t444;
t338 = pkin(3) * t472 - pkin(8) * t559 + t385;
t386 = -pkin(7) * t560 + t469 * t444;
t561 = t469 * t474;
t358 = -pkin(8) * t561 + t386;
t188 = t635 * t338 + t636 * t358;
t172 = t472 * qJ(5) + t188;
t519 = -t338 * t636 + t635 * t358;
t173 = -t472 * pkin(4) + t519;
t402 = t496 * t474;
t404 = t427 * t474;
t279 = t402 * t473 + t404 * t471;
t229 = mrSges(7,2) * t472 + mrSges(7,3) * t279;
t282 = t402 * t471 - t404 * t473;
t231 = -mrSges(7,1) * t472 - mrSges(7,3) * t282;
t361 = -mrSges(6,2) * t402 + mrSges(6,3) * t472;
t577 = t472 * mrSges(6,1);
t583 = t404 * mrSges(6,2);
t368 = -t577 - t583;
t440 = -t471 * qJ(5) - t473 * t680;
t441 = t473 * qJ(5) - t471 * t680;
t122 = t404 * pkin(9) - t472 * t680 + t519;
t137 = pkin(9) * t402 + t172;
t66 = t122 * t473 - t137 * t471;
t67 = t122 * t471 + t137 * t473;
t662 = t282 / 0.2e1;
t665 = t279 / 0.2e1;
t698 = Ifges(7,5) * t662 + Ifges(7,6) * t665;
t716 = -t472 / 0.2e1;
t518 = Ifges(7,3) * t716 - t67 * mrSges(7,2) / 0.2e1 + t66 * mrSges(7,1) / 0.2e1 + t698;
t647 = -t441 / 0.2e1;
t648 = -t440 / 0.2e1;
t682 = -mrSges(5,1) / 0.2e1;
t686 = -m(7) / 0.2e1;
t688 = -m(6) / 0.2e1;
t739 = (-pkin(4) * t173 + qJ(5) * t172) * t688 + (t440 * t66 + t441 * t67) * t686 + pkin(4) * t368 / 0.2e1 - qJ(5) * t361 / 0.2e1 - t172 * mrSges(6,3) / 0.2e1 + t173 * mrSges(6,1) / 0.2e1 - t519 * t682 + t188 * mrSges(5,2) / 0.2e1 + t231 * t648 + t229 * t647 + t518;
t535 = -t723 / 0.2e1;
t708 = Ifges(6,4) + Ifges(5,5);
t736 = Ifges(5,6) - Ifges(6,6);
t619 = Ifges(7,4) * t505;
t117 = Ifges(7,2) * t278 + t474 * Ifges(7,6) + t619;
t134 = Ifges(7,1) * t278 - t619;
t733 = -t134 + t117;
t618 = Ifges(7,4) * t504;
t179 = Ifges(7,2) * t334 + t618;
t180 = Ifges(7,1) * t334 - t618;
t732 = t179 - t180;
t500 = t727 / 0.2e1 - t703 / 0.2e1;
t730 = t140 * t471 + t473 * t695;
t575 = t474 * mrSges(7,1);
t699 = t505 * mrSges(7,3);
t230 = t575 - t699;
t672 = t230 / 0.2e1;
t27 = (t724 / 0.2e1 - t228 / 0.2e1 + t574 / 0.2e1) * t473 + (t699 / 0.2e1 + t672 + t575 / 0.2e1) * t471;
t729 = -t278 / 0.2e1;
t658 = -t334 / 0.2e1;
t728 = t334 / 0.2e1;
t536 = -t700 / 0.2e1;
t318 = pkin(4) * t427 - t503;
t345 = mrSges(6,1) * t427 - mrSges(6,3) * t496;
t720 = -m(6) * t318 - t345;
t719 = -t385 * t469 + t386 * t470;
t576 = t473 * mrSges(7,2);
t578 = t471 * mrSges(7,1);
t516 = t576 + t578;
t718 = qJD(6) * t516;
t651 = t402 / 0.2e1;
t650 = -t404 / 0.2e1;
t642 = t472 / 0.2e1;
t657 = -t504 / 0.2e1;
t715 = t504 / 0.2e1;
t713 = t505 / 0.2e1;
t710 = -mrSges(5,1) - mrSges(6,1);
t709 = mrSges(5,2) - mrSges(6,3);
t707 = Ifges(6,2) + Ifges(5,3);
t439 = (pkin(3) * t469 + pkin(7)) * t474;
t702 = t439 * mrSges(5,2);
t697 = -t427 * t708 - t496 * t736;
t696 = -t401 * t708 + t403 * t736;
t694 = t439 * mrSges(5,1) + Ifges(5,4) * t404 / 0.2e1 + Ifges(5,6) * t716 + Ifges(6,5) * t650 + Ifges(6,6) * t642 + (Ifges(5,2) + Ifges(6,3)) * t651;
t468 = t470 ^ 2;
t692 = -0.2e1 * t403;
t691 = 2 * qJD(4);
t690 = m(4) / 0.2e1;
t689 = m(5) / 0.2e1;
t687 = m(6) / 0.2e1;
t685 = -m(7) / 0.4e1;
t684 = m(7) / 0.2e1;
t683 = m(4) * pkin(7);
t681 = mrSges(7,3) / 0.2e1;
t677 = -t695 / 0.2e1;
t675 = t176 / 0.2e1;
t666 = t278 / 0.2e1;
t570 = qJ(5) * t427;
t513 = -pkin(4) * t496 - t570;
t655 = -t513 / 0.2e1;
t654 = -t359 / 0.2e1;
t653 = -t360 / 0.2e1;
t652 = -t402 / 0.2e1;
t646 = -t469 / 0.2e1;
t645 = t470 / 0.2e1;
t644 = -t471 / 0.2e1;
t643 = t471 / 0.2e1;
t641 = -t473 / 0.2e1;
t387 = t401 * qJ(5);
t290 = -pkin(4) * t403 + t387;
t634 = m(6) * t290;
t633 = pkin(7) * t474;
t631 = t57 * mrSges(7,2);
t58 = t115 * t471 + t127 * t473;
t630 = t58 * mrSges(7,1);
t72 = -t144 * t471 + t145 * t473;
t627 = t72 * mrSges(7,1);
t626 = t73 * mrSges(7,2);
t623 = Ifges(4,4) * t469;
t622 = Ifges(4,4) * t470;
t621 = Ifges(5,4) * t403;
t620 = Ifges(5,4) * t496;
t617 = Ifges(4,5) * t470;
t616 = Ifges(6,5) * t401;
t615 = Ifges(6,5) * t427;
t611 = Ifges(4,2) * t469;
t610 = Ifges(4,6) * t469;
t598 = t279 * mrSges(7,1);
t593 = t282 * mrSges(7,2);
t118 = Ifges(7,4) * t282 + Ifges(7,2) * t279 - Ifges(7,6) * t472;
t120 = Ifges(7,1) * t282 + Ifges(7,4) * t279 - t472 * Ifges(7,5);
t130 = -mrSges(7,1) * t278 + mrSges(7,2) * t505;
t131 = t593 - t598;
t169 = t466 + t550;
t218 = t402 * pkin(4) + t404 * qJ(5) + t439;
t171 = pkin(5) * t402 + t218;
t217 = pkin(4) * t401 - t520;
t273 = -Ifges(6,1) * t404 + t472 * Ifges(6,4) + Ifges(6,5) * t402;
t275 = -Ifges(5,1) * t404 - Ifges(5,4) * t402 + t472 * Ifges(5,5);
t293 = mrSges(6,1) * t401 + mrSges(6,3) * t403;
t582 = t404 * mrSges(6,3);
t585 = t402 * mrSges(6,1);
t294 = t582 + t585;
t584 = t404 * mrSges(5,2);
t586 = t402 * mrSges(5,1);
t295 = -t584 + t586;
t573 = t474 * mrSges(6,3);
t587 = t401 * mrSges(6,2);
t362 = -t573 - t587;
t363 = mrSges(5,2) * t474 - t401 * mrSges(5,3);
t364 = -mrSges(5,2) * t472 - mrSges(5,3) * t402;
t365 = -mrSges(5,1) * t474 + t403 * mrSges(5,3);
t366 = mrSges(6,1) * t474 - t403 * mrSges(6,2);
t367 = mrSges(5,1) * t472 + mrSges(5,3) * t404;
t372 = -pkin(7) * t561 + t425;
t398 = t472 * Ifges(4,6) + (-t611 + t622) * t474;
t399 = Ifges(4,5) * t472 + (Ifges(4,1) * t470 - t623) * t474;
t579 = t470 * mrSges(4,2);
t580 = t469 * mrSges(4,1);
t413 = (t579 + t580) * t474;
t434 = mrSges(4,2) * t474 - mrSges(4,3) * t562;
t435 = -mrSges(4,2) * t472 - mrSges(4,3) * t561;
t436 = -mrSges(4,1) * t474 - mrSges(4,3) * t560;
t437 = mrSges(4,1) * t472 - mrSges(4,3) * t559;
t537 = -Ifges(6,6) / 0.2e1 + Ifges(5,6) / 0.2e1;
t538 = -Ifges(5,5) / 0.2e1 - Ifges(6,4) / 0.2e1;
t489 = t402 * t537 - t404 * t538;
t272 = -Ifges(6,1) * t403 - t474 * Ifges(6,4) + t616;
t396 = Ifges(5,4) * t401;
t274 = -Ifges(5,1) * t403 - t474 * Ifges(5,5) - t396;
t530 = t272 / 0.2e1 + t274 / 0.2e1;
t393 = Ifges(6,5) * t403;
t268 = -t474 * Ifges(6,6) + Ifges(6,3) * t401 - t393;
t270 = -Ifges(5,2) * t401 - t474 * Ifges(5,6) - t621;
t531 = t268 / 0.2e1 - t270 / 0.2e1;
t3 = -t550 * t367 + m(5) * (t186 * t188 + t438 * t439 + t519 * t550) - t519 * t365 + t694 * t401 + t119 * t662 + t117 * t665 + t118 * t666 + m(4) * (t372 * t385 + t373 * t386) + t386 * t434 + t373 * t435 + t385 * t436 + t372 * t437 + t438 * t295 + t166 * t361 + t172 * t362 + t188 * t363 + t186 * t364 + t173 * t366 + t169 * t368 + t218 * t293 + t217 * t294 + t66 * t230 + t57 * t231 + t67 * t228 + t58 * t229 + t170 * t131 - t171 * t130 + t120 * t713 + m(6) * (t166 * t172 + t169 * t173 + t217 * t218) + m(7) * (-t170 * t171 + t57 * t66 + t58 * t67) + (-Ifges(7,5) * t505 / 0.2e1 + Ifges(7,6) * t729 + t399 * t645 + t398 * t646 + pkin(7) * t413 - pkin(1) * mrSges(3,1) + (t617 / 0.2e1 - t610 / 0.2e1 - Ifges(3,4)) * t472 + t538 * t403 - t537 * t401) * t472 - t530 * t404 + t531 * t402 - (t702 + t273 / 0.2e1 + t275 / 0.2e1) * t403 + (-pkin(1) * mrSges(3,2) + (Ifges(3,1) - Ifges(3,2) - Ifges(7,3) + t468 * Ifges(4,1) / 0.2e1 - Ifges(4,3) + (t579 + t683) * pkin(7) + (pkin(7) * mrSges(4,1) - t622 + t611 / 0.2e1) * t469 - t707) * t472 + t489 + (Ifges(3,4) + t610 - t617) * t474 + t698) * t474;
t592 = t3 * qJD(1);
t194 = t403 * t680 - t387;
t291 = -t403 * mrSges(6,1) + t401 * mrSges(6,3);
t292 = -t403 * mrSges(5,1) - t401 * mrSges(5,2);
t296 = -Ifges(6,3) * t403 - t616;
t297 = Ifges(5,2) * t403 - t396;
t298 = -Ifges(6,1) * t401 - t393;
t299 = -Ifges(5,1) * t401 + t621;
t548 = -t365 + t366;
t549 = t362 + t363;
t4 = t438 * t292 + t217 * t291 + t290 * t293 + t119 * t729 + t133 * t666 + t72 * t230 + t73 * t228 + t194 * t130 + t548 * t186 - t549 * t550 + (t278 * t57 + t505 * t58) * mrSges(7,3) + m(6) * (-t166 * t550 + t169 * t186 + t217 * t290) + m(7) * (t170 * t194 + t57 * t72 + t58 * t73) - (t298 / 0.2e1 + t299 / 0.2e1 - t186 * mrSges(5,3) - t166 * mrSges(6,2) + t531) * t403 + (t296 / 0.2e1 - t297 / 0.2e1 - t550 * mrSges(5,3) - t169 * mrSges(6,2) - t530) * t401 - t696 * t474 / 0.2e1 + t733 * t713 - t743;
t588 = t4 * qJD(1);
t581 = t427 * mrSges(6,2);
t8 = t57 * t228 - t58 * t230 + (-t58 * mrSges(7,3) + t134 / 0.2e1 - t117 / 0.2e1) * t505 + (-t57 * mrSges(7,3) + t119 / 0.2e1 - t133 / 0.2e1) * t278 + t743;
t572 = t8 * qJD(1);
t16 = -t278 * t228 + t505 * t230 - t548 * t403 - t549 * t401 + m(7) * (-t278 * t58 + t505 * t57) + m(6) * (-t166 * t401 - t169 * t403) + m(5) * (-t186 * t401 - t403 * t550) + (-t434 * t469 - t436 * t470 + m(4) * (-t372 * t470 - t373 * t469)) * t472;
t568 = qJD(1) * t16;
t556 = t473 * t228;
t558 = t471 * t230;
t23 = -(t130 - t293) * t403 + (-t362 - t556 + t558) * t474 + m(7) * (-t170 * t403 + (t471 * t57 - t473 * t58) * t474) + m(6) * (-t166 * t474 + t217 * t403);
t567 = qJD(1) * t23;
t557 = t471 * t504;
t555 = t473 * t334;
t553 = -t558 / 0.2e1 + t556 / 0.2e1;
t162 = t441 * mrSges(7,1) + mrSges(7,2) * t440;
t542 = qJD(6) * t162;
t541 = t685 + t688;
t539 = -mrSges(5,3) / 0.2e1 - mrSges(6,2) / 0.2e1;
t532 = -t557 / 0.2e1;
t528 = -t362 / 0.2e1 - t363 / 0.2e1;
t527 = -t365 / 0.2e1 + t366 / 0.2e1;
t343 = mrSges(6,1) * t496 + t427 * mrSges(6,3);
t344 = mrSges(5,1) * t496 - t427 * mrSges(5,2);
t524 = -t359 * t403 - t360 * t401;
t267 = -t496 * t680 - t570;
t346 = Ifges(6,3) * t496 - t615;
t418 = Ifges(6,5) * t496;
t347 = Ifges(6,3) * t427 + t418;
t421 = Ifges(5,4) * t427;
t348 = -Ifges(5,2) * t496 - t421;
t349 = -Ifges(5,2) * t427 + t620;
t350 = -Ifges(6,1) * t427 + t418;
t351 = Ifges(6,1) * t496 + t615;
t352 = -Ifges(5,1) * t427 - t620;
t353 = Ifges(5,1) * t496 - t421;
t476 = (-t217 * t513 + t290 * t318 + (-t166 + t186) * t359 + (t169 - t550) * t360) * t687 - (t270 / 0.4e1 - t299 / 0.4e1 - t298 / 0.4e1 - t268 / 0.4e1) * t496 + (t58 * t715 + t657 * t72 + t728 * t741) * mrSges(7,3) + (t713 * mrSges(7,3) + t684 * t741 + t672) * t695 + t733 * t504 / 0.4e1 - t697 * t474 / 0.4e1 + (t170 * t267 + t194 * t232 + (t58 - t72) * t140) * t684 + t194 * t675 + t293 * t655 + t732 * t505 / 0.4e1 + t457 * t292 / 0.2e1 + t438 * t344 / 0.2e1 + t217 * t343 / 0.2e1 + t290 * t345 / 0.2e1 + t318 * t291 / 0.2e1 + t267 * t130 / 0.2e1 - (-t349 / 0.4e1 + t352 / 0.4e1 + t350 / 0.4e1 + t347 / 0.4e1 + mrSges(5,3) * t653) * t403 + (t346 / 0.4e1 - t353 / 0.4e1 - t351 / 0.4e1 - t348 / 0.4e1 + mrSges(5,3) * t654) * t401 + t527 * t360 + t528 * t359 + (-t403 * t653 + t401 * t654 - (-t186 / 0.2e1 + t166 / 0.2e1) * t496 + (t550 / 0.2e1 - t169 / 0.2e1) * t427) * mrSges(6,2) - t748 + (-t297 / 0.4e1 - t274 / 0.4e1 - t272 / 0.4e1 + t296 / 0.4e1) * t427;
t2 = t476 + (-Ifges(6,2) / 0.2e1 - Ifges(5,3) / 0.2e1) * t472 + t489 + t739;
t7 = t318 * t343 + t457 * t344 + t178 * t728 + t182 * t658 + t747 - (-t347 / 0.2e1 - t350 / 0.2e1 - t352 / 0.2e1 + t349 / 0.2e1) * t496 + (t346 / 0.2e1 - t348 / 0.2e1 - t351 / 0.2e1 - t353 / 0.2e1) * t427 + t720 * t513 + t732 * t715 + t740 * t267;
t512 = t2 * qJD(1) + t7 * qJD(2);
t477 = (-t278 * t728 + t505 * t657) * mrSges(7,3) + (-t401 * t539 + t528) * t427 - (-t403 * t539 - t527) * t496 + (-t372 * t469 + t373 * t470) * t690 + (-t186 * t427 + t496 * t550 + t524) * t689 + (-t166 * t427 + t169 * t496 + t524) * t687 + (-t140 * t505 - t278 * t695 - t334 * t58 + t504 * t57) * t684 + t230 * t715 + t228 * t658 + t436 * t646 + t434 * t645;
t485 = -m(5) * t439 / 0.2e1 + t218 * t688 + t171 * t686 - t598 / 0.2e1 + t593 / 0.2e1;
t10 = t477 + (-t683 / 0.2e1 - t579 / 0.2e1 - t580 / 0.2e1) * t474 - (-mrSges(5,2) / 0.2e1 + mrSges(6,3) / 0.2e1) * t404 + (-mrSges(6,1) / 0.2e1 + t682) * t402 + t485;
t26 = (-t334 ^ 2 - t504 ^ 2) * mrSges(7,3) + m(7) * (-t140 * t504 - t334 * t695) + (t427 ^ 2 + t496 ^ 2) * (mrSges(6,2) + mrSges(5,3)) + (m(4) * qJ(3) + mrSges(4,3)) * (t469 ^ 2 + t468) + (m(6) + m(5)) * (t359 * t496 - t360 * t427);
t511 = -qJD(1) * t10 - qJD(2) * t26;
t478 = -(-t130 / 0.2e1 + t293 / 0.2e1) * t496 - (t675 - t345 / 0.2e1) * t403 + (t581 / 0.2e1 + (t532 - t555 / 0.2e1) * mrSges(7,3)) * t474 + (-t217 * t496 + t318 * t403 - t360 * t474) * t687 + (t170 * t496 - t232 * t403 - t474 * t730) * t684;
t482 = t173 * t688 + (t471 * t67 + t473 * t66) * t686 + t583 / 0.2e1 + t229 * t644 + t577 / 0.2e1 + t231 * t641;
t14 = t478 + t482;
t47 = (-t720 - t740) * t496;
t510 = qJD(1) * t14 - qJD(2) * t47;
t486 = -t634 / 0.2e1 + (-t278 * t441 + t440 * t505) * t684;
t502 = t634 / 0.2e1 + t194 * t686;
t32 = -t128 - t291 - t292 + t486 - t502;
t484 = t513 * t687 + (-t334 * t441 + t440 * t504) * t684;
t501 = m(6) * t655 + t267 * t686;
t41 = -t175 + t343 + t344 - t484 + t501;
t509 = qJD(1) * t32 - qJD(2) * t41;
t70 = 0.2e1 * mrSges(7,2) * t729 + 0.2e1 * t536;
t80 = 0.2e1 * t657 * mrSges(7,1) + 0.2e1 * t535;
t508 = qJD(1) * t70 + qJD(2) * t80;
t497 = m(7) * (-t278 * t471 + t473 * t505);
t75 = -t497 / 0.2e1 + t541 * t692;
t494 = (-t334 * t471 + t473 * t504) * t684;
t98 = -0.2e1 * t496 * t541 + t494;
t507 = qJD(1) * t75 - qJD(2) * t98;
t498 = m(7) * t730;
t495 = t730 * t684;
t19 = -t747 + (t180 / 0.2e1 - t179 / 0.2e1) * t504 + (t182 / 0.2e1 - t178 / 0.2e1) * t334;
t305 = t557 * t681;
t44 = mrSges(7,3) * t532 + t305;
t479 = (t134 / 0.4e1 - t117 / 0.4e1) * t504 + (mrSges(7,3) * t677 + t180 / 0.4e1 - t179 / 0.4e1) * t505 + t230 * t677 + t748;
t6 = t479 - t518;
t493 = t6 * qJD(1) + t19 * qJD(2) + t44 * qJD(5);
t483 = (t278 * t648 + t505 * t647) * mrSges(7,3) + t440 * t673 + t230 * t647 - t500;
t12 = (t57 / 0.2e1 + t73 / 0.2e1) * mrSges(7,2) + (t58 / 0.2e1 - t72 / 0.2e1) * mrSges(7,1) + t483 + t500;
t492 = t12 * qJD(1) + t162 * qJD(4);
t481 = t573 + (t186 - 0.2e1 * t554) * t688 + ((-t441 * t474 + t58) * t473 + (t440 * t474 - t57) * t471) * t686;
t488 = t186 * t687 + (t471 * t73 + t473 * t72) * t684;
t17 = t481 + t488 + t27;
t259 = mrSges(6,3) + m(6) * qJ(5) + m(7) * (-t440 * t471 + t441 * t473) + t516;
t30 = t495 - t498 / 0.2e1 + ((t728 + t658) * t473 + (t715 + t657) * t471) * mrSges(7,3);
t491 = -t17 * qJD(1) - t30 * qJD(2) + t259 * qJD(4);
t490 = t27 * qJD(1) - t44 * qJD(2) - qJD(4) * t516;
t99 = -t496 * t684 + t494;
t81 = t535 + t723 / 0.2e1;
t74 = -t403 * t687 + t497 / 0.2e1 + (t685 - m(6) / 0.4e1) * t692;
t71 = t536 + t700 / 0.2e1;
t54 = t484 + t501;
t45 = t486 + t502;
t42 = t44 * qJD(6);
t29 = t498 / 0.2e1 + t555 * t681 + t305 - t581 + (-t334 * t641 + t504 * t643) * mrSges(7,3) + m(6) * t360 + t495;
t28 = t699 * t644 + t724 * t641 + (t578 / 0.2e1 + t576 / 0.2e1) * t474 + t553;
t18 = t575 * t644 + t574 * t641 - t587 + (-t278 * t641 + t505 * t643) * mrSges(7,3) - t481 + t488 + t553;
t13 = t478 - t482;
t11 = t630 / 0.2e1 + t631 / 0.2e1 + t627 / 0.2e1 - t626 / 0.2e1 + t483 - t500;
t9 = t477 - t584 / 0.2e1 + t582 / 0.2e1 + t585 / 0.2e1 + t586 / 0.2e1 + t633 * t690 + mrSges(4,2) * t559 / 0.2e1 + mrSges(4,1) * t561 / 0.2e1 - t485;
t5 = t479 + t518;
t1 = Ifges(5,6) * t652 + Ifges(6,6) * t651 + t642 * t707 + t650 * t708 + t476 - t739;
t15 = [qJD(2) * t3 + qJD(3) * t16 + qJD(4) * t4 + qJD(5) * t23 + qJD(6) * t8, t9 * qJD(3) + t1 * qJD(4) + t13 * qJD(5) + t5 * qJD(6) + t592 + (t695 * t229 + 0.2e1 * (t188 * t360 + t359 * t519 + t439 * t457) * t689 + (t275 + t273) * t496 / 0.2e1 + 0.2e1 * (-pkin(2) * t633 + qJ(3) * t719) * t690 + t719 * mrSges(4,3) + 0.2e1 * (t172 * t360 + t173 * t359 + t218 * t318) * t687 + t182 * t662 + t179 * t665 + t347 * t651 + t349 * t652 + (Ifges(4,5) * t469 + Ifges(4,6) * t470) * t642 + t398 * t645 + mrSges(3,2) * t465 + (Ifges(7,5) * t504 + Ifges(7,6) * t334) * t716 - Ifges(3,6) * t472 + t469 * t399 / 0.2e1 + t457 * t295 - pkin(2) * t413 + t359 * t368 + t360 * t361 + t360 * t364 - t359 * t367 + t218 * t345 + t318 * t294 + t232 * t131 + (t334 * t67 - t504 * t66) * mrSges(7,3) - t171 * t176 - (-t173 * mrSges(6,2) - mrSges(5,3) * t519 - t642 * t708 - t702) * t496 + t120 * t715 + t118 * t728 - t140 * t231 + 0.2e1 * (-t140 * t66 - t171 * t232 + t67 * t695) * t684 - t172 * t581 + t470 * qJ(3) * t435 - t469 * qJ(3) * t437 + (Ifges(3,5) + (Ifges(4,2) * t470 + t623) * t646 + (Ifges(4,1) * t469 + t622) * t645 + (-mrSges(4,1) * t470 + mrSges(4,2) * t469 - mrSges(3,1)) * pkin(7)) * t474 + (-t188 * mrSges(5,3) - t642 * t736 + t694) * t427 + (t353 + t351) * t650) * qJD(2), qJD(2) * t9 + qJD(4) * t45 + qJD(5) * t74 + qJD(6) * t71 + t568, t588 + t1 * qJD(2) + t45 * qJD(3) + (mrSges(6,2) * t571 + pkin(4) * t587 + t186 * t710 + t440 * t724 + t441 * t699 + t550 * t709 + t552 + t626 - t627 + t696) * qJD(4) + t18 * qJD(5) + t11 * qJD(6) + ((t440 * t72 + t441 * t73) * t684 + (-pkin(4) * t186 - qJ(5) * t550) * t687) * t691, qJD(2) * t13 + qJD(3) * t74 + qJD(4) * t18 + qJD(6) * t28 + t567, t572 + t5 * qJD(2) + t71 * qJD(3) + t11 * qJD(4) + t28 * qJD(5) + (t552 - t630 - t631) * qJD(6); qJD(3) * t10 + qJD(4) * t2 + qJD(5) * t14 + qJD(6) * t6 - t592, qJD(3) * t26 + qJD(4) * t7 - qJD(5) * t47 + qJD(6) * t19, qJD(4) * t54 + qJD(5) * t99 + qJD(6) * t81 - t511, t54 * qJD(3) + (-mrSges(6,2) * t569 + pkin(4) * t581 + t710 * t360 + t709 * t359 + (t334 * t440 + t441 * t504) * mrSges(7,3) + t697 + t746) * qJD(4) + t29 * qJD(5) - t749 + ((t140 * t441 + t440 * t695) * t684 + (-pkin(4) * t360 - qJ(5) * t359) * t687) * t691 + t512, qJD(3) * t99 + qJD(4) * t29 + t42 + t510, t81 * qJD(3) - qJD(4) * t746 + t493 + t749; -qJD(2) * t10 - qJD(4) * t32 + qJD(5) * t75 + qJD(6) * t70 - t568, qJD(4) * t41 - qJD(5) * t98 + qJD(6) * t80 + t511, 0, -t509, t507, t508; -qJD(2) * t2 + qJD(3) * t32 - qJD(5) * t17 + qJD(6) * t12 - t588, -qJD(3) * t41 - qJD(5) * t30 - t512, t509, qJD(5) * t259 + t542, t491, t492 - t542; -qJD(2) * t14 - qJD(3) * t75 + qJD(4) * t17 - qJD(6) * t27 - t567, qJD(3) * t98 + qJD(4) * t30 + t42 - t510, -t507, -t491 + t718, 0, -t490 - t718; -qJD(2) * t6 - qJD(3) * t70 - qJD(4) * t12 + qJD(5) * t27 - t572, -qJD(3) * t80 - t493, -t508, -qJD(5) * t516 - t492, t490, 0;];
Cq  = t15;
