% Calculate matrix of centrifugal and coriolis load on the joints for
% S6RPRRPR6
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% qJD [6x1]
%   Generalized joint velocities
% pkin [11x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d3,d4,d6,theta2,theta5]';
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
% Datum: 2019-03-09 05:18
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Cq = S6RPRRPR6_coriolismatJ_fixb_slag_vp2(qJ, qJD, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(11,1),zeros(7,1),zeros(7,3),zeros(7,6)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPRRPR6_coriolismatJ_fixb_slag_vp2: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RPRRPR6_coriolismatJ_fixb_slag_vp2: qJD has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [11 1]), ...
  'S6RPRRPR6_coriolismatJ_fixb_slag_vp2: pkin has to be [11x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6RPRRPR6_coriolismatJ_fixb_slag_vp2: m has to be [7x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [7,3]), ...
  'S6RPRRPR6_coriolismatJ_fixb_slag_vp2: mrSges has to be [7x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [7 6]), ...
  'S6RPRRPR6_coriolismatJ_fixb_slag_vp2: Ifges has to be [7x6] (double)'); 

%% Symbolic Calculation
% From coriolismat_joint_fixb_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 05:15:33
% EndTime: 2019-03-09 05:16:02
% DurationCPUTime: 15.80s
% Computational Cost: add. (44213->705), mult. (88785->977), div. (0->0), fcn. (105119->10), ass. (0->360)
t452 = sin(qJ(6));
t455 = cos(qJ(6));
t449 = sin(pkin(11));
t453 = sin(qJ(4));
t456 = cos(qJ(4));
t578 = cos(pkin(11));
t481 = t449 * t453 - t578 * t456;
t482 = t449 * t456 + t453 * t578;
t362 = -t452 * t481 + t455 * t482;
t506 = -t452 * t482 - t455 * t481;
t549 = Ifges(7,5) * t506 - Ifges(7,6) * t362;
t615 = -qJ(5) - pkin(8);
t432 = t615 * t453;
t433 = t615 * t456;
t386 = t449 * t432 - t578 * t433;
t313 = pkin(9) * t481 - t386;
t705 = t578 * t432 + t449 * t433;
t719 = -pkin(9) * t482 + t705;
t202 = t313 * t455 - t452 * t719;
t728 = t313 * t452 + t455 * t719;
t752 = t202 * mrSges(7,1) - t728 * mrSges(7,2);
t29 = t549 + t752;
t754 = t29 * qJD(6);
t753 = t362 * mrSges(7,1) + t506 * mrSges(7,2);
t747 = qJD(6) * t753;
t450 = sin(pkin(10));
t451 = cos(pkin(10));
t454 = sin(qJ(3));
t633 = cos(qJ(3));
t423 = t450 * t454 - t451 * t633;
t425 = t450 * t633 + t454 * t451;
t324 = t481 * t425;
t326 = t482 * t425;
t507 = t324 * t452 - t455 * t326;
t230 = t324 * t455 + t326 * t452;
t608 = Ifges(7,4) * t230;
t107 = Ifges(7,2) * t507 + t423 * Ifges(7,6) - t608;
t224 = Ifges(7,4) * t507;
t109 = -Ifges(7,1) * t230 + t423 * Ifges(7,5) + t224;
t177 = -mrSges(7,2) * t423 + mrSges(7,3) * t507;
t179 = mrSges(7,1) * t423 + t230 * mrSges(7,3);
t616 = pkin(7) + qJ(2);
t431 = t616 * t451;
t513 = t616 * t450;
t377 = t431 * t454 + t633 * t513;
t563 = t425 * t453;
t304 = pkin(4) * t563 + t377;
t237 = pkin(5) * t326 + t304;
t444 = -pkin(4) * t456 - pkin(3);
t390 = pkin(5) * t481 + t444;
t679 = t202 / 0.2e1;
t680 = t728 / 0.2e1;
t718 = Ifges(7,2) * t230 + t224;
t722 = t230 * mrSges(7,1);
t731 = t507 * mrSges(7,2) - t722;
t733 = Ifges(7,1) * t507 + t608;
t751 = t179 * t679 + t177 * t680 + (t718 / 0.4e1 + t109 / 0.4e1) * t506 + t731 * t390 / 0.2e1 + t753 * t237 / 0.2e1 + (-t107 / 0.4e1 + t733 / 0.4e1) * t362;
t625 = pkin(9) * t324;
t524 = -pkin(2) * t451 - pkin(1);
t628 = pkin(8) * t425;
t352 = pkin(3) * t423 + t524 - t628;
t378 = t431 * t633 - t454 * t513;
t268 = t456 * t352 - t453 * t378;
t562 = t425 * t456;
t225 = -qJ(5) * t562 + t268;
t180 = t423 * pkin(4) + t225;
t269 = t453 * t352 + t378 * t456;
t226 = -qJ(5) * t563 + t269;
t196 = t449 * t226;
t99 = t578 * t180 - t196;
t72 = pkin(5) * t423 + t625 + t99;
t511 = t578 * t226;
t100 = t449 * t180 + t511;
t624 = pkin(9) * t326;
t75 = t100 - t624;
t48 = t452 * t72 + t455 * t75;
t114 = -t225 * t449 - t511;
t81 = t114 + t624;
t115 = t578 * t225 - t196;
t82 = t115 + t625;
t58 = -t452 * t82 + t455 * t81;
t750 = t48 + t58;
t744 = t115 - t99;
t743 = t237 * t731;
t741 = t390 * t753;
t629 = pkin(8) * t423;
t632 = pkin(3) * t425;
t376 = t629 + t632;
t274 = t456 * t376 + t453 * t377;
t564 = t423 * t456;
t189 = t425 * pkin(4) + qJ(5) * t564 + t274;
t275 = t453 * t376 - t456 * t377;
t565 = t423 * t453;
t246 = qJ(5) * t565 + t275;
t112 = t578 * t189 - t246 * t449;
t113 = t449 * t189 + t578 * t246;
t325 = t482 * t423;
t327 = t481 * t423;
t231 = t325 * t455 - t327 * t452;
t176 = -mrSges(7,2) * t425 + t231 * mrSges(7,3);
t234 = t325 * t452 + t327 * t455;
t178 = mrSges(7,1) * t425 - t234 * mrSges(7,3);
t281 = -mrSges(6,2) * t425 + t325 * mrSges(6,3);
t283 = mrSges(6,1) * t425 - mrSges(6,3) * t327;
t523 = t578 * pkin(4);
t443 = t523 + pkin(5);
t631 = pkin(4) * t449;
t408 = t443 * t455 - t452 * t631;
t409 = t443 * t452 + t455 * t631;
t662 = t327 / 0.2e1;
t664 = t325 / 0.2e1;
t487 = Ifges(6,5) * t662 + Ifges(6,6) * t664;
t672 = t234 / 0.2e1;
t676 = t231 / 0.2e1;
t486 = Ifges(7,5) * t672 + Ifges(7,6) * t676;
t74 = pkin(5) * t425 - pkin(9) * t327 + t112;
t76 = pkin(9) * t325 + t113;
t54 = -t452 * t76 + t455 * t74;
t55 = t452 * t74 + t455 * t76;
t639 = t425 / 0.2e1;
t497 = Ifges(7,3) * t639 - t55 * mrSges(7,2) / 0.2e1 + t54 * mrSges(7,1) / 0.2e1 + t486;
t518 = t565 / 0.2e1;
t445 = Ifges(5,5) * t456;
t534 = -t445 / 0.2e1;
t687 = m(6) * pkin(4);
t541 = t687 / 0.2e1;
t646 = t408 / 0.2e1;
t689 = m(7) / 0.2e1;
t714 = Ifges(5,3) + Ifges(6,3);
t740 = (t408 * t54 + t409 * t55) * t689 + t112 * mrSges(6,1) / 0.2e1 - t113 * mrSges(6,2) / 0.2e1 + t274 * mrSges(5,1) / 0.2e1 - t275 * mrSges(5,2) / 0.2e1 + t178 * t646 + t409 * t176 / 0.2e1 + (t112 * t578 + t113 * t449) * t541 + t423 * t534 + Ifges(5,6) * t518 + t281 * t631 / 0.2e1 + t283 * t523 / 0.2e1 + t487 + t714 * t639 + t497;
t607 = Ifges(7,4) * t362;
t264 = Ifges(7,2) * t506 + t607;
t734 = Ifges(7,1) * t506 - t607;
t739 = t264 / 0.4e1 - t734 / 0.4e1;
t678 = t722 / 0.2e1;
t585 = t362 * mrSges(7,3);
t706 = t100 + t114;
t282 = -mrSges(6,2) * t423 - t326 * mrSges(6,3);
t284 = mrSges(6,1) * t423 + t324 * mrSges(6,3);
t642 = t482 / 0.2e1;
t655 = t362 / 0.2e1;
t658 = -t506 / 0.2e1;
t472 = t179 * t655 + t177 * t658 + t481 * t282 / 0.2e1 + t284 * t642;
t47 = -t452 * t75 + t455 * t72;
t494 = -t362 * t47 + t48 * t506;
t366 = -t423 * mrSges(5,2) - mrSges(5,3) * t563;
t557 = t456 * t366;
t368 = t423 * mrSges(5,1) - mrSges(5,3) * t562;
t559 = t453 * t368;
t59 = t452 * t81 + t455 * t82;
t691 = m(6) / 0.2e1;
t737 = t472 - (-t706 * t481 + t482 * t744) * t691 - (t362 * t59 + t506 * t58 + t494) * t689 + t559 / 0.2e1 - t557 / 0.2e1;
t736 = t472 - (-t100 * t481 + t324 * t705 - t326 * t386 - t482 * t99) * t691 - (-t202 * t507 + t230 * t728 + t494) * t689;
t594 = t234 * mrSges(7,2);
t597 = t231 * mrSges(7,1);
t556 = t597 / 0.2e1 - t594 / 0.2e1;
t590 = t327 * mrSges(6,2);
t591 = t325 * mrSges(6,1);
t735 = t556 - t590 / 0.2e1 + t591 / 0.2e1;
t710 = Ifges(7,5) * t507;
t724 = Ifges(7,6) * t230;
t555 = t710 + t724;
t358 = Ifges(7,4) * t506;
t267 = Ifges(7,1) * t362 + t358;
t717 = -Ifges(7,2) * t362 + t358;
t730 = t717 / 0.4e1 + t267 / 0.4e1;
t516 = -t724 / 0.2e1 - t710 / 0.2e1;
t580 = t456 * mrSges(5,1);
t581 = t453 * mrSges(5,2);
t496 = t580 - t581;
t347 = t496 * t425;
t727 = -t347 / 0.2e1;
t726 = -t362 / 0.2e1;
t725 = t705 / 0.2e1;
t582 = t482 * mrSges(6,3);
t720 = -t274 * t453 + t275 * t456;
t715 = -t507 / 0.2e1;
t603 = Ifges(5,6) * t453;
t707 = Ifges(4,4) + t534 + t603 / 0.2e1;
t446 = Ifges(5,4) * t456;
t704 = -Ifges(5,2) * t453 + t446;
t437 = Ifges(5,1) * t453 + t446;
t703 = t386 * t324 / 0.2e1 + t326 * t725;
t702 = -Ifges(6,5) * t481 - Ifges(6,6) * t482 + t549;
t701 = -Ifges(6,5) * t326 + Ifges(6,6) * t324 + t555;
t694 = t456 ^ 2;
t695 = t453 ^ 2;
t700 = t694 + t695;
t698 = -mrSges(5,3) * t700 / 0.2e1;
t124 = -mrSges(7,1) * t507 - mrSges(7,2) * t230;
t610 = Ifges(6,4) * t324;
t186 = -Ifges(6,2) * t326 + t423 * Ifges(6,6) - t610;
t241 = -t324 * mrSges(6,1) - t326 * mrSges(6,2);
t245 = -Ifges(6,1) * t326 + t610;
t261 = -mrSges(7,1) * t506 + mrSges(7,2) * t362;
t539 = pkin(4) * t562;
t280 = -t324 * pkin(5) + t539;
t370 = mrSges(6,1) * t482 - mrSges(6,2) * t481;
t418 = Ifges(6,4) * t481;
t372 = -Ifges(6,2) * t482 - t418;
t609 = Ifges(6,4) * t482;
t373 = -Ifges(6,2) * t481 + t609;
t374 = -Ifges(6,1) * t481 - t609;
t375 = Ifges(6,1) * t482 - t418;
t630 = pkin(4) * t453;
t393 = pkin(5) * t482 + t630;
t434 = t453 * mrSges(5,1) + t456 * mrSges(5,2);
t636 = t444 / 0.2e1;
t640 = t423 / 0.4e1;
t647 = t393 / 0.2e1;
t657 = t506 / 0.2e1;
t666 = t280 / 0.2e1;
t677 = t230 / 0.2e1;
t685 = -t47 / 0.2e1;
t696 = t377 * t434 / 0.2e1 + t304 * t370 / 0.2e1 + t702 * t640 - (t186 / 0.4e1 - t245 / 0.4e1) * t482 + t282 * t725 + pkin(3) * t727 - t386 * t284 / 0.2e1 + t730 * t507 + t261 * t666 + t124 * t647 + t241 * t636 - (t372 / 0.4e1 + t375 / 0.4e1) * t326 + (t373 / 0.4e1 - t374 / 0.4e1) * t324 + t739 * t230 + (t237 * t393 + t280 * t390 + t750 * t728 + (t47 - t59) * t202) * t689 + (-t202 * t677 + t506 * t685 + t59 * t657 + t728 * t715 + t726 * t750) * mrSges(7,3) + t751;
t693 = m(5) / 0.2e1;
t692 = -m(6) / 0.2e1;
t690 = -m(7) / 0.2e1;
t688 = -pkin(8) / 0.2e1;
t684 = t109 / 0.2e1;
t681 = -t728 / 0.2e1;
t674 = t507 / 0.2e1;
t671 = -t230 / 0.2e1;
t668 = t267 / 0.2e1;
t665 = -t324 / 0.2e1;
t663 = -t326 / 0.2e1;
t611 = Ifges(5,4) * t453;
t435 = Ifges(5,2) * t456 + t611;
t350 = t425 * t435;
t661 = -t350 / 0.4e1;
t351 = t425 * t437;
t660 = -t351 / 0.4e1;
t653 = t373 / 0.2e1;
t650 = t375 / 0.2e1;
t645 = -t409 / 0.2e1;
t644 = -t481 / 0.2e1;
t641 = t423 / 0.2e1;
t638 = -t435 / 0.4e1;
t438 = Ifges(5,1) * t456 - t611;
t637 = t438 / 0.4e1;
t635 = t453 / 0.2e1;
t634 = t456 / 0.2e1;
t627 = pkin(8) * t453;
t626 = pkin(8) * t456;
t622 = t47 * mrSges(7,2);
t621 = t48 * mrSges(7,1);
t618 = t58 * mrSges(7,1);
t617 = t59 * mrSges(7,2);
t614 = mrSges(7,3) * t408;
t613 = mrSges(7,3) * t409;
t606 = Ifges(5,5) * t423;
t604 = Ifges(5,6) * t423;
t106 = Ifges(7,4) * t234 + Ifges(7,2) * t231 + t425 * Ifges(7,6);
t108 = Ifges(7,1) * t234 + Ifges(7,4) * t231 + t425 * Ifges(7,5);
t123 = t594 - t597;
t185 = Ifges(6,4) * t327 + Ifges(6,2) * t325 + t425 * Ifges(6,6);
t187 = Ifges(6,1) * t327 + Ifges(6,4) * t325 + t425 * Ifges(6,5);
t320 = Ifges(6,4) * t326;
t188 = -Ifges(6,1) * t324 + t423 * Ifges(6,5) - t320;
t305 = -pkin(4) * t565 + t378;
t238 = -t325 * pkin(5) + t305;
t242 = t590 - t591;
t243 = mrSges(6,1) * t326 - mrSges(6,2) * t324;
t296 = Ifges(5,6) * t425 - t423 * t704;
t298 = Ifges(5,5) * t425 - t423 * t438;
t348 = t434 * t423;
t349 = t434 * t425;
t365 = -t425 * mrSges(5,2) + mrSges(5,3) * t565;
t367 = t425 * mrSges(5,1) + mrSges(5,3) * t564;
t415 = t425 * mrSges(4,1);
t299 = t425 * t438 + t606;
t558 = t456 * t299;
t297 = t425 * t704 + t604;
t560 = t453 * t297;
t3 = m(5) * (t268 * t274 + t269 * t275 + t377 * t378) + m(6) * (t100 * t113 + t112 * t99 + t304 * t305) + m(7) * (t237 * t238 + t47 * t54 + t48 * t55) - t377 * t348 + t378 * t349 + t269 * t365 + t275 * t366 + t268 * t367 + t274 * t368 + t304 * t242 + t305 * t243 + t100 * t281 + t113 * t282 + t99 * t283 + t112 * t284 + t237 * t123 + t238 * t124 + t48 * t176 + t55 * t177 + t47 * t178 + t54 * t179 + t188 * t662 + t185 * t663 + t186 * t664 + t187 * t665 + t108 * t671 + t109 * t672 + t106 * t674 + t107 * t676 + t524 * t415 + (-t524 * mrSges(4,2) - t558 / 0.2e1 + t560 / 0.2e1 + t486 + t487 + t707 * t423) * t423 + (Ifges(6,5) * t665 + Ifges(6,6) * t663 + Ifges(7,5) * t671 + Ifges(7,6) * t674 + t298 * t634 - t453 * t296 / 0.2e1 - t707 * t425 + (Ifges(4,2) + Ifges(7,3) - Ifges(4,1) + t714) * t423) * t425;
t592 = t3 * qJD(1);
t588 = t506 * mrSges(7,3);
t244 = Ifges(6,2) * t324 - t320;
t4 = t377 * t347 + t268 * t366 - t269 * t368 + t304 * t241 + t280 * t124 + t115 * t282 + t114 * t284 + t743 + t718 * t674 + t507 * t684 + t733 * t671 + t107 * t677 + t59 * t177 + t58 * t179 + (t230 * t48 - t47 * t507) * mrSges(7,3) + m(7) * (t237 * t280 + t47 * t58 + t48 * t59) + m(6) * (t100 * t115 + t99 * t114) - (t244 / 0.2e1 + t188 / 0.2e1 - t99 * mrSges(6,3)) * t326 + (t186 / 0.2e1 - t245 / 0.2e1 + t100 * mrSges(6,3)) * t324 + ((-t606 / 0.2e1 + t350 / 0.2e1 - t299 / 0.2e1 + t268 * mrSges(5,3)) * t453 + (-t604 / 0.2e1 - t351 / 0.2e1 - t297 / 0.2e1 - t269 * mrSges(5,3) + (m(6) * t304 + t243) * pkin(4)) * t456) * t425 + t701 * t641;
t584 = t4 * qJD(1);
t583 = t481 * mrSges(6,3);
t7 = t47 * t177 - t48 * t179 + t555 * t641 + t743 - (-t48 * mrSges(7,3) + t733 / 0.2e1 - t107 / 0.2e1) * t230 + (-t47 * mrSges(7,3) + t684 + t718 / 0.2e1) * t507;
t579 = t7 * qJD(1);
t371 = mrSges(6,1) * t481 + mrSges(6,2) * t482;
t460 = (t274 * t456 + t453 * t275) * t693 + (-t112 * t481 + t113 * t482) * t691 + (t362 * t55 + t506 * t54) * t689 + t178 * t657 + t176 * t655 + t283 * t644 + t281 * t642 + t365 * t635 + t367 * t634;
t461 = -m(5) * (-t629 * t700 - t632) / 0.2e1 + (t325 * t705 + t327 * t386 + t425 * t444) * t692 + (-t202 * t234 + t231 * t728 + t390 * t425) * t690;
t12 = (-mrSges(4,2) + (t694 / 0.2e1 + t695 / 0.2e1) * mrSges(5,3)) * t423 + (t231 * t655 + t234 * t658) * mrSges(7,3) + (t325 * t642 + t481 * t662) * mrSges(6,3) + t460 + t415 + (-t261 / 0.2e1 + t580 / 0.2e1 - t581 / 0.2e1 - t371 / 0.2e1) * t425 + t461;
t577 = qJD(1) * t12;
t21 = t507 * t177 + t230 * t179 - t326 * t282 + t324 * t284 + m(7) * (t230 * t47 + t48 * t507) + m(6) * (-t100 * t326 + t324 * t99);
t576 = qJD(1) * t21;
t568 = t377 * t425;
t14 = t234 * t177 + t231 * t179 + t327 * t282 + t325 * t284 + (mrSges(4,3) * t423 - t557 + t559) * t423 + (mrSges(4,3) * t425 + t124 + t243 + t349) * t425 + m(7) * (t231 * t47 + t234 * t48 + t237 * t425) + m(6) * (t100 * t327 + t304 * t425 + t325 * t99) + m(5) * (t568 + (t268 * t453 - t269 * t456) * t423) + m(4) * (-t378 * t423 + t568) + (m(3) * qJ(2) + mrSges(3,3)) * (t450 ^ 2 + t451 ^ 2);
t575 = t14 * qJD(1);
t470 = (-t230 * t726 + t507 * t658) * mrSges(7,3) + t177 * t657 + t179 * t726;
t18 = t470 - t556;
t574 = t18 * qJD(1);
t571 = t304 * t453;
t167 = t409 * mrSges(7,1) + mrSges(7,2) * t408;
t543 = t167 * qJD(6);
t542 = mrSges(6,3) * t631;
t538 = t630 / 0.2e1;
t537 = -t627 / 0.2e1;
t536 = -t626 / 0.2e1;
t533 = t588 / 0.2e1;
t532 = -t585 / 0.2e1;
t531 = t583 / 0.2e1;
t530 = -t583 / 0.2e1;
t529 = -t582 / 0.2e1;
t508 = t445 - t603;
t504 = t539 / 0.2e1;
t503 = mrSges(6,3) * t523;
t501 = -t326 * t531 - t582 * t665 + t585 * t677 + t588 * t715;
t500 = t230 * t532 + t324 * t529 - t326 * t530 + t507 * t533;
t498 = 0.2e1 * (m(7) / 0.4e1 + m(6) / 0.4e1) * t425;
t27 = t741 + (t734 / 0.2e1 - t264 / 0.2e1) * t362 + (t668 + t717 / 0.2e1) * t506;
t459 = (mrSges(7,3) * t681 + t730) * t507 - (mrSges(7,3) * t679 - t739) * t230 + t549 * t640 + t751;
t6 = t459 - t497;
t493 = t6 * qJD(1) + t27 * qJD(3);
t467 = (t230 * t408 + t409 * t507) * t689 + (t324 * t578 - t326 * t449) * t541;
t479 = m(6) * t504 + m(7) * t666;
t36 = t731 + t241 - t467 + t479;
t475 = (-t449 * t481 - t482 * t578) * t687;
t484 = m(7) * (-t362 * t408 + t409 * t506);
t466 = t484 / 0.2e1 + t475 / 0.2e1;
t480 = t753 + t370;
t483 = m(6) * t538 + m(7) * t647;
t71 = -t466 + t480 + t483;
t492 = qJD(1) * t36 + qJD(3) * t71;
t468 = (-t324 * t481 - t326 * t482) * t692 + (t230 * t506 + t362 * t507) * t690;
t63 = t498 + t468;
t491 = qJD(1) * t63;
t69 = 0.2e1 * t715 * mrSges(7,2) + 0.2e1 * t678;
t489 = qJD(1) * t69 - qJD(3) * t753;
t474 = t386 * t744 + t706 * t705;
t2 = t696 - t560 / 0.4e1 + t558 / 0.4e1 + t371 * t504 + t115 * t530 + t99 * t531 + t368 * t536 + t366 * t537 + t703 * mrSges(6,3) - (t437 + t704) * t563 / 0.4e1 - (t244 + t188) * t481 / 0.4e1 + t628 * t698 + t243 * t538 + ((t444 * t562 + t571) * pkin(4) + t474) * t691 + t453 * t660 + t456 * t661 + t562 * t637 + t562 * t638 + t508 * t640 - t740 + t706 * t529;
t23 = -pkin(3) * t434 + t741 + t734 * t655 + t717 * t657 + t264 * t726 + (t704 / 0.2e1 + t437 / 0.2e1) * t456 - (t653 - t374 / 0.2e1) * t482 - (t372 / 0.2e1 + t650) * t481 + (t438 / 0.2e1 - t435 / 0.2e1 + pkin(4) * t371) * t453 + t506 * t668 + (m(6) * t630 + t370) * t444 + (m(7) * t390 + t261) * t393;
t478 = t2 * qJD(1) + t23 * qJD(3);
t462 = (t231 * t408 + t234 * t409) * t689 + (t325 * t578 + t327 * t449) * t541 + mrSges(5,1) * t518 + mrSges(5,2) * t564 / 0.2e1 + t735;
t10 = t462 + t500 + t737;
t477 = t10 * qJD(1);
t471 = t238 * t689 + t305 * t691 - t735;
t15 = t471 + t501 + t736;
t42 = (t362 ^ 2 + t506 ^ 2) * mrSges(7,3) + (t481 ^ 2 + t482 ^ 2) * mrSges(6,3) + m(7) * (-t202 * t506 - t362 * t728) + m(6) * (-t386 * t481 - t482 * t705);
t476 = -qJD(1) * t15 + qJD(3) * t42;
t28 = (t681 + t680) * mrSges(7,2) + (t679 - t202 / 0.2e1) * mrSges(7,1);
t465 = (-t230 * t645 + t408 * t715) * mrSges(7,3) + t177 * t646 + t179 * t645 - t516;
t9 = (t685 + t59 / 0.2e1) * mrSges(7,2) + (-t48 / 0.2e1 - t58 / 0.2e1) * mrSges(7,1) + t465 + t516;
t473 = t9 * qJD(1) + t28 * qJD(3) - t167 * qJD(4);
t120 = t466 + t483;
t70 = -t722 / 0.2e1 + t678;
t66 = t467 + t479;
t64 = t498 - t468;
t17 = t470 + t556;
t16 = t471 + t500 - t736;
t13 = t727 + t327 * t530 + t325 * t529 + t234 * t533 + t231 * t532 + t460 - t461 + (t261 + t371) * t639 + t423 * t698;
t11 = t462 + t501 - t737;
t8 = -t621 / 0.2e1 - t622 / 0.2e1 + t618 / 0.2e1 - t617 / 0.2e1 + t465 - t516;
t5 = t459 + t497;
t1 = (-(t114 / 0.2e1 + t100 / 0.2e1) * t482 - (t115 / 0.2e1 - t99 / 0.2e1) * t481 + t703) * mrSges(6,3) + (-t604 / 0.4e1 + t660 - t297 / 0.4e1 + t366 * t688 + pkin(4) * t243 / 0.2e1) * t453 + (pkin(4) * t571 + t474) * t691 + (t661 + t299 / 0.4e1 + t368 * t688) * t456 - (t244 / 0.4e1 + t188 / 0.4e1) * t481 + ((-t437 / 0.4e1 - t704 / 0.4e1 + mrSges(5,3) * t537) * t453 + (t637 + t638 + mrSges(5,3) * t536 + (m(6) * t636 + t371 / 0.2e1) * pkin(4)) * t456) * t425 + t445 * t640 + t696 + t740;
t19 = [qJD(2) * t14 + qJD(3) * t3 + qJD(4) * t4 + qJD(5) * t21 + qJD(6) * t7, t575 + 0.2e1 * ((t231 * t506 + t234 * t362) * t689 + (-t325 * t481 + t327 * t482) * t691) * qJD(2) + t13 * qJD(3) + t11 * qJD(4) + t64 * qJD(5) + t17 * qJD(6), t13 * qJD(2) + t1 * qJD(4) + t16 * qJD(5) + t5 * qJD(6) + t592 + (t444 * t242 - Ifges(4,6) * t425 + t386 * t281 + t390 * t123 + t377 * mrSges(4,2) - t378 * mrSges(4,1) + t305 * t371 + pkin(3) * t348 + t238 * t261 - t378 * t496 + (Ifges(5,5) * t453 + Ifges(6,5) * t482 + Ifges(7,5) * t362 + Ifges(5,6) * t456 - Ifges(6,6) * t481 + Ifges(7,6) * t506) * t639 + t705 * t283 + 0.2e1 * (t112 * t705 + t113 * t386 + t305 * t444) * t691 + 0.2e1 * (-pkin(3) * t378 + pkin(8) * t720) * t693 + t720 * mrSges(5,3) + (-Ifges(4,5) - t456 * t437 / 0.2e1 + t435 * t635) * t423 - t367 * t627 - t54 * t585 + t234 * t668 + t264 * t676 + t187 * t642 + t185 * t644 + t327 * t650 + t325 * t653 + t108 * t655 + t106 * t657 + t365 * t626 + t296 * t634 + t298 * t635 - t202 * t176 + 0.2e1 * (-t202 * t55 + t238 * t390 + t54 * t728) * t689 - t112 * t582 - t113 * t583 + t55 * t588 + t728 * t178) * qJD(3), t584 + t11 * qJD(2) + t1 * qJD(3) + (-Ifges(5,5) * t563 - Ifges(5,6) * t562 + m(7) * (t408 * t58 + t409 * t59) + t230 * t613 - t507 * t614 + t618 - t617 + (t114 * t578 + t115 * t449) * t687 + t114 * mrSges(6,1) - t115 * mrSges(6,2) - t268 * mrSges(5,2) - t269 * mrSges(5,1) + t326 * t503 + t324 * t542 + t701) * qJD(4) + t66 * qJD(5) + t8 * qJD(6), qJD(2) * t64 + qJD(3) * t16 + qJD(4) * t66 + qJD(6) * t70 + t576, t579 + t17 * qJD(2) + t5 * qJD(3) + t8 * qJD(4) + t70 * qJD(5) + (t555 - t621 - t622) * qJD(6); qJD(3) * t12 - qJD(4) * t10 - qJD(5) * t63 + qJD(6) * t18 - t575, 0, t577 (-t434 + t475 - t480 + t484) * qJD(4) - t747 - t477, -t491, -qJD(4) * t753 + t574 - t747; -qJD(2) * t12 + qJD(4) * t2 - qJD(5) * t15 + qJD(6) * t6 - t592, -t577, qJD(4) * t23 + qJD(5) * t42 + qJD(6) * t27 (-t362 * t613 - t506 * t614 + t481 * t503 - t482 * t542 + m(7) * (t202 * t408 + t409 * t728) + (-t386 * t578 + t449 * t705) * t687 - t386 * mrSges(6,1) - t705 * mrSges(6,2) + t508 - t496 * pkin(8) + t702 + t752) * qJD(4) + t120 * qJD(5) + t754 + t478, qJD(4) * t120 + t476, t29 * qJD(4) + t493 + t754; qJD(2) * t10 - qJD(3) * t2 - qJD(5) * t36 + qJD(6) * t9 - t584, t477, -qJD(5) * t71 + qJD(6) * t28 - t478, -t543, -t492, t473 - t543; qJD(2) * t63 + qJD(3) * t15 + qJD(4) * t36 - qJD(6) * t69 - t576, t491, qJD(4) * t71 - t476 + t747, t492, 0, -t489; -qJD(2) * t18 - qJD(3) * t6 - qJD(4) * t9 + qJD(5) * t69 - t579, -t574, -t28 * qJD(4) - qJD(5) * t753 - t493, -t473, t489, 0;];
Cq  = t19;
