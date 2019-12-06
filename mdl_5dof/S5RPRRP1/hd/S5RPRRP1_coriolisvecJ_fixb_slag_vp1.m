% Calculate vector of centrifugal and Coriolis load on the joints for
% S5RPRRP1
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% qJD [5x1]
%   Generalized joint velocities
% pkin [7x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d1,d3,d4]';
% m_mdh [6x1]
%   mass of all robot links (including the base)
% rSges [6x3]
%   center of mass of all robot links (in body frames)
%   rows: links of the robot (starting with base)
%   columns: x-, y-, z-coordinates
% Icges [6x6]
%   inertia of all robot links about their respective center of mass, in body frames
%   rows: links of the robot (starting with base)
%   columns: xx, yy, zz, xy, xz, yz (see inertiavector2matrix.m)
% 
% Output:
% tauc [5x1]
%   joint torques required to compensate Coriolis and centrifugal load

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-05 18:00
% Revision: 77da58f92bca3eff71542919beafa37024070d86 (2019-12-05)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauc = S5RPRRP1_coriolisvecJ_fixb_slag_vp1(qJ, qJD, ...
  pkin, m, rSges, Icges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(7,1),zeros(6,1),zeros(6,3),zeros(6,6)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RPRRP1_coriolisvecJ_fixb_slag_vp1: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RPRRP1_coriolisvecJ_fixb_slag_vp1: qJD has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [7 1]), ...
  'S5RPRRP1_coriolisvecJ_fixb_slag_vp1: pkin has to be [7x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5RPRRP1_coriolisvecJ_fixb_slag_vp1: m has to be [6x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [6,3]), ...
  'S5RPRRP1_coriolisvecJ_fixb_slag_vp1: rSges has to be [6x3] (double)');
assert(isreal(Icges) && all(size(Icges) == [6 6]), ...
  'S5RPRRP1_coriolisvecJ_fixb_slag_vp1: Icges has to be [6x6] (double)'); 

%% Symbolic Calculation
% From coriolisvec_joint_fixb_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-05 17:59:14
% EndTime: 2019-12-05 17:59:51
% DurationCPUTime: 24.55s
% Computational Cost: add. (11021->641), mult. (16065->807), div. (0->0), fcn. (12434->6), ass. (0->354)
t736 = Icges(5,4) + Icges(6,4);
t735 = Icges(5,1) + Icges(6,1);
t734 = Icges(5,2) + Icges(6,2);
t368 = qJ(3) + qJ(4);
t344 = cos(t368);
t733 = t736 * t344;
t343 = sin(t368);
t732 = t736 * t343;
t731 = Icges(5,5) + Icges(6,5);
t730 = Icges(5,6) + Icges(6,6);
t729 = t734 * t344 + t732;
t723 = t735 * t343 + t733;
t727 = -t734 * t343 + t733;
t726 = t735 * t344 - t732;
t728 = Icges(5,3) + Icges(6,3);
t370 = sin(qJ(1));
t372 = cos(qJ(1));
t708 = t729 * t370 + t730 * t372;
t712 = -t730 * t370 + t729 * t372;
t721 = -t731 * t370 + t723 * t372;
t558 = t344 * t370;
t725 = t736 * t558;
t724 = t731 * t343 + t730 * t344;
t722 = t731 * t372;
t560 = t343 * t370;
t707 = t735 * t560 + t722 + t725;
t720 = -t730 * t343 + t731 * t344;
t719 = -t726 + t729;
t718 = t723 + t727;
t717 = t727 * t372 + t721;
t716 = -t726 * t372 + t712;
t715 = -t726 * t370 + t708;
t714 = t724 * t370 + t728 * t372;
t713 = -t728 * t370 + t724 * t372;
t710 = t726 * t343 + t727 * t344;
t632 = t721 * t343 + t712 * t344;
t662 = -t713 * t370 + t372 * t632;
t706 = t720 * t370;
t705 = t707 * t343 + t708 * t344;
t365 = qJD(3) + qJD(4);
t703 = t718 * t365;
t702 = t719 * t365;
t701 = -t721 * qJD(1) + t715 * t365;
t700 = t716 * t365 + (t723 * t370 + t722) * qJD(1);
t290 = t365 * t370;
t699 = -t712 * qJD(1) - t727 * t290 - t707 * t365;
t698 = -t708 * qJD(1) + t717 * t365;
t697 = t720 * t372;
t664 = -t713 * t372 - t712 * t558 - t560 * t721;
t696 = t710 * t370 + t697;
t695 = t710 * t372 - t706;
t694 = t714 * t370;
t692 = t713 * qJD(1);
t691 = t714 * qJD(1);
t665 = t714 * t372 + t708 * t558 + t707 * t560;
t663 = -t705 * t372 + t694;
t597 = pkin(4) * t343;
t369 = sin(qJ(3));
t599 = pkin(3) * t369;
t284 = t597 + t599;
t434 = rSges(6,1) * t343 + rSges(6,2) * t344;
t689 = -t284 - t434;
t688 = t710 * qJD(1) - t365 * t724;
t687 = t705 * qJD(1) + t706 * t365 + t692;
t686 = -t632 * qJD(1) - t697 * t365 + t691;
t685 = t700 * t343 - t698 * t344 + t692;
t684 = t701 * t343 + t699 * t344 + t691;
t683 = qJD(1) * t720 + t703 * t343 + t702 * t344;
t291 = t365 * t372;
t650 = -qJD(1) * t718 + t290 * t716 - t291 * t715;
t651 = (t734 * t560 - t707 - t725) * t291 + t717 * t290 + t719 * qJD(1);
t682 = -t650 * t343 + t651 * t344;
t371 = cos(qJ(3));
t681 = rSges(4,2) * t371;
t373 = -pkin(7) - pkin(6);
t364 = -qJ(5) + t373;
t678 = t372 * t284 + (-rSges(6,3) + t364) * t370;
t677 = t695 * qJD(1);
t676 = t696 * qJD(1) + t290 * t664;
t598 = pkin(3) * t371;
t675 = t687 * t370 + t684 * t372;
t674 = t686 * t370 - t685 * t372;
t673 = -t684 * t370 + t687 * t372;
t672 = t685 * t370 + t686 * t372;
t671 = t665 * t291 + t676;
t670 = t662 * t290 + t663 * t291 - t677;
t669 = t688 * t370 + t683 * t372;
t668 = -t683 * t370 + t688 * t372;
t667 = t699 * t343 - t701 * t344;
t666 = t698 * t343 + t700 * t344;
t661 = t708 * t343 - t344 * t707;
t660 = t712 * t343 - t344 * t721;
t270 = rSges(6,1) * t344 - rSges(6,2) * t343;
t232 = t270 * t372;
t502 = qJD(3) * t371;
t492 = pkin(3) * t502;
t559 = t344 * t365;
t243 = pkin(4) * t559 + t492;
t360 = t372 * rSges(6,3);
t500 = qJD(5) * t370;
t501 = qJD(3) * t372;
t322 = t501 * t598;
t505 = qJD(1) * t373;
t328 = t372 * t505;
t521 = t322 + t328;
t577 = -t365 * t232 - t243 * t372 + t500 + t521 + (-t364 * t372 + t360 + (-t599 - t689) * t370) * qJD(1);
t496 = t372 * t599;
t552 = t370 * t373;
t548 = -t434 * t372 + t496 + t552 - t678;
t649 = -qJD(1) * t724 - t697 * t290 + t706 * t291;
t506 = qJD(1) * t372;
t646 = t344 * t290 + t343 * t506;
t553 = t370 * t371;
t329 = Icges(4,4) * t553;
t557 = t369 * t370;
t569 = Icges(4,5) * t372;
t211 = Icges(4,1) * t557 + t329 + t569;
t574 = Icges(4,4) * t371;
t433 = Icges(4,1) * t369 + t574;
t212 = -Icges(4,5) * t370 + t372 * t433;
t297 = -Icges(4,2) * t369 + t574;
t248 = t297 * t372;
t397 = t370 * (t212 + t248) - t372 * (-Icges(4,2) * t557 + t211 + t329);
t575 = Icges(4,4) * t369;
t430 = Icges(4,2) * t371 + t575;
t209 = Icges(4,6) * t372 + t370 * t430;
t210 = -Icges(4,6) * t370 + t372 * t430;
t299 = Icges(4,1) * t371 - t575;
t249 = t299 * t370;
t250 = t299 * t372;
t398 = t370 * (t210 - t250) - t372 * (t209 - t249);
t645 = -t398 * t369 + t397 * t371;
t525 = t297 + t433;
t526 = -t430 + t299;
t644 = (t369 * t525 - t371 * t526) * qJD(1);
t643 = 0.2e1 * qJD(3);
t271 = rSges(5,1) * t344 - rSges(5,2) * t343;
t231 = t271 * t370;
t435 = rSges(5,1) * t343 + rSges(5,2) * t344;
t361 = t372 * rSges(5,3);
t187 = rSges(5,1) * t560 + rSges(5,2) * t558 + t361;
t333 = pkin(3) * t557;
t590 = pkin(6) + t373;
t237 = -t372 * t590 + t333;
t307 = t372 * pkin(1) + t370 * qJ(2);
t594 = pkin(6) * t372;
t466 = t307 + t594;
t442 = t237 + t466;
t342 = qJD(2) * t372;
t520 = t322 + t342;
t73 = -t271 * t291 + (t187 + t442) * qJD(1) - t520;
t642 = (qJD(1) * t231 + t291 * t435) * t73;
t417 = t210 * t371 + t212 * t369;
t637 = t417 * t372;
t236 = t370 * t590 + t496;
t346 = t372 * qJ(2);
t303 = pkin(1) * t370 - t346;
t283 = qJD(1) * t303;
t636 = qJD(1) * t236 - t283;
t635 = rSges(6,1) * t560 + rSges(6,2) * t558 + t370 * t284 + t360;
t595 = pkin(6) * t370;
t467 = -t303 - t595;
t464 = -rSges(3,2) * t372 + t370 * rSges(3,3);
t634 = t307 + t464;
t503 = qJD(3) * t370;
t476 = t369 * t503;
t633 = pkin(3) * t476;
t340 = qJD(5) * t372;
t480 = t344 * t506;
t507 = qJD(1) * t370;
t631 = t646 * rSges(6,1) + rSges(6,2) * t480 + t370 * t243 + t284 * t506 + t364 * t507 + t340;
t315 = pkin(4) * t558;
t292 = qJD(1) * t315;
t561 = t343 * t365;
t495 = pkin(4) * t561;
t554 = t370 * t270;
t630 = t270 * t507 + t372 * t495 + t292 - qJD(1) * t554 - (t434 + t597) * t291;
t205 = t434 * t365;
t293 = pkin(4) * t480;
t490 = t290 * t343;
t629 = pkin(4) * t490 - t370 * t205 + t270 * t506 + t290 * t434 + t293;
t628 = t291 * t232 + t577 * t372;
t124 = t210 * t369 - t212 * t371;
t129 = qJD(1) * t209 - qJD(3) * t248;
t131 = -qJD(3) * t250 + (t370 * t433 + t569) * qJD(1);
t427 = Icges(4,5) * t369 + Icges(4,6) * t371;
t208 = -Icges(4,3) * t370 + t372 * t427;
t511 = qJD(1) * t208;
t627 = qJD(3) * t124 + t129 * t371 + t131 * t369 + t511;
t276 = t430 * qJD(3);
t277 = t433 * qJD(3);
t295 = Icges(4,5) * t371 - Icges(4,6) * t369;
t626 = qJD(1) * t295 + qJD(3) * (t297 * t369 - t299 * t371) + t276 * t371 + t277 * t369;
t130 = qJD(1) * t210 + t297 * t503;
t132 = qJD(1) * t212 + qJD(3) * t249;
t418 = t209 * t369 - t211 * t371;
t207 = Icges(4,3) * t372 + t370 * t427;
t512 = qJD(1) * t207;
t625 = qJD(3) * t418 - t130 * t371 - t132 * t369 + t512;
t623 = t665 - t662;
t366 = t370 ^ 2;
t273 = t365 * t507;
t612 = -t273 / 0.2e1;
t274 = qJD(1) * t291;
t610 = t274 / 0.2e1;
t609 = -t290 / 0.2e1;
t608 = t290 / 0.2e1;
t607 = -t291 / 0.2e1;
t606 = t291 / 0.2e1;
t605 = t370 / 0.2e1;
t604 = -t372 / 0.2e1;
t603 = t372 / 0.2e1;
t602 = rSges(3,2) - pkin(1);
t601 = -rSges(5,3) - pkin(1);
t600 = -rSges(6,3) - pkin(1);
t596 = pkin(4) * t344;
t593 = pkin(6) * qJD(1) ^ 2;
t592 = -qJD(1) / 0.2e1;
t591 = qJD(1) / 0.2e1;
t588 = rSges(3,3) * t372;
t233 = t271 * t372;
t118 = -t365 * t233 + (t370 * t435 + t361) * qJD(1);
t206 = t435 * t365;
t374 = qJD(3) ^ 2;
t498 = qJD(1) * qJD(2);
t336 = t372 * t498;
t444 = -t372 * t593 + t336;
t447 = qJD(1) * t492;
t385 = -t333 * t374 + t372 * t447 + t444;
t154 = (t333 - t594) * qJD(1) - t521;
t238 = qJD(1) * t307 - t342;
t546 = -t154 - t238;
t53 = -t206 * t290 + t271 * t274 + (-t118 + t546) * qJD(1) + t385;
t585 = t370 * t53;
t306 = rSges(4,1) * t371 - rSges(4,2) * t369;
t252 = t306 * t372;
t362 = t372 * rSges(4,3);
t436 = rSges(4,1) * t369 + t681;
t135 = -qJD(3) * t252 + (t370 * t436 + t362) * qJD(1);
t280 = t436 * qJD(3);
t477 = t306 * t501;
t69 = -t280 * t503 + (-t135 - t238 + t477) * qJD(1) + t444;
t584 = t370 * t69;
t484 = t646 * rSges(5,1) + rSges(5,2) * t480;
t120 = (-rSges(5,2) * t561 - rSges(5,3) * qJD(1)) * t370 + t484;
t475 = t370 * t502;
t320 = pkin(3) * t475;
t479 = t369 * t506;
t482 = pkin(3) * t479 + t370 * t505 + t320;
t153 = pkin(6) * t507 + t482;
t341 = qJD(2) * t370;
t519 = qJ(2) * t506 + t341;
t532 = qJD(1) * (-pkin(1) * t507 + t519) + t370 * t498;
t412 = -t370 * t593 + t532;
t396 = qJD(1) * t153 + t370 * t447 + t374 * t496 + t412;
t52 = qJD(1) * t120 + t206 * t291 + t271 * t273 + t396;
t583 = t372 * t52;
t483 = t506 * t681 + (t475 + t479) * rSges(4,1);
t504 = qJD(3) * t369;
t136 = (-rSges(4,2) * t504 - rSges(4,3) * qJD(1)) * t370 + t483;
t255 = t306 * t503;
t68 = t280 * t501 + (t136 + t255) * qJD(1) + t412;
t582 = t372 * t68;
t357 = t370 * rSges(5,3);
t189 = t372 * t435 - t357;
t523 = t320 + t341;
t439 = t290 * t271 + t523;
t443 = t236 + t467;
t72 = (t189 + t443) * qJD(1) + t439;
t581 = t372 * t72;
t448 = t342 - t500;
t468 = t270 + t596;
t547 = t333 - (-t364 + t373) * t372 - t635;
t62 = -t322 - t468 * t291 + (t442 - t547) * qJD(1) - t448;
t580 = t62 * t205;
t579 = t72 * t271;
t494 = rSges(6,2) * t561;
t576 = (-rSges(6,3) * qJD(1) - t494) * t370 - t482 + t631;
t438 = -t236 * t501 - t237 * t503;
t43 = t290 * t547 + t291 * t548 + t438;
t566 = qJD(1) * t43;
t358 = t370 * rSges(4,3);
t216 = t436 * t372 - t358;
t101 = t255 + t341 + (t216 + t467) * qJD(1);
t565 = t101 * t372;
t562 = t295 * t372;
t245 = t370 * t295;
t549 = t548 * t372;
t537 = -t187 - t237;
t524 = t315 + t554;
t518 = rSges(3,2) * t507 + rSges(3,3) * t506;
t517 = t341 - t283;
t508 = qJD(1) * t427;
t414 = t371 * t297 + t369 * t299;
t126 = t372 * t414 - t245;
t499 = t126 * qJD(1);
t497 = -rSges(4,3) - pkin(1) - pkin(6);
t334 = pkin(3) * t553;
t493 = pkin(3) * t504;
t488 = -t237 + t547;
t83 = t372 * t207 + t209 * t553 + t211 * t557;
t84 = -t372 * t208 - t210 * t553 - t212 * t557;
t215 = rSges(4,1) * t557 + rSges(4,2) * t553 + t362;
t474 = -t507 / 0.2e1;
t473 = t506 / 0.2e1;
t472 = -t503 / 0.2e1;
t470 = -t501 / 0.2e1;
t455 = -t231 * t290 - t291 * t233;
t449 = qJD(1) * t233 - t290 * t435;
t285 = t596 + t598;
t102 = -t477 - t342 + (t215 + t466) * qJD(1);
t424 = t101 * t370 - t102 * t372;
t419 = t209 * t371 + t211 * t369;
t411 = (-t372 ^ 2 - t366) * t492;
t406 = (t370 * t84 + t372 * t83) * qJD(3);
t190 = t370 * t207;
t85 = -t419 * t372 + t190;
t86 = -t370 * t208 + t637;
t405 = (t370 * t86 + t372 * t85) * qJD(3);
t121 = (-t215 * t370 - t216 * t372) * qJD(3);
t404 = t435 + t599;
t391 = -qJD(1) * t417 - qJD(3) * t562 + t512;
t390 = qJD(1) * t419 + qJD(3) * t245 + t511;
t387 = t414 * qJD(1) - t427 * qJD(3);
t386 = t290 * t468 + t340 + t523;
t379 = t154 * t501 + (-t153 * t370 + (t236 * t370 - t237 * t372) * qJD(1)) * qJD(3);
t376 = (t664 * t370 + t665 * t372) * t612 + (t662 * t370 + t663 * t372) * t610 + (t649 * t370 + t682 * t372) * t609 + (t675 * t372 + t674 * t370 + (-t663 * t370 + t662 * t372) * qJD(1)) * t608 + (-t682 * t370 + t649 * t372) * t607 + (t673 * t372 + t672 * t370 + (-t665 * t370 + t664 * t372) * qJD(1)) * t606 + (t669 * qJD(1) - t663 * t273 + t662 * t274 + t674 * t290 + t675 * t291) * t605 + (t668 * qJD(1) - t665 * t273 + t664 * t274 + t672 * t290 + t673 * t291) * t603 + (t651 * t343 + t650 * t344) * t592 + (t667 * t372 + t666 * t370 + (t661 * t370 + t660 * t372) * qJD(1)) * t591 + t671 * t474 + t670 * t473;
t304 = rSges(3,2) * t370 + t588;
t251 = t306 * t370;
t196 = t372 * t236;
t195 = (-t285 + t598) * t372;
t194 = t285 * t370 - t334;
t166 = t372 * t189;
t158 = qJD(1) * t634 - t342;
t157 = t341 + (-t303 + t304) * qJD(1);
t143 = t372 * t154;
t134 = t336 + (-qJD(1) * t464 - t238) * qJD(1);
t133 = qJD(1) * t518 + t532;
t125 = t370 * t414 + t562;
t122 = t125 * qJD(1);
t104 = t372 * t118;
t67 = -t187 * t290 - t189 * t291 + t438;
t66 = -t370 * t626 + t387 * t372;
t65 = t387 * t370 + t372 * t626;
t64 = t417 * qJD(3) - t129 * t369 + t131 * t371;
t63 = -qJD(3) * t419 - t130 * t369 + t132 * t371;
t61 = (t443 - t548) * qJD(1) + t386;
t42 = t405 - t499;
t41 = t122 + t406;
t28 = -t205 * t290 + t270 * t274 + (t274 * t344 - t290 * t561) * pkin(4) + (-t500 + t546 - t577) * qJD(1) + t385;
t27 = t205 * t291 + t270 * t273 + (t273 * t344 + t291 * t561) * pkin(4) + (t340 + t576) * qJD(1) + t396;
t26 = t118 * t291 - t120 * t290 - t187 * t274 + t189 * t273 + t379;
t9 = -t273 * t548 + t274 * t547 - t290 * t576 + t291 * t577 + t379;
t1 = [(t122 + ((-t85 + t190 + t84) * t370 + (t86 - t637 + (t208 - t419) * t370 + t83) * t372) * qJD(3)) * t472 - t695 * t274 / 0.2e1 + t660 * t610 + ((t623 + t662) * t291 + t676) * t609 + (t499 + (t366 * t208 + (-t190 + t84 + (t208 + t419) * t372) * t372) * qJD(3) + t42) * t470 + (-qJD(3) * t414 + t276 * t369 - t277 * t371 + t702 * t343 - t703 * t344) * qJD(1) + (t28 * (-t303 + t678) + t61 * t448 + t27 * (t307 + t635) + t62 * (-rSges(6,2) * t490 + t519 + t631) + (t28 * t434 + t61 * (rSges(6,1) * t559 + t243 - t494) - t27 * t364) * t372 + (t61 * (t364 + t600) * t372 + (t61 * (-qJ(2) + t689) + t62 * t600) * t370) * qJD(1) - (-t61 + (-t548 - t595) * qJD(1) + t386 + t636) * t62) * m(6) + (t53 * (-t303 - t357 + t552) + t72 * (t328 + t520) + t52 * (t333 + t187 + t307) + t73 * (-rSges(5,2) * t490 + t482 + t484 + t519) + (t365 * t579 - t52 * t373 + t404 * t53) * t372 + (t601 * t581 + (t72 * (-qJ(2) - t404) + t73 * t601) * t370) * qJD(1) - (-t72 + (t189 - t595) * qJD(1) + t439 + t636) * t73) * m(5) + (t69 * (-t358 + t467) + t101 * t342 + t68 * (t215 + t307) + t102 * (-rSges(4,2) * t476 + t483 + t519) + (qJD(3) * t101 * t306 + t68 * pkin(6) + t436 * t69) * t372 + (t497 * t565 + (t101 * (-qJ(2) - t436) + t102 * t497) * t370) * qJD(1) - (-t101 + t255 + (t216 - t595) * qJD(1) + t517) * t102) * m(4) + (t134 * (t370 * t602 + t346 + t588) + t157 * t342 + t133 * t634 + t158 * (t518 + t519) + (t157 * t602 * t372 + (t157 * (-rSges(3,3) - qJ(2)) - t158 * pkin(1)) * t370) * qJD(1) - (qJD(1) * t304 - t157 + t517) * t158) * m(3) + (t41 + t64 + t65) * t503 / 0.2e1 + (qJD(1) * t124 + t63 + t66) * t501 / 0.2e1 + (-t661 + t696) * t612 + (((t705 + t713) * t372 + t632 * t370 + t664 - t694) * t291 + (t623 - t665) * t290 + t670 + t677) * t607 + (t667 + t668) * t606 + (t666 + t669 + t671) * t608 + (t372 * t126 + (-t418 + t125) * t370) * qJD(3) * t592; 0.2e1 * (t27 * t604 + t28 * t605) * m(6) + 0.2e1 * (t585 / 0.2e1 - t583 / 0.2e1) * m(5) + 0.2e1 * (t584 / 0.2e1 - t582 / 0.2e1) * m(4) + 0.2e1 * (t133 * t604 + t134 * t605) * m(3); ((-t369 * t526 - t371 * t525) * qJD(1) + (t369 * t397 + t371 * t398) * qJD(3)) * t592 + t376 + ((t245 * t501 - t508) * t372 + (-t644 + (-t372 * t562 - t645) * qJD(3)) * t370) * t470 + ((-t503 * t562 - t508) * t370 + (t644 + (t370 * t245 + t645) * qJD(3)) * t372) * t472 + (t370 * t64 + t372 * t63 + (t124 * t372 + t370 * t418) * qJD(1)) * t591 + (qJD(1) * t65 + ((t390 * t370 + t372 * t625) * t372 + (t391 * t370 - t372 * t627) * t370 + (-t85 * t370 + t86 * t372) * qJD(1)) * t643) * t605 + (qJD(1) * t66 + ((-t370 * t625 + t390 * t372) * t372 + (t370 * t627 + t391 * t372) * t370 + (-t83 * t370 + t84 * t372) * qJD(1)) * t643) * t603 + (t41 + t406) * t474 + (t42 + t405) * t473 + (t28 * (t334 + t524) + t9 * (-t196 + t549) + (t27 * (-t270 - t285) + t580 + t488 * t566) * t372 + (t9 * t488 + (t236 - t548) * t566) * t370 + (-qJD(1) * t194 + t630) * t62 + ((-t493 - t495) * t370 - (-t195 + t232) * qJD(1) + t629 + t633) * t61 + (t143 + (-t153 - t576) * t370 - t195 * t291 - (-t194 - t554) * t290 - t411 + t628) * t43) * m(6) + (-t642 - t67 * (t411 + t455) + t53 * t334 + t26 * (-t166 - t196) + t67 * (t104 + t143) + (t52 * (-t271 - t598) + t73 * t206 + (t537 * t67 + t579) * qJD(1)) * t372 + (t53 * t271 + t26 * t537 + t67 * (-t120 - t153) + (t73 * t271 + t67 * (t189 + t236)) * qJD(1)) * t370 + (-t449 + (-t206 - t493) * t370 + t633) * t72) * m(5) + (-(t101 * t252 + t102 * t251) * qJD(1) - (t121 * (-t251 * t370 - t252 * t372) - t424 * t436) * qJD(3) + 0.2e1 * t121 * (t135 * t372 - t136 * t370 + (-t215 * t372 + t216 * t370) * qJD(1)) - t424 * t280 + (t584 - t582 + (t102 * t370 + t565) * qJD(1)) * t306) * m(4); t376 + (t28 * t524 + (-t27 * t468 + t547 * t566 + t580) * t372 + (t370 * t547 + t549) * t9 + (-t292 + t630) * t62 + (-qJD(1) * t232 - t370 * t495 - t293 + t629) * t61 + (t554 * t290 - (-t290 * t370 - t291 * t372) * t596 + (-qJD(1) * t548 - t576) * t370 + t628) * t43) * m(6) + (-t449 * t72 - t642 + t26 * (-t187 * t370 - t166) - (t370 * t72 - t372 * t73) * t206 + (t585 - t583 + (t370 * t73 + t581) * qJD(1)) * t271 + (-t455 - t120 * t370 + t104 + (-t187 * t372 + t189 * t370) * qJD(1)) * t67) * m(5); 0.2e1 * (t27 * t605 + t28 * t603 - t43 * (-t290 * t372 + t291 * t370) / 0.2e1) * m(6);];
tauc = t1(:);
