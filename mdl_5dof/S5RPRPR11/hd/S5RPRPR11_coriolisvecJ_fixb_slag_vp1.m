% Calculate vector of centrifugal and Coriolis load on the joints for
% S5RPRPR11
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% qJD [5x1]
%   Generalized joint velocities
% pkin [8x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d1,d3,d5,theta2]';
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
% Datum: 2019-12-31 18:28
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauc = S5RPRPR11_coriolisvecJ_fixb_slag_vp1(qJ, qJD, ...
  pkin, m, rSges, Icges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(8,1),zeros(6,1),zeros(6,3),zeros(6,6)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RPRPR11_coriolisvecJ_fixb_slag_vp1: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RPRPR11_coriolisvecJ_fixb_slag_vp1: qJD has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S5RPRPR11_coriolisvecJ_fixb_slag_vp1: pkin has to be [8x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5RPRPR11_coriolisvecJ_fixb_slag_vp1: m has to be [6x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [6,3]), ...
  'S5RPRPR11_coriolisvecJ_fixb_slag_vp1: rSges has to be [6x3] (double)');
assert(isreal(Icges) && all(size(Icges) == [6 6]), ...
  'S5RPRPR11_coriolisvecJ_fixb_slag_vp1: Icges has to be [6x6] (double)'); 

%% Symbolic Calculation
% From coriolisvec_joint_fixb_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 18:27:11
% EndTime: 2019-12-31 18:27:45
% DurationCPUTime: 27.28s
% Computational Cost: add. (15343->740), mult. (21390->943), div. (0->0), fcn. (19371->8), ass. (0->385)
t733 = Icges(5,4) + Icges(4,5);
t389 = sin(qJ(1));
t390 = cos(qJ(1));
t382 = pkin(8) + qJ(3);
t362 = sin(t382);
t352 = Icges(5,5) * t362;
t363 = cos(t382);
t455 = Icges(5,1) * t363 + t352;
t732 = -t389 * t455 + t733 * t390;
t571 = t362 * t389;
t334 = Icges(4,4) * t571;
t569 = t363 * t389;
t722 = Icges(4,1) * t569 - t334 - t732;
t452 = Icges(5,3) * t363 - t352;
t595 = Icges(4,4) * t362;
t731 = Icges(4,2) * t363 + t452 + t595;
t590 = Icges(5,5) * t363;
t294 = Icges(5,1) * t362 - t590;
t353 = Icges(4,4) * t363;
t730 = Icges(4,1) * t362 + t294 + t353;
t289 = Icges(4,5) * t363 - Icges(4,6) * t362;
t214 = Icges(4,3) * t389 + t289 * t390;
t291 = Icges(5,4) * t363 + Icges(5,6) * t362;
t216 = Icges(5,2) * t389 + t291 * t390;
t729 = t214 + t216;
t220 = Icges(5,4) * t389 + t390 * t455;
t297 = Icges(4,1) * t363 - t595;
t222 = Icges(4,5) * t389 + t297 * t390;
t721 = t220 + t222;
t287 = Icges(5,3) * t362 + t590;
t453 = -Icges(4,2) * t362 + t353;
t728 = t287 - t453;
t720 = (Icges(4,6) - Icges(5,6)) * t363 + t733 * t362;
t727 = t297 + t455;
t581 = Icges(4,3) * t390;
t213 = Icges(4,5) * t569 - Icges(4,6) * t571 - t581;
t211 = -Icges(5,6) * t390 + t287 * t389;
t584 = Icges(4,6) * t390;
t217 = Icges(4,4) * t569 - Icges(4,2) * t571 - t584;
t577 = t217 * t362;
t693 = -t211 * t362 - t722 * t363 + t577;
t726 = -t213 * t390 - t389 * t693;
t724 = -t211 + t217;
t568 = t363 * t390;
t333 = Icges(5,5) * t568;
t570 = t362 * t390;
t583 = Icges(5,6) * t389;
t212 = Icges(5,3) * t570 + t333 + t583;
t218 = Icges(4,6) * t389 + t390 * t453;
t723 = -t212 + t218;
t719 = t731 * qJD(3);
t718 = t730 * qJD(3);
t710 = -t731 * t362 + t730 * t363;
t717 = t212 * t570 + t729 * t389 + t721 * t568;
t215 = -Icges(5,2) * t390 + t291 * t389;
t197 = t389 * t215;
t716 = -t211 * t570 - t389 * t213 - t722 * t568 - t197;
t715 = t728 * qJD(3);
t714 = t727 * qJD(3);
t713 = -t289 - t291;
t650 = t720 * t390;
t651 = t720 * t389;
t508 = qJD(3) - qJD(5);
t321 = t508 * t389;
t516 = qJD(3) * t390;
t322 = -qJD(5) * t390 + t516;
t388 = sin(qJ(5));
t617 = cos(qJ(5));
t494 = t362 * t617;
t239 = t388 * t569 - t389 * t494;
t435 = t362 * t388 + t363 * t617;
t240 = t435 * t389;
t118 = Icges(6,5) * t240 - Icges(6,6) * t239 + Icges(6,3) * t390;
t592 = Icges(6,4) * t240;
t121 = -Icges(6,2) * t239 + Icges(6,6) * t390 + t592;
t593 = Icges(6,4) * t239;
t125 = -Icges(6,1) * t240 - Icges(6,5) * t390 + t593;
t241 = t388 * t568 - t390 * t494;
t242 = t435 * t390;
t43 = -t118 * t389 - t241 * t121 - t125 * t242;
t120 = Icges(6,5) * t242 - Icges(6,6) * t241 - Icges(6,3) * t389;
t210 = Icges(6,4) * t242;
t123 = -Icges(6,2) * t241 - Icges(6,6) * t389 + t210;
t209 = Icges(6,4) * t241;
t126 = Icges(6,1) * t242 - Icges(6,5) * t389 - t209;
t469 = t120 * t389 + t241 * t123 - t242 * t126;
t282 = -t363 * t388 + t494;
t173 = Icges(6,5) * t282 - Icges(6,6) * t435;
t272 = Icges(6,4) * t282;
t176 = -Icges(6,2) * t435 + t272;
t271 = Icges(6,4) * t435;
t179 = Icges(6,1) * t282 - t271;
t55 = -t173 * t389 - t176 * t241 + t179 * t242;
t12 = t55 * qJD(1) - t321 * t469 - t322 * t43;
t520 = qJD(1) * t390;
t169 = t508 * t435;
t677 = t508 * t282;
t85 = Icges(6,5) * t169 + Icges(6,6) * t677;
t86 = Icges(6,4) * t169 + Icges(6,2) * t677;
t87 = Icges(6,1) * t169 + Icges(6,4) * t677;
t521 = qJD(1) * t389;
t97 = t169 * t390 - t282 * t521;
t98 = -t390 * t677 - t435 * t521;
t13 = -t173 * t520 + t176 * t97 + t179 * t98 - t241 * t86 + t242 * t87 - t389 * t85;
t100 = qJD(1) * t242 - t389 * t677;
t99 = t169 * t389 + t282 * t520;
t14 = t100 * t179 - t173 * t521 + t176 * t99 - t239 * t86 + t240 * t87 + t390 * t85;
t59 = Icges(6,4) * t100 + Icges(6,2) * t99 - Icges(6,6) * t521;
t61 = Icges(6,1) * t100 + Icges(6,4) * t99 - Icges(6,5) * t521;
t17 = t121 * t677 - t125 * t169 + t282 * t61 - t435 * t59;
t58 = Icges(6,4) * t98 + Icges(6,2) * t97 - Icges(6,6) * t520;
t60 = Icges(6,1) * t98 + Icges(6,4) * t97 - Icges(6,5) * t520;
t18 = t123 * t677 + t126 * t169 + t282 * t60 - t435 * t58;
t472 = qJD(1) * t508;
t306 = t389 * t472;
t307 = t390 * t472;
t409 = qJD(1) * (Icges(6,1) * t435 + t176 + t272) + t321 * (Icges(6,1) * t241 + t123 + t210) - t322 * (Icges(6,1) * t239 + t121 + t592);
t41 = t118 * t390 - t121 * t239 - t125 * t240;
t419 = qJD(1) * (-Icges(6,5) * t435 - Icges(6,6) * t282) + (Icges(6,5) * t239 + Icges(6,6) * t240) * t322 - (Icges(6,5) * t241 + Icges(6,6) * t242) * t321;
t42 = t390 * t120 - t239 * t123 + t240 * t126;
t54 = t173 * t390 - t176 * t239 + t179 * t240;
t57 = Icges(6,5) * t100 + Icges(6,6) * t99 - Icges(6,3) * t521;
t6 = -t118 * t520 + t121 * t97 - t125 * t98 - t241 * t59 + t242 * t61 - t389 * t57;
t614 = -qJD(1) / 0.2e1;
t622 = t322 / 0.2e1;
t66 = -t121 * t435 - t125 * t282;
t559 = Icges(6,2) * t282 - t179 + t271;
t631 = -t321 * (-Icges(6,2) * t242 + t126 - t209) + t322 * (-Icges(6,2) * t240 - t125 - t593);
t666 = -qJD(1) * t559 - t631;
t668 = t43 / 0.2e1 + t42 / 0.2e1;
t67 = -t123 * t435 + t126 * t282;
t56 = Icges(6,5) * t98 + Icges(6,6) * t97 - Icges(6,3) * t520;
t7 = -t120 * t520 + t123 * t97 + t126 * t98 - t241 * t58 + t242 * t60 - t389 * t56;
t8 = -t100 * t125 - t118 * t521 + t121 * t99 - t239 * t59 + t240 * t61 + t390 * t57;
t9 = t100 * t126 - t120 * t521 + t123 * t99 - t239 * t58 + t240 * t60 + t390 * t56;
t709 = (-t17 * t390 + t18 * t389 + (t389 * t66 + t390 * t67) * qJD(1)) * t614 - t389 * (qJD(1) * t13 + t321 * t7 - t322 * t6) / 0.2e1 + t390 * (qJD(1) * t14 + t321 * t9 - t322 * t8) / 0.2e1 - (qJD(1) * t54 + t321 * t42 - t322 * t41) * t521 / 0.2e1 - t12 * t520 / 0.2e1 + (t389 * t469 + t390 * t668) * t307 + (-t389 * t668 + t390 * t41) * t306 - (-t241 * t666 - t242 * t409 + (-t469 * qJD(1) - t6) * t390 + (t43 * qJD(1) - t419 + t7) * t389) * t321 / 0.2e1 + (-t239 * t666 - t240 * t409 + (t41 * qJD(1) + t9) * t389 + (t42 * qJD(1) + t419 - t8) * t390) * t622;
t567 = t390 * t215;
t684 = -t567 + t726;
t683 = -t217 * t570 - t716;
t682 = -t218 * t570 + t717;
t466 = -t212 * t571 + t216 * t390 - t220 * t569;
t187 = t222 * t569;
t475 = t390 * t214 - t187;
t78 = -t218 * t571 - t475;
t708 = -t466 + t78;
t707 = t409 * t282;
t706 = t710 * t389 - t650;
t705 = t710 * t390 + t651;
t704 = t719 * t390 + (t389 * t453 - t211 - t584) * qJD(1);
t703 = t719 * t389 + (t287 * t390 - t218 + t583) * qJD(1);
t702 = -t718 * t390 + (-t297 * t389 + t732) * qJD(1);
t701 = -t721 * qJD(1) + t718 * t389;
t700 = t720 * qJD(3);
t576 = t218 * t362;
t699 = t212 * t362 + t721 * t363 - t576;
t654 = t721 * t362 + t723 * t363;
t655 = t722 * t362 + t724 * t363;
t696 = t389 * (-t390 * t731 + t721) - t390 * (-Icges(4,2) * t569 - t452 * t389 - t334 + t722);
t695 = t714 * t363 + t715 * t362 + (-t362 * t730 - t363 * t731) * qJD(3) + t720 * qJD(1);
t694 = t724 * t390 + (-Icges(5,1) * t570 + t294 * t390 + t333 - t723) * t389;
t692 = t729 * qJD(1);
t691 = t730 - t728;
t690 = -t731 + t727;
t687 = t710 * qJD(1) + t713 * qJD(3);
t181 = -rSges(6,1) * t435 - rSges(6,2) * t282;
t679 = t181 * t321;
t678 = t181 * t322;
t676 = t705 * qJD(1);
t460 = rSges(6,1) * t240 - rSges(6,2) * t239;
t127 = rSges(6,3) * t390 + t460;
t298 = pkin(4) * t569 + pkin(7) * t390;
t562 = t127 + t298;
t675 = -qJD(3) * t654 + t362 * t704 + t363 * t702 + t692;
t644 = qJD(1) * t215;
t674 = -qJD(1) * t213 + qJD(3) * t655 - t362 * t703 + t363 * t701 - t644;
t673 = (t389 * t682 - t390 * t683) * qJD(3);
t672 = (t708 * t389 - t390 * t684) * qJD(3);
t671 = t706 * qJD(1);
t670 = qJD(1) * t693 - t389 * t700 + t692;
t669 = -t644 - t700 * t390 + (-t289 * t389 + t581 - t699) * qJD(1);
t162 = rSges(6,1) * t239 + rSges(6,2) * t240;
t164 = rSges(6,1) * t241 + rSges(6,2) * t242;
t665 = -t162 * t321 - t164 * t322;
t664 = 0.2e1 * qJD(3);
t663 = t671 + t672;
t662 = t673 + t676;
t661 = -t687 * t389 + t390 * t695;
t660 = t389 * t695 + t687 * t390;
t659 = t693 * qJD(3) + t362 * t701 + t363 * t703;
t658 = qJD(3) * t699 + t362 * t702 - t363 * t704;
t605 = t362 * rSges(5,1);
t301 = -rSges(5,3) * t363 + t605;
t351 = qJD(4) * t362;
t656 = (qJD(3) * t301 - t351) * t389;
t368 = t390 * qJ(2);
t326 = pkin(1) * t389 - t368;
t386 = cos(pkin(8));
t354 = pkin(2) * t386 + pkin(1);
t387 = -pkin(6) - qJ(2);
t357 = t390 * t387;
t531 = -t389 * t354 - t357;
t203 = t326 + t531;
t309 = qJD(1) * t326;
t653 = qJD(1) * t203 - t309;
t375 = t389 * rSges(4,3);
t236 = rSges(4,1) * t568 - rSges(4,2) * t570 + t375;
t336 = t390 * t354;
t473 = -t387 * t389 + t336;
t652 = t236 + t473;
t367 = t389 * qJ(2);
t328 = t390 * pkin(1) + t367;
t608 = rSges(3,2) * sin(pkin(8));
t610 = rSges(3,1) * t386;
t440 = t389 * rSges(3,3) + (-t608 + t610) * t390;
t649 = t328 + t440;
t341 = pkin(4) * t568;
t299 = -pkin(7) * t389 + t341;
t648 = -t362 * t696 + t694 * t363;
t647 = (-t362 * t691 + t363 * t690) * qJD(1);
t646 = t567 + t717;
t645 = t713 * qJD(1);
t303 = pkin(3) * t363 + qJ(4) * t362;
t265 = t303 * t389;
t243 = qJD(1) * t265;
t638 = -t243 + t653;
t629 = m(5) / 0.2e1;
t628 = m(6) / 0.2e1;
t621 = t389 / 0.2e1;
t620 = -t390 / 0.2e1;
t619 = -rSges(5,1) - pkin(3);
t618 = rSges(6,3) + pkin(7);
t616 = pkin(4) * t363;
t613 = qJD(1) / 0.2e1;
t612 = pkin(1) - t354;
t611 = t98 * rSges(6,1) + t97 * rSges(6,2);
t609 = rSges(4,1) * t363;
t302 = rSges(4,1) * t362 + rSges(4,2) * t363;
t268 = t302 * t390;
t366 = qJD(2) * t390;
t517 = qJD(3) * t389;
t492 = t302 * t517;
t84 = qJD(1) * t652 - t366 - t492;
t607 = t268 * t84;
t377 = t389 * rSges(5,2);
t234 = rSges(4,1) * t569 - rSges(4,2) * t571 - t390 * rSges(4,3);
t365 = qJD(2) * t389;
t491 = t302 * t516;
t464 = t365 - t491;
t555 = t203 - t326;
t83 = (-t234 + t555) * qJD(1) + t464;
t604 = t389 * t83;
t603 = -rSges(5,3) - qJ(4);
t421 = -t362 * t516 - t363 * t521;
t199 = pkin(4) * t421 - pkin(7) * t520;
t62 = -rSges(6,3) * t520 + t611;
t602 = t199 + t62;
t490 = t362 * t517;
t318 = pkin(4) * t490;
t200 = qJD(1) * t299 - t318;
t465 = rSges(6,1) * t100 + rSges(6,2) * t99;
t63 = -rSges(6,3) * t521 + t465;
t601 = t200 + t63;
t549 = t242 * rSges(6,1) - t241 * rSges(6,2);
t129 = -rSges(6,3) * t389 + t549;
t561 = t129 + t299;
t283 = qJD(1) * t328 - t366;
t348 = t387 * t521;
t556 = t348 - (-t390 * t612 - t367) * qJD(1) - t283;
t515 = qJD(4) * t363;
t232 = qJD(3) * t303 - t515;
t304 = rSges(5,1) * t363 + rSges(5,3) * t362;
t547 = -t304 * qJD(3) - t232;
t235 = rSges(5,1) * t568 + rSges(5,3) * t570 + t377;
t342 = pkin(3) * t568;
t269 = qJ(4) * t570 + t342;
t546 = -t235 - t269;
t545 = t389 * t265 + t390 * t269;
t300 = pkin(3) * t362 - qJ(4) * t363;
t266 = t300 * t390;
t514 = qJD(4) * t389;
t544 = -qJD(1) * t266 + t363 * t514;
t543 = -t269 - t299;
t510 = qJD(1) * qJD(2);
t358 = qJ(2) * t520;
t527 = t358 + t365;
t542 = qJD(1) * (-pkin(1) * t521 + t527) + t389 * t510;
t273 = t300 * t517;
t541 = -t273 - t366;
t536 = -t300 - t301;
t535 = -t303 - t304;
t489 = t363 * t516;
t534 = rSges(5,2) * t520 + rSges(5,3) * t489;
t493 = t362 * t521;
t533 = rSges(4,2) * t493 + rSges(4,3) * t520;
t513 = qJD(4) * t390;
t324 = t362 * t513;
t532 = t324 + t365;
t349 = t389 * t608;
t530 = rSges(3,3) * t520 + qJD(1) * t349;
t529 = t348 + t366;
t528 = t390 * rSges(3,3) + t349;
t526 = t389 ^ 2 + t390 ^ 2;
t519 = qJD(3) * t362;
t518 = qJD(3) * t363;
t509 = qJD(3) * qJD(4);
t507 = -t387 - t618;
t506 = t389 * t610;
t313 = qJ(4) * t489;
t130 = pkin(3) * t421 - qJ(4) * t493 + t313 + t324;
t247 = t362 * t520 + t363 * t517;
t319 = pkin(3) * t490;
t487 = t362 * t514;
t131 = qJ(4) * t247 + qJD(1) * t342 - t319 + t487;
t504 = t390 * t130 + t389 * t131 + t265 * t520;
t503 = -t131 + t556;
t502 = qJD(1) * (-t358 + (t389 * t612 - t357) * qJD(1)) + t542;
t501 = -t265 + t555;
t500 = t269 + t473;
t262 = t300 * t389;
t499 = -t262 * t517 - t266 * t516 + t351;
t356 = t390 * t510;
t485 = t363 * t509;
t498 = qJD(1) * t273 + t390 * t485 + t356;
t497 = t313 + t532;
t496 = t319 + t529;
t495 = t336 + t269;
t486 = -pkin(1) - t610;
t482 = -t517 / 0.2e1;
t479 = t516 / 0.2e1;
t478 = -pkin(4) * t362 - t300;
t476 = t390 * t536;
t474 = -t213 + t576;
t471 = t131 * t517 + t362 * t509 + (t130 + t243) * t516;
t182 = rSges(6,1) * t282 - rSges(6,2) * t435;
t470 = -t182 + t478;
t463 = t265 * t517 + t269 * t516 - t515;
t461 = -rSges(4,2) * t362 + t609;
t436 = (-t616 * qJD(3) - t232) * qJD(3);
t437 = t389 * t485 + t502 + (t130 + t324) * qJD(1);
t438 = t478 * t516;
t88 = rSges(6,1) * t169 + rSges(6,2) * t677;
t15 = -t182 * t307 - t321 * t88 + t436 * t389 + (t438 + t602) * qJD(1) + t437;
t16 = t182 * t306 - t322 * t88 + t436 * t390 + ((pkin(4) * qJD(3) - qJD(4)) * t571 + t503 - t601) * qJD(1) + t498;
t459 = t15 * t389 + t16 * t390;
t412 = -t322 * t182 + t438 + t532;
t38 = (t501 - t562) * qJD(1) + t412;
t39 = t487 - t182 * t321 - t318 + (t500 + t561) * qJD(1) + t541;
t458 = -t38 * t390 - t389 * t39;
t457 = -t389 * t84 - t390 * t83;
t442 = -pkin(4) * t518 - t232 - t88;
t264 = t302 * t389;
t263 = t301 * t389;
t113 = (t234 * t389 + t236 * t390) * qJD(3);
t423 = -t303 - t616;
t422 = qJD(3) * t476 + t532;
t418 = -t354 + t423;
t413 = -t354 + t535;
t380 = t390 * rSges(5,2);
t325 = t363 * t513;
t285 = t461 * qJD(3);
t274 = t300 * t521;
t267 = t301 * t390;
t249 = t506 - t528;
t248 = t489 - t493;
t246 = t526 * t519;
t233 = t304 * t389 - t380;
t171 = qJD(1) * t649 - t366;
t170 = t365 + (-t249 - t326) * qJD(1);
t161 = -qJD(3) * t264 + (t390 * t461 + t375) * qJD(1);
t160 = -qJD(3) * t263 + (t304 * t390 + t377) * qJD(1);
t159 = rSges(4,1) * t421 - rSges(4,2) * t489 + t533;
t158 = rSges(5,1) * t421 - rSges(5,3) * t493 + t534;
t133 = t356 + (-qJD(1) * t440 - t283) * qJD(1);
t132 = qJD(1) * (-qJD(1) * t506 + t530) + t542;
t72 = (t233 * t389 + t235 * t390) * qJD(3) + t463;
t71 = -t656 + (t235 + t500) * qJD(1) + t541;
t70 = (-t233 + t501) * qJD(1) + t422;
t65 = -t285 * t516 + t356 + (-t161 + t492 + t556) * qJD(1);
t64 = -t285 * t517 + (t159 - t491) * qJD(1) + t502;
t40 = t127 * t321 + t129 * t322 + (t298 * t389 + t299 * t390) * qJD(3) + t463;
t31 = t547 * t516 + (-t160 + t503 + t656) * qJD(1) + t498;
t30 = qJD(1) * t158 + (qJD(1) * t476 + t389 * t547) * qJD(3) + t437;
t21 = (t158 * t390 + t160 * t389 + (t233 * t390 + t389 * t546) * qJD(1)) * qJD(3) + t471;
t10 = t127 * t307 - t129 * t306 + t321 * t63 + t322 * t62 + (t199 * t390 + t200 * t389 + (t298 * t390 + t389 * t543) * qJD(1)) * qJD(3) + t471;
t1 = [t12 * t622 + (t66 + t54) * t306 / 0.2e1 + (t67 + t55) * t307 / 0.2e1 + (t13 + t18) * t321 / 0.2e1 + (((t78 - t187 + (t214 + t577) * t390 + t716) * t390 + (t646 + t684 - t726) * t389) * qJD(3) + t676) * t479 + (t710 * qJD(3) + t169 * t179 + t677 * t176 + t282 * t87 + t714 * t362 - t715 * t363 - t435 * t86) * qJD(1) + (-(-t562 * qJD(1) - t38 + t412 + t638) * t39 + t16 * (-t460 + t531) + t38 * (t318 - t465 + t496) + t15 * (t341 + t495 + t549) + t39 * (t497 + t611) + (-t16 * t618 + t39 * (-pkin(3) - pkin(4)) * t519) * t390 + (t16 * t423 + t38 * (-qJ(4) * t518 - t351) + t15 * t507) * t389 + ((t38 * t618 + t39 * t418) * t389 + (t38 * t418 + t39 * t507) * t390) * qJD(1)) * m(6) + (-(-qJD(1) * t233 + t422 + t638 - t70) * t71 + t31 * (t380 + t531) + t70 * t496 + t30 * (t235 + t495) + t71 * (t497 + t534) + (t71 * t619 * t519 + (-t71 * t387 + t413 * t70) * qJD(1)) * t390 + (-t30 * t387 + t31 * t619 * t363 + (-t70 * qJD(4) + t31 * t603) * t362 + t70 * (t363 * t603 + t605) * qJD(3) + (-t70 * rSges(5,2) + t413 * t71) * qJD(1)) * t389) * m(5) + (t65 * (-t234 + t531) + t83 * t529 + t64 * t652 + t84 * (t365 + t533) + (t302 * t604 - t607) * qJD(3) + ((-t83 * rSges(4,3) + t84 * (-t354 - t609)) * t389 + (t83 * (-t354 - t461) - t84 * t387) * t390) * qJD(1) - (-qJD(1) * t234 + t464 + t653 - t83) * t84) * m(4) + (t133 * (t389 * t486 + t368 + t528) + t170 * t366 + t132 * t649 + t171 * (t527 + t530) + (t170 * (t486 + t608) * t390 + (t170 * (-rSges(3,3) - qJ(2)) + t171 * t486) * t389) * qJD(1) - (-qJD(1) * t249 - t170 - t309 + t365) * t171) * m(3) - (t14 + t12 + t17) * t322 / 0.2e1 + (t658 + t661) * t517 / 0.2e1 + (((t390 * t474 - t646 + t682) * t390 + (t389 * t474 - t197 + t466 + t475 + t683) * t389) * qJD(3) + t663 - t671) * t482 - (-t659 + t660 + t662) * t516 / 0.2e1 + ((t655 + t706) * t389 + (t654 + t705) * t390) * qJD(3) * t613; 0.2e1 * (t15 * t620 + t16 * t621) * m(6) + 0.2e1 * (t30 * t620 + t31 * t621) * m(5) + 0.2e1 * (t620 * t64 + t621 * t65) * m(4) + 0.2e1 * (t132 * t620 + t133 * t621) * m(3); (t659 * t390 + t658 * t389 + (t389 * t655 + t390 * t654) * qJD(1)) * t613 + ((-t517 * t650 - t645) * t389 + ((t389 * t651 + t648) * qJD(3) + t647) * t390) * t482 + ((-t516 * t651 + t645) * t390 + ((t390 * t650 + t648) * qJD(3) + t647) * t389) * t479 + (-t38 * (t325 + t678) - t39 * (t544 + t679) - t40 * (t499 - t665) - (t38 * (-t162 + t262) + t39 * (-pkin(4) * t570 + t164)) * qJD(1) - (t458 * t303 + (-t362 * t40 * t526 + t363 * t458) * pkin(4)) * qJD(3) + t38 * t274 + t10 * t545 + t40 * t504 + (t15 * t470 + t39 * t442 + t10 * t562 + t40 * t601 + (t38 * t182 + t40 * (-t129 + t543)) * qJD(1)) * t389 + (t16 * t470 + t38 * t442 + t10 * t561 + t40 * t602 + (t39 * t470 + t40 * t562) * qJD(1)) * t390) * m(6) + (t70 * t274 + t21 * t545 + t72 * t504 + (t31 * t536 + t70 * t547 + t21 * t235 + t72 * t158 + (t72 * t233 + t536 * t71) * qJD(1)) * t390 + (t30 * t536 + t71 * t547 + t21 * t233 + t72 * t160 + (t70 * t301 + t546 * t72) * qJD(1)) * t389 - t70 * (t325 + (t262 + t263) * qJD(1)) - t71 * (-qJD(1) * t267 + t544) - t72 * t499 - ((-t72 * t267 + t535 * t70) * t390 + (-t72 * t263 + t535 * t71) * t389) * qJD(3)) * m(5) + (-(t264 * t83 - t607) * qJD(1) - (t113 * (-t264 * t389 - t268 * t390) + t457 * t461) * qJD(3) + 0.2e1 * t113 * (t159 * t390 + t161 * t389 + (t234 * t390 - t236 * t389) * qJD(1)) + t457 * t285 + (-t64 * t389 - t65 * t390 + (-t390 * t84 + t604) * qJD(1)) * t302) * m(4) + (t661 * qJD(1) + ((t682 * qJD(1) + t674 * t390) * t390 + (t669 * t389 + t683 * qJD(1) + (-t670 + t675) * t390) * t389) * t664) * t621 + (t660 * qJD(1) + ((t708 * qJD(1) + t670 * t390) * t390 + (t675 * t389 + t684 * qJD(1) + (-t669 + t674) * t390) * t389) * t664) * t620 + (t707 - t631 * t435 + (t694 * t362 + t363 * t696) * qJD(3) + (t362 * t690 + t363 * t691 - t559 * t435) * qJD(1)) * t614 + (t663 + t672) * t521 / 0.2e1 + (t662 + t673) * t520 / 0.2e1 - t709; -m(5) * (t246 * t72 + t247 * t71 + t248 * t70) - m(6) * (t246 * t40 + t247 * t39 + t248 * t38) + 0.2e1 * ((t516 * t70 + t517 * t71 - t21) * t629 + (t38 * t516 + t39 * t517 - t10) * t628) * t363 + 0.2e1 * ((qJD(3) * t72 + t30 * t389 + t31 * t390 + t520 * t71 - t521 * t70) * t629 + (qJD(3) * t40 - t38 * t521 + t39 * t520 + t459) * t628) * t362; (-t435 * t666 - t707) * t614 + ((t38 * t88 - t10 * t129 + t40 * (-qJD(1) * t127 - t62)) * t390 + (t39 * t88 - t10 * t127 + t40 * (qJD(1) * t129 - t63)) * t389 + ((-t38 * t389 + t39 * t390) * qJD(1) + t459) * t182 - t38 * (qJD(1) * t162 - t678) - t39 * (-qJD(1) * t164 - t679) - t40 * t665) * m(6) + t709;];
tauc = t1(:);
