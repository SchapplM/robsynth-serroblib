% Calculate vector of centrifugal and Coriolis load on the joints for
% S5RPRPR16
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% qJD [5x1]
%   Generalized joint velocities
% pkin [7x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d1,d3,d5]';
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
% Datum: 2019-12-31 18:39
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauc = S5RPRPR16_coriolisvecJ_fixb_slag_vp1(qJ, qJD, ...
  pkin, m, rSges, Icges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(7,1),zeros(6,1),zeros(6,3),zeros(6,6)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RPRPR16_coriolisvecJ_fixb_slag_vp1: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RPRPR16_coriolisvecJ_fixb_slag_vp1: qJD has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [7 1]), ...
  'S5RPRPR16_coriolisvecJ_fixb_slag_vp1: pkin has to be [7x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5RPRPR16_coriolisvecJ_fixb_slag_vp1: m has to be [6x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [6,3]), ...
  'S5RPRPR16_coriolisvecJ_fixb_slag_vp1: rSges has to be [6x3] (double)');
assert(isreal(Icges) && all(size(Icges) == [6 6]), ...
  'S5RPRPR16_coriolisvecJ_fixb_slag_vp1: Icges has to be [6x6] (double)'); 

%% Symbolic Calculation
% From coriolisvec_joint_fixb_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 18:38:45
% EndTime: 2019-12-31 18:39:24
% DurationCPUTime: 33.15s
% Computational Cost: add. (8888->887), mult. (23435->1155), div. (0->0), fcn. (21032->6), ass. (0->442)
t362 = sin(qJ(3));
t365 = cos(qJ(3));
t575 = Icges(5,6) * t365;
t281 = Icges(5,2) * t362 + t575;
t584 = Icges(4,4) * t365;
t445 = Icges(4,1) * t362 + t584;
t721 = -t281 - t445;
t349 = Icges(5,6) * t362;
t282 = -Icges(5,2) * t365 + t349;
t585 = Icges(4,4) * t362;
t290 = Icges(4,1) * t365 - t585;
t720 = t282 - t290;
t288 = -Icges(4,2) * t362 + t584;
t437 = -Icges(5,3) * t362 + t575;
t705 = -t288 - t437;
t363 = sin(qJ(1));
t366 = cos(qJ(1));
t443 = Icges(4,2) * t365 + t585;
t193 = -Icges(4,6) * t363 + t366 * t443;
t197 = -Icges(4,5) * t363 + t366 * t445;
t436 = Icges(5,3) * t365 + t349;
t198 = Icges(5,5) * t363 + t366 * t436;
t200 = Icges(5,4) * t363 + t281 * t366;
t428 = t198 * t365 + t200 * t362;
t691 = t193 * t365 + t197 * t362 + t428;
t719 = t436 + t443;
t284 = Icges(4,5) * t365 - Icges(4,6) * t362;
t442 = Icges(5,4) * t365 - Icges(5,5) * t362;
t718 = t442 - t284;
t440 = Icges(4,5) * t362 + Icges(4,6) * t365;
t189 = -Icges(4,3) * t363 + t366 * t440;
t285 = Icges(5,4) * t362 + Icges(5,5) * t365;
t202 = Icges(5,1) * t363 + t285 * t366;
t717 = t189 - t202;
t424 = t288 * t365 + t290 * t362;
t426 = t282 * t362 - t365 * t437;
t704 = t424 - t426;
t192 = Icges(4,6) * t366 + t363 * t443;
t199 = Icges(5,5) * t366 - t363 * t436;
t716 = t192 - t199;
t715 = t193 + t198;
t563 = t363 * t365;
t327 = Icges(4,4) * t563;
t567 = t362 * t363;
t579 = Icges(4,5) * t366;
t196 = Icges(4,1) * t567 + t327 + t579;
t201 = Icges(5,4) * t366 - t281 * t363;
t714 = -t196 + t201;
t713 = t197 + t200;
t712 = t705 * t366;
t711 = t720 * t363;
t710 = t720 * t366;
t709 = t719 * qJD(3);
t708 = t721 * qJD(3);
t706 = t285 - t440;
t659 = t718 * t366;
t703 = t691 * t366;
t188 = Icges(4,3) * t366 + t363 * t440;
t72 = t366 * t188 + t192 * t563 + t196 * t567;
t203 = Icges(5,1) * t366 - t285 * t363;
t181 = t366 * t203;
t427 = t199 * t365 + t201 * t362;
t79 = -t363 * t427 + t181;
t702 = t72 + t79;
t180 = t366 * t202;
t671 = -t366 * t189 - t193 * t563 - t197 * t567 - t363 * t428 + t180;
t176 = t363 * t188;
t431 = t192 * t365 + t196 * t362;
t74 = -t431 * t366 + t176;
t561 = t365 * t366;
t566 = t362 * t366;
t77 = t199 * t561 + t201 * t566 + t363 * t203;
t701 = t74 + t77;
t670 = -t717 * t363 + t703;
t700 = t704 * t363 - t659;
t669 = t362 * t715 - t365 * t713;
t666 = t362 * t716 + t714 * t365;
t699 = qJD(1) * t716 + t712 * qJD(3);
t232 = t437 * t363;
t516 = qJD(3) * t363;
t698 = -qJD(1) * t715 - qJD(3) * t232 - t288 * t516;
t697 = t710 * qJD(3) + (t363 * t445 - t201 + t579) * qJD(1);
t696 = qJD(1) * t713 - qJD(3) * t711;
t237 = t363 * t284;
t240 = t442 * t363;
t664 = t237 - t240;
t695 = t704 * qJD(1) + qJD(3) * t706;
t694 = (Icges(4,2) * t567 - t232 - t327 + t714) * t366 + (-t712 + t713) * t363;
t693 = (-t711 - t716) * t366 + (t710 + t715) * t363;
t692 = t427 - t431;
t347 = t366 * qJ(2);
t294 = pkin(1) * t363 - t347;
t278 = qJD(1) * t294;
t344 = qJD(2) * t363;
t518 = qJD(1) * t366;
t530 = qJ(2) * t518 + t344;
t690 = t530 - t344 + t278;
t689 = -t705 - t721;
t688 = t720 + t719;
t687 = -t709 * t365 + t708 * t362 + (t362 * t705 - t365 * t720) * qJD(3) + t718 * qJD(1);
t686 = t700 * qJD(1);
t361 = sin(qJ(5));
t364 = cos(qJ(5));
t439 = Icges(6,5) * t361 + Icges(6,6) * t364;
t186 = Icges(6,3) * t365 + t362 * t439;
t187 = -Icges(6,3) * t362 + t365 * t439;
t511 = qJD(5) * t362;
t514 = qJD(3) * t366;
t265 = t363 * t511 + t514;
t510 = qJD(5) * t365;
t332 = qJD(1) + t510;
t509 = qJD(5) * t366;
t410 = t362 * t509 - t516;
t581 = Icges(6,4) * t361;
t441 = Icges(6,2) * t364 + t581;
t190 = Icges(6,6) * t365 + t362 * t441;
t580 = Icges(6,4) * t364;
t326 = t362 * t580;
t569 = t361 * t362;
t577 = Icges(6,5) * t365;
t194 = Icges(6,1) * t569 + t326 + t577;
t432 = t190 * t364 + t194 * t361;
t568 = t361 * t366;
t230 = -t364 * t563 - t568;
t562 = t364 * t366;
t231 = -t361 * t563 + t562;
t582 = Icges(6,4) * t231;
t117 = Icges(6,2) * t230 + Icges(6,6) * t567 + t582;
t216 = Icges(6,4) * t230;
t120 = Icges(6,1) * t231 + Icges(6,5) * t567 + t216;
t434 = t117 * t364 + t120 * t361;
t564 = t363 * t364;
t229 = t361 * t561 + t564;
t215 = Icges(6,4) * t229;
t228 = t361 * t363 - t364 * t561;
t115 = -Icges(6,2) * t228 - Icges(6,6) * t566 + t215;
t214 = Icges(6,4) * t228;
t119 = -Icges(6,1) * t229 + Icges(6,5) * t566 + t214;
t435 = t115 * t364 - t119 * t361;
t627 = t410 * (-t186 * t366 + t435) - t265 * (t186 * t363 + t434) - t332 * (t187 + t432);
t685 = t627 * t362;
t684 = (t671 * t363 + t366 * t702) * qJD(3);
t683 = (t363 * t670 + t366 * t701) * qJD(3);
t107 = t366 * t424 - t237;
t108 = t366 * t426 - t240;
t682 = (-t107 + t108) * qJD(1);
t681 = (t188 + t203) * qJD(1);
t680 = t717 * qJD(1);
t679 = (t362 * t689 + t365 * t688) * qJD(1);
t113 = -Icges(6,5) * t229 + Icges(6,6) * t228 + Icges(6,3) * t566;
t678 = -t230 * t115 + t119 * t231;
t39 = -t113 * t567 - t678;
t114 = Icges(6,5) * t231 + Icges(6,6) * t230 + Icges(6,3) * t567;
t40 = t114 * t567 + t230 * t117 + t231 * t120;
t67 = t186 * t567 + t190 * t230 + t194 * t231;
t11 = t265 * t40 + t67 * t332 - t39 * t410;
t43 = -t365 * t113 + t362 * t435;
t37 = t113 * t566 - t115 * t228 - t119 * t229;
t677 = t684 + t686;
t676 = t682 + t683;
t675 = qJD(3) * t692 + t362 * t698 + t365 * t696;
t674 = qJD(3) * t691 - t362 * t699 + t365 * t697;
t673 = t363 * t695 - t366 * t687;
t672 = t363 * t687 + t366 * t695;
t354 = t366 * rSges(4,3);
t207 = rSges(4,1) * t567 + rSges(4,2) * t563 + t354;
t300 = t366 * pkin(1) + t363 * qJ(2);
t356 = t366 * pkin(6);
t639 = t356 + t300;
t665 = t207 + t639;
t663 = -qJD(1) * t691 + t659 * qJD(3) + t681;
t662 = -qJD(1) * t692 + t664 * qJD(3) + t680;
t486 = t362 * t514;
t519 = qJD(1) * t365;
t493 = t363 * t519;
t661 = t486 + t493;
t515 = qJD(3) * t365;
t490 = t363 * t515;
t494 = t362 * t518;
t400 = t490 + t494;
t660 = t706 * qJD(1);
t657 = qJD(3) * t666 - t362 * t696 + t365 * t698 + t681;
t656 = qJD(3) * t669 + t697 * t362 + t699 * t365 + t680;
t123 = rSges(6,1) * t229 - rSges(6,2) * t228 - rSges(6,3) * t566;
t314 = pkin(7) * t490;
t655 = -t123 * t332 + t314;
t379 = t186 * t566 + t190 * t228 - t194 * t229;
t654 = t379 * t332 + t410 * t37;
t653 = -t362 * t693 + t365 * t694;
t652 = 0.2e1 * qJD(3);
t651 = pkin(7) * t365;
t650 = t228 * t117 - t229 * t120;
t649 = pkin(7) * qJD(3);
t343 = qJD(4) * t365;
t454 = rSges(5,2) * t365 - rSges(5,3) * t362;
t642 = (-qJD(3) * t454 - t343) * t366;
t609 = pkin(6) * t363;
t473 = -t294 - t609;
t469 = -rSges(3,2) * t366 + t363 * rSges(3,3);
t640 = t300 + t469;
t602 = rSges(5,2) * t362;
t292 = rSges(5,3) * t365 + t602;
t355 = t366 * rSges(5,1);
t211 = -t292 * t363 + t355;
t277 = t366 * pkin(4) + pkin(7) * t567;
t571 = qJ(4) * t365;
t610 = pkin(3) * t362;
t291 = t571 - t610;
t637 = qJD(3) * t291;
t210 = rSges(5,1) * t363 + t292 * t366;
t444 = Icges(6,1) * t361 + t580;
t634 = t362 * t444 + t577;
t455 = rSges(6,1) * t361 + rSges(6,2) * t364;
t206 = rSges(6,3) * t365 + t362 * t455;
t626 = t410 * (-Icges(6,2) * t229 - t119 - t214) - t265 * (-Icges(6,2) * t231 + t120 + t216) - t332 * (-Icges(6,2) * t569 + t194 + t326);
t359 = t363 ^ 2;
t625 = m(5) / 0.2e1;
t624 = m(6) / 0.2e1;
t623 = -pkin(1) - pkin(6);
t503 = qJD(3) * qJD(5);
t481 = t365 * t503;
t174 = qJD(1) * t410 + t363 * t481;
t622 = t174 / 0.2e1;
t175 = qJD(1) * t265 - t366 * t481;
t621 = t175 / 0.2e1;
t620 = t410 / 0.2e1;
t619 = -t410 / 0.2e1;
t618 = -t265 / 0.2e1;
t617 = t265 / 0.2e1;
t616 = -t332 / 0.2e1;
t615 = t332 / 0.2e1;
t614 = t363 / 0.2e1;
t613 = t365 / 0.2e1;
t611 = rSges(3,2) - pkin(1);
t468 = -qJD(5) - t519;
t378 = t363 * t468 - t486;
t104 = -t332 * t568 + t364 * t378;
t105 = t332 * t562 + t361 * t378;
t489 = t365 * t514;
t520 = qJD(1) * t363;
t495 = t362 * t520;
t399 = -t489 + t495;
t54 = Icges(6,5) * t105 + Icges(6,6) * t104 + Icges(6,3) * t399;
t56 = Icges(6,4) * t105 + Icges(6,2) * t104 + Icges(6,6) * t399;
t58 = Icges(6,1) * t105 + Icges(6,4) * t104 + Icges(6,5) * t399;
t8 = (qJD(3) * t435 + t54) * t365 + (qJD(3) * t113 + t361 * t58 + t364 * t56 + (-t115 * t361 - t119 * t364) * qJD(5)) * t362;
t608 = t8 * t410;
t421 = t468 * t366;
t517 = qJD(3) * t362;
t102 = t364 * t421 + (t332 * t361 + t364 * t517) * t363;
t485 = t362 * t516;
t103 = -t332 * t564 + (t421 + t485) * t361;
t53 = Icges(6,5) * t103 + Icges(6,6) * t102 + Icges(6,3) * t400;
t55 = Icges(6,4) * t103 + Icges(6,2) * t102 + Icges(6,6) * t400;
t57 = Icges(6,1) * t103 + Icges(6,4) * t102 + Icges(6,5) * t400;
t9 = (qJD(3) * t434 + t53) * t365 + (-qJD(3) * t114 + t361 * t57 + t364 * t55 + (-t117 * t361 + t120 * t364) * qJD(5)) * t362;
t607 = t9 * t265;
t601 = rSges(3,3) * t366;
t125 = t231 * rSges(6,1) + t230 * rSges(6,2) + rSges(6,3) * t567;
t208 = -rSges(6,3) * t362 + t365 * t455;
t247 = (rSges(6,1) * t364 - rSges(6,2) * t361) * t362;
t143 = qJD(3) * t208 + qJD(5) * t247;
t532 = pkin(7) * t494 + t314;
t182 = -pkin(4) * t520 + t532;
t342 = qJD(4) * t362;
t217 = t342 + t637;
t335 = pkin(7) * t566;
t367 = qJD(3) ^ 2;
t499 = pkin(3) * t400 + qJ(4) * t485;
t513 = qJD(4) * t363;
t122 = (-qJ(4) * t518 - t513) * t365 + t499;
t368 = qJD(1) ^ 2;
t505 = qJD(1) * qJD(2);
t541 = qJD(1) * (-pkin(1) * t520 + t530) + t363 * t505;
t422 = -t368 * t609 + t541;
t297 = pkin(3) * t365 + qJ(4) * t362;
t504 = qJD(1) * qJD(3);
t483 = t297 * t504;
t402 = qJD(1) * t122 + t363 * t483 + t422;
t488 = t363 * t343;
t512 = qJD(4) * t366;
t59 = t103 * rSges(6,1) + t102 * rSges(6,2) + rSges(6,3) * t400;
t12 = t367 * t335 - t143 * t265 - t174 * t206 + t332 * t59 + (t182 - t488) * qJD(1) + (pkin(7) * t493 - t217 * t366 + (-qJD(5) * t125 - t512) * t362) * qJD(3) + t402;
t599 = t12 * t366;
t317 = pkin(7) * t489;
t183 = qJD(1) * t277 - t317;
t338 = t366 * t505;
t466 = -t356 * t368 + t338;
t401 = qJD(4) * t485 + t217 * t516 + t366 * t483 + t466;
t487 = t365 * t512;
t498 = pkin(3) * t489 + qJ(4) * t661;
t404 = t487 - t498;
t121 = pkin(3) * t495 + t404;
t345 = qJD(2) * t366;
t218 = qJD(1) * t300 - t345;
t557 = -t121 - t218;
t456 = rSges(6,1) * t105 + rSges(6,2) * t104;
t60 = rSges(6,3) * t399 + t456;
t13 = -t143 * t410 + t175 * t206 - t332 * t60 + (-pkin(7) * t363 * t367 + t123 * t503) * t362 + (-t183 + (-qJD(4) + t649) * t561 + t557) * qJD(1) + t401;
t598 = t13 * t363;
t308 = rSges(5,3) * t485;
t146 = -rSges(5,2) * t490 - qJD(1) * t210 + t308;
t273 = t292 * qJD(3);
t540 = -t217 - t273;
t32 = (t146 - t488) * qJD(1) + (-t454 * t520 + (-t342 + t540) * t366) * qJD(3) + t402;
t595 = t32 * t366;
t254 = t454 * t366;
t147 = qJD(1) * t211 + qJD(3) * t254;
t33 = t273 * t516 + (-t147 + t557 + t642) * qJD(1) + t401;
t594 = t33 * t363;
t334 = pkin(3) * t567;
t248 = -qJ(4) * t563 + t334;
t465 = t248 + t639;
t537 = -t297 * t514 - t345;
t42 = t487 + t125 * t332 - t206 * t265 - t317 + (t277 + t465) * qJD(1) + t537;
t593 = t366 * t42;
t492 = t365 * t518;
t497 = rSges(4,1) * t400 + rSges(4,2) * t492;
t145 = (-rSges(4,2) * t517 - rSges(4,3) * qJD(1)) * t363 + t497;
t299 = rSges(4,1) * t365 - rSges(4,2) * t362;
t259 = t299 * t516;
t457 = rSges(4,1) * t362 + rSges(4,2) * t365;
t274 = t457 * qJD(3);
t61 = t274 * t514 + (t145 + t259) * qJD(1) + t422;
t592 = t366 * t61;
t352 = t363 * rSges(4,3);
t209 = t457 * t366 - t352;
t88 = t259 + t344 + (t209 + t473) * qJD(1);
t591 = t366 * t88;
t590 = t42 * t206;
t589 = t43 * t175;
t44 = t114 * t365 + t362 * t434;
t588 = t44 * t174;
t255 = t299 * t366;
t144 = -qJD(3) * t255 + (t363 * t457 + t354) * qJD(1);
t491 = t299 * t514;
t62 = -t274 * t516 + (-t144 - t218 + t491) * qJD(1) + t466;
t587 = t62 * t363;
t565 = t363 * t297;
t560 = qJD(3) * t343 + t121 * t514;
t556 = -t122 - t146;
t555 = -t122 - t182;
t554 = t363 * t217 + t297 * t518;
t253 = t297 * t366;
t544 = -t253 * t514 + t343;
t325 = qJ(4) * t561;
t252 = pkin(3) * t566 - t325;
t543 = -t210 + t252;
t542 = -t211 - t248;
t539 = -t248 - t277;
t276 = pkin(4) * t363 - t335;
t538 = t252 - t276;
t529 = rSges(3,2) * t520 + rSges(3,3) * t518;
t528 = -t252 * t514 + t342;
t526 = t366 ^ 2 + t359;
t506 = -pkin(4) + t623;
t502 = -rSges(5,1) + t623;
t501 = -rSges(4,3) + t623;
t500 = -t125 + t539;
t484 = t567 / 0.2e1;
t482 = t362 * t503;
t479 = t518 / 0.2e1;
t478 = -t517 / 0.2e1;
t477 = -t516 / 0.2e1;
t476 = t516 / 0.2e1;
t475 = -t514 / 0.2e1;
t474 = t514 / 0.2e1;
t471 = t206 + t651;
t467 = t334 + t639;
t464 = t252 + t473;
t461 = qJD(5) * t478;
t460 = qJD(1) * t565 - t362 * t512;
t459 = qJD(1) * t253 + t291 * t516 + t362 * t513;
t38 = -t114 * t566 - t650;
t453 = t363 * t38 - t366 * t37;
t452 = t363 * t37 + t366 * t38;
t451 = t363 * t40 - t366 * t39;
t450 = t363 * t39 + t366 * t40;
t449 = t363 * t44 - t366 * t43;
t448 = t363 * t43 + t366 * t44;
t89 = qJD(1) * t665 - t345 - t491;
t447 = t363 * t88 - t366 * t89;
t446 = t499 + t530;
t433 = t123 * t363 + t125 * t366;
t236 = (Icges(6,5) * t364 - Icges(6,6) * t361) * t362;
t126 = qJD(3) * t187 + qJD(5) * t236;
t191 = -Icges(6,6) * t362 + t365 * t441;
t129 = (-Icges(6,2) * t361 + t580) * t511 + t191 * qJD(3);
t195 = -Icges(6,5) * t362 + t365 * t444;
t244 = (Icges(6,1) * t364 - t581) * t362;
t132 = qJD(3) * t195 + qJD(5) * t244;
t19 = (qJD(3) * t432 + t126) * t365 + (-qJD(3) * t186 + t129 * t364 + t132 * t361 + (-t190 * t361 + t194 * t364) * qJD(5)) * t362;
t71 = t186 * t365 + t362 * t432;
t419 = t19 * t332 - t482 * t71;
t418 = t297 * t516 + t344 - t488;
t412 = -t202 - t427;
t90 = (-t207 * t363 - t209 * t366) * qJD(3);
t403 = -t292 + t610;
t398 = -t113 * t410 - t114 * t265 - t186 * t332;
t397 = -(-Icges(6,5) * t228 - Icges(6,6) * t229) * t410 + (Icges(6,5) * t230 - Icges(6,6) * t231) * t265 + t236 * t332;
t396 = qJD(1) * t252 - t278 + t418;
t377 = (Icges(6,1) * t230 - t117 - t582) * t265 - (-Icges(6,1) * t228 - t115 - t215) * t410 + (-t190 + t244) * t332;
t34 = t123 * t265 + t125 * t410 + (t276 * t366 + t363 * t539) * qJD(3) + t528;
t41 = -t206 * t410 + (-t276 + t464) * qJD(1) + t418 + t655;
t370 = t34 * t433 + (-t363 * t42 - t366 * t41) * t206;
t295 = rSges(3,2) * t363 + t601;
t261 = t297 * t520;
t258 = t454 * t516;
t251 = t299 * t363;
t250 = t454 * t363;
t226 = t526 * t515;
t225 = t485 - t492;
t213 = t366 * t252;
t167 = t206 * t366;
t166 = t206 * t363;
t165 = t634 * t366;
t164 = t634 * t363;
t163 = t190 * t366;
t162 = t190 * t363;
t159 = qJD(1) * t640 - t345;
t158 = t344 + (-t294 + t295) * qJD(1);
t155 = rSges(6,1) * t230 - rSges(6,2) * t231;
t154 = -rSges(6,1) * t228 - rSges(6,2) * t229;
t142 = t338 + (-qJD(1) * t469 - t218) * qJD(1);
t141 = qJD(1) * t529 + t541;
t110 = t366 * t121;
t70 = (t210 * t366 + t363 * t542) * qJD(3) + t528;
t69 = -t642 + (t211 + t465) * qJD(1) + t537;
t68 = -t258 + (-t210 + t464) * qJD(1) + t418;
t17 = (t147 * t366 + t556 * t363 + (t363 * t543 + t366 * t542) * qJD(1)) * qJD(3) + t560;
t16 = t104 * t190 + t105 * t194 - t126 * t566 - t129 * t228 + t132 * t229 + t186 * t399;
t15 = t102 * t190 + t103 * t194 + t126 * t567 + t129 * t230 + t132 * t231 + t186 * t400;
t14 = t265 * t44 + t332 * t71 - t410 * t43;
t10 = t265 * t38 - t654;
t7 = t123 * t174 - t125 * t175 + t410 * t59 + t265 * t60 + (t183 * t366 + t555 * t363 + (t363 * t538 + t366 * t539) * qJD(1)) * qJD(3) + t560;
t6 = t104 * t117 + t105 * t120 + t114 * t399 - t228 * t55 + t229 * t57 - t53 * t566;
t5 = t104 * t115 - t105 * t119 - t113 * t399 - t228 * t56 + t229 * t58 - t54 * t566;
t4 = t102 * t117 + t103 * t120 + t114 * t400 + t230 * t55 + t231 * t57 + t53 * t567;
t3 = t102 * t115 - t103 * t119 - t113 * t400 + t230 * t56 + t231 * t58 + t54 * t567;
t2 = t16 * t332 + t174 * t38 + t175 * t37 + t265 * t6 + t379 * t482 - t410 * t5;
t1 = t15 * t332 + t174 * t40 + t175 * t39 + t265 * t4 - t3 * t410 - t482 * t67;
t18 = [t15 * t617 - t379 * t621 + t67 * t622 - t608 / 0.2e1 + t607 / 0.2e1 + t588 / 0.2e1 + t589 / 0.2e1 + t11 * t620 + t419 + (t16 + t11) * t619 + ((t39 + (t113 * t363 + t114 * t366) * t362 + t650 + t678) * t265 + t10 + t654) * t618 + (((t181 + t72 + t670 - t703) * t366 + (-t74 + t176 + (t412 + t189 - t431) * t366 + t671) * t363) * qJD(3) + t686) * t477 + (-(-t41 + t396 + t655) * t42 + t590 * t410 + t13 * (-t123 + t252 + t335 + t347) + t41 * (rSges(6,3) * t489 + t317 + t345 - t404 - t456) + t12 * (t467 + t125 + t277) + t42 * (t446 + t59 + t532) + (t13 * t506 + (-qJ(4) * t12 - qJD(4) * t42) * t365) * t363) * m(6) + (-(-t258 - t68 + t396) * t69 + t33 * (-t325 + t347) + t68 * (t345 + t498) + t32 * (t355 + t467) + t69 * (t308 + t446) + (t33 * t403 + t68 * (-rSges(5,2) * t515 + rSges(5,3) * t517 - t343)) * t366 + (t33 * t502 - t32 * t602 + (t32 * (-rSges(5,3) - qJ(4)) + t69 * (-rSges(5,2) * qJD(3) - qJD(4))) * t365) * t363) * m(5) + (t62 * (-t352 + t473) + t88 * t345 + t61 * t665 + (qJD(3) * t299 * t88 + t457 * t62) * t366 + (-rSges(4,2) * t485 - t259 + t497 + t690 + t88) * t89) * m(4) + (t142 * (t363 * t611 + t347 + t601) + t158 * t345 + t141 * t640 + (t529 + t158 + t690) * t159) * m(3) + (((t363 * t412 + t181 - t79) * t363 + t189 * t359 + (-t77 - t180 - t176 + (t189 + t431) * t366 + t671) * t366) * qJD(3) + t676 - t682) * t475 + (t672 + t675) * t474 + (t673 + t674 + t677) * t476 - (t366 * t107 + (-t666 + t700) * t363) * t504 / 0.2e1 + (t708 * t365 + t709 * t362 - t704 * qJD(3) + (-(-t276 - t609) * t42 + (t41 * t506 - t42 * t571) * t366 + (t42 * t506 + (-qJ(2) + (-rSges(6,3) - pkin(3) - pkin(7)) * t362) * t41) * t363) * m(6) + (-(-t210 - t609) * t69 + (t68 * t502 + t69 * (-t292 - t571)) * t366 + (t68 * (-qJ(2) - t403) + t69 * t502) * t363) * m(5) + (-(t209 - t609) * t89 + t501 * t591 + (t88 * (-qJ(2) - t457) + t89 * t501) * t363) * m(4) + (t158 * t611 * t366 + (t158 * (-rSges(3,3) - qJ(2)) - t159 * pkin(1)) * t363 - t295 * t159) * m(3) + (t108 + t669) * t474) * qJD(1); 0.2e1 * (-t599 / 0.2e1 + t598 / 0.2e1) * m(6) + 0.2e1 * (-t595 / 0.2e1 + t594 / 0.2e1) * m(5) + 0.2e1 * (t587 / 0.2e1 - t592 / 0.2e1) * m(4) + 0.2e1 * (-t141 * t366 / 0.2e1 + t142 * t614) * m(3); (-qJD(1) * t449 + t363 * t8 + t366 * t9) * t615 + (((t162 * t364 + t164 * t361 - t114) * t265 - (-t163 * t364 - t165 * t361 + t113) * t410 + (t191 * t364 + t195 * t361 - t186) * t332 - t71 * qJD(5)) * t362 + (qJD(5) * t449 - t627) * t365) * t616 + (-qJD(1) * t451 + t3 * t363 + t366 * t4) * t617 + ((t162 * t230 + t164 * t231) * t265 - (-t163 * t230 - t165 * t231) * t410 + (t191 * t230 + t195 * t231) * t332 + (-t362 * t67 - t39 * t561) * qJD(5) + ((qJD(5) * t40 - t398) * t365 - t685) * t363) * t618 + (-qJD(1) * t453 + t363 * t5 + t366 * t6) * t619 + ((-t162 * t228 + t164 * t229) * t265 - (t163 * t228 - t165 * t229) * t410 + (-t191 * t228 + t195 * t229) * t332 + (t362 * t379 + t38 * t563) * qJD(5) + ((-qJD(5) * t37 + t398) * t365 + t685) * t366) * t620 + t452 * t621 + t450 * t622 - t363 * t11 * t510 / 0.2e1 + t14 * t511 / 0.2e1 + t10 * t509 * t613 + t448 * t461 - ((t362 * t694 + t365 * t693) * qJD(3) + (t362 * t688 - t365 * t689) * qJD(1)) * qJD(1) / 0.2e1 + (t675 * t366 + t674 * t363 + (t666 * t363 + t669 * t366) * qJD(1)) * qJD(1) / 0.2e1 + ((t516 * t659 + t660) * t363 + ((t363 * t664 + t653) * qJD(3) + t679) * t366) * t477 + ((t514 * t664 + t660) * t366 + ((t366 * t659 - t653) * qJD(3) - t679) * t363) * t475 + (t13 * t565 - t7 * t213 + (t12 * (-t297 - t471) + t7 * (t123 + t276)) * t366 + (t590 * qJD(1) + t13 * t471 + t7 * t500) * t363 + t593 * t637 - t370 * t510 + (t554 + t206 * t518 + (-pkin(7) * t517 + t143) * t363 - t167 * t332 + t208 * t410 - t459 + t567 * t649 - t123 * t511) * t41 + (t261 + (-t143 - t217) * t366 - t166 * t332 + t208 * t265 - t460 + t125 * t511) * t42 + (t110 + (qJD(1) * t500 + t183 + t60) * t366 + (-t59 + t555 + (-t123 + t538) * qJD(1)) * t363 - t166 * t410 + t167 * t265 - t544 - (-t363 * t565 - t526 * t651) * qJD(3)) * t34) * m(6) + (t33 * t565 + t68 * t554 + t69 * t261 - t17 * t213 + t70 * t110 + (-t33 * t454 + t68 * t273 + t17 * t542 + t70 * t556 + (-t454 * t69 + t543 * t70) * qJD(1)) * t363 + (t32 * (-t297 + t454) + t69 * t540 + t17 * t210 + t70 * t147 + (-t454 * t68 + t542 * t70) * qJD(1)) * t366 - t68 * (-qJD(1) * t254 + t459) - t69 * (-qJD(1) * t250 + t460) - t70 * t544 - ((t69 * (-t291 - t292) + t70 * t254) * t366 + (t68 * t292 + t70 * (-t565 + t250)) * t363) * qJD(3)) * m(5) + (-(t251 * t89 + t255 * t88) * qJD(1) - (t90 * (-t251 * t363 - t255 * t366) - t447 * t457) * qJD(3) + 0.2e1 * t90 * (t144 * t366 - t145 * t363 + (-t207 * t366 + t209 * t363) * qJD(1)) - t447 * t274 + (t587 - t592 + (t363 * t89 + t591) * qJD(1)) * t299) * m(4) + (t2 + t673 * qJD(1) + ((t670 * qJD(1) + t657 * t366) * t366 + (t663 * t363 - t701 * qJD(1) + (-t656 + t662) * t366) * t363) * t652) * t614 + (t1 + t672 * qJD(1) + ((t671 * qJD(1) + t662 * t366) * t366 + (t656 * t363 - t702 * qJD(1) + (-t657 + t663) * t366) * t363) * t652) * t366 / 0.2e1 - (t11 + t677 + t684) * t520 / 0.2e1 + (t10 + t676 + t683) * t479; -m(5) * (t225 * t68 + t226 * t70 - t661 * t69) - m(6) * (t225 * t41 + t226 * t34 - t42 * t661) + 0.2e1 * ((-t514 * t69 + t516 * t68 + t17) * t625 + (t41 * t516 - t42 * t514 + t7) * t624) * t362 + 0.2e1 * ((qJD(3) * t70 - t518 * t68 - t520 * t69 - t594 + t595) * t625 + (qJD(3) * t34 - t41 * t518 - t42 * t520 - t598 + t599) * t624) * t365; t1 * t484 + (t362 * t451 + t365 * t67) * t622 + ((qJD(3) * t451 + t15) * t365 + (qJD(1) * t450 - qJD(3) * t67 - t3 * t366 + t363 * t4) * t362) * t617 - t2 * t566 / 0.2e1 + (t362 * t453 - t365 * t379) * t621 + ((qJD(3) * t453 + t16) * t365 + (qJD(1) * t452 + qJD(3) * t379 + t363 * t6 - t366 * t5) * t362) * t619 + t14 * t478 + (t419 + t588 + t589 + t607 - t608) * t613 + (t362 * t449 + t365 * t71) * t461 + ((qJD(3) * t449 + t19) * t365 + (qJD(1) * t448 - qJD(3) * t71 + t363 * t9 - t366 * t8) * t362) * t615 + (-t230 * t626 + t231 * t377 + t397 * t567) * t618 + (t228 * t626 + t377 * t229 - t397 * t566) * t620 + (t397 * t365 + (t361 * t377 - t364 * t626) * t362) * t616 + (t362 * t479 + t365 * t476) * t11 + (qJD(1) * t484 + t365 * t475) * t10 + ((qJD(3) * t370 + t12 * t125 - t13 * t123 - t41 * t60 + t42 * t59) * t365 + (t41 * (qJD(3) * t123 - t143 * t366) + t42 * (-qJD(3) * t125 - t143 * t363) + t7 * t433 + t34 * (t123 * t518 - t125 * t520 + t363 * t60 + t366 * t59) + (-t12 * t363 - t13 * t366 + (t363 * t41 - t593) * qJD(1)) * t206) * t362 - t41 * (-t154 * t332 - t247 * t410) - t42 * (t155 * t332 - t247 * t265) - t34 * (t154 * t265 + t155 * t410)) * m(6);];
tauc = t18(:);
