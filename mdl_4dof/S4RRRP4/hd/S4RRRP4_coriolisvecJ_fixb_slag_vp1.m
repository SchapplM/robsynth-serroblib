% Calculate vector of centrifugal and Coriolis load on the joints for
% S4RRRP4
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [4x1]
%   Generalized joint coordinates (joint angles)
% qJD [4x1]
%   Generalized joint velocities
% pkin [6x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,d1,d2,d3]';
% m_mdh [5x1]
%   mass of all robot links (including the base)
% rSges [5x3]
%   center of mass of all robot links (in body frames)
%   rows: links of the robot (starting with base)
%   columns: x-, y-, z-coordinates
% Icges [5x6]
%   inertia of all robot links about their respective center of mass, in body frames
%   rows: links of the robot (starting with base)
%   columns: xx, yy, zz, xy, xz, yz (see inertiavector2matrix.m)
% 
% Output:
% tauc [4x1]
%   joint torques required to compensate Coriolis and centrifugal load

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-31 17:15
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauc = S4RRRP4_coriolisvecJ_fixb_slag_vp1(qJ, qJD, ...
  pkin, m, rSges, Icges)
%% Coder Information
%#codegen
%$cgargs {zeros(4,1),zeros(4,1),zeros(6,1),zeros(5,1),zeros(5,3),zeros(5,6)}
assert(isreal(qJ) && all(size(qJ) == [4 1]), ...
  'S4RRRP4_coriolisvecJ_fixb_slag_vp1: qJ has to be [4x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [4 1]), ...
  'S4RRRP4_coriolisvecJ_fixb_slag_vp1: qJD has to be [4x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [6 1]), ...
  'S4RRRP4_coriolisvecJ_fixb_slag_vp1: pkin has to be [6x1] (double)');
assert(isreal(m) && all(size(m) == [5 1]), ...
  'S4RRRP4_coriolisvecJ_fixb_slag_vp1: m has to be [5x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [5,3]), ...
  'S4RRRP4_coriolisvecJ_fixb_slag_vp1: rSges has to be [5x3] (double)');
assert(isreal(Icges) && all(size(Icges) == [5 6]), ...
  'S4RRRP4_coriolisvecJ_fixb_slag_vp1: Icges has to be [5x6] (double)'); 

%% Symbolic Calculation
% From coriolisvec_joint_fixb_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 17:15:10
% EndTime: 2019-12-31 17:15:37
% DurationCPUTime: 23.00s
% Computational Cost: add. (10642->590), mult. (14973->782), div. (0->0), fcn. (11796->6), ass. (0->344)
t705 = Icges(4,4) + Icges(5,4);
t706 = Icges(4,1) + Icges(5,1);
t702 = Icges(4,2) + Icges(5,2);
t699 = Icges(4,5) + Icges(5,5);
t698 = Icges(4,6) + Icges(5,6);
t351 = qJ(2) + qJ(3);
t326 = cos(t351);
t704 = t705 * t326;
t325 = sin(t351);
t703 = t705 * t325;
t701 = t326 * t706 - t703;
t700 = -t325 * t702 + t704;
t697 = Icges(4,3) + Icges(5,3);
t692 = t702 * t326 + t703;
t696 = t706 * t325 + t704;
t353 = sin(qJ(1));
t538 = t325 * t353;
t695 = t705 * t538;
t355 = cos(qJ(1));
t694 = t698 * t355;
t693 = t699 * t355;
t535 = t326 * t353;
t611 = t705 * t535 - t702 * t538 - t694;
t610 = t698 * t353 + t700 * t355;
t677 = t699 * t353 + t701 * t355;
t657 = -t706 * t535 + t693 + t695;
t676 = -t698 * t325 + t699 * t326;
t688 = t699 * t325 + t698 * t326;
t687 = t697 * t355;
t686 = t700 + t696;
t685 = t701 - t692;
t684 = -t692 * t355 + t677;
t683 = -t696 * t355 - t610;
t682 = t696 * t353 + t611;
t678 = -t699 * t535 + t698 * t538 + t687;
t681 = t697 * t353 + t676 * t355;
t680 = t692 * t325 - t326 * t696;
t675 = t611 * t325 + t657 * t326;
t679 = t677 * t535;
t348 = qJD(2) + qJD(3);
t674 = t685 * t348;
t673 = t686 * t348;
t672 = -t677 * qJD(1) + t682 * t348;
t671 = t683 * t348 + (-t353 * t701 + t693) * qJD(1);
t277 = t348 * t353;
t670 = t610 * qJD(1) - t692 * t277 - t657 * t348;
t669 = t684 * t348 + (-t353 * t700 + t694) * qJD(1);
t668 = t688 * t355;
t667 = t688 * t353;
t666 = t610 * t325;
t665 = t353 * t680 + t668;
t664 = -t355 * t680 + t667;
t356 = -pkin(6) - pkin(5);
t347 = -qJ(4) + t356;
t663 = rSges(5,3) - t347;
t662 = t675 * t353;
t661 = t355 * t681 - t679;
t660 = t681 * qJD(1);
t534 = t326 * t355;
t607 = -t353 * t681 - t534 * t677;
t659 = t353 * t678 + t657 * t534;
t658 = t678 * t355;
t633 = t658 - t662;
t632 = -t538 * t610 - t661;
t537 = t325 * t355;
t631 = -t537 * t611 - t659;
t630 = -t537 * t610 - t607;
t656 = t688 * qJD(1) - t673 * t325 + t674 * t326;
t655 = -t325 * t669 + t326 * t671 + t660;
t654 = qJD(1) * t678 + t325 * t670 + t326 * t672;
t278 = t348 * t355;
t653 = t685 * qJD(1) + t683 * t277 + t682 * t278;
t652 = (t702 * t535 + t657 + t695) * t278 + t684 * t277 + t686 * qJD(1);
t651 = qJD(1) * t680 + t676 * t348;
t650 = qJD(1) * t675 - t348 * t667 + t660;
t649 = -t668 * t348 + (-t326 * t677 - t353 * t676 + t666 + t687) * qJD(1);
t648 = t664 * qJD(1);
t354 = cos(qJ(2));
t317 = t354 * pkin(2) + pkin(1);
t271 = pkin(3) * t326 + t317;
t647 = -rSges(5,1) * t535 + rSges(5,2) * t538 - t353 * t271 + t355 * t663;
t337 = t353 * rSges(5,3);
t646 = rSges(5,1) * t534 - rSges(5,2) * t537 + t355 * t271 + t337;
t352 = sin(qJ(2));
t576 = pkin(2) * qJD(2);
t474 = t352 * t576;
t584 = pkin(3) * t325;
t236 = -t348 * t584 - t474;
t323 = qJD(4) * t353;
t480 = qJD(1) * t353;
t467 = t325 * t480;
t479 = qJD(1) * t355;
t645 = rSges(5,2) * t467 + rSges(5,3) * t479 + t355 * t236 + t323;
t644 = t665 * qJD(1);
t643 = -t650 * t353 + t654 * t355;
t642 = t649 * t353 + t655 * t355;
t641 = t654 * t353 + t650 * t355;
t640 = t655 * t353 - t649 * t355;
t639 = t632 * t277 - t633 * t278 - t644;
t638 = t630 * t277 - t631 * t278 + t648;
t637 = t651 * t353 + t656 * t355;
t636 = t656 * t353 - t651 * t355;
t635 = t325 * t672 - t326 * t670;
t634 = t325 * t671 + t326 * t669;
t629 = -t657 * t325 + t611 * t326;
t628 = t325 * t677 + t610 * t326;
t385 = -t278 * t325 - t326 * t480;
t477 = qJD(2) * t355;
t462 = t352 * t477;
t424 = pkin(2) * t462;
t471 = t278 * t326;
t485 = -t347 + t356;
t494 = t271 - t317;
t574 = rSges(5,1) * t385 - rSges(5,2) * t471 + t424 + (-t353 * t494 + t355 * t485) * qJD(1) + t645;
t260 = rSges(5,1) * t325 + rSges(5,2) * t326;
t230 = t260 * t353;
t578 = rSges(5,1) * t326;
t262 = -rSges(5,2) * t325 + t578;
t324 = qJD(4) * t355;
t427 = -t236 * t353 + t324;
t306 = t353 * t474;
t487 = t356 * t480 + t306;
t533 = t347 * t353;
t573 = -t348 * t230 - t427 + t487 + (t337 - t533 + (t262 + t494) * t355) * qJD(1);
t321 = t355 * t356;
t489 = t353 * t317 + t321;
t520 = t489 + t647;
t300 = t355 * t317;
t519 = t353 * t485 - t300 + t646;
t617 = -t652 * t325 + t653 * t326;
t616 = qJD(1) * t676 - t277 * t668 + t278 * t667;
t615 = 0.2e1 * qJD(2);
t614 = rSges(3,2) * t352;
t207 = t262 * t348;
t423 = qJD(1) * t348;
t264 = t353 * t423;
t528 = t354 * qJD(2) ^ 2;
t397 = -pkin(2) * t355 * t528 + qJD(1) * t306;
t343 = t353 * pkin(5);
t581 = pkin(1) - t317;
t142 = (-t355 * t581 - t343) * qJD(1) - t487;
t292 = t355 * pkin(1) + t343;
t272 = t292 * qJD(1);
t521 = -t142 - t272;
t536 = t326 * t348;
t28 = -t207 * t278 + t260 * t264 + (t264 * t325 - t278 * t536) * pkin(3) + (t324 + t521 - t573) * qJD(1) + t397;
t450 = -t260 - t584;
t426 = -t353 * t356 + t300;
t178 = t426 - t292;
t468 = -t178 - t519;
t61 = -t306 - t324 + t450 * t277 + (t292 - t468) * qJD(1);
t613 = qJD(1) * t61 + t28;
t345 = t355 * pkin(5);
t291 = pkin(1) * t353 - t345;
t177 = t291 - t489;
t273 = qJD(1) * t291;
t612 = qJD(1) * t177 - t273;
t609 = t323 - t424;
t336 = Icges(3,4) * t354;
t410 = -Icges(3,2) * t352 + t336;
t284 = Icges(3,1) * t352 + t336;
t608 = t666 + t678;
t261 = rSges(4,1) * t325 + rSges(4,2) * t326;
t231 = t261 * t353;
t233 = t261 * t355;
t579 = rSges(4,1) * t326;
t263 = -rSges(4,2) * t325 + t579;
t192 = rSges(4,1) * t535 - rSges(4,2) * t538 - t355 * rSges(4,3);
t338 = t353 * rSges(4,3);
t194 = rSges(4,1) * t534 - rSges(4,2) * t537 + t338;
t478 = qJD(2) * t353;
t513 = -t177 * t478 + t178 * t477;
t66 = t192 * t277 + t194 * t278 + t513;
t386 = -t261 * t278 - t424;
t508 = t177 - t291;
t67 = (-t192 + t508) * qJD(1) + t386;
t507 = -t178 - t194;
t68 = -t261 * t277 - t306 + (t292 - t507) * qJD(1);
t606 = -t67 * (qJD(1) * t231 - t278 * t263) - t66 * (-t277 * t231 - t233 * t278) - t68 * (-qJD(1) * t233 - t263 * t277);
t530 = t353 * t354;
t532 = t352 * t353;
t557 = Icges(3,3) * t355;
t209 = Icges(3,5) * t530 - Icges(3,6) * t532 - t557;
t311 = Icges(3,4) * t532;
t566 = Icges(3,5) * t355;
t213 = Icges(3,1) * t530 - t311 - t566;
t560 = Icges(3,6) * t355;
t211 = Icges(3,4) * t530 - Icges(3,2) * t532 - t560;
t547 = t211 * t352;
t402 = -t213 * t354 + t547;
t83 = -t355 * t209 - t353 * t402;
t281 = Icges(3,5) * t354 - Icges(3,6) * t352;
t280 = Icges(3,5) * t352 + Icges(3,6) * t354;
t387 = qJD(2) * t280;
t569 = Icges(3,4) * t352;
t285 = Icges(3,1) * t354 - t569;
t214 = Icges(3,5) * t353 + t285 * t355;
t212 = Icges(3,6) * t353 + t355 * t410;
t546 = t212 * t352;
t401 = -t214 * t354 + t546;
t601 = -t355 * t387 + (-t281 * t353 + t401 + t557) * qJD(1);
t210 = Icges(3,3) * t353 + t281 * t355;
t482 = qJD(1) * t210;
t600 = qJD(1) * t402 - t353 * t387 + t482;
t282 = Icges(3,2) * t354 + t569;
t398 = t352 * t282 - t284 * t354;
t597 = t398 * qJD(1) + t281 * qJD(2);
t596 = t353 * (-t282 * t355 + t214) - t355 * (-Icges(3,2) * t530 + t213 - t311);
t593 = t264 / 0.2e1;
t265 = t355 * t423;
t592 = t265 / 0.2e1;
t591 = -t277 / 0.2e1;
t590 = t277 / 0.2e1;
t589 = -t278 / 0.2e1;
t588 = t278 / 0.2e1;
t587 = t353 / 0.2e1;
t586 = -t355 / 0.2e1;
t585 = pkin(2) * t352;
t583 = -qJD(1) / 0.2e1;
t582 = qJD(1) / 0.2e1;
t580 = rSges(3,1) * t354;
t577 = rSges(3,2) * t354;
t339 = t353 * rSges(3,3);
t575 = t67 * t261;
t41 = -t277 * t520 + t278 * t519 + t513;
t554 = qJD(1) * t41;
t486 = rSges(3,2) * t532 + t355 * rSges(3,3);
t216 = rSges(3,1) * t530 - t486;
t287 = rSges(3,1) * t352 + t577;
t463 = t287 * t477;
t119 = -t463 + (-t216 - t291) * qJD(1);
t552 = t119 * t353;
t551 = t119 * t355;
t464 = t287 * t478;
t529 = t354 * t355;
t531 = t352 * t355;
t217 = rSges(3,1) * t529 - rSges(3,2) * t531 + t339;
t499 = t217 + t292;
t120 = qJD(1) * t499 - t464;
t246 = t287 * t355;
t550 = t120 * t246;
t541 = t277 * t326;
t540 = t280 * t353;
t539 = t280 * t355;
t322 = pkin(5) * t479;
t141 = -t424 - t322 + (t353 * t581 - t321) * qJD(1);
t247 = qJD(1) * (-pkin(1) * t480 + t322);
t522 = qJD(1) * t141 + t247;
t512 = -t353 * t177 + t355 * t178;
t511 = t353 * t192 + t355 * t194;
t510 = -t353 * t209 - t213 * t529;
t509 = t353 * t210 + t214 * t529;
t279 = pkin(3) * t467;
t498 = t260 * t480 + t279;
t492 = rSges(4,2) * t467 + rSges(4,3) * t479;
t491 = -t282 + t285;
t490 = t284 + t410;
t488 = rSges(3,3) * t479 + t480 * t614;
t481 = qJD(1) * t281;
t123 = -t353 * t398 - t539;
t476 = t123 * qJD(1);
t475 = qJD(1) * qJD(2);
t473 = t354 * t576;
t114 = rSges(4,1) * t385 - rSges(4,2) * t471 + t492;
t116 = -t348 * t231 + (t263 * t355 + t338) * qJD(1);
t472 = t355 * t114 + t353 * t116 + t192 * t479;
t470 = t353 * t528;
t469 = t355 * t141 + t353 * t142 - t177 * t479;
t465 = t352 * t479;
t460 = -pkin(1) - t580;
t459 = t355 * t475;
t458 = t480 / 0.2e1;
t457 = t479 / 0.2e1;
t456 = -t478 / 0.2e1;
t453 = t477 / 0.2e1;
t452 = -t261 - t585;
t274 = -t584 - t585;
t451 = t274 + t585;
t448 = t352 * (-t353 ^ 2 - t355 ^ 2);
t232 = t260 * t355;
t437 = -t277 * t230 - t232 * t278;
t172 = t214 * t530;
t436 = t355 * t210 - t172;
t429 = -qJD(1) * t232 - t262 * t277;
t428 = -t209 + t546;
t422 = -t520 * t353 + t519 * t355;
t421 = -pkin(3) * t536 - t207;
t208 = t263 * t348;
t418 = -t208 - t473;
t417 = -pkin(3) * t471 + qJD(1) * t230 - t278 * t262;
t415 = t580 - t614;
t414 = -t353 * t68 - t355 * t67;
t407 = -t120 * t353 - t551;
t121 = t211 * t354 + t213 * t352;
t122 = t212 * t354 + t214 * t352;
t396 = t573 * t353 + t574 * t355 - t520 * t479;
t395 = -t260 + t274;
t245 = t287 * t353;
t84 = -t212 * t532 - t436;
t393 = (t353 * t84 - t355 * t83) * qJD(2);
t85 = -t211 * t531 - t510;
t86 = -t212 * t531 + t509;
t392 = (t353 * t86 - t355 * t85) * qJD(2);
t389 = qJD(2) * t284;
t388 = qJD(2) * t282;
t117 = (t216 * t353 + t217 * t355) * qJD(2);
t384 = t421 - t473;
t381 = -t178 * t353 * t475 + t141 * t477 + t142 * t478 - t177 * t459;
t380 = t211 * t355 - t212 * t353;
t379 = (-t352 * t490 + t354 * t491) * qJD(1);
t378 = t278 * t450 + t609;
t135 = qJD(1) * t212 - t353 * t388;
t137 = qJD(1) * t214 - t353 * t389;
t365 = qJD(1) * t209 - qJD(2) * t121 - t135 * t352 + t137 * t354;
t134 = -t355 * t388 + (-t353 * t410 + t560) * qJD(1);
t136 = -t355 * t389 + (-t285 * t353 + t566) * qJD(1);
t364 = -qJD(2) * t122 - t134 * t352 + t136 * t354 + t482;
t267 = t410 * qJD(2);
t268 = t285 * qJD(2);
t363 = qJD(1) * t280 - t267 * t352 + t268 * t354 + (-t282 * t354 - t284 * t352) * qJD(2);
t362 = (t632 * t353 - t633 * t355) * t593 + (t630 * t353 - t631 * t355) * t592 + (t616 * t353 + t617 * t355) * t591 + (t643 * t355 + t642 * t353 + (t631 * t353 + t630 * t355) * qJD(1)) * t590 + (t641 * t355 + t640 * t353 + (t633 * t353 + t632 * t355) * qJD(1)) * t589 + (t617 * t353 - t616 * t355) * t588 + (t637 * qJD(1) + t631 * t264 + t630 * t265 + t642 * t277 + t643 * t278) * t587 + (t636 * qJD(1) + t633 * t264 + t632 * t265 + t640 * t277 + t641 * t278) * t586 + (t653 * t325 + t652 * t326) * t583 + (t635 * t355 + t634 * t353 + (t629 * t353 + t628 * t355) * qJD(1)) * t582 + t639 * t458 + t638 * t457;
t361 = -t352 * t596 + t380 * t354;
t269 = t415 * qJD(2);
t198 = t451 * t355;
t197 = t451 * t353;
t139 = -qJD(2) * t245 + (t355 * t415 + t339) * qJD(1);
t138 = -t477 * t577 + (-t354 * t480 - t462) * rSges(3,1) + t488;
t124 = -t355 * t398 + t540;
t118 = t124 * qJD(1);
t78 = -t269 * t477 + (-t139 - t272 + t464) * qJD(1);
t77 = -t269 * t478 + t247 + (t138 - t463) * qJD(1);
t65 = t363 * t353 - t355 * t597;
t64 = t353 * t597 + t363 * t355;
t63 = -qJD(2) * t401 + t134 * t354 + t136 * t352;
t62 = -t402 * qJD(2) + t135 * t354 + t137 * t352;
t60 = (t508 + t520) * qJD(1) + t378;
t53 = -t208 * t278 + t261 * t264 + (-t116 + t521) * qJD(1) + t397;
t52 = qJD(1) * t114 - t208 * t277 - t261 * t265 + (-t352 * t459 - t470) * pkin(2) + t522;
t43 = t118 + t392;
t42 = t393 + t476;
t27 = -pkin(2) * t470 - t207 * t277 - t260 * t265 + (-t265 * t325 - t277 * t536) * pkin(3) + (t574 + t609) * qJD(1) + t522;
t16 = t114 * t278 + t116 * t277 + t192 * t265 - t194 * t264 + t381;
t9 = -t264 * t519 - t265 * t520 + t277 * t573 + t278 * t574 + t381;
t1 = [(t118 + ((t84 - t172 + (t210 + t547) * t355 + t510) * t355 + t509 * t353) * qJD(2)) * t453 + (((t353 * t610 + t355 * t611) * t325 + t632 + t659 + t661) * t278 + (t657 * t535 + (t353 * t611 - t355 * t610) * t325 - t607 + t633 - t658) * t277 + t648) * t588 + (-t476 + ((t355 * t428 - t509 + t86) * t355 + (t353 * t428 + t436 + t85) * t353) * qJD(2) + t42) * t456 + (t63 + t64) * t478 / 0.2e1 + (-qJD(2) * t398 + t267 * t354 + t268 * t352 + t674 * t325 + t673 * t326) * qJD(1) + (-(qJD(1) * t520 + t378 - t60 + t612) * t61 + t28 * t647 + t60 * t427 + t27 * (-t533 + t646) + t61 * t645 + (t230 * t60 - t232 * t61) * t348 + ((t60 * (-t262 - t271) - t61 * t347) * t355 + (-t60 * t663 + t61 * (-t271 - t578)) * t353) * qJD(1)) * m(5) + (-(-qJD(1) * t192 + t386 + t612 - t67) * t68 + t53 * (-t192 - t489) + t67 * t487 + t52 * (t194 + t426) + t68 * (-t424 + t492) + (-t233 * t68 + t353 * t575) * t348 + ((-t67 * rSges(4,3) + t68 * (-t317 - t579)) * t353 + (t67 * (-t263 - t317) - t68 * t356) * t355) * qJD(1)) * m(4) + (t78 * (t353 * t460 + t345 + t486) + t77 * t499 + t120 * (t322 + t488) + (t287 * t552 - t550) * qJD(2) + ((-pkin(1) - t415) * t551 + (t119 * (-rSges(3,3) - pkin(5)) + t120 * t460) * t353) * qJD(1) - (-qJD(1) * t216 - t119 - t273 - t463) * t120) * m(3) - (t62 + t65 + t43) * t477 / 0.2e1 + (t629 - t665) * t593 + (t628 + t664) * t592 + ((t355 * t608 + t607 + t630 - t662) * t278 + ((t675 + t681) * t355 + t608 * t353 + t631 - t679) * t277 + t639 + t644) * t591 + (t634 + t637) * t590 + (-t635 + t636 + t638) * t589 + ((t121 + t123) * t353 + (t122 + t124) * t355) * t475 / 0.2e1; t362 + ((-t477 * t540 - t481) * t355 + (t379 + (t355 * t539 + t361) * qJD(2)) * t353) * t453 + ((-t478 * t539 + t481) * t353 + (t379 + (t353 * t540 + t361) * qJD(2)) * t355) * t456 + (t353 * t63 - t355 * t62 + (t121 * t353 + t122 * t355) * qJD(1)) * t582 + ((t352 * t491 + t354 * t490) * qJD(1) + (t380 * t352 + t354 * t596) * qJD(2)) * t583 + (qJD(1) * t64 + (-(t353 * t600 + t365 * t355) * t355 + (t353 * t601 + t364 * t355) * t353 + (t85 * t353 + t86 * t355) * qJD(1)) * t615) * t587 + (qJD(1) * t65 + (-(t365 * t353 - t355 * t600) * t355 + t353 * (t364 * t353 - t355 * t601) + (t83 * t353 + t84 * t355) * qJD(1)) * t615) * t586 + (t393 + t42) * t458 + (t43 + t392) * t457 + (t60 * t498 + t9 * (t422 + t512) + t41 * (t396 + t469) + (t384 * t60 + t395 * t613) * t355 + (t27 * t395 + t384 * t61 + t468 * t554) * t353 - t60 * (-qJD(1) * t197 + t417) - t61 * (-pkin(3) * t541 + qJD(1) * t198 + t429) - t41 * (t197 * t277 + t198 * t278 + t437) - (-t61 * t465 + ((-t353 * t61 - t355 * t60) * t354 + t41 * t448) * qJD(2)) * pkin(2)) * m(5) + (t16 * (t511 + t512) + t66 * (t469 + t472) + (t418 * t67 + (qJD(1) * t68 + t53) * t452) * t355 + (t52 * t452 + t68 * t418 + (t507 * t66 + t575) * qJD(1)) * t353 - (-t68 * t465 + (t354 * t414 + t448 * t66) * qJD(2)) * pkin(2) + t606) * m(4) + (-(t119 * t245 - t550) * qJD(1) - (t117 * (-t245 * t353 - t246 * t355) + t407 * t415) * qJD(2) + 0.2e1 * t117 * (t355 * t138 + t353 * t139 + (t216 * t355 - t217 * t353) * qJD(1)) + t407 * t269 + (-t77 * t353 - t78 * t355 + (-t120 * t355 + t552) * qJD(1)) * t287) * m(3); t362 + (t9 * t422 + t41 * t396 + (t421 * t61 - t519 * t554) * t353 - t61 * t429 - t41 * t437 - (-t61 * t541 + (-t61 * t479 + t41 * (-t277 * t353 - t278 * t355)) * t325) * pkin(3) + (t27 * t353 + t355 * t613) * t450 + (t355 * t421 - t279 - t417 + t498) * t60) * m(5) + (t16 * t511 + t66 * (-t194 * t480 + t472) + t414 * t208 + (-t52 * t353 - t53 * t355 + (t353 * t67 - t355 * t68) * qJD(1)) * t261 + t606) * m(4); 0.2e1 * (t27 * t586 + t28 * t587 - t41 * (-t277 * t355 + t278 * t353) / 0.2e1) * m(5);];
tauc = t1(:);
