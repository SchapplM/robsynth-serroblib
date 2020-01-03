% Calculate vector of centrifugal and Coriolis load on the joints for
% S5RPRPP4
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% qJD [5x1]
%   Generalized joint velocities
% pkin [7x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d1,d3,theta4]';
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
% Datum: 2019-12-31 18:15
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauc = S5RPRPP4_coriolisvecJ_fixb_slag_vp1(qJ, qJD, ...
  pkin, m, rSges, Icges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(7,1),zeros(6,1),zeros(6,3),zeros(6,6)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RPRPP4_coriolisvecJ_fixb_slag_vp1: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RPRPP4_coriolisvecJ_fixb_slag_vp1: qJD has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [7 1]), ...
  'S5RPRPP4_coriolisvecJ_fixb_slag_vp1: pkin has to be [7x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5RPRPP4_coriolisvecJ_fixb_slag_vp1: m has to be [6x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [6,3]), ...
  'S5RPRPP4_coriolisvecJ_fixb_slag_vp1: rSges has to be [6x3] (double)');
assert(isreal(Icges) && all(size(Icges) == [6 6]), ...
  'S5RPRPP4_coriolisvecJ_fixb_slag_vp1: Icges has to be [6x6] (double)'); 

%% Symbolic Calculation
% From coriolisvec_joint_fixb_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 18:14:18
% EndTime: 2019-12-31 18:14:49
% DurationCPUTime: 25.68s
% Computational Cost: add. (7427->583), mult. (12070->697), div. (0->0), fcn. (9170->6), ass. (0->325)
t687 = -Icges(6,4) - Icges(5,5);
t686 = Icges(5,6) - Icges(6,6);
t328 = sin(qJ(1));
t330 = cos(qJ(1));
t323 = qJ(3) + pkin(7);
t300 = cos(t323);
t299 = sin(t323);
t536 = Icges(5,4) * t299;
t404 = Icges(5,2) * t300 + t536;
t685 = t686 * t328 - t330 * t404;
t535 = Icges(5,4) * t300;
t407 = Icges(5,1) * t299 + t535;
t684 = t687 * t328 + t330 * t407;
t683 = -Icges(4,3) - Icges(5,3);
t517 = t300 * t330;
t519 = t299 * t330;
t675 = Icges(6,5) * t519 - Icges(6,3) * t517 + t685;
t267 = Icges(6,5) * t517;
t682 = Icges(6,1) * t519 - t267 + t684;
t327 = sin(qJ(3));
t329 = cos(qJ(3));
t680 = Icges(4,5) * t327 + Icges(5,5) * t299 + Icges(4,6) * t329 + Icges(5,6) * t300;
t403 = Icges(6,4) * t299 - Icges(6,6) * t300;
t140 = Icges(6,2) * t330 + t328 * t403;
t674 = -t680 * t328 + t683 * t330;
t681 = t140 - t674;
t531 = Icges(6,5) * t299;
t222 = Icges(6,1) * t300 + t531;
t224 = Icges(5,1) * t300 - t536;
t664 = t222 + t224;
t530 = Icges(6,2) * t328;
t141 = Icges(6,4) * t519 - Icges(6,6) * t517 - t530;
t623 = t683 * t328 + t680 * t330;
t669 = -t141 - t623;
t399 = -Icges(6,3) * t300 + t531;
t679 = t399 - t404;
t292 = Icges(6,5) * t300;
t406 = Icges(6,1) * t299 - t292;
t678 = -t406 - t407;
t677 = -Icges(4,5) * t329 + Icges(4,6) * t327 + t686 * t299 + t687 * t300;
t538 = Icges(4,4) * t327;
t405 = Icges(4,2) * t329 + t538;
t163 = -Icges(4,6) * t328 + t330 * t405;
t537 = Icges(4,4) * t329;
t408 = Icges(4,1) * t327 + t537;
t165 = -Icges(4,5) * t328 + t330 * t408;
t651 = t163 * t329 + t165 * t327 + t682 * t299 - t675 * t300;
t520 = t299 * t328;
t266 = Icges(6,5) * t520;
t518 = t300 * t328;
t528 = Icges(6,6) * t330;
t136 = -Icges(6,3) * t518 + t266 + t528;
t142 = Icges(5,6) * t330 + t328 * t404;
t676 = -t136 + t142;
t144 = Icges(6,4) * t330 + t328 * t406;
t268 = Icges(5,4) * t518;
t532 = Icges(5,5) * t330;
t146 = Icges(5,1) * t520 + t268 + t532;
t673 = -t144 - t146;
t671 = t664 * t330;
t220 = -Icges(5,2) * t299 + t535;
t251 = -Icges(4,2) * t327 + t537;
t253 = Icges(4,1) * t329 - t538;
t670 = t220 * t300 + t224 * t299 + t251 * t329 + t253 * t327;
t162 = Icges(4,6) * t330 + t328 * t405;
t515 = t328 * t329;
t286 = Icges(4,4) * t515;
t516 = t327 * t328;
t533 = Icges(4,5) * t330;
t164 = Icges(4,1) * t516 + t286 + t533;
t648 = -t142 * t300 - t146 * t299 - t162 * t329 - t164 * t327;
t608 = t144 * t299 - t648;
t587 = Icges(6,3) * t299 + t292;
t662 = t222 * t299 - t300 * t587 + t670;
t511 = t330 * t140 + t144 * t520;
t50 = -t136 * t518 + t511;
t641 = -t142 * t518 - t146 * t520 - t162 * t515 - t164 * t516 + t674 * t330;
t668 = t50 - t641;
t667 = t679 * qJD(3);
t666 = t678 * qJD(3);
t665 = -t220 + t587;
t525 = t136 * t300;
t663 = t525 - t608;
t610 = -t136 * t517 - t681 * t328;
t609 = t677 * t328;
t661 = -t403 - t680;
t606 = t677 * t330;
t660 = t651 * t330;
t173 = t220 * t330;
t463 = qJD(3) * t330;
t659 = qJD(3) * t173 - t587 * t463 + (t328 * t399 - t142 + t528) * qJD(1);
t166 = t587 * t328;
t464 = qJD(3) * t328;
t658 = qJD(3) * t166 - t220 * t464 + (t330 * t399 + t685) * qJD(1);
t657 = -t671 * qJD(3) + (t328 * t407 + t144 + t532) * qJD(1);
t176 = t224 * t328;
t656 = qJD(3) * t176 + t222 * t464 + (t330 * t406 + t684) * qJD(1);
t655 = (Icges(6,1) * t518 + t176 + t266 - t676) * t330 + (-t671 - t675) * t328;
t613 = -t163 * t515 - t165 * t516 + t669 * t330 + t675 * t518 - t520 * t682;
t654 = -t144 * t519 + t648 * t330 - t610;
t638 = t669 * t328 + t660;
t653 = t662 * t328 - t606;
t652 = t222 * t519 + t670 * t330 - t517 * t587 + t609;
t612 = t163 * t327 - t165 * t329 - t675 * t299 - t300 * t682;
t611 = t162 * t327 - t164 * t329 + t676 * t299 + t673 * t300;
t626 = rSges(6,1) + pkin(4);
t540 = rSges(6,3) + qJ(5);
t650 = t540 * t299;
t649 = (Icges(5,2) * t520 + t166 - t268 + t673) * t330 + (-Icges(6,3) * t519 + t173 - t267 + t682) * t328;
t647 = t623 * qJD(1);
t646 = -t665 - t678;
t645 = -t664 - t679;
t644 = t662 * qJD(1) + qJD(3) * t661;
t643 = t681 * qJD(1);
t234 = t405 * qJD(3);
t235 = t408 * qJD(3);
t642 = -t234 * t329 - t235 * t327 + t667 * t300 + t666 * t299 + (-t251 * t327 + t253 * t329 + t299 * t665 + t300 * t664) * qJD(3) + t677 * qJD(1);
t289 = qJD(5) * t299;
t547 = rSges(6,3) * t300;
t552 = rSges(6,1) * t299;
t410 = -t547 + t552;
t526 = qJ(5) * t300;
t563 = pkin(4) * t299;
t622 = t526 - t563 - t410;
t501 = qJD(3) * t622 + t289;
t564 = pkin(3) * t329;
t484 = t251 + t408;
t485 = -t405 + t253;
t640 = (t299 * t646 + t300 * t645 + t327 * t484 - t329 * t485) * qJD(1);
t637 = t653 * qJD(1);
t636 = t647 - t609 * qJD(3) + (t330 * t403 - t530 - t663) * qJD(1);
t635 = -t651 * qJD(1) + t606 * qJD(3) + t643;
t634 = (t638 * t328 + t654 * t330) * qJD(3);
t633 = (t613 * t328 + t330 * t668) * qJD(3);
t632 = t652 * qJD(1);
t103 = qJD(1) * t163 + t251 * t464;
t209 = t253 * t328;
t105 = qJD(1) * t165 + qJD(3) * t209;
t631 = t611 * qJD(3) - t103 * t329 - t105 * t327 - t656 * t299 + t658 * t300 + t643;
t208 = t251 * t330;
t102 = qJD(1) * t162 - qJD(3) * t208;
t210 = t253 * t330;
t104 = -qJD(3) * t210 + (t328 * t408 + t533) * qJD(1);
t630 = qJD(1) * t141 + t612 * qJD(3) + t102 * t329 + t104 * t327 + t657 * t299 - t659 * t300 + t647;
t347 = t328 * (t165 + t208) - t330 * (-Icges(4,2) * t516 + t164 + t286);
t348 = t328 * (t163 - t210) - t330 * (t162 - t209);
t629 = -t655 * t299 + t300 * t649 - t348 * t327 + t347 * t329;
t628 = 0.2e1 * qJD(3);
t569 = t328 / 0.2e1;
t627 = -t330 / 0.2e1;
t625 = rSges(4,2) * t329;
t549 = rSges(5,2) * t300;
t411 = rSges(5,1) * t299 + t549;
t624 = t330 * t411;
t229 = pkin(4) * t300 + qJ(5) * t299;
t230 = rSges(6,1) * t300 + rSges(6,3) * t299;
t486 = t229 + t230;
t321 = t330 * rSges(6,2);
t621 = t520 * t626 + t321;
t565 = pkin(3) * t327;
t368 = t411 + t565;
t466 = qJD(1) * t330;
t620 = t299 * t466 + t300 * t464;
t619 = t633 + t637;
t618 = -t632 + t634;
t617 = t328 * t644 - t330 * t642;
t616 = t328 * t642 + t330 * t644;
t615 = qJD(3) * t663 - t103 * t327 + t105 * t329 + t658 * t299 + t656 * t300;
t614 = t651 * qJD(3) - t102 * t327 + t104 * t329 + t659 * t299 + t657 * t300;
t607 = t661 * qJD(1);
t326 = -qJ(4) - pkin(6);
t442 = t327 * t466;
t467 = qJD(1) * t328;
t301 = qJD(4) * t330;
t436 = t329 * t464;
t481 = pkin(3) * t436 + t301;
t421 = pkin(3) * t442 + t326 * t467 + t481;
t113 = pkin(6) * t467 + t421;
t550 = rSges(5,2) * t299;
t453 = qJD(3) * t550;
t380 = -rSges(5,3) * qJD(1) - t453;
t448 = rSges(5,1) * t620 + t466 * t549;
t93 = t328 * t380 + t448;
t603 = -t113 - t93;
t461 = qJD(5) * t300;
t370 = qJD(3) * t486 - t461;
t601 = t330 * t370;
t285 = t326 * t466;
t290 = pkin(3) * t516;
t279 = t463 * t564;
t462 = qJD(4) * t328;
t424 = -t279 + t462;
t561 = pkin(6) * t330;
t114 = -t285 + (t290 - t561) * qJD(1) + t424;
t261 = t330 * pkin(1) + t328 * qJ(2);
t303 = qJD(2) * t330;
t202 = qJD(1) * t261 - t303;
t592 = -t114 - t202;
t591 = t540 * t518 - t621;
t556 = pkin(6) + t326;
t192 = t328 * t556 + t330 * t565;
t305 = t330 * qJ(2);
t257 = pkin(1) * t328 - t305;
t238 = qJD(1) * t257;
t590 = qJD(1) * t192 - t238;
t562 = pkin(6) * t328;
t429 = -t257 - t562;
t426 = -rSges(3,2) * t330 + t328 * rSges(3,3);
t589 = t261 + t426;
t482 = -t328 * rSges(6,2) - rSges(6,3) * t517;
t588 = -qJ(5) * t517 + t482;
t324 = t328 ^ 2;
t586 = t330 ^ 2 + t324;
t585 = t464 * t650 + t626 * t620;
t459 = qJD(1) * qJD(2);
t302 = qJD(2) * t328;
t478 = qJ(2) * t466 + t302;
t492 = qJD(1) * (-pkin(1) * t467 + t478) + t328 * t459;
t560 = pkin(6) * qJD(1) ^ 2;
t383 = -t328 * t560 + t492;
t423 = qJD(3) * qJD(1) * t564;
t457 = qJD(3) ^ 2 * t565;
t340 = t328 * t423 + t330 * t457 + t383 + (t113 + t301) * qJD(1);
t414 = t289 + t501;
t460 = qJD(5) * t328;
t554 = -(-qJ(5) * t466 - t460) * t300 - t482 * qJD(1) - t585;
t10 = -t414 * t463 + (t328 * t370 - t554) * qJD(1) + t340;
t295 = t330 * t459;
t420 = -t330 * t560 + t295;
t382 = t330 * t423 + t420;
t184 = t230 * t330;
t555 = (pkin(4) * t467 - qJ(5) * t463) * t299 + (-qJ(5) * t467 + (-pkin(4) * qJD(3) + qJD(5)) * t330) * t300 - qJD(3) * t184 + (t328 * t410 + t321) * qJD(1);
t11 = (qJD(3) * t414 - t457) * t328 + (-t462 - t555 + t592 + t601) * qJD(1) + t382;
t584 = t10 * t627 + t11 * t569;
t583 = t586 * t564;
t567 = rSges(3,2) - pkin(1);
t566 = -rSges(6,2) - pkin(1);
t559 = -qJD(1) / 0.2e1;
t557 = -pkin(1) + t326;
t553 = rSges(5,1) * t300;
t548 = rSges(3,3) * t330;
t231 = -t550 + t553;
t191 = t231 * t464;
t204 = t411 * qJD(3);
t26 = t204 * t463 + (t93 + t191) * qJD(1) + t340;
t544 = t26 * t330;
t320 = t330 * rSges(4,3);
t319 = t330 * rSges(5,3);
t316 = t328 * rSges(4,3);
t412 = rSges(4,1) * t327 + t625;
t187 = t412 * t330 - t316;
t260 = rSges(4,1) * t329 - rSges(4,2) * t327;
t225 = t260 * t464;
t74 = t225 + t302 + (t187 + t429) * qJD(1);
t543 = t330 * t74;
t447 = t466 * t625 + (t436 + t442) * rSges(4,1);
t465 = qJD(3) * t327;
t109 = (-rSges(4,2) * t465 - rSges(4,3) * qJD(1)) * t328 + t447;
t236 = t412 * qJD(3);
t42 = t236 * t463 + (t109 + t225) * qJD(1) + t383;
t542 = t42 * t330;
t212 = t260 * t330;
t108 = -qJD(3) * t212 + (t328 * t412 + t320) * qJD(1);
t440 = t260 * t463;
t43 = -t236 * t464 + (-t108 - t202 + t440) * qJD(1) + t420;
t541 = t43 * t328;
t181 = t231 * t328;
t150 = rSges(5,1) * t520 + rSges(5,2) * t518 + t319;
t193 = -t330 * t556 + t290;
t500 = -t150 - t193;
t499 = t519 * t626 + t588;
t494 = t486 * t328;
t493 = -t229 * t330 - t184;
t479 = t279 + t303;
t477 = rSges(3,2) * t467 + rSges(3,3) * t466;
t476 = t302 - t238;
t458 = -rSges(4,3) - pkin(1) - pkin(6);
t456 = -t113 + t554;
t452 = -t193 + t591;
t451 = t192 + t499;
t446 = t302 + t481;
t445 = t285 + t479;
t186 = rSges(4,1) * t516 + rSges(4,2) * t515 + t320;
t444 = t290 + t261;
t437 = t327 * t464;
t433 = -t464 / 0.2e1;
t431 = -t463 / 0.2e1;
t428 = t261 + t561;
t425 = t141 - t525;
t422 = t191 + t446;
t419 = t192 + t429;
t418 = t193 + t428;
t75 = -t440 - t303 + (t186 + t428) * qJD(1);
t409 = t328 * t74 - t330 * t75;
t381 = t421 + t478;
t185 = t231 * t330;
t379 = -t231 * t463 + t462;
t315 = t328 * rSges(5,3);
t152 = -t315 + t624;
t371 = -t152 * t330 + t328 * t500;
t94 = (-t186 * t328 - t187 * t330) * qJD(3);
t369 = t552 + t563 + t565;
t158 = t330 * t192;
t363 = -t158 + t371;
t341 = -t300 * t460 + t464 * t486 + t446;
t291 = pkin(3) * t515;
t258 = rSges(3,2) * t328 + t548;
t211 = t260 * t328;
t157 = t192 * t463;
t125 = qJD(1) * t589 - t303;
t124 = t302 + (-t257 + t258) * qJD(1);
t111 = t330 * t114;
t110 = t114 * t463;
t107 = t295 + (-qJD(1) * t426 - t202) * qJD(1);
t106 = qJD(1) * t477 + t492;
t91 = -qJD(3) * t185 + (t328 * t411 + t319) * qJD(1);
t49 = (t150 + t418) * qJD(1) + t379 - t479;
t48 = (t152 + t419) * qJD(1) + t422;
t47 = qJD(3) * t371 - t157;
t37 = -t303 - t601 + (t418 - t591) * qJD(1) + t424;
t36 = (t419 + t499) * qJD(1) + t341;
t32 = -t157 + t289 + (t328 * t452 - t330 * t499) * qJD(3);
t27 = (-qJD(3) * t204 - t457) * t328 + (-t379 - t91 + t592) * qJD(1) + t382;
t1 = t110 + (t461 + t555 * t330 + t456 * t328 + (t328 * t451 + t330 * t452) * qJD(1)) * qJD(3);
t2 = [(-t662 * qJD(3) + t234 * t327 - t235 * t329 - t667 * t299 + t666 * t300) * qJD(1) + (t11 * (t305 + t588) + t36 * t445 + t10 * (t444 + t621) + t37 * (t381 + t585) + (t11 * t557 - t36 * qJD(4) + (-t37 * qJD(5) - t10 * t540) * t300) * t328 + (-t10 * t326 + t11 * t369 + (-t461 + (t300 * t626 + t650) * qJD(3)) * t36) * t330 + ((-t300 * t37 * t540 + t36 * t566) * t330 + (t36 * (-qJ(2) - t369 + t526 + t547) + t37 * t566) * t328) * qJD(1) - (-t36 + (t499 - t562) * qJD(1) + t341 + t590) * t37) * m(6) + ((-t326 * t330 + t150 + t444) * t26 + (t328 * t557 + t330 * t368 + t305 - t315) * t27 + (t445 + (-pkin(1) * qJD(1) + qJD(3) * t553 + t380) * t330 + (-qJD(4) + (-qJ(2) - t368) * qJD(1)) * t328) * t48 + (t48 - (t152 - t562) * qJD(1) - t422 - t590 + t381 + t448 + (-t453 + (-rSges(5,3) - pkin(1)) * qJD(1)) * t328) * t49) * m(5) + (t43 * (-t316 + t429) + t74 * t303 + t42 * (t186 + t261) + t75 * (-rSges(4,2) * t437 + t447 + t478) + (qJD(3) * t260 * t74 + t42 * pkin(6) + t412 * t43) * t330 + (t458 * t543 + (t74 * (-qJ(2) - t412) + t75 * t458) * t328) * qJD(1) - (t225 - t74 + (t187 - t562) * qJD(1) + t476) * t75) * m(4) + (t107 * (t328 * t567 + t305 + t548) + t124 * t303 + t106 * t589 + t125 * (t477 + t478) + (t124 * t567 * t330 + (t124 * (-rSges(3,3) - qJ(2)) - t125 * pkin(1)) * t328) * qJD(1) - (qJD(1) * t258 - t124 + t476) * t125) * m(3) + (((t511 + t638 - t641 - t660) * t330 + ((t425 - t608 + t623) * t330 - t610 + t613 - t654) * t328) * qJD(3) + t637) * t433 + (((t328 * t425 - t50 + t511) * t328 + t623 * t324 + ((t608 - t669) * t330 + t610 + t613) * t330) * qJD(3) + t618 + t632) * t431 + (t614 + t617 + t619) * t464 / 0.2e1 + (t612 * qJD(1) + t615 + t616) * t463 / 0.2e1 + (t652 * t330 + (-t611 + t653) * t328) * qJD(3) * t559; 0.2e1 * t584 * m(6) + 0.2e1 * (-t544 / 0.2e1 + t27 * t569) * m(5) + 0.2e1 * (t541 / 0.2e1 - t542 / 0.2e1) * m(4) + 0.2e1 * (t106 * t627 + t107 * t569) * m(3); (t11 * t291 - t1 * t158 - t36 * t299 * t460 - (-t36 * t493 + t37 * t494) * qJD(1) + (t1 * t452 + (qJD(1) * t37 + t11) * t486 + (-pkin(3) * t465 + t501 - (t622 - t565) * qJD(3)) * t36) * t328 + (t10 * (-t486 - t564) - t1 * t499 + t36 * t486 * qJD(1)) * t330 + (-t461 - (-t494 * t328 + t493 * t330 - t583) * qJD(3) + t111 + (qJD(1) * t452 + t555) * t330 + (qJD(1) * t451 + t456) * t328) * t32) * m(6) + (-(t49 * t624 + t47 * (-t185 * t330 - t583) + (-t47 * t181 - t48 * t368) * t328) * qJD(3) + t27 * (t291 + t181) + t48 * (-pkin(3) * t437 - t204 * t328) + (-t231 - t564) * t544 + t49 * t204 * t330 + (t463 * t91 + t603 * t464 + t110) * t363 + t47 * (t603 * t328 + t330 * t91 + t111) + (-t48 * t185 - t49 * t181 + (t328 * t49 + t330 * t48) * t231 + (qJD(3) * t363 + t47) * (t500 * t330 + (t152 + t192) * t328)) * qJD(1)) * m(5) + (-(t211 * t75 + t212 * t74) * qJD(1) - (t94 * (-t211 * t328 - t212 * t330) - t409 * t412) * qJD(3) + 0.2e1 * t94 * (t108 * t330 - t109 * t328 + (-t186 * t330 + t187 * t328) * qJD(1)) - t409 * t236 + (t541 - t542 + (t328 * t75 + t543) * qJD(1)) * t260) * m(4) + ((t649 * t299 + t655 * t300 + t327 * t347 + t329 * t348) * qJD(3) + (t299 * t645 - t300 * t646 - t327 * t485 - t329 * t484) * qJD(1)) * t559 + (t615 * t330 + t614 * t328 + (t611 * t328 + t612 * t330) * qJD(1)) * qJD(1) / 0.2e1 + ((t464 * t606 + t607) * t328 + ((-t328 * t609 + t629) * qJD(3) + t640) * t330) * t433 + ((-t463 * t609 + t607) * t330 + ((t330 * t606 - t629) * qJD(3) - t640) * t328) * t431 + (t617 * qJD(1) + ((t638 * qJD(1) + t631 * t330) * t330 + (t635 * t328 - t654 * qJD(1) + (-t630 + t636) * t330) * t328) * t628) * t569 + (t616 * qJD(1) + ((t613 * qJD(1) + t636 * t330) * t330 + (t630 * t328 - t668 * qJD(1) + (-t631 + t635) * t330) * t328) * t628) * t330 / 0.2e1 - (t619 + t633) * t467 / 0.2e1 + (t618 + t634) * t466 / 0.2e1; m(5) * (t26 * t328 + t27 * t330) + m(6) * (t10 * t328 + t11 * t330); (t1 * t299 + 0.2e1 * (-t584 + (0.1e1 / 0.2e1 - t586 / 0.2e1) * qJD(3) * t32) * t300) * m(6);];
tauc = t2(:);
