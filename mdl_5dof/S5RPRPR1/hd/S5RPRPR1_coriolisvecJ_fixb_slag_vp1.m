% Calculate vector of centrifugal and Coriolis load on the joints for
% S5RPRPR1
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% qJD [5x1]
%   Generalized joint velocities
% pkin [8x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d1,d3,d5,theta4]';
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
% Datum: 2019-12-05 17:48
% Revision: 77da58f92bca3eff71542919beafa37024070d86 (2019-12-05)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauc = S5RPRPR1_coriolisvecJ_fixb_slag_vp1(qJ, qJD, ...
  pkin, m, rSges, Icges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(8,1),zeros(6,1),zeros(6,3),zeros(6,6)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RPRPR1_coriolisvecJ_fixb_slag_vp1: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RPRPR1_coriolisvecJ_fixb_slag_vp1: qJD has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S5RPRPR1_coriolisvecJ_fixb_slag_vp1: pkin has to be [8x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5RPRPR1_coriolisvecJ_fixb_slag_vp1: m has to be [6x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [6,3]), ...
  'S5RPRPR1_coriolisvecJ_fixb_slag_vp1: rSges has to be [6x3] (double)');
assert(isreal(Icges) && all(size(Icges) == [6 6]), ...
  'S5RPRPR1_coriolisvecJ_fixb_slag_vp1: Icges has to be [6x6] (double)'); 

%% Symbolic Calculation
% From coriolisvec_joint_fixb_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-05 17:47:03
% EndTime: 2019-12-05 17:47:39
% DurationCPUTime: 25.52s
% Computational Cost: add. (10875->690), mult. (13896->861), div. (0->0), fcn. (10666->8), ass. (0->384)
t664 = Icges(4,3) + Icges(5,3);
t345 = qJ(3) + pkin(8);
t318 = sin(t345);
t319 = cos(t345);
t349 = sin(qJ(3));
t351 = cos(qJ(3));
t663 = Icges(4,5) * t349 + Icges(5,5) * t318 + Icges(4,6) * t351 + Icges(5,6) * t319;
t350 = sin(qJ(1));
t352 = cos(qJ(1));
t660 = t663 * t350 + t664 * t352;
t662 = Icges(4,5) * t351 + Icges(5,5) * t319 - Icges(4,6) * t349 - Icges(5,6) * t318;
t546 = Icges(5,4) * t318;
t417 = Icges(5,2) * t319 + t546;
t166 = -Icges(5,6) * t350 + t352 * t417;
t545 = Icges(5,4) * t319;
t420 = Icges(5,1) * t318 + t545;
t168 = -Icges(5,5) * t350 + t352 * t420;
t548 = Icges(4,4) * t349;
t418 = Icges(4,2) * t351 + t548;
t196 = -Icges(4,6) * t350 + t352 * t418;
t547 = Icges(4,4) * t351;
t421 = Icges(4,1) * t349 + t547;
t198 = -Icges(4,5) * t350 + t352 * t421;
t652 = t166 * t319 + t168 * t318 + t196 * t351 + t198 * t349;
t661 = t652 * t352;
t604 = -t664 * t350 + t663 * t352;
t242 = -Icges(5,2) * t318 + t545;
t244 = Icges(5,1) * t319 - t546;
t276 = -Icges(4,2) * t349 + t547;
t278 = Icges(4,1) * t351 - t548;
t659 = t242 * t319 + t244 * t318 + t276 * t351 + t278 * t349;
t165 = Icges(5,6) * t352 + t350 * t417;
t533 = t319 * t350;
t288 = Icges(5,4) * t533;
t534 = t318 * t350;
t541 = Icges(5,5) * t352;
t167 = Icges(5,1) * t534 + t288 + t541;
t195 = Icges(4,6) * t352 + t350 * t418;
t530 = t350 * t351;
t304 = Icges(4,4) * t530;
t532 = t349 * t350;
t542 = Icges(4,5) * t352;
t197 = Icges(4,1) * t532 + t304 + t542;
t601 = t165 * t319 + t167 * t318 + t195 * t351 + t197 * t349;
t605 = t660 * t350;
t622 = t662 * t350;
t617 = t662 * t352;
t657 = t165 * t533 + t167 * t534 + t195 * t530 + t197 * t532 + t660 * t352;
t614 = -t166 * t533 - t168 * t534 - t196 * t530 - t198 * t532 - t604 * t352;
t656 = -t601 * t352 + t605;
t646 = -t604 * t350 + t661;
t627 = t166 * t318 - t168 * t319 + t196 * t349 - t198 * t351;
t655 = t659 * t350 + t617;
t654 = t659 * t352 - t622;
t620 = t165 * t318 - t167 * t319 + t195 * t349 - t197 * t351;
t653 = t659 * qJD(1) - qJD(3) * t663;
t651 = t604 * qJD(1);
t650 = t660 * qJD(1);
t216 = t417 * qJD(3);
t217 = t420 * qJD(3);
t254 = t418 * qJD(3);
t255 = t421 * qJD(3);
t649 = t216 * t319 + t217 * t318 + t254 * t351 + t255 * t349 + (t242 * t318 - t244 * t319 + t276 * t349 - t278 * t351) * qJD(3) + t662 * qJD(1);
t575 = pkin(3) * t351;
t645 = t655 * qJD(1);
t644 = t601 * qJD(1) + t622 * qJD(3) + t651;
t643 = -t652 * qJD(1) - t617 * qJD(3) + t650;
t642 = (t646 * t350 + t656 * t352) * qJD(3);
t641 = (t614 * t350 + t657 * t352) * qJD(3);
t640 = t654 * qJD(1);
t483 = qJD(3) * t350;
t102 = qJD(1) * t166 + t242 * t483;
t203 = t244 * t350;
t104 = qJD(1) * t168 + qJD(3) * t203;
t116 = qJD(1) * t196 + t276 * t483;
t233 = t278 * t350;
t118 = qJD(1) * t198 + qJD(3) * t233;
t639 = t620 * qJD(3) - t102 * t319 - t104 * t318 - t116 * t351 - t118 * t349 + t650;
t202 = t242 * t352;
t101 = qJD(1) * t165 - qJD(3) * t202;
t204 = t244 * t352;
t103 = -qJD(3) * t204 + (t350 * t420 + t541) * qJD(1);
t232 = t276 * t352;
t115 = qJD(1) * t195 - qJD(3) * t232;
t234 = t278 * t352;
t117 = -qJD(3) * t234 + (t350 * t421 + t542) * qJD(1);
t638 = t627 * qJD(3) + t101 * t319 + t103 * t318 + t115 * t351 + t117 * t349 + t651;
t504 = t276 + t421;
t505 = -t418 + t278;
t506 = t242 + t420;
t507 = -t417 + t244;
t637 = (t318 * t506 - t319 * t507 + t349 * t504 - t351 * t505) * qJD(1);
t320 = qJ(5) + t345;
t310 = sin(t320);
t311 = cos(t320);
t535 = t311 * t350;
t266 = Icges(6,4) * t535;
t536 = t310 * t350;
t540 = Icges(6,5) * t352;
t155 = Icges(6,1) * t536 + t266 + t540;
t344 = qJD(3) + qJD(5);
t271 = t344 * t350;
t272 = t344 * t352;
t544 = Icges(6,4) * t310;
t225 = Icges(6,1) * t311 - t544;
t416 = Icges(6,2) * t311 + t544;
t619 = -t416 + t225;
t543 = Icges(6,4) * t311;
t419 = Icges(6,1) * t310 + t543;
t156 = -Icges(6,5) * t350 + t352 * t419;
t223 = -Icges(6,2) * t310 + t543;
t623 = t223 * t352 + t156;
t587 = qJD(1) * t619 - t271 * t623 + t272 * (-Icges(6,2) * t536 + t155 + t266);
t621 = t223 + t419;
t154 = -Icges(6,6) * t350 + t352 * t416;
t624 = -t225 * t352 + t154;
t153 = Icges(6,6) * t352 + t350 * t416;
t625 = -t225 * t350 + t153;
t588 = qJD(1) * t621 - t271 * t624 + t272 * t625;
t636 = t310 * t588 - t311 * t587;
t635 = 0.2e1 * qJD(3);
t573 = pkin(4) * t319;
t634 = t641 + t645;
t633 = -t640 + t642;
t632 = t653 * t350 + t649 * t352;
t631 = -t649 * t350 + t653 * t352;
t630 = -t601 * qJD(3) - t102 * t318 + t104 * t319 - t116 * t349 + t118 * t351;
t629 = t652 * qJD(3) - t101 * t318 + t103 * t319 - t115 * t349 + t117 * t351;
t628 = rSges(4,2) * t351;
t425 = rSges(5,1) * t318 + rSges(5,2) * t319;
t626 = t352 * t425;
t576 = pkin(3) * t349;
t386 = t425 + t576;
t618 = t663 * qJD(1);
t424 = rSges(6,1) * t310 + rSges(6,2) * t311;
t370 = t350 * (t198 + t232) - t352 * (-Icges(4,2) * t532 + t197 + t304);
t371 = t350 * (t196 - t234) - t352 * (t195 - t233);
t373 = t350 * (t168 + t202) - t352 * (-Icges(5,2) * t534 + t167 + t288);
t374 = t350 * (t166 - t204) - t352 * (t165 - t203);
t615 = -t374 * t318 + t373 * t319 - t371 * t349 + t370 * t351;
t413 = Icges(6,5) * t310 + Icges(6,6) * t311;
t152 = -Icges(6,3) * t350 + t352 * t413;
t60 = -t352 * t152 - t154 * t535 - t156 * t536;
t404 = t223 * t311 + t225 * t310;
t221 = Icges(6,5) * t311 - Icges(6,6) * t310;
t539 = t221 * t352;
t75 = t350 * t404 + t539;
t613 = t75 * qJD(1) + t271 * t60;
t411 = t154 * t311 + t156 * t310;
t612 = t352 * t411;
t485 = qJD(3) * t318;
t470 = rSges(5,2) * t485;
t396 = -rSges(5,3) * qJD(1) - t470;
t458 = t319 * t483;
t486 = qJD(1) * t352;
t466 = rSges(5,1) * t458 + t425 * t486;
t106 = t350 * t396 + t466;
t321 = qJD(4) * t352;
t456 = t351 * t483;
t295 = pkin(3) * t456;
t348 = -qJ(4) - pkin(6);
t461 = t349 * t486;
t487 = qJD(1) * t350;
t463 = -pkin(3) * t461 - t348 * t487 - t295;
t434 = t321 - t463;
t132 = pkin(6) * t487 + t434;
t607 = -t106 - t132;
t308 = pkin(3) * t532;
t481 = qJD(4) * t350;
t482 = qJD(3) * t352;
t297 = t482 * t575;
t303 = t348 * t486;
t502 = t297 + t303;
t571 = pkin(6) * t352;
t133 = t481 + (t308 - t571) * qJD(1) - t502;
t286 = t352 * pkin(1) + t350 * qJ(2);
t323 = qJD(2) * t352;
t218 = qJD(1) * t286 - t323;
t606 = -t133 - t218;
t477 = t352 * t576;
t567 = pkin(6) + t348;
t211 = t350 * t567 + t477;
t325 = t352 * qJ(2);
t282 = pkin(1) * t350 - t325;
t260 = qJD(1) * t282;
t603 = qJD(1) * t211 - t260;
t572 = pkin(6) * t350;
t448 = -t282 - t572;
t445 = -rSges(3,2) * t352 + t350 * rSges(3,3);
t602 = t286 + t445;
t346 = t350 ^ 2;
t600 = (-t352 ^ 2 - t346) * t575;
t439 = t621 * t344;
t440 = t619 * t344;
t591 = qJD(1) * t221 + t310 * t439 - t311 * t440;
t441 = -qJD(1) * t153 + t344 * t623;
t443 = (t350 * t419 + t540) * qJD(1) + t624 * t344;
t495 = qJD(1) * t152;
t590 = t310 * t443 - t311 * t441 + t495;
t442 = qJD(1) * t154 + t155 * t344 + t223 * t271;
t444 = -qJD(1) * t156 + t344 * t625;
t151 = Icges(6,3) * t352 + t350 * t413;
t496 = qJD(1) * t151;
t589 = t310 * t444 - t311 * t442 + t496;
t251 = t344 * t487;
t586 = -t251 / 0.2e1;
t252 = qJD(1) * t272;
t585 = t252 / 0.2e1;
t584 = -t271 / 0.2e1;
t583 = t271 / 0.2e1;
t582 = -t272 / 0.2e1;
t581 = t272 / 0.2e1;
t580 = t350 / 0.2e1;
t579 = t352 / 0.2e1;
t578 = rSges(3,2) - pkin(1);
t577 = -rSges(6,3) - pkin(1);
t574 = pkin(4) * t318;
t570 = pkin(6) * qJD(1) ^ 2;
t569 = -qJD(1) / 0.2e1;
t568 = qJD(1) / 0.2e1;
t565 = rSges(5,1) * t319;
t563 = rSges(6,1) * t311;
t560 = rSges(6,2) * t310;
t558 = rSges(3,3) * t352;
t175 = t424 * t344;
t227 = -t560 + t563;
t268 = pkin(4) * t458;
t353 = qJD(3) ^ 2;
t479 = qJD(1) * qJD(2);
t315 = qJ(2) * t486;
t322 = qJD(2) * t350;
t500 = t315 + t322;
t511 = qJD(1) * (-pkin(1) * t487 + t500) + t350 * t479;
t399 = -t350 * t570 + t511;
t436 = qJD(3) * qJD(1) * t575;
t476 = t353 * t576;
t362 = t350 * t436 + t352 * t476 + t399 + (t132 + t321) * qJD(1);
t262 = t573 + t575;
t245 = t262 * qJD(3);
t261 = t574 + t576;
t343 = -pkin(7) + t348;
t468 = t350 * t245 + t261 * t486 + t343 * t487;
t79 = t463 + t468;
t475 = t344 * t563;
t467 = t350 * t475 + t424 * t486;
t474 = t344 * t560;
t96 = (-rSges(6,3) * qJD(1) - t474) * t350 + t467;
t20 = t352 * t353 * t574 + t175 * t272 + t227 * t251 + (t79 + t96 + t268) * qJD(1) + t362;
t557 = t20 * t352;
t237 = t350 * t261;
t387 = -t482 * t573 + t481;
t314 = t352 * t479;
t433 = -t352 * t570 + t314;
t398 = t352 * t436 + t433;
t80 = -t245 * t352 + (-t343 * t352 + (t261 - t576) * t350) * qJD(1) + t502;
t189 = t227 * t352;
t339 = t352 * rSges(6,3);
t95 = -t344 * t189 + (t350 * t424 + t339) * qJD(1);
t21 = -t175 * t271 + t227 * t252 - t353 * t237 + (-t387 - t80 - t95 + t606) * qJD(1) + t398;
t556 = t21 * t350;
t341 = t352 * rSges(4,3);
t340 = t352 * rSges(5,3);
t498 = t321 + t322;
t462 = t295 + t498;
t389 = t271 * t227 + t268 + t462;
t432 = t211 + t448;
t508 = -t352 * t261 - t350 * t343;
t130 = t348 * t350 + t477 + t508;
t335 = t350 * rSges(6,3);
t158 = t352 * t424 - t335;
t527 = -t130 + t158;
t48 = (t432 + t527) * qJD(1) + t389;
t554 = t352 * t48;
t337 = t350 * rSges(4,3);
t426 = rSges(4,1) * t349 + t628;
t208 = t352 * t426 - t337;
t285 = rSges(4,1) * t351 - rSges(4,2) * t349;
t246 = t285 * t483;
t97 = t246 + t322 + (t208 + t448) * qJD(1);
t553 = t352 * t97;
t248 = -rSges(5,2) * t318 + t565;
t210 = t248 * t483;
t219 = t425 * qJD(3);
t41 = t219 * t482 + (t106 + t210) * qJD(1) + t362;
t552 = t41 * t352;
t464 = t486 * t628 + (t456 + t461) * rSges(4,1);
t484 = qJD(3) * t349;
t122 = (-rSges(4,2) * t484 - rSges(4,3) * qJD(1)) * t350 + t464;
t256 = t426 * qJD(3);
t54 = t256 * t482 + (t122 + t246) * qJD(1) + t399;
t551 = t54 * t352;
t236 = t285 * t352;
t121 = -qJD(3) * t236 + (t350 * t426 + t341) * qJD(1);
t459 = t285 * t482;
t55 = -t256 * t483 + (-t121 - t218 + t459) * qJD(1) + t433;
t550 = t55 * t350;
t549 = -t132 - t79;
t205 = t248 * t350;
t531 = t350 * t152;
t182 = t350 * t221;
t76 = t352 * t404 - t182;
t529 = t76 * qJD(1);
t526 = -t130 + t211;
t131 = t237 - t308 + (-t343 + t348) * t352;
t212 = -t352 * t567 + t308;
t525 = -t131 - t212;
t173 = rSges(5,1) * t534 + rSges(5,2) * t533 + t340;
t516 = -t173 - t212;
t501 = t297 + t323;
t499 = rSges(3,2) * t487 + rSges(3,3) * t486;
t497 = t322 - t260;
t478 = -rSges(4,3) - pkin(1) - pkin(6);
t59 = t352 * t151 + t153 * t535 + t155 * t536;
t157 = rSges(6,1) * t536 + rSges(6,2) * t535 + t339;
t469 = -t157 + t525;
t207 = rSges(4,1) * t532 + rSges(4,2) * t530 + t341;
t457 = t349 * t483;
t455 = -t487 / 0.2e1;
t454 = t486 / 0.2e1;
t453 = -t483 / 0.2e1;
t451 = -t482 / 0.2e1;
t449 = t227 + t573;
t447 = t286 + t571;
t188 = t227 * t350;
t438 = qJD(1) * t188 + t272 * t424;
t437 = qJD(1) * t189 - t271 * t424;
t435 = t210 + t462;
t431 = t212 + t447;
t98 = -t459 - t323 + (t207 + t447) * qJD(1);
t423 = t350 * t97 - t352 * t98;
t179 = t211 * t482;
t38 = -t157 * t271 - t158 * t272 - t179 + (t130 * t352 + t350 * t525) * qJD(3);
t422 = t38 * (-t188 * t271 - t272 * t189);
t412 = t153 * t311 + t155 * t310;
t78 = t154 * t310 - t156 * t311;
t397 = t59 + t531;
t206 = t248 * t352;
t395 = -t248 * t482 + t481;
t336 = t350 * rSges(5,3);
t174 = -t336 + t626;
t388 = -t352 * t174 + t350 * t516;
t107 = (-t207 * t350 - t208 * t352) * qJD(3);
t190 = t352 * t211;
t385 = -t190 + t388;
t380 = -qJD(1) * t413 + t182 * t272 - t271 * t539;
t141 = t350 * t151;
t61 = -t352 * t412 + t141;
t377 = -qJD(1) * t411 - t344 * t539 + t496;
t376 = qJD(1) * t412 + t182 * t344 + t495;
t366 = qJD(1) * t404 - t413 * t344;
t10 = t377 * t350 - t352 * t590;
t11 = -t350 * t589 + t376 * t352;
t12 = t350 * t590 + t377 * t352;
t22 = t272 * t59 + t613;
t62 = -t531 + t612;
t23 = t271 * t62 + t272 * t61 - t529;
t30 = t366 * t350 + t352 * t591;
t31 = -t350 * t591 + t366 * t352;
t36 = -t310 * t442 - t311 * t444;
t37 = t310 * t441 + t311 * t443;
t77 = -t153 * t310 + t155 * t311;
t9 = t376 * t350 + t352 * t589;
t363 = (qJD(1) * t30 + t10 * t271 - t251 * t61 + t252 * t62 + t272 * t9) * t580 + (-t310 * t587 - t311 * t588) * t569 + t22 * t455 + t23 * t454 + (qJD(1) * t31 + t11 * t272 + t12 * t271 - t251 * t59 + t252 * t60) * t579 + (t10 * t350 + t352 * t9 + (-t350 * t61 + t352 * t62) * qJD(1)) * t583 + (t350 * t60 + t352 * t59) * t586 + (t350 * t62 + t352 * t61) * t585 + (t11 * t352 + t12 * t350 + (-t350 * t59 + t352 * t60) * qJD(1)) * t581 + (t350 * t37 + t352 * t36 + (-t350 * t77 + t352 * t78) * qJD(1)) * t568 + (t380 * t350 + t352 * t636) * t584 + (-t350 * t636 + t380 * t352) * t582;
t309 = pkin(3) * t530;
t283 = rSges(3,2) * t350 + t558;
t235 = t285 * t350;
t181 = (-t262 + t575) * t352;
t180 = t262 * t350 - t309;
t146 = qJD(1) * t602 - t323;
t145 = t322 + (-t282 + t283) * qJD(1);
t144 = t352 * t158;
t124 = t352 * t133;
t123 = t133 * t482;
t120 = t314 + (-qJD(1) * t445 - t218) * qJD(1);
t119 = qJD(1) * t499 + t511;
t105 = -qJD(3) * t206 + (t350 * t425 + t340) * qJD(1);
t86 = t352 * t95;
t65 = (t173 + t431) * qJD(1) + t395 - t501;
t64 = (t174 + t432) * qJD(1) + t435;
t63 = qJD(3) * t388 - t179;
t49 = -t227 * t272 + (t131 + t157 + t431) * qJD(1) + t387 - t501;
t42 = (-qJD(3) * t219 - t476) * t350 + (-t105 - t395 + t606) * qJD(1) + t398;
t5 = -t157 * t252 + t158 * t251 - t271 * t96 + t272 * t95 + t123 + (t352 * t80 + t549 * t350 + (t350 * t526 + t352 * t525) * qJD(1)) * qJD(3);
t1 = [((t397 + t62 - t612) * t272 + t613) * t584 + t78 * t585 - t252 * t76 / 0.2e1 + (t77 + t75) * t586 + (t529 + (t350 * t411 - t141 + t60) * t272 + (t397 - t59) * t271 + ((t152 + t412) * t272 - t411 * t271) * t352 + t23) * t582 + (t36 + t31) * t581 + (((t605 + t614 - t656) * t350 + (-t661 + (-t601 + t604) * t350 + t646 + t657) * t352) * qJD(3) + t645) * t453 + (-t659 * qJD(3) + t216 * t318 - t217 * t319 + t254 * t349 - t255 * t351 - t310 * t440 - t311 * t439) * qJD(1) + (-(-t48 + (t527 - t572) * qJD(1) + t389 + t603) * t49 + t21 * (-t282 - t335 - t508) + t48 * (t323 - t481) + t20 * (t237 + t157 + t286) + t49 * (-t350 * t474 + t315 + t467 + t468 + t498) + (t21 * t424 + t48 * (t245 - t474 + t475) - t20 * t343) * t352 + ((t343 + t577) * t554 + (t48 * (-qJ(2) - t261 - t424) + t49 * t577) * t350) * qJD(1)) * m(6) + ((-t348 * t352 + t173 + t286 + t308) * t41 + (t325 - t336 + t386 * t352 + (-pkin(1) + t348) * t350) * t42 + (t303 + t501 + (-pkin(1) * qJD(1) + qJD(3) * t565 + t396) * t352 + (-qJD(4) + (-qJ(2) - t386) * qJD(1)) * t350) * t64 + (t64 - (t174 - t572) * qJD(1) - t435 + t434 + t466 + t500 + (-t470 + (-rSges(5,3) - pkin(1)) * qJD(1)) * t350 - t603) * t65) * m(5) + (-(t246 - t97 + (t208 - t572) * qJD(1) + t497) * t98 + t55 * (-t337 + t448) + t97 * t323 + t54 * (t207 + t286) + t98 * (-rSges(4,2) * t457 + t464 + t500) + (qJD(3) * t285 * t97 + t54 * pkin(6) + t426 * t55) * t352 + (t478 * t553 + (t97 * (-qJ(2) - t426) + t98 * t478) * t350) * qJD(1)) * m(4) + (-(qJD(1) * t283 - t145 + t497) * t146 + t120 * (t350 * t578 + t325 + t558) + t145 * t323 + t119 * t602 + t146 * (t499 + t500) + (t145 * t578 * t352 + (t145 * (-rSges(3,3) - qJ(2)) - t146 * pkin(1)) * t350) * qJD(1)) * m(3) + (t37 + t30 + t22) * t583 + ((t604 * t346 + ((t601 + t604) * t352 - t605 + t614) * t352) * qJD(3) + t633 + t640) * t451 + (t629 + t632 + t634) * t483 / 0.2e1 + (qJD(1) * t627 + t630 + t631) * t482 / 0.2e1 + (t654 * t352 + (-t620 + t655) * t350) * qJD(3) * t569; 0.2e1 * (-t557 / 0.2e1 + t556 / 0.2e1) * m(6) + 0.2e1 * (t42 * t580 - t552 / 0.2e1) * m(5) + 0.2e1 * (t550 / 0.2e1 - t551 / 0.2e1) * m(4) + 0.2e1 * (-t119 * t352 / 0.2e1 + t120 * t580) * m(3); t363 + ((t318 * t373 + t319 * t374 + t349 * t370 + t351 * t371) * qJD(3) + (-t318 * t507 - t319 * t506 - t349 * t505 - t351 * t504) * qJD(1)) * t569 + (t630 * t352 + t629 * t350 + (t350 * t620 + t352 * t627) * qJD(1)) * t568 + ((-t483 * t617 - t618) * t350 + ((t350 * t622 + t615) * qJD(3) + t637) * t352) * t453 + ((t482 * t622 - t618) * t352 + ((-t352 * t617 - t615) * qJD(3) - t637) * t350) * t451 + (-t48 * (-qJD(1) * t181 + t437) - t49 * (qJD(1) * t180 + t438) - t422 - (t38 * (t181 * t352 + t600) + (-t38 * t180 - t261 * t48) * t350) * qJD(3) + t21 * t309 + t5 * (-t144 - t190) + t38 * (t124 + t86) + (t20 * (-t227 - t262) + t49 * t175 + t5 * t130 + t38 * t80 + (t38 * t469 + t449 * t48) * qJD(1)) * t352 + (t21 * t449 + t48 * (-pkin(3) * t484 - pkin(4) * t485 - t175) + t5 * t469 + t38 * (-t96 + t549) + (t49 * t449 + t38 * (t158 + t526)) * qJD(1)) * t350) * m(6) + (-(t65 * t626 + t63 * (-t206 * t352 + t600) + (-t63 * t205 - t386 * t64) * t350) * qJD(3) + t42 * (t309 + t205) - t64 * pkin(3) * t457 + (-t248 - t575) * t552 + (t105 * t482 + t607 * t483 + t123) * t385 + t63 * (t105 * t352 + t607 * t350 + t124) + (-t64 * t350 + t65 * t352) * t219 + (-t64 * t206 - t65 * t205 + (t350 * t65 + t352 * t64) * t248 + (qJD(3) * t385 + t63) * (t516 * t352 + (t174 + t211) * t350)) * qJD(1)) * m(5) + (-(t235 * t98 + t236 * t97) * qJD(1) - (t107 * (-t235 * t350 - t236 * t352) - t423 * t426) * qJD(3) + 0.2e1 * t107 * (t352 * t121 - t350 * t122 + (-t207 * t352 + t208 * t350) * qJD(1)) - t423 * t256 + (t550 - t551 + (t350 * t98 + t553) * qJD(1)) * t285) * m(4) + (t632 * qJD(1) + ((t646 * qJD(1) + t639 * t352) * t352 + (t643 * t350 - t656 * qJD(1) + (-t638 + t644) * t352) * t350) * t635) * t580 + (t631 * qJD(1) + ((t614 * qJD(1) + t644 * t352) * t352 + (t638 * t350 - t657 * qJD(1) + (-t639 + t643) * t352) * t350) * t635) * t579 + (t634 + t641) * t455 + (t633 + t642) * t454; m(5) * (t350 * t41 + t352 * t42) + m(6) * (t20 * t350 + t21 * t352); t363 + (t5 * (-t157 * t350 - t144) + t38 * (-t350 * t96 + t86 + (-t157 * t352 + t158 * t350) * qJD(1)) - (t350 * t48 - t352 * t49) * t175 + (-t557 + t556 + (t350 * t49 + t554) * qJD(1)) * t227 - t437 * t48 - t438 * t49 - t422) * m(6);];
tauc = t1(:);
