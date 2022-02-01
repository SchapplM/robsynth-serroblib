% Calculate vector of centrifugal and Coriolis load on the joints for
% S5RPRPR4
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% qJD [5x1]
%   Generalized joint velocities
% pkin [9x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d1,d3,d5,theta2,theta4]';
% m [6x1]
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
% Datum: 2022-01-23 09:23
% Revision: fd3771346c4aea32fdeb66112c511235427c26a7 (2022-01-20)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauc = S5RPRPR4_coriolisvecJ_fixb_slag_vp1(qJ, qJD, ...
  pkin, m, rSges, Icges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(9,1),zeros(6,1),zeros(6,3),zeros(6,6)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RPRPR4_coriolisvecJ_fixb_slag_vp1: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RPRPR4_coriolisvecJ_fixb_slag_vp1: qJD has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S5RPRPR4_coriolisvecJ_fixb_slag_vp1: pkin has to be [9x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5RPRPR4_coriolisvecJ_fixb_slag_vp1: m has to be [6x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [6,3]), ...
  'S5RPRPR4_coriolisvecJ_fixb_slag_vp1: rSges has to be [6x3] (double)');
assert(isreal(Icges) && all(size(Icges) == [6 6]), ...
  'S5RPRPR4_coriolisvecJ_fixb_slag_vp1: Icges has to be [6x6] (double)'); 

%% Symbolic Calculation
% From coriolisvec_joint_fixb_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2022-01-23 09:22:25
% EndTime: 2022-01-23 09:22:59
% DurationCPUTime: 29.01s
% Computational Cost: add. (16678->702), mult. (13398->892), div. (0->0), fcn. (10356->10), ass. (0->411)
t673 = Icges(4,3) + Icges(5,3);
t341 = qJ(3) + pkin(9);
t331 = sin(t341);
t333 = cos(t341);
t241 = Icges(5,5) * t333 - Icges(5,6) * t331;
t344 = sin(qJ(3));
t346 = cos(qJ(3));
t279 = Icges(4,5) * t346 - Icges(4,6) * t344;
t669 = t241 + t279;
t342 = qJ(1) + pkin(8);
t334 = cos(t342);
t672 = t673 * t334;
t332 = sin(t342);
t534 = t332 * t346;
t535 = t332 * t344;
t536 = t332 * t333;
t538 = t331 * t332;
t654 = -Icges(4,5) * t534 - Icges(5,5) * t536 + Icges(4,6) * t535 + Icges(5,6) * t538 + t672;
t653 = t673 * t332 + t669 * t334;
t561 = Icges(5,6) * t334;
t172 = Icges(5,4) * t536 - Icges(5,2) * t538 - t561;
t562 = Icges(4,6) * t334;
t189 = Icges(4,4) * t534 - Icges(4,2) * t535 - t562;
t671 = t172 * t331 + t189 * t344;
t570 = Icges(5,4) * t331;
t245 = Icges(5,1) * t333 - t570;
t175 = Icges(5,5) * t332 + t245 * t334;
t571 = Icges(4,4) * t344;
t283 = Icges(4,1) * t346 - t571;
t192 = Icges(4,5) * t332 + t283 * t334;
t670 = -t175 * t536 - t192 * t534;
t242 = Icges(5,2) * t333 + t570;
t315 = Icges(5,4) * t333;
t244 = Icges(5,1) * t331 + t315;
t280 = Icges(4,2) * t346 + t571;
t336 = Icges(4,4) * t346;
t282 = Icges(4,1) * t344 + t336;
t667 = t242 * t331 - t244 * t333 + t280 * t344 - t282 * t346;
t268 = Icges(5,4) * t538;
t567 = Icges(5,5) * t334;
t174 = Icges(5,1) * t536 - t268 - t567;
t291 = Icges(4,4) * t535;
t568 = Icges(4,5) * t334;
t191 = Icges(4,1) * t534 - t291 - t568;
t657 = -t174 * t333 - t191 * t346 + t671;
t240 = Icges(5,5) * t331 + Icges(5,6) * t333;
t278 = Icges(4,5) * t344 + Icges(4,6) * t346;
t668 = t278 + t240;
t666 = t334 * t653 + t670;
t529 = t334 * t346;
t533 = t333 * t334;
t665 = -t174 * t533 - t191 * t529 + t332 * t654;
t635 = t175 * t533 + t192 * t529 + t332 * t653;
t664 = -t657 * t332 + t334 * t654;
t411 = -Icges(5,2) * t331 + t315;
t173 = Icges(5,6) * t332 + t334 * t411;
t412 = -Icges(4,2) * t344 + t336;
t190 = Icges(4,6) * t332 + t334 * t412;
t643 = -t173 * t538 - t190 * t535 - t666;
t530 = t334 * t344;
t537 = t331 * t334;
t642 = -t172 * t537 - t189 * t530 - t665;
t641 = -t173 * t537 - t190 * t530 + t635;
t543 = t278 * t334;
t545 = t240 * t334;
t663 = -t332 * t667 - t543 - t545;
t544 = t278 * t332;
t546 = t240 * t332;
t662 = -t334 * t667 + t544 + t546;
t661 = t173 * t331 + t190 * t344;
t225 = t411 * qJD(3);
t226 = t245 * qJD(3);
t260 = t412 * qJD(3);
t261 = t283 * qJD(3);
t660 = -t225 * t331 + t226 * t333 - t260 * t344 + t261 * t346 + (-t242 * t333 - t244 * t331 - t280 * t346 - t282 * t344) * qJD(3) + t668 * qJD(1);
t475 = rSges(5,1) * t536;
t343 = -qJ(4) - pkin(6);
t298 = t334 * t343;
t337 = t346 * pkin(3);
t328 = t337 + pkin(2);
t499 = t332 * t328 + t298;
t345 = sin(qJ(1));
t593 = pkin(1) * t345;
t659 = -t475 - t593 - t499;
t658 = t175 * t333 + t192 * t346 - t661;
t656 = t667 * qJD(1) + qJD(3) * t669;
t655 = t662 * qJD(1);
t652 = (t641 * t332 - t642 * t334) * qJD(3);
t651 = (t643 * t332 - t664 * t334) * qJD(3);
t650 = t663 * qJD(1);
t649 = t650 + t651;
t648 = t652 + t655;
t380 = qJD(3) * t242;
t106 = qJD(1) * t173 - t332 * t380;
t382 = qJD(3) * t244;
t108 = qJD(1) * t175 - t332 * t382;
t381 = qJD(3) * t280;
t120 = qJD(1) * t190 - t332 * t381;
t383 = qJD(3) * t282;
t122 = qJD(1) * t192 - t332 * t383;
t647 = t657 * qJD(3) - t106 * t333 - t108 * t331 - t120 * t346 - t122 * t344;
t105 = -t334 * t380 + (-t332 * t411 + t561) * qJD(1);
t107 = -t334 * t382 + (-t245 * t332 + t567) * qJD(1);
t119 = -t334 * t381 + (-t332 * t412 + t562) * qJD(1);
t121 = -t334 * t383 + (-t283 * t332 + t568) * qJD(1);
t646 = t658 * qJD(3) + t105 * t333 + t107 * t331 + t119 * t346 + t121 * t344;
t645 = t656 * t332 + t660 * t334;
t644 = t660 * t332 - t656 * t334;
t640 = t172 * t333 + t174 * t331 + t189 * t346 + t191 * t344;
t639 = t173 * t333 + t175 * t331 + t190 * t346 + t192 * t344;
t324 = t334 * pkin(6);
t251 = pkin(2) * t332 - t324;
t157 = t251 - t499;
t153 = qJD(1) * t157;
t237 = qJD(1) * t251;
t304 = qJD(4) * t332;
t638 = t304 - t153 + t237;
t637 = t668 * qJD(3);
t636 = t654 + t661;
t634 = t653 * qJD(1);
t349 = qJD(1) ^ 2;
t633 = rSges(4,2) * t344;
t453 = -t251 - t593;
t427 = t157 + t453;
t335 = qJ(5) + t341;
t327 = cos(t335);
t581 = rSges(6,2) * t327;
t326 = sin(t335);
t585 = rSges(6,1) * t326;
t235 = t581 + t585;
t340 = qJD(3) + qJD(5);
t239 = t334 * t340;
t433 = -t235 * t239 + t304;
t591 = pkin(4) * t331;
t592 = pkin(3) * t344;
t271 = -t591 - t592;
t246 = t271 * qJD(3);
t532 = t334 * t246;
t590 = pkin(4) * t333;
t262 = t328 + t590;
t339 = pkin(7) - t343;
t531 = t334 * t339;
t128 = t262 * t332 - t499 - t531;
t542 = t326 * t332;
t257 = rSges(6,2) * t542;
t501 = -t334 * rSges(6,3) - t257;
t540 = t327 * t332;
t165 = rSges(6,1) * t540 + t501;
t630 = -t128 - t165;
t48 = t532 + (t427 + t630) * qJD(1) + t433;
t632 = t332 * t48;
t238 = t332 * t340;
t483 = qJD(3) * t344;
t276 = t332 * pkin(3) * t483;
t305 = qJD(4) * t334;
t497 = t276 + t305;
t421 = t235 * t238 + t497;
t485 = qJD(3) * t332;
t431 = t485 * t591;
t316 = t332 * rSges(6,3);
t539 = t327 * t334;
t541 = t326 * t334;
t166 = rSges(6,1) * t539 - rSges(6,2) * t541 + t316;
t323 = t332 * pkin(6);
t252 = t334 * pkin(2) + t323;
t493 = -t334 * t328 + t332 * t343;
t158 = -t252 - t493;
t347 = cos(qJ(1));
t338 = t347 * pkin(1);
t452 = t252 + t338;
t426 = t158 + t452;
t294 = t339 * t332;
t628 = t334 * t262 + t294;
t616 = t493 + t628;
t621 = t616 + t166 + t426;
t49 = t621 * qJD(1) - t421 - t431;
t631 = t334 * t49;
t569 = Icges(6,4) * t326;
t233 = Icges(6,1) * t327 - t569;
t164 = Icges(6,5) * t332 + t233 * t334;
t230 = Icges(6,2) * t327 + t569;
t629 = -t230 * t334 + t164;
t627 = -t230 + t233;
t302 = Icges(6,4) * t327;
t232 = Icges(6,1) * t326 + t302;
t410 = -Icges(6,2) * t326 + t302;
t626 = t232 + t410;
t625 = -qJD(3) * t639 - t105 * t331 + t107 * t333 - t119 * t344 + t121 * t346 + t634;
t624 = qJD(1) * t654 + t640 * qJD(3) + t106 * t331 - t108 * t333 + t120 * t344 - t122 * t346;
t623 = t657 * qJD(1) - t637 * t332 + t634;
t622 = -t637 * t334 + (-t332 * t669 - t658 + t672) * qJD(1);
t620 = 0.2e1 * qJD(3);
t318 = t332 * rSges(4,3);
t196 = rSges(4,1) * t529 - rSges(4,2) * t530 + t318;
t618 = t196 + t452;
t269 = rSges(5,2) * t538;
t617 = -t334 * rSges(5,3) - t269;
t446 = t334 * rSges(3,1) - rSges(3,2) * t332;
t615 = t338 + t446;
t229 = Icges(6,5) * t327 - Icges(6,6) * t326;
t228 = Icges(6,5) * t326 + Icges(6,6) * t327;
t547 = t228 * t334;
t162 = Icges(6,6) * t332 + t334 * t410;
t555 = t162 * t326;
t557 = Icges(6,3) * t334;
t613 = -t340 * t547 + (-t164 * t327 - t229 * t332 + t555 + t557) * qJD(1);
t560 = Icges(6,6) * t334;
t161 = Icges(6,4) * t540 - Icges(6,2) * t542 - t560;
t256 = Icges(6,4) * t542;
t566 = Icges(6,5) * t334;
t163 = Icges(6,1) * t540 - t256 - t566;
t409 = t161 * t326 - t163 * t327;
t160 = Icges(6,3) * t332 + t229 * t334;
t492 = qJD(1) * t160;
t548 = t228 * t332;
t612 = qJD(1) * t409 - t340 * t548 + t492;
t402 = t230 * t326 - t232 * t327;
t607 = qJD(1) * t402 + t229 * t340;
t604 = t332 * (-t242 * t334 + t175) - t334 * (-Icges(5,2) * t536 + t174 - t268);
t508 = -Icges(4,2) * t534 + t191 - t291;
t510 = t282 * t332 + t189;
t603 = -t344 * t508 - t346 * t510;
t602 = qJD(1) * t626 + t238 * t629 - t239 * (-Icges(6,2) * t540 + t163 - t256);
t430 = qJD(1) * t340;
t222 = t332 * t430;
t601 = t222 / 0.2e1;
t223 = t334 * t430;
t600 = t223 / 0.2e1;
t599 = -t238 / 0.2e1;
t598 = t238 / 0.2e1;
t597 = -t239 / 0.2e1;
t596 = t239 / 0.2e1;
t595 = t332 / 0.2e1;
t594 = -t334 / 0.2e1;
t589 = -qJD(1) / 0.2e1;
t588 = qJD(1) / 0.2e1;
t587 = pkin(2) - t328;
t586 = rSges(4,1) * t346;
t584 = rSges(6,1) * t327;
t582 = rSges(5,2) * t333;
t285 = rSges(4,1) * t344 + rSges(4,2) * t346;
t219 = t285 * t334;
t465 = t285 * t485;
t98 = qJD(1) * t618 - t465;
t580 = t219 * t98;
t317 = t332 * rSges(5,3);
t494 = rSges(4,2) * t535 + t334 * rSges(4,3);
t195 = rSges(4,1) * t534 - t494;
t484 = qJD(3) * t334;
t463 = t285 * t484;
t97 = -t463 + (-t195 + t453) * qJD(1);
t579 = t332 * t97;
t578 = t334 * t48;
t577 = t334 * t97;
t248 = rSges(5,1) * t331 + t582;
t176 = t475 + t617;
t451 = -t248 - t592;
t395 = t451 * t484;
t375 = t304 + t395;
t62 = (-t176 + t427) * qJD(1) + t375;
t576 = t62 * t248;
t159 = Icges(6,5) * t540 - Icges(6,6) * t542 - t557;
t556 = t159 * t334;
t79 = -t332 * t402 - t547;
t528 = t79 * qJD(1);
t487 = qJD(1) * t332;
t284 = t343 * t487;
t468 = t284 + t497;
t115 = (-t334 * t587 - t323) * qJD(1) - t468;
t234 = t252 * qJD(1);
t526 = -t115 - t234;
t525 = -t616 - t158;
t524 = -t332 * t159 - t163 * t539;
t523 = t332 * t160 + t164 * t539;
t519 = -t332 * t157 + t334 * t158;
t518 = t332 * t165 + t334 * t166;
t177 = rSges(5,1) * t533 - rSges(5,2) * t537 + t317;
t515 = -t158 - t177;
t509 = -t282 * t334 - t190;
t507 = -t280 * t334 + t192;
t486 = qJD(1) * t334;
t506 = t339 * t486 + t532;
t504 = -t242 + t245;
t503 = t244 + t411;
t502 = rSges(5,3) * t486 + qJD(1) * t269;
t500 = t262 - t328;
t498 = rSges(4,3) * t486 + t487 * t633;
t496 = -t280 + t283;
t495 = t282 + t412;
t489 = qJD(1) * t241;
t488 = qJD(1) * t279;
t482 = qJD(3) * t346;
t480 = qJD(1) * qJD(4);
t479 = pkin(3) * t530;
t478 = t349 * t593;
t477 = t349 * t338;
t474 = t340 * t581;
t391 = rSges(6,3) * t486 + qJD(1) * t257 - t334 * t474;
t471 = t326 * t239;
t94 = (-t327 * t487 - t471) * rSges(6,1) + t391;
t193 = t235 * t332;
t236 = -rSges(6,2) * t326 + t584;
t95 = -t340 * t193 + (t236 * t334 + t316) * qJD(1);
t476 = t165 * t486 + t332 * t95 + t334 * t94;
t473 = pkin(3) * t482;
t303 = pkin(6) * t486;
t462 = t334 * t483;
t432 = pkin(3) * t462;
t114 = -t432 - t303 + t304 + (t332 * t587 - t298) * qJD(1);
t472 = t115 * t485 + (t114 - t153) * t484;
t470 = t334 * t114 + t332 * t115 - t157 * t486;
t466 = t248 * t485;
t461 = -t157 * t485 + t158 * t484 + qJD(2);
t460 = -pkin(2) - t586;
t459 = t487 / 0.2e1;
t458 = t486 / 0.2e1;
t457 = -t485 / 0.2e1;
t454 = t484 / 0.2e1;
t250 = rSges(5,1) * t333 - rSges(5,2) * t331;
t450 = -t250 - t337;
t449 = t271 + t592;
t447 = -t262 - t584;
t385 = t232 * t340;
t445 = qJD(1) * t164 - t161 * t340 - t332 * t385;
t444 = -t162 * t340 - t334 * t385 + (-t233 * t332 + t566) * qJD(1);
t443 = qJD(1) * t162 + t163 * t340 - t230 * t238;
t442 = (-t332 * t410 + t560) * qJD(1) + t629 * t340;
t131 = t164 * t540;
t441 = t160 * t334 - t131;
t438 = -t159 + t555;
t435 = t626 * t340;
t434 = t627 * t340;
t429 = (-t332 ^ 2 - t334 ^ 2) * t592;
t428 = qJD(1) * (-pkin(2) * t487 + t303) - t478;
t227 = t250 * qJD(3);
t423 = -t227 - t473;
t422 = -t337 - t590;
t63 = -t466 + (t177 + t426) * qJD(1) - t497;
t420 = t63 * t451;
t249 = rSges(3,1) * t332 + rSges(3,2) * t334;
t417 = t586 - t633;
t416 = -t332 * t98 - t577;
t77 = t161 * t327 + t163 * t326;
t403 = t195 * t332 + t196 * t334;
t399 = qJD(1) * t276 + t334 * t480 - t477;
t397 = -t235 + t271;
t348 = qJD(3) ^ 2;
t394 = -qJD(3) * t227 - t337 * t348;
t393 = t422 * t348;
t392 = qJD(1) * t114 + t332 * t480 + t428;
t218 = t285 * t332;
t203 = t248 * t332;
t384 = t409 * t332;
t208 = t236 * t340;
t374 = -qJD(3) * t590 - t208 - t473;
t373 = qJD(1) * t229 - t238 * t547 + t239 * t548;
t372 = t172 * t334 - t173 * t332;
t371 = -t344 * t507 + t346 * t509;
t359 = -t326 * t442 + t327 * t444 + t492;
t10 = t332 * t613 + t359 * t334;
t360 = qJD(1) * t159 - t326 * t443 + t327 * t445;
t11 = t360 * t332 - t334 * t612;
t12 = t359 * t332 - t334 * t613;
t58 = -t384 - t556;
t59 = -t162 * t542 - t441;
t21 = t238 * t59 - t239 * t58 + t528;
t60 = -t161 * t541 - t524;
t61 = -t162 * t541 + t523;
t80 = -t334 * t402 + t548;
t76 = t80 * qJD(1);
t22 = t238 * t61 - t239 * t60 + t76;
t365 = (-t232 * t334 - t162) * t238 - (-t232 * t332 - t161) * t239 + t627 * qJD(1);
t350 = -t326 * t602 + t365 * t327;
t358 = qJD(1) * t228 - t326 * t435 + t327 * t434;
t36 = t332 * t607 + t358 * t334;
t37 = t358 * t332 - t334 * t607;
t38 = t326 * t445 + t327 * t443;
t39 = t326 * t444 + t327 * t442;
t78 = t162 * t327 + t164 * t326;
t9 = t332 * t612 + t360 * t334;
t370 = (qJD(1) * t36 + t10 * t238 + t222 * t60 + t223 * t61 - t239 * t9) * t595 + (t365 * t326 + t327 * t602) * t589 + t21 * t459 + (qJD(1) * t37 - t11 * t239 + t12 * t238 + t222 * t58 + t223 * t59) * t594 + t22 * t458 + (t10 * t332 - t334 * t9 + (t332 * t60 + t334 * t61) * qJD(1)) * t598 + (t332 * t59 - t334 * t58) * t601 + (t332 * t61 - t334 * t60) * t600 + (-t11 * t334 + t12 * t332 + (t332 * t58 + t334 * t59) * qJD(1)) * t597 + (t332 * t39 - t334 * t38 + (t332 * t77 + t334 * t78) * qJD(1)) * t588 + (t332 * t373 + t334 * t350) * t599 + (t332 * t350 - t334 * t373) * t596;
t369 = (-t331 * t503 + t333 * t504) * qJD(1);
t368 = (-t344 * t495 + t346 * t496) * qJD(1);
t194 = t235 * t334;
t33 = t165 * t238 + t166 * t239 + (t128 * t332 + t334 * t616) * qJD(3) + t461;
t367 = t33 * (-t238 * t193 - t194 * t239) + t49 * (-qJD(1) * t194 - t236 * t238);
t125 = -rSges(4,2) * t334 * t482 + (-t346 * t487 - t462) * rSges(4,1) + t498;
t126 = -qJD(3) * t218 + (t334 * t417 + t318) * qJD(1);
t366 = t125 * t334 + t126 * t332 + (t195 * t334 - t196 * t332) * qJD(1);
t351 = -t331 * t604 + t372 * t333;
t263 = t417 * qJD(3);
t204 = t248 * t334;
t180 = t449 * t334;
t179 = t449 * t332;
t169 = t239 * t236;
t110 = -qJD(3) * t203 + (t250 * t334 + t317) * qJD(1);
t109 = -t484 * t582 + (-t331 * t484 - t333 * t487) * rSges(5,1) + t502;
t96 = qJD(3) * t403 + qJD(2);
t75 = t246 * t332 + t276 + t284 + (t334 * t500 + t294) * qJD(1);
t74 = t432 + (-t332 * t500 + t298) * qJD(1) + t506;
t69 = -t477 - t263 * t484 + (-t126 - t234 + t465) * qJD(1);
t68 = -t263 * t485 + (t125 - t463) * qJD(1) + t428;
t53 = (t176 * t332 + t177 * t334) * qJD(3) + t461;
t50 = t366 * qJD(3);
t41 = t394 * t334 + (-t110 + t466 + t526) * qJD(1) + t399;
t40 = t394 * t332 + (t109 + t395) * qJD(1) + t392;
t24 = -t208 * t239 + t222 * t235 + t334 * t393 + (-t75 - t95 + t431 + t526) * qJD(1) + t399;
t23 = -t208 * t238 - t223 * t235 + t332 * t393 + (t74 + t94 + t532) * qJD(1) + t392;
t14 = (t109 * t334 + t110 * t332 + (t176 * t334 + t332 * t515) * qJD(1)) * qJD(3) + t472;
t5 = t165 * t223 - t166 * t222 + t238 * t95 + t239 * t94 + (t332 * t75 + t334 * t74 + (t128 * t334 + t332 * t525) * qJD(1)) * qJD(3) + t472;
t1 = [m(3) * ((-t249 * t349 - t478) * t615 + (-t477 + (-0.2e1 * t446 - t338 + t615) * t349) * (-t249 - t593)) + (t76 + (t59 + (t161 * t334 + t162 * t332) * t326 + t441 + t524) * t239 + (-t163 * t540 + t556 + t58 + (t161 * t332 - t162 * t334) * t326 + t523) * t238) * t596 + (t79 + t77) * t601 + (t78 + t80) * t600 + (-t528 + (t61 - t384 - t523) * t239 + (t438 * t332 - t131 + t60) * t238 + ((t160 + t409) * t238 + t438 * t239) * t334 + t21) * t599 + (t39 + t36) * t598 + ((t635 * t332 + ((t653 + t671) * t334 + t643 + t665 + t670) * t334) * qJD(3) + t655) * t454 + (t23 * (t166 + t338 + t628) - (t271 * t631 + t591 * t632) * qJD(3) + (-rSges(6,1) * t471 + t391 - t433 + t506 + t638) * t49 + (t447 * t332 - t501 + t531 - t593) * t24 + (t305 + (t340 * t585 - t246 + t474) * t332 - t421) * t48) * m(6) + (t225 * t333 + t226 * t331 + t260 * t346 + t261 * t344 + t326 * t434 + t327 * t435 - t667 * qJD(3) + ((-t345 * t49 - t347 * t48) * pkin(1) + (-t236 - t262) * t578 + (t48 * (-rSges(6,3) - t339) + t49 * t447) * t332 + t48 * t621 - t49 * (-t593 + t630)) * m(6)) * qJD(1) + (t41 * (-t617 + t659) + t40 * (t177 + t338 - t493) + (t332 * t576 + t334 * t420) * qJD(3) + (t468 + (-t317 - t338 + (-t250 - t328) * t334) * qJD(1)) * t62 + (t62 - t375 + t502 + (t176 + t593 + t659) * qJD(1) + t638) * t63) * m(5) + (-(-t463 - t237 - t97 + (-t195 - t593) * qJD(1)) * t98 + t69 * (t332 * t460 + t324 + t494 - t593) + t68 * t618 + t98 * (t303 + t498) + (t285 * t579 - t580) * qJD(3) + ((-t345 * t98 - t347 * t97) * pkin(1) + (-pkin(2) - t417) * t577 + (t97 * (-rSges(4,3) - pkin(6)) + t98 * t460) * t332) * qJD(1)) * m(4) + (t38 + t37 + t22) * t597 + (((t334 * t636 - t635 + t641) * t334 + (t332 * t636 + t642 + t666) * t332) * qJD(3) + t649 - t650) * t457 + (t645 + t646) * t485 / 0.2e1 - (t644 - t647 + t648) * t484 / 0.2e1 + ((t640 + t663) * t332 + (t639 + t662) * t334) * qJD(3) * t588; m(4) * t50 + m(5) * t14 + m(6) * t5; t370 + (((t332 * t507 - t334 * t508) * t346 + (t332 * t509 + t334 * t510) * t344 + t604 * t333 + t372 * t331) * qJD(3) + (t331 * t504 + t333 * t503 + t344 * t496 + t346 * t495) * qJD(1)) * t589 + (t647 * t334 + t646 * t332 + (t332 * t640 + t334 * t639) * qJD(1)) * t588 + ((-t485 * t545 + t489) * t332 + (t369 + (t332 * t546 + t351) * qJD(3)) * t334 + (-t485 * t543 + t488) * t332 + (t368 + (-t603 * t334 + (t544 + t371) * t332) * qJD(3)) * t334) * t457 + ((-t484 * t546 - t489) * t334 + (t369 + (t334 * t545 + t351) * qJD(3)) * t332 + (-t484 * t544 - t488) * t334 + (t368 + (t371 * t332 + (t543 - t603) * t334) * qJD(3)) * t332) * t454 + (t5 * (t518 + t519) + t33 * (t470 + t476) + (t24 * t397 + t48 * t374 + t5 * t616 + t33 * t74 + (t33 * t128 + t397 * t49) * qJD(1)) * t334 + (t23 * t397 + t49 * t374 + t5 * t128 + t33 * t75 + (t48 * (t235 + t591) + t33 * (-t166 + t525)) * qJD(1)) * t332 + t48 * t169 - (t48 * (-t179 + t193) + t49 * (t180 - t479)) * qJD(1) - (t33 * t429 + (t33 * t180 + t422 * t48) * t334 + (t33 * t179 + t422 * t49) * t332) * qJD(3) - t367) * m(6) + (-(t218 * t97 - t580) * qJD(1) - (t96 * (-t218 * t332 - t219 * t334) + t416 * t417) * qJD(3) + t50 * t403 + t96 * t366 + t416 * t263 + (-t68 * t332 - t69 * t334 + (-t334 * t98 + t579) * qJD(1)) * t285) * m(4) + (t645 * qJD(1) + ((t641 * qJD(1) + t624 * t334) * t334 + (t622 * t332 + t642 * qJD(1) + (-t623 + t625) * t334) * t332) * t620) * t595 + (t644 * qJD(1) + ((t643 * qJD(1) + t623 * t334) * t334 + (t625 * t332 + t664 * qJD(1) + (-t622 + t624) * t334) * t332) * t620) * t594 + (-(t62 * t203 + t63 * (-t204 - t479)) * qJD(1) - (t53 * t429 + (-t53 * t204 + t450 * t62) * t334 + (-t53 * t203 + t450 * t63) * t332) * qJD(3) + t14 * t519 + t53 * t470 + (t41 * t451 + t62 * t423 + t14 * t177 + t53 * t109 + (t53 * t176 + t420) * qJD(1)) * t334 + (t40 * t451 + t63 * t423 + t14 * t176 + t53 * t110 + (t515 * t53 + t576) * qJD(1)) * t332) * m(5) + (t649 + t651) * t459 + (t648 + t652) * t458; 0.2e1 * (t23 * t594 + t24 * t595) * m(6) + 0.2e1 * (t40 * t594 + t41 * t595) * m(5); t370 + (t5 * t518 + t33 * (-t166 * t487 + t476) + (-t332 * t49 - t578) * t208 + (-t23 * t332 - t24 * t334 + (-t631 + t632) * qJD(1)) * t235 - t48 * (qJD(1) * t193 - t169) - t367) * m(6);];
tauc = t1(:);
