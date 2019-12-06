% Calculate vector of centrifugal and Coriolis load on the joints for
% S5RRPRR5
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% qJD [5x1]
%   Generalized joint velocities
% pkin [9x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d1,d2,d4,d5,theta3]';
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
% Datum: 2019-12-05 18:34
% Revision: 77da58f92bca3eff71542919beafa37024070d86 (2019-12-05)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauc = S5RRPRR5_coriolisvecJ_fixb_slag_vp1(qJ, qJD, ...
  pkin, m, rSges, Icges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(9,1),zeros(6,1),zeros(6,3),zeros(6,6)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RRPRR5_coriolisvecJ_fixb_slag_vp1: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RRPRR5_coriolisvecJ_fixb_slag_vp1: qJD has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S5RRPRR5_coriolisvecJ_fixb_slag_vp1: pkin has to be [9x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5RRPRR5_coriolisvecJ_fixb_slag_vp1: m has to be [6x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [6,3]), ...
  'S5RRPRR5_coriolisvecJ_fixb_slag_vp1: rSges has to be [6x3] (double)');
assert(isreal(Icges) && all(size(Icges) == [6 6]), ...
  'S5RRPRR5_coriolisvecJ_fixb_slag_vp1: Icges has to be [6x6] (double)'); 

%% Symbolic Calculation
% From coriolisvec_joint_fixb_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-05 18:33:37
% EndTime: 2019-12-05 18:33:52
% DurationCPUTime: 9.62s
% Computational Cost: add. (18232->576), mult. (12208->730), div. (0->0), fcn. (9444->10), ass. (0->350)
t280 = pkin(9) + qJ(4);
t275 = qJ(5) + t280;
t270 = sin(t275);
t271 = cos(t275);
t209 = rSges(6,1) * t270 + rSges(6,2) * t271;
t283 = qJ(1) + qJ(2);
t276 = sin(t283);
t451 = t270 * t276;
t243 = Icges(6,4) * t451;
t277 = cos(t283);
t449 = t271 * t276;
t148 = -Icges(6,1) * t449 + Icges(6,5) * t277 + t243;
t450 = t270 * t277;
t244 = Icges(6,4) * t450;
t448 = t271 * t277;
t149 = Icges(6,1) * t448 + Icges(6,5) * t276 - t244;
t281 = qJD(4) + qJD(5);
t230 = t276 * t281;
t231 = t277 * t281;
t282 = qJD(1) + qJD(2);
t265 = Icges(6,4) * t271;
t346 = -Icges(6,2) * t270 + t265;
t535 = Icges(6,1) * t270 + t265;
t543 = t535 + t346;
t300 = t230 * (-Icges(6,2) * t448 + t149 - t244) + t231 * (Icges(6,2) * t449 + t148 + t243) + t282 * t543;
t146 = Icges(6,6) * t277 - t276 * t346;
t147 = Icges(6,4) * t448 - Icges(6,2) * t450 + Icges(6,6) * t276;
t474 = Icges(6,4) * t270;
t205 = Icges(6,2) * t271 + t474;
t208 = Icges(6,1) * t271 - t474;
t511 = t230 * (t277 * t535 + t147) + t231 * (-t276 * t535 + t146) + t282 * (t205 - t208);
t551 = t300 * t270 + t511 * t271;
t273 = sin(t280);
t274 = cos(t280);
t447 = t273 * t276;
t255 = Icges(5,4) * t447;
t445 = t274 * t276;
t159 = -Icges(5,1) * t445 + Icges(5,5) * t277 + t255;
t446 = t273 * t277;
t256 = Icges(5,4) * t446;
t444 = t274 * t277;
t160 = Icges(5,1) * t444 + Icges(5,5) * t276 - t256;
t307 = t276 * (-Icges(5,2) * t444 + t160 - t256) + t277 * (Icges(5,2) * t445 + t159 + t255);
t266 = Icges(5,4) * t274;
t347 = -Icges(5,2) * t273 + t266;
t157 = Icges(5,6) * t277 - t276 * t347;
t158 = Icges(5,4) * t444 - Icges(5,2) * t446 + Icges(5,6) * t276;
t515 = t277 * t157 + t158 * t276;
t550 = -t307 * t273 - t515 * t274;
t549 = 2 * qJD(4);
t173 = t209 * t276;
t455 = t209 * t277;
t441 = t276 * t282;
t491 = rSges(6,1) * t271;
t544 = t209 * t281;
t393 = -t544 * t277 - t441 * t491;
t481 = t277 * rSges(6,3);
t97 = (rSges(6,2) * t451 + t481) * t282 + t393;
t548 = t173 * t230 + t231 * t455 + t277 * t97;
t534 = Icges(5,1) * t273 + t266;
t415 = t534 + t347;
t475 = Icges(5,4) * t273;
t224 = Icges(5,2) * t274 + t475;
t227 = Icges(5,1) * t274 - t475;
t416 = t224 - t227;
t547 = (t273 * t415 + t274 * t416) * t282;
t145 = Icges(6,5) * t448 - Icges(6,6) * t450 + Icges(6,3) * t276;
t546 = t277 * t145 - t149 * t449;
t353 = rSges(3,1) * t276 + rSges(3,2) * t277;
t193 = t353 * t282;
t286 = sin(qJ(1));
t485 = pkin(1) * qJD(1);
t405 = t286 * t485;
t175 = t193 + t405;
t488 = rSges(6,2) * t270;
t210 = -t488 + t491;
t545 = t282 * t173 - t209 * t441 - t210 * t231;
t180 = t210 * t281;
t439 = t277 * t282;
t285 = cos(pkin(9));
t390 = pkin(3) * t285 + pkin(2);
t498 = pkin(4) * t274;
t236 = t390 + t498;
t195 = t277 * t236;
t495 = pkin(7) + qJ(3);
t264 = t276 * t495;
t279 = -pkin(8) - t495;
t121 = t276 * t279 + t277 * t390 - t195 + t264;
t247 = rSges(6,2) * t450;
t333 = -rSges(6,1) * t448 - rSges(6,3) * t276;
t150 = -t247 - t333;
t411 = qJD(4) * t276;
t388 = t273 * t411;
t245 = pkin(4) * t388;
t268 = qJD(3) * t277;
t414 = t245 + t268;
t357 = t209 * t230 + t414;
t287 = cos(qJ(1));
t403 = t287 * t485;
t438 = t277 * t285;
t443 = t276 * qJ(3);
t332 = -pkin(3) * t438 + t443;
t143 = -t264 + t332;
t232 = t277 * pkin(2) + t443;
t432 = t143 - t232;
t59 = -t403 + (t121 - t150 + t432) * t282 + t357;
t542 = (t180 * t276 + t209 * t439 - t230 * t210 - t282 * t455) * t59;
t437 = t279 * t282;
t248 = t276 * t437;
t381 = t282 * t495;
t251 = t276 * t381;
t363 = -t236 + t390;
t101 = t363 * t439 + t245 + t248 + t251;
t394 = t282 * t247 + t544 * t276;
t98 = t282 * t333 + t394;
t541 = -t101 - t98;
t343 = t158 * t273 - t160 * t274;
t538 = t343 * t277;
t490 = rSges(4,2) * sin(pkin(9));
t263 = t277 * t490;
t335 = -rSges(4,1) * t438 - rSges(4,3) * t276;
t172 = -t263 - t335;
t235 = t282 * t263;
t537 = -t282 * t172 - t235;
t440 = t276 * t285;
t234 = t282 * rSges(4,1) * t440;
t260 = pkin(2) * t441;
t536 = t234 + t260;
t198 = t282 * t232;
t412 = t268 - t198;
t303 = (-t279 - t495) * t277 + t363 * t276;
t118 = t282 * t303;
t119 = t282 * t121;
t533 = -t277 * t118 + t276 * t119;
t135 = t282 * t143;
t137 = t282 * t150;
t532 = t119 + t135 - t137 - t198 + t357 - t248 - t394 - t414;
t204 = Icges(6,5) * t271 - Icges(6,6) * t270;
t320 = t204 * t276;
t144 = Icges(6,3) * t277 - t320;
t435 = -t276 * t144 - t148 * t448;
t531 = t435 - t546;
t489 = rSges(5,2) * t273;
t493 = rSges(5,1) * t274;
t351 = -t489 + t493;
t202 = t351 * qJD(4);
t228 = rSges(5,1) * t273 + rSges(5,2) * t274;
t530 = -t202 * t277 + t228 * t441;
t456 = t205 * t281;
t528 = -Icges(6,6) * t282 + t456;
t203 = Icges(6,5) * t270 + Icges(6,6) * t271;
t527 = -Icges(6,3) * t282 + t203 * t281;
t526 = t210 * t276 - t481;
t482 = t277 * rSges(5,3);
t525 = t276 * t351 - t482;
t524 = -Icges(6,5) * t282 + t281 * t535;
t103 = t157 * t274 + t159 * t273;
t323 = t347 * t282;
t520 = -Icges(5,6) * t282 + qJD(4) * t224;
t111 = t520 * t276 - t277 * t323;
t325 = t227 * t282;
t518 = -Icges(5,5) * t282 + qJD(4) * t534;
t113 = t518 * t276 - t277 * t325;
t223 = Icges(5,5) * t274 - Icges(5,6) * t273;
t155 = Icges(5,3) * t277 - t223 * t276;
t523 = qJD(4) * t103 + t111 * t273 - t113 * t274 - t155 * t282;
t104 = t158 * t274 + t160 * t273;
t110 = -t276 * t323 - t520 * t277;
t112 = -t276 * t325 - t518 * t277;
t156 = Icges(5,5) * t444 - Icges(5,6) * t446 + Icges(5,3) * t276;
t522 = qJD(4) * t104 + t110 * t273 - t112 * t274 - t156 * t282;
t222 = Icges(5,5) * t273 + Icges(5,6) * t274;
t521 = -Icges(5,3) * t282 + qJD(4) * t222;
t200 = t347 * qJD(4);
t201 = t227 * qJD(4);
t519 = t200 * t273 - t201 * t274 - t222 * t282 + (t224 * t274 + t273 * t534) * qJD(4);
t257 = rSges(5,2) * t446;
t334 = -rSges(5,1) * t444 - t276 * rSges(5,3);
t164 = -t257 - t334;
t192 = t228 * t411;
t391 = rSges(5,1) * t388 + (t273 * t439 + t274 * t411) * rSges(5,2);
t517 = -t282 * t164 + t135 + t192 + t251 - t268 - t391;
t366 = -t208 * t281 + t456;
t367 = t543 * t281;
t514 = -t203 * t282 + t270 * t367 + t271 * t366;
t322 = t346 * t282;
t372 = t149 * t281 - t276 * t322 - t528 * t277;
t324 = t208 * t282;
t374 = t147 * t281 + t276 * t324 + t524 * t277;
t513 = -t145 * t282 + t270 * t372 + t271 * t374;
t373 = t148 * t281 + t528 * t276 - t277 * t322;
t375 = t146 * t281 - t524 * t276 + t277 * t324;
t512 = -t144 * t282 + t270 * t373 + t271 * t375;
t510 = t277 ^ 2;
t190 = t281 * t441;
t509 = -t190 / 0.2e1;
t191 = t282 * t231;
t508 = t191 / 0.2e1;
t507 = -t230 / 0.2e1;
t506 = t230 / 0.2e1;
t505 = -t231 / 0.2e1;
t504 = t231 / 0.2e1;
t503 = t276 / 0.2e1;
t502 = t277 / 0.2e1;
t501 = -t282 / 0.2e1;
t500 = t282 / 0.2e1;
t497 = t286 * pkin(1);
t496 = t287 * pkin(1);
t494 = rSges(4,1) * t285;
t486 = rSges(4,3) * t277;
t484 = pkin(4) * qJD(4);
t359 = -t268 + t403;
t68 = t192 + (-t164 + t432) * t282 - t359;
t483 = t276 * t68;
t166 = t203 * t276;
t342 = t205 * t270 - t271 * t535;
t88 = -t277 * t342 + t166;
t480 = t88 * t282;
t479 = rSges(4,3) + qJ(3);
t478 = rSges(6,3) - t279;
t181 = t222 * t276;
t340 = t224 * t273 - t274 * t534;
t106 = -t277 * t340 + t181;
t465 = t106 * t282;
t462 = t157 * t273;
t461 = t159 * t274;
t458 = t203 * t277;
t454 = t222 * t277;
t453 = t223 * t282;
t187 = t228 * t276;
t452 = t228 * t277;
t442 = t276 * t145;
t436 = t277 * t144 + t146 * t451;
t434 = t277 * t155 + t157 * t447;
t433 = t276 * t155 + t159 * t444;
t423 = t236 * t441 + t277 * t437;
t422 = -t195 + t247;
t196 = t282 * (-t276 * pkin(2) + t277 * qJ(3));
t267 = qJD(3) * t276;
t421 = -t196 - t267;
t418 = -t277 * t381 + t390 * t441;
t413 = t257 - t264;
t410 = qJD(4) * t277;
t409 = (qJD(4) ^ 2) * t498;
t289 = qJD(1) ^ 2;
t408 = t289 * t496;
t404 = t273 * t484;
t402 = t146 * t450;
t395 = -t282 * (-pkin(3) * t440 + t277 * pkin(7)) + t421;
t385 = t274 * t410;
t387 = t273 * t410;
t392 = -rSges(5,1) * t387 - rSges(5,2) * t385 - t441 * t493;
t389 = t228 * t410;
t384 = -t441 / 0.2e1;
t383 = t439 / 0.2e1;
t382 = -pkin(2) - t494;
t380 = -t411 / 0.2e1;
t378 = -t410 / 0.2e1;
t377 = t410 / 0.2e1;
t376 = pkin(4) * t273 + t209;
t233 = rSges(3,1) * t277 - t276 * rSges(3,2);
t369 = -t156 - t461;
t358 = qJ(3) * t439 - t260;
t365 = -0.2e1 * t267 - t358;
t364 = pkin(4) * t387;
t194 = -rSges(3,1) * t439 + rSges(3,2) * t441;
t360 = -t267 + t405;
t356 = t358 + t418 + t365;
t352 = t490 - t494;
t90 = t147 * t271 + t149 * t270;
t345 = t147 * t270 - t149 * t271;
t344 = -t461 + t462;
t339 = -t408 + (t268 + t412) * t282;
t338 = -t390 - t493;
t337 = t198 + t359;
t71 = t277 * t156 + t158 * t447 - t160 * t445;
t331 = -t267 - t393 + t423;
t329 = -t267 - t392 + t418;
t176 = -t233 * t282 - t403;
t328 = t282 * (t282 * t332 - t251) + t339;
t70 = -t159 * t445 + t434;
t327 = (t276 * t71 + t277 * t70) * qJD(4);
t72 = -t157 * t446 + t433;
t73 = t156 * t276 - t538;
t326 = (t276 * t73 + t277 * t72) * qJD(4);
t317 = t166 * t231 + t204 * t282 - t230 * t458;
t315 = t282 * t525 + t389 + t395;
t314 = -t527 * t277 + (-t320 + t345) * t282;
t313 = -t204 * t439 + t527 * t276 + (t146 * t270 - t148 * t271) * t282;
t312 = -t453 * t276 - t521 * t277 + t282 * t343;
t311 = t521 * t276 - t453 * t277 + t282 * t344;
t310 = t204 * t281 + t282 * t342;
t309 = t223 * qJD(4) + t282 * t340;
t13 = t313 * t276 - t512 * t277;
t14 = t314 * t276 - t513 * t277;
t15 = t512 * t276 + t313 * t277;
t16 = t513 * t276 + t314 * t277;
t63 = -t148 * t449 + t436;
t124 = t147 * t451;
t64 = t124 + t546;
t87 = t276 * t342 + t458;
t83 = t87 * t282;
t28 = t230 * t64 + t231 * t63 + t83;
t65 = -t402 - t435;
t66 = -t277 * t345 + t442;
t29 = t230 * t66 + t231 * t65 + t480;
t42 = t310 * t276 - t514 * t277;
t43 = t514 * t276 + t310 * t277;
t44 = -t270 * t375 + t271 * t373;
t45 = -t270 * t374 + t271 * t372;
t89 = t146 * t271 + t148 * t270;
t304 = (t13 * t231 + t14 * t230 - t190 * t65 + t191 * t66 + t282 * t42) * t503 + (t317 * t276 - t551 * t277) * t507 + (t551 * t276 + t317 * t277) * t505 + (t15 * t231 + t16 * t230 - t190 * t63 + t191 * t64 + t282 * t43) * t502 + (-t270 * t511 + t271 * t300) * t501 + t28 * t384 + t29 * t383 + ((t282 * t66 + t13) * t277 + (-t282 * t65 + t14) * t276) * t506 + (t276 * t64 + t277 * t63) * t509 + (t276 * t66 + t277 * t65) * t508 + ((t282 * t64 + t15) * t277 + (-t282 * t63 + t16) * t276) * t504 + ((t282 * t90 + t44) * t277 + (-t282 * t89 + t45) * t276) * t500;
t136 = t282 * t526;
t302 = t209 * t231 - t118 + t136 + t364 + t395;
t299 = t276 * t303;
t295 = t147 * t450 - t442 + (-t148 * t276 - t149 * t277) * t271 + t436;
t102 = (t277 * t164 + t525 * t276) * qJD(4);
t105 = t276 * t340 + t454;
t99 = t105 * t282;
t34 = t99 + t327;
t35 = t326 + t465;
t48 = -qJD(4) * t344 + t111 * t274 + t113 * t273;
t49 = -qJD(4) * t343 + t110 * t274 + t112 * t273;
t54 = t309 * t276 - t519 * t277;
t55 = t519 * t276 + t309 * t277;
t293 = (t83 + (t66 + t295) * t231 + (t124 - t65 - t402 - t531) * t230) * t507 + (t99 + ((t434 + t73 + t538) * t277 + (-t72 + (t369 - t462) * t277 + t71 + t433) * t276) * qJD(4)) * t380 + (t89 + t87) * t509 + (t90 + t88) * t508 + (t29 - t480 + (t64 + (t146 * t277 - t147 * t276) * t270 + t531) * t231 + (-t63 + t295) * t230) * t505 + (t44 + t43) * t504 + (t35 - t465 + ((t71 + (-t156 + t462) * t277 - t433) * t277 + (t276 * t369 + t434 - t70) * t276) * qJD(4)) * t378 + (t48 + t55) * t377 + (t45 + t42 + t28) * t506 + (t49 + t54 + t34) * t411 / 0.2e1 + ((t103 + t105) * t380 + (t104 + t106) * t377 - qJD(4) * t340 + t200 * t274 + t201 * t273 - t270 * t366 + t271 * t367) * t282;
t100 = -t364 + t418 - t423;
t272 = t289 * t497;
t30 = t276 * t409 + t180 * t230 + t191 * t209 + t272 + (-t100 + t356 - t97 + t364) * t282;
t31 = -t277 * t409 - t180 * t231 + t190 * t209 + (t245 - t541) * t282 + t328;
t58 = t302 + t405;
t292 = (-t30 * t478 + t31 * (-t210 - t236) + (rSges(6,3) * t58 - t488 * t59) * t282) * t276 + (-t30 * t491 + t59 * t404 + t31 * t478 + (-t59 * rSges(6,3) - t58 * (-t236 - t491)) * t282) * t277;
t114 = (rSges(5,2) * t447 + t482) * t282 + t392;
t52 = t202 * t411 + t272 + (-t114 + t356 + t389) * t282;
t115 = t282 * t334 + t391;
t53 = t530 * qJD(4) + t115 * t282 + t328;
t67 = t315 + t405;
t291 = (t52 * t338 + t53 * (rSges(5,3) + t495) + (-t68 * rSges(5,3) - t338 * t67) * t282) * t277 + (-t68 * t282 * t489 + t53 * (t338 + t489) + (t282 * t67 - t52) * rSges(5,3)) * t276;
t161 = t282 * (t276 * t352 + t486);
t116 = -t161 - t196 + t360;
t117 = (-t172 - t232) * t282 - t359;
t79 = t272 + (t234 + (-t276 * t490 - t486) * t282 + t365) * t282;
t80 = (t282 * t335 + t235) * t282 + t339;
t290 = (-t79 * t479 + t80 * (-pkin(2) + t352)) * t276 + (t79 * t382 + t80 * t479) * t277 + ((t116 * t479 - t117 * t490) * t276 + (-t116 * t382 - t117 * t479) * t277) * t282;
t154 = t194 * t282 - t408;
t153 = t193 * t282 + t272;
t133 = t277 * t150;
t57 = -qJD(4) * t299 - t121 * t410 + t231 * t150 + t230 * t526;
t12 = t533 * qJD(4) + t100 * t410 - t101 * t411 - t190 * t150 + t191 * t526 - t230 * t98 + t231 * t97;
t1 = [m(3) * (t153 * (-t233 - t496) + t154 * (-t353 - t497) + (t176 - t194 + t403) * t175) + t293 + (t30 * (t422 - t496) + t59 * (t331 + t405) - t31 * t497 + t292 + (-t59 + t532) * t58) * m(6) + (t52 * (t413 - t496) + t68 * (t329 + t405) - t53 * t497 + t291 + (-t337 - t68 + t403 + t517) * t67) * m(5) + (t79 * (t263 - t496) + t117 * (t360 + t536) - t80 * t497 + t290 + (-t117 - t337 + t359 + t537) * t116) * m(4); t293 + (t30 * t422 + t292 + (-t302 + t331) * t59 + t532 * t58) * m(6) + (t52 * t413 + t291 + (-t315 + t329) * t68 + (t412 + t517) * t67) * m(5) + (t79 * t263 + t290 + (t161 - t421 - t267 + t536) * t117 + (t412 - t268 + t537) * t116) * m(4) + (-t153 * t233 - t154 * t353 - t175 * t194 + t176 * t193 - (t175 * t233 + t176 * t353) * t282) * m(3); m(4) * (t276 * t80 + t79 * t277) + m(5) * (t276 * t53 + t52 * t277) + m(6) * (t276 * t31 + t30 * t277); ((t104 * t282 + t48) * t277 + (-t103 * t282 + t49) * t276) * t500 + ((-t411 * t454 + t453) * t276 + (-t547 + (t276 * t181 + t550) * qJD(4)) * t277) * t380 + ((t181 * t410 + t453) * t277 + (t547 + (-t277 * t454 - t550) * qJD(4)) * t276) * t378 + ((-t273 * t416 + t274 * t415) * t282 + (-t273 * t515 + t274 * t307) * qJD(4)) * t501 + t304 + (t282 * t54 + ((t311 * t276 - t523 * t277 + t282 * t73) * t277 + (t312 * t276 - t522 * t277 - t282 * t72) * t276) * t549) * t503 + (t282 * t55 + ((t523 * t276 + t311 * t277 + t282 * t71) * t277 + (t522 * t276 + t312 * t277 - t282 * t70) * t276) * t549) * t502 + (t327 + t34) * t384 + (t326 + t35) * t383 + (t12 * (t133 - t299) + (-t12 * t121 - t31 * t376) * t277 + (t12 * t526 + t30 * t376) * t276 + (-(-t276 ^ 2 - t510) * t404 + t533 + (t100 + t136) * t277 + (-t137 + t541) * t276 + t548) * t57 + t542 + (-pkin(4) * t385 + (t274 * t484 + t180) * t277 + t545) * t58) * m(6) + (-(-t187 * t67 + t452 * t68) * t282 - (t102 * (-t187 * t276 - t277 * t452) + (t277 * t67 + t483) * t351) * qJD(4) + t228 * t439 * t68 + 0.2e1 * t102 * ((-rSges(5,3) * t439 + t114) * t277 + (-t115 + (t277 * t351 - t164) * t282) * t276) + t187 * t52 + t202 * t483 - t452 * t53 - t530 * t67) * m(5); t304 + (t12 * (t526 * t276 + t133) + t30 * t173 - t31 * t455 + (-t510 * t282 * rSges(6,3) + (-t98 + (t210 * t277 - t150) * t282) * t276 + t548) * t57 + t542 + (t277 * t180 + t545) * t58) * m(6);];
tauc = t1(:);
