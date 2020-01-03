% Calculate vector of centrifugal and Coriolis load on the joints for
% S5RRRRR5
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% qJD [5x1]
%   Generalized joint velocities
% pkin [9x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d1,d2,d3,d4,d5]';
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
% Datum: 2020-01-03 12:13
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauc = S5RRRRR5_coriolisvecJ_fixb_slag_vp1(qJ, qJD, ...
  pkin, m, rSges, Icges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(9,1),zeros(6,1),zeros(6,3),zeros(6,6)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RRRRR5_coriolisvecJ_fixb_slag_vp1: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RRRRR5_coriolisvecJ_fixb_slag_vp1: qJD has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S5RRRRR5_coriolisvecJ_fixb_slag_vp1: pkin has to be [9x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5RRRRR5_coriolisvecJ_fixb_slag_vp1: m has to be [6x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [6,3]), ...
  'S5RRRRR5_coriolisvecJ_fixb_slag_vp1: rSges has to be [6x3] (double)');
assert(isreal(Icges) && all(size(Icges) == [6 6]), ...
  'S5RRRRR5_coriolisvecJ_fixb_slag_vp1: Icges has to be [6x6] (double)'); 

%% Symbolic Calculation
% From coriolisvec_joint_fixb_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2020-01-03 12:13:03
% EndTime: 2020-01-03 12:13:18
% DurationCPUTime: 9.60s
% Computational Cost: add. (22097->569), mult. (13212->745), div. (0->0), fcn. (10256->10), ass. (0->351)
t276 = sin(qJ(1));
t471 = pkin(1) * qJD(1);
t386 = t276 * t471;
t274 = qJ(1) + qJ(2);
t262 = sin(t274);
t272 = qJD(1) + qJD(2);
t427 = t262 * t272;
t394 = pkin(2) * t427;
t312 = t386 + t394;
t260 = qJD(3) + t272;
t266 = qJ(3) + t274;
t255 = sin(t266);
t256 = cos(t266);
t513 = t255 * rSges(4,1) + t256 * rSges(4,2);
t518 = t513 * t260;
t125 = t312 + t518;
t277 = cos(qJ(4));
t265 = Icges(5,4) * t277;
t275 = sin(qJ(4));
t336 = -Icges(5,2) * t275 + t265;
t316 = t336 * t260;
t465 = Icges(5,4) * t275;
t234 = Icges(5,2) * t277 + t465;
t500 = -Icges(5,6) * t260 + qJD(4) * t234;
t100 = -t255 * t500 + t256 * t316;
t237 = Icges(5,1) * t277 - t465;
t318 = t237 * t260;
t510 = Icges(5,1) * t275 + t265;
t497 = -Icges(5,5) * t260 + qJD(4) * t510;
t102 = -t255 * t497 + t256 * t318;
t432 = t256 * t277;
t433 = t256 * t275;
t147 = Icges(5,4) * t432 - Icges(5,2) * t433 + Icges(5,6) * t255;
t229 = Icges(5,4) * t433;
t149 = Icges(5,1) * t432 + Icges(5,5) * t255 - t229;
t332 = t147 * t277 + t149 * t275;
t148 = -Icges(5,5) * t256 + t237 * t255;
t447 = t148 * t277;
t146 = -Icges(5,6) * t256 + t255 * t336;
t449 = t146 * t275;
t333 = t447 - t449;
t526 = -qJD(4) * t333 - t100 * t277 - t102 * t275 + t260 * t332;
t264 = cos(t274);
t511 = t262 * rSges(3,1) + t264 * rSges(3,2);
t182 = t511 * t272;
t166 = t386 + t182;
t241 = rSges(5,1) * t275 + rSges(5,2) * t277;
t322 = qJD(4) * t241;
t475 = rSges(5,1) * t277;
t367 = pkin(3) + t475;
t430 = t260 * t275;
t389 = rSges(5,2) * t430;
t231 = rSges(5,2) * t433;
t392 = rSges(5,1) * t432;
t151 = rSges(5,3) * t255 - t231 + t392;
t478 = pkin(3) * t256;
t193 = pkin(8) * t255 + t478;
t399 = qJD(4) * t255;
t375 = t241 * t399;
t425 = t264 * t272;
t240 = pkin(2) * t425;
t278 = cos(qJ(1));
t259 = t278 * t471;
t400 = t240 + t259;
t75 = -t375 + (t151 + t193) * t260 + t400;
t469 = t260 * t75;
t473 = rSges(5,2) * t275;
t480 = rSges(5,3) + pkin(8);
t397 = qJD(4) * t275;
t371 = t256 * t397;
t396 = qJD(4) * t277;
t436 = t256 * t260;
t412 = rSges(5,3) * t436 + t255 * t389;
t437 = t255 * t277;
t104 = rSges(5,2) * t256 * t396 + (t260 * t437 + t371) * rSges(5,1) - t412;
t217 = pkin(8) * t436;
t441 = t255 * t260;
t156 = pkin(3) * t441 - t217;
t341 = -t473 + t475;
t214 = t341 * qJD(4);
t268 = t276 * pkin(1);
t281 = qJD(1) ^ 2;
t395 = t281 * t268;
t479 = pkin(2) * t272 ^ 2;
t319 = -t262 * t479 - t395;
t398 = qJD(4) * t256;
t374 = t241 * t398;
t55 = -t214 * t399 + (-t104 - t156 - t374) * t260 + t319;
t372 = t255 * t396;
t373 = t255 * t397;
t411 = rSges(5,3) * t441 + t260 * t392;
t105 = -rSges(5,1) * t373 + (-t256 * t430 - t372) * rSges(5,2) + t411;
t269 = t278 * pkin(1);
t258 = t281 * t269;
t401 = t264 * t479 + t258;
t408 = pkin(3) * t436 + pkin(8) * t441;
t382 = t260 * t408 + t401;
t56 = t105 * t260 + (t214 * t256 - t241 * t441) * qJD(4) + t382;
t230 = rSges(5,1) * t437;
t438 = t255 * t275;
t150 = -rSges(5,2) * t438 - rSges(5,3) * t256 + t230;
t249 = t255 * pkin(3);
t192 = -pkin(8) * t256 + t249;
t305 = -t260 * (t150 + t192) - t374;
t74 = -t305 + t312;
t283 = (-t322 * t74 - t367 * t469 - t473 * t56 + t480 * t55) * t255 + (-t322 * t75 + t367 * t55 - t389 * t74 - t480 * t56) * t256;
t380 = t217 + t412;
t129 = t260 * t151;
t177 = t260 * t193;
t508 = -t129 - t177 + t375 + t408 + t411;
t525 = t283 + (-t305 + t380) * t75 + t508 * t74;
t477 = pkin(4) * t277;
t257 = pkin(3) + t477;
t279 = -pkin(9) - pkin(8);
t385 = pkin(4) * t397;
t273 = qJ(4) + qJ(5);
t263 = cos(t273);
t271 = qJD(4) + qJD(5);
t426 = t263 * t271;
t387 = rSges(6,2) * t426;
t261 = sin(t273);
t428 = t261 * t271;
t390 = rSges(6,1) * t428;
t292 = -t385 - t387 - t390;
t189 = t256 * t271;
t153 = t260 * t189;
t472 = rSges(6,2) * t261;
t474 = rSges(6,1) * t263;
t204 = -t472 + t474;
t181 = t204 * t271;
t188 = t255 * t271;
t202 = rSges(6,1) * t261 + rSges(6,2) * t263;
t350 = pkin(4) * t371;
t393 = qJD(4) ^ 2 * t477;
t431 = t260 * t261;
t388 = rSges(6,2) * t431;
t414 = rSges(6,3) * t436 + t255 * t388;
t439 = t255 * t263;
t85 = t256 * t387 + (t256 * t428 + t260 * t439) * rSges(6,1) - t414;
t245 = t256 * t279;
t95 = t350 + t217 + (t245 + (-pkin(3) + t257) * t255) * t260;
t476 = -t85 - t95;
t32 = -t255 * t393 - t153 * t202 - t181 * t188 + (-t156 - t350 + t476) * t260 + t319;
t152 = t260 * t188;
t351 = pkin(4) * t373;
t434 = t256 * t263;
t391 = rSges(6,1) * t434;
t413 = rSges(6,3) * t441 + t260 * t391;
t86 = -t255 * t390 + (-t255 * t426 - t256 * t431) * rSges(6,2) + t413;
t185 = t257 * t436;
t429 = t260 * t279;
t96 = t185 + (-t385 - t429) * t255 - t408;
t33 = t256 * t393 - t152 * t202 + t181 * t189 + (t86 + t96 - t351) * t260 + t382;
t326 = -t189 * t202 - t350;
t296 = -t326 + t394;
t407 = t255 * t257 + t245;
t127 = -t192 + t407;
t224 = rSges(6,1) * t439;
t440 = t255 * t261;
t136 = -rSges(6,2) * t440 - rSges(6,3) * t256 + t224;
t383 = t127 + t136 + t192;
t57 = t260 * t383 + t296 + t386;
t310 = -t188 * t202 - t351;
t223 = t256 * t257;
t128 = t478 - t223 + (pkin(8) + t279) * t255;
t435 = t256 * t261;
t225 = rSges(6,2) * t435;
t137 = rSges(6,3) * t255 - t225 + t391;
t424 = -t128 + t137;
t58 = (t193 + t424) * t260 + t310 + t400;
t282 = (t32 * (rSges(6,3) - t279) - t33 * t472 + t57 * t292 + (t58 * (-t257 - t474) - t57 * t279) * t260) * t255 + (t32 * t474 + t58 * (t292 - t429) - t33 * rSges(6,3) - t57 * t388) * t256;
t470 = t260 * t58;
t118 = t260 * t128;
t124 = t260 * t137;
t507 = t118 - t124 - t177 - t310 + t185 + t413;
t524 = t383 * t470 + t507 * t57 + t282;
t464 = Icges(6,4) * t261;
t201 = Icges(6,1) * t263 - t464;
t134 = -Icges(6,5) * t256 + t201 * t255;
t221 = Icges(6,4) * t435;
t135 = Icges(6,1) * t434 + Icges(6,5) * t255 - t221;
t198 = Icges(6,2) * t263 + t464;
t250 = Icges(6,4) * t263;
t335 = -Icges(6,2) * t261 + t250;
t512 = Icges(6,1) * t261 + t250;
t517 = t512 + t335;
t287 = t188 * (-Icges(6,2) * t434 + t135 - t221) - t189 * (-t198 * t255 + t134) + t260 * t517;
t315 = t335 * t255;
t132 = -Icges(6,6) * t256 + t315;
t133 = Icges(6,4) * t434 - Icges(6,2) * t435 + Icges(6,6) * t255;
t491 = t188 * (t256 * t512 + t133) - t189 * (t255 * t512 + t132) + t260 * (t198 - t201);
t523 = t287 * t261 + t263 * t491;
t402 = t510 + t336;
t403 = t234 - t237;
t519 = (t275 * t402 + t277 * t403) * t260;
t516 = 0.2e1 * qJD(4);
t197 = Icges(6,5) * t263 - Icges(6,6) * t261;
t313 = t197 * t255;
t130 = -Icges(6,3) * t256 + t313;
t451 = t133 * t261;
t334 = -t135 * t263 + t451;
t323 = -t130 + t334;
t515 = t189 * t323;
t191 = t256 * rSges(4,1) - rSges(4,2) * t255;
t174 = t260 * t191;
t126 = t400 + t174;
t443 = t198 * t271;
t506 = -Icges(6,6) * t260 + t443;
t163 = t202 * t255;
t164 = t202 * t256;
t50 = t136 * t188 + t137 * t189 + (t127 * t255 - t128 * t256) * qJD(4);
t505 = -t57 * (-t260 * t163 + t189 * t204) - t50 * (-t163 * t188 - t189 * t164) - t58 * (-t164 * t260 - t188 * t204);
t196 = Icges(6,5) * t261 + Icges(6,6) * t263;
t504 = -Icges(6,3) * t260 + t196 * t271;
t503 = -Icges(6,5) * t260 + t271 * t512;
t106 = t277 * t146 + t148 * t275;
t233 = Icges(5,5) * t277 - Icges(5,6) * t275;
t144 = -Icges(5,3) * t256 + t233 * t255;
t502 = qJD(4) * t106 + t100 * t275 - t102 * t277 - t144 * t260;
t232 = Icges(5,5) * t275 + Icges(5,6) * t277;
t501 = -Icges(5,3) * t260 + qJD(4) * t232;
t207 = t336 * qJD(4);
t208 = t237 * qJD(4);
t499 = qJD(4) * (t234 * t277 + t275 * t510) + t207 * t275 - t208 * t277 - t232 * t260;
t101 = t255 * t318 + t256 * t497;
t145 = Icges(5,5) * t432 - Icges(5,6) * t433 + Icges(5,3) * t255;
t99 = t255 * t316 + t256 * t500;
t498 = qJD(4) * t332 + t101 * t277 - t145 * t260 - t275 * t99;
t352 = -t201 * t271 + t443;
t353 = t517 * t271;
t494 = -t196 * t260 + t261 * t353 + t263 * t352;
t131 = Icges(6,5) * t434 - Icges(6,6) * t435 + Icges(6,3) * t255;
t357 = t135 * t271 - t256 * t506 - t260 * t315;
t317 = t201 * t260;
t359 = t133 * t271 + t255 * t317 + t256 * t503;
t493 = -t131 * t260 + t261 * t357 + t263 * t359;
t358 = t134 * t271 - t255 * t506 + t335 * t436;
t360 = t132 * t271 + t255 * t503 - t256 * t317;
t492 = -t130 * t260 + t261 * t358 + t263 * t360;
t490 = t152 / 0.2e1;
t489 = -t153 / 0.2e1;
t488 = -t188 / 0.2e1;
t487 = t188 / 0.2e1;
t486 = t189 / 0.2e1;
t485 = -t189 / 0.2e1;
t484 = -t255 / 0.2e1;
t483 = -t256 / 0.2e1;
t482 = -t260 / 0.2e1;
t481 = t260 / 0.2e1;
t445 = t196 * t255;
t92 = -t198 * t435 + t434 * t512 + t445;
t468 = t92 * t260;
t442 = t232 * t255;
t110 = -t234 * t433 + t432 * t510 + t442;
t455 = t110 * t260;
t452 = t132 * t261;
t450 = t134 * t263;
t448 = t147 * t275;
t158 = t196 * t256;
t169 = t232 * t256;
t314 = t233 * t260;
t419 = t255 * t510 + t146;
t418 = t256 * t510 + t147;
t417 = -t234 * t255 + t148;
t416 = -Icges(5,2) * t432 + t149 - t229;
t406 = t223 - t225;
t405 = t230 + t249;
t254 = pkin(2) * t264;
t404 = -t231 + t254;
t384 = t136 * t436 + (-t124 + t86) * t255;
t379 = t224 + t407;
t378 = t254 + t406;
t253 = pkin(2) * t262;
t377 = t253 + t405;
t376 = t253 + t513;
t369 = t441 / 0.2e1;
t368 = -t436 / 0.2e1;
t365 = t399 / 0.2e1;
t363 = t398 / 0.2e1;
t361 = pkin(4) * t275 + t202;
t356 = -t131 - t452;
t355 = -t131 + t450;
t205 = rSges(3,1) * t264 - rSges(3,2) * t262;
t167 = t205 * t272 + t259;
t347 = t253 + t379;
t183 = rSges(3,1) * t425 - rSges(3,2) * t427;
t155 = rSges(4,1) * t436 - rSges(4,2) * t441;
t344 = t191 + t254;
t340 = -t255 * t75 + t256 * t74;
t88 = -t133 * t263 - t135 * t261;
t331 = -t149 * t277 + t448;
t330 = -t198 * t261 + t263 * t512;
t328 = -t234 * t275 + t277 * t510;
t325 = t155 + t240;
t121 = t148 * t437;
t65 = -t144 * t256 - t146 * t438 + t121;
t122 = t149 * t437;
t66 = t145 * t256 + t147 * t438 - t122;
t321 = (-t255 * t66 - t256 * t65) * qJD(4);
t123 = t146 * t433;
t67 = -t144 * t255 - t148 * t432 + t123;
t68 = t255 * t145 - t331 * t256;
t320 = (-t255 * t68 - t256 * t67) * qJD(4);
t90 = (t150 * t255 + t151 * t256) * qJD(4);
t307 = t158 * t188 - t189 * t445 - t197 * t260;
t304 = -t256 * t504 + (-t313 + t334) * t260;
t303 = -t197 * t436 + t504 * t255 + (t450 - t452) * t260;
t302 = -t255 * t314 - t256 * t501 + t260 * t331;
t301 = t255 * t501 - t256 * t314 + t260 * t333;
t300 = -t197 * t271 + t260 * t330;
t299 = -t233 * qJD(4) + t260 * t328;
t298 = t275 * t417 + t277 * t419;
t297 = t275 * t416 + t277 * t418;
t13 = t303 * t255 + t256 * t492;
t14 = t304 * t255 - t256 * t493;
t15 = -t255 * t492 + t303 * t256;
t16 = t255 * t493 + t304 * t256;
t113 = t134 * t439;
t61 = -t130 * t256 - t132 * t440 + t113;
t114 = t135 * t439;
t62 = t131 * t256 + t133 * t440 - t114;
t91 = t255 * t330 - t158;
t89 = t91 * t260;
t28 = -t188 * t62 - t189 * t61 + t89;
t115 = t132 * t435;
t63 = -t130 * t255 - t134 * t434 + t115;
t64 = t131 * t255 - t256 * t334;
t29 = -t188 * t64 - t189 * t63 - t468;
t40 = -t261 * t360 + t263 * t358;
t41 = t261 * t359 - t263 * t357;
t44 = t300 * t255 + t256 * t494;
t45 = -t255 * t494 + t300 * t256;
t87 = t132 * t263 + t134 * t261;
t295 = (-t13 * t189 - t14 * t188 + t152 * t63 - t153 * t64 + t260 * t44) * t484 + (t307 * t255 + t523 * t256) * t487 + (-t523 * t255 + t307 * t256) * t486 + (-t15 * t189 + t152 * t61 - t153 * t62 - t16 * t188 + t260 * t45) * t483 + (-t261 * t491 + t263 * t287) * t482 + t28 * t369 + t29 * t368 + ((-t260 * t64 - t13) * t256 + (t260 * t63 - t14) * t255) * t488 + (-t255 * t62 - t256 * t61) * t490 + (-t255 * t64 - t256 * t63) * t489 + ((-t260 * t62 - t15) * t256 + (t260 * t61 - t16) * t255) * t485 + ((-t260 * t88 - t40) * t256 + (t260 * t87 - t41) * t255) * t481;
t109 = t255 * t328 - t169;
t108 = t109 * t260;
t34 = t108 + t321;
t35 = t320 - t455;
t49 = qJD(4) * t331 + t101 * t275 + t277 * t99;
t53 = t299 * t255 + t256 * t499;
t54 = -t255 * t499 + t299 * t256;
t284 = (t89 - (t255 * t356 + t113 + t64) * t189 + (t114 - t115 + t63 + (t130 - t451) * t255) * t188 + (t188 * t355 - t515) * t256) * t487 + (t108 + ((t67 + t122 - t123 + (t144 - t448) * t255) * t255 + (-t121 - t68 + (t144 - t331) * t256 + (t447 + t449) * t255) * t256) * qJD(4)) * t365 + t88 * t489 + t153 * t92 / 0.2e1 + (t87 + t91) * t490 + (t468 - (-t115 + t62) * t189 + (-t113 + t61) * t188 + (-t188 * t323 - t189 * t355) * t256 + (-t188 * t356 + t515) * t255 + t29) * t486 + (t40 + t45) * t485 + (t35 + t455 + ((t123 - t66 + (t145 - t447) * t256) * t256 + (-t121 + t65 + (t145 + t449) * t255) * t255) * qJD(4)) * t363 + (qJD(4) * t328 + t207 * t277 + t208 * t275 - t261 * t352 + t263 * t353) * t260 + (t41 + t44 + t28) * t488 - (t49 + t53 + t34) * t399 / 0.2e1 - (t54 - t526) * t398 / 0.2e1 + (t256 * t110 + (t106 + t109) * t255) * qJD(4) * t481;
t176 = t241 * t256;
t175 = t241 * t255;
t139 = t183 * t272 + t258;
t138 = -t182 * t272 - t395;
t120 = t255 * t136;
t112 = t155 * t260 + t401;
t111 = -t260 * t518 + t319;
t12 = t136 * t153 - t137 * t152 + t188 * t86 - t189 * t85 + ((t127 * t260 - t95) * t256 + (t96 + t118) * t255) * qJD(4);
t1 = [m(3) * (t138 * (t205 + t269) + t139 * (t268 + t511) + (-t167 + t183 + t259) * t166) + t284 + (t32 * (t269 + t378) + t58 * (-t312 + t414) + t33 * (t268 + t347) + t282 + (t58 + t507) * t57) * m(6) + (t55 * (t269 + t404) + t75 * (-t312 + t380) + t56 * (t268 + t377) + t283 + (t75 + t508) * t74) * m(5) + m(4) * (t111 * (t269 + t344) + t112 * (t268 + t376) + (-t126 + t259 + t325) * t125); t284 + (t32 * t378 + t33 * t347 + (-t394 + t414 + t296) * t58 + t524) * m(6) + (t56 * t377 + t55 * t404 + t525) * m(5) + (t111 * t344 + t112 * t376 + (t325 - t240 - t174) * t125) * m(4) + (t138 * t205 + t139 * t511 + t166 * t183 - t167 * t182 - (t166 * t205 - t167 * t511) * t272) * m(3); t284 + (t32 * t406 + t33 * t379 + (-t326 + t414) * t58 + t524) * m(6) + (-t55 * t231 + t56 * t405 + t525) * m(5) + (t111 * t191 + t112 * t513 + t125 * t155 - t126 * t518 - (t125 * t191 - t126 * t513) * t260) * m(4); t295 + (t526 * t256 + (t106 * t260 - t49) * t255) * t481 + ((-t398 * t442 - t314) * t256 + (-t519 + (-t297 * t255 + (t169 + t298) * t256) * qJD(4)) * t255) * t363 + ((-t403 * t275 + t402 * t277) * t260 + ((t255 * t416 - t256 * t417) * t277 + (-t255 * t418 + t256 * t419) * t275) * qJD(4)) * t482 + ((t169 * t399 - t314) * t255 + (t519 + (-t298 * t256 + (-t442 + t297) * t255) * qJD(4)) * t256) * t365 + (t260 * t53 + ((-t301 * t255 - t256 * t502 - t260 * t68) * t256 + (-t302 * t255 + t256 * t498 + t260 * t67) * t255) * t516) * t484 + (t260 * t54 + ((t255 * t502 - t301 * t256 - t260 * t66) * t256 + (-t255 * t498 - t302 * t256 + t260 * t65) * t255) * t516) * t483 + (t321 + t34) * t369 + (t320 + t35) * t368 + (t12 * t120 + t50 * t384 + (t12 * t127 + t50 * t96 - t32 * t361 + t58 * (-pkin(4) * t396 - t181) + (t50 * t128 - t361 * t57) * t260) * t255 + (t12 * t424 + t50 * t476 + t33 * t361 + t57 * t181 + (t50 * t127 - t361 * t58) * t260) * t256 - (-t58 * t372 + ((-t255 * t57 - t256 * t58) * t260 + t50 * (-t255 ^ 2 - t256 ^ 2) * qJD(4)) * t275) * pkin(4) + t505) * m(6) + (0.2e1 * t90 * ((t150 * t260 - t104) * t256 + (t105 - t129) * t255) + t340 * t214 + ((t56 - t469) * t256 + (-t260 * t74 - t55) * t255) * t241 - (-t175 * t74 - t176 * t75) * t260 - (t90 * (-t175 * t255 - t176 * t256) + t340 * t341) * qJD(4)) * m(5); t295 + (t12 * (t137 * t256 + t120) + t50 * (-t256 * t85 + t384) + (-t255 * t58 + t256 * t57) * t181 + ((t33 - t470) * t256 + (-t260 * t57 - t32) * t255) * t202 + t505) * m(6);];
tauc = t1(:);
