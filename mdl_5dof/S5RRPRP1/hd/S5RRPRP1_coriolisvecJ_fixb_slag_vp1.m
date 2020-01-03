% Calculate vector of centrifugal and Coriolis load on the joints for
% S5RRPRP1
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% qJD [5x1]
%   Generalized joint velocities
% pkin [8x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d1,d2,d4,theta3]';
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
% Datum: 2020-01-03 11:59
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauc = S5RRPRP1_coriolisvecJ_fixb_slag_vp1(qJ, qJD, ...
  pkin, m, rSges, Icges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(8,1),zeros(6,1),zeros(6,3),zeros(6,6)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RRPRP1_coriolisvecJ_fixb_slag_vp1: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RRPRP1_coriolisvecJ_fixb_slag_vp1: qJD has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S5RRPRP1_coriolisvecJ_fixb_slag_vp1: pkin has to be [8x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5RRPRP1_coriolisvecJ_fixb_slag_vp1: m has to be [6x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [6,3]), ...
  'S5RRPRP1_coriolisvecJ_fixb_slag_vp1: rSges has to be [6x3] (double)');
assert(isreal(Icges) && all(size(Icges) == [6 6]), ...
  'S5RRPRP1_coriolisvecJ_fixb_slag_vp1: Icges has to be [6x6] (double)'); 

%% Symbolic Calculation
% From coriolisvec_joint_fixb_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2020-01-03 11:58:54
% EndTime: 2020-01-03 11:59:09
% DurationCPUTime: 11.31s
% Computational Cost: add. (12262->444), mult. (9083->534), div. (0->0), fcn. (6938->8), ass. (0->273)
t246 = qJ(1) + qJ(2);
t236 = pkin(8) + t246;
t229 = sin(t236);
t230 = cos(t236);
t250 = cos(qJ(4));
t392 = t230 * t250;
t248 = sin(qJ(4));
t393 = t230 * t248;
t122 = Icges(6,4) * t392 - Icges(6,2) * t393 + Icges(6,6) * t229;
t124 = Icges(5,4) * t392 - Icges(5,2) * t393 + Icges(5,6) * t229;
t523 = t122 + t124;
t531 = Icges(5,5) + Icges(6,5);
t530 = -Icges(5,6) - Icges(6,6);
t529 = Icges(5,3) + Icges(6,3);
t194 = Icges(6,4) * t393;
t126 = Icges(6,1) * t392 + Icges(6,5) * t229 - t194;
t195 = Icges(5,4) * t393;
t128 = Icges(5,1) * t392 + Icges(5,5) * t229 - t195;
t527 = t126 + t128;
t201 = Icges(6,5) * t250 - Icges(6,6) * t248;
t203 = Icges(5,5) * t250 - Icges(5,6) * t248;
t528 = t201 + t203;
t239 = Icges(6,4) * t250;
t311 = -Icges(6,2) * t248 + t239;
t121 = -Icges(6,6) * t230 + t229 * t311;
t240 = Icges(5,4) * t250;
t312 = -Icges(5,2) * t248 + t240;
t123 = -Icges(5,6) * t230 + t229 * t312;
t524 = t121 + t123;
t414 = Icges(6,4) * t248;
t209 = Icges(6,1) * t250 - t414;
t290 = t209 * t229;
t125 = -Icges(6,5) * t230 + t290;
t415 = Icges(5,4) * t248;
t211 = Icges(5,1) * t250 - t415;
t127 = -Icges(5,5) * t230 + t211 * t229;
t522 = t125 + t127;
t526 = t311 + t312;
t525 = t523 * t248;
t245 = qJD(1) + qJD(2);
t387 = t245 * t250;
t477 = t229 * t528 - t230 * t529;
t476 = t229 * t529 + t392 * t531 + t393 * t530;
t204 = Icges(6,2) * t250 + t414;
t206 = Icges(5,2) * t250 + t415;
t521 = t204 + t206;
t456 = Icges(6,1) * t248 + t239;
t302 = -t204 * t248 + t250 * t456;
t455 = Icges(5,1) * t248 + t240;
t514 = -t206 * t248 + t250 * t455 + t302;
t520 = -t250 * t527 + t525;
t519 = t455 + t456;
t388 = t245 * t248;
t395 = t229 * t250;
t479 = t522 * t395;
t518 = t527 * t395;
t478 = t524 * t393;
t517 = t526 * qJD(4);
t516 = (t209 + t211) * qJD(4);
t200 = Icges(6,5) * t248 + Icges(6,6) * t250;
t202 = Icges(5,5) * t248 + Icges(5,6) * t250;
t515 = t200 + t202;
t513 = (t519 * qJD(4) - t245 * t531) * t248 + (t521 * qJD(4) + t245 * t530) * t250;
t396 = t229 * t248;
t486 = -t230 * t477 - t396 * t524 + t479;
t485 = -t230 * t476 - t396 * t523 + t518;
t484 = -t229 * t477 - t392 * t522 + t478;
t483 = -t229 * t476 + t230 * t520;
t138 = t200 * t230;
t140 = t202 * t230;
t512 = t229 * t514 - t138 - t140;
t398 = t202 * t229;
t399 = t200 * t229;
t511 = -t206 * t393 + t230 * t302 + t392 * t455 + t398 + t399;
t247 = -qJ(5) - pkin(7);
t510 = rSges(6,3) - t247;
t468 = t522 * t250;
t467 = t524 * t248;
t507 = -t211 * t388 - t387 * t526;
t431 = pkin(4) * t250;
t233 = pkin(3) + t431;
t506 = -rSges(6,2) * t393 + t230 * t233;
t505 = -qJD(4) * t528 + t514 * t245;
t286 = t201 * t245;
t287 = t203 * t245;
t504 = t286 + t287;
t502 = t467 - t468;
t501 = -t516 * t250 + t517 * t248 - t515 * t245 + (t248 * t519 + t250 * t521) * qJD(4);
t500 = t515 * qJD(4) - t245 * t529;
t499 = t512 * t245;
t394 = t230 * t245;
t432 = pkin(3) * t230;
t383 = -rSges(6,1) * t392 + t432 - t506 + (pkin(7) - t510) * t229;
t498 = t245 * t383;
t497 = (-t229 * t500 + t230 * t504 + t245 * t502) * t230;
t162 = pkin(7) * t229 + t432;
t238 = cos(t246);
t390 = t238 * t245;
t214 = pkin(2) * t390;
t496 = -t162 * t245 + t214;
t495 = (t229 * t483 - t230 * t484) * qJD(4);
t494 = (t229 * t485 - t230 * t486) * qJD(4);
t237 = sin(t246);
t458 = rSges(3,1) * t237 + rSges(3,2) * t238;
t154 = t458 * t245;
t249 = sin(qJ(1));
t423 = pkin(1) * qJD(1);
t354 = t249 * t423;
t134 = t354 + t154;
t493 = t511 * t245;
t492 = t494 + t499;
t491 = -t493 + t495;
t490 = -t209 * t394 * t248 + t502 * qJD(4) + t229 * t513 + t507 * t230;
t489 = -qJD(4) * t520 + t229 * t507 - t230 * t513 - t290 * t388;
t488 = t229 * t505 + t230 * t501;
t487 = -t229 * t501 + t230 * t505;
t482 = t248 * t522 + t250 * t524;
t366 = t456 + t311;
t367 = t204 - t209;
t481 = (t248 * t366 + t250 * t367) * t245;
t364 = t455 + t312;
t365 = t206 - t211;
t480 = (t248 * t364 + t250 * t365) * t245;
t426 = rSges(4,2) * t229;
t160 = rSges(4,1) * t230 - t426;
t153 = t245 * t160;
t186 = rSges(4,1) * t394;
t475 = t186 - t153;
t473 = t229 * t504 + t230 * t500 - t245 * t520;
t199 = rSges(5,2) * t393;
t132 = rSges(5,1) * t392 + rSges(5,3) * t229 - t199;
t114 = t245 * t132;
t217 = rSges(5,1) * t248 + rSges(5,2) * t250;
t361 = qJD(4) * t229;
t321 = -t217 * t361 + t214;
t397 = t229 * t245;
t369 = pkin(3) * t394 + pkin(7) * t397;
t352 = t230 * t387;
t371 = rSges(5,1) * t352 + rSges(5,3) * t397;
t472 = -t114 - t321 + t369 + t371 + t496;
t471 = t248 * t527 + t250 * t523;
t231 = pkin(2) * t237;
t319 = -rSges(4,1) * t229 - rSges(4,2) * t230;
t470 = t319 - t231;
t462 = rSges(6,1) * t352 + rSges(6,3) * t397 + t233 * t394;
t469 = t462 + t496 + t498;
t464 = 0.2e1 * qJD(4);
t463 = t245 * t470;
t188 = pkin(7) * t394;
t353 = t229 * t388;
t373 = -rSges(5,2) * t353 - rSges(5,3) * t394;
t461 = t188 - t373;
t215 = t230 * t247;
t460 = rSges(6,1) * t395 + t229 * t233 + t215;
t222 = qJD(5) * t229;
t459 = rSges(6,2) * t353 + rSges(6,3) * t394 + t222;
t226 = t229 * pkin(3);
t161 = -pkin(7) * t230 + t226;
t334 = t161 + t231;
t384 = -rSges(6,2) * t396 - rSges(6,3) * t230 - t161 + t460;
t451 = -t245 * (t334 + t384) + t222;
t244 = t245 ^ 2;
t435 = t245 / 0.2e1;
t434 = rSges(5,3) + pkin(7);
t433 = pkin(2) * t244;
t242 = t249 * pkin(1);
t251 = cos(qJ(1));
t243 = t251 * pkin(1);
t359 = qJD(4) * t248;
t345 = t230 * t359;
t283 = t229 * t387 + t345;
t358 = qJD(4) * t250;
t344 = t230 * t358;
t430 = -pkin(4) * t345 - t188 - (t215 + (-pkin(3) + t233) * t229) * t245 - rSges(6,1) * t283 - rSges(6,2) * t344 + t459;
t282 = -t229 * t358 - t230 * t388;
t346 = t229 * t359;
t357 = qJD(5) * t230;
t389 = t245 * t247;
t429 = -t357 + (-pkin(4) * t359 - t389) * t229 - t369 - rSges(6,1) * t346 + rSges(6,2) * t282 + t462;
t428 = rSges(5,1) * t250;
t427 = rSges(6,1) * t250;
t425 = rSges(5,2) * t248;
t424 = rSges(6,2) * t248;
t235 = t251 * t423;
t60 = t235 + (t132 + t162) * t245 + t321;
t422 = t245 * t60;
t421 = t250 * rSges(6,2);
t391 = t237 * t245;
t382 = t229 * t456 + t121;
t381 = t230 * t456 + t122;
t380 = t229 * t455 + t123;
t379 = t230 * t455 + t124;
t378 = -t204 * t229 + t125;
t377 = -Icges(6,2) * t392 + t126 - t194;
t376 = -t206 * t229 + t127;
t375 = -Icges(5,2) * t392 + t128 - t195;
t232 = pkin(2) * t238;
t368 = -t199 + t232;
t253 = qJD(1) ^ 2;
t234 = t253 * t243;
t363 = t238 * t433 + t234;
t362 = t214 + t235;
t360 = qJD(4) * t230;
t356 = pkin(2) * t391;
t355 = t253 * t242;
t351 = t245 * t369 + t363;
t350 = t232 + t506;
t197 = rSges(5,1) * t395;
t349 = t197 + t226 + t231;
t347 = t217 * t360;
t341 = pkin(3) + t428;
t339 = t361 / 0.2e1;
t338 = -t360 / 0.2e1;
t337 = t360 / 0.2e1;
t216 = rSges(6,1) * t248 + t421;
t333 = pkin(4) * t248 + t216;
t169 = rSges(3,1) * t238 - rSges(3,2) * t237;
t135 = t169 * t245 + t235;
t331 = t231 + t460;
t327 = -pkin(4) * t396 - t216 * t229;
t326 = -pkin(4) * t393 - t216 * t230;
t130 = -rSges(5,2) * t396 - rSges(5,3) * t230 + t197;
t325 = t130 + t334;
t155 = rSges(3,1) * t390 - rSges(3,2) * t391;
t324 = t160 + t232;
t49 = t333 * t360 + t354 - t451;
t323 = t49 * t333;
t316 = qJD(4) * t333;
t264 = -t229 * t316 - t357;
t262 = t264 + t362;
t50 = (t162 - t383) * t245 + t262;
t322 = t50 * t333;
t318 = -t425 + t428;
t219 = -t424 + t427;
t59 = t245 * t325 + t347 + t354;
t317 = -t229 * t60 + t230 * t59;
t304 = t130 * t229 + t132 * t230;
t178 = t219 * qJD(4);
t298 = (qJD(4) * t431 + t178) * qJD(4);
t297 = qJD(4) * t217;
t292 = -t237 * t433 - t355;
t285 = -t354 - t356;
t275 = qJD(4) * (-t421 + (-rSges(6,1) - pkin(4)) * t248);
t268 = t248 * t378 + t250 * t382;
t267 = t248 * t377 + t250 * t381;
t266 = t248 * t376 + t250 * t380;
t265 = t248 * t375 + t250 * t379;
t91 = rSges(5,1) * t283 + rSges(5,2) * t344 + t373;
t93 = -rSges(5,1) * t346 + rSges(5,2) * t282 + t371;
t263 = (t130 * t245 - t91) * t230 + (t93 - t114) * t229;
t111 = t354 - t463;
t112 = t362 + t153;
t257 = (-t111 * t426 + t112 * t470) * t245;
t256 = ((((t477 - t520) * t230 - t479 + t483) * t230 + ((t477 - t525) * t229 + (t467 + t468) * t230 - t478 + t484 + t518) * t229) * qJD(4) + t499) * t339 + (t487 - t490) * t338 + ((((-t468 + t476) * t230 + t478 + t485) * t230 + ((t467 + t476) * t229 - t479 + t486) * t229) * qJD(4) + t491 + t493) * t337 - (t488 - t489 + t492) * t361 / 0.2e1 + (t511 * t230 + (t482 + t512) * t229) * qJD(4) * t435 + (qJD(4) * t514 + t248 * t516 + t250 * t517 - t471 * t338) * t245;
t136 = pkin(3) * t397 - t188;
t179 = t318 * qJD(4);
t47 = -t179 * t361 + (-t136 - t91 - t347) * t245 + t292;
t48 = t245 * t93 + (t179 * t230 - t217 * t397) * qJD(4) + t351;
t255 = (-t297 * t59 - t341 * t422 - t425 * t48 + t434 * t47) * t229 + (-rSges(5,2) * t388 * t59 - t297 * t60 + t341 * t47 - t434 * t48) * t230;
t23 = -t298 * t229 + (-t230 * t316 - t136 + t222 + t430) * t245 + t292;
t24 = t298 * t230 + (t264 + t429) * t245 + t351;
t254 = (t23 * t510 - t24 * t424 + (t50 * (-t233 - t427) - t49 * t247) * t245 + t49 * t275) * t229 + (t23 * t427 - t24 * rSges(6,3) + t49 * (-rSges(6,2) * t388 - qJD(5)) + (-t389 + t275) * t50) * t230;
t152 = t217 * t230;
t150 = t217 * t229;
t116 = t155 * t245 + t234;
t115 = -t154 * t245 - t355;
t101 = t245 * (-rSges(4,2) * t397 + t186) + t363;
t100 = t244 * t319 + t292;
t65 = qJD(4) * t304 + qJD(3);
t42 = qJD(3) + (t229 * t384 - t230 * t383) * qJD(4);
t25 = t263 * qJD(4);
t5 = ((t245 * t384 + t430) * t230 + (t429 + t498) * t229) * qJD(4);
t1 = [m(3) * (t115 * (t169 + t243) + t116 * (t242 + t458) + (-t135 + t155 + t235) * t134) + t256 + (t23 * (t243 + t350) + t50 * (t285 + t459) + t24 * (t242 + t331) + t254 + (t235 + t50 - t262 + t469) * t49) * m(6) + (t47 * (t243 + t368) + t60 * (t285 + t461) + t48 * (t242 + t349) + t255 + (t60 + t472) * t59) * m(5) + (t100 * (t243 + t324) - t112 * t354 + t101 * (t242 - t470) + t257 + (t112 + t475) * t111) * m(4); t256 + (t23 * t350 + t24 * t331 + t254 - (-t229 * t323 - t230 * t322) * qJD(4) + (-t356 + t459 - t451) * t50 + (-t214 + t357 + t469) * t49) * m(6) + (t325 * t422 + t48 * t349 + t47 * t368 + t255 + (t347 - t356 + t461) * t60 + t472 * t59) * m(5) + (t100 * t324 - t101 * t470 + t475 * t111 - t112 * t463 + t257) * m(4) + (t115 * t169 + t116 * t458 + t134 * t155 - t135 * t154 - (t134 * t169 - t135 * t458) * t245) * m(3); m(5) * t25 + m(6) * t5; -(((t364 + t366) * t250 + (-t365 - t367) * t248) * t245 + (((-t376 - t378) * t230 + (t375 + t377) * t229) * t250 + ((t380 + t382) * t230 + (-t379 - t381) * t229) * t248) * qJD(4)) * t245 / 0.2e1 + ((t245 * t471 + t490) * t230 + (t245 * t482 + t489) * t229) * t435 + ((t140 * t361 - t287) * t229 + (t480 + (-t266 * t230 + (-t398 + t265) * t229) * qJD(4)) * t230 + (t138 * t361 - t286) * t229 + (t481 + (-t268 * t230 + (-t399 + t267) * t229) * qJD(4)) * t230) * t339 + ((-t360 * t398 - t287) * t230 + (-t480 + (-t265 * t229 + (t140 + t266) * t230) * qJD(4)) * t229 + (-t360 * t399 - t286) * t230 + (-t481 + (-t267 * t229 + (t138 + t268) * t230) * qJD(4)) * t229) * t337 + (t25 * t304 + t65 * t263 + t317 * t179 + ((t48 - t422) * t230 + (-t245 * t59 - t47) * t229) * t217 - (-t150 * t59 - t152 * t60) * t245 - (t65 * (-t150 * t229 - t152 * t230) + t317 * t318) * qJD(4)) * m(5) - (t488 * t245 + (t483 * t394 + (t473 * t229 + t484 * t245 + t497) * t229) * t464) * t229 / 0.2e1 - (t487 * t245 + ((t485 * t245 + t497) * t230 + (t230 * t473 + t486 * t245) * t229) * t464) * t230 / 0.2e1 + (-(t326 * t50 + t327 * t49) * t245 - ((t219 * t49 + t326 * t42) * t230 + (t50 * (-t219 - t431) + t327 * t42) * t229) * qJD(4) + (-t5 * t383 + t42 * t430 + t24 * t333 + t49 * t178 + (t384 * t42 - t322) * t245) * t230 + (t5 * t384 + t42 * t429 - t23 * t333 + t50 * (-pkin(4) * t358 - t178) + (t383 * t42 - t323) * t245) * t229) * m(6) + (t492 + t494) * t397 / 0.2e1 - (t491 + t495) * t394 / 0.2e1; m(6) * (-t229 * t24 - t23 * t230);];
tauc = t1(:);
