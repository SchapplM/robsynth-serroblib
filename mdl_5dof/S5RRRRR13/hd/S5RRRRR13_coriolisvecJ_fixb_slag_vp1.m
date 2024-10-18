% Calculate vector of centrifugal and Coriolis load on the joints for
% S5RRRRR13
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% qJD [5x1]
%   Generalized joint velocities
% pkin [10x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,alpha4,d1,d2,d3,d4,d5]';
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
% Datum: 2024-09-27 17:33
% Revision: 4128c0dfad36ceb08eaa7a9fbc8f797c888f59be (2024-07-29)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauc = S5RRRRR13_coriolisvecJ_fixb_slag_vp1(qJ, qJD, ...
  pkin, m, rSges, Icges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(10,1),zeros(6,1),zeros(6,3),zeros(6,6)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RRRRR13_coriolisvecJ_fixb_slag_vp1: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RRRRR13_coriolisvecJ_fixb_slag_vp1: qJD has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S5RRRRR13_coriolisvecJ_fixb_slag_vp1: pkin has to be [10x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5RRRRR13_coriolisvecJ_fixb_slag_vp1: m has to be [6x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [6,3]), ...
  'S5RRRRR13_coriolisvecJ_fixb_slag_vp1: rSges has to be [6x3] (double)');
assert(isreal(Icges) && all(size(Icges) == [6 6]), ...
  'S5RRRRR13_coriolisvecJ_fixb_slag_vp1: Icges has to be [6x6] (double)'); 

%% Symbolic Calculation
% From coriolisvec_joint_fixb_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2024-09-27 17:30:24
% EndTime: 2024-09-27 17:30:39
% DurationCPUTime: 10.70s
% Computational Cost: add. (49409->692), mult. (29421->912), div. (0->0), fcn. (25705->16), ass. (0->350)
t332 = qJ(1) + qJ(2);
t325 = qJ(3) + t332;
t316 = sin(t325);
t317 = cos(t325);
t275 = rSges(4,1) * t316 + rSges(4,2) * t317;
t330 = qJD(1) + qJD(2);
t320 = qJD(3) + t330;
t255 = t320 * t275;
t336 = sin(qJ(1));
t495 = pkin(1) * qJD(1);
t449 = t336 * t495;
t322 = sin(t332);
t472 = t322 * t330;
t452 = pkin(2) * t472;
t208 = -t449 - t255 - t452;
t337 = cos(qJ(4));
t334 = cos(pkin(5));
t335 = sin(qJ(4));
t468 = t334 * t335;
t252 = t316 * t337 + t317 * t468;
t467 = t334 * t337;
t390 = -t316 * t335 + t317 * t467;
t333 = sin(pkin(5));
t475 = t317 * t333;
t169 = -Icges(5,5) * t252 - Icges(5,6) * t390 + Icges(5,3) * t475;
t239 = Icges(5,4) * t252;
t171 = Icges(5,2) * t390 - Icges(5,6) * t475 + t239;
t238 = Icges(5,4) * t390;
t175 = -Icges(5,1) * t252 + Icges(5,5) * t475 - t238;
t84 = -t334 * t169 + (t171 * t337 - t175 * t335) * t333;
t331 = qJ(4) + qJ(5);
t321 = sin(t331);
t435 = pkin(5) + t331;
t399 = cos(t435) / 0.2e1;
t436 = pkin(5) - t331;
t412 = cos(t436);
t359 = t412 / 0.2e1 + t399;
t219 = t316 * t321 - t317 * t359;
t411 = sin(t436);
t502 = sin(t435) / 0.2e1;
t283 = t502 - t411 / 0.2e1;
t323 = cos(t331);
t220 = t283 * t317 + t316 * t323;
t535 = Icges(6,4) * t220;
t147 = Icges(6,2) * t219 + Icges(6,6) * t475 - t535;
t222 = -t316 * t359 - t317 * t321;
t477 = t316 * t333;
t224 = t283 * t316 - t317 * t323;
t534 = Icges(6,4) * t224;
t148 = Icges(6,2) * t222 + Icges(6,6) * t477 - t534;
t398 = t411 / 0.2e1;
t282 = t502 + t398;
t284 = t399 - t412 / 0.2e1;
t483 = Icges(6,4) * t284;
t205 = Icges(6,2) * t282 + Icges(6,6) * t334 - t483;
t329 = qJD(4) + qJD(5);
t426 = t329 * t333;
t249 = t316 * t426;
t250 = t317 * t426;
t303 = qJD(4) * t334 + t320;
t281 = qJD(5) * t334 + t303;
t536 = (-Icges(6,1) * t222 + t148 - t534) * t249 + t250 * (-Icges(6,1) * t219 + t147 - t535) + (-Icges(6,1) * t282 + t205 - t483) * t281;
t404 = t171 * t390 - t175 * t252;
t72 = t169 * t475 + t404;
t544 = t317 * t72;
t144 = -Icges(6,5) * t220 + Icges(6,6) * t219 + Icges(6,3) * t475;
t210 = Icges(6,4) * t219;
t150 = -Icges(6,1) * t220 + Icges(6,5) * t475 + t210;
t540 = t222 * t147 - t224 * t150;
t56 = -t144 * t477 - t540;
t145 = -Icges(6,5) * t224 + Icges(6,6) * t222 + Icges(6,3) * t477;
t211 = Icges(6,4) * t222;
t151 = -Icges(6,1) * t224 + Icges(6,5) * t477 + t211;
t57 = t145 * t477 + t148 * t222 - t151 * t224;
t204 = -Icges(6,5) * t284 + Icges(6,6) * t282 + Icges(6,3) * t334;
t271 = Icges(6,4) * t282;
t206 = -Icges(6,1) * t284 + Icges(6,5) * t334 + t271;
t95 = t204 * t477 + t205 * t222 - t206 * t224;
t25 = t249 * t57 - t250 * t56 + t281 * t95;
t54 = t144 * t475 + t147 * t219 - t150 * t220;
t68 = -t144 * t334 - t147 * t282 + t150 * t284;
t253 = -t316 * t467 - t317 * t335;
t389 = t316 * t468 - t317 * t337;
t539 = -t253 * t171 - t175 * t389;
t537 = (-Icges(6,5) * t219 - Icges(6,6) * t220) * t250 - (Icges(6,5) * t222 + Icges(6,6) * t224) * t249 - (Icges(6,5) * t282 + Icges(6,6) * t284) * t281;
t257 = Icges(5,3) * t334 + (Icges(5,5) * t335 + Icges(5,6) * t337) * t333;
t485 = Icges(5,4) * t335;
t258 = Icges(5,6) * t334 + (Icges(5,2) * t337 + t485) * t333;
t484 = Icges(5,4) * t337;
t310 = t333 * t484;
t469 = t333 * t335;
t259 = Icges(5,1) * t469 + Icges(5,5) * t334 + t310;
t367 = -t252 * t259 + t257 * t475 - t258 * t390;
t530 = t367 * t303;
t152 = rSges(6,1) * t220 - rSges(6,2) * t219 - rSges(6,3) * t475;
t268 = pkin(3) * t316 - pkin(9) * t475;
t339 = pkin(9) + pkin(10);
t285 = pkin(4) * t468 - t333 * t339;
t497 = pkin(4) * t337;
t318 = pkin(3) + t497;
t456 = -t285 * t317 - t316 * t318;
t188 = t268 + t456;
t235 = t320 * t268;
t513 = t152 * t281 - t188 * t303 + t235;
t494 = pkin(4) * qJD(4);
t424 = t494 * t467;
t514 = -t318 * t320 - t424;
t344 = -t320 * t359 - t323 * t329;
t295 = t329 * t398;
t522 = t329 * t502;
t429 = t320 * t321 - t295 + t522;
t131 = t316 * t429 + t317 * t344;
t402 = -t283 * t320 - t321 * t329;
t294 = t329 * t399;
t396 = t329 * t412;
t430 = t320 * t323 + t294 + t396 / 0.2e1;
t132 = -t316 * t430 + t317 * t402;
t473 = t320 * t333;
t444 = t317 * t473;
t92 = rSges(6,1) * t132 + rSges(6,2) * t131 + rSges(6,3) * t444;
t528 = t316 * t514 + t513 + t92;
t183 = qJD(4) * t389 - t320 * t390;
t184 = qJD(4) * t253 - t252 * t320;
t108 = rSges(5,1) * t184 + rSges(5,2) * t183 + rSges(5,3) * t444;
t280 = pkin(9) * t444;
t177 = rSges(5,1) * t252 + rSges(5,2) * t390 - rSges(5,3) * t475;
t516 = t177 * t303 + t235;
t527 = t280 + t108 + t516;
t324 = cos(t332);
t286 = rSges(3,1) * t322 + rSges(3,2) * t324;
t479 = t286 * t330;
t236 = -t449 - t479;
t264 = pkin(4) * t469 + (-pkin(9) + t339) * t334;
t455 = qJD(4) * t333;
t443 = t316 * t455;
t269 = pkin(3) * t317 + pkin(9) * t477;
t518 = t269 * t320;
t526 = t264 * t443 - t518;
t368 = t204 * t475 + t205 * t219 - t206 * t220;
t525 = t250 * t54 + t368 * t281;
t447 = t335 * t494;
t394 = -t285 * t320 - t447;
t207 = -rSges(6,1) * t284 + rSges(6,2) * t282 + rSges(6,3) * t334;
t442 = t317 * t455;
t421 = t264 * t442;
t378 = -t207 * t250 - t421;
t360 = t378 - t452;
t347 = t360 - t449;
t65 = t347 - t513;
t154 = -rSges(6,1) * t224 + rSges(6,2) * t222 + rSges(6,3) * t477;
t428 = -t285 * t316 + t317 * t318;
t189 = t428 - t269;
t364 = -t154 * t281 - t189 * t303 + t207 * t249 + t526;
t471 = t324 * t330;
t451 = pkin(2) * t471;
t349 = t364 - t451;
t338 = cos(qJ(1));
t448 = t338 * t495;
t66 = -t349 + t448;
t521 = (t394 * t66 + t514 * t65) * t317;
t263 = rSges(5,3) * t334 + (rSges(5,1) * t335 + rSges(5,2) * t337) * t333;
t422 = t263 * t442;
t372 = -t422 - t452;
t356 = t372 - t449;
t96 = t356 - t516;
t489 = t317 * t96;
t179 = -rSges(5,1) * t389 + rSges(5,2) * t253 + rSges(5,3) * t477;
t388 = -t179 * t303 + t263 * t443 - t518;
t365 = t388 - t451;
t97 = -t365 + t448;
t520 = (-pkin(3) * t489 + (-t97 * pkin(3) + t96 * (-rSges(5,3) - pkin(9)) * t333) * t316) * t320;
t517 = t148 * t219 - t151 * t220;
t486 = Icges(5,4) * t389;
t173 = Icges(5,2) * t253 + Icges(5,6) * t477 - t486;
t240 = Icges(5,4) * t253;
t176 = -Icges(5,1) * t389 + Icges(5,5) * t477 + t240;
t466 = -t173 * t390 - t176 * t252;
t512 = -t452 + t527;
t511 = (pkin(3) - t318) * t320 - t424;
t510 = -t452 + t528;
t369 = t316 * (Icges(5,2) * t389 + t176 + t240) - t317 * (-Icges(5,2) * t252 - t175 + t238);
t350 = t249 * (Icges(6,2) * t224 + t151 + t211) - t250 * (-Icges(6,2) * t220 - t150 - t210) + t281 * (Icges(6,2) * t284 + t206 + t271);
t225 = t320 * t249;
t509 = t225 / 0.2e1;
t226 = t320 * t250;
t508 = t226 / 0.2e1;
t507 = -t249 / 0.2e1;
t506 = t249 / 0.2e1;
t505 = -t250 / 0.2e1;
t504 = t250 / 0.2e1;
t501 = t334 / 0.2e1;
t500 = pkin(1) * t336;
t499 = pkin(2) * t322;
t498 = pkin(2) * t330 ^ 2;
t326 = t338 * pkin(1);
t245 = t294 - t396 / 0.2e1;
t247 = t295 + t522;
t133 = t316 * t344 - t317 * t429;
t134 = t316 * t402 + t317 * t430;
t445 = t316 * t473;
t87 = Icges(6,5) * t134 + Icges(6,6) * t133 + Icges(6,3) * t445;
t89 = Icges(6,4) * t134 + Icges(6,2) * t133 + Icges(6,6) * t445;
t91 = Icges(6,1) * t134 + Icges(6,4) * t133 + Icges(6,5) * t445;
t30 = -t147 * t245 - t150 * t247 + t282 * t89 - t284 * t91 + t334 * t87;
t491 = t30 * t250;
t86 = Icges(6,5) * t132 + Icges(6,6) * t131 + Icges(6,3) * t444;
t88 = Icges(6,4) * t132 + Icges(6,2) * t131 + Icges(6,6) * t444;
t90 = Icges(6,1) * t132 + Icges(6,4) * t131 + Icges(6,5) * t444;
t31 = t148 * t245 + t151 * t247 + t282 * t88 - t284 * t90 + t334 * t86;
t490 = t31 * t249;
t488 = t68 * t225;
t69 = t145 * t334 + t148 * t282 - t151 * t284;
t487 = t69 * t226;
t482 = t169 * t316;
t170 = -Icges(5,5) * t389 + Icges(5,6) * t253 + Icges(5,3) * t477;
t481 = t170 * t317;
t478 = t316 * t320;
t476 = t317 * t320;
t460 = -t207 - t264;
t459 = -t285 * t478 - t316 * t447;
t274 = (Icges(5,1) * t337 - t485) * t333;
t458 = -t258 + t274;
t457 = -Icges(5,2) * t469 + t259 + t310;
t341 = qJD(1) ^ 2;
t454 = t341 * t500;
t453 = t341 * t326;
t93 = rSges(6,1) * t134 + rSges(6,2) * t133 + rSges(6,3) * t445;
t450 = t152 * t444 + t475 * t92 + t477 * t93;
t74 = -t169 * t477 - t539;
t75 = t170 * t477 + t173 * t253 - t176 * t389;
t441 = t477 / 0.2e1;
t440 = -t475 / 0.2e1;
t439 = t473 / 0.2e1;
t438 = -t455 / 0.2e1;
t437 = t455 / 0.2e1;
t120 = t316 * t511 + t317 * t394 - t280;
t121 = -pkin(9) * t445 - t317 * t511 + t459;
t19 = t152 * t226 - t154 * t225 + t249 * t93 + t250 * t92 + ((-t188 * t320 + t120) * t317 + (-t189 * t320 + t121) * t316) * t455;
t434 = t19 * (t152 * t477 + t154 * t475);
t287 = rSges(3,1) * t324 - rSges(3,2) * t322;
t276 = rSges(4,1) * t317 - rSges(4,2) * t316;
t161 = -rSges(6,1) * t219 - rSges(6,2) * t220;
t162 = rSges(6,1) * t222 + rSges(6,2) * t224;
t433 = t161 * t249 + t162 * t250;
t218 = rSges(6,1) * t282 + rSges(6,2) * t284;
t432 = t162 * t281 - t218 * t249;
t431 = -t161 * t281 - t218 * t250;
t427 = t333 ^ 2 * qJD(4) ^ 2 * t497;
t420 = t316 * t439;
t419 = t317 * t439;
t418 = t316 * t438;
t417 = t316 * t437;
t416 = t317 * t438;
t415 = t317 * t437;
t414 = t320 * t437;
t262 = rSges(3,1) * t471 - rSges(3,2) * t472;
t233 = rSges(4,1) * t476 - rSges(4,2) * t478;
t315 = pkin(2) * t324;
t413 = t276 + t315;
t185 = -qJD(4) * t252 + t253 * t320;
t186 = qJD(4) * t390 - t320 * t389;
t407 = rSges(5,1) * t186 + rSges(5,2) * t185;
t406 = -t316 * t97 - t489;
t405 = t179 + t269;
t403 = (Icges(5,5) * t390 - Icges(5,6) * t252) * t317 - (Icges(5,5) * t253 + Icges(5,6) * t389) * t316;
t401 = t316 * t414;
t400 = t317 * t414;
t55 = -t145 * t475 - t517;
t397 = -t276 * t320 - t451;
t395 = t315 + t405;
t387 = (rSges(5,1) * t337 - rSges(5,2) * t335) * t333;
t73 = -t170 * t475 - t466;
t386 = (t316 * t73 - t544) * t333;
t385 = (t316 * t75 - t317 * t74) * t333;
t383 = -t322 * t498 - t454;
t382 = -t324 * t498 - t453;
t381 = t428 + t154;
t272 = (Icges(5,5) * t337 - Icges(5,6) * t335) * t333;
t379 = -t275 - t499;
t376 = t320 * (-pkin(3) * t478 + t280) + t383;
t375 = t315 + t381;
t374 = -t233 - t451;
t99 = (t177 * t316 + t179 * t317) * t455;
t373 = -t152 + t456;
t371 = -t407 - t451;
t370 = (Icges(5,1) * t253 - t173 + t486) * t316 - (Icges(5,1) * t390 - t171 - t239) * t317;
t366 = -t177 - t268;
t363 = t373 - t499;
t15 = -t131 * t147 - t132 * t150 + t222 * t89 - t224 * t91 + (-t144 * t476 + t316 * t87) * t333;
t16 = t131 * t148 + t132 * t151 + t222 * t88 - t224 * t90 + (t145 * t476 + t316 * t86) * t333;
t17 = -t133 * t147 - t134 * t150 - t219 * t89 + t220 * t91 + (-t144 * t478 - t317 * t87) * t333;
t18 = t133 * t148 + t134 * t151 - t219 * t88 + t220 * t90 + (t145 * t478 - t317 * t86) * t333;
t24 = t249 * t55 - t525;
t190 = Icges(6,5) * t247 + Icges(6,6) * t245;
t191 = Icges(6,4) * t247 + Icges(6,2) * t245;
t192 = Icges(6,1) * t247 + Icges(6,4) * t245;
t49 = t131 * t205 + t132 * t206 + t191 * t222 - t192 * t224 + (t190 * t316 + t204 * t476) * t333;
t50 = t133 * t205 + t134 * t206 - t191 * t219 + t192 * t220 + (-t190 * t317 + t204 * t478) * t333;
t62 = t190 * t334 + t191 * t282 - t192 * t284 + t205 * t245 + t206 * t247;
t53 = t62 * t281;
t361 = (-t15 * t250 + t16 * t249 + t225 * t56 + t226 * t57 + t281 * t49) * t441 + (t350 * t222 + t224 * t536 - t477 * t537) * t507 + (-t350 * t219 - t220 * t536 + t537 * t475) * t504 - (t350 * t282 + t284 * t536 - t334 * t537) * t281 / 0.2e1 + (-t17 * t250 + t18 * t249 + t225 * t54 + t226 * t55 + t281 * t50) * t440 + t24 * t420 + t25 * t419 + (-t334 * t368 + (t316 * t55 - t317 * t54) * t333) * t509 + (t334 * t95 + (t316 * t57 - t317 * t56) * t333) * t508 + (t334 * t49 + ((t320 * t57 - t15) * t317 + (t320 * t56 + t16) * t316) * t333) * t506 + (t334 * t50 + ((t320 * t55 - t17) * t317 + (t320 * t54 + t18) * t316) * t333) * t505 + (t487 + t488 + t490 + t53 - t491) * t501 + t281 * (t334 * t62 + ((t320 * t69 - t30) * t317 + (t320 * t68 + t31) * t316) * t333) / 0.2e1;
t358 = -t93 - t459;
t355 = t366 - t499;
t102 = Icges(5,5) * t184 + Icges(5,6) * t183 + Icges(5,3) * t444;
t103 = Icges(5,5) * t186 + Icges(5,6) * t185 + Icges(5,3) * t445;
t104 = Icges(5,4) * t184 + Icges(5,2) * t183 + Icges(5,6) * t444;
t105 = Icges(5,4) * t186 + Icges(5,2) * t185 + Icges(5,6) * t445;
t106 = Icges(5,1) * t184 + Icges(5,4) * t183 + Icges(5,5) * t444;
t107 = Icges(5,1) * t186 + Icges(5,4) * t185 + Icges(5,5) * t445;
t353 = ((t320 * t75 - t105 * t253 + t107 * t389 - t171 * t183 + t175 * t184 - (t103 * t316 - t169 * t476) * t333) * t317 + (t320 * t74 + t104 * t253 - t106 * t389 + t173 * t183 + t176 * t184 + (t102 * t316 + t170 * t476) * t333) * t316) * t333;
t352 = ((t320 * t73 - t105 * t390 - t107 * t252 - t171 * t185 + t175 * t186 - (-t103 * t317 - t169 * t478) * t333) * t317 + (t320 * t72 + t104 * t390 + t106 * t252 + t173 * t185 + t176 * t186 + (-t102 * t317 + t170 * t478) * t333) * t316) * t333;
t47 = t103 * t334 + (t105 * t337 + t107 * t335 + (-t171 * t335 - t175 * t337) * qJD(4)) * t333;
t48 = t102 * t334 + (t104 * t337 + t106 * t335 + (-t173 * t335 + t176 * t337) * qJD(4)) * t333;
t85 = t170 * t334 + (t173 * t337 + t176 * t335) * t333;
t351 = ((t320 * t85 - t47) * t317 + (t320 * t84 + t48) * t316) * t333;
t346 = t358 - t451;
t113 = t253 * t258 + t257 * t477 - t259 * t389;
t111 = t113 * t303;
t39 = qJD(4) * t386 - t530;
t40 = qJD(4) * t385 + t111;
t265 = qJD(4) * t272;
t266 = (-Icges(5,2) * t335 + t484) * t455;
t267 = qJD(4) * t274;
t63 = t183 * t258 + t184 * t259 + t253 * t266 - t389 * t267 + (t257 * t476 + t265 * t316) * t333;
t64 = t185 * t258 + t186 * t259 + t390 * t266 + t252 * t267 + (t257 * t478 - t265 * t317) * t333;
t110 = t265 * t334 + (t266 * t337 + t267 * t335 + (-t258 * t335 + t259 * t337) * qJD(4)) * t333;
t98 = t110 * t303;
t342 = -t368 * t509 + t490 / 0.2e1 + t487 / 0.2e1 + t488 / 0.2e1 + t98 + t53 + t25 * t504 - t491 / 0.2e1 + (t111 + ((-t404 + t72 + t75) * t316 + (t73 + (t481 - t482) * t333 - t74 + t466) * t317) * t455) * t415 + t49 * t506 + t95 * t508 + ((t56 + (t144 * t316 + t145 * t317) * t333 + t517 + t540) * t249 + t525 + t24) * t507 + (t50 + t25) * t505 + (t48 + t63) * t417 + (t84 - t367) * t401 + (t113 + t85) * t400 + (t39 + t530 + (t544 + (t466 + t74 + (t481 + t482) * t333 + t539) * t316) * t455) * t418 + (t47 + t64 + t40) * t416;
t270 = qJD(4) * t387;
t242 = t253 * pkin(4);
t241 = t390 * pkin(4);
t237 = t287 * t330 + t448;
t229 = -t262 * t330 - t453;
t228 = -t330 * t479 - t454;
t209 = -t397 + t448;
t203 = -t233 * t320 + t382;
t202 = -t255 * t320 + t383;
t201 = rSges(5,1) * t253 + rSges(5,2) * t389;
t200 = rSges(5,1) * t390 - rSges(5,2) * t252;
t193 = rSges(6,1) * t247 + rSges(6,2) * t245;
t187 = t207 * t445;
t136 = t334 * t154;
t109 = rSges(5,3) * t445 + t407;
t83 = t334 * t92;
t71 = -t109 * t303 - t518 * t320 + (t263 * t478 - t270 * t317) * t455 + t382;
t70 = t108 * t303 + (-t263 * t476 - t270 * t316) * t455 + t376;
t52 = t152 * t249 + t154 * t250 + (-t188 * t316 + t189 * t317) * t455;
t44 = -t121 * t303 - t193 * t250 + t207 * t225 - t281 * t93 - t317 * t427 + t320 * t526 + t382;
t43 = t120 * t303 - t193 * t249 - t207 * t226 + t281 * t92 - t316 * t427 - t320 * t421 + t376;
t1 = [t342 + m(3) * (t229 * (-t286 - t500) + t228 * (t287 + t326) + (-t262 - t448 + t237) * t236) + (t44 * (t363 - t500) + t65 * (t346 - t448) + t43 * (t326 + t375) + (-t449 - t347 + t65 + t510) * t66 + t521) * m(6) + (t71 * (t355 - t500) + t96 * (t371 - t448) + t70 * (t326 + t395) + (-t356 + t96 - t449 + t512) * t97 + t520) * m(5) + m(4) * (t203 * (t379 - t500) + t202 * (t326 + t413) + (t374 - t448 + t209) * t208); t342 + (t44 * t363 + t43 * t375 + (-t360 + t510) * t66 + (-t349 + t346) * t65 + t521) * m(6) + (t71 * t355 + t70 * t395 + (-t372 + t512) * t97 + (-t365 + t371) * t96 + t520) * m(5) + (t202 * t413 + t203 * t379 + (t374 - t397) * t208) * m(4) + (t228 * t287 - t229 * t286 - t236 * t262 - t237 * t479 - (-t236 * t287 - t237 * t286) * t330) * m(3); t342 + (t44 * t373 + t43 * t381 + (-t378 + t528) * t66 + (t358 - t364) * t65 + t521) * m(6) + (t71 * t366 + t70 * t405 + (t422 + t527) * t97 + (-t388 - t407) * t96 + t520) * m(5) + (-(-t208 * t276 - t209 * t275) * t320 + t202 * t276 - t203 * t275 - t208 * t233 - t209 * t255) * m(4); t361 + ((t253 * t457 + t272 * t477 - t389 * t458) * t303 + (t253 * t369 - t370 * t389 - t403 * t477) * t455) * t418 + ((t252 * t458 - t272 * t475 + t390 * t457) * t303 + (t370 * t252 + t369 * t390 + t403 * t475) * t455) * t415 + (t334 * t64 + t352) * t416 + (t334 * t63 + t353) * t417 + t40 * t419 + t39 * t420 + (qJD(4) * t352 + t303 * t64) * t440 + (t113 * t334 + t385) * t400 + (qJD(4) * t351 + t98) * t501 + t303 * (t110 * t334 + t351) / 0.2e1 + (qJD(4) * t353 + t303 * t63) * t441 + (-t334 * t367 + t386) * t401 - t303 * (t334 * t272 * t303 + ((t335 * t458 + t337 * t457) * t303 + ((t335 * t370 + t337 * t369) * t333 - t403 * t334) * qJD(4)) * t333) / 0.2e1 + (t65 * t187 + t43 * t136 + t66 * t83 + t434 + t52 * t450 + (t44 * (-t152 + t188) + t65 * (-t121 - t93) + t43 * t189 + t66 * t120) * t334 + ((t44 * t460 - t65 * t193 + t19 * t189 + t52 * t120 + (-t188 * t52 + t460 * t66) * t320) * t317 + (t43 * t460 - t66 * t193 - t19 * t188 + t52 * t121 + (t65 * t264 + t52 * (-t154 - t189)) * t320) * t316) * t333 - t65 * (-t241 * t303 + t431) - t66 * (t242 * t303 + t432) - t52 * ((t241 * t316 + t242 * t317) * t455 + t433)) * m(6) + ((t108 * t97 - t109 * t96 - t177 * t71 + t179 * t70) * t334 + (0.2e1 * t99 * ((t177 * t320 + t108) * t317 + (-t179 * t320 + t109) * t316) + t406 * t270 + ((-t320 * t97 - t71) * t317 + (t320 * t96 - t70) * t316) * t263) * t333 - (-t200 * t96 + t201 * t97) * t303 - (t99 * (t200 * t316 + t201 * t317) + t406 * t387) * t455) * m(5); t361 + (t44 * (-t152 * t334 - t207 * t475) + t43 * (-t207 * t477 + t136) + t434 + (-t432 + t83 + (-t193 * t316 - t207 * t476) * t333) * t66 + (-t193 * t475 - t334 * t93 + t187 - t431) * t65 + (-t154 * t445 - t433 + t450) * t52) * m(6);];
tauc = t1(:);
