% Calculate vector of centrifugal and Coriolis load on the joints for
% S5RRPRR4
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
% Datum: 2019-12-05 18:32
% Revision: 77da58f92bca3eff71542919beafa37024070d86 (2019-12-05)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauc = S5RRPRR4_coriolisvecJ_fixb_slag_vp1(qJ, qJD, ...
  pkin, m, rSges, Icges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(9,1),zeros(6,1),zeros(6,3),zeros(6,6)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RRPRR4_coriolisvecJ_fixb_slag_vp1: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RRPRR4_coriolisvecJ_fixb_slag_vp1: qJD has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S5RRPRR4_coriolisvecJ_fixb_slag_vp1: pkin has to be [9x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5RRPRR4_coriolisvecJ_fixb_slag_vp1: m has to be [6x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [6,3]), ...
  'S5RRPRR4_coriolisvecJ_fixb_slag_vp1: rSges has to be [6x3] (double)');
assert(isreal(Icges) && all(size(Icges) == [6 6]), ...
  'S5RRPRR4_coriolisvecJ_fixb_slag_vp1: Icges has to be [6x6] (double)'); 

%% Symbolic Calculation
% From coriolisvec_joint_fixb_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-05 18:31:49
% EndTime: 2019-12-05 18:32:01
% DurationCPUTime: 8.65s
% Computational Cost: add. (17736->546), mult. (11158->695), div. (0->0), fcn. (8628->10), ass. (0->352)
t275 = qJ(1) + qJ(2);
t264 = pkin(9) + t275;
t261 = cos(t264);
t276 = sin(qJ(4));
t440 = t261 * t276;
t242 = rSges(5,2) * t440;
t278 = cos(qJ(4));
t439 = t261 * t278;
t401 = rSges(5,1) * t439;
t260 = sin(t264);
t484 = rSges(5,3) * t260;
t328 = -t401 - t484;
t154 = -t242 - t328;
t273 = qJD(1) + qJD(2);
t142 = t273 * t154;
t251 = rSges(5,1) * t276 + rSges(5,2) * t278;
t410 = qJD(4) * t260;
t187 = t251 * t410;
t493 = pkin(7) * t260;
t495 = pkin(3) * t261;
t192 = t493 + t495;
t268 = cos(t275);
t437 = t268 * t273;
t405 = pkin(2) * t437;
t350 = -t273 * t192 - t405;
t408 = qJD(4) * t276;
t384 = t260 * t408;
t407 = qJD(4) * t278;
t386 = rSges(5,1) * t384 + (t260 * t407 + t273 * t440) * rSges(5,2);
t556 = t350 + t187 - t142 - t386;
t274 = qJ(4) + qJ(5);
t265 = sin(t274);
t267 = cos(t274);
t210 = rSges(6,1) * t265 + rSges(6,2) * t267;
t255 = t260 * rSges(4,2);
t489 = rSges(4,1) * t261;
t191 = -t255 + t489;
t447 = t260 * t273;
t233 = rSges(4,2) * t447;
t497 = pkin(2) * t268;
t349 = -t489 - t497;
t555 = -t405 - t233 + (-t191 - t349) * t273;
t449 = t260 * t265;
t220 = Icges(6,4) * t449;
t448 = t260 * t267;
t139 = -Icges(6,1) * t448 + Icges(6,5) * t261 + t220;
t443 = t261 * t265;
t221 = Icges(6,4) * t443;
t442 = t261 * t267;
t140 = Icges(6,1) * t442 + Icges(6,5) * t260 - t221;
t272 = qJD(4) + qJD(5);
t189 = t260 * t272;
t190 = t261 * t272;
t258 = Icges(6,4) * t267;
t339 = -Icges(6,2) * t265 + t258;
t535 = Icges(6,1) * t265 + t258;
t548 = t535 + t339;
t291 = t189 * (-Icges(6,2) * t442 + t140 - t221) + t190 * (Icges(6,2) * t448 + t139 + t220) + t273 * t548;
t137 = Icges(6,6) * t261 - t260 * t339;
t138 = Icges(6,4) * t442 - Icges(6,2) * t443 + Icges(6,6) * t260;
t475 = Icges(6,4) * t265;
t206 = Icges(6,2) * t267 + t475;
t209 = Icges(6,1) * t267 - t475;
t513 = t189 * (t261 * t535 + t138) + t190 * (-t260 * t535 + t137) + t273 * (t206 - t209);
t554 = t291 * t265 + t267 * t513;
t553 = 0.2e1 * qJD(4);
t162 = t210 * t260;
t454 = t210 * t261;
t400 = rSges(6,1) * t448;
t549 = t210 * t272;
t388 = -t261 * t549 - t273 * t400;
t536 = rSges(6,2) * t449 + t261 * rSges(6,3);
t89 = t273 * t536 + t388;
t552 = t162 * t189 + t190 * t454 + t261 * t89;
t269 = Icges(5,4) * t278;
t340 = -Icges(5,2) * t276 + t269;
t534 = Icges(5,1) * t276 + t269;
t412 = t534 + t340;
t476 = Icges(5,4) * t276;
t245 = Icges(5,2) * t278 + t476;
t248 = Icges(5,1) * t278 - t476;
t413 = t245 - t248;
t551 = (t276 * t412 + t278 * t413) * t273;
t494 = pkin(4) * t278;
t262 = pkin(3) + t494;
t222 = t261 * t262;
t280 = -pkin(8) - pkin(7);
t490 = pkin(7) + t280;
t127 = t260 * t490 - t222 + t495;
t224 = rSges(6,2) * t443;
t394 = rSges(6,1) * t442;
t327 = rSges(6,3) * t260 + t394;
t143 = -t224 + t327;
t430 = t127 - t143;
t266 = sin(t275);
t346 = rSges(3,1) * t266 + rSges(3,2) * t268;
t183 = t346 * t273;
t277 = sin(qJ(1));
t483 = pkin(1) * qJD(1);
t397 = t277 * t483;
t166 = t183 + t397;
t486 = rSges(6,1) * t267;
t211 = -rSges(6,2) * t265 + t486;
t550 = t273 * t162 - t190 * t211 - t210 * t447;
t498 = pkin(2) * t266;
t547 = -t536 + t498;
t491 = pkin(3) - t262;
t539 = t490 * t261;
t546 = -t260 * t491 + t539;
t182 = t211 * t272;
t441 = t261 * t273;
t231 = pkin(4) * t384;
t357 = t189 * t210 + t231;
t370 = -t192 - t497;
t279 = cos(qJ(1));
t396 = t279 * t483;
t57 = -t396 + (t370 + t430) * t273 + t357;
t545 = (t182 * t260 - t189 * t211 + t210 * t441 - t273 * t454) * t57;
t136 = Icges(6,5) * t442 - Icges(6,6) * t443 + Icges(6,3) * t260;
t63 = t261 * t136 + t138 * t449 - t140 * t448;
t335 = t206 * t265 - t267 * t535;
t204 = Icges(6,5) * t265 + Icges(6,6) * t267;
t458 = t204 * t261;
t92 = t260 * t335 + t458;
t544 = t189 * t63 + t92 * t273;
t150 = Icges(5,4) * t439 - Icges(5,2) * t440 + Icges(5,6) * t260;
t240 = Icges(5,4) * t440;
t152 = Icges(5,1) * t439 + Icges(5,5) * t260 - t240;
t336 = t150 * t276 - t152 * t278;
t542 = t261 * t336;
t338 = t138 * t265 - t140 * t267;
t540 = t338 * t261;
t537 = -t222 - t497;
t446 = t260 * t276;
t241 = rSges(5,2) * t446;
t414 = t261 * rSges(5,3) + t241;
t389 = t273 * t224 + t260 * t549;
t435 = t273 * t280;
t417 = t260 * t435 + t231;
t533 = t273 * t430 + t350 + t357 - t389 - t417;
t488 = rSges(5,1) * t278;
t344 = -rSges(5,2) * t276 + t488;
t225 = t344 * qJD(4);
t532 = -t225 * t261 + t251 * t447;
t455 = t206 * t272;
t529 = -Icges(6,6) * t273 + t455;
t319 = t340 * t273;
t522 = -Icges(5,6) * t273 + qJD(4) * t245;
t102 = t260 * t522 - t261 * t319;
t320 = t248 * t273;
t520 = -Icges(5,5) * t273 + qJD(4) * t534;
t104 = t260 * t520 - t261 * t320;
t244 = Icges(5,5) * t278 - Icges(5,6) * t276;
t147 = Icges(5,3) * t261 - t244 * t260;
t149 = Icges(5,6) * t261 - t260 * t340;
t239 = Icges(5,4) * t446;
t445 = t260 * t278;
t151 = -Icges(5,1) * t445 + Icges(5,5) * t261 + t239;
t96 = t149 * t278 + t151 * t276;
t527 = qJD(4) * t96 + t102 * t276 - t104 * t278 - t147 * t273;
t101 = -t260 * t319 - t261 * t522;
t103 = -t260 * t320 - t261 * t520;
t148 = Icges(5,5) * t439 - Icges(5,6) * t440 + Icges(5,3) * t260;
t97 = t150 * t278 + t152 * t276;
t526 = qJD(4) * t97 + t101 * t276 - t103 * t278 - t148 * t273;
t525 = -Icges(6,3) * t273 + t204 * t272;
t524 = -Icges(6,5) * t273 + t272 * t535;
t243 = Icges(5,5) * t276 + Icges(5,6) * t278;
t523 = -Icges(5,3) * t273 + qJD(4) * t243;
t216 = t340 * qJD(4);
t217 = t248 * qJD(4);
t521 = qJD(4) * (t245 * t278 + t276 * t534) + t216 * t276 - t217 * t278 - t243 * t273;
t422 = -Icges(5,2) * t439 + t152 - t240;
t424 = t261 * t534 + t150;
t518 = t276 * t422 + t278 * t424;
t423 = Icges(5,2) * t445 + t151 + t239;
t425 = -t260 * t534 + t149;
t517 = -t276 * t423 - t278 * t425;
t358 = -t209 * t272 + t455;
t359 = t548 * t272;
t516 = -t204 * t273 + t265 * t359 + t267 * t358;
t318 = t339 * t273;
t364 = t140 * t272 - t260 * t318 - t261 * t529;
t321 = t273 * t209;
t366 = t138 * t272 + t260 * t321 + t261 * t524;
t515 = -t136 * t273 + t265 * t364 + t267 * t366;
t205 = Icges(6,5) * t267 - Icges(6,6) * t265;
t135 = Icges(6,3) * t261 - t205 * t260;
t365 = t139 * t272 + t260 * t529 - t261 * t318;
t367 = t137 * t272 - t260 * t524 + t261 * t321;
t514 = -t135 * t273 + t265 * t365 + t267 * t367;
t271 = t273 ^ 2;
t164 = t272 * t447;
t512 = -t164 / 0.2e1;
t165 = t273 * t190;
t511 = t165 / 0.2e1;
t510 = -t189 / 0.2e1;
t509 = t189 / 0.2e1;
t508 = -t190 / 0.2e1;
t507 = t190 / 0.2e1;
t506 = t260 / 0.2e1;
t505 = t261 / 0.2e1;
t504 = -t273 / 0.2e1;
t503 = t273 / 0.2e1;
t502 = -rSges(5,3) - pkin(7);
t501 = pkin(1) * t277;
t500 = pkin(1) * t279;
t499 = pkin(1) * qJD(1) ^ 2;
t496 = pkin(2) * t271;
t492 = t261 * pkin(7);
t76 = -t396 + t187 + (-t154 + t370) * t273;
t481 = t260 * t76;
t156 = t204 * t260;
t93 = -t261 * t335 + t156;
t480 = t93 * t273;
t377 = t261 * t491;
t108 = (t377 + t493) * t273 + t417;
t90 = -t273 * t327 + t389;
t479 = -t108 - t90;
t169 = t243 * t260;
t333 = t276 * t245 - t278 * t534;
t111 = -t261 * t333 + t169;
t466 = t111 * t273;
t463 = t137 * t265;
t462 = t139 * t267;
t461 = t149 * t276;
t460 = t151 * t278;
t456 = t205 * t273;
t452 = t243 * t261;
t451 = t244 * t273;
t175 = t251 * t260;
t450 = t251 * t261;
t444 = t261 * t127;
t438 = t266 * t273;
t434 = t261 * t135 + t137 * t449;
t433 = -t260 * t135 - t139 * t442;
t432 = t261 * t147 + t149 * t446;
t431 = t260 * t147 + t151 * t439;
t421 = -t261 * t435 - t262 * t447;
t415 = -rSges(4,1) * t447 - rSges(4,2) * t441;
t263 = t277 * t499;
t411 = t266 * t496 + t263;
t409 = qJD(4) * t261;
t406 = pkin(2) * t438;
t404 = qJD(4) ^ 2 * t494;
t403 = t279 * t499;
t402 = rSges(5,1) * t445;
t395 = pkin(4) * t408;
t381 = t261 * t407;
t382 = t261 * t408;
t387 = -rSges(5,1) * t382 - rSges(5,2) * t381 - t273 * t402;
t385 = t251 * t409;
t380 = -t447 / 0.2e1;
t379 = t441 / 0.2e1;
t378 = -pkin(3) - t488;
t374 = -t410 / 0.2e1;
t372 = -t409 / 0.2e1;
t371 = t409 / 0.2e1;
t212 = rSges(3,1) * t268 - t266 * rSges(3,2);
t368 = -t262 - t486;
t361 = -t148 - t460;
t356 = pkin(4) * t382;
t236 = pkin(3) * t447;
t355 = t236 - t387;
t354 = -pkin(3) - t368;
t351 = -t273 * (-pkin(3) * t260 + t492) + t406;
t184 = -rSges(3,1) * t437 + rSges(3,2) * t438;
t345 = -rSges(4,1) * t260 - rSges(4,2) * t261;
t81 = t138 * t267 + t140 * t265;
t337 = -t460 + t461;
t332 = t402 - t414;
t331 = t400 - t536;
t64 = -t137 * t443 - t433;
t67 = t261 * t148 + t150 * t446 - t152 * t445;
t330 = t255 + t349;
t329 = t539 - t536;
t167 = -t212 * t273 - t396;
t66 = -t151 * t445 + t432;
t325 = (t260 * t67 + t261 * t66) * qJD(4);
t68 = -t149 * t440 + t431;
t69 = t148 * t260 - t542;
t324 = (t260 * t69 + t261 * t68) * qJD(4);
t323 = -t268 * t496 - t403;
t322 = t414 + t492 - t498;
t315 = t397 + t406;
t313 = -t261 * t280 - t547;
t141 = t273 * t332;
t312 = t345 - t498;
t311 = -t271 * t192 + t323;
t308 = t156 * t190 - t189 * t458 + t456;
t306 = t224 - t394 + t537;
t304 = -t260 * t456 - t261 * t525 + t273 * t338;
t303 = -t261 * t456 + t525 * t260 + (-t462 + t463) * t273;
t302 = -t451 * t260 - t261 * t523 + t273 * t336;
t301 = t260 * t523 - t451 * t261 + t273 * t337;
t300 = t205 * t272 + t273 * t335;
t299 = t244 * qJD(4) + t273 * t333;
t298 = t242 - t401 - t495 - t497;
t297 = t141 + t351 + t385;
t296 = t356 - t388 - t421;
t13 = t303 * t260 - t261 * t514;
t14 = t304 * t260 - t261 * t515;
t15 = t260 * t514 + t303 * t261;
t16 = t260 * t515 + t304 * t261;
t62 = -t139 * t448 + t434;
t28 = t190 * t62 + t544;
t65 = t136 * t260 - t540;
t29 = t189 * t65 + t190 * t64 + t480;
t40 = -t265 * t367 + t267 * t365;
t41 = -t265 * t366 + t267 * t364;
t44 = t300 * t260 - t261 * t516;
t45 = t260 * t516 + t300 * t261;
t80 = t137 * t267 + t139 * t265;
t295 = (t13 * t190 + t14 * t189 - t164 * t64 + t165 * t65 + t273 * t44) * t506 + (t308 * t260 - t554 * t261) * t510 + (t554 * t260 + t308 * t261) * t508 + (t15 * t190 + t16 * t189 - t164 * t62 + t165 * t63 + t273 * t45) * t505 + (-t265 * t513 + t267 * t291) * t504 + t28 * t380 + t29 * t379 + ((t273 * t65 + t13) * t261 + (-t273 * t64 + t14) * t260) * t509 + (t260 * t63 + t261 * t62) * t512 + (t260 * t65 + t261 * t64) * t511 + ((t273 * t63 + t15) * t261 + (-t273 * t62 + t16) * t260) * t507 + ((t273 * t81 + t40) * t261 + (-t273 * t80 + t41) * t260) * t503;
t293 = t261 * t154 + t260 * t332;
t288 = t190 * t210 + t351 + t356 + (t331 + t546) * t273;
t287 = (-t136 - t462) * t260 + t540 + t434;
t105 = t273 * t414 + t387;
t106 = t273 * t328 + t386;
t286 = (-t106 - t142) * t260 + (t105 + t141) * t261;
t110 = t260 * t333 + t452;
t109 = t110 * t273;
t32 = t109 + t325;
t33 = t324 + t466;
t49 = -qJD(4) * t337 + t102 * t278 + t104 * t276;
t50 = -qJD(4) * t336 + t101 * t278 + t103 * t276;
t54 = t299 * t260 - t261 * t521;
t55 = t260 * t521 + t299 * t261;
t285 = ((t65 + t287) * t190 + t544) * t510 + (t109 + ((t432 + t69 + t542) * t261 + (-t68 + (t361 - t461) * t261 + t67 + t431) * t260) * qJD(4)) * t374 + (t80 + t92) * t512 + (t93 + t81) * t511 + (-t480 + (t63 + (-t136 + t463) * t261 - t338 * t260 + t433) * t190 + (-t62 + t287) * t189 + t29) * t508 + (t40 + t45) * t507 + (t33 - t466 + ((t67 + (-t148 + t461) * t261 - t431) * t261 + (t260 * t361 + t432 - t66) * t260) * qJD(4)) * t372 + (t49 + t55) * t371 + (t41 + t44 + t28) * t509 + (t50 + t54 + t32) * t410 / 0.2e1 + ((t110 + t96) * t374 + (t111 + t97) * t371 - qJD(4) * t333 + t216 * t278 + t217 * t276 - t265 * t358 + t267 * t359) * t273;
t107 = t236 + (-pkin(7) * t273 - t395) * t261 + t421;
t168 = pkin(7) * t441 - t236;
t36 = t260 * t404 + t165 * t210 + t182 * t189 + (-t107 - t168 - t89 + t356) * t273 + t411;
t37 = -t261 * t404 + t164 * t210 - t182 * t190 + (t231 - t479) * t273 + t311;
t56 = t288 + t397;
t284 = (t36 * (-rSges(6,3) + t280) + t37 * t368) * t260 + (t57 * t547 - t56 * (-t327 + t537)) * t273;
t60 = t225 * t410 + (-t105 - t168 + t385) * t273 + t411;
t61 = qJD(4) * t532 + t106 * t273 + t311;
t75 = t297 + t397;
t283 = (t61 * t378 + t60 * t502) * t260 + (t76 * (-t241 + t498) - t75 * (-t484 - t493 - t497) + (-t75 * t378 + t76 * t502) * t261) * t273;
t180 = t273 * t345;
t146 = t184 * t273 - t403;
t145 = t183 * t273 + t263;
t134 = -t396 + (-t191 - t497) * t273;
t133 = -t180 + t315;
t123 = t261 * t143;
t116 = t273 * (-rSges(4,1) * t441 + t233) + t323;
t115 = -t273 * t415 + t411;
t82 = qJD(4) * t293 + qJD(3);
t51 = qJD(3) + t190 * t143 + t189 * t331 + (t260 * t546 - t444) * qJD(4);
t46 = t286 * qJD(4);
t12 = -t164 * t143 + t190 * t89 + t165 * t331 - t189 * t90 + ((t273 * t539 + t107) * t261 + (-t108 + (t127 - t377) * t273) * t260) * qJD(4);
t1 = [t285 + m(3) * (t145 * (-t212 - t500) + t146 * (-t346 - t501) + (t167 - t184 + t396) * t166) + (t36 * (t306 - t500) + t57 * (t296 + t397) + t37 * (t313 - t501) + t284 + (-t57 + t533) * t56) * m(6) + (t60 * (t298 - t500) + t76 * (t355 + t397) + t61 * (t322 - t501) + t283 + (-t76 + t556) * t75) * m(5) + (t115 * (t330 - t500) + t134 * (t315 - t415) + t116 * (t312 - t501) + (-t134 + t555) * t133) * m(4); t285 + (t306 * t36 + t313 * t37 + t284 + (-t288 + t296) * t57 + t533 * t56) * m(6) + (t298 * t60 + t322 * t61 + t283 + (t355 - t297) * t76 + t556 * t75) * m(5) + (t115 * t330 + t116 * t312 + t555 * t133 + (t180 - t415) * t134) * m(4) + (-t145 * t212 - t146 * t346 - t166 * t184 + t167 * t183 - (t166 * t212 + t167 * t346) * t273) * m(3); m(5) * t46 + m(6) * t12; t295 + ((-t276 * t413 + t278 * t412) * t273 + ((t260 * t422 + t261 * t423) * t278 + (-t260 * t424 - t261 * t425) * t276) * qJD(4)) * t504 + ((t169 * t409 + t451) * t261 + (t551 + (t518 * t260 + (-t452 - t517) * t261) * qJD(4)) * t260) * t372 + ((t273 * t97 + t49) * t261 + (-t273 * t96 + t50) * t260) * t503 + ((-t410 * t452 + t451) * t260 + (-t551 + (t517 * t261 + (t169 - t518) * t260) * qJD(4)) * t261) * t374 + (t273 * t54 + ((t301 * t260 - t261 * t527 + t273 * t69) * t261 + (t302 * t260 - t261 * t526 - t273 * t68) * t260) * t553) * t506 + (t273 * t55 + ((t260 * t527 + t301 * t261 + t273 * t67) * t261 + (t260 * t526 + t302 * t261 - t273 * t66) * t260) * t553) * t505 + (t325 + t32) * t380 + (t324 + t33) * t379 + (t12 * (-t444 + t123 + (t260 * t354 + t329) * t260) + (t260 * t36 - t261 * t37) * (pkin(4) * t276 + t210) + (t261 * t107 + t479 * t260 + (t329 * t261 + (t261 * t354 + t430) * t260) * t273 - (-t260 ^ 2 - t261 ^ 2) * t395 + t552) * t51 + t545 + (-(-pkin(4) * t407 - t182) * t261 - pkin(4) * t381 + t550) * t56) * m(6) + (t251 * t441 * t76 + t175 * t60 + t225 * t481 + t286 * t82 + t293 * t46 - t450 * t61 - t532 * t75 - (-t175 * t75 + t450 * t76) * t273 - (t82 * (-t175 * t260 - t261 * t450) + (t261 * t75 + t481) * t344) * qJD(4)) * m(5); t295 + (t12 * (t260 * t331 + t123) + t36 * t162 - t37 * t454 + (-t260 * t90 + (-t260 * t143 + t261 * t331) * t273 + t552) * t51 + t545 + (t261 * t182 + t550) * t56) * m(6);];
tauc = t1(:);
