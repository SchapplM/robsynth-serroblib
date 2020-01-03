% Calculate vector of centrifugal and Coriolis load on the joints for
% S4RRRP3
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
% Datum: 2019-12-31 17:14
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauc = S4RRRP3_coriolisvecJ_fixb_slag_vp1(qJ, qJD, ...
  pkin, m, rSges, Icges)
%% Coder Information
%#codegen
%$cgargs {zeros(4,1),zeros(4,1),zeros(6,1),zeros(5,1),zeros(5,3),zeros(5,6)}
assert(isreal(qJ) && all(size(qJ) == [4 1]), ...
  'S4RRRP3_coriolisvecJ_fixb_slag_vp1: qJ has to be [4x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [4 1]), ...
  'S4RRRP3_coriolisvecJ_fixb_slag_vp1: qJD has to be [4x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [6 1]), ...
  'S4RRRP3_coriolisvecJ_fixb_slag_vp1: pkin has to be [6x1] (double)');
assert(isreal(m) && all(size(m) == [5 1]), ...
  'S4RRRP3_coriolisvecJ_fixb_slag_vp1: m has to be [5x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [5,3]), ...
  'S4RRRP3_coriolisvecJ_fixb_slag_vp1: rSges has to be [5x3] (double)');
assert(isreal(Icges) && all(size(Icges) == [5 6]), ...
  'S4RRRP3_coriolisvecJ_fixb_slag_vp1: Icges has to be [5x6] (double)'); 

%% Symbolic Calculation
% From coriolisvec_joint_fixb_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 17:14:01
% EndTime: 2019-12-31 17:14:13
% DurationCPUTime: 9.71s
% Computational Cost: add. (7957->413), mult. (8975->513), div. (0->0), fcn. (6997->6), ass. (0->244)
t242 = cos(qJ(3));
t240 = sin(qJ(3));
t389 = Icges(4,4) * t240;
t190 = Icges(4,2) * t242 + t389;
t235 = Icges(5,5) * t240;
t298 = Icges(5,3) * t242 - t235;
t477 = t190 + t298;
t388 = Icges(5,5) * t242;
t192 = Icges(5,1) * t240 - t388;
t236 = Icges(4,4) * t242;
t194 = Icges(4,1) * t240 + t236;
t482 = t192 + t194;
t239 = qJ(1) + qJ(2);
t233 = sin(t239);
t187 = Icges(4,5) * t242 - Icges(4,6) * t240;
t234 = cos(t239);
t276 = t187 * t234;
t118 = Icges(4,3) * t233 + t276;
t189 = Icges(5,4) * t242 + Icges(5,6) * t240;
t277 = t189 * t234;
t120 = Icges(5,2) * t233 + t277;
t481 = t118 + t120;
t300 = Icges(5,1) * t242 + t235;
t123 = -Icges(5,4) * t234 + t233 * t300;
t375 = t233 * t240;
t213 = Icges(4,4) * t375;
t374 = t233 * t242;
t125 = Icges(4,1) * t374 - Icges(4,5) * t234 - t213;
t480 = t123 + t125;
t279 = t300 * t234;
t124 = Icges(5,4) * t233 + t279;
t195 = Icges(4,1) * t242 - t389;
t280 = t195 * t234;
t126 = Icges(4,5) * t233 + t280;
t479 = t124 + t126;
t185 = Icges(5,3) * t240 + t388;
t299 = -Icges(4,2) * t240 + t236;
t478 = t185 - t299;
t475 = t195 + t300;
t465 = -t477 * t240 + t482 * t242;
t369 = t234 * t242;
t212 = Icges(5,5) * t369;
t370 = t234 * t240;
t116 = Icges(5,6) * t233 + Icges(5,3) * t370 + t212;
t474 = t116 * t370 + t481 * t233 + t479 * t369;
t119 = -Icges(5,2) * t234 + t189 * t233;
t111 = t233 * t119;
t115 = -Icges(5,6) * t234 + t185 * t233;
t117 = Icges(4,5) * t374 - Icges(4,6) * t375 - Icges(4,3) * t234;
t473 = -t115 * t370 - t233 * t117 - t480 * t369 - t111;
t471 = t478 * qJD(3);
t470 = t475 * qJD(3);
t186 = Icges(4,5) * t240 + Icges(4,6) * t242;
t188 = Icges(5,4) * t240 - Icges(5,6) * t242;
t469 = t186 + t188;
t468 = -t187 - t189;
t238 = qJD(1) + qJD(2);
t467 = (-Icges(4,6) + Icges(5,6)) * t238 + t477 * qJD(3);
t466 = (Icges(5,4) + Icges(4,5)) * t238 - t482 * qJD(3);
t121 = Icges(4,4) * t374 - Icges(4,2) * t375 - Icges(4,6) * t234;
t385 = t121 * t240;
t293 = -t125 * t242 + t385;
t372 = t234 * t119;
t296 = t115 * t240 + t123 * t242;
t421 = t233 * t296;
t49 = -t372 + t421;
t464 = -t234 * t117 - t233 * t293 + t49;
t437 = -t121 * t370 - t473;
t278 = t299 * t234;
t122 = Icges(4,6) * t233 + t278;
t436 = -t122 * t370 + t474;
t378 = t188 * t234;
t381 = t186 * t234;
t463 = t465 * t233 - t378 - t381;
t379 = t188 * t233;
t382 = t186 * t233;
t462 = t465 * t234 + t379 + t382;
t307 = -t116 * t375 + t120 * t234 - t124 * t374;
t100 = t126 * t374;
t312 = t234 * t118 - t100;
t52 = -t122 * t375 - t312;
t461 = -t307 + t52;
t205 = rSges(5,1) * t242 + rSges(5,3) * t240;
t460 = pkin(3) * t242 + qJ(4) * t240 + t205;
t444 = rSges(5,1) + pkin(3);
t376 = t233 * t238;
t459 = t467 * t234 - t478 * t376;
t371 = t234 * t238;
t458 = t185 * t371 + t467 * t233 - t238 * t278;
t457 = t466 * t234 - t475 * t376;
t456 = (-t279 - t280) * t238 - t466 * t233;
t435 = rSges(5,3) + qJ(4);
t455 = (t116 - t122) * t242 - t479 * t240;
t433 = (-t115 + t121) * t242 + t480 * t240;
t454 = t470 * t242 + t471 * t240 + t469 * t238 + (-t240 * t482 - t477 * t242) * qJD(3);
t453 = (Icges(5,2) + Icges(4,3)) * t238 - t469 * qJD(3);
t384 = t122 * t240;
t452 = -t116 * t240 - t242 * t479 + t384;
t451 = t293 - t296;
t450 = t468 * qJD(3) + t465 * t238;
t449 = t462 * t238;
t228 = t234 * rSges(5,2);
t357 = t233 * t460 - t228;
t448 = t357 * t238;
t447 = (t233 * t436 - t234 * t437) * qJD(3);
t446 = (t461 * t233 - t234 * t464) * qJD(3);
t445 = t463 * t238;
t443 = t445 + t446;
t442 = t447 + t449;
t441 = t451 * qJD(3) + t456 * t240 + t458 * t242;
t440 = -t452 * qJD(3) + t457 * t240 - t459 * t242;
t439 = -t450 * t233 + t454 * t234;
t438 = t454 * t233 + t450 * t234;
t345 = t444 * t240 - t435 * t242;
t313 = t234 * t345;
t231 = t234 * pkin(6);
t165 = pkin(2) * t233 - t231;
t199 = pkin(6) * t371;
t432 = t238 * t165 + t199;
t431 = t233 * rSges(5,2) + pkin(3) * t369;
t430 = t455 * qJD(3) + t481 * t238 + t459 * t240 + t457 * t242;
t429 = t456 * t242 - t458 * t240 + (-t117 - t119) * t238 + t433 * qJD(3);
t428 = t372 + t474;
t338 = qJD(3) * t242;
t325 = t233 * t338;
t368 = t238 * t240;
t427 = t234 * t368 + t325;
t241 = sin(qJ(1));
t397 = pkin(1) * qJD(1);
t332 = t241 * t397;
t163 = rSges(3,1) * t233 + rSges(3,2) * t234;
t383 = t163 * t238;
t132 = -t332 - t383;
t426 = (t276 + t277 + t451) * t238 + t453 * t233;
t425 = t453 * t234 + t452 * t238 + t468 * t376;
t424 = 0.2e1 * qJD(3);
t423 = t233 / 0.2e1;
t336 = qJD(4) * t242;
t356 = rSges(5,1) * t369 + t370 * t435 + t431;
t46 = -t336 + (t357 * t233 + t356 * t234) * qJD(3);
t422 = qJD(3) * t46;
t166 = t234 * pkin(2) + t233 * pkin(6);
t420 = t166 + t356;
t339 = qJD(3) * t240;
t326 = t233 * t339;
t351 = t444 * t326;
t419 = rSges(4,1) * t369 + t233 * rSges(4,3);
t337 = qJD(4) * t240;
t208 = t234 * t337;
t323 = t234 * t338;
t418 = rSges(5,2) * t371 + t435 * t323 + t208;
t343 = rSges(4,2) * t375 + t234 * rSges(4,3);
t128 = rSges(4,1) * t374 - t343;
t114 = t238 * t128;
t331 = t233 * t368;
t281 = rSges(4,3) * t371 + (-t323 + t331) * rSges(4,2);
t324 = t234 * t339;
t417 = -rSges(4,1) * t324 + t114 + t281 + t432;
t202 = rSges(4,1) * t240 + rSges(4,2) * t242;
t341 = qJD(3) * t233;
t162 = t202 * t341;
t130 = -rSges(4,2) * t370 + t419;
t272 = t130 + t166;
t412 = t238 * t272 - t162;
t409 = t418 + t432 + t448;
t359 = -Icges(4,2) * t374 + t125 - t213;
t363 = t194 * t233 + t121;
t408 = -t240 * t359 - t242 * t363;
t361 = -t298 * t233 + t123;
t365 = -t192 * t233 + t115;
t407 = -t240 * t361 + t242 * t365;
t403 = t238 / 0.2e1;
t401 = pkin(1) * t241;
t243 = cos(qJ(1));
t237 = t243 * pkin(1);
t273 = -t238 * t374 - t324;
t400 = t273 * t444 - t435 * t331 + t418;
t322 = t233 * t337;
t399 = rSges(5,3) * t325 + t322 + t427 * qJ(4) - t351 + (t205 * t234 + t431) * t238;
t398 = rSges(4,1) * t242;
t340 = qJD(3) * t234;
t327 = t202 * t340;
t274 = -t327 - t332;
t61 = (-t128 - t165) * t238 + t274;
t396 = t234 * t61;
t287 = -qJD(3) * t313 + t208;
t268 = t287 - t332;
t47 = (-t165 - t357) * t238 + t268;
t395 = t238 * t47;
t380 = t187 * t238;
t377 = t189 * t238;
t364 = -Icges(5,1) * t370 + t116 + t212;
t362 = -t194 * t234 - t122;
t360 = -t298 * t234 + t124;
t358 = -t190 * t234 + t126;
t355 = t345 * t233;
t353 = -qJD(3) * t460 + t336;
t349 = -t298 + t300;
t348 = t185 - t192;
t347 = -t190 + t195;
t346 = t194 + t299;
t342 = t228 + t231;
t244 = qJD(1) ^ 2;
t335 = t244 * t401;
t334 = t244 * t237;
t333 = t243 * t397;
t328 = rSges(4,1) * t326 + rSges(4,2) * t427;
t319 = -pkin(2) - t398;
t318 = -t341 / 0.2e1;
t315 = t340 / 0.2e1;
t164 = t234 * rSges(3,1) - rSges(3,2) * t233;
t311 = -t117 + t384;
t308 = t238 * (-pkin(2) * t376 + t199) - t335;
t136 = rSges(3,1) * t371 - rSges(3,2) * t376;
t306 = t336 + t353;
t304 = -rSges(4,2) * t240 + t398;
t62 = t333 + t412;
t303 = -t233 * t62 - t396;
t288 = -t341 * t345 + t322;
t65 = (t128 * t233 + t130 * t234) * qJD(3);
t271 = -t240 * t360 + t242 * t364;
t270 = -t240 * t358 + t242 * t362;
t269 = t233 * t319 + t231 + t343;
t267 = (t240 * t348 + t242 * t349) * t238;
t266 = (-t240 * t346 + t242 * t347) * t238;
t247 = (t319 * t396 + (t61 * (-rSges(4,3) - pkin(6)) + t62 * t319) * t233) * t238;
t246 = (((t52 - t100 + (t118 + t385) * t234 + t473) * t234 + (t49 - t421 + t428) * t233) * qJD(3) + t449) * t315 + (t465 * qJD(3) + t470 * t240 - t471 * t242) * t238 + (t439 + t440) * t341 / 0.2e1 + (((t234 * t311 - t428 + t436) * t234 + (t233 * t311 - t111 + t307 + t312 + t437) * t233) * qJD(3) + t443 - t445) * t318 - (t438 - t441 + t442) * t340 / 0.2e1 + ((t433 + t463) * t233 + (-t455 + t462) * t234) * qJD(3) * t403;
t137 = t166 * t238;
t23 = -t334 + t306 * t340 + (-t137 + (qJD(3) * t345 - t337) * t233 - t399) * t238;
t48 = t238 * t420 + t288 + t333;
t245 = (-t48 * t444 * t339 + (-t240 * t435 - t242 * t444 - pkin(2)) * t395) * t234 + (-t23 * pkin(2) + (-t47 * qJD(4) - t23 * t435) * t240 + (-qJD(3) * t435 * t47 - t23 * t444) * t242 + (t47 * (-rSges(5,2) - pkin(6)) + t48 * (-pkin(2) - t460)) * t238) * t233;
t176 = t304 * qJD(3);
t156 = t202 * t234;
t152 = t202 * t233;
t133 = t164 * t238 + t333;
t108 = -t136 * t238 - t334;
t107 = -t238 * t383 - t335;
t89 = t238 * t419 - t328;
t87 = rSges(4,1) * t273 + t281;
t44 = -t334 - t176 * t340 + (-t137 - t89 + t162) * t238;
t43 = t238 * t87 + (-t176 * t233 - t202 * t371) * qJD(3) + t308;
t22 = (t208 + t400) * t238 + (t233 * t306 - t238 * t313) * qJD(3) + t308;
t5 = (t337 + (t400 + t448) * t234 + (-t356 * t238 + t399) * t233) * qJD(3);
t1 = [m(3) * (t108 * (-t163 - t401) + t107 * (t164 + t237) + (-t136 - t333 + t133) * t132) + t246 + (t23 * (t342 - t401) + t47 * (-t333 + t351) + t22 * (t237 + t420) + t245 + (-t332 + t47 - t268 + t409) * t48) * m(5) + (t44 * (t269 - t401) + t61 * (t328 - t333) + t43 * (t237 + t272) + t247 + (-t274 + t61 - t332 + t417) * t62) * m(4); t246 + (t23 * t342 + t245 + (-t287 + t409) * t48 + (t351 + t288) * t47 + (t22 + t395) * t420) * m(5) + (t269 * t44 + t272 * t43 + t247 + (t327 + t417) * t62 + (t328 + t412) * t61) * m(4) + (t107 * t164 - t108 * t163 - t132 * t136 - t133 * t383 - (-t132 * t164 - t133 * t163) * t238) * m(3); -(((t346 - t348) * t242 + (t347 + t349) * t240) * t238 + (((-t359 - t361) * t234 + (t358 + t360) * t233) * t242 + ((t363 - t365) * t234 + (t362 + t364) * t233) * t240) * qJD(3)) * t238 / 0.2e1 + ((-t238 * t455 + t441) * t234 + (t238 * t433 + t440) * t233) * t403 + ((-t341 * t378 + t377) * t233 + (t267 + (-t407 * t234 + (t379 + t271) * t233) * qJD(3)) * t234 + (-t341 * t381 + t380) * t233 + (t266 + (-t408 * t234 + (t382 + t270) * t233) * qJD(3)) * t234) * t318 + ((-t340 * t379 - t377) * t234 + (t267 + (t271 * t233 + (t378 - t407) * t234) * qJD(3)) * t233 + (-t340 * t382 - t380) * t234 + (t266 + (t270 * t233 + (t381 - t408) * t234) * qJD(3)) * t233) * t315 + (-(t240 * t46 + (t233 * t48 + t234 * t47) * t242) * qJD(4) - (-t313 * t48 + t355 * t47) * t238 - ((-t313 * t46 - t460 * t47) * t234 + (-t355 * t46 - t460 * t48) * t233) * qJD(3) + (-t23 * t345 + t47 * t353 + t5 * t356 + t46 * t400 + (-t345 * t48 + t357 * t46) * t238) * t234 + (-t22 * t345 + t48 * t353 + t5 * t357 + t46 * t399 + (t345 * t47 - t356 * t46) * t238) * t233) * m(5) + (0.2e1 * t65 * ((t87 + t114) * t234 + (-t130 * t238 + t89) * t233) + t303 * t176 + ((-t238 * t62 - t44) * t234 + (t238 * t61 - t43) * t233) * t202 - (t152 * t61 - t156 * t62) * t238 - (t65 * (-t152 * t233 - t156 * t234) + t303 * t304) * qJD(3)) * m(4) + (t439 * t238 + ((t429 * t234 + t436 * t238) * t234 + (t425 * t233 + t437 * t238 + (-t426 + t430) * t234) * t233) * t424) * t423 - (t438 * t238 + ((t426 * t234 + t461 * t238) * t234 + (t430 * t233 + t464 * t238 + (-t425 + t429) * t234) * t233) * t424) * t234 / 0.2e1 + (t443 + t446) * t376 / 0.2e1 + (t442 + t447) * t371 / 0.2e1; (-t242 * t5 + 0.2e1 * (t422 / 0.2e1 + t22 * t423 + t23 * t234 / 0.2e1 - (t233 ^ 2 + t234 ^ 2) * t422 / 0.2e1) * t240) * m(5);];
tauc = t1(:);
