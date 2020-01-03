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
% Datum: 2020-01-03 12:02
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
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
% StartTime: 2020-01-03 12:01:48
% EndTime: 2020-01-03 12:02:00
% DurationCPUTime: 8.14s
% Computational Cost: add. (17736->561), mult. (11158->740), div. (0->0), fcn. (8628->10), ass. (0->351)
t272 = qJ(1) + qJ(2);
t259 = pkin(9) + t272;
t252 = sin(t259);
t253 = cos(t259);
t270 = qJD(1) + qJD(2);
t275 = cos(qJ(4));
t264 = Icges(5,4) * t275;
t273 = sin(qJ(4));
t330 = -Icges(5,2) * t273 + t264;
t311 = t330 * t270;
t454 = Icges(5,4) * t273;
t233 = Icges(5,2) * t275 + t454;
t489 = -Icges(5,6) * t270 + qJD(4) * t233;
t102 = -t252 * t489 + t253 * t311;
t236 = Icges(5,1) * t275 - t454;
t313 = t236 * t270;
t498 = Icges(5,1) * t273 + t264;
t486 = -Icges(5,5) * t270 + qJD(4) * t498;
t104 = -t252 * t486 + t253 * t313;
t422 = t253 * t275;
t423 = t253 * t273;
t146 = Icges(5,4) * t422 - Icges(5,2) * t423 + Icges(5,6) * t252;
t228 = Icges(5,4) * t423;
t148 = Icges(5,1) * t422 + Icges(5,5) * t252 - t228;
t326 = t146 * t275 + t148 * t273;
t147 = -Icges(5,5) * t253 + t236 * t252;
t436 = t147 * t275;
t145 = -Icges(5,6) * t253 + t252 * t330;
t438 = t145 * t273;
t327 = t436 - t438;
t513 = -qJD(4) * t327 - t102 * t275 - t104 * t273 + t270 * t326;
t186 = t252 * rSges(4,1) + t253 * rSges(4,2);
t261 = sin(t272);
t254 = pkin(2) * t261;
t499 = t254 + t186;
t263 = cos(t272);
t500 = t261 * rSges(3,1) + t263 * rSges(3,2);
t180 = t500 * t270;
t274 = sin(qJ(1));
t461 = pkin(1) * qJD(1);
t378 = t274 * t461;
t164 = t378 + t180;
t271 = qJ(4) + qJ(5);
t260 = sin(t271);
t262 = cos(t271);
t453 = Icges(6,4) * t260;
t200 = Icges(6,1) * t262 - t453;
t312 = t200 * t252;
t134 = -Icges(6,5) * t253 + t312;
t426 = t253 * t260;
t210 = Icges(6,4) * t426;
t425 = t253 * t262;
t135 = Icges(6,1) * t425 + Icges(6,5) * t252 - t210;
t269 = qJD(4) + qJD(5);
t184 = t252 * t269;
t185 = t253 * t269;
t197 = Icges(6,2) * t262 + t453;
t249 = Icges(6,4) * t262;
t329 = -Icges(6,2) * t260 + t249;
t501 = Icges(6,1) * t260 + t249;
t508 = t501 + t329;
t286 = t184 * (-Icges(6,2) * t425 + t135 - t210) - t185 * (-t197 * t252 + t134) + t270 * t508;
t310 = t329 * t252;
t132 = -Icges(6,6) * t253 + t310;
t133 = Icges(6,4) * t425 - Icges(6,2) * t426 + Icges(6,6) * t252;
t481 = t184 * (t253 * t501 + t133) - t185 * (t252 * t501 + t132) + t270 * (t197 - t200);
t512 = t286 * t260 + t262 * t481;
t394 = t498 + t330;
t395 = t233 - t236;
t510 = (t273 * t394 + t275 * t395) * t270;
t468 = pkin(3) * t253;
t189 = pkin(7) * t252 + t468;
t417 = t263 * t270;
t239 = pkin(2) * t417;
t509 = -t270 * t189 + t239;
t507 = 0.2e1 * qJD(4);
t196 = Icges(6,5) * t262 - Icges(6,6) * t260;
t130 = -Icges(6,3) * t253 + t196 * t252;
t440 = t133 * t260;
t328 = -t135 * t262 + t440;
t318 = -t130 + t328;
t506 = t185 * t318;
t505 = t270 * t499;
t464 = rSges(4,2) * t252;
t187 = t253 * rSges(4,1) - t464;
t178 = t270 * t187;
t424 = t253 * t270;
t222 = rSges(4,1) * t424;
t503 = -t178 + t222;
t224 = pkin(7) * t424;
t416 = t270 * t273;
t381 = rSges(5,2) * t416;
t402 = -rSges(5,3) * t424 - t252 * t381;
t502 = t224 - t402;
t230 = rSges(5,2) * t423;
t384 = rSges(5,1) * t422;
t151 = rSges(5,3) * t252 - t230 + t384;
t136 = t270 * t151;
t240 = rSges(5,1) * t273 + rSges(5,2) * t275;
t391 = qJD(4) * t252;
t338 = -t240 * t391 + t239;
t429 = t252 * t270;
t397 = pkin(3) * t424 + pkin(7) * t429;
t401 = rSges(5,3) * t429 + t270 * t384;
t497 = t397 + t401 - t136 - t338 + t509;
t467 = pkin(4) * t275;
t256 = pkin(3) + t467;
t212 = t253 * t256;
t277 = -pkin(8) - pkin(7);
t126 = t468 - t212 + (pkin(7) + t277) * t252;
t120 = t270 * t126;
t214 = rSges(6,2) * t426;
t383 = rSges(6,1) * t425;
t138 = rSges(6,3) * t252 - t214 + t383;
t127 = t270 * t138;
t190 = t256 * t424;
t201 = rSges(6,1) * t260 + rSges(6,2) * t262;
t389 = qJD(4) * t273;
t370 = t252 * t389;
t346 = pkin(4) * t370;
t301 = -t184 * t201 + t239 - t346;
t403 = rSges(6,3) * t429 + t270 * t383;
t496 = t190 + t403 + t120 - t127 - t301 + t509;
t433 = t197 * t269;
t495 = -Icges(6,6) * t270 + t433;
t159 = t201 * t252;
t160 = t201 * t253;
t462 = rSges(6,2) * t260;
t465 = rSges(6,1) * t262;
t203 = -t462 + t465;
t248 = t252 * pkin(3);
t188 = -pkin(7) * t253 + t248;
t244 = t253 * t277;
t398 = t252 * t256 + t244;
t125 = -t188 + t398;
t430 = t252 * t262;
t213 = rSges(6,1) * t430;
t431 = t252 * t260;
t137 = -rSges(6,2) * t431 - rSges(6,3) * t253 + t213;
t51 = t137 * t184 + t138 * t185 + qJD(3) + (t125 * t252 - t126 * t253) * qJD(4);
t368 = t253 * t389;
t345 = pkin(4) * t368;
t319 = -t185 * t201 - t345;
t357 = t188 + t254;
t320 = t125 + t137 + t357;
t56 = t270 * t320 - t319 + t378;
t276 = cos(qJ(1));
t258 = t276 * t461;
t413 = -t126 + t138;
t57 = t258 + (t189 + t413) * t270 + t301;
t494 = -(-t270 * t159 + t185 * t203) * t56 - t51 * (-t159 * t184 - t185 * t160) - t57 * (-t160 * t270 - t184 * t203);
t232 = Icges(5,5) * t275 - Icges(5,6) * t273;
t143 = -Icges(5,3) * t253 + t232 * t252;
t96 = t145 * t275 + t147 * t273;
t493 = qJD(4) * t96 + t102 * t273 - t104 * t275 - t143 * t270;
t195 = Icges(6,5) * t260 + Icges(6,6) * t262;
t492 = -Icges(6,3) * t270 + t195 * t269;
t491 = -Icges(6,5) * t270 + t269 * t501;
t231 = Icges(5,5) * t273 + Icges(5,6) * t275;
t490 = -Icges(5,3) * t270 + qJD(4) * t231;
t206 = t330 * qJD(4);
t207 = t236 * qJD(4);
t488 = qJD(4) * (t233 * t275 + t273 * t498) + t206 * t273 - t207 * t275 - t231 * t270;
t101 = t252 * t311 + t253 * t489;
t103 = t252 * t313 + t253 * t486;
t144 = Icges(5,5) * t422 - Icges(5,6) * t423 + Icges(5,3) * t252;
t487 = qJD(4) * t326 - t101 * t273 + t103 * t275 - t144 * t270;
t347 = -t200 * t269 + t433;
t348 = t508 * t269;
t484 = -t195 * t270 + t260 * t348 + t262 * t347;
t131 = Icges(6,5) * t425 - Icges(6,6) * t426 + Icges(6,3) * t252;
t352 = t135 * t269 - t253 * t495 - t270 * t310;
t354 = t133 * t269 + t253 * t491 + t270 * t312;
t483 = -t131 * t270 + t260 * t352 + t262 * t354;
t353 = t134 * t269 - t252 * t495 + t329 * t424;
t355 = t132 * t269 - t200 * t424 + t252 * t491;
t482 = -t130 * t270 + t260 * t353 + t262 * t355;
t268 = t270 ^ 2;
t162 = t270 * t184;
t480 = t162 / 0.2e1;
t163 = t270 * t185;
t479 = -t163 / 0.2e1;
t478 = -t184 / 0.2e1;
t477 = t184 / 0.2e1;
t476 = t185 / 0.2e1;
t475 = -t185 / 0.2e1;
t474 = -t252 / 0.2e1;
t473 = -t253 / 0.2e1;
t472 = -t270 / 0.2e1;
t471 = t270 / 0.2e1;
t470 = rSges(5,3) + pkin(7);
t469 = pkin(2) * t268;
t266 = t274 * pkin(1);
t267 = t276 * pkin(1);
t466 = rSges(5,1) * t275;
t463 = rSges(5,2) * t273;
t460 = t270 * t57;
t76 = t258 + (t151 + t189) * t270 + t338;
t459 = t270 * t76;
t435 = t195 * t252;
t93 = -t197 * t426 + t425 * t501 + t435;
t458 = t93 * t270;
t107 = t345 + t224 + (t244 + (-pkin(3) + t256) * t252) * t270;
t418 = t262 * t269;
t379 = rSges(6,2) * t418;
t420 = t260 * t270;
t380 = rSges(6,2) * t420;
t404 = -rSges(6,3) * t424 - t252 * t380;
t421 = t260 * t269;
t89 = t253 * t379 + (t253 * t421 + t262 * t429) * rSges(6,1) + t404;
t457 = -t107 - t89;
t432 = t231 * t252;
t111 = -t233 * t423 + t422 * t498 + t432;
t444 = t111 * t270;
t441 = t132 * t260;
t439 = t134 * t262;
t437 = t146 * t273;
t154 = t195 * t253;
t308 = t196 * t270;
t168 = t231 * t253;
t309 = t232 * t270;
t428 = t252 * t273;
t427 = t252 * t275;
t419 = t261 * t270;
t415 = t270 * t277;
t414 = t275 * t148;
t408 = t252 * t498 + t145;
t407 = t253 * t498 + t146;
t406 = -t233 * t252 + t147;
t405 = -Icges(5,2) * t422 + t148 - t228;
t255 = pkin(2) * t263;
t396 = -t230 + t255;
t279 = qJD(1) ^ 2;
t257 = t279 * t267;
t393 = t263 * t469 + t257;
t390 = qJD(4) * t253;
t388 = qJD(4) * t275;
t387 = pkin(2) * t419;
t386 = qJD(4) ^ 2 * t467;
t385 = t279 * t266;
t382 = rSges(6,1) * t421;
t377 = pkin(4) * t389;
t90 = -t252 * t382 + (-t252 * t418 - t253 * t420) * rSges(6,2) + t403;
t376 = t137 * t424 + (-t127 + t90) * t252;
t375 = t270 * t397 + t393;
t374 = t212 - t214 + t255;
t229 = rSges(5,1) * t427;
t373 = t229 + t248 + t254;
t371 = t240 * t390;
t369 = t252 * t388;
t366 = t429 / 0.2e1;
t365 = -t424 / 0.2e1;
t364 = pkin(3) + t466;
t362 = t391 / 0.2e1;
t360 = t390 / 0.2e1;
t356 = pkin(4) * t273 + t201;
t351 = -t131 - t441;
t350 = -t131 + t439;
t204 = rSges(3,1) * t263 - rSges(3,2) * t261;
t165 = t204 * t270 + t258;
t343 = t213 + t254 + t398;
t150 = -rSges(5,2) * t428 - rSges(5,3) * t253 + t229;
t340 = t150 + t357;
t181 = rSges(3,1) * t417 - rSges(3,2) * t419;
t339 = t187 + t255;
t335 = -t463 + t466;
t75 = t270 * t340 + t371 + t378;
t334 = -t252 * t76 + t253 * t75;
t81 = -t133 * t262 - t135 * t260;
t325 = -t414 + t437;
t324 = t150 * t252 + t151 * t253;
t323 = -t197 * t260 + t262 * t501;
t321 = -t233 * t273 + t275 * t498;
t317 = qJD(4) * t240;
t122 = t147 * t427;
t66 = -t143 * t253 - t145 * t428 + t122;
t123 = t252 * t414;
t67 = t144 * t253 + t146 * t428 - t123;
t316 = (-t252 * t67 - t253 * t66) * qJD(4);
t124 = t145 * t423;
t68 = -t143 * t252 - t147 * t422 + t124;
t69 = t252 * t144 - t253 * t325;
t315 = (-t252 * t69 - t253 * t68) * qJD(4);
t314 = -t261 * t469 - t385;
t307 = -t378 - t387;
t303 = t154 * t184 - t185 * t435 - t308;
t300 = -t252 * t308 - t253 * t492 + t270 * t328;
t299 = -t253 * t308 + t492 * t252 + (t439 - t441) * t270;
t298 = -t252 * t309 - t253 * t490 + t270 * t325;
t297 = t252 * t490 - t253 * t309 + t270 * t327;
t296 = -t196 * t269 + t270 * t323;
t295 = -t232 * qJD(4) + t270 * t321;
t294 = t273 * t406 + t275 * t408;
t293 = t273 * t405 + t275 * t407;
t13 = t299 * t252 + t253 * t482;
t14 = t300 * t252 - t253 * t483;
t15 = -t252 * t482 + t299 * t253;
t16 = t252 * t483 + t300 * t253;
t112 = t134 * t430;
t62 = -t130 * t253 - t132 * t431 + t112;
t113 = t135 * t430;
t63 = t131 * t253 + t133 * t431 - t113;
t92 = t252 * t323 - t154;
t91 = t92 * t270;
t28 = -t184 * t63 - t185 * t62 + t91;
t114 = t132 * t426;
t64 = -t130 * t252 - t134 * t425 + t114;
t65 = t131 * t252 - t253 * t328;
t29 = -t184 * t65 - t185 * t64 - t458;
t40 = -t260 * t355 + t262 * t353;
t41 = t260 * t354 - t262 * t352;
t44 = t296 * t252 + t253 * t484;
t45 = -t252 * t484 + t296 * t253;
t80 = t132 * t262 + t134 * t260;
t292 = (-t13 * t185 - t14 * t184 + t162 * t64 - t163 * t65 + t270 * t44) * t474 + (t303 * t252 + t253 * t512) * t477 + (-t252 * t512 + t303 * t253) * t476 + (-t15 * t185 - t16 * t184 + t162 * t62 - t163 * t63 + t270 * t45) * t473 + (-t260 * t481 + t262 * t286) * t472 + t28 * t366 + t29 * t365 + ((-t270 * t65 - t13) * t253 + (t270 * t64 - t14) * t252) * t478 + (-t252 * t63 - t253 * t62) * t480 + (-t252 * t65 - t253 * t64) * t479 + ((-t270 * t63 - t15) * t253 + (t270 * t62 - t16) * t252) * t475 + ((-t270 * t81 - t40) * t253 + (t270 * t80 - t41) * t252) * t471;
t291 = -t377 - t379 - t382;
t105 = rSges(5,2) * t253 * t388 + (t270 * t427 + t368) * rSges(5,1) + t402;
t106 = -rSges(5,1) * t370 + (-t253 * t416 - t369) * rSges(5,2) + t401;
t288 = (t150 * t270 - t105) * t253 + (t106 - t136) * t252;
t128 = t378 + t505;
t129 = t258 + t239 + t178;
t283 = (-t128 * t464 - t129 * t499) * t270;
t110 = t252 * t321 - t168;
t109 = t110 * t270;
t32 = t109 + t316;
t33 = t315 - t444;
t50 = qJD(4) * t325 + t101 * t275 + t103 * t273;
t54 = t295 * t252 + t253 * t488;
t55 = -t252 * t488 + t295 * t253;
t282 = t163 * t93 / 0.2e1 + (t109 + ((t68 + t123 - t124 + (t143 - t437) * t252) * t252 + (-t122 - t69 + (t143 - t325) * t253 + (t436 + t438) * t252) * t253) * qJD(4)) * t362 + (t91 - (t252 * t351 + t112 + t65) * t185 + (t113 - t114 + t64 + (t130 - t440) * t252) * t184 + (t184 * t350 - t506) * t253) * t477 + t81 * t479 + (t80 + t92) * t480 + (t458 - (-t114 + t63) * t185 + (-t112 + t62) * t184 + (-t184 * t318 - t185 * t350) * t253 + (-t184 * t351 + t506) * t252 + t29) * t476 + (t40 + t45) * t475 + (t444 + ((t124 - t67 + (t144 - t436) * t253) * t253 + (-t122 + t66 + (t144 + t438) * t252) * t252) * qJD(4) + t33) * t360 + (qJD(4) * t321 + t206 * t275 + t207 * t273 - t260 * t347 + t262 * t348) * t270 + (t41 + t44 + t28) * t478 - (t50 + t54 + t32) * t391 / 0.2e1 - (t55 - t513) * t390 / 0.2e1 + (t253 * t111 + (t110 + t96) * t252) * qJD(4) * t471;
t166 = pkin(3) * t429 - t224;
t215 = t335 * qJD(4);
t60 = -t215 * t391 + (-t105 - t166 - t371) * t270 + t314;
t61 = t106 * t270 + (t215 * t253 - t240 * t429) * qJD(4) + t375;
t281 = (-t317 * t75 - t364 * t459 - t463 * t61 + t470 * t60) * t252 + (-t317 * t76 + t364 * t60 - t381 * t75 - t470 * t61) * t253;
t179 = t203 * t269;
t36 = -t252 * t386 - t163 * t201 - t179 * t184 + (-t166 - t345 + t457) * t270 + t314;
t108 = t190 + (-t377 - t415) * t252 - t397;
t37 = t253 * t386 - t162 * t201 + t179 * t185 + (t108 + t90 - t346) * t270 + t375;
t280 = (t36 * (rSges(6,3) - t277) - t37 * t462 + t56 * t291 + (t57 * (-t256 - t465) - t56 * t277) * t270) * t252 + (t36 * t465 + t57 * (t291 - t415) - t37 * rSges(6,3) - t56 * t380) * t253;
t174 = t240 * t253;
t173 = t240 * t252;
t142 = t181 * t270 + t257;
t141 = -t180 * t270 - t385;
t121 = t252 * t137;
t118 = t270 * (-rSges(4,2) * t429 + t222) + t393;
t117 = -t186 * t268 + t314;
t82 = qJD(4) * t324 + qJD(3);
t46 = t288 * qJD(4);
t12 = t137 * t163 - t138 * t162 + t184 * t90 - t185 * t89 + ((t125 * t270 - t107) * t253 + (t108 + t120) * t252) * qJD(4);
t1 = [t282 + m(3) * (t141 * (t204 + t267) + t142 * (t266 + t500) + (-t165 + t181 + t258) * t164) + (t36 * (t267 + t374) + t57 * (t307 - t404) + t37 * (t266 + t343) + t280 + (t57 + t496) * t56) * m(6) + (t60 * (t267 + t396) + t76 * (t307 + t502) + t61 * (t266 + t373) + t281 + (t76 + t497) * t75) * m(5) + (t117 * (t267 + t339) - t129 * t378 + t118 * (t266 + t499) + t283 + (t129 + t503) * t128) * m(4); t282 + (t320 * t460 + t37 * t343 + t36 * t374 + t280 + (-t387 - t404 - t319) * t57 + t496 * t56) * m(6) + (t340 * t459 + t61 * t373 + t60 * t396 + t281 + (-t387 + t371 + t502) * t76 + t497 * t75) * m(5) + (t117 * t339 + t118 * t499 + t128 * t503 + t129 * t505 + t283) * m(4) + (-(t164 * t204 - t165 * t500) * t270 + t141 * t204 + t142 * t500 + t164 * t181 - t165 * t180) * m(3); m(5) * t46 + m(6) * t12; ((t168 * t391 - t309) * t252 + (t510 + (-t294 * t253 + (-t432 + t293) * t252) * qJD(4)) * t253) * t362 + ((-t273 * t395 + t275 * t394) * t270 + ((t252 * t405 - t253 * t406) * t275 + (-t252 * t407 + t253 * t408) * t273) * qJD(4)) * t472 + (t513 * t253 + (t270 * t96 - t50) * t252) * t471 + t292 + ((-t390 * t432 - t309) * t253 + (-t510 + (-t293 * t252 + (t168 + t294) * t253) * qJD(4)) * t252) * t360 + (t270 * t54 + ((-t297 * t252 - t253 * t493 - t270 * t69) * t253 + (-t298 * t252 + t253 * t487 + t270 * t68) * t252) * t507) * t474 + (t270 * t55 + ((t252 * t493 - t297 * t253 - t270 * t67) * t253 + (-t252 * t487 - t298 * t253 + t270 * t66) * t252) * t507) * t473 + (t32 + t316) * t366 + (t33 + t315) * t365 + (-(-t57 * t369 + ((-t252 * t56 - t253 * t57) * t270 + t51 * (-t252 ^ 2 - t253 ^ 2) * qJD(4)) * t273) * pkin(4) + t12 * t121 + t51 * t376 + (t12 * t125 + t51 * t108 - t36 * t356 + t57 * (-pkin(4) * t388 - t179) + (t51 * t126 - t356 * t56) * t270) * t252 + (t12 * t413 + t51 * t457 + t37 * t356 + t56 * t179 + (t51 * t125 - t356 * t57) * t270) * t253 + t494) * m(6) + (t46 * t324 + t82 * t288 + t334 * t215 + ((t61 - t459) * t253 + (-t270 * t75 - t60) * t252) * t240 - (-t173 * t75 - t174 * t76) * t270 - (t82 * (-t173 * t252 - t174 * t253) + t334 * t335) * qJD(4)) * m(5); t292 + (t12 * (t138 * t253 + t121) + t51 * (-t253 * t89 + t376) + (-t252 * t57 + t253 * t56) * t179 + ((t37 - t460) * t253 + (-t270 * t56 - t36) * t252) * t201 + t494) * m(6);];
tauc = t1(:);
