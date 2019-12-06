% Calculate vector of centrifugal and Coriolis load on the joints for
% S5PRRRR4
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% qJD [5x1]
%   Generalized joint velocities
% pkin [9x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d2,d3,d4,d5,theta1]';
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
% Datum: 2019-12-05 17:08
% Revision: 77da58f92bca3eff71542919beafa37024070d86 (2019-12-05)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauc = S5PRRRR4_coriolisvecJ_fixb_slag_vp1(qJ, qJD, ...
  pkin, m, rSges, Icges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(9,1),zeros(6,1),zeros(6,3),zeros(6,6)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5PRRRR4_coriolisvecJ_fixb_slag_vp1: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5PRRRR4_coriolisvecJ_fixb_slag_vp1: qJD has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S5PRRRR4_coriolisvecJ_fixb_slag_vp1: pkin has to be [9x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5PRRRR4_coriolisvecJ_fixb_slag_vp1: m has to be [6x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [6,3]), ...
  'S5PRRRR4_coriolisvecJ_fixb_slag_vp1: rSges has to be [6x3] (double)');
assert(isreal(Icges) && all(size(Icges) == [6 6]), ...
  'S5PRRRR4_coriolisvecJ_fixb_slag_vp1: Icges has to be [6x6] (double)'); 

%% Symbolic Calculation
% From coriolisvec_joint_fixb_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-05 17:07:40
% EndTime: 2019-12-05 17:07:51
% DurationCPUTime: 6.09s
% Computational Cost: add. (17428->514), mult. (10814->675), div. (0->0), fcn. (8468->8), ass. (0->319)
t268 = qJ(4) + qJ(5);
t263 = cos(t268);
t255 = Icges(6,4) * t263;
t262 = sin(t268);
t203 = Icges(6,1) * t262 + t255;
t327 = -Icges(6,2) * t262 + t255;
t476 = t203 + t327;
t265 = pkin(9) + qJ(2);
t261 = qJ(3) + t265;
t256 = sin(t261);
t419 = t256 * t262;
t213 = Icges(6,4) * t419;
t257 = cos(t261);
t418 = t256 * t263;
t137 = Icges(6,1) * t418 - Icges(6,5) * t257 - t213;
t135 = Icges(6,4) * t418 - Icges(6,2) * t419 - Icges(6,6) * t257;
t430 = t135 * t262;
t326 = -t137 * t263 + t430;
t308 = t326 * t256;
t200 = Icges(6,5) * t263 - Icges(6,6) * t262;
t309 = t200 * t257;
t134 = Icges(6,3) * t256 + t309;
t434 = Icges(6,4) * t262;
t204 = Icges(6,1) * t263 - t434;
t313 = t204 * t257;
t138 = Icges(6,5) * t256 + t313;
t411 = t257 * t263;
t404 = t256 * t134 + t138 * t411;
t475 = -t308 - t404;
t259 = sin(t265);
t441 = pkin(2) * qJD(2);
t374 = t259 * t441;
t186 = rSges(4,1) * t256 + rSges(4,2) * t257;
t267 = qJD(2) + qJD(3);
t426 = t186 * t267;
t153 = -t374 - t426;
t266 = qJD(4) + qJD(5);
t442 = rSges(6,2) * t263;
t375 = t266 * t442;
t407 = t262 * t266;
t474 = -rSges(6,1) * t407 - t375;
t473 = 2 * qJD(4);
t389 = rSges(6,1) * t411 + t256 * rSges(6,3);
t270 = cos(qJ(4));
t408 = t257 * t270;
t472 = rSges(5,1) * t408 + t256 * rSges(5,3);
t251 = t256 * pkin(7);
t189 = t257 * pkin(3) + t251;
t264 = Icges(5,4) * t270;
t269 = sin(qJ(4));
t328 = -Icges(5,2) * t269 + t264;
t235 = Icges(5,1) * t269 + t264;
t416 = t256 * t269;
t388 = rSges(5,2) * t416 + t257 * rSges(5,3);
t415 = t256 * t270;
t151 = rSges(5,1) * t415 - t388;
t139 = t267 * t151;
t252 = t257 * pkin(7);
t188 = pkin(3) * t256 - t252;
t181 = t267 * t188;
t410 = t257 * t267;
t225 = pkin(7) * t410;
t382 = qJD(4) * t270;
t372 = rSges(5,2) * t382;
t406 = t267 * t269;
t315 = rSges(5,2) * t256 * t406 + rSges(5,3) * t410 - t257 * t372;
t383 = qJD(4) * t269;
t363 = t257 * t383;
t471 = -rSges(5,1) * t363 + t139 + t181 + t225 + t315;
t199 = Icges(6,5) * t262 + Icges(6,6) * t263;
t295 = Icges(6,3) * t267 - t199 * t266;
t417 = t256 * t267;
t311 = t327 * t257;
t136 = Icges(6,6) * t256 + t311;
t429 = t136 * t262;
t470 = -t200 * t417 + t257 * t295 + t267 * (-t138 * t263 + t429);
t469 = t256 * t295 + (t309 + t326) * t267;
t232 = Icges(5,5) * t270 - Icges(5,6) * t269;
t231 = Icges(5,5) * t269 + Icges(5,6) * t270;
t292 = Icges(5,3) * t267 - qJD(4) * t231;
t435 = Icges(5,4) * t269;
t236 = Icges(5,1) * t270 - t435;
t314 = t236 * t257;
t149 = Icges(5,5) * t256 + t314;
t312 = t328 * t257;
t147 = Icges(5,6) * t256 + t312;
t427 = t147 * t269;
t323 = -t149 * t270 + t427;
t468 = -t232 * t417 + t257 * t292 + t267 * t323;
t310 = t232 * t257;
t228 = Icges(5,4) * t416;
t148 = Icges(5,1) * t415 - Icges(5,5) * t257 - t228;
t146 = Icges(5,4) * t415 - Icges(5,2) * t416 - Icges(5,6) * t257;
t428 = t146 * t269;
t324 = -t148 * t270 + t428;
t467 = t256 * t292 + (t310 + t324) * t267;
t201 = Icges(6,2) * t263 + t434;
t321 = t201 * t262 - t203 * t263;
t466 = t200 * t266 + t267 * t321;
t237 = rSges(5,1) * t269 + rSges(5,2) * t270;
t385 = qJD(4) * t256;
t182 = t237 * t385;
t409 = t257 * t269;
t152 = -rSges(5,2) * t409 + t472;
t304 = t152 + t189;
t465 = t267 * t304 - t182;
t233 = Icges(5,2) * t270 + t435;
t320 = t269 * t233 - t270 * t235;
t464 = t232 * qJD(4) + t267 * t320;
t144 = Icges(5,5) * t415 - Icges(5,6) * t416 - Icges(5,3) * t257;
t66 = -t257 * t144 - t256 * t324;
t271 = -pkin(8) - pkin(7);
t239 = t257 * t271;
t447 = pkin(4) * t270;
t258 = pkin(3) + t447;
t390 = -t256 * t258 - t239;
t126 = t188 + t390;
t118 = t267 * t126;
t216 = rSges(6,2) * t419;
t140 = rSges(6,1) * t418 - t257 * rSges(6,3) - t216;
t128 = t267 * t140;
t371 = t263 * t417;
t393 = rSges(6,3) * t410 + t267 * t216;
t463 = -rSges(6,1) * t371 - t258 * t417 - t118 + t128 + t181 + t393;
t340 = pkin(4) * t363;
t446 = pkin(3) - t258;
t108 = -t340 - t225 + (t256 * t446 - t239) * t267;
t364 = t256 * t383;
t221 = pkin(4) * t364;
t414 = t256 * t271;
t391 = -t267 * t414 - t221;
t109 = (-t257 * t446 - t251) * t267 + t391;
t342 = t257 * t258 - t414;
t127 = t342 - t189;
t412 = t257 * t262;
t376 = rSges(6,2) * t412;
t141 = -t376 + t389;
t341 = t266 * t267;
t164 = t256 * t341;
t165 = t257 * t341;
t183 = t256 * t266;
t184 = t257 * t266;
t90 = -t257 * t375 + (-t257 * t407 - t371) * rSges(6,1) + t393;
t367 = t256 * t474 - t267 * t376;
t91 = t389 * t267 + t367;
t12 = t140 * t165 - t141 * t164 + t183 * t91 + t184 * t90 + ((t108 - t118) * t257 + (-t127 * t267 + t109) * t256) * qJD(4);
t205 = rSges(6,1) * t262 + t442;
t162 = t205 * t256;
t163 = t205 * t257;
t443 = rSges(6,2) * t262;
t444 = rSges(6,1) * t263;
t206 = -t443 + t444;
t51 = t140 * t183 + t141 * t184 + qJD(1) + (-t126 * t256 + t127 * t257) * qJD(4);
t307 = -t184 * t205 - t340;
t290 = t307 - t374;
t58 = (t126 - t140 - t188) * t267 + t290;
t401 = -t127 - t141;
t368 = t189 - t401;
t260 = cos(t265);
t373 = t260 * t441;
t398 = t183 * t205 + t221;
t59 = t267 * t368 + t373 - t398;
t462 = -(t162 * t267 - t184 * t206) * t58 - t59 * (-t267 * t163 - t183 * t206) - t51 * (-t183 * t162 - t163 * t184) + t12 * (t256 * t140 + t257 * t141);
t395 = -Icges(5,2) * t415 + t148 - t228;
t397 = t235 * t256 + t146;
t461 = -t269 * t395 - t270 * t397;
t460 = t183 * (-t201 * t257 + t138) - t184 * (-Icges(6,2) * t418 + t137 - t213) + t267 * t476;
t459 = t164 / 0.2e1;
t458 = t165 / 0.2e1;
t457 = -t183 / 0.2e1;
t456 = t183 / 0.2e1;
t455 = -t184 / 0.2e1;
t454 = t184 / 0.2e1;
t453 = t256 / 0.2e1;
t452 = -t257 / 0.2e1;
t451 = -t267 / 0.2e1;
t450 = t267 / 0.2e1;
t449 = pkin(2) * t259;
t448 = pkin(2) * qJD(2) ^ 2;
t445 = rSges(5,1) * t270;
t384 = qJD(4) * t257;
t365 = t237 * t384;
t306 = -t365 - t374;
t79 = (-t151 - t188) * t267 + t306;
t440 = t257 * t79;
t439 = t267 * t58;
t424 = t199 * t257;
t93 = -t256 * t321 - t424;
t438 = t93 * t267;
t421 = t231 * t257;
t111 = -t256 * t320 - t421;
t431 = t111 * t267;
t425 = t199 * t256;
t423 = t201 * t266;
t422 = t231 * t256;
t420 = t232 * t267;
t133 = Icges(6,5) * t418 - Icges(6,6) * t419 - Icges(6,3) * t257;
t405 = -t256 * t133 - t137 * t411;
t403 = -t256 * t144 - t148 * t408;
t145 = Icges(5,3) * t256 + t310;
t402 = t256 * t145 + t149 * t408;
t396 = -t235 * t257 - t147;
t394 = -t233 * t257 + t149;
t387 = -t233 + t236;
t386 = t235 + t328;
t381 = t259 * t448;
t380 = t260 * t448;
t379 = (qJD(4) ^ 2) * t447;
t378 = t256 * t91 + (t128 + t90) * t257;
t369 = t257 * t406;
t366 = rSges(5,1) * t364 + rSges(5,2) * t369 + t256 * t372;
t362 = t417 / 0.2e1;
t361 = t410 / 0.2e1;
t360 = -pkin(3) - t445;
t359 = -t385 / 0.2e1;
t356 = t384 / 0.2e1;
t354 = -pkin(4) * t269 - t205;
t187 = t257 * rSges(4,1) - rSges(4,2) * t256;
t297 = Icges(6,5) * t267 - t203 * t266;
t352 = -t135 * t266 + t256 * t297 + t267 * t313;
t351 = -t136 * t266 - t204 * t417 + t257 * t297;
t296 = Icges(6,6) * t267 - t423;
t350 = t137 * t266 + t256 * t296 + t267 * t311;
t349 = t138 * t266 + t257 * t296 - t327 * t417;
t123 = t149 * t415;
t348 = t257 * t145 - t123;
t347 = -t133 + t429;
t345 = -t144 + t427;
t344 = t476 * t266;
t343 = t204 * t266 - t423;
t337 = t267 * (-pkin(3) * t417 + t225) - t381;
t167 = rSges(4,1) * t410 - rSges(4,2) * t417;
t180 = t206 * t266;
t336 = -pkin(4) * t382 - t180;
t113 = t138 * t418;
t335 = t136 * t419 - t113;
t333 = -rSges(5,2) * t269 + t445;
t332 = -t256 * t59 - t257 * t58;
t80 = t373 + t465;
t331 = -t256 * t80 - t440;
t81 = t135 * t263 + t137 * t262;
t97 = t146 * t270 + t148 * t269;
t98 = t147 * t270 + t269 * t149;
t322 = t151 * t256 + t152 * t257;
t319 = t342 + t389;
t67 = -t147 * t416 - t348;
t317 = (t256 * t67 - t257 * t66) * qJD(4);
t68 = -t146 * t409 - t403;
t69 = -t147 * t409 + t402;
t316 = (t256 * t69 - t257 * t68) * qJD(4);
t305 = -t140 + t390;
t303 = t183 * t424 - t184 * t425 - t200 * t267;
t302 = -t269 * t394 + t270 * t396;
t301 = t256 * t360 + t252 + t388;
t300 = -rSges(6,3) * t417 - t367 - t391;
t283 = t133 * t267 - t262 * t350 + t263 * t352;
t13 = t469 * t256 + t283 * t257;
t282 = t134 * t267 - t262 * t349 + t263 * t351;
t14 = t470 * t256 + t282 * t257;
t15 = t283 * t256 - t469 * t257;
t16 = t282 * t256 - t470 * t257;
t286 = (-t203 * t257 - t136) * t183 - (-t203 * t256 - t135) * t184 + (-t201 + t204) * t267;
t274 = -t262 * t460 + t286 * t263;
t62 = -t133 * t257 - t308;
t63 = -t134 * t257 - t335;
t28 = t183 * t63 - t184 * t62 + t438;
t64 = -t135 * t412 - t405;
t65 = -t136 * t412 + t404;
t94 = -t257 * t321 + t425;
t92 = t94 * t267;
t29 = t183 * t65 - t184 * t64 + t92;
t40 = t262 * t352 + t263 * t350;
t41 = t262 * t351 + t263 * t349;
t281 = t199 * t267 - t262 * t344 + t263 * t343;
t44 = t466 * t256 + t281 * t257;
t45 = t281 * t256 - t466 * t257;
t82 = t136 * t263 + t138 * t262;
t299 = (-t13 * t184 + t14 * t183 + t164 * t64 + t165 * t65 + t267 * t44) * t453 + (-t256 * t303 + t257 * t274) * t457 + (t256 * t274 + t257 * t303) * t454 + (-t15 * t184 + t16 * t183 + t164 * t62 + t165 * t63 + t267 * t45) * t452 + (t286 * t262 + t263 * t460) * t451 + t28 * t362 + t29 * t361 + ((t267 * t65 - t13) * t257 + (t267 * t64 + t14) * t256) * t456 + (t256 * t63 - t257 * t62) * t459 + (t256 * t65 - t257 * t64) * t458 + ((t267 * t63 - t15) * t257 + (t267 * t62 + t16) * t256) * t455 + ((t267 * t82 - t40) * t257 + (t267 * t81 + t41) * t256) * t450;
t298 = (-t269 * t386 + t270 * t387) * t267;
t294 = Icges(5,5) * t267 - qJD(4) * t235;
t293 = Icges(5,6) * t267 - qJD(4) * t233;
t106 = (-t267 * t415 - t363) * rSges(5,1) + t315;
t107 = t472 * t267 - t366;
t288 = (t106 + t139) * t257 + (-t152 * t267 + t107) * t256;
t102 = t257 * t293 - t328 * t417;
t104 = -t236 * t417 + t257 * t294;
t280 = -qJD(4) * t98 - t102 * t269 + t104 * t270 + t145 * t267;
t103 = t256 * t293 + t267 * t312;
t105 = t256 * t294 + t267 * t314;
t279 = -qJD(4) * t97 - t103 * t269 + t105 * t270 + t144 * t267;
t209 = t328 * qJD(4);
t210 = t236 * qJD(4);
t278 = -t209 * t269 + t210 * t270 + t231 * t267 + (-t233 * t270 - t235 * t269) * qJD(4);
t277 = (t360 * t440 + (t79 * (-rSges(5,3) - pkin(7)) + t80 * t360) * t256) * t267;
t112 = -t257 * t320 + t422;
t110 = t112 * t267;
t32 = t317 + t431;
t33 = t110 + t316;
t49 = -qJD(4) * t324 + t103 * t270 + t105 * t269;
t50 = -qJD(4) * t323 + t102 * t270 + t104 * t269;
t54 = t464 * t256 + t278 * t257;
t55 = t278 * t256 - t464 * t257;
t276 = (t110 + ((t67 - t123 + (t145 + t428) * t257 + t403) * t257 + t402 * t256) * qJD(4)) * t356 + (t92 + (t63 + (t134 + t430) * t257 + t335 + t405) * t184 + (-t257 * t347 - t475 + t62) * t183) * t454 + (t81 + t93) * t459 + (t82 + t94) * t458 + (t28 - t438 + (t65 + t475) * t184 + (t256 * t347 - t113 + t64) * t183 + ((t134 + t326) * t183 + t347 * t184) * t257) * t457 + (t41 + t44) * t456 + (-t431 + ((t257 * t345 - t402 + t69) * t257 + (t256 * t345 + t348 + t68) * t256) * qJD(4) + t32) * t359 + (t50 + t54) * t385 / 0.2e1 + (-qJD(4) * t320 + t209 * t270 + t210 * t269 + t262 * t343 + t263 * t344) * t267 + (t40 + t45 + t29) * t455 - (t49 + t55 + t33) * t384 / 0.2e1 + ((t111 + t97) * t256 + (t112 + t98) * t257) * qJD(4) * t450;
t36 = -t256 * t379 - t165 * t205 - t180 * t183 + (t108 + t90 - t340) * t267 + t337;
t275 = (-t36 * t443 + t59 * (-pkin(4) * t383 + t474) + (t58 * (-t258 - t444) - t59 * t271) * t267) * t257;
t254 = pkin(2) * t260;
t218 = t333 * qJD(4);
t176 = t237 * t257;
t175 = t237 * t256;
t168 = t189 * t267;
t154 = t187 * t267 + t373;
t132 = -t167 * t267 - t380;
t131 = -t267 * t426 - t381;
t83 = qJD(4) * t322 + qJD(1);
t61 = -t380 - t218 * t384 + (-t107 - t168 + t182) * t267;
t60 = t106 * t267 + (-t218 * t256 - t237 * t410) * qJD(4) + t337;
t46 = t288 * qJD(4);
t37 = -t257 * t379 - t380 + t164 * t205 - t180 * t184 + (-t109 - t168 - t91 + t221) * t267;
t1 = [m(5) * t46 + m(6) * t12; m(4) * (t132 * (-t186 - t449) + t131 * (t187 + t254) + (-t167 - t373 + t154) * t153) + t276 + (t37 * (t305 - t449) + t58 * (t300 - t373) + t36 * (t254 + t319) + t275 + (-t290 + t58 - t374 + t463) * t59) * m(6) + (t61 * (t301 - t449) + t79 * (t366 - t373) + t60 * (t254 + t304) + t277 + (-t306 + t79 - t374 + t471) * t80) * m(5); t276 + (t37 * t305 + t36 * t319 + t368 * t439 + t275 + (-t307 + t463) * t59 + (-t398 + t300) * t58) * m(6) + (t61 * t301 + t60 * t304 + t277 + (t365 + t471) * t80 + (t366 + t465) * t79) * m(5) + (-(-t153 * t187 - t154 * t186) * t267 + t131 * t187 - t132 * t186 - t153 * t167 - t154 * t426) * m(4); ((t269 * t387 + t270 * t386) * t267 + ((t256 * t394 - t257 * t395) * t270 + (t256 * t396 + t257 * t397) * t269) * qJD(4)) * t451 + t299 + ((-t384 * t422 - t420) * t257 + (t298 + (t302 * t256 + (t421 - t461) * t257) * qJD(4)) * t256) * t356 + ((-t385 * t421 + t420) * t256 + (t298 + (-t461 * t257 + (t422 + t302) * t256) * qJD(4)) * t257) * t359 + ((t267 * t98 - t49) * t257 + (t267 * t97 + t50) * t256) * t450 + (t54 * t267 + ((-t467 * t256 - t279 * t257 + t267 * t69) * t257 + (t468 * t256 + t280 * t257 + t267 * t68) * t256) * t473) * t453 + (t55 * t267 + ((-t279 * t256 + t467 * t257 + t267 * t67) * t257 + (t280 * t256 - t468 * t257 + t267 * t66) * t256) * t473) * t452 + (t317 + t32) * t362 + (t316 + t33) * t361 + (t51 * t378 + (t37 * t354 + t58 * t336 + t12 * t127 + t51 * t108 + (-t51 * t126 + t354 * t59) * t267) * t257 + (t36 * t354 + t59 * t336 - t12 * t126 + t51 * t109 + (t58 * t205 + t401 * t51) * t267) * t256 - (-t59 * t369 + (t332 * t270 + t51 * (-t256 ^ 2 - t257 ^ 2) * t269) * qJD(4)) * pkin(4) + t462) * m(6) + (-(t175 * t79 - t176 * t80) * t267 - (t83 * (-t175 * t256 - t176 * t257) + t331 * t333) * qJD(4) + t46 * t322 + t83 * t288 + t331 * t218 + ((-t267 * t80 - t61) * t257 + (t267 * t79 - t60) * t256) * t237) * m(5); t299 + (t51 * (-t141 * t417 + t378) + t332 * t180 + ((-t267 * t59 - t37) * t257 + (-t36 + t439) * t256) * t205 + t462) * m(6);];
tauc = t1(:);
