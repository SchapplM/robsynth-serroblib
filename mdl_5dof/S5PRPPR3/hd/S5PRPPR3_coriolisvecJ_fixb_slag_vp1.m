% Calculate vector of centrifugal and Coriolis load on the joints for
% S5PRPPR3
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% qJD [5x1]
%   Generalized joint velocities
% pkin [8x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d2,d5,theta1,theta3]';
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
% Datum: 2019-12-05 15:27
% Revision: 77da58f92bca3eff71542919beafa37024070d86 (2019-12-05)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauc = S5PRPPR3_coriolisvecJ_fixb_slag_vp1(qJ, qJD, ...
  pkin, m, rSges, Icges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(8,1),zeros(6,1),zeros(6,3),zeros(6,6)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5PRPPR3_coriolisvecJ_fixb_slag_vp1: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5PRPPR3_coriolisvecJ_fixb_slag_vp1: qJD has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S5PRPPR3_coriolisvecJ_fixb_slag_vp1: pkin has to be [8x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5PRPPR3_coriolisvecJ_fixb_slag_vp1: m has to be [6x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [6,3]), ...
  'S5PRPPR3_coriolisvecJ_fixb_slag_vp1: rSges has to be [6x3] (double)');
assert(isreal(Icges) && all(size(Icges) == [6 6]), ...
  'S5PRPPR3_coriolisvecJ_fixb_slag_vp1: Icges has to be [6x6] (double)'); 

%% Symbolic Calculation
% From coriolisvec_joint_fixb_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-05 15:26:18
% EndTime: 2019-12-05 15:26:41
% DurationCPUTime: 16.80s
% Computational Cost: add. (8196->533), mult. (12596->823), div. (0->0), fcn. (11595->8), ass. (0->294)
t238 = qJ(2) + pkin(8);
t234 = sin(t238);
t235 = cos(t238);
t243 = sin(qJ(2));
t245 = cos(qJ(2));
t452 = -Icges(3,5) * t243 - Icges(3,6) * t245 + (Icges(5,5) - Icges(4,6)) * t235 + (Icges(5,4) - Icges(4,5)) * t234;
t240 = cos(pkin(7));
t367 = qJD(5) * t235;
t239 = sin(pkin(7));
t375 = qJD(2) * t239;
t206 = t240 * t367 + t375;
t242 = sin(qJ(5));
t244 = cos(qJ(5));
t331 = rSges(6,1) * t242 + rSges(6,2) * t244;
t266 = -rSges(6,3) * t234 + t235 * t331;
t455 = t206 * t266;
t374 = qJD(2) * t240;
t207 = t239 * t367 - t374;
t454 = t207 * t266;
t236 = t239 ^ 2;
t237 = t240 ^ 2;
t378 = t236 + t237;
t453 = t452 * qJD(2);
t407 = Icges(3,4) * t243;
t314 = -Icges(3,2) * t245 - t407;
t406 = Icges(3,4) * t245;
t319 = -Icges(3,1) * t243 - t406;
t376 = qJD(2) * t235;
t377 = qJD(2) * t234;
t398 = Icges(5,6) * t235;
t404 = Icges(4,4) * t235;
t405 = Icges(4,4) * t234;
t451 = -(-t243 * t314 + t245 * t319) * qJD(2) + (-Icges(5,6) * t234 - t405 + (-Icges(4,2) - Icges(5,3)) * t235) * t377 + (t398 + t404 + (Icges(4,1) + Icges(5,2)) * t234) * t376;
t305 = Icges(5,3) * t234 - t398;
t124 = Icges(5,5) * t239 + t240 * t305;
t125 = -Icges(5,5) * t240 + t239 * t305;
t393 = t234 * t240;
t227 = Icges(5,6) * t393;
t391 = t235 * t240;
t126 = Icges(5,4) * t239 - Icges(5,2) * t391 + t227;
t394 = t234 * t239;
t226 = Icges(5,6) * t394;
t392 = t235 * t239;
t127 = -Icges(5,4) * t240 - Icges(5,2) * t392 + t226;
t313 = -Icges(4,2) * t234 + t404;
t318 = Icges(4,1) * t235 - t405;
t438 = -(Icges(4,6) * t239 + t240 * t313) * t235 - (Icges(4,5) * t239 + t240 * t318) * t234;
t439 = (-Icges(4,6) * t240 + t239 * t313) * t235 + (-Icges(4,5) * t240 + t239 * t318) * t234;
t450 = (t124 * t239 - t125 * t240) * t235 + ((Icges(5,3) * t391 + t126 + t227) * t239 - (Icges(5,3) * t392 + t127 + t226) * t240) * t234 + t239 * t438 + t240 * t439;
t370 = qJD(4) * t239;
t222 = t234 * t370;
t212 = pkin(3) * t234 - qJ(4) * t235;
t284 = qJD(2) * t212;
t102 = -t239 * t284 + t222;
t369 = qJD(4) * t240;
t224 = t234 * t369;
t103 = -t240 * t284 + t224;
t412 = pkin(2) * qJD(2);
t363 = t243 * t412;
t372 = qJD(3) * t240;
t210 = -t239 * t363 - t372;
t233 = qJD(3) * t239;
t211 = -t240 * t363 + t233;
t380 = t239 * t210 + t240 * t211;
t449 = -qJD(4) * t234 + t239 * t102 + t240 * t103 + t284 * t378 + t380;
t445 = t452 * t239;
t444 = t453 * t239;
t443 = t453 * t240;
t442 = t452 * t240;
t315 = -Icges(3,2) * t243 + t406;
t156 = Icges(3,6) * t239 + t240 * t315;
t320 = Icges(3,1) * t245 - t407;
t158 = Icges(3,5) * t239 + t240 * t320;
t441 = t451 * t240 + (-t124 * t235 - t126 * t234 + t156 * t245 + t158 * t243 - t438) * qJD(2);
t155 = -Icges(3,6) * t240 + t239 * t315;
t157 = -Icges(3,5) * t240 + t239 * t320;
t440 = t451 * t239 + (-t125 * t235 - t127 * t234 + t155 * t245 + t157 * t243 + t439) * qJD(2);
t330 = rSges(5,2) * t234 + rSges(5,3) * t235;
t437 = t378 * qJD(2) * t330;
t425 = t206 / 0.2e1;
t423 = t207 / 0.2e1;
t353 = pkin(6) * t378;
t228 = rSges(3,1) * t243 + rSges(3,2) * t245;
t287 = qJD(2) * t228;
t307 = Icges(6,5) * t242 + Icges(6,6) * t244;
t263 = -Icges(6,3) * t234 + t235 * t307;
t401 = Icges(6,4) * t242;
t310 = Icges(6,2) * t244 + t401;
t264 = -Icges(6,6) * t234 + t235 * t310;
t400 = Icges(6,4) * t244;
t316 = Icges(6,1) * t242 + t400;
t265 = -Icges(6,5) * t234 + t235 * t316;
t429 = -t243 * (t314 * t239 + t157) - t245 * (-t319 * t239 + t155);
t387 = t240 * t244;
t390 = t239 * t242;
t192 = t234 * t387 - t390;
t174 = Icges(6,4) * t192;
t388 = t240 * t242;
t389 = t239 * t244;
t194 = t234 * t389 + t388;
t175 = Icges(6,4) * t194;
t179 = (Icges(6,2) * t242 - t400) * t235;
t193 = t234 * t388 + t389;
t195 = -t234 * t390 + t387;
t368 = qJD(5) * t234;
t71 = Icges(6,1) * t193 + Icges(6,5) * t391 + t174;
t72 = -Icges(6,1) * t195 + Icges(6,5) * t392 + t175;
t254 = t206 * (-Icges(6,2) * t193 + t174 + t71) + t207 * (Icges(6,2) * t195 + t175 + t72) + t368 * (-t265 + t179);
t180 = (-Icges(6,1) * t244 + t401) * t235;
t402 = Icges(6,4) * t195;
t403 = Icges(6,4) * t193;
t69 = Icges(6,2) * t192 + Icges(6,6) * t391 + t403;
t70 = Icges(6,2) * t194 + Icges(6,6) * t392 - t402;
t255 = t206 * (-Icges(6,1) * t192 + t403 + t69) + t207 * (-Icges(6,1) * t194 - t402 + t70) + t368 * (-t264 - t180);
t246 = qJD(2) ^ 2;
t110 = Icges(6,3) * t235 + t234 * t307;
t178 = (-Icges(6,5) * t244 + Icges(6,6) * t242) * t235;
t63 = qJD(2) * t110 + qJD(5) * t178;
t291 = t235 * t63 + t263 * t377;
t67 = Icges(6,5) * t193 + Icges(6,6) * t192 + Icges(6,3) * t391;
t19 = t192 * t69 + t193 * t71 + t391 * t67;
t68 = -Icges(6,5) * t195 + Icges(6,6) * t194 + Icges(6,3) * t392;
t20 = t192 * t70 + t193 * t72 + t391 * t68;
t411 = t20 * t239;
t328 = t19 * t240 + t411;
t43 = -t192 * t264 - t193 * t265 - t263 * t391;
t409 = t43 * t235;
t112 = Icges(6,6) * t235 + t234 * t310;
t64 = qJD(2) * t112 + qJD(5) * t179;
t114 = Icges(6,5) * t235 + t234 * t316;
t65 = qJD(2) * t114 + qJD(5) * t180;
t356 = t242 * t376;
t96 = qJD(5) * t192 + t240 * t356;
t355 = t244 * t376;
t97 = -qJD(5) * t193 + t240 * t355;
t253 = (t192 * t64 + t193 * t65 + t240 * t291 - t264 * t97 - t265 * t96) * t234 + (-t234 * t328 + t409) * qJD(2);
t357 = t234 * t374;
t48 = Icges(6,5) * t96 + Icges(6,6) * t97 - Icges(6,3) * t357;
t294 = t235 * t48 - t377 * t67;
t50 = Icges(6,4) * t96 + Icges(6,2) * t97 - Icges(6,6) * t357;
t52 = Icges(6,1) * t96 + Icges(6,4) * t97 - Icges(6,5) * t357;
t7 = t192 * t50 + t193 * t52 + t240 * t294 + t69 * t97 + t71 * t96;
t358 = t234 * t375;
t98 = qJD(5) * t194 + t239 * t356;
t99 = qJD(5) * t195 + t239 * t355;
t49 = Icges(6,5) * t98 + Icges(6,6) * t99 - Icges(6,3) * t358;
t293 = t235 * t49 - t377 * t68;
t51 = Icges(6,4) * t98 + Icges(6,2) * t99 - Icges(6,6) * t358;
t53 = Icges(6,1) * t98 + Icges(6,4) * t99 - Icges(6,5) * t358;
t8 = t192 * t51 + t193 * t53 + t240 * t293 + t70 * t97 + t72 * t96;
t428 = qJD(5) * t253 / 0.2e1 + t7 * t425 + t8 * t423;
t427 = -t378 * t377 / 0.2e1;
t426 = -t206 / 0.2e1;
t424 = -t207 / 0.2e1;
t422 = -t235 / 0.2e1;
t421 = pkin(2) * t243;
t420 = pkin(2) * t245;
t21 = t194 * t69 - t195 * t71 + t392 * t67;
t410 = t21 * t240;
t44 = -t194 * t264 + t195 * t265 - t263 * t392;
t408 = t44 * t235;
t395 = t263 * t234;
t118 = -qJ(3) * t240 + t239 * t420;
t119 = qJ(3) * t239 + t240 * t420;
t386 = t239 * t118 + t240 * t119;
t381 = t210 * t375 + t211 * t374;
t379 = t224 + t233;
t371 = qJD(4) * t235;
t365 = qJD(2) * qJD(4);
t364 = t246 * t420;
t362 = t245 * t412;
t354 = t118 * t375 + t119 * t374 + qJD(1);
t352 = t235 * t365;
t351 = t376 / 0.2e1;
t350 = -t375 / 0.2e1;
t348 = -t368 / 0.2e1;
t347 = t368 / 0.2e1;
t346 = -t212 - t421;
t214 = rSges(4,1) * t234 + rSges(4,2) * t235;
t345 = -t214 - t421;
t217 = rSges(4,1) * t235 - rSges(4,2) * t234;
t344 = -t217 - t420;
t342 = t222 - t372;
t215 = pkin(3) * t235 + qJ(4) * t234;
t176 = t215 * t239;
t177 = t215 * t240;
t341 = t239 * t176 + t240 * t177 + t386;
t340 = t378 * t421;
t339 = t330 + t346;
t216 = -rSges(5,2) * t235 + rSges(5,3) * t234;
t338 = -t215 - t216 - t420;
t337 = qJD(2) * t348;
t336 = qJD(5) * t351;
t136 = qJD(2) * t215 - t371;
t335 = -t136 - t362;
t199 = t217 * qJD(2);
t334 = -t199 - t362;
t333 = -pkin(6) * t235 - t420;
t332 = t102 * t375 + t103 * t374 + t234 * t365 + t381;
t229 = rSges(3,1) * t245 - rSges(3,2) * t243;
t219 = t239 * t352;
t292 = t333 * t246;
t54 = rSges(6,1) * t96 + rSges(6,2) * t97 - rSges(6,3) * t357;
t116 = rSges(6,3) * t235 + t234 * t331;
t183 = (-rSges(6,1) * t244 + rSges(6,2) * t242) * t235;
t66 = qJD(2) * t116 + qJD(5) * t183;
t73 = rSges(6,1) * t193 + rSges(6,2) * t192 + rSges(6,3) * t391;
t16 = t54 * t368 - t206 * t66 + t219 + t239 * t292 + (-t136 * t239 + (t235 * t73 - t266 * t393) * qJD(5)) * qJD(2);
t220 = t240 * t352;
t55 = rSges(6,1) * t98 + rSges(6,2) * t99 - rSges(6,3) * t358;
t74 = -rSges(6,1) * t195 + rSges(6,2) * t194 + rSges(6,3) * t392;
t17 = -t55 * t368 + t207 * t66 + t220 + t240 * t292 + (-t136 * t240 + (-t235 * t74 + t266 * t394) * qJD(5)) * qJD(2);
t329 = -t16 * t240 + t17 * t239;
t327 = -t67 * t206 - t68 * t207;
t22 = t194 * t70 - t195 * t72 + t392 * t68;
t326 = t22 * t239 + t410;
t323 = t242 * t71 + t244 * t69;
t25 = t234 * t67 - t235 * t323;
t322 = t242 * t72 + t244 * t70;
t26 = t234 * t68 - t235 * t322;
t325 = t26 * t239 + t25 * t240;
t324 = t239 * t73 - t240 * t74;
t321 = qJD(2) * t345;
t304 = -t242 * t265 - t244 * t264;
t303 = t378 * t229;
t302 = t378 * t287;
t301 = t239 * t337;
t300 = t240 * t337;
t198 = t216 * qJD(2);
t299 = -t198 + t335;
t298 = -pkin(6) * t234 + t346;
t296 = qJD(2) * t339;
t295 = -qJD(2) * t199 - t364;
t208 = pkin(4) * t239 + pkin(6) * t391;
t209 = -pkin(4) * t240 + pkin(6) * t392;
t270 = t176 * t375 + t177 * t374 + t354 - t371;
t18 = t206 * t74 - t207 * t73 + (t208 * t240 + t209 * t239) * qJD(2) + t270;
t290 = t18 * t324;
t134 = rSges(5,1) * t239 + t216 * t240;
t135 = -rSges(5,1) * t240 + t216 * t239;
t23 = (t134 * t240 + t135 * t239) * qJD(2) + t270;
t289 = t23 * t330;
t132 = -rSges(4,3) * t240 + t217 * t239;
t133 = rSges(4,3) * t239 + t217 * t240;
t45 = (t132 * t239 + t133 * t240) * qJD(2) + t354;
t288 = t45 * t214;
t286 = qJD(2) * t214;
t283 = t266 + t298;
t273 = qJD(2) * t298;
t272 = -t364 + (-t136 - t198) * qJD(2);
t271 = (t110 + t304) * t234;
t269 = (t319 * t240 - t156) * t245 + (-t314 * t240 - t158) * t243;
t268 = t178 * t368 + t206 * (Icges(6,5) * t192 - Icges(6,6) * t193) + t207 * (Icges(6,5) * t194 + Icges(6,6) * t195);
t267 = -pkin(6) * t376 + t335 - t66;
t262 = t235 * t268;
t252 = (t194 * t64 - t195 * t65 + t239 * t291 - t264 * t99 - t265 * t98) * t234 + (-t234 * t326 + t408) * qJD(2);
t46 = -t235 * t304 - t395;
t251 = ((qJD(2) * t304 + t63) * t234 + (-qJD(2) * t263 - t242 * t65 - t244 * t64 + (-t242 * t264 + t244 * t265) * qJD(5)) * t235) * t234 + (-t234 * t325 + t46 * t235) * qJD(2);
t250 = (t263 * t240 + t323) * t206 + (t263 * t239 + t322) * t207;
t247 = (qJD(5) * t271 + t250) * t235;
t225 = t235 * t369;
t223 = t235 * t370;
t154 = t240 * t286;
t152 = t239 * t286;
t105 = t295 * t240;
t104 = t295 * t239;
t101 = t240 * t321 + t233;
t100 = t239 * t321 - t372;
t93 = t265 * t240;
t92 = t265 * t239;
t91 = t264 * t240;
t90 = t264 * t239;
t83 = rSges(6,1) * t194 + rSges(6,2) * t195;
t82 = rSges(6,1) * t192 - rSges(6,2) * t193;
t75 = t302 * qJD(2);
t62 = qJD(2) * t303 + qJD(1);
t61 = t240 * t296 + t379;
t60 = t239 * t296 + t342;
t59 = t240 * t272 + t220;
t58 = t239 * t272 + t219;
t47 = (-t152 * t239 - t154 * t240) * qJD(2) + t381;
t38 = t240 * t273 - t368 * t74 + t379 - t454;
t37 = t239 * t273 + t368 * t73 + t342 + t455;
t24 = qJD(2) * t437 + t332;
t12 = t206 * t55 - t207 * t54 + (qJD(2) * qJD(5) * t324 - t246 * t353) * t234 + t332;
t11 = t206 * t25 + t207 * t26 + t368 * t46;
t10 = t194 * t51 - t195 * t53 + t239 * t293 + t70 * t99 + t72 * t98;
t9 = t194 * t50 - t195 * t52 + t239 * t294 + t69 * t99 + t71 * t98;
t6 = t206 * t21 + t207 * t22 + t368 * t44;
t5 = t19 * t206 + t20 * t207 + t368 * t43;
t4 = (qJD(2) * t322 + t49) * t234 + (qJD(2) * t68 - t242 * t53 - t244 * t51 + (t242 * t70 - t244 * t72) * qJD(5)) * t235;
t3 = (qJD(2) * t323 + t48) * t234 + (qJD(2) * t67 - t242 * t52 - t244 * t50 + (t242 * t69 - t244 * t71) * qJD(5)) * t235;
t2 = qJD(5) * t252 + t10 * t207 + t206 * t9;
t1 = [-m(3) * t75 + m(4) * t47 + m(5) * t24 + m(6) * t12; (t19 * t239 - t20 * t240) * t300 + (t21 * t239 - t22 * t240) * t301 - t240 * t2 / 0.2e1 + (-t10 * t240 + t239 * t9) * t423 + ((t194 * t91 - t195 * t93) * t206 + (t194 * t90 - t195 * t92) * t207 + (t408 + (t112 * t194 - t114 * t195 - t410) * t234) * qJD(5) + (((-t22 + t395) * qJD(5) + t327) * t234 + t247) * t239) * t424 + (t239 * t7 - t240 * t8) * t425 + ((t192 * t91 + t193 * t93) * t206 + (t192 * t90 + t193 * t92) * t207 + (t409 + (t112 * t192 + t114 * t193 - t411) * t234) * qJD(5) + (((-t19 + t395) * qJD(5) + t327) * t234 + t247) * t240) * t426 + t239 * t428 + (t239 * t25 - t240 * t26) * t336 - t11 * t367 / 0.2e1 + (((-t242 * t93 - t244 * t91 + t67) * t206 + (-t242 * t92 - t244 * t90 + t68) * t207 + t46 * qJD(5)) * t235 + ((t271 + (-t112 * t244 - t114 * t242 - t263) * t235 - t325) * qJD(5) + t250) * t234) * t348 + (-t38 * (t116 * t207 + t225) - t37 * (-t116 * t206 + t223) - ((t37 * t73 - t38 * t74) * t235 + t290 * t234) * qJD(5) - (t37 * t239 + t38 * t240) * (-t215 + t333) * qJD(2) + t12 * t341 + (t17 * t283 + t38 * t267 + t12 * (t208 + t73)) * t240 + (t16 * t283 + t37 * t267 + t12 * (t209 + t74)) * t239 + (-(-t234 * t353 - t340) * qJD(2) - t353 * t377 + t449 + (t454 + t54) * t240 + (-t455 + t55) * t239) * t18) * m(6) + (t24 * t341 + (t24 * t134 + t299 * t61 + t339 * t59) * t240 + (t24 * t135 + t299 * t60 + t339 * t58) * t239 - t61 * t225 - t60 * t223 - ((t240 * t289 + t338 * t61) * t240 + (t239 * t289 + t338 * t60) * t239) * qJD(2) + (t340 * qJD(2) + t437 + t449) * t23) * m(5) + (t47 * t386 + t45 * t380 + (t101 * t334 + t105 * t345 + t47 * t133 - t45 * t154) * t240 + (t100 * t334 + t104 * t345 + t47 * t132 - t45 * t152) * t239 - (-t45 * t340 + (t101 * t344 - t240 * t288) * t240 + (t100 * t344 - t239 * t288) * t239) * qJD(2)) * m(4) + (t440 * t237 + (t443 * t239 + (-t441 - t444) * t240) * t239) * t375 + (-t444 * t237 + (t441 * t239 + (-t440 + t443) * t240) * t239) * t374 + ((-t429 * t240 + (t269 - t445) * t239 + t450) * t374 + t442 * t236 * qJD(2)) * t350 + ((t269 * t239 + (-t429 - t442) * t240 + t450) * t375 + t445 * qJD(2) * t237) * t374 / 0.2e1 + ((-t4 + t5) * t240 + (t3 + t6) * t239) * t347 + (-t302 * t62 - t303 * t75 + (t228 * t229 * t246 + t287 * t62) * t378) * m(3); m(4) * (-t104 * t240 + t105 * t239) + m(5) * (t239 * t59 - t240 * t58) + m(6) * t329; 0.2e1 * (t12 * t422 + t18 * t427) * m(6) + 0.2e1 * (t23 * t427 + t24 * t422) * m(5) + 0.2e1 * (m(5) * (qJD(2) * t23 + t239 * t58 + t240 * t59) / 0.2e1 + m(6) * (qJD(2) * t18 + t16 * t239 + t17 * t240) / 0.2e1) * t234; -t5 * t357 / 0.2e1 + t391 * t428 + (t234 * t43 + t235 * t328) * t300 + ((t239 * t8 + t240 * t7) * t235 + t253) * t425 + t234 * t6 * t350 + t2 * t392 / 0.2e1 + (t234 * t44 + t235 * t326) * t301 + ((t10 * t239 + t240 * t9) * t235 + t252) * t423 + t11 * t351 + t234 * (qJD(5) * t251 + t206 * t3 + t207 * t4) / 0.2e1 + (t234 * t46 + t235 * t325) * t336 + ((t239 * t4 + t240 * t3) * t235 + t251) * t347 + (t254 * t192 - t193 * t255 + t240 * t262) * t426 + (t194 * t254 + t195 * t255 + t239 * t262) * t424 + (t268 * t234 + (t255 * t242 - t244 * t254) * t235) * t348 + ((t16 * t73 - t17 * t74 + t37 * t54 - t38 * t55 + (t290 - (-t239 * t38 + t240 * t37) * t266) * qJD(2)) * t234 + (t38 * (-qJD(2) * t74 + t239 * t66) + t37 * (qJD(2) * t73 - t240 * t66) - t12 * t324 + t18 * (-t239 * t54 + t240 * t55) - t329 * t266) * t235 - t38 * (t183 * t207 - t368 * t83) - t37 * (-t183 * t206 + t368 * t82) - t18 * (t206 * t83 - t207 * t82)) * m(6);];
tauc = t1(:);
