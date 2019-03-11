% Calculate joint inertia matrix for
% S6RRRRPR8
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% pkin [10x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d2,d3,d4,d6]';
% m_mdh [7x1]
%   mass of all robot links (including the base)
% rSges [7x3]
%   center of mass of all robot links (in body frames)
%   rows: links of the robot (starting with base)
%   columns: x-, y-, z-coordinates
% Icges [7x6]
%   inertia of all robot links about their respective center of mass, in body frames
%   rows: links of the robot (starting with base)
%   columns: xx, yy, zz, xy, xz, yz (see inertiavector2matrix.m)
% 
% Output:
% Mq [6x6]
%   inertia matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-09 22:46
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Mq = S6RRRRPR8_inertiaJ_slag_vp1(qJ, ...
  pkin, m, rSges, Icges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(10,1),zeros(7,1),zeros(7,3),zeros(7,6)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRRRPR8_inertiaJ_slag_vp1: qJ has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S6RRRRPR8_inertiaJ_slag_vp1: pkin has to be [10x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6RRRRPR8_inertiaJ_slag_vp1: m has to be [7x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [7,3]), ...
  'S6RRRRPR8_inertiaJ_slag_vp1: rSges has to be [7x3] (double)');
assert(isreal(Icges) && all(size(Icges) == [7 6]), ...
  'S6RRRRPR8_inertiaJ_slag_vp1: Icges has to be [7x6] (double)'); 

%% Symbolic Calculation
% From inertia_joint_joint_floatb_twist_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 22:39:03
% EndTime: 2019-03-09 22:39:17
% DurationCPUTime: 7.10s
% Computational Cost: add. (17885->568), mult. (24541->782), div. (0->0), fcn. (28698->10), ass. (0->286)
t292 = qJ(3) + qJ(4);
t286 = sin(t292);
t287 = cos(t292);
t300 = cos(qJ(1));
t296 = sin(qJ(1));
t299 = cos(qJ(2));
t362 = t296 * t299;
t252 = t286 * t362 + t287 * t300;
t253 = -t286 * t300 + t287 * t362;
t295 = sin(qJ(2));
t365 = t295 * t296;
t163 = Icges(6,5) * t253 + Icges(6,6) * t365 + Icges(6,3) * t252;
t169 = Icges(5,4) * t253 - Icges(5,2) * t252 + Icges(5,6) * t365;
t418 = t163 - t169;
t361 = t299 * t300;
t254 = t286 * t361 - t296 * t287;
t255 = t286 * t296 + t287 * t361;
t364 = t295 * t300;
t164 = Icges(6,5) * t255 + Icges(6,6) * t364 + Icges(6,3) * t254;
t170 = Icges(5,4) * t255 - Icges(5,2) * t254 + Icges(5,6) * t364;
t417 = t164 - t170;
t165 = Icges(5,5) * t253 - Icges(5,6) * t252 + Icges(5,3) * t365;
t167 = Icges(6,4) * t253 + Icges(6,2) * t365 + Icges(6,6) * t252;
t416 = t165 + t167;
t166 = Icges(5,5) * t255 - Icges(5,6) * t254 + Icges(5,3) * t364;
t168 = Icges(6,4) * t255 + Icges(6,2) * t364 + Icges(6,6) * t254;
t415 = t166 + t168;
t171 = Icges(6,1) * t253 + Icges(6,4) * t365 + Icges(6,5) * t252;
t173 = Icges(5,1) * t253 - Icges(5,4) * t252 + Icges(5,5) * t365;
t414 = t171 + t173;
t172 = Icges(6,1) * t255 + Icges(6,4) * t364 + Icges(6,5) * t254;
t174 = Icges(5,1) * t255 - Icges(5,4) * t254 + Icges(5,5) * t364;
t413 = t172 + t174;
t412 = t418 * t252 + t414 * t253 + t416 * t365;
t411 = t417 * t252 + t413 * t253 + t415 * t365;
t410 = t418 * t254 + t414 * t255 + t416 * t364;
t409 = t417 * t254 + t413 * t255 + t415 * t364;
t227 = -Icges(5,6) * t299 + (Icges(5,4) * t287 - Icges(5,2) * t286) * t295;
t370 = t286 * t295;
t224 = -Icges(6,6) * t299 + (Icges(6,5) * t287 + Icges(6,3) * t286) * t295;
t228 = -Icges(6,4) * t299 + (Icges(6,1) * t287 + Icges(6,5) * t286) * t295;
t229 = -Icges(5,5) * t299 + (Icges(5,1) * t287 - Icges(5,4) * t286) * t295;
t369 = t287 * t295;
t403 = t224 * t370 + (t228 + t229) * t369;
t225 = -Icges(5,3) * t299 + (Icges(5,5) * t287 - Icges(5,6) * t286) * t295;
t226 = -Icges(6,2) * t299 + (Icges(6,4) * t287 + Icges(6,6) * t286) * t295;
t405 = -t225 - t226;
t395 = -t227 * t370 + t299 * t405 + t403;
t408 = t395 * t299;
t118 = t224 * t252 + t226 * t365 + t228 * t253;
t119 = t225 * t365 - t227 * t252 + t229 * t253;
t407 = -t118 - t119;
t120 = t224 * t254 + t226 * t364 + t228 * t255;
t121 = t225 * t364 - t227 * t254 + t229 * t255;
t406 = -t120 - t121;
t404 = Icges(3,5) * t295;
t402 = t404 / 0.2e1;
t401 = -t296 / 0.2e1;
t400 = t299 / 0.2e1;
t383 = t300 / 0.2e1;
t399 = t407 * t299 + (t296 * t412 + t411 * t300) * t295;
t398 = t406 * t299 + (t296 * t410 + t300 * t409) * t295;
t93 = -t167 * t299 + (t163 * t286 + t171 * t287) * t295;
t95 = -t165 * t299 + (-t169 * t286 + t173 * t287) * t295;
t397 = -t93 - t95;
t94 = -t168 * t299 + (t164 * t286 + t172 * t287) * t295;
t96 = -t166 * t299 + (-t170 * t286 + t174 * t287) * t295;
t396 = t94 + t96;
t293 = sin(qJ(6));
t297 = cos(qJ(6));
t188 = t252 * t297 - t253 * t293;
t189 = t252 * t293 + t253 * t297;
t125 = Icges(7,5) * t189 + Icges(7,6) * t188 - Icges(7,3) * t365;
t127 = Icges(7,4) * t189 + Icges(7,2) * t188 - Icges(7,6) * t365;
t129 = Icges(7,1) * t189 + Icges(7,4) * t188 - Icges(7,5) * t365;
t38 = -t125 * t365 + t127 * t188 + t129 * t189;
t190 = t254 * t297 - t255 * t293;
t191 = t254 * t293 + t255 * t297;
t126 = Icges(7,5) * t191 + Icges(7,6) * t190 - Icges(7,3) * t364;
t128 = Icges(7,4) * t191 + Icges(7,2) * t190 - Icges(7,6) * t364;
t130 = Icges(7,1) * t191 + Icges(7,4) * t190 - Icges(7,5) * t364;
t39 = -t126 * t365 + t128 * t188 + t130 * t189;
t16 = t296 * t39 - t300 * t38;
t394 = t411 * t296 - t300 * t412 + t16;
t40 = -t125 * t364 + t127 * t190 + t129 * t191;
t41 = -t126 * t364 + t128 * t190 + t130 * t191;
t17 = t296 * t41 - t300 * t40;
t393 = t296 * t409 - t300 * t410 + t17;
t232 = (t286 * t297 - t287 * t293) * t295;
t233 = (t286 * t293 + t287 * t297) * t295;
t159 = Icges(7,5) * t233 + Icges(7,6) * t232 + Icges(7,3) * t299;
t160 = Icges(7,4) * t233 + Icges(7,2) * t232 + Icges(7,6) * t299;
t161 = Icges(7,1) * t233 + Icges(7,4) * t232 + Icges(7,5) * t299;
t70 = -t159 * t364 + t160 * t190 + t161 * t191;
t10 = -t70 * t299 + (t296 * t40 + t300 * t41) * t295;
t69 = -t159 * t365 + t160 * t188 + t161 * t189;
t9 = -t69 * t299 + (t296 * t38 + t300 * t39) * t295;
t392 = -t10 * t364 - t9 * t365;
t55 = t125 * t299 + t127 * t232 + t129 * t233;
t56 = t126 * t299 + t128 * t232 + t130 * t233;
t391 = (t56 * t296 - t55 * t300) * t400 + t10 * t401 + t9 * t383;
t390 = t296 ^ 2;
t389 = t300 ^ 2;
t386 = t296 / 0.2e1;
t385 = -t299 / 0.2e1;
t384 = -t300 / 0.2e1;
t382 = rSges(7,3) + pkin(10);
t381 = pkin(2) * t299;
t380 = pkin(8) * t295;
t298 = cos(qJ(3));
t285 = pkin(3) * t298 + pkin(2);
t378 = -pkin(2) + t285;
t377 = rSges(6,3) * t252;
t336 = t299 * t159 + t232 * t160 + t233 * t161;
t75 = t336 * t299;
t13 = -t75 + (t55 * t296 + t56 * t300) * t295;
t375 = t299 * t13;
t373 = t300 * rSges(3,3);
t371 = Icges(3,4) * t299;
t294 = sin(qJ(3));
t244 = -Icges(4,6) * t299 + (Icges(4,4) * t298 - Icges(4,2) * t294) * t295;
t368 = t294 * t244;
t367 = t294 * t296;
t366 = t294 * t300;
t301 = -pkin(9) - pkin(8);
t363 = t295 * t301;
t355 = t191 * rSges(7,1) + t190 * rSges(7,2);
t132 = -rSges(7,3) * t364 + t355;
t239 = t255 * pkin(5);
t360 = -pkin(10) * t364 + t132 + t239;
t315 = -rSges(7,1) * t189 - rSges(7,2) * t188;
t131 = -rSges(7,3) * t365 - t315;
t162 = rSges(7,1) * t233 + rSges(7,2) * t232 + rSges(7,3) * t299;
t101 = -t299 * t131 - t162 * t365;
t175 = rSges(6,1) * t253 + rSges(6,2) * t365 + t377;
t234 = t252 * qJ(5);
t195 = pkin(4) * t253 + t234;
t179 = t195 * t364;
t358 = t175 * t364 + t179;
t177 = t255 * rSges(6,1) + rSges(6,2) * t364 + t254 * rSges(6,3);
t196 = t255 * pkin(4) + t254 * qJ(5);
t357 = -t177 - t196;
t178 = t255 * rSges(5,1) - t254 * rSges(5,2) + rSges(5,3) * t364;
t340 = t300 * t363;
t347 = pkin(3) * t367 + t285 * t361;
t307 = -t340 + t347;
t345 = pkin(2) * t361 + pkin(8) * t364;
t206 = t307 - t345;
t356 = -t178 - t206;
t256 = (pkin(4) * t287 + qJ(5) * t286) * t295;
t354 = t299 * t195 + t256 * t365;
t346 = -pkin(3) * t366 - t296 * t363;
t205 = (t378 * t299 - t380) * t296 + t346;
t223 = (pkin(8) + t301) * t299 + t378 * t295;
t353 = t299 * t205 + t223 * t365;
t316 = -rSges(5,1) * t253 + rSges(5,2) * t252;
t176 = rSges(5,3) * t365 - t316;
t231 = -rSges(5,3) * t299 + (rSges(5,1) * t287 - rSges(5,2) * t286) * t295;
t142 = t299 * t176 + t231 * t365;
t351 = -t223 - t231;
t230 = -rSges(6,2) * t299 + (rSges(6,1) * t287 + rSges(6,3) * t286) * t295;
t350 = -t230 - t256;
t251 = -rSges(4,3) * t299 + (rSges(4,1) * t298 - rSges(4,2) * t294) * t295;
t274 = t295 * pkin(2) - t299 * pkin(8);
t349 = -t251 - t274;
t348 = t390 * (t380 + t381) + t300 * t345;
t344 = t300 * pkin(1) + t296 * pkin(7);
t343 = -t13 + t408 + (t296 * t397 - t300 * t396) * t295;
t342 = -t56 / 0.2e1 - t70 / 0.2e1;
t341 = -t69 / 0.2e1 - t55 / 0.2e1;
t123 = t131 * t364;
t217 = pkin(5) * t253 - pkin(10) * t365;
t338 = t217 * t364 + t123 + t179;
t337 = -t196 - t360;
t265 = pkin(5) * t369 + pkin(10) * t299;
t335 = -t162 - t256 - t265;
t334 = -t206 + t357;
t333 = -t223 + t350;
t332 = -t274 + t351;
t263 = -t294 * t361 + t296 * t298;
t264 = t298 * t361 + t367;
t204 = t264 * rSges(4,1) + t263 * rSges(4,2) + rSges(4,3) * t364;
t290 = t300 * pkin(7);
t331 = t290 - t346;
t330 = t365 / 0.2e1;
t329 = t364 / 0.2e1;
t198 = Icges(4,5) * t264 + Icges(4,6) * t263 + Icges(4,3) * t364;
t200 = Icges(4,4) * t264 + Icges(4,2) * t263 + Icges(4,6) * t364;
t202 = Icges(4,1) * t264 + Icges(4,4) * t263 + Icges(4,5) * t364;
t108 = -t198 * t299 + (-t200 * t294 + t202 * t298) * t295;
t241 = -Icges(4,3) * t299 + (Icges(4,5) * t298 - Icges(4,6) * t294) * t295;
t247 = -Icges(4,5) * t299 + (Icges(4,1) * t298 - Icges(4,4) * t294) * t295;
t134 = t241 * t364 + t244 * t263 + t247 * t264;
t328 = t108 / 0.2e1 + t134 / 0.2e1;
t261 = -t294 * t362 - t298 * t300;
t262 = t298 * t362 - t366;
t197 = Icges(4,5) * t262 + Icges(4,6) * t261 + Icges(4,3) * t365;
t199 = Icges(4,4) * t262 + Icges(4,2) * t261 + Icges(4,6) * t365;
t201 = Icges(4,1) * t262 + Icges(4,4) * t261 + Icges(4,5) * t365;
t107 = -t197 * t299 + (-t199 * t294 + t201 * t298) * t295;
t133 = t241 * t365 + t244 * t261 + t247 * t262;
t327 = t133 / 0.2e1 + t107 / 0.2e1;
t326 = -t285 * t299 - pkin(1);
t325 = -t206 + t337;
t324 = -t223 + t335;
t323 = t296 * t205 + t300 * t206 + t348;
t105 = t299 * t175 + t230 * t365 + t354;
t322 = -t274 + t333;
t321 = -t234 + t331;
t320 = t375 + t392;
t319 = t364 * t398 + t365 * t399 - t392;
t318 = rSges(3,1) * t299 - rSges(3,2) * t295;
t317 = -rSges(4,1) * t262 - rSges(4,2) * t261;
t314 = -t274 + t324;
t312 = -Icges(3,2) * t295 + t371;
t311 = Icges(3,5) * t299 - Icges(3,6) * t295;
t308 = rSges(3,1) * t361 - rSges(3,2) * t364 + t296 * rSges(3,3);
t306 = t296 * t195 + t300 * t196 + t323;
t63 = t299 * t217 + t265 * t365 - t101 + t354;
t305 = t196 + t344 + t347;
t304 = t343 * t299 + t319;
t303 = -t75 + (t55 + t69 - t397 - t407) * t330 + (t56 + t70 + t396 - t406) * t329;
t302 = -t391 + t398 * t386 + (t296 * t396 + t300 * t397) * t385 + t399 * t384 + t394 * t330 + t393 * t329;
t273 = rSges(2,1) * t300 - rSges(2,2) * t296;
t272 = -rSges(2,1) * t296 - rSges(2,2) * t300;
t271 = rSges(3,1) * t295 + rSges(3,2) * t299;
t268 = Icges(3,6) * t299 + t404;
t243 = Icges(3,3) * t296 + t311 * t300;
t242 = -Icges(3,3) * t300 + t311 * t296;
t221 = t295 * t298 * t247;
t220 = t308 + t344;
t219 = t373 + t290 + (-pkin(1) - t318) * t296;
t208 = t349 * t300;
t207 = t349 * t296;
t203 = rSges(4,3) * t365 - t317;
t187 = t300 * t308 + (t318 * t296 - t373) * t296;
t181 = t205 * t364;
t156 = t176 * t364;
t152 = t204 + t344 + t345;
t151 = t290 + (-t381 - pkin(1) + (-rSges(4,3) - pkin(8)) * t295) * t296 + t317;
t150 = t332 * t300;
t149 = t332 * t296;
t148 = -t204 * t299 - t251 * t364;
t147 = t203 * t299 + t251 * t365;
t146 = -t299 * t241 - t295 * t368 + t221;
t143 = -t178 * t299 - t231 * t364;
t141 = t307 + t178 + t344;
t140 = (-rSges(5,3) * t295 + t326) * t296 + t316 + t331;
t139 = t322 * t300;
t138 = t322 * t296;
t135 = (t203 * t300 - t204 * t296) * t295;
t122 = -t178 * t365 + t156;
t113 = t203 * t296 + t204 * t300 + t348;
t112 = t177 + t305 - t340;
t111 = -t377 + (-rSges(6,1) - pkin(4)) * t253 + (-rSges(6,2) * t295 + t326) * t296 + t321;
t110 = t314 * t300;
t109 = t314 * t296;
t106 = t357 * t299 + t350 * t364;
t104 = t356 * t299 + t351 * t364;
t103 = t142 + t353;
t102 = t132 * t299 + t162 * t364;
t100 = t198 * t364 + t200 * t263 + t202 * t264;
t99 = t197 * t364 + t199 * t263 + t201 * t264;
t98 = t198 * t365 + t200 * t261 + t202 * t262;
t97 = t197 * t365 + t199 * t261 + t201 * t262;
t80 = t356 * t365 + t156 + t181;
t78 = t239 + (-t301 - t382) * t364 + t305 + t355;
t77 = (-pkin(4) - pkin(5)) * t253 + (t382 * t295 + t326) * t296 + t315 + t321;
t76 = t357 * t365 + t358;
t74 = t132 * t365 - t123;
t73 = t334 * t299 + t333 * t364;
t72 = t105 + t353;
t71 = t176 * t296 + t178 * t300 + t323;
t64 = t337 * t299 + t335 * t364;
t62 = t334 * t365 + t181 + t358;
t61 = t175 * t296 + t177 * t300 + t306;
t60 = t100 * t296 - t300 * t99;
t59 = t296 * t98 - t300 * t97;
t54 = t325 * t299 + t324 * t364;
t53 = t63 + t353;
t50 = t337 * t365 + t338;
t37 = -t134 * t299 + (t100 * t300 + t296 * t99) * t295;
t36 = -t133 * t299 + (t296 * t97 + t300 * t98) * t295;
t33 = t325 * t365 + t181 + t338;
t32 = t360 * t300 + (t131 + t217) * t296 + t306;
t1 = [Icges(2,3) + t221 + (Icges(3,1) * t295 - t286 * t227 - t368 + t371) * t295 + (Icges(3,4) * t295 + Icges(3,2) * t299 - t241 + t405) * t299 + m(7) * (t77 ^ 2 + t78 ^ 2) + m(6) * (t111 ^ 2 + t112 ^ 2) + m(5) * (t140 ^ 2 + t141 ^ 2) + m(4) * (t151 ^ 2 + t152 ^ 2) + m(3) * (t219 ^ 2 + t220 ^ 2) + m(2) * (t272 ^ 2 + t273 ^ 2) + t336 + t403; (-t95 / 0.2e1 - t93 / 0.2e1 - t118 / 0.2e1 - t119 / 0.2e1 + (-Icges(3,6) * t300 + t312 * t296) * t385 + t300 * t402 + t268 * t383 - t327 + t341) * t300 + (t96 / 0.2e1 + t94 / 0.2e1 + t120 / 0.2e1 + t121 / 0.2e1 + (Icges(3,6) * t296 + t312 * t300) * t400 + t296 * t402 + t268 * t386 + t328 - t342) * t296 + m(4) * (t151 * t208 + t152 * t207) + m(7) * (t109 * t78 + t110 * t77) + m(6) * (t111 * t139 + t112 * t138) + m(5) * (t140 * t150 + t141 * t149) + m(3) * (-t219 * t300 - t220 * t296) * t271; m(7) * (t109 ^ 2 + t110 ^ 2 + t32 ^ 2) + m(6) * (t138 ^ 2 + t139 ^ 2 + t61 ^ 2) + m(5) * (t149 ^ 2 + t150 ^ 2 + t71 ^ 2) + m(4) * (t113 ^ 2 + t207 ^ 2 + t208 ^ 2) + m(3) * (t187 ^ 2 + (t389 + t390) * t271 ^ 2) + (-t389 * t242 - t394 - t59) * t300 + (t390 * t243 + t60 + (-t296 * t242 + t300 * t243) * t300 + t393) * t296; t303 + (t327 * t296 + t328 * t300) * t295 + m(7) * (t53 * t77 + t54 * t78) + m(6) * (t111 * t72 + t112 * t73) + m(5) * (t103 * t140 + t104 * t141) + m(4) * (t147 * t151 + t148 * t152) + (-t146 - t395) * t299; t302 + (t60 * t383 + t59 * t386) * t295 + m(7) * (t109 * t54 + t110 * t53 + t33 * t32) + m(6) * (t138 * t73 + t139 * t72 + t62 * t61) + m(5) * (t103 * t150 + t104 * t149 + t71 * t80) + m(4) * (t113 * t135 + t147 * t208 + t148 * t207) + t37 * t386 + t36 * t384 + (-t107 * t300 + t108 * t296) * t385; (t296 * t36 + t300 * t37) * t295 + (t146 * t299 + (-t107 * t296 - t108 * t300) * t295 + t343) * t299 + m(7) * (t33 ^ 2 + t53 ^ 2 + t54 ^ 2) + m(6) * (t62 ^ 2 + t72 ^ 2 + t73 ^ 2) + m(5) * (t103 ^ 2 + t104 ^ 2 + t80 ^ 2) + m(4) * (t135 ^ 2 + t147 ^ 2 + t148 ^ 2) + t319; t303 + m(7) * (t63 * t77 + t64 * t78) + m(6) * (t105 * t111 + t106 * t112) + m(5) * (t140 * t142 + t141 * t143) - t408; m(7) * (t109 * t64 + t110 * t63 + t50 * t32) + m(6) * (t105 * t139 + t106 * t138 + t61 * t76) + m(5) * (t122 * t71 + t142 * t150 + t143 * t149) + t302; m(7) * (t33 * t50 + t53 * t63 + t54 * t64) + m(6) * (t105 * t72 + t106 * t73 + t62 * t76) + m(5) * (t103 * t142 + t104 * t143 + t122 * t80) + t304; m(7) * (t50 ^ 2 + t63 ^ 2 + t64 ^ 2) + m(6) * (t105 ^ 2 + t106 ^ 2 + t76 ^ 2) + m(5) * (t122 ^ 2 + t142 ^ 2 + t143 ^ 2) + t304; m(7) * (t252 * t78 + t254 * t77) + m(6) * (t111 * t254 + t112 * t252); m(7) * (t109 * t252 + t110 * t254 + t32 * t370) + m(6) * (t138 * t252 + t139 * t254 + t61 * t370); m(7) * (t252 * t54 + t254 * t53 + t33 * t370) + m(6) * (t252 * t73 + t254 * t72 + t62 * t370); m(7) * (t252 * t64 + t254 * t63 + t50 * t370) + m(6) * (t105 * t254 + t106 * t252 + t76 * t370); 0.2e1 * (m(6) / 0.2e1 + m(7) / 0.2e1) * (t286 ^ 2 * t295 ^ 2 + t252 ^ 2 + t254 ^ 2); m(7) * (t101 * t77 + t102 * t78) + t75 + (t341 * t296 + t342 * t300) * t295; m(7) * (t101 * t110 + t102 * t109 + t32 * t74) + (t16 * t401 + t17 * t384) * t295 + t391; m(7) * (t101 * t53 + t102 * t54 + t33 * t74) + t320; m(7) * (t101 * t63 + t102 * t64 + t50 * t74) + t320; m(7) * (t101 * t254 + t102 * t252 + t74 * t370); -t375 + m(7) * (t101 ^ 2 + t102 ^ 2 + t74 ^ 2) + (t10 * t300 + t296 * t9) * t295;];
%% Postprocessing: Reshape Output
% From vec2symmat_6_matlab.m
res = [t1(1) t1(2) t1(4) t1(7) t1(11) t1(16); t1(2) t1(3) t1(5) t1(8) t1(12) t1(17); t1(4) t1(5) t1(6) t1(9) t1(13) t1(18); t1(7) t1(8) t1(9) t1(10) t1(14) t1(19); t1(11) t1(12) t1(13) t1(14) t1(15) t1(20); t1(16) t1(17) t1(18) t1(19) t1(20) t1(21);];
Mq  = res;
