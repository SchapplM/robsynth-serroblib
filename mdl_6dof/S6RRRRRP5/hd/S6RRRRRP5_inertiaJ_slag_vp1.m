% Calculate joint inertia matrix for
% S6RRRRRP5
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% pkin [10x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d2,d3,d4,d5]';
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
% Datum: 2019-03-10 01:25
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Mq = S6RRRRRP5_inertiaJ_slag_vp1(qJ, ...
  pkin, m, rSges, Icges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(10,1),zeros(7,1),zeros(7,3),zeros(7,6)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRRRRP5_inertiaJ_slag_vp1: qJ has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S6RRRRRP5_inertiaJ_slag_vp1: pkin has to be [10x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6RRRRRP5_inertiaJ_slag_vp1: m has to be [7x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [7,3]), ...
  'S6RRRRRP5_inertiaJ_slag_vp1: rSges has to be [7x3] (double)');
assert(isreal(Icges) && all(size(Icges) == [7 6]), ...
  'S6RRRRRP5_inertiaJ_slag_vp1: Icges has to be [7x6] (double)'); 

%% Symbolic Calculation
% From inertia_joint_joint_floatb_twist_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-10 01:18:12
% EndTime: 2019-03-10 01:18:28
% DurationCPUTime: 7.58s
% Computational Cost: add. (21822->582), mult. (22285->788), div. (0->0), fcn. (24059->10), ass. (0->288)
t298 = qJ(3) + qJ(4);
t287 = qJ(5) + t298;
t281 = sin(t287);
t282 = cos(t287);
t304 = cos(qJ(1));
t301 = sin(qJ(1));
t303 = cos(qJ(2));
t370 = t301 * t303;
t224 = -t281 * t370 - t282 * t304;
t225 = -t281 * t304 + t282 * t370;
t300 = sin(qJ(2));
t373 = t300 * t301;
t150 = Icges(7,5) * t225 + Icges(7,6) * t224 + Icges(7,3) * t373;
t152 = Icges(6,5) * t225 + Icges(6,6) * t224 + Icges(6,3) * t373;
t420 = t150 + t152;
t369 = t303 * t304;
t226 = -t281 * t369 + t282 * t301;
t227 = t281 * t301 + t282 * t369;
t372 = t300 * t304;
t151 = Icges(7,5) * t227 + Icges(7,6) * t226 + Icges(7,3) * t372;
t153 = Icges(6,5) * t227 + Icges(6,6) * t226 + Icges(6,3) * t372;
t419 = t151 + t153;
t154 = Icges(7,4) * t225 + Icges(7,2) * t224 + Icges(7,6) * t373;
t156 = Icges(6,4) * t225 + Icges(6,2) * t224 + Icges(6,6) * t373;
t418 = t154 + t156;
t155 = Icges(7,4) * t227 + Icges(7,2) * t226 + Icges(7,6) * t372;
t157 = Icges(6,4) * t227 + Icges(6,2) * t226 + Icges(6,6) * t372;
t417 = t155 + t157;
t158 = Icges(7,1) * t225 + Icges(7,4) * t224 + Icges(7,5) * t373;
t160 = Icges(6,1) * t225 + Icges(6,4) * t224 + Icges(6,5) * t373;
t416 = t158 + t160;
t159 = Icges(7,1) * t227 + Icges(7,4) * t226 + Icges(7,5) * t372;
t161 = Icges(6,1) * t227 + Icges(6,4) * t226 + Icges(6,5) * t372;
t415 = t159 + t161;
t414 = t224 * t418 + t225 * t416 + t373 * t420;
t413 = t224 * t417 + t225 * t415 + t373 * t419;
t412 = t226 * t418 + t227 * t416 + t372 * t420;
t411 = t226 * t417 + t227 * t415 + t372 * t419;
t204 = -Icges(7,3) * t303 + (Icges(7,5) * t282 - Icges(7,6) * t281) * t300;
t206 = -Icges(7,6) * t303 + (Icges(7,4) * t282 - Icges(7,2) * t281) * t300;
t208 = -Icges(7,5) * t303 + (Icges(7,1) * t282 - Icges(7,4) * t281) * t300;
t98 = t204 * t373 + t206 * t224 + t208 * t225;
t205 = -Icges(6,3) * t303 + (Icges(6,5) * t282 - Icges(6,6) * t281) * t300;
t207 = -Icges(6,6) * t303 + (Icges(6,4) * t282 - Icges(6,2) * t281) * t300;
t209 = -Icges(6,5) * t303 + (Icges(6,1) * t282 - Icges(6,4) * t281) * t300;
t99 = t205 * t373 + t207 * t224 + t209 * t225;
t410 = -t98 - t99;
t404 = (t208 + t209) * t282 * t300;
t406 = (-t206 - t207) * t281;
t407 = -t204 - t205;
t395 = t300 * t406 + t303 * t407 + t404;
t409 = t395 * t303;
t100 = t204 * t372 + t206 * t226 + t208 * t227;
t101 = t205 * t372 + t207 * t226 + t209 * t227;
t408 = -t100 - t101;
t405 = Icges(3,5) * t300;
t284 = sin(t298);
t299 = sin(qJ(3));
t259 = pkin(3) * t299 + pkin(4) * t284;
t245 = pkin(5) * t281 + t259;
t403 = -rSges(7,1) * t225 - rSges(7,2) * t224 + t245 * t304;
t402 = t405 / 0.2e1;
t401 = t410 * t303 + (t301 * t414 + t304 * t413) * t300;
t400 = t408 * t303 + (t301 * t412 + t304 * t411) * t300;
t399 = t301 * t413 - t304 * t414;
t398 = t301 * t411 - t304 * t412;
t76 = -t303 * t150 + (-t154 * t281 + t158 * t282) * t300;
t78 = -t303 * t152 + (-t156 * t281 + t160 * t282) * t300;
t397 = -t76 - t78;
t77 = -t303 * t151 + (-t155 * t281 + t159 * t282) * t300;
t79 = -t303 * t153 + (-t157 * t281 + t161 * t282) * t300;
t396 = t77 + t79;
t305 = -pkin(9) - pkin(8);
t294 = -pkin(10) + t305;
t286 = -qJ(6) + t294;
t352 = t259 * t304 + t294 * t373;
t302 = cos(qJ(3));
t283 = pkin(3) * t302 + pkin(2);
t285 = cos(t298);
t257 = pkin(4) * t285 + t283;
t239 = pkin(5) * t282 + t257;
t356 = t239 - t257;
t394 = (-t286 * t300 + t303 * t356) * t301 + t352 + rSges(7,3) * t373 - t403;
t347 = t286 - t294;
t393 = (t347 - rSges(7,3)) * t303 + (rSges(7,1) * t282 - rSges(7,2) * t281 + t356) * t300;
t392 = rSges(7,1) * t227 + rSges(7,2) * t226 + rSges(7,3) * t372 + t239 * t369 + t245 * t301;
t296 = t301 ^ 2;
t297 = t304 ^ 2;
t391 = t301 / 0.2e1;
t390 = -t303 / 0.2e1;
t389 = -t304 / 0.2e1;
t388 = t304 / 0.2e1;
t387 = pkin(2) * t303;
t386 = pkin(8) * t300;
t385 = -pkin(2) + t283;
t384 = t409 + (t301 * t397 - t304 * t396) * t300;
t383 = t304 * rSges(3,3);
t381 = Icges(3,4) * t303;
t216 = -Icges(5,6) * t303 + (Icges(5,4) * t285 - Icges(5,2) * t284) * t300;
t377 = t284 * t216;
t233 = -Icges(4,6) * t303 + (Icges(4,4) * t302 - Icges(4,2) * t299) * t300;
t376 = t299 * t233;
t375 = t299 * t301;
t374 = t299 * t304;
t371 = t300 * t305;
t368 = t394 * t372;
t354 = t257 * t369 + t259 * t301;
t366 = -t347 * t372 - t354 + t392;
t349 = pkin(3) * t374 + t301 * t371;
t351 = t257 - t283;
t138 = t351 * t370 + t349 - t352;
t131 = t138 * t372;
t322 = -rSges(6,1) * t225 - rSges(6,2) * t224;
t163 = rSges(6,3) * t373 - t322;
t145 = t163 * t372;
t365 = t131 + t145;
t345 = t294 - t305;
t192 = t300 * t351 + t303 * t345;
t364 = t138 * t303 + t192 * t373;
t350 = -pkin(3) * t375 - t283 * t369;
t139 = -t345 * t372 + t350 + t354;
t165 = rSges(6,1) * t227 + rSges(6,2) * t226 + rSges(6,3) * t372;
t363 = -t139 - t165;
t243 = -t284 * t369 + t285 * t301;
t244 = t284 * t301 + t285 * t369;
t175 = rSges(5,1) * t244 + rSges(5,2) * t243 + rSges(5,3) * t372;
t313 = -t304 * t371 - t350;
t348 = pkin(2) * t369 + pkin(8) * t372;
t191 = t313 - t348;
t361 = -t175 - t191;
t190 = (t303 * t385 - t386) * t301 - t349;
t214 = (pkin(8) + t305) * t303 + t385 * t300;
t360 = t190 * t303 + t214 * t373;
t212 = -t303 * rSges(6,3) + (rSges(6,1) * t282 - rSges(6,2) * t281) * t300;
t359 = -t192 - t212;
t125 = t163 * t303 + t212 * t373;
t241 = -t284 * t370 - t285 * t304;
t242 = -t284 * t304 + t285 * t370;
t323 = -rSges(5,1) * t242 - rSges(5,2) * t241;
t174 = rSges(5,3) * t373 - t323;
t222 = -t303 * rSges(5,3) + (rSges(5,1) * t285 - rSges(5,2) * t284) * t300;
t129 = t174 * t303 + t222 * t373;
t357 = -t214 - t222;
t240 = -t303 * rSges(4,3) + (rSges(4,1) * t302 - rSges(4,2) * t299) * t300;
t268 = t300 * pkin(2) - t303 * pkin(8);
t355 = -t240 - t268;
t353 = t296 * (t386 + t387) + t304 * t348;
t346 = pkin(1) * t304 + pkin(7) * t301;
t344 = t296 + t297;
t217 = -Icges(5,5) * t303 + (Icges(5,1) * t285 - Icges(5,4) * t284) * t300;
t197 = t300 * t285 * t217;
t215 = -Icges(5,3) * t303 + (Icges(5,5) * t285 - Icges(5,6) * t284) * t300;
t124 = -t303 * t215 - t300 * t377 + t197;
t168 = Icges(5,5) * t242 + Icges(5,6) * t241 + Icges(5,3) * t373;
t170 = Icges(5,4) * t242 + Icges(5,2) * t241 + Icges(5,6) * t373;
t172 = Icges(5,1) * t242 + Icges(5,4) * t241 + Icges(5,5) * t373;
t82 = -t303 * t168 + (-t170 * t284 + t172 * t285) * t300;
t169 = Icges(5,5) * t244 + Icges(5,6) * t243 + Icges(5,3) * t372;
t171 = Icges(5,4) * t244 + Icges(5,2) * t243 + Icges(5,6) * t372;
t173 = Icges(5,1) * t244 + Icges(5,4) * t243 + Icges(5,5) * t372;
t83 = -t303 * t169 + (-t171 * t284 + t173 * t285) * t300;
t343 = t124 * t303 - (t301 * t82 + t304 * t83) * t300 + t384;
t342 = t372 * t400 + t373 * t401;
t230 = -Icges(4,3) * t303 + (Icges(4,5) * t302 - Icges(4,6) * t299) * t300;
t236 = -Icges(4,5) * t303 + (Icges(4,1) * t302 - Icges(4,4) * t299) * t300;
t253 = -t299 * t369 + t301 * t302;
t254 = t302 * t369 + t375;
t122 = t230 * t372 + t233 * t253 + t236 * t254;
t183 = Icges(4,5) * t254 + Icges(4,6) * t253 + Icges(4,3) * t372;
t185 = Icges(4,4) * t254 + Icges(4,2) * t253 + Icges(4,6) * t372;
t187 = Icges(4,1) * t254 + Icges(4,4) * t253 + Icges(4,5) * t372;
t93 = -t303 * t183 + (-t185 * t299 + t187 * t302) * t300;
t341 = t122 / 0.2e1 + t93 / 0.2e1;
t251 = -t299 * t370 - t302 * t304;
t252 = t302 * t370 - t374;
t121 = t230 * t373 + t233 * t251 + t236 * t252;
t182 = Icges(4,5) * t252 + Icges(4,6) * t251 + Icges(4,3) * t373;
t184 = Icges(4,4) * t252 + Icges(4,2) * t251 + Icges(4,6) * t373;
t186 = Icges(4,1) * t252 + Icges(4,4) * t251 + Icges(4,5) * t373;
t92 = -t303 * t182 + (-t184 * t299 + t186 * t302) * t300;
t340 = t92 / 0.2e1 + t121 / 0.2e1;
t339 = t131 + t368;
t338 = -t124 - t395;
t337 = -t139 - t366;
t336 = -t191 + t363;
t335 = -t192 - t393;
t334 = -t214 + t359;
t333 = -t268 + t357;
t189 = rSges(4,1) * t254 + rSges(4,2) * t253 + rSges(4,3) * t372;
t332 = t373 / 0.2e1;
t331 = t372 / 0.2e1;
t330 = -t191 + t337;
t54 = t303 * t394 + t373 * t393;
t329 = -t214 + t335;
t328 = t190 * t301 + t191 * t304 + t353;
t66 = t125 + t364;
t327 = -t268 + t334;
t108 = t215 * t373 + t216 * t241 + t217 * t242;
t68 = t168 * t373 + t170 * t241 + t172 * t242;
t69 = t169 * t373 + t171 * t241 + t173 * t242;
t19 = -t108 * t303 + (t301 * t68 + t304 * t69) * t300;
t109 = t215 * t372 + t216 * t243 + t217 * t244;
t70 = t168 * t372 + t170 * t243 + t172 * t244;
t71 = t169 * t372 + t171 * t243 + t173 * t244;
t20 = -t109 * t303 + (t301 * t70 + t304 * t71) * t300;
t326 = t19 * t373 + t20 * t372 + t342;
t325 = rSges(3,1) * t303 - rSges(3,2) * t300;
t324 = -rSges(4,1) * t252 - rSges(4,2) * t251;
t320 = -t268 + t329;
t318 = -Icges(3,2) * t300 + t381;
t317 = Icges(3,5) * t303 - Icges(3,6) * t300;
t314 = rSges(3,1) * t369 - rSges(3,2) * t372 + rSges(3,3) * t301;
t312 = t138 * t301 + t139 * t304 + t328;
t44 = t54 + t364;
t311 = t303 * t384 + t342;
t310 = (-t397 - t410) * t332 + (t396 - t408) * t331;
t309 = t400 * t391 + (t301 * t396 + t304 * t397) * t390 + t401 * t389 + t399 * t332 + t398 * t331;
t308 = t303 * t343 + t326;
t307 = t310 + (t108 + t82) * t332 + (t109 + t83) * t331;
t39 = t301 * t69 - t304 * t68;
t40 = t301 * t71 - t304 * t70;
t306 = t19 * t389 + t20 * t391 + t40 * t331 + t39 * t332 + t309 + (t83 * t301 - t82 * t304) * t390;
t292 = t304 * pkin(7);
t266 = rSges(2,1) * t304 - rSges(2,2) * t301;
t265 = -rSges(2,1) * t301 - rSges(2,2) * t304;
t264 = rSges(3,1) * t300 + rSges(3,2) * t303;
t261 = Icges(3,6) * t303 + t405;
t232 = Icges(3,3) * t301 + t304 * t317;
t231 = -Icges(3,3) * t304 + t301 * t317;
t210 = t300 * t302 * t236;
t203 = t314 + t346;
t202 = t383 + t292 + (-pkin(1) - t325) * t301;
t194 = t355 * t304;
t193 = t355 * t301;
t188 = rSges(4,3) * t373 - t324;
t177 = t304 * t314 + (t301 * t325 - t383) * t301;
t176 = t190 * t372;
t149 = t174 * t372;
t143 = t189 + t346 + t348;
t142 = t292 + (-t387 - pkin(1) + (-rSges(4,3) - pkin(8)) * t300) * t301 + t324;
t141 = t333 * t304;
t140 = t333 * t301;
t137 = -t189 * t303 - t240 * t372;
t136 = t188 * t303 + t240 * t373;
t132 = -t303 * t230 - t300 * t376 + t210;
t130 = -t175 * t303 - t222 * t372;
t128 = t313 + t175 + t346;
t127 = t292 + (-rSges(5,3) * t300 - t283 * t303 - pkin(1)) * t301 + t323 + t349;
t126 = -t165 * t303 - t212 * t372;
t123 = (t188 * t304 - t189 * t301) * t300;
t115 = -t294 * t372 + t165 + t346 + t354;
t114 = t292 + (-rSges(6,3) * t300 - t257 * t303 - pkin(1)) * t301 + t322 + t352;
t112 = t327 * t304;
t111 = t327 * t301;
t110 = -t175 * t373 + t149;
t107 = -t286 * t372 + t346 + t392;
t106 = t292 + (-t239 * t303 - pkin(1) + (-rSges(7,3) + t286) * t300) * t301 + t403;
t105 = -t165 * t373 + t145;
t102 = t188 * t301 + t189 * t304 + t353;
t91 = t303 * t361 + t357 * t372;
t90 = t129 + t360;
t89 = t183 * t372 + t185 * t253 + t187 * t254;
t88 = t182 * t372 + t184 * t253 + t186 * t254;
t87 = t183 * t373 + t185 * t251 + t187 * t252;
t86 = t182 * t373 + t184 * t251 + t186 * t252;
t85 = t320 * t304;
t84 = t320 * t301;
t67 = t303 * t363 + t359 * t372;
t65 = t361 * t373 + t149 + t176;
t56 = t174 * t301 + t175 * t304 + t328;
t55 = -t303 * t366 - t372 * t393;
t53 = t363 * t373 + t365;
t52 = t303 * t336 + t334 * t372;
t51 = t66 + t360;
t50 = -t366 * t373 + t368;
t49 = t301 * t89 - t304 * t88;
t48 = t301 * t87 - t304 * t86;
t46 = t336 * t373 + t176 + t365;
t45 = t303 * t337 + t335 * t372;
t41 = t163 * t301 + t165 * t304 + t312;
t32 = t303 * t330 + t329 * t372;
t31 = t44 + t360;
t30 = -t122 * t303 + (t301 * t88 + t304 * t89) * t300;
t29 = -t121 * t303 + (t301 * t86 + t304 * t87) * t300;
t23 = t337 * t373 + t339;
t10 = t330 * t373 + t176 + t339;
t1 = t301 * t394 + t304 * t366 + t312;
t2 = [Icges(2,3) + t197 + t210 + (Icges(3,4) * t300 + Icges(3,2) * t303 - t215 - t230 + t407) * t303 + (Icges(3,1) * t300 - t376 - t377 + t381 + t406) * t300 + m(7) * (t106 ^ 2 + t107 ^ 2) + m(6) * (t114 ^ 2 + t115 ^ 2) + m(5) * (t127 ^ 2 + t128 ^ 2) + m(4) * (t142 ^ 2 + t143 ^ 2) + m(3) * (t202 ^ 2 + t203 ^ 2) + m(2) * (t265 ^ 2 + t266 ^ 2) + t404; (-t82 / 0.2e1 - t78 / 0.2e1 - t76 / 0.2e1 - t98 / 0.2e1 - t99 / 0.2e1 - t108 / 0.2e1 + (-Icges(3,6) * t304 + t301 * t318) * t390 + t304 * t402 + t261 * t388 - t340) * t304 + (t83 / 0.2e1 + t79 / 0.2e1 + t77 / 0.2e1 + t100 / 0.2e1 + t101 / 0.2e1 + t109 / 0.2e1 + t303 * (Icges(3,6) * t301 + t304 * t318) / 0.2e1 + t301 * t402 + t261 * t391 + t341) * t301 + m(7) * (t106 * t85 + t107 * t84) + m(6) * (t111 * t115 + t112 * t114) + m(5) * (t127 * t141 + t128 * t140) + m(4) * (t142 * t194 + t143 * t193) + m(3) * (-t202 * t304 - t203 * t301) * t264; m(7) * (t1 ^ 2 + t84 ^ 2 + t85 ^ 2) + m(6) * (t111 ^ 2 + t112 ^ 2 + t41 ^ 2) + m(5) * (t140 ^ 2 + t141 ^ 2 + t56 ^ 2) + m(4) * (t102 ^ 2 + t193 ^ 2 + t194 ^ 2) + m(3) * (t264 ^ 2 * t344 + t177 ^ 2) + (-t297 * t231 - t39 - t399 - t48) * t304 + (t296 * t232 + t40 + t49 + (-t301 * t231 + t304 * t232) * t304 + t398) * t301; (t301 * t340 + t304 * t341) * t300 + m(7) * (t106 * t31 + t107 * t32) + m(6) * (t114 * t51 + t115 * t52) + m(5) * (t127 * t90 + t128 * t91) + m(4) * (t136 * t142 + t137 * t143) + (-t132 + t338) * t303 + t307; t306 + (t388 * t49 + t391 * t48) * t300 + m(7) * (t1 * t10 + t31 * t85 + t32 * t84) + m(6) * (t111 * t52 + t112 * t51 + t41 * t46) + m(5) * (t140 * t91 + t141 * t90 + t56 * t65) + m(4) * (t102 * t123 + t136 * t194 + t137 * t193) + t30 * t391 + (t93 * t301 - t92 * t304) * t390 + t29 * t389; (t301 * t29 + t304 * t30) * t300 + (t132 * t303 + (-t301 * t92 - t304 * t93) * t300 + t343) * t303 + m(5) * (t65 ^ 2 + t90 ^ 2 + t91 ^ 2) + m(4) * (t123 ^ 2 + t136 ^ 2 + t137 ^ 2) + m(7) * (t10 ^ 2 + t31 ^ 2 + t32 ^ 2) + m(6) * (t46 ^ 2 + t51 ^ 2 + t52 ^ 2) + t326; m(7) * (t106 * t44 + t107 * t45) + m(6) * (t114 * t66 + t115 * t67) + m(5) * (t127 * t129 + t128 * t130) + t338 * t303 + t307; t306 + m(7) * (t1 * t23 + t44 * t85 + t45 * t84) + m(6) * (t111 * t67 + t112 * t66 + t41 * t53) + m(5) * (t110 * t56 + t129 * t141 + t130 * t140); m(7) * (t10 * t23 + t31 * t44 + t32 * t45) + m(6) * (t46 * t53 + t51 * t66 + t52 * t67) + m(5) * (t110 * t65 + t129 * t90 + t130 * t91) + t308; m(7) * (t23 ^ 2 + t44 ^ 2 + t45 ^ 2) + m(6) * (t53 ^ 2 + t66 ^ 2 + t67 ^ 2) + m(5) * (t110 ^ 2 + t129 ^ 2 + t130 ^ 2) + t308; -t409 + m(7) * (t106 * t54 + t107 * t55) + m(6) * (t114 * t125 + t115 * t126) + t310; m(7) * (t1 * t50 + t54 * t85 + t55 * t84) + m(6) * (t105 * t41 + t111 * t126 + t112 * t125) + t309; m(7) * (t10 * t50 + t31 * t54 + t32 * t55) + m(6) * (t105 * t46 + t125 * t51 + t126 * t52) + t311; m(7) * (t23 * t50 + t44 * t54 + t45 * t55) + m(6) * (t105 * t53 + t125 * t66 + t126 * t67) + t311; m(7) * (t50 ^ 2 + t54 ^ 2 + t55 ^ 2) + m(6) * (t105 ^ 2 + t125 ^ 2 + t126 ^ 2) + t311; m(7) * (t106 * t304 + t107 * t301) * t300; m(7) * (-t303 * t1 + (t301 * t84 + t304 * t85) * t300); m(7) * (-t303 * t10 + (t301 * t32 + t304 * t31) * t300); m(7) * (-t303 * t23 + (t301 * t45 + t304 * t44) * t300); m(7) * (-t303 * t50 + (t301 * t55 + t304 * t54) * t300); m(7) * (t300 ^ 2 * t344 + t303 ^ 2);];
%% Postprocessing: Reshape Output
% From vec2symmat_6_matlab.m
res = [t2(1) t2(2) t2(4) t2(7) t2(11) t2(16); t2(2) t2(3) t2(5) t2(8) t2(12) t2(17); t2(4) t2(5) t2(6) t2(9) t2(13) t2(18); t2(7) t2(8) t2(9) t2(10) t2(14) t2(19); t2(11) t2(12) t2(13) t2(14) t2(15) t2(20); t2(16) t2(17) t2(18) t2(19) t2(20) t2(21);];
Mq  = res;
