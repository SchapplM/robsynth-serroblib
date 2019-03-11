% Calculate joint inertia matrix for
% S6RRPRRR4
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% pkin [12x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,alpha2,d1,d2,d4,d5,d6,theta3]';
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
% Datum: 2019-03-09 13:37
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Mq = S6RRPRRR4_inertiaJ_slag_vp1(qJ, ...
  pkin, m, rSges, Icges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(12,1),zeros(7,1),zeros(7,3),zeros(7,6)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRPRRR4_inertiaJ_slag_vp1: qJ has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [12 1]), ...
  'S6RRPRRR4_inertiaJ_slag_vp1: pkin has to be [12x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6RRPRRR4_inertiaJ_slag_vp1: m has to be [7x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [7,3]), ...
  'S6RRPRRR4_inertiaJ_slag_vp1: rSges has to be [7x3] (double)');
assert(isreal(Icges) && all(size(Icges) == [7 6]), ...
  'S6RRPRRR4_inertiaJ_slag_vp1: Icges has to be [7x6] (double)'); 

%% Symbolic Calculation
% From inertia_joint_joint_floatb_twist_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 13:27:52
% EndTime: 2019-03-09 13:28:10
% DurationCPUTime: 7.93s
% Computational Cost: add. (35892->639), mult. (73167->886), div. (0->0), fcn. (95636->14), ass. (0->312)
t316 = sin(pkin(6));
t388 = sin(pkin(12));
t389 = cos(pkin(12));
t397 = sin(qJ(2));
t398 = cos(qJ(2));
t326 = t388 * t398 + t389 * t397;
t288 = t326 * t316;
t381 = qJ(4) + qJ(5);
t314 = sin(t381);
t317 = cos(pkin(6));
t348 = cos(t381);
t275 = t288 * t314 - t317 * t348;
t276 = t288 * t348 + t317 * t314;
t299 = -t397 * t388 + t398 * t389;
t287 = t299 * t316;
t199 = Icges(6,5) * t276 - Icges(6,6) * t275 - Icges(6,3) * t287;
t200 = Icges(6,4) * t276 - Icges(6,2) * t275 - Icges(6,6) * t287;
t201 = Icges(6,1) * t276 - Icges(6,4) * t275 - Icges(6,5) * t287;
t289 = t326 * t317;
t320 = sin(qJ(1));
t323 = cos(qJ(1));
t272 = t289 * t323 + t299 * t320;
t339 = t316 * t348;
t232 = t272 * t314 + t323 * t339;
t385 = t316 * t323;
t233 = t272 * t348 - t314 * t385;
t325 = t317 * t299;
t271 = -t320 * t326 + t323 * t325;
t101 = -t199 * t271 - t200 * t232 + t201 * t233;
t273 = -t320 * t325 - t323 * t326;
t318 = sin(qJ(6));
t321 = cos(qJ(6));
t191 = -t233 * t318 - t271 * t321;
t192 = t233 * t321 - t271 * t318;
t124 = Icges(7,5) * t192 + Icges(7,6) * t191 + Icges(7,3) * t232;
t126 = Icges(7,4) * t192 + Icges(7,2) * t191 + Icges(7,6) * t232;
t128 = Icges(7,1) * t192 + Icges(7,4) * t191 + Icges(7,5) * t232;
t52 = t124 * t232 + t126 * t191 + t128 * t192;
t274 = -t289 * t320 + t299 * t323;
t386 = t316 * t320;
t235 = t274 * t348 + t314 * t386;
t193 = -t235 * t318 - t273 * t321;
t194 = t235 * t321 - t273 * t318;
t234 = t274 * t314 - t320 * t339;
t125 = Icges(7,5) * t194 + Icges(7,6) * t193 + Icges(7,3) * t234;
t127 = Icges(7,4) * t194 + Icges(7,2) * t193 + Icges(7,6) * t234;
t129 = Icges(7,1) * t194 + Icges(7,4) * t193 + Icges(7,5) * t234;
t53 = t125 * t232 + t127 * t191 + t129 * t192;
t226 = -t276 * t318 - t287 * t321;
t227 = t276 * t321 - t287 * t318;
t152 = Icges(7,5) * t227 + Icges(7,6) * t226 + Icges(7,3) * t275;
t153 = Icges(7,4) * t227 + Icges(7,2) * t226 + Icges(7,6) * t275;
t154 = Icges(7,1) * t227 + Icges(7,4) * t226 + Icges(7,5) * t275;
t67 = t152 * t232 + t153 * t191 + t154 * t192;
t11 = -t271 * t52 - t273 * t53 - t287 * t67;
t155 = Icges(6,5) * t233 - Icges(6,6) * t232 - Icges(6,3) * t271;
t157 = Icges(6,4) * t233 - Icges(6,2) * t232 - Icges(6,6) * t271;
t159 = Icges(6,1) * t233 - Icges(6,4) * t232 - Icges(6,5) * t271;
t78 = -t155 * t271 - t157 * t232 + t159 * t233;
t156 = Icges(6,5) * t235 - Icges(6,6) * t234 - Icges(6,3) * t273;
t158 = Icges(6,4) * t235 - Icges(6,2) * t234 - Icges(6,6) * t273;
t160 = Icges(6,1) * t235 - Icges(6,4) * t234 - Icges(6,5) * t273;
t79 = -t156 * t271 - t158 * t232 + t160 * t233;
t419 = -t101 * t287 - t271 * t78 - t273 * t79 + t11;
t102 = -t199 * t273 - t200 * t234 + t201 * t235;
t54 = t124 * t234 + t126 * t193 + t128 * t194;
t55 = t125 * t234 + t127 * t193 + t129 * t194;
t68 = t152 * t234 + t153 * t193 + t154 * t194;
t12 = -t271 * t54 - t273 * t55 - t287 * t68;
t80 = -t155 * t273 - t157 * t234 + t159 * t235;
t81 = -t156 * t273 - t158 * t234 + t160 * t235;
t418 = -t102 * t287 - t271 * t80 - t273 * t81 + t12;
t15 = t317 * t67 + (t320 * t53 - t323 * t52) * t316;
t417 = t15 + t101 * t317 + (t320 * t79 - t323 * t78) * t316;
t16 = t317 * t68 + (t320 * t55 - t323 * t54) * t316;
t416 = t16 + t102 * t317 + (t320 * t81 - t323 * t80) * t316;
t109 = -t287 * t199 - t275 * t200 + t276 * t201;
t106 = t109 * t287;
t60 = t125 * t275 + t127 * t226 + t129 * t227;
t392 = t60 * t273;
t59 = t124 * t275 + t126 * t226 + t128 * t227;
t393 = t59 * t271;
t77 = t275 * t152 + t226 * t153 + t227 * t154;
t73 = t77 * t287;
t22 = -t392 - t73 - t393;
t90 = -t156 * t287 - t158 * t275 + t160 * t276;
t390 = t90 * t273;
t89 = -t155 * t287 - t157 * t275 + t159 * t276;
t391 = t89 * t271;
t415 = t22 - t106 - t390 - t391;
t108 = t109 * t317;
t76 = t77 * t317;
t24 = t76 + (t60 * t320 - t59 * t323) * t316;
t414 = t24 + t108 + (t90 * t320 - t89 * t323) * t316;
t335 = -t192 * rSges(7,1) - t191 * rSges(7,2);
t130 = t232 * rSges(7,3) - t335;
t395 = t233 * pkin(5);
t378 = t232 * pkin(11) + t130 + t395;
t131 = t194 * rSges(7,1) + t193 * rSges(7,2) + t234 * rSges(7,3);
t377 = t235 * pkin(5) + pkin(11) * t234 + t131;
t161 = rSges(7,1) * t227 + rSges(7,2) * t226 + rSges(7,3) * t275;
t374 = pkin(5) * t276 + pkin(11) * t275 + t161;
t247 = Icges(4,5) * t288 + Icges(4,6) * t287 + Icges(4,3) * t317;
t248 = Icges(4,4) * t288 + Icges(4,2) * t287 + Icges(4,6) * t317;
t249 = Icges(4,1) * t288 + Icges(4,4) * t287 + Icges(4,5) * t317;
t283 = Icges(3,3) * t317 + (Icges(3,5) * t397 + Icges(3,6) * t398) * t316;
t284 = Icges(3,6) * t317 + (Icges(3,4) * t397 + Icges(3,2) * t398) * t316;
t285 = Icges(3,5) * t317 + (Icges(3,1) * t397 + Icges(3,4) * t398) * t316;
t355 = t316 * t397;
t413 = t316 * t398 * t284 + t287 * t248 + t288 * t249 + t285 * t355 + (t247 + t283) * t317;
t412 = t316 ^ 2;
t411 = m(4) / 0.2e1;
t410 = m(5) / 0.2e1;
t409 = m(6) / 0.2e1;
t408 = m(7) / 0.2e1;
t407 = t232 / 0.2e1;
t406 = t234 / 0.2e1;
t405 = -t271 / 0.2e1;
t404 = -t273 / 0.2e1;
t403 = t275 / 0.2e1;
t402 = -t287 / 0.2e1;
t401 = t317 / 0.2e1;
t400 = t320 / 0.2e1;
t399 = -t323 / 0.2e1;
t396 = pkin(1) * t323;
t322 = cos(qJ(4));
t312 = pkin(4) * t322 + pkin(3);
t394 = -pkin(3) + t312;
t324 = -pkin(10) - pkin(9);
t387 = t271 * t324;
t319 = sin(qJ(4));
t384 = t317 * t319;
t313 = pkin(2) * t398 + pkin(1);
t383 = t320 * t313;
t290 = t317 * t397 * pkin(2) + (-pkin(8) - qJ(3)) * t316;
t382 = t323 * t290;
t380 = t378 * t273;
t379 = t377 * t287;
t376 = t374 * t271;
t375 = t413 * t317;
t204 = Icges(4,5) * t274 + Icges(4,6) * t273 + Icges(4,3) * t386;
t351 = t323 * t397;
t354 = t320 * t398;
t296 = -t317 * t354 - t351;
t352 = t323 * t398;
t353 = t320 * t397;
t297 = -t317 * t353 + t352;
t254 = Icges(3,5) * t297 + Icges(3,6) * t296 + Icges(3,3) * t386;
t373 = t204 + t254;
t307 = t323 * t313;
t264 = -t396 + t307 + (-t316 * pkin(8) - t290) * t320;
t252 = t317 * t264;
t347 = -t274 * pkin(3) + pkin(9) * t273;
t372 = -t317 * t347 + t252;
t269 = t271 * pkin(9);
t218 = t272 * pkin(3) - t269;
t311 = pkin(8) * t385;
t263 = t382 + t311 + (-pkin(1) + t313) * t320;
t371 = -t218 - t263;
t370 = t263 * t386 + t264 * t385;
t203 = Icges(4,5) * t272 + Icges(4,6) * t271 - Icges(4,3) * t385;
t294 = t317 * t352 - t353;
t295 = t317 * t351 + t354;
t253 = Icges(3,5) * t295 + Icges(3,6) * t294 - Icges(3,3) * t385;
t369 = -t253 - t203;
t300 = pkin(2) * t355 + t317 * qJ(3);
t368 = -t288 * pkin(3) + t287 * pkin(9) - t300;
t367 = t59 / 0.2e1 + t67 / 0.2e1;
t366 = t60 / 0.2e1 + t68 / 0.2e1;
t365 = t319 * t386;
t364 = t319 * t385;
t277 = -t288 * t319 + t317 * t322;
t278 = t288 * t322 + t384;
t214 = Icges(5,5) * t278 + Icges(5,6) * t277 - Icges(5,3) * t287;
t215 = Icges(5,4) * t278 + Icges(5,2) * t277 - Icges(5,6) * t287;
t216 = Icges(5,1) * t278 + Icges(5,4) * t277 - Icges(5,5) * t287;
t242 = -t274 * t319 + t322 * t386;
t243 = t274 * t322 + t365;
t104 = -t214 * t273 + t215 * t242 + t216 * t243;
t166 = Icges(5,5) * t243 + Icges(5,6) * t242 - Icges(5,3) * t273;
t168 = Icges(5,4) * t243 + Icges(5,2) * t242 - Icges(5,6) * t273;
t170 = Icges(5,1) * t243 + Icges(5,4) * t242 - Icges(5,5) * t273;
t92 = -t166 * t287 + t168 * t277 + t170 * t278;
t363 = -t104 / 0.2e1 - t92 / 0.2e1;
t240 = -t272 * t319 - t322 * t385;
t241 = t272 * t322 - t364;
t103 = -t214 * t271 + t215 * t240 + t216 * t241;
t165 = Icges(5,5) * t241 + Icges(5,6) * t240 - Icges(5,3) * t271;
t167 = Icges(5,4) * t241 + Icges(5,2) * t240 - Icges(5,6) * t271;
t169 = Icges(5,1) * t241 + Icges(5,4) * t240 - Icges(5,5) * t271;
t91 = -t165 * t287 + t167 * t277 + t169 * t278;
t362 = -t91 / 0.2e1 - t103 / 0.2e1;
t357 = pkin(4) * t365 + t273 * t324 + t274 * t312;
t150 = t347 + t357;
t361 = t317 * t150 + t372;
t303 = pkin(4) * t364;
t149 = t272 * t394 + t269 - t303 + t387;
t360 = -t149 + t371;
t112 = -t287 * t214 + t277 * t215 + t278 * t216;
t198 = pkin(4) * t384 + t394 * t288 - (-pkin(9) - t324) * t287;
t359 = -t198 + t368;
t163 = t235 * rSges(6,1) - t234 * rSges(6,2) - t273 * rSges(6,3);
t172 = t243 * rSges(5,1) + t242 * rSges(5,2) - t273 * rSges(5,3);
t210 = t274 * rSges(4,1) + t273 * rSges(4,2) + rSges(4,3) * t386;
t261 = t297 * rSges(3,1) + t296 * rSges(3,2) + rSges(3,3) * t386;
t346 = t316 * (-rSges(4,1) * t288 - rSges(4,2) * t287 - rSges(4,3) * t317 - t300);
t345 = -t290 * t320 + t307;
t344 = t218 * t386 - t347 * t385 + t370;
t217 = rSges(5,1) * t278 + rSges(5,2) * t277 - rSges(5,3) * t287;
t343 = t316 * (-t217 + t368);
t71 = t77 * t275;
t18 = t59 * t232 + t60 * t234 + t71;
t5 = t232 * t52 + t234 * t53 + t275 * t67;
t6 = t232 * t54 + t234 * t55 + t275 * t68;
t340 = t11 * t407 + t12 * t406 + t18 * t402 + t22 * t403 + t6 * t404 + t5 * t405;
t338 = -t271 * t419 - t418 * t273 - t415 * t287;
t337 = -rSges(4,1) * t272 - rSges(4,2) * t271;
t336 = -t233 * rSges(6,1) + t232 * rSges(6,2);
t334 = -t382 - t383;
t202 = rSges(6,1) * t276 - rSges(6,2) * t275 - rSges(6,3) * t287;
t333 = t316 * (-t202 + t359);
t332 = t149 * t386 + t150 * t385 + t344;
t331 = t316 * (t359 - t374);
t330 = t345 + t357;
t171 = rSges(5,1) * t241 + rSges(5,2) * t240 - rSges(5,3) * t271;
t260 = rSges(3,1) * t295 + rSges(3,2) * t294 - rSges(3,3) * t385;
t329 = -t106 - t393 / 0.2e1 - t392 / 0.2e1 - t73 - t391 / 0.2e1 - t390 / 0.2e1 + (t67 + t101) * t405 + (t68 + t102) * t404;
t328 = -t272 * t312 + t303 + t334;
t327 = t417 * t405 + t416 * t404 + t414 * t402 + t415 * t401 + t418 * t386 / 0.2e1 - t419 * t385 / 0.2e1;
t305 = rSges(2,1) * t323 - rSges(2,2) * t320;
t304 = -rSges(2,1) * t320 - rSges(2,2) * t323;
t286 = t317 * rSges(3,3) + (rSges(3,1) * t397 + rSges(3,2) * t398) * t316;
t258 = Icges(3,1) * t297 + Icges(3,4) * t296 + Icges(3,5) * t386;
t257 = Icges(3,1) * t295 + Icges(3,4) * t294 - Icges(3,5) * t385;
t256 = Icges(3,4) * t297 + Icges(3,2) * t296 + Icges(3,6) * t386;
t255 = Icges(3,4) * t295 + Icges(3,2) * t294 - Icges(3,6) * t385;
t239 = pkin(8) * t386 + t261 + t396;
t238 = -pkin(1) * t320 - t260 + t311;
t222 = -t260 * t317 - t286 * t385;
t221 = t261 * t317 - t286 * t386;
t209 = -rSges(4,3) * t385 - t337;
t208 = Icges(4,1) * t274 + Icges(4,4) * t273 + Icges(4,5) * t386;
t207 = Icges(4,1) * t272 + Icges(4,4) * t271 - Icges(4,5) * t385;
t206 = Icges(4,4) * t274 + Icges(4,2) * t273 + Icges(4,6) * t386;
t205 = Icges(4,4) * t272 + Icges(4,2) * t271 - Icges(4,6) * t385;
t197 = (t260 * t320 + t261 * t323) * t316;
t196 = t283 * t386 + t284 * t296 + t285 * t297;
t195 = -t283 * t385 + t284 * t294 + t285 * t295;
t181 = t345 + t210;
t180 = -t383 + (rSges(4,3) * t316 - t290) * t323 + t337;
t177 = t271 * t202;
t175 = t271 * t198;
t174 = t317 * t254 + (t256 * t398 + t258 * t397) * t316;
t173 = t317 * t253 + (t255 * t398 + t257 * t397) * t316;
t162 = -t271 * rSges(6,3) - t336;
t151 = t287 * t163;
t143 = t273 * t162;
t141 = t287 * t150;
t140 = (-t209 - t263) * t317 + t323 * t346;
t139 = t210 * t317 + t320 * t346 + t252;
t136 = t273 * t149;
t135 = t247 * t386 + t248 * t273 + t249 * t274;
t134 = -t247 * t385 + t248 * t271 + t249 * t272;
t133 = t345 - t347 + t172;
t132 = -t171 - t218 + t334;
t123 = t330 + t163;
t122 = (rSges(6,3) - t324) * t271 + t328 + t336;
t120 = (t209 * t320 + t210 * t323) * t316 + t370;
t119 = -t172 * t287 + t217 * t273;
t118 = t171 * t287 - t217 * t271;
t116 = t202 * t273 - t151;
t115 = t162 * t287 - t177;
t114 = t204 * t317 + t206 * t287 + t208 * t288;
t113 = t203 * t317 + t205 * t287 + t207 * t288;
t111 = t112 * t317;
t110 = t112 * t287;
t107 = -t171 * t273 + t172 * t271;
t105 = t163 * t271 - t143;
t100 = (-t171 + t371) * t317 + t323 * t343;
t99 = t172 * t317 + t320 * t343 + t372;
t96 = t330 + t377;
t95 = -t395 - t387 + (-rSges(7,3) - pkin(11)) * t232 + t328 + t335;
t94 = t131 * t275 - t161 * t234;
t93 = -t130 * t275 + t161 * t232;
t88 = (t171 * t320 + t172 * t323) * t316 + t344;
t87 = -t166 * t273 + t168 * t242 + t170 * t243;
t86 = -t165 * t273 + t167 * t242 + t169 * t243;
t85 = -t166 * t271 + t168 * t240 + t170 * t241;
t84 = -t165 * t271 + t167 * t240 + t169 * t241;
t75 = -t141 - t151 + (t198 + t202) * t273;
t74 = -t175 - t177 - (-t149 - t162) * t287;
t72 = t130 * t234 - t131 * t232;
t70 = t273 * t374 - t379;
t69 = t287 * t378 - t376;
t66 = (-t162 + t360) * t317 + t323 * t333;
t65 = t163 * t317 + t320 * t333 + t361;
t62 = -t136 - t143 + (t150 + t163) * t271;
t61 = t271 * t377 - t380;
t58 = (t162 * t320 + t163 * t323) * t316 + t332;
t51 = -t141 + (t198 + t374) * t273 - t379;
t50 = -t175 - (-t149 - t378) * t287 - t376;
t49 = (t360 - t378) * t317 + t323 * t331;
t48 = t317 * t377 + t320 * t331 + t361;
t47 = -t136 + (t150 + t377) * t271 - t380;
t46 = t111 + (t92 * t320 - t91 * t323) * t316;
t45 = (t320 * t378 + t323 * t377) * t316 + t332;
t44 = -t91 * t271 - t92 * t273 - t110;
t38 = t104 * t317 + (t320 * t87 - t323 * t86) * t316;
t37 = t103 * t317 + (t320 * t85 - t323 * t84) * t316;
t36 = -t104 * t287 - t271 * t86 - t273 * t87;
t35 = -t103 * t287 - t271 * t84 - t273 * t85;
t1 = [m(7) * (t95 ^ 2 + t96 ^ 2) + m(6) * (t122 ^ 2 + t123 ^ 2) + m(5) * (t132 ^ 2 + t133 ^ 2) + m(4) * (t180 ^ 2 + t181 ^ 2) + m(3) * (t238 ^ 2 + t239 ^ 2) + m(2) * (t304 ^ 2 + t305 ^ 2) + Icges(2,3) + t112 + t109 + t77 + t413; t76 + t111 + t108 + m(7) * (t48 * t96 + t49 * t95) + m(6) * (t122 * t66 + t123 * t65) + m(5) * (t100 * t132 + t133 * t99) + m(4) * (t139 * t181 + t140 * t180) + m(3) * (t221 * t239 + t222 * t238) + ((-t173 / 0.2e1 - t113 / 0.2e1 - t101 / 0.2e1 - t134 / 0.2e1 - t195 / 0.2e1 - t89 / 0.2e1 + t362 - t367) * t323 + (t90 / 0.2e1 + t174 / 0.2e1 + t114 / 0.2e1 + t102 / 0.2e1 + t135 / 0.2e1 + t196 / 0.2e1 - t363 + t366) * t320) * t316 + t375; (t46 + t375 + t414) * t317 + m(7) * (t45 ^ 2 + t48 ^ 2 + t49 ^ 2) + m(6) * (t58 ^ 2 + t65 ^ 2 + t66 ^ 2) + m(5) * (t100 ^ 2 + t88 ^ 2 + t99 ^ 2) + m(4) * (t120 ^ 2 + t139 ^ 2 + t140 ^ 2) + m(3) * (t197 ^ 2 + t221 ^ 2 + t222 ^ 2) + ((-t37 + ((t205 * t271 + t207 * t272 + t255 * t294 + t257 * t295) * t316 + t369 * t412 * t323) * t323 + (-t113 - t134 - t173 - t195) * t317 - t417) * t323 + (t38 + ((t206 * t273 + t208 * t274 + t256 * t296 + t258 * t297) * t316 + t373 * t412 * t320) * t320 + (t135 + t196 + t174 + t114) * t317 + (-t205 * t273 - t206 * t271 - t207 * t274 - t208 * t272 - t255 * t296 - t256 * t294 - t257 * t297 - t258 * t295 + (t320 * t369 + t323 * t373) * t316) * t385 + t416) * t320) * t316; 0.2e1 * ((t320 * t95 - t323 * t96) * t408 + (t122 * t320 - t123 * t323) * t409 + (t132 * t320 - t133 * t323) * t410 + (t180 * t320 - t181 * t323) * t411) * t316; m(7) * (t317 * t45 + (t320 * t49 - t323 * t48) * t316) + m(6) * (t317 * t58 + (t320 * t66 - t323 * t65) * t316) + m(5) * (t317 * t88 + (t100 * t320 - t323 * t99) * t316) + m(4) * (t120 * t317 + (-t139 * t323 + t140 * t320) * t316); 0.2e1 * (t411 + t410 + t409 + t408) * (t317 ^ 2 + (t320 ^ 2 + t323 ^ 2) * t412); t362 * t271 + m(7) * (t50 * t95 + t51 * t96) + m(6) * (t122 * t74 + t123 * t75) + m(5) * (t118 * t132 + t119 * t133) + t363 * t273 - t110 + t329; t327 + (t35 * t399 + t36 * t400) * t316 + m(7) * (t45 * t47 + t48 * t51 + t49 * t50) + m(6) * (t58 * t62 + t65 * t75 + t66 * t74) + m(5) * (t100 * t118 + t107 * t88 + t119 * t99) + t37 * t405 + t38 * t404 + t46 * t402 + t44 * t401; m(5) * (t107 * t317 + (t118 * t320 - t119 * t323) * t316) + m(6) * (t317 * t62 + (t320 * t74 - t323 * t75) * t316) + m(7) * (t317 * t47 + (t320 * t50 - t323 * t51) * t316); -t271 * t35 - t273 * t36 - t287 * t44 + m(7) * (t47 ^ 2 + t50 ^ 2 + t51 ^ 2) + m(6) * (t62 ^ 2 + t74 ^ 2 + t75 ^ 2) + m(5) * (t107 ^ 2 + t118 ^ 2 + t119 ^ 2) + t338; m(7) * (t69 * t95 + t70 * t96) + m(6) * (t115 * t122 + t116 * t123) + t329; t327 + m(7) * (t45 * t61 + t48 * t70 + t49 * t69) + m(6) * (t105 * t58 + t115 * t66 + t116 * t65); m(6) * (t105 * t317 + (t115 * t320 - t116 * t323) * t316) + m(7) * (t317 * t61 + (t320 * t69 - t323 * t70) * t316); m(7) * (t47 * t61 + t50 * t69 + t51 * t70) + m(6) * (t105 * t62 + t115 * t74 + t116 * t75) + t338; m(7) * (t61 ^ 2 + t69 ^ 2 + t70 ^ 2) + m(6) * (t105 ^ 2 + t115 ^ 2 + t116 ^ 2) + t338; m(7) * (t93 * t95 + t94 * t96) + t71 + t366 * t234 + t367 * t232; t18 * t401 + t16 * t406 + t24 * t403 + m(7) * (t45 * t72 + t48 * t94 + t49 * t93) + t15 * t407 + (t399 * t5 + t400 * t6) * t316; m(7) * (t317 * t72 + (t320 * t93 - t323 * t94) * t316); m(7) * (t47 * t72 + t50 * t93 + t51 * t94) + t340; m(7) * (t61 * t72 + t69 * t93 + t70 * t94) + t340; t234 * t6 + t232 * t5 + t275 * t18 + m(7) * (t72 ^ 2 + t93 ^ 2 + t94 ^ 2);];
%% Postprocessing: Reshape Output
% From vec2symmat_6_matlab.m
res = [t1(1) t1(2) t1(4) t1(7) t1(11) t1(16); t1(2) t1(3) t1(5) t1(8) t1(12) t1(17); t1(4) t1(5) t1(6) t1(9) t1(13) t1(18); t1(7) t1(8) t1(9) t1(10) t1(14) t1(19); t1(11) t1(12) t1(13) t1(14) t1(15) t1(20); t1(16) t1(17) t1(18) t1(19) t1(20) t1(21);];
Mq  = res;
