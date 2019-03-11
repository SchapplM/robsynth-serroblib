% Calculate joint inertia matrix for
% S6RRRPRR7
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% pkin [12x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,alpha2,d1,d2,d3,d5,d6,theta4]';
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
% Datum: 2019-03-09 18:43
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Mq = S6RRRPRR7_inertiaJ_slag_vp1(qJ, ...
  pkin, m, rSges, Icges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(12,1),zeros(7,1),zeros(7,3),zeros(7,6)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRRPRR7_inertiaJ_slag_vp1: qJ has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [12 1]), ...
  'S6RRRPRR7_inertiaJ_slag_vp1: pkin has to be [12x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6RRRPRR7_inertiaJ_slag_vp1: m has to be [7x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [7,3]), ...
  'S6RRRPRR7_inertiaJ_slag_vp1: rSges has to be [7x3] (double)');
assert(isreal(Icges) && all(size(Icges) == [7 6]), ...
  'S6RRRPRR7_inertiaJ_slag_vp1: Icges has to be [7x6] (double)'); 

%% Symbolic Calculation
% From inertia_joint_joint_floatb_twist_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 18:32:56
% EndTime: 2019-03-09 18:33:13
% DurationCPUTime: 7.31s
% Computational Cost: add. (34850->679), mult. (46780->926), div. (0->0), fcn. (58142->14), ass. (0->327)
t331 = cos(pkin(6));
t339 = cos(qJ(2));
t340 = cos(qJ(1));
t400 = t339 * t340;
t335 = sin(qJ(2));
t336 = sin(qJ(1));
t403 = t335 * t336;
t305 = -t331 * t400 + t403;
t401 = t336 * t339;
t402 = t335 * t340;
t307 = t331 * t401 + t402;
t330 = sin(pkin(6));
t406 = t330 * t339;
t306 = t331 * t402 + t401;
t329 = qJ(3) + pkin(12);
t364 = qJ(5) + t329;
t322 = sin(t364);
t357 = cos(t364);
t405 = t330 * t340;
t262 = t306 * t357 - t322 * t405;
t333 = sin(qJ(6));
t337 = cos(qJ(6));
t223 = -t262 * t333 + t305 * t337;
t224 = t262 * t337 + t305 * t333;
t349 = t330 * t357;
t261 = t306 * t322 + t340 * t349;
t136 = Icges(7,5) * t224 + Icges(7,6) * t223 + Icges(7,3) * t261;
t138 = Icges(7,4) * t224 + Icges(7,2) * t223 + Icges(7,6) * t261;
t140 = Icges(7,1) * t224 + Icges(7,4) * t223 + Icges(7,5) * t261;
t58 = t136 * t261 + t138 * t223 + t140 * t224;
t308 = -t331 * t403 + t400;
t408 = t330 * t336;
t264 = t308 * t357 + t322 * t408;
t225 = -t264 * t333 + t307 * t337;
t226 = t264 * t337 + t307 * t333;
t263 = t308 * t322 - t336 * t349;
t137 = Icges(7,5) * t226 + Icges(7,6) * t225 + Icges(7,3) * t263;
t139 = Icges(7,4) * t226 + Icges(7,2) * t225 + Icges(7,6) * t263;
t141 = Icges(7,1) * t226 + Icges(7,4) * t225 + Icges(7,5) * t263;
t59 = t137 * t261 + t139 * t223 + t141 * t224;
t284 = t331 * t322 + t335 * t349;
t259 = -t284 * t333 - t337 * t406;
t260 = t284 * t337 - t333 * t406;
t409 = t330 * t335;
t283 = t322 * t409 - t331 * t357;
t167 = Icges(7,5) * t260 + Icges(7,6) * t259 + Icges(7,3) * t283;
t168 = Icges(7,4) * t260 + Icges(7,2) * t259 + Icges(7,6) * t283;
t169 = Icges(7,1) * t260 + Icges(7,4) * t259 + Icges(7,5) * t283;
t73 = t167 * t261 + t168 * t223 + t169 * t224;
t11 = t305 * t58 + t307 * t59 - t406 * t73;
t227 = Icges(6,5) * t284 - Icges(6,6) * t283 - Icges(6,3) * t406;
t228 = Icges(6,4) * t284 - Icges(6,2) * t283 - Icges(6,6) * t406;
t229 = Icges(6,1) * t284 - Icges(6,4) * t283 - Icges(6,5) * t406;
t113 = t227 * t305 - t228 * t261 + t229 * t262;
t173 = Icges(6,5) * t262 - Icges(6,6) * t261 + Icges(6,3) * t305;
t175 = Icges(6,4) * t262 - Icges(6,2) * t261 + Icges(6,6) * t305;
t177 = Icges(6,1) * t262 - Icges(6,4) * t261 + Icges(6,5) * t305;
t83 = t173 * t305 - t175 * t261 + t177 * t262;
t174 = Icges(6,5) * t264 - Icges(6,6) * t263 + Icges(6,3) * t307;
t176 = Icges(6,4) * t264 - Icges(6,2) * t263 + Icges(6,6) * t307;
t178 = Icges(6,1) * t264 - Icges(6,4) * t263 + Icges(6,5) * t307;
t84 = t174 * t305 - t176 * t261 + t178 * t262;
t433 = -t113 * t406 + t305 * t83 + t307 * t84 + t11;
t114 = t227 * t307 - t228 * t263 + t229 * t264;
t60 = t136 * t263 + t138 * t225 + t140 * t226;
t61 = t137 * t263 + t139 * t225 + t141 * t226;
t74 = t167 * t263 + t168 * t225 + t169 * t226;
t12 = t305 * t60 + t307 * t61 - t406 * t74;
t85 = t173 * t307 - t175 * t263 + t177 * t264;
t86 = t174 * t307 - t176 * t263 + t178 * t264;
t432 = -t114 * t406 + t305 * t85 + t307 * t86 + t12;
t15 = t73 * t331 + (t336 * t59 - t340 * t58) * t330;
t431 = t15 + t113 * t331 + (t336 * t84 - t340 * t83) * t330;
t16 = t74 * t331 + (t336 * t61 - t340 * t60) * t330;
t430 = t16 + t114 * t331 + (t336 * t86 - t340 * t85) * t330;
t389 = -t283 * t228 + t284 * t229;
t119 = -t227 * t406 + t389;
t65 = t137 * t283 + t139 * t259 + t141 * t260;
t414 = t65 * t307;
t64 = t136 * t283 + t138 * t259 + t140 * t260;
t415 = t64 * t305;
t81 = t283 * t167 + t259 * t168 + t260 * t169;
t21 = -t406 * t81 + t414 + t415;
t96 = -t174 * t406 - t176 * t283 + t178 * t284;
t412 = t96 * t307;
t95 = -t173 * t406 - t175 * t283 + t177 * t284;
t413 = t95 * t305;
t429 = -t119 * t406 + t21 + t412 + t413;
t118 = t119 * t331;
t80 = t81 * t331;
t23 = t80 + (t65 * t336 - t64 * t340) * t330;
t428 = t23 + t118 + (t96 * t336 - t95 * t340) * t330;
t352 = -t224 * rSges(7,1) - t223 * rSges(7,2);
t142 = rSges(7,3) * t261 - t352;
t418 = t262 * pkin(5);
t398 = pkin(11) * t261 + t142 + t418;
t170 = rSges(7,1) * t260 + rSges(7,2) * t259 + rSges(7,3) * t283;
t427 = pkin(5) * t284 + pkin(11) * t283 + t170;
t426 = t330 ^ 2;
t425 = t261 / 0.2e1;
t424 = t263 / 0.2e1;
t423 = t283 / 0.2e1;
t422 = t305 / 0.2e1;
t421 = t307 / 0.2e1;
t420 = t331 / 0.2e1;
t334 = sin(qJ(3));
t419 = pkin(3) * t334;
t338 = cos(qJ(3));
t323 = t338 * pkin(3) + pkin(2);
t417 = -pkin(2) + t323;
t332 = -qJ(4) - pkin(9);
t411 = -t119 - t81;
t243 = Icges(3,5) * t306 - Icges(3,6) * t305 - Icges(3,3) * t405;
t410 = t243 * t340;
t407 = t330 * t338;
t404 = t331 * t334;
t399 = t398 * t307;
t143 = t226 * rSges(7,1) + t225 * rSges(7,2) + t263 * rSges(7,3);
t397 = t264 * pkin(5) + pkin(11) * t263 + t143;
t324 = sin(t329);
t312 = pkin(4) * t324 + t419;
t295 = t312 * t405;
t376 = t334 * t405;
t314 = pkin(3) * t376;
t328 = -pkin(10) + t332;
t380 = t328 - t332;
t325 = cos(t329);
t311 = pkin(4) * t325 + t323;
t382 = t311 - t323;
t160 = -t305 * t380 + t306 * t382 - t295 + t314;
t301 = t305 * pkin(9);
t195 = -t305 * t332 + t306 * t417 - t301 - t314;
t164 = t307 * t195;
t396 = t307 * t160 + t164;
t377 = t334 * t408;
t369 = -pkin(3) * t377 + t307 * t332 - t308 * t323;
t371 = -t307 * t328 + t308 * t311 + t312 * t408;
t161 = t369 + t371;
t266 = t308 * pkin(2) + pkin(9) * t307;
t196 = -t266 - t369;
t395 = -t161 - t196;
t241 = pkin(3) * t404 + ((pkin(9) + t332) * t339 + t417 * t335) * t330;
t393 = t195 * t406 + t305 * t241;
t254 = t331 * t266;
t392 = t331 * t196 + t254;
t271 = -t308 * t324 + t325 * t408;
t272 = t308 * t325 + t324 * t408;
t194 = t272 * rSges(5,1) + t271 * rSges(5,2) + t307 * rSges(5,3);
t391 = -t194 - t196;
t265 = pkin(2) * t306 + t301;
t390 = -t195 - t265;
t275 = -t306 * t334 - t338 * t405;
t276 = t306 * t338 - t376;
t207 = rSges(4,1) * t276 + rSges(4,2) * t275 + rSges(4,3) * t305;
t388 = -t207 - t265;
t289 = -t324 * t409 + t325 * t331;
t290 = t324 * t331 + t325 * t409;
t234 = Icges(5,4) * t290 + Icges(5,2) * t289 - Icges(5,6) * t406;
t235 = Icges(5,1) * t290 + Icges(5,4) * t289 - Icges(5,5) * t406;
t387 = t289 * t234 + t290 * t235;
t353 = -t262 * rSges(6,1) + t261 * rSges(6,2);
t179 = rSges(6,3) * t305 - t353;
t230 = rSges(6,1) * t284 - rSges(6,2) * t283 - rSges(6,3) * t406;
t130 = t179 * t406 + t305 * t230;
t216 = (t312 - t419) * t331 + (t335 * t382 + t339 * t380) * t330;
t386 = -t216 - t241;
t303 = t331 * t338 - t334 * t409;
t304 = t335 * t407 + t404;
t239 = Icges(4,4) * t304 + Icges(4,2) * t303 - Icges(4,6) * t406;
t240 = Icges(4,1) * t304 + Icges(4,4) * t303 - Icges(4,5) * t406;
t385 = t303 * t239 + t304 * t240;
t236 = rSges(5,1) * t290 + rSges(5,2) * t289 - rSges(5,3) * t406;
t384 = -t236 - t241;
t383 = t265 * t408 + t266 * t405;
t381 = t340 * pkin(1) + pkin(8) * t408;
t379 = t64 / 0.2e1 + t73 / 0.2e1;
t378 = t65 / 0.2e1 + t74 / 0.2e1;
t375 = t331 * t161 + t392;
t374 = -t160 + t390;
t180 = t264 * rSges(6,1) - t263 * rSges(6,2) + t307 * rSges(6,3);
t373 = -t180 + t395;
t372 = -t230 + t386;
t277 = -t308 * t334 + t336 * t407;
t278 = t308 * t338 + t377;
t208 = t278 * rSges(4,1) + t277 * rSges(4,2) + t307 * rSges(4,3);
t286 = Icges(3,3) * t331 + (Icges(3,5) * t335 + Icges(3,6) * t339) * t330;
t287 = Icges(3,6) * t331 + (Icges(3,4) * t335 + Icges(3,2) * t339) * t330;
t288 = Icges(3,5) * t331 + (Icges(3,1) * t335 + Icges(3,4) * t339) * t330;
t370 = t331 * t286 + t287 * t406 + t288 * t409;
t250 = t308 * rSges(3,1) - t307 * rSges(3,2) + rSges(3,3) * t408;
t367 = -t406 / 0.2e1;
t365 = t305 * t433 + t432 * t307;
t363 = -t336 * pkin(1) + pkin(8) * t405;
t242 = rSges(4,1) * t304 + rSges(4,2) * t303 - rSges(4,3) * t406;
t309 = (pkin(2) * t335 - pkin(9) * t339) * t330;
t362 = t330 * (-t242 - t309);
t361 = t395 - t397;
t360 = t160 * t406 + t305 * t216 + t393;
t359 = t386 - t427;
t358 = t195 * t408 + t196 * t405 + t383;
t76 = t305 * t427 + t398 * t406;
t356 = t330 * (-t309 + t384);
t78 = t81 * t283;
t18 = t64 * t261 + t65 * t263 + t78;
t3 = t261 * t58 + t263 * t59 + t283 * t73;
t4 = t261 * t60 + t263 * t61 + t283 * t74;
t355 = t11 * t425 + t12 * t424 + t18 * t367 + t21 * t423 + t3 * t422 + t4 * t421;
t269 = -t306 * t324 - t325 * t405;
t270 = t306 * t325 - t324 * t405;
t354 = -t270 * rSges(5,1) - t269 * rSges(5,2);
t351 = t371 + t381;
t350 = t330 * (-t309 + t372);
t348 = t160 * t408 + t161 * t405 + t358;
t188 = Icges(5,5) * t272 + Icges(5,6) * t271 + Icges(5,3) * t307;
t190 = Icges(5,4) * t272 + Icges(5,2) * t271 + Icges(5,6) * t307;
t192 = Icges(5,1) * t272 + Icges(5,4) * t271 + Icges(5,5) * t307;
t102 = -t188 * t406 + t190 * t289 + t192 * t290;
t202 = Icges(4,5) * t278 + Icges(4,6) * t277 + Icges(4,3) * t307;
t204 = Icges(4,4) * t278 + Icges(4,2) * t277 + Icges(4,6) * t307;
t206 = Icges(4,1) * t278 + Icges(4,4) * t277 + Icges(4,5) * t307;
t112 = -t202 * t406 + t204 * t303 + t206 * t304;
t233 = Icges(5,5) * t290 + Icges(5,6) * t289 - Icges(5,3) * t406;
t117 = t233 * t307 + t234 * t271 + t235 * t272;
t238 = Icges(4,5) * t304 + Icges(4,6) * t303 - Icges(4,3) * t406;
t124 = t238 * t307 + t239 * t277 + t240 * t278;
t347 = t102 / 0.2e1 + t117 / 0.2e1 + t124 / 0.2e1 + t112 / 0.2e1;
t187 = Icges(5,5) * t270 + Icges(5,6) * t269 + Icges(5,3) * t305;
t189 = Icges(5,4) * t270 + Icges(5,2) * t269 + Icges(5,6) * t305;
t191 = Icges(5,1) * t270 + Icges(5,4) * t269 + Icges(5,5) * t305;
t101 = -t187 * t406 + t189 * t289 + t191 * t290;
t201 = Icges(4,5) * t276 + Icges(4,6) * t275 + Icges(4,3) * t305;
t203 = Icges(4,4) * t276 + Icges(4,2) * t275 + Icges(4,6) * t305;
t205 = Icges(4,1) * t276 + Icges(4,4) * t275 + Icges(4,5) * t305;
t111 = -t201 * t406 + t203 * t303 + t205 * t304;
t116 = t233 * t305 + t234 * t269 + t235 * t270;
t123 = t238 * t305 + t239 * t275 + t240 * t276;
t346 = t111 / 0.2e1 + t101 / 0.2e1 + t123 / 0.2e1 + t116 / 0.2e1;
t345 = t330 * (-t309 + t359);
t344 = -t306 * t311 + t295 + t363;
t343 = t415 / 0.2e1 + t414 / 0.2e1 + t413 / 0.2e1 + t412 / 0.2e1 + (t113 + t73) * t422 + (t114 + t74) * t421;
t342 = -t406 * t429 + t365;
t249 = t306 * rSges(3,1) - t305 * rSges(3,2) - rSges(3,3) * t405;
t341 = t431 * t422 + t430 * t421 + t429 * t420 + t432 * t408 / 0.2e1 + t428 * t367 - t433 * t405 / 0.2e1;
t317 = rSges(2,1) * t340 - t336 * rSges(2,2);
t316 = -t336 * rSges(2,1) - rSges(2,2) * t340;
t291 = t331 * rSges(3,3) + (rSges(3,1) * t335 + rSges(3,2) * t339) * t330;
t248 = Icges(3,1) * t308 - Icges(3,4) * t307 + Icges(3,5) * t408;
t247 = Icges(3,1) * t306 - Icges(3,4) * t305 - Icges(3,5) * t405;
t246 = Icges(3,4) * t308 - Icges(3,2) * t307 + Icges(3,6) * t408;
t245 = Icges(3,4) * t306 - Icges(3,2) * t305 - Icges(3,6) * t405;
t244 = Icges(3,5) * t308 - Icges(3,6) * t307 + Icges(3,3) * t408;
t232 = t250 + t381;
t231 = -t249 + t363;
t215 = -t331 * t249 - t291 * t405;
t214 = t250 * t331 - t291 * t408;
t197 = t370 * t331;
t193 = rSges(5,3) * t305 - t354;
t172 = (t249 * t336 + t250 * t340) * t330;
t166 = t286 * t408 - t287 * t307 + t288 * t308;
t165 = -t286 * t405 - t305 * t287 + t306 * t288;
t163 = t307 * t179;
t155 = t266 + t208 + t381;
t154 = t363 + t388;
t149 = -t208 * t406 - t242 * t307;
t148 = t207 * t406 + t242 * t305;
t147 = t331 * t244 + (t246 * t339 + t248 * t335) * t330;
t146 = t331 * t243 + (t245 * t339 + t247 * t335) * t330;
t145 = -t369 + t194 + t381;
t144 = -t306 * t323 + t314 + (-rSges(5,3) + t332) * t305 + t354 + t363;
t134 = t351 + t180;
t133 = (-rSges(6,3) + t328) * t305 + t344 + t353;
t132 = -t238 * t406 + t385;
t131 = -t180 * t406 - t307 * t230;
t128 = t132 * t331;
t127 = t207 * t307 - t208 * t305;
t126 = t331 * t388 + t340 * t362;
t125 = t331 * t208 + t336 * t362 + t254;
t122 = -t233 * t406 + t387;
t121 = t122 * t331;
t120 = -t180 * t305 + t163;
t115 = (t207 * t336 + t208 * t340) * t330 + t383;
t108 = t351 + t397;
t107 = -t418 + t305 * t328 + (-rSges(7,3) - pkin(11)) * t261 + t344 + t352;
t106 = t202 * t307 + t204 * t277 + t206 * t278;
t105 = t201 * t307 + t203 * t277 + t205 * t278;
t104 = t202 * t305 + t204 * t275 + t206 * t276;
t103 = t201 * t305 + t203 * t275 + t205 * t276;
t100 = t307 * t384 + t391 * t406;
t99 = t193 * t406 + t236 * t305 + t393;
t98 = t143 * t283 - t170 * t263;
t97 = -t142 * t283 + t170 * t261;
t94 = t188 * t307 + t190 * t271 + t192 * t272;
t93 = t187 * t307 + t189 * t271 + t191 * t272;
t92 = t188 * t305 + t190 * t269 + t192 * t270;
t91 = t187 * t305 + t189 * t269 + t191 * t270;
t90 = (-t193 + t390) * t331 + t340 * t356;
t89 = t331 * t194 + t336 * t356 + t392;
t82 = t142 * t263 - t143 * t261;
t79 = t307 * t193 + t305 * t391 + t164;
t77 = -t307 * t427 - t397 * t406;
t75 = (t193 * t336 + t194 * t340) * t330 + t358;
t70 = -t305 * t397 + t399;
t69 = t307 * t372 + t373 * t406;
t68 = t360 + t130;
t67 = (-t179 + t374) * t331 + t340 * t350;
t66 = t331 * t180 + t336 * t350 + t375;
t57 = t305 * t373 + t163 + t396;
t56 = (t179 * t336 + t180 * t340) * t330 + t348;
t55 = t128 + (-t111 * t340 + t112 * t336) * t330;
t54 = t111 * t305 + t112 * t307 - t132 * t406;
t53 = t307 * t359 + t361 * t406;
t52 = t76 + t360;
t51 = (t374 - t398) * t331 + t340 * t345;
t50 = t331 * t397 + t336 * t345 + t375;
t49 = t124 * t331 + (-t105 * t340 + t106 * t336) * t330;
t48 = t123 * t331 + (-t103 * t340 + t104 * t336) * t330;
t47 = t121 + (-t101 * t340 + t102 * t336) * t330;
t46 = t105 * t305 + t106 * t307 - t124 * t406;
t45 = t103 * t305 + t104 * t307 - t123 * t406;
t44 = t101 * t305 + t102 * t307 - t122 * t406;
t39 = t117 * t331 + (t336 * t94 - t340 * t93) * t330;
t38 = t116 * t331 + (t336 * t92 - t340 * t91) * t330;
t37 = -t117 * t406 + t305 * t93 + t307 * t94;
t36 = -t116 * t406 + t305 * t91 + t307 * t92;
t35 = t305 * t361 + t396 + t399;
t34 = (t336 * t398 + t340 * t397) * t330 + t348;
t1 = [Icges(2,3) + (-t227 - t233 - t238) * t406 + m(7) * (t107 ^ 2 + t108 ^ 2) + m(6) * (t133 ^ 2 + t134 ^ 2) + m(5) * (t144 ^ 2 + t145 ^ 2) + m(4) * (t154 ^ 2 + t155 ^ 2) + m(3) * (t231 ^ 2 + t232 ^ 2) + m(2) * (t316 ^ 2 + t317 ^ 2) + t370 + t81 + t385 + t387 + t389; t118 + t121 + t197 + t80 + t128 + m(7) * (t107 * t51 + t108 * t50) + m(6) * (t133 * t67 + t134 * t66) + m(5) * (t144 * t90 + t145 * t89) + m(4) * (t125 * t155 + t126 * t154) + m(3) * (t214 * t232 + t215 * t231) + ((-t113 / 0.2e1 - t165 / 0.2e1 - t95 / 0.2e1 - t146 / 0.2e1 - t346 - t379) * t340 + (t147 / 0.2e1 + t114 / 0.2e1 + t166 / 0.2e1 + t96 / 0.2e1 + t347 + t378) * t336) * t330; (t47 + t55 + t197 + t428) * t331 + m(7) * (t34 ^ 2 + t50 ^ 2 + t51 ^ 2) + m(6) * (t56 ^ 2 + t66 ^ 2 + t67 ^ 2) + m(5) * (t75 ^ 2 + t89 ^ 2 + t90 ^ 2) + m(4) * (t115 ^ 2 + t125 ^ 2 + t126 ^ 2) + m(3) * (t172 ^ 2 + t214 ^ 2 + t215 ^ 2) + ((-t48 - t38 + ((-t305 * t245 + t306 * t247) * t330 - t426 * t410) * t340 - t431) * t340 + (t49 + t39 + ((-t246 * t307 + t248 * t308 + (t244 * t336 - t410) * t330) * t336 + (t244 * t405 + t245 * t307 + t305 * t246 - t247 * t308 - t306 * t248) * t340) * t330 + t430) * t336 + ((-t146 - t165) * t340 + (t147 + t166) * t336) * t331) * t330; (-t122 - t132 + t411) * t406 + t347 * t307 + t346 * t305 + m(7) * (t107 * t52 + t108 * t53) + m(6) * (t133 * t68 + t134 * t69) + m(5) * (t100 * t145 + t144 * t99) + m(4) * (t148 * t154 + t149 * t155) + t343; ((-t36 / 0.2e1 - t45 / 0.2e1) * t340 + (-t47 / 0.2e1 - t55 / 0.2e1) * t339 + (t37 / 0.2e1 + t46 / 0.2e1) * t336) * t330 + m(7) * (t34 * t35 + t50 * t53 + t51 * t52) + m(6) * (t56 * t57 + t66 * t69 + t67 * t68) + m(5) * (t100 * t89 + t75 * t79 + t90 * t99) + m(4) * (t115 * t127 + t125 * t149 + t126 * t148) + (t44 / 0.2e1 + t54 / 0.2e1) * t331 + (t39 / 0.2e1 + t49 / 0.2e1) * t307 + (t38 / 0.2e1 + t48 / 0.2e1) * t305 + t341; (t37 + t46) * t307 + (t36 + t45) * t305 + (-t44 - t54 - t429) * t406 + m(7) * (t35 ^ 2 + t52 ^ 2 + t53 ^ 2) + m(6) * (t57 ^ 2 + t68 ^ 2 + t69 ^ 2) + m(5) * (t100 ^ 2 + t79 ^ 2 + t99 ^ 2) + m(4) * (t127 ^ 2 + t148 ^ 2 + t149 ^ 2) + t365; m(7) * (t107 * t307 + t108 * t305) + m(6) * (t133 * t307 + t134 * t305) + m(5) * (t144 * t307 + t145 * t305); m(7) * (t305 * t50 + t307 * t51 - t34 * t406) + m(6) * (t305 * t66 + t307 * t67 - t406 * t56) + m(5) * (t305 * t89 + t307 * t90 - t406 * t75); m(7) * (t305 * t53 + t307 * t52 - t35 * t406) + m(6) * (t305 * t69 + t307 * t68 - t406 * t57) + m(5) * (t100 * t305 + t307 * t99 - t406 * t79); 0.2e1 * (m(5) / 0.2e1 + m(6) / 0.2e1 + m(7) / 0.2e1) * (t339 ^ 2 * t426 + t305 ^ 2 + t307 ^ 2); t411 * t406 + m(7) * (t107 * t76 + t108 * t77) + m(6) * (t130 * t133 + t131 * t134) + t343; m(7) * (t34 * t70 + t50 * t77 + t51 * t76) + m(6) * (t120 * t56 + t130 * t67 + t131 * t66) + t341; m(7) * (t35 * t70 + t52 * t76 + t53 * t77) + m(6) * (t120 * t57 + t130 * t68 + t131 * t69) + t342; m(6) * (-t120 * t406 + t130 * t307 + t131 * t305) + m(7) * (t305 * t77 + t307 * t76 - t406 * t70); m(7) * (t70 ^ 2 + t76 ^ 2 + t77 ^ 2) + m(6) * (t120 ^ 2 + t130 ^ 2 + t131 ^ 2) + t342; t78 + m(7) * (t107 * t97 + t108 * t98) + t378 * t263 + t379 * t261; t15 * t425 + m(7) * (t34 * t82 + t50 * t98 + t51 * t97) + t23 * t423 + t16 * t424 + t18 * t420 + (t336 * t4 / 0.2e1 - t340 * t3 / 0.2e1) * t330; m(7) * (t35 * t82 + t52 * t97 + t53 * t98) + t355; m(7) * (t305 * t98 + t307 * t97 - t406 * t82); m(7) * (t70 * t82 + t76 * t97 + t77 * t98) + t355; t263 * t4 + t261 * t3 + t283 * t18 + m(7) * (t82 ^ 2 + t97 ^ 2 + t98 ^ 2);];
%% Postprocessing: Reshape Output
% From vec2symmat_6_matlab.m
res = [t1(1) t1(2) t1(4) t1(7) t1(11) t1(16); t1(2) t1(3) t1(5) t1(8) t1(12) t1(17); t1(4) t1(5) t1(6) t1(9) t1(13) t1(18); t1(7) t1(8) t1(9) t1(10) t1(14) t1(19); t1(11) t1(12) t1(13) t1(14) t1(15) t1(20); t1(16) t1(17) t1(18) t1(19) t1(20) t1(21);];
Mq  = res;
