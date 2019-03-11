% Calculate joint inertia matrix for
% S6RRRRPP4
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% pkin [10x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d2,d3,d4,theta5]';
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
% Datum: 2019-03-09 21:04
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Mq = S6RRRRPP4_inertiaJ_slag_vp1(qJ, ...
  pkin, m, rSges, Icges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(10,1),zeros(7,1),zeros(7,3),zeros(7,6)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRRRPP4_inertiaJ_slag_vp1: qJ has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S6RRRRPP4_inertiaJ_slag_vp1: pkin has to be [10x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6RRRRPP4_inertiaJ_slag_vp1: m has to be [7x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [7,3]), ...
  'S6RRRRPP4_inertiaJ_slag_vp1: rSges has to be [7x3] (double)');
assert(isreal(Icges) && all(size(Icges) == [7 6]), ...
  'S6RRRRPP4_inertiaJ_slag_vp1: Icges has to be [7x6] (double)'); 

%% Symbolic Calculation
% From inertia_joint_joint_floatb_twist_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 20:58:30
% EndTime: 2019-03-09 20:58:44
% DurationCPUTime: 6.94s
% Computational Cost: add. (15617->527), mult. (17205->743), div. (0->0), fcn. (18609->10), ass. (0->266)
t288 = qJ(3) + qJ(4);
t276 = pkin(10) + t288;
t273 = sin(t276);
t274 = cos(t276);
t294 = cos(qJ(1));
t291 = sin(qJ(1));
t293 = cos(qJ(2));
t353 = t291 * t293;
t217 = t273 * t353 + t274 * t294;
t218 = -t273 * t294 + t274 * t353;
t290 = sin(qJ(2));
t356 = t290 * t291;
t138 = Icges(7,5) * t218 + Icges(7,6) * t356 + Icges(7,3) * t217;
t144 = Icges(6,4) * t218 - Icges(6,2) * t217 + Icges(6,6) * t356;
t403 = t138 - t144;
t352 = t293 * t294;
t219 = t273 * t352 - t291 * t274;
t220 = t273 * t291 + t274 * t352;
t355 = t290 * t294;
t139 = Icges(7,5) * t220 + Icges(7,6) * t355 + Icges(7,3) * t219;
t145 = Icges(6,4) * t220 - Icges(6,2) * t219 + Icges(6,6) * t355;
t402 = t139 - t145;
t146 = Icges(7,1) * t218 + Icges(7,4) * t356 + Icges(7,5) * t217;
t148 = Icges(6,1) * t218 - Icges(6,4) * t217 + Icges(6,5) * t356;
t401 = t146 + t148;
t147 = Icges(7,1) * t220 + Icges(7,4) * t355 + Icges(7,5) * t219;
t149 = Icges(6,1) * t220 - Icges(6,4) * t219 + Icges(6,5) * t355;
t400 = t147 + t149;
t140 = Icges(6,5) * t218 - Icges(6,6) * t217 + Icges(6,3) * t356;
t142 = Icges(7,4) * t218 + Icges(7,2) * t356 + Icges(7,6) * t217;
t277 = sin(t288);
t278 = cos(t288);
t234 = -t277 * t353 - t278 * t294;
t235 = -t277 * t294 + t278 * t353;
t158 = Icges(5,5) * t235 + Icges(5,6) * t234 + Icges(5,3) * t356;
t399 = t140 + t142 + t158;
t141 = Icges(6,5) * t220 - Icges(6,6) * t219 + Icges(6,3) * t355;
t143 = Icges(7,4) * t220 + Icges(7,2) * t355 + Icges(7,6) * t219;
t236 = -t277 * t352 + t278 * t291;
t237 = t277 * t291 + t278 * t352;
t159 = Icges(5,5) * t237 + Icges(5,6) * t236 + Icges(5,3) * t355;
t398 = t141 + t143 + t159;
t384 = rSges(7,3) + qJ(6);
t385 = rSges(7,1) + pkin(5);
t397 = -t384 * t217 - t385 * t218;
t201 = -Icges(6,6) * t293 + (Icges(6,4) * t274 - Icges(6,2) * t273) * t290;
t360 = t273 * t290;
t214 = -Icges(5,6) * t293 + (Icges(5,4) * t278 - Icges(5,2) * t277) * t290;
t362 = t214 * t277;
t198 = -Icges(7,6) * t293 + (Icges(7,5) * t274 + Icges(7,3) * t273) * t290;
t202 = -Icges(7,4) * t293 + (Icges(7,1) * t274 + Icges(7,5) * t273) * t290;
t203 = -Icges(6,5) * t293 + (Icges(6,1) * t274 - Icges(6,4) * t273) * t290;
t215 = -Icges(5,5) * t293 + (Icges(5,1) * t278 - Icges(5,4) * t277) * t290;
t387 = t198 * t360 + (t278 * t215 + (t202 + t203) * t274) * t290;
t199 = -Icges(6,3) * t293 + (Icges(6,5) * t274 - Icges(6,6) * t273) * t290;
t200 = -Icges(7,2) * t293 + (Icges(7,4) * t274 + Icges(7,6) * t273) * t290;
t213 = -Icges(5,3) * t293 + (Icges(5,5) * t278 - Icges(5,6) * t277) * t290;
t389 = -t199 - t200 - t213;
t375 = -t201 * t360 - t290 * t362 + t293 * t389 + t387;
t396 = t375 * t293;
t160 = Icges(5,4) * t235 + Icges(5,2) * t234 + Icges(5,6) * t356;
t162 = Icges(5,1) * t235 + Icges(5,4) * t234 + Icges(5,5) * t356;
t395 = t160 * t234 + t162 * t235 + t217 * t403 + t401 * t218 + t399 * t356;
t161 = Icges(5,4) * t237 + Icges(5,2) * t236 + Icges(5,6) * t355;
t163 = Icges(5,1) * t237 + Icges(5,4) * t236 + Icges(5,5) * t355;
t394 = t161 * t234 + t163 * t235 + t217 * t402 + t218 * t400 + t356 * t398;
t393 = t160 * t236 + t162 * t237 + t219 * t403 + t401 * t220 + t399 * t355;
t392 = t161 * t236 + t163 * t237 + t219 * t402 + t220 * t400 + t355 * t398;
t104 = t213 * t356 + t214 * t234 + t215 * t235;
t97 = t198 * t217 + t200 * t356 + t202 * t218;
t98 = t199 * t356 - t201 * t217 + t203 * t218;
t391 = -t98 - t104 - t97;
t100 = t199 * t355 - t201 * t219 + t203 * t220;
t105 = t213 * t355 + t214 * t236 + t215 * t237;
t99 = t198 * t219 + t200 * t355 + t202 * t220;
t390 = -t99 - t100 - t105;
t388 = Icges(3,5) * t290;
t386 = t388 / 0.2e1;
t383 = rSges(7,2) * t356 - t397;
t382 = -rSges(7,2) * t293 + (t384 * t273 + t274 * t385) * t290;
t381 = t391 * t293 + (t291 * t395 + t294 * t394) * t290;
t380 = t390 * t293 + (t291 * t393 + t294 * t392) * t290;
t379 = t291 * t394 - t294 * t395;
t378 = t291 * t392 - t294 * t393;
t73 = -t142 * t293 + (t138 * t273 + t146 * t274) * t290;
t75 = -t140 * t293 + (-t144 * t273 + t148 * t274) * t290;
t79 = -t158 * t293 + (-t160 * t277 + t162 * t278) * t290;
t377 = -t73 - t75 - t79;
t74 = -t143 * t293 + (t139 * t273 + t147 * t274) * t290;
t76 = -t141 * t293 + (-t145 * t273 + t149 * t274) * t290;
t80 = -t159 * t293 + (-t161 * t277 + t163 * t278) * t290;
t376 = t74 + t76 + t80;
t286 = t291 ^ 2;
t287 = t294 ^ 2;
t374 = m(6) / 0.2e1;
t373 = m(7) / 0.2e1;
t295 = -pkin(9) - pkin(8);
t372 = t291 / 0.2e1;
t371 = -t293 / 0.2e1;
t370 = -t294 / 0.2e1;
t369 = t294 / 0.2e1;
t368 = pkin(2) * t293;
t367 = pkin(8) * t290;
t292 = cos(qJ(3));
t275 = t292 * pkin(3) + pkin(2);
t366 = -pkin(2) + t275;
t365 = t294 * rSges(3,3);
t363 = Icges(3,4) * t293;
t289 = sin(qJ(3));
t227 = -Icges(4,6) * t293 + (Icges(4,4) * t292 - Icges(4,2) * t289) * t290;
t361 = t227 * t289;
t358 = t289 * t291;
t357 = t289 * t294;
t354 = t290 * t295;
t336 = pkin(3) * t357 + t291 * t354;
t250 = pkin(4) * t278 + t275;
t338 = t250 - t275;
t252 = pkin(3) * t289 + pkin(4) * t277;
t284 = -qJ(5) + t295;
t339 = -t294 * t252 - t284 * t356;
t128 = t338 * t353 + t336 + t339;
t121 = t128 * t355;
t309 = -rSges(6,1) * t218 + rSges(6,2) * t217;
t151 = rSges(6,3) * t356 - t309;
t351 = t151 * t355 + t121;
t333 = t284 - t295;
t184 = t338 * t290 + t333 * t293;
t350 = t293 * t128 + t184 * t356;
t337 = -pkin(3) * t358 - t275 * t352;
t341 = t250 * t352 + t291 * t252;
t129 = -t333 * t355 + t337 + t341;
t153 = t220 * rSges(6,1) - t219 * rSges(6,2) + rSges(6,3) * t355;
t349 = -t129 - t153;
t348 = rSges(7,2) * t355 + t384 * t219 + t385 * t220;
t165 = t237 * rSges(5,1) + t236 * rSges(5,2) + rSges(5,3) * t355;
t301 = -t294 * t354 - t337;
t335 = pkin(2) * t352 + pkin(8) * t355;
t183 = t301 - t335;
t347 = -t165 - t183;
t182 = (t366 * t293 - t367) * t291 - t336;
t207 = (pkin(8) + t295) * t293 + t366 * t290;
t346 = t293 * t182 + t207 * t356;
t205 = -rSges(6,3) * t293 + (rSges(6,1) * t274 - rSges(6,2) * t273) * t290;
t345 = -t184 - t205;
t310 = -rSges(5,1) * t235 - rSges(5,2) * t234;
t164 = rSges(5,3) * t356 - t310;
t216 = -rSges(5,3) * t293 + (rSges(5,1) * t278 - rSges(5,2) * t277) * t290;
t119 = t293 * t164 + t216 * t356;
t343 = -t207 - t216;
t233 = -rSges(4,3) * t293 + (rSges(4,1) * t292 - rSges(4,2) * t289) * t290;
t262 = t290 * pkin(2) - t293 * pkin(8);
t342 = -t233 - t262;
t340 = t286 * (t367 + t368) + t294 * t335;
t334 = t294 * pkin(1) + t291 * pkin(7);
t332 = t286 + t287;
t331 = t396 + (t377 * t291 - t376 * t294) * t290;
t224 = -Icges(4,3) * t293 + (Icges(4,5) * t292 - Icges(4,6) * t289) * t290;
t230 = -Icges(4,5) * t293 + (Icges(4,1) * t292 - Icges(4,4) * t289) * t290;
t244 = -t289 * t353 - t292 * t294;
t245 = t292 * t353 - t357;
t113 = t224 * t356 + t227 * t244 + t230 * t245;
t174 = Icges(4,5) * t245 + Icges(4,6) * t244 + Icges(4,3) * t356;
t176 = Icges(4,4) * t245 + Icges(4,2) * t244 + Icges(4,6) * t356;
t178 = Icges(4,1) * t245 + Icges(4,4) * t244 + Icges(4,5) * t356;
t89 = -t174 * t293 + (-t176 * t289 + t178 * t292) * t290;
t330 = t113 / 0.2e1 + t89 / 0.2e1;
t246 = -t289 * t352 + t291 * t292;
t247 = t292 * t352 + t358;
t114 = t224 * t355 + t227 * t246 + t230 * t247;
t175 = Icges(4,5) * t247 + Icges(4,6) * t246 + Icges(4,3) * t355;
t177 = Icges(4,4) * t247 + Icges(4,2) * t246 + Icges(4,6) * t355;
t179 = Icges(4,1) * t247 + Icges(4,4) * t246 + Icges(4,5) * t355;
t90 = -t175 * t293 + (-t177 * t289 + t179 * t292) * t290;
t329 = t90 / 0.2e1 + t114 / 0.2e1;
t327 = t383 * t355 + t121;
t326 = -t129 - t348;
t325 = -t183 + t349;
t324 = -t184 - t382;
t323 = -t207 + t345;
t322 = -t262 + t343;
t181 = t247 * rSges(4,1) + t246 * rSges(4,2) + rSges(4,3) * t355;
t282 = t294 * pkin(7);
t321 = t282 - t339;
t320 = t356 / 0.2e1;
t319 = t355 / 0.2e1;
t318 = -t250 * t293 - pkin(1);
t317 = -t183 + t326;
t316 = t291 * t182 + t294 * t183 + t340;
t62 = t293 * t151 + t205 * t356 + t350;
t315 = -t207 + t324;
t314 = -t262 + t323;
t313 = t380 * t355 + t381 * t356;
t312 = rSges(3,1) * t293 - rSges(3,2) * t290;
t311 = -rSges(4,1) * t245 - rSges(4,2) * t244;
t308 = -t262 + t315;
t306 = -Icges(3,2) * t290 + t363;
t305 = Icges(3,5) * t293 - Icges(3,6) * t290;
t302 = rSges(3,1) * t352 - rSges(3,2) * t355 + t291 * rSges(3,3);
t300 = t291 * t128 + t294 * t129 + t316;
t48 = t383 * t293 + t382 * t356 + t350;
t299 = -t284 * t355 + t334 + t341;
t298 = t331 * t293 + t313;
t297 = (-t377 - t391) * t320 + (t376 - t390) * t319;
t296 = t380 * t372 + (t376 * t291 + t377 * t294) * t371 + t381 * t370 + t379 * t320 + t378 * t319;
t285 = t290 ^ 2;
t260 = rSges(2,1) * t294 - rSges(2,2) * t291;
t259 = -rSges(2,1) * t291 - rSges(2,2) * t294;
t258 = rSges(3,1) * t290 + rSges(3,2) * t293;
t254 = Icges(3,6) * t293 + t388;
t226 = Icges(3,3) * t291 + t305 * t294;
t225 = -Icges(3,3) * t294 + t305 * t291;
t206 = t290 * t292 * t230;
t197 = t302 + t334;
t196 = t365 + t282 + (-pkin(1) - t312) * t291;
t189 = t342 * t294;
t188 = t342 * t291;
t180 = rSges(4,3) * t356 - t311;
t169 = t294 * t302 + (t312 * t291 - t365) * t291;
t168 = t182 * t355;
t154 = t164 * t355;
t133 = t181 + t334 + t335;
t132 = t282 + (-t368 - pkin(1) + (-rSges(4,3) - pkin(8)) * t290) * t291 + t311;
t131 = t322 * t294;
t130 = t322 * t291;
t127 = -t181 * t293 - t233 * t355;
t126 = t180 * t293 + t233 * t356;
t122 = -t293 * t224 - t290 * t361 + t206;
t120 = -t165 * t293 - t216 * t355;
t118 = t301 + t165 + t334;
t117 = t282 + (-rSges(5,3) * t290 - t275 * t293 - pkin(1)) * t291 + t310 + t336;
t115 = (t180 * t294 - t181 * t291) * t290;
t110 = t299 + t153;
t109 = (-rSges(6,3) * t290 + t318) * t291 + t309 + t321;
t108 = t314 * t294;
t107 = t314 * t291;
t106 = -t165 * t356 + t154;
t101 = t180 * t291 + t181 * t294 + t340;
t96 = t308 * t294;
t95 = t308 * t291;
t88 = t299 + t348;
t87 = (-rSges(7,2) * t290 + t318) * t291 + t321 + t397;
t86 = t347 * t293 + t343 * t355;
t85 = t119 + t346;
t84 = t175 * t355 + t177 * t246 + t179 * t247;
t83 = t174 * t355 + t176 * t246 + t178 * t247;
t82 = t175 * t356 + t177 * t244 + t179 * t245;
t81 = t174 * t356 + t176 * t244 + t178 * t245;
t64 = t347 * t356 + t154 + t168;
t63 = t349 * t293 + t345 * t355;
t53 = t164 * t291 + t165 * t294 + t316;
t52 = t349 * t356 + t351;
t51 = t325 * t293 + t323 * t355;
t50 = t62 + t346;
t49 = t326 * t293 + t324 * t355;
t47 = t291 * t84 - t294 * t83;
t46 = t291 * t82 - t294 * t81;
t44 = t293 * t317 + t315 * t355;
t43 = t48 + t346;
t42 = t325 * t356 + t168 + t351;
t41 = t326 * t356 + t327;
t40 = t151 * t291 + t153 * t294 + t300;
t33 = -t114 * t293 + (t291 * t83 + t294 * t84) * t290;
t32 = -t113 * t293 + (t291 * t81 + t294 * t82) * t290;
t22 = t317 * t356 + t168 + t327;
t21 = t383 * t291 + t348 * t294 + t300;
t1 = [Icges(2,3) + t206 + (Icges(3,1) * t290 - t201 * t273 - t361 - t362 + t363) * t290 + (Icges(3,4) * t290 + Icges(3,2) * t293 - t224 + t389) * t293 + m(6) * (t109 ^ 2 + t110 ^ 2) + m(7) * (t87 ^ 2 + t88 ^ 2) + m(5) * (t117 ^ 2 + t118 ^ 2) + m(4) * (t132 ^ 2 + t133 ^ 2) + m(3) * (t196 ^ 2 + t197 ^ 2) + m(2) * (t259 ^ 2 + t260 ^ 2) + t387; ((-Icges(3,6) * t294 + t306 * t291) * t371 + t294 * t386 - t97 / 0.2e1 - t98 / 0.2e1 - t104 / 0.2e1 - t79 / 0.2e1 - t75 / 0.2e1 - t73 / 0.2e1 + t254 * t369 - t330) * t294 + ((Icges(3,6) * t291 + t306 * t294) * t293 / 0.2e1 + t291 * t386 + t80 / 0.2e1 + t76 / 0.2e1 + t74 / 0.2e1 + t99 / 0.2e1 + t100 / 0.2e1 + t105 / 0.2e1 + t254 * t372 + t329) * t291 + m(5) * (t117 * t131 + t118 * t130) + m(6) * (t107 * t110 + t108 * t109) + m(7) * (t87 * t96 + t88 * t95) + m(4) * (t132 * t189 + t133 * t188) + m(3) * (-t196 * t294 - t197 * t291) * t258; m(7) * (t21 ^ 2 + t95 ^ 2 + t96 ^ 2) + m(6) * (t107 ^ 2 + t108 ^ 2 + t40 ^ 2) + m(5) * (t130 ^ 2 + t131 ^ 2 + t53 ^ 2) + m(4) * (t101 ^ 2 + t188 ^ 2 + t189 ^ 2) + m(3) * (t258 ^ 2 * t332 + t169 ^ 2) + (-t287 * t225 - t379 - t46) * t294 + (t286 * t226 + t47 + (-t291 * t225 + t294 * t226) * t294 + t378) * t291; (t291 * t330 + t294 * t329) * t290 + m(6) * (t109 * t50 + t110 * t51) + m(7) * (t43 * t87 + t44 * t88) + m(5) * (t117 * t85 + t118 * t86) + m(4) * (t126 * t132 + t127 * t133) + (-t122 - t375) * t293 + t297; t296 + m(7) * (t21 * t22 + t43 * t96 + t44 * t95) + m(6) * (t107 * t51 + t108 * t50 + t40 * t42) + m(5) * (t130 * t86 + t131 * t85 + t53 * t64) + m(4) * (t101 * t115 + t126 * t189 + t127 * t188) + (t369 * t47 + t372 * t46) * t290 + (t90 * t291 - t89 * t294) * t371 + t32 * t370 + t33 * t372; (t291 * t32 + t294 * t33) * t290 + (t122 * t293 + (-t291 * t89 - t294 * t90) * t290 + t331) * t293 + m(5) * (t64 ^ 2 + t85 ^ 2 + t86 ^ 2) + m(4) * (t115 ^ 2 + t126 ^ 2 + t127 ^ 2) + m(6) * (t42 ^ 2 + t50 ^ 2 + t51 ^ 2) + m(7) * (t22 ^ 2 + t43 ^ 2 + t44 ^ 2) + t313; -t396 + m(6) * (t109 * t62 + t110 * t63) + m(7) * (t48 * t87 + t49 * t88) + m(5) * (t117 * t119 + t118 * t120) + t297; m(7) * (t21 * t41 + t48 * t96 + t49 * t95) + m(6) * (t107 * t63 + t108 * t62 + t40 * t52) + m(5) * (t106 * t53 + t119 * t131 + t120 * t130) + t296; m(6) * (t42 * t52 + t50 * t62 + t51 * t63) + m(7) * (t22 * t41 + t43 * t48 + t44 * t49) + m(5) * (t106 * t64 + t119 * t85 + t120 * t86) + t298; m(7) * (t41 ^ 2 + t48 ^ 2 + t49 ^ 2) + m(6) * (t52 ^ 2 + t62 ^ 2 + t63 ^ 2) + m(5) * (t106 ^ 2 + t119 ^ 2 + t120 ^ 2) + t298; 0.2e1 * ((t109 * t294 + t110 * t291) * t374 + (t291 * t88 + t294 * t87) * t373) * t290; m(7) * (-t21 * t293 + (t291 * t95 + t294 * t96) * t290) + m(6) * (-t293 * t40 + (t107 * t291 + t108 * t294) * t290); m(6) * (-t293 * t42 + (t291 * t51 + t294 * t50) * t290) + m(7) * (-t22 * t293 + (t291 * t44 + t294 * t43) * t290); m(7) * (-t293 * t41 + (t291 * t49 + t294 * t48) * t290) + m(6) * (-t293 * t52 + (t291 * t63 + t294 * t62) * t290); 0.2e1 * (t374 + t373) * (t285 * t332 + t293 ^ 2); m(7) * (t217 * t88 + t219 * t87); m(7) * (t21 * t360 + t217 * t95 + t219 * t96); m(7) * (t217 * t44 + t219 * t43 + t22 * t360); m(7) * (t217 * t49 + t219 * t48 + t360 * t41); m(7) * (t217 * t291 + t219 * t294 - t273 * t293) * t290; m(7) * (t273 ^ 2 * t285 + t217 ^ 2 + t219 ^ 2);];
%% Postprocessing: Reshape Output
% From vec2symmat_6_matlab.m
res = [t1(1) t1(2) t1(4) t1(7) t1(11) t1(16); t1(2) t1(3) t1(5) t1(8) t1(12) t1(17); t1(4) t1(5) t1(6) t1(9) t1(13) t1(18); t1(7) t1(8) t1(9) t1(10) t1(14) t1(19); t1(11) t1(12) t1(13) t1(14) t1(15) t1(20); t1(16) t1(17) t1(18) t1(19) t1(20) t1(21);];
Mq  = res;
