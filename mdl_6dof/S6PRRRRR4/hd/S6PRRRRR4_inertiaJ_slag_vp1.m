% Calculate joint inertia matrix for
% S6PRRRRR4
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% pkin [13x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,alpha2,alpha3,d2,d3,d4,d5,d6,theta1]';
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
% Datum: 2019-03-09 01:02
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Mq = S6PRRRRR4_inertiaJ_slag_vp1(qJ, ...
  pkin, m, rSges, Icges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(13,1),zeros(7,1),zeros(7,3),zeros(7,6)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6PRRRRR4_inertiaJ_slag_vp1: qJ has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [13 1]), ...
  'S6PRRRRR4_inertiaJ_slag_vp1: pkin has to be [13x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6PRRRRR4_inertiaJ_slag_vp1: m has to be [7x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [7,3]), ...
  'S6PRRRRR4_inertiaJ_slag_vp1: rSges has to be [7x3] (double)');
assert(isreal(Icges) && all(size(Icges) == [7 6]), ...
  'S6PRRRRR4_inertiaJ_slag_vp1: Icges has to be [7x6] (double)'); 

%% Symbolic Calculation
% From inertia_joint_joint_floatb_twist_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 00:55:00
% EndTime: 2019-03-09 00:55:18
% DurationCPUTime: 8.58s
% Computational Cost: add. (69714->680), mult. (166316->964), div. (0->0), fcn. (218766->16), ass. (0->333)
t312 = sin(pkin(13));
t314 = cos(pkin(13));
t319 = sin(qJ(2));
t315 = cos(pkin(6));
t322 = cos(qJ(2));
t382 = t315 * t322;
t303 = -t312 * t319 + t314 * t382;
t383 = t315 * t319;
t304 = t312 * t322 + t314 * t383;
t318 = sin(qJ(3));
t313 = sin(pkin(6));
t391 = sin(pkin(7));
t349 = t313 * t391;
t392 = cos(pkin(7));
t395 = cos(qJ(3));
t272 = t304 * t395 + (t303 * t392 - t314 * t349) * t318;
t350 = t313 * t392;
t292 = -t303 * t391 - t314 * t350;
t381 = qJ(4) + qJ(5);
t311 = sin(t381);
t351 = cos(t381);
t249 = t272 * t311 - t292 * t351;
t250 = t272 * t351 + t292 * t311;
t335 = t395 * t391;
t331 = t313 * t335;
t336 = t392 * t395;
t271 = -t303 * t336 + t304 * t318 + t314 * t331;
t187 = Icges(6,5) * t250 - Icges(6,6) * t249 + Icges(6,3) * t271;
t189 = Icges(6,4) * t250 - Icges(6,2) * t249 + Icges(6,6) * t271;
t191 = Icges(6,1) * t250 - Icges(6,4) * t249 + Icges(6,5) * t271;
t102 = t187 * t271 - t189 * t249 + t191 * t250;
t305 = -t312 * t382 - t314 * t319;
t306 = -t312 * t383 + t314 * t322;
t274 = t306 * t395 + (t305 * t392 + t312 * t349) * t318;
t293 = -t305 * t391 + t312 * t350;
t251 = t274 * t311 - t293 * t351;
t252 = t274 * t351 + t293 * t311;
t273 = -t305 * t336 + t306 * t318 - t312 * t331;
t188 = Icges(6,5) * t252 - Icges(6,6) * t251 + Icges(6,3) * t273;
t190 = Icges(6,4) * t252 - Icges(6,2) * t251 + Icges(6,6) * t273;
t192 = Icges(6,1) * t252 - Icges(6,4) * t251 + Icges(6,5) * t273;
t103 = t188 * t271 - t190 * t249 + t192 * t250;
t291 = t315 * t391 * t318 + (t318 * t322 * t392 + t319 * t395) * t313;
t302 = t315 * t392 - t322 * t349;
t265 = t291 * t311 - t302 * t351;
t266 = t291 * t351 + t302 * t311;
t384 = t313 * t319;
t290 = -t313 * t322 * t336 - t315 * t335 + t318 * t384;
t218 = Icges(6,5) * t266 - Icges(6,6) * t265 + Icges(6,3) * t290;
t219 = Icges(6,4) * t266 - Icges(6,2) * t265 + Icges(6,6) * t290;
t220 = Icges(6,1) * t266 - Icges(6,4) * t265 + Icges(6,5) * t290;
t120 = t218 * t271 - t219 * t249 + t220 * t250;
t316 = sin(qJ(6));
t320 = cos(qJ(6));
t213 = -t250 * t316 + t271 * t320;
t214 = t250 * t320 + t271 * t316;
t146 = Icges(7,5) * t214 + Icges(7,6) * t213 + Icges(7,3) * t249;
t148 = Icges(7,4) * t214 + Icges(7,2) * t213 + Icges(7,6) * t249;
t150 = Icges(7,1) * t214 + Icges(7,4) * t213 + Icges(7,5) * t249;
t76 = t146 * t249 + t148 * t213 + t150 * t214;
t215 = -t252 * t316 + t273 * t320;
t216 = t252 * t320 + t273 * t316;
t147 = Icges(7,5) * t216 + Icges(7,6) * t215 + Icges(7,3) * t251;
t149 = Icges(7,4) * t216 + Icges(7,2) * t215 + Icges(7,6) * t251;
t151 = Icges(7,1) * t216 + Icges(7,4) * t215 + Icges(7,5) * t251;
t77 = t147 * t249 + t149 * t213 + t151 * t214;
t247 = -t266 * t316 + t290 * t320;
t248 = t266 * t320 + t290 * t316;
t183 = Icges(7,5) * t248 + Icges(7,6) * t247 + Icges(7,3) * t265;
t184 = Icges(7,4) * t248 + Icges(7,2) * t247 + Icges(7,6) * t265;
t185 = Icges(7,1) * t248 + Icges(7,4) * t247 + Icges(7,5) * t265;
t92 = t183 * t249 + t184 * t213 + t185 * t214;
t13 = t271 * t76 + t273 * t77 + t290 * t92;
t423 = t102 * t271 + t103 * t273 + t120 * t290 + t13;
t104 = t187 * t273 - t189 * t251 + t191 * t252;
t105 = t188 * t273 - t190 * t251 + t192 * t252;
t121 = t218 * t273 - t219 * t251 + t220 * t252;
t78 = t146 * t251 + t148 * t215 + t150 * t216;
t79 = t147 * t251 + t149 * t215 + t151 * t216;
t93 = t183 * t251 + t184 * t215 + t185 * t216;
t14 = t271 * t78 + t273 * t79 + t290 * t93;
t422 = t104 * t271 + t105 * t273 + t121 * t290 + t14;
t17 = t292 * t76 + t293 * t77 + t302 * t92;
t44 = t102 * t292 + t103 * t293 + t120 * t302;
t421 = t17 + t44;
t18 = t292 * t78 + t293 * t79 + t302 * t93;
t45 = t104 * t292 + t105 * t293 + t121 * t302;
t420 = t18 + t45;
t21 = t315 * t92 + (t312 * t77 - t314 * t76) * t313;
t50 = t120 * t315 + (-t102 * t314 + t103 * t312) * t313;
t419 = t21 + t50;
t22 = t315 * t93 + (t312 * t79 - t314 * t78) * t313;
t51 = t121 * t315 + (-t104 * t314 + t105 * t312) * t313;
t418 = t22 + t51;
t112 = t187 * t290 - t189 * t265 + t191 * t266;
t113 = t188 * t290 - t190 * t265 + t192 * t266;
t129 = t218 * t290 - t219 * t265 + t220 * t266;
t101 = t183 * t265 + t184 * t247 + t185 * t248;
t84 = t146 * t265 + t148 * t247 + t150 * t248;
t85 = t147 * t265 + t149 * t247 + t151 * t248;
t29 = t101 * t290 + t271 * t84 + t273 * t85;
t417 = t112 * t271 + t113 * t273 + t129 * t290 + t29;
t31 = t101 * t302 + t292 * t84 + t293 * t85;
t61 = t112 * t292 + t113 * t293 + t129 * t302;
t416 = t31 + t61;
t33 = t101 * t315 + (t312 * t85 - t314 * t84) * t313;
t64 = t129 * t315 + (-t112 * t314 + t113 * t312) * t313;
t415 = t33 + t64;
t153 = rSges(7,1) * t216 + rSges(7,2) * t215 + rSges(7,3) * t251;
t377 = pkin(5) * t252 + pkin(12) * t251 + t153;
t186 = rSges(7,1) * t248 + rSges(7,2) * t247 + rSges(7,3) * t265;
t369 = pkin(5) * t266 + pkin(12) * t265 + t186;
t193 = rSges(6,1) * t250 - rSges(6,2) * t249 + rSges(6,3) * t271;
t167 = t273 * t193;
t194 = rSges(6,1) * t252 - rSges(6,2) * t251 + rSges(6,3) * t273;
t128 = -t271 * t194 + t167;
t278 = t304 * pkin(2) + pkin(9) * t292;
t279 = t306 * pkin(2) + pkin(9) * t293;
t385 = t313 * t314;
t386 = t312 * t313;
t361 = t278 * t386 + t279 * t385;
t408 = m(7) / 0.2e1;
t409 = m(6) / 0.2e1;
t414 = t409 + t408;
t317 = sin(qJ(4));
t321 = cos(qJ(4));
t256 = -t274 * t317 + t293 * t321;
t388 = t293 * t317;
t257 = t274 * t321 + t388;
t205 = rSges(5,1) * t257 + rSges(5,2) * t256 + rSges(5,3) * t273;
t232 = rSges(4,1) * t274 - rSges(4,2) * t273 + rSges(4,3) * t293;
t413 = -m(4) * t232 - m(5) * t205 - m(6) * t194 - m(7) * t377;
t412 = -0.2e1 * t271;
t411 = -0.2e1 * t292;
t410 = m(5) / 0.2e1;
t407 = t249 / 0.2e1;
t406 = t251 / 0.2e1;
t405 = t265 / 0.2e1;
t404 = t271 / 0.2e1;
t403 = t273 / 0.2e1;
t402 = t290 / 0.2e1;
t401 = t292 / 0.2e1;
t400 = t293 / 0.2e1;
t399 = t302 / 0.2e1;
t398 = t312 / 0.2e1;
t397 = -t314 / 0.2e1;
t396 = t315 / 0.2e1;
t394 = pkin(4) * t321;
t389 = t292 * t317;
t387 = t302 * t317;
t152 = rSges(7,1) * t214 + rSges(7,2) * t213 + rSges(7,3) * t249;
t144 = t273 * t152;
t209 = pkin(5) * t250 + pkin(12) * t249;
t197 = t273 * t209;
t380 = t144 + t197;
t379 = t377 * t290;
t378 = t152 + t209;
t181 = pkin(4) * t389 + pkin(11) * t271 + t272 * t394;
t163 = t273 * t181;
t182 = pkin(4) * t388 + pkin(11) * t273 + t274 * t394;
t376 = t182 * t412 + 0.2e1 * t163;
t375 = t369 * t271;
t374 = 0.2e1 * t128;
t171 = t293 * t181;
t244 = t272 * pkin(3) + t271 * pkin(10);
t236 = t293 * t244;
t373 = t171 + t236;
t245 = t274 * pkin(3) + t273 * pkin(10);
t237 = t302 * t245;
t372 = t302 * t182 + t237;
t371 = -t181 - t193;
t370 = -t182 - t194;
t254 = -t272 * t317 + t292 * t321;
t255 = t272 * t321 + t389;
t204 = rSges(5,1) * t255 + rSges(5,2) * t254 + rSges(5,3) * t271;
t368 = -t204 - t244;
t217 = pkin(4) * t387 + pkin(11) * t290 + t291 * t394;
t264 = t291 * pkin(3) + t290 * pkin(10);
t253 = t292 * t264;
t367 = t292 * t217 + t253;
t221 = rSges(6,1) * t266 - rSges(6,2) * t265 + rSges(6,3) * t290;
t366 = -t217 - t221;
t275 = -t291 * t317 + t302 * t321;
t276 = t291 * t321 + t387;
t233 = rSges(5,1) * t276 + rSges(5,2) * t275 + rSges(5,3) * t290;
t365 = -t233 - t264;
t364 = t245 * t411 + 0.2e1 * t236;
t277 = t315 * t279;
t363 = t315 * t245 + t277;
t362 = 0.2e1 * t361;
t359 = -t181 - t378;
t358 = -t182 - t377;
t357 = t315 * t182 + t363;
t356 = -t244 + t371;
t355 = -t217 - t369;
t354 = -t264 + t366;
t261 = rSges(4,1) * t291 - rSges(4,2) * t290 + rSges(4,3) * t302;
t294 = pkin(2) * t384 + pkin(9) * t302;
t344 = (-t261 - t294) * t313;
t343 = -t244 + t359;
t341 = -t264 + t355;
t340 = t377 * t412 + 0.2e1 * t144 + 0.2e1 * t197;
t241 = t244 * t386;
t242 = t245 * t385;
t339 = 0.2e1 * t241 + 0.2e1 * t242 + t362;
t338 = t241 + t242 + t361;
t337 = (-t294 + t365) * t313;
t24 = t101 * t265 + t249 * t84 + t251 * t85;
t3 = t249 * t76 + t251 * t77 + t265 * t92;
t4 = t249 * t78 + t251 * t79 + t265 * t93;
t334 = t13 * t407 + t14 * t406 + t24 * t402 + t29 * t405 + t3 * t404 + t4 * t403;
t333 = t271 * t423 + t422 * t273 + t417 * t290;
t100 = t152 * t251 - t153 * t249;
t130 = t204 * t273 - t205 * t271;
t332 = (-t294 + t354) * t313;
t178 = t181 * t386;
t179 = t182 * t385;
t329 = t178 + t179 + t338;
t328 = (-t294 + t341) * t313;
t326 = t417 * t399 + t422 * t400 + t401 * t423 + t416 * t402 + t420 * t403 + t421 * t404;
t325 = t419 * t404 + t418 * t403 + t415 * t402 + t417 * t396 + t422 * t386 / 0.2e1 - t423 * t385 / 0.2e1;
t231 = rSges(4,1) * t272 - rSges(4,2) * t271 + rSges(4,3) * t292;
t324 = m(4) * t231 + m(5) * t204 + m(6) * t193 + m(7) * t378;
t301 = rSges(3,3) * t315 + (rSges(3,1) * t319 + rSges(3,2) * t322) * t313;
t300 = Icges(3,5) * t315 + (Icges(3,1) * t319 + Icges(3,4) * t322) * t313;
t299 = Icges(3,6) * t315 + (Icges(3,4) * t319 + Icges(3,2) * t322) * t313;
t298 = Icges(3,3) * t315 + (Icges(3,5) * t319 + Icges(3,6) * t322) * t313;
t287 = rSges(3,1) * t306 + rSges(3,2) * t305 + rSges(3,3) * t386;
t286 = rSges(3,1) * t304 + rSges(3,2) * t303 - rSges(3,3) * t385;
t285 = Icges(3,1) * t306 + Icges(3,4) * t305 + Icges(3,5) * t386;
t284 = Icges(3,1) * t304 + Icges(3,4) * t303 - Icges(3,5) * t385;
t283 = Icges(3,4) * t306 + Icges(3,2) * t305 + Icges(3,6) * t386;
t282 = Icges(3,4) * t304 + Icges(3,2) * t303 - Icges(3,6) * t385;
t281 = Icges(3,5) * t306 + Icges(3,6) * t305 + Icges(3,3) * t386;
t280 = Icges(3,5) * t304 + Icges(3,6) * t303 - Icges(3,3) * t385;
t263 = -t286 * t315 - t301 * t385;
t262 = t287 * t315 - t301 * t386;
t260 = Icges(4,1) * t291 - Icges(4,4) * t290 + Icges(4,5) * t302;
t259 = Icges(4,4) * t291 - Icges(4,2) * t290 + Icges(4,6) * t302;
t258 = Icges(4,5) * t291 - Icges(4,6) * t290 + Icges(4,3) * t302;
t246 = (t286 * t312 + t287 * t314) * t313;
t230 = Icges(5,1) * t276 + Icges(5,4) * t275 + Icges(5,5) * t290;
t229 = Icges(5,4) * t276 + Icges(5,2) * t275 + Icges(5,6) * t290;
t228 = Icges(5,5) * t276 + Icges(5,6) * t275 + Icges(5,3) * t290;
t227 = Icges(4,1) * t274 - Icges(4,4) * t273 + Icges(4,5) * t293;
t226 = Icges(4,1) * t272 - Icges(4,4) * t271 + Icges(4,5) * t292;
t225 = Icges(4,4) * t274 - Icges(4,2) * t273 + Icges(4,6) * t293;
t224 = Icges(4,4) * t272 - Icges(4,2) * t271 + Icges(4,6) * t292;
t223 = Icges(4,5) * t274 - Icges(4,6) * t273 + Icges(4,3) * t293;
t222 = Icges(4,5) * t272 - Icges(4,6) * t271 + Icges(4,3) * t292;
t208 = t271 * t221;
t207 = t271 * t217;
t203 = Icges(5,1) * t257 + Icges(5,4) * t256 + Icges(5,5) * t273;
t202 = Icges(5,1) * t255 + Icges(5,4) * t254 + Icges(5,5) * t271;
t201 = Icges(5,4) * t257 + Icges(5,2) * t256 + Icges(5,6) * t273;
t200 = Icges(5,4) * t255 + Icges(5,2) * t254 + Icges(5,6) * t271;
t199 = Icges(5,5) * t257 + Icges(5,6) * t256 + Icges(5,3) * t273;
t198 = Icges(5,5) * t255 + Icges(5,6) * t254 + Icges(5,3) * t271;
t175 = t232 * t302 - t261 * t293;
t174 = -t231 * t302 + t261 * t292;
t172 = t290 * t194;
t170 = t290 * t182;
t160 = (-t231 - t278) * t315 + t314 * t344;
t159 = t232 * t315 + t312 * t344 + t277;
t158 = t258 * t302 - t259 * t290 + t260 * t291;
t157 = t231 * t293 - t232 * t292;
t156 = t258 * t293 - t259 * t273 + t260 * t274;
t155 = t258 * t292 - t259 * t271 + t260 * t272;
t154 = (t231 * t312 + t232 * t314) * t313 + t361;
t141 = t205 * t290 - t233 * t273;
t140 = -t204 * t290 + t233 * t271;
t139 = -t221 * t273 + t172;
t138 = -t193 * t290 + t208;
t137 = t223 * t302 - t225 * t290 + t227 * t291;
t136 = t222 * t302 - t224 * t290 + t226 * t291;
t135 = t228 * t290 + t229 * t275 + t230 * t276;
t134 = t223 * t293 - t225 * t273 + t227 * t274;
t133 = t222 * t293 - t224 * t273 + t226 * t274;
t132 = t223 * t292 - t225 * t271 + t227 * t272;
t131 = t222 * t292 - t224 * t271 + t226 * t272;
t127 = t205 * t302 + t293 * t365 + t237;
t126 = t233 * t292 + t302 * t368 + t253;
t125 = (-t278 + t368) * t315 + t314 * t337;
t124 = t205 * t315 + t312 * t337 + t363;
t123 = t228 * t273 + t229 * t256 + t230 * t257;
t122 = t228 * t271 + t229 * t254 + t230 * t255;
t119 = t153 * t265 - t186 * t251;
t118 = -t152 * t265 + t186 * t249;
t117 = t204 * t293 + t236 + (-t205 - t245) * t292;
t116 = (t204 * t312 + t205 * t314) * t313 + t338;
t115 = t199 * t290 + t201 * t275 + t203 * t276;
t114 = t198 * t290 + t200 * t275 + t202 * t276;
t111 = t199 * t273 + t201 * t256 + t203 * t257;
t110 = t198 * t273 + t200 * t256 + t202 * t257;
t109 = t199 * t271 + t201 * t254 + t203 * t255;
t108 = t198 * t271 + t200 * t254 + t202 * t255;
t107 = t273 * t366 + t170 + t172;
t106 = t290 * t371 + t207 + t208;
t99 = (-t278 + t356) * t315 + t314 * t332;
t98 = t194 * t315 + t312 * t332 + t357;
t97 = t194 * t302 + t293 * t354 + t372;
t96 = t221 * t292 + t302 * t356 + t367;
t95 = -t273 * t369 + t379;
t94 = -t290 * t378 + t375;
t91 = t271 * t370 + t163 + t167;
t90 = t158 * t315 + (-t136 * t314 + t137 * t312) * t313;
t89 = (t193 * t312 + t194 * t314) * t313 + t329;
t88 = t136 * t292 + t137 * t293 + t158 * t302;
t87 = -t271 * t377 + t380;
t86 = t193 * t293 + (-t245 + t370) * t292 + t373;
t83 = t156 * t315 + (-t133 * t314 + t134 * t312) * t313;
t82 = t155 * t315 + (-t131 * t314 + t132 * t312) * t313;
t81 = t133 * t292 + t134 * t293 + t156 * t302;
t80 = t131 * t292 + t132 * t293 + t155 * t302;
t75 = t273 * t355 + t170 + t379;
t74 = t290 * t359 + t207 + t375;
t73 = (-t278 + t343) * t315 + t314 * t328;
t72 = t312 * t328 + t315 * t377 + t357;
t71 = t293 * t341 + t302 * t377 + t372;
t70 = t292 * t369 + t302 * t343 + t367;
t69 = t271 * t358 + t163 + t380;
t68 = (t312 * t378 + t314 * t377) * t313 + t329;
t67 = t378 * t293 + (-t245 + t358) * t292 + t373;
t66 = t135 * t315 + (-t114 * t314 + t115 * t312) * t313;
t65 = t114 * t292 + t115 * t293 + t135 * t302;
t62 = t114 * t271 + t115 * t273 + t135 * t290;
t55 = t123 * t315 + (-t110 * t314 + t111 * t312) * t313;
t54 = t122 * t315 + (-t108 * t314 + t109 * t312) * t313;
t53 = t110 * t292 + t111 * t293 + t123 * t302;
t52 = t108 * t292 + t109 * t293 + t122 * t302;
t49 = t110 * t271 + t111 * t273 + t123 * t290;
t48 = t108 * t271 + t109 * t273 + t122 * t290;
t1 = [m(3) + m(4) + m(5) + m(6) + m(7) + m(2); m(4) * t362 / 0.2e1 + t339 * t410 + (m(3) * t287 - t413) * t385 + (m(3) * t286 + t324) * t386 + t414 * (0.2e1 * t178 + 0.2e1 * t179 + t339); (t68 ^ 2 + t72 ^ 2 + t73 ^ 2) * m(7) + (t89 ^ 2 + t98 ^ 2 + t99 ^ 2) * m(6) + (t116 ^ 2 + t124 ^ 2 + t125 ^ 2) * m(5) + (t154 ^ 2 + t159 ^ 2 + t160 ^ 2) * m(4) + m(3) * (t246 ^ 2 + t262 ^ 2 + t263 ^ 2) + (t55 + t83 + (t281 * t386 + t283 * t305 + t285 * t306) * t386 + t418) * t386 + (-t54 - t82 + (-t280 * t385 + t282 * t303 + t284 * t304) * t385 + (-t280 * t386 + t281 * t385 - t282 * t305 - t283 * t303 - t284 * t306 - t285 * t304) * t386 - t419) * t385 + ((t298 * t386 + t299 * t305 + t300 * t306) * t386 - (-t298 * t385 + t299 * t303 + t300 * t304) * t385 + t66 + t90 + ((t283 * t322 + t285 * t319) * t312 - (t282 * t322 + t284 * t319) * t314) * t313 ^ 2 + ((-t280 * t314 + t281 * t312 + t299 * t322 + t300 * t319) * t313 + t315 * t298) * t315 + t415) * t315; t364 * t410 + t324 * t293 + t413 * t292 + t414 * (t182 * t411 + 0.2e1 * t171 + t364); (t154 * t157 + t159 * t175 + t160 * t174) * m(4) + (t116 * t117 + t124 * t127 + t125 * t126) * m(5) + (t67 * t68 + t70 * t73 + t71 * t72) * m(7) + (t86 * t89 + t96 * t99 + t97 * t98) * m(6) + (t61 / 0.2e1 + t31 / 0.2e1 + t65 / 0.2e1 + t88 / 0.2e1) * t315 + (t33 / 0.2e1 + t64 / 0.2e1 + t66 / 0.2e1 + t90 / 0.2e1) * t302 + (t55 / 0.2e1 + t83 / 0.2e1 + t51 / 0.2e1 + t22 / 0.2e1) * t293 + (t54 / 0.2e1 + t82 / 0.2e1 + t50 / 0.2e1 + t21 / 0.2e1) * t292 + ((-t44 / 0.2e1 - t52 / 0.2e1 - t80 / 0.2e1 - t17 / 0.2e1) * t314 + (t45 / 0.2e1 + t53 / 0.2e1 + t81 / 0.2e1 + t18 / 0.2e1) * t312) * t313; (t67 ^ 2 + t70 ^ 2 + t71 ^ 2) * m(7) + (t86 ^ 2 + t96 ^ 2 + t97 ^ 2) * m(6) + (t117 ^ 2 + t126 ^ 2 + t127 ^ 2) * m(5) + (t157 ^ 2 + t174 ^ 2 + t175 ^ 2) * m(4) + (t65 + t88 + t416) * t302 + (t53 + t81 + t420) * t293 + (t52 + t80 + t421) * t292; t130 * m(5) + (t374 + t376) * t409 + (t340 + t376) * t408; (t68 * t69 + t72 * t75 + t73 * t74) * m(7) + (t106 * t99 + t107 * t98 + t89 * t91) * m(6) + (t116 * t130 + t124 * t141 + t125 * t140) * m(5) + t66 * t402 + t62 * t396 + t55 * t403 + t54 * t404 + t325 + (t397 * t48 + t398 * t49) * t313; (t67 * t69 + t70 * t74 + t71 * t75) * m(7) + (t106 * t96 + t107 * t97 + t86 * t91) * m(6) + (t117 * t130 + t126 * t140 + t127 * t141) * m(5) + t326 + t49 * t400 + t48 * t401 + t65 * t402 + t62 * t399 + t53 * t403 + t52 * t404; (t69 ^ 2 + t74 ^ 2 + t75 ^ 2) * m(7) + (t106 ^ 2 + t107 ^ 2 + t91 ^ 2) * m(6) + t273 * t49 + t271 * t48 + (t130 ^ 2 + t140 ^ 2 + t141 ^ 2) * m(5) + t290 * t62 + t333; t340 * t408 + t374 * t409; (t128 * t89 + t138 * t99 + t139 * t98) * m(6) + (t68 * t87 + t72 * t95 + t73 * t94) * m(7) + t325; (t128 * t86 + t138 * t96 + t139 * t97) * m(6) + (t67 * t87 + t70 * t94 + t71 * t95) * m(7) + t326; (t69 * t87 + t74 * t94 + t75 * t95) * m(7) + (t106 * t138 + t107 * t139 + t128 * t91) * m(6) + t333; (t87 ^ 2 + t94 ^ 2 + t95 ^ 2) * m(7) + (t128 ^ 2 + t138 ^ 2 + t139 ^ 2) * m(6) + t333; t100 * m(7); (t100 * t68 + t118 * t73 + t119 * t72) * m(7) + t33 * t405 + t22 * t406 + t21 * t407 + t24 * t396 + (t3 * t397 + t398 * t4) * t313; (t100 * t67 + t118 * t70 + t119 * t71) * m(7) + t24 * t399 + t31 * t405 + t4 * t400 + t3 * t401 + t17 * t407 + t18 * t406; (t100 * t69 + t118 * t74 + t119 * t75) * m(7) + t334; (t100 * t87 + t118 * t94 + t119 * t95) * m(7) + t334; t251 * t4 + t249 * t3 + t265 * t24 + (t100 ^ 2 + t118 ^ 2 + t119 ^ 2) * m(7);];
%% Postprocessing: Reshape Output
% From vec2symmat_6_matlab.m
res = [t1(1) t1(2) t1(4) t1(7) t1(11) t1(16); t1(2) t1(3) t1(5) t1(8) t1(12) t1(17); t1(4) t1(5) t1(6) t1(9) t1(13) t1(18); t1(7) t1(8) t1(9) t1(10) t1(14) t1(19); t1(11) t1(12) t1(13) t1(14) t1(15) t1(20); t1(16) t1(17) t1(18) t1(19) t1(20) t1(21);];
Mq  = res;
