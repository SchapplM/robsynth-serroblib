% Calculate joint inertia matrix for
% S6RRRPRR12
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
% Datum: 2019-03-09 19:50
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Mq = S6RRRPRR12_inertiaJ_slag_vp1(qJ, ...
  pkin, m, rSges, Icges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(12,1),zeros(7,1),zeros(7,3),zeros(7,6)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRRPRR12_inertiaJ_slag_vp1: qJ has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [12 1]), ...
  'S6RRRPRR12_inertiaJ_slag_vp1: pkin has to be [12x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6RRRPRR12_inertiaJ_slag_vp1: m has to be [7x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [7,3]), ...
  'S6RRRPRR12_inertiaJ_slag_vp1: rSges has to be [7x3] (double)');
assert(isreal(Icges) && all(size(Icges) == [7 6]), ...
  'S6RRRPRR12_inertiaJ_slag_vp1: Icges has to be [7x6] (double)'); 

%% Symbolic Calculation
% From inertia_joint_joint_floatb_twist_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 19:38:14
% EndTime: 2019-03-09 19:38:29
% DurationCPUTime: 5.65s
% Computational Cost: add. (32139->693), mult. (58727->948), div. (0->0), fcn. (74786->14), ass. (0->324)
t321 = sin(pkin(6));
t329 = cos(qJ(1));
t383 = t321 * t329;
t323 = cos(pkin(6));
t325 = sin(qJ(3));
t326 = sin(qJ(2));
t394 = cos(qJ(3));
t351 = t321 * t394;
t294 = t323 * t325 + t326 * t351;
t319 = pkin(12) + qJ(5);
t313 = sin(t319);
t314 = cos(t319);
t328 = cos(qJ(2));
t384 = t321 * t328;
t262 = -t294 * t313 - t314 * t384;
t263 = t294 * t314 - t313 * t384;
t386 = t321 * t326;
t293 = -t323 * t394 + t325 * t386;
t189 = Icges(6,5) * t263 + Icges(6,6) * t262 + Icges(6,3) * t293;
t190 = Icges(6,4) * t263 + Icges(6,2) * t262 + Icges(6,6) * t293;
t191 = Icges(6,1) * t263 + Icges(6,4) * t262 + Icges(6,5) * t293;
t101 = t293 * t189 + t262 * t190 + t263 * t191;
t320 = sin(pkin(12));
t322 = cos(pkin(12));
t273 = -t294 * t320 - t322 * t384;
t359 = t320 * t384;
t274 = t294 * t322 - t359;
t195 = Icges(5,5) * t274 + Icges(5,6) * t273 + Icges(5,3) * t293;
t196 = Icges(5,4) * t274 + Icges(5,2) * t273 + Icges(5,6) * t293;
t197 = Icges(5,1) * t274 + Icges(5,4) * t273 + Icges(5,5) * t293;
t104 = t293 * t195 + t273 * t196 + t274 * t197;
t240 = Icges(4,5) * t294 - Icges(4,6) * t293 - Icges(4,3) * t384;
t241 = Icges(4,4) * t294 - Icges(4,2) * t293 - Icges(4,6) * t384;
t242 = Icges(4,1) * t294 - Icges(4,4) * t293 - Icges(4,5) * t384;
t131 = -t240 * t384 - t293 * t241 + t294 * t242;
t315 = qJ(6) + t319;
t310 = sin(t315);
t311 = cos(t315);
t257 = -t294 * t310 - t311 * t384;
t258 = t294 * t311 - t310 * t384;
t185 = Icges(7,5) * t258 + Icges(7,6) * t257 + Icges(7,3) * t293;
t186 = Icges(7,4) * t258 + Icges(7,2) * t257 + Icges(7,6) * t293;
t187 = Icges(7,1) * t258 + Icges(7,4) * t257 + Icges(7,5) * t293;
t95 = t293 * t185 + t257 * t186 + t258 * t187;
t401 = -t101 - t104 - t131 - t95;
t327 = sin(qJ(1));
t380 = t327 * t328;
t381 = t326 * t329;
t296 = t323 * t381 + t380;
t275 = t296 * t325 + t329 * t351;
t400 = t275 / 0.2e1;
t379 = t328 * t329;
t382 = t326 * t327;
t298 = -t323 * t382 + t379;
t277 = t298 * t325 - t327 * t351;
t399 = t277 / 0.2e1;
t398 = t293 / 0.2e1;
t295 = -t323 * t379 + t382;
t397 = t295 / 0.2e1;
t297 = t323 * t380 + t381;
t396 = t297 / 0.2e1;
t395 = t323 / 0.2e1;
t393 = pkin(4) * t320;
t312 = t322 * pkin(4) + pkin(3);
t392 = -pkin(3) + t312;
t324 = -pkin(10) - qJ(4);
t276 = t296 * t394 - t325 * t383;
t221 = -t276 * t310 + t295 * t311;
t222 = t276 * t311 + t295 * t310;
t138 = Icges(7,5) * t222 + Icges(7,6) * t221 + Icges(7,3) * t275;
t140 = Icges(7,4) * t222 + Icges(7,2) * t221 + Icges(7,6) * t275;
t142 = Icges(7,1) * t222 + Icges(7,4) * t221 + Icges(7,5) * t275;
t75 = t138 * t293 + t140 * t257 + t142 * t258;
t391 = t75 * t275;
t385 = t321 * t327;
t278 = t298 * t394 + t325 * t385;
t223 = -t278 * t310 + t297 * t311;
t224 = t278 * t311 + t297 * t310;
t139 = Icges(7,5) * t224 + Icges(7,6) * t223 + Icges(7,3) * t277;
t141 = Icges(7,4) * t224 + Icges(7,2) * t223 + Icges(7,6) * t277;
t143 = Icges(7,1) * t224 + Icges(7,4) * t223 + Icges(7,5) * t277;
t76 = t139 * t293 + t141 * t257 + t143 * t258;
t390 = t76 * t277;
t244 = Icges(3,5) * t296 - Icges(3,6) * t295 - Icges(3,3) * t383;
t389 = t244 * t329;
t388 = t295 * t320;
t387 = t297 * t320;
t301 = pkin(5) * t313 + t393;
t348 = -t301 + t393;
t318 = -pkin(11) + t324;
t364 = -t318 + t324;
t300 = pkin(5) * t314 + t312;
t366 = t300 - t312;
t126 = t275 * t364 + t276 * t366 - t295 * t348;
t337 = -rSges(7,1) * t222 - rSges(7,2) * t221;
t144 = rSges(7,3) * t275 - t337;
t378 = t126 + t144;
t353 = -pkin(4) * t387 + t277 * t324 - t278 * t312;
t354 = -t277 * t318 + t278 * t300 + t297 * t301;
t127 = t353 + t354;
t145 = t224 * rSges(7,1) + t223 * rSges(7,2) + t277 * rSges(7,3);
t377 = t127 + t145;
t266 = t275 * qJ(4);
t362 = pkin(4) * t388;
t146 = -t275 * t324 + t276 * t392 - t266 + t362;
t215 = pkin(3) * t276 + t266;
t207 = t297 * t215;
t376 = t297 * t146 + t207;
t216 = t278 * pkin(3) + qJ(4) * t277;
t147 = -t216 - t353;
t375 = -t147 - t216;
t238 = -t278 * t320 + t297 * t322;
t239 = t278 * t322 + t387;
t167 = t239 * rSges(5,1) + t238 * rSges(5,2) + t277 * rSges(5,3);
t374 = -t167 - t216;
t175 = t293 * t364 + t294 * t366 + t348 * t384;
t188 = rSges(7,1) * t258 + rSges(7,2) * t257 + rSges(7,3) * t293;
t373 = t175 + t188;
t193 = -pkin(4) * t359 + t392 * t294 + (-qJ(4) - t324) * t293;
t256 = t294 * pkin(3) + t293 * qJ(4);
t372 = -t193 - t256;
t198 = rSges(5,1) * t274 + rSges(5,2) * t273 + rSges(5,3) * t293;
t371 = -t198 - t256;
t370 = t215 * t384 + t295 * t256;
t261 = t298 * pkin(2) + pkin(9) * t297;
t255 = t323 * t261;
t369 = t323 * t216 + t255;
t260 = pkin(2) * t296 + t295 * pkin(9);
t368 = -t215 - t260;
t367 = t260 * t385 + t261 * t383;
t365 = t329 * pkin(1) + pkin(8) * t385;
t91 = t95 * t293;
t26 = t390 + t91 + t391;
t61 = t138 * t275 + t140 * t221 + t142 * t222;
t62 = t139 * t275 + t141 * t221 + t143 * t222;
t85 = t185 * t275 + t186 * t221 + t187 * t222;
t7 = t275 * t61 + t277 * t62 + t293 * t85;
t63 = t138 * t277 + t140 * t223 + t142 * t224;
t64 = t139 * t277 + t141 * t223 + t143 * t224;
t86 = t185 * t277 + t186 * t223 + t187 * t224;
t8 = t275 * t63 + t277 * t64 + t293 * t86;
t363 = t293 * t26 + t275 * t7 + t277 * t8;
t227 = -t276 * t313 + t295 * t314;
t228 = t276 * t314 + t295 * t313;
t148 = Icges(6,5) * t228 + Icges(6,6) * t227 + Icges(6,3) * t275;
t150 = Icges(6,4) * t228 + Icges(6,2) * t227 + Icges(6,6) * t275;
t152 = Icges(6,1) * t228 + Icges(6,4) * t227 + Icges(6,5) * t275;
t77 = t148 * t293 + t150 * t262 + t152 * t263;
t87 = t189 * t275 + t190 * t227 + t191 * t228;
t361 = t87 / 0.2e1 + t77 / 0.2e1;
t229 = -t278 * t313 + t297 * t314;
t230 = t278 * t314 + t297 * t313;
t149 = Icges(6,5) * t230 + Icges(6,6) * t229 + Icges(6,3) * t277;
t151 = Icges(6,4) * t230 + Icges(6,2) * t229 + Icges(6,6) * t277;
t153 = Icges(6,1) * t230 + Icges(6,4) * t229 + Icges(6,5) * t277;
t78 = t149 * t293 + t151 * t262 + t153 * t263;
t88 = t189 * t277 + t190 * t229 + t191 * t230;
t360 = t88 / 0.2e1 + t78 / 0.2e1;
t358 = t323 * t147 + t369;
t357 = -t146 + t368;
t157 = t230 * rSges(6,1) + t229 * rSges(6,2) + t277 * rSges(6,3);
t356 = -t157 + t375;
t192 = rSges(6,1) * t263 + rSges(6,2) * t262 + rSges(6,3) * t293;
t355 = -t192 + t372;
t206 = t278 * rSges(4,1) - t277 * rSges(4,2) + t297 * rSges(4,3);
t283 = Icges(3,3) * t323 + (Icges(3,5) * t326 + Icges(3,6) * t328) * t321;
t284 = Icges(3,6) * t323 + (Icges(3,4) * t326 + Icges(3,2) * t328) * t321;
t285 = Icges(3,5) * t323 + (Icges(3,1) * t326 + Icges(3,4) * t328) * t321;
t352 = t323 * t283 + t284 * t384 + t285 * t386;
t251 = t298 * rSges(3,1) - t297 * rSges(3,2) + rSges(3,3) * t385;
t350 = -t384 / 0.2e1;
t349 = -t327 * pkin(1) + pkin(8) * t383;
t243 = rSges(4,1) * t294 - rSges(4,2) * t293 - rSges(4,3) * t384;
t299 = (pkin(2) * t326 - pkin(9) * t328) * t321;
t347 = t321 * (-t243 - t299);
t346 = t375 - t377;
t345 = t146 * t384 + t295 * t193 + t370;
t344 = t372 - t373;
t343 = t215 * t385 + t216 * t383 + t367;
t342 = t391 / 0.2e1 + t390 / 0.2e1 + t85 * t400 + t86 * t399 + t91;
t341 = t321 * (-t299 + t371);
t13 = t295 * t61 + t297 * t62 - t384 * t85;
t14 = t295 * t63 + t297 * t64 - t384 * t86;
t30 = t75 * t295 + t76 * t297 - t384 * t95;
t340 = t13 * t400 + t14 * t399 + t26 * t350 + t30 * t398 + t8 * t396 + t7 * t397;
t15 = t85 * t323 + (t327 * t62 - t329 * t61) * t321;
t16 = t86 * t323 + (t327 * t64 - t329 * t63) * t321;
t94 = t95 * t323;
t34 = t94 + (t76 * t327 - t75 * t329) * t321;
t339 = t15 * t400 + t16 * t399 + t26 * t395 + t34 * t398 + t8 * t385 / 0.2e1 - t7 * t383 / 0.2e1;
t338 = -rSges(6,1) * t228 - rSges(6,2) * t227;
t336 = t261 + t365;
t335 = t321 * (-t299 + t355);
t334 = t146 * t385 + t147 * t383 + t343;
t333 = t321 * (-t299 + t344);
t332 = -t260 + t349;
t205 = rSges(4,1) * t276 - rSges(4,2) * t275 + rSges(4,3) * t295;
t236 = -t276 * t320 + t295 * t322;
t237 = t276 * t322 + t388;
t166 = rSges(5,1) * t237 + rSges(5,2) * t236 + rSges(5,3) * t275;
t250 = t296 * rSges(3,1) - t295 * rSges(3,2) - rSges(3,3) * t383;
t199 = Icges(4,5) * t276 - Icges(4,6) * t275 + Icges(4,3) * t295;
t201 = Icges(4,4) * t276 - Icges(4,2) * t275 + Icges(4,6) * t295;
t203 = Icges(4,1) * t276 - Icges(4,4) * t275 + Icges(4,5) * t295;
t117 = -t199 * t384 - t201 * t293 + t203 * t294;
t122 = t240 * t295 - t241 * t275 + t242 * t276;
t160 = Icges(5,5) * t237 + Icges(5,6) * t236 + Icges(5,3) * t275;
t162 = Icges(5,4) * t237 + Icges(5,2) * t236 + Icges(5,6) * t275;
t164 = Icges(5,1) * t237 + Icges(5,4) * t236 + Icges(5,5) * t275;
t80 = t160 * t293 + t162 * t273 + t164 * t274;
t92 = t195 * t275 + t196 * t236 + t197 * t237;
t331 = t75 / 0.2e1 + t85 / 0.2e1 + t92 / 0.2e1 + t122 / 0.2e1 + t117 / 0.2e1 + t80 / 0.2e1 + t361;
t200 = Icges(4,5) * t278 - Icges(4,6) * t277 + Icges(4,3) * t297;
t202 = Icges(4,4) * t278 - Icges(4,2) * t277 + Icges(4,6) * t297;
t204 = Icges(4,1) * t278 - Icges(4,4) * t277 + Icges(4,5) * t297;
t118 = -t200 * t384 - t202 * t293 + t204 * t294;
t123 = t240 * t297 - t241 * t277 + t242 * t278;
t161 = Icges(5,5) * t239 + Icges(5,6) * t238 + Icges(5,3) * t277;
t163 = Icges(5,4) * t239 + Icges(5,2) * t238 + Icges(5,6) * t277;
t165 = Icges(5,1) * t239 + Icges(5,4) * t238 + Icges(5,5) * t277;
t81 = t161 * t293 + t163 * t273 + t165 * t274;
t93 = t195 * t277 + t196 * t238 + t197 * t239;
t330 = t76 / 0.2e1 + t86 / 0.2e1 + t123 / 0.2e1 + t93 / 0.2e1 + t118 / 0.2e1 + t81 / 0.2e1 + t360;
t303 = rSges(2,1) * t329 - t327 * rSges(2,2);
t302 = -t327 * rSges(2,1) - rSges(2,2) * t329;
t286 = t323 * rSges(3,3) + (rSges(3,1) * t326 + rSges(3,2) * t328) * t321;
t249 = Icges(3,1) * t298 - Icges(3,4) * t297 + Icges(3,5) * t385;
t248 = Icges(3,1) * t296 - Icges(3,4) * t295 - Icges(3,5) * t383;
t247 = Icges(3,4) * t298 - Icges(3,2) * t297 + Icges(3,6) * t385;
t246 = Icges(3,4) * t296 - Icges(3,2) * t295 - Icges(3,6) * t383;
t245 = Icges(3,5) * t298 - Icges(3,6) * t297 + Icges(3,3) * t385;
t232 = t251 + t365;
t231 = -t250 + t349;
t210 = -t323 * t250 - t286 * t383;
t209 = t251 * t323 - t286 * t385;
t194 = t352 * t323;
t183 = (t250 * t327 + t251 * t329) * t321;
t182 = t283 * t385 - t284 * t297 + t285 * t298;
t181 = -t283 * t383 - t295 * t284 + t296 * t285;
t174 = t275 * t188;
t171 = t336 + t206;
t170 = -t205 + t332;
t159 = -t206 * t384 - t243 * t297;
t158 = t205 * t384 + t243 * t295;
t156 = rSges(6,3) * t275 - t338;
t155 = t323 * t245 + (t247 * t328 + t249 * t326) * t321;
t154 = t323 * t244 + (t246 * t328 + t248 * t326) * t321;
t132 = t293 * t145;
t130 = t131 * t323;
t129 = t277 * t144;
t128 = t205 * t297 - t206 * t295;
t125 = (-t205 - t260) * t323 + t329 * t347;
t124 = t323 * t206 + t327 * t347 + t255;
t121 = t336 - t374;
t120 = -t166 - t215 + t332;
t119 = (t205 * t327 + t206 * t329) * t321 + t367;
t116 = t336 - t353 + t157;
t115 = -t362 - t276 * t312 + (-rSges(6,3) + t324) * t275 + t332 + t338;
t114 = t157 * t293 - t192 * t277;
t113 = -t156 * t293 + t192 * t275;
t112 = -t188 * t277 + t132;
t111 = -t144 * t293 + t174;
t110 = t336 + t354 + t145;
t109 = -t276 * t300 - t295 * t301 + (-rSges(7,3) + t318) * t275 + t332 + t337;
t108 = t200 * t297 - t202 * t277 + t204 * t278;
t107 = t199 * t297 - t201 * t277 + t203 * t278;
t106 = t200 * t295 - t202 * t275 + t204 * t276;
t105 = t199 * t295 - t201 * t275 + t203 * t276;
t103 = t104 * t323;
t102 = t156 * t277 - t157 * t275;
t100 = t101 * t323;
t99 = -t145 * t275 + t129;
t98 = t297 * t371 + t374 * t384;
t97 = t166 * t384 + t198 * t295 + t370;
t96 = t101 * t293;
t90 = (-t166 + t368) * t323 + t329 * t341;
t89 = t323 * t167 + t327 * t341 + t369;
t82 = t297 * t166 + t295 * t374 + t207;
t79 = (t166 * t327 + t167 * t329) * t321 + t343;
t74 = t161 * t277 + t163 * t238 + t165 * t239;
t73 = t160 * t277 + t162 * t238 + t164 * t239;
t72 = t161 * t275 + t163 * t236 + t165 * t237;
t71 = t160 * t275 + t162 * t236 + t164 * t237;
t68 = t149 * t277 + t151 * t229 + t153 * t230;
t67 = t148 * t277 + t150 * t229 + t152 * t230;
t66 = t149 * t275 + t151 * t227 + t153 * t228;
t65 = t148 * t275 + t150 * t227 + t152 * t228;
t60 = t297 * t355 + t356 * t384;
t59 = t156 * t384 + t192 * t295 + t345;
t58 = t293 * t127 - t277 * t373 + t132;
t57 = t275 * t175 - t293 * t378 + t174;
t56 = (-t156 + t357) * t323 + t329 * t335;
t55 = t157 * t323 + t327 * t335 + t358;
t54 = t130 + (-t117 * t329 + t118 * t327) * t321;
t53 = t117 * t295 + t118 * t297 - t131 * t384;
t52 = t277 * t126 - t275 * t377 + t129;
t51 = t297 * t156 + t295 * t356 + t376;
t50 = (t156 * t327 + t157 * t329) * t321 + t334;
t49 = t123 * t323 + (-t107 * t329 + t108 * t327) * t321;
t48 = t122 * t323 + (-t105 * t329 + t106 * t327) * t321;
t47 = t107 * t295 + t108 * t297 - t123 * t384;
t46 = t105 * t295 + t106 * t297 - t122 * t384;
t45 = t297 * t344 + t346 * t384;
t44 = t295 * t373 + t378 * t384 + t345;
t43 = (t357 - t378) * t323 + t329 * t333;
t42 = t323 * t377 + t327 * t333 + t358;
t41 = t295 * t346 + t297 * t378 + t376;
t40 = (t327 * t378 + t329 * t377) * t321 + t334;
t39 = t103 + (t81 * t327 - t80 * t329) * t321;
t38 = -t104 * t384 + t80 * t295 + t81 * t297;
t37 = t100 + (t78 * t327 - t77 * t329) * t321;
t36 = -t101 * t384 + t77 * t295 + t78 * t297;
t35 = t77 * t275 + t78 * t277 + t96;
t32 = t93 * t323 + (t327 * t74 - t329 * t73) * t321;
t31 = t92 * t323 + (t327 * t72 - t329 * t71) * t321;
t28 = t295 * t73 + t297 * t74 - t384 * t93;
t27 = t295 * t71 + t297 * t72 - t384 * t92;
t22 = t88 * t323 + (t327 * t68 - t329 * t67) * t321;
t21 = t87 * t323 + (t327 * t66 - t329 * t65) * t321;
t20 = t295 * t67 + t297 * t68 - t384 * t88;
t19 = t295 * t65 + t297 * t66 - t384 * t87;
t18 = t275 * t67 + t277 * t68 + t293 * t88;
t17 = t275 * t65 + t277 * t66 + t293 * t87;
t1 = [m(7) * (t109 ^ 2 + t110 ^ 2) + m(6) * (t115 ^ 2 + t116 ^ 2) + m(5) * (t120 ^ 2 + t121 ^ 2) + m(4) * (t170 ^ 2 + t171 ^ 2) + m(3) * (t231 ^ 2 + t232 ^ 2) + m(2) * (t302 ^ 2 + t303 ^ 2) + Icges(2,3) + t352 - t401; t130 + t103 + t100 + t94 + t194 + m(3) * (t209 * t232 + t210 * t231) + m(4) * (t124 * t171 + t125 * t170) + m(7) * (t109 * t43 + t110 * t42) + m(6) * (t115 * t56 + t116 * t55) + m(5) * (t120 * t90 + t121 * t89) + ((-t154 / 0.2e1 - t181 / 0.2e1 - t331) * t329 + (t182 / 0.2e1 + t155 / 0.2e1 + t330) * t327) * t321; (t34 + t37 + t39 + t54 + t194) * t323 + m(7) * (t40 ^ 2 + t42 ^ 2 + t43 ^ 2) + m(6) * (t50 ^ 2 + t55 ^ 2 + t56 ^ 2) + m(5) * (t79 ^ 2 + t89 ^ 2 + t90 ^ 2) + m(4) * (t119 ^ 2 + t124 ^ 2 + t125 ^ 2) + m(3) * (t183 ^ 2 + t209 ^ 2 + t210 ^ 2) + ((-t15 - t21 - t48 - t31 + (-t295 * t246 + t296 * t248 - t321 * t389) * t383) * t329 + (t16 + t22 + t49 + t32 + ((-t247 * t297 + t249 * t298 + (t245 * t327 - t389) * t321) * t327 + (t245 * t383 + t246 * t297 + t295 * t247 - t248 * t298 - t296 * t249) * t329) * t321) * t327 + ((-t154 - t181) * t329 + (t155 + t182) * t327) * t323) * t321; t401 * t384 + m(7) * (t109 * t44 + t110 * t45) + m(6) * (t115 * t59 + t116 * t60) + m(5) * (t120 * t97 + t121 * t98) + m(4) * (t158 * t170 + t159 * t171) + t330 * t297 + t331 * t295; (t53 / 0.2e1 + t30 / 0.2e1 + t36 / 0.2e1 + t38 / 0.2e1) * t323 + (t49 / 0.2e1 + t16 / 0.2e1 + t22 / 0.2e1 + t32 / 0.2e1) * t297 + (t48 / 0.2e1 + t15 / 0.2e1 + t21 / 0.2e1 + t31 / 0.2e1) * t295 + m(4) * (t119 * t128 + t124 * t159 + t125 * t158) + m(7) * (t40 * t41 + t42 * t45 + t43 * t44) + m(6) * (t50 * t51 + t55 * t60 + t56 * t59) + m(5) * (t79 * t82 + t89 * t98 + t90 * t97) + ((-t13 / 0.2e1 - t19 / 0.2e1 - t27 / 0.2e1 - t46 / 0.2e1) * t329 + (-t34 / 0.2e1 - t37 / 0.2e1 - t39 / 0.2e1 - t54 / 0.2e1) * t328 + (t14 / 0.2e1 + t20 / 0.2e1 + t28 / 0.2e1 + t47 / 0.2e1) * t327) * t321; (-t30 - t36 - t38 - t53) * t384 + (t14 + t20 + t47 + t28) * t297 + (t13 + t19 + t46 + t27) * t295 + m(7) * (t41 ^ 2 + t44 ^ 2 + t45 ^ 2) + m(6) * (t51 ^ 2 + t59 ^ 2 + t60 ^ 2) + m(5) * (t82 ^ 2 + t97 ^ 2 + t98 ^ 2) + m(4) * (t128 ^ 2 + t158 ^ 2 + t159 ^ 2); m(7) * (t109 * t277 + t110 * t275) + m(6) * (t115 * t277 + t116 * t275) + m(5) * (t120 * t277 + t121 * t275); m(7) * (t275 * t42 + t277 * t43 + t293 * t40) + m(6) * (t275 * t55 + t277 * t56 + t293 * t50) + m(5) * (t275 * t89 + t277 * t90 + t293 * t79); m(7) * (t275 * t45 + t277 * t44 + t293 * t41) + m(6) * (t275 * t60 + t277 * t59 + t293 * t51) + m(5) * (t275 * t98 + t277 * t97 + t293 * t82); 0.2e1 * (m(5) / 0.2e1 + m(6) / 0.2e1 + m(7) / 0.2e1) * (t275 ^ 2 + t277 ^ 2 + t293 ^ 2); t96 + t360 * t277 + t361 * t275 + m(7) * (t109 * t57 + t110 * t58) + m(6) * (t113 * t115 + t114 * t116) + t342; t22 * t399 + t21 * t400 + t37 * t398 + t35 * t395 + (t327 * t18 / 0.2e1 - t329 * t17 / 0.2e1) * t321 + m(7) * (t40 * t52 + t42 * t58 + t43 * t57) + m(6) * (t102 * t50 + t113 * t56 + t114 * t55) + t339; t35 * t350 + t18 * t396 + t20 * t399 + t19 * t400 + t17 * t397 + t36 * t398 + m(7) * (t41 * t52 + t44 * t57 + t45 * t58) + m(6) * (t102 * t51 + t113 * t59 + t114 * t60) + t340; m(6) * (t102 * t293 + t113 * t277 + t114 * t275) + m(7) * (t275 * t58 + t277 * t57 + t293 * t52); t275 * t17 + t277 * t18 + t293 * t35 + m(7) * (t52 ^ 2 + t57 ^ 2 + t58 ^ 2) + m(6) * (t102 ^ 2 + t113 ^ 2 + t114 ^ 2) + t363; m(7) * (t109 * t111 + t110 * t112) + t342; m(7) * (t111 * t43 + t112 * t42 + t40 * t99) + t339; m(7) * (t111 * t44 + t112 * t45 + t41 * t99) + t340; m(7) * (t111 * t277 + t112 * t275 + t293 * t99); m(7) * (t111 * t57 + t112 * t58 + t52 * t99) + t363; m(7) * (t111 ^ 2 + t112 ^ 2 + t99 ^ 2) + t363;];
%% Postprocessing: Reshape Output
% From vec2symmat_6_matlab.m
res = [t1(1) t1(2) t1(4) t1(7) t1(11) t1(16); t1(2) t1(3) t1(5) t1(8) t1(12) t1(17); t1(4) t1(5) t1(6) t1(9) t1(13) t1(18); t1(7) t1(8) t1(9) t1(10) t1(14) t1(19); t1(11) t1(12) t1(13) t1(14) t1(15) t1(20); t1(16) t1(17) t1(18) t1(19) t1(20) t1(21);];
Mq  = res;
