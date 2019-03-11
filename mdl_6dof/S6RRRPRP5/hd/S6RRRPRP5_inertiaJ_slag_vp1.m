% Calculate joint inertia matrix for
% S6RRRPRP5
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% pkin [10x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d2,d3,d5,theta4]';
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
% Datum: 2019-03-09 16:52
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Mq = S6RRRPRP5_inertiaJ_slag_vp1(qJ, ...
  pkin, m, rSges, Icges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(10,1),zeros(7,1),zeros(7,3),zeros(7,6)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRRPRP5_inertiaJ_slag_vp1: qJ has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S6RRRPRP5_inertiaJ_slag_vp1: pkin has to be [10x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6RRRPRP5_inertiaJ_slag_vp1: m has to be [7x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [7,3]), ...
  'S6RRRPRP5_inertiaJ_slag_vp1: rSges has to be [7x3] (double)');
assert(isreal(Icges) && all(size(Icges) == [7 6]), ...
  'S6RRRPRP5_inertiaJ_slag_vp1: Icges has to be [7x6] (double)'); 

%% Symbolic Calculation
% From inertia_joint_joint_floatb_twist_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 16:47:47
% EndTime: 2019-03-09 16:48:02
% DurationCPUTime: 6.81s
% Computational Cost: add. (13641->533), mult. (14842->749), div. (0->0), fcn. (16044->10), ass. (0->267)
t267 = qJ(3) + pkin(10);
t260 = qJ(5) + t267;
t255 = sin(t260);
t256 = cos(t260);
t277 = cos(qJ(1));
t274 = sin(qJ(1));
t276 = cos(qJ(2));
t333 = t274 * t276;
t199 = t255 * t333 + t256 * t277;
t200 = -t255 * t277 + t256 * t333;
t273 = sin(qJ(2));
t336 = t273 * t274;
t123 = Icges(7,5) * t200 + Icges(7,6) * t336 + Icges(7,3) * t199;
t129 = Icges(6,4) * t200 - Icges(6,2) * t199 + Icges(6,6) * t336;
t391 = t123 - t129;
t332 = t276 * t277;
t201 = t255 * t332 - t274 * t256;
t202 = t274 * t255 + t256 * t332;
t335 = t273 * t277;
t124 = Icges(7,5) * t202 + Icges(7,6) * t335 + Icges(7,3) * t201;
t130 = Icges(6,4) * t202 - Icges(6,2) * t201 + Icges(6,6) * t335;
t390 = t124 - t130;
t125 = Icges(6,5) * t200 - Icges(6,6) * t199 + Icges(6,3) * t336;
t127 = Icges(7,4) * t200 + Icges(7,2) * t336 + Icges(7,6) * t199;
t389 = t125 + t127;
t126 = Icges(6,5) * t202 - Icges(6,6) * t201 + Icges(6,3) * t335;
t128 = Icges(7,4) * t202 + Icges(7,2) * t335 + Icges(7,6) * t201;
t388 = t126 + t128;
t131 = Icges(7,1) * t200 + Icges(7,4) * t336 + Icges(7,5) * t199;
t133 = Icges(6,1) * t200 - Icges(6,4) * t199 + Icges(6,5) * t336;
t387 = t131 + t133;
t132 = Icges(7,1) * t202 + Icges(7,4) * t335 + Icges(7,5) * t201;
t134 = Icges(6,1) * t202 - Icges(6,4) * t201 + Icges(6,5) * t335;
t386 = t132 + t134;
t365 = rSges(7,3) + qJ(6);
t373 = rSges(7,1) + pkin(5);
t385 = -t365 * t199 - t373 * t200;
t384 = t199 * t391 + t387 * t200 + t389 * t336;
t383 = t199 * t390 + t200 * t386 + t336 * t388;
t382 = t201 * t391 + t387 * t202 + t389 * t335;
t381 = t201 * t390 + t202 * t386 + t335 * t388;
t180 = -Icges(7,6) * t276 + (Icges(7,5) * t256 + Icges(7,3) * t255) * t273;
t182 = -Icges(7,2) * t276 + (Icges(7,4) * t256 + Icges(7,6) * t255) * t273;
t184 = -Icges(7,4) * t276 + (Icges(7,1) * t256 + Icges(7,5) * t255) * t273;
t84 = t180 * t199 + t182 * t336 + t184 * t200;
t181 = -Icges(6,3) * t276 + (Icges(6,5) * t256 - Icges(6,6) * t255) * t273;
t183 = -Icges(6,6) * t276 + (Icges(6,4) * t256 - Icges(6,2) * t255) * t273;
t185 = -Icges(6,5) * t276 + (Icges(6,1) * t256 - Icges(6,4) * t255) * t273;
t85 = t181 * t336 - t183 * t199 + t185 * t200;
t380 = -t85 - t84;
t86 = t201 * t180 + t182 * t335 + t202 * t184;
t87 = t181 * t335 - t201 * t183 + t202 * t185;
t379 = -t86 - t87;
t339 = t255 * t273;
t375 = t180 * t339 + (t184 + t185) * t256 * t273;
t377 = -t181 - t182;
t366 = -t183 * t339 + t276 * t377 + t375;
t378 = t366 * t276;
t376 = Icges(3,5) * t273;
t374 = t376 / 0.2e1;
t372 = t380 * t276 + (t274 * t384 + t277 * t383) * t273;
t371 = t379 * t276 + (t274 * t382 + t277 * t381) * t273;
t370 = t274 * t383 - t277 * t384;
t369 = t274 * t381 - t277 * t382;
t60 = -t276 * t127 + (t123 * t255 + t131 * t256) * t273;
t62 = -t276 * t125 + (-t129 * t255 + t133 * t256) * t273;
t368 = -t60 - t62;
t61 = -t276 * t128 + (t124 * t255 + t132 * t256) * t273;
t63 = -t276 * t126 + (-t130 * t255 + t134 * t256) * t273;
t367 = t61 + t63;
t364 = rSges(7,2) * t336 - t385;
t363 = -t276 * rSges(7,2) + (t365 * t255 + t256 * t373) * t273;
t258 = sin(t267);
t259 = cos(t267);
t190 = -Icges(5,3) * t276 + (Icges(5,5) * t259 - Icges(5,6) * t258) * t273;
t272 = sin(qJ(3));
t275 = cos(qJ(3));
t210 = -Icges(4,3) * t276 + (Icges(4,5) * t275 - Icges(4,6) * t272) * t273;
t362 = -t190 - t210;
t191 = -Icges(5,6) * t276 + (Icges(5,4) * t259 - Icges(5,2) * t258) * t273;
t213 = -Icges(4,6) * t276 + (Icges(4,4) * t275 - Icges(4,2) * t272) * t273;
t361 = -t191 * t258 - t213 * t272;
t206 = -t258 * t333 - t259 * t277;
t207 = -t258 * t277 + t259 * t333;
t141 = Icges(5,5) * t207 + Icges(5,6) * t206 + Icges(5,3) * t336;
t143 = Icges(5,4) * t207 + Icges(5,2) * t206 + Icges(5,6) * t336;
t145 = Icges(5,1) * t207 + Icges(5,4) * t206 + Icges(5,5) * t336;
t52 = t141 * t336 + t143 * t206 + t145 * t207;
t208 = -t258 * t332 + t274 * t259;
t209 = t274 * t258 + t259 * t332;
t142 = Icges(5,5) * t209 + Icges(5,6) * t208 + Icges(5,3) * t335;
t144 = Icges(5,4) * t209 + Icges(5,2) * t208 + Icges(5,6) * t335;
t146 = Icges(5,1) * t209 + Icges(5,4) * t208 + Icges(5,5) * t335;
t53 = t142 * t336 + t144 * t206 + t146 * t207;
t226 = -t272 * t333 - t275 * t277;
t337 = t272 * t277;
t227 = t275 * t333 - t337;
t157 = Icges(4,5) * t227 + Icges(4,6) * t226 + Icges(4,3) * t336;
t159 = Icges(4,4) * t227 + Icges(4,2) * t226 + Icges(4,6) * t336;
t161 = Icges(4,1) * t227 + Icges(4,4) * t226 + Icges(4,5) * t336;
t68 = t157 * t336 + t159 * t226 + t161 * t227;
t228 = -t272 * t332 + t274 * t275;
t334 = t274 * t272;
t229 = t275 * t332 + t334;
t158 = Icges(4,5) * t229 + Icges(4,6) * t228 + Icges(4,3) * t335;
t160 = Icges(4,4) * t229 + Icges(4,2) * t228 + Icges(4,6) * t335;
t162 = Icges(4,1) * t229 + Icges(4,4) * t228 + Icges(4,5) * t335;
t69 = t158 * t336 + t160 * t226 + t162 * t227;
t192 = -Icges(5,5) * t276 + (Icges(5,1) * t259 - Icges(5,4) * t258) * t273;
t90 = t190 * t336 + t191 * t206 + t192 * t207;
t216 = -Icges(4,5) * t276 + (Icges(4,1) * t275 - Icges(4,4) * t272) * t273;
t98 = t210 * t336 + t213 * t226 + t216 * t227;
t360 = (-t90 - t98) * t276 + ((t53 + t69) * t277 + (t52 + t68) * t274) * t273;
t54 = t141 * t335 + t208 * t143 + t209 * t145;
t55 = t142 * t335 + t208 * t144 + t209 * t146;
t70 = t157 * t335 + t228 * t159 + t229 * t161;
t71 = t158 * t335 + t228 * t160 + t229 * t162;
t91 = t190 * t335 + t208 * t191 + t209 * t192;
t99 = t210 * t335 + t228 * t213 + t229 * t216;
t359 = (-t91 - t99) * t276 + ((t55 + t71) * t277 + (t54 + t70) * t274) * t273;
t64 = -t276 * t141 + (-t143 * t258 + t145 * t259) * t273;
t76 = -t276 * t157 + (-t159 * t272 + t161 * t275) * t273;
t358 = -t64 - t76;
t65 = -t276 * t142 + (-t144 * t258 + t146 * t259) * t273;
t77 = -t276 * t158 + (-t160 * t272 + t162 * t275) * t273;
t357 = t65 + t77;
t356 = (t192 * t259 + t216 * t275) * t273;
t269 = t274 ^ 2;
t270 = t277 ^ 2;
t355 = m(5) / 0.2e1;
t354 = m(6) / 0.2e1;
t353 = m(7) / 0.2e1;
t352 = t274 / 0.2e1;
t351 = -t276 / 0.2e1;
t350 = -t277 / 0.2e1;
t349 = pkin(2) * t276;
t348 = pkin(8) * t273;
t257 = t275 * pkin(3) + pkin(2);
t347 = -pkin(2) + t257;
t271 = -qJ(4) - pkin(8);
t346 = t378 + (t274 * t368 - t277 * t367) * t273;
t344 = t277 * rSges(3,3);
t342 = Icges(3,4) * t276;
t331 = t273 * t361 + t276 * t362 + t356;
t314 = pkin(3) * t337 + t271 * t336;
t232 = pkin(4) * t259 + t257;
t316 = t232 - t257;
t234 = pkin(3) * t272 + pkin(4) * t258;
t266 = -pkin(9) + t271;
t317 = -t277 * t234 - t266 * t336;
t113 = t316 * t333 + t314 + t317;
t163 = (t276 * t347 - t348) * t274 - t314;
t151 = t163 * t335;
t330 = t113 * t335 + t151;
t311 = t266 - t271;
t315 = -pkin(3) * t334 - t257 * t332;
t319 = t232 * t332 + t274 * t234;
t114 = -t311 * t335 + t315 + t319;
t285 = -t271 * t335 - t315;
t313 = pkin(2) * t332 + pkin(8) * t335;
t164 = t285 - t313;
t329 = -t114 - t164;
t328 = t364 * t335;
t327 = rSges(7,2) * t335 + t365 * t201 + t202 * t373;
t148 = t209 * rSges(5,1) + t208 * rSges(5,2) + rSges(5,3) * t335;
t326 = -t148 - t164;
t189 = (pkin(8) + t271) * t276 + t347 * t273;
t325 = t276 * t163 + t189 * t336;
t167 = t273 * t316 + t276 * t311;
t324 = -t167 - t189;
t293 = -t200 * rSges(6,1) + t199 * rSges(6,2);
t136 = rSges(6,3) * t336 - t293;
t187 = -t276 * rSges(6,3) + (rSges(6,1) * t256 - rSges(6,2) * t255) * t273;
t102 = t276 * t136 + t187 * t336;
t198 = -t276 * rSges(5,3) + (rSges(5,1) * t259 - rSges(5,2) * t258) * t273;
t321 = -t189 - t198;
t219 = -t276 * rSges(4,3) + (rSges(4,1) * t275 - rSges(4,2) * t272) * t273;
t244 = t273 * pkin(2) - t276 * pkin(8);
t320 = -t219 - t244;
t318 = t269 * (t348 + t349) + t277 * t313;
t312 = t277 * pkin(1) + t274 * pkin(7);
t310 = t269 + t270;
t309 = t335 * t371 + t336 * t372;
t138 = t202 * rSges(6,1) - t201 * rSges(6,2) + rSges(6,3) * t335;
t308 = -t138 + t329;
t307 = -t187 + t324;
t306 = -t244 + t321;
t166 = t229 * rSges(4,1) + t228 * rSges(4,2) + rSges(4,3) * t335;
t264 = t277 * pkin(7);
t305 = t264 - t317;
t304 = t336 / 0.2e1;
t303 = t335 / 0.2e1;
t302 = -t232 * t276 - pkin(1);
t301 = t276 * t113 + t167 * t336 + t325;
t300 = -t327 + t329;
t299 = t274 * t163 + t277 * t164 + t318;
t298 = -t363 + t324;
t297 = -t244 + t307;
t66 = t276 * t364 + t336 * t363;
t296 = rSges(3,1) * t276 - rSges(3,2) * t273;
t295 = -t227 * rSges(4,1) - t226 * rSges(4,2);
t294 = -t207 * rSges(5,1) - t206 * rSges(5,2);
t292 = -t244 + t298;
t290 = -Icges(3,2) * t273 + t342;
t289 = Icges(3,5) * t276 - Icges(3,6) * t273;
t286 = rSges(3,1) * t332 - rSges(3,2) * t335 + t274 * rSges(3,3);
t284 = t64 / 0.2e1 + t98 / 0.2e1 + t90 / 0.2e1 + t76 / 0.2e1;
t283 = t91 / 0.2e1 + t77 / 0.2e1 + t65 / 0.2e1 + t99 / 0.2e1;
t282 = t274 * t113 + t277 * t114 + t299;
t281 = t276 * t346 + t309;
t280 = (-t368 - t380) * t304 + (t367 - t379) * t303;
t279 = -t266 * t335 + t312 + t319;
t278 = t371 * t352 + (t274 * t367 + t277 * t368) * t351 + t372 * t350 + t370 * t304 + t369 * t303;
t268 = t273 ^ 2;
t242 = rSges(2,1) * t277 - t274 * rSges(2,2);
t241 = -t274 * rSges(2,1) - rSges(2,2) * t277;
t240 = rSges(3,1) * t273 + rSges(3,2) * t276;
t236 = Icges(3,6) * t276 + t376;
t212 = Icges(3,3) * t274 + t277 * t289;
t211 = -Icges(3,3) * t277 + t274 * t289;
t179 = t286 + t312;
t178 = t344 + t264 + (-pkin(1) - t296) * t274;
t172 = t320 * t277;
t171 = t320 * t274;
t165 = rSges(4,3) * t336 - t295;
t152 = t277 * t286 + (t274 * t296 - t344) * t274;
t147 = rSges(5,3) * t336 - t294;
t120 = t136 * t335;
t118 = t166 + t312 + t313;
t117 = t264 + (-t349 - pkin(1) + (-rSges(4,3) - pkin(8)) * t273) * t274 + t295;
t116 = t306 * t277;
t115 = t306 * t274;
t112 = -t276 * t166 - t219 * t335;
t111 = t165 * t276 + t219 * t336;
t105 = t285 + t148 + t312;
t104 = t264 + (-rSges(5,3) * t273 - t257 * t276 - pkin(1)) * t274 + t294 + t314;
t103 = -t276 * t138 - t187 * t335;
t100 = (t165 * t277 - t166 * t274) * t273;
t95 = t279 + t138;
t94 = (-rSges(6,3) * t273 + t302) * t274 + t293 + t305;
t93 = t297 * t277;
t92 = t297 * t274;
t89 = -t138 * t336 + t120;
t88 = t274 * t165 + t166 * t277 + t318;
t83 = t292 * t277;
t82 = t292 * t274;
t75 = t279 + t327;
t74 = (-rSges(7,2) * t273 + t302) * t274 + t305 + t385;
t73 = t276 * t326 + t321 * t335;
t72 = t147 * t276 + t198 * t336 + t325;
t67 = -t276 * t327 - t335 * t363;
t51 = t151 + (t147 * t277 + t274 * t326) * t273;
t42 = -t327 * t336 + t328;
t41 = t274 * t147 + t148 * t277 + t299;
t40 = t276 * t308 + t307 * t335;
t39 = t301 + t102;
t38 = t71 * t274 - t277 * t70;
t37 = t69 * t274 - t277 * t68;
t36 = t276 * t300 + t298 * t335;
t35 = t66 + t301;
t34 = t308 * t336 + t120 + t330;
t33 = t274 * t136 + t138 * t277 + t282;
t30 = t55 * t274 - t277 * t54;
t29 = t53 * t274 - t277 * t52;
t18 = t300 * t336 + t328 + t330;
t17 = t274 * t364 + t327 * t277 + t282;
t1 = [Icges(2,3) + (Icges(3,1) * t273 - t183 * t255 + t342 + t361) * t273 + (Icges(3,4) * t273 + Icges(3,2) * t276 + t362 + t377) * t276 + m(7) * (t74 ^ 2 + t75 ^ 2) + m(6) * (t94 ^ 2 + t95 ^ 2) + m(5) * (t104 ^ 2 + t105 ^ 2) + m(4) * (t117 ^ 2 + t118 ^ 2) + m(3) * (t178 ^ 2 + t179 ^ 2) + m(2) * (t241 ^ 2 + t242 ^ 2) + t356 + t375; (-t62 / 0.2e1 - t60 / 0.2e1 - t84 / 0.2e1 - t85 / 0.2e1 + t274 * t290 * t351 - t284 + (-Icges(3,6) * t351 + t374 + t236 / 0.2e1) * t277) * t277 + (t63 / 0.2e1 + t61 / 0.2e1 + t86 / 0.2e1 + t87 / 0.2e1 + t276 * (Icges(3,6) * t274 + t277 * t290) / 0.2e1 + t274 * t374 + t236 * t352 + t283) * t274 + m(7) * (t74 * t83 + t75 * t82) + m(6) * (t92 * t95 + t93 * t94) + m(5) * (t104 * t116 + t105 * t115) + m(4) * (t117 * t172 + t118 * t171) + m(3) * (-t178 * t277 - t179 * t274) * t240; m(7) * (t17 ^ 2 + t82 ^ 2 + t83 ^ 2) + m(6) * (t33 ^ 2 + t92 ^ 2 + t93 ^ 2) + m(5) * (t115 ^ 2 + t116 ^ 2 + t41 ^ 2) + m(4) * (t171 ^ 2 + t172 ^ 2 + t88 ^ 2) + m(3) * (t240 ^ 2 * t310 + t152 ^ 2) + (-t270 * t211 - t29 - t37 - t370) * t277 + (t269 * t212 + t30 + t38 + (-t274 * t211 + t277 * t212) * t277 + t369) * t274; (-t331 - t366) * t276 + m(7) * (t35 * t74 + t36 * t75) + m(6) * (t39 * t94 + t40 * t95) + m(5) * (t104 * t72 + t105 * t73) + m(4) * (t111 * t117 + t112 * t118) + (t274 * t284 + t277 * t283) * t273 + t280; ((t38 / 0.2e1 + t30 / 0.2e1) * t277 + (t37 / 0.2e1 + t29 / 0.2e1) * t274) * t273 + m(7) * (t17 * t18 + t35 * t83 + t36 * t82) + m(6) * (t33 * t34 + t39 * t93 + t40 * t92) + m(5) * (t115 * t73 + t116 * t72 + t41 * t51) + m(4) * (t100 * t88 + t111 * t172 + t112 * t171) + t278 + t359 * t352 + (t357 * t274 + t358 * t277) * t351 + t360 * t350; (t331 * t276 + t346) * t276 + m(5) * (t51 ^ 2 + t72 ^ 2 + t73 ^ 2) + m(4) * (t100 ^ 2 + t111 ^ 2 + t112 ^ 2) + m(7) * (t18 ^ 2 + t35 ^ 2 + t36 ^ 2) + m(6) * (t34 ^ 2 + t39 ^ 2 + t40 ^ 2) + ((-t357 * t276 + t359) * t277 + (t358 * t276 + t360) * t274) * t273 + t309; 0.2e1 * ((t274 * t75 + t277 * t74) * t353 + (t274 * t95 + t277 * t94) * t354 + (t104 * t277 + t105 * t274) * t355) * t273; m(7) * (-t276 * t17 + (t274 * t82 + t277 * t83) * t273) + m(6) * (-t276 * t33 + (t274 * t92 + t277 * t93) * t273) + m(5) * (-t276 * t41 + (t115 * t274 + t116 * t277) * t273); m(7) * (-t276 * t18 + (t274 * t36 + t277 * t35) * t273) + m(6) * (-t276 * t34 + (t274 * t40 + t277 * t39) * t273) + m(5) * (-t276 * t51 + (t274 * t73 + t277 * t72) * t273); 0.2e1 * (t355 + t354 + t353) * (t268 * t310 + t276 ^ 2); -t378 + m(7) * (t66 * t74 + t67 * t75) + m(6) * (t102 * t94 + t103 * t95) + t280; m(7) * (t17 * t42 + t66 * t83 + t67 * t82) + m(6) * (t102 * t93 + t103 * t92 + t33 * t89) + t278; m(7) * (t18 * t42 + t35 * t66 + t36 * t67) + m(6) * (t102 * t39 + t103 * t40 + t34 * t89) + t281; m(6) * (-t89 * t276 + (t102 * t277 + t103 * t274) * t273) + m(7) * (-t42 * t276 + (t274 * t67 + t277 * t66) * t273); m(7) * (t42 ^ 2 + t66 ^ 2 + t67 ^ 2) + m(6) * (t102 ^ 2 + t103 ^ 2 + t89 ^ 2) + t281; m(7) * (t199 * t75 + t201 * t74); m(7) * (t17 * t339 + t199 * t82 + t201 * t83); m(7) * (t18 * t339 + t199 * t36 + t201 * t35); m(7) * (t199 * t274 + t201 * t277 - t255 * t276) * t273; m(7) * (t199 * t67 + t201 * t66 + t339 * t42); m(7) * (t255 ^ 2 * t268 + t199 ^ 2 + t201 ^ 2);];
%% Postprocessing: Reshape Output
% From vec2symmat_6_matlab.m
res = [t1(1) t1(2) t1(4) t1(7) t1(11) t1(16); t1(2) t1(3) t1(5) t1(8) t1(12) t1(17); t1(4) t1(5) t1(6) t1(9) t1(13) t1(18); t1(7) t1(8) t1(9) t1(10) t1(14) t1(19); t1(11) t1(12) t1(13) t1(14) t1(15) t1(20); t1(16) t1(17) t1(18) t1(19) t1(20) t1(21);];
Mq  = res;
