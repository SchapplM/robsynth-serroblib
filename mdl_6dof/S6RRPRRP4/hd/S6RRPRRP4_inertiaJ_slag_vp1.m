% Calculate joint inertia matrix for
% S6RRPRRP4
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% pkin [10x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d2,d4,d5,theta3]';
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
% Datum: 2019-03-09 11:56
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Mq = S6RRPRRP4_inertiaJ_slag_vp1(qJ, ...
  pkin, m, rSges, Icges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(10,1),zeros(7,1),zeros(7,3),zeros(7,6)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRPRRP4_inertiaJ_slag_vp1: qJ has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S6RRPRRP4_inertiaJ_slag_vp1: pkin has to be [10x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6RRPRRP4_inertiaJ_slag_vp1: m has to be [7x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [7,3]), ...
  'S6RRPRRP4_inertiaJ_slag_vp1: rSges has to be [7x3] (double)');
assert(isreal(Icges) && all(size(Icges) == [7 6]), ...
  'S6RRPRRP4_inertiaJ_slag_vp1: Icges has to be [7x6] (double)'); 

%% Symbolic Calculation
% From inertia_joint_joint_floatb_twist_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 11:52:36
% EndTime: 2019-03-09 11:52:47
% DurationCPUTime: 5.33s
% Computational Cost: add. (11835->467), mult. (11428->664), div. (0->0), fcn. (12390->10), ass. (0->236)
t228 = qJ(2) + pkin(10);
t222 = cos(t228);
t231 = qJ(4) + qJ(5);
t224 = cos(t231);
t238 = cos(qJ(1));
t291 = t238 * t224;
t223 = sin(t231);
t235 = sin(qJ(1));
t296 = t235 * t223;
t177 = t222 * t296 + t291;
t292 = t238 * t223;
t295 = t235 * t224;
t178 = t222 * t295 - t292;
t221 = sin(t228);
t300 = t221 * t235;
t104 = Icges(6,4) * t178 - Icges(6,2) * t177 + Icges(6,6) * t300;
t98 = Icges(7,5) * t178 + Icges(7,6) * t300 + Icges(7,3) * t177;
t355 = -t104 + t98;
t179 = t222 * t292 - t295;
t180 = t222 * t291 + t296;
t299 = t221 * t238;
t105 = Icges(6,4) * t180 - Icges(6,2) * t179 + Icges(6,6) * t299;
t99 = Icges(7,5) * t180 + Icges(7,6) * t299 + Icges(7,3) * t179;
t354 = -t105 + t99;
t100 = Icges(6,5) * t178 - Icges(6,6) * t177 + Icges(6,3) * t300;
t102 = Icges(7,4) * t178 + Icges(7,2) * t300 + Icges(7,6) * t177;
t353 = t100 + t102;
t101 = Icges(6,5) * t180 - Icges(6,6) * t179 + Icges(6,3) * t299;
t103 = Icges(7,4) * t180 + Icges(7,2) * t299 + Icges(7,6) * t179;
t352 = t101 + t103;
t106 = Icges(7,1) * t178 + Icges(7,4) * t300 + Icges(7,5) * t177;
t108 = Icges(6,1) * t178 - Icges(6,4) * t177 + Icges(6,5) * t300;
t351 = t106 + t108;
t107 = Icges(7,1) * t180 + Icges(7,4) * t299 + Icges(7,5) * t179;
t109 = Icges(6,1) * t180 - Icges(6,4) * t179 + Icges(6,5) * t299;
t350 = t107 + t109;
t329 = rSges(7,3) + qJ(6);
t337 = rSges(7,1) + pkin(5);
t349 = -t329 * t177 - t337 * t178;
t348 = t177 * t355 + t351 * t178 + t353 * t300;
t347 = t177 * t354 + t178 * t350 + t300 * t352;
t346 = t179 * t355 + t351 * t180 + t353 * t299;
t345 = t179 * t354 + t180 * t350 + t299 * t352;
t145 = -Icges(7,6) * t222 + (Icges(7,5) * t224 + Icges(7,3) * t223) * t221;
t147 = -Icges(7,2) * t222 + (Icges(7,4) * t224 + Icges(7,6) * t223) * t221;
t149 = -Icges(7,4) * t222 + (Icges(7,1) * t224 + Icges(7,5) * t223) * t221;
t69 = t145 * t177 + t147 * t300 + t149 * t178;
t146 = -Icges(6,3) * t222 + (Icges(6,5) * t224 - Icges(6,6) * t223) * t221;
t148 = -Icges(6,6) * t222 + (Icges(6,4) * t224 - Icges(6,2) * t223) * t221;
t150 = -Icges(6,5) * t222 + (Icges(6,1) * t224 - Icges(6,4) * t223) * t221;
t70 = t146 * t300 - t148 * t177 + t150 * t178;
t344 = -t70 - t69;
t71 = t145 * t179 + t147 * t299 + t149 * t180;
t72 = t146 * t299 - t148 * t179 + t150 * t180;
t343 = -t71 - t72;
t302 = t221 * t223;
t339 = t145 * t302 + (t149 + t150) * t221 * t224;
t341 = -t146 - t147;
t330 = -t148 * t302 + t222 * t341 + t339;
t342 = t330 * t222;
t340 = Icges(3,3) + Icges(4,3);
t234 = sin(qJ(2));
t237 = cos(qJ(2));
t338 = Icges(3,5) * t237 + Icges(4,5) * t222 - Icges(3,6) * t234 - Icges(4,6) * t221;
t336 = t344 * t222 + (t235 * t348 + t238 * t347) * t221;
t335 = t343 * t222 + (t235 * t346 + t238 * t345) * t221;
t334 = t235 * t347 - t238 * t348;
t333 = t235 * t345 - t238 * t346;
t51 = -t222 * t102 + (t106 * t224 + t223 * t98) * t221;
t53 = -t222 * t100 + (-t104 * t223 + t108 * t224) * t221;
t332 = -t51 - t53;
t52 = -t222 * t103 + (t107 * t224 + t223 * t99) * t221;
t54 = -t222 * t101 + (-t105 * t223 + t109 * t224) * t221;
t331 = t52 + t54;
t328 = rSges(7,2) * t300 - t349;
t327 = -t222 * rSges(7,2) + (t329 * t223 + t224 * t337) * t221;
t326 = t221 / 0.2e1;
t325 = t222 / 0.2e1;
t324 = t234 / 0.2e1;
t323 = t237 / 0.2e1;
t322 = -t338 * t235 + t340 * t238;
t321 = t340 * t235 + t338 * t238;
t229 = t235 ^ 2;
t230 = t238 ^ 2;
t320 = -t222 / 0.2e1;
t319 = t235 / 0.2e1;
t318 = -t238 / 0.2e1;
t317 = pkin(2) * t234;
t316 = pkin(3) * t222;
t315 = pkin(8) * t221;
t236 = cos(qJ(4));
t219 = t236 * pkin(4) + pkin(3);
t314 = -pkin(3) + t219;
t313 = t342 + (t332 * t235 - t331 * t238) * t221;
t311 = rSges(3,1) * t237;
t310 = rSges(3,2) * t234;
t309 = t238 * rSges(3,3);
t308 = t328 * t299;
t259 = -t178 * rSges(6,1) + t177 * rSges(6,2);
t111 = rSges(6,3) * t300 - t259;
t152 = -t222 * rSges(6,3) + (rSges(6,1) * t224 - rSges(6,2) * t223) * t221;
t85 = t222 * t111 + t152 * t300;
t307 = Icges(3,4) * t234;
t306 = Icges(3,4) * t237;
t305 = Icges(4,4) * t221;
t304 = Icges(4,4) * t222;
t233 = sin(qJ(4));
t158 = -Icges(5,6) * t222 + (Icges(5,4) * t236 - Icges(5,2) * t233) * t221;
t303 = t158 * t233;
t239 = -pkin(9) - pkin(8);
t298 = t221 * t239;
t297 = t222 * t238;
t294 = t235 * t233;
t293 = t235 * t236;
t232 = -qJ(3) - pkin(7);
t290 = t238 * t232;
t289 = t238 * t233;
t288 = t238 * t236;
t287 = rSges(7,2) * t299 + t329 * t179 + t337 * t180;
t113 = t180 * rSges(6,1) - t179 * rSges(6,2) + rSges(6,3) * t299;
t246 = pkin(4) * t294 + t219 * t297 - t238 * t298;
t280 = pkin(3) * t297 + pkin(8) * t299;
t122 = t246 - t280;
t286 = -t113 - t122;
t281 = -pkin(4) * t289 - t235 * t298;
t121 = (t222 * t314 - t315) * t235 + t281;
t143 = (pkin(8) + t239) * t222 + t314 * t221;
t285 = t222 * t121 + t143 * t300;
t220 = t237 * pkin(2) + pkin(1);
t213 = t238 * t220;
t227 = t238 * pkin(7);
t282 = t235 * (t290 + t227 + (-pkin(1) + t220) * t235) + t238 * (-t238 * pkin(1) + t213 + (-pkin(7) - t232) * t235);
t279 = t235 * rSges(3,3) + t238 * t311;
t278 = t229 + t230;
t189 = -t222 * t294 - t288;
t190 = t222 * t293 - t289;
t123 = Icges(5,5) * t190 + Icges(5,6) * t189 + Icges(5,3) * t300;
t125 = Icges(5,4) * t190 + Icges(5,2) * t189 + Icges(5,6) * t300;
t127 = Icges(5,1) * t190 + Icges(5,4) * t189 + Icges(5,5) * t300;
t61 = -t222 * t123 + (-t125 * t233 + t127 * t236) * t221;
t157 = -Icges(5,3) * t222 + (Icges(5,5) * t236 - Icges(5,6) * t233) * t221;
t159 = -Icges(5,5) * t222 + (Icges(5,1) * t236 - Icges(5,4) * t233) * t221;
t74 = t157 * t300 + t158 * t189 + t159 * t190;
t277 = t61 / 0.2e1 + t74 / 0.2e1;
t191 = -t222 * t289 + t293;
t192 = t222 * t288 + t294;
t124 = Icges(5,5) * t192 + Icges(5,6) * t191 + Icges(5,3) * t299;
t126 = Icges(5,4) * t192 + Icges(5,2) * t191 + Icges(5,6) * t299;
t128 = Icges(5,1) * t192 + Icges(5,4) * t191 + Icges(5,5) * t299;
t62 = -t222 * t124 + (-t126 * t233 + t128 * t236) * t221;
t75 = t157 * t299 + t158 * t191 + t159 * t192;
t276 = t62 / 0.2e1 + t75 / 0.2e1;
t275 = t335 * t299 + t336 * t300;
t274 = -t122 - t287;
t130 = t192 * rSges(5,1) + t191 * rSges(5,2) + rSges(5,3) * t299;
t273 = t300 / 0.2e1;
t272 = t299 / 0.2e1;
t271 = Icges(3,5) * t324 + Icges(4,5) * t326 + Icges(3,6) * t323 + Icges(4,6) * t325;
t270 = -t221 * rSges(4,1) - t222 * rSges(4,2) - t317;
t269 = -t221 * pkin(3) + t222 * pkin(8) - t317;
t268 = -t235 * t232 + t213;
t267 = -t219 * t222 - t220;
t59 = t222 * t328 + t327 * t300;
t266 = t229 * (t315 + t316) + t238 * t280 + t282;
t265 = -t143 + t269;
t160 = -t222 * rSges(5,3) + (rSges(5,1) * t236 - rSges(5,2) * t233) * t221;
t264 = -t160 + t269;
t263 = -t281 - t290;
t262 = -t310 + t311;
t261 = rSges(4,1) * t222 - rSges(4,2) * t221;
t260 = -t190 * rSges(5,1) - t189 * rSges(5,2);
t258 = Icges(3,1) * t237 - t307;
t257 = Icges(4,1) * t222 - t305;
t256 = -Icges(3,2) * t234 + t306;
t255 = -Icges(4,2) * t221 + t304;
t248 = -t152 + t265;
t247 = rSges(4,1) * t297 - rSges(4,2) * t299 + t235 * rSges(4,3);
t245 = t235 * t121 + t238 * t122 + t266;
t244 = t265 - t327;
t243 = t222 * t313 + t275;
t242 = (-t332 - t344) * t273 + (t331 - t343) * t272;
t241 = (t331 * t235 + t332 * t238) * t320 + t335 * t319 + t336 * t318 + t334 * t273 + t333 * t272;
t240 = t246 + t268;
t206 = t238 * rSges(2,1) - t235 * rSges(2,2);
t205 = -t235 * rSges(2,1) - t238 * rSges(2,2);
t204 = t234 * rSges(3,1) + t237 * rSges(3,2);
t162 = t270 * t238;
t161 = t270 * t235;
t156 = t235 * pkin(7) + (pkin(1) - t310) * t238 + t279;
t155 = t309 + t227 + (-pkin(1) - t262) * t235;
t142 = t221 * t236 * t159;
t139 = t247 + t268;
t138 = (rSges(4,3) - t232) * t238 + (-t220 - t261) * t235;
t133 = t238 * (-t238 * t310 + t279) + (t235 * t262 - t309) * t235;
t129 = rSges(5,3) * t300 - t260;
t117 = t264 * t238;
t116 = t264 * t235;
t97 = t121 * t299;
t94 = t111 * t299;
t92 = t268 + t130 + t280;
t91 = -t290 + (-t316 - t220 + (-rSges(5,3) - pkin(8)) * t221) * t235 + t260;
t90 = -t130 * t222 - t160 * t299;
t89 = t129 * t222 + t160 * t300;
t88 = t248 * t238;
t87 = t248 * t235;
t86 = -t222 * t113 - t152 * t299;
t84 = -t222 * t157 - t221 * t303 + t142;
t83 = t240 + t113;
t82 = (-rSges(6,3) * t221 + t267) * t235 + t259 + t263;
t81 = t238 * t247 + (-t238 * rSges(4,3) + t235 * t261) * t235 + t282;
t80 = (t129 * t238 - t130 * t235) * t221;
t77 = t244 * t238;
t76 = t244 * t235;
t73 = -t113 * t300 + t94;
t68 = t240 + t287;
t67 = (-rSges(7,2) * t221 + t267) * t235 + t263 + t349;
t60 = -t222 * t287 - t299 * t327;
t58 = t124 * t299 + t126 * t191 + t128 * t192;
t57 = t123 * t299 + t125 * t191 + t127 * t192;
t56 = t124 * t300 + t126 * t189 + t128 * t190;
t55 = t123 * t300 + t125 * t189 + t127 * t190;
t50 = t286 * t222 + (-t143 - t152) * t299;
t49 = t285 + t85;
t48 = t129 * t235 + t130 * t238 + t266;
t35 = -t287 * t300 + t308;
t34 = t286 * t300 + t94 + t97;
t33 = t274 * t222 + (-t143 - t327) * t299;
t32 = t59 + t285;
t31 = t111 * t235 + t113 * t238 + t245;
t30 = t274 * t300 + t308 + t97;
t29 = t235 * t58 - t238 * t57;
t28 = t235 * t56 - t238 * t55;
t25 = t235 * t328 + t287 * t238 + t245;
t16 = -t75 * t222 + (t235 * t57 + t238 * t58) * t221;
t15 = -t74 * t222 + (t235 * t55 + t238 * t56) * t221;
t1 = [t237 * (Icges(3,2) * t237 + t307) + t234 * (Icges(3,1) * t234 + t306) + Icges(2,3) + t142 + (Icges(4,1) * t221 - t148 * t223 - t303 + t304) * t221 + (Icges(4,2) * t222 - t157 + t305 + t341) * t222 + m(7) * (t67 ^ 2 + t68 ^ 2) + m(6) * (t82 ^ 2 + t83 ^ 2) + m(5) * (t91 ^ 2 + t92 ^ 2) + m(4) * (t138 ^ 2 + t139 ^ 2) + m(3) * (t155 ^ 2 + t156 ^ 2) + m(2) * (t205 ^ 2 + t206 ^ 2) + t339; (-t237 * (-Icges(3,6) * t238 + t235 * t256) / 0.2e1 - t234 * (-Icges(3,5) * t238 + t235 * t258) / 0.2e1 - t69 / 0.2e1 - t70 / 0.2e1 + (-Icges(4,6) * t238 + t235 * t255) * t320 - t221 * (-Icges(4,5) * t238 + t235 * t257) / 0.2e1 - t53 / 0.2e1 - t51 / 0.2e1 + t271 * t238 - t277) * t238 + ((Icges(3,6) * t235 + t238 * t256) * t323 + (Icges(3,5) * t235 + t238 * t258) * t324 + t71 / 0.2e1 + t72 / 0.2e1 + (Icges(4,6) * t235 + t238 * t255) * t325 + (Icges(4,5) * t235 + t238 * t257) * t326 + t54 / 0.2e1 + t52 / 0.2e1 + t271 * t235 + t276) * t235 + m(7) * (t67 * t77 + t68 * t76) + m(6) * (t82 * t88 + t83 * t87) + m(5) * (t116 * t92 + t117 * t91) + m(4) * (t138 * t162 + t139 * t161) + m(3) * (-t155 * t238 - t156 * t235) * t204; m(7) * (t25 ^ 2 + t76 ^ 2 + t77 ^ 2) + m(6) * (t31 ^ 2 + t87 ^ 2 + t88 ^ 2) + m(5) * (t116 ^ 2 + t117 ^ 2 + t48 ^ 2) + m(4) * (t161 ^ 2 + t162 ^ 2 + t81 ^ 2) + m(3) * (t204 ^ 2 * t278 + t133 ^ 2) + (t230 * t322 - t28 - t334) * t238 + (t29 + t321 * t229 + (t235 * t322 + t238 * t321) * t238 + t333) * t235; m(7) * (t235 * t67 - t238 * t68) + m(6) * (t235 * t82 - t238 * t83) + m(5) * (t235 * t91 - t238 * t92) + m(4) * (t138 * t235 - t139 * t238); m(7) * (t235 * t77 - t238 * t76) + m(6) * (t235 * t88 - t238 * t87) + m(5) * (-t116 * t238 + t117 * t235) + m(4) * (-t161 * t238 + t162 * t235); 0.2e1 * (m(4) / 0.2e1 + m(5) / 0.2e1 + m(6) / 0.2e1 + m(7) / 0.2e1) * t278; (-t84 - t330) * t222 + m(7) * (t32 * t67 + t33 * t68) + m(6) * (t49 * t82 + t50 * t83) + m(5) * (t89 * t91 + t90 * t92) + (t235 * t277 + t238 * t276) * t221 + t242; t241 + m(7) * (t30 * t25 + t32 * t77 + t33 * t76) + m(6) * (t34 * t31 + t49 * t88 + t50 * t87) + m(5) * (t116 * t90 + t117 * t89 + t48 * t80) + (t28 * t319 + t238 * t29 / 0.2e1) * t221 + (t62 * t235 - t61 * t238) * t320 + t15 * t318 + t16 * t319; m(5) * (t235 * t89 - t238 * t90) + m(6) * (t235 * t49 - t238 * t50) + m(7) * (t235 * t32 - t238 * t33); (t84 * t222 + t313) * t222 + (t235 * t15 + t238 * t16 - t222 * (t235 * t61 + t238 * t62)) * t221 + m(7) * (t30 ^ 2 + t32 ^ 2 + t33 ^ 2) + m(6) * (t34 ^ 2 + t49 ^ 2 + t50 ^ 2) + m(5) * (t80 ^ 2 + t89 ^ 2 + t90 ^ 2) + t275; -t342 + m(7) * (t59 * t67 + t60 * t68) + m(6) * (t82 * t85 + t83 * t86) + t242; m(7) * (t35 * t25 + t59 * t77 + t60 * t76) + m(6) * (t31 * t73 + t85 * t88 + t86 * t87) + t241; m(6) * (t235 * t85 - t238 * t86) + m(7) * (t235 * t59 - t238 * t60); m(7) * (t30 * t35 + t32 * t59 + t33 * t60) + m(6) * (t34 * t73 + t49 * t85 + t50 * t86) + t243; m(7) * (t35 ^ 2 + t59 ^ 2 + t60 ^ 2) + m(6) * (t73 ^ 2 + t85 ^ 2 + t86 ^ 2) + t243; m(7) * (t177 * t68 + t179 * t67); m(7) * (t177 * t76 + t179 * t77 + t25 * t302); m(7) * (-t177 * t238 + t179 * t235); m(7) * (t177 * t33 + t179 * t32 + t30 * t302); m(7) * (t177 * t60 + t179 * t59 + t302 * t35); m(7) * (t221 ^ 2 * t223 ^ 2 + t177 ^ 2 + t179 ^ 2);];
%% Postprocessing: Reshape Output
% From vec2symmat_6_matlab.m
res = [t1(1) t1(2) t1(4) t1(7) t1(11) t1(16); t1(2) t1(3) t1(5) t1(8) t1(12) t1(17); t1(4) t1(5) t1(6) t1(9) t1(13) t1(18); t1(7) t1(8) t1(9) t1(10) t1(14) t1(19); t1(11) t1(12) t1(13) t1(14) t1(15) t1(20); t1(16) t1(17) t1(18) t1(19) t1(20) t1(21);];
Mq  = res;
