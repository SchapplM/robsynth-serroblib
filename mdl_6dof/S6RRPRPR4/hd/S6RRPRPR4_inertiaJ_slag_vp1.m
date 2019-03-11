% Calculate joint inertia matrix for
% S6RRPRPR4
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% pkin [12x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,alpha2,d1,d2,d4,d6,theta3,theta5]';
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
% Datum: 2019-03-09 10:27
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Mq = S6RRPRPR4_inertiaJ_slag_vp1(qJ, ...
  pkin, m, rSges, Icges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(12,1),zeros(7,1),zeros(7,3),zeros(7,6)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRPRPR4_inertiaJ_slag_vp1: qJ has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [12 1]), ...
  'S6RRPRPR4_inertiaJ_slag_vp1: pkin has to be [12x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6RRPRPR4_inertiaJ_slag_vp1: m has to be [7x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [7,3]), ...
  'S6RRPRPR4_inertiaJ_slag_vp1: rSges has to be [7x3] (double)');
assert(isreal(Icges) && all(size(Icges) == [7 6]), ...
  'S6RRPRPR4_inertiaJ_slag_vp1: Icges has to be [7x6] (double)'); 

%% Symbolic Calculation
% From inertia_joint_joint_floatb_twist_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 10:21:22
% EndTime: 2019-03-09 10:21:40
% DurationCPUTime: 7.11s
% Computational Cost: add. (25458->623), mult. (52771->871), div. (0->0), fcn. (68534->14), ass. (0->290)
t335 = sin(pkin(11));
t336 = cos(pkin(11));
t342 = sin(qJ(2));
t343 = cos(qJ(2));
t256 = -t342 * t335 + t343 * t336;
t273 = sin(pkin(6));
t244 = t256 * t273;
t283 = t335 * t343 + t336 * t342;
t245 = t283 * t273;
t274 = cos(pkin(6));
t202 = Icges(4,5) * t245 + Icges(4,6) * t244 + Icges(4,3) * t274;
t203 = Icges(4,4) * t245 + Icges(4,2) * t244 + Icges(4,6) * t274;
t204 = Icges(4,1) * t245 + Icges(4,4) * t244 + Icges(4,5) * t274;
t240 = Icges(3,3) * t274 + (Icges(3,5) * t342 + Icges(3,6) * t343) * t273;
t241 = Icges(3,6) * t274 + (Icges(3,4) * t342 + Icges(3,2) * t343) * t273;
t242 = Icges(3,5) * t274 + (Icges(3,1) * t342 + Icges(3,4) * t343) * t273;
t308 = t273 * t342;
t355 = t273 * t343 * t241 + t244 * t203 + t245 * t204 + t242 * t308 + (t202 + t240) * t274;
t349 = m(7) / 0.2e1;
t350 = m(6) / 0.2e1;
t319 = t350 + t349;
t354 = 0.2e1 * t319;
t353 = t273 ^ 2;
t352 = m(4) / 0.2e1;
t351 = m(5) / 0.2e1;
t246 = t283 * t274;
t278 = sin(qJ(1));
t281 = cos(qJ(1));
t228 = t246 * t281 + t278 * t256;
t320 = qJ(4) + pkin(12);
t271 = sin(t320);
t302 = cos(t320);
t292 = t273 * t302;
t187 = t228 * t271 + t281 * t292;
t231 = -t278 * t246 + t256 * t281;
t189 = t231 * t271 - t278 * t292;
t232 = t245 * t271 - t274 * t302;
t331 = t273 * t281;
t188 = t228 * t302 - t271 * t331;
t282 = t274 * t256;
t226 = -t278 * t283 + t281 * t282;
t276 = sin(qJ(6));
t279 = cos(qJ(6));
t145 = -t188 * t276 - t226 * t279;
t146 = t188 * t279 - t226 * t276;
t84 = Icges(7,5) * t146 + Icges(7,6) * t145 + Icges(7,3) * t187;
t86 = Icges(7,4) * t146 + Icges(7,2) * t145 + Icges(7,6) * t187;
t88 = Icges(7,1) * t146 + Icges(7,4) * t145 + Icges(7,5) * t187;
t26 = t145 * t86 + t146 * t88 + t187 * t84;
t332 = t273 * t278;
t190 = t231 * t302 + t271 * t332;
t229 = -t278 * t282 - t281 * t283;
t147 = -t190 * t276 - t229 * t279;
t148 = t190 * t279 - t229 * t276;
t85 = Icges(7,5) * t148 + Icges(7,6) * t147 + Icges(7,3) * t189;
t87 = Icges(7,4) * t148 + Icges(7,2) * t147 + Icges(7,6) * t189;
t89 = Icges(7,1) * t148 + Icges(7,4) * t147 + Icges(7,5) * t189;
t27 = t145 * t87 + t146 * t89 + t187 * t85;
t233 = t245 * t302 + t274 * t271;
t181 = -t233 * t276 - t244 * t279;
t182 = t233 * t279 - t244 * t276;
t109 = Icges(7,5) * t182 + Icges(7,6) * t181 + Icges(7,3) * t232;
t110 = Icges(7,4) * t182 + Icges(7,2) * t181 + Icges(7,6) * t232;
t111 = Icges(7,1) * t182 + Icges(7,4) * t181 + Icges(7,5) * t232;
t38 = t109 * t187 + t110 * t145 + t111 * t146;
t1 = t187 * t26 + t189 * t27 + t232 * t38;
t348 = -t1 / 0.2e1;
t33 = t181 * t86 + t182 * t88 + t232 * t84;
t34 = t181 * t87 + t182 * t89 + t232 * t85;
t46 = t232 * t109 + t181 * t110 + t182 * t111;
t40 = t46 * t232;
t7 = t33 * t187 + t34 * t189 + t40;
t347 = t7 / 0.2e1;
t346 = t187 / 0.2e1;
t345 = t189 / 0.2e1;
t344 = t232 / 0.2e1;
t341 = pkin(1) * t281;
t340 = t188 * pkin(5);
t280 = cos(qJ(4));
t269 = pkin(4) * t280 + pkin(3);
t339 = -pkin(3) + t269;
t293 = -t146 * rSges(7,1) - t145 * rSges(7,2);
t90 = t187 * rSges(7,3) - t293;
t338 = t187 * pkin(10) + t340 + t90;
t91 = t148 * rSges(7,1) + t147 * rSges(7,2) + t189 * rSges(7,3);
t337 = t190 * pkin(5) + pkin(10) * t189 + t91;
t275 = -qJ(5) - pkin(9);
t334 = t226 * t275;
t247 = t274 * t342 * pkin(2) + (-pkin(8) - qJ(3)) * t273;
t333 = t247 * t281;
t277 = sin(qJ(4));
t330 = t274 * t277;
t270 = pkin(2) * t343 + pkin(1);
t329 = t278 * t270;
t328 = t355 * t274;
t112 = rSges(7,1) * t182 + rSges(7,2) * t181 + rSges(7,3) * t232;
t327 = -pkin(5) * t233 - pkin(10) * t232 - t112;
t158 = Icges(4,5) * t228 + Icges(4,6) * t226 - Icges(4,3) * t331;
t305 = t281 * t343;
t306 = t278 * t342;
t251 = t274 * t305 - t306;
t304 = t281 * t342;
t307 = t278 * t343;
t252 = t274 * t304 + t307;
t209 = Icges(3,5) * t252 + Icges(3,6) * t251 - Icges(3,3) * t331;
t326 = -t158 - t209;
t159 = Icges(4,5) * t231 + Icges(4,6) * t229 + Icges(4,3) * t332;
t253 = -t274 * t307 - t304;
t254 = -t274 * t306 + t305;
t210 = Icges(3,5) * t254 + Icges(3,6) * t253 + Icges(3,3) * t332;
t325 = t159 + t210;
t264 = t281 * t270;
t219 = -t341 + t264 + (-t273 * pkin(8) - t247) * t278;
t207 = t274 * t219;
t303 = -t231 * pkin(3) + pkin(9) * t229;
t324 = -t274 * t303 + t207;
t224 = t226 * pkin(9);
t173 = t228 * pkin(3) - t224;
t268 = pkin(8) * t331;
t218 = t333 + t268 + (-pkin(1) + t270) * t278;
t323 = -t173 - t218;
t322 = t218 * t332 + t219 * t331;
t257 = pkin(2) * t308 + t274 * qJ(3);
t321 = -pkin(3) * t245 + pkin(9) * t244 - t257;
t318 = t33 / 0.2e1 + t38 / 0.2e1;
t39 = t109 * t189 + t110 * t147 + t111 * t148;
t317 = t39 / 0.2e1 + t34 / 0.2e1;
t316 = t277 * t332;
t315 = t277 * t331;
t310 = pkin(4) * t316 + t229 * t275 + t231 * t269;
t108 = t303 + t310;
t314 = t274 * t108 + t324;
t260 = pkin(4) * t315;
t107 = t228 * t339 + t224 - t260 + t334;
t313 = -t107 + t323;
t154 = Icges(6,5) * t233 - Icges(6,6) * t232 - Icges(6,3) * t244;
t155 = Icges(6,4) * t233 - Icges(6,2) * t232 - Icges(6,6) * t244;
t156 = Icges(6,1) * t233 - Icges(6,4) * t232 - Icges(6,5) * t244;
t73 = -t244 * t154 - t232 * t155 + t233 * t156;
t234 = -t245 * t277 + t274 * t280;
t235 = t245 * t280 + t330;
t169 = Icges(5,5) * t235 + Icges(5,6) * t234 - Icges(5,3) * t244;
t170 = Icges(5,4) * t235 + Icges(5,2) * t234 - Icges(5,6) * t244;
t171 = Icges(5,1) * t235 + Icges(5,4) * t234 - Icges(5,5) * t244;
t76 = -t244 * t169 + t234 * t170 + t235 * t171;
t153 = pkin(4) * t330 + t339 * t245 - (-pkin(9) - t275) * t244;
t312 = -t153 + t321;
t120 = t190 * rSges(6,1) - t189 * rSges(6,2) - t229 * rSges(6,3);
t197 = -t231 * t277 + t280 * t332;
t198 = t231 * t280 + t316;
t128 = t198 * rSges(5,1) + t197 * rSges(5,2) - t229 * rSges(5,3);
t165 = t231 * rSges(4,1) + t229 * rSges(4,2) + rSges(4,3) * t332;
t216 = t254 * rSges(3,1) + t253 * rSges(3,2) + rSges(3,3) * t332;
t301 = t273 * (-rSges(4,1) * t245 - rSges(4,2) * t244 - rSges(4,3) * t274 - t257);
t300 = -t278 * t247 + t264;
t299 = t173 * t332 - t303 * t331 + t322;
t172 = rSges(5,1) * t235 + rSges(5,2) * t234 - rSges(5,3) * t244;
t298 = t273 * (-t172 + t321);
t295 = -t228 * rSges(4,1) - t226 * rSges(4,2);
t294 = -t188 * rSges(6,1) + t187 * rSges(6,2);
t291 = -t329 - t333;
t157 = rSges(6,1) * t233 - rSges(6,2) * t232 - rSges(6,3) * t244;
t290 = t273 * (-t157 + t312);
t289 = t107 * t332 + t108 * t331 + t299;
t288 = t273 * (t312 + t327);
t287 = t300 + t310;
t195 = -t228 * t277 - t280 * t331;
t196 = t228 * t280 - t315;
t127 = t196 * rSges(5,1) + t195 * rSges(5,2) - t226 * rSges(5,3);
t215 = t252 * rSges(3,1) + t251 * rSges(3,2) - rSges(3,3) * t331;
t113 = Icges(6,5) * t188 - Icges(6,6) * t187 - Icges(6,3) * t226;
t115 = Icges(6,4) * t188 - Icges(6,2) * t187 - Icges(6,6) * t226;
t117 = Icges(6,1) * t188 - Icges(6,4) * t187 - Icges(6,5) * t226;
t56 = -t113 * t244 - t115 * t232 + t117 * t233;
t121 = Icges(5,5) * t196 + Icges(5,6) * t195 - Icges(5,3) * t226;
t123 = Icges(5,4) * t196 + Icges(5,2) * t195 - Icges(5,6) * t226;
t125 = Icges(5,1) * t196 + Icges(5,4) * t195 - Icges(5,5) * t226;
t58 = -t121 * t244 + t123 * t234 + t125 * t235;
t66 = -t154 * t226 - t155 * t187 + t156 * t188;
t68 = -t169 * t226 + t170 * t195 + t171 * t196;
t286 = -t68 / 0.2e1 - t66 / 0.2e1 - t58 / 0.2e1 - t56 / 0.2e1 - t318;
t114 = Icges(6,5) * t190 - Icges(6,6) * t189 - Icges(6,3) * t229;
t116 = Icges(6,4) * t190 - Icges(6,2) * t189 - Icges(6,6) * t229;
t118 = Icges(6,1) * t190 - Icges(6,4) * t189 - Icges(6,5) * t229;
t57 = -t114 * t244 - t116 * t232 + t118 * t233;
t122 = Icges(5,5) * t198 + Icges(5,6) * t197 - Icges(5,3) * t229;
t124 = Icges(5,4) * t198 + Icges(5,2) * t197 - Icges(5,6) * t229;
t126 = Icges(5,1) * t198 + Icges(5,4) * t197 - Icges(5,5) * t229;
t59 = -t122 * t244 + t124 * t234 + t126 * t235;
t67 = -t154 * t229 - t155 * t189 + t156 * t190;
t69 = -t169 * t229 + t170 * t197 + t171 * t198;
t285 = -t59 / 0.2e1 - t57 / 0.2e1 - t69 / 0.2e1 - t67 / 0.2e1 - t317;
t284 = -t228 * t269 + t260 + t291;
t262 = rSges(2,1) * t281 - t278 * rSges(2,2);
t261 = -t278 * rSges(2,1) - rSges(2,2) * t281;
t243 = t274 * rSges(3,3) + (rSges(3,1) * t342 + rSges(3,2) * t343) * t273;
t214 = Icges(3,1) * t254 + Icges(3,4) * t253 + Icges(3,5) * t332;
t213 = Icges(3,1) * t252 + Icges(3,4) * t251 - Icges(3,5) * t331;
t212 = Icges(3,4) * t254 + Icges(3,2) * t253 + Icges(3,6) * t332;
t211 = Icges(3,4) * t252 + Icges(3,2) * t251 - Icges(3,6) * t331;
t194 = pkin(8) * t332 + t216 + t341;
t193 = -t278 * pkin(1) - t215 + t268;
t177 = -t274 * t215 - t243 * t331;
t176 = t216 * t274 - t243 * t332;
t164 = -rSges(4,3) * t331 - t295;
t163 = Icges(4,1) * t231 + Icges(4,4) * t229 + Icges(4,5) * t332;
t162 = Icges(4,1) * t228 + Icges(4,4) * t226 - Icges(4,5) * t331;
t161 = Icges(4,4) * t231 + Icges(4,2) * t229 + Icges(4,6) * t332;
t160 = Icges(4,4) * t228 + Icges(4,2) * t226 - Icges(4,6) * t331;
t152 = (t215 * t278 + t216 * t281) * t273;
t151 = t240 * t332 + t241 * t253 + t242 * t254;
t150 = -t240 * t331 + t251 * t241 + t252 * t242;
t135 = t300 + t165;
t134 = -t329 + (rSges(4,3) * t273 - t247) * t281 + t295;
t131 = t226 * t153;
t130 = t274 * t210 + (t212 * t343 + t214 * t342) * t273;
t129 = t274 * t209 + (t211 * t343 + t213 * t342) * t273;
t119 = -t226 * rSges(6,3) - t294;
t101 = t244 * t108;
t100 = (-t164 - t218) * t274 + t281 * t301;
t99 = t165 * t274 + t278 * t301 + t207;
t96 = t229 * t107;
t95 = t202 * t332 + t203 * t229 + t204 * t231;
t94 = -t202 * t331 + t226 * t203 + t228 * t204;
t93 = t300 - t303 + t128;
t92 = -t127 - t173 + t291;
t83 = t287 + t120;
t82 = (rSges(6,3) - t275) * t226 + t284 + t294;
t81 = (t164 * t278 + t165 * t281) * t273 + t322;
t80 = -t128 * t244 + t172 * t229;
t79 = t127 * t244 - t172 * t226;
t78 = t159 * t274 + t161 * t244 + t163 * t245;
t77 = t158 * t274 + t160 * t244 + t162 * t245;
t75 = t76 * t274;
t74 = t76 * t244;
t72 = t73 * t274;
t71 = -t127 * t229 + t128 * t226;
t70 = t73 * t244;
t65 = (-t127 + t323) * t274 + t281 * t298;
t64 = t128 * t274 + t278 * t298 + t324;
t63 = t287 + t337;
t62 = -t340 - t334 + (-rSges(7,3) - pkin(10)) * t187 + t284 + t293;
t61 = -t112 * t189 + t232 * t91;
t60 = t112 * t187 - t232 * t90;
t55 = (t127 * t278 + t128 * t281) * t273 + t299;
t54 = -t122 * t229 + t124 * t197 + t126 * t198;
t53 = -t121 * t229 + t123 * t197 + t125 * t198;
t52 = -t122 * t226 + t124 * t195 + t126 * t196;
t51 = -t121 * t226 + t123 * t195 + t125 * t196;
t50 = -t114 * t229 - t116 * t189 + t118 * t190;
t49 = -t113 * t229 - t115 * t189 + t117 * t190;
t48 = -t114 * t226 - t116 * t187 + t118 * t188;
t47 = -t113 * t226 - t115 * t187 + t117 * t188;
t45 = t46 * t274;
t44 = -t120 * t244 - t101 + (t153 + t157) * t229;
t43 = -t157 * t226 - t131 - (-t107 - t119) * t244;
t42 = -t187 * t91 + t189 * t90;
t41 = t46 * t244;
t37 = (-t119 + t313) * t274 + t281 * t290;
t36 = t120 * t274 + t278 * t290 + t314;
t35 = -t119 * t229 - t96 + (t108 + t120) * t226;
t32 = (t119 * t278 + t120 * t281) * t273 + t289;
t31 = -t101 - t337 * t244 + (t153 - t327) * t229;
t30 = -t131 + t327 * t226 - (-t107 - t338) * t244;
t29 = t147 * t87 + t148 * t89 + t189 * t85;
t28 = t147 * t86 + t148 * t88 + t189 * t84;
t25 = (t313 - t338) * t274 + t281 * t288;
t24 = t274 * t337 + t278 * t288 + t314;
t23 = -t96 - t338 * t229 + (t108 + t337) * t226;
t22 = t75 + (t59 * t278 - t58 * t281) * t273;
t21 = (t278 * t338 + t281 * t337) * t273 + t289;
t20 = -t58 * t226 - t59 * t229 - t74;
t19 = t72 + (t57 * t278 - t56 * t281) * t273;
t18 = -t56 * t226 - t57 * t229 - t70;
t17 = t69 * t274 + (t278 * t54 - t281 * t53) * t273;
t16 = t68 * t274 + (t278 * t52 - t281 * t51) * t273;
t15 = -t226 * t53 - t229 * t54 - t244 * t69;
t14 = -t226 * t51 - t229 * t52 - t244 * t68;
t13 = t67 * t274 + (t278 * t50 - t281 * t49) * t273;
t12 = t66 * t274 + (t278 * t48 - t281 * t47) * t273;
t11 = -t226 * t49 - t229 * t50 - t244 * t67;
t10 = -t226 * t47 - t229 * t48 - t244 * t66;
t9 = t45 + (t34 * t278 - t33 * t281) * t273;
t8 = -t33 * t226 - t34 * t229 - t41;
t6 = t39 * t274 + (t278 * t29 - t28 * t281) * t273;
t5 = t38 * t274 + (-t26 * t281 + t27 * t278) * t273;
t4 = -t226 * t28 - t229 * t29 - t244 * t39;
t3 = -t226 * t26 - t229 * t27 - t244 * t38;
t2 = t187 * t28 + t189 * t29 + t232 * t39;
t97 = [m(7) * (t62 ^ 2 + t63 ^ 2) + m(6) * (t82 ^ 2 + t83 ^ 2) + m(5) * (t92 ^ 2 + t93 ^ 2) + m(4) * (t134 ^ 2 + t135 ^ 2) + m(3) * (t193 ^ 2 + t194 ^ 2) + m(2) * (t261 ^ 2 + t262 ^ 2) + Icges(2,3) + t76 + t73 + t46 + t355; t75 + t72 + t45 + m(7) * (t24 * t63 + t25 * t62) + m(6) * (t36 * t83 + t37 * t82) + m(5) * (t64 * t93 + t65 * t92) + m(4) * (t100 * t134 + t135 * t99) + m(3) * (t176 * t194 + t177 * t193) + ((-t129 / 0.2e1 - t77 / 0.2e1 - t94 / 0.2e1 - t150 / 0.2e1 + t286) * t281 + (t130 / 0.2e1 + t78 / 0.2e1 + t95 / 0.2e1 + t151 / 0.2e1 - t285) * t278) * t273 + t328; (t9 + t22 + t19 + t328) * t274 + m(7) * (t21 ^ 2 + t24 ^ 2 + t25 ^ 2) + m(6) * (t32 ^ 2 + t36 ^ 2 + t37 ^ 2) + m(5) * (t55 ^ 2 + t64 ^ 2 + t65 ^ 2) + m(4) * (t100 ^ 2 + t81 ^ 2 + t99 ^ 2) + m(3) * (t152 ^ 2 + t176 ^ 2 + t177 ^ 2) + ((-t12 - t16 - t5 + ((t226 * t160 + t228 * t162 + t251 * t211 + t252 * t213) * t273 + t326 * t353 * t281) * t281 + (-t129 - t150 - t77 - t94) * t274) * t281 + (t6 + t13 + t17 + ((t229 * t161 + t231 * t163 + t253 * t212 + t254 * t214) * t273 + t325 * t353 * t278) * t278 + (t95 + t151 + t130 + t78) * t274 + (-t229 * t160 - t226 * t161 - t231 * t162 - t228 * t163 - t253 * t211 - t251 * t212 - t254 * t213 - t252 * t214 + (t278 * t326 + t281 * t325) * t273) * t331) * t278) * t273; 0.2e1 * ((t278 * t62 - t281 * t63) * t349 + (t278 * t82 - t281 * t83) * t350 + (t278 * t92 - t281 * t93) * t351 + (t134 * t278 - t135 * t281) * t352) * t273; m(7) * (t274 * t21 + (-t24 * t281 + t25 * t278) * t273) + m(6) * (t274 * t32 + (t278 * t37 - t281 * t36) * t273) + m(5) * (t274 * t55 + (t278 * t65 - t281 * t64) * t273) + m(4) * (t274 * t81 + (t100 * t278 - t281 * t99) * t273); 0.2e1 * (t352 + t351 + t319) * (t274 ^ 2 + (t278 ^ 2 + t281 ^ 2) * t353); -t41 - t74 - t70 + m(7) * (t30 * t62 + t31 * t63) + m(6) * (t43 * t82 + t44 * t83) + m(5) * (t79 * t92 + t80 * t93) + t285 * t229 + t286 * t226; (t8 / 0.2e1 + t20 / 0.2e1 + t18 / 0.2e1) * t274 - (t9 / 0.2e1 + t22 / 0.2e1 + t19 / 0.2e1) * t244 + (-t6 / 0.2e1 - t17 / 0.2e1 - t13 / 0.2e1) * t229 + (-t5 / 0.2e1 - t12 / 0.2e1 - t16 / 0.2e1) * t226 + m(7) * (t21 * t23 + t24 * t31 + t25 * t30) + m(6) * (t32 * t35 + t36 * t44 + t37 * t43) + m(5) * (t55 * t71 + t64 * t80 + t65 * t79) + ((-t3 / 0.2e1 - t10 / 0.2e1 - t14 / 0.2e1) * t281 + (t4 / 0.2e1 + t11 / 0.2e1 + t15 / 0.2e1) * t278) * t273; m(5) * (t71 * t274 + (t278 * t79 - t281 * t80) * t273) + m(6) * (t35 * t274 + (t278 * t43 - t281 * t44) * t273) + m(7) * (t23 * t274 + (t278 * t30 - t281 * t31) * t273); -(t8 + t20 + t18) * t244 + (-t4 - t11 - t15) * t229 + (-t3 - t10 - t14) * t226 + m(7) * (t23 ^ 2 + t30 ^ 2 + t31 ^ 2) + m(6) * (t35 ^ 2 + t43 ^ 2 + t44 ^ 2) + m(5) * (t71 ^ 2 + t79 ^ 2 + t80 ^ 2); m(7) * (-t226 * t63 - t229 * t62) + m(6) * (-t226 * t83 - t229 * t82); m(7) * (-t21 * t244 - t226 * t24 - t229 * t25) + m(6) * (-t226 * t36 - t229 * t37 - t244 * t32); (-t244 * t274 + (t226 * t281 - t229 * t278) * t273) * t354; m(7) * (-t226 * t31 - t229 * t30 - t23 * t244) + m(6) * (-t226 * t44 - t229 * t43 - t244 * t35); (t226 ^ 2 + t229 ^ 2 + t244 ^ 2) * t354; m(7) * (t60 * t62 + t61 * t63) + t40 + t317 * t189 + t318 * t187; m(7) * (t21 * t42 + t24 * t61 + t25 * t60) + t6 * t345 + t274 * t347 + t5 * t346 + t9 * t344 + (t278 * t2 / 0.2e1 + t281 * t348) * t273; m(7) * (t42 * t274 + (t278 * t60 - t281 * t61) * t273); -t244 * t347 - t229 * t2 / 0.2e1 + t4 * t345 + t226 * t348 + m(7) * (t23 * t42 + t30 * t60 + t31 * t61) + t8 * t344 + t3 * t346; m(7) * (-t226 * t61 - t229 * t60 - t244 * t42); t189 * t2 + t187 * t1 + t232 * t7 + m(7) * (t42 ^ 2 + t60 ^ 2 + t61 ^ 2);];
%% Postprocessing: Reshape Output
% From vec2symmat_6_matlab.m
res = [t97(1) t97(2) t97(4) t97(7) t97(11) t97(16); t97(2) t97(3) t97(5) t97(8) t97(12) t97(17); t97(4) t97(5) t97(6) t97(9) t97(13) t97(18); t97(7) t97(8) t97(9) t97(10) t97(14) t97(19); t97(11) t97(12) t97(13) t97(14) t97(15) t97(20); t97(16) t97(17) t97(18) t97(19) t97(20) t97(21);];
Mq  = res;
