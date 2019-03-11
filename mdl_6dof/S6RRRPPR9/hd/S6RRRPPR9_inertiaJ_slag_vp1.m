% Calculate joint inertia matrix for
% S6RRRPPR9
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% pkin [11x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,alpha2,d1,d2,d3,d6,theta4]';
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
% Datum: 2019-03-09 16:20
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Mq = S6RRRPPR9_inertiaJ_slag_vp1(qJ, ...
  pkin, m, rSges, Icges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(11,1),zeros(7,1),zeros(7,3),zeros(7,6)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRRPPR9_inertiaJ_slag_vp1: qJ has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [11 1]), ...
  'S6RRRPPR9_inertiaJ_slag_vp1: pkin has to be [11x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6RRRPPR9_inertiaJ_slag_vp1: m has to be [7x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [7,3]), ...
  'S6RRRPPR9_inertiaJ_slag_vp1: rSges has to be [7x3] (double)');
assert(isreal(Icges) && all(size(Icges) == [7 6]), ...
  'S6RRRPPR9_inertiaJ_slag_vp1: Icges has to be [7x6] (double)'); 

%% Symbolic Calculation
% From inertia_joint_joint_floatb_twist_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 16:11:41
% EndTime: 2019-03-09 16:11:54
% DurationCPUTime: 5.03s
% Computational Cost: add. (20473->625), mult. (52749->875), div. (0->0), fcn. (68517->12), ass. (0->282)
t310 = m(6) / 0.2e1 + m(7) / 0.2e1;
t342 = 0.2e1 * t310;
t274 = cos(pkin(6));
t280 = cos(qJ(2));
t281 = cos(qJ(1));
t322 = t280 * t281;
t277 = sin(qJ(2));
t278 = sin(qJ(1));
t325 = t277 * t278;
t261 = -t274 * t325 + t322;
t276 = sin(qJ(3));
t273 = sin(pkin(6));
t328 = t273 * t278;
t335 = cos(qJ(3));
t241 = t261 * t335 + t276 * t328;
t323 = t278 * t280;
t324 = t277 * t281;
t260 = t274 * t323 + t324;
t272 = sin(pkin(11));
t332 = cos(pkin(11));
t209 = t241 * t332 + t260 * t272;
t208 = t241 * t272 - t260 * t332;
t275 = sin(qJ(6));
t279 = cos(qJ(6));
t152 = t208 * t279 - t209 * t275;
t153 = t208 * t275 + t209 * t279;
t302 = t273 * t335;
t240 = t261 * t276 - t278 * t302;
t97 = t153 * rSges(7,1) + t152 * rSges(7,2) - t240 * rSges(7,3);
t341 = t209 * pkin(5) + t97;
t326 = t273 * t281;
t329 = t273 * t277;
t256 = -t274 * t335 + t276 * t329;
t257 = t274 * t276 + t277 * t302;
t327 = t273 * t280;
t211 = Icges(4,5) * t257 - Icges(4,6) * t256 - Icges(4,3) * t327;
t212 = Icges(4,4) * t257 - Icges(4,2) * t256 - Icges(4,6) * t327;
t213 = Icges(4,1) * t257 - Icges(4,4) * t256 - Icges(4,5) * t327;
t110 = -t211 * t327 - t256 * t212 + t257 * t213;
t236 = t257 * t272 + t332 * t327;
t237 = t257 * t332 - t272 * t327;
t188 = t236 * t279 - t237 * t275;
t189 = t236 * t275 + t237 * t279;
t112 = Icges(7,5) * t189 + Icges(7,6) * t188 - Icges(7,3) * t256;
t113 = Icges(7,4) * t189 + Icges(7,2) * t188 - Icges(7,6) * t256;
t114 = Icges(7,1) * t189 + Icges(7,4) * t188 - Icges(7,5) * t256;
t46 = -t256 * t112 + t188 * t113 + t189 * t114;
t164 = Icges(6,5) * t237 + Icges(6,6) * t256 + Icges(6,3) * t236;
t166 = Icges(6,4) * t237 + Icges(6,2) * t256 + Icges(6,6) * t236;
t168 = Icges(6,1) * t237 + Icges(6,4) * t256 + Icges(6,5) * t236;
t82 = t236 * t164 + t256 * t166 + t237 * t168;
t165 = Icges(5,5) * t237 - Icges(5,6) * t236 + Icges(5,3) * t256;
t167 = Icges(5,4) * t237 - Icges(5,2) * t236 + Icges(5,6) * t256;
t169 = Icges(5,1) * t237 - Icges(5,4) * t236 + Icges(5,5) * t256;
t83 = t256 * t165 - t236 * t167 + t237 * t169;
t340 = -t110 - t46 - t82 - t83;
t259 = t274 * t324 + t323;
t238 = t259 * t276 + t281 * t302;
t239 = t259 * t335 - t276 * t326;
t258 = -t274 * t322 + t325;
t206 = t239 * t272 - t258 * t332;
t207 = t239 * t332 + t258 * t272;
t150 = t206 * t279 - t207 * t275;
t151 = t206 * t275 + t207 * t279;
t90 = Icges(7,5) * t151 + Icges(7,6) * t150 - Icges(7,3) * t238;
t92 = Icges(7,4) * t151 + Icges(7,2) * t150 - Icges(7,6) * t238;
t94 = Icges(7,1) * t151 + Icges(7,4) * t150 - Icges(7,5) * t238;
t26 = t152 * t92 + t153 * t94 - t240 * t90;
t91 = Icges(7,5) * t153 + Icges(7,6) * t152 - Icges(7,3) * t240;
t93 = Icges(7,4) * t153 + Icges(7,2) * t152 - Icges(7,6) * t240;
t95 = Icges(7,1) * t153 + Icges(7,4) * t152 - Icges(7,5) * t240;
t27 = t152 * t93 + t153 * t95 - t240 * t91;
t42 = -t112 * t240 + t113 * t152 + t114 * t153;
t2 = -t238 * t26 - t240 * t27 - t256 * t42;
t339 = t2 / 0.2e1;
t338 = -t238 / 0.2e1;
t337 = -t240 / 0.2e1;
t336 = -t256 / 0.2e1;
t294 = -t151 * rSges(7,1) - t150 * rSges(7,2);
t96 = -t238 * rSges(7,3) - t294;
t334 = t207 * pkin(5) - t238 * pkin(10) + t96;
t333 = -pkin(10) * t240 + t341;
t331 = qJ(4) * t240;
t215 = Icges(3,5) * t259 - Icges(3,6) * t258 - Icges(3,3) * t326;
t330 = t215 * t281;
t115 = rSges(7,1) * t189 + rSges(7,2) * t188 - rSges(7,3) * t256;
t321 = pkin(5) * t237 - pkin(10) * t256 + t115;
t135 = t209 * rSges(5,1) - t208 * rSges(5,2) + t240 * rSges(5,3);
t235 = t241 * pkin(3);
t193 = t235 + t331;
t320 = -t135 - t193;
t198 = t206 * qJ(5);
t154 = t207 * pkin(4) + t198;
t192 = t239 * pkin(3) + t238 * qJ(4);
t181 = t260 * t192;
t319 = t260 * t154 + t181;
t155 = t209 * pkin(4) + t208 * qJ(5);
t318 = -t155 - t193;
t171 = rSges(5,1) * t237 - rSges(5,2) * t236 + rSges(5,3) * t256;
t226 = pkin(3) * t257 + qJ(4) * t256;
t317 = -t171 - t226;
t316 = t192 * t327 + t258 * t226;
t228 = t261 * pkin(2) + t260 * pkin(9);
t225 = t274 * t228;
t315 = t274 * t193 + t225;
t191 = pkin(4) * t237 + qJ(5) * t236;
t314 = -t191 - t226;
t227 = t259 * pkin(2) + t258 * pkin(9);
t313 = -t192 - t227;
t312 = t227 * t328 + t228 * t326;
t311 = t281 * pkin(1) + pkin(8) * t328;
t28 = t188 * t92 + t189 * t94 - t256 * t90;
t41 = -t112 * t238 + t113 * t150 + t114 * t151;
t309 = -t28 / 0.2e1 - t41 / 0.2e1;
t29 = t188 * t93 + t189 * t95 - t256 * t91;
t308 = -t42 / 0.2e1 - t29 / 0.2e1;
t134 = t209 * rSges(6,1) + t240 * rSges(6,2) + t208 * rSges(6,3);
t307 = -t134 + t318;
t306 = t274 * t155 + t315;
t305 = -t154 + t313;
t170 = rSges(6,1) * t237 + rSges(6,2) * t256 + rSges(6,3) * t236;
t304 = -t170 + t314;
t180 = t241 * rSges(4,1) - t240 * rSges(4,2) + t260 * rSges(4,3);
t245 = Icges(3,3) * t274 + (Icges(3,5) * t277 + Icges(3,6) * t280) * t273;
t246 = Icges(3,6) * t274 + (Icges(3,4) * t277 + Icges(3,2) * t280) * t273;
t247 = Icges(3,5) * t274 + (Icges(3,1) * t277 + Icges(3,4) * t280) * t273;
t303 = t274 * t245 + t246 * t327 + t247 * t329;
t222 = t261 * rSges(3,1) - t260 * rSges(3,2) + rSges(3,3) * t328;
t301 = -t278 * pkin(1) + pkin(8) * t326;
t214 = rSges(4,1) * t257 - rSges(4,2) * t256 - rSges(4,3) * t327;
t262 = (pkin(2) * t277 - pkin(9) * t280) * t273;
t300 = t273 * (-t214 - t262);
t299 = t318 - t333;
t298 = t314 - t321;
t297 = t154 * t327 + t258 * t191 + t316;
t296 = t192 * t328 + t193 * t326 + t312;
t295 = t273 * (-t262 + t317);
t293 = -t238 * rSges(6,2) - t206 * rSges(6,3);
t292 = t228 + t311;
t291 = t273 * (-t262 + t304);
t290 = t154 * t328 + t155 * t326 + t296;
t289 = t235 + t292;
t288 = t273 * (-t262 + t298);
t287 = -t227 + t301;
t179 = rSges(4,1) * t239 - rSges(4,2) * t238 + rSges(4,3) * t258;
t133 = rSges(5,1) * t207 - rSges(5,2) * t206 + rSges(5,3) * t238;
t221 = t259 * rSges(3,1) - t258 * rSges(3,2) - rSges(3,3) * t326;
t286 = -t192 + t287;
t285 = t155 + t289;
t284 = -t198 + t286;
t104 = t211 * t258 - t212 * t238 + t213 * t239;
t120 = Icges(6,5) * t207 + Icges(6,6) * t238 + Icges(6,3) * t206;
t124 = Icges(6,4) * t207 + Icges(6,2) * t238 + Icges(6,6) * t206;
t128 = Icges(6,1) * t207 + Icges(6,4) * t238 + Icges(6,5) * t206;
t63 = t120 * t236 + t124 * t256 + t128 * t237;
t122 = Icges(5,5) * t207 - Icges(5,6) * t206 + Icges(5,3) * t238;
t126 = Icges(5,4) * t207 - Icges(5,2) * t206 + Icges(5,6) * t238;
t130 = Icges(5,1) * t207 - Icges(5,4) * t206 + Icges(5,5) * t238;
t65 = t122 * t256 - t126 * t236 + t130 * t237;
t72 = t164 * t206 + t166 * t238 + t168 * t207;
t73 = t165 * t238 - t167 * t206 + t169 * t207;
t173 = Icges(4,5) * t239 - Icges(4,6) * t238 + Icges(4,3) * t258;
t175 = Icges(4,4) * t239 - Icges(4,2) * t238 + Icges(4,6) * t258;
t177 = Icges(4,1) * t239 - Icges(4,4) * t238 + Icges(4,5) * t258;
t88 = -t173 * t327 - t175 * t256 + t177 * t257;
t283 = t88 / 0.2e1 + t65 / 0.2e1 + t63 / 0.2e1 + t104 / 0.2e1 + t73 / 0.2e1 + t72 / 0.2e1 - t309;
t105 = t211 * t260 - t212 * t240 + t213 * t241;
t121 = Icges(6,5) * t209 + Icges(6,6) * t240 + Icges(6,3) * t208;
t125 = Icges(6,4) * t209 + Icges(6,2) * t240 + Icges(6,6) * t208;
t129 = Icges(6,1) * t209 + Icges(6,4) * t240 + Icges(6,5) * t208;
t64 = t121 * t236 + t125 * t256 + t129 * t237;
t123 = Icges(5,5) * t209 - Icges(5,6) * t208 + Icges(5,3) * t240;
t127 = Icges(5,4) * t209 - Icges(5,2) * t208 + Icges(5,6) * t240;
t131 = Icges(5,1) * t209 - Icges(5,4) * t208 + Icges(5,5) * t240;
t66 = t123 * t256 - t127 * t236 + t131 * t237;
t74 = t164 * t208 + t166 * t240 + t168 * t209;
t75 = t165 * t240 - t167 * t208 + t169 * t209;
t174 = Icges(4,5) * t241 - Icges(4,6) * t240 + Icges(4,3) * t260;
t176 = Icges(4,4) * t241 - Icges(4,2) * t240 + Icges(4,6) * t260;
t178 = Icges(4,1) * t241 - Icges(4,4) * t240 + Icges(4,5) * t260;
t89 = -t174 * t327 - t176 * t256 + t178 * t257;
t282 = t74 / 0.2e1 + t64 / 0.2e1 + t105 / 0.2e1 + t75 / 0.2e1 + t89 / 0.2e1 + t66 / 0.2e1 - t308;
t264 = rSges(2,1) * t281 - t278 * rSges(2,2);
t263 = -t278 * rSges(2,1) - rSges(2,2) * t281;
t248 = t274 * rSges(3,3) + (rSges(3,1) * t277 + rSges(3,2) * t280) * t273;
t220 = Icges(3,1) * t261 - Icges(3,4) * t260 + Icges(3,5) * t328;
t219 = Icges(3,1) * t259 - Icges(3,4) * t258 - Icges(3,5) * t326;
t218 = Icges(3,4) * t261 - Icges(3,2) * t260 + Icges(3,6) * t328;
t217 = Icges(3,4) * t259 - Icges(3,2) * t258 - Icges(3,6) * t326;
t216 = Icges(3,5) * t261 - Icges(3,6) * t260 + Icges(3,3) * t328;
t197 = t222 + t311;
t196 = -t221 + t301;
t184 = -t274 * t221 - t248 * t326;
t183 = t222 * t274 - t248 * t328;
t163 = t303 * t274;
t158 = (t221 * t278 + t222 * t281) * t273;
t157 = t245 * t328 - t246 * t260 + t247 * t261;
t156 = -t245 * t326 - t258 * t246 + t259 * t247;
t139 = t292 + t180;
t138 = -t179 + t287;
t132 = rSges(6,1) * t207 - t293;
t119 = -t180 * t327 - t214 * t260;
t118 = t179 * t327 + t214 * t258;
t117 = t274 * t216 + (t218 * t280 + t220 * t277) * t273;
t116 = t274 * t215 + (t217 * t280 + t219 * t277) * t273;
t109 = t110 * t274;
t108 = t179 * t260 - t180 * t258;
t107 = (-t179 - t227) * t274 + t281 * t300;
t106 = t274 * t180 + t278 * t300 + t225;
t100 = t135 + t289 + t331;
t99 = -t133 + t286;
t98 = (t179 * t278 + t180 * t281) * t273 + t312;
t87 = t174 * t260 - t176 * t240 + t178 * t241;
t86 = t173 * t260 - t175 * t240 + t177 * t241;
t85 = t174 * t258 - t176 * t238 + t178 * t239;
t84 = t173 * t258 - t175 * t238 + t177 * t239;
t81 = t83 * t274;
t80 = t82 * t274;
t79 = t134 + t285 + t331;
t78 = (-rSges(6,1) - pkin(4)) * t207 + t284 + t293;
t77 = t317 * t260 + t320 * t327;
t76 = t133 * t327 + t171 * t258 + t316;
t71 = (-t133 + t313) * t274 + t281 * t295;
t70 = t274 * t135 + t278 * t295 + t315;
t69 = t115 * t240 - t256 * t97;
t68 = -t115 * t238 + t256 * t96;
t67 = t260 * t133 + t320 * t258 + t181;
t62 = (t133 * t278 + t135 * t281) * t273 + t296;
t61 = t123 * t240 - t127 * t208 + t131 * t209;
t60 = t122 * t240 - t126 * t208 + t130 * t209;
t59 = t121 * t208 + t125 * t240 + t129 * t209;
t58 = t120 * t208 + t124 * t240 + t128 * t209;
t57 = t123 * t238 - t127 * t206 + t131 * t207;
t56 = t122 * t238 - t126 * t206 + t130 * t207;
t55 = t121 * t206 + t125 * t238 + t129 * t207;
t54 = t120 * t206 + t124 * t238 + t128 * t207;
t53 = (-pkin(10) + qJ(4)) * t240 + t285 + t341;
t52 = (rSges(7,3) + pkin(10)) * t238 + (-pkin(4) - pkin(5)) * t207 + t284 + t294;
t51 = t304 * t260 + t307 * t327;
t50 = t132 * t327 + t170 * t258 + t297;
t49 = t238 * t97 - t240 * t96;
t48 = (-t132 + t305) * t274 + t281 * t291;
t47 = t274 * t134 + t278 * t291 + t306;
t45 = t46 * t274;
t44 = t46 * t256;
t43 = t260 * t132 + t307 * t258 + t319;
t40 = (t132 * t278 + t134 * t281) * t273 + t290;
t39 = t109 + (t89 * t278 - t88 * t281) * t273;
t38 = -t110 * t327 + t88 * t258 + t89 * t260;
t37 = t105 * t274 + (t278 * t87 - t281 * t86) * t273;
t36 = t104 * t274 + (t278 * t85 - t281 * t84) * t273;
t35 = -t105 * t327 + t258 * t86 + t260 * t87;
t34 = -t104 * t327 + t258 * t84 + t260 * t85;
t33 = t298 * t260 + t299 * t327;
t32 = t321 * t258 + t334 * t327 + t297;
t31 = (t305 - t334) * t274 + t281 * t288;
t30 = t333 * t274 + t278 * t288 + t306;
t25 = t150 * t93 + t151 * t95 - t238 * t91;
t24 = t150 * t92 + t151 * t94 - t238 * t90;
t23 = t299 * t258 + t334 * t260 + t319;
t22 = (t334 * t278 + t333 * t281) * t273 + t290;
t21 = t81 + (t66 * t278 - t65 * t281) * t273;
t20 = t80 + (t64 * t278 - t63 * t281) * t273;
t19 = t65 * t258 + t66 * t260 - t83 * t327;
t18 = t63 * t258 + t64 * t260 - t82 * t327;
t17 = t75 * t274 + (t278 * t61 - t281 * t60) * t273;
t16 = t74 * t274 + (t278 * t59 - t281 * t58) * t273;
t15 = t73 * t274 + (t278 * t57 - t281 * t56) * t273;
t14 = t72 * t274 + (t278 * t55 - t281 * t54) * t273;
t13 = t258 * t60 + t260 * t61 - t75 * t327;
t12 = t258 * t58 + t260 * t59 - t74 * t327;
t11 = t258 * t56 + t260 * t57 - t73 * t327;
t10 = t258 * t54 + t260 * t55 - t72 * t327;
t9 = t45 + (t29 * t278 - t28 * t281) * t273;
t8 = t28 * t258 + t29 * t260 - t46 * t327;
t7 = -t28 * t238 - t29 * t240 - t44;
t6 = t42 * t274 + (-t26 * t281 + t27 * t278) * t273;
t5 = t41 * t274 + (-t24 * t281 + t25 * t278) * t273;
t4 = t258 * t26 + t260 * t27 - t42 * t327;
t3 = t24 * t258 + t25 * t260 - t41 * t327;
t1 = -t238 * t24 - t240 * t25 - t256 * t41;
t101 = [m(4) * (t138 ^ 2 + t139 ^ 2) + m(3) * (t196 ^ 2 + t197 ^ 2) + m(2) * (t263 ^ 2 + t264 ^ 2) + m(7) * (t52 ^ 2 + t53 ^ 2) + m(6) * (t78 ^ 2 + t79 ^ 2) + m(5) * (t100 ^ 2 + t99 ^ 2) + Icges(2,3) + t303 - t340; t81 + t109 + t80 + t45 + t163 + m(3) * (t183 * t197 + t184 * t196) + m(4) * (t106 * t139 + t107 * t138) + m(7) * (t30 * t53 + t31 * t52) + m(6) * (t47 * t79 + t48 * t78) + m(5) * (t100 * t70 + t71 * t99) + ((-t116 / 0.2e1 - t156 / 0.2e1 - t283) * t281 + (t117 / 0.2e1 + t157 / 0.2e1 + t282) * t278) * t273; (t9 + t21 + t20 + t39 + t163) * t274 + m(7) * (t22 ^ 2 + t30 ^ 2 + t31 ^ 2) + m(6) * (t40 ^ 2 + t47 ^ 2 + t48 ^ 2) + m(5) * (t62 ^ 2 + t70 ^ 2 + t71 ^ 2) + m(4) * (t106 ^ 2 + t107 ^ 2 + t98 ^ 2) + m(3) * (t158 ^ 2 + t183 ^ 2 + t184 ^ 2) + ((-t5 - t36 - t15 - t14 + (-t258 * t217 + t259 * t219 - t330 * t273) * t326) * t281 + (t6 + t37 + t17 + t16 + ((-t218 * t260 + t220 * t261 + (t216 * t278 - t330) * t273) * t278 + (t216 * t326 + t217 * t260 + t258 * t218 - t219 * t261 - t259 * t220) * t281) * t273) * t278 + ((-t116 - t156) * t281 + (t117 + t157) * t278) * t274) * t273; t340 * t327 + m(7) * (t32 * t52 + t33 * t53) + m(6) * (t50 * t78 + t51 * t79) + m(5) * (t100 * t77 + t76 * t99) + m(4) * (t118 * t138 + t119 * t139) + t282 * t260 + t283 * t258; (t8 / 0.2e1 + t18 / 0.2e1 + t19 / 0.2e1 + t38 / 0.2e1) * t274 + (t6 / 0.2e1 + t16 / 0.2e1 + t17 / 0.2e1 + t37 / 0.2e1) * t260 + (t5 / 0.2e1 + t14 / 0.2e1 + t15 / 0.2e1 + t36 / 0.2e1) * t258 + m(7) * (t22 * t23 + t30 * t33 + t31 * t32) + m(6) * (t40 * t43 + t47 * t51 + t48 * t50) + m(5) * (t62 * t67 + t70 * t77 + t71 * t76) + m(4) * (t106 * t119 + t107 * t118 + t108 * t98) + ((-t3 / 0.2e1 - t10 / 0.2e1 - t11 / 0.2e1 - t34 / 0.2e1) * t281 + (-t9 / 0.2e1 - t20 / 0.2e1 - t21 / 0.2e1 - t39 / 0.2e1) * t280 + (t4 / 0.2e1 + t12 / 0.2e1 + t13 / 0.2e1 + t35 / 0.2e1) * t278) * t273; (-t18 - t19 - t38 - t8) * t327 + (t4 + t12 + t13 + t35) * t260 + (t3 + t11 + t10 + t34) * t258 + m(7) * (t23 ^ 2 + t32 ^ 2 + t33 ^ 2) + m(6) * (t43 ^ 2 + t50 ^ 2 + t51 ^ 2) + m(5) * (t67 ^ 2 + t76 ^ 2 + t77 ^ 2) + m(4) * (t108 ^ 2 + t118 ^ 2 + t119 ^ 2); m(7) * (t238 * t53 + t240 * t52) + m(6) * (t238 * t79 + t240 * t78) + m(5) * (t100 * t238 + t240 * t99); m(7) * (t22 * t256 + t238 * t30 + t240 * t31) + m(6) * (t238 * t47 + t240 * t48 + t256 * t40) + m(5) * (t238 * t70 + t240 * t71 + t256 * t62); m(7) * (t23 * t256 + t238 * t33 + t240 * t32) + m(6) * (t238 * t51 + t240 * t50 + t256 * t43) + m(5) * (t238 * t77 + t240 * t76 + t256 * t67); 0.2e1 * (m(5) / 0.2e1 + t310) * (t238 ^ 2 + t240 ^ 2 + t256 ^ 2); m(7) * (t206 * t53 + t208 * t52) + m(6) * (t206 * t79 + t208 * t78); m(7) * (t206 * t30 + t208 * t31 + t22 * t236) + m(6) * (t206 * t47 + t208 * t48 + t236 * t40); m(7) * (t206 * t33 + t208 * t32 + t23 * t236) + m(6) * (t206 * t51 + t208 * t50 + t236 * t43); (t206 * t238 + t208 * t240 + t236 * t256) * t342; (t206 ^ 2 + t208 ^ 2 + t236 ^ 2) * t342; m(7) * (t52 * t68 + t53 * t69) - t44 + t308 * t240 + t309 * t238; m(7) * (t49 * t22 + t30 * t69 + t31 * t68) + t9 * t336 + t6 * t337 + t5 * t338 + t274 * t7 / 0.2e1 + (t278 * t339 - t281 * t1 / 0.2e1) * t273; -t7 * t327 / 0.2e1 + m(7) * (t49 * t23 + t32 * t68 + t33 * t69) + t260 * t339 + t258 * t1 / 0.2e1 + t3 * t338 + t4 * t337 + t8 * t336; m(7) * (t238 * t69 + t240 * t68 + t256 * t49); m(7) * (t206 * t69 + t208 * t68 + t236 * t49); -t240 * t2 - t238 * t1 - t256 * t7 + m(7) * (t49 ^ 2 + t68 ^ 2 + t69 ^ 2);];
%% Postprocessing: Reshape Output
% From vec2symmat_6_matlab.m
res = [t101(1) t101(2) t101(4) t101(7) t101(11) t101(16); t101(2) t101(3) t101(5) t101(8) t101(12) t101(17); t101(4) t101(5) t101(6) t101(9) t101(13) t101(18); t101(7) t101(8) t101(9) t101(10) t101(14) t101(19); t101(11) t101(12) t101(13) t101(14) t101(15) t101(20); t101(16) t101(17) t101(18) t101(19) t101(20) t101(21);];
Mq  = res;
