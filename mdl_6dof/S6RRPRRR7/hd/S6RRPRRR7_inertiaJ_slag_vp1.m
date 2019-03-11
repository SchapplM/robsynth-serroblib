% Calculate joint inertia matrix for
% S6RRPRRR7
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% pkin [10x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d2,d4,d5,d6]';
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
% Datum: 2019-03-09 14:01
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Mq = S6RRPRRR7_inertiaJ_slag_vp1(qJ, ...
  pkin, m, rSges, Icges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(10,1),zeros(7,1),zeros(7,3),zeros(7,6)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRPRRR7_inertiaJ_slag_vp1: qJ has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S6RRPRRR7_inertiaJ_slag_vp1: pkin has to be [10x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6RRPRRR7_inertiaJ_slag_vp1: m has to be [7x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [7,3]), ...
  'S6RRPRRR7_inertiaJ_slag_vp1: rSges has to be [7x3] (double)');
assert(isreal(Icges) && all(size(Icges) == [7 6]), ...
  'S6RRPRRR7_inertiaJ_slag_vp1: Icges has to be [7x6] (double)'); 

%% Symbolic Calculation
% From inertia_joint_joint_floatb_twist_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 13:55:46
% EndTime: 2019-03-09 13:55:58
% DurationCPUTime: 5.03s
% Computational Cost: add. (10658->474), mult. (21597->683), div. (0->0), fcn. (26116->10), ass. (0->241)
t349 = Icges(3,1) + Icges(4,1);
t347 = Icges(4,4) + Icges(3,5);
t215 = sin(qJ(2));
t348 = (Icges(3,4) - Icges(4,5)) * t215;
t346 = Icges(4,2) + Icges(3,3);
t214 = sin(qJ(4));
t218 = cos(qJ(2));
t326 = cos(qJ(4));
t266 = t215 * t326;
t180 = -t218 * t214 + t266;
t216 = sin(qJ(1));
t172 = t180 * t216;
t219 = cos(qJ(1));
t284 = t218 * t219;
t174 = t214 * t284 - t219 * t266;
t179 = t215 * t214 + t218 * t326;
t175 = t179 * t219;
t213 = sin(qJ(5));
t217 = cos(qJ(5));
t144 = -t175 * t213 - t216 * t217;
t288 = t213 * t216;
t145 = t175 * t217 - t288;
t173 = t179 * t216;
t142 = -t173 * t213 + t217 * t219;
t287 = t213 * t219;
t143 = t173 * t217 + t287;
t230 = Icges(6,5) * t143 + Icges(6,6) * t142 - Icges(6,3) * t172;
t293 = Icges(6,6) * t172;
t296 = Icges(6,2) * t142;
t304 = Icges(6,4) * t143;
t232 = -t293 + t296 + t304;
t299 = Icges(6,5) * t172;
t309 = Icges(6,1) * t143;
t235 = Icges(6,4) * t142 - t299 + t309;
t221 = t144 * t232 + t145 * t235 + t174 * t230;
t91 = t145 * Icges(6,5) + t144 * Icges(6,6) + t174 * Icges(6,3);
t92 = t145 * Icges(6,4) + t144 * Icges(6,2) + t174 * Icges(6,6);
t93 = t145 * Icges(6,1) + t144 * Icges(6,4) + t174 * Icges(6,5);
t41 = t144 * t92 + t145 * t93 + t174 * t91;
t110 = t179 * Icges(6,3) + (t217 * Icges(6,5) - t213 * Icges(6,6)) * t180;
t111 = t179 * Icges(6,6) + (t217 * Icges(6,4) - t213 * Icges(6,2)) * t180;
t112 = t179 * Icges(6,5) + (t217 * Icges(6,1) - t213 * Icges(6,4)) * t180;
t52 = t110 * t174 + t111 * t144 + t112 * t145;
t10 = -t172 * t221 + t41 * t174 + t52 * t179;
t333 = t172 ^ 2;
t224 = Icges(6,3) * t333 + (-0.2e1 * t299 + t309) * t143 + (-0.2e1 * t293 + t296 + 0.2e1 * t304) * t142;
t40 = t142 * t92 + t143 * t93 - t172 * t91;
t20 = t216 * t40 - t219 * t224;
t21 = t216 * t41 - t219 * t221;
t212 = qJ(5) + qJ(6);
t202 = sin(t212);
t203 = cos(t212);
t130 = -t173 * t202 + t203 * t219;
t131 = t173 * t203 + t202 * t219;
t292 = Icges(7,6) * t172;
t295 = Icges(7,2) * t130;
t298 = Icges(7,5) * t172;
t303 = Icges(7,4) * t131;
t308 = Icges(7,1) * t131;
t223 = Icges(7,3) * t333 + (-0.2e1 * t298 + t308) * t131 + (-0.2e1 * t292 + t295 + 0.2e1 * t303) * t130;
t132 = -t175 * t202 - t203 * t216;
t133 = t175 * t203 - t202 * t216;
t83 = t133 * Icges(7,5) + t132 * Icges(7,6) + t174 * Icges(7,3);
t84 = t133 * Icges(7,4) + t132 * Icges(7,2) + t174 * Icges(7,6);
t85 = t133 * Icges(7,1) + t132 * Icges(7,4) + t174 * Icges(7,5);
t33 = t130 * t84 + t131 * t85 - t172 * t83;
t16 = t216 * t33 - t219 * t223;
t229 = Icges(7,5) * t131 + Icges(7,6) * t130 - Icges(7,3) * t172;
t231 = -t292 + t295 + t303;
t234 = Icges(7,4) * t130 - t298 + t308;
t222 = t132 * t231 + t133 * t234 + t174 * t229;
t34 = t132 * t84 + t133 * t85 + t174 * t83;
t17 = t216 * t34 - t219 * t222;
t328 = t174 / 0.2e1;
t329 = -t172 / 0.2e1;
t337 = t216 / 0.2e1;
t339 = t179 / 0.2e1;
t341 = -t219 / 0.2e1;
t42 = t179 * t229 + (-t202 * t231 + t203 * t234) * t180;
t43 = t179 * t83 + (-t202 * t84 + t203 * t85) * t180;
t104 = t179 * Icges(7,3) + (t203 * Icges(7,5) - t202 * Icges(7,6)) * t180;
t105 = t179 * Icges(7,6) + (t203 * Icges(7,4) - t202 * Icges(7,2)) * t180;
t106 = t179 * Icges(7,5) + (t203 * Icges(7,1) - t202 * Icges(7,4)) * t180;
t49 = -t104 * t172 + t105 * t130 + t106 * t131;
t7 = -t172 * t223 + t33 * t174 + t49 * t179;
t50 = t104 * t174 + t105 * t132 + t106 * t133;
t8 = -t172 * t222 + t34 * t174 + t50 * t179;
t259 = t16 * t329 + t17 * t328 + (t43 * t216 - t42 * t219) * t339 + t8 * t337 + t7 * t341;
t44 = t179 * t230 + (-t213 * t232 + t217 * t235) * t180;
t45 = t179 * t91 + (-t213 * t92 + t217 * t93) * t180;
t51 = -t110 * t172 + t111 * t142 + t112 * t143;
t9 = -t172 * t224 + t40 * t174 + t51 * t179;
t345 = t20 * t329 + t21 * t328 + (t45 * t216 - t44 * t219) * t339 + t10 * t337 + t9 * t341 + t259;
t344 = t347 * t218 + (-Icges(3,6) + Icges(4,6)) * t215;
t343 = t218 * t349 - t348;
t342 = -t216 / 0.2e1;
t340 = t219 / 0.2e1;
t336 = -t344 * t216 + t346 * t219;
t335 = t346 * t216 + t344 * t219;
t115 = t175 * Icges(5,5) - t174 * Icges(5,6) - t216 * Icges(5,3);
t116 = t175 * Icges(5,4) - t174 * Icges(5,2) - t216 * Icges(5,6);
t117 = t175 * Icges(5,1) - t174 * Icges(5,4) - t216 * Icges(5,5);
t211 = t219 ^ 2;
t294 = Icges(5,6) * t219;
t297 = Icges(5,2) * t172;
t300 = Icges(5,5) * t219;
t305 = Icges(5,4) * t173;
t310 = Icges(5,1) * t173;
t272 = t16 + t20 + (t115 * t219 + t116 * t172 + t117 * t173) * t216 - (Icges(5,3) * t211 + (0.2e1 * t300 + t310) * t173 - (-0.2e1 * t294 - t297 - 0.2e1 * t305) * t172) * t219;
t233 = t294 + t297 + t305;
t236 = Icges(5,4) * t172 + t300 + t310;
t271 = t17 + t21 + (-t115 * t216 - t174 * t116 + t117 * t175) * t216 - (t175 * t236 - t174 * t233 - t216 * (Icges(5,5) * t173 + Icges(5,6) * t172 + Icges(5,3) * t219)) * t219;
t334 = t271 * t216 - t272 * t219;
t262 = t175 * pkin(4) + pkin(9) * t174;
t168 = t172 * pkin(9);
t263 = t173 * pkin(4) - t168;
t281 = t216 * t263 + t219 * t262;
t201 = pkin(5) * t217 + pkin(4);
t220 = -pkin(10) - pkin(9);
t280 = -t174 * t220 + t175 * t201;
t82 = -pkin(5) * t288 - t262 + t280;
t87 = t133 * rSges(7,1) + t132 * rSges(7,2) + t174 * rSges(7,3);
t319 = t82 + t87;
t321 = -pkin(4) + t201;
t81 = pkin(5) * t287 + t172 * t220 + t173 * t321 + t168;
t252 = -t131 * rSges(7,1) - t130 * rSges(7,2);
t86 = -t172 * rSges(7,3) - t252;
t320 = t81 + t86;
t25 = -t216 * t320 - t219 * t319 - t281;
t210 = t216 ^ 2;
t332 = m(4) / 0.2e1;
t331 = m(6) / 0.2e1;
t330 = m(7) / 0.2e1;
t327 = -rSges(5,3) - pkin(8);
t325 = pkin(8) * t216;
t323 = t219 * pkin(8);
t318 = rSges(4,2) * t219;
t314 = t219 * rSges(3,3);
t313 = t42 * t172;
t312 = t43 * t174;
t311 = t180 * t203 * t106 + t179 * t104;
t306 = Icges(3,4) * t218;
t301 = Icges(4,5) * t218;
t291 = qJ(3) * t215;
t290 = t105 * t202;
t289 = t111 * t213;
t286 = t215 * t219;
t283 = t180 * t217 * t112 + t179 * t110;
t102 = t321 * t180 + (-pkin(9) - t220) * t179;
t109 = rSges(7,3) * t179 + (rSges(7,1) * t203 - rSges(7,2) * t202) * t180;
t282 = t102 + t109;
t279 = t175 * rSges(5,1) - t174 * rSges(5,2);
t276 = pkin(2) * t284 + qJ(3) * t286;
t278 = t210 * (pkin(2) * t218 + t291) + t219 * t276;
t188 = pkin(2) * t215 - qJ(3) * t218;
t277 = -rSges(4,1) * t215 + rSges(4,3) * t218 - t188;
t275 = t219 * pkin(1) + t216 * pkin(7);
t274 = t210 + t211;
t53 = (-t180 * t290 + t311) * t179;
t273 = -t172 * t7 + t174 * t8 + t179 * (t312 + t53 - t313);
t270 = t44 / 0.2e1 + t51 / 0.2e1;
t269 = t52 / 0.2e1 + t45 / 0.2e1;
t95 = t145 * rSges(6,1) + t144 * rSges(6,2) + t174 * rSges(6,3);
t268 = rSges(4,1) * t284 + t216 * rSges(4,2) + rSges(4,3) * t286;
t267 = -pkin(5) * t213 - pkin(8);
t265 = t347 * t215 / 0.2e1 + (-Icges(4,6) / 0.2e1 + Icges(3,6) / 0.2e1) * t218;
t264 = -pkin(3) * t215 - t188;
t199 = pkin(3) * t284;
t261 = t216 * (t216 * t218 * pkin(3) + t323) + t219 * (t199 - t325) + t278;
t260 = t275 + t276;
t257 = -t313 / 0.2e1 + t312 / 0.2e1 + t49 * t329 + t50 * t328 + t53;
t113 = rSges(6,3) * t179 + (rSges(6,1) * t217 - rSges(6,2) * t213) * t180;
t256 = -t113 + t264;
t140 = rSges(5,1) * t180 - rSges(5,2) * t179;
t255 = -t140 + t264;
t254 = rSges(3,1) * t218 - rSges(3,2) * t215;
t253 = -rSges(5,1) * t173 - rSges(5,2) * t172;
t251 = t199 + t260;
t248 = -Icges(3,2) * t215 + t306;
t245 = Icges(4,3) * t215 + t301;
t107 = t255 * t216;
t108 = t255 * t219;
t244 = t107 * t216 + t108 * t219;
t80 = t216 * t253 - t219 * t279;
t239 = t264 - t282;
t238 = rSges(3,1) * t284 - rSges(3,2) * t286 + t216 * rSges(3,3);
t207 = t219 * pkin(7);
t226 = t207 + (-t291 - pkin(1) + (-pkin(2) - pkin(3)) * t218) * t216;
t96 = t219 * t327 + t226 + t253;
t97 = t216 * t327 + t251 + t279;
t237 = m(5) * (t216 * t97 + t219 * t96);
t94 = rSges(6,1) * t143 + rSges(6,2) * t142 - rSges(6,3) * t172;
t48 = -t216 * t94 - t219 * t95 - t281;
t135 = Icges(5,5) * t180 - Icges(5,6) * t179;
t136 = Icges(5,4) * t180 - Icges(5,2) * t179;
t137 = Icges(5,1) * t180 - Icges(5,4) * t179;
t228 = t42 / 0.2e1 + t49 / 0.2e1 + t135 * t340 + t136 * t172 / 0.2e1 + t137 * t173 / 0.2e1 - t179 * t233 / 0.2e1 + t180 * t236 / 0.2e1 + t270;
t227 = -t43 / 0.2e1 - t50 / 0.2e1 + t135 * t337 + t136 * t328 - t137 * t175 / 0.2e1 + t116 * t339 - t117 * t180 / 0.2e1 - t269;
t192 = rSges(2,1) * t219 - rSges(2,2) * t216;
t191 = -rSges(2,1) * t216 - rSges(2,2) * t219;
t190 = rSges(3,1) * t215 + rSges(3,2) * t218;
t149 = t277 * t219;
t148 = t277 * t216;
t147 = t238 + t275;
t146 = t314 + t207 + (-pkin(1) - t254) * t216;
t141 = t180 * pkin(4) + t179 * pkin(9);
t129 = t219 * t141;
t128 = t216 * t141;
t123 = t260 + t268;
t122 = t318 + t207 + (-pkin(1) + (-rSges(4,1) - pkin(2)) * t218 + (-rSges(4,3) - qJ(3)) * t215) * t216;
t114 = t219 * t238 + (t216 * t254 - t314) * t216;
t98 = t172 * t109;
t90 = t219 * t268 + (-t318 + (rSges(4,1) * t218 + rSges(4,3) * t215) * t216) * t216 + t278;
t89 = t113 * t219 + t129;
t88 = t113 * t216 + t128;
t79 = t179 * t87;
t78 = t219 * t256 - t129;
t77 = t216 * t256 - t128;
t76 = t174 * t86;
t71 = t251 + t262 + t95 - t325;
t70 = t226 - t263 - t94 - t323;
t69 = t219 * t282 + t129;
t68 = t216 * t282 + t128;
t67 = -t113 * t174 + t179 * t95;
t66 = -t113 * t172 - t179 * t94;
t65 = -t80 + t261;
t64 = -t109 * t174 + t79;
t63 = -t179 * t86 - t98;
t62 = t216 * t267 + t251 + t280 + t87;
t61 = -t173 * t201 + t267 * t219 - (-rSges(7,3) + t220) * t172 + t226 + t252;
t60 = t219 * t239 - t129;
t59 = t216 * t239 - t128;
t56 = t172 * t95 + t174 * t94;
t55 = (-t180 * t289 + t283) * t179;
t54 = t172 * t87 + t76;
t39 = -t48 + t261;
t28 = -t174 * t282 + t179 * t82 + t79;
t27 = -t102 * t172 - t179 * t320 - t98;
t26 = t172 * t319 + t174 * t81 + t76;
t24 = -t25 + t261;
t1 = [-t179 * t136 + Icges(2,3) + (t137 - t289 - t290) * t180 + m(7) * (t61 ^ 2 + t62 ^ 2) + m(6) * (t70 ^ 2 + t71 ^ 2) + m(5) * (t96 ^ 2 + t97 ^ 2) + m(4) * (t122 ^ 2 + t123 ^ 2) + m(3) * (t146 ^ 2 + t147 ^ 2) + m(2) * (t191 ^ 2 + t192 ^ 2) + t283 + t311 + ((Icges(3,2) + Icges(4,3)) * t218 + t348) * t218 + (t215 * t349 - t301 + t306) * t215; (t265 * t219 + (Icges(3,6) * t340 + Icges(4,6) * t341 + t245 * t337 + t248 * t342) * t218 + (t347 * t340 + t343 * t342) * t215 - t228) * t219 + (t265 * t216 + (Icges(3,6) * t337 + Icges(4,6) * t342 + t245 * t341 + t248 * t340) * t218 + (t347 * t337 + t343 * t340) * t215 - t227) * t216 + m(7) * (t59 * t62 + t60 * t61) + m(6) * (t70 * t78 + t71 * t77) + m(5) * (t107 * t97 + t108 * t96) + m(4) * (t122 * t149 + t123 * t148) + m(3) * (-t146 * t219 - t147 * t216) * t190; m(7) * (t24 ^ 2 + t59 ^ 2 + t60 ^ 2) + m(6) * (t39 ^ 2 + t77 ^ 2 + t78 ^ 2) + m(5) * (t107 ^ 2 + t108 ^ 2 + t65 ^ 2) + m(4) * (t148 ^ 2 + t149 ^ 2 + t90 ^ 2) + m(3) * (t190 ^ 2 * t274 + t114 ^ 2) + (t211 * t336 - t272) * t219 + (t335 * t210 + (t216 * t336 + t219 * t335) * t219 + t271) * t216; 0.2e1 * ((t216 * t62 + t219 * t61) * t330 + (t216 * t71 + t219 * t70) * t331 + t237 / 0.2e1 + (t122 * t219 + t123 * t216) * t332) * t215; m(7) * (-t218 * t24 + (t216 * t59 + t219 * t60) * t215) + m(6) * (-t218 * t39 + (t216 * t77 + t219 * t78) * t215) + m(5) * (t215 * t244 - t218 * t65) + m(4) * (-t218 * t90 + (t148 * t216 + t149 * t219) * t215); 0.2e1 * (t332 + m(5) / 0.2e1 + t331 + t330) * (t215 ^ 2 * t274 + t218 ^ 2); t228 * t219 + t227 * t216 + m(7) * (t61 * t69 + t62 * t68) + m(6) * (t70 * t89 + t71 * t88) + t140 * t237; m(7) * (t25 * t24 + t59 * t68 + t60 * t69) + m(6) * (t48 * t39 + t77 * t88 + t78 * t89) + m(5) * (t140 * t244 + t65 * t80) - t334; m(5) * (t140 * t215 * t274 - t218 * t80) + m(6) * (-t218 * t48 + (t216 * t88 + t219 * t89) * t215) + m(7) * (-t218 * t25 + (t216 * t68 + t219 * t69) * t215); m(7) * (t25 ^ 2 + t68 ^ 2 + t69 ^ 2) + m(6) * (t48 ^ 2 + t88 ^ 2 + t89 ^ 2) + m(5) * (t140 ^ 2 * t274 + t80 ^ 2) + t334; t55 + t269 * t174 - t270 * t172 + m(7) * (t27 * t61 + t28 * t62) + m(6) * (t66 * t70 + t67 * t71) + t257; m(7) * (t24 * t26 + t27 * t60 + t28 * t59) + m(6) * (t56 * t39 + t66 * t78 + t67 * t77) + t345; m(6) * (-t218 * t56 + (t216 * t67 + t219 * t66) * t215) + m(7) * (-t218 * t26 + (t216 * t28 + t219 * t27) * t215); m(7) * (t26 * t25 + t27 * t69 + t28 * t68) + m(6) * (t56 * t48 + t66 * t89 + t67 * t88) - t345; t174 * t10 - t172 * t9 + t179 * (-t44 * t172 + t45 * t174 + t55) + m(7) * (t26 ^ 2 + t27 ^ 2 + t28 ^ 2) + m(6) * (t56 ^ 2 + t66 ^ 2 + t67 ^ 2) + t273; m(7) * (t61 * t63 + t62 * t64) + t257; m(7) * (t24 * t54 + t59 * t64 + t60 * t63) + t259; m(7) * (-t218 * t54 + (t216 * t64 + t219 * t63) * t215); m(7) * (t54 * t25 + t63 * t69 + t64 * t68) - t259; m(7) * (t26 * t54 + t27 * t63 + t28 * t64) + t273; m(7) * (t54 ^ 2 + t63 ^ 2 + t64 ^ 2) + t273;];
%% Postprocessing: Reshape Output
% From vec2symmat_6_matlab.m
res = [t1(1) t1(2) t1(4) t1(7) t1(11) t1(16); t1(2) t1(3) t1(5) t1(8) t1(12) t1(17); t1(4) t1(5) t1(6) t1(9) t1(13) t1(18); t1(7) t1(8) t1(9) t1(10) t1(14) t1(19); t1(11) t1(12) t1(13) t1(14) t1(15) t1(20); t1(16) t1(17) t1(18) t1(19) t1(20) t1(21);];
Mq  = res;
