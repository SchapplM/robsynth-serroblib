% Calculate joint inertia matrix for
% S6RRPRPR3
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% pkin [11x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d2,d4,d6,theta3,theta5]';
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
% Datum: 2019-03-09 10:20
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Mq = S6RRPRPR3_inertiaJ_slag_vp1(qJ, ...
  pkin, m, rSges, Icges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(11,1),zeros(7,1),zeros(7,3),zeros(7,6)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRPRPR3_inertiaJ_slag_vp1: qJ has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [11 1]), ...
  'S6RRPRPR3_inertiaJ_slag_vp1: pkin has to be [11x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6RRPRPR3_inertiaJ_slag_vp1: m has to be [7x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [7,3]), ...
  'S6RRPRPR3_inertiaJ_slag_vp1: rSges has to be [7x3] (double)');
assert(isreal(Icges) && all(size(Icges) == [7 6]), ...
  'S6RRPRPR3_inertiaJ_slag_vp1: Icges has to be [7x6] (double)'); 

%% Symbolic Calculation
% From inertia_joint_joint_floatb_twist_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 10:16:02
% EndTime: 2019-03-09 10:16:12
% DurationCPUTime: 4.53s
% Computational Cost: add. (11439->494), mult. (9848->703), div. (0->0), fcn. (10477->12), ass. (0->245)
t329 = Icges(3,3) + Icges(4,3);
t220 = qJ(2) + pkin(10);
t210 = sin(t220);
t212 = cos(t220);
t226 = sin(qJ(2));
t229 = cos(qJ(2));
t328 = Icges(3,5) * t229 + Icges(4,5) * t212 - Icges(3,6) * t226 - Icges(4,6) * t210;
t219 = qJ(4) + pkin(11);
t209 = sin(t219);
t211 = cos(t219);
t128 = -Icges(6,3) * t212 + (Icges(6,5) * t211 - Icges(6,6) * t209) * t210;
t225 = sin(qJ(4));
t228 = cos(qJ(4));
t136 = -Icges(5,3) * t212 + (Icges(5,5) * t228 - Icges(5,6) * t225) * t210;
t327 = -t128 - t136;
t129 = -Icges(6,6) * t212 + (Icges(6,4) * t211 - Icges(6,2) * t209) * t210;
t137 = -Icges(5,6) * t212 + (Icges(5,4) * t228 - Icges(5,2) * t225) * t210;
t326 = -t209 * t129 - t225 * t137;
t325 = t210 / 0.2e1;
t324 = t212 / 0.2e1;
t323 = t226 / 0.2e1;
t322 = t229 / 0.2e1;
t227 = sin(qJ(1));
t230 = cos(qJ(1));
t280 = t227 * t209;
t159 = -t211 * t230 - t212 * t280;
t279 = t227 * t211;
t160 = -t209 * t230 + t212 * t279;
t288 = t210 * t227;
t91 = Icges(6,5) * t160 + Icges(6,6) * t159 + Icges(6,3) * t288;
t93 = Icges(6,4) * t160 + Icges(6,2) * t159 + Icges(6,6) * t288;
t95 = Icges(6,1) * t160 + Icges(6,4) * t159 + Icges(6,5) * t288;
t31 = t159 * t93 + t160 * t95 + t91 * t288;
t286 = t212 * t230;
t161 = -t209 * t286 + t279;
t162 = t211 * t286 + t280;
t287 = t210 * t230;
t92 = Icges(6,5) * t162 + Icges(6,6) * t161 + Icges(6,3) * t287;
t94 = Icges(6,4) * t162 + Icges(6,2) * t161 + Icges(6,6) * t287;
t96 = Icges(6,1) * t162 + Icges(6,4) * t161 + Icges(6,5) * t287;
t32 = t159 * t94 + t160 * t96 + t92 * t288;
t276 = t228 * t230;
t278 = t227 * t225;
t172 = -t212 * t278 - t276;
t277 = t227 * t228;
t283 = t225 * t230;
t173 = t212 * t277 - t283;
t106 = Icges(5,5) * t173 + Icges(5,6) * t172 + Icges(5,3) * t288;
t108 = Icges(5,4) * t173 + Icges(5,2) * t172 + Icges(5,6) * t288;
t110 = Icges(5,1) * t173 + Icges(5,4) * t172 + Icges(5,5) * t288;
t44 = t106 * t288 + t108 * t172 + t110 * t173;
t174 = -t212 * t283 + t277;
t175 = t212 * t276 + t278;
t107 = Icges(5,5) * t175 + Icges(5,6) * t174 + Icges(5,3) * t287;
t109 = Icges(5,4) * t175 + Icges(5,2) * t174 + Icges(5,6) * t287;
t111 = Icges(5,1) * t175 + Icges(5,4) * t174 + Icges(5,5) * t287;
t45 = t107 * t288 + t109 * t172 + t111 * t173;
t130 = -Icges(6,5) * t212 + (Icges(6,1) * t211 - Icges(6,4) * t209) * t210;
t56 = t128 * t288 + t129 * t159 + t130 * t160;
t138 = -Icges(5,5) * t212 + (Icges(5,1) * t228 - Icges(5,4) * t225) * t210;
t62 = t136 * t288 + t137 * t172 + t138 * t173;
t321 = (-t62 - t56) * t212 + ((t32 + t45) * t230 + (t31 + t44) * t227) * t210;
t33 = t161 * t93 + t162 * t95 + t91 * t287;
t34 = t161 * t94 + t162 * t96 + t92 * t287;
t46 = t106 * t287 + t174 * t108 + t175 * t110;
t47 = t107 * t287 + t174 * t109 + t175 * t111;
t57 = t128 * t287 + t161 * t129 + t162 * t130;
t63 = t136 * t287 + t174 * t137 + t175 * t138;
t320 = (-t63 - t57) * t212 + ((t34 + t47) * t230 + (t33 + t46) * t227) * t210;
t40 = -t212 * t91 + (-t209 * t93 + t211 * t95) * t210;
t48 = -t106 * t212 + (-t108 * t225 + t110 * t228) * t210;
t319 = -t40 - t48;
t41 = -t212 * t92 + (-t209 * t94 + t211 * t96) * t210;
t49 = -t107 * t212 + (-t109 * t225 + t111 * t228) * t210;
t318 = t41 + t49;
t206 = t228 * pkin(4) + pkin(3);
t183 = pkin(5) * t211 + t206;
t184 = pkin(4) * t225 + pkin(5) * t209;
t213 = qJ(6) + t219;
t204 = sin(t213);
t205 = cos(t213);
t281 = t227 * t205;
t146 = -t204 * t286 + t281;
t282 = t227 * t204;
t147 = t205 * t286 + t282;
t89 = t147 * rSges(7,1) + t146 * rSges(7,2) + rSges(7,3) * t287;
t317 = t183 * t286 + t227 * t184 + t89;
t223 = -qJ(5) - pkin(8);
t218 = -pkin(9) + t223;
t268 = t218 - t223;
t272 = t183 - t206;
t114 = t272 * t210 + t268 * t212;
t127 = -rSges(7,3) * t212 + (rSges(7,1) * t205 - rSges(7,2) * t204) * t210;
t316 = -t114 - t127;
t315 = (t130 * t211 + t138 * t228) * t210;
t314 = -t328 * t227 + t329 * t230;
t313 = t329 * t227 + t328 * t230;
t221 = t227 ^ 2;
t222 = t230 ^ 2;
t312 = m(6) / 0.2e1;
t311 = m(7) / 0.2e1;
t144 = -t205 * t230 - t212 * t282;
t145 = -t204 * t230 + t212 * t281;
t82 = Icges(7,5) * t145 + Icges(7,6) * t144 + Icges(7,3) * t288;
t84 = Icges(7,4) * t145 + Icges(7,2) * t144 + Icges(7,6) * t288;
t86 = Icges(7,1) * t145 + Icges(7,4) * t144 + Icges(7,5) * t288;
t27 = t144 * t84 + t145 * t86 + t82 * t288;
t83 = Icges(7,5) * t147 + Icges(7,6) * t146 + Icges(7,3) * t287;
t85 = Icges(7,4) * t147 + Icges(7,2) * t146 + Icges(7,6) * t287;
t87 = Icges(7,1) * t147 + Icges(7,4) * t146 + Icges(7,5) * t287;
t28 = t144 * t85 + t145 * t87 + t83 * t288;
t124 = -Icges(7,3) * t212 + (Icges(7,5) * t205 - Icges(7,6) * t204) * t210;
t125 = -Icges(7,6) * t212 + (Icges(7,4) * t205 - Icges(7,2) * t204) * t210;
t126 = -Icges(7,5) * t212 + (Icges(7,1) * t205 - Icges(7,4) * t204) * t210;
t52 = t124 * t288 + t125 * t144 + t126 * t145;
t5 = -t52 * t212 + (t227 * t27 + t230 * t28) * t210;
t29 = t146 * t84 + t147 * t86 + t82 * t287;
t30 = t146 * t85 + t147 * t87 + t83 * t287;
t53 = t124 * t287 + t146 * t125 + t147 * t126;
t6 = -t53 * t212 + (t227 * t29 + t230 * t30) * t210;
t310 = t6 * t287 + t5 * t288;
t309 = -t212 / 0.2e1;
t308 = t227 / 0.2e1;
t307 = -t230 / 0.2e1;
t306 = pkin(2) * t226;
t305 = pkin(3) * t212;
t304 = pkin(8) * t210;
t303 = -pkin(3) + t206;
t302 = t326 * t210 + t327 * t212 + t315;
t273 = -pkin(4) * t278 - t206 * t286;
t301 = -t268 * t287 + t273 + t317;
t300 = rSges(3,1) * t229;
t299 = rSges(3,2) * t226;
t298 = t230 * rSges(3,3);
t116 = t210 * t205 * t126;
t290 = t204 * t125;
t59 = -t212 * t124 - t210 * t290 + t116;
t297 = t59 * t212;
t235 = -t223 * t287 - t273;
t270 = pkin(3) * t286 + pkin(8) * t287;
t105 = t235 - t270;
t98 = t162 * rSges(6,1) + t161 * rSges(6,2) + rSges(6,3) * t287;
t296 = -t105 - t98;
t271 = pkin(4) * t283 + t223 * t288;
t104 = (t303 * t212 - t304) * t227 - t271;
t123 = (pkin(8) + t223) * t212 + t303 * t210;
t295 = t212 * t104 + t123 * t288;
t248 = -t145 * rSges(7,1) - t144 * rSges(7,2);
t88 = rSges(7,3) * t288 - t248;
t67 = t127 * t288 + t212 * t88;
t294 = Icges(3,4) * t226;
t293 = Icges(3,4) * t229;
t292 = Icges(4,4) * t210;
t291 = Icges(4,4) * t212;
t224 = -qJ(3) - pkin(7);
t285 = t224 * t230;
t207 = pkin(2) * t229 + pkin(1);
t199 = t230 * t207;
t217 = t230 * pkin(7);
t275 = t227 * (t285 + t217 + (-pkin(1) + t207) * t227) + t230 * (-pkin(1) * t230 + t199 + (-pkin(7) - t224) * t227);
t269 = t227 * rSges(3,3) + t230 * t300;
t267 = t221 + t222;
t266 = t312 + t311;
t265 = -t105 - t301;
t113 = t175 * rSges(5,1) + t174 * rSges(5,2) + rSges(5,3) * t287;
t264 = t288 / 0.2e1;
t263 = t287 / 0.2e1;
t262 = Icges(3,5) * t323 + Icges(4,5) * t325 + Icges(3,6) * t322 + Icges(4,6) * t324;
t261 = -rSges(4,1) * t210 - rSges(4,2) * t212 - t306;
t260 = -pkin(3) * t210 + pkin(8) * t212 - t306;
t37 = -t212 * t82 + (-t204 * t84 + t205 * t86) * t210;
t38 = -t212 * t83 + (-t204 * t85 + t205 * t87) * t210;
t259 = (t37 + t52) * t264 + (t38 + t53) * t263;
t258 = -t227 * t224 + t199;
t9 = -t297 + (t227 * t37 + t230 * t38) * t210;
t257 = -t212 * t9 + t310;
t256 = t221 * (t304 + t305) + t230 * t270 + t275;
t15 = t28 * t227 - t230 * t27;
t16 = t30 * t227 - t230 * t29;
t255 = t15 * t264 + t16 * t263 + t5 * t307 + t6 * t308 + (t38 * t227 - t37 * t230) * t309;
t254 = -t123 + t260;
t141 = -rSges(5,3) * t212 + (rSges(5,1) * t228 - rSges(5,2) * t225) * t210;
t253 = -t141 + t260;
t252 = -t299 + t300;
t251 = rSges(4,1) * t212 - rSges(4,2) * t210;
t250 = -t173 * rSges(5,1) - t172 * rSges(5,2);
t249 = -t160 * rSges(6,1) - t159 * rSges(6,2);
t247 = Icges(3,1) * t229 - t294;
t246 = Icges(4,1) * t212 - t292;
t245 = -Icges(3,2) * t226 + t293;
t244 = -Icges(4,2) * t210 + t291;
t131 = -rSges(6,3) * t212 + (rSges(6,1) * t211 - rSges(6,2) * t209) * t210;
t237 = -t131 + t254;
t236 = rSges(4,1) * t286 - rSges(4,2) * t287 + t227 * rSges(4,3);
t234 = t57 / 0.2e1 + t49 / 0.2e1 + t41 / 0.2e1 + t63 / 0.2e1;
t233 = t62 / 0.2e1 + t56 / 0.2e1 + t48 / 0.2e1 + t40 / 0.2e1;
t232 = t227 * t104 + t230 * t105 + t256;
t231 = t254 + t316;
t192 = rSges(2,1) * t230 - t227 * rSges(2,2);
t191 = -t227 * rSges(2,1) - rSges(2,2) * t230;
t190 = rSges(3,1) * t226 + rSges(3,2) * t229;
t143 = t261 * t230;
t142 = t261 * t227;
t135 = t227 * pkin(7) + (pkin(1) - t299) * t230 + t269;
t134 = t298 + t217 + (-pkin(1) - t252) * t227;
t121 = t236 + t258;
t120 = (rSges(4,3) - t224) * t230 + (-t207 - t251) * t227;
t115 = t230 * (-t230 * t299 + t269) + (t252 * t227 - t298) * t227;
t112 = rSges(5,3) * t288 - t250;
t103 = t253 * t230;
t102 = t253 * t227;
t97 = rSges(6,3) * t288 - t249;
t90 = t104 * t287;
t80 = t88 * t287;
t78 = -t184 * t230 + (-t210 * t218 + t272 * t212) * t227 + t271;
t77 = t258 + t113 + t270;
t76 = -t285 + (-t305 - t207 + (-rSges(5,3) - pkin(8)) * t210) * t227 + t250;
t75 = -t212 * t113 - t141 * t287;
t74 = t112 * t212 + t141 * t288;
t73 = t237 * t230;
t72 = t237 * t227;
t70 = t235 + t258 + t98;
t69 = -t285 + (-rSges(6,3) * t210 - t206 * t212 - t207) * t227 + t249 + t271;
t68 = -t127 * t287 - t212 * t89;
t66 = t230 * t236 + (-t230 * rSges(4,3) + t251 * t227) * t227 + t275;
t65 = (t112 * t230 - t113 * t227) * t210;
t61 = -t218 * t287 + t258 + t317;
t60 = (t184 - t224) * t230 + (-t183 * t212 - t207 + (-rSges(7,3) + t218) * t210) * t227 + t248;
t58 = -t89 * t288 + t80;
t55 = t231 * t230;
t54 = t231 * t227;
t43 = t296 * t212 + (-t123 - t131) * t287;
t42 = t131 * t288 + t212 * t97 + t295;
t39 = t227 * t112 + t113 * t230 + t256;
t26 = t90 + (t296 * t227 + t230 * t97) * t210;
t25 = t227 * t97 + t230 * t98 + t232;
t24 = t265 * t212 + (-t123 + t316) * t287;
t23 = t114 * t288 + t212 * t78 + t295 + t67;
t22 = t47 * t227 - t230 * t46;
t21 = t45 * t227 - t230 * t44;
t20 = t80 + t90 + (t265 * t227 + t230 * t78) * t210;
t18 = t34 * t227 - t230 * t33;
t17 = t32 * t227 - t230 * t31;
t14 = t301 * t230 + (t78 + t88) * t227 + t232;
t1 = [t229 * (Icges(3,2) * t229 + t294) + t226 * (Icges(3,1) * t226 + t293) + Icges(2,3) + t116 + (Icges(4,2) * t212 - t124 + t292 + t327) * t212 + (Icges(4,1) * t210 - t290 + t291 + t326) * t210 + m(7) * (t60 ^ 2 + t61 ^ 2) + m(6) * (t69 ^ 2 + t70 ^ 2) + m(5) * (t76 ^ 2 + t77 ^ 2) + m(4) * (t120 ^ 2 + t121 ^ 2) + m(3) * (t134 ^ 2 + t135 ^ 2) + m(2) * (t191 ^ 2 + t192 ^ 2) + t315; (-(-Icges(3,6) * t230 + t245 * t227) * t229 / 0.2e1 - (-Icges(3,5) * t230 + t247 * t227) * t226 / 0.2e1 - t52 / 0.2e1 + (-Icges(4,6) * t230 + t244 * t227) * t309 - (-Icges(4,5) * t230 + t246 * t227) * t210 / 0.2e1 - t37 / 0.2e1 + t262 * t230 - t233) * t230 + (t53 / 0.2e1 + (Icges(4,6) * t227 + t244 * t230) * t324 + (Icges(4,5) * t227 + t246 * t230) * t325 + (Icges(3,6) * t227 + t245 * t230) * t322 + (Icges(3,5) * t227 + t247 * t230) * t323 + t38 / 0.2e1 + t262 * t227 + t234) * t227 + m(7) * (t54 * t61 + t55 * t60) + m(6) * (t69 * t73 + t70 * t72) + m(5) * (t102 * t77 + t103 * t76) + m(4) * (t120 * t143 + t121 * t142) + m(3) * (-t134 * t230 - t135 * t227) * t190; m(7) * (t14 ^ 2 + t54 ^ 2 + t55 ^ 2) + m(6) * (t25 ^ 2 + t72 ^ 2 + t73 ^ 2) + m(5) * (t102 ^ 2 + t103 ^ 2 + t39 ^ 2) + m(4) * (t142 ^ 2 + t143 ^ 2 + t66 ^ 2) + m(3) * (t267 * t190 ^ 2 + t115 ^ 2) + (t222 * t314 - t15 - t17 - t21) * t230 + (t16 + t18 + t22 + t313 * t221 + (t227 * t314 + t230 * t313) * t230) * t227; m(7) * (t227 * t60 - t230 * t61) + m(6) * (t227 * t69 - t230 * t70) + m(5) * (t227 * t76 - t230 * t77) + m(4) * (t227 * t120 - t121 * t230); m(7) * (t227 * t55 - t230 * t54) + m(6) * (t227 * t73 - t230 * t72) + m(5) * (-t102 * t230 + t227 * t103) + m(4) * (-t142 * t230 + t227 * t143); 0.2e1 * (m(4) / 0.2e1 + m(5) / 0.2e1 + t266) * t267; (-t59 - t302) * t212 + m(7) * (t23 * t60 + t24 * t61) + m(6) * (t42 * t69 + t43 * t70) + m(5) * (t74 * t76 + t75 * t77) + (t233 * t227 + t234 * t230) * t210 + t259; m(7) * (t14 * t20 + t23 * t55 + t24 * t54) + m(6) * (t26 * t25 + t42 * t73 + t43 * t72) + m(5) * (t102 * t75 + t103 * t74 + t65 * t39) + ((t22 / 0.2e1 + t18 / 0.2e1) * t230 + (t21 / 0.2e1 + t17 / 0.2e1) * t227) * t210 + t255 + (t318 * t227 + t319 * t230) * t309 + t320 * t308 + t321 * t307; m(5) * (t74 * t227 - t230 * t75) + m(6) * (t42 * t227 - t230 * t43) + m(7) * (t23 * t227 - t230 * t24); m(7) * (t20 ^ 2 + t23 ^ 2 + t24 ^ 2) + m(6) * (t26 ^ 2 + t42 ^ 2 + t43 ^ 2) + m(5) * (t65 ^ 2 + t74 ^ 2 + t75 ^ 2) + (t302 * t212 - t9) * t212 + (t320 * t230 + t321 * t227 + (t319 * t227 - t318 * t230) * t212) * t210 + t310; 0.2e1 * ((t227 * t61 + t230 * t60) * t311 + (t227 * t70 + t230 * t69) * t312) * t210; m(7) * (-t212 * t14 + (t227 * t54 + t230 * t55) * t210) + m(6) * (-t212 * t25 + (t227 * t72 + t230 * t73) * t210); 0; m(7) * (-t212 * t20 + (t227 * t24 + t23 * t230) * t210) + m(6) * (-t212 * t26 + (t227 * t43 + t230 * t42) * t210); 0.2e1 * t266 * (t267 * t210 ^ 2 + t212 ^ 2); m(7) * (t67 * t60 + t61 * t68) - t297 + t259; m(7) * (t58 * t14 + t54 * t68 + t67 * t55) + t255; m(7) * (t67 * t227 - t230 * t68); m(7) * (t58 * t20 + t67 * t23 + t24 * t68) + t257; m(7) * (-t58 * t212 + (t227 * t68 + t230 * t67) * t210); m(7) * (t58 ^ 2 + t67 ^ 2 + t68 ^ 2) + t257;];
%% Postprocessing: Reshape Output
% From vec2symmat_6_matlab.m
res = [t1(1) t1(2) t1(4) t1(7) t1(11) t1(16); t1(2) t1(3) t1(5) t1(8) t1(12) t1(17); t1(4) t1(5) t1(6) t1(9) t1(13) t1(18); t1(7) t1(8) t1(9) t1(10) t1(14) t1(19); t1(11) t1(12) t1(13) t1(14) t1(15) t1(20); t1(16) t1(17) t1(18) t1(19) t1(20) t1(21);];
Mq  = res;
