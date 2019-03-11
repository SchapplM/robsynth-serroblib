% Calculate joint inertia matrix for
% S6RRRPRP9
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% pkin [9x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d2,d3,d5]';
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
% Datum: 2019-03-09 17:28
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Mq = S6RRRPRP9_inertiaJ_slag_vp1(qJ, ...
  pkin, m, rSges, Icges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(9,1),zeros(7,1),zeros(7,3),zeros(7,6)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRRPRP9_inertiaJ_slag_vp1: qJ has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S6RRRPRP9_inertiaJ_slag_vp1: pkin has to be [9x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6RRRPRP9_inertiaJ_slag_vp1: m has to be [7x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [7,3]), ...
  'S6RRRPRP9_inertiaJ_slag_vp1: rSges has to be [7x3] (double)');
assert(isreal(Icges) && all(size(Icges) == [7 6]), ...
  'S6RRRPRP9_inertiaJ_slag_vp1: Icges has to be [7x6] (double)'); 

%% Symbolic Calculation
% From inertia_joint_joint_floatb_twist_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 17:21:54
% EndTime: 2019-03-09 17:22:05
% DurationCPUTime: 5.24s
% Computational Cost: add. (9108->543), mult. (22817->768), div. (0->0), fcn. (27538->8), ass. (0->255)
t333 = rSges(7,1) + pkin(5);
t330 = rSges(7,3) + qJ(6);
t231 = sin(qJ(5));
t233 = sin(qJ(2));
t235 = cos(qJ(3));
t299 = t233 * t235;
t232 = sin(qJ(3));
t301 = t232 * t233;
t309 = cos(qJ(5));
t193 = t231 * t299 - t301 * t309;
t194 = (t231 * t232 + t235 * t309) * t233;
t236 = cos(qJ(2));
t122 = Icges(7,5) * t194 + Icges(7,6) * t236 + Icges(7,3) * t193;
t123 = Icges(6,5) * t194 - Icges(6,6) * t193 + Icges(6,3) * t236;
t124 = Icges(7,4) * t194 + Icges(7,2) * t236 + Icges(7,6) * t193;
t125 = Icges(6,4) * t194 - Icges(6,2) * t193 + Icges(6,6) * t236;
t126 = Icges(7,1) * t194 + Icges(7,4) * t236 + Icges(7,5) * t193;
t127 = Icges(6,1) * t194 - Icges(6,4) * t193 + Icges(6,5) * t236;
t335 = (t123 + t124) * t236 + (t126 + t127) * t194 + (t122 - t125) * t193;
t237 = cos(qJ(1));
t294 = t237 * t235;
t234 = sin(qJ(1));
t297 = t234 * t236;
t205 = t232 * t297 + t294;
t295 = t237 * t232;
t206 = t235 * t297 - t295;
t157 = -t205 * t309 + t206 * t231;
t158 = t205 * t231 + t206 * t309;
t334 = -t330 * t157 - t333 * t158;
t300 = t233 * t234;
t89 = Icges(7,5) * t158 - Icges(7,6) * t300 + Icges(7,3) * t157;
t93 = Icges(7,4) * t158 - Icges(7,2) * t300 + Icges(7,6) * t157;
t97 = Icges(7,1) * t158 - Icges(7,4) * t300 + Icges(7,5) * t157;
t36 = t193 * t89 + t194 * t97 + t236 * t93;
t91 = Icges(6,5) * t158 - Icges(6,6) * t157 - Icges(6,3) * t300;
t95 = Icges(6,4) * t158 - Icges(6,2) * t157 - Icges(6,6) * t300;
t99 = Icges(6,1) * t158 - Icges(6,4) * t157 - Icges(6,5) * t300;
t38 = -t193 * t95 + t194 * t99 + t236 * t91;
t332 = t36 + t38;
t207 = -t234 * t235 + t236 * t295;
t208 = t234 * t232 + t236 * t294;
t159 = -t207 * t309 + t208 * t231;
t160 = t207 * t231 + t208 * t309;
t298 = t233 * t237;
t90 = Icges(7,5) * t160 - Icges(7,6) * t298 + Icges(7,3) * t159;
t94 = Icges(7,4) * t160 - Icges(7,2) * t298 + Icges(7,6) * t159;
t98 = Icges(7,1) * t160 - Icges(7,4) * t298 + Icges(7,5) * t159;
t37 = t193 * t90 + t194 * t98 + t236 * t94;
t100 = Icges(6,1) * t160 - Icges(6,4) * t159 - Icges(6,5) * t298;
t92 = Icges(6,5) * t160 - Icges(6,6) * t159 - Icges(6,3) * t298;
t96 = Icges(6,4) * t160 - Icges(6,2) * t159 - Icges(6,6) * t298;
t39 = t194 * t100 - t193 * t96 + t236 * t92;
t331 = t37 + t39;
t178 = -Icges(4,3) * t236 + (Icges(4,5) * t235 - Icges(4,6) * t232) * t233;
t181 = -Icges(5,2) * t236 + (Icges(5,4) * t235 + Icges(5,6) * t232) * t233;
t329 = -t178 - t181;
t307 = t335 * t236;
t311 = -t307 + (t332 * t234 + t331 * t237) * t233;
t27 = t159 * t89 + t160 * t97 - t298 * t93;
t28 = t159 * t90 + t160 * t98 - t298 * t94;
t52 = t159 * t122 - t124 * t298 + t160 * t126;
t3 = -t52 * t236 + (t234 * t27 + t237 * t28) * t233;
t29 = -t159 * t95 + t160 * t99 - t298 * t91;
t30 = t160 * t100 - t159 * t96 - t298 * t92;
t53 = -t123 * t298 - t159 * t125 + t160 * t127;
t4 = -t53 * t236 + (t234 * t29 + t237 * t30) * t233;
t312 = t3 + t4;
t23 = t157 * t89 + t158 * t97 - t300 * t93;
t24 = t157 * t90 + t158 * t98 - t300 * t94;
t50 = t157 * t122 - t124 * t300 + t158 * t126;
t1 = -t50 * t236 + (t23 * t234 + t237 * t24) * t233;
t25 = -t157 * t95 + t158 * t99 - t300 * t91;
t26 = t158 * t100 - t157 * t96 - t300 * t92;
t51 = -t123 * t300 - t157 * t125 + t158 * t127;
t2 = -t51 * t236 + (t234 * t25 + t237 * t26) * t233;
t313 = t1 + t2;
t328 = (t234 * t313 + t312 * t237) * t233 - t311 * t236;
t327 = t233 / 0.2e1;
t325 = t236 / 0.2e1;
t293 = -rSges(7,2) * t300 - t334;
t323 = t236 * t293;
t322 = t237 * t293;
t321 = (t331 * t234 / 0.2e1 - t332 * t237 / 0.2e1) * t236;
t182 = -Icges(4,6) * t236 + (Icges(4,4) * t235 - Icges(4,2) * t232) * t233;
t177 = -Icges(5,6) * t236 + (Icges(5,5) * t235 + Icges(5,3) * t232) * t233;
t185 = -Icges(5,4) * t236 + (Icges(5,1) * t235 + Icges(5,5) * t232) * t233;
t186 = -Icges(4,5) * t236 + (Icges(4,1) * t235 - Icges(4,4) * t232) * t233;
t317 = t177 * t301 + (t185 + t186) * t299;
t320 = (-t182 * t301 + t329 * t236 + t317) * t236;
t318 = t330 * t159 + t160 * t333;
t315 = t234 ^ 2;
t314 = t237 ^ 2;
t310 = Icges(3,5) * t327 + Icges(3,6) * t325;
t308 = pkin(2) * t236;
t306 = t205 * rSges(5,3);
t305 = t237 * rSges(3,3);
t304 = Icges(3,4) * t233;
t303 = Icges(3,4) * t236;
t253 = -t158 * rSges(6,1) + t157 * rSges(6,2);
t102 = -rSges(6,3) * t300 - t253;
t302 = t102 * t237;
t296 = t236 * t237;
t292 = -rSges(7,2) * t298 + t318;
t290 = t236 * rSges(7,2) + t330 * t193 + t333 * t194;
t144 = t208 * rSges(5,1) + rSges(5,2) * t298 + t207 * rSges(5,3);
t164 = t208 * pkin(3) + t207 * qJ(4);
t289 = -t144 - t164;
t198 = t205 * qJ(4);
t163 = t206 * pkin(3) + t198;
t146 = t163 * t298;
t223 = pkin(9) * t300;
t174 = t206 * pkin(4) - t223;
t288 = t174 * t298 + t146;
t286 = t160 * rSges(6,1) - t159 * rSges(6,2);
t209 = (pkin(3) * t235 + qJ(4) * t232) * t233;
t285 = t236 * t163 + t209 * t300;
t203 = t208 * pkin(4);
t175 = -pkin(9) * t298 + t203;
t284 = -t164 - t175;
t189 = -t236 * rSges(5,2) + (rSges(5,1) * t235 + rSges(5,3) * t232) * t233;
t282 = -t189 - t209;
t190 = -t236 * rSges(4,3) + (rSges(4,1) * t235 - rSges(4,2) * t232) * t233;
t218 = t233 * pkin(2) - t236 * pkin(8);
t281 = -t190 - t218;
t278 = pkin(2) * t296 + pkin(8) * t298;
t280 = t315 * (pkin(8) * t233 + t308) + t237 * t278;
t210 = pkin(4) * t299 + t236 * pkin(9);
t279 = -t209 - t210;
t277 = t237 * pkin(1) + t234 * pkin(7);
t229 = t237 * pkin(7);
t276 = t229 - t198;
t275 = t1 / 0.2e1 + t2 / 0.2e1;
t274 = -t3 / 0.2e1 - t4 / 0.2e1;
t7 = -t23 * t237 + t24 * t234;
t8 = t26 * t234 - t25 * t237;
t273 = -t8 / 0.2e1 - t7 / 0.2e1;
t10 = t30 * t234 - t29 * t237;
t9 = t28 * t234 - t27 * t237;
t272 = -t9 / 0.2e1 - t10 / 0.2e1;
t130 = Icges(5,5) * t206 + Icges(5,6) * t300 + Icges(5,3) * t205;
t134 = Icges(5,4) * t206 + Icges(5,2) * t300 + Icges(5,6) * t205;
t138 = Icges(5,1) * t206 + Icges(5,4) * t300 + Icges(5,5) * t205;
t69 = -t236 * t134 + (t130 * t232 + t138 * t235) * t233;
t132 = Icges(4,5) * t206 - Icges(4,6) * t205 + Icges(4,3) * t300;
t136 = Icges(4,4) * t206 - Icges(4,2) * t205 + Icges(4,6) * t300;
t140 = Icges(4,1) * t206 - Icges(4,4) * t205 + Icges(4,5) * t300;
t71 = -t236 * t132 + (-t136 * t232 + t140 * t235) * t233;
t270 = t71 / 0.2e1 + t69 / 0.2e1;
t131 = Icges(5,5) * t208 + Icges(5,6) * t298 + Icges(5,3) * t207;
t135 = Icges(5,4) * t208 + Icges(5,2) * t298 + Icges(5,6) * t207;
t139 = Icges(5,1) * t208 + Icges(5,4) * t298 + Icges(5,5) * t207;
t70 = -t236 * t135 + (t131 * t232 + t139 * t235) * t233;
t133 = Icges(4,5) * t208 - Icges(4,6) * t207 + Icges(4,3) * t298;
t137 = Icges(4,4) * t208 - Icges(4,2) * t207 + Icges(4,6) * t298;
t141 = Icges(4,1) * t208 - Icges(4,4) * t207 + Icges(4,5) * t298;
t72 = -t236 * t133 + (-t137 * t232 + t141 * t235) * t233;
t269 = -t72 / 0.2e1 - t70 / 0.2e1;
t104 = -rSges(6,3) * t298 + t286;
t268 = -t104 + t284;
t129 = t194 * rSges(6,1) - t193 * rSges(6,2) + t236 * rSges(6,3);
t265 = -t129 + t279;
t264 = -t218 + t282;
t145 = t208 * rSges(4,1) - t207 * rSges(4,2) + rSges(4,3) * t298;
t263 = -pkin(1) - t308;
t262 = t233 * t290;
t261 = t284 - t292;
t260 = t279 - t290;
t259 = -t218 + t265;
t258 = t234 * t163 + t237 * t164 + t280;
t257 = t236 * t174 + t210 * t300 + t285;
t256 = t277 + t278;
t255 = rSges(3,1) * t236 - rSges(3,2) * t233;
t254 = -t206 * rSges(4,1) + t205 * rSges(4,2);
t252 = -t218 + t260;
t251 = Icges(3,1) * t236 - t304;
t250 = -Icges(3,2) * t233 + t303;
t249 = Icges(3,5) * t236 - Icges(3,6) * t233;
t246 = rSges(3,1) * t296 - rSges(3,2) * t298 + t234 * rSges(3,3);
t245 = -t36 / 0.2e1 - t51 / 0.2e1 - t50 / 0.2e1 - t38 / 0.2e1;
t244 = -t39 / 0.2e1 - t37 / 0.2e1 - t53 / 0.2e1 - t52 / 0.2e1;
t243 = t234 * t174 + t237 * t175 + t258;
t73 = -t236 * t102 - t129 * t300;
t242 = t223 + (-pkin(3) - pkin(4)) * t206 + t276;
t241 = t164 + t256;
t240 = t203 + t241;
t82 = t205 * t177 + t181 * t300 + t206 * t185;
t83 = t178 * t300 - t205 * t182 + t206 * t186;
t239 = t82 / 0.2e1 + t83 / 0.2e1 - t245 + t270;
t84 = t207 * t177 + t181 * t298 + t208 * t185;
t85 = t178 * t298 - t207 * t182 + t208 * t186;
t238 = t85 / 0.2e1 + t84 / 0.2e1 - t244 - t269;
t217 = t237 * rSges(2,1) - t234 * rSges(2,2);
t216 = -t234 * rSges(2,1) - t237 * rSges(2,2);
t215 = t233 * rSges(3,1) + t236 * rSges(3,2);
t180 = Icges(3,3) * t234 + t237 * t249;
t179 = -Icges(3,3) * t237 + t234 * t249;
t170 = t246 + t277;
t169 = t305 + t229 + (-pkin(1) - t255) * t234;
t162 = t281 * t237;
t161 = t281 * t234;
t143 = rSges(4,3) * t300 - t254;
t142 = t206 * rSges(5,1) + rSges(5,2) * t300 + t306;
t121 = t237 * t246 + (t234 * t255 - t305) * t234;
t118 = t264 * t237;
t117 = t264 * t234;
t112 = t256 + t145;
t111 = t229 + ((-rSges(4,3) - pkin(8)) * t233 + t263) * t234 + t254;
t110 = -t236 * t145 - t190 * t298;
t109 = t236 * t143 + t190 * t300;
t88 = t259 * t237;
t87 = t259 * t234;
t86 = (t143 * t237 - t145 * t234) * t233;
t81 = t241 + t144;
t80 = -t306 + (-rSges(5,1) - pkin(3)) * t206 + ((-rSges(5,2) - pkin(8)) * t233 + t263) * t234 + t276;
t79 = t234 * t143 + t237 * t145 + t280;
t78 = t236 * t289 + t282 * t298;
t77 = t236 * t142 + t189 * t300 + t285;
t76 = t252 * t237;
t75 = t252 * t234;
t74 = t236 * t104 + t129 * t298;
t68 = t133 * t298 - t207 * t137 + t208 * t141;
t67 = t132 * t298 - t207 * t136 + t208 * t140;
t66 = t207 * t131 + t135 * t298 + t208 * t139;
t65 = t207 * t130 + t134 * t298 + t208 * t138;
t64 = t133 * t300 - t205 * t137 + t206 * t141;
t63 = t132 * t300 - t205 * t136 + t206 * t140;
t62 = t205 * t131 + t135 * t300 + t206 * t139;
t61 = t205 * t130 + t134 * t300 + t206 * t138;
t60 = (-rSges(6,3) - pkin(9)) * t298 + t240 + t286;
t59 = ((rSges(6,3) - pkin(8)) * t233 + t263) * t234 + t242 + t253;
t58 = t146 + (t142 * t237 + t234 * t289) * t233;
t55 = (t104 * t234 - t302) * t233;
t54 = t234 * t142 + t237 * t144 + t258;
t45 = t236 * t268 + t265 * t298;
t44 = -t73 + t257;
t43 = (-rSges(7,2) - pkin(9)) * t298 + t240 + t318;
t42 = ((rSges(7,2) - pkin(8)) * t233 + t263) * t234 + t242 + t334;
t41 = t236 * t292 + t237 * t262;
t40 = -t290 * t300 - t323;
t35 = (t234 * t268 + t302) * t233 + t288;
t34 = t68 * t234 - t67 * t237;
t33 = t66 * t234 - t65 * t237;
t32 = t64 * t234 - t63 * t237;
t31 = t62 * t234 - t61 * t237;
t22 = t234 * t102 + t237 * t104 + t243;
t21 = (t234 * t292 - t322) * t233;
t20 = t236 * t261 + t260 * t298;
t19 = t234 * t262 + t257 + t323;
t18 = -t85 * t236 + (t234 * t67 + t237 * t68) * t233;
t17 = -t84 * t236 + (t234 * t65 + t237 * t66) * t233;
t16 = -t83 * t236 + (t234 * t63 + t237 * t64) * t233;
t15 = -t82 * t236 + (t234 * t61 + t237 * t62) * t233;
t14 = (t234 * t261 + t322) * t233 + t288;
t13 = t234 * t293 + t237 * t292 + t243;
t5 = [Icges(2,3) + (Icges(3,1) * t233 - t232 * t182 + t303) * t233 + (Icges(3,2) * t236 + t304 + t329) * t236 + m(7) * (t42 ^ 2 + t43 ^ 2) + m(6) * (t59 ^ 2 + t60 ^ 2) + m(4) * (t111 ^ 2 + t112 ^ 2) + m(5) * (t80 ^ 2 + t81 ^ 2) + m(3) * (t169 ^ 2 + t170 ^ 2) + m(2) * (t216 ^ 2 + t217 ^ 2) + t317 + t335; (-t236 * (-Icges(3,6) * t237 + t234 * t250) / 0.2e1 - t233 * (-Icges(3,5) * t237 + t234 * t251) / 0.2e1 + t237 * t310 - t239) * t237 + ((Icges(3,6) * t234 + t237 * t250) * t325 + (Icges(3,5) * t234 + t237 * t251) * t327 + t234 * t310 + t238) * t234 + m(7) * (t42 * t76 + t43 * t75) + m(6) * (t59 * t88 + t60 * t87) + m(5) * (t117 * t81 + t118 * t80) + m(4) * (t111 * t162 + t112 * t161) + m(3) * (-t169 * t237 - t170 * t234) * t215; m(4) * (t161 ^ 2 + t162 ^ 2 + t79 ^ 2) + m(3) * (t121 ^ 2 + (t314 + t315) * t215 ^ 2) + m(7) * (t13 ^ 2 + t75 ^ 2 + t76 ^ 2) + m(6) * (t22 ^ 2 + t87 ^ 2 + t88 ^ 2) + m(5) * (t117 ^ 2 + t118 ^ 2 + t54 ^ 2) + (-t314 * t179 - t31 - t32 - t7 - t8) * t237 + (t315 * t180 + t10 + t33 + t34 + t9 + (-t234 * t179 + t237 * t180) * t237) * t234; -t320 + m(7) * (t19 * t42 + t20 * t43) + m(6) * (t44 * t59 + t45 * t60) + m(4) * (t109 * t111 + t110 * t112) + m(5) * (t77 * t80 + t78 * t81) + (t234 * t239 + t237 * t238) * t233 - t307; -t321 + (-t15 / 0.2e1 - t16 / 0.2e1 + t270 * t236 - t275) * t237 + (t18 / 0.2e1 + t17 / 0.2e1 + t269 * t236 - t274) * t234 + m(7) * (t14 * t13 + t19 * t76 + t20 * t75) + m(6) * (t35 * t22 + t44 * t88 + t45 * t87) + m(5) * (t117 * t78 + t118 * t77 + t58 * t54) + m(4) * (t109 * t162 + t110 * t161 + t79 * t86) + ((t34 / 0.2e1 + t33 / 0.2e1 - t272) * t237 + (t32 / 0.2e1 + t31 / 0.2e1 - t273) * t234) * t233; (-t311 + t320) * t236 + m(6) * (t35 ^ 2 + t44 ^ 2 + t45 ^ 2) + m(7) * (t14 ^ 2 + t19 ^ 2 + t20 ^ 2) + m(5) * (t58 ^ 2 + t77 ^ 2 + t78 ^ 2) + m(4) * (t109 ^ 2 + t110 ^ 2 + t86 ^ 2) + ((t17 + t18 + (-t70 - t72) * t236 + t312) * t237 + (t15 + t16 + (-t69 - t71) * t236 + t313) * t234) * t233; m(7) * (t205 * t43 + t207 * t42) + m(6) * (t205 * t60 + t207 * t59) + m(5) * (t205 * t81 + t207 * t80); m(7) * (t13 * t301 + t205 * t75 + t207 * t76) + m(6) * (t205 * t87 + t207 * t88 + t22 * t301) + m(5) * (t205 * t117 + t207 * t118 + t301 * t54); m(6) * (t205 * t45 + t207 * t44 + t301 * t35) + m(7) * (t14 * t301 + t207 * t19 + t205 * t20) + m(5) * (t205 * t78 + t207 * t77 + t301 * t58); 0.2e1 * (m(5) / 0.2e1 + m(6) / 0.2e1 + m(7) / 0.2e1) * (t233 ^ 2 * t232 ^ 2 + t205 ^ 2 + t207 ^ 2); m(7) * (t40 * t42 + t41 * t43) + m(6) * (t59 * t73 + t60 * t74) + (t234 * t245 + t237 * t244) * t233 + t307; t275 * t237 + t321 + t274 * t234 + m(7) * (t21 * t13 + t40 * t76 + t41 * t75) + m(6) * (t55 * t22 + t73 * t88 + t74 * t87) + (t234 * t273 + t237 * t272) * t233; m(6) * (t55 * t35 + t44 * t73 + t45 * t74) + m(7) * (t14 * t21 + t19 * t40 + t20 * t41) - t328; m(6) * (t74 * t205 + t73 * t207 + t301 * t55) + m(7) * (t41 * t205 + t40 * t207 + t21 * t301); m(7) * (t21 ^ 2 + t40 ^ 2 + t41 ^ 2) + m(6) * (t55 ^ 2 + t73 ^ 2 + t74 ^ 2) + t328; m(7) * (t157 * t43 + t159 * t42); m(7) * (t13 * t193 + t157 * t75 + t159 * t76); m(7) * (t14 * t193 + t157 * t20 + t159 * t19); m(7) * (t157 * t205 + t159 * t207 + t193 * t301); m(7) * (t157 * t41 + t159 * t40 + t193 * t21); m(7) * (t157 ^ 2 + t159 ^ 2 + t193 ^ 2);];
%% Postprocessing: Reshape Output
% From vec2symmat_6_matlab.m
res = [t5(1) t5(2) t5(4) t5(7) t5(11) t5(16); t5(2) t5(3) t5(5) t5(8) t5(12) t5(17); t5(4) t5(5) t5(6) t5(9) t5(13) t5(18); t5(7) t5(8) t5(9) t5(10) t5(14) t5(19); t5(11) t5(12) t5(13) t5(14) t5(15) t5(20); t5(16) t5(17) t5(18) t5(19) t5(20) t5(21);];
Mq  = res;
