% Calculate joint inertia matrix for
% S6RPRRRP12
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% pkin [12x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,alpha2,alpha3,d1,d3,d4,d5,theta2]';
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
% Datum: 2019-03-09 06:53
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Mq = S6RPRRRP12_inertiaJ_slag_vp1(qJ, ...
  pkin, m, rSges, Icges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(12,1),zeros(7,1),zeros(7,3),zeros(7,6)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPRRRP12_inertiaJ_slag_vp1: qJ has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [12 1]), ...
  'S6RPRRRP12_inertiaJ_slag_vp1: pkin has to be [12x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6RPRRRP12_inertiaJ_slag_vp1: m has to be [7x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [7,3]), ...
  'S6RPRRRP12_inertiaJ_slag_vp1: rSges has to be [7x3] (double)');
assert(isreal(Icges) && all(size(Icges) == [7 6]), ...
  'S6RPRRRP12_inertiaJ_slag_vp1: Icges has to be [7x6] (double)'); 

%% Symbolic Calculation
% From inertia_joint_joint_floatb_twist_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 06:44:07
% EndTime: 2019-03-09 06:44:23
% DurationCPUTime: 6.92s
% Computational Cost: add. (42020->568), mult. (116549->788), div. (0->0), fcn. (154074->14), ass. (0->255)
t244 = cos(pkin(6));
t248 = sin(qJ(1));
t305 = cos(pkin(12));
t272 = t248 * t305;
t242 = sin(pkin(12));
t249 = cos(qJ(1));
t299 = t249 * t242;
t229 = t244 * t299 + t272;
t247 = sin(qJ(3));
t271 = t249 * t305;
t300 = t248 * t242;
t263 = -t244 * t271 + t300;
t306 = cos(pkin(7));
t258 = t263 * t306;
t243 = sin(pkin(6));
t304 = sin(pkin(7));
t273 = t243 * t304;
t310 = cos(qJ(3));
t211 = t229 * t310 + (-t249 * t273 - t258) * t247;
t274 = t243 * t306;
t220 = -t249 * t274 + t263 * t304;
t246 = sin(qJ(4));
t309 = cos(qJ(4));
t191 = t211 * t309 + t220 * t246;
t268 = t310 * t304;
t266 = t243 * t268;
t210 = t229 * t247 + t249 * t266 + t310 * t258;
t245 = sin(qJ(5));
t308 = cos(qJ(5));
t156 = t191 * t245 - t210 * t308;
t157 = t191 * t308 + t210 * t245;
t190 = t211 * t246 - t220 * t309;
t316 = rSges(7,3) + qJ(6);
t317 = rSges(7,1) + pkin(5);
t298 = rSges(7,2) * t190 + t316 * t156 + t317 * t157;
t262 = t244 * t272 + t299;
t221 = t248 * t274 + t262 * t304;
t315 = m(3) / 0.2e1;
t314 = m(4) / 0.2e1;
t313 = m(5) / 0.2e1;
t312 = m(6) / 0.2e1;
t311 = m(7) / 0.2e1;
t303 = t242 * t243;
t302 = t243 * t248;
t301 = t243 * t249;
t105 = rSges(6,1) * t157 - rSges(6,2) * t156 + rSges(6,3) * t190;
t147 = pkin(4) * t191 + t190 * pkin(11);
t297 = -t105 - t147;
t230 = -t244 * t300 + t271;
t256 = t262 * t306;
t213 = t230 * t310 + (t248 * t273 - t256) * t247;
t193 = t213 * t309 + t221 * t246;
t212 = t230 * t247 - t248 * t266 + t310 * t256;
t158 = t193 * t245 - t212 * t308;
t159 = t193 * t308 + t212 * t245;
t192 = t213 * t246 - t221 * t309;
t296 = t192 * rSges(7,2) + t316 * t158 + t317 * t159;
t107 = t159 * rSges(6,1) - t158 * rSges(6,2) + t192 * rSges(6,3);
t148 = t193 * pkin(4) + t192 * pkin(11);
t295 = -t107 - t148;
t267 = t306 * t305;
t219 = t244 * t304 * t247 + (t310 * t242 + t247 * t267) * t243;
t228 = t244 * t306 - t305 * t273;
t209 = t219 * t309 + t228 * t246;
t218 = -t243 * t310 * t267 - t244 * t268 + t247 * t303;
t182 = t209 * t245 - t218 * t308;
t183 = t209 * t308 + t218 * t245;
t208 = t219 * t246 - t228 * t309;
t294 = rSges(7,2) * t208 + t316 * t182 + t317 * t183;
t129 = rSges(6,1) * t183 - rSges(6,2) * t182 + rSges(6,3) * t208;
t177 = pkin(4) * t209 + pkin(11) * t208;
t293 = -t129 - t177;
t178 = pkin(3) * t211 + t210 * pkin(10);
t173 = t221 * t178;
t292 = t221 * t147 + t173;
t179 = t213 * pkin(3) + t212 * pkin(10);
t174 = t228 * t179;
t291 = t228 * t148 + t174;
t198 = pkin(3) * t219 + pkin(10) * t218;
t181 = t220 * t198;
t290 = t220 * t177 + t181;
t289 = t249 * pkin(1) + qJ(2) * t302;
t100 = Icges(7,1) * t157 + Icges(7,4) * t190 + Icges(7,5) * t156;
t92 = Icges(7,5) * t157 + Icges(7,6) * t190 + Icges(7,3) * t156;
t96 = Icges(7,4) * t157 + Icges(7,2) * t190 + Icges(7,6) * t156;
t32 = t100 * t157 + t156 * t92 + t190 * t96;
t101 = Icges(7,1) * t159 + Icges(7,4) * t192 + Icges(7,5) * t158;
t93 = Icges(7,5) * t159 + Icges(7,6) * t192 + Icges(7,3) * t158;
t97 = Icges(7,4) * t159 + Icges(7,2) * t192 + Icges(7,6) * t158;
t33 = t101 * t157 + t156 * t93 + t190 * t97;
t122 = Icges(7,5) * t183 + Icges(7,6) * t208 + Icges(7,3) * t182;
t124 = Icges(7,4) * t183 + Icges(7,2) * t208 + Icges(7,6) * t182;
t126 = Icges(7,1) * t183 + Icges(7,4) * t208 + Icges(7,5) * t182;
t50 = t122 * t156 + t124 * t190 + t126 * t157;
t1 = t190 * t32 + t192 * t33 + t208 * t50;
t102 = Icges(6,1) * t157 - Icges(6,4) * t156 + Icges(6,5) * t190;
t94 = Icges(6,5) * t157 - Icges(6,6) * t156 + Icges(6,3) * t190;
t98 = Icges(6,4) * t157 - Icges(6,2) * t156 + Icges(6,6) * t190;
t34 = t102 * t157 - t156 * t98 + t190 * t94;
t103 = Icges(6,1) * t159 - Icges(6,4) * t158 + Icges(6,5) * t192;
t95 = Icges(6,5) * t159 - Icges(6,6) * t158 + Icges(6,3) * t192;
t99 = Icges(6,4) * t159 - Icges(6,2) * t158 + Icges(6,6) * t192;
t35 = t103 * t157 - t156 * t99 + t190 * t95;
t123 = Icges(6,5) * t183 - Icges(6,6) * t182 + Icges(6,3) * t208;
t125 = Icges(6,4) * t183 - Icges(6,2) * t182 + Icges(6,6) * t208;
t127 = Icges(6,1) * t183 - Icges(6,4) * t182 + Icges(6,5) * t208;
t51 = t123 * t190 - t125 * t156 + t127 * t157;
t2 = t190 * t34 + t192 * t35 + t208 * t51;
t288 = t1 / 0.2e1 + t2 / 0.2e1;
t36 = t100 * t159 + t158 * t92 + t192 * t96;
t37 = t101 * t159 + t158 * t93 + t192 * t97;
t52 = t122 * t158 + t124 * t192 + t126 * t159;
t3 = t190 * t36 + t192 * t37 + t208 * t52;
t38 = t102 * t159 - t158 * t98 + t192 * t94;
t39 = t103 * t159 - t158 * t99 + t192 * t95;
t53 = t123 * t192 - t125 * t158 + t127 * t159;
t4 = t190 * t38 + t192 * t39 + t208 * t53;
t287 = t3 / 0.2e1 + t4 / 0.2e1;
t5 = t210 * t32 + t212 * t33 + t218 * t50;
t6 = t210 * t34 + t212 * t35 + t218 * t51;
t286 = t6 / 0.2e1 + t5 / 0.2e1;
t7 = t210 * t36 + t212 * t37 + t218 * t52;
t8 = t210 * t38 + t212 * t39 + t218 * t53;
t285 = t7 / 0.2e1 + t8 / 0.2e1;
t10 = t220 * t34 + t221 * t35 + t228 * t51;
t9 = t220 * t32 + t221 * t33 + t228 * t50;
t284 = t10 / 0.2e1 + t9 / 0.2e1;
t11 = t220 * t36 + t221 * t37 + t228 * t52;
t12 = t220 * t38 + t221 * t39 + t228 * t53;
t283 = t12 / 0.2e1 + t11 / 0.2e1;
t41 = t100 * t183 + t182 * t92 + t208 * t96;
t42 = t101 * t183 + t182 * t93 + t208 * t97;
t64 = t182 * t122 + t208 * t124 + t183 * t126;
t58 = t64 * t208;
t13 = t41 * t190 + t42 * t192 + t58;
t43 = t102 * t183 - t182 * t98 + t208 * t94;
t44 = t103 * t183 - t182 * t99 + t208 * t95;
t65 = t208 * t123 - t182 * t125 + t183 * t127;
t59 = t65 * t208;
t14 = t43 * t190 + t44 * t192 + t59;
t282 = t13 / 0.2e1 + t14 / 0.2e1;
t60 = t64 * t218;
t15 = t41 * t210 + t42 * t212 + t60;
t61 = t65 * t218;
t16 = t43 * t210 + t44 * t212 + t61;
t281 = t15 / 0.2e1 + t16 / 0.2e1;
t62 = t64 * t228;
t17 = t41 * t220 + t42 * t221 + t62;
t63 = t65 * t228;
t18 = t43 * t220 + t44 * t221 + t63;
t280 = t17 / 0.2e1 + t18 / 0.2e1;
t279 = -t147 - t298;
t278 = -t148 - t296;
t277 = -t177 - t294;
t160 = Icges(5,5) * t209 - Icges(5,6) * t208 + Icges(5,3) * t218;
t161 = Icges(5,4) * t209 - Icges(5,2) * t208 + Icges(5,6) * t218;
t162 = Icges(5,1) * t209 - Icges(5,4) * t208 + Icges(5,5) * t218;
t85 = t218 * t160 - t208 * t161 + t209 * t162;
t194 = Icges(4,5) * t219 - Icges(4,6) * t218 + Icges(4,3) * t228;
t195 = Icges(4,4) * t219 - Icges(4,2) * t218 + Icges(4,6) * t228;
t196 = Icges(4,1) * t219 - Icges(4,4) * t218 + Icges(4,5) * t228;
t276 = t228 * t194 - t218 * t195 + t219 * t196;
t137 = t193 * rSges(5,1) - t192 * rSges(5,2) + t212 * rSges(5,3);
t172 = t213 * rSges(4,1) - t212 * rSges(4,2) + t221 * rSges(4,3);
t275 = -t248 * pkin(1) + qJ(2) * t301;
t265 = t42 / 0.2e1 + t44 / 0.2e1 + t53 / 0.2e1 + t52 / 0.2e1;
t264 = t43 / 0.2e1 + t50 / 0.2e1 + t51 / 0.2e1 + t41 / 0.2e1;
t171 = rSges(4,1) * t211 - rSges(4,2) * t210 + rSges(4,3) * t220;
t136 = rSges(5,1) * t191 - rSges(5,2) * t190 + rSges(5,3) * t210;
t261 = -pkin(2) * t229 - t220 * pkin(9) + t275;
t131 = Icges(5,5) * t193 - Icges(5,6) * t192 + Icges(5,3) * t212;
t133 = Icges(5,4) * t193 - Icges(5,2) * t192 + Icges(5,6) * t212;
t135 = Icges(5,1) * t193 - Icges(5,4) * t192 + Icges(5,5) * t212;
t72 = t131 * t218 - t133 * t208 + t135 * t209;
t79 = t160 * t212 - t161 * t192 + t162 * t193;
t260 = t79 / 0.2e1 + t72 / 0.2e1 + t265;
t130 = Icges(5,5) * t191 - Icges(5,6) * t190 + Icges(5,3) * t210;
t132 = Icges(5,4) * t191 - Icges(5,2) * t190 + Icges(5,6) * t210;
t134 = Icges(5,1) * t191 - Icges(5,4) * t190 + Icges(5,5) * t210;
t71 = t130 * t218 - t132 * t208 + t134 * t209;
t78 = t160 * t210 - t161 * t190 + t162 * t191;
t259 = t78 / 0.2e1 + t71 / 0.2e1 + t264;
t254 = -t178 + t261;
t253 = -t147 + t254;
t252 = t230 * pkin(2) + t221 * pkin(9) + t289;
t251 = t179 + t252;
t250 = t148 + t251;
t236 = rSges(2,1) * t249 - t248 * rSges(2,2);
t235 = -t248 * rSges(2,1) - rSges(2,2) * t249;
t207 = t230 * rSges(3,1) - t262 * rSges(3,2) + rSges(3,3) * t302 + t289;
t206 = -t229 * rSges(3,1) + t263 * rSges(3,2) + rSges(3,3) * t301 + t275;
t197 = rSges(4,1) * t219 - rSges(4,2) * t218 + rSges(4,3) * t228;
t170 = Icges(4,1) * t213 - Icges(4,4) * t212 + Icges(4,5) * t221;
t169 = Icges(4,1) * t211 - Icges(4,4) * t210 + Icges(4,5) * t220;
t168 = Icges(4,4) * t213 - Icges(4,2) * t212 + Icges(4,6) * t221;
t167 = Icges(4,4) * t211 - Icges(4,2) * t210 + Icges(4,6) * t220;
t166 = Icges(4,5) * t213 - Icges(4,6) * t212 + Icges(4,3) * t221;
t165 = Icges(4,5) * t211 - Icges(4,6) * t210 + Icges(4,3) * t220;
t163 = rSges(5,1) * t209 - rSges(5,2) * t208 + rSges(5,3) * t218;
t150 = t210 * t177;
t144 = t252 + t172;
t143 = -t171 + t261;
t139 = t218 * t148;
t138 = t212 * t147;
t121 = t172 * t228 - t197 * t221;
t120 = -t171 * t228 + t197 * t220;
t111 = t171 * t221 - t172 * t220;
t110 = t276 * t228;
t109 = t194 * t221 - t195 * t212 + t196 * t213;
t108 = t194 * t220 - t195 * t210 + t196 * t211;
t91 = t251 + t137;
t90 = -t136 + t254;
t89 = t137 * t218 - t163 * t212;
t88 = -t136 * t218 + t163 * t210;
t87 = t166 * t228 - t168 * t218 + t170 * t219;
t86 = t165 * t228 - t167 * t218 + t169 * t219;
t84 = t136 * t212 - t137 * t210;
t83 = t85 * t228;
t82 = t85 * t218;
t81 = t137 * t228 + t174 + (-t163 - t198) * t221;
t80 = t163 * t220 + t181 + (-t136 - t178) * t228;
t77 = t250 + t107;
t76 = -t105 + t253;
t75 = t107 * t208 - t129 * t192;
t74 = -t105 * t208 + t129 * t190;
t73 = t136 * t221 + t173 + (-t137 - t179) * t220;
t70 = t131 * t212 - t133 * t192 + t135 * t193;
t69 = t130 * t212 - t132 * t192 + t134 * t193;
t68 = t131 * t210 - t133 * t190 + t135 * t191;
t67 = t130 * t210 - t132 * t190 + t134 * t191;
t66 = t105 * t192 - t107 * t190;
t57 = t250 + t296;
t56 = t253 - t298;
t55 = t107 * t218 + t293 * t212 + t139;
t54 = t129 * t210 + t297 * t218 + t150;
t49 = t107 * t228 + (-t198 + t293) * t221 + t291;
t48 = t129 * t220 + (-t178 + t297) * t228 + t290;
t47 = t105 * t212 + t295 * t210 + t138;
t46 = -t294 * t192 + t296 * t208;
t45 = t294 * t190 - t298 * t208;
t40 = t105 * t221 + (-t179 + t295) * t220 + t292;
t31 = t277 * t212 + t296 * t218 + t139;
t30 = t294 * t210 + t279 * t218 + t150;
t29 = -t296 * t190 + t298 * t192;
t28 = t296 * t228 + (-t198 + t277) * t221 + t291;
t27 = t294 * t220 + (-t178 + t279) * t228 + t290;
t26 = t278 * t210 + t298 * t212 + t138;
t25 = t298 * t221 + (-t179 + t278) * t220 + t292;
t24 = t71 * t220 + t72 * t221 + t83;
t23 = t71 * t210 + t72 * t212 + t82;
t22 = t220 * t69 + t221 * t70 + t228 * t79;
t21 = t220 * t67 + t221 * t68 + t228 * t78;
t20 = t210 * t69 + t212 * t70 + t218 * t79;
t19 = t210 * t67 + t212 * t68 + t218 * t78;
t104 = [m(7) * (t56 ^ 2 + t57 ^ 2) + m(6) * (t76 ^ 2 + t77 ^ 2) + m(5) * (t90 ^ 2 + t91 ^ 2) + m(4) * (t143 ^ 2 + t144 ^ 2) + m(3) * (t206 ^ 2 + t207 ^ 2) + m(2) * (t235 ^ 2 + t236 ^ 2) + t244 * (Icges(3,3) * t244 + (Icges(3,5) * t242 + t305 * Icges(3,6)) * t243) + (Icges(3,5) * t244 + (Icges(3,1) * t242 + t305 * Icges(3,4)) * t243) * t303 + t276 + Icges(2,3) + t243 * t305 * (Icges(3,6) * t244 + (Icges(3,4) * t242 + t305 * Icges(3,2)) * t243) + t85 + t65 + t64; 0.2e1 * ((t248 * t56 - t249 * t57) * t311 + (t248 * t76 - t249 * t77) * t312 + (t248 * t90 - t249 * t91) * t313 + (t143 * t248 - t144 * t249) * t314 + (t206 * t248 - t207 * t249) * t315) * t243; 0.2e1 * (t315 + t314 + t313 + t312 + t311) * (t244 ^ 2 + (t248 ^ 2 + t249 ^ 2) * t243 ^ 2); t62 + t63 + t83 + t110 + m(7) * (t27 * t56 + t28 * t57) + m(6) * (t48 * t76 + t49 * t77) + m(5) * (t80 * t90 + t81 * t91) + m(4) * (t120 * t143 + t121 * t144) + (t109 / 0.2e1 + t87 / 0.2e1 + t260) * t221 + (t108 / 0.2e1 + t86 / 0.2e1 + t259) * t220; m(4) * (t111 * t244 + (t120 * t248 - t121 * t249) * t243) + m(5) * (t73 * t244 + (t248 * t80 - t249 * t81) * t243) + m(6) * (t40 * t244 + (t248 * t48 - t249 * t49) * t243) + m(7) * (t25 * t244 + (t248 * t27 - t249 * t28) * t243); (t110 + t17 + t18 + t24) * t228 + (t11 + t12 + t22 + (t166 * t221 - t168 * t212 + t170 * t213) * t221 + (t109 + t87) * t228) * t221 + (t10 + t9 + t21 + (t165 * t220 - t167 * t210 + t169 * t211) * t220 + (t86 + t108) * t228 + (t165 * t221 + t166 * t220 - t167 * t212 - t168 * t210 + t169 * t213 + t170 * t211) * t221) * t220 + m(7) * (t25 ^ 2 + t27 ^ 2 + t28 ^ 2) + m(6) * (t40 ^ 2 + t48 ^ 2 + t49 ^ 2) + m(5) * (t73 ^ 2 + t80 ^ 2 + t81 ^ 2) + m(4) * (t111 ^ 2 + t120 ^ 2 + t121 ^ 2); t60 + t61 + t82 + m(7) * (t30 * t56 + t31 * t57) + m(6) * (t54 * t76 + t55 * t77) + m(5) * (t88 * t90 + t89 * t91) + t260 * t212 + t259 * t210; m(5) * (t84 * t244 + (t248 * t88 - t249 * t89) * t243) + m(6) * (t47 * t244 + (t248 * t54 - t249 * t55) * t243) + m(7) * (t26 * t244 + (t248 * t30 - t249 * t31) * t243); (t23 / 0.2e1 + t281) * t228 + (t20 / 0.2e1 + t285) * t221 + (t19 / 0.2e1 + t286) * t220 + (t24 / 0.2e1 + t280) * t218 + (t22 / 0.2e1 + t283) * t212 + (t21 / 0.2e1 + t284) * t210 + m(7) * (t25 * t26 + t27 * t30 + t28 * t31) + m(6) * (t40 * t47 + t48 * t54 + t49 * t55) + m(5) * (t73 * t84 + t80 * t88 + t81 * t89); (t15 + t16 + t23) * t218 + (t7 + t8 + t20) * t212 + (t6 + t5 + t19) * t210 + m(7) * (t26 ^ 2 + t30 ^ 2 + t31 ^ 2) + m(6) * (t47 ^ 2 + t54 ^ 2 + t55 ^ 2) + m(5) * (t84 ^ 2 + t88 ^ 2 + t89 ^ 2); t58 + t59 + m(7) * (t45 * t56 + t46 * t57) + m(6) * (t74 * t76 + t75 * t77) + t265 * t192 + t264 * t190; m(6) * (t66 * t244 + (t248 * t74 - t249 * t75) * t243) + m(7) * (t29 * t244 + (t248 * t45 - t249 * t46) * t243); t282 * t228 + t287 * t221 + t288 * t220 + t280 * t208 + t283 * t192 + t284 * t190 + m(7) * (t25 * t29 + t27 * t45 + t28 * t46) + m(6) * (t66 * t40 + t48 * t74 + t49 * t75); t282 * t218 + t287 * t212 + t288 * t210 + t281 * t208 + t285 * t192 + t286 * t190 + m(7) * (t26 * t29 + t30 * t45 + t31 * t46) + m(6) * (t66 * t47 + t54 * t74 + t55 * t75); (t13 + t14) * t208 + (t3 + t4) * t192 + (t1 + t2) * t190 + m(7) * (t29 ^ 2 + t45 ^ 2 + t46 ^ 2) + m(6) * (t66 ^ 2 + t74 ^ 2 + t75 ^ 2); m(7) * (t156 * t57 + t158 * t56); m(7) * (t182 * t244 + (-t156 * t249 + t158 * t248) * t243); m(7) * (t156 * t28 + t158 * t27 + t182 * t25); m(7) * (t156 * t31 + t158 * t30 + t182 * t26); m(7) * (t156 * t46 + t158 * t45 + t182 * t29); m(7) * (t156 ^ 2 + t158 ^ 2 + t182 ^ 2);];
%% Postprocessing: Reshape Output
% From vec2symmat_6_matlab.m
res = [t104(1) t104(2) t104(4) t104(7) t104(11) t104(16); t104(2) t104(3) t104(5) t104(8) t104(12) t104(17); t104(4) t104(5) t104(6) t104(9) t104(13) t104(18); t104(7) t104(8) t104(9) t104(10) t104(14) t104(19); t104(11) t104(12) t104(13) t104(14) t104(15) t104(20); t104(16) t104(17) t104(18) t104(19) t104(20) t104(21);];
Mq  = res;
