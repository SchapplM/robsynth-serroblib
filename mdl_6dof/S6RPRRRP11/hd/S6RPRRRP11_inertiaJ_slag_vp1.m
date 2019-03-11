% Calculate joint inertia matrix for
% S6RPRRRP11
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
% Datum: 2019-03-09 06:42
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Mq = S6RPRRRP11_inertiaJ_slag_vp1(qJ, ...
  pkin, m, rSges, Icges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(12,1),zeros(7,1),zeros(7,3),zeros(7,6)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPRRRP11_inertiaJ_slag_vp1: qJ has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [12 1]), ...
  'S6RPRRRP11_inertiaJ_slag_vp1: pkin has to be [12x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6RPRRRP11_inertiaJ_slag_vp1: m has to be [7x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [7,3]), ...
  'S6RPRRRP11_inertiaJ_slag_vp1: rSges has to be [7x3] (double)');
assert(isreal(Icges) && all(size(Icges) == [7 6]), ...
  'S6RPRRRP11_inertiaJ_slag_vp1: Icges has to be [7x6] (double)'); 

%% Symbolic Calculation
% From inertia_joint_joint_floatb_twist_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 06:34:33
% EndTime: 2019-03-09 06:34:48
% DurationCPUTime: 7.14s
% Computational Cost: add. (42899->577), mult. (118238->796), div. (0->0), fcn. (156120->14), ass. (0->260)
t323 = rSges(7,3) + qJ(6) + pkin(11);
t244 = cos(pkin(6));
t249 = sin(qJ(1));
t309 = cos(pkin(12));
t273 = t249 * t309;
t242 = sin(pkin(12));
t251 = cos(qJ(1));
t300 = t251 * t242;
t228 = t244 * t300 + t273;
t248 = sin(qJ(3));
t272 = t251 * t309;
t301 = t249 * t242;
t263 = -t244 * t272 + t301;
t310 = cos(pkin(7));
t258 = t263 * t310;
t243 = sin(pkin(6));
t308 = sin(pkin(7));
t274 = t243 * t308;
t315 = cos(qJ(3));
t211 = t228 * t315 + (-t251 * t274 - t258) * t248;
t275 = t243 * t310;
t219 = -t251 * t275 + t263 * t308;
t247 = sin(qJ(4));
t314 = cos(qJ(4));
t192 = t211 * t314 + t219 * t247;
t269 = t315 * t308;
t266 = t243 * t269;
t210 = t228 * t248 + t251 * t266 + t258 * t315;
t246 = sin(qJ(5));
t250 = cos(qJ(5));
t155 = -t192 * t246 + t210 * t250;
t307 = t210 * t246;
t156 = t192 * t250 + t307;
t191 = t211 * t247 - t219 * t314;
t324 = t156 * rSges(7,1) + t155 * rSges(7,2) + pkin(5) * t307 + t323 * t191;
t262 = t244 * t273 + t300;
t220 = t249 * t275 + t262 * t308;
t229 = -t244 * t301 + t272;
t256 = t262 * t310;
t213 = t229 * t315 + (t249 * t274 - t256) * t248;
t194 = t213 * t314 + t220 * t247;
t212 = t229 * t248 - t249 * t266 + t256 * t315;
t157 = -t194 * t246 + t212 * t250;
t306 = t212 * t246;
t158 = t194 * t250 + t306;
t193 = t213 * t247 - t220 * t314;
t239 = pkin(5) * t250 + pkin(4);
t322 = t158 * rSges(7,1) + t157 * rSges(7,2) + pkin(5) * t306 + t323 * t193 + t194 * t239;
t320 = m(3) / 0.2e1;
t319 = m(4) / 0.2e1;
t318 = m(5) / 0.2e1;
t317 = m(6) / 0.2e1;
t316 = m(7) / 0.2e1;
t313 = -pkin(4) + t239;
t189 = t191 * pkin(11);
t312 = t192 * t313 - t189 + t324;
t148 = t194 * pkin(4) + t193 * pkin(11);
t311 = -t148 + t322;
t268 = t310 * t309;
t304 = t242 * t243;
t217 = -t243 * t268 * t315 - t244 * t269 + t248 * t304;
t305 = t217 * t246;
t303 = t243 * t249;
t302 = t243 * t251;
t107 = rSges(6,1) * t156 + rSges(6,2) * t155 + rSges(6,3) * t191;
t147 = pkin(4) * t192 + t189;
t299 = -t107 - t147;
t109 = t158 * rSges(6,1) + t157 * rSges(6,2) + t193 * rSges(6,3);
t298 = -t109 - t148;
t218 = t244 * t308 * t248 + (t242 * t315 + t248 * t268) * t243;
t227 = t244 * t310 - t274 * t309;
t209 = t218 * t314 + t227 * t247;
t183 = -t209 * t246 + t217 * t250;
t184 = t209 * t250 + t305;
t208 = t218 * t247 - t227 * t314;
t297 = rSges(7,1) * t184 + rSges(7,2) * t183 + pkin(5) * t305 + t313 * t209 + (-pkin(11) + t323) * t208;
t130 = rSges(6,1) * t184 + rSges(6,2) * t183 + rSges(6,3) * t208;
t176 = t209 * pkin(4) + t208 * pkin(11);
t296 = -t130 - t176;
t177 = t211 * pkin(3) + t210 * pkin(10);
t172 = t220 * t177;
t295 = t220 * t147 + t172;
t178 = t213 * pkin(3) + t212 * pkin(10);
t173 = t227 * t178;
t294 = t227 * t148 + t173;
t199 = pkin(3) * t218 + pkin(10) * t217;
t181 = t219 * t199;
t293 = t219 * t176 + t181;
t292 = t251 * pkin(1) + qJ(2) * t303;
t102 = Icges(7,1) * t156 + Icges(7,4) * t155 + Icges(7,5) * t191;
t94 = Icges(7,5) * t156 + Icges(7,6) * t155 + Icges(7,3) * t191;
t98 = Icges(7,4) * t156 + Icges(7,2) * t155 + Icges(7,6) * t191;
t34 = t102 * t156 + t155 * t98 + t191 * t94;
t103 = Icges(7,1) * t158 + Icges(7,4) * t157 + Icges(7,5) * t193;
t95 = Icges(7,5) * t158 + Icges(7,6) * t157 + Icges(7,3) * t193;
t99 = Icges(7,4) * t158 + Icges(7,2) * t157 + Icges(7,6) * t193;
t35 = t103 * t156 + t155 * t99 + t191 * t95;
t123 = Icges(7,5) * t184 + Icges(7,6) * t183 + Icges(7,3) * t208;
t125 = Icges(7,4) * t184 + Icges(7,2) * t183 + Icges(7,6) * t208;
t127 = Icges(7,1) * t184 + Icges(7,4) * t183 + Icges(7,5) * t208;
t50 = t123 * t191 + t125 * t155 + t127 * t156;
t1 = t191 * t34 + t193 * t35 + t208 * t50;
t100 = Icges(6,4) * t156 + Icges(6,2) * t155 + Icges(6,6) * t191;
t104 = Icges(6,1) * t156 + Icges(6,4) * t155 + Icges(6,5) * t191;
t96 = Icges(6,5) * t156 + Icges(6,6) * t155 + Icges(6,3) * t191;
t36 = t100 * t155 + t104 * t156 + t191 * t96;
t101 = Icges(6,4) * t158 + Icges(6,2) * t157 + Icges(6,6) * t193;
t105 = Icges(6,1) * t158 + Icges(6,4) * t157 + Icges(6,5) * t193;
t97 = Icges(6,5) * t158 + Icges(6,6) * t157 + Icges(6,3) * t193;
t37 = t101 * t155 + t105 * t156 + t191 * t97;
t124 = Icges(6,5) * t184 + Icges(6,6) * t183 + Icges(6,3) * t208;
t126 = Icges(6,4) * t184 + Icges(6,2) * t183 + Icges(6,6) * t208;
t128 = Icges(6,1) * t184 + Icges(6,4) * t183 + Icges(6,5) * t208;
t51 = t124 * t191 + t126 * t155 + t128 * t156;
t2 = t191 * t36 + t193 * t37 + t208 * t51;
t291 = t1 / 0.2e1 + t2 / 0.2e1;
t38 = t102 * t158 + t157 * t98 + t193 * t94;
t39 = t103 * t158 + t157 * t99 + t193 * t95;
t52 = t123 * t193 + t125 * t157 + t127 * t158;
t3 = t191 * t38 + t193 * t39 + t208 * t52;
t40 = t100 * t157 + t104 * t158 + t193 * t96;
t41 = t101 * t157 + t105 * t158 + t193 * t97;
t53 = t124 * t193 + t126 * t157 + t128 * t158;
t4 = t191 * t40 + t193 * t41 + t208 * t53;
t290 = t3 / 0.2e1 + t4 / 0.2e1;
t5 = t210 * t34 + t212 * t35 + t217 * t50;
t6 = t210 * t36 + t212 * t37 + t217 * t51;
t289 = t6 / 0.2e1 + t5 / 0.2e1;
t7 = t210 * t38 + t212 * t39 + t217 * t52;
t8 = t210 * t40 + t212 * t41 + t217 * t53;
t288 = t8 / 0.2e1 + t7 / 0.2e1;
t10 = t219 * t36 + t220 * t37 + t227 * t51;
t9 = t219 * t34 + t220 * t35 + t227 * t50;
t286 = t10 / 0.2e1 + t9 / 0.2e1;
t11 = t219 * t38 + t220 * t39 + t227 * t52;
t12 = t219 * t40 + t220 * t41 + t227 * t53;
t285 = t12 / 0.2e1 + t11 / 0.2e1;
t43 = t102 * t184 + t183 * t98 + t208 * t94;
t44 = t103 * t184 + t183 * t99 + t208 * t95;
t62 = t208 * t123 + t183 * t125 + t184 * t127;
t56 = t62 * t208;
t13 = t43 * t191 + t44 * t193 + t56;
t45 = t100 * t183 + t104 * t184 + t208 * t96;
t46 = t101 * t183 + t105 * t184 + t208 * t97;
t63 = t208 * t124 + t183 * t126 + t184 * t128;
t57 = t63 * t208;
t14 = t45 * t191 + t46 * t193 + t57;
t284 = t13 / 0.2e1 + t14 / 0.2e1;
t58 = t62 * t217;
t15 = t43 * t210 + t44 * t212 + t58;
t59 = t63 * t217;
t16 = t45 * t210 + t46 * t212 + t59;
t283 = t16 / 0.2e1 + t15 / 0.2e1;
t60 = t62 * t227;
t17 = t43 * t219 + t44 * t220 + t60;
t61 = t63 * t227;
t18 = t45 * t219 + t46 * t220 + t61;
t282 = t18 / 0.2e1 + t17 / 0.2e1;
t281 = -t147 - t312;
t280 = -t148 - t311;
t279 = -t176 - t297;
t159 = Icges(5,5) * t209 - Icges(5,6) * t208 + Icges(5,3) * t217;
t160 = Icges(5,4) * t209 - Icges(5,2) * t208 + Icges(5,6) * t217;
t161 = Icges(5,1) * t209 - Icges(5,4) * t208 + Icges(5,5) * t217;
t85 = t217 * t159 - t208 * t160 + t209 * t161;
t195 = Icges(4,5) * t218 - Icges(4,6) * t217 + Icges(4,3) * t227;
t196 = Icges(4,4) * t218 - Icges(4,2) * t217 + Icges(4,6) * t227;
t197 = Icges(4,1) * t218 - Icges(4,4) * t217 + Icges(4,5) * t227;
t278 = t227 * t195 - t217 * t196 + t218 * t197;
t138 = t194 * rSges(5,1) - t193 * rSges(5,2) + t212 * rSges(5,3);
t171 = t213 * rSges(4,1) - t212 * rSges(4,2) + t220 * rSges(4,3);
t276 = -t249 * pkin(1) + qJ(2) * t302;
t265 = t43 / 0.2e1 + t45 / 0.2e1 + t50 / 0.2e1 + t51 / 0.2e1;
t264 = t46 / 0.2e1 + t53 / 0.2e1 + t52 / 0.2e1 + t44 / 0.2e1;
t170 = rSges(4,1) * t211 - rSges(4,2) * t210 + rSges(4,3) * t219;
t137 = rSges(5,1) * t192 - rSges(5,2) * t191 + rSges(5,3) * t210;
t261 = -t228 * pkin(2) - t219 * pkin(9) + t276;
t131 = Icges(5,5) * t192 - Icges(5,6) * t191 + Icges(5,3) * t210;
t133 = Icges(5,4) * t192 - Icges(5,2) * t191 + Icges(5,6) * t210;
t135 = Icges(5,1) * t192 - Icges(5,4) * t191 + Icges(5,5) * t210;
t69 = t131 * t217 - t133 * t208 + t135 * t209;
t78 = t159 * t210 - t160 * t191 + t161 * t192;
t260 = t78 / 0.2e1 + t69 / 0.2e1 + t265;
t132 = Icges(5,5) * t194 - Icges(5,6) * t193 + Icges(5,3) * t212;
t134 = Icges(5,4) * t194 - Icges(5,2) * t193 + Icges(5,6) * t212;
t136 = Icges(5,1) * t194 - Icges(5,4) * t193 + Icges(5,5) * t212;
t70 = t132 * t217 - t134 * t208 + t136 * t209;
t79 = t159 * t212 - t160 * t193 + t161 * t194;
t259 = t79 / 0.2e1 + t70 / 0.2e1 + t264;
t254 = -t177 + t261;
t253 = t229 * pkin(2) + t220 * pkin(9) + t292;
t252 = t178 + t253;
t235 = rSges(2,1) * t251 - t249 * rSges(2,2);
t234 = -t249 * rSges(2,1) - rSges(2,2) * t251;
t207 = t229 * rSges(3,1) - rSges(3,2) * t262 + rSges(3,3) * t303 + t292;
t206 = -t228 * rSges(3,1) + rSges(3,2) * t263 + rSges(3,3) * t302 + t276;
t198 = rSges(4,1) * t218 - rSges(4,2) * t217 + rSges(4,3) * t227;
t169 = Icges(4,1) * t213 - Icges(4,4) * t212 + Icges(4,5) * t220;
t168 = Icges(4,1) * t211 - Icges(4,4) * t210 + Icges(4,5) * t219;
t167 = Icges(4,4) * t213 - Icges(4,2) * t212 + Icges(4,6) * t220;
t166 = Icges(4,4) * t211 - Icges(4,2) * t210 + Icges(4,6) * t219;
t165 = Icges(4,5) * t213 - Icges(4,6) * t212 + Icges(4,3) * t220;
t164 = Icges(4,5) * t211 - Icges(4,6) * t210 + Icges(4,3) * t219;
t162 = rSges(5,1) * t209 - rSges(5,2) * t208 + rSges(5,3) * t217;
t150 = t210 * t176;
t145 = t253 + t171;
t144 = -t170 + t261;
t140 = t217 * t148;
t139 = t212 * t147;
t121 = t171 * t227 - t198 * t220;
t120 = -t170 * t227 + t198 * t219;
t113 = t170 * t220 - t171 * t219;
t112 = t278 * t227;
t111 = t195 * t220 - t196 * t212 + t197 * t213;
t110 = t195 * t219 - t196 * t210 + t197 * t211;
t93 = t252 + t138;
t92 = -t137 + t254;
t89 = t138 * t217 - t162 * t212;
t88 = -t137 * t217 + t162 * t210;
t87 = t165 * t227 - t167 * t217 + t169 * t218;
t86 = t164 * t227 - t166 * t217 + t168 * t218;
t84 = t137 * t212 - t138 * t210;
t83 = t85 * t227;
t82 = t85 * t217;
t81 = t227 * t138 + t173 + (-t162 - t199) * t220;
t80 = t219 * t162 + t181 + (-t137 - t177) * t227;
t77 = t252 - t298;
t76 = t254 + t299;
t75 = t109 * t208 - t130 * t193;
t74 = -t107 * t208 + t130 * t191;
t73 = t252 + t322;
t72 = -t192 * t239 + t254 - t324;
t71 = t220 * t137 + t172 + (-t138 - t178) * t219;
t68 = t132 * t212 - t134 * t193 + t136 * t194;
t67 = t131 * t212 - t133 * t193 + t135 * t194;
t66 = t132 * t210 - t134 * t191 + t136 * t192;
t65 = t131 * t210 - t133 * t191 + t135 * t192;
t64 = t107 * t193 - t109 * t191;
t55 = t217 * t109 + t212 * t296 + t140;
t54 = t210 * t130 + t217 * t299 + t150;
t49 = t227 * t109 + (-t199 + t296) * t220 + t294;
t48 = t219 * t130 + (-t177 + t299) * t227 + t293;
t47 = t212 * t107 + t210 * t298 + t139;
t42 = t220 * t107 + (-t178 + t298) * t219 + t295;
t33 = -t193 * t297 + t208 * t311;
t32 = t191 * t297 - t208 * t312;
t31 = t212 * t279 + t217 * t311 + t140;
t30 = t210 * t297 + t217 * t281 + t150;
t29 = t311 * t227 + (-t199 + t279) * t220 + t294;
t28 = t297 * t219 + (-t177 + t281) * t227 + t293;
t27 = -t191 * t311 + t193 * t312;
t26 = t210 * t280 + t212 * t312 + t139;
t25 = t69 * t219 + t70 * t220 + t83;
t24 = t69 * t210 + t70 * t212 + t82;
t23 = t312 * t220 + (-t178 + t280) * t219 + t295;
t22 = t219 * t67 + t220 * t68 + t227 * t79;
t21 = t219 * t65 + t220 * t66 + t227 * t78;
t20 = t210 * t67 + t212 * t68 + t217 * t79;
t19 = t210 * t65 + t212 * t66 + t217 * t78;
t90 = [t244 * (Icges(3,3) * t244 + (Icges(3,5) * t242 + Icges(3,6) * t309) * t243) + m(6) * (t76 ^ 2 + t77 ^ 2) + m(5) * (t92 ^ 2 + t93 ^ 2) + m(4) * (t144 ^ 2 + t145 ^ 2) + m(3) * (t206 ^ 2 + t207 ^ 2) + m(2) * (t234 ^ 2 + t235 ^ 2) + t243 * t309 * (Icges(3,6) * t244 + (Icges(3,4) * t242 + Icges(3,2) * t309) * t243) + m(7) * (t72 ^ 2 + t73 ^ 2) + Icges(2,3) + t278 + (Icges(3,5) * t244 + (Icges(3,1) * t242 + Icges(3,4) * t309) * t243) * t304 + t85 + t63 + t62; 0.2e1 * ((t249 * t72 - t251 * t73) * t316 + (t249 * t76 - t251 * t77) * t317 + (t249 * t92 - t251 * t93) * t318 + (t144 * t249 - t145 * t251) * t319 + (t206 * t249 - t207 * t251) * t320) * t243; 0.2e1 * (t320 + t319 + t318 + t317 + t316) * (t244 ^ 2 + (t249 ^ 2 + t251 ^ 2) * t243 ^ 2); t61 + t60 + t83 + t112 + m(7) * (t28 * t72 + t29 * t73) + m(6) * (t48 * t76 + t49 * t77) + m(5) * (t80 * t92 + t81 * t93) + m(4) * (t120 * t144 + t121 * t145) + (t87 / 0.2e1 + t111 / 0.2e1 + t259) * t220 + (t110 / 0.2e1 + t86 / 0.2e1 + t260) * t219; m(4) * (t113 * t244 + (t120 * t249 - t121 * t251) * t243) + m(5) * (t71 * t244 + (t249 * t80 - t251 * t81) * t243) + m(6) * (t42 * t244 + (t249 * t48 - t251 * t49) * t243) + m(7) * (t23 * t244 + (t249 * t28 - t251 * t29) * t243); (t112 + t17 + t18 + t25) * t227 + (t11 + t12 + t22 + (t165 * t220 - t167 * t212 + t169 * t213) * t220 + (t111 + t87) * t227) * t220 + (t9 + t10 + t21 + (t164 * t219 - t166 * t210 + t168 * t211) * t219 + (t86 + t110) * t227 + (t164 * t220 + t165 * t219 - t166 * t212 - t167 * t210 + t168 * t213 + t169 * t211) * t220) * t219 + m(7) * (t23 ^ 2 + t28 ^ 2 + t29 ^ 2) + m(6) * (t42 ^ 2 + t48 ^ 2 + t49 ^ 2) + m(5) * (t71 ^ 2 + t80 ^ 2 + t81 ^ 2) + m(4) * (t113 ^ 2 + t120 ^ 2 + t121 ^ 2); t59 + t58 + t82 + m(7) * (t30 * t72 + t31 * t73) + m(6) * (t54 * t76 + t55 * t77) + m(5) * (t88 * t92 + t89 * t93) + t259 * t212 + t260 * t210; m(5) * (t84 * t244 + (t249 * t88 - t251 * t89) * t243) + m(6) * (t47 * t244 + (t249 * t54 - t251 * t55) * t243) + m(7) * (t26 * t244 + (t249 * t30 - t251 * t31) * t243); (t24 / 0.2e1 + t283) * t227 + (t20 / 0.2e1 + t288) * t220 + (t19 / 0.2e1 + t289) * t219 + (t25 / 0.2e1 + t282) * t217 + (t22 / 0.2e1 + t285) * t212 + (t21 / 0.2e1 + t286) * t210 + m(7) * (t23 * t26 + t28 * t30 + t29 * t31) + m(6) * (t42 * t47 + t48 * t54 + t49 * t55) + m(5) * (t71 * t84 + t80 * t88 + t81 * t89); (t16 + t15 + t24) * t217 + (t8 + t7 + t20) * t212 + (t5 + t6 + t19) * t210 + m(7) * (t26 ^ 2 + t30 ^ 2 + t31 ^ 2) + m(6) * (t47 ^ 2 + t54 ^ 2 + t55 ^ 2) + m(5) * (t84 ^ 2 + t88 ^ 2 + t89 ^ 2); t57 + t56 + m(7) * (t32 * t72 + t33 * t73) + m(6) * (t74 * t76 + t75 * t77) + t264 * t193 + t265 * t191; m(6) * (t64 * t244 + (t249 * t74 - t251 * t75) * t243) + m(7) * (t27 * t244 + (t249 * t32 - t251 * t33) * t243); t284 * t227 + t290 * t220 + t291 * t219 + t282 * t208 + t285 * t193 + t286 * t191 + m(7) * (t23 * t27 + t28 * t32 + t29 * t33) + m(6) * (t64 * t42 + t48 * t74 + t49 * t75); t284 * t217 + t290 * t212 + t291 * t210 + t283 * t208 + t288 * t193 + t289 * t191 + m(7) * (t26 * t27 + t30 * t32 + t31 * t33) + m(6) * (t64 * t47 + t54 * t74 + t55 * t75); (t14 + t13) * t208 + (t4 + t3) * t193 + (t1 + t2) * t191 + m(7) * (t27 ^ 2 + t32 ^ 2 + t33 ^ 2) + m(6) * (t64 ^ 2 + t74 ^ 2 + t75 ^ 2); m(7) * (t191 * t73 + t193 * t72); m(7) * (t208 * t244 + (-t191 * t251 + t193 * t249) * t243); m(7) * (t191 * t29 + t193 * t28 + t208 * t23); m(7) * (t191 * t31 + t193 * t30 + t208 * t26); m(7) * (t191 * t33 + t193 * t32 + t208 * t27); m(7) * (t191 ^ 2 + t193 ^ 2 + t208 ^ 2);];
%% Postprocessing: Reshape Output
% From vec2symmat_6_matlab.m
res = [t90(1) t90(2) t90(4) t90(7) t90(11) t90(16); t90(2) t90(3) t90(5) t90(8) t90(12) t90(17); t90(4) t90(5) t90(6) t90(9) t90(13) t90(18); t90(7) t90(8) t90(9) t90(10) t90(14) t90(19); t90(11) t90(12) t90(13) t90(14) t90(15) t90(20); t90(16) t90(17) t90(18) t90(19) t90(20) t90(21);];
Mq  = res;
