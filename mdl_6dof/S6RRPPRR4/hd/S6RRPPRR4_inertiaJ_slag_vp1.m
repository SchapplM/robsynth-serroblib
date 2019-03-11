% Calculate joint inertia matrix for
% S6RRPPRR4
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% pkin [11x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,alpha2,d1,d2,d5,d6,theta3]';
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
% Datum: 2019-03-09 09:06
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Mq = S6RRPPRR4_inertiaJ_slag_vp1(qJ, ...
  pkin, m, rSges, Icges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(11,1),zeros(7,1),zeros(7,3),zeros(7,6)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRPPRR4_inertiaJ_slag_vp1: qJ has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [11 1]), ...
  'S6RRPPRR4_inertiaJ_slag_vp1: pkin has to be [11x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6RRPPRR4_inertiaJ_slag_vp1: m has to be [7x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [7,3]), ...
  'S6RRPPRR4_inertiaJ_slag_vp1: rSges has to be [7x3] (double)');
assert(isreal(Icges) && all(size(Icges) == [7 6]), ...
  'S6RRPPRR4_inertiaJ_slag_vp1: Icges has to be [7x6] (double)'); 

%% Symbolic Calculation
% From inertia_joint_joint_floatb_twist_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 09:01:15
% EndTime: 2019-03-09 09:01:28
% DurationCPUTime: 5.52s
% Computational Cost: add. (16997->563), mult. (43314->797), div. (0->0), fcn. (56244->12), ass. (0->258)
t300 = sin(pkin(11));
t301 = cos(pkin(11));
t309 = sin(qJ(2));
t311 = cos(qJ(2));
t232 = -t309 * t300 + t311 * t301;
t246 = sin(pkin(6));
t221 = t232 * t246;
t256 = t300 * t311 + t301 * t309;
t222 = t256 * t246;
t247 = cos(pkin(6));
t173 = Icges(5,5) * t247 - Icges(5,6) * t222 - Icges(5,3) * t221;
t174 = Icges(4,5) * t222 + Icges(4,6) * t221 + Icges(4,3) * t247;
t175 = Icges(5,4) * t247 - Icges(5,2) * t222 - Icges(5,6) * t221;
t176 = Icges(4,4) * t222 + Icges(4,2) * t221 + Icges(4,6) * t247;
t177 = Icges(5,1) * t247 - Icges(5,4) * t222 - Icges(5,5) * t221;
t178 = Icges(4,1) * t222 + Icges(4,4) * t221 + Icges(4,5) * t247;
t217 = Icges(3,3) * t247 + (Icges(3,5) * t309 + Icges(3,6) * t311) * t246;
t218 = Icges(3,6) * t247 + (Icges(3,4) * t309 + Icges(3,2) * t311) * t246;
t219 = Icges(3,5) * t247 + (Icges(3,1) * t309 + Icges(3,4) * t311) * t246;
t279 = t246 * t309;
t327 = t246 * t311 * t218 + t219 * t279 + (-t175 + t178) * t222 + (-t173 + t176) * t221 + (t174 + t177 + t217) * t247;
t317 = m(7) / 0.2e1;
t318 = m(6) / 0.2e1;
t319 = m(5) / 0.2e1;
t271 = t319 + t318 + t317;
t326 = 0.2e1 * t271;
t250 = sin(qJ(1));
t252 = cos(qJ(1));
t254 = t247 * t232;
t205 = -t250 * t254 - t252 * t256;
t253 = t247 * t256;
t207 = t232 * t252 - t250 * t253;
t298 = t246 * t250;
t119 = Icges(5,5) * t298 - Icges(5,6) * t207 - Icges(5,3) * t205;
t128 = Icges(4,4) * t207 + Icges(4,2) * t205 + Icges(4,6) * t298;
t325 = t119 - t128;
t202 = -t250 * t256 + t252 * t254;
t204 = t250 * t232 + t252 * t253;
t297 = t246 * t252;
t120 = -Icges(5,5) * t297 - Icges(5,6) * t204 - Icges(5,3) * t202;
t127 = Icges(4,4) * t204 + Icges(4,2) * t202 - Icges(4,6) * t297;
t324 = t120 - t127;
t121 = Icges(5,4) * t298 - Icges(5,2) * t207 - Icges(5,6) * t205;
t130 = Icges(4,1) * t207 + Icges(4,4) * t205 + Icges(4,5) * t298;
t323 = t121 - t130;
t122 = -Icges(5,4) * t297 - Icges(5,2) * t204 - Icges(5,6) * t202;
t129 = Icges(4,1) * t204 + Icges(4,4) * t202 - Icges(4,5) * t297;
t322 = t122 - t129;
t321 = t246 ^ 2;
t320 = m(4) / 0.2e1;
t249 = sin(qJ(5));
t310 = cos(qJ(5));
t163 = t205 * t310 + t249 * t298;
t165 = -t202 * t310 + t249 * t297;
t280 = t246 * t310;
t164 = -t205 * t249 + t250 * t280;
t248 = sin(qJ(6));
t251 = cos(qJ(6));
t114 = -t164 * t248 + t207 * t251;
t115 = t164 * t251 + t207 * t248;
t66 = Icges(7,5) * t115 + Icges(7,6) * t114 + Icges(7,3) * t163;
t68 = Icges(7,4) * t115 + Icges(7,2) * t114 + Icges(7,6) * t163;
t70 = Icges(7,1) * t115 + Icges(7,4) * t114 + Icges(7,5) * t163;
t17 = t114 * t68 + t115 * t70 + t163 * t66;
t166 = -t202 * t249 - t252 * t280;
t116 = -t166 * t248 + t204 * t251;
t117 = t166 * t251 + t204 * t248;
t67 = Icges(7,5) * t117 + Icges(7,6) * t116 - Icges(7,3) * t165;
t69 = Icges(7,4) * t117 + Icges(7,2) * t116 - Icges(7,6) * t165;
t71 = Icges(7,1) * t117 + Icges(7,4) * t116 - Icges(7,5) * t165;
t18 = t114 * t69 + t115 * t71 + t163 * t67;
t208 = t221 * t310 + t247 * t249;
t209 = -t221 * t249 + t247 * t310;
t152 = -t209 * t248 + t222 * t251;
t153 = t209 * t251 + t222 * t248;
t95 = Icges(7,5) * t153 + Icges(7,6) * t152 + Icges(7,3) * t208;
t96 = Icges(7,4) * t153 + Icges(7,2) * t152 + Icges(7,6) * t208;
t97 = Icges(7,1) * t153 + Icges(7,4) * t152 + Icges(7,5) * t208;
t29 = t114 * t96 + t115 * t97 + t163 * t95;
t1 = t163 * t17 - t165 * t18 + t208 * t29;
t316 = t1 / 0.2e1;
t21 = t152 * t68 + t153 * t70 + t208 * t66;
t22 = t152 * t69 + t153 * t71 + t208 * t67;
t39 = t152 * t96 + t153 * t97 + t208 * t95;
t32 = t39 * t208;
t7 = t21 * t163 - t22 * t165 + t32;
t315 = t7 / 0.2e1;
t314 = t163 / 0.2e1;
t313 = -t165 / 0.2e1;
t312 = t208 / 0.2e1;
t308 = pkin(1) * t252;
t307 = t166 * pkin(5);
t306 = t204 * pkin(3);
t305 = t202 * rSges(5,3);
t72 = t115 * rSges(7,1) + t114 * rSges(7,2) + t163 * rSges(7,3);
t304 = t164 * pkin(5) + pkin(10) * t163 + t72;
t264 = -t117 * rSges(7,1) - t116 * rSges(7,2);
t73 = -t165 * rSges(7,3) - t264;
t303 = -t165 * pkin(10) + t307 + t73;
t98 = rSges(7,1) * t153 + rSges(7,2) * t152 + rSges(7,3) * t208;
t302 = pkin(5) * t209 + pkin(10) * t208 + t98;
t223 = t247 * t309 * pkin(2) + (-pkin(8) - qJ(3)) * t246;
t299 = t223 * t252;
t244 = pkin(2) * t311 + pkin(1);
t296 = t250 * t244;
t143 = t207 * pkin(3) - qJ(4) * t205;
t237 = t252 * t244;
t193 = -t308 + t237 + (-t246 * pkin(8) - t223) * t250;
t181 = t247 * t193;
t295 = t247 * t143 + t181;
t194 = t202 * qJ(4);
t142 = -t194 + t306;
t242 = pkin(8) * t297;
t192 = t299 + t242 + (-pkin(1) + t244) * t250;
t294 = -t142 - t192;
t293 = t192 * t298 + t193 * t297;
t233 = pkin(2) * t279 + t247 * qJ(3);
t292 = -pkin(3) * t222 + qJ(4) * t221 - t233;
t170 = -pkin(4) * t297 + t204 * pkin(9);
t291 = t327 * t247;
t290 = t21 / 0.2e1 + t29 / 0.2e1;
t30 = t116 * t96 + t117 * t97 - t165 * t95;
t289 = -t22 / 0.2e1 - t30 / 0.2e1;
t138 = Icges(6,5) * t209 - Icges(6,6) * t208 + Icges(6,3) * t222;
t139 = Icges(6,4) * t209 - Icges(6,2) * t208 + Icges(6,6) * t222;
t140 = Icges(6,1) * t209 - Icges(6,4) * t208 + Icges(6,5) * t222;
t56 = t222 * t138 - t208 * t139 + t209 * t140;
t123 = Icges(5,1) * t298 - Icges(5,4) * t207 - Icges(5,5) * t205;
t126 = Icges(4,5) * t207 + Icges(4,6) * t205 + Icges(4,3) * t298;
t275 = t252 * t309;
t278 = t250 * t311;
t229 = -t247 * t278 - t275;
t276 = t252 * t311;
t277 = t250 * t309;
t230 = -t247 * t277 + t276;
t183 = Icges(3,5) * t230 + Icges(3,6) * t229 + Icges(3,3) * t298;
t288 = t123 + t126 + t183;
t124 = -Icges(5,1) * t297 - Icges(5,4) * t204 - Icges(5,5) * t202;
t125 = Icges(4,5) * t204 + Icges(4,6) * t202 - Icges(4,3) * t297;
t227 = t247 * t276 - t277;
t228 = t247 * t275 + t278;
t182 = Icges(3,5) * t228 + Icges(3,6) * t227 - Icges(3,3) * t297;
t287 = -t125 - t182 - t124;
t169 = pkin(4) * t298 + pkin(9) * t207;
t286 = t247 * t169 + t295;
t285 = -t170 + t294;
t93 = t164 * rSges(6,1) - t163 * rSges(6,2) + t207 * rSges(6,3);
t282 = -pkin(4) * t247 - pkin(9) * t222 + t292;
t134 = t207 * rSges(4,1) + t205 * rSges(4,2) + rSges(4,3) * t298;
t190 = t230 * rSges(3,1) + t229 * rSges(3,2) + rSges(3,3) * t298;
t131 = rSges(5,1) * t298 - t207 * rSges(5,2) - t205 * rSges(5,3);
t274 = t246 * (-rSges(4,1) * t222 - rSges(4,2) * t221 - rSges(4,3) * t247 - t233);
t273 = t194 - t296;
t272 = -t250 * t223 + t237;
t270 = t142 * t298 + t143 * t297 + t293;
t269 = t246 * (-rSges(5,1) * t247 + rSges(5,2) * t222 + rSges(5,3) * t221 + t292);
t266 = -t204 * rSges(4,1) - t202 * rSges(4,2);
t265 = -t166 * rSges(6,1) - t165 * rSges(6,2);
t141 = rSges(6,1) * t209 - rSges(6,2) * t208 + rSges(6,3) * t222;
t263 = t246 * (-t141 + t282);
t87 = Icges(6,5) * t164 - Icges(6,6) * t163 + Icges(6,3) * t207;
t89 = Icges(6,4) * t164 - Icges(6,2) * t163 + Icges(6,6) * t207;
t91 = Icges(6,1) * t164 - Icges(6,4) * t163 + Icges(6,5) * t207;
t40 = -t208 * t89 + t209 * t91 + t222 * t87;
t49 = t138 * t207 - t139 * t163 + t140 * t164;
t262 = t49 / 0.2e1 + t40 / 0.2e1 + t290;
t88 = Icges(6,5) * t166 + Icges(6,6) * t165 + Icges(6,3) * t204;
t90 = Icges(6,4) * t166 + Icges(6,2) * t165 + Icges(6,6) * t204;
t92 = Icges(6,1) * t166 + Icges(6,4) * t165 + Icges(6,5) * t204;
t41 = -t208 * t90 + t209 * t92 + t222 * t88;
t50 = t138 * t204 + t139 * t165 + t140 * t166;
t261 = t50 / 0.2e1 + t41 / 0.2e1 - t289;
t260 = t169 * t297 + t170 * t298 + t270;
t259 = t246 * (t282 - t302);
t258 = t143 + t272;
t257 = -t170 + t273 - t299;
t189 = t228 * rSges(3,1) + t227 * rSges(3,2) - rSges(3,3) * t297;
t255 = t169 + t258;
t235 = rSges(2,1) * t252 - t250 * rSges(2,2);
t234 = -t250 * rSges(2,1) - rSges(2,2) * t252;
t220 = t247 * rSges(3,3) + (rSges(3,1) * t309 + rSges(3,2) * t311) * t246;
t187 = Icges(3,1) * t230 + Icges(3,4) * t229 + Icges(3,5) * t298;
t186 = Icges(3,1) * t228 + Icges(3,4) * t227 - Icges(3,5) * t297;
t185 = Icges(3,4) * t230 + Icges(3,2) * t229 + Icges(3,6) * t298;
t184 = Icges(3,4) * t228 + Icges(3,2) * t227 - Icges(3,6) * t297;
t161 = pkin(8) * t298 + t190 + t308;
t160 = -t250 * pkin(1) - t189 + t242;
t146 = -t247 * t189 - t220 * t297;
t145 = t190 * t247 - t220 * t298;
t133 = -rSges(4,3) * t297 - t266;
t132 = -rSges(5,1) * t297 - t204 * rSges(5,2) - t305;
t118 = (t189 * t250 + t190 * t252) * t246;
t113 = t217 * t298 + t218 * t229 + t219 * t230;
t112 = -t217 * t297 + t227 * t218 + t228 * t219;
t102 = t272 + t134;
t101 = -t296 + (rSges(4,3) * t246 - t223) * t252 + t266;
t100 = t247 * t183 + (t185 * t311 + t187 * t309) * t246;
t99 = t247 * t182 + (t184 * t311 + t186 * t309) * t246;
t94 = t204 * rSges(6,3) - t265;
t83 = t258 + t131;
t82 = t305 + (rSges(5,1) * t246 - t223) * t252 + (rSges(5,2) - pkin(3)) * t204 + t273;
t79 = (-t133 - t192) * t247 + t252 * t274;
t78 = t134 * t247 + t250 * t274 + t181;
t77 = t174 * t298 + t176 * t205 + t178 * t207;
t76 = -t174 * t297 + t202 * t176 + t204 * t178;
t75 = -t202 * t173 - t204 * t175 - t177 * t297;
t74 = -t173 * t205 - t175 * t207 + t177 * t298;
t65 = (t133 * t250 + t134 * t252) * t246 + t293;
t64 = -t141 * t207 + t222 * t93;
t63 = t141 * t204 - t222 * t94;
t62 = t255 + t93;
t61 = (-rSges(6,3) - pkin(3)) * t204 + t257 + t265;
t60 = t126 * t247 + t128 * t221 + t130 * t222;
t59 = t125 * t247 + t127 * t221 + t129 * t222;
t58 = -t120 * t221 - t122 * t222 + t124 * t247;
t57 = -t119 * t221 - t121 * t222 + t123 * t247;
t55 = t56 * t247;
t54 = (-t132 + t294) * t247 + t252 * t269;
t53 = t131 * t247 + t250 * t269 + t295;
t52 = t56 * t222;
t51 = -t204 * t93 + t207 * t94;
t48 = (t131 * t252 + t132 * t250) * t246 + t270;
t47 = -t165 * t98 - t208 * t73;
t46 = -t163 * t98 + t208 * t72;
t45 = t255 + t304;
t44 = -t306 - t307 + (rSges(7,3) + pkin(10)) * t165 + t257 + t264;
t43 = (-t94 + t285) * t247 + t252 * t263;
t42 = t247 * t93 + t250 * t263 + t286;
t38 = t39 * t247;
t37 = t165 * t90 + t166 * t92 + t204 * t88;
t36 = t165 * t89 + t166 * t91 + t204 * t87;
t35 = -t163 * t90 + t164 * t92 + t207 * t88;
t34 = -t163 * t89 + t164 * t91 + t207 * t87;
t33 = t39 * t222;
t31 = t163 * t73 + t165 * t72;
t28 = -t207 * t302 + t222 * t304;
t27 = t204 * t302 - t222 * t303;
t26 = (t250 * t94 + t252 * t93) * t246 + t260;
t25 = (t285 - t303) * t247 + t252 * t259;
t24 = t247 * t304 + t250 * t259 + t286;
t23 = -t204 * t304 + t207 * t303;
t20 = t116 * t69 + t117 * t71 - t165 * t67;
t19 = t116 * t68 + t117 * t70 - t165 * t66;
t16 = (t250 * t303 + t252 * t304) * t246 + t260;
t15 = t55 + (t40 * t250 - t41 * t252) * t246;
t14 = t41 * t204 + t40 * t207 + t52;
t13 = t50 * t247 + (t250 * t36 - t252 * t37) * t246;
t12 = t49 * t247 + (t250 * t34 - t252 * t35) * t246;
t11 = t204 * t37 + t207 * t36 + t222 * t50;
t10 = t204 * t35 + t207 * t34 + t222 * t49;
t9 = t38 + (t21 * t250 - t22 * t252) * t246;
t8 = t22 * t204 + t21 * t207 + t33;
t6 = t30 * t247 + (t19 * t250 - t20 * t252) * t246;
t5 = t29 * t247 + (t17 * t250 - t18 * t252) * t246;
t4 = t19 * t207 + t20 * t204 + t222 * t30;
t3 = t17 * t207 + t18 * t204 + t222 * t29;
t2 = t163 * t19 - t165 * t20 + t208 * t30;
t80 = [m(7) * (t44 ^ 2 + t45 ^ 2) + m(6) * (t61 ^ 2 + t62 ^ 2) + m(5) * (t82 ^ 2 + t83 ^ 2) + m(4) * (t101 ^ 2 + t102 ^ 2) + m(3) * (t160 ^ 2 + t161 ^ 2) + m(2) * (t234 ^ 2 + t235 ^ 2) + Icges(2,3) + t56 + t39 + t327; t55 + t38 + m(6) * (t42 * t62 + t43 * t61) + m(5) * (t53 * t83 + t54 * t82) + m(4) * (t101 * t79 + t102 * t78) + m(3) * (t145 * t161 + t146 * t160) + m(7) * (t24 * t45 + t25 * t44) + ((-t99 / 0.2e1 - t59 / 0.2e1 - t58 / 0.2e1 - t75 / 0.2e1 - t76 / 0.2e1 - t112 / 0.2e1 - t261) * t252 + (t100 / 0.2e1 + t60 / 0.2e1 + t57 / 0.2e1 + t74 / 0.2e1 + t77 / 0.2e1 + t113 / 0.2e1 + t262) * t250) * t246 + t291; (t9 + t15 + t291) * t247 + m(7) * (t16 ^ 2 + t24 ^ 2 + t25 ^ 2) + m(6) * (t26 ^ 2 + t42 ^ 2 + t43 ^ 2) + m(5) * (t48 ^ 2 + t53 ^ 2 + t54 ^ 2) + m(4) * (t65 ^ 2 + t78 ^ 2 + t79 ^ 2) + m(3) * (t118 ^ 2 + t145 ^ 2 + t146 ^ 2) + ((-t13 - t6 + ((t227 * t184 + t228 * t186 - t202 * t324 - t204 * t322) * t246 + t287 * t321 * t252) * t252 + (-t112 - t58 - t59 - t75 - t76 - t99) * t247) * t252 + (t5 + t12 + ((t185 * t229 + t187 * t230 - t205 * t325 - t207 * t323) * t246 + t288 * t321 * t250) * t250 + (t74 + t77 + t113 + t100 + t60 + t57) * t247 + (-t184 * t229 - t227 * t185 - t186 * t230 - t228 * t187 + (t250 * t287 + t252 * t288) * t246 + t322 * t207 + t324 * t205 + t323 * t204 + t325 * t202) * t297) * t250) * t246; 0.2e1 * ((t250 * t44 - t252 * t45) * t317 + (t250 * t61 - t252 * t62) * t318 + (t250 * t82 - t252 * t83) * t319 + (t101 * t250 - t102 * t252) * t320) * t246; m(7) * (t247 * t16 + (-t24 * t252 + t25 * t250) * t246) + m(6) * (t247 * t26 + (t250 * t43 - t252 * t42) * t246) + m(5) * (t247 * t48 + (t250 * t54 - t252 * t53) * t246) + m(4) * (t247 * t65 + (t250 * t79 - t252 * t78) * t246); 0.2e1 * (t320 + t271) * (t247 ^ 2 + (t250 ^ 2 + t252 ^ 2) * t321); m(7) * (-t202 * t45 - t205 * t44) + m(6) * (-t202 * t62 - t205 * t61) + m(5) * (-t202 * t83 - t205 * t82); m(7) * (-t16 * t221 - t202 * t24 - t205 * t25) + m(6) * (-t202 * t42 - t205 * t43 - t221 * t26) + m(5) * (-t202 * t53 - t205 * t54 - t221 * t48); (-t221 * t247 + (t202 * t252 - t205 * t250) * t246) * t326; (t202 ^ 2 + t205 ^ 2 + t221 ^ 2) * t326; t33 + t52 + m(7) * (t27 * t44 + t28 * t45) + m(6) * (t61 * t63 + t62 * t64) + t262 * t207 + t261 * t204; (t8 / 0.2e1 + t14 / 0.2e1) * t247 + (t9 / 0.2e1 + t15 / 0.2e1) * t222 + (t5 / 0.2e1 + t12 / 0.2e1) * t207 + (t6 / 0.2e1 + t13 / 0.2e1) * t204 + m(7) * (t16 * t23 + t24 * t28 + t25 * t27) + m(6) * (t26 * t51 + t42 * t64 + t43 * t63) + ((-t4 / 0.2e1 - t11 / 0.2e1) * t252 + (t3 / 0.2e1 + t10 / 0.2e1) * t250) * t246; m(6) * (t51 * t247 + (t250 * t63 - t252 * t64) * t246) + m(7) * (t23 * t247 + (t250 * t27 - t252 * t28) * t246); m(6) * (-t202 * t64 - t205 * t63 - t221 * t51) + m(7) * (-t202 * t28 - t205 * t27 - t221 * t23); (t8 + t14) * t222 + (t3 + t10) * t207 + (t4 + t11) * t204 + m(7) * (t23 ^ 2 + t27 ^ 2 + t28 ^ 2) + m(6) * (t51 ^ 2 + t63 ^ 2 + t64 ^ 2); m(7) * (t44 * t47 + t45 * t46) + t32 + t289 * t165 + t290 * t163; t6 * t313 + t5 * t314 + m(7) * (t16 * t31 + t24 * t46 + t25 * t47) + t247 * t315 + t9 * t312 + (t250 * t316 - t252 * t2 / 0.2e1) * t246; m(7) * (t31 * t247 + (t250 * t47 - t252 * t46) * t246); m(7) * (-t202 * t46 - t205 * t47 - t221 * t31); m(7) * (t23 * t31 + t27 * t47 + t28 * t46) + t4 * t313 + t3 * t314 + t222 * t315 + t207 * t316 + t8 * t312 + t204 * t2 / 0.2e1; t163 * t1 - t165 * t2 + t208 * t7 + m(7) * (t31 ^ 2 + t46 ^ 2 + t47 ^ 2);];
%% Postprocessing: Reshape Output
% From vec2symmat_6_matlab.m
res = [t80(1) t80(2) t80(4) t80(7) t80(11) t80(16); t80(2) t80(3) t80(5) t80(8) t80(12) t80(17); t80(4) t80(5) t80(6) t80(9) t80(13) t80(18); t80(7) t80(8) t80(9) t80(10) t80(14) t80(19); t80(11) t80(12) t80(13) t80(14) t80(15) t80(20); t80(16) t80(17) t80(18) t80(19) t80(20) t80(21);];
Mq  = res;
