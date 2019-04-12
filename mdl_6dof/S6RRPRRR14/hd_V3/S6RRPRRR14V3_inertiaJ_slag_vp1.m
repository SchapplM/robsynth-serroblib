% Calculate joint inertia matrix for
% S6RRPRRR14V3
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% pkin [1x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[dummy]';
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
% Datum: 2019-04-12 15:12
% Revision: b693519ea345eb34ae9622239e7f1167217e9d53 (2019-04-09)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Mq = S6RRPRRR14V3_inertiaJ_slag_vp1(qJ, ...
  pkin, m, rSges, Icges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(1,1),zeros(7,1),zeros(7,3),zeros(7,6)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRPRRR14V3_inertiaJ_slag_vp1: qJ has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [1 1]), ...
  'S6RRPRRR14V3_inertiaJ_slag_vp1: pkin has to be [1x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6RRPRRR14V3_inertiaJ_slag_vp1: m has to be [7x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [7,3]), ...
  'S6RRPRRR14V3_inertiaJ_slag_vp1: rSges has to be [7x3] (double)');
assert(isreal(Icges) && all(size(Icges) == [7 6]), ...
  'S6RRPRRR14V3_inertiaJ_slag_vp1: Icges has to be [7x6] (double)'); 

%% Symbolic Calculation
% From inertia_joint_joint_floatb_twist_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2019-04-12 15:03:28
% EndTime: 2019-04-12 15:03:44
% DurationCPUTime: 5.45s
% Computational Cost: add. (11285->501), mult. (29084->759), div. (0->0), fcn. (36550->10), ass. (0->224)
t283 = Icges(3,1) + Icges(4,1);
t281 = Icges(4,4) + Icges(3,5);
t209 = sin(qJ(2));
t282 = (Icges(3,4) - Icges(4,5)) * t209;
t280 = Icges(4,2) + Icges(3,3);
t213 = cos(qJ(2));
t279 = t281 * t213 + (-Icges(3,6) + Icges(4,6)) * t209;
t278 = t283 * t213 - t282;
t210 = sin(qJ(1));
t277 = -t210 / 0.2e1;
t261 = t210 / 0.2e1;
t214 = cos(qJ(1));
t273 = -t214 / 0.2e1;
t276 = t214 / 0.2e1;
t204 = t210 ^ 2;
t205 = t214 ^ 2;
t243 = t204 + t205;
t274 = -t213 / 0.2e1;
t207 = sin(qJ(5));
t212 = cos(qJ(4));
t260 = cos(qJ(5));
t231 = t209 * t260;
t173 = -t213 * t207 + t212 * t231;
t206 = sin(qJ(6));
t211 = cos(qJ(6));
t208 = sin(qJ(4));
t250 = t208 * t209;
t138 = -t173 * t206 + t211 * t250;
t139 = t173 * t211 + t206 * t250;
t248 = t209 * t212;
t172 = t207 * t248 + t213 * t260;
t91 = Icges(7,5) * t139 + Icges(7,6) * t138 + Icges(7,3) * t172;
t94 = Icges(7,4) * t139 + Icges(7,2) * t138 + Icges(7,6) * t172;
t97 = Icges(7,1) * t139 + Icges(7,4) * t138 + Icges(7,5) * t172;
t37 = t138 * t94 + t139 * t97 + t172 * t91;
t121 = Icges(6,5) * t173 - Icges(6,6) * t172 + Icges(6,3) * t250;
t124 = Icges(6,4) * t173 - Icges(6,2) * t172 + Icges(6,6) * t250;
t127 = Icges(6,1) * t173 - Icges(6,4) * t172 + Icges(6,5) * t250;
t59 = t121 * t250 - t172 * t124 + t173 * t127;
t272 = -t37 - t59;
t271 = t210 * t280 + t214 * t279;
t270 = -t210 * t279 + t214 * t280;
t269 = m(4) / 0.2e1;
t268 = m(5) / 0.2e1;
t267 = m(6) / 0.2e1;
t266 = m(7) / 0.2e1;
t246 = t210 * t213;
t175 = -t208 * t214 + t212 * t246;
t140 = t175 * t207 - t210 * t231;
t245 = t213 * t214;
t177 = t210 * t208 + t212 * t245;
t142 = t177 * t207 - t214 * t231;
t247 = t209 * t214;
t143 = t177 * t260 + t207 * t247;
t176 = t208 * t245 - t210 * t212;
t115 = -t143 * t206 + t176 * t211;
t116 = t143 * t211 + t176 * t206;
t249 = t209 * t210;
t141 = t175 * t260 + t207 * t249;
t174 = t208 * t246 + t212 * t214;
t113 = -t141 * t206 + t174 * t211;
t114 = t141 * t211 + t174 * t206;
t70 = Icges(7,5) * t114 + Icges(7,6) * t113 + Icges(7,3) * t140;
t72 = Icges(7,4) * t114 + Icges(7,2) * t113 + Icges(7,6) * t140;
t74 = Icges(7,1) * t114 + Icges(7,4) * t113 + Icges(7,5) * t140;
t24 = t115 * t72 + t116 * t74 + t142 * t70;
t71 = Icges(7,5) * t116 + Icges(7,6) * t115 + Icges(7,3) * t142;
t73 = Icges(7,4) * t116 + Icges(7,2) * t115 + Icges(7,6) * t142;
t75 = Icges(7,1) * t116 + Icges(7,4) * t115 + Icges(7,5) * t142;
t25 = t115 * t73 + t116 * t75 + t142 * t71;
t31 = t115 * t94 + t116 * t97 + t142 * t91;
t2 = t140 * t24 + t142 * t25 + t172 * t31;
t265 = t2 / 0.2e1;
t264 = t140 / 0.2e1;
t263 = t142 / 0.2e1;
t262 = t172 / 0.2e1;
t259 = rSges(4,1) * t213;
t122 = Icges(5,5) * t175 - Icges(5,6) * t174 + Icges(5,3) * t249;
t125 = Icges(5,4) * t175 - Icges(5,2) * t174 + Icges(5,6) * t249;
t128 = Icges(5,1) * t175 - Icges(5,4) * t174 + Icges(5,5) * t249;
t64 = -t122 * t213 + (-t125 * t208 + t128 * t212) * t209;
t258 = t213 * t64;
t123 = Icges(5,5) * t177 - Icges(5,6) * t176 + Icges(5,3) * t247;
t126 = Icges(5,4) * t177 - Icges(5,2) * t176 + Icges(5,6) * t247;
t129 = Icges(5,1) * t177 - Icges(5,4) * t176 + Icges(5,5) * t247;
t65 = -t123 * t213 + (-t126 * t208 + t129 * t212) * t209;
t257 = t213 * t65;
t256 = t214 * rSges(4,2);
t254 = Icges(3,4) * t213;
t252 = Icges(4,5) * t213;
t244 = t243 * qJ(3) * t209;
t92 = Icges(6,5) * t141 - Icges(6,6) * t140 + Icges(6,3) * t174;
t95 = Icges(6,4) * t141 - Icges(6,2) * t140 + Icges(6,6) * t174;
t98 = Icges(6,1) * t141 - Icges(6,4) * t140 + Icges(6,5) * t174;
t38 = -t140 * t95 + t141 * t98 + t174 * t92;
t93 = Icges(6,5) * t143 - Icges(6,6) * t142 + Icges(6,3) * t176;
t96 = Icges(6,4) * t143 - Icges(6,2) * t142 + Icges(6,6) * t176;
t99 = Icges(6,1) * t143 - Icges(6,4) * t142 + Icges(6,5) * t176;
t39 = -t140 * t96 + t141 * t99 + t174 * t93;
t53 = t121 * t174 - t124 * t140 + t127 * t141;
t13 = t174 * t38 + t176 * t39 + t250 * t53;
t22 = t113 * t72 + t114 * t74 + t140 * t70;
t23 = t113 * t73 + t114 * t75 + t140 * t71;
t30 = t113 * t94 + t114 * t97 + t140 * t91;
t3 = t174 * t22 + t176 * t23 + t250 * t30;
t242 = t3 / 0.2e1 + t13 / 0.2e1;
t40 = -t142 * t95 + t143 * t98 + t176 * t92;
t41 = -t142 * t96 + t143 * t99 + t176 * t93;
t54 = t121 * t176 - t124 * t142 + t127 * t143;
t14 = t174 * t40 + t176 * t41 + t250 * t54;
t4 = t174 * t24 + t176 * t25 + t250 * t31;
t241 = t4 / 0.2e1 + t14 / 0.2e1;
t15 = -t53 * t213 + (t210 * t38 + t214 * t39) * t209;
t5 = -t30 * t213 + (t210 * t22 + t214 * t23) * t209;
t240 = t5 / 0.2e1 + t15 / 0.2e1;
t16 = -t54 * t213 + (t210 * t40 + t214 * t41) * t209;
t6 = -t31 * t213 + (t210 * t24 + t214 * t25) * t209;
t239 = t6 / 0.2e1 + t16 / 0.2e1;
t10 = t23 * t210 - t214 * t22;
t19 = t39 * t210 - t214 * t38;
t238 = t10 / 0.2e1 + t19 / 0.2e1;
t11 = t25 * t210 - t214 * t24;
t20 = t41 * t210 - t214 * t40;
t237 = t11 / 0.2e1 + t20 / 0.2e1;
t26 = t138 * t72 + t139 * t74 + t172 * t70;
t27 = t138 * t73 + t139 * t75 + t172 * t71;
t12 = t27 * t210 - t26 * t214;
t45 = -t172 * t95 + t173 * t98 + t250 * t92;
t46 = -t172 * t96 + t173 * t99 + t250 * t93;
t236 = t12 / 0.2e1 + t46 * t261 + t45 * t273;
t235 = t30 / 0.2e1 + t26 / 0.2e1;
t234 = t31 / 0.2e1 + t27 / 0.2e1;
t233 = qJ(3) * t249;
t77 = t116 * rSges(7,1) + t115 * rSges(7,2) + t142 * rSges(7,3);
t102 = t143 * rSges(6,1) - t142 * rSges(6,2) + t176 * rSges(6,3);
t132 = t177 * rSges(5,1) - t176 * rSges(5,2) + rSges(5,3) * t247;
t232 = rSges(4,1) * t245 + t210 * rSges(4,2) + rSges(4,3) * t247;
t230 = Icges(4,6) * t274 + Icges(3,6) * t213 / 0.2e1 + t281 * t209 / 0.2e1;
t229 = -rSges(5,1) * t175 + rSges(5,2) * t174;
t226 = -Icges(3,2) * t209 + t254;
t223 = Icges(4,3) * t209 + t252;
t218 = t53 / 0.2e1 + t45 / 0.2e1 + t235;
t217 = t46 / 0.2e1 + t54 / 0.2e1 + t234;
t101 = rSges(6,1) * t141 - rSges(6,2) * t140 + rSges(6,3) * t174;
t76 = rSges(7,1) * t114 + rSges(7,2) * t113 + rSges(7,3) * t140;
t149 = -Icges(5,3) * t213 + (Icges(5,5) * t212 - Icges(5,6) * t208) * t209;
t154 = -Icges(5,6) * t213 + (Icges(5,4) * t212 - Icges(5,2) * t208) * t209;
t159 = -Icges(5,5) * t213 + (Icges(5,1) * t212 - Icges(5,4) * t208) * t209;
t83 = t149 * t249 - t154 * t174 + t159 * t175;
t216 = t83 / 0.2e1 + t64 / 0.2e1 + t218;
t84 = t149 * t247 - t176 * t154 + t177 * t159;
t215 = t65 / 0.2e1 + t84 / 0.2e1 + t217;
t195 = qJ(3) * t245;
t194 = qJ(3) * t247;
t192 = qJ(3) * t246;
t188 = rSges(2,1) * t214 - t210 * rSges(2,2);
t187 = -t210 * rSges(2,1) - rSges(2,2) * t214;
t186 = rSges(3,1) * t209 + rSges(3,2) * t213;
t185 = rSges(4,1) * t209 - rSges(4,3) * t213;
t167 = t210 * rSges(3,3) + (rSges(3,1) * t213 - rSges(3,2) * t209) * t214;
t165 = rSges(3,1) * t246 - rSges(3,2) * t249 - t214 * rSges(3,3);
t164 = -rSges(5,3) * t213 + (rSges(5,1) * t212 - rSges(5,2) * t208) * t209;
t148 = -t185 * t214 + t195;
t147 = -t185 * t210 + t192;
t146 = t159 * t248;
t145 = t194 + t232;
t144 = t256 + (-t259 + (-rSges(4,3) - qJ(3)) * t209) * t210;
t134 = -t164 * t214 + t195;
t133 = -t164 * t210 + t192;
t131 = rSges(5,3) * t249 - t229;
t130 = rSges(6,1) * t173 - rSges(6,2) * t172 + rSges(6,3) * t250;
t120 = t210 * t165 + t167 * t214;
t118 = t194 + t132;
t117 = (-rSges(5,3) - qJ(3)) * t249 + t229;
t110 = -t130 * t214 + t195;
t109 = -t130 * t210 + t192;
t106 = t214 * t232 + (-t256 + (rSges(4,3) * t209 + t259) * t210) * t210 + t244;
t105 = -t213 * t132 - t164 * t247;
t104 = t131 * t213 + t164 * t249;
t103 = -t213 * t149 - t154 * t250 + t146;
t100 = rSges(7,1) * t139 + rSges(7,2) * t138 + rSges(7,3) * t172;
t90 = t194 + t102;
t89 = -t101 - t233;
t88 = -t100 * t214 + t195;
t87 = -t100 * t210 + t192;
t85 = (t131 * t214 - t132 * t210) * t209;
t82 = t210 * t131 + t132 * t214 + t244;
t79 = -t213 * t102 - t130 * t247;
t78 = t101 * t213 + t130 * t249;
t69 = t102 * t250 - t130 * t176;
t68 = -t101 * t250 + t130 * t174;
t67 = t194 + t77;
t66 = -t76 - t233;
t63 = t123 * t247 - t176 * t126 + t177 * t129;
t62 = t122 * t247 - t176 * t125 + t177 * t128;
t61 = t123 * t249 - t126 * t174 + t129 * t175;
t60 = t122 * t249 - t125 * t174 + t128 * t175;
t58 = (t101 * t214 - t102 * t210) * t209;
t57 = t59 * t250;
t56 = t210 * t101 + t102 * t214 + t244;
t55 = t101 * t176 - t102 * t174;
t52 = -t100 * t247 - t213 * t77;
t51 = t100 * t249 + t213 * t76;
t50 = -t100 * t176 + t250 * t77;
t49 = t100 * t174 - t250 * t76;
t48 = -t100 * t142 + t172 * t77;
t47 = t100 * t140 - t172 * t76;
t44 = (-t210 * t77 + t214 * t76) * t209;
t43 = t210 * t76 + t214 * t77 + t244;
t42 = -t174 * t77 + t176 * t76;
t36 = t37 * t250;
t35 = t37 * t172;
t34 = -t140 * t77 + t142 * t76;
t33 = t63 * t210 - t214 * t62;
t32 = t61 * t210 - t214 * t60;
t29 = -t84 * t213 + (t210 * t62 + t214 * t63) * t209;
t28 = -t83 * t213 + (t210 * t60 + t214 * t61) * t209;
t18 = -t59 * t213 + (t45 * t210 + t46 * t214) * t209;
t17 = t45 * t174 + t46 * t176 + t57;
t9 = -t37 * t213 + (t26 * t210 + t27 * t214) * t209;
t8 = t26 * t174 + t27 * t176 + t36;
t7 = t26 * t140 + t27 * t142 + t35;
t1 = t140 * t22 + t142 * t23 + t172 * t30;
t21 = [Icges(2,3) + t146 + m(7) * (t66 ^ 2 + t67 ^ 2) + m(6) * (t89 ^ 2 + t90 ^ 2) + m(5) * (t117 ^ 2 + t118 ^ 2) + m(4) * (t144 ^ 2 + t145 ^ 2) + m(3) * (t165 ^ 2 + t167 ^ 2) + m(2) * (t187 ^ 2 + t188 ^ 2) + (-t149 + (Icges(3,2) + Icges(4,3)) * t213 + t282) * t213 + (-t208 * t154 + t283 * t209 - t252 + t254) * t209 - t272; (t230 * t214 + (Icges(3,6) * t276 + Icges(4,6) * t273 + t223 * t261 + t226 * t277) * t213 + (t276 * t281 + t278 * t277) * t209 - t216) * t214 + (t230 * t210 + (Icges(3,6) * t261 + Icges(4,6) * t277 + t223 * t273 + t226 * t276) * t213 + (t261 * t281 + t278 * t276) * t209 + t215) * t210 + m(7) * (t66 * t88 + t67 * t87) + m(6) * (t109 * t90 + t110 * t89) + m(5) * (t117 * t134 + t118 * t133) + m(4) * (t144 * t148 + t145 * t147) + m(3) * (t165 * t214 - t167 * t210) * t186; m(7) * (t43 ^ 2 + t87 ^ 2 + t88 ^ 2) + m(6) * (t109 ^ 2 + t110 ^ 2 + t56 ^ 2) + m(5) * (t133 ^ 2 + t134 ^ 2 + t82 ^ 2) + m(4) * (t106 ^ 2 + t147 ^ 2 + t148 ^ 2) + m(3) * (t186 ^ 2 * t243 + t120 ^ 2) + (t205 * t270 - t10 - t19 - t32) * t214 + (t11 + t20 + t33 + t271 * t204 + (t210 * t270 + t214 * t271) * t214) * t210; 0.2e1 * ((t210 * t67 + t214 * t66) * t266 + (t210 * t90 + t214 * t89) * t267 + (t117 * t214 + t118 * t210) * t268 + (t144 * t214 + t145 * t210) * t269) * t209; m(7) * (-t213 * t43 + (t210 * t87 + t214 * t88) * t209) + m(6) * (-t213 * t56 + (t109 * t210 + t110 * t214) * t209) + m(5) * (-t213 * t82 + (t133 * t210 + t134 * t214) * t209) + m(4) * (-t213 * t106 + (t147 * t210 + t148 * t214) * t209); 0.2e1 * (t269 + t268 + t267 + t266) * (t209 ^ 2 * t243 + t213 ^ 2); (-t103 + t272) * t213 + m(7) * (t51 * t66 + t52 * t67) + m(6) * (t78 * t89 + t79 * t90) + m(5) * (t104 * t117 + t105 * t118) + (t210 * t216 + t214 * t215) * t209; -t236 * t213 + (-t28 / 0.2e1 + t258 / 0.2e1 - t240) * t214 + (-t257 / 0.2e1 + t29 / 0.2e1 + t239) * t210 + m(7) * (t44 * t43 + t51 * t88 + t52 * t87) + m(6) * (t109 * t79 + t110 * t78 + t58 * t56) + m(5) * (t104 * t134 + t105 * t133 + t82 * t85) + ((t33 / 0.2e1 + t237) * t214 + (t32 / 0.2e1 + t238) * t210) * t209; m(5) * (-t85 * t213 + (t104 * t214 + t105 * t210) * t209) + m(6) * (-t58 * t213 + (t210 * t79 + t214 * t78) * t209) + m(7) * (-t44 * t213 + (t210 * t52 + t214 * t51) * t209); (t103 * t213 - t18 - t9) * t213 + m(7) * (t44 ^ 2 + t51 ^ 2 + t52 ^ 2) + m(6) * (t58 ^ 2 + t78 ^ 2 + t79 ^ 2) + m(5) * (t104 ^ 2 + t105 ^ 2 + t85 ^ 2) + ((t16 + t29 + t6 - t257) * t214 + (t15 + t28 + t5 - t258) * t210) * t209; t36 + t57 + m(7) * (t49 * t66 + t50 * t67) + m(6) * (t68 * t89 + t69 * t90) + t217 * t176 + t218 * t174; -t242 * t214 + t241 * t210 + t236 * t250 + t237 * t176 + t238 * t174 + m(7) * (t42 * t43 + t49 * t88 + t50 * t87) + m(6) * (t109 * t69 + t110 * t68 + t55 * t56); m(6) * (-t55 * t213 + (t210 * t69 + t214 * t68) * t209) + m(7) * (-t42 * t213 + (t210 * t50 + t214 * t49) * t209); (-t8 / 0.2e1 - t17 / 0.2e1) * t213 + t239 * t176 + t240 * t174 + m(7) * (t42 * t44 + t49 * t51 + t50 * t52) + m(6) * (t55 * t58 + t68 * t78 + t69 * t79) + (t241 * t214 + t242 * t210 + (t9 / 0.2e1 + t18 / 0.2e1) * t208) * t209; (t17 + t8) * t250 + (t4 + t14) * t176 + (t3 + t13) * t174 + m(7) * (t42 ^ 2 + t49 ^ 2 + t50 ^ 2) + m(6) * (t55 ^ 2 + t68 ^ 2 + t69 ^ 2); t35 + m(7) * (t47 * t66 + t48 * t67) + t234 * t142 + t235 * t140; t10 * t264 + t12 * t262 + t1 * t273 + m(7) * (t34 * t43 + t47 * t88 + t48 * t87) + t2 * t261 + t11 * t263; m(7) * (-t34 * t213 + (t210 * t48 + t214 * t47) * t209); m(7) * (t34 * t44 + t47 * t51 + t48 * t52) + t6 * t263 + t9 * t262 + t5 * t264 + t7 * t274 + (t1 * t261 + t214 * t265) * t209; t7 * t250 / 0.2e1 + m(7) * (t34 * t42 + t47 * t49 + t48 * t50) + t176 * t265 + t3 * t264 + t4 * t263 + t8 * t262 + t174 * t1 / 0.2e1; t142 * t2 + t140 * t1 + t172 * t7 + m(7) * (t34 ^ 2 + t47 ^ 2 + t48 ^ 2);];
%% Postprocessing: Reshape Output
% From vec2symmat_6_matlab.m
res = [t21(1) t21(2) t21(4) t21(7) t21(11) t21(16); t21(2) t21(3) t21(5) t21(8) t21(12) t21(17); t21(4) t21(5) t21(6) t21(9) t21(13) t21(18); t21(7) t21(8) t21(9) t21(10) t21(14) t21(19); t21(11) t21(12) t21(13) t21(14) t21(15) t21(20); t21(16) t21(17) t21(18) t21(19) t21(20) t21(21);];
Mq  = res;
