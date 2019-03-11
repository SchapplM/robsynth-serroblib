% Calculate joint inertia matrix for
% S6RPRPRR13
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% pkin [12x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,alpha2,alpha3,d1,d3,d5,d6,theta2]';
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
% Datum: 2019-03-09 04:26
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Mq = S6RPRPRR13_inertiaJ_slag_vp1(qJ, ...
  pkin, m, rSges, Icges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(12,1),zeros(7,1),zeros(7,3),zeros(7,6)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPRPRR13_inertiaJ_slag_vp1: qJ has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [12 1]), ...
  'S6RPRPRR13_inertiaJ_slag_vp1: pkin has to be [12x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6RPRPRR13_inertiaJ_slag_vp1: m has to be [7x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [7,3]), ...
  'S6RPRPRR13_inertiaJ_slag_vp1: rSges has to be [7x3] (double)');
assert(isreal(Icges) && all(size(Icges) == [7 6]), ...
  'S6RPRPRR13_inertiaJ_slag_vp1: Icges has to be [7x6] (double)'); 

%% Symbolic Calculation
% From inertia_joint_joint_floatb_twist_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 04:21:46
% EndTime: 2019-03-09 04:21:56
% DurationCPUTime: 5.17s
% Computational Cost: add. (25430->515), mult. (70403->717), div. (0->0), fcn. (92444->14), ass. (0->240)
t222 = sin(qJ(1));
t219 = cos(pkin(6));
t274 = cos(pkin(12));
t251 = t222 * t274;
t217 = sin(pkin(12));
t224 = cos(qJ(1));
t268 = t224 * t217;
t235 = t219 * t251 + t268;
t218 = sin(pkin(6));
t275 = cos(pkin(7));
t253 = t218 * t275;
t273 = sin(pkin(7));
t196 = t222 * t253 + t235 * t273;
t250 = t224 * t274;
t269 = t222 * t217;
t236 = -t219 * t250 + t269;
t195 = -t224 * t253 + t236 * t273;
t244 = t275 * t274;
t284 = cos(qJ(3));
t246 = t284 * t273;
t272 = t217 * t218;
t282 = sin(qJ(3));
t192 = -t218 * t244 * t284 - t219 * t246 + t272 * t282;
t245 = t273 * t282;
t193 = t219 * t245 + (t217 * t284 + t244 * t282) * t218;
t252 = t218 * t274;
t202 = t219 * t275 - t252 * t273;
t160 = Icges(5,5) * t202 - Icges(5,6) * t193 + Icges(5,3) * t192;
t161 = Icges(4,5) * t193 - Icges(4,6) * t192 + Icges(4,3) * t202;
t162 = Icges(5,4) * t202 - Icges(5,2) * t193 + Icges(5,6) * t192;
t163 = Icges(4,4) * t193 - Icges(4,2) * t192 + Icges(4,6) * t202;
t164 = Icges(5,1) * t202 - Icges(5,4) * t193 + Icges(5,5) * t192;
t165 = Icges(4,1) * t193 - Icges(4,4) * t192 + Icges(4,5) * t202;
t301 = (t161 + t164) * t202 + (-t162 + t165) * t193 + (t160 - t163) * t192;
t291 = m(7) / 0.2e1;
t292 = m(6) / 0.2e1;
t293 = m(5) / 0.2e1;
t249 = t293 + t292 + t291;
t300 = 0.2e1 * t249;
t203 = t219 * t268 + t251;
t232 = t236 * t275;
t240 = t218 * t246;
t184 = t203 * t282 + t224 * t240 + t232 * t284;
t239 = t218 * t245;
t185 = t203 * t284 - t224 * t239 - t232 * t282;
t117 = Icges(5,5) * t195 - Icges(5,6) * t185 + Icges(5,3) * t184;
t123 = Icges(4,4) * t185 - Icges(4,2) * t184 + Icges(4,6) * t195;
t299 = t117 - t123;
t204 = -t219 * t269 + t250;
t230 = t235 * t275;
t186 = t204 * t282 - t222 * t240 + t230 * t284;
t187 = t204 * t284 + t222 * t239 - t230 * t282;
t118 = Icges(5,5) * t196 - Icges(5,6) * t187 + Icges(5,3) * t186;
t124 = Icges(4,4) * t187 - Icges(4,2) * t186 + Icges(4,6) * t196;
t298 = t118 - t124;
t121 = Icges(5,4) * t195 - Icges(5,2) * t185 + Icges(5,6) * t184;
t127 = Icges(4,1) * t185 - Icges(4,4) * t184 + Icges(4,5) * t195;
t297 = -t121 + t127;
t122 = Icges(5,4) * t196 - Icges(5,2) * t187 + Icges(5,6) * t186;
t128 = Icges(4,1) * t187 - Icges(4,4) * t186 + Icges(4,5) * t196;
t296 = -t122 + t128;
t295 = m(3) / 0.2e1;
t294 = m(4) / 0.2e1;
t221 = sin(qJ(5));
t283 = cos(qJ(5));
t154 = -t184 * t283 + t195 * t221;
t156 = -t186 * t283 + t196 * t221;
t155 = t184 * t221 + t195 * t283;
t220 = sin(qJ(6));
t223 = cos(qJ(6));
t109 = -t155 * t220 + t185 * t223;
t110 = t155 * t223 + t185 * t220;
t65 = Icges(7,5) * t110 + Icges(7,6) * t109 + Icges(7,3) * t154;
t67 = Icges(7,4) * t110 + Icges(7,2) * t109 + Icges(7,6) * t154;
t69 = Icges(7,1) * t110 + Icges(7,4) * t109 + Icges(7,5) * t154;
t17 = t109 * t67 + t110 * t69 + t154 * t65;
t157 = t186 * t221 + t196 * t283;
t111 = -t157 * t220 + t187 * t223;
t112 = t157 * t223 + t187 * t220;
t66 = Icges(7,5) * t112 + Icges(7,6) * t111 + Icges(7,3) * t156;
t68 = Icges(7,4) * t112 + Icges(7,2) * t111 + Icges(7,6) * t156;
t70 = Icges(7,1) * t112 + Icges(7,4) * t111 + Icges(7,5) * t156;
t18 = t109 * t68 + t110 * t70 + t154 * t66;
t182 = -t192 * t283 + t202 * t221;
t183 = t192 * t221 + t202 * t283;
t148 = -t183 * t220 + t193 * t223;
t149 = t183 * t223 + t193 * t220;
t87 = Icges(7,5) * t149 + Icges(7,6) * t148 + Icges(7,3) * t182;
t88 = Icges(7,4) * t149 + Icges(7,2) * t148 + Icges(7,6) * t182;
t89 = Icges(7,1) * t149 + Icges(7,4) * t148 + Icges(7,5) * t182;
t26 = t109 * t88 + t110 * t89 + t154 * t87;
t1 = t154 * t17 + t156 * t18 + t182 * t26;
t290 = t1 / 0.2e1;
t19 = t111 * t67 + t112 * t69 + t156 * t65;
t20 = t111 * t68 + t112 * t70 + t156 * t66;
t27 = t111 * t88 + t112 * t89 + t156 * t87;
t2 = t154 * t19 + t156 * t20 + t182 * t27;
t289 = t2 / 0.2e1;
t23 = t148 * t67 + t149 * t69 + t182 * t65;
t24 = t148 * t68 + t149 * t70 + t182 * t66;
t34 = t148 * t88 + t149 * t89 + t182 * t87;
t31 = t34 * t182;
t7 = t23 * t154 + t24 * t156 + t31;
t288 = t7 / 0.2e1;
t287 = t154 / 0.2e1;
t286 = t156 / 0.2e1;
t285 = t182 / 0.2e1;
t281 = pkin(3) * t185;
t280 = pkin(5) * t155;
t279 = t301 * t202;
t241 = -rSges(7,1) * t110 - rSges(7,2) * t109;
t71 = rSges(7,3) * t154 - t241;
t278 = pkin(11) * t154 + t280 + t71;
t72 = rSges(7,1) * t112 + rSges(7,2) * t111 + rSges(7,3) * t156;
t277 = pkin(5) * t157 + pkin(11) * t156 + t72;
t90 = rSges(7,1) * t149 + rSges(7,2) * t148 + rSges(7,3) * t182;
t276 = pkin(5) * t183 + pkin(11) * t182 + t90;
t271 = t218 * t222;
t270 = t218 * t224;
t120 = Icges(4,5) * t187 - Icges(4,6) * t186 + Icges(4,3) * t196;
t126 = Icges(5,1) * t196 - Icges(5,4) * t187 + Icges(5,5) * t186;
t267 = t120 + t126;
t119 = Icges(4,5) * t185 - Icges(4,6) * t184 + Icges(4,3) * t195;
t125 = Icges(5,1) * t195 - Icges(5,4) * t185 + Icges(5,5) * t184;
t266 = t125 + t119;
t172 = t184 * qJ(4);
t142 = t172 + t281;
t133 = t196 * t142;
t158 = pkin(4) * t195 + pkin(10) * t185;
t265 = t158 * t196 + t133;
t143 = pkin(3) * t187 + qJ(4) * t186;
t135 = t202 * t143;
t159 = pkin(4) * t196 + pkin(10) * t187;
t264 = t159 * t202 + t135;
t263 = -t142 - t158;
t262 = -t143 - t159;
t169 = pkin(3) * t193 + qJ(4) * t192;
t147 = t195 * t169;
t188 = pkin(4) * t202 + pkin(10) * t193;
t261 = t188 * t195 + t147;
t260 = -t169 - t188;
t259 = pkin(1) * t224 + qJ(2) * t271;
t258 = t23 / 0.2e1 + t26 / 0.2e1;
t257 = t24 / 0.2e1 + t27 / 0.2e1;
t113 = Icges(6,5) * t183 - Icges(6,6) * t182 + Icges(6,3) * t193;
t114 = Icges(6,4) * t183 - Icges(6,2) * t182 + Icges(6,6) * t193;
t115 = Icges(6,1) * t183 - Icges(6,4) * t182 + Icges(6,5) * t193;
t54 = t113 * t193 - t114 * t182 + t115 * t183;
t98 = rSges(6,1) * t157 - rSges(6,2) * t156 + rSges(6,3) * t187;
t130 = rSges(4,1) * t187 - rSges(4,2) * t186 + rSges(4,3) * t196;
t132 = rSges(5,1) * t196 - rSges(5,2) * t187 + rSges(5,3) * t186;
t254 = -t222 * pkin(1) + qJ(2) * t270;
t243 = -rSges(5,1) * t195 - rSges(5,3) * t184;
t242 = -rSges(6,1) * t155 + rSges(6,2) * t154;
t91 = Icges(6,5) * t155 - Icges(6,6) * t154 + Icges(6,3) * t185;
t93 = Icges(6,4) * t155 - Icges(6,2) * t154 + Icges(6,6) * t185;
t95 = Icges(6,1) * t155 - Icges(6,4) * t154 + Icges(6,5) * t185;
t40 = -t182 * t93 + t183 * t95 + t193 * t91;
t48 = t113 * t185 - t114 * t154 + t115 * t155;
t238 = t48 / 0.2e1 + t40 / 0.2e1 + t258;
t92 = Icges(6,5) * t157 - Icges(6,6) * t156 + Icges(6,3) * t187;
t94 = Icges(6,4) * t157 - Icges(6,2) * t156 + Icges(6,6) * t187;
t96 = Icges(6,1) * t157 - Icges(6,4) * t156 + Icges(6,5) * t187;
t41 = -t182 * t94 + t183 * t96 + t193 * t92;
t49 = t113 * t187 - t114 * t156 + t115 * t157;
t237 = t41 / 0.2e1 + t49 / 0.2e1 + t257;
t129 = rSges(4,1) * t185 - rSges(4,2) * t184 + rSges(4,3) * t195;
t234 = -pkin(2) * t203 - pkin(9) * t195 + t254;
t233 = -t172 + t234;
t228 = -t158 + t233;
t227 = t204 * pkin(2) + pkin(9) * t196 + t259;
t226 = t143 + t227;
t225 = t159 + t226;
t211 = rSges(2,1) * t224 - rSges(2,2) * t222;
t210 = -rSges(2,1) * t222 - rSges(2,2) * t224;
t181 = t204 * rSges(3,1) - rSges(3,2) * t235 + rSges(3,3) * t271 + t259;
t180 = -t203 * rSges(3,1) + rSges(3,2) * t236 + rSges(3,3) * t270 + t254;
t167 = rSges(5,1) * t202 - rSges(5,2) * t193 + rSges(5,3) * t192;
t166 = rSges(4,1) * t193 - rSges(4,2) * t192 + rSges(4,3) * t202;
t131 = -rSges(5,2) * t185 - t243;
t116 = rSges(6,1) * t183 - rSges(6,2) * t182 + rSges(6,3) * t193;
t102 = t227 + t130;
t101 = -t129 + t234;
t97 = rSges(6,3) * t185 - t242;
t86 = t130 * t202 - t166 * t196;
t85 = -t129 * t202 + t166 * t195;
t83 = t226 + t132;
t82 = (rSges(5,2) - pkin(3)) * t185 + t233 + t243;
t79 = t129 * t196 - t130 * t195;
t76 = t160 * t186 - t162 * t187 + t164 * t196;
t75 = t160 * t184 - t162 * t185 + t164 * t195;
t74 = t161 * t196 - t163 * t186 + t165 * t187;
t73 = t161 * t195 - t163 * t184 + t165 * t185;
t64 = -t116 * t187 + t193 * t98;
t63 = t116 * t185 - t193 * t97;
t62 = t225 + t98;
t61 = (-rSges(6,3) - pkin(3)) * t185 + t228 + t242;
t60 = t132 * t202 + t135 + (-t167 - t169) * t196;
t59 = t167 * t195 + t147 + (-t131 - t142) * t202;
t58 = t118 * t192 - t122 * t193 + t126 * t202;
t57 = t117 * t192 - t121 * t193 + t125 * t202;
t56 = t120 * t202 - t124 * t192 + t128 * t193;
t55 = t119 * t202 - t123 * t192 + t127 * t193;
t53 = -t185 * t98 + t187 * t97;
t52 = t54 * t202;
t51 = t54 * t193;
t50 = t131 * t196 + t133 + (-t132 - t143) * t195;
t47 = -t156 * t90 + t182 * t72;
t46 = t154 * t90 - t182 * t71;
t45 = t225 + t277;
t44 = -t281 - t280 + (-rSges(7,3) - pkin(11)) * t154 + t228 + t241;
t43 = t202 * t98 + (-t116 + t260) * t196 + t264;
t42 = t116 * t195 + (-t97 + t263) * t202 + t261;
t39 = -t156 * t94 + t157 * t96 + t187 * t92;
t38 = -t156 * t93 + t157 * t95 + t187 * t91;
t37 = -t154 * t94 + t155 * t96 + t185 * t92;
t36 = -t154 * t93 + t155 * t95 + t185 * t91;
t35 = -t154 * t72 + t156 * t71;
t33 = t34 * t202;
t32 = t34 * t193;
t30 = t196 * t97 + (-t98 + t262) * t195 + t265;
t29 = -t187 * t276 + t193 * t277;
t28 = t185 * t276 - t193 * t278;
t25 = -t185 * t277 + t187 * t278;
t22 = t277 * t202 + (t260 - t276) * t196 + t264;
t21 = t276 * t195 + (t263 - t278) * t202 + t261;
t16 = t278 * t196 + (t262 - t277) * t195 + t265;
t15 = t40 * t195 + t41 * t196 + t52;
t14 = t40 * t185 + t41 * t187 + t51;
t13 = t195 * t38 + t196 * t39 + t202 * t49;
t12 = t195 * t36 + t196 * t37 + t202 * t48;
t11 = t185 * t38 + t187 * t39 + t193 * t49;
t10 = t185 * t36 + t187 * t37 + t193 * t48;
t9 = t23 * t195 + t24 * t196 + t33;
t8 = t23 * t185 + t24 * t187 + t32;
t6 = t19 * t195 + t196 * t20 + t202 * t27;
t5 = t17 * t195 + t18 * t196 + t202 * t26;
t4 = t185 * t19 + t187 * t20 + t193 * t27;
t3 = t17 * t185 + t18 * t187 + t193 * t26;
t77 = [m(7) * (t44 ^ 2 + t45 ^ 2) + m(6) * (t61 ^ 2 + t62 ^ 2) + m(5) * (t82 ^ 2 + t83 ^ 2) + m(4) * (t101 ^ 2 + t102 ^ 2) + m(3) * (t180 ^ 2 + t181 ^ 2) + m(2) * (t210 ^ 2 + t211 ^ 2) + t219 * (Icges(3,3) * t219 + (Icges(3,5) * t217 + Icges(3,6) * t274) * t218) + (Icges(3,5) * t219 + (Icges(3,1) * t217 + Icges(3,4) * t274) * t218) * t272 + Icges(2,3) + (Icges(3,6) * t219 + (Icges(3,4) * t217 + Icges(3,2) * t274) * t218) * t252 + t54 + t34 + t301; 0.2e1 * ((t222 * t44 - t224 * t45) * t291 + (t222 * t61 - t224 * t62) * t292 + (t222 * t82 - t224 * t83) * t293 + (t101 * t222 - t102 * t224) * t294 + (t180 * t222 - t181 * t224) * t295) * t218; 0.2e1 * (t295 + t294 + t249) * (t219 ^ 2 + (t222 ^ 2 + t224 ^ 2) * t218 ^ 2); t33 + t52 + m(7) * (t21 * t44 + t22 * t45) + m(6) * (t42 * t61 + t43 * t62) + m(5) * (t59 * t82 + t60 * t83) + m(4) * (t101 * t85 + t102 * t86) + (t56 / 0.2e1 + t58 / 0.2e1 + t74 / 0.2e1 + t76 / 0.2e1 + t237) * t196 + (t75 / 0.2e1 + t73 / 0.2e1 + t55 / 0.2e1 + t57 / 0.2e1 + t238) * t195 + t279; m(4) * (t79 * t219 + (t222 * t85 - t224 * t86) * t218) + m(5) * (t50 * t219 + (t222 * t59 - t224 * t60) * t218) + m(6) * (t30 * t219 + (t222 * t42 - t224 * t43) * t218) + m(7) * (t16 * t219 + (t21 * t222 - t22 * t224) * t218); (t9 + t15 + t279) * t202 + m(7) * (t16 ^ 2 + t21 ^ 2 + t22 ^ 2) + m(6) * (t30 ^ 2 + t42 ^ 2 + t43 ^ 2) + m(5) * (t50 ^ 2 + t59 ^ 2 + t60 ^ 2) + m(4) * (t79 ^ 2 + t85 ^ 2 + t86 ^ 2) + (t13 + t6 + (t186 * t298 + t187 * t296 + t267 * t196) * t196 + (t56 + t58 + t74 + t76) * t202) * t196 + (t5 + t12 + (t184 * t299 + t185 * t297 + t266 * t195) * t195 + (t55 + t57 + t75 + t73) * t202 + (t184 * t298 + t185 * t296 + t186 * t299 + t187 * t297 + t267 * t195 + t266 * t196) * t196) * t195; m(7) * (t184 * t45 + t186 * t44) + m(6) * (t184 * t62 + t186 * t61) + m(5) * (t184 * t83 + t186 * t82); (t192 * t219 + (-t184 * t224 + t186 * t222) * t218) * t300; m(7) * (t16 * t192 + t184 * t22 + t186 * t21) + m(6) * (t184 * t43 + t186 * t42 + t192 * t30) + m(5) * (t184 * t60 + t186 * t59 + t192 * t50); (t184 ^ 2 + t186 ^ 2 + t192 ^ 2) * t300; t32 + t51 + m(7) * (t28 * t44 + t29 * t45) + m(6) * (t61 * t63 + t62 * t64) + t237 * t187 + t238 * t185; m(6) * (t53 * t219 + (t222 * t63 - t224 * t64) * t218) + m(7) * (t25 * t219 + (t222 * t28 - t224 * t29) * t218); (t8 / 0.2e1 + t14 / 0.2e1) * t202 + (t4 / 0.2e1 + t11 / 0.2e1) * t196 + (t3 / 0.2e1 + t10 / 0.2e1) * t195 + (t9 / 0.2e1 + t15 / 0.2e1) * t193 + (t6 / 0.2e1 + t13 / 0.2e1) * t187 + (t5 / 0.2e1 + t12 / 0.2e1) * t185 + m(7) * (t16 * t25 + t21 * t28 + t22 * t29) + m(6) * (t30 * t53 + t42 * t63 + t43 * t64); m(6) * (t184 * t64 + t186 * t63 + t192 * t53) + m(7) * (t184 * t29 + t186 * t28 + t192 * t25); (t8 + t14) * t193 + (t4 + t11) * t187 + (t3 + t10) * t185 + m(7) * (t25 ^ 2 + t28 ^ 2 + t29 ^ 2) + m(6) * (t53 ^ 2 + t63 ^ 2 + t64 ^ 2); m(7) * (t44 * t46 + t45 * t47) + t31 + t257 * t156 + t258 * t154; m(7) * (t35 * t219 + (t222 * t46 - t224 * t47) * t218); m(7) * (t16 * t35 + t21 * t46 + t22 * t47) + t195 * t290 + t202 * t288 + t6 * t286 + t9 * t285 + t5 * t287 + t196 * t289; m(7) * (t184 * t47 + t186 * t46 + t192 * t35); m(7) * (t25 * t35 + t28 * t46 + t29 * t47) + t8 * t285 + t4 * t286 + t185 * t290 + t193 * t288 + t3 * t287 + t187 * t289; t156 * t2 + t154 * t1 + t182 * t7 + m(7) * (t35 ^ 2 + t46 ^ 2 + t47 ^ 2);];
%% Postprocessing: Reshape Output
% From vec2symmat_6_matlab.m
res = [t77(1) t77(2) t77(4) t77(7) t77(11) t77(16); t77(2) t77(3) t77(5) t77(8) t77(12) t77(17); t77(4) t77(5) t77(6) t77(9) t77(13) t77(18); t77(7) t77(8) t77(9) t77(10) t77(14) t77(19); t77(11) t77(12) t77(13) t77(14) t77(15) t77(20); t77(16) t77(17) t77(18) t77(19) t77(20) t77(21);];
Mq  = res;
