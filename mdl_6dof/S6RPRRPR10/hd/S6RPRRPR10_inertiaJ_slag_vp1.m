% Calculate joint inertia matrix for
% S6RPRRPR10
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% pkin [9x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d3,d4,d6]';
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
% Datum: 2019-03-09 05:38
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Mq = S6RPRRPR10_inertiaJ_slag_vp1(qJ, ...
  pkin, m, rSges, Icges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(9,1),zeros(7,1),zeros(7,3),zeros(7,6)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPRRPR10_inertiaJ_slag_vp1: qJ has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S6RPRRPR10_inertiaJ_slag_vp1: pkin has to be [9x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6RPRRPR10_inertiaJ_slag_vp1: m has to be [7x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [7,3]), ...
  'S6RPRRPR10_inertiaJ_slag_vp1: rSges has to be [7x3] (double)');
assert(isreal(Icges) && all(size(Icges) == [7 6]), ...
  'S6RPRRPR10_inertiaJ_slag_vp1: Icges has to be [7x6] (double)'); 

%% Symbolic Calculation
% From inertia_joint_joint_floatb_twist_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 05:34:29
% EndTime: 2019-03-09 05:34:36
% DurationCPUTime: 3.66s
% Computational Cost: add. (5249->437), mult. (13027->631), div. (0->0), fcn. (15286->8), ass. (0->210)
t184 = sin(qJ(3));
t275 = Icges(4,6) * t184;
t188 = cos(qJ(3));
t260 = t188 / 0.2e1;
t274 = Icges(4,5) * t260 - t275 / 0.2e1;
t183 = sin(qJ(4));
t187 = cos(qJ(4));
t123 = Icges(6,6) * t184 + (Icges(6,5) * t187 + Icges(6,3) * t183) * t188;
t124 = Icges(5,3) * t184 + (Icges(5,5) * t187 - Icges(5,6) * t183) * t188;
t127 = Icges(6,2) * t184 + (Icges(6,4) * t187 + Icges(6,6) * t183) * t188;
t131 = Icges(6,4) * t184 + (Icges(6,1) * t187 + Icges(6,5) * t183) * t188;
t132 = Icges(5,5) * t184 + (Icges(5,1) * t187 - Icges(5,4) * t183) * t188;
t235 = t187 * t188;
t241 = t183 * t188;
t273 = t123 * t241 + (t131 + t132) * t235 + (t124 + t127) * t184;
t222 = m(6) / 0.2e1 + m(7) / 0.2e1;
t272 = 0.2e1 * t222;
t185 = sin(qJ(1));
t189 = cos(qJ(1));
t233 = t189 * t187;
t238 = t185 * t183;
t150 = t184 * t238 - t233;
t237 = t185 * t187;
t151 = t183 * t189 + t184 * t237;
t236 = t185 * t188;
t82 = Icges(6,5) * t151 - Icges(6,6) * t236 + Icges(6,3) * t150;
t86 = Icges(6,4) * t151 - Icges(6,2) * t236 + Icges(6,6) * t150;
t90 = Icges(6,1) * t151 - Icges(6,4) * t236 + Icges(6,5) * t150;
t36 = t150 * t82 + t151 * t90 - t236 * t86;
t239 = t184 * t189;
t152 = t183 * t239 + t237;
t154 = -t184 * t233 + t238;
t234 = t188 * t189;
t83 = Icges(6,5) * t154 + Icges(6,6) * t234 - Icges(6,3) * t152;
t87 = Icges(6,4) * t154 + Icges(6,2) * t234 - Icges(6,6) * t152;
t91 = Icges(6,1) * t154 + Icges(6,4) * t234 - Icges(6,5) * t152;
t37 = t150 * t83 + t151 * t91 - t236 * t87;
t84 = Icges(5,5) * t151 - Icges(5,6) * t150 - Icges(5,3) * t236;
t88 = Icges(5,4) * t151 - Icges(5,2) * t150 - Icges(5,6) * t236;
t92 = Icges(5,1) * t151 - Icges(5,4) * t150 - Icges(5,5) * t236;
t38 = -t150 * t88 + t151 * t92 - t236 * t84;
t85 = Icges(5,5) * t154 + Icges(5,6) * t152 + Icges(5,3) * t234;
t89 = Icges(5,4) * t154 + Icges(5,2) * t152 + Icges(5,6) * t234;
t93 = Icges(5,1) * t154 + Icges(5,4) * t152 + Icges(5,5) * t234;
t39 = -t150 * t89 + t151 * t93 - t236 * t85;
t55 = t123 * t150 - t127 * t236 + t131 * t151;
t128 = Icges(5,6) * t184 + (Icges(5,4) * t187 - Icges(5,2) * t183) * t188;
t56 = -t124 * t236 - t128 * t150 + t132 * t151;
t271 = ((t37 + t39) * t189 + (-t36 - t38) * t185) * t188 + (t55 + t56) * t184;
t40 = -t152 * t82 + t154 * t90 + t234 * t86;
t41 = -t152 * t83 + t154 * t91 + t234 * t87;
t42 = t152 * t88 + t154 * t92 + t234 * t84;
t43 = t152 * t89 + t154 * t93 + t234 * t85;
t57 = -t152 * t123 + t127 * t234 + t154 * t131;
t58 = t124 * t234 + t152 * t128 + t154 * t132;
t270 = ((t41 + t43) * t189 + (-t40 - t42) * t185) * t188 + (t57 + t58) * t184;
t44 = t184 * t86 + (t183 * t82 + t187 * t90) * t188;
t46 = t184 * t84 + (-t183 * t88 + t187 * t92) * t188;
t269 = t44 + t46;
t45 = t184 * t87 + (t183 * t83 + t187 * t91) * t188;
t47 = t184 * t85 + (-t183 * t89 + t187 * t93) * t188;
t268 = t45 + t47;
t267 = (rSges(4,1) * t184 + rSges(4,2) * t188) * t189;
t182 = sin(qJ(6));
t186 = cos(qJ(6));
t104 = t150 * t186 - t151 * t182;
t105 = t150 * t182 + t151 * t186;
t62 = Icges(7,5) * t105 + Icges(7,6) * t104 + Icges(7,3) * t236;
t63 = Icges(7,4) * t105 + Icges(7,2) * t104 + Icges(7,6) * t236;
t64 = Icges(7,1) * t105 + Icges(7,4) * t104 + Icges(7,5) * t236;
t14 = t104 * t63 + t105 * t64 + t236 * t62;
t106 = -t152 * t186 - t154 * t182;
t107 = -t152 * t182 + t154 * t186;
t192 = Icges(7,5) * t107 + Icges(7,6) * t106 - Icges(7,3) * t234;
t216 = Icges(7,6) * t234;
t242 = Icges(7,2) * t106;
t243 = Icges(7,4) * t107;
t193 = -t216 + t242 + t243;
t217 = Icges(7,5) * t234;
t246 = Icges(7,1) * t107;
t194 = Icges(7,4) * t106 - t217 + t246;
t190 = t104 * t193 + t105 * t194 + t192 * t236;
t140 = (-t182 * t187 + t183 * t186) * t188;
t141 = (t182 * t183 + t186 * t187) * t188;
t78 = Icges(7,5) * t141 + Icges(7,6) * t140 - Icges(7,3) * t184;
t79 = Icges(7,4) * t141 + Icges(7,2) * t140 - Icges(7,6) * t184;
t80 = Icges(7,1) * t141 + Icges(7,4) * t140 - Icges(7,5) * t184;
t28 = t104 * t79 + t105 * t80 + t236 * t78;
t1 = t14 * t236 - t28 * t184 - t190 * t234;
t15 = t106 * t63 + t107 * t64 - t234 * t62;
t181 = t189 ^ 2;
t265 = t188 ^ 2;
t191 = t265 * t181 * Icges(7,3) + (-0.2e1 * t217 + t246) * t107 + (-0.2e1 * t216 + t242 + 0.2e1 * t243) * t106;
t29 = t106 * t79 + t107 * t80 - t234 * t78;
t2 = t15 * t236 - t29 * t184 - t191 * t234;
t23 = t140 * t193 + t141 * t194 - t184 * t192;
t21 = t23 * t234;
t22 = t140 * t63 + t141 * t64 - t184 * t62;
t253 = t140 * t79 + t141 * t80;
t32 = (-t184 * t78 + t253) * t184;
t3 = t22 * t236 - t21 - t32;
t266 = (t185 * t1 - t189 * t2) * t188 - t184 * t3;
t180 = t185 ^ 2;
t264 = -pkin(1) - pkin(7);
t261 = t185 / 0.2e1;
t259 = -t189 / 0.2e1;
t258 = t189 / 0.2e1;
t164 = rSges(4,1) * t188 - rSges(4,2) * t184;
t257 = m(4) * t164;
t255 = t184 * (t23 * t185 + t22 * t189);
t254 = (-t128 * t241 + t273) * t184;
t252 = t152 * rSges(6,3);
t110 = t151 * pkin(4) + qJ(5) * t150;
t228 = t151 * rSges(6,1) + t150 * rSges(6,3);
t94 = -rSges(6,2) * t236 + t228;
t251 = -t110 - t94;
t143 = t152 * qJ(5);
t111 = t154 * pkin(4) - t143;
t96 = t154 * rSges(6,1) + rSges(6,2) * t234 - t252;
t250 = -t111 - t96;
t65 = t105 * rSges(7,1) + t104 * rSges(7,2) + rSges(7,3) * t236;
t249 = t151 * pkin(5) + pkin(9) * t236 + t65;
t205 = -t107 * rSges(7,1) - t106 * rSges(7,2);
t66 = -rSges(7,3) * t234 - t205;
t248 = t154 * pkin(5) - pkin(9) * t234 + t66;
t81 = rSges(7,1) * t141 + rSges(7,2) * t140 - rSges(7,3) * t184;
t247 = pkin(5) * t235 - pkin(9) * t184 + t81;
t240 = t184 * t185;
t155 = (pkin(4) * t187 + qJ(5) * t183) * t188;
t232 = t184 * t110 + t155 * t236;
t173 = pkin(3) * t239;
t142 = t189 * (pkin(8) * t234 - t173);
t231 = t189 * t111 + t142;
t166 = pkin(3) * t188 + pkin(8) * t184;
t157 = t185 * t166;
t229 = t185 * t155 + t157;
t227 = t151 * rSges(5,1) - t150 * rSges(5,2);
t226 = -t155 - t166;
t225 = t189 * pkin(1) + t185 * qJ(2);
t224 = t180 + t181;
t223 = -t21 / 0.2e1 - t32;
t221 = pkin(8) * t236;
t220 = t28 / 0.2e1 + t22 / 0.2e1;
t219 = -t110 - t249;
t218 = -t111 - t248;
t214 = rSges(4,1) * t240 + rSges(4,2) * t236 + t189 * rSges(4,3);
t213 = t189 * pkin(7) + t225;
t212 = (-rSges(6,2) - pkin(8)) * t188;
t211 = (-rSges(5,3) - pkin(8)) * t188;
t210 = t247 * t188;
t171 = pkin(3) * t240;
t209 = t171 + t213;
t206 = -t154 * rSges(5,1) - t152 * rSges(5,2);
t202 = Icges(4,5) * t184 + Icges(4,6) * t188;
t176 = t189 * qJ(2);
t199 = t185 * t264 + t173 + t176;
t198 = t143 + t199;
t197 = t29 / 0.2e1 + t45 / 0.2e1 + t58 / 0.2e1 + t57 / 0.2e1 + t47 / 0.2e1;
t196 = t110 + t209;
t195 = -t46 / 0.2e1 - t44 / 0.2e1 - t56 / 0.2e1 - t55 / 0.2e1 - t220;
t165 = rSges(2,1) * t189 - t185 * rSges(2,2);
t163 = -t185 * rSges(2,1) - rSges(2,2) * t189;
t160 = Icges(4,5) * t188 - t275;
t156 = t171 - t221;
t138 = -rSges(3,2) * t189 + t185 * rSges(3,3) + t225;
t137 = rSges(3,3) * t189 + t176 + (rSges(3,2) - pkin(1)) * t185;
t136 = rSges(5,3) * t184 + (rSges(5,1) * t187 - rSges(5,2) * t183) * t188;
t135 = rSges(6,2) * t184 + (rSges(6,1) * t187 + rSges(6,3) * t183) * t188;
t126 = Icges(4,3) * t185 - t189 * t202;
t125 = Icges(4,3) * t189 + t185 * t202;
t122 = t155 * t234;
t113 = t213 + t214;
t112 = t176 + t267 + (-rSges(4,3) + t264) * t185;
t109 = (-t136 - t166) * t189;
t108 = t136 * t185 + t157;
t97 = rSges(5,3) * t234 - t206;
t95 = -rSges(5,3) * t236 + t227;
t77 = -t185 * t214 + (t185 * rSges(4,3) - t267) * t189;
t76 = (-t135 + t226) * t189;
t75 = t135 * t185 + t229;
t72 = t185 * t211 + t209 + t227;
t71 = t189 * t211 + t199 + t206;
t70 = t136 * t234 - t184 * t97;
t69 = t136 * t236 + t184 * t95;
t61 = (t226 - t247) * t189;
t60 = t185 * t247 + t229;
t59 = (-t185 * t97 - t189 * t95) * t188;
t54 = t185 * t212 + t196 + t228;
t53 = t252 + t189 * t212 + (-rSges(6,1) - pkin(4)) * t154 + t198;
t52 = t189 * t97 + t142 + (-t156 - t95) * t185;
t51 = t135 * t234 + t184 * t250 + t122;
t50 = t135 * t236 + t184 * t94 + t232;
t49 = t184 * t66 - t234 * t81;
t48 = -t184 * t65 - t236 * t81;
t35 = t196 - t221 + t249;
t34 = (-pkin(4) - pkin(5)) * t154 + (rSges(7,3) - pkin(8) + pkin(9)) * t234 + t198 + t205;
t33 = (t185 * t250 + t189 * t251) * t188;
t31 = (t185 * t66 + t189 * t65) * t188;
t30 = t189 * t96 + (-t156 + t251) * t185 + t231;
t25 = t184 * t218 + t189 * t210 + t122;
t24 = t184 * t249 + t185 * t210 + t232;
t20 = (t185 * t218 + t189 * t219) * t188;
t19 = t43 * t185 + t189 * t42;
t18 = t41 * t185 + t189 * t40;
t17 = t39 * t185 + t189 * t38;
t16 = t37 * t185 + t189 * t36;
t13 = t248 * t189 + (-t156 + t219) * t185 + t231;
t5 = t15 * t189 + t185 * t191;
t4 = t14 * t189 + t185 * t190;
t6 = [Icges(3,1) + Icges(2,3) + (Icges(4,1) * t188 - t183 * t128) * t188 + m(7) * (t34 ^ 2 + t35 ^ 2) + m(6) * (t53 ^ 2 + t54 ^ 2) + m(5) * (t71 ^ 2 + t72 ^ 2) + m(4) * (t112 ^ 2 + t113 ^ 2) + m(3) * (t137 ^ 2 + t138 ^ 2) + m(2) * (t163 ^ 2 + t165 ^ 2) + t253 + t273 + (-0.2e1 * Icges(4,4) * t188 + Icges(4,2) * t184 - t78) * t184; m(7) * (t185 * t34 - t189 * t35) + m(6) * (t185 * t53 - t189 * t54) + m(5) * (t185 * t71 - t189 * t72) + m(4) * (t185 * t112 - t113 * t189) + m(3) * (t185 * t137 - t138 * t189); 0.2e1 * (m(3) / 0.2e1 + m(4) / 0.2e1 + m(5) / 0.2e1 + t222) * t224; m(7) * (t34 * t60 + t35 * t61) + m(6) * (t53 * t75 + t54 * t76) + m(5) * (t108 * t71 + t109 * t72) + (-t113 * t257 + t160 * t258 + t189 * t274 - t195) * t189 + (t23 / 0.2e1 + t112 * t257 + t160 * t261 + t197 + t274 * t185) * t185; m(5) * (t108 * t185 - t109 * t189) + m(6) * (t75 * t185 - t189 * t76) + m(7) * (t60 * t185 - t189 * t61) + t224 * t257; m(7) * (t13 ^ 2 + t60 ^ 2 + t61 ^ 2) + m(6) * (t30 ^ 2 + t75 ^ 2 + t76 ^ 2) + m(5) * (t108 ^ 2 + t109 ^ 2 + t52 ^ 2) + m(4) * (t164 ^ 2 * t224 + t77 ^ 2) + (t181 * t125 + t16 + t17 + t4) * t189 + (t180 * t126 + t18 + t19 + t5 + (t185 * t125 + t189 * t126) * t189) * t185; m(7) * (t24 * t35 + t25 * t34) + m(6) * (t50 * t54 + t51 * t53) + m(5) * (t69 * t72 + t70 * t71) + (t185 * t195 + t189 * t197) * t188 - t223 + t254; m(5) * (t70 * t185 - t189 * t69) + m(6) * (t51 * t185 - t189 * t50) + m(7) * (t25 * t185 - t189 * t24); -t185 * t2 / 0.2e1 + t1 * t259 + t255 / 0.2e1 + m(7) * (t13 * t20 + t24 * t61 + t25 * t60) + m(6) * (t33 * t30 + t50 * t76 + t51 * t75) + m(5) * (t108 * t70 + t109 * t69 + t59 * t52) + ((t5 / 0.2e1 + t19 / 0.2e1 + t18 / 0.2e1) * t189 + (-t4 / 0.2e1 - t17 / 0.2e1 - t16 / 0.2e1) * t185) * t188 + (t268 * t185 + t269 * t189) * t184 / 0.2e1 + t270 * t261 + t271 * t258; (-t3 + t254) * t184 + m(7) * (t20 ^ 2 + t24 ^ 2 + t25 ^ 2) + m(6) * (t33 ^ 2 + t50 ^ 2 + t51 ^ 2) + m(5) * (t59 ^ 2 + t69 ^ 2 + t70 ^ 2) + ((t268 * t184 - t2 + t270) * t189 + (-t269 * t184 + t1 - t271) * t185) * t188; m(7) * (t150 * t34 - t152 * t35) + m(6) * (t150 * t53 - t152 * t54); (t150 * t185 + t152 * t189) * t272; m(7) * (t13 * t241 + t150 * t60 - t152 * t61) + m(6) * (t150 * t75 - t152 * t76 + t241 * t30); m(7) * (t150 * t25 - t152 * t24 + t20 * t241) + m(6) * (t150 * t51 - t152 * t50 + t241 * t33); (t183 ^ 2 * t265 + t150 ^ 2 + t152 ^ 2) * t272; m(7) * (t34 * t49 + t35 * t48) + (t185 * t220 + t259 * t29) * t188 + t223; m(7) * (t49 * t185 - t189 * t48); -t255 / 0.2e1 + m(7) * (t13 * t31 + t48 * t61 + t49 * t60) + (t1 / 0.2e1 - t188 * t5 / 0.2e1) * t189 + (t2 / 0.2e1 + t4 * t260) * t185; m(7) * (t20 * t31 + t24 * t48 + t25 * t49) - t266; m(7) * (t150 * t49 - t152 * t48 + t241 * t31); m(7) * (t31 ^ 2 + t48 ^ 2 + t49 ^ 2) + t266;];
%% Postprocessing: Reshape Output
% From vec2symmat_6_matlab.m
res = [t6(1) t6(2) t6(4) t6(7) t6(11) t6(16); t6(2) t6(3) t6(5) t6(8) t6(12) t6(17); t6(4) t6(5) t6(6) t6(9) t6(13) t6(18); t6(7) t6(8) t6(9) t6(10) t6(14) t6(19); t6(11) t6(12) t6(13) t6(14) t6(15) t6(20); t6(16) t6(17) t6(18) t6(19) t6(20) t6(21);];
Mq  = res;
