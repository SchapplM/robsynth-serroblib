% Calculate joint inertia matrix for
% S6RRPRPP1
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% pkin [10x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d2,d4,theta3,theta5]';
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
% Datum: 2019-03-09 09:48
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Mq = S6RRPRPP1_inertiaJ_slag_vp1(qJ, ...
  pkin, m, rSges, Icges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(10,1),zeros(7,1),zeros(7,3),zeros(7,6)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRPRPP1_inertiaJ_slag_vp1: qJ has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S6RRPRPP1_inertiaJ_slag_vp1: pkin has to be [10x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6RRPRPP1_inertiaJ_slag_vp1: m has to be [7x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [7,3]), ...
  'S6RRPRPP1_inertiaJ_slag_vp1: rSges has to be [7x3] (double)');
assert(isreal(Icges) && all(size(Icges) == [7 6]), ...
  'S6RRPRPP1_inertiaJ_slag_vp1: Icges has to be [7x6] (double)'); 

%% Symbolic Calculation
% From inertia_joint_joint_floatb_twist_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 09:45:10
% EndTime: 2019-03-09 09:45:19
% DurationCPUTime: 3.63s
% Computational Cost: add. (8042->448), mult. (8092->649), div. (0->0), fcn. (8605->10), ass. (0->214)
t294 = Icges(3,3) + Icges(4,3);
t192 = qJ(2) + pkin(9);
t185 = sin(t192);
t187 = cos(t192);
t198 = sin(qJ(2));
t201 = cos(qJ(2));
t293 = Icges(3,5) * t201 + Icges(4,5) * t187 - Icges(3,6) * t198 - Icges(4,6) * t185;
t191 = qJ(4) + pkin(10);
t186 = cos(t191);
t202 = cos(qJ(1));
t184 = sin(t191);
t199 = sin(qJ(1));
t247 = t199 * t184;
t139 = t186 * t202 + t187 * t247;
t246 = t199 * t186;
t140 = -t184 * t202 + t187 * t246;
t290 = rSges(7,3) + qJ(6);
t291 = rSges(7,1) + pkin(5);
t292 = -t290 * t139 - t291 * t140;
t107 = -Icges(6,3) * t187 + (Icges(6,5) * t186 - Icges(6,6) * t184) * t185;
t108 = -Icges(7,2) * t187 + (Icges(7,4) * t186 + Icges(7,6) * t184) * t185;
t197 = sin(qJ(4));
t200 = cos(qJ(4));
t118 = -Icges(5,3) * t187 + (Icges(5,5) * t200 - Icges(5,6) * t197) * t185;
t289 = -t107 - t108 - t118;
t288 = t185 / 0.2e1;
t287 = t187 / 0.2e1;
t286 = t198 / 0.2e1;
t285 = t201 / 0.2e1;
t284 = -t293 * t199 + t202 * t294;
t283 = t199 * t294 + t293 * t202;
t252 = t185 * t199;
t64 = Icges(7,5) * t140 + Icges(7,6) * t252 + Icges(7,3) * t139;
t68 = Icges(7,4) * t140 + Icges(7,2) * t252 + Icges(7,6) * t139;
t72 = Icges(7,1) * t140 + Icges(7,4) * t252 + Icges(7,5) * t139;
t19 = t139 * t64 + t140 * t72 + t252 * t68;
t250 = t187 * t202;
t141 = t184 * t250 - t246;
t142 = t186 * t250 + t247;
t251 = t185 * t202;
t65 = Icges(7,5) * t142 + Icges(7,6) * t251 + Icges(7,3) * t141;
t69 = Icges(7,4) * t142 + Icges(7,2) * t251 + Icges(7,6) * t141;
t73 = Icges(7,1) * t142 + Icges(7,4) * t251 + Icges(7,5) * t141;
t20 = t139 * t65 + t140 * t73 + t252 * t69;
t66 = Icges(6,5) * t140 - Icges(6,6) * t139 + Icges(6,3) * t252;
t70 = Icges(6,4) * t140 - Icges(6,2) * t139 + Icges(6,6) * t252;
t74 = Icges(6,1) * t140 - Icges(6,4) * t139 + Icges(6,5) * t252;
t21 = -t139 * t70 + t140 * t74 + t252 * t66;
t67 = Icges(6,5) * t142 - Icges(6,6) * t141 + Icges(6,3) * t251;
t71 = Icges(6,4) * t142 - Icges(6,2) * t141 + Icges(6,6) * t251;
t75 = Icges(6,1) * t142 - Icges(6,4) * t141 + Icges(6,5) * t251;
t22 = -t139 * t71 + t140 * t75 + t252 * t67;
t243 = t200 * t202;
t245 = t199 * t197;
t151 = -t187 * t245 - t243;
t244 = t199 * t200;
t248 = t197 * t202;
t152 = t187 * t244 - t248;
t87 = Icges(5,5) * t152 + Icges(5,6) * t151 + Icges(5,3) * t252;
t89 = Icges(5,4) * t152 + Icges(5,2) * t151 + Icges(5,6) * t252;
t91 = Icges(5,1) * t152 + Icges(5,4) * t151 + Icges(5,5) * t252;
t34 = t151 * t89 + t152 * t91 + t252 * t87;
t153 = -t187 * t248 + t244;
t154 = t187 * t243 + t245;
t88 = Icges(5,5) * t154 + Icges(5,6) * t153 + Icges(5,3) * t251;
t90 = Icges(5,4) * t154 + Icges(5,2) * t153 + Icges(5,6) * t251;
t92 = Icges(5,1) * t154 + Icges(5,4) * t153 + Icges(5,5) * t251;
t35 = t151 * t90 + t152 * t92 + t252 * t88;
t106 = -Icges(7,6) * t187 + (Icges(7,5) * t186 + Icges(7,3) * t184) * t185;
t110 = -Icges(7,4) * t187 + (Icges(7,1) * t186 + Icges(7,5) * t184) * t185;
t42 = t106 * t139 + t108 * t252 + t110 * t140;
t109 = -Icges(6,6) * t187 + (Icges(6,4) * t186 - Icges(6,2) * t184) * t185;
t111 = -Icges(6,5) * t187 + (Icges(6,1) * t186 - Icges(6,4) * t184) * t185;
t43 = t107 * t252 - t109 * t139 + t111 * t140;
t119 = -Icges(5,6) * t187 + (Icges(5,4) * t200 - Icges(5,2) * t197) * t185;
t120 = -Icges(5,5) * t187 + (Icges(5,1) * t200 - Icges(5,4) * t197) * t185;
t46 = t118 * t252 + t119 * t151 + t120 * t152;
t282 = (-t42 - t43 - t46) * t187 + ((t20 + t22 + t35) * t202 + (t19 + t21 + t34) * t199) * t185;
t23 = t141 * t64 + t142 * t72 + t251 * t68;
t24 = t141 * t65 + t142 * t73 + t251 * t69;
t25 = -t141 * t70 + t142 * t74 + t251 * t66;
t26 = -t141 * t71 + t142 * t75 + t251 * t67;
t36 = t153 * t89 + t154 * t91 + t251 * t87;
t37 = t153 * t90 + t154 * t92 + t251 * t88;
t44 = t141 * t106 + t108 * t251 + t142 * t110;
t45 = t107 * t251 - t141 * t109 + t142 * t111;
t47 = t118 * t251 + t153 * t119 + t154 * t120;
t281 = (-t44 - t45 - t47) * t187 + ((t24 + t26 + t37) * t202 + (t23 + t25 + t36) * t199) * t185;
t28 = -t187 * t68 + (t184 * t64 + t186 * t72) * t185;
t30 = -t187 * t66 + (-t184 * t70 + t186 * t74) * t185;
t38 = -t187 * t87 + (-t197 * t89 + t200 * t91) * t185;
t280 = -t28 - t30 - t38;
t29 = -t187 * t69 + (t184 * t65 + t186 * t73) * t185;
t31 = -t187 * t67 + (-t184 * t71 + t186 * t75) * t185;
t39 = -t187 * t88 + (-t197 * t90 + t200 * t92) * t185;
t279 = t29 + t31 + t39;
t254 = t184 * t185;
t278 = t106 * t254 + (t200 * t120 + (t110 + t111) * t186) * t185;
t277 = t187 ^ 2;
t193 = t199 ^ 2;
t194 = t202 ^ 2;
t276 = m(6) / 0.2e1;
t275 = m(7) / 0.2e1;
t274 = -t187 / 0.2e1;
t271 = pkin(2) * t198;
t270 = pkin(3) * t187;
t269 = pkin(8) * t185;
t181 = pkin(4) * t200 + pkin(3);
t268 = -pkin(3) + t181;
t267 = rSges(7,2) * t252 - t292;
t266 = rSges(7,2) * t251 + t290 * t141 + t291 * t142;
t79 = t142 * rSges(6,1) - t141 * rSges(6,2) + rSges(6,3) * t251;
t195 = -qJ(5) - pkin(8);
t208 = pkin(4) * t245 + t181 * t250 - t195 * t251;
t239 = pkin(3) * t250 + pkin(8) * t251;
t86 = t208 - t239;
t265 = -t79 - t86;
t105 = (pkin(8) + t195) * t187 + t268 * t185;
t240 = -pkin(4) * t248 - t195 * t252;
t85 = (t187 * t268 - t269) * t199 + t240;
t264 = t105 * t252 + t187 * t85;
t263 = rSges(3,1) * t201;
t262 = rSges(3,2) * t198;
t261 = t202 * rSges(3,3);
t259 = Icges(3,4) * t198;
t258 = Icges(3,4) * t201;
t257 = Icges(4,4) * t185;
t256 = Icges(4,4) * t187;
t255 = t119 * t197;
t196 = -qJ(3) - pkin(7);
t249 = t196 * t202;
t242 = -t187 * rSges(7,2) + (t290 * t184 + t291 * t186) * t185;
t182 = pkin(2) * t201 + pkin(1);
t175 = t202 * t182;
t190 = t202 * pkin(7);
t241 = t199 * (t249 + t190 + (-pkin(1) + t182) * t199) + t202 * (-pkin(1) * t202 + t175 + (-pkin(7) - t196) * t199);
t238 = t199 * rSges(3,3) + t202 * t263;
t237 = t193 + t194;
t236 = t276 + t275;
t235 = -t109 * t254 - t185 * t255 + t289 * t187 + t278;
t234 = -t86 - t266;
t94 = t154 * rSges(5,1) + t153 * rSges(5,2) + rSges(5,3) * t251;
t233 = Icges(3,5) * t286 + Icges(4,5) * t288 + Icges(3,6) * t285 + Icges(4,6) * t287;
t232 = -rSges(4,1) * t185 - rSges(4,2) * t187 - t271;
t231 = -pkin(3) * t185 + pkin(8) * t187 - t271;
t230 = -t199 * t196 + t175;
t229 = -t181 * t187 - t182;
t228 = t193 * (t269 + t270) + t202 * t239 + t241;
t227 = -t105 + t231;
t121 = -t187 * rSges(5,3) + (rSges(5,1) * t200 - rSges(5,2) * t197) * t185;
t226 = -t121 + t231;
t225 = -t240 - t249;
t224 = -t262 + t263;
t223 = rSges(4,1) * t187 - rSges(4,2) * t185;
t222 = -t152 * rSges(5,1) - t151 * rSges(5,2);
t221 = -t140 * rSges(6,1) + t139 * rSges(6,2);
t220 = Icges(3,1) * t201 - t259;
t219 = Icges(4,1) * t187 - t257;
t218 = -Icges(3,2) * t198 + t258;
t217 = -Icges(4,2) * t185 + t256;
t113 = -t187 * rSges(6,3) + (rSges(6,1) * t186 - rSges(6,2) * t184) * t185;
t210 = -t113 + t227;
t209 = rSges(4,1) * t250 - rSges(4,2) * t251 + t199 * rSges(4,3);
t207 = t199 * t85 + t202 * t86 + t228;
t206 = t227 - t242;
t205 = t30 / 0.2e1 + t28 / 0.2e1 + t46 / 0.2e1 + t43 / 0.2e1 + t42 / 0.2e1 + t38 / 0.2e1;
t204 = t45 / 0.2e1 + t44 / 0.2e1 + t39 / 0.2e1 + t31 / 0.2e1 + t29 / 0.2e1 + t47 / 0.2e1;
t203 = t208 + t230;
t183 = t185 ^ 2;
t168 = rSges(2,1) * t202 - t199 * rSges(2,2);
t167 = -t199 * rSges(2,1) - rSges(2,2) * t202;
t166 = rSges(3,1) * t198 + rSges(3,2) * t201;
t123 = t232 * t202;
t122 = t232 * t199;
t117 = t199 * pkin(7) + (pkin(1) - t262) * t202 + t238;
t116 = t261 + t190 + (-pkin(1) - t224) * t199;
t103 = t209 + t230;
t102 = (rSges(4,3) - t196) * t202 + (-t182 - t223) * t199;
t97 = t202 * (-t202 * t262 + t238) + (t199 * t224 - t261) * t199;
t93 = rSges(5,3) * t252 - t222;
t84 = t226 * t202;
t83 = t226 * t199;
t77 = rSges(6,3) * t252 - t221;
t63 = t85 * t251;
t62 = t230 + t94 + t239;
t61 = -t249 + (-t270 - t182 + (-rSges(5,3) - pkin(8)) * t185) * t199 + t222;
t60 = -t121 * t251 - t187 * t94;
t59 = t121 * t252 + t187 * t93;
t58 = t210 * t202;
t57 = t210 * t199;
t55 = t203 + t79;
t54 = (-rSges(6,3) * t185 + t229) * t199 + t221 + t225;
t53 = t202 * t209 + (-t202 * rSges(4,3) + t199 * t223) * t199 + t241;
t52 = (-t199 * t94 + t202 * t93) * t185;
t49 = t206 * t202;
t48 = t206 * t199;
t41 = t203 + t266;
t40 = (-rSges(7,2) * t185 + t229) * t199 + t225 + t292;
t33 = t265 * t187 + (-t105 - t113) * t251;
t32 = t113 * t252 + t187 * t77 + t264;
t27 = t199 * t93 + t202 * t94 + t228;
t18 = t63 + (t199 * t265 + t202 * t77) * t185;
t17 = t234 * t187 + (-t105 - t242) * t251;
t16 = t187 * t267 + t242 * t252 + t264;
t15 = t199 * t77 + t202 * t79 + t207;
t14 = t63 + (t199 * t234 + t202 * t267) * t185;
t13 = t37 * t199 - t202 * t36;
t12 = t35 * t199 - t202 * t34;
t11 = t199 * t267 + t202 * t266 + t207;
t10 = t26 * t199 - t202 * t25;
t9 = t24 * t199 - t202 * t23;
t8 = t22 * t199 - t202 * t21;
t7 = -t19 * t202 + t20 * t199;
t1 = [t201 * (Icges(3,2) * t201 + t259) + t198 * (Icges(3,1) * t198 + t258) + Icges(2,3) + (Icges(4,1) * t185 - t109 * t184 - t255 + t256) * t185 + (Icges(4,2) * t187 + t257 + t289) * t187 + m(7) * (t40 ^ 2 + t41 ^ 2) + m(6) * (t54 ^ 2 + t55 ^ 2) + m(5) * (t61 ^ 2 + t62 ^ 2) + m(4) * (t102 ^ 2 + t103 ^ 2) + m(3) * (t116 ^ 2 + t117 ^ 2) + m(2) * (t167 ^ 2 + t168 ^ 2) + t278; ((-Icges(4,6) * t202 + t199 * t217) * t274 - t185 * (-Icges(4,5) * t202 + t199 * t219) / 0.2e1 - t201 * (-Icges(3,6) * t202 + t199 * t218) / 0.2e1 - t198 * (-Icges(3,5) * t202 + t199 * t220) / 0.2e1 + t233 * t202 - t205) * t202 + ((Icges(4,6) * t199 + t202 * t217) * t287 + (Icges(4,5) * t199 + t202 * t219) * t288 + (Icges(3,6) * t199 + t202 * t218) * t285 + (Icges(3,5) * t199 + t202 * t220) * t286 + t233 * t199 + t204) * t199 + m(7) * (t40 * t49 + t41 * t48) + m(6) * (t54 * t58 + t55 * t57) + m(5) * (t61 * t84 + t62 * t83) + m(4) * (t102 * t123 + t103 * t122) + m(3) * (-t116 * t202 - t117 * t199) * t166; m(6) * (t15 ^ 2 + t57 ^ 2 + t58 ^ 2) + m(7) * (t11 ^ 2 + t48 ^ 2 + t49 ^ 2) + m(5) * (t27 ^ 2 + t83 ^ 2 + t84 ^ 2) + m(4) * (t122 ^ 2 + t123 ^ 2 + t53 ^ 2) + m(3) * (t166 ^ 2 * t237 + t97 ^ 2) + (t194 * t284 - t12 - t7 - t8) * t202 + (t10 + t13 + t9 + t283 * t193 + (t199 * t284 + t202 * t283) * t202) * t199; m(7) * (t199 * t40 - t202 * t41) + m(6) * (t199 * t54 - t202 * t55) + m(5) * (t199 * t61 - t202 * t62) + m(4) * (t199 * t102 - t103 * t202); m(6) * (t199 * t58 - t202 * t57) + m(7) * (t199 * t49 - t202 * t48) + m(5) * (t199 * t84 - t202 * t83) + m(4) * (-t122 * t202 + t199 * t123); 0.2e1 * (m(4) / 0.2e1 + m(5) / 0.2e1 + t236) * t237; -t235 * t187 + m(7) * (t16 * t40 + t17 * t41) + m(6) * (t32 * t54 + t33 * t55) + m(5) * (t59 * t61 + t60 * t62) + (t199 * t205 + t202 * t204) * t185; m(6) * (t15 * t18 + t32 * t58 + t33 * t57) + m(7) * (t11 * t14 + t16 * t49 + t17 * t48) + m(5) * (t27 * t52 + t59 * t84 + t60 * t83) + ((t13 / 0.2e1 + t10 / 0.2e1 + t9 / 0.2e1) * t202 + (t8 / 0.2e1 + t12 / 0.2e1 + t7 / 0.2e1) * t199) * t185 + (t199 * t279 + t202 * t280) * t274 + t281 * t199 / 0.2e1 - t282 * t202 / 0.2e1; m(5) * (t59 * t199 - t202 * t60) + m(6) * (t32 * t199 - t202 * t33) + m(7) * (t16 * t199 - t17 * t202); m(7) * (t14 ^ 2 + t16 ^ 2 + t17 ^ 2) + m(6) * (t18 ^ 2 + t32 ^ 2 + t33 ^ 2) + m(5) * (t52 ^ 2 + t59 ^ 2 + t60 ^ 2) + t235 * t277 + (t281 * t202 + t282 * t199 + (t199 * t280 - t202 * t279) * t187) * t185; 0.2e1 * ((t199 * t41 + t202 * t40) * t275 + (t199 * t55 + t202 * t54) * t276) * t185; m(6) * (-t187 * t15 + (t199 * t57 + t202 * t58) * t185) + m(7) * (-t187 * t11 + (t199 * t48 + t202 * t49) * t185); 0; m(7) * (-t187 * t14 + (t16 * t202 + t17 * t199) * t185) + m(6) * (-t187 * t18 + (t199 * t33 + t202 * t32) * t185); 0.2e1 * t236 * (t183 * t237 + t277); m(7) * (t139 * t41 + t141 * t40); m(7) * (t11 * t254 + t139 * t48 + t141 * t49); m(7) * (-t139 * t202 + t141 * t199); m(7) * (t139 * t17 + t14 * t254 + t141 * t16); m(7) * (t139 * t199 + t141 * t202 - t184 * t187) * t185; m(7) * (t183 * t184 ^ 2 + t139 ^ 2 + t141 ^ 2);];
%% Postprocessing: Reshape Output
% From vec2symmat_6_matlab.m
res = [t1(1) t1(2) t1(4) t1(7) t1(11) t1(16); t1(2) t1(3) t1(5) t1(8) t1(12) t1(17); t1(4) t1(5) t1(6) t1(9) t1(13) t1(18); t1(7) t1(8) t1(9) t1(10) t1(14) t1(19); t1(11) t1(12) t1(13) t1(14) t1(15) t1(20); t1(16) t1(17) t1(18) t1(19) t1(20) t1(21);];
Mq  = res;
