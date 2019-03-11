% Calculate joint inertia matrix for
% S6RRPRPP3
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% pkin [9x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d2,d4,theta3]';
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
% Datum: 2019-03-09 09:57
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Mq = S6RRPRPP3_inertiaJ_slag_vp1(qJ, ...
  pkin, m, rSges, Icges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(9,1),zeros(7,1),zeros(7,3),zeros(7,6)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRPRPP3_inertiaJ_slag_vp1: qJ has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S6RRPRPP3_inertiaJ_slag_vp1: pkin has to be [9x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6RRPRPP3_inertiaJ_slag_vp1: m has to be [7x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [7,3]), ...
  'S6RRPRPP3_inertiaJ_slag_vp1: rSges has to be [7x3] (double)');
assert(isreal(Icges) && all(size(Icges) == [7 6]), ...
  'S6RRPRPP3_inertiaJ_slag_vp1: Icges has to be [7x6] (double)'); 

%% Symbolic Calculation
% From inertia_joint_joint_floatb_twist_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 09:54:16
% EndTime: 2019-03-09 09:54:25
% DurationCPUTime: 3.50s
% Computational Cost: add. (6430->483), mult. (8958->702), div. (0->0), fcn. (9687->8), ass. (0->212)
t202 = sin(qJ(1));
t266 = t202 / 0.2e1;
t204 = cos(qJ(1));
t281 = t204 / 0.2e1;
t280 = rSges(7,1) + pkin(5);
t279 = rSges(7,3) + qJ(6);
t194 = pkin(9) + qJ(4);
t188 = sin(t194);
t189 = cos(t194);
t201 = sin(qJ(2));
t203 = cos(qJ(2));
t120 = -Icges(5,3) * t203 + (Icges(5,5) * t189 - Icges(5,6) * t188) * t201;
t127 = -Icges(7,1) * t203 + (Icges(7,4) * t188 + Icges(7,5) * t189) * t201;
t128 = -Icges(6,1) * t203 + (-Icges(6,4) * t189 + Icges(6,5) * t188) * t201;
t278 = -t120 - t127 - t128;
t244 = t202 * t203;
t144 = t188 * t244 + t189 * t204;
t242 = t204 * t188;
t145 = t189 * t244 - t242;
t247 = t201 * t202;
t68 = Icges(7,5) * t247 + Icges(7,6) * t144 + Icges(7,3) * t145;
t72 = Icges(7,4) * t247 + Icges(7,2) * t144 + Icges(7,6) * t145;
t76 = Icges(7,1) * t247 + Icges(7,4) * t144 + Icges(7,5) * t145;
t20 = t144 * t72 + t145 * t68 + t247 * t76;
t146 = -t202 * t189 + t203 * t242;
t243 = t203 * t204;
t147 = t202 * t188 + t189 * t243;
t246 = t201 * t204;
t69 = Icges(7,5) * t246 + Icges(7,6) * t146 + Icges(7,3) * t147;
t73 = Icges(7,4) * t246 + Icges(7,2) * t146 + Icges(7,6) * t147;
t77 = Icges(7,1) * t246 + Icges(7,4) * t146 + Icges(7,5) * t147;
t21 = t144 * t73 + t145 * t69 + t247 * t77;
t70 = Icges(6,5) * t247 - Icges(6,6) * t145 + Icges(6,3) * t144;
t74 = Icges(6,4) * t247 - Icges(6,2) * t145 + Icges(6,6) * t144;
t78 = Icges(6,1) * t247 - Icges(6,4) * t145 + Icges(6,5) * t144;
t22 = t144 * t70 - t145 * t74 + t247 * t78;
t71 = Icges(6,5) * t246 - Icges(6,6) * t147 + Icges(6,3) * t146;
t75 = Icges(6,4) * t246 - Icges(6,2) * t147 + Icges(6,6) * t146;
t79 = Icges(6,1) * t246 - Icges(6,4) * t147 + Icges(6,5) * t146;
t23 = t144 * t71 - t145 * t75 + t247 * t79;
t80 = Icges(5,5) * t145 - Icges(5,6) * t144 + Icges(5,3) * t247;
t82 = Icges(5,4) * t145 - Icges(5,2) * t144 + Icges(5,6) * t247;
t84 = Icges(5,1) * t145 - Icges(5,4) * t144 + Icges(5,5) * t247;
t28 = -t144 * t82 + t145 * t84 + t247 * t80;
t81 = Icges(5,5) * t147 - Icges(5,6) * t146 + Icges(5,3) * t246;
t83 = Icges(5,4) * t147 - Icges(5,2) * t146 + Icges(5,6) * t246;
t85 = Icges(5,1) * t147 - Icges(5,4) * t146 + Icges(5,5) * t246;
t29 = -t144 * t83 + t145 * t85 + t247 * t81;
t123 = -Icges(7,5) * t203 + (Icges(7,6) * t188 + Icges(7,3) * t189) * t201;
t125 = -Icges(7,4) * t203 + (Icges(7,2) * t188 + Icges(7,6) * t189) * t201;
t45 = t123 * t145 + t125 * t144 + t127 * t247;
t124 = -Icges(6,5) * t203 + (-Icges(6,6) * t189 + Icges(6,3) * t188) * t201;
t126 = -Icges(6,4) * t203 + (-Icges(6,2) * t189 + Icges(6,6) * t188) * t201;
t46 = t124 * t144 - t126 * t145 + t128 * t247;
t121 = -Icges(5,6) * t203 + (Icges(5,4) * t189 - Icges(5,2) * t188) * t201;
t122 = -Icges(5,5) * t203 + (Icges(5,1) * t189 - Icges(5,4) * t188) * t201;
t49 = t120 * t247 - t121 * t144 + t122 * t145;
t277 = (-t45 - t46 - t49) * t203 + ((t21 + t23 + t29) * t204 + (t20 + t22 + t28) * t202) * t201;
t24 = t146 * t72 + t147 * t68 + t246 * t76;
t25 = t146 * t73 + t147 * t69 + t246 * t77;
t26 = t146 * t70 - t147 * t74 + t246 * t78;
t27 = t146 * t71 - t147 * t75 + t246 * t79;
t30 = -t146 * t82 + t147 * t84 + t246 * t80;
t31 = -t146 * t83 + t147 * t85 + t246 * t81;
t47 = t147 * t123 + t146 * t125 + t127 * t246;
t48 = t146 * t124 - t147 * t126 + t128 * t246;
t50 = t120 * t246 - t146 * t121 + t147 * t122;
t276 = (-t47 - t48 - t50) * t203 + ((t25 + t27 + t31) * t204 + (t24 + t26 + t30) * t202) * t201;
t32 = -t203 * t80 + (-t188 * t82 + t189 * t84) * t201;
t34 = -t203 * t76 + (t188 * t72 + t189 * t68) * t201;
t36 = -t203 * t78 + (t188 * t70 - t189 * t74) * t201;
t275 = -t32 - t34 - t36;
t33 = -t203 * t81 + (-t188 * t83 + t189 * t85) * t201;
t35 = -t203 * t77 + (t188 * t73 + t189 * t69) * t201;
t37 = -t203 * t79 + (t188 * t71 - t189 * t75) * t201;
t274 = t33 + t35 + t37;
t249 = t189 * t201;
t250 = t188 * t201;
t273 = (t124 + t125) * t250 + (t122 + t123) * t249;
t196 = t202 ^ 2;
t272 = t203 ^ 2;
t197 = t204 ^ 2;
t271 = 0.2e1 * t201;
t270 = m(4) / 0.2e1;
t269 = m(5) / 0.2e1;
t268 = m(6) / 0.2e1;
t267 = m(7) / 0.2e1;
t264 = -t204 / 0.2e1;
t171 = rSges(3,1) * t201 + rSges(3,2) * t203;
t263 = m(3) * t171;
t262 = pkin(2) * t203;
t199 = cos(pkin(9));
t187 = pkin(3) * t199 + pkin(2);
t261 = -pkin(2) + t187;
t260 = rSges(7,2) * t144;
t259 = rSges(6,3) * t144;
t258 = t204 * rSges(3,3);
t101 = t147 * pkin(4) + t146 * qJ(5);
t89 = rSges(6,1) * t246 - t147 * rSges(6,2) + t146 * rSges(6,3);
t257 = -t101 - t89;
t256 = t279 * t145 + t280 * t247 + t260;
t255 = t146 * rSges(7,2) + t279 * t147 + t280 * t246;
t132 = t144 * qJ(5);
t100 = pkin(4) * t145 + t132;
t154 = (pkin(4) * t189 + qJ(5) * t188) * t201;
t254 = t203 * t100 + t154 * t247;
t253 = Icges(3,4) * t201;
t252 = Icges(3,4) * t203;
t251 = qJ(3) * t201;
t198 = sin(pkin(9));
t248 = t198 * t204;
t245 = t202 * t198;
t170 = pkin(2) * t201 - qJ(3) * t203;
t200 = -pkin(8) - qJ(3);
t240 = -(qJ(3) + t200) * t203 - t261 * t201 - t170;
t239 = (rSges(7,2) * t188 + rSges(7,3) * t189) * t201 + qJ(6) * t249 - t280 * t203;
t131 = -t203 * rSges(6,1) + (-rSges(6,2) * t189 + rSges(6,3) * t188) * t201;
t238 = -t131 - t154;
t237 = t203 * rSges(4,3) - (rSges(4,1) * t199 - rSges(4,2) * t198) * t201 - t170;
t234 = pkin(2) * t243 + qJ(3) * t246;
t236 = t196 * (t251 + t262) + t204 * t234;
t235 = -pkin(3) * t248 - t200 * t247;
t233 = t204 * pkin(1) + t202 * pkin(7);
t232 = t196 + t197;
t231 = t268 + t267;
t230 = -t121 * t250 - t126 * t249 + t278 * t203 + t273;
t229 = -t101 - t255;
t129 = -t203 * rSges(5,3) + (rSges(5,1) * t189 - rSges(5,2) * t188) * t201;
t228 = -t129 + t240;
t227 = -t154 - t239;
t91 = t147 * rSges(5,1) - t146 * rSges(5,2) + rSges(5,3) * t246;
t162 = -t198 * t243 + t202 * t199;
t163 = t199 * t243 + t245;
t226 = t163 * rSges(4,1) + t162 * rSges(4,2) + rSges(4,3) * t246;
t192 = t204 * pkin(7);
t225 = t192 - t235;
t224 = -t187 * t203 - pkin(1);
t209 = pkin(3) * t245 + t187 * t243 - t200 * t246;
t223 = t202 * ((t203 * t261 - t251) * t202 + t235) + t204 * (t209 - t234) + t236;
t222 = t238 + t240;
t221 = -t132 + t225;
t220 = rSges(3,1) * t203 - rSges(3,2) * t201;
t160 = -t198 * t244 - t199 * t204;
t161 = t199 * t244 - t248;
t219 = -rSges(4,1) * t161 - rSges(4,2) * t160;
t218 = -rSges(5,1) * t145 + rSges(5,2) * t144;
t217 = t227 + t240;
t216 = Icges(3,1) * t203 - t253;
t215 = -Icges(3,2) * t201 + t252;
t214 = Icges(3,5) * t203 - Icges(3,6) * t201;
t211 = t202 * t100 + t204 * t101 + t223;
t210 = rSges(3,1) * t243 - rSges(3,2) * t246 + t202 * rSges(3,3);
t208 = t209 + t233;
t207 = t34 / 0.2e1 + t32 / 0.2e1 + t36 / 0.2e1 + t49 / 0.2e1 + t46 / 0.2e1 + t45 / 0.2e1;
t206 = t35 / 0.2e1 + t33 / 0.2e1 + t37 / 0.2e1 + t50 / 0.2e1 + t48 / 0.2e1 + t47 / 0.2e1;
t205 = t101 + t208;
t195 = t201 ^ 2;
t173 = rSges(2,1) * t204 - t202 * rSges(2,2);
t172 = -t202 * rSges(2,1) - rSges(2,2) * t204;
t167 = Icges(3,5) * t201 + Icges(3,6) * t203;
t149 = Icges(3,3) * t202 + t204 * t214;
t148 = -Icges(3,3) * t204 + t202 * t214;
t142 = -Icges(4,5) * t203 + (Icges(4,1) * t199 - Icges(4,4) * t198) * t201;
t141 = -Icges(4,6) * t203 + (Icges(4,4) * t199 - Icges(4,2) * t198) * t201;
t117 = t210 + t233;
t116 = t258 + t192 + (-pkin(1) - t220) * t202;
t109 = t237 * t204;
t108 = t237 * t202;
t107 = Icges(4,1) * t163 + Icges(4,4) * t162 + Icges(4,5) * t246;
t106 = Icges(4,1) * t161 + Icges(4,4) * t160 + Icges(4,5) * t247;
t105 = Icges(4,4) * t163 + Icges(4,2) * t162 + Icges(4,6) * t246;
t104 = Icges(4,4) * t161 + Icges(4,2) * t160 + Icges(4,6) * t247;
t103 = Icges(4,5) * t163 + Icges(4,6) * t162 + Icges(4,3) * t246;
t102 = Icges(4,5) * t161 + Icges(4,6) * t160 + Icges(4,3) * t247;
t99 = t204 * t210 + (t202 * t220 - t258) * t202;
t93 = t100 * t246;
t90 = rSges(5,3) * t247 - t218;
t87 = rSges(6,1) * t247 - rSges(6,2) * t145 + t259;
t66 = t226 + t233 + t234;
t65 = t192 + (-t262 - pkin(1) + (-rSges(4,3) - qJ(3)) * t201) * t202 + t219;
t64 = t228 * t204;
t63 = t228 * t202;
t62 = -t129 * t246 - t203 * t91;
t61 = t129 * t247 + t203 * t90;
t60 = t208 + t91;
t59 = (-rSges(5,3) * t201 + t224) * t202 + t218 + t225;
t58 = t222 * t204;
t57 = t222 * t202;
t53 = t217 * t204;
t52 = t217 * t202;
t51 = (-t202 * t91 + t204 * t90) * t201;
t44 = t202 * (rSges(4,3) * t247 - t219) + t204 * t226 + t236;
t43 = t205 + t89;
t42 = -t259 + (rSges(6,2) - pkin(4)) * t145 + (-rSges(6,1) * t201 + t224) * t202 + t221;
t41 = t203 * t257 + t238 * t246;
t40 = t131 * t247 + t203 * t87 + t254;
t39 = t205 + t255;
t38 = -t260 + (-pkin(4) - t279) * t145 + (-t280 * t201 + t224) * t202 + t221;
t19 = t203 * t229 + t227 * t246;
t18 = t203 * t256 + t239 * t247 + t254;
t17 = t93 + (t202 * t257 + t204 * t87) * t201;
t16 = t202 * t90 + t204 * t91 + t223;
t15 = t93 + (t202 * t229 + t204 * t256) * t201;
t14 = t202 * t87 + t204 * t89 + t211;
t13 = t202 * t256 + t204 * t255 + t211;
t12 = t31 * t202 - t204 * t30;
t11 = t29 * t202 - t204 * t28;
t10 = t27 * t202 - t204 * t26;
t9 = t25 * t202 - t204 * t24;
t8 = t23 * t202 - t204 * t22;
t7 = -t20 * t204 + t21 * t202;
t1 = [Icges(2,3) + (t253 - (Icges(4,5) * t199 - Icges(4,6) * t198) * t201 + (Icges(3,2) + Icges(4,3)) * t203 + t278) * t203 + (Icges(3,1) * t201 - t121 * t188 - t126 * t189 - t141 * t198 + t142 * t199 + t252) * t201 + m(7) * (t38 ^ 2 + t39 ^ 2) + m(6) * (t42 ^ 2 + t43 ^ 2) + m(5) * (t59 ^ 2 + t60 ^ 2) + m(4) * (t65 ^ 2 + t66 ^ 2) + m(3) * (t116 ^ 2 + t117 ^ 2) + m(2) * (t172 ^ 2 + t173 ^ 2) + t273; m(4) * (t108 * t66 + t109 * t65) + m(6) * (t42 * t58 + t43 * t57) + m(5) * (t59 * t64 + t60 * t63) + m(7) * (t38 * t53 + t39 * t52) + (-t160 * t141 / 0.2e1 - t161 * t142 / 0.2e1 - t116 * t263 + t167 * t281 + (Icges(3,6) * t281 - t202 * t215 / 0.2e1 + t102 / 0.2e1) * t203 - t207) * t204 + (t162 * t141 / 0.2e1 + t163 * t142 / 0.2e1 - t117 * t263 + t167 * t266 + (Icges(3,6) * t266 + t215 * t281 - t103 / 0.2e1) * t203 + t206) * t202 + ((Icges(3,5) * t202 - t105 * t198 + t107 * t199 + t204 * t216) * t266 + (-Icges(3,5) * t204 - t104 * t198 + t106 * t199 + t202 * t216) * t264) * t201; m(7) * (t13 ^ 2 + t52 ^ 2 + t53 ^ 2) + m(6) * (t14 ^ 2 + t57 ^ 2 + t58 ^ 2) + m(5) * (t16 ^ 2 + t63 ^ 2 + t64 ^ 2) + m(4) * (t108 ^ 2 + t109 ^ 2 + t44 ^ 2) + m(3) * (t171 ^ 2 * t232 + t99 ^ 2) + (-t197 * t148 - t11 - t7 - t8 + (t102 * t247 + t160 * t104 + t161 * t106) * t204) * t204 + (t12 + t10 + t9 + t196 * t149 + (t103 * t246 + t162 * t105 + t163 * t107) * t202 + (-t102 * t246 - t103 * t247 - t162 * t104 - t105 * t160 - t163 * t106 - t107 * t161 - t202 * t148 + t204 * t149) * t204) * t202; ((t202 * t39 + t204 * t38) * t267 + (t202 * t43 + t204 * t42) * t268 + (t202 * t60 + t204 * t59) * t269 + (t202 * t66 + t204 * t65) * t270) * t271; m(7) * (-t203 * t13 + (t202 * t52 + t204 * t53) * t201) + m(6) * (-t203 * t14 + (t202 * t57 + t204 * t58) * t201) + m(5) * (-t203 * t16 + (t202 * t63 + t204 * t64) * t201) + m(4) * (-t203 * t44 + (t108 * t202 + t109 * t204) * t201); 0.2e1 * (t270 + t269 + t231) * (t195 * t232 + t272); -t230 * t203 + m(7) * (t18 * t38 + t19 * t39) + m(6) * (t40 * t42 + t41 * t43) + m(5) * (t59 * t61 + t60 * t62) + (t202 * t207 + t204 * t206) * t201; m(7) * (t13 * t15 + t18 * t53 + t19 * t52) + m(6) * (t14 * t17 + t40 * t58 + t41 * t57) + m(5) * (t16 * t51 + t61 * t64 + t62 * t63) + ((t12 / 0.2e1 + t10 / 0.2e1 + t9 / 0.2e1) * t204 + (t8 / 0.2e1 + t7 / 0.2e1 + t11 / 0.2e1) * t202) * t201 + t276 * t266 - (t202 * t274 + t204 * t275) * t203 / 0.2e1 + t277 * t264; m(5) * (-t51 * t203 + (t202 * t62 + t204 * t61) * t201) + m(6) * (-t17 * t203 + (t202 * t41 + t204 * t40) * t201) + m(7) * (-t15 * t203 + (t18 * t204 + t19 * t202) * t201); m(6) * (t17 ^ 2 + t40 ^ 2 + t41 ^ 2) + m(7) * (t15 ^ 2 + t18 ^ 2 + t19 ^ 2) + m(5) * (t51 ^ 2 + t61 ^ 2 + t62 ^ 2) + t230 * t272 + ((-t274 * t203 + t276) * t204 + (t275 * t203 + t277) * t202) * t201; m(7) * (t144 * t39 + t146 * t38) + m(6) * (t144 * t43 + t146 * t42); m(7) * (t13 * t250 + t144 * t52 + t146 * t53) + m(6) * (t14 * t250 + t144 * t57 + t146 * t58); t231 * (t144 * t202 + t146 * t204 - t188 * t203) * t271; m(6) * (t144 * t41 + t146 * t40 + t17 * t250) + m(7) * (t144 * t19 + t146 * t18 + t15 * t250); 0.2e1 * t231 * (t188 ^ 2 * t195 + t144 ^ 2 + t146 ^ 2); m(7) * (t145 * t39 + t147 * t38); m(7) * (t13 * t249 + t145 * t52 + t147 * t53); m(7) * (t145 * t202 + t147 * t204 - t189 * t203) * t201; m(7) * (t145 * t19 + t147 * t18 + t15 * t249); m(7) * (t188 * t189 * t195 + t144 * t145 + t146 * t147); m(7) * (t189 ^ 2 * t195 + t145 ^ 2 + t147 ^ 2);];
%% Postprocessing: Reshape Output
% From vec2symmat_6_matlab.m
res = [t1(1) t1(2) t1(4) t1(7) t1(11) t1(16); t1(2) t1(3) t1(5) t1(8) t1(12) t1(17); t1(4) t1(5) t1(6) t1(9) t1(13) t1(18); t1(7) t1(8) t1(9) t1(10) t1(14) t1(19); t1(11) t1(12) t1(13) t1(14) t1(15) t1(20); t1(16) t1(17) t1(18) t1(19) t1(20) t1(21);];
Mq  = res;
