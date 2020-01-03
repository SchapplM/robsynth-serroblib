% Calculate joint inertia matrix for
% S5RRRRP10
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% pkin [9x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,alpha2,d1,d2,d3,d4]';
% m_mdh [6x1]
%   mass of all robot links (including the base)
% rSges [6x3]
%   center of mass of all robot links (in body frames)
%   rows: links of the robot (starting with base)
%   columns: x-, y-, z-coordinates
% Icges [6x6]
%   inertia of all robot links about their respective center of mass, in body frames
%   rows: links of the robot (starting with base)
%   columns: xx, yy, zz, xy, xz, yz (see inertiavector2matrix.m)
% 
% Output:
% Mq [5x5]
%   inertia matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-31 22:13
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Mq = S5RRRRP10_inertiaJ_slag_vp1(qJ, ...
  pkin, m, rSges, Icges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(9,1),zeros(6,1),zeros(6,3),zeros(6,6)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RRRRP10_inertiaJ_slag_vp1: qJ has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S5RRRRP10_inertiaJ_slag_vp1: pkin has to be [9x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5RRRRP10_inertiaJ_slag_vp1: m has to be [6x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [6,3]), ...
  'S5RRRRP10_inertiaJ_slag_vp1: rSges has to be [6x3] (double)');
assert(isreal(Icges) && all(size(Icges) == [6 6]), ...
  'S5RRRRP10_inertiaJ_slag_vp1: Icges has to be [6x6] (double)'); 

%% Symbolic Calculation
% From inertia_joint_joint_floatb_twist_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 22:08:11
% EndTime: 2019-12-31 22:08:23
% DurationCPUTime: 4.17s
% Computational Cost: add. (14254->490), mult. (36082->699), div. (0->0), fcn. (46014->10), ass. (0->231)
t282 = rSges(6,3) + qJ(5) + pkin(9);
t219 = cos(pkin(5));
t224 = sin(qJ(1));
t226 = cos(qJ(2));
t267 = t224 * t226;
t223 = sin(qJ(2));
t227 = cos(qJ(1));
t268 = t223 * t227;
t204 = t219 * t268 + t267;
t222 = sin(qJ(3));
t218 = sin(pkin(5));
t270 = t218 * t227;
t279 = cos(qJ(3));
t185 = t204 * t279 - t222 * t270;
t265 = t226 * t227;
t269 = t223 * t224;
t203 = -t219 * t265 + t269;
t221 = sin(qJ(4));
t225 = cos(qJ(4));
t152 = -t185 * t221 + t203 * t225;
t275 = t203 * t221;
t153 = t185 * t225 + t275;
t240 = t218 * t279;
t184 = t204 * t222 + t227 * t240;
t283 = rSges(6,1) * t153 + rSges(6,2) * t152 + pkin(4) * t275 + t282 * t184;
t206 = -t219 * t269 + t265;
t272 = t218 * t224;
t187 = t206 * t279 + t222 * t272;
t205 = t219 * t267 + t268;
t154 = -t187 * t221 + t205 * t225;
t274 = t205 * t221;
t155 = t187 * t225 + t274;
t186 = t206 * t222 - t224 * t240;
t216 = pkin(4) * t225 + pkin(3);
t281 = t155 * rSges(6,1) + t154 * rSges(6,2) + pkin(4) * t274 + t282 * t186 + t187 * t216;
t278 = -pkin(3) + t216;
t180 = t184 * pkin(9);
t277 = t185 * t278 - t180 + t283;
t142 = t187 * pkin(3) + pkin(9) * t186;
t276 = -t142 + t281;
t273 = t218 * t223;
t271 = t218 * t226;
t266 = t224 * t227;
t104 = t155 * rSges(5,1) + t154 * rSges(5,2) + t186 * rSges(5,3);
t264 = -t104 - t142;
t202 = t219 * t222 + t223 * t240;
t182 = -t202 * t221 - t225 * t271;
t244 = t221 * t271;
t183 = t202 * t225 - t244;
t201 = -t219 * t279 + t222 * t273;
t263 = rSges(6,1) * t183 + rSges(6,2) * t182 - pkin(4) * t244 + t278 * t202 + (-pkin(9) + t282) * t201;
t125 = rSges(5,1) * t183 + rSges(5,2) * t182 + rSges(5,3) * t201;
t171 = pkin(3) * t202 + pkin(9) * t201;
t262 = -t125 - t171;
t141 = pkin(3) * t185 + t180;
t261 = t141 * t271 + t203 * t171;
t173 = t206 * pkin(2) + pkin(8) * t205;
t170 = t219 * t173;
t260 = t219 * t142 + t170;
t172 = pkin(2) * t204 + t203 * pkin(8);
t259 = -t141 - t172;
t258 = t172 * t272 + t173 * t270;
t257 = t227 * pkin(1) + pkin(7) * t272;
t89 = Icges(6,5) * t153 + Icges(6,6) * t152 + Icges(6,3) * t184;
t93 = Icges(6,4) * t153 + Icges(6,2) * t152 + Icges(6,6) * t184;
t97 = Icges(6,1) * t153 + Icges(6,4) * t152 + Icges(6,5) * t184;
t32 = t152 * t93 + t153 * t97 + t184 * t89;
t90 = Icges(6,5) * t155 + Icges(6,6) * t154 + Icges(6,3) * t186;
t94 = Icges(6,4) * t155 + Icges(6,2) * t154 + Icges(6,6) * t186;
t98 = Icges(6,1) * t155 + Icges(6,4) * t154 + Icges(6,5) * t186;
t33 = t152 * t94 + t153 * t98 + t184 * t90;
t118 = Icges(6,5) * t183 + Icges(6,6) * t182 + Icges(6,3) * t201;
t120 = Icges(6,4) * t183 + Icges(6,2) * t182 + Icges(6,6) * t201;
t122 = Icges(6,1) * t183 + Icges(6,4) * t182 + Icges(6,5) * t201;
t50 = t118 * t184 + t120 * t152 + t122 * t153;
t1 = t184 * t32 + t186 * t33 + t201 * t50;
t91 = Icges(5,5) * t153 + Icges(5,6) * t152 + Icges(5,3) * t184;
t95 = Icges(5,4) * t153 + Icges(5,2) * t152 + Icges(5,6) * t184;
t99 = Icges(5,1) * t153 + Icges(5,4) * t152 + Icges(5,5) * t184;
t34 = t152 * t95 + t153 * t99 + t184 * t91;
t100 = Icges(5,1) * t155 + Icges(5,4) * t154 + Icges(5,5) * t186;
t92 = Icges(5,5) * t155 + Icges(5,6) * t154 + Icges(5,3) * t186;
t96 = Icges(5,4) * t155 + Icges(5,2) * t154 + Icges(5,6) * t186;
t35 = t100 * t153 + t152 * t96 + t184 * t92;
t119 = Icges(5,5) * t183 + Icges(5,6) * t182 + Icges(5,3) * t201;
t121 = Icges(5,4) * t183 + Icges(5,2) * t182 + Icges(5,6) * t201;
t123 = Icges(5,1) * t183 + Icges(5,4) * t182 + Icges(5,5) * t201;
t51 = t119 * t184 + t121 * t152 + t123 * t153;
t2 = t184 * t34 + t186 * t35 + t201 * t51;
t256 = t1 / 0.2e1 + t2 / 0.2e1;
t36 = t154 * t93 + t155 * t97 + t186 * t89;
t37 = t154 * t94 + t155 * t98 + t186 * t90;
t52 = t118 * t186 + t120 * t154 + t122 * t155;
t3 = t184 * t36 + t186 * t37 + t201 * t52;
t38 = t154 * t95 + t155 * t99 + t186 * t91;
t39 = t100 * t155 + t154 * t96 + t186 * t92;
t53 = t119 * t186 + t121 * t154 + t123 * t155;
t4 = t184 * t38 + t186 * t39 + t201 * t53;
t255 = t4 / 0.2e1 + t3 / 0.2e1;
t5 = t203 * t32 + t205 * t33 - t271 * t50;
t6 = t203 * t34 + t205 * t35 - t271 * t51;
t254 = t5 / 0.2e1 + t6 / 0.2e1;
t7 = t203 * t36 + t205 * t37 - t271 * t52;
t8 = t203 * t38 + t205 * t39 - t271 * t53;
t253 = t8 / 0.2e1 + t7 / 0.2e1;
t61 = t201 * t118 + t182 * t120 + t183 * t122;
t62 = t201 * t119 + t182 * t121 + t183 * t123;
t156 = Icges(4,5) * t202 - Icges(4,6) * t201 - Icges(4,3) * t271;
t157 = Icges(4,4) * t202 - Icges(4,2) * t201 - Icges(4,6) * t271;
t158 = Icges(4,1) * t202 - Icges(4,4) * t201 - Icges(4,5) * t271;
t82 = -t156 * t271 - t201 * t157 + t202 * t158;
t252 = -t61 - t62 - t82;
t10 = t51 * t219 + (t224 * t35 - t227 * t34) * t218;
t9 = t50 * t219 + (t224 * t33 - t227 * t32) * t218;
t250 = t9 / 0.2e1 + t10 / 0.2e1;
t11 = t52 * t219 + (t224 * t37 - t227 * t36) * t218;
t12 = t53 * t219 + (t224 * t39 - t227 * t38) * t218;
t249 = t12 / 0.2e1 + t11 / 0.2e1;
t43 = t182 * t93 + t183 * t97 + t201 * t89;
t44 = t182 * t94 + t183 * t98 + t201 * t90;
t57 = t61 * t201;
t13 = t43 * t184 + t44 * t186 + t57;
t45 = t182 * t95 + t183 * t99 + t201 * t91;
t46 = t100 * t183 + t182 * t96 + t201 * t92;
t58 = t62 * t201;
t14 = t45 * t184 + t46 * t186 + t58;
t248 = t13 / 0.2e1 + t14 / 0.2e1;
t15 = t43 * t203 + t44 * t205 - t271 * t61;
t16 = t45 * t203 + t46 * t205 - t271 * t62;
t247 = t15 / 0.2e1 + t16 / 0.2e1;
t59 = t61 * t219;
t17 = t59 + (t44 * t224 - t43 * t227) * t218;
t60 = t62 * t219;
t18 = t60 + (t46 * t224 - t45 * t227) * t218;
t246 = t18 / 0.2e1 + t17 / 0.2e1;
t245 = -t142 - t276;
t243 = -t171 - t263;
t133 = t187 * rSges(4,1) - t186 * rSges(4,2) + t205 * rSges(4,3);
t191 = Icges(3,3) * t219 + (Icges(3,5) * t223 + Icges(3,6) * t226) * t218;
t192 = Icges(3,6) * t219 + (Icges(3,4) * t223 + Icges(3,2) * t226) * t218;
t193 = Icges(3,5) * t219 + (Icges(3,1) * t223 + Icges(3,4) * t226) * t218;
t241 = t219 * t191 + t192 * t271 + t193 * t273;
t167 = t206 * rSges(3,1) - t205 * rSges(3,2) + rSges(3,3) * t272;
t239 = -t224 * pkin(1) + pkin(7) * t270;
t159 = rSges(4,1) * t202 - rSges(4,2) * t201 - rSges(4,3) * t271;
t207 = (pkin(2) * t223 - pkin(8) * t226) * t218;
t238 = t218 * (-t159 - t207);
t237 = t141 * t272 + t142 * t270 + t258;
t236 = t218 * (-t207 + t262);
t234 = t173 + t257;
t233 = t218 * (-t207 + t243);
t232 = t45 / 0.2e1 + t43 / 0.2e1 + t51 / 0.2e1 + t50 / 0.2e1;
t231 = t53 / 0.2e1 + t52 / 0.2e1 + t46 / 0.2e1 + t44 / 0.2e1;
t230 = -t172 + t239;
t132 = rSges(4,1) * t185 - rSges(4,2) * t184 + rSges(4,3) * t203;
t102 = rSges(5,1) * t153 + rSges(5,2) * t152 + rSges(5,3) * t184;
t166 = t204 * rSges(3,1) - t203 * rSges(3,2) - rSges(3,3) * t270;
t126 = Icges(4,5) * t185 - Icges(4,6) * t184 + Icges(4,3) * t203;
t128 = Icges(4,4) * t185 - Icges(4,2) * t184 + Icges(4,6) * t203;
t130 = Icges(4,1) * t185 - Icges(4,4) * t184 + Icges(4,5) * t203;
t67 = -t126 * t271 - t128 * t201 + t130 * t202;
t76 = t156 * t203 - t157 * t184 + t158 * t185;
t229 = t67 / 0.2e1 + t76 / 0.2e1 + t232;
t127 = Icges(4,5) * t187 - Icges(4,6) * t186 + Icges(4,3) * t205;
t129 = Icges(4,4) * t187 - Icges(4,2) * t186 + Icges(4,6) * t205;
t131 = Icges(4,1) * t187 - Icges(4,4) * t186 + Icges(4,5) * t205;
t68 = -t127 * t271 - t129 * t201 + t131 * t202;
t77 = t156 * t205 - t157 * t186 + t158 * t187;
t228 = t68 / 0.2e1 + t77 / 0.2e1 + t231;
t209 = rSges(2,1) * t227 - t224 * rSges(2,2);
t208 = -t224 * rSges(2,1) - rSges(2,2) * t227;
t194 = rSges(3,3) * t219 + (rSges(3,1) * t223 + rSges(3,2) * t226) * t218;
t165 = Icges(3,1) * t206 - Icges(3,4) * t205 + Icges(3,5) * t272;
t164 = Icges(3,1) * t204 - Icges(3,4) * t203 - Icges(3,5) * t270;
t163 = Icges(3,4) * t206 - Icges(3,2) * t205 + Icges(3,6) * t272;
t162 = Icges(3,4) * t204 - Icges(3,2) * t203 - Icges(3,6) * t270;
t161 = Icges(3,5) * t206 - Icges(3,6) * t205 + Icges(3,3) * t272;
t160 = Icges(3,5) * t204 - Icges(3,6) * t203 - Icges(3,3) * t270;
t146 = t167 + t257;
t145 = -t166 + t239;
t136 = -t219 * t166 - t194 * t270;
t135 = t167 * t219 - t194 * t272;
t134 = t205 * t141;
t117 = t241 * t219;
t113 = (t166 * t224 + t167 * t227) * t218;
t112 = t191 * t272 - t192 * t205 + t193 * t206;
t111 = -t191 * t270 - t203 * t192 + t204 * t193;
t106 = t234 + t133;
t105 = -t132 + t230;
t88 = -t133 * t271 - t159 * t205;
t87 = t132 * t271 + t159 * t203;
t86 = t161 * t219 + (t163 * t226 + t165 * t223) * t218;
t85 = t160 * t219 + (t162 * t226 + t164 * t223) * t218;
t81 = t82 * t219;
t80 = t132 * t205 - t133 * t203;
t79 = (-t132 - t172) * t219 + t227 * t238;
t78 = t133 * t219 + t224 * t238 + t170;
t75 = t234 - t264;
t74 = -t102 - t141 + t230;
t73 = (t132 * t224 + t133 * t227) * t218 + t258;
t72 = t104 * t201 - t125 * t186;
t71 = -t102 * t201 + t125 * t184;
t70 = t234 + t281;
t69 = -t185 * t216 + t230 - t283;
t66 = t127 * t205 - t129 * t186 + t131 * t187;
t65 = t126 * t205 - t128 * t186 + t130 * t187;
t64 = t127 * t203 - t129 * t184 + t131 * t185;
t63 = t126 * t203 - t128 * t184 + t130 * t185;
t56 = t102 * t186 - t104 * t184;
t55 = t205 * t262 + t264 * t271;
t54 = t102 * t271 + t125 * t203 + t261;
t49 = (-t102 + t259) * t219 + t227 * t236;
t48 = t104 * t219 + t224 * t236 + t260;
t47 = t102 * t205 + t203 * t264 + t134;
t42 = (t102 * t224 + t104 * t227) * t218 + t237;
t41 = -t186 * t263 + t201 * t276;
t40 = t184 * t263 - t201 * t277;
t31 = t205 * t243 + t245 * t271;
t30 = t203 * t263 + t271 * t277 + t261;
t29 = (t259 - t277) * t219 + t227 * t233;
t28 = t219 * t276 + t224 * t233 + t260;
t27 = -t184 * t276 + t186 * t277;
t26 = t81 + (t68 * t224 - t67 * t227) * t218;
t25 = t67 * t203 + t68 * t205 - t271 * t82;
t24 = t203 * t245 + t205 * t277 + t134;
t23 = (t224 * t277 + t227 * t276) * t218 + t237;
t22 = t77 * t219 + (t224 * t66 - t227 * t65) * t218;
t21 = t76 * t219 + (t224 * t64 - t227 * t63) * t218;
t20 = t203 * t65 + t205 * t66 - t271 * t77;
t19 = t203 * t63 + t205 * t64 - t271 * t76;
t83 = [Icges(2,3) + m(6) * (t69 ^ 2 + t70 ^ 2) + m(5) * (t74 ^ 2 + t75 ^ 2) + m(4) * (t105 ^ 2 + t106 ^ 2) + m(3) * (t145 ^ 2 + t146 ^ 2) + m(2) * (t208 ^ 2 + t209 ^ 2) + t241 - t252; t60 + t59 + t81 + t117 + m(6) * (t28 * t70 + t29 * t69) + m(5) * (t48 * t75 + t49 * t74) + m(4) * (t105 * t79 + t106 * t78) + m(3) * (t135 * t146 + t136 * t145) + ((-t111 / 0.2e1 - t85 / 0.2e1 - t229) * t227 + (t112 / 0.2e1 + t86 / 0.2e1 + t228) * t224) * t218; (t17 + t18 + t26 + t117) * t219 + m(6) * (t23 ^ 2 + t28 ^ 2 + t29 ^ 2) + m(5) * (t42 ^ 2 + t48 ^ 2 + t49 ^ 2) + m(4) * (t73 ^ 2 + t78 ^ 2 + t79 ^ 2) + m(3) * (t113 ^ 2 + t135 ^ 2 + t136 ^ 2) + (-t227 * t9 + t224 * t12 + t224 * t11 - t227 * t10 + t224 * t22 - t227 * t21 + (t224 * ((-t163 * t205 + t165 * t206) * t224 - (-t162 * t205 + t164 * t206) * t227) - t227 * ((-t203 * t163 + t204 * t165) * t224 - (-t203 * t162 + t204 * t164) * t227) + (t224 * (t161 * t224 ^ 2 - t160 * t266) - t227 * (t160 * t227 ^ 2 - t161 * t266)) * t218) * t218 + ((-t111 - t85) * t227 + (t112 + t86) * t224) * t219) * t218; t252 * t271 + m(6) * (t30 * t69 + t31 * t70) + m(5) * (t54 * t74 + t55 * t75) + m(4) * (t105 * t87 + t106 * t88) + t228 * t205 + t229 * t203; (t25 / 0.2e1 + t247) * t219 + (t22 / 0.2e1 + t249) * t205 + (t21 / 0.2e1 + t250) * t203 + m(6) * (t23 * t24 + t28 * t31 + t29 * t30) + m(5) * (t42 * t47 + t48 * t55 + t49 * t54) + m(4) * (t73 * t80 + t78 * t88 + t79 * t87) + ((-t19 / 0.2e1 - t254) * t227 + (-t26 / 0.2e1 - t246) * t226 + (t20 / 0.2e1 + t253) * t224) * t218; (-t15 - t16 - t25) * t271 + (t8 + t7 + t20) * t205 + (t5 + t6 + t19) * t203 + m(6) * (t24 ^ 2 + t30 ^ 2 + t31 ^ 2) + m(5) * (t47 ^ 2 + t54 ^ 2 + t55 ^ 2) + m(4) * (t80 ^ 2 + t87 ^ 2 + t88 ^ 2); t58 + t57 + m(6) * (t40 * t69 + t41 * t70) + m(5) * (t71 * t74 + t72 * t75) + t231 * t186 + t232 * t184; t248 * t219 + t246 * t201 + t249 * t186 + t250 * t184 + m(6) * (t23 * t27 + t28 * t41 + t29 * t40) + m(5) * (t42 * t56 + t48 * t72 + t49 * t71) + (t224 * t255 - t227 * t256) * t218; -t248 * t271 + t255 * t205 + t256 * t203 + t247 * t201 + t253 * t186 + t254 * t184 + m(6) * (t24 * t27 + t30 * t40 + t31 * t41) + m(5) * (t47 * t56 + t54 * t71 + t55 * t72); (t13 + t14) * t201 + (t3 + t4) * t186 + (t1 + t2) * t184 + m(6) * (t27 ^ 2 + t40 ^ 2 + t41 ^ 2) + m(5) * (t56 ^ 2 + t71 ^ 2 + t72 ^ 2); m(6) * (t184 * t70 + t186 * t69); m(6) * (t184 * t28 + t186 * t29 + t201 * t23); m(6) * (t184 * t31 + t186 * t30 + t201 * t24); m(6) * (t184 * t41 + t186 * t40 + t201 * t27); m(6) * (t184 ^ 2 + t186 ^ 2 + t201 ^ 2);];
%% Postprocessing: Reshape Output
% From vec2symmat_5_matlab.m
res = [t83(1), t83(2), t83(4), t83(7), t83(11); t83(2), t83(3), t83(5), t83(8), t83(12); t83(4), t83(5), t83(6), t83(9), t83(13); t83(7), t83(8), t83(9), t83(10), t83(14); t83(11), t83(12), t83(13), t83(14), t83(15);];
Mq = res;
