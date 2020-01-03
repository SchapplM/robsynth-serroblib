% Calculate joint inertia matrix for
% S5RRRRP11
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
% Datum: 2019-12-31 22:20
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Mq = S5RRRRP11_inertiaJ_slag_vp1(qJ, ...
  pkin, m, rSges, Icges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(9,1),zeros(6,1),zeros(6,3),zeros(6,6)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RRRRP11_inertiaJ_slag_vp1: qJ has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S5RRRRP11_inertiaJ_slag_vp1: pkin has to be [9x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5RRRRP11_inertiaJ_slag_vp1: m has to be [6x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [6,3]), ...
  'S5RRRRP11_inertiaJ_slag_vp1: rSges has to be [6x3] (double)');
assert(isreal(Icges) && all(size(Icges) == [6 6]), ...
  'S5RRRRP11_inertiaJ_slag_vp1: Icges has to be [6x6] (double)'); 

%% Symbolic Calculation
% From inertia_joint_joint_floatb_twist_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 22:14:28
% EndTime: 2019-12-31 22:14:39
% DurationCPUTime: 4.16s
% Computational Cost: add. (13864->481), mult. (35604->691), div. (0->0), fcn. (45576->10), ass. (0->226)
t218 = cos(pkin(5));
t222 = sin(qJ(1));
t223 = cos(qJ(2));
t263 = t222 * t223;
t221 = sin(qJ(2));
t224 = cos(qJ(1));
t264 = t221 * t224;
t204 = t218 * t264 + t263;
t220 = sin(qJ(3));
t217 = sin(pkin(5));
t266 = t217 * t224;
t273 = cos(qJ(3));
t184 = t204 * t273 - t220 * t266;
t261 = t223 * t224;
t265 = t221 * t222;
t203 = -t218 * t261 + t265;
t219 = sin(qJ(4));
t272 = cos(qJ(4));
t153 = t184 * t219 - t203 * t272;
t154 = t184 * t272 + t203 * t219;
t238 = t217 * t273;
t183 = t204 * t220 + t224 * t238;
t274 = rSges(6,3) + qJ(5);
t275 = rSges(6,1) + pkin(4);
t270 = rSges(6,2) * t183 + t274 * t153 + t275 * t154;
t269 = t217 * t221;
t268 = t217 * t222;
t267 = t217 * t223;
t262 = t222 * t224;
t206 = -t218 * t265 + t261;
t186 = t206 * t273 + t220 * t268;
t205 = t218 * t263 + t264;
t155 = t186 * t219 - t205 * t272;
t156 = t186 * t272 + t205 * t219;
t185 = t206 * t220 - t222 * t238;
t260 = t185 * rSges(6,2) + t274 * t155 + t275 * t156;
t102 = t156 * rSges(5,1) - t155 * rSges(5,2) + t185 * rSges(5,3);
t142 = t186 * pkin(3) + pkin(9) * t185;
t259 = -t102 - t142;
t202 = t218 * t220 + t221 * t238;
t181 = t202 * t219 + t272 * t267;
t182 = t202 * t272 - t219 * t267;
t201 = -t218 * t273 + t220 * t269;
t258 = rSges(6,2) * t201 + t274 * t181 + t275 * t182;
t124 = rSges(5,1) * t182 - rSges(5,2) * t181 + rSges(5,3) * t201;
t172 = pkin(3) * t202 + pkin(9) * t201;
t257 = -t124 - t172;
t141 = pkin(3) * t184 + t183 * pkin(9);
t256 = t141 * t267 + t203 * t172;
t174 = t206 * pkin(2) + pkin(8) * t205;
t171 = t218 * t174;
t255 = t218 * t142 + t171;
t173 = pkin(2) * t204 + t203 * pkin(8);
t254 = -t141 - t173;
t253 = t173 * t268 + t174 * t266;
t252 = t224 * pkin(1) + pkin(7) * t268;
t87 = Icges(6,5) * t154 + Icges(6,6) * t183 + Icges(6,3) * t153;
t91 = Icges(6,4) * t154 + Icges(6,2) * t183 + Icges(6,6) * t153;
t95 = Icges(6,1) * t154 + Icges(6,4) * t183 + Icges(6,5) * t153;
t32 = t153 * t87 + t154 * t95 + t183 * t91;
t88 = Icges(6,5) * t156 + Icges(6,6) * t185 + Icges(6,3) * t155;
t92 = Icges(6,4) * t156 + Icges(6,2) * t185 + Icges(6,6) * t155;
t96 = Icges(6,1) * t156 + Icges(6,4) * t185 + Icges(6,5) * t155;
t33 = t153 * t88 + t154 * t96 + t183 * t92;
t117 = Icges(6,5) * t182 + Icges(6,6) * t201 + Icges(6,3) * t181;
t119 = Icges(6,4) * t182 + Icges(6,2) * t201 + Icges(6,6) * t181;
t121 = Icges(6,1) * t182 + Icges(6,4) * t201 + Icges(6,5) * t181;
t50 = t117 * t153 + t119 * t183 + t121 * t154;
t1 = t183 * t32 + t185 * t33 + t201 * t50;
t89 = Icges(5,5) * t154 - Icges(5,6) * t153 + Icges(5,3) * t183;
t93 = Icges(5,4) * t154 - Icges(5,2) * t153 + Icges(5,6) * t183;
t97 = Icges(5,1) * t154 - Icges(5,4) * t153 + Icges(5,5) * t183;
t34 = -t153 * t93 + t154 * t97 + t183 * t89;
t90 = Icges(5,5) * t156 - Icges(5,6) * t155 + Icges(5,3) * t185;
t94 = Icges(5,4) * t156 - Icges(5,2) * t155 + Icges(5,6) * t185;
t98 = Icges(5,1) * t156 - Icges(5,4) * t155 + Icges(5,5) * t185;
t35 = -t153 * t94 + t154 * t98 + t183 * t90;
t118 = Icges(5,5) * t182 - Icges(5,6) * t181 + Icges(5,3) * t201;
t120 = Icges(5,4) * t182 - Icges(5,2) * t181 + Icges(5,6) * t201;
t122 = Icges(5,1) * t182 - Icges(5,4) * t181 + Icges(5,5) * t201;
t51 = t118 * t183 - t120 * t153 + t122 * t154;
t2 = t183 * t34 + t185 * t35 + t201 * t51;
t251 = t2 / 0.2e1 + t1 / 0.2e1;
t36 = t155 * t87 + t156 * t95 + t185 * t91;
t37 = t155 * t88 + t156 * t96 + t185 * t92;
t52 = t117 * t155 + t119 * t185 + t121 * t156;
t3 = t183 * t36 + t185 * t37 + t201 * t52;
t38 = -t155 * t93 + t156 * t97 + t185 * t89;
t39 = -t155 * t94 + t156 * t98 + t185 * t90;
t53 = t118 * t185 - t120 * t155 + t122 * t156;
t4 = t183 * t38 + t185 * t39 + t201 * t53;
t250 = t4 / 0.2e1 + t3 / 0.2e1;
t5 = t203 * t32 + t205 * t33 - t50 * t267;
t6 = t203 * t34 + t205 * t35 - t51 * t267;
t249 = t6 / 0.2e1 + t5 / 0.2e1;
t7 = t203 * t36 + t205 * t37 - t52 * t267;
t8 = t203 * t38 + t205 * t39 - t53 * t267;
t248 = t8 / 0.2e1 + t7 / 0.2e1;
t63 = t181 * t117 + t201 * t119 + t182 * t121;
t64 = t201 * t118 - t181 * t120 + t182 * t122;
t157 = Icges(4,5) * t202 - Icges(4,6) * t201 - Icges(4,3) * t267;
t158 = Icges(4,4) * t202 - Icges(4,2) * t201 - Icges(4,6) * t267;
t159 = Icges(4,1) * t202 - Icges(4,4) * t201 - Icges(4,5) * t267;
t82 = -t157 * t267 - t201 * t158 + t202 * t159;
t247 = -t63 - t64 - t82;
t10 = t51 * t218 + (t222 * t35 - t224 * t34) * t217;
t9 = t50 * t218 + (t222 * t33 - t224 * t32) * t217;
t246 = t10 / 0.2e1 + t9 / 0.2e1;
t11 = t52 * t218 + (t222 * t37 - t224 * t36) * t217;
t12 = t53 * t218 + (t222 * t39 - t224 * t38) * t217;
t245 = t12 / 0.2e1 + t11 / 0.2e1;
t41 = t181 * t87 + t182 * t95 + t201 * t91;
t42 = t181 * t88 + t182 * t96 + t201 * t92;
t59 = t63 * t201;
t13 = t41 * t183 + t42 * t185 + t59;
t43 = -t181 * t93 + t182 * t97 + t201 * t89;
t44 = -t181 * t94 + t182 * t98 + t201 * t90;
t60 = t64 * t201;
t14 = t43 * t183 + t44 * t185 + t60;
t244 = t14 / 0.2e1 + t13 / 0.2e1;
t15 = t41 * t203 + t42 * t205 - t63 * t267;
t16 = t43 * t203 + t44 * t205 - t64 * t267;
t243 = t16 / 0.2e1 + t15 / 0.2e1;
t61 = t63 * t218;
t17 = t61 + (t42 * t222 - t41 * t224) * t217;
t62 = t64 * t218;
t18 = t62 + (t44 * t222 - t43 * t224) * t217;
t242 = t18 / 0.2e1 + t17 / 0.2e1;
t241 = -t142 - t260;
t240 = -t172 - t258;
t132 = t186 * rSges(4,1) - t185 * rSges(4,2) + t205 * rSges(4,3);
t190 = Icges(3,3) * t218 + (Icges(3,5) * t221 + Icges(3,6) * t223) * t217;
t191 = Icges(3,6) * t218 + (Icges(3,4) * t221 + Icges(3,2) * t223) * t217;
t192 = Icges(3,5) * t218 + (Icges(3,1) * t221 + Icges(3,4) * t223) * t217;
t239 = t218 * t190 + t191 * t267 + t192 * t269;
t168 = t206 * rSges(3,1) - t205 * rSges(3,2) + rSges(3,3) * t268;
t237 = -t222 * pkin(1) + pkin(7) * t266;
t160 = rSges(4,1) * t202 - rSges(4,2) * t201 - rSges(4,3) * t267;
t207 = (pkin(2) * t221 - pkin(8) * t223) * t217;
t236 = t217 * (-t160 - t207);
t235 = t141 * t268 + t142 * t266 + t253;
t234 = t217 * (-t207 + t257);
t233 = t174 + t252;
t232 = t217 * (-t207 + t240);
t231 = t41 / 0.2e1 + t43 / 0.2e1 + t51 / 0.2e1 + t50 / 0.2e1;
t230 = t52 / 0.2e1 + t42 / 0.2e1 + t44 / 0.2e1 + t53 / 0.2e1;
t229 = -t173 + t237;
t131 = rSges(4,1) * t184 - rSges(4,2) * t183 + rSges(4,3) * t203;
t100 = rSges(5,1) * t154 - rSges(5,2) * t153 + rSges(5,3) * t183;
t167 = t204 * rSges(3,1) - t203 * rSges(3,2) - rSges(3,3) * t266;
t125 = Icges(4,5) * t184 - Icges(4,6) * t183 + Icges(4,3) * t203;
t127 = Icges(4,4) * t184 - Icges(4,2) * t183 + Icges(4,6) * t203;
t129 = Icges(4,1) * t184 - Icges(4,4) * t183 + Icges(4,5) * t203;
t69 = -t125 * t267 - t127 * t201 + t129 * t202;
t76 = t157 * t203 - t158 * t183 + t159 * t184;
t228 = t69 / 0.2e1 + t76 / 0.2e1 + t231;
t126 = Icges(4,5) * t186 - Icges(4,6) * t185 + Icges(4,3) * t205;
t128 = Icges(4,4) * t186 - Icges(4,2) * t185 + Icges(4,6) * t205;
t130 = Icges(4,1) * t186 - Icges(4,4) * t185 + Icges(4,5) * t205;
t70 = -t126 * t267 - t128 * t201 + t130 * t202;
t77 = t157 * t205 - t158 * t185 + t159 * t186;
t227 = t70 / 0.2e1 + t77 / 0.2e1 + t230;
t226 = t142 + t233;
t225 = -t141 + t229;
t209 = rSges(2,1) * t224 - t222 * rSges(2,2);
t208 = -t222 * rSges(2,1) - rSges(2,2) * t224;
t193 = rSges(3,3) * t218 + (rSges(3,1) * t221 + rSges(3,2) * t223) * t217;
t166 = Icges(3,1) * t206 - Icges(3,4) * t205 + Icges(3,5) * t268;
t165 = Icges(3,1) * t204 - Icges(3,4) * t203 - Icges(3,5) * t266;
t164 = Icges(3,4) * t206 - Icges(3,2) * t205 + Icges(3,6) * t268;
t163 = Icges(3,4) * t204 - Icges(3,2) * t203 - Icges(3,6) * t266;
t162 = Icges(3,5) * t206 - Icges(3,6) * t205 + Icges(3,3) * t268;
t161 = Icges(3,5) * t204 - Icges(3,6) * t203 - Icges(3,3) * t266;
t146 = t168 + t252;
t145 = -t167 + t237;
t135 = -t218 * t167 - t193 * t266;
t134 = t168 * t218 - t193 * t268;
t133 = t205 * t141;
t116 = t239 * t218;
t113 = (t167 * t222 + t168 * t224) * t217;
t112 = t190 * t268 - t191 * t205 + t192 * t206;
t111 = -t190 * t266 - t203 * t191 + t204 * t192;
t104 = t233 + t132;
t103 = -t131 + t229;
t86 = -t132 * t267 - t160 * t205;
t85 = t131 * t267 + t160 * t203;
t84 = t162 * t218 + (t164 * t223 + t166 * t221) * t217;
t83 = t161 * t218 + (t163 * t223 + t165 * t221) * t217;
t81 = t82 * t218;
t80 = t131 * t205 - t132 * t203;
t79 = (-t131 - t173) * t218 + t224 * t236;
t78 = t132 * t218 + t222 * t236 + t171;
t75 = t226 + t102;
t74 = -t100 + t225;
t73 = (t131 * t222 + t132 * t224) * t217 + t253;
t72 = t102 * t201 - t124 * t185;
t71 = -t100 * t201 + t124 * t183;
t68 = t126 * t205 - t128 * t185 + t130 * t186;
t67 = t125 * t205 - t127 * t185 + t129 * t186;
t66 = t126 * t203 - t128 * t183 + t130 * t184;
t65 = t125 * t203 - t127 * t183 + t129 * t184;
t58 = t100 * t185 - t102 * t183;
t57 = t226 + t260;
t56 = t225 - t270;
t55 = t257 * t205 + t259 * t267;
t54 = t100 * t267 + t124 * t203 + t256;
t49 = (-t100 + t254) * t218 + t224 * t234;
t48 = t102 * t218 + t222 * t234 + t255;
t47 = t100 * t205 + t259 * t203 + t133;
t46 = -t258 * t185 + t260 * t201;
t45 = t258 * t183 - t270 * t201;
t40 = (t100 * t222 + t102 * t224) * t217 + t235;
t31 = t240 * t205 + t241 * t267;
t30 = t258 * t203 + t270 * t267 + t256;
t29 = (t254 - t270) * t218 + t224 * t232;
t28 = t260 * t218 + t222 * t232 + t255;
t27 = -t260 * t183 + t270 * t185;
t26 = t241 * t203 + t270 * t205 + t133;
t25 = (t270 * t222 + t260 * t224) * t217 + t235;
t24 = t81 + (t70 * t222 - t69 * t224) * t217;
t23 = t69 * t203 + t70 * t205 - t82 * t267;
t22 = t77 * t218 + (t222 * t68 - t224 * t67) * t217;
t21 = t76 * t218 + (t222 * t66 - t224 * t65) * t217;
t20 = t203 * t67 + t205 * t68 - t77 * t267;
t19 = t203 * t65 + t205 * t66 - t76 * t267;
t99 = [Icges(2,3) + m(6) * (t56 ^ 2 + t57 ^ 2) + m(5) * (t74 ^ 2 + t75 ^ 2) + m(4) * (t103 ^ 2 + t104 ^ 2) + m(3) * (t145 ^ 2 + t146 ^ 2) + m(2) * (t208 ^ 2 + t209 ^ 2) + t239 - t247; t62 + t61 + t81 + t116 + m(6) * (t28 * t57 + t29 * t56) + m(5) * (t48 * t75 + t49 * t74) + m(4) * (t103 * t79 + t104 * t78) + m(3) * (t134 * t146 + t135 * t145) + ((-t83 / 0.2e1 - t111 / 0.2e1 - t228) * t224 + (t84 / 0.2e1 + t112 / 0.2e1 + t227) * t222) * t217; (t18 + t17 + t24 + t116) * t218 + m(6) * (t25 ^ 2 + t28 ^ 2 + t29 ^ 2) + m(5) * (t40 ^ 2 + t48 ^ 2 + t49 ^ 2) + m(4) * (t73 ^ 2 + t78 ^ 2 + t79 ^ 2) + m(3) * (t113 ^ 2 + t134 ^ 2 + t135 ^ 2) + (t222 * t12 + t222 * t11 - t224 * t10 - t224 * t9 - t224 * t21 + t222 * t22 + (t222 * ((-t164 * t205 + t166 * t206) * t222 - (-t163 * t205 + t165 * t206) * t224) - t224 * ((-t203 * t164 + t204 * t166) * t222 - (-t203 * t163 + t204 * t165) * t224) + (t222 * (t162 * t222 ^ 2 - t161 * t262) - t224 * (t161 * t224 ^ 2 - t162 * t262)) * t217) * t217 + ((-t111 - t83) * t224 + (t112 + t84) * t222) * t218) * t217; t247 * t267 + m(6) * (t30 * t56 + t31 * t57) + m(5) * (t54 * t74 + t55 * t75) + m(4) * (t103 * t85 + t104 * t86) + t227 * t205 + t228 * t203; (t23 / 0.2e1 + t243) * t218 + (t22 / 0.2e1 + t245) * t205 + (t21 / 0.2e1 + t246) * t203 + m(6) * (t25 * t26 + t28 * t31 + t29 * t30) + m(5) * (t40 * t47 + t48 * t55 + t49 * t54) + m(4) * (t73 * t80 + t78 * t86 + t79 * t85) + ((-t19 / 0.2e1 - t249) * t224 + (-t24 / 0.2e1 - t242) * t223 + (t20 / 0.2e1 + t248) * t222) * t217; (-t15 - t16 - t23) * t267 + (t8 + t7 + t20) * t205 + (t6 + t5 + t19) * t203 + m(6) * (t26 ^ 2 + t30 ^ 2 + t31 ^ 2) + m(5) * (t47 ^ 2 + t54 ^ 2 + t55 ^ 2) + m(4) * (t80 ^ 2 + t85 ^ 2 + t86 ^ 2); t60 + t59 + m(6) * (t45 * t56 + t46 * t57) + m(5) * (t71 * t74 + t72 * t75) + t230 * t185 + t231 * t183; t244 * t218 + t242 * t201 + t245 * t185 + t246 * t183 + m(6) * (t25 * t27 + t28 * t46 + t29 * t45) + m(5) * (t58 * t40 + t48 * t72 + t49 * t71) + (t250 * t222 - t251 * t224) * t217; -t244 * t267 + t250 * t205 + t251 * t203 + t243 * t201 + t248 * t185 + t249 * t183 + m(6) * (t26 * t27 + t30 * t45 + t31 * t46) + m(5) * (t58 * t47 + t54 * t71 + t55 * t72); (t13 + t14) * t201 + (t4 + t3) * t185 + (t1 + t2) * t183 + m(6) * (t27 ^ 2 + t45 ^ 2 + t46 ^ 2) + m(5) * (t58 ^ 2 + t71 ^ 2 + t72 ^ 2); m(6) * (t153 * t57 + t155 * t56); m(6) * (t153 * t28 + t155 * t29 + t181 * t25); m(6) * (t153 * t31 + t155 * t30 + t181 * t26); m(6) * (t153 * t46 + t155 * t45 + t181 * t27); m(6) * (t153 ^ 2 + t155 ^ 2 + t181 ^ 2);];
%% Postprocessing: Reshape Output
% From vec2symmat_5_matlab.m
res = [t99(1), t99(2), t99(4), t99(7), t99(11); t99(2), t99(3), t99(5), t99(8), t99(12); t99(4), t99(5), t99(6), t99(9), t99(13); t99(7), t99(8), t99(9), t99(10), t99(14); t99(11), t99(12), t99(13), t99(14), t99(15);];
Mq = res;
