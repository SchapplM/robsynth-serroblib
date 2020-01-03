% Calculate joint inertia matrix for
% S5RRPRR10
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% pkin [10x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,alpha2,d1,d2,d4,d5,theta3]';
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
% Datum: 2019-12-31 20:28
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Mq = S5RRPRR10_inertiaJ_slag_vp1(qJ, ...
  pkin, m, rSges, Icges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(10,1),zeros(6,1),zeros(6,3),zeros(6,6)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RRPRR10_inertiaJ_slag_vp1: qJ has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S5RRPRR10_inertiaJ_slag_vp1: pkin has to be [10x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5RRPRR10_inertiaJ_slag_vp1: m has to be [6x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [6,3]), ...
  'S5RRPRR10_inertiaJ_slag_vp1: rSges has to be [6x3] (double)');
assert(isreal(Icges) && all(size(Icges) == [6 6]), ...
  'S5RRPRR10_inertiaJ_slag_vp1: Icges has to be [6x6] (double)'); 

%% Symbolic Calculation
% From inertia_joint_joint_floatb_twist_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 20:23:50
% EndTime: 2019-12-31 20:24:01
% DurationCPUTime: 4.26s
% Computational Cost: add. (14569->454), mult. (37297->662), div. (0->0), fcn. (48768->12), ass. (0->218)
t247 = sin(pkin(10));
t248 = cos(pkin(10));
t255 = sin(qJ(2));
t257 = cos(qJ(2));
t192 = -t255 * t247 + t257 * t248;
t205 = sin(pkin(5));
t180 = t192 * t205;
t213 = t247 * t257 + t248 * t255;
t181 = t213 * t205;
t206 = cos(pkin(5));
t145 = Icges(4,5) * t181 + Icges(4,6) * t180 + Icges(4,3) * t206;
t146 = Icges(4,4) * t181 + Icges(4,2) * t180 + Icges(4,6) * t206;
t147 = Icges(4,1) * t181 + Icges(4,4) * t180 + Icges(4,5) * t206;
t176 = Icges(3,3) * t206 + (Icges(3,5) * t255 + Icges(3,6) * t257) * t205;
t177 = Icges(3,6) * t206 + (Icges(3,4) * t255 + Icges(3,2) * t257) * t205;
t178 = Icges(3,5) * t206 + (Icges(3,1) * t255 + Icges(3,4) * t257) * t205;
t231 = t205 * t255;
t267 = t205 * t257 * t177 + t180 * t146 + t181 * t147 + t178 * t231 + (t145 + t176) * t206;
t266 = t205 ^ 2;
t265 = m(4) / 0.2e1;
t264 = m(5) / 0.2e1;
t263 = m(6) / 0.2e1;
t182 = t213 * t206;
t209 = sin(qJ(1));
t211 = cos(qJ(1));
t167 = t182 * t211 + t192 * t209;
t208 = sin(qJ(4));
t256 = cos(qJ(4));
t232 = t205 * t256;
t138 = t167 * t208 + t211 * t232;
t169 = -t182 * t209 + t192 * t211;
t140 = t169 * t208 - t209 * t232;
t246 = t205 * t211;
t139 = t167 * t256 - t208 * t246;
t212 = t206 * t192;
t166 = -t209 * t213 + t211 * t212;
t207 = sin(qJ(5));
t210 = cos(qJ(5));
t102 = -t139 * t207 - t166 * t210;
t103 = t139 * t210 - t166 * t207;
t61 = Icges(6,5) * t103 + Icges(6,6) * t102 + Icges(6,3) * t138;
t63 = Icges(6,4) * t103 + Icges(6,2) * t102 + Icges(6,6) * t138;
t65 = Icges(6,1) * t103 + Icges(6,4) * t102 + Icges(6,5) * t138;
t17 = t102 * t63 + t103 * t65 + t138 * t61;
t170 = t181 * t208 - t206 * t256;
t244 = t209 * t205;
t141 = t169 * t256 + t208 * t244;
t168 = -t209 * t212 - t211 * t213;
t104 = -t141 * t207 - t168 * t210;
t105 = t141 * t210 - t168 * t207;
t62 = Icges(6,5) * t105 + Icges(6,6) * t104 + Icges(6,3) * t140;
t64 = Icges(6,4) * t105 + Icges(6,2) * t104 + Icges(6,6) * t140;
t66 = Icges(6,1) * t105 + Icges(6,4) * t104 + Icges(6,5) * t140;
t18 = t102 * t64 + t103 * t66 + t138 * t62;
t171 = t181 * t256 + t206 * t208;
t130 = -t171 * t207 - t180 * t210;
t131 = t171 * t210 - t180 * t207;
t85 = Icges(6,5) * t131 + Icges(6,6) * t130 + Icges(6,3) * t170;
t86 = Icges(6,4) * t131 + Icges(6,2) * t130 + Icges(6,6) * t170;
t87 = Icges(6,1) * t131 + Icges(6,4) * t130 + Icges(6,5) * t170;
t28 = t102 * t86 + t103 * t87 + t138 * t85;
t1 = t138 * t17 + t140 * t18 + t170 * t28;
t262 = -t1 / 0.2e1;
t21 = t130 * t63 + t131 * t65 + t170 * t61;
t22 = t130 * t64 + t131 * t66 + t170 * t62;
t39 = t130 * t86 + t131 * t87 + t170 * t85;
t31 = t39 * t170;
t7 = t21 * t138 + t22 * t140 + t31;
t261 = t7 / 0.2e1;
t260 = t138 / 0.2e1;
t259 = t140 / 0.2e1;
t258 = t170 / 0.2e1;
t254 = t139 * pkin(4);
t253 = t211 * pkin(1);
t219 = -t103 * rSges(6,1) - t102 * rSges(6,2);
t67 = rSges(6,3) * t138 - t219;
t252 = pkin(9) * t138 + t254 + t67;
t68 = rSges(6,1) * t105 + rSges(6,2) * t104 + rSges(6,3) * t140;
t251 = pkin(4) * t141 + pkin(9) * t140 + t68;
t250 = t267 * t206;
t88 = rSges(6,1) * t131 + rSges(6,2) * t130 + rSges(6,3) * t170;
t249 = pkin(4) * t171 + pkin(9) * t170 + t88;
t203 = pkin(2) * t257 + pkin(1);
t245 = t209 * t203;
t183 = t206 * t255 * pkin(2) + (-pkin(7) - qJ(3)) * t205;
t243 = t211 * t183;
t123 = t169 * pkin(3) - pkin(8) * t168;
t197 = t211 * t203;
t160 = -t253 + t197 + (-pkin(7) * t205 - t183) * t209;
t149 = t206 * t160;
t242 = t123 * t206 + t149;
t122 = t167 * pkin(3) - t166 * pkin(8);
t201 = pkin(7) * t246;
t159 = t243 + t201 + (-pkin(1) + t203) * t209;
t241 = -t122 - t159;
t240 = t159 * t244 + t160 * t246;
t107 = Icges(4,5) * t167 + Icges(4,6) * t166 - Icges(4,3) * t246;
t228 = t211 * t257;
t229 = t209 * t255;
t187 = t206 * t228 - t229;
t227 = t211 * t255;
t230 = t209 * t257;
t188 = t206 * t227 + t230;
t150 = Icges(3,5) * t188 + Icges(3,6) * t187 - Icges(3,3) * t246;
t239 = -t150 - t107;
t108 = Icges(4,5) * t169 + Icges(4,6) * t168 + Icges(4,3) * t244;
t189 = -t206 * t230 - t227;
t190 = -t206 * t229 + t228;
t151 = Icges(3,5) * t190 + Icges(3,6) * t189 + Icges(3,3) * t244;
t238 = t151 + t108;
t193 = pkin(2) * t231 + qJ(3) * t206;
t237 = -pkin(3) * t181 + pkin(8) * t180 - t193;
t118 = Icges(5,5) * t171 - Icges(5,6) * t170 - Icges(5,3) * t180;
t119 = Icges(5,4) * t171 - Icges(5,2) * t170 - Icges(5,6) * t180;
t120 = Icges(5,1) * t171 - Icges(5,4) * t170 - Icges(5,5) * t180;
t53 = -t118 * t180 - t119 * t170 + t120 * t171;
t236 = t28 / 0.2e1 + t21 / 0.2e1;
t29 = t104 * t86 + t105 * t87 + t140 * t85;
t235 = t29 / 0.2e1 + t22 / 0.2e1;
t84 = rSges(5,1) * t141 - rSges(5,2) * t140 - t168 * rSges(5,3);
t114 = t169 * rSges(4,1) + t168 * rSges(4,2) + rSges(4,3) * t244;
t157 = t190 * rSges(3,1) + t189 * rSges(3,2) + rSges(3,3) * t244;
t226 = t205 * (-rSges(4,1) * t181 - rSges(4,2) * t180 - rSges(4,3) * t206 - t193);
t225 = -t183 * t209 + t197;
t224 = t122 * t244 + t123 * t246 + t240;
t121 = rSges(5,1) * t171 - rSges(5,2) * t170 - rSges(5,3) * t180;
t223 = t205 * (-t121 + t237);
t220 = -t167 * rSges(4,1) - t166 * rSges(4,2);
t218 = t205 * (t237 - t249);
t77 = Icges(5,5) * t139 - Icges(5,6) * t138 - Icges(5,3) * t166;
t79 = Icges(5,4) * t139 - Icges(5,2) * t138 - Icges(5,6) * t166;
t81 = Icges(5,1) * t139 - Icges(5,4) * t138 - Icges(5,5) * t166;
t40 = -t170 * t79 + t171 * t81 - t180 * t77;
t48 = -t118 * t166 - t119 * t138 + t120 * t139;
t217 = -t48 / 0.2e1 - t40 / 0.2e1 - t236;
t78 = Icges(5,5) * t141 - Icges(5,6) * t140 - Icges(5,3) * t168;
t80 = Icges(5,4) * t141 - Icges(5,2) * t140 - Icges(5,6) * t168;
t82 = Icges(5,1) * t141 - Icges(5,4) * t140 - Icges(5,5) * t168;
t41 = -t170 * t80 + t171 * t82 - t180 * t78;
t49 = -t118 * t168 - t119 * t140 + t120 * t141;
t216 = -t41 / 0.2e1 - t49 / 0.2e1 - t235;
t215 = t123 + t225;
t83 = rSges(5,1) * t139 - rSges(5,2) * t138 - rSges(5,3) * t166;
t156 = rSges(3,1) * t188 + rSges(3,2) * t187 - rSges(3,3) * t246;
t214 = -t122 - t243 - t245;
t195 = rSges(2,1) * t211 - rSges(2,2) * t209;
t194 = -rSges(2,1) * t209 - rSges(2,2) * t211;
t179 = t206 * rSges(3,3) + (rSges(3,1) * t255 + rSges(3,2) * t257) * t205;
t155 = Icges(3,1) * t190 + Icges(3,4) * t189 + Icges(3,5) * t244;
t154 = Icges(3,1) * t188 + Icges(3,4) * t187 - Icges(3,5) * t246;
t153 = Icges(3,4) * t190 + Icges(3,2) * t189 + Icges(3,6) * t244;
t152 = Icges(3,4) * t188 + Icges(3,2) * t187 - Icges(3,6) * t246;
t137 = pkin(7) * t244 + t157 + t253;
t136 = -pkin(1) * t209 - t156 + t201;
t126 = -t156 * t206 - t179 * t246;
t125 = t157 * t206 - t179 * t244;
t113 = -rSges(4,3) * t246 - t220;
t112 = Icges(4,1) * t169 + Icges(4,4) * t168 + Icges(4,5) * t244;
t111 = Icges(4,1) * t167 + Icges(4,4) * t166 - Icges(4,5) * t246;
t110 = Icges(4,4) * t169 + Icges(4,2) * t168 + Icges(4,6) * t244;
t109 = Icges(4,4) * t167 + Icges(4,2) * t166 - Icges(4,6) * t246;
t106 = (t156 * t209 + t157 * t211) * t205;
t101 = t176 * t244 + t177 * t189 + t178 * t190;
t100 = -t176 * t246 + t177 * t187 + t178 * t188;
t92 = t225 + t114;
t91 = -t245 + (rSges(4,3) * t205 - t183) * t211 + t220;
t90 = t206 * t151 + (t153 * t257 + t155 * t255) * t205;
t89 = t206 * t150 + (t152 * t257 + t154 * t255) * t205;
t72 = (-t113 - t159) * t206 + t211 * t226;
t71 = t206 * t114 + t209 * t226 + t149;
t70 = t145 * t244 + t146 * t168 + t147 * t169;
t69 = -t145 * t246 + t146 * t166 + t147 * t167;
t60 = t215 + t84;
t59 = t214 - t83;
t58 = (t113 * t209 + t114 * t211) * t205 + t240;
t57 = t121 * t168 - t180 * t84;
t56 = -t121 * t166 + t180 * t83;
t55 = t108 * t206 + t110 * t180 + t112 * t181;
t54 = t107 * t206 + t109 * t180 + t111 * t181;
t52 = t53 * t206;
t51 = t53 * t180;
t50 = t166 * t84 - t168 * t83;
t47 = (-t83 + t241) * t206 + t211 * t223;
t46 = t206 * t84 + t209 * t223 + t242;
t45 = t215 + t251;
t44 = -t254 + (-rSges(6,3) - pkin(9)) * t138 + t214 + t219;
t43 = -t140 * t88 + t170 * t68;
t42 = t138 * t88 - t170 * t67;
t38 = (t209 * t83 + t211 * t84) * t205 + t224;
t37 = t39 * t206;
t36 = -t140 * t80 + t141 * t82 - t168 * t78;
t35 = -t140 * t79 + t141 * t81 - t168 * t77;
t34 = -t138 * t80 + t139 * t82 - t166 * t78;
t33 = -t138 * t79 + t139 * t81 - t166 * t77;
t32 = t39 * t180;
t30 = -t138 * t68 + t140 * t67;
t27 = t168 * t249 - t180 * t251;
t26 = -t166 * t249 + t180 * t252;
t25 = (t241 - t252) * t206 + t211 * t218;
t24 = t206 * t251 + t209 * t218 + t242;
t23 = t166 * t251 - t168 * t252;
t20 = t104 * t64 + t105 * t66 + t140 * t62;
t19 = t104 * t63 + t105 * t65 + t140 * t61;
t16 = (t209 * t252 + t211 * t251) * t205 + t224;
t15 = t52 + (t41 * t209 - t40 * t211) * t205;
t14 = -t40 * t166 - t41 * t168 - t51;
t13 = t49 * t206 + (t209 * t36 - t211 * t35) * t205;
t12 = t48 * t206 + (t209 * t34 - t211 * t33) * t205;
t11 = -t166 * t35 - t168 * t36 - t180 * t49;
t10 = -t166 * t33 - t168 * t34 - t180 * t48;
t9 = t37 + (t22 * t209 - t21 * t211) * t205;
t8 = -t21 * t166 - t22 * t168 - t32;
t6 = t29 * t206 + (-t19 * t211 + t20 * t209) * t205;
t5 = t28 * t206 + (-t17 * t211 + t18 * t209) * t205;
t4 = -t166 * t19 - t168 * t20 - t180 * t29;
t3 = -t166 * t17 - t168 * t18 - t180 * t28;
t2 = t138 * t19 + t140 * t20 + t170 * t29;
t73 = [Icges(2,3) + m(6) * (t44 ^ 2 + t45 ^ 2) + m(5) * (t59 ^ 2 + t60 ^ 2) + m(4) * (t91 ^ 2 + t92 ^ 2) + m(3) * (t136 ^ 2 + t137 ^ 2) + m(2) * (t194 ^ 2 + t195 ^ 2) + t53 + t39 + t267; t37 + t52 + m(6) * (t24 * t45 + t25 * t44) + m(5) * (t46 * t60 + t47 * t59) + m(4) * (t71 * t92 + t72 * t91) + m(3) * (t125 * t137 + t126 * t136) + ((-t54 / 0.2e1 - t100 / 0.2e1 - t69 / 0.2e1 - t89 / 0.2e1 + t217) * t211 + (t90 / 0.2e1 + t55 / 0.2e1 + t101 / 0.2e1 + t70 / 0.2e1 - t216) * t209) * t205 + t250; (t9 + t15 + t250) * t206 + m(6) * (t16 ^ 2 + t24 ^ 2 + t25 ^ 2) + m(5) * (t38 ^ 2 + t46 ^ 2 + t47 ^ 2) + m(4) * (t58 ^ 2 + t71 ^ 2 + t72 ^ 2) + m(3) * (t106 ^ 2 + t125 ^ 2 + t126 ^ 2) + ((-t12 - t5 + ((t166 * t109 + t167 * t111 + t187 * t152 + t188 * t154) * t205 + t239 * t266 * t211) * t211 + (-t100 - t54 - t69 - t89) * t206) * t211 + (t6 + t13 + ((t168 * t110 + t169 * t112 + t189 * t153 + t190 * t155) * t205 + t238 * t266 * t209) * t209 + (t101 + t70 + t90 + t55) * t206 + (-t168 * t109 - t166 * t110 - t169 * t111 - t167 * t112 - t189 * t152 - t187 * t153 - t190 * t154 - t188 * t155 + (t209 * t239 + t211 * t238) * t205) * t246) * t209) * t205; 0.2e1 * ((t209 * t44 - t211 * t45) * t263 + (t209 * t59 - t211 * t60) * t264 + (t209 * t91 - t211 * t92) * t265) * t205; m(6) * (t206 * t16 + (t209 * t25 - t211 * t24) * t205) + m(5) * (t206 * t38 + (t209 * t47 - t211 * t46) * t205) + m(4) * (t206 * t58 + (t209 * t72 - t211 * t71) * t205); 0.2e1 * (t265 + t264 + t263) * (t206 ^ 2 + (t209 ^ 2 + t211 ^ 2) * t266); -t32 - t51 + m(6) * (t26 * t44 + t27 * t45) + m(5) * (t56 * t59 + t57 * t60) + t216 * t168 + t217 * t166; (t8 / 0.2e1 + t14 / 0.2e1) * t206 - (t9 / 0.2e1 + t15 / 0.2e1) * t180 + (-t6 / 0.2e1 - t13 / 0.2e1) * t168 + (-t5 / 0.2e1 - t12 / 0.2e1) * t166 + m(6) * (t16 * t23 + t24 * t27 + t25 * t26) + m(5) * (t38 * t50 + t46 * t57 + t47 * t56) + ((-t3 / 0.2e1 - t10 / 0.2e1) * t211 + (t4 / 0.2e1 + t11 / 0.2e1) * t209) * t205; m(5) * (t50 * t206 + (t209 * t56 - t211 * t57) * t205) + m(6) * (t23 * t206 + (t209 * t26 - t211 * t27) * t205); -(t8 + t14) * t180 + (-t4 - t11) * t168 + (-t3 - t10) * t166 + m(6) * (t23 ^ 2 + t26 ^ 2 + t27 ^ 2) + m(5) * (t50 ^ 2 + t56 ^ 2 + t57 ^ 2); t31 + m(6) * (t42 * t44 + t43 * t45) + t235 * t140 + t236 * t138; t206 * t261 + m(6) * (t16 * t30 + t24 * t43 + t25 * t42) + t9 * t258 + t5 * t260 + t6 * t259 + (t211 * t262 + t209 * t2 / 0.2e1) * t205; m(6) * (t30 * t206 + (t209 * t42 - t211 * t43) * t205); t166 * t262 + m(6) * (t23 * t30 + t26 * t42 + t27 * t43) + t4 * t259 + t3 * t260 + t8 * t258 - t168 * t2 / 0.2e1 - t180 * t261; m(6) * (t30 ^ 2 + t42 ^ 2 + t43 ^ 2) + t140 * t2 + t138 * t1 + t170 * t7;];
%% Postprocessing: Reshape Output
% From vec2symmat_5_matlab.m
res = [t73(1), t73(2), t73(4), t73(7), t73(11); t73(2), t73(3), t73(5), t73(8), t73(12); t73(4), t73(5), t73(6), t73(9), t73(13); t73(7), t73(8), t73(9), t73(10), t73(14); t73(11), t73(12), t73(13), t73(14), t73(15);];
Mq = res;
