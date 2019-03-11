% Calculate joint inertia matrix for
% S6PPRRPR2
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% pkin [12x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,alpha2,alpha3,d3,d4,d6,theta1,theta2]';
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
% Datum: 2019-03-08 18:51
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Mq = S6PPRRPR2_inertiaJ_slag_vp1(qJ, ...
  pkin, m, rSges, Icges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(12,1),zeros(7,1),zeros(7,3),zeros(7,6)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6PPRRPR2_inertiaJ_slag_vp1: qJ has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [12 1]), ...
  'S6PPRRPR2_inertiaJ_slag_vp1: pkin has to be [12x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6PPRRPR2_inertiaJ_slag_vp1: m has to be [7x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [7,3]), ...
  'S6PPRRPR2_inertiaJ_slag_vp1: rSges has to be [7x3] (double)');
assert(isreal(Icges) && all(size(Icges) == [7 6]), ...
  'S6PPRRPR2_inertiaJ_slag_vp1: Icges has to be [7x6] (double)'); 

%% Symbolic Calculation
% From inertia_joint_joint_floatb_twist_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-08 18:49:10
% EndTime: 2019-03-08 18:49:17
% DurationCPUTime: 3.80s
% Computational Cost: add. (24646->453), mult. (68455->659), div. (0->0), fcn. (90213->14), ass. (0->202)
t197 = m(7) / 0.2e1 + m(6) / 0.2e1;
t219 = 0.2e1 * t197;
t174 = sin(pkin(11));
t176 = cos(pkin(11));
t177 = cos(pkin(6));
t203 = sin(pkin(12));
t190 = t177 * t203;
t205 = cos(pkin(12));
t168 = t174 * t205 + t176 * t190;
t180 = sin(qJ(3));
t191 = t177 * t205;
t185 = t174 * t203 - t176 * t191;
t206 = cos(pkin(7));
t183 = t185 * t206;
t175 = sin(pkin(6));
t204 = sin(pkin(7));
t192 = t175 * t204;
t212 = cos(qJ(3));
t152 = t168 * t212 + (-t176 * t192 - t183) * t180;
t193 = t175 * t206;
t161 = -t176 * t193 + t185 * t204;
t179 = sin(qJ(4));
t211 = cos(qJ(4));
t140 = t152 * t211 + t161 * t179;
t169 = -t174 * t190 + t176 * t205;
t184 = t174 * t191 + t176 * t203;
t182 = t184 * t206;
t154 = t169 * t212 + (t174 * t192 - t182) * t180;
t162 = t174 * t193 + t184 * t204;
t142 = t154 * t211 + t162 * t179;
t187 = t206 * t205;
t160 = t177 * t204 * t180 + (t180 * t187 + t203 * t212) * t175;
t167 = t177 * t206 - t205 * t192;
t156 = t160 * t211 + t167 * t179;
t139 = t152 * t179 - t161 * t211;
t188 = t212 * t204;
t186 = t175 * t188;
t151 = t168 * t180 + t176 * t186 + t212 * t183;
t178 = sin(qJ(6));
t181 = cos(qJ(6));
t110 = t139 * t181 - t151 * t178;
t111 = t139 * t178 + t151 * t181;
t70 = Icges(7,5) * t111 + Icges(7,6) * t110 + Icges(7,3) * t140;
t72 = Icges(7,4) * t111 + Icges(7,2) * t110 + Icges(7,6) * t140;
t74 = Icges(7,1) * t111 + Icges(7,4) * t110 + Icges(7,5) * t140;
t24 = t110 * t72 + t111 * t74 + t140 * t70;
t141 = t154 * t179 - t162 * t211;
t153 = t169 * t180 - t174 * t186 + t212 * t182;
t112 = t141 * t181 - t153 * t178;
t113 = t141 * t178 + t153 * t181;
t71 = Icges(7,5) * t113 + Icges(7,6) * t112 + Icges(7,3) * t142;
t73 = Icges(7,4) * t113 + Icges(7,2) * t112 + Icges(7,6) * t142;
t75 = Icges(7,1) * t113 + Icges(7,4) * t112 + Icges(7,5) * t142;
t25 = t110 * t73 + t111 * t75 + t140 * t71;
t155 = t160 * t179 - t167 * t211;
t159 = -t177 * t188 + (t180 * t203 - t187 * t212) * t175;
t143 = t155 * t181 - t159 * t178;
t144 = t155 * t178 + t159 * t181;
t100 = Icges(7,1) * t144 + Icges(7,4) * t143 + Icges(7,5) * t156;
t98 = Icges(7,5) * t144 + Icges(7,6) * t143 + Icges(7,3) * t156;
t99 = Icges(7,4) * t144 + Icges(7,2) * t143 + Icges(7,6) * t156;
t35 = t100 * t111 + t110 * t99 + t140 * t98;
t1 = t140 * t24 + t142 * t25 + t156 * t35;
t218 = t1 / 0.2e1;
t26 = t112 * t72 + t113 * t74 + t142 * t70;
t27 = t112 * t73 + t113 * t75 + t142 * t71;
t36 = t100 * t113 + t112 * t99 + t142 * t98;
t2 = t140 * t26 + t142 * t27 + t156 * t36;
t217 = t2 / 0.2e1;
t30 = t143 * t72 + t144 * t74 + t156 * t70;
t31 = t143 * t73 + t144 * t75 + t156 * t71;
t49 = t100 * t144 + t143 * t99 + t156 * t98;
t7 = t140 * t30 + t142 * t31 + t156 * t49;
t216 = t7 / 0.2e1;
t215 = t140 / 0.2e1;
t214 = t142 / 0.2e1;
t213 = t156 / 0.2e1;
t107 = pkin(4) * t140 + qJ(5) * t139;
t94 = rSges(6,1) * t151 - rSges(6,2) * t140 + rSges(6,3) * t139;
t210 = -t107 - t94;
t108 = pkin(4) * t142 + qJ(5) * t141;
t95 = rSges(6,1) * t153 - rSges(6,2) * t142 + rSges(6,3) * t141;
t209 = -t108 - t95;
t76 = rSges(7,1) * t111 + rSges(7,2) * t110 + rSges(7,3) * t140;
t208 = pkin(5) * t151 + pkin(10) * t140 + t76;
t77 = rSges(7,1) * t113 + rSges(7,2) * t112 + rSges(7,3) * t142;
t207 = pkin(5) * t153 + pkin(10) * t142 + t77;
t101 = rSges(7,1) * t144 + rSges(7,2) * t143 + rSges(7,3) * t156;
t202 = pkin(5) * t159 + pkin(10) * t156 + t101;
t135 = pkin(3) * t152 + pkin(9) * t151;
t132 = t162 * t135;
t201 = t162 * t107 + t132;
t136 = pkin(3) * t154 + pkin(9) * t153;
t134 = t167 * t136;
t200 = t167 * t108 + t134;
t130 = rSges(6,1) * t159 - rSges(6,2) * t156 + rSges(6,3) * t155;
t137 = pkin(4) * t156 + qJ(5) * t155;
t199 = -t130 - t137;
t150 = pkin(3) * t160 + pkin(9) * t159;
t138 = t161 * t150;
t198 = t161 * t137 + t138;
t196 = -t107 - t208;
t195 = -t108 - t207;
t194 = -t137 - t202;
t189 = m(3) + m(4) + m(5) + m(6) + m(7);
t149 = rSges(4,1) * t160 - rSges(4,2) * t159 + rSges(4,3) * t167;
t148 = Icges(4,1) * t160 - Icges(4,4) * t159 + Icges(4,5) * t167;
t147 = Icges(4,4) * t160 - Icges(4,2) * t159 + Icges(4,6) * t167;
t146 = Icges(4,5) * t160 - Icges(4,6) * t159 + Icges(4,3) * t167;
t131 = rSges(5,1) * t156 - rSges(5,2) * t155 + rSges(5,3) * t159;
t129 = Icges(5,1) * t156 - Icges(5,4) * t155 + Icges(5,5) * t159;
t128 = Icges(6,1) * t159 - Icges(6,4) * t156 + Icges(6,5) * t155;
t127 = Icges(5,4) * t156 - Icges(5,2) * t155 + Icges(5,6) * t159;
t126 = Icges(6,4) * t159 - Icges(6,2) * t156 + Icges(6,6) * t155;
t125 = Icges(5,5) * t156 - Icges(5,6) * t155 + Icges(5,3) * t159;
t124 = Icges(6,5) * t159 - Icges(6,6) * t156 + Icges(6,3) * t155;
t123 = rSges(4,1) * t154 - rSges(4,2) * t153 + rSges(4,3) * t162;
t122 = rSges(4,1) * t152 - rSges(4,2) * t151 + rSges(4,3) * t161;
t121 = Icges(4,1) * t154 - Icges(4,4) * t153 + Icges(4,5) * t162;
t120 = Icges(4,1) * t152 - Icges(4,4) * t151 + Icges(4,5) * t161;
t119 = Icges(4,4) * t154 - Icges(4,2) * t153 + Icges(4,6) * t162;
t118 = Icges(4,4) * t152 - Icges(4,2) * t151 + Icges(4,6) * t161;
t117 = Icges(4,5) * t154 - Icges(4,6) * t153 + Icges(4,3) * t162;
t116 = Icges(4,5) * t152 - Icges(4,6) * t151 + Icges(4,3) * t161;
t109 = t151 * t137;
t104 = t159 * t108;
t102 = t153 * t107;
t97 = rSges(5,1) * t142 - rSges(5,2) * t141 + rSges(5,3) * t153;
t96 = rSges(5,1) * t140 - rSges(5,2) * t139 + rSges(5,3) * t151;
t93 = Icges(5,1) * t142 - Icges(5,4) * t141 + Icges(5,5) * t153;
t92 = Icges(5,1) * t140 - Icges(5,4) * t139 + Icges(5,5) * t151;
t91 = Icges(6,1) * t153 - Icges(6,4) * t142 + Icges(6,5) * t141;
t90 = Icges(6,1) * t151 - Icges(6,4) * t140 + Icges(6,5) * t139;
t89 = Icges(5,4) * t142 - Icges(5,2) * t141 + Icges(5,6) * t153;
t88 = Icges(5,4) * t140 - Icges(5,2) * t139 + Icges(5,6) * t151;
t87 = Icges(6,4) * t153 - Icges(6,2) * t142 + Icges(6,6) * t141;
t86 = Icges(6,4) * t151 - Icges(6,2) * t140 + Icges(6,6) * t139;
t85 = Icges(5,5) * t142 - Icges(5,6) * t141 + Icges(5,3) * t153;
t84 = Icges(5,5) * t140 - Icges(5,6) * t139 + Icges(5,3) * t151;
t83 = Icges(6,5) * t153 - Icges(6,6) * t142 + Icges(6,3) * t141;
t82 = Icges(6,5) * t151 - Icges(6,6) * t140 + Icges(6,3) * t139;
t80 = t123 * t167 - t149 * t162;
t79 = -t122 * t167 + t149 * t161;
t78 = t122 * t162 - t123 * t161;
t69 = -t131 * t153 + t159 * t97;
t68 = t131 * t151 - t159 * t96;
t67 = t125 * t159 - t127 * t155 + t129 * t156;
t66 = t124 * t155 - t126 * t156 + t128 * t159;
t65 = -t151 * t97 + t153 * t96;
t64 = t167 * t97 + t134 + (-t131 - t150) * t162;
t63 = t131 * t161 + t138 + (-t135 - t96) * t167;
t62 = t125 * t153 - t127 * t141 + t129 * t142;
t61 = t125 * t151 - t127 * t139 + t129 * t140;
t60 = t124 * t141 - t126 * t142 + t128 * t153;
t59 = t124 * t139 - t126 * t140 + t128 * t151;
t58 = -t101 * t142 + t156 * t77;
t57 = t101 * t140 - t156 * t76;
t56 = t162 * t96 + t132 + (-t136 - t97) * t161;
t55 = t199 * t153 + t159 * t95 + t104;
t54 = t130 * t151 + t210 * t159 + t109;
t53 = -t155 * t89 + t156 * t93 + t159 * t85;
t52 = -t155 * t88 + t156 * t92 + t159 * t84;
t51 = t155 * t83 - t156 * t87 + t159 * t91;
t50 = t155 * t82 - t156 * t86 + t159 * t90;
t48 = -t141 * t89 + t142 * t93 + t153 * t85;
t47 = -t141 * t88 + t142 * t92 + t153 * t84;
t46 = -t139 * t89 + t140 * t93 + t151 * t85;
t45 = -t139 * t88 + t140 * t92 + t151 * t84;
t44 = t141 * t83 - t142 * t87 + t153 * t91;
t43 = t141 * t82 - t142 * t86 + t153 * t90;
t42 = t139 * t83 - t140 * t87 + t151 * t91;
t41 = t139 * t82 - t140 * t86 + t151 * t90;
t40 = -t140 * t77 + t142 * t76;
t39 = t167 * t95 + (-t150 + t199) * t162 + t200;
t38 = t130 * t161 + (-t135 + t210) * t167 + t198;
t37 = t209 * t151 + t153 * t94 + t102;
t34 = t162 * t94 + (-t136 + t209) * t161 + t201;
t33 = t194 * t153 + t207 * t159 + t104;
t32 = t202 * t151 + t196 * t159 + t109;
t29 = t207 * t167 + (-t150 + t194) * t162 + t200;
t28 = t202 * t161 + (-t135 + t196) * t167 + t198;
t23 = t195 * t151 + t208 * t153 + t102;
t22 = t208 * t162 + (-t136 + t195) * t161 + t201;
t21 = t161 * t52 + t162 * t53 + t167 * t67;
t20 = t161 * t50 + t162 * t51 + t167 * t66;
t19 = t151 * t52 + t153 * t53 + t159 * t67;
t18 = t151 * t50 + t153 * t51 + t159 * t66;
t17 = t161 * t47 + t162 * t48 + t167 * t62;
t16 = t161 * t45 + t162 * t46 + t167 * t61;
t15 = t161 * t43 + t162 * t44 + t167 * t60;
t14 = t161 * t41 + t162 * t42 + t167 * t59;
t13 = t151 * t47 + t153 * t48 + t159 * t62;
t12 = t151 * t45 + t153 * t46 + t159 * t61;
t11 = t151 * t43 + t153 * t44 + t159 * t60;
t10 = t151 * t41 + t153 * t42 + t159 * t59;
t9 = t161 * t30 + t162 * t31 + t167 * t49;
t8 = t151 * t30 + t153 * t31 + t159 * t49;
t6 = t161 * t26 + t162 * t27 + t167 * t36;
t5 = t161 * t24 + t162 * t25 + t167 * t35;
t4 = t151 * t26 + t153 * t27 + t159 * t36;
t3 = t151 * t24 + t153 * t25 + t159 * t35;
t81 = [m(2) + t189; t189 * t177; 0.2e1 * (m(3) / 0.2e1 + m(4) / 0.2e1 + m(5) / 0.2e1 + t197) * (t177 ^ 2 + (t174 ^ 2 + t176 ^ 2) * t175 ^ 2); m(4) * t78 + m(5) * t56 + m(6) * t34 + m(7) * t22; m(7) * (t177 * t22 + (t174 * t28 - t176 * t29) * t175) + m(4) * (t177 * t78 + (t174 * t79 - t176 * t80) * t175) + m(5) * (t177 * t56 + (t174 * t63 - t176 * t64) * t175) + m(6) * (t177 * t34 + (t174 * t38 - t176 * t39) * t175); (t9 + t21 + t20 + (t146 * t167 - t147 * t159 + t148 * t160) * t167) * t167 + (t6 + t17 + t15 + (t117 * t162 - t119 * t153 + t121 * t154) * t162 + (t117 * t167 - t119 * t159 + t121 * t160 + t146 * t162 - t147 * t153 + t148 * t154) * t167) * t162 + (t5 + t14 + t16 + (t116 * t161 - t118 * t151 + t120 * t152) * t161 + (t116 * t167 - t118 * t159 + t120 * t160 + t146 * t161 - t147 * t151 + t148 * t152) * t167 + (t116 * t162 + t117 * t161 - t153 * t118 - t119 * t151 + t154 * t120 + t121 * t152) * t162) * t161 + m(7) * (t22 ^ 2 + t28 ^ 2 + t29 ^ 2) + m(6) * (t34 ^ 2 + t38 ^ 2 + t39 ^ 2) + m(5) * (t56 ^ 2 + t63 ^ 2 + t64 ^ 2) + m(4) * (t78 ^ 2 + t79 ^ 2 + t80 ^ 2); m(5) * t65 + m(6) * t37 + m(7) * t23; m(7) * (t177 * t23 + (t174 * t32 - t176 * t33) * t175) + m(5) * (t177 * t65 + (t174 * t68 - t176 * t69) * t175) + m(6) * (t177 * t37 + (t174 * t54 - t176 * t55) * t175); (t8 / 0.2e1 + t19 / 0.2e1 + t18 / 0.2e1) * t167 + (t4 / 0.2e1 + t13 / 0.2e1 + t11 / 0.2e1) * t162 + (t3 / 0.2e1 + t10 / 0.2e1 + t12 / 0.2e1) * t161 + (t9 / 0.2e1 + t21 / 0.2e1 + t20 / 0.2e1) * t159 + (t6 / 0.2e1 + t17 / 0.2e1 + t15 / 0.2e1) * t153 + (t5 / 0.2e1 + t14 / 0.2e1 + t16 / 0.2e1) * t151 + m(7) * (t22 * t23 + t28 * t32 + t29 * t33) + m(6) * (t34 * t37 + t38 * t54 + t39 * t55) + m(5) * (t56 * t65 + t63 * t68 + t64 * t69); (t8 + t19 + t18) * t159 + (t4 + t11 + t13) * t153 + (t3 + t12 + t10) * t151 + m(7) * (t23 ^ 2 + t32 ^ 2 + t33 ^ 2) + m(6) * (t37 ^ 2 + t54 ^ 2 + t55 ^ 2) + m(5) * (t65 ^ 2 + t68 ^ 2 + t69 ^ 2); t155 * t219; (t155 * t177 + (-t139 * t176 + t141 * t174) * t175) * t219; m(7) * (t139 * t29 + t141 * t28 + t155 * t22) + m(6) * (t139 * t39 + t141 * t38 + t155 * t34); m(7) * (t139 * t33 + t141 * t32 + t155 * t23) + m(6) * (t139 * t55 + t141 * t54 + t155 * t37); (t139 ^ 2 + t141 ^ 2 + t155 ^ 2) * t219; m(7) * t40; m(7) * (t177 * t40 + (t174 * t57 - t176 * t58) * t175); m(7) * (t22 * t40 + t28 * t57 + t29 * t58) + t6 * t214 + t9 * t213 + t162 * t217 + t161 * t218 + t167 * t216 + t5 * t215; m(7) * (t23 * t40 + t32 * t57 + t33 * t58) + t153 * t217 + t3 * t215 + t4 * t214 + t151 * t218 + t8 * t213 + t159 * t216; m(7) * (t139 * t58 + t141 * t57 + t155 * t40); t142 * t2 + t140 * t1 + t156 * t7 + m(7) * (t40 ^ 2 + t57 ^ 2 + t58 ^ 2);];
%% Postprocessing: Reshape Output
% From vec2symmat_6_matlab.m
res = [t81(1) t81(2) t81(4) t81(7) t81(11) t81(16); t81(2) t81(3) t81(5) t81(8) t81(12) t81(17); t81(4) t81(5) t81(6) t81(9) t81(13) t81(18); t81(7) t81(8) t81(9) t81(10) t81(14) t81(19); t81(11) t81(12) t81(13) t81(14) t81(15) t81(20); t81(16) t81(17) t81(18) t81(19) t81(20) t81(21);];
Mq  = res;
