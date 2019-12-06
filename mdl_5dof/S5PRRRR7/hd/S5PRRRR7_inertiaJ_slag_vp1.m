% Calculate joint inertia matrix for
% S5PRRRR7
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% pkin [9x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d2,d3,d4,d5,theta1]';
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
% Datum: 2019-12-05 17:13
% Revision: 77da58f92bca3eff71542919beafa37024070d86 (2019-12-05)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Mq = S5PRRRR7_inertiaJ_slag_vp1(qJ, ...
  pkin, m, rSges, Icges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(9,1),zeros(6,1),zeros(6,3),zeros(6,6)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5PRRRR7_inertiaJ_slag_vp1: qJ has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S5PRRRR7_inertiaJ_slag_vp1: pkin has to be [9x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5PRRRR7_inertiaJ_slag_vp1: m has to be [6x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [6,3]), ...
  'S5PRRRR7_inertiaJ_slag_vp1: rSges has to be [6x3] (double)');
assert(isreal(Icges) && all(size(Icges) == [6 6]), ...
  'S5PRRRR7_inertiaJ_slag_vp1: Icges has to be [6x6] (double)'); 

%% Symbolic Calculation
% From inertia_joint_joint_floatb_twist_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-05 17:11:48
% EndTime: 2019-12-05 17:11:55
% DurationCPUTime: 2.18s
% Computational Cost: add. (9164->327), mult. (10544->503), div. (0->0), fcn. (11580->10), ass. (0->178)
t166 = cos(pkin(9));
t222 = t166 ^ 2;
t165 = sin(pkin(9));
t223 = t165 ^ 2;
t224 = t222 + t223;
t164 = qJ(3) + qJ(4);
t159 = sin(t164);
t160 = cos(t164);
t168 = sin(qJ(2));
t170 = cos(qJ(2));
t124 = -Icges(5,6) * t170 + (Icges(5,4) * t160 - Icges(5,2) * t159) * t168;
t125 = -Icges(5,5) * t170 + (Icges(5,1) * t160 - Icges(5,4) * t159) * t168;
t208 = t165 * t170;
t137 = -t159 * t208 - t166 * t160;
t138 = -t166 * t159 + t160 * t208;
t123 = -Icges(5,3) * t170 + (Icges(5,5) * t160 - Icges(5,6) * t159) * t168;
t201 = t170 * t123;
t209 = t165 * t168;
t88 = Icges(5,5) * t138 + Icges(5,6) * t137 + Icges(5,3) * t209;
t90 = Icges(5,4) * t138 + Icges(5,2) * t137 + Icges(5,6) * t209;
t92 = Icges(5,1) * t138 + Icges(5,4) * t137 + Icges(5,5) * t209;
t42 = t137 * t90 + t138 * t92 + t88 * t209;
t205 = t166 * t170;
t139 = -t159 * t205 + t165 * t160;
t140 = t165 * t159 + t160 * t205;
t206 = t166 * t168;
t89 = Icges(5,5) * t140 + Icges(5,6) * t139 + Icges(5,3) * t206;
t91 = Icges(5,4) * t140 + Icges(5,2) * t139 + Icges(5,6) * t206;
t93 = Icges(5,1) * t140 + Icges(5,4) * t139 + Icges(5,5) * t206;
t43 = t137 * t91 + t138 * t93 + t89 * t209;
t11 = -(t137 * t124 + t138 * t125) * t170 + (t43 * t166 + (t42 - t201) * t165) * t168;
t44 = t139 * t90 + t140 * t92 + t88 * t206;
t45 = t139 * t91 + t140 * t93 + t89 * t206;
t12 = -(t139 * t124 + t140 * t125) * t170 + (t44 * t165 + (t45 - t201) * t166) * t168;
t225 = t11 * t209 + t12 * t206;
t221 = t170 ^ 2;
t161 = qJ(5) + t164;
t156 = sin(t161);
t157 = cos(t161);
t119 = -Icges(6,6) * t170 + (Icges(6,4) * t157 - Icges(6,2) * t156) * t168;
t120 = -Icges(6,5) * t170 + (Icges(6,1) * t157 - Icges(6,4) * t156) * t168;
t127 = -t156 * t208 - t166 * t157;
t128 = -t166 * t156 + t157 * t208;
t118 = -Icges(6,3) * t170 + (Icges(6,5) * t157 - Icges(6,6) * t156) * t168;
t202 = t170 * t118;
t79 = Icges(6,5) * t128 + Icges(6,6) * t127 + Icges(6,3) * t209;
t81 = Icges(6,4) * t128 + Icges(6,2) * t127 + Icges(6,6) * t209;
t83 = Icges(6,1) * t128 + Icges(6,4) * t127 + Icges(6,5) * t209;
t35 = t127 * t81 + t128 * t83 + t79 * t209;
t129 = -t156 * t205 + t165 * t157;
t130 = t165 * t156 + t157 * t205;
t80 = Icges(6,5) * t130 + Icges(6,6) * t129 + Icges(6,3) * t206;
t82 = Icges(6,4) * t130 + Icges(6,2) * t129 + Icges(6,6) * t206;
t84 = Icges(6,1) * t130 + Icges(6,4) * t129 + Icges(6,5) * t206;
t36 = t127 * t82 + t128 * t84 + t80 * t209;
t5 = -(t127 * t119 + t128 * t120) * t170 + (t36 * t166 + (t35 - t202) * t165) * t168;
t37 = t129 * t81 + t130 * t83 + t79 * t206;
t38 = t129 * t82 + t130 * t84 + t80 * t206;
t6 = -(t129 * t119 + t130 * t120) * t170 + (t37 * t165 + (t38 - t202) * t166) * t168;
t220 = t6 * t206 + t5 * t209;
t219 = t165 / 0.2e1;
t218 = -t166 / 0.2e1;
t217 = -t170 / 0.2e1;
t169 = cos(qJ(3));
t215 = t169 * pkin(3);
t194 = pkin(4) * t160;
t172 = pkin(8) * t168 + t194 * t170;
t187 = pkin(4) * t159;
t72 = t172 * t165 - t187 * t166;
t85 = t128 * rSges(6,1) + t127 * rSges(6,2) + rSges(6,3) * t209;
t76 = t85 * t206;
t213 = t72 * t206 + t76;
t86 = t130 * rSges(6,1) + t129 * rSges(6,2) + rSges(6,3) * t206;
t212 = t187 * t165 + t172 * t166 + t86;
t174 = pkin(7) * t168 + t215 * t170;
t167 = sin(qJ(3));
t210 = t165 * t167;
t111 = pkin(3) * t210 + t174 * t166;
t95 = t140 * rSges(5,1) + t139 * rSges(5,2) + rSges(5,3) * t206;
t211 = -t111 - t95;
t121 = -t170 * rSges(6,3) + (rSges(6,1) * t157 - rSges(6,2) * t156) * t168;
t64 = t121 * t209 + t170 * t85;
t126 = -t170 * rSges(5,3) + (rSges(5,1) * t160 - rSges(5,2) * t159) * t168;
t94 = t138 * rSges(5,1) + t137 * rSges(5,2) + rSges(5,3) * t209;
t66 = t126 * t209 + t170 * t94;
t207 = t166 * t167;
t204 = t167 * t170;
t203 = t169 * t170;
t141 = -Icges(4,3) * t170 + (Icges(4,5) * t169 - Icges(4,6) * t167) * t168;
t200 = t170 * t141;
t110 = -pkin(3) * t207 + t174 * t165;
t122 = -pkin(7) * t170 + t215 * t168;
t199 = t170 * t110 + t122 * t209;
t112 = -pkin(8) * t170 + t194 * t168;
t198 = -t112 - t121;
t197 = -t122 - t126;
t144 = -t170 * rSges(4,3) + (rSges(4,1) * t169 - rSges(4,2) * t167) * t168;
t155 = t168 * pkin(2) - t170 * pkin(6);
t196 = -t144 - t155;
t195 = t224 * (pkin(2) * t170 + pkin(6) * t168);
t192 = -t111 - t212;
t191 = -t122 + t198;
t190 = -t155 + t197;
t189 = t209 / 0.2e1;
t188 = t206 / 0.2e1;
t40 = t112 * t209 + t170 * t72 + t64;
t186 = t165 * t110 + t166 * t111 + t195;
t46 = -t170 * t79 + (-t156 * t81 + t157 * t83) * t168;
t47 = -t170 * t80 + (-t156 * t82 + t157 * t84) * t168;
t13 = t221 * t118 + (t47 * t166 + t46 * t165 - (-t119 * t156 + t120 * t157) * t170) * t168;
t185 = -t170 * t13 + t220;
t184 = -t155 + t191;
t19 = t36 * t165 - t35 * t166;
t20 = t38 * t165 - t37 * t166;
t183 = t20 * t188 + t19 * t189 + t5 * t218 + t6 * t219 + (t47 * t165 - t46 * t166) * t217;
t178 = Icges(3,5) * t170 - Icges(3,6) * t168;
t48 = -t170 * t88 + (-t159 * t90 + t160 * t92) * t168;
t49 = -t170 * t89 + (-t159 * t91 + t160 * t93) * t168;
t18 = t221 * t123 + (t49 * t166 + t48 * t165 - (-t124 * t159 + t125 * t160) * t170) * t168;
t175 = (-t13 - t18) * t170 + t220 + t225;
t23 = t43 * t165 - t42 * t166;
t24 = t45 * t165 - t44 * t166;
t173 = t11 * t218 + t12 * t219 + t24 * t188 + t23 * t189 + t183 + (t49 * t165 - t48 * t166) * t217;
t154 = t168 * rSges(3,1) + t170 * rSges(3,2);
t150 = t166 * t203 + t210;
t149 = t165 * t169 - t166 * t204;
t148 = t165 * t203 - t207;
t147 = -t165 * t204 - t166 * t169;
t143 = -Icges(4,5) * t170 + (Icges(4,1) * t169 - Icges(4,4) * t167) * t168;
t142 = -Icges(4,6) * t170 + (Icges(4,4) * t169 - Icges(4,2) * t167) * t168;
t132 = Icges(3,3) * t165 + t178 * t166;
t131 = -Icges(3,3) * t166 + t178 * t165;
t114 = t196 * t166;
t113 = t196 * t165;
t109 = t150 * rSges(4,1) + t149 * rSges(4,2) + rSges(4,3) * t206;
t108 = t148 * rSges(4,1) + t147 * rSges(4,2) + rSges(4,3) * t209;
t106 = Icges(4,1) * t150 + Icges(4,4) * t149 + Icges(4,5) * t206;
t105 = Icges(4,1) * t148 + Icges(4,4) * t147 + Icges(4,5) * t209;
t104 = Icges(4,4) * t150 + Icges(4,2) * t149 + Icges(4,6) * t206;
t103 = Icges(4,4) * t148 + Icges(4,2) * t147 + Icges(4,6) * t209;
t102 = Icges(4,5) * t150 + Icges(4,6) * t149 + Icges(4,3) * t206;
t101 = Icges(4,5) * t148 + Icges(4,6) * t147 + Icges(4,3) * t209;
t97 = t224 * (rSges(3,1) * t170 - rSges(3,2) * t168);
t96 = t110 * t206;
t78 = t94 * t206;
t75 = t190 * t166;
t74 = t190 * t165;
t71 = -t170 * t109 - t144 * t206;
t70 = t170 * t108 + t144 * t209;
t67 = -t126 * t206 - t170 * t95;
t65 = -t121 * t206 - t170 * t86;
t63 = (t108 * t166 - t109 * t165) * t168;
t62 = t184 * t166;
t61 = t184 * t165;
t60 = -t95 * t209 + t78;
t59 = -t86 * t209 + t76;
t58 = t165 * t108 + t166 * t109 + t195;
t57 = -t170 * t102 + (-t104 * t167 + t106 * t169) * t168;
t56 = -t170 * t101 + (-t103 * t167 + t105 * t169) * t168;
t55 = t211 * t170 + t197 * t206;
t54 = t199 + t66;
t53 = t102 * t206 + t149 * t104 + t150 * t106;
t52 = t101 * t206 + t149 * t103 + t150 * t105;
t51 = t102 * t209 + t147 * t104 + t148 * t106;
t50 = t101 * t209 + t147 * t103 + t148 * t105;
t41 = -t212 * t170 + t198 * t206;
t39 = t211 * t209 + t78 + t96;
t34 = t165 * t94 + t166 * t95 + t186;
t33 = -t212 * t209 + t213;
t32 = t192 * t170 + t191 * t206;
t31 = t40 + t199;
t30 = t53 * t165 - t52 * t166;
t29 = t51 * t165 - t50 * t166;
t27 = t192 * t209 + t213 + t96;
t25 = t212 * t166 + (t72 + t85) * t165 + t186;
t17 = -(t149 * t142 + t150 * t143) * t170 + (t52 * t165 + (t53 - t200) * t166) * t168;
t16 = -(t147 * t142 + t148 * t143) * t170 + (t51 * t166 + (t50 - t200) * t165) * t168;
t1 = [m(2) + m(3) + m(4) + m(5) + m(6); m(3) * t97 + m(4) * t58 + m(5) * t34 + m(6) * t25; m(6) * (t25 ^ 2 + t61 ^ 2 + t62 ^ 2) + m(5) * (t34 ^ 2 + t74 ^ 2 + t75 ^ 2) + m(4) * (t113 ^ 2 + t114 ^ 2 + t58 ^ 2) + m(3) * (t224 * t154 ^ 2 + t97 ^ 2) + (-t222 * t131 - t19 - t23 - t29) * t166 + (t223 * t132 + t20 + t24 + t30 + (-t165 * t131 + t166 * t132) * t166) * t165; m(4) * t63 + m(5) * t39 + m(6) * t27; m(6) * (t27 * t25 + t31 * t62 + t32 * t61) + m(5) * (t39 * t34 + t54 * t75 + t55 * t74) + m(4) * (t71 * t113 + t70 * t114 + t63 * t58) + (t29 * t219 + t166 * t30 / 0.2e1) * t168 + t173 + t17 * t219 + t16 * t218 + (t57 * t165 - t56 * t166) * t217; t17 * t206 + t16 * t209 + m(6) * (t27 ^ 2 + t31 ^ 2 + t32 ^ 2) + m(5) * (t39 ^ 2 + t54 ^ 2 + t55 ^ 2) - t170 * t18 + m(4) * (t63 ^ 2 + t70 ^ 2 + t71 ^ 2) - t170 * (t221 * t141 + (t57 * t166 + t56 * t165 - (-t142 * t167 + t143 * t169) * t170) * t168) + t185 + t225; m(5) * t60 + m(6) * t33; m(6) * (t33 * t25 + t40 * t62 + t41 * t61) + m(5) * (t60 * t34 + t66 * t75 + t67 * t74) + t173; m(6) * (t33 * t27 + t40 * t31 + t41 * t32) + m(5) * (t60 * t39 + t66 * t54 + t67 * t55) + t175; m(6) * (t33 ^ 2 + t40 ^ 2 + t41 ^ 2) + m(5) * (t60 ^ 2 + t66 ^ 2 + t67 ^ 2) + t175; m(6) * t59; m(6) * (t59 * t25 + t65 * t61 + t64 * t62) + t183; m(6) * (t59 * t27 + t64 * t31 + t65 * t32) + t185; m(6) * (t59 * t33 + t64 * t40 + t65 * t41) + t185; m(6) * (t59 ^ 2 + t64 ^ 2 + t65 ^ 2) + t185;];
%% Postprocessing: Reshape Output
% From vec2symmat_5_matlab.m
res = [t1(1), t1(2), t1(4), t1(7), t1(11); t1(2), t1(3), t1(5), t1(8), t1(12); t1(4), t1(5), t1(6), t1(9), t1(13); t1(7), t1(8), t1(9), t1(10), t1(14); t1(11), t1(12), t1(13), t1(14), t1(15);];
Mq = res;
