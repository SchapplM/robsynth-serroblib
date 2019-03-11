% Calculate joint inertia matrix for
% S6PPRRRP2
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% pkin [12x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,alpha2,alpha3,d3,d4,d5,theta1,theta2]';
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
% Datum: 2019-03-08 18:58
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Mq = S6PPRRRP2_inertiaJ_slag_vp1(qJ, ...
  pkin, m, rSges, Icges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(12,1),zeros(7,1),zeros(7,3),zeros(7,6)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6PPRRRP2_inertiaJ_slag_vp1: qJ has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [12 1]), ...
  'S6PPRRRP2_inertiaJ_slag_vp1: pkin has to be [12x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6PPRRRP2_inertiaJ_slag_vp1: m has to be [7x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [7,3]), ...
  'S6PPRRRP2_inertiaJ_slag_vp1: rSges has to be [7x3] (double)');
assert(isreal(Icges) && all(size(Icges) == [7 6]), ...
  'S6PPRRRP2_inertiaJ_slag_vp1: Icges has to be [7x6] (double)'); 

%% Symbolic Calculation
% From inertia_joint_joint_floatb_twist_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-08 18:56:06
% EndTime: 2019-03-08 18:56:16
% DurationCPUTime: 5.01s
% Computational Cost: add. (34930->471), mult. (97067->677), div. (0->0), fcn. (128548->14), ass. (0->208)
t229 = rSges(7,1) + pkin(5);
t228 = rSges(7,3) + qJ(6);
t227 = cos(qJ(3));
t226 = cos(qJ(4));
t225 = cos(qJ(5));
t181 = sin(pkin(11));
t183 = cos(pkin(11));
t184 = cos(pkin(6));
t217 = sin(pkin(12));
t196 = t184 * t217;
t219 = cos(pkin(12));
t175 = t181 * t219 + t183 * t196;
t187 = sin(qJ(3));
t197 = t184 * t219;
t191 = t181 * t217 - t183 * t197;
t220 = cos(pkin(7));
t189 = t191 * t220;
t182 = sin(pkin(6));
t218 = sin(pkin(7));
t198 = t182 * t218;
t158 = t175 * t227 + (-t183 * t198 - t189) * t187;
t199 = t182 * t220;
t168 = -t183 * t199 + t191 * t218;
t186 = sin(qJ(4));
t145 = t158 * t226 + t168 * t186;
t194 = t227 * t218;
t192 = t182 * t194;
t157 = t175 * t187 + t183 * t192 + t189 * t227;
t185 = sin(qJ(5));
t121 = t145 * t185 - t157 * t225;
t122 = t145 * t225 + t157 * t185;
t144 = t158 * t186 - t168 * t226;
t224 = rSges(7,2) * t144 + t121 * t228 + t122 * t229;
t176 = -t181 * t196 + t183 * t219;
t190 = t181 * t197 + t183 * t217;
t188 = t190 * t220;
t160 = t176 * t227 + (t181 * t198 - t188) * t187;
t169 = t181 * t199 + t190 * t218;
t147 = t160 * t226 + t169 * t186;
t159 = t176 * t187 - t181 * t192 + t188 * t227;
t123 = t147 * t185 - t159 * t225;
t124 = t147 * t225 + t159 * t185;
t146 = t160 * t186 - t169 * t226;
t223 = rSges(7,2) * t146 + t123 * t228 + t124 * t229;
t117 = pkin(4) * t145 + pkin(10) * t144;
t89 = rSges(6,1) * t122 - rSges(6,2) * t121 + rSges(6,3) * t144;
t222 = -t117 - t89;
t118 = pkin(4) * t147 + pkin(10) * t146;
t91 = rSges(6,1) * t124 - rSges(6,2) * t123 + rSges(6,3) * t146;
t221 = -t118 - t91;
t193 = t220 * t219;
t167 = t184 * t218 * t187 + (t187 * t193 + t217 * t227) * t182;
t174 = t184 * t220 - t198 * t219;
t162 = t167 * t226 + t174 * t186;
t166 = -t184 * t194 + (t187 * t217 - t193 * t227) * t182;
t148 = t162 * t185 - t166 * t225;
t149 = t162 * t225 + t166 * t185;
t161 = t167 * t186 - t174 * t226;
t216 = rSges(7,2) * t161 + t148 * t228 + t149 * t229;
t112 = rSges(6,1) * t149 - rSges(6,2) * t148 + rSges(6,3) * t161;
t142 = pkin(4) * t162 + pkin(10) * t161;
t215 = -t112 - t142;
t140 = pkin(3) * t158 + pkin(9) * t157;
t137 = t169 * t140;
t214 = t117 * t169 + t137;
t141 = pkin(3) * t160 + pkin(9) * t159;
t139 = t174 * t141;
t213 = t118 * t174 + t139;
t154 = pkin(3) * t167 + pkin(9) * t166;
t143 = t168 * t154;
t212 = t142 * t168 + t143;
t76 = Icges(7,5) * t122 + Icges(7,6) * t144 + Icges(7,3) * t121;
t80 = Icges(7,4) * t122 + Icges(7,2) * t144 + Icges(7,6) * t121;
t84 = Icges(7,1) * t122 + Icges(7,4) * t144 + Icges(7,5) * t121;
t32 = t121 * t76 + t122 * t84 + t144 * t80;
t77 = Icges(7,5) * t124 + Icges(7,6) * t146 + Icges(7,3) * t123;
t81 = Icges(7,4) * t124 + Icges(7,2) * t146 + Icges(7,6) * t123;
t85 = Icges(7,1) * t124 + Icges(7,4) * t146 + Icges(7,5) * t123;
t33 = t121 * t77 + t122 * t85 + t144 * t81;
t105 = Icges(7,5) * t149 + Icges(7,6) * t161 + Icges(7,3) * t148;
t107 = Icges(7,4) * t149 + Icges(7,2) * t161 + Icges(7,6) * t148;
t109 = Icges(7,1) * t149 + Icges(7,4) * t161 + Icges(7,5) * t148;
t50 = t105 * t121 + t107 * t144 + t109 * t122;
t1 = t144 * t32 + t146 * t33 + t161 * t50;
t78 = Icges(6,5) * t122 - Icges(6,6) * t121 + Icges(6,3) * t144;
t82 = Icges(6,4) * t122 - Icges(6,2) * t121 + Icges(6,6) * t144;
t86 = Icges(6,1) * t122 - Icges(6,4) * t121 + Icges(6,5) * t144;
t34 = -t121 * t82 + t122 * t86 + t144 * t78;
t79 = Icges(6,5) * t124 - Icges(6,6) * t123 + Icges(6,3) * t146;
t83 = Icges(6,4) * t124 - Icges(6,2) * t123 + Icges(6,6) * t146;
t87 = Icges(6,1) * t124 - Icges(6,4) * t123 + Icges(6,5) * t146;
t35 = -t121 * t83 + t122 * t87 + t144 * t79;
t106 = Icges(6,5) * t149 - Icges(6,6) * t148 + Icges(6,3) * t161;
t108 = Icges(6,4) * t149 - Icges(6,2) * t148 + Icges(6,6) * t161;
t110 = Icges(6,1) * t149 - Icges(6,4) * t148 + Icges(6,5) * t161;
t51 = t106 * t144 - t108 * t121 + t110 * t122;
t2 = t144 * t34 + t146 * t35 + t161 * t51;
t211 = t2 / 0.2e1 + t1 / 0.2e1;
t36 = t123 * t76 + t124 * t84 + t146 * t80;
t37 = t123 * t77 + t124 * t85 + t146 * t81;
t52 = t105 * t123 + t107 * t146 + t109 * t124;
t3 = t144 * t36 + t146 * t37 + t161 * t52;
t38 = -t123 * t82 + t124 * t86 + t146 * t78;
t39 = -t123 * t83 + t124 * t87 + t146 * t79;
t53 = t106 * t146 - t108 * t123 + t110 * t124;
t4 = t144 * t38 + t146 * t39 + t161 * t53;
t210 = t4 / 0.2e1 + t3 / 0.2e1;
t5 = t157 * t32 + t159 * t33 + t166 * t50;
t6 = t157 * t34 + t159 * t35 + t166 * t51;
t209 = t5 / 0.2e1 + t6 / 0.2e1;
t7 = t157 * t36 + t159 * t37 + t166 * t52;
t8 = t157 * t38 + t159 * t39 + t166 * t53;
t208 = t8 / 0.2e1 + t7 / 0.2e1;
t10 = t168 * t34 + t169 * t35 + t174 * t51;
t9 = t168 * t32 + t169 * t33 + t174 * t50;
t207 = t10 / 0.2e1 + t9 / 0.2e1;
t206 = -t117 - t224;
t205 = -t118 - t223;
t11 = t168 * t36 + t169 * t37 + t174 * t52;
t12 = t168 * t38 + t169 * t39 + t174 * t53;
t204 = t12 / 0.2e1 + t11 / 0.2e1;
t41 = t148 * t76 + t149 * t84 + t161 * t80;
t42 = t148 * t77 + t149 * t85 + t161 * t81;
t61 = t105 * t148 + t107 * t161 + t109 * t149;
t13 = t144 * t41 + t146 * t42 + t161 * t61;
t43 = -t148 * t82 + t149 * t86 + t161 * t78;
t44 = -t148 * t83 + t149 * t87 + t161 * t79;
t62 = t106 * t161 - t108 * t148 + t110 * t149;
t14 = t144 * t43 + t146 * t44 + t161 * t62;
t203 = t13 / 0.2e1 + t14 / 0.2e1;
t15 = t157 * t41 + t159 * t42 + t166 * t61;
t16 = t157 * t43 + t159 * t44 + t166 * t62;
t202 = t15 / 0.2e1 + t16 / 0.2e1;
t17 = t168 * t41 + t169 * t42 + t174 * t61;
t18 = t168 * t43 + t169 * t44 + t174 * t62;
t201 = t17 / 0.2e1 + t18 / 0.2e1;
t200 = -t142 - t216;
t195 = m(3) + m(4) + m(5) + m(6) + m(7);
t153 = rSges(4,1) * t167 - rSges(4,2) * t166 + rSges(4,3) * t174;
t152 = Icges(4,1) * t167 - Icges(4,4) * t166 + Icges(4,5) * t174;
t151 = Icges(4,4) * t167 - Icges(4,2) * t166 + Icges(4,6) * t174;
t150 = Icges(4,5) * t167 - Icges(4,6) * t166 + Icges(4,3) * t174;
t136 = rSges(5,1) * t162 - rSges(5,2) * t161 + rSges(5,3) * t166;
t135 = Icges(5,1) * t162 - Icges(5,4) * t161 + Icges(5,5) * t166;
t134 = Icges(5,4) * t162 - Icges(5,2) * t161 + Icges(5,6) * t166;
t133 = Icges(5,5) * t162 - Icges(5,6) * t161 + Icges(5,3) * t166;
t132 = rSges(4,1) * t160 - rSges(4,2) * t159 + rSges(4,3) * t169;
t131 = rSges(4,1) * t158 - rSges(4,2) * t157 + rSges(4,3) * t168;
t130 = Icges(4,1) * t160 - Icges(4,4) * t159 + Icges(4,5) * t169;
t129 = Icges(4,1) * t158 - Icges(4,4) * t157 + Icges(4,5) * t168;
t128 = Icges(4,4) * t160 - Icges(4,2) * t159 + Icges(4,6) * t169;
t127 = Icges(4,4) * t158 - Icges(4,2) * t157 + Icges(4,6) * t168;
t126 = Icges(4,5) * t160 - Icges(4,6) * t159 + Icges(4,3) * t169;
t125 = Icges(4,5) * t158 - Icges(4,6) * t157 + Icges(4,3) * t168;
t120 = t157 * t142;
t114 = t166 * t118;
t113 = t159 * t117;
t104 = rSges(5,1) * t147 - rSges(5,2) * t146 + rSges(5,3) * t159;
t103 = rSges(5,1) * t145 - rSges(5,2) * t144 + rSges(5,3) * t157;
t102 = Icges(5,1) * t147 - Icges(5,4) * t146 + Icges(5,5) * t159;
t101 = Icges(5,1) * t145 - Icges(5,4) * t144 + Icges(5,5) * t157;
t100 = Icges(5,4) * t147 - Icges(5,2) * t146 + Icges(5,6) * t159;
t99 = Icges(5,4) * t145 - Icges(5,2) * t144 + Icges(5,6) * t157;
t98 = Icges(5,5) * t147 - Icges(5,6) * t146 + Icges(5,3) * t159;
t97 = Icges(5,5) * t145 - Icges(5,6) * t144 + Icges(5,3) * t157;
t96 = t132 * t174 - t153 * t169;
t95 = -t131 * t174 + t153 * t168;
t92 = t131 * t169 - t132 * t168;
t75 = t104 * t166 - t136 * t159;
t74 = -t103 * t166 + t136 * t157;
t73 = t133 * t166 - t134 * t161 + t135 * t162;
t72 = t103 * t159 - t104 * t157;
t71 = t174 * t104 + t139 + (-t136 - t154) * t169;
t70 = t168 * t136 + t143 + (-t103 - t140) * t174;
t69 = t133 * t159 - t134 * t146 + t135 * t147;
t68 = t133 * t157 - t134 * t144 + t135 * t145;
t67 = -t112 * t146 + t161 * t91;
t66 = t112 * t144 - t161 * t89;
t65 = t169 * t103 + t137 + (-t104 - t141) * t168;
t64 = -t100 * t161 + t102 * t162 + t166 * t98;
t63 = t101 * t162 - t161 * t99 + t166 * t97;
t60 = -t100 * t146 + t102 * t147 + t159 * t98;
t59 = t101 * t147 - t146 * t99 + t159 * t97;
t58 = -t100 * t144 + t102 * t145 + t157 * t98;
t57 = t101 * t145 - t144 * t99 + t157 * t97;
t56 = -t144 * t91 + t146 * t89;
t55 = t159 * t215 + t166 * t91 + t114;
t54 = t157 * t112 + t166 * t222 + t120;
t49 = t174 * t91 + (-t154 + t215) * t169 + t213;
t48 = t168 * t112 + (-t140 + t222) * t174 + t212;
t47 = t157 * t221 + t159 * t89 + t113;
t46 = -t146 * t216 + t161 * t223;
t45 = t144 * t216 - t161 * t224;
t40 = t169 * t89 + (-t141 + t221) * t168 + t214;
t31 = t159 * t200 + t166 * t223 + t114;
t30 = t157 * t216 + t166 * t206 + t120;
t29 = -t144 * t223 + t146 * t224;
t28 = t223 * t174 + (-t154 + t200) * t169 + t213;
t27 = t216 * t168 + (-t140 + t206) * t174 + t212;
t26 = t157 * t205 + t159 * t224 + t113;
t25 = t224 * t169 + (-t141 + t205) * t168 + t214;
t24 = t168 * t63 + t169 * t64 + t174 * t73;
t23 = t157 * t63 + t159 * t64 + t166 * t73;
t22 = t168 * t59 + t169 * t60 + t174 * t69;
t21 = t168 * t57 + t169 * t58 + t174 * t68;
t20 = t157 * t59 + t159 * t60 + t166 * t69;
t19 = t157 * t57 + t159 * t58 + t166 * t68;
t88 = [m(2) + t195; t195 * t184; 0.2e1 * (m(7) / 0.2e1 + m(6) / 0.2e1 + m(5) / 0.2e1 + m(4) / 0.2e1 + m(3) / 0.2e1) * (t184 ^ 2 + (t181 ^ 2 + t183 ^ 2) * t182 ^ 2); m(4) * t92 + m(5) * t65 + m(6) * t40 + m(7) * t25; m(7) * (t25 * t184 + (t181 * t27 - t183 * t28) * t182) + m(6) * (t40 * t184 + (t181 * t48 - t183 * t49) * t182) + m(5) * (t65 * t184 + (t181 * t70 - t183 * t71) * t182) + m(4) * (t92 * t184 + (t181 * t95 - t183 * t96) * t182); (t17 + t18 + t24 + (t150 * t174 - t151 * t166 + t152 * t167) * t174) * t174 + (t11 + t12 + t22 + (t126 * t169 - t128 * t159 + t130 * t160) * t169 + (t126 * t174 - t128 * t166 + t130 * t167 + t150 * t169 - t151 * t159 + t152 * t160) * t174) * t169 + (t10 + t9 + t21 + (t125 * t168 - t127 * t157 + t129 * t158) * t168 + (t125 * t174 - t127 * t166 + t129 * t167 + t150 * t168 - t151 * t157 + t152 * t158) * t174 + (t125 * t169 + t126 * t168 - t127 * t159 - t128 * t157 + t129 * t160 + t130 * t158) * t169) * t168 + m(7) * (t25 ^ 2 + t27 ^ 2 + t28 ^ 2) + m(6) * (t40 ^ 2 + t48 ^ 2 + t49 ^ 2) + m(5) * (t65 ^ 2 + t70 ^ 2 + t71 ^ 2) + m(4) * (t92 ^ 2 + t95 ^ 2 + t96 ^ 2); m(5) * t72 + m(6) * t47 + m(7) * t26; m(7) * (t26 * t184 + (t181 * t30 - t183 * t31) * t182) + m(6) * (t47 * t184 + (t181 * t54 - t183 * t55) * t182) + m(5) * (t72 * t184 + (t181 * t74 - t183 * t75) * t182); (t23 / 0.2e1 + t202) * t174 + (t20 / 0.2e1 + t208) * t169 + (t19 / 0.2e1 + t209) * t168 + (t24 / 0.2e1 + t201) * t166 + (t22 / 0.2e1 + t204) * t159 + (t21 / 0.2e1 + t207) * t157 + m(7) * (t25 * t26 + t27 * t30 + t28 * t31) + m(6) * (t40 * t47 + t48 * t54 + t49 * t55) + m(5) * (t65 * t72 + t70 * t74 + t71 * t75); (t15 + t16 + t23) * t166 + (t8 + t7 + t20) * t159 + (t5 + t6 + t19) * t157 + m(7) * (t26 ^ 2 + t30 ^ 2 + t31 ^ 2) + m(6) * (t47 ^ 2 + t54 ^ 2 + t55 ^ 2) + m(5) * (t72 ^ 2 + t74 ^ 2 + t75 ^ 2); m(6) * t56 + m(7) * t29; m(7) * (t29 * t184 + (t181 * t45 - t183 * t46) * t182) + m(6) * (t56 * t184 + (t181 * t66 - t183 * t67) * t182); t203 * t174 + t210 * t169 + t211 * t168 + t201 * t161 + t204 * t146 + t207 * t144 + m(7) * (t25 * t29 + t27 * t45 + t28 * t46) + m(6) * (t40 * t56 + t48 * t66 + t49 * t67); t203 * t166 + t202 * t161 + t210 * t159 + t211 * t157 + t208 * t146 + t209 * t144 + m(7) * (t26 * t29 + t30 * t45 + t31 * t46) + m(6) * (t47 * t56 + t54 * t66 + t55 * t67); (t13 + t14) * t161 + (t3 + t4) * t146 + (t1 + t2) * t144 + m(7) * (t29 ^ 2 + t45 ^ 2 + t46 ^ 2) + m(6) * (t56 ^ 2 + t66 ^ 2 + t67 ^ 2); m(7) * t148; m(7) * (t148 * t184 + (-t121 * t183 + t123 * t181) * t182); m(7) * (t121 * t28 + t123 * t27 + t148 * t25); m(7) * (t121 * t31 + t123 * t30 + t148 * t26); m(7) * (t121 * t46 + t123 * t45 + t148 * t29); m(7) * (t121 ^ 2 + t123 ^ 2 + t148 ^ 2);];
%% Postprocessing: Reshape Output
% From vec2symmat_6_matlab.m
res = [t88(1) t88(2) t88(4) t88(7) t88(11) t88(16); t88(2) t88(3) t88(5) t88(8) t88(12) t88(17); t88(4) t88(5) t88(6) t88(9) t88(13) t88(18); t88(7) t88(8) t88(9) t88(10) t88(14) t88(19); t88(11) t88(12) t88(13) t88(14) t88(15) t88(20); t88(16) t88(17) t88(18) t88(19) t88(20) t88(21);];
Mq  = res;
