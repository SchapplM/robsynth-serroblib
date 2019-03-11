% Calculate joint inertia matrix for
% S6RPPRRR1
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% pkin [11x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d4,d5,d6,theta2,theta3]';
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
% Datum: 2019-03-09 02:19
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Mq = S6RPPRRR1_inertiaJ_slag_vp1(qJ, ...
  pkin, m, rSges, Icges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(11,1),zeros(7,1),zeros(7,3),zeros(7,6)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPPRRR1_inertiaJ_slag_vp1: qJ has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [11 1]), ...
  'S6RPPRRR1_inertiaJ_slag_vp1: pkin has to be [11x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6RPPRRR1_inertiaJ_slag_vp1: m has to be [7x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [7,3]), ...
  'S6RPPRRR1_inertiaJ_slag_vp1: rSges has to be [7x3] (double)');
assert(isreal(Icges) && all(size(Icges) == [7 6]), ...
  'S6RPPRRR1_inertiaJ_slag_vp1: Icges has to be [7x6] (double)'); 

%% Symbolic Calculation
% From inertia_joint_joint_floatb_twist_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 02:18:02
% EndTime: 2019-03-09 02:18:05
% DurationCPUTime: 1.39s
% Computational Cost: add. (6623->280), mult. (4343->413), div. (0->0), fcn. (4513->12), ass. (0->142)
t130 = qJ(1) + pkin(10);
t122 = sin(t130);
t119 = t122 ^ 2;
t129 = pkin(11) + qJ(4);
t125 = qJ(5) + t129;
t117 = cos(t125);
t124 = cos(t130);
t168 = t117 * t124;
t116 = sin(t125);
t169 = t116 * t124;
t134 = sin(qJ(6));
t165 = t124 * t134;
t136 = cos(qJ(6));
t166 = t122 * t136;
t89 = -t117 * t165 + t166;
t164 = t124 * t136;
t167 = t122 * t134;
t90 = t117 * t164 + t167;
t53 = t90 * rSges(7,1) + t89 * rSges(7,2) + rSges(7,3) * t169;
t190 = pkin(5) * t168 + pkin(9) * t169 + t53;
t120 = t124 ^ 2;
t171 = Icges(6,4) * t117;
t144 = -Icges(6,2) * t116 + t171;
t70 = Icges(6,6) * t122 + t144 * t124;
t172 = Icges(6,4) * t116;
t146 = Icges(6,1) * t117 - t172;
t72 = Icges(6,5) * t122 + t146 * t124;
t151 = -t116 * t70 + t117 * t72;
t69 = -Icges(6,6) * t124 + t144 * t122;
t71 = -Icges(6,5) * t124 + t146 * t122;
t152 = t116 * t69 - t117 * t71;
t142 = Icges(6,5) * t117 - Icges(6,6) * t116;
t67 = -Icges(6,3) * t124 + t142 * t122;
t68 = Icges(6,3) * t122 + t142 * t124;
t170 = t116 * t122;
t87 = -t117 * t167 - t164;
t88 = t117 * t166 - t165;
t46 = Icges(7,5) * t88 + Icges(7,6) * t87 + Icges(7,3) * t170;
t48 = Icges(7,4) * t88 + Icges(7,2) * t87 + Icges(7,6) * t170;
t50 = Icges(7,1) * t88 + Icges(7,4) * t87 + Icges(7,5) * t170;
t14 = t46 * t170 + t48 * t87 + t50 * t88;
t47 = Icges(7,5) * t90 + Icges(7,6) * t89 + Icges(7,3) * t169;
t49 = Icges(7,4) * t90 + Icges(7,2) * t89 + Icges(7,6) * t169;
t51 = Icges(7,1) * t90 + Icges(7,4) * t89 + Icges(7,5) * t169;
t15 = t47 * t170 + t49 * t87 + t51 * t88;
t8 = t122 * t15 - t124 * t14;
t189 = -t120 * t67 - (t151 * t122 + (t152 - t68) * t124) * t122 - t8;
t188 = t122 / 0.2e1;
t187 = -t124 / 0.2e1;
t16 = t46 * t169 + t48 * t89 + t50 * t90;
t17 = t47 * t169 + t49 * t89 + t51 * t90;
t9 = t122 * t17 - t124 * t16;
t186 = (t119 * t68 + t9 + (t152 * t124 + (t151 - t67) * t122) * t124) * t122;
t135 = sin(qJ(1));
t185 = pkin(1) * t135;
t121 = sin(t129);
t184 = pkin(4) * t121;
t183 = pkin(5) * t117;
t133 = -pkin(7) - qJ(3);
t132 = cos(pkin(11));
t118 = t132 * pkin(3) + pkin(2);
t123 = cos(t129);
t107 = pkin(4) * t123 + t118;
t93 = t124 * t107;
t182 = t124 * (-t118 * t124 + t93) + (t107 - t118) * t119;
t140 = rSges(6,1) * t168 - rSges(6,2) * t169 + t122 * rSges(6,3);
t153 = rSges(6,1) * t117 - rSges(6,2) * t116;
t38 = t122 * (-rSges(6,3) * t124 + t153 * t122) + t124 * t140;
t84 = -rSges(7,3) * t117 + (rSges(7,1) * t136 - rSges(7,2) * t134) * t116;
t181 = -pkin(5) * t116 + pkin(9) * t117 - t84;
t180 = rSges(5,1) * t123;
t179 = rSges(5,2) * t121;
t82 = -Icges(7,6) * t117 + (Icges(7,4) * t136 - Icges(7,2) * t134) * t116;
t178 = t134 * t82;
t20 = -t117 * t46 + (-t134 * t48 + t136 * t50) * t116;
t177 = t20 * t124;
t21 = -t117 * t47 + (-t134 * t49 + t136 * t51) * t116;
t176 = t21 * t122;
t175 = rSges(4,3) + qJ(3);
t174 = Icges(5,4) * t121;
t173 = Icges(5,4) * t123;
t162 = t122 * rSges(5,3) + t124 * t180;
t161 = t119 + t120;
t98 = rSges(6,1) * t116 + rSges(6,2) * t117;
t159 = -t98 - t184;
t155 = -rSges(7,1) * t88 - rSges(7,2) * t87;
t52 = rSges(7,3) * t170 - t155;
t22 = t122 * t52 + t119 * (pkin(9) * t116 + t183) + t190 * t124;
t81 = -Icges(7,3) * t117 + (Icges(7,5) * t136 - Icges(7,6) * t134) * t116;
t83 = -Icges(7,5) * t117 + (Icges(7,1) * t136 - Icges(7,4) * t134) * t116;
t27 = t81 * t170 + t82 * t87 + t83 * t88;
t3 = -t117 * t27 + (t122 * t14 + t124 * t15) * t116;
t28 = t81 * t169 + t82 * t89 + t83 * t90;
t4 = -t117 * t28 + (t122 * t16 + t124 * t17) * t116;
t158 = t3 * t187 + t4 * t188 - t117 * (t176 - t177) / 0.2e1 + t8 * t170 / 0.2e1 + t9 * t169 / 0.2e1;
t157 = t181 - t184;
t137 = cos(qJ(1));
t127 = t137 * pkin(1);
t128 = -pkin(8) + t133;
t156 = -t122 * t128 + t127 + t93;
t154 = -t179 + t180;
t96 = Icges(6,2) * t117 + t172;
t97 = Icges(6,1) * t116 + t171;
t150 = -t116 * t96 + t117 * t97;
t147 = Icges(5,1) * t123 - t174;
t145 = -Icges(5,2) * t121 + t173;
t143 = Icges(5,5) * t123 - Icges(5,6) * t121;
t141 = t189 * t124 + t186;
t131 = sin(pkin(11));
t139 = rSges(4,1) * t132 - rSges(4,2) * t131 + pkin(2);
t95 = Icges(6,5) * t116 + Icges(6,6) * t117;
t138 = t176 / 0.2e1 - t177 / 0.2e1 + (t116 * t72 + t117 * t70 + t122 * t95 + t150 * t124 + t28) * t188 + (t116 * t71 + t117 * t69 + t150 * t122 - t124 * t95 + t27) * t187;
t113 = rSges(2,1) * t137 - t135 * rSges(2,2);
t112 = -t135 * rSges(2,1) - rSges(2,2) * t137;
t104 = rSges(5,1) * t121 + rSges(5,2) * t123;
t92 = rSges(3,1) * t124 - rSges(3,2) * t122 + t127;
t91 = -rSges(3,1) * t122 - rSges(3,2) * t124 - t185;
t76 = Icges(5,3) * t122 + t143 * t124;
t75 = -Icges(5,3) * t124 + t143 * t122;
t66 = t159 * t124;
t65 = t159 * t122;
t62 = t116 * t136 * t83;
t61 = t175 * t122 + t139 * t124 + t127;
t60 = -t139 * t122 + t175 * t124 - t185;
t57 = -t122 * t133 + t127 + (t118 - t179) * t124 + t162;
t56 = -t185 + (rSges(5,3) - t133) * t124 + (-t118 - t154) * t122;
t55 = t181 * t124;
t54 = t181 * t122;
t45 = t140 + t156;
t44 = -t185 + (rSges(6,3) - t128) * t124 + (-t107 - t153) * t122;
t43 = t157 * t124;
t42 = t157 * t122;
t39 = t124 * (-t124 * t179 + t162) + (-t124 * rSges(5,3) + t154 * t122) * t122;
t33 = -t117 * t53 - t84 * t169;
t32 = t117 * t52 + t84 * t170;
t31 = -t116 * t178 - t117 * t81 + t62;
t30 = t156 + t190;
t29 = -t185 - t124 * t128 + (-t183 - t107 + (-rSges(7,3) - pkin(9)) * t116) * t122 + t155;
t26 = (-t122 * t53 + t124 * t52) * t116;
t23 = t38 + t182;
t11 = t22 + t182;
t1 = [Icges(4,2) * t132 ^ 2 + t123 * (Icges(5,2) * t123 + t174) + t121 * (Icges(5,1) * t121 + t173) + Icges(2,3) + Icges(3,3) + t62 + (Icges(4,1) * t131 + 0.2e1 * Icges(4,4) * t132) * t131 + (-t81 + t96) * t117 + (t97 - t178) * t116 + m(7) * (t29 ^ 2 + t30 ^ 2) + m(6) * (t44 ^ 2 + t45 ^ 2) + m(5) * (t56 ^ 2 + t57 ^ 2) + m(4) * (t60 ^ 2 + t61 ^ 2) + m(3) * (t91 ^ 2 + t92 ^ 2) + m(2) * (t112 ^ 2 + t113 ^ 2); 0; m(3) + m(4) + m(5) + m(6) + m(7); m(7) * (t122 * t29 - t124 * t30) + m(6) * (t122 * t44 - t124 * t45) + m(5) * (t122 * t56 - t124 * t57) + m(4) * (t122 * t60 - t124 * t61); 0; 0.2e1 * (m(4) / 0.2e1 + m(5) / 0.2e1 + m(6) / 0.2e1 + m(7) / 0.2e1) * t161; (t121 * (-Icges(5,5) * t124 + t147 * t122) + t123 * (-Icges(5,6) * t124 + t145 * t122)) * t187 + (t121 * (Icges(5,5) * t122 + t147 * t124) + t123 * (Icges(5,6) * t122 + t145 * t124)) * t188 + m(7) * (t29 * t43 + t30 * t42) + m(6) * (t44 * t66 + t45 * t65) + m(5) * (-t122 * t57 - t124 * t56) * t104 + (t120 / 0.2e1 + t119 / 0.2e1) * (Icges(5,5) * t121 + Icges(5,6) * t123) + t138; m(5) * t39 + m(6) * t23 + m(7) * t11; m(6) * (t122 * t66 - t124 * t65) + m(7) * (t122 * t43 - t124 * t42); m(7) * (t11 ^ 2 + t42 ^ 2 + t43 ^ 2) + m(6) * (t23 ^ 2 + t65 ^ 2 + t66 ^ 2) + t122 * t119 * t76 + m(5) * (t161 * t104 ^ 2 + t39 ^ 2) + t186 + (-t120 * t75 + (-t122 * t75 + t124 * t76) * t122 + t189) * t124; m(7) * (t29 * t55 + t30 * t54) + m(6) * (-t122 * t45 - t124 * t44) * t98 + t138; m(6) * t38 + m(7) * t22; m(7) * (t122 * t55 - t124 * t54); m(7) * (t11 * t22 + t42 * t54 + t43 * t55) + m(6) * (t23 * t38 + (-t122 * t65 - t124 * t66) * t98) + t141; m(6) * (t161 * t98 ^ 2 + t38 ^ 2) + m(7) * (t22 ^ 2 + t54 ^ 2 + t55 ^ 2) + t141; -t31 * t117 + m(7) * (t29 * t32 + t30 * t33) + ((t21 / 0.2e1 + t28 / 0.2e1) * t124 + (t27 / 0.2e1 + t20 / 0.2e1) * t122) * t116; m(7) * t26; m(7) * (t122 * t32 - t124 * t33); m(7) * (t11 * t26 + t32 * t43 + t33 * t42) + t158; m(7) * (t22 * t26 + t32 * t55 + t33 * t54) + t158; t117 ^ 2 * t31 + m(7) * (t26 ^ 2 + t32 ^ 2 + t33 ^ 2) + (t124 * t4 + t122 * t3 - t117 * (t122 * t20 + t124 * t21)) * t116;];
%% Postprocessing: Reshape Output
% From vec2symmat_6_matlab.m
res = [t1(1) t1(2) t1(4) t1(7) t1(11) t1(16); t1(2) t1(3) t1(5) t1(8) t1(12) t1(17); t1(4) t1(5) t1(6) t1(9) t1(13) t1(18); t1(7) t1(8) t1(9) t1(10) t1(14) t1(19); t1(11) t1(12) t1(13) t1(14) t1(15) t1(20); t1(16) t1(17) t1(18) t1(19) t1(20) t1(21);];
Mq  = res;
