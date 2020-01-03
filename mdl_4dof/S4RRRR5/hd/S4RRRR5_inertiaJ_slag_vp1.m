% Calculate joint inertia matrix for
% S4RRRR5
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [4x1]
%   Generalized joint coordinates (joint angles)
% pkin [7x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,d1,d2,d3,d4]';
% m_mdh [5x1]
%   mass of all robot links (including the base)
% rSges [5x3]
%   center of mass of all robot links (in body frames)
%   rows: links of the robot (starting with base)
%   columns: x-, y-, z-coordinates
% Icges [5x6]
%   inertia of all robot links about their respective center of mass, in body frames
%   rows: links of the robot (starting with base)
%   columns: xx, yy, zz, xy, xz, yz (see inertiavector2matrix.m)
% 
% Output:
% Mq [4x4]
%   inertia matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-31 17:28
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Mq = S4RRRR5_inertiaJ_slag_vp1(qJ, ...
  pkin, m, rSges, Icges)
%% Coder Information
%#codegen
%$cgargs {zeros(4,1),zeros(7,1),zeros(5,1),zeros(5,3),zeros(5,6)}
assert(isreal(qJ) && all(size(qJ) == [4 1]), ...
  'S4RRRR5_inertiaJ_slag_vp1: qJ has to be [4x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [7 1]), ...
  'S4RRRR5_inertiaJ_slag_vp1: pkin has to be [7x1] (double)');
assert(isreal(m) && all(size(m) == [5 1]), ...
  'S4RRRR5_inertiaJ_slag_vp1: m has to be [5x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [5,3]), ...
  'S4RRRR5_inertiaJ_slag_vp1: rSges has to be [5x3] (double)');
assert(isreal(Icges) && all(size(Icges) == [5 6]), ...
  'S4RRRR5_inertiaJ_slag_vp1: Icges has to be [5x6] (double)'); 

%% Symbolic Calculation
% From inertia_joint_joint_floatb_twist_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 17:27:23
% EndTime: 2019-12-31 17:27:27
% DurationCPUTime: 1.57s
% Computational Cost: add. (3532->281), mult. (5587->412), div. (0->0), fcn. (6106->8), ass. (0->152)
t137 = sin(qJ(2));
t193 = Icges(3,5) * t137;
t192 = t193 / 0.2e1;
t139 = cos(qJ(3));
t128 = pkin(3) * t139 + pkin(2);
t142 = -pkin(7) - pkin(6);
t140 = cos(qJ(2));
t141 = cos(qJ(1));
t166 = t140 * t141;
t168 = t137 * t141;
t136 = sin(qJ(3));
t138 = sin(qJ(1));
t171 = t136 * t138;
t135 = qJ(3) + qJ(4);
t129 = sin(t135);
t130 = cos(t135);
t103 = -t129 * t166 + t130 * t138;
t104 = t129 * t138 + t130 * t166;
t65 = t104 * rSges(5,1) + t103 * rSges(5,2) + rSges(5,3) * t168;
t191 = pkin(3) * t171 + t128 * t166 - t142 * t168 + t65;
t190 = t138 ^ 2;
t189 = t141 ^ 2;
t169 = t137 * t138;
t167 = t138 * t140;
t101 = -t129 * t167 - t130 * t141;
t102 = -t129 * t141 + t130 * t167;
t58 = Icges(5,5) * t102 + Icges(5,6) * t101 + Icges(5,3) * t169;
t60 = Icges(5,4) * t102 + Icges(5,2) * t101 + Icges(5,6) * t169;
t62 = Icges(5,1) * t102 + Icges(5,4) * t101 + Icges(5,5) * t169;
t19 = t101 * t60 + t102 * t62 + t58 * t169;
t59 = Icges(5,5) * t104 + Icges(5,6) * t103 + Icges(5,3) * t168;
t61 = Icges(5,4) * t104 + Icges(5,2) * t103 + Icges(5,6) * t168;
t63 = Icges(5,1) * t104 + Icges(5,4) * t103 + Icges(5,5) * t168;
t20 = t101 * t61 + t102 * t63 + t59 * t169;
t85 = -Icges(5,3) * t140 + (Icges(5,5) * t130 - Icges(5,6) * t129) * t137;
t86 = -Icges(5,6) * t140 + (Icges(5,4) * t130 - Icges(5,2) * t129) * t137;
t87 = -Icges(5,5) * t140 + (Icges(5,1) * t130 - Icges(5,4) * t129) * t137;
t38 = t101 * t86 + t102 * t87 + t85 * t169;
t5 = -t140 * t38 + (t138 * t19 + t141 * t20) * t137;
t21 = t103 * t60 + t104 * t62 + t58 * t168;
t22 = t103 * t61 + t104 * t63 + t59 * t168;
t39 = t103 * t86 + t104 * t87 + t85 * t168;
t6 = -t140 * t39 + (t138 * t21 + t141 * t22) * t137;
t188 = t6 * t168 + t5 * t169;
t187 = t138 / 0.2e1;
t186 = -t140 / 0.2e1;
t185 = -t141 / 0.2e1;
t184 = t141 / 0.2e1;
t118 = rSges(3,1) * t137 + rSges(3,2) * t140;
t183 = m(3) * t118;
t182 = pkin(2) * t140;
t181 = -pkin(2) + t128;
t180 = pkin(6) + t142;
t163 = pkin(2) * t166 + pkin(6) * t168;
t179 = -t163 + t191;
t150 = -t102 * rSges(5,1) - t101 * rSges(5,2);
t64 = rSges(5,3) * t169 - t150;
t88 = -rSges(5,3) * t140 + (rSges(5,1) * t130 - rSges(5,2) * t129) * t137;
t47 = t140 * t64 + t88 * t169;
t84 = t181 * t137 + t180 * t140;
t178 = -t84 - t88;
t177 = t129 * t86;
t94 = -Icges(4,6) * t140 + (Icges(4,4) * t139 - Icges(4,2) * t136) * t137;
t176 = t136 * t94;
t175 = t141 * rSges(3,3);
t79 = t137 * t130 * t87;
t44 = -t137 * t177 - t140 * t85 + t79;
t174 = t44 * t140;
t172 = Icges(3,4) * t140;
t170 = t136 * t141;
t100 = -rSges(4,3) * t140 + (rSges(4,1) * t139 - rSges(4,2) * t136) * t137;
t121 = t137 * pkin(2) - t140 * pkin(6);
t165 = -t100 - t121;
t164 = t190 * (pkin(6) * t137 + t182) + t141 * t163;
t162 = t141 * pkin(1) + t138 * pkin(5);
t161 = pkin(3) * t170;
t160 = -t121 + t178;
t111 = -t136 * t166 + t138 * t139;
t112 = t139 * t166 + t171;
t68 = Icges(4,5) * t112 + Icges(4,6) * t111 + Icges(4,3) * t168;
t70 = Icges(4,4) * t112 + Icges(4,2) * t111 + Icges(4,6) * t168;
t72 = Icges(4,1) * t112 + Icges(4,4) * t111 + Icges(4,5) * t168;
t34 = -t140 * t68 + (-t136 * t70 + t139 * t72) * t137;
t91 = -Icges(4,3) * t140 + (Icges(4,5) * t139 - Icges(4,6) * t136) * t137;
t97 = -Icges(4,5) * t140 + (Icges(4,1) * t139 - Icges(4,4) * t136) * t137;
t42 = t111 * t94 + t112 * t97 + t91 * t168;
t159 = t34 / 0.2e1 + t42 / 0.2e1;
t109 = -t136 * t167 - t139 * t141;
t110 = t139 * t167 - t170;
t67 = Icges(4,5) * t110 + Icges(4,6) * t109 + Icges(4,3) * t169;
t69 = Icges(4,4) * t110 + Icges(4,2) * t109 + Icges(4,6) * t169;
t71 = Icges(4,1) * t110 + Icges(4,4) * t109 + Icges(4,5) * t169;
t33 = -t140 * t67 + (-t136 * t69 + t139 * t71) * t137;
t41 = t109 * t94 + t110 * t97 + t91 * t169;
t158 = t41 / 0.2e1 + t33 / 0.2e1;
t74 = t112 * rSges(4,1) + t111 * rSges(4,2) + rSges(4,3) * t168;
t157 = t169 / 0.2e1;
t156 = t168 / 0.2e1;
t25 = -t140 * t58 + (-t129 * t60 + t130 * t62) * t137;
t26 = -t140 * t59 + (-t129 * t61 + t130 * t63) * t137;
t155 = (t25 + t38) * t157 + (t26 + t39) * t156;
t7 = -t174 + (t138 * t25 + t141 * t26) * t137;
t154 = -t140 * t7 + t188;
t12 = t138 * t20 - t141 * t19;
t13 = t138 * t22 - t141 * t21;
t153 = t12 * t157 + t13 * t156 + t5 * t185 + t6 * t187 + (t26 * t138 - t25 * t141) * t186;
t152 = rSges(3,1) * t140 - rSges(3,2) * t137;
t151 = -rSges(4,1) * t110 - rSges(4,2) * t109;
t146 = -Icges(3,2) * t137 + t172;
t145 = Icges(3,5) * t140 - Icges(3,6) * t137;
t144 = rSges(3,1) * t166 - rSges(3,2) * t168 + t138 * rSges(3,3);
t133 = t141 * pkin(5);
t120 = rSges(2,1) * t141 - rSges(2,2) * t138;
t119 = -rSges(2,1) * t138 - rSges(2,2) * t141;
t115 = Icges(3,6) * t140 + t193;
t93 = Icges(3,3) * t138 + t145 * t141;
t92 = -Icges(3,3) * t141 + t145 * t138;
t83 = t137 * t139 * t97;
t82 = t144 + t162;
t81 = t175 + t133 + (-pkin(1) - t152) * t138;
t78 = t165 * t141;
t77 = t165 * t138;
t75 = -t161 + (-t180 * t137 + t181 * t140) * t138;
t73 = rSges(4,3) * t169 - t151;
t66 = t141 * t144 + (t152 * t138 - t175) * t138;
t56 = t64 * t168;
t55 = t74 + t162 + t163;
t54 = t133 + (-t182 - pkin(1) + (-rSges(4,3) - pkin(6)) * t137) * t138 + t151;
t53 = t160 * t141;
t52 = t160 * t138;
t51 = -t100 * t168 - t140 * t74;
t50 = t100 * t169 + t140 * t73;
t49 = -t137 * t176 - t140 * t91 + t83;
t48 = -t140 * t65 - t88 * t168;
t46 = t162 + t191;
t45 = t161 + t133 + (-t128 * t140 - pkin(1) + (-rSges(5,3) + t142) * t137) * t138 + t150;
t43 = (-t138 * t74 + t141 * t73) * t137;
t40 = -t65 * t169 + t56;
t35 = t138 * t73 + t141 * t74 + t164;
t32 = -t179 * t140 + t178 * t168;
t31 = t140 * t75 + t84 * t169 + t47;
t30 = t111 * t70 + t112 * t72 + t68 * t168;
t29 = t111 * t69 + t112 * t71 + t67 * t168;
t28 = t109 * t70 + t110 * t72 + t68 * t169;
t27 = t109 * t69 + t110 * t71 + t67 * t169;
t18 = t56 + (-t179 * t138 + t141 * t75) * t137;
t17 = t179 * t141 + (t64 + t75) * t138 + t164;
t16 = t138 * t30 - t141 * t29;
t15 = t138 * t28 - t141 * t27;
t9 = -t140 * t42 + (t138 * t29 + t141 * t30) * t137;
t8 = -t140 * t41 + (t138 * t27 + t141 * t28) * t137;
t1 = [Icges(2,3) + t79 + t83 + (Icges(3,4) * t137 + Icges(3,2) * t140 - t85 - t91) * t140 + (Icges(3,1) * t137 + t172 - t176 - t177) * t137 + m(5) * (t45 ^ 2 + t46 ^ 2) + m(4) * (t54 ^ 2 + t55 ^ 2) + m(3) * (t81 ^ 2 + t82 ^ 2) + m(2) * (t119 ^ 2 + t120 ^ 2); m(5) * (t45 * t53 + t46 * t52) + m(4) * (t54 * t78 + t55 * t77) + (-t38 / 0.2e1 - t25 / 0.2e1 + t141 * t192 + (-Icges(3,6) * t141 + t146 * t138) * t186 - t81 * t183 + t115 * t184 - t158) * t141 + (t39 / 0.2e1 + t26 / 0.2e1 + t138 * t192 + t140 * (Icges(3,6) * t138 + t146 * t141) / 0.2e1 - t82 * t183 + t115 * t187 + t159) * t138; m(5) * (t17 ^ 2 + t52 ^ 2 + t53 ^ 2) + m(4) * (t35 ^ 2 + t77 ^ 2 + t78 ^ 2) + m(3) * (t66 ^ 2 + (t189 + t190) * t118 ^ 2) + (t190 * t93 + t13 + t16) * t138 + (-t189 * t92 - t12 - t15 + (-t138 * t92 + t141 * t93) * t138) * t141; (-t44 - t49) * t140 + m(5) * (t31 * t45 + t32 * t46) + m(4) * (t50 * t54 + t51 * t55) + (t158 * t138 + t159 * t141) * t137 + t155; t9 * t187 + t8 * t185 + (t34 * t138 - t33 * t141) * t186 + (t15 * t187 + t16 * t184) * t137 + m(5) * (t17 * t18 + t31 * t53 + t32 * t52) + m(4) * (t43 * t35 + t50 * t78 + t51 * t77) + t153; (t49 * t140 - t7) * t140 + m(5) * (t18 ^ 2 + t31 ^ 2 + t32 ^ 2) + m(4) * (t43 ^ 2 + t50 ^ 2 + t51 ^ 2) + (t141 * t9 + t138 * t8 - t140 * (t138 * t33 + t141 * t34)) * t137 + t188; m(5) * (t45 * t47 + t46 * t48) - t174 + t155; m(5) * (t17 * t40 + t47 * t53 + t48 * t52) + t153; m(5) * (t18 * t40 + t31 * t47 + t32 * t48) + t154; m(5) * (t40 ^ 2 + t47 ^ 2 + t48 ^ 2) + t154;];
%% Postprocessing: Reshape Output
% From vec2symmat_4_matlab.m
res = [t1(1), t1(2), t1(4), t1(7); t1(2), t1(3), t1(5), t1(8); t1(4), t1(5), t1(6), t1(9); t1(7), t1(8), t1(9), t1(10);];
Mq = res;
