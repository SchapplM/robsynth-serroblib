% Calculate joint inertia matrix for
% S5RRRRR6
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% pkin [9x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d1,d2,d3,d4,d5]';
% m [6x1]
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
% Datum: 2022-01-20 12:09
% Revision: fd3771346c4aea32fdeb66112c511235427c26a7 (2022-01-20)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Mq = S5RRRRR6_inertiaJ_slag_vp1(qJ, ...
  pkin, m, rSges, Icges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(9,1),zeros(6,1),zeros(6,3),zeros(6,6)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RRRRR6_inertiaJ_slag_vp1: qJ has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S5RRRRR6_inertiaJ_slag_vp1: pkin has to be [9x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5RRRRR6_inertiaJ_slag_vp1: m has to be [6x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [6,3]), ...
  'S5RRRRR6_inertiaJ_slag_vp1: rSges has to be [6x3] (double)');
assert(isreal(Icges) && all(size(Icges) == [6 6]), ...
  'S5RRRRR6_inertiaJ_slag_vp1: Icges has to be [6x6] (double)'); 

%% Symbolic Calculation
% From inertia_joint_joint_floatb_twist_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2022-01-20 12:07:38
% EndTime: 2022-01-20 12:07:39
% DurationCPUTime: 0.97s
% Computational Cost: add. (4193->241), mult. (2794->340), div. (0->0), fcn. (2510->10), ass. (0->138)
t128 = qJ(1) + qJ(2);
t120 = sin(t128);
t117 = t120 ^ 2;
t122 = cos(t128);
t118 = t122 ^ 2;
t133 = -pkin(8) - pkin(7);
t129 = sin(qJ(3));
t131 = cos(qJ(3));
t98 = t129 * rSges(4,1) + t131 * rSges(4,2);
t191 = m(4) * t98;
t127 = qJ(3) + qJ(4);
t119 = sin(t127);
t121 = cos(t127);
t83 = t119 * rSges(5,1) + t121 * rSges(5,2);
t190 = m(5) * t83;
t123 = qJ(5) + t127;
t114 = sin(t123);
t115 = cos(t123);
t78 = t114 * rSges(6,1) + t115 * rSges(6,2);
t189 = m(6) * t78;
t188 = t120 / 0.2e1;
t187 = -t122 / 0.2e1;
t186 = pkin(3) * t129;
t130 = sin(qJ(1));
t185 = t130 * pkin(1);
t116 = t131 * pkin(3) + pkin(2);
t112 = t122 * pkin(7);
t160 = -t122 * t116 + t120 * t133;
t166 = -t122 * pkin(2) - t120 * pkin(7);
t168 = t122 * t133;
t184 = t120 * (t168 + t112 + (-pkin(2) + t116) * t120) + t122 * (-t160 + t166);
t178 = rSges(6,2) * t114;
t181 = rSges(6,1) * t115;
t139 = t120 * rSges(6,3) + (-t178 + t181) * t122;
t176 = t122 * rSges(6,3) + t120 * t178;
t14 = t120 * (t120 * t181 - t176) + t122 * t139;
t179 = rSges(5,2) * t119;
t182 = rSges(5,1) * t121;
t140 = t120 * rSges(5,3) + (-t179 + t182) * t122;
t175 = t122 * rSges(5,3) + t120 * t179;
t17 = t120 * (t120 * t182 - t175) + t122 * t140;
t183 = rSges(4,1) * t131;
t180 = rSges(4,2) * t129;
t126 = pkin(9) - t133;
t90 = pkin(4) * t121 + t116;
t177 = t126 * t120 + t122 * t90;
t174 = Icges(4,4) * t129;
t173 = Icges(4,4) * t131;
t172 = Icges(5,4) * t119;
t171 = Icges(5,4) * t121;
t170 = Icges(6,4) * t114;
t169 = Icges(6,4) * t115;
t167 = t122 * rSges(4,3) + t120 * t180;
t165 = t117 + t118;
t144 = -Icges(6,2) * t114 + t169;
t50 = Icges(6,6) * t120 + t144 * t122;
t147 = Icges(6,1) * t115 - t170;
t52 = Icges(6,5) * t120 + t147 * t122;
t158 = -t114 * t50 + t115 * t52;
t49 = -Icges(6,6) * t122 + t144 * t120;
t51 = -Icges(6,5) * t122 + t147 * t120;
t159 = t114 * t49 - t115 * t51;
t141 = Icges(6,5) * t115 - Icges(6,6) * t114;
t47 = -Icges(6,3) * t122 + t141 * t120;
t48 = Icges(6,3) * t120 + t141 * t122;
t1 = t120 * (t117 * t48 + (t159 * t122 + (t158 - t47) * t120) * t122);
t2 = t118 * t47 + (t158 * t120 + (t159 - t48) * t122) * t120;
t164 = -t122 * t2 + t1;
t76 = Icges(6,2) * t115 + t170;
t77 = Icges(6,1) * t114 + t169;
t157 = -t114 * t76 + t115 * t77;
t75 = Icges(6,5) * t114 + Icges(6,6) * t115;
t163 = (t114 * t52 + t115 * t50 + t120 * t75 + t157 * t122) * t188 + (t114 * t51 + t115 * t49 + t157 * t120 - t122 * t75) * t187;
t162 = -t83 - t186;
t161 = -pkin(4) * t119 - t78;
t6 = t120 * ((-t126 - t133) * t122 + (-t116 + t90) * t120) + t122 * (t160 + t177) + t14;
t85 = t122 * rSges(3,1) - t120 * rSges(3,2);
t84 = -t120 * rSges(3,1) - t122 * rSges(3,2);
t145 = -Icges(5,2) * t119 + t171;
t55 = -Icges(5,6) * t122 + t145 * t120;
t148 = Icges(5,1) * t121 - t172;
t57 = -Icges(5,5) * t122 + t148 * t120;
t156 = t119 * t55 - t121 * t57;
t56 = Icges(5,6) * t120 + t145 * t122;
t58 = Icges(5,5) * t120 + t148 * t122;
t155 = -t119 * t56 + t121 * t58;
t81 = Icges(5,2) * t121 + t172;
t82 = Icges(5,1) * t119 + t171;
t154 = -t119 * t81 + t121 * t82;
t96 = Icges(4,2) * t131 + t174;
t97 = Icges(4,1) * t129 + t173;
t151 = -t129 * t96 + t131 * t97;
t142 = Icges(5,5) * t121 - Icges(5,6) * t119;
t53 = -Icges(5,3) * t122 + t142 * t120;
t54 = Icges(5,3) * t120 + t142 * t122;
t4 = t120 * (t117 * t54 + (t156 * t122 + (t155 - t53) * t120) * t122);
t5 = t118 * t53 + (t155 * t120 + (t156 - t54) * t122) * t120;
t150 = t1 + t4 + (-t5 - t2) * t122;
t149 = Icges(4,1) * t131 - t174;
t146 = -Icges(4,2) * t129 + t173;
t143 = Icges(4,5) * t131 - Icges(4,6) * t129;
t138 = t120 * rSges(4,3) + (-t180 + t183) * t122;
t137 = t161 - t186;
t136 = t114 * t77 + t115 * t76 + t119 * t82 + t121 * t81 + t129 * t97 + t131 * t96 + Icges(3,3);
t80 = Icges(5,5) * t119 + Icges(5,6) * t121;
t135 = t163 + (t119 * t58 + t120 * t80 + t121 * t56 + t154 * t122) * t188 + (t119 * t57 + t154 * t120 + t121 * t55 - t122 * t80) * t187;
t26 = t139 + t177;
t38 = t138 - t166;
t37 = t112 + (-pkin(2) - t183) * t120 + t167;
t32 = t140 - t160;
t25 = t122 * t126 + (-t90 - t181) * t120 + t176;
t31 = -t168 + (-t116 - t182) * t120 + t175;
t95 = Icges(4,5) * t129 + Icges(4,6) * t131;
t134 = t135 + (t120 * t95 + t151 * t122 + t129 * (Icges(4,5) * t120 + t149 * t122) + t131 * (Icges(4,6) * t120 + t146 * t122)) * t188 + (t151 * t120 - t122 * t95 + t129 * (-Icges(4,5) * t122 + t149 * t120) + t131 * (-Icges(4,6) * t122 + t146 * t120)) * t187;
t132 = cos(qJ(1));
t125 = t132 * pkin(1);
t100 = t132 * rSges(2,1) - t130 * rSges(2,2);
t99 = -t130 * rSges(2,1) - t132 * rSges(2,2);
t72 = t125 + t85;
t71 = t84 - t185;
t62 = Icges(4,3) * t120 + t143 * t122;
t61 = -Icges(4,3) * t122 + t143 * t120;
t60 = t162 * t122;
t59 = t162 * t120;
t46 = t161 * t122;
t45 = t161 * t120;
t36 = t137 * t122;
t35 = t137 * t120;
t34 = t125 + t38;
t33 = t37 - t185;
t28 = t125 + t32;
t27 = t31 - t185;
t22 = t120 * (t120 * t183 - t167) + t122 * t138;
t21 = t125 + t26;
t20 = t25 - t185;
t7 = t17 + t184;
t3 = t6 + t184;
t8 = [Icges(2,3) + m(6) * (t20 ^ 2 + t21 ^ 2) + m(5) * (t27 ^ 2 + t28 ^ 2) + m(4) * (t33 ^ 2 + t34 ^ 2) + m(3) * (t71 ^ 2 + t72 ^ 2) + m(2) * (t100 ^ 2 + t99 ^ 2) + t136; m(6) * (t25 * t20 + t26 * t21) + m(5) * (t31 * t27 + t32 * t28) + m(4) * (t37 * t33 + t38 * t34) + m(3) * (t84 * t71 + t85 * t72) + t136; m(6) * (t25 ^ 2 + t26 ^ 2) + m(5) * (t31 ^ 2 + t32 ^ 2) + m(4) * (t37 ^ 2 + t38 ^ 2) + m(3) * (t84 ^ 2 + t85 ^ 2) + t136; m(6) * (t36 * t20 + t35 * t21) + m(5) * (t60 * t27 + t59 * t28) + (-t120 * t34 - t122 * t33) * t191 + t134; m(6) * (t36 * t25 + t35 * t26) + m(5) * (t60 * t31 + t59 * t32) + (-t120 * t38 - t122 * t37) * t191 + t134; m(6) * (t3 ^ 2 + t35 ^ 2 + t36 ^ 2) + t4 + m(5) * (t59 ^ 2 + t60 ^ 2 + t7 ^ 2) + m(4) * (t165 * t98 ^ 2 + t22 ^ 2) + t120 * t117 * t62 + t164 + (-t118 * t61 - t5 + (-t120 * t61 + t122 * t62) * t120) * t122; m(6) * (t46 * t20 + t45 * t21) + (-t120 * t28 - t122 * t27) * t190 + t135; m(6) * (t46 * t25 + t45 * t26) + (-t120 * t32 - t122 * t31) * t190 + t135; m(6) * (t6 * t3 + t45 * t35 + t46 * t36) + m(5) * (t17 * t7 + (-t120 * t59 - t122 * t60) * t83) + t150; m(5) * (t165 * t83 ^ 2 + t17 ^ 2) + m(6) * (t45 ^ 2 + t46 ^ 2 + t6 ^ 2) + t150; (-t120 * t21 - t122 * t20) * t189 + t163; (-t120 * t26 - t122 * t25) * t189 + t163; m(6) * (t14 * t3 + (-t120 * t35 - t122 * t36) * t78) + t164; m(6) * (t14 * t6 + (-t120 * t45 - t122 * t46) * t78) + t164; m(6) * (t165 * t78 ^ 2 + t14 ^ 2) + t164;];
%% Postprocessing: Reshape Output
% From vec2symmat_5_matlab.m
res = [t8(1), t8(2), t8(4), t8(7), t8(11); t8(2), t8(3), t8(5), t8(8), t8(12); t8(4), t8(5), t8(6), t8(9), t8(13); t8(7), t8(8), t8(9), t8(10), t8(14); t8(11), t8(12), t8(13), t8(14), t8(15);];
Mq = res;
