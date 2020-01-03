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
% Datum: 2020-01-03 12:16
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
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
% StartTime: 2020-01-03 12:14:47
% EndTime: 2020-01-03 12:14:51
% DurationCPUTime: 1.25s
% Computational Cost: add. (4193->243), mult. (2794->342), div. (0->0), fcn. (2510->10), ass. (0->143)
t133 = qJ(3) + qJ(4);
t128 = qJ(5) + t133;
t119 = sin(t128);
t120 = cos(t128);
t195 = -rSges(6,1) * t120 + rSges(6,2) * t119;
t134 = qJ(1) + qJ(2);
t127 = cos(t134);
t125 = sin(t134);
t51 = -t125 * rSges(6,3) + t127 * t195;
t137 = cos(qJ(3));
t121 = t137 * pkin(3) + pkin(2);
t126 = cos(t133);
t92 = pkin(4) * t126 + t121;
t194 = t127 * t92 - t51;
t122 = t125 ^ 2;
t123 = t127 ^ 2;
t139 = -pkin(8) - pkin(7);
t124 = sin(t133);
t85 = rSges(5,1) * t124 + rSges(5,2) * t126;
t193 = m(5) * t85;
t80 = rSges(6,1) * t119 + rSges(6,2) * t120;
t192 = m(6) * t80;
t191 = -t125 / 0.2e1;
t190 = -t127 / 0.2e1;
t135 = sin(qJ(3));
t102 = rSges(4,1) * t135 + rSges(4,2) * t137;
t189 = m(4) * t102;
t188 = pkin(3) * t135;
t180 = t125 * t121 + t127 * t139;
t132 = -pkin(9) + t139;
t181 = t125 * t92 + t127 * t132;
t142 = -t127 * rSges(6,3) - t125 * t195;
t40 = t125 * t142;
t187 = t125 * (-t180 + t181) + t40;
t94 = t127 * t121;
t186 = -t94 - (t132 - t139) * t125 + t194;
t173 = t124 * t127;
t43 = pkin(4) * t173 + t127 * t80;
t185 = rSges(4,1) * t137;
t184 = rSges(5,1) * t126;
t179 = Icges(4,4) * t135;
t178 = Icges(4,4) * t137;
t177 = Icges(5,4) * t124;
t176 = Icges(5,4) * t126;
t175 = Icges(6,4) * t119;
t174 = Icges(6,4) * t120;
t172 = t127 * t135;
t86 = t125 * rSges(3,1) + t127 * rSges(3,2);
t171 = t127 * pkin(2) + t125 * pkin(7);
t170 = t122 + t123;
t78 = Icges(6,2) * t120 + t175;
t79 = Icges(6,1) * t119 + t174;
t162 = t119 * t78 - t120 * t79;
t151 = -Icges(6,2) * t119 + t174;
t47 = -Icges(6,6) * t127 + t125 * t151;
t48 = -Icges(6,6) * t125 - t127 * t151;
t154 = Icges(6,1) * t120 - t175;
t49 = -Icges(6,5) * t127 + t125 * t154;
t50 = -Icges(6,5) * t125 - t127 * t154;
t77 = Icges(6,5) * t119 + Icges(6,6) * t120;
t169 = (t119 * t50 + t120 * t48 - t125 * t77 + t127 * t162) * t191 + (t119 * t49 + t120 * t47 - t125 * t162 - t127 * t77) * t190;
t168 = -pkin(4) * t124 - t80;
t87 = t127 * rSges(3,1) - rSges(3,2) * t125;
t167 = t125 * t139 - t94;
t166 = (-rSges(4,2) * t135 + t185) * t125;
t163 = t119 * t48 - t120 * t50;
t164 = -t119 * t47 + t120 * t49;
t148 = Icges(6,5) * t120 - Icges(6,6) * t119;
t45 = -Icges(6,3) * t127 + t125 * t148;
t46 = -Icges(6,3) * t125 - t127 * t148;
t1 = t123 * t45 + (t163 * t125 + (-t164 + t46) * t127) * t125;
t2 = t122 * t46 + (t164 * t127 + (-t163 + t45) * t125) * t127;
t165 = -t127 * t1 - t125 * t2;
t152 = -Icges(5,2) * t124 + t176;
t54 = -Icges(5,6) * t127 + t125 * t152;
t155 = Icges(5,1) * t126 - t177;
t56 = -Icges(5,5) * t127 + t125 * t155;
t161 = -t124 * t54 + t126 * t56;
t55 = -Icges(5,6) * t125 - t127 * t152;
t57 = -Icges(5,5) * t125 - t127 * t155;
t160 = t124 * t55 - t126 * t57;
t83 = Icges(5,2) * t126 + t177;
t84 = Icges(5,1) * t124 + t176;
t159 = t124 * t83 - t126 * t84;
t156 = Icges(4,1) * t137 - t179;
t153 = -Icges(4,2) * t135 + t178;
t150 = Icges(4,5) * t137 - Icges(4,6) * t135;
t149 = Icges(5,5) * t126 - Icges(5,6) * t124;
t100 = Icges(4,2) * t137 + t179;
t101 = Icges(4,1) * t135 + t178;
t147 = t100 * t135 - t101 * t137;
t58 = rSges(5,2) * t173 - t125 * rSges(5,3) - t127 * t184;
t146 = rSges(4,2) * t172 - t125 * rSges(4,3) - t127 * t185;
t145 = t137 * t100 + t135 * t101 + t119 * t79 + t120 * t78 + t124 * t84 + t126 * t83 + Icges(3,3);
t82 = Icges(5,5) * t124 + Icges(5,6) * t126;
t144 = t169 + (t124 * t57 - t125 * t82 + t126 * t55 + t127 * t159) * t191 + (t124 * t56 - t125 * t159 + t126 * t54 - t127 * t82) * t190;
t143 = -t127 * rSges(5,3) + (-rSges(5,2) * t124 + t184) * t125;
t52 = -Icges(5,3) * t127 + t125 * t149;
t53 = -Icges(5,3) * t125 - t127 * t149;
t4 = t123 * t52 + (t160 * t125 + (-t161 + t53) * t127) * t125;
t5 = t122 * t53 + (t161 * t127 + (-t160 + t52) * t125) * t127;
t141 = (-t1 - t4) * t127 + (-t2 - t5) * t125;
t38 = -t146 + t171;
t30 = t143 + t180;
t25 = t142 + t181;
t31 = -t167 - t58;
t26 = -t125 * t132 + t194;
t117 = t125 * pkin(2);
t37 = t117 + (-rSges(4,3) - pkin(7)) * t127 + t166;
t99 = Icges(4,5) * t135 + Icges(4,6) * t137;
t140 = t144 + (-t125 * t99 + t147 * t127 + t135 * (-Icges(4,5) * t125 - t127 * t156) + t137 * (-Icges(4,6) * t125 - t127 * t153)) * t191 + (-t147 * t125 - t127 * t99 + t135 * (-Icges(4,5) * t127 + t125 * t156) + t137 * (-Icges(4,6) * t127 + t125 * t153)) * t190;
t138 = cos(qJ(1));
t136 = sin(qJ(1));
t131 = t138 * pkin(1);
t129 = t136 * pkin(1);
t107 = pkin(3) * t172;
t104 = rSges(2,1) * t138 - rSges(2,2) * t136;
t103 = rSges(2,1) * t136 + rSges(2,2) * t138;
t73 = t131 + t87;
t72 = t129 + t86;
t62 = -Icges(4,3) * t125 - t127 * t150;
t61 = -Icges(4,3) * t127 + t125 * t150;
t60 = t127 * t85 + t107;
t59 = (-t85 - t188) * t125;
t44 = t167 + t171;
t42 = t168 * t125;
t41 = t125 * t143;
t39 = t125 * (pkin(7) * t127 - t117 + t180);
t36 = t107 + t43;
t35 = (t168 - t188) * t125;
t34 = t131 + t38;
t33 = t129 + t37;
t28 = t131 + t31;
t27 = t129 + t30;
t22 = -t127 * t146 + t125 * (-rSges(4,3) * t127 + t166);
t21 = t131 + t26;
t20 = t129 + t25;
t17 = -t127 * t58 + t41;
t14 = -t127 * t51 + t40;
t7 = t39 + t41 + (-t44 - t58) * t127;
t6 = t127 * t186 + t187;
t3 = t39 + (-t44 + t186) * t127 + t187;
t8 = [Icges(2,3) + m(6) * (t20 ^ 2 + t21 ^ 2) + m(5) * (t27 ^ 2 + t28 ^ 2) + m(4) * (t33 ^ 2 + t34 ^ 2) + m(2) * (t103 ^ 2 + t104 ^ 2) + m(3) * (t72 ^ 2 + t73 ^ 2) + t145; m(6) * (t20 * t25 + t21 * t26) + m(5) * (t27 * t30 + t28 * t31) + m(4) * (t33 * t37 + t34 * t38) + m(3) * (t72 * t86 + t73 * t87) + t145; m(6) * (t25 ^ 2 + t26 ^ 2) + m(5) * (t30 ^ 2 + t31 ^ 2) + m(4) * (t37 ^ 2 + t38 ^ 2) + m(3) * (t86 ^ 2 + t87 ^ 2) + t145; t140 + m(6) * (t20 * t36 + t21 * t35) + m(5) * (t27 * t60 + t28 * t59) + (-t125 * t34 + t127 * t33) * t189; t140 + m(6) * (t25 * t36 + t26 * t35) + m(5) * (t30 * t60 + t31 * t59) + (-t125 * t38 + t127 * t37) * t189; m(6) * (t3 ^ 2 + t35 ^ 2 + t36 ^ 2) + m(5) * (t59 ^ 2 + t60 ^ 2 + t7 ^ 2) + m(4) * (t102 ^ 2 * t170 + t22 ^ 2) + t165 + (-t123 * t61 - t4) * t127 + (-t122 * t62 - t5 + (-t125 * t61 - t127 * t62) * t127) * t125; m(6) * (t20 * t43 + t21 * t42) + (-t125 * t28 + t127 * t27) * t193 + t144; m(6) * (t25 * t43 + t26 * t42) + (-t125 * t31 + t127 * t30) * t193 + t144; m(6) * (t3 * t6 + t35 * t42 + t36 * t43) + m(5) * (t17 * t7 + (-t125 * t59 + t127 * t60) * t85) + t141; m(5) * (t170 * t85 ^ 2 + t17 ^ 2) + m(6) * (t42 ^ 2 + t43 ^ 2 + t6 ^ 2) + t141; (-t125 * t21 + t127 * t20) * t192 + t169; (-t125 * t26 + t127 * t25) * t192 + t169; m(6) * (t14 * t3 + (-t125 * t35 + t127 * t36) * t80) + t165; m(6) * (t14 * t6 + (-t125 * t42 + t127 * t43) * t80) + t165; m(6) * (t170 * t80 ^ 2 + t14 ^ 2) + t165;];
%% Postprocessing: Reshape Output
% From vec2symmat_5_matlab.m
res = [t8(1), t8(2), t8(4), t8(7), t8(11); t8(2), t8(3), t8(5), t8(8), t8(12); t8(4), t8(5), t8(6), t8(9), t8(13); t8(7), t8(8), t8(9), t8(10), t8(14); t8(11), t8(12), t8(13), t8(14), t8(15);];
Mq = res;
