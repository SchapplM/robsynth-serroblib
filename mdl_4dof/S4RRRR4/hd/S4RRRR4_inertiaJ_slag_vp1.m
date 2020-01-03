% Calculate joint inertia matrix for
% S4RRRR4
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
% Datum: 2019-12-31 17:26
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Mq = S4RRRR4_inertiaJ_slag_vp1(qJ, ...
  pkin, m, rSges, Icges)
%% Coder Information
%#codegen
%$cgargs {zeros(4,1),zeros(7,1),zeros(5,1),zeros(5,3),zeros(5,6)}
assert(isreal(qJ) && all(size(qJ) == [4 1]), ...
  'S4RRRR4_inertiaJ_slag_vp1: qJ has to be [4x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [7 1]), ...
  'S4RRRR4_inertiaJ_slag_vp1: pkin has to be [7x1] (double)');
assert(isreal(m) && all(size(m) == [5 1]), ...
  'S4RRRR4_inertiaJ_slag_vp1: m has to be [5x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [5,3]), ...
  'S4RRRR4_inertiaJ_slag_vp1: rSges has to be [5x3] (double)');
assert(isreal(Icges) && all(size(Icges) == [5 6]), ...
  'S4RRRR4_inertiaJ_slag_vp1: Icges has to be [5x6] (double)'); 

%% Symbolic Calculation
% From inertia_joint_joint_floatb_twist_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 17:25:41
% EndTime: 2019-12-31 17:25:44
% DurationCPUTime: 1.12s
% Computational Cost: add. (2843->234), mult. (3809->354), div. (0->0), fcn. (4051->8), ass. (0->129)
t116 = qJ(2) + qJ(3);
t110 = cos(t116);
t122 = cos(qJ(1));
t153 = t110 * t122;
t109 = sin(t116);
t154 = t109 * t122;
t117 = sin(qJ(4));
t150 = t122 * t117;
t119 = sin(qJ(1));
t120 = cos(qJ(4));
t151 = t119 * t120;
t87 = -t110 * t150 + t151;
t149 = t122 * t120;
t152 = t119 * t117;
t88 = t110 * t149 + t152;
t50 = t88 * rSges(5,1) + t87 * rSges(5,2) + rSges(5,3) * t154;
t174 = pkin(3) * t153 + pkin(7) * t154 + t50;
t115 = t122 ^ 2;
t156 = Icges(4,4) * t110;
t129 = -Icges(4,2) * t109 + t156;
t72 = Icges(4,6) * t119 + t122 * t129;
t157 = Icges(4,4) * t109;
t131 = Icges(4,1) * t110 - t157;
t74 = Icges(4,5) * t119 + t122 * t131;
t136 = -t109 * t72 + t110 * t74;
t71 = -Icges(4,6) * t122 + t119 * t129;
t73 = -Icges(4,5) * t122 + t119 * t131;
t137 = t109 * t71 - t110 * t73;
t127 = Icges(4,5) * t110 - Icges(4,6) * t109;
t69 = -Icges(4,3) * t122 + t119 * t127;
t70 = Icges(4,3) * t119 + t122 * t127;
t155 = t109 * t119;
t85 = -t110 * t152 - t149;
t86 = t110 * t151 - t150;
t43 = Icges(5,5) * t86 + Icges(5,6) * t85 + Icges(5,3) * t155;
t45 = Icges(5,4) * t86 + Icges(5,2) * t85 + Icges(5,6) * t155;
t47 = Icges(5,1) * t86 + Icges(5,4) * t85 + Icges(5,5) * t155;
t13 = t155 * t43 + t45 * t85 + t47 * t86;
t44 = Icges(5,5) * t88 + Icges(5,6) * t87 + Icges(5,3) * t154;
t46 = Icges(5,4) * t88 + Icges(5,2) * t87 + Icges(5,6) * t154;
t48 = Icges(5,1) * t88 + Icges(5,4) * t87 + Icges(5,5) * t154;
t14 = t155 * t44 + t46 * t85 + t48 * t86;
t8 = t119 * t14 - t122 * t13;
t173 = -t115 * t69 - (t136 * t119 + (t137 - t70) * t122) * t119 - t8;
t114 = t119 ^ 2;
t172 = t119 / 0.2e1;
t171 = -t122 / 0.2e1;
t15 = t154 * t43 + t45 * t87 + t47 * t88;
t16 = t154 * t44 + t46 * t87 + t48 * t88;
t9 = t119 * t16 - t122 * t15;
t170 = (t114 * t70 + t9 + (t137 * t122 + (t136 - t69) * t119) * t122) * t119;
t118 = sin(qJ(2));
t169 = pkin(2) * t118;
t168 = pkin(3) * t110;
t121 = cos(qJ(2));
t108 = pkin(2) * t121 + pkin(1);
t102 = t122 * t108;
t113 = t122 * pkin(5);
t123 = -pkin(6) - pkin(5);
t148 = t122 * t123;
t167 = t119 * (t148 + t113 + (-pkin(1) + t108) * t119) + t122 * (-t122 * pkin(1) + t102 + (-pkin(5) - t123) * t119);
t125 = rSges(4,1) * t153 - rSges(4,2) * t154 + t119 * rSges(4,3);
t138 = rSges(4,1) * t110 - rSges(4,2) * t109;
t38 = t119 * (-t122 * rSges(4,3) + t119 * t138) + t122 * t125;
t66 = -t110 * rSges(5,3) + (rSges(5,1) * t120 - rSges(5,2) * t117) * t109;
t166 = -pkin(3) * t109 + pkin(7) * t110 - t66;
t165 = rSges(3,1) * t121;
t164 = rSges(3,2) * t118;
t64 = -Icges(5,6) * t110 + (Icges(5,4) * t120 - Icges(5,2) * t117) * t109;
t163 = t117 * t64;
t162 = t122 * rSges(3,3);
t20 = -t110 * t43 + (-t117 * t45 + t120 * t47) * t109;
t161 = t20 * t122;
t21 = -t110 * t44 + (-t117 * t46 + t120 * t48) * t109;
t160 = t21 * t119;
t159 = Icges(3,4) * t118;
t158 = Icges(3,4) * t121;
t146 = t119 * rSges(3,3) + t122 * t165;
t145 = t114 + t115;
t93 = rSges(4,1) * t109 + rSges(4,2) * t110;
t144 = -t93 - t169;
t140 = -t86 * rSges(5,1) - t85 * rSges(5,2);
t49 = rSges(5,3) * t155 - t140;
t22 = t119 * t49 + t114 * (pkin(7) * t109 + t168) + t174 * t122;
t143 = -t119 * t123 + t102;
t63 = -Icges(5,3) * t110 + (Icges(5,5) * t120 - Icges(5,6) * t117) * t109;
t65 = -Icges(5,5) * t110 + (Icges(5,1) * t120 - Icges(5,4) * t117) * t109;
t25 = t155 * t63 + t64 * t85 + t65 * t86;
t3 = -t25 * t110 + (t119 * t13 + t122 * t14) * t109;
t26 = t154 * t63 + t64 * t87 + t65 * t88;
t4 = -t26 * t110 + (t119 * t15 + t122 * t16) * t109;
t142 = t8 * t155 / 0.2e1 + t3 * t171 + t4 * t172 - t110 * (t160 - t161) / 0.2e1 + t9 * t154 / 0.2e1;
t141 = t166 - t169;
t139 = -t164 + t165;
t91 = Icges(4,2) * t110 + t157;
t92 = Icges(4,1) * t109 + t156;
t135 = -t109 * t91 + t110 * t92;
t132 = Icges(3,1) * t121 - t159;
t130 = -Icges(3,2) * t118 + t158;
t128 = Icges(3,5) * t121 - Icges(3,6) * t118;
t126 = t122 * t173 + t170;
t90 = Icges(4,5) * t109 + Icges(4,6) * t110;
t124 = t160 / 0.2e1 - t161 / 0.2e1 + (t109 * t74 + t110 * t72 + t119 * t90 + t122 * t135 + t26) * t172 + (t109 * t73 + t110 * t71 + t119 * t135 - t122 * t90 + t25) * t171;
t101 = rSges(2,1) * t122 - rSges(2,2) * t119;
t100 = -rSges(2,1) * t119 - rSges(2,2) * t122;
t99 = rSges(3,1) * t118 + rSges(3,2) * t121;
t78 = Icges(3,3) * t119 + t122 * t128;
t77 = -Icges(3,3) * t122 + t119 * t128;
t68 = t144 * t122;
t67 = t144 * t119;
t58 = t119 * pkin(5) + (pkin(1) - t164) * t122 + t146;
t57 = t162 + t113 + (-pkin(1) - t139) * t119;
t56 = t109 * t120 * t65;
t55 = t125 + t143;
t54 = (rSges(4,3) - t123) * t122 + (-t108 - t138) * t119;
t53 = t166 * t122;
t52 = t166 * t119;
t51 = t122 * (-t122 * t164 + t146) + (t119 * t139 - t162) * t119;
t40 = t141 * t122;
t39 = t141 * t119;
t33 = t143 + t174;
t32 = -t148 + (-t168 - t108 + (-rSges(5,3) - pkin(7)) * t109) * t119 + t140;
t31 = -t110 * t50 - t154 * t66;
t30 = t110 * t49 + t155 * t66;
t29 = -t109 * t163 - t110 * t63 + t56;
t28 = t38 + t167;
t27 = (-t119 * t50 + t122 * t49) * t109;
t11 = t22 + t167;
t1 = [t118 * (Icges(3,1) * t118 + t158) + t121 * (Icges(3,2) * t121 + t159) + Icges(2,3) + t56 + (-t63 + t91) * t110 + (t92 - t163) * t109 + m(5) * (t32 ^ 2 + t33 ^ 2) + m(4) * (t54 ^ 2 + t55 ^ 2) + m(3) * (t57 ^ 2 + t58 ^ 2) + m(2) * (t100 ^ 2 + t101 ^ 2); (t118 * (Icges(3,5) * t119 + t122 * t132) + t121 * (Icges(3,6) * t119 + t122 * t130)) * t172 + (t118 * (-Icges(3,5) * t122 + t119 * t132) + t121 * (-Icges(3,6) * t122 + t119 * t130)) * t171 + m(5) * (t32 * t40 + t33 * t39) + m(4) * (t54 * t68 + t55 * t67) + m(3) * (-t119 * t58 - t122 * t57) * t99 + (t114 / 0.2e1 + t115 / 0.2e1) * (Icges(3,5) * t118 + Icges(3,6) * t121) + t124; m(5) * (t11 ^ 2 + t39 ^ 2 + t40 ^ 2) + m(4) * (t28 ^ 2 + t67 ^ 2 + t68 ^ 2) + m(3) * (t145 * t99 ^ 2 + t51 ^ 2) + t119 * t114 * t78 + t170 + (-t115 * t77 + (-t119 * t77 + t122 * t78) * t119 + t173) * t122; m(5) * (t32 * t53 + t33 * t52) + m(4) * (-t119 * t55 - t122 * t54) * t93 + t124; m(5) * (t11 * t22 + t39 * t52 + t40 * t53) + m(4) * (t38 * t28 + (-t119 * t67 - t122 * t68) * t93) + t126; m(4) * (t145 * t93 ^ 2 + t38 ^ 2) + m(5) * (t22 ^ 2 + t52 ^ 2 + t53 ^ 2) + t126; m(5) * (t30 * t32 + t31 * t33) - t29 * t110 + ((t21 / 0.2e1 + t26 / 0.2e1) * t122 + (t20 / 0.2e1 + t25 / 0.2e1) * t119) * t109; m(5) * (t11 * t27 + t30 * t40 + t31 * t39) + t142; m(5) * (t22 * t27 + t30 * t53 + t31 * t52) + t142; m(5) * (t27 ^ 2 + t30 ^ 2 + t31 ^ 2) + t110 ^ 2 * t29 + (t122 * t4 + t119 * t3 - t110 * (t119 * t20 + t122 * t21)) * t109;];
%% Postprocessing: Reshape Output
% From vec2symmat_4_matlab.m
res = [t1(1), t1(2), t1(4), t1(7); t1(2), t1(3), t1(5), t1(8); t1(4), t1(5), t1(6), t1(9); t1(7), t1(8), t1(9), t1(10);];
Mq = res;
