% Calculate joint inertia matrix for
% S5RPRRR9
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% pkin [9x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d1,d3,d4,d5,theta2]';
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
% Datum: 2019-12-31 19:08
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Mq = S5RPRRR9_inertiaJ_slag_vp1(qJ, ...
  pkin, m, rSges, Icges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(9,1),zeros(6,1),zeros(6,3),zeros(6,6)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RPRRR9_inertiaJ_slag_vp1: qJ has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S5RPRRR9_inertiaJ_slag_vp1: pkin has to be [9x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5RPRRR9_inertiaJ_slag_vp1: m has to be [6x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [6,3]), ...
  'S5RPRRR9_inertiaJ_slag_vp1: rSges has to be [6x3] (double)');
assert(isreal(Icges) && all(size(Icges) == [6 6]), ...
  'S5RPRRR9_inertiaJ_slag_vp1: Icges has to be [6x6] (double)'); 

%% Symbolic Calculation
% From inertia_joint_joint_floatb_twist_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 19:07:19
% EndTime: 2019-12-31 19:07:23
% DurationCPUTime: 1.35s
% Computational Cost: add. (4513->258), mult. (4115->398), div. (0->0), fcn. (4311->10), ass. (0->135)
t129 = sin(qJ(1));
t123 = t129 ^ 2;
t122 = pkin(9) + qJ(3);
t117 = qJ(4) + t122;
t113 = cos(t117);
t131 = cos(qJ(1));
t162 = t113 * t131;
t112 = sin(t117);
t163 = t112 * t131;
t130 = cos(qJ(5));
t159 = t129 * t130;
t128 = sin(qJ(5));
t161 = t128 * t131;
t89 = -t113 * t161 + t159;
t158 = t130 * t131;
t160 = t129 * t128;
t90 = t113 * t158 + t160;
t53 = t90 * rSges(6,1) + t89 * rSges(6,2) + rSges(6,3) * t163;
t183 = pkin(4) * t162 + pkin(8) * t163 + t53;
t124 = t131 ^ 2;
t165 = Icges(5,4) * t113;
t138 = -Icges(5,2) * t112 + t165;
t74 = Icges(5,6) * t129 + t138 * t131;
t166 = Icges(5,4) * t112;
t140 = Icges(5,1) * t113 - t166;
t76 = Icges(5,5) * t129 + t140 * t131;
t145 = -t112 * t74 + t113 * t76;
t73 = -Icges(5,6) * t131 + t138 * t129;
t75 = -Icges(5,5) * t131 + t140 * t129;
t146 = t112 * t73 - t113 * t75;
t136 = Icges(5,5) * t113 - Icges(5,6) * t112;
t71 = -Icges(5,3) * t131 + t136 * t129;
t72 = Icges(5,3) * t129 + t136 * t131;
t164 = t112 * t129;
t87 = -t113 * t160 - t158;
t88 = t113 * t159 - t161;
t46 = Icges(6,5) * t88 + Icges(6,6) * t87 + Icges(6,3) * t164;
t48 = Icges(6,4) * t88 + Icges(6,2) * t87 + Icges(6,6) * t164;
t50 = Icges(6,1) * t88 + Icges(6,4) * t87 + Icges(6,5) * t164;
t14 = t46 * t164 + t48 * t87 + t50 * t88;
t47 = Icges(6,5) * t90 + Icges(6,6) * t89 + Icges(6,3) * t163;
t49 = Icges(6,4) * t90 + Icges(6,2) * t89 + Icges(6,6) * t163;
t51 = Icges(6,1) * t90 + Icges(6,4) * t89 + Icges(6,5) * t163;
t15 = t47 * t164 + t49 * t87 + t51 * t88;
t8 = t15 * t129 - t131 * t14;
t182 = -t124 * t71 - (t145 * t129 + (t146 - t72) * t131) * t129 - t8;
t181 = t129 / 0.2e1;
t180 = -t131 / 0.2e1;
t16 = t46 * t163 + t89 * t48 + t90 * t50;
t17 = t47 * t163 + t89 * t49 + t90 * t51;
t9 = t17 * t129 - t131 * t16;
t179 = (t123 * t72 + t9 + (t146 * t131 + (t145 - t71) * t129) * t131) * t129;
t115 = sin(t122);
t178 = pkin(3) * t115;
t177 = pkin(4) * t113;
t127 = -pkin(6) - qJ(2);
t126 = cos(pkin(9));
t114 = t126 * pkin(2) + pkin(1);
t116 = cos(t122);
t103 = pkin(3) * t116 + t114;
t97 = t131 * t103;
t176 = t131 * (-t114 * t131 + t97) + (t103 - t114) * t123;
t64 = -t113 * rSges(6,3) + (rSges(6,1) * t130 - rSges(6,2) * t128) * t112;
t175 = -pkin(4) * t112 + pkin(8) * t113 - t64;
t134 = rSges(5,1) * t162 - rSges(5,2) * t163 + t129 * rSges(5,3);
t147 = rSges(5,1) * t113 - rSges(5,2) * t112;
t40 = t129 * (-t131 * rSges(5,3) + t147 * t129) + t131 * t134;
t174 = rSges(4,1) * t116;
t173 = rSges(4,2) * t115;
t62 = -Icges(6,6) * t113 + (Icges(6,4) * t130 - Icges(6,2) * t128) * t112;
t172 = t128 * t62;
t20 = -t113 * t46 + (-t128 * t48 + t130 * t50) * t112;
t171 = t20 * t131;
t21 = -t113 * t47 + (-t128 * t49 + t130 * t51) * t112;
t170 = t21 * t129;
t169 = rSges(3,3) + qJ(2);
t168 = Icges(4,4) * t115;
t167 = Icges(4,4) * t116;
t156 = t129 * rSges(4,3) + t131 * t174;
t154 = t123 + t124;
t95 = rSges(5,1) * t112 + rSges(5,2) * t113;
t153 = -t95 - t178;
t149 = -t88 * rSges(6,1) - t87 * rSges(6,2);
t52 = rSges(6,3) * t164 - t149;
t22 = t129 * t52 + t123 * (pkin(8) * t112 + t177) + t183 * t131;
t121 = -pkin(7) + t127;
t152 = -t129 * t121 + t97;
t61 = -Icges(6,3) * t113 + (Icges(6,5) * t130 - Icges(6,6) * t128) * t112;
t63 = -Icges(6,5) * t113 + (Icges(6,1) * t130 - Icges(6,4) * t128) * t112;
t26 = t61 * t164 + t62 * t87 + t63 * t88;
t3 = -t26 * t113 + (t129 * t14 + t131 * t15) * t112;
t27 = t61 * t163 + t89 * t62 + t90 * t63;
t4 = -t27 * t113 + (t129 * t16 + t131 * t17) * t112;
t151 = t8 * t164 / 0.2e1 + t3 * t180 + t4 * t181 - t113 * (t170 - t171) / 0.2e1 + t9 * t163 / 0.2e1;
t150 = t175 - t178;
t148 = -t173 + t174;
t93 = Icges(5,2) * t113 + t166;
t94 = Icges(5,1) * t112 + t165;
t144 = -t112 * t93 + t113 * t94;
t141 = Icges(4,1) * t116 - t168;
t139 = -Icges(4,2) * t115 + t167;
t137 = Icges(4,5) * t116 - Icges(4,6) * t115;
t135 = t182 * t131 + t179;
t125 = sin(pkin(9));
t133 = rSges(3,1) * t126 - rSges(3,2) * t125 + pkin(1);
t92 = Icges(5,5) * t112 + Icges(5,6) * t113;
t132 = -t171 / 0.2e1 + t170 / 0.2e1 + (t112 * t76 + t113 * t74 + t129 * t92 + t144 * t131 + t27) * t181 + (t112 * t75 + t113 * t73 + t144 * t129 - t131 * t92 + t26) * t180;
t109 = rSges(2,1) * t131 - t129 * rSges(2,2);
t108 = -t129 * rSges(2,1) - rSges(2,2) * t131;
t102 = rSges(4,1) * t115 + rSges(4,2) * t116;
t80 = Icges(4,3) * t129 + t137 * t131;
t79 = -Icges(4,3) * t131 + t137 * t129;
t70 = t169 * t129 + t133 * t131;
t69 = -t133 * t129 + t169 * t131;
t66 = t153 * t131;
t65 = t153 * t129;
t60 = -t129 * t127 + (t114 - t173) * t131 + t156;
t59 = (rSges(4,3) - t127) * t131 + (-t114 - t148) * t129;
t58 = t112 * t130 * t63;
t55 = t134 + t152;
t54 = (rSges(5,3) - t121) * t131 + (-t103 - t147) * t129;
t45 = t175 * t131;
t44 = t175 * t129;
t43 = t131 * (-t131 * t173 + t156) + (-t131 * rSges(4,3) + t148 * t129) * t129;
t39 = t150 * t131;
t38 = t150 * t129;
t33 = t152 + t183;
t32 = -t121 * t131 + (-t177 - t103 + (-rSges(6,3) - pkin(8)) * t112) * t129 + t149;
t31 = -t113 * t53 - t64 * t163;
t30 = t113 * t52 + t64 * t164;
t29 = -t112 * t172 - t113 * t61 + t58;
t28 = (-t129 * t53 + t131 * t52) * t112;
t23 = t40 + t176;
t11 = t22 + t176;
t1 = [Icges(3,2) * t126 ^ 2 + t116 * (Icges(4,2) * t116 + t168) + t115 * (Icges(4,1) * t115 + t167) + Icges(2,3) + t58 + (Icges(3,1) * t125 + 0.2e1 * Icges(3,4) * t126) * t125 + (-t61 + t93) * t113 + (t94 - t172) * t112 + m(6) * (t32 ^ 2 + t33 ^ 2) + m(5) * (t54 ^ 2 + t55 ^ 2) + m(4) * (t59 ^ 2 + t60 ^ 2) + m(3) * (t69 ^ 2 + t70 ^ 2) + m(2) * (t108 ^ 2 + t109 ^ 2); m(6) * (t129 * t32 - t131 * t33) + m(5) * (t129 * t54 - t131 * t55) + m(4) * (t129 * t59 - t131 * t60) + m(3) * (t129 * t69 - t131 * t70); 0.2e1 * (m(3) / 0.2e1 + m(4) / 0.2e1 + m(5) / 0.2e1 + m(6) / 0.2e1) * t154; (t115 * (-Icges(4,5) * t131 + t141 * t129) + t116 * (-Icges(4,6) * t131 + t139 * t129)) * t180 + (t115 * (Icges(4,5) * t129 + t141 * t131) + t116 * (Icges(4,6) * t129 + t139 * t131)) * t181 + m(6) * (t32 * t39 + t33 * t38) + m(5) * (t54 * t66 + t55 * t65) + m(4) * (-t129 * t60 - t131 * t59) * t102 + (t123 / 0.2e1 + t124 / 0.2e1) * (Icges(4,5) * t115 + Icges(4,6) * t116) + t132; m(5) * (t66 * t129 - t131 * t65) + m(6) * (t39 * t129 - t131 * t38); m(6) * (t11 ^ 2 + t38 ^ 2 + t39 ^ 2) + m(5) * (t23 ^ 2 + t65 ^ 2 + t66 ^ 2) + t129 * t123 * t80 + m(4) * (t154 * t102 ^ 2 + t43 ^ 2) + t179 + (-t124 * t79 + (-t129 * t79 + t131 * t80) * t129 + t182) * t131; m(6) * (t32 * t45 + t33 * t44) + m(5) * (-t129 * t55 - t131 * t54) * t95 + t132; m(6) * (t45 * t129 - t131 * t44); m(6) * (t11 * t22 + t38 * t44 + t39 * t45) + m(5) * (t40 * t23 + (-t129 * t65 - t131 * t66) * t95) + t135; m(5) * (t154 * t95 ^ 2 + t40 ^ 2) + m(6) * (t22 ^ 2 + t44 ^ 2 + t45 ^ 2) + t135; m(6) * (t30 * t32 + t31 * t33) - t29 * t113 + ((t21 / 0.2e1 + t27 / 0.2e1) * t131 + (t26 / 0.2e1 + t20 / 0.2e1) * t129) * t112; m(6) * (t30 * t129 - t131 * t31); m(6) * (t11 * t28 + t30 * t39 + t31 * t38) + t151; m(6) * (t22 * t28 + t30 * t45 + t31 * t44) + t151; m(6) * (t28 ^ 2 + t30 ^ 2 + t31 ^ 2) + t113 ^ 2 * t29 + (t131 * t4 + t129 * t3 - t113 * (t129 * t20 + t131 * t21)) * t112;];
%% Postprocessing: Reshape Output
% From vec2symmat_5_matlab.m
res = [t1(1), t1(2), t1(4), t1(7), t1(11); t1(2), t1(3), t1(5), t1(8), t1(12); t1(4), t1(5), t1(6), t1(9), t1(13); t1(7), t1(8), t1(9), t1(10), t1(14); t1(11), t1(12), t1(13), t1(14), t1(15);];
Mq = res;
