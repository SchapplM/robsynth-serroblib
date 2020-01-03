% Calculate joint inertia matrix for
% S5RPRRR12
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% pkin [8x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d1,d3,d4,d5]';
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
% Datum: 2019-12-31 19:13
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Mq = S5RPRRR12_inertiaJ_slag_vp1(qJ, ...
  pkin, m, rSges, Icges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(8,1),zeros(6,1),zeros(6,3),zeros(6,6)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RPRRR12_inertiaJ_slag_vp1: qJ has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S5RPRRR12_inertiaJ_slag_vp1: pkin has to be [8x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5RPRRR12_inertiaJ_slag_vp1: m has to be [6x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [6,3]), ...
  'S5RPRRR12_inertiaJ_slag_vp1: rSges has to be [6x3] (double)');
assert(isreal(Icges) && all(size(Icges) == [6 6]), ...
  'S5RPRRR12_inertiaJ_slag_vp1: Icges has to be [6x6] (double)'); 

%% Symbolic Calculation
% From inertia_joint_joint_floatb_twist_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 19:12:30
% EndTime: 2019-12-31 19:12:33
% DurationCPUTime: 1.38s
% Computational Cost: add. (3023->255), mult. (4109->385), div. (0->0), fcn. (4293->8), ass. (0->132)
t126 = qJ(3) + qJ(4);
t117 = sin(t126);
t118 = cos(t126);
t127 = sin(qJ(5));
t130 = cos(qJ(5));
t65 = t117 * rSges(6,3) + (rSges(6,1) * t130 - rSges(6,2) * t127) * t118;
t189 = t118 * pkin(4) + t117 * pkin(8) + t65;
t129 = sin(qJ(1));
t132 = cos(qJ(1));
t188 = t129 * t132;
t187 = (rSges(5,1) * t117 + rSges(5,2) * t118) * t132;
t128 = sin(qJ(3));
t131 = cos(qJ(3));
t186 = (rSges(4,1) * t128 + rSges(4,2) * t131) * t132;
t124 = t129 ^ 2;
t125 = t132 ^ 2;
t185 = t129 / 0.2e1;
t184 = t132 / 0.2e1;
t183 = pkin(3) * t131;
t182 = pkin(4) * t117;
t162 = t132 * t127;
t164 = t129 * t130;
t91 = t117 * t162 + t164;
t161 = t132 * t130;
t165 = t129 * t127;
t92 = -t117 * t161 + t165;
t152 = -t92 * rSges(6,1) - t91 * rSges(6,2);
t167 = t118 * t132;
t50 = rSges(6,3) * t167 - t152;
t181 = t132 * t50 + t125 * (pkin(8) * t118 - t182);
t169 = t117 * t129;
t109 = pkin(4) * t169;
t168 = t118 * t129;
t89 = -t117 * t165 + t161;
t90 = t117 * t164 + t162;
t177 = t90 * rSges(6,1) + t89 * rSges(6,2);
t49 = -rSges(6,3) * t168 + t177;
t180 = pkin(8) * t168 - t109 - t49;
t62 = Icges(6,3) * t117 + (Icges(6,5) * t130 - Icges(6,6) * t127) * t118;
t64 = Icges(6,5) * t117 + (Icges(6,1) * t130 - Icges(6,4) * t127) * t118;
t179 = t118 * t130 * t64 + t117 * t62;
t52 = t189 * t129;
t63 = Icges(6,6) * t117 + (Icges(6,4) * t130 - Icges(6,2) * t127) * t118;
t176 = t127 * t63;
t43 = Icges(6,5) * t90 + Icges(6,6) * t89 - Icges(6,3) * t168;
t45 = Icges(6,4) * t90 + Icges(6,2) * t89 - Icges(6,6) * t168;
t47 = Icges(6,1) * t90 + Icges(6,4) * t89 - Icges(6,5) * t168;
t21 = t117 * t43 + (-t127 * t45 + t130 * t47) * t118;
t175 = t21 * t132;
t44 = Icges(6,5) * t92 + Icges(6,6) * t91 + Icges(6,3) * t167;
t46 = Icges(6,4) * t92 + Icges(6,2) * t91 + Icges(6,6) * t167;
t48 = Icges(6,1) * t92 + Icges(6,4) * t91 + Icges(6,5) * t167;
t22 = t117 * t44 + (-t127 * t46 + t130 * t48) * t118;
t174 = t22 * t129;
t173 = Icges(4,4) * t128;
t172 = Icges(4,4) * t131;
t171 = Icges(5,4) * t117;
t170 = Icges(5,4) * t118;
t166 = t128 * t129;
t163 = t129 * t131;
t133 = -pkin(7) - pkin(6);
t160 = t132 * t128 * pkin(3) + t129 * t133;
t159 = t132 * pkin(1) + t129 * qJ(2);
t158 = t124 + t125;
t74 = rSges(5,1) * t169 + rSges(5,2) * t168 + t132 * rSges(5,3);
t157 = rSges(4,1) * t166 + rSges(4,2) * t163 + t132 * rSges(4,3);
t120 = t132 * qJ(2);
t156 = t120 + t160;
t155 = (-rSges(6,3) - pkin(8)) * t118;
t16 = t43 * t167 + t91 * t45 + t92 * t47;
t17 = t44 * t167 + t91 * t46 + t92 * t48;
t10 = t17 * t129 + t16 * t132;
t138 = Icges(5,5) * t117 + Icges(5,6) * t118;
t68 = Icges(5,3) * t132 + t138 * t129;
t69 = Icges(5,3) * t129 - t138 * t132;
t14 = -t43 * t168 + t89 * t45 + t90 * t47;
t15 = -t44 * t168 + t89 * t46 + t90 * t48;
t9 = t15 * t129 + t14 * t132;
t154 = (t125 * t68 + t9) * t132 + (t124 * t69 + t10 + (t129 * t68 + t132 * t69) * t132) * t129;
t26 = -t62 * t168 + t89 * t63 + t90 * t64;
t3 = t26 * t117 + (-t129 * t14 + t132 * t15) * t118;
t27 = t62 * t167 + t91 * t63 + t92 * t64;
t4 = t27 * t117 + (-t129 * t16 + t132 * t17) * t118;
t153 = -t9 * t168 / 0.2e1 + t3 * t184 + t4 * t185 + t117 * (t174 + t175) / 0.2e1 + t10 * t167 / 0.2e1;
t96 = -Icges(5,2) * t117 + t170;
t97 = Icges(5,1) * t118 - t171;
t147 = t117 * t97 + t118 * t96;
t114 = pkin(3) * t163;
t98 = t118 * rSges(5,1) - t117 * rSges(5,2);
t66 = t129 * t98 + t114;
t67 = (-t98 - t183) * t132;
t144 = t66 * t129 - t67 * t132;
t143 = Icges(4,1) * t128 + t172;
t142 = Icges(5,1) * t117 + t170;
t141 = Icges(4,2) * t131 + t173;
t140 = Icges(5,2) * t118 + t171;
t139 = Icges(4,5) * t128 + Icges(4,6) * t131;
t113 = pkin(3) * t166;
t137 = -t132 * t133 + t113 + t159;
t57 = t120 + t186 + (-rSges(4,3) - pkin(1) - pkin(6)) * t129;
t58 = t132 * pkin(6) + t157 + t159;
t136 = m(4) * (t129 * t57 - t132 * t58);
t54 = t187 + (-rSges(5,3) - pkin(1)) * t129 + t156;
t55 = t137 + t74;
t135 = m(5) * (t129 * t54 - t132 * t55);
t95 = Icges(5,5) * t118 - Icges(5,6) * t117;
t134 = t174 / 0.2e1 + t175 / 0.2e1 + (-t117 * (Icges(5,6) * t129 - t140 * t132) + t118 * (Icges(5,5) * t129 - t142 * t132) + t129 * t95 - t147 * t132 + t27) * t185 + (-t117 * (Icges(5,6) * t132 + t140 * t129) + t118 * (Icges(5,5) * t132 + t142 * t129) + t147 * t129 + t132 * t95 + t26) * t184;
t106 = t132 * rSges(2,1) - t129 * rSges(2,2);
t105 = t131 * rSges(4,1) - t128 * rSges(4,2);
t104 = -t129 * rSges(2,1) - t132 * rSges(2,2);
t88 = t113 + (-pkin(6) - t133) * t132;
t86 = -t132 * rSges(3,2) + t129 * rSges(3,3) + t159;
t85 = t132 * rSges(3,3) + t120 + (rSges(3,2) - pkin(1)) * t129;
t78 = Icges(4,3) * t129 - t139 * t132;
t77 = Icges(4,3) * t132 + t139 * t129;
t76 = t132 * (-t129 * pkin(6) - t160);
t61 = t132 * (t129 * rSges(5,3) - t187);
t53 = t189 * t132;
t51 = -t129 * t157 + (t129 * rSges(4,3) - t186) * t132;
t41 = (-t189 - t183) * t132;
t40 = t114 + t52;
t39 = -t129 * t74 + t61;
t34 = -t117 * t50 + t65 * t167;
t33 = t117 * t49 + t65 * t168;
t32 = t129 * t155 + t109 + t137 + t177;
t31 = -t129 * pkin(1) + (t155 + t182) * t132 + t152 + t156;
t30 = t61 + t76 + (-t74 - t88) * t129;
t29 = (-t118 * t176 + t179) * t117;
t28 = (-t129 * t50 - t132 * t49) * t118;
t23 = t180 * t129 + t181;
t18 = t76 + (-t88 + t180) * t129 + t181;
t1 = [-t128 * (-Icges(4,2) * t128 + t172) + t131 * (Icges(4,1) * t131 - t173) - t117 * t96 + Icges(3,1) + Icges(2,3) + (t97 - t176) * t118 + m(6) * (t31 ^ 2 + t32 ^ 2) + m(5) * (t54 ^ 2 + t55 ^ 2) + m(4) * (t57 ^ 2 + t58 ^ 2) + m(3) * (t85 ^ 2 + t86 ^ 2) + m(2) * (t104 ^ 2 + t106 ^ 2) + t179; m(6) * (t129 * t31 - t132 * t32) + t135 + t136 + m(3) * (t129 * t85 - t132 * t86); 0.2e1 * (m(3) / 0.2e1 + m(4) / 0.2e1 + m(5) / 0.2e1 + m(6) / 0.2e1) * t158; (-t128 * (Icges(4,6) * t129 - t141 * t132) + t131 * (Icges(4,5) * t129 - t143 * t132)) * t185 + (-t128 * (Icges(4,6) * t132 + t141 * t129) + t131 * (Icges(4,5) * t132 + t143 * t129)) * t184 + m(6) * (t40 * t31 + t41 * t32) + m(5) * (t66 * t54 + t67 * t55) + t105 * t136 + (t124 / 0.2e1 + t125 / 0.2e1) * (Icges(4,5) * t131 - Icges(4,6) * t128) + t134; m(5) * t144 + m(6) * (t40 * t129 - t41 * t132) + m(4) * t158 * t105; m(6) * (t18 ^ 2 + t40 ^ 2 + t41 ^ 2) + m(5) * (t30 ^ 2 + t66 ^ 2 + t67 ^ 2) + m(4) * (t158 * t105 ^ 2 + t51 ^ 2) + t132 * (t125 * t77 + t78 * t188) + t129 * (t124 * t78 + t77 * t188) + t154; m(6) * (t52 * t31 - t53 * t32) + t98 * t135 + t134; m(6) * (t52 * t129 + t53 * t132) + m(5) * t158 * t98; m(6) * (t23 * t18 + t52 * t40 - t53 * t41) + m(5) * (t144 * t98 + t39 * t30) + t154; m(5) * (t158 * t98 ^ 2 + t39 ^ 2) + m(6) * (t23 ^ 2 + t52 ^ 2 + t53 ^ 2) + t154; m(6) * (t34 * t31 + t33 * t32) + t29 + ((t27 / 0.2e1 + t22 / 0.2e1) * t132 + (-t26 / 0.2e1 - t21 / 0.2e1) * t129) * t118; m(6) * (t34 * t129 - t33 * t132); m(6) * (t28 * t18 + t33 * t41 + t34 * t40) + t153; m(6) * (t28 * t23 - t33 * t53 + t34 * t52) + t153; m(6) * (t28 ^ 2 + t33 ^ 2 + t34 ^ 2) + t117 * t29 + (-t129 * t3 + t132 * t4 + t117 * (-t129 * t21 + t132 * t22)) * t118;];
%% Postprocessing: Reshape Output
% From vec2symmat_5_matlab.m
res = [t1(1), t1(2), t1(4), t1(7), t1(11); t1(2), t1(3), t1(5), t1(8), t1(12); t1(4), t1(5), t1(6), t1(9), t1(13); t1(7), t1(8), t1(9), t1(10), t1(14); t1(11), t1(12), t1(13), t1(14), t1(15);];
Mq = res;
