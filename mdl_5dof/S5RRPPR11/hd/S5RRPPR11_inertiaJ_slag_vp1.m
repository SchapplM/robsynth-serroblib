% Calculate joint inertia matrix for
% S5RRPPR11
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% pkin [8x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d1,d2,d5,theta4]';
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
% Datum: 2019-12-31 19:48
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Mq = S5RRPPR11_inertiaJ_slag_vp1(qJ, ...
  pkin, m, rSges, Icges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(8,1),zeros(6,1),zeros(6,3),zeros(6,6)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RRPPR11_inertiaJ_slag_vp1: qJ has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S5RRPPR11_inertiaJ_slag_vp1: pkin has to be [8x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5RRPPR11_inertiaJ_slag_vp1: m has to be [6x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [6,3]), ...
  'S5RRPPR11_inertiaJ_slag_vp1: rSges has to be [6x3] (double)');
assert(isreal(Icges) && all(size(Icges) == [6 6]), ...
  'S5RRPPR11_inertiaJ_slag_vp1: Icges has to be [6x6] (double)'); 

%% Symbolic Calculation
% From inertia_joint_joint_floatb_twist_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 19:46:27
% EndTime: 2019-12-31 19:46:32
% DurationCPUTime: 2.01s
% Computational Cost: add. (2216->305), mult. (3641->453), div. (0->0), fcn. (3748->8), ass. (0->147)
t201 = Icges(4,1) + Icges(3,3);
t129 = sin(qJ(2));
t131 = cos(qJ(2));
t200 = (-Icges(4,4) + Icges(3,5)) * t131 + (Icges(4,5) - Icges(3,6)) * t129;
t130 = sin(qJ(1));
t199 = -t130 / 0.2e1;
t187 = t130 / 0.2e1;
t132 = cos(qJ(1));
t198 = -t132 / 0.2e1;
t197 = t132 / 0.2e1;
t196 = t129 / 0.2e1;
t195 = t201 * t130 + t200 * t132;
t194 = -t200 * t130 + t201 * t132;
t127 = cos(pkin(8));
t111 = pkin(4) * t127 + pkin(3);
t128 = -pkin(7) - qJ(4);
t126 = sin(pkin(8));
t176 = t129 * t132;
t159 = t126 * t176;
t170 = t131 * t132;
t121 = pkin(8) + qJ(5);
t113 = cos(t121);
t112 = sin(t121);
t175 = t130 * t112;
t61 = t113 * t176 - t175;
t174 = t130 * t113;
t62 = t112 * t176 + t174;
t32 = t62 * rSges(6,1) + t61 * rSges(6,2) + rSges(6,3) * t170;
t193 = pkin(4) * t159 + t130 * t111 - t128 * t170 + t32;
t192 = -qJ(4) - t128;
t123 = t130 ^ 2;
t125 = t132 ^ 2;
t191 = 0.2e1 * t131;
t190 = m(4) / 0.2e1;
t189 = m(5) / 0.2e1;
t188 = m(6) / 0.2e1;
t186 = pkin(4) * t126;
t168 = pkin(2) * t170 + qJ(3) * t176;
t177 = qJ(3) * t129;
t185 = t123 * (pkin(2) * t131 + t177) + t132 * t168;
t98 = pkin(2) * t129 - qJ(3) * t131;
t184 = rSges(4,2) * t129 + rSges(4,3) * t131 - t98;
t183 = t132 * rSges(4,1);
t182 = t132 * rSges(3,3);
t181 = Icges(3,4) * t129;
t180 = Icges(3,4) * t131;
t179 = Icges(4,6) * t129;
t178 = Icges(4,6) * t131;
t173 = t130 * t126;
t172 = t130 * t127;
t171 = t130 * t131;
t169 = t132 * t111;
t167 = t132 * pkin(1) + t130 * pkin(6);
t166 = t130 * pkin(3) + qJ(4) * t170;
t165 = t123 + t125;
t164 = t189 + t188;
t84 = t127 * t176 - t173;
t85 = t159 + t172;
t163 = t85 * rSges(5,1) + t84 * rSges(5,2) + rSges(5,3) * t170;
t63 = t112 * t132 + t129 * t174;
t64 = -t113 * t132 + t129 * t175;
t27 = Icges(6,5) * t64 + Icges(6,6) * t63 + Icges(6,3) * t171;
t29 = Icges(6,4) * t64 + Icges(6,2) * t63 + Icges(6,6) * t171;
t31 = Icges(6,1) * t64 + Icges(6,4) * t63 + Icges(6,5) * t171;
t11 = t129 * t27 + (-t112 * t31 - t113 * t29) * t131;
t50 = Icges(6,3) * t129 + (-Icges(6,5) * t112 - Icges(6,6) * t113) * t131;
t51 = Icges(6,6) * t129 + (-Icges(6,4) * t112 - Icges(6,2) * t113) * t131;
t52 = Icges(6,5) * t129 + (-Icges(6,1) * t112 - Icges(6,4) * t113) * t131;
t14 = t50 * t171 + t51 * t63 + t52 * t64;
t162 = t11 / 0.2e1 + t14 / 0.2e1;
t26 = Icges(6,5) * t62 + Icges(6,6) * t61 + Icges(6,3) * t170;
t28 = Icges(6,4) * t62 + Icges(6,2) * t61 + Icges(6,6) * t170;
t30 = Icges(6,1) * t62 + Icges(6,4) * t61 + Icges(6,5) * t170;
t10 = t129 * t26 + (-t112 * t30 - t113 * t28) * t131;
t13 = t50 * t170 + t61 * t51 + t62 * t52;
t161 = t13 / 0.2e1 + t10 / 0.2e1;
t160 = -Icges(4,4) * t129 / 0.2e1 + Icges(3,5) * t196 + (-Icges(4,5) / 0.2e1 + Icges(3,6) / 0.2e1) * t131;
t119 = t132 * pkin(3);
t158 = t130 * (qJ(4) * t171 - t119) + t132 * t166 + t185;
t157 = -qJ(4) * t129 - t98;
t156 = t167 + t168;
t155 = t157 - rSges(5,3) * t129 - (-rSges(5,1) * t126 - rSges(5,2) * t127) * t131;
t86 = t126 * t132 + t129 * t172;
t87 = -t127 * t132 + t129 * t173;
t154 = -rSges(5,1) * t87 - rSges(5,2) * t86;
t153 = -t64 * rSges(6,1) - t63 * rSges(6,2);
t152 = rSges(3,1) * t131 - rSges(3,2) * t129;
t151 = -t112 * t52 - t113 * t51;
t33 = rSges(6,3) * t171 - t153;
t54 = rSges(6,3) * t129 + (-rSges(6,1) * t112 - rSges(6,2) * t113) * t131;
t19 = -t129 * t33 + t54 * t171;
t20 = t129 * t32 - t54 * t170;
t146 = t130 * t20 + t132 * t19;
t137 = -t192 * t129 + t131 * t186 + t157 - t54;
t24 = t137 * t130;
t25 = t137 * t132;
t145 = t130 * t24 + t132 * t25;
t35 = t155 * t130;
t36 = t155 * t132;
t144 = t130 * t35 + t132 * t36;
t143 = Icges(3,1) * t131 - t181;
t142 = -Icges(3,2) * t129 + t180;
t139 = -Icges(4,2) * t131 + t179;
t138 = Icges(4,3) * t129 - t178;
t136 = rSges(3,1) * t170 - rSges(3,2) * t176 + t130 * rSges(3,3);
t135 = t130 * rSges(4,1) - rSges(4,2) * t170 + rSges(4,3) * t176;
t118 = t132 * pkin(6);
t16 = t169 + t118 + (-pkin(1) + (-qJ(3) - t186) * t129 + (-rSges(6,3) - pkin(2) + t128) * t131) * t130 + t153;
t17 = t156 + t193;
t22 = t118 + t119 + (-t177 - pkin(1) + (-rSges(5,3) - pkin(2) - qJ(4)) * t131) * t130 + t154;
t23 = t156 + t163 + t166;
t133 = (t130 * t17 + t132 * t16) * t188 + (t130 * t23 + t132 * t22) * t189;
t124 = t131 ^ 2;
t122 = t129 ^ 2;
t102 = rSges(2,1) * t132 - t130 * rSges(2,2);
t101 = -t130 * rSges(2,1) - rSges(2,2) * t132;
t100 = rSges(3,1) * t129 + rSges(3,2) * t131;
t59 = Icges(5,5) * t129 + (-Icges(5,1) * t126 - Icges(5,4) * t127) * t131;
t58 = Icges(5,6) * t129 + (-Icges(5,4) * t126 - Icges(5,2) * t127) * t131;
t49 = t184 * t132;
t48 = t184 * t130;
t47 = t129 * t50;
t46 = t136 + t167;
t45 = t182 + t118 + (-pkin(1) - t152) * t130;
t44 = t135 + t156;
t43 = t183 + t118 + (-pkin(1) + (rSges(4,2) - pkin(2)) * t131 + (-rSges(4,3) - qJ(3)) * t129) * t130;
t42 = Icges(5,1) * t87 + Icges(5,4) * t86 + Icges(5,5) * t171;
t41 = Icges(5,1) * t85 + Icges(5,4) * t84 + Icges(5,5) * t170;
t40 = Icges(5,4) * t87 + Icges(5,2) * t86 + Icges(5,6) * t171;
t39 = Icges(5,4) * t85 + Icges(5,2) * t84 + Icges(5,6) * t170;
t38 = Icges(5,5) * t87 + Icges(5,6) * t86 + Icges(5,3) * t171;
t37 = Icges(5,5) * t85 + Icges(5,6) * t84 + Icges(5,3) * t170;
t34 = t132 * t136 + (t152 * t130 - t182) * t130;
t21 = t132 * t135 + (-t183 + (-rSges(4,2) * t131 + rSges(4,3) * t129) * t130) * t130 + t185;
t18 = (t151 * t131 + t47) * t129;
t15 = (-t130 * t32 + t132 * t33) * t131;
t12 = t130 * (rSges(5,3) * t171 - t154) + t132 * t163 + t158;
t9 = t27 * t171 + t29 * t63 + t31 * t64;
t8 = t26 * t171 + t28 * t63 + t30 * t64;
t7 = t27 * t170 + t61 * t29 + t62 * t31;
t6 = t26 * t170 + t61 * t28 + t62 * t30;
t5 = (-t166 + t193) * t132 + (-t169 + t119 + t33 + (t129 * t186 + t192 * t131) * t130) * t130 + t158;
t4 = t8 * t130 - t132 * t9;
t3 = t6 * t130 - t132 * t7;
t2 = t14 * t129 + (t130 * t9 + t132 * t8) * t131;
t1 = t13 * t129 + (t130 * t7 + t132 * t6) * t131;
t53 = [Icges(2,3) + t47 + m(6) * (t16 ^ 2 + t17 ^ 2) + m(4) * (t43 ^ 2 + t44 ^ 2) + m(5) * (t22 ^ 2 + t23 ^ 2) + m(3) * (t45 ^ 2 + t46 ^ 2) + m(2) * (t101 ^ 2 + t102 ^ 2) + (-t126 * t59 - t127 * t58 + t151 + t179 + t181 + (Icges(3,2) + Icges(4,3)) * t131) * t131 + (t180 + t178 + (-Icges(5,5) * t126 - Icges(5,6) * t127) * t131 + (Icges(3,1) + Icges(4,2) + Icges(5,3)) * t129) * t129; (-t86 * t58 / 0.2e1 - t87 * t59 / 0.2e1 + t160 * t132 - t162) * t132 + (t84 * t58 / 0.2e1 + t85 * t59 / 0.2e1 + t160 * t130 + t161) * t130 + m(6) * (t16 * t25 + t17 * t24) + m(4) * (t43 * t49 + t44 * t48) + m(5) * (t22 * t36 + t23 * t35) + m(3) * (-t130 * t46 - t132 * t45) * t100 + ((Icges(3,5) * t197 + t143 * t199 + Icges(4,4) * t198 + t139 * t187 - t38 / 0.2e1) * t132 + (Icges(3,5) * t187 + t143 * t197 + Icges(4,4) * t199 + t139 * t198 + t37 / 0.2e1) * t130) * t129 + ((Icges(3,6) * t197 + t142 * t199 + Icges(4,5) * t198 + t138 * t187 + t126 * t42 / 0.2e1 + t127 * t40 / 0.2e1) * t132 + (Icges(3,6) * t187 + t142 * t197 + Icges(4,5) * t199 + t138 * t198 - t126 * t41 / 0.2e1 - t127 * t39 / 0.2e1) * t130) * t131; m(6) * (t24 ^ 2 + t25 ^ 2 + t5 ^ 2) + m(5) * (t12 ^ 2 + t35 ^ 2 + t36 ^ 2) + m(4) * (t21 ^ 2 + t48 ^ 2 + t49 ^ 2) + m(3) * (t165 * t100 ^ 2 + t34 ^ 2) + (-t4 + (t38 * t171 + t86 * t40 + t87 * t42) * t132 + t194 * t125) * t132 + (t3 + (t37 * t170 + t84 * t39 + t85 * t41) * t130 + t195 * t123 + (t194 * t130 + t195 * t132 - t38 * t170 - t37 * t171 - t86 * t39 - t84 * t40 - t87 * t41 - t85 * t42) * t132) * t130; 0.2e1 * ((t130 * t44 + t132 * t43) * t190 + t133) * t129; m(6) * (t145 * t129 - t131 * t5) + m(5) * (-t131 * t12 + t144 * t129) + m(4) * (-t131 * t21 + (t130 * t48 + t132 * t49) * t129); 0.2e1 * (t190 + t164) * (t165 * t122 + t124); t133 * t191; m(6) * (t129 * t5 + t145 * t131) + m(5) * (t129 * t12 + t144 * t131); t164 * (-0.1e1 + t165) * t129 * t191; 0.2e1 * t164 * (t165 * t124 + t122); m(6) * (t16 * t19 + t17 * t20) + t18 + (t162 * t130 + t161 * t132) * t131; m(6) * (t15 * t5 + t19 * t25 + t20 * t24) + t2 * t198 + t1 * t187 + (t10 * t130 - t11 * t132) * t196 + (t4 * t187 + t3 * t197) * t131; m(6) * (t146 * t129 - t15 * t131); m(6) * (t15 * t129 + t146 * t131); m(6) * (t15 ^ 2 + t19 ^ 2 + t20 ^ 2) + t129 * t18 + (t132 * t1 + t130 * t2 + t129 * (t10 * t132 + t11 * t130)) * t131;];
%% Postprocessing: Reshape Output
% From vec2symmat_5_matlab.m
res = [t53(1), t53(2), t53(4), t53(7), t53(11); t53(2), t53(3), t53(5), t53(8), t53(12); t53(4), t53(5), t53(6), t53(9), t53(13); t53(7), t53(8), t53(9), t53(10), t53(14); t53(11), t53(12), t53(13), t53(14), t53(15);];
Mq = res;
