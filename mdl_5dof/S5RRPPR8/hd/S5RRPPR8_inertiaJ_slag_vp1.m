% Calculate joint inertia matrix for
% S5RRPPR8
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
% Datum: 2019-12-31 19:39
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Mq = S5RRPPR8_inertiaJ_slag_vp1(qJ, ...
  pkin, m, rSges, Icges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(8,1),zeros(6,1),zeros(6,3),zeros(6,6)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RRPPR8_inertiaJ_slag_vp1: qJ has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S5RRPPR8_inertiaJ_slag_vp1: pkin has to be [8x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5RRPPR8_inertiaJ_slag_vp1: m has to be [6x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [6,3]), ...
  'S5RRPPR8_inertiaJ_slag_vp1: rSges has to be [6x3] (double)');
assert(isreal(Icges) && all(size(Icges) == [6 6]), ...
  'S5RRPPR8_inertiaJ_slag_vp1: Icges has to be [6x6] (double)'); 

%% Symbolic Calculation
% From inertia_joint_joint_floatb_twist_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 19:38:08
% EndTime: 2019-12-31 19:38:13
% DurationCPUTime: 1.79s
% Computational Cost: add. (1884->262), mult. (3138->391), div. (0->0), fcn. (3299->8), ass. (0->132)
t201 = Icges(3,1) + Icges(4,1);
t199 = Icges(3,5) + Icges(4,4);
t119 = sin(qJ(2));
t200 = (Icges(3,4) - Icges(4,5)) * t119;
t198 = Icges(4,2) + Icges(3,3);
t121 = cos(qJ(2));
t197 = t199 * t121 + (-Icges(3,6) + Icges(4,6)) * t119;
t196 = t201 * t121 - t200;
t120 = sin(qJ(1));
t114 = t120 ^ 2;
t112 = pkin(8) + qJ(5);
t105 = sin(t112);
t106 = cos(t112);
t76 = -t105 * t121 + t106 * t119;
t195 = t76 / 0.2e1;
t194 = -t120 / 0.2e1;
t193 = t120 / 0.2e1;
t122 = cos(qJ(1));
t192 = -t122 / 0.2e1;
t191 = t122 / 0.2e1;
t128 = t105 * t119 + t106 * t121;
t190 = -t128 / 0.2e1;
t189 = Icges(6,5) * t195 + Icges(6,6) * t190;
t117 = cos(pkin(8));
t116 = sin(pkin(8));
t161 = t116 * t119;
t127 = t117 * t121 + t161;
t188 = t127 / 0.2e1;
t160 = t116 * t121;
t81 = t117 * t119 - t160;
t187 = -t81 / 0.2e1;
t185 = -t197 * t120 + t198 * t122;
t184 = t198 * t120 + t197 * t122;
t103 = pkin(4) * t117 + pkin(3);
t118 = -pkin(7) - qJ(4);
t153 = pkin(4) * t161;
t158 = t121 * t122;
t52 = t76 * t122;
t53 = t128 * t122;
t25 = t53 * rSges(6,1) + t52 * rSges(6,2) - rSges(6,3) * t120;
t183 = t103 * t158 + t120 * t118 + t122 * t153 + t25;
t115 = t122 ^ 2;
t182 = m(4) / 0.2e1;
t181 = m(5) / 0.2e1;
t93 = rSges(3,1) * t119 + rSges(3,2) * t121;
t180 = m(3) * t93;
t179 = pkin(3) * t121;
t71 = t81 * t122;
t72 = t127 * t122;
t178 = t72 * rSges(5,1) + t71 * rSges(5,2);
t162 = qJ(3) * t119;
t159 = t119 * t122;
t168 = pkin(2) * t158 + qJ(3) * t159;
t177 = t114 * (pkin(2) * t121 + t162) + t122 * t168;
t91 = pkin(2) * t119 - qJ(3) * t121;
t176 = -rSges(4,1) * t119 + rSges(4,3) * t121 - t91;
t51 = t128 * t120;
t175 = Icges(6,4) * t51;
t174 = Icges(6,5) * t51;
t50 = t76 * t120;
t173 = Icges(6,2) * t50;
t172 = Icges(6,6) * t50;
t171 = t122 * rSges(4,2);
t170 = t122 * rSges(3,3);
t169 = -rSges(5,3) - qJ(4);
t166 = Icges(3,4) * t121;
t164 = Icges(4,5) * t121;
t163 = Icges(6,3) * t122;
t157 = t122 * pkin(1) + t120 * pkin(6);
t156 = t114 + t115;
t155 = t181 + m(6) / 0.2e1;
t123 = Icges(6,6) * t122 + t173 + t175;
t124 = Icges(6,1) * t51 + Icges(6,4) * t50 + Icges(6,5) * t122;
t37 = Icges(6,4) * t76 - Icges(6,2) * t128;
t38 = Icges(6,1) * t76 - Icges(6,4) * t128;
t154 = t122 * t189 + t50 * t37 / 0.2e1 + t51 * t38 / 0.2e1 + t123 * t190 + t124 * t195;
t22 = Icges(6,4) * t53 + Icges(6,2) * t52 - Icges(6,6) * t120;
t23 = Icges(6,1) * t53 + Icges(6,4) * t52 - Icges(6,5) * t120;
t152 = t120 * t189 - t37 * t52 / 0.2e1 - t38 * t53 / 0.2e1 + t22 * t128 / 0.2e1 - t23 * t76 / 0.2e1;
t150 = rSges(4,1) * t158 + t120 * rSges(4,2) + rSges(4,3) * t159;
t149 = -pkin(3) * t119 - t91;
t101 = pkin(3) * t158;
t148 = t122 * t101 + t114 * t179 + t177;
t147 = t157 + t168;
t146 = -rSges(5,1) * t81 + rSges(5,2) * t127 + t149;
t145 = Icges(5,5) * t187 + Icges(5,6) * t188 + t199 * t119 / 0.2e1 + (Icges(3,6) / 0.2e1 - Icges(4,6) / 0.2e1) * t121;
t69 = t81 * t120;
t70 = t127 * t120;
t144 = t70 * rSges(5,1) + t69 * rSges(5,2);
t143 = -t51 * rSges(6,1) - t50 * rSges(6,2);
t21 = Icges(6,5) * t53 + Icges(6,6) * t52 - Icges(6,3) * t120;
t142 = t122 * ((t122 * t21 + t50 * t22 + t51 * t23) * t120 - (Icges(6,1) * t51 ^ 2 + (t173 + 0.2e1 * t175) * t50 + (t163 + 0.2e1 * t172 + 0.2e1 * t174) * t122) * t122) - t120 * ((-t120 * t21 + t22 * t52 + t23 * t53) * t120 - (t53 * t124 + t52 * t123 - t120 * (t163 + t172 + t174)) * t122);
t141 = rSges(3,1) * t121 - rSges(3,2) * t119;
t39 = rSges(6,1) * t76 - rSges(6,2) * t128;
t135 = t149 - t39 + pkin(4) * t160 - (-pkin(3) + t103) * t119;
t14 = t135 * t120;
t15 = t135 * t122;
t136 = t120 * t14 + t122 * t15;
t132 = -Icges(3,2) * t119 + t166;
t129 = Icges(4,3) * t119 + t164;
t126 = rSges(3,1) * t158 - rSges(3,2) * t159 + t120 * rSges(3,3);
t110 = t122 * pkin(6);
t12 = t110 + (-rSges(6,3) + t118) * t122 + (-pkin(1) + (-pkin(2) - t103) * t121 + (-pkin(4) * t116 - qJ(3)) * t119) * t120 + t143;
t13 = t147 + t183;
t125 = m(6) * (t12 * t122 + t120 * t13);
t95 = rSges(2,1) * t122 - t120 * rSges(2,2);
t94 = -t120 * rSges(2,1) - rSges(2,2) * t122;
t47 = t176 * t122;
t46 = t176 * t120;
t45 = t126 + t157;
t44 = t170 + t110 + (-pkin(1) - t141) * t120;
t42 = Icges(5,1) * t81 - Icges(5,4) * t127;
t41 = Icges(5,4) * t81 - Icges(5,2) * t127;
t34 = t147 + t150;
t33 = t171 + t110 + (-pkin(1) + (-rSges(4,1) - pkin(2)) * t121 + (-rSges(4,3) - qJ(3)) * t119) * t120;
t32 = t122 * t126 + (t120 * t141 - t170) * t120;
t31 = Icges(5,1) * t72 + Icges(5,4) * t71 - Icges(5,5) * t120;
t30 = Icges(5,1) * t70 + Icges(5,4) * t69 + Icges(5,5) * t122;
t29 = Icges(5,4) * t72 + Icges(5,2) * t71 - Icges(5,6) * t120;
t28 = Icges(5,4) * t70 + Icges(5,2) * t69 + Icges(5,6) * t122;
t27 = Icges(5,5) * t72 + Icges(5,6) * t71 - Icges(5,3) * t120;
t26 = Icges(5,5) * t70 + Icges(5,6) * t69 + Icges(5,3) * t122;
t24 = rSges(6,3) * t122 - t143;
t20 = t146 * t122;
t19 = t146 * t120;
t18 = t120 * t169 + t101 + t147 + t178;
t17 = t110 + t169 * t122 + (-t162 - pkin(1) + (-pkin(2) - pkin(3)) * t121) * t120 - t144;
t16 = t122 * t150 + (-t171 + (rSges(4,1) * t121 + rSges(4,3) * t119) * t120) * t120 + t177;
t11 = -t120 * t24 - t122 * t25;
t8 = t120 * t144 + t122 * t178 + t148;
t3 = (-t101 + t183) * t122 + (-t122 * t118 + t24 + (t103 * t121 + t153 - t179) * t120) * t120 + t148;
t1 = [-t128 * t37 + t76 * t38 - t127 * t41 + t81 * t42 + Icges(2,3) + m(6) * (t12 ^ 2 + t13 ^ 2) + m(4) * (t33 ^ 2 + t34 ^ 2) + m(5) * (t17 ^ 2 + t18 ^ 2) + m(3) * (t44 ^ 2 + t45 ^ 2) + m(2) * (t94 ^ 2 + t95 ^ 2) + ((Icges(3,2) + Icges(4,3)) * t121 + t200) * t121 + (t201 * t119 - t164 + t166) * t119; m(6) * (t12 * t15 + t13 * t14) + m(4) * (t33 * t47 + t34 * t46) + m(5) * (t17 * t20 + t18 * t19) + (-t44 * t180 - t69 * t41 / 0.2e1 - t70 * t42 / 0.2e1 + t28 * t188 + t30 * t187 + t145 * t122 + (Icges(3,6) * t191 + Icges(4,6) * t192 + t129 * t193 + t132 * t194) * t121 - t154) * t122 + (-t45 * t180 + t71 * t41 / 0.2e1 + t72 * t42 / 0.2e1 - t127 * t29 / 0.2e1 + t81 * t31 / 0.2e1 + t145 * t120 + (Icges(3,6) * t193 + Icges(4,6) * t194 + t129 * t192 + t132 * t191) * t121 - t152) * t120 + (t196 * t122 * t194 + t199 * t120 * t193 + (t196 * t120 + t199 * t122) * t191) * t119; m(6) * (t14 ^ 2 + t15 ^ 2 + t3 ^ 2) + m(5) * (t19 ^ 2 + t20 ^ 2 + t8 ^ 2) + m(4) * (t16 ^ 2 + t46 ^ 2 + t47 ^ 2) + m(3) * (t156 * t93 ^ 2 + t32 ^ 2) - t142 + ((t122 * t26 + t69 * t28 + t70 * t30) * t122 + t185 * t115) * t122 + ((-t120 * t27 + t71 * t29 + t72 * t31) * t120 + t184 * t114 + (-t69 * t29 - t70 * t31 - t71 * t28 - t72 * t30 + (t26 + t185) * t120 + (-t27 + t184) * t122) * t122) * t120; 0.2e1 * (t125 / 0.2e1 + (t120 * t34 + t122 * t33) * t182 + (t120 * t18 + t122 * t17) * t181) * t119; m(6) * (t119 * t136 - t121 * t3) + m(5) * (-t121 * t8 + (t120 * t19 + t122 * t20) * t119) + m(4) * (-t121 * t16 + (t120 * t46 + t122 * t47) * t119); 0.2e1 * (t182 + t155) * (t119 ^ 2 * t156 + t121 ^ 2); m(6) * (-t120 * t12 + t122 * t13) + m(5) * (-t120 * t17 + t122 * t18); m(6) * (-t120 * t15 + t122 * t14) + m(5) * (-t120 * t20 + t122 * t19); 0; 0.2e1 * t155 * t156; t120 * t152 + t122 * t154 + t125 * t39; m(6) * (t11 * t3 + t136 * t39) + t142; m(6) * (t119 * t156 * t39 - t11 * t121); 0; m(6) * (t156 * t39 ^ 2 + t11 ^ 2) - t142;];
%% Postprocessing: Reshape Output
% From vec2symmat_5_matlab.m
res = [t1(1), t1(2), t1(4), t1(7), t1(11); t1(2), t1(3), t1(5), t1(8), t1(12); t1(4), t1(5), t1(6), t1(9), t1(13); t1(7), t1(8), t1(9), t1(10), t1(14); t1(11), t1(12), t1(13), t1(14), t1(15);];
Mq = res;
