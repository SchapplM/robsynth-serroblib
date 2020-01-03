% Calculate joint inertia matrix for
% S5RRPPR10
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% pkin [8x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d1,d2,d5,theta3]';
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
% Datum: 2019-12-31 19:45
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Mq = S5RRPPR10_inertiaJ_slag_vp1(qJ, ...
  pkin, m, rSges, Icges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(8,1),zeros(6,1),zeros(6,3),zeros(6,6)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RRPPR10_inertiaJ_slag_vp1: qJ has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S5RRPPR10_inertiaJ_slag_vp1: pkin has to be [8x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5RRPPR10_inertiaJ_slag_vp1: m has to be [6x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [6,3]), ...
  'S5RRPPR10_inertiaJ_slag_vp1: rSges has to be [6x3] (double)');
assert(isreal(Icges) && all(size(Icges) == [6 6]), ...
  'S5RRPPR10_inertiaJ_slag_vp1: Icges has to be [6x6] (double)'); 

%% Symbolic Calculation
% From inertia_joint_joint_floatb_twist_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 19:43:03
% EndTime: 2019-12-31 19:43:09
% DurationCPUTime: 2.14s
% Computational Cost: add. (2355->312), mult. (5814->462), div. (0->0), fcn. (6751->8), ass. (0->146)
t194 = Icges(4,1) + Icges(5,1);
t193 = -Icges(4,4) + Icges(5,5);
t192 = Icges(5,4) + Icges(4,5);
t191 = Icges(4,2) + Icges(5,3);
t190 = Icges(4,6) - Icges(5,6);
t129 = sin(qJ(1));
t189 = -t129 / 0.2e1;
t177 = t129 / 0.2e1;
t132 = cos(qJ(1));
t188 = t132 / 0.2e1;
t128 = sin(qJ(2));
t164 = t128 * t129;
t125 = sin(pkin(8));
t126 = cos(pkin(8));
t131 = cos(qJ(2));
t162 = t129 * t131;
t98 = t125 * t162 + t126 * t132;
t99 = -t125 * t132 + t126 * t162;
t185 = -t164 * t190 + t191 * t98 + t193 * t99;
t49 = Icges(4,5) * t99 - Icges(4,6) * t98 + Icges(4,3) * t164;
t51 = Icges(5,4) * t99 + Icges(5,2) * t164 + Icges(5,6) * t98;
t187 = t49 + t51;
t161 = t131 * t132;
t100 = t125 * t161 - t126 * t129;
t101 = t125 * t129 + t126 * t161;
t163 = t128 * t132;
t50 = Icges(4,5) * t101 - Icges(4,6) * t100 + Icges(4,3) * t163;
t52 = Icges(5,4) * t101 + Icges(5,2) * t163 + Icges(5,6) * t100;
t186 = t50 + t52;
t184 = t100 * t191 + t101 * t193 - t163 * t190;
t183 = t164 * t192 + t193 * t98 + t194 * t99;
t182 = t100 * t193 + t101 * t194 + t163 * t192;
t123 = t129 ^ 2;
t124 = t132 ^ 2;
t181 = 0.2e1 * t128;
t180 = m(4) / 0.2e1;
t179 = m(5) / 0.2e1;
t178 = m(6) / 0.2e1;
t176 = -t132 / 0.2e1;
t175 = -rSges(6,3) - pkin(7);
t174 = rSges(5,3) * t98;
t173 = pkin(2) * t131;
t127 = sin(qJ(5));
t130 = cos(qJ(5));
t66 = t100 * t130 - t101 * t127;
t67 = t100 * t127 + t101 * t130;
t172 = rSges(6,1) * t67 + rSges(6,2) * t66;
t159 = pkin(2) * t161 + qJ(3) * t163;
t171 = t123 * (qJ(3) * t128 + t173) + t132 * t159;
t170 = t132 * rSges(3,3);
t108 = pkin(2) * t128 - qJ(3) * t131;
t169 = -t108 + rSges(4,3) * t131 - (rSges(4,1) * t126 - rSges(4,2) * t125) * t128;
t120 = t132 * pkin(6);
t90 = t98 * qJ(4);
t168 = t120 - t90;
t167 = Icges(3,4) * t128;
t166 = Icges(3,4) * t131;
t165 = t125 * t128;
t160 = -(pkin(3) * t126 + qJ(4) * t125) * t128 - t108;
t158 = t132 * pkin(1) + t129 * pkin(6);
t157 = t123 + t124;
t156 = t179 + t178;
t86 = (t125 * t130 - t126 * t127) * t128;
t87 = (t125 * t127 + t126 * t130) * t128;
t42 = Icges(6,5) * t87 + Icges(6,6) * t86 + Icges(6,3) * t131;
t43 = Icges(6,4) * t87 + Icges(6,2) * t86 + Icges(6,6) * t131;
t44 = Icges(6,1) * t87 + Icges(6,4) * t86 + Icges(6,5) * t131;
t155 = t131 * t42 + t43 * t86 + t44 * t87;
t154 = t101 * rSges(5,1) + rSges(5,2) * t163 + t100 * rSges(5,3);
t153 = t101 * rSges(4,1) - t100 * rSges(4,2) + rSges(4,3) * t163;
t64 = -t127 * t99 + t130 * t98;
t65 = t127 * t98 + t130 * t99;
t26 = Icges(6,5) * t65 + Icges(6,6) * t64 - Icges(6,3) * t164;
t28 = Icges(6,4) * t65 + Icges(6,2) * t64 - Icges(6,6) * t164;
t30 = Icges(6,1) * t65 + Icges(6,4) * t64 - Icges(6,5) * t164;
t10 = t131 * t26 + t28 * t86 + t30 * t87;
t12 = -t164 * t42 + t43 * t64 + t44 * t65;
t152 = -t12 / 0.2e1 - t10 / 0.2e1;
t27 = Icges(6,5) * t67 + Icges(6,6) * t66 - Icges(6,3) * t163;
t29 = Icges(6,4) * t67 + Icges(6,2) * t66 - Icges(6,6) * t163;
t31 = Icges(6,1) * t67 + Icges(6,4) * t66 - Icges(6,5) * t163;
t11 = t131 * t27 + t29 * t86 + t31 * t87;
t13 = -t163 * t42 + t43 * t66 + t44 * t67;
t151 = -t13 / 0.2e1 - t11 / 0.2e1;
t72 = -Icges(5,6) * t131 + (Icges(5,5) * t126 + Icges(5,3) * t125) * t128;
t75 = -Icges(4,6) * t131 + (Icges(4,4) * t126 - Icges(4,2) * t125) * t128;
t150 = t75 / 0.2e1 - t72 / 0.2e1;
t76 = -Icges(5,4) * t131 + (Icges(5,1) * t126 + Icges(5,5) * t125) * t128;
t77 = -Icges(4,5) * t131 + (Icges(4,1) * t126 - Icges(4,4) * t125) * t128;
t149 = t77 / 0.2e1 + t76 / 0.2e1;
t148 = rSges(5,2) * t131 - (rSges(5,1) * t126 + rSges(5,3) * t125) * t128 + t160;
t147 = -pkin(1) - t173;
t145 = t101 * pkin(3) + t100 * qJ(4);
t146 = t129 * (pkin(3) * t99 + t90) + t132 * t145 + t171;
t45 = rSges(6,1) * t87 + rSges(6,2) * t86 + rSges(6,3) * t131;
t144 = -pkin(4) * t126 * t128 - pkin(7) * t131 + t160 - t45;
t143 = t158 + t159;
t142 = -rSges(4,1) * t99 + rSges(4,2) * t98;
t141 = -rSges(6,1) * t65 - rSges(6,2) * t64;
t140 = rSges(3,1) * t131 - rSges(3,2) * t128;
t137 = Icges(3,1) * t131 - t167;
t136 = -Icges(3,2) * t128 + t166;
t135 = Icges(3,5) * t131 - Icges(3,6) * t128;
t134 = rSges(3,1) * t161 - rSges(3,2) * t163 + t129 * rSges(3,3);
t133 = t143 + t145;
t122 = t128 ^ 2;
t111 = rSges(2,1) * t132 - rSges(2,2) * t129;
t110 = -rSges(2,1) * t129 - rSges(2,2) * t132;
t109 = rSges(3,1) * t128 + rSges(3,2) * t131;
t105 = Icges(3,5) * t128 + Icges(3,6) * t131;
t95 = t101 * pkin(4);
t81 = Icges(3,3) * t129 + t132 * t135;
t80 = -Icges(3,3) * t132 + t129 * t135;
t71 = t134 + t158;
t70 = t170 + t120 + (-pkin(1) - t140) * t129;
t69 = t169 * t132;
t68 = t169 * t129;
t46 = t132 * t134 + (t129 * t140 - t170) * t129;
t39 = t148 * t132;
t38 = t148 * t129;
t35 = t143 + t153;
t34 = t120 + ((-rSges(4,3) - qJ(3)) * t128 + t147) * t129 + t142;
t33 = -rSges(6,3) * t163 + t172;
t32 = -rSges(6,3) * t164 - t141;
t25 = t144 * t132;
t24 = t144 * t129;
t23 = t133 + t154;
t22 = -t174 + (-rSges(5,1) - pkin(3)) * t99 + ((-rSges(5,2) - qJ(3)) * t128 + t147) * t129 + t168;
t21 = t129 * (rSges(4,3) * t164 - t142) + t132 * t153 + t171;
t20 = t131 * t33 + t163 * t45;
t19 = -t131 * t32 - t164 * t45;
t18 = t163 * t175 + t133 + t172 + t95;
t17 = (-pkin(3) - pkin(4)) * t99 + ((-qJ(3) - t175) * t128 + t147) * t129 + t141 + t168;
t16 = t155 * t131;
t15 = (t129 * t33 - t132 * t32) * t128;
t14 = t129 * (rSges(5,1) * t99 + rSges(5,2) * t164 + t174) + t132 * t154 + t146;
t9 = -t163 * t27 + t29 * t66 + t31 * t67;
t8 = -t163 * t26 + t28 * t66 + t30 * t67;
t7 = -t164 * t27 + t29 * t64 + t31 * t65;
t6 = -t164 * t26 + t28 * t64 + t30 * t65;
t5 = (-pkin(7) * t163 + t33 + t95) * t132 + (pkin(4) * t99 - pkin(7) * t164 + t32) * t129 + t146;
t4 = t129 * t9 - t132 * t8;
t3 = t129 * t7 - t132 * t6;
t2 = t13 * t131 + (-t129 * t8 - t132 * t9) * t128;
t1 = t12 * t131 + (-t129 * t6 - t132 * t7) * t128;
t36 = [Icges(2,3) + (t167 + (t125 * t190 - t126 * t192) * t128 + (Icges(4,3) + Icges(5,2) + Icges(3,2)) * t131) * t131 + (Icges(3,1) * t128 + t166 + (t76 + t77) * t126 + (t72 - t75) * t125) * t128 + m(6) * (t17 ^ 2 + t18 ^ 2) + m(5) * (t22 ^ 2 + t23 ^ 2) + m(4) * (t34 ^ 2 + t35 ^ 2) + m(3) * (t70 ^ 2 + t71 ^ 2) + m(2) * (t110 ^ 2 + t111 ^ 2) + t155; (t105 * t188 - t149 * t99 + t150 * t98 + t152) * t132 + (-t150 * t100 + t149 * t101 + t105 * t177 - t151) * t129 + m(6) * (t17 * t25 + t18 * t24) + m(5) * (t22 * t39 + t23 * t38) + m(4) * (t34 * t69 + t35 * t68) + m(3) * (-t129 * t71 - t132 * t70) * t109 + ((t49 / 0.2e1 + t51 / 0.2e1 + Icges(3,6) * t188 + t136 * t189) * t132 + (Icges(3,6) * t177 + t136 * t188 - t50 / 0.2e1 - t52 / 0.2e1) * t129) * t131 + ((Icges(3,5) * t129 + t125 * t184 + t126 * t182 + t132 * t137) * t177 + (-Icges(3,5) * t132 + t125 * t185 + t126 * t183 + t129 * t137) * t176) * t128; m(6) * (t24 ^ 2 + t25 ^ 2 + t5 ^ 2) + m(5) * (t14 ^ 2 + t38 ^ 2 + t39 ^ 2) + m(4) * (t21 ^ 2 + t68 ^ 2 + t69 ^ 2) + m(3) * (t109 ^ 2 * t157 + t46 ^ 2) + (-t124 * t80 - t3 + (t164 * t187 + t183 * t99 + t185 * t98) * t132) * t132 + (t4 + t123 * t81 + (t100 * t184 + t101 * t182 + t163 * t186) * t129 + (-t185 * t100 - t183 * t101 - t129 * t80 + t132 * t81 - t163 * t187 - t164 * t186 - t182 * t99 - t184 * t98) * t132) * t129; ((t129 * t18 + t132 * t17) * t178 + (t129 * t23 + t132 * t22) * t179 + (t129 * t35 + t132 * t34) * t180) * t181; m(6) * (-t131 * t5 + (t129 * t24 + t132 * t25) * t128) + m(5) * (-t131 * t14 + (t129 * t38 + t132 * t39) * t128) + m(4) * (-t131 * t21 + (t129 * t68 + t132 * t69) * t128); 0.2e1 * (t180 + t156) * (t122 * t157 + t131 ^ 2); m(6) * (t100 * t17 + t18 * t98) + m(5) * (t100 * t22 + t23 * t98); m(6) * (t100 * t25 + t165 * t5 + t24 * t98) + m(5) * (t100 * t39 + t14 * t165 + t38 * t98); t156 * (t100 * t132 - t125 * t131 + t129 * t98) * t181; 0.2e1 * t156 * (t122 * t125 ^ 2 + t100 ^ 2 + t98 ^ 2); m(6) * (t17 * t19 + t18 * t20) + t16 + (t129 * t152 + t132 * t151) * t128; m(6) * (t15 * t5 + t19 * t25 + t20 * t24) + t2 * t177 + t1 * t176 + t131 * (-t10 * t132 + t11 * t129) / 0.2e1 + (t4 * t176 + t189 * t3) * t128; m(6) * (-t15 * t131 + (t129 * t20 + t132 * t19) * t128); m(6) * (t100 * t19 + t15 * t165 + t20 * t98); m(6) * (t15 ^ 2 + t19 ^ 2 + t20 ^ 2) + t131 * t16 + (-t132 * t2 - t129 * t1 + t131 * (-t10 * t129 - t11 * t132)) * t128;];
%% Postprocessing: Reshape Output
% From vec2symmat_5_matlab.m
res = [t36(1), t36(2), t36(4), t36(7), t36(11); t36(2), t36(3), t36(5), t36(8), t36(12); t36(4), t36(5), t36(6), t36(9), t36(13); t36(7), t36(8), t36(9), t36(10), t36(14); t36(11), t36(12), t36(13), t36(14), t36(15);];
Mq = res;
