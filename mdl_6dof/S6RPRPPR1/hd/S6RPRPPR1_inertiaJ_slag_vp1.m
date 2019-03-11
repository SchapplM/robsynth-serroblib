% Calculate joint inertia matrix for
% S6RPRPPR1
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% pkin [11x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d3,d6,theta2,theta4,theta5]';
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
% Datum: 2019-03-09 02:40
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Mq = S6RPRPPR1_inertiaJ_slag_vp1(qJ, ...
  pkin, m, rSges, Icges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(11,1),zeros(7,1),zeros(7,3),zeros(7,6)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPRPPR1_inertiaJ_slag_vp1: qJ has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [11 1]), ...
  'S6RPRPPR1_inertiaJ_slag_vp1: pkin has to be [11x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6RPRPPR1_inertiaJ_slag_vp1: m has to be [7x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [7,3]), ...
  'S6RPRPPR1_inertiaJ_slag_vp1: rSges has to be [7x3] (double)');
assert(isreal(Icges) && all(size(Icges) == [7 6]), ...
  'S6RPRPPR1_inertiaJ_slag_vp1: Icges has to be [7x6] (double)'); 

%% Symbolic Calculation
% From inertia_joint_joint_floatb_twist_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 02:38:00
% EndTime: 2019-03-09 02:38:05
% DurationCPUTime: 1.99s
% Computational Cost: add. (4987->329), mult. (3643->478), div. (0->0), fcn. (3733->12), ass. (0->155)
t206 = Icges(4,3) + Icges(5,3);
t128 = qJ(3) + pkin(10);
t121 = sin(t128);
t124 = cos(t128);
t134 = sin(qJ(3));
t136 = cos(qJ(3));
t205 = Icges(4,5) * t136 + Icges(5,5) * t124 - Icges(4,6) * t134 - Icges(5,6) * t121;
t129 = qJ(1) + pkin(9);
t122 = sin(t129);
t194 = t122 / 0.2e1;
t125 = cos(t129);
t204 = t125 / 0.2e1;
t203 = t134 / 0.2e1;
t202 = t136 / 0.2e1;
t201 = -t205 * t122 + t206 * t125;
t200 = t206 * t122 + t205 * t125;
t131 = cos(pkin(11));
t115 = pkin(5) * t131 + pkin(4);
t171 = t124 * t125;
t130 = sin(pkin(11));
t173 = t122 * t130;
t133 = -pkin(8) - qJ(5);
t175 = t121 * t133;
t176 = t121 * t125;
t127 = pkin(11) + qJ(6);
t120 = sin(t127);
t123 = cos(t127);
t74 = -t120 * t171 + t122 * t123;
t75 = t120 * t122 + t123 * t171;
t33 = t75 * rSges(7,1) + t74 * rSges(7,2) + rSges(7,3) * t176;
t199 = pkin(5) * t173 + t115 * t171 - t125 * t175 + t33;
t118 = t122 ^ 2;
t198 = t124 ^ 2;
t119 = t125 ^ 2;
t197 = m(6) / 0.2e1;
t196 = m(7) / 0.2e1;
t195 = -m(6) - m(7);
t193 = -t125 / 0.2e1;
t107 = rSges(4,1) * t134 + rSges(4,2) * t136;
t192 = m(4) * t107;
t135 = sin(qJ(1));
t191 = pkin(1) * t135;
t190 = pkin(3) * t134;
t189 = pkin(4) * t124;
t116 = pkin(3) * t136 + pkin(2);
t102 = t125 * t116;
t114 = t125 * pkin(7);
t132 = -qJ(4) - pkin(7);
t168 = t125 * t132;
t188 = t122 * (t168 + t114 + (-pkin(2) + t116) * t122) + t125 * (-pkin(2) * t125 + t102 + (-pkin(7) - t132) * t122);
t187 = rSges(4,1) * t136;
t186 = rSges(4,2) * t134;
t54 = -Icges(7,6) * t124 + (Icges(7,4) * t123 - Icges(7,2) * t120) * t121;
t185 = t120 * t54;
t184 = t125 * rSges(4,3);
t183 = pkin(4) * t171 + qJ(5) * t176;
t182 = Icges(4,4) * t134;
t181 = Icges(4,4) * t136;
t180 = Icges(5,4) * t121;
t179 = Icges(5,4) * t124;
t178 = t115 * t124;
t177 = t121 * t122;
t174 = t122 * t124;
t172 = t122 * t131;
t170 = t125 * t130;
t169 = t125 * t131;
t167 = t122 * rSges(4,3) + t125 * t187;
t166 = t118 + t119;
t165 = t197 + t196;
t86 = -t124 * t170 + t172;
t87 = t124 * t169 + t173;
t164 = t87 * rSges(6,1) + t86 * rSges(6,2) + rSges(6,3) * t176;
t72 = -t120 * t174 - t123 * t125;
t73 = -t120 * t125 + t123 * t174;
t26 = Icges(7,5) * t73 + Icges(7,6) * t72 + Icges(7,3) * t177;
t28 = Icges(7,4) * t73 + Icges(7,2) * t72 + Icges(7,6) * t177;
t30 = Icges(7,1) * t73 + Icges(7,4) * t72 + Icges(7,5) * t177;
t11 = -t124 * t26 + (-t120 * t28 + t123 * t30) * t121;
t51 = -Icges(7,3) * t124 + (Icges(7,5) * t123 - Icges(7,6) * t120) * t121;
t57 = -Icges(7,5) * t124 + (Icges(7,1) * t123 - Icges(7,4) * t120) * t121;
t14 = t177 * t51 + t54 * t72 + t57 * t73;
t163 = t11 / 0.2e1 + t14 / 0.2e1;
t27 = Icges(7,5) * t75 + Icges(7,6) * t74 + Icges(7,3) * t176;
t29 = Icges(7,4) * t75 + Icges(7,2) * t74 + Icges(7,6) * t176;
t31 = Icges(7,1) * t75 + Icges(7,4) * t74 + Icges(7,5) * t176;
t12 = -t124 * t27 + (-t120 * t29 + t123 * t31) * t121;
t15 = t176 * t51 + t54 * t74 + t57 * t75;
t162 = t12 / 0.2e1 + t15 / 0.2e1;
t161 = Icges(4,5) * t203 + Icges(4,6) * t202 + Icges(5,5) * t121 / 0.2e1 + Icges(5,6) * t124 / 0.2e1;
t160 = -pkin(4) * t121 + qJ(5) * t124 - t190;
t159 = -rSges(5,1) * t121 - rSges(5,2) * t124 - t190;
t151 = qJ(5) * t121 + t189;
t158 = t118 * t151 + t125 * t183 + t188;
t157 = t160 + rSges(6,3) * t124 - (rSges(6,1) * t131 - rSges(6,2) * t130) * t121;
t84 = -t124 * t173 - t169;
t85 = t124 * t172 - t170;
t156 = -rSges(6,1) * t85 - rSges(6,2) * t84;
t155 = -rSges(7,1) * t73 - rSges(7,2) * t72;
t137 = cos(qJ(1));
t126 = t137 * pkin(1);
t154 = -t122 * t132 + t102 + t126;
t153 = -t186 + t187;
t152 = rSges(5,1) * t124 - rSges(5,2) * t121;
t60 = -rSges(7,3) * t124 + (rSges(7,1) * t123 - rSges(7,2) * t120) * t121;
t146 = t160 - (qJ(5) + t133) * t124 - (-pkin(4) + t115) * t121 - t60;
t145 = Icges(4,1) * t136 - t182;
t144 = Icges(5,1) * t124 - t180;
t143 = -Icges(4,2) * t134 + t181;
t142 = -Icges(5,2) * t121 + t179;
t139 = rSges(5,1) * t171 - rSges(5,2) * t176 + t122 * rSges(5,3);
t109 = rSges(2,1) * t137 - t135 * rSges(2,2);
t108 = -t135 * rSges(2,1) - rSges(2,2) * t137;
t89 = rSges(3,1) * t125 - rSges(3,2) * t122 + t126;
t88 = -rSges(3,1) * t122 - rSges(3,2) * t125 - t191;
t70 = -Icges(6,5) * t124 + (Icges(6,1) * t131 - Icges(6,4) * t130) * t121;
t69 = -Icges(6,6) * t124 + (Icges(6,4) * t131 - Icges(6,2) * t130) * t121;
t62 = t159 * t125;
t61 = t159 * t122;
t47 = t121 * t123 * t57;
t46 = pkin(7) * t122 + t126 + (pkin(2) - t186) * t125 + t167;
t45 = t184 - t191 + t114 + (-pkin(2) - t153) * t122;
t44 = t139 + t154;
t43 = -t191 + (rSges(5,3) - t132) * t125 + (-t116 - t152) * t122;
t42 = t157 * t125;
t41 = t157 * t122;
t40 = Icges(6,1) * t87 + Icges(6,4) * t86 + Icges(6,5) * t176;
t39 = Icges(6,1) * t85 + Icges(6,4) * t84 + Icges(6,5) * t177;
t38 = Icges(6,4) * t87 + Icges(6,2) * t86 + Icges(6,6) * t176;
t37 = Icges(6,4) * t85 + Icges(6,2) * t84 + Icges(6,6) * t177;
t36 = Icges(6,5) * t87 + Icges(6,6) * t86 + Icges(6,3) * t176;
t35 = Icges(6,5) * t85 + Icges(6,6) * t84 + Icges(6,3) * t177;
t34 = t125 * (-t125 * t186 + t167) + (t122 * t153 - t184) * t122;
t32 = rSges(7,3) * t177 - t155;
t25 = t146 * t125;
t24 = t146 * t122;
t23 = t154 + t164 + t183;
t22 = -t191 - t168 + (-t189 - t116 + (-rSges(6,3) - qJ(5)) * t121) * t122 + t156;
t21 = -t124 * t33 - t176 * t60;
t20 = t124 * t32 + t177 * t60;
t19 = -t121 * t185 - t124 * t51 + t47;
t18 = t154 + t199;
t17 = -t191 + (pkin(5) * t130 - t132) * t125 + (-t178 - t116 + (-rSges(7,3) + t133) * t121) * t122 + t155;
t16 = t125 * t139 + (-t125 * rSges(5,3) + t122 * t152) * t122 + t188;
t13 = (-t122 * t33 + t125 * t32) * t121;
t10 = t122 * (rSges(6,3) * t177 - t156) + t125 * t164 + t158;
t9 = t176 * t27 + t29 * t74 + t31 * t75;
t8 = t176 * t26 + t28 * t74 + t30 * t75;
t7 = t177 * t27 + t29 * t72 + t31 * t73;
t6 = t177 * t26 + t28 * t72 + t30 * t73;
t5 = (-t183 + t199) * t125 + (-pkin(5) * t170 + t32 + (-t151 - t175 + t178) * t122) * t122 + t158;
t4 = t122 * t9 - t125 * t8;
t3 = t122 * t7 - t125 * t6;
t2 = -t124 * t15 + (t122 * t8 + t125 * t9) * t121;
t1 = -t124 * t14 + (t122 * t6 + t125 * t7) * t121;
t48 = [t136 * (Icges(4,2) * t136 + t182) + t134 * (Icges(4,1) * t134 + t181) + Icges(2,3) + Icges(3,3) + t47 + (-t51 + t180 - (Icges(6,5) * t131 - Icges(6,6) * t130) * t121 + (Icges(5,2) + Icges(6,3)) * t124) * t124 + (Icges(5,1) * t121 - t130 * t69 + t131 * t70 + t179 - t185) * t121 + m(7) * (t17 ^ 2 + t18 ^ 2) + m(6) * (t22 ^ 2 + t23 ^ 2) + m(4) * (t45 ^ 2 + t46 ^ 2) + m(5) * (t43 ^ 2 + t44 ^ 2) + m(3) * (t88 ^ 2 + t89 ^ 2) + m(2) * (t108 ^ 2 + t109 ^ 2); 0; m(3) + m(4) + m(5) - t195; m(7) * (t17 * t25 + t18 * t24) + m(6) * (t22 * t42 + t23 * t41) + m(5) * (t43 * t62 + t44 * t61) + (-t45 * t192 - t69 * t84 / 0.2e1 - t70 * t85 / 0.2e1 - t134 * (-Icges(4,5) * t125 + t122 * t145) / 0.2e1 - t136 * (-Icges(4,6) * t125 + t122 * t143) / 0.2e1 + t161 * t125 - t163) * t125 + (-t46 * t192 + t69 * t86 / 0.2e1 + t70 * t87 / 0.2e1 + (Icges(4,5) * t122 + t125 * t145) * t203 + (Icges(4,6) * t122 + t125 * t143) * t202 + t161 * t122 + t162) * t122 + ((Icges(5,6) * t204 - t122 * t142 / 0.2e1 + t35 / 0.2e1) * t125 + (Icges(5,6) * t194 + t142 * t204 - t36 / 0.2e1) * t122) * t124 + ((Icges(5,5) * t122 + t125 * t144 - t130 * t38 + t131 * t40) * t194 + (-Icges(5,5) * t125 + t122 * t144 - t130 * t37 + t131 * t39) * t193) * t121; m(4) * t34 + m(5) * t16 + m(6) * t10 + m(7) * t5; m(7) * (t24 ^ 2 + t25 ^ 2 + t5 ^ 2) + m(6) * (t10 ^ 2 + t41 ^ 2 + t42 ^ 2) + m(5) * (t16 ^ 2 + t61 ^ 2 + t62 ^ 2) + m(4) * (t107 ^ 2 * t166 + t34 ^ 2) + (-t3 + (t177 * t35 + t84 * t37 + t85 * t39) * t125 + t201 * t119) * t125 + (t4 + (t176 * t36 + t86 * t38 + t87 * t40) * t122 + t200 * t118 + (t201 * t122 + t200 * t125 - t176 * t35 - t177 * t36 - t37 * t86 - t38 * t84 - t39 * t87 - t40 * t85) * t125) * t122; m(7) * (t122 * t17 - t125 * t18) + m(6) * (t122 * t22 - t125 * t23) + m(5) * (t122 * t43 - t125 * t44); 0; m(7) * (t122 * t25 - t125 * t24) + m(6) * (t122 * t42 - t125 * t41) + m(5) * (t122 * t62 - t125 * t61); 0.2e1 * (m(5) / 0.2e1 + t165) * t166; 0.2e1 * ((t122 * t18 + t125 * t17) * t196 + (t122 * t23 + t125 * t22) * t197) * t121; t195 * t124; m(7) * (-t124 * t5 + (t122 * t24 + t125 * t25) * t121) + m(6) * (-t10 * t124 + (t122 * t41 + t125 * t42) * t121); 0; 0.2e1 * t165 * (t121 ^ 2 * t166 + t198); m(7) * (t17 * t20 + t18 * t21) - t19 * t124 + (t122 * t163 + t125 * t162) * t121; m(7) * t13; t2 * t194 - t124 * (-t11 * t125 + t12 * t122) / 0.2e1 + m(7) * (t13 * t5 + t20 * t25 + t21 * t24) + t1 * t193 + (t3 * t194 + t4 * t204) * t121; m(7) * (t122 * t20 - t125 * t21); m(7) * (-t124 * t13 + (t122 * t21 + t125 * t20) * t121); t198 * t19 + m(7) * (t13 ^ 2 + t20 ^ 2 + t21 ^ 2) + (t125 * t2 + t122 * t1 - t124 * (t11 * t122 + t12 * t125)) * t121;];
%% Postprocessing: Reshape Output
% From vec2symmat_6_matlab.m
res = [t48(1) t48(2) t48(4) t48(7) t48(11) t48(16); t48(2) t48(3) t48(5) t48(8) t48(12) t48(17); t48(4) t48(5) t48(6) t48(9) t48(13) t48(18); t48(7) t48(8) t48(9) t48(10) t48(14) t48(19); t48(11) t48(12) t48(13) t48(14) t48(15) t48(20); t48(16) t48(17) t48(18) t48(19) t48(20) t48(21);];
Mq  = res;
