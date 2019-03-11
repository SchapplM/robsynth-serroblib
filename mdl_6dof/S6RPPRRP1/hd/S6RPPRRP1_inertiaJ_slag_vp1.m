% Calculate joint inertia matrix for
% S6RPPRRP1
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% pkin [10x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d4,d5,theta2,theta3]';
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
% Datum: 2019-03-09 01:59
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Mq = S6RPPRRP1_inertiaJ_slag_vp1(qJ, ...
  pkin, m, rSges, Icges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(10,1),zeros(7,1),zeros(7,3),zeros(7,6)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPPRRP1_inertiaJ_slag_vp1: qJ has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S6RPPRRP1_inertiaJ_slag_vp1: pkin has to be [10x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6RPPRRP1_inertiaJ_slag_vp1: m has to be [7x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [7,3]), ...
  'S6RPPRRP1_inertiaJ_slag_vp1: rSges has to be [7x3] (double)');
assert(isreal(Icges) && all(size(Icges) == [7 6]), ...
  'S6RPPRRP1_inertiaJ_slag_vp1: Icges has to be [7x6] (double)'); 

%% Symbolic Calculation
% From inertia_joint_joint_floatb_twist_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 01:57:37
% EndTime: 2019-03-09 01:57:40
% DurationCPUTime: 1.72s
% Computational Cost: add. (5686->302), mult. (4741->433), div. (0->0), fcn. (5002->10), ass. (0->142)
t127 = -qJ(6) - pkin(8);
t188 = rSges(7,3) - t127;
t123 = pkin(10) + qJ(4);
t118 = sin(t123);
t187 = Icges(5,5) * t118;
t129 = sin(qJ(5));
t120 = cos(t123);
t131 = cos(qJ(5));
t79 = -Icges(7,6) * t120 + (Icges(7,4) * t131 - Icges(7,2) * t129) * t118;
t80 = -Icges(6,6) * t120 + (Icges(6,4) * t131 - Icges(6,2) * t129) * t118;
t186 = (t79 + t80) * t129;
t185 = t187 / 0.2e1;
t77 = -Icges(7,3) * t120 + (Icges(7,5) * t131 - Icges(7,6) * t129) * t118;
t78 = -Icges(6,3) * t120 + (Icges(6,5) * t131 - Icges(6,6) * t129) * t118;
t184 = t77 + t78;
t124 = qJ(1) + pkin(9);
t119 = sin(t124);
t158 = t118 * t119;
t121 = cos(t124);
t150 = t121 * t131;
t154 = t119 * t129;
t89 = -t120 * t154 - t150;
t151 = t121 * t129;
t153 = t119 * t131;
t90 = t120 * t153 - t151;
t44 = Icges(7,5) * t90 + Icges(7,6) * t89 + Icges(7,3) * t158;
t48 = Icges(7,4) * t90 + Icges(7,2) * t89 + Icges(7,6) * t158;
t52 = Icges(7,1) * t90 + Icges(7,4) * t89 + Icges(7,5) * t158;
t11 = t44 * t158 + t48 * t89 + t52 * t90;
t157 = t118 * t121;
t91 = -t120 * t151 + t153;
t92 = t120 * t150 + t154;
t45 = Icges(7,5) * t92 + Icges(7,6) * t91 + Icges(7,3) * t157;
t49 = Icges(7,4) * t92 + Icges(7,2) * t91 + Icges(7,6) * t157;
t53 = Icges(7,1) * t92 + Icges(7,4) * t91 + Icges(7,5) * t157;
t12 = t45 * t158 + t49 * t89 + t53 * t90;
t46 = Icges(6,5) * t90 + Icges(6,6) * t89 + Icges(6,3) * t158;
t50 = Icges(6,4) * t90 + Icges(6,2) * t89 + Icges(6,6) * t158;
t54 = Icges(6,1) * t90 + Icges(6,4) * t89 + Icges(6,5) * t158;
t13 = t46 * t158 + t50 * t89 + t54 * t90;
t47 = Icges(6,5) * t92 + Icges(6,6) * t91 + Icges(6,3) * t157;
t51 = Icges(6,4) * t92 + Icges(6,2) * t91 + Icges(6,6) * t157;
t55 = Icges(6,1) * t92 + Icges(6,4) * t91 + Icges(6,5) * t157;
t14 = t47 * t158 + t51 * t89 + t55 * t90;
t81 = -Icges(7,5) * t120 + (Icges(7,1) * t131 - Icges(7,4) * t129) * t118;
t27 = t77 * t158 + t79 * t89 + t81 * t90;
t82 = -Icges(6,5) * t120 + (Icges(6,1) * t131 - Icges(6,4) * t129) * t118;
t28 = t78 * t158 + t80 * t89 + t82 * t90;
t182 = (-t27 - t28) * t120 + ((t12 + t14) * t121 + (t11 + t13) * t119) * t118;
t15 = t44 * t157 + t48 * t91 + t52 * t92;
t16 = t45 * t157 + t49 * t91 + t53 * t92;
t17 = t46 * t157 + t50 * t91 + t54 * t92;
t18 = t47 * t157 + t51 * t91 + t55 * t92;
t29 = t77 * t157 + t79 * t91 + t81 * t92;
t30 = t78 * t157 + t80 * t91 + t82 * t92;
t181 = (-t29 - t30) * t120 + ((t16 + t18) * t121 + (t15 + t17) * t119) * t118;
t21 = -t120 * t44 + (-t129 * t48 + t131 * t52) * t118;
t23 = -t120 * t46 + (-t129 * t50 + t131 * t54) * t118;
t180 = -t21 - t23;
t22 = -t120 * t45 + (-t129 * t49 + t131 * t53) * t118;
t24 = -t120 * t47 + (-t129 * t51 + t131 * t55) * t118;
t179 = t22 + t24;
t178 = (t81 + t82) * t118 * t131;
t114 = pkin(5) * t131 + pkin(4);
t152 = t120 * t121;
t177 = t92 * rSges(7,1) + t91 * rSges(7,2) + pkin(5) * t154 + t114 * t152 + t188 * t157;
t116 = t119 ^ 2;
t176 = t120 ^ 2;
t117 = t121 ^ 2;
t175 = t119 / 0.2e1;
t174 = -t120 / 0.2e1;
t100 = rSges(5,1) * t118 + rSges(5,2) * t120;
t172 = m(5) * t100;
t130 = sin(qJ(1));
t171 = pkin(1) * t130;
t170 = pkin(4) * t120;
t169 = -pkin(4) + t114;
t168 = pkin(8) + t127;
t167 = t118 * t186 + t184 * t120 - t178;
t145 = -rSges(7,1) * t90 - rSges(7,2) * t89;
t166 = -pkin(5) * t151 + (-t168 * t118 + t169 * t120) * t119 + rSges(7,3) * t158 - t145;
t149 = pkin(4) * t152 + pkin(8) * t157;
t165 = -t149 + t177;
t164 = (t168 - rSges(7,3)) * t120 + (rSges(7,1) * t131 - rSges(7,2) * t129 + t169) * t118;
t163 = t116 * (pkin(8) * t118 + t170) + t121 * t149;
t162 = rSges(4,3) + qJ(3);
t101 = pkin(4) * t118 - pkin(8) * t120;
t84 = -rSges(6,3) * t120 + (rSges(6,1) * t131 - rSges(6,2) * t129) * t118;
t161 = -t101 - t84;
t159 = Icges(5,4) * t120;
t148 = t116 + t117;
t147 = -t101 - t164;
t59 = t92 * rSges(6,1) + t91 * rSges(6,2) + rSges(6,3) * t157;
t146 = -rSges(6,1) * t90 - rSges(6,2) * t89;
t126 = cos(pkin(10));
t113 = pkin(3) * t126 + pkin(2);
t132 = cos(qJ(1));
t122 = t132 * pkin(1);
t128 = -pkin(7) - qJ(3);
t144 = t121 * t113 - t119 * t128 + t122;
t143 = rSges(5,1) * t120 - rSges(5,2) * t118;
t139 = -Icges(5,2) * t118 + t159;
t138 = Icges(5,5) * t120 - Icges(5,6) * t118;
t137 = rSges(5,1) * t152 - rSges(5,2) * t157 + t119 * rSges(5,3);
t125 = sin(pkin(10));
t135 = rSges(4,1) * t126 - rSges(4,2) * t125 + pkin(2);
t134 = t21 / 0.2e1 + t28 / 0.2e1 + t27 / 0.2e1 + t23 / 0.2e1;
t133 = t24 / 0.2e1 + t22 / 0.2e1 + t30 / 0.2e1 + t29 / 0.2e1;
t110 = rSges(2,1) * t132 - t130 * rSges(2,2);
t109 = -t130 * rSges(2,1) - rSges(2,2) * t132;
t97 = Icges(5,6) * t120 + t187;
t94 = rSges(3,1) * t121 - rSges(3,2) * t119 + t122;
t93 = -rSges(3,1) * t119 - rSges(3,2) * t121 - t171;
t70 = Icges(5,3) * t119 + t138 * t121;
t69 = -Icges(5,3) * t121 + t138 * t119;
t65 = t162 * t119 + t135 * t121 + t122;
t64 = -t135 * t119 + t162 * t121 - t171;
t63 = t161 * t121;
t62 = t161 * t119;
t61 = t137 + t144;
t60 = -t171 + (rSges(5,3) - t128) * t121 + (-t113 - t143) * t119;
t57 = rSges(6,3) * t158 - t146;
t41 = t121 * t137 + (-t121 * rSges(5,3) + t143 * t119) * t119;
t40 = t147 * t121;
t39 = t147 * t119;
t36 = -t120 * t59 - t84 * t157;
t35 = t120 * t57 + t84 * t158;
t34 = t144 + t59 + t149;
t33 = -t171 - t121 * t128 + (-t170 - t113 + (-rSges(6,3) - pkin(8)) * t118) * t119 + t146;
t32 = t144 + t177;
t31 = -t171 + (pkin(5) * t129 - t128) * t121 + (-t114 * t120 - t118 * t188 - t113) * t119 + t145;
t26 = (-t119 * t59 + t121 * t57) * t118;
t25 = t119 * t57 + t121 * t59 + t163;
t20 = -t165 * t120 - t164 * t157;
t19 = t166 * t120 + t164 * t158;
t10 = (-t165 * t119 + t166 * t121) * t118;
t9 = t166 * t119 + t165 * t121 + t163;
t8 = t119 * t18 - t121 * t17;
t7 = t119 * t16 - t121 * t15;
t6 = t119 * t14 - t121 * t13;
t5 = -t11 * t121 + t119 * t12;
t1 = [Icges(4,2) * t126 ^ 2 + Icges(2,3) + Icges(3,3) + (Icges(4,1) * t125 + 0.2e1 * Icges(4,4) * t126) * t125 + (Icges(5,4) * t118 + Icges(5,2) * t120 - t184) * t120 + (Icges(5,1) * t118 + t159 - t186) * t118 + m(7) * (t31 ^ 2 + t32 ^ 2) + m(6) * (t33 ^ 2 + t34 ^ 2) + m(5) * (t60 ^ 2 + t61 ^ 2) + m(4) * (t64 ^ 2 + t65 ^ 2) + m(3) * (t93 ^ 2 + t94 ^ 2) + m(2) * (t109 ^ 2 + t110 ^ 2) + t178; 0; m(3) + m(4) + m(5) + m(6) + m(7); m(7) * (t119 * t31 - t121 * t32) + m(6) * (t119 * t33 - t121 * t34) + m(5) * (t119 * t60 - t121 * t61) + m(4) * (t119 * t64 - t121 * t65); 0; 0.2e1 * (m(4) / 0.2e1 + m(5) / 0.2e1 + m(6) / 0.2e1 + m(7) / 0.2e1) * t148; m(7) * (t31 * t40 + t32 * t39) + m(6) * (t33 * t63 + t34 * t62) + (t139 * t119 * t174 - t60 * t172 - t134 + (t185 - Icges(5,6) * t174 + t97 / 0.2e1) * t121) * t121 + (t119 * t185 + t120 * (Icges(5,6) * t119 + t139 * t121) / 0.2e1 - t61 * t172 + t97 * t175 + t133) * t119; m(5) * t41 + m(6) * t25 + m(7) * t9; m(6) * (t119 * t63 - t121 * t62) + m(7) * (t119 * t40 - t121 * t39); m(7) * (t39 ^ 2 + t40 ^ 2 + t9 ^ 2) + m(6) * (t25 ^ 2 + t62 ^ 2 + t63 ^ 2) + m(5) * (t148 * t100 ^ 2 + t41 ^ 2) + (t116 * t70 + t7 + t8) * t119 + (-t117 * t69 - t5 - t6 + (-t119 * t69 + t121 * t70) * t119) * t121; t167 * t120 + m(7) * (t19 * t31 + t20 * t32) + m(6) * (t33 * t35 + t34 * t36) + (t134 * t119 + t133 * t121) * t118; m(6) * t26 + m(7) * t10; m(6) * (t119 * t35 - t121 * t36) + m(7) * (t119 * t19 - t121 * t20); m(7) * (t10 * t9 + t19 * t40 + t20 * t39) + m(6) * (t25 * t26 + t35 * t63 + t36 * t62) + ((t8 / 0.2e1 + t7 / 0.2e1) * t121 + (t5 / 0.2e1 + t6 / 0.2e1) * t119) * t118 + t181 * t175 + (t179 * t119 + t180 * t121) * t174 - t182 * t121 / 0.2e1; m(7) * (t10 ^ 2 + t19 ^ 2 + t20 ^ 2) + m(6) * (t26 ^ 2 + t35 ^ 2 + t36 ^ 2) - t167 * t176 + ((-t179 * t120 + t181) * t121 + (t180 * t120 + t182) * t119) * t118; m(7) * (t119 * t32 + t121 * t31) * t118; -m(7) * t120; 0; m(7) * (-t120 * t9 + (t119 * t39 + t121 * t40) * t118); m(7) * (-t10 * t120 + (t119 * t20 + t121 * t19) * t118); m(7) * (t148 * t118 ^ 2 + t176);];
%% Postprocessing: Reshape Output
% From vec2symmat_6_matlab.m
res = [t1(1) t1(2) t1(4) t1(7) t1(11) t1(16); t1(2) t1(3) t1(5) t1(8) t1(12) t1(17); t1(4) t1(5) t1(6) t1(9) t1(13) t1(18); t1(7) t1(8) t1(9) t1(10) t1(14) t1(19); t1(11) t1(12) t1(13) t1(14) t1(15) t1(20); t1(16) t1(17) t1(18) t1(19) t1(20) t1(21);];
Mq  = res;
