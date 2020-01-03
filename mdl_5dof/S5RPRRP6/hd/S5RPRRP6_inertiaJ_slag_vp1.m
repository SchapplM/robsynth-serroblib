% Calculate joint inertia matrix for
% S5RPRRP6
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% pkin [8x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d1,d3,d4,theta2]';
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
% Datum: 2019-12-31 18:43
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Mq = S5RPRRP6_inertiaJ_slag_vp1(qJ, ...
  pkin, m, rSges, Icges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(8,1),zeros(6,1),zeros(6,3),zeros(6,6)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RPRRP6_inertiaJ_slag_vp1: qJ has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S5RPRRP6_inertiaJ_slag_vp1: pkin has to be [8x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5RPRRP6_inertiaJ_slag_vp1: m has to be [6x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [6,3]), ...
  'S5RPRRP6_inertiaJ_slag_vp1: rSges has to be [6x3] (double)');
assert(isreal(Icges) && all(size(Icges) == [6 6]), ...
  'S5RPRRP6_inertiaJ_slag_vp1: Icges has to be [6x6] (double)'); 

%% Symbolic Calculation
% From inertia_joint_joint_floatb_twist_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 18:42:09
% EndTime: 2019-12-31 18:42:13
% DurationCPUTime: 1.43s
% Computational Cost: add. (3659->267), mult. (4353->386), div. (0->0), fcn. (4652->8), ass. (0->135)
t120 = -qJ(5) - pkin(7);
t184 = rSges(6,3) - t120;
t122 = sin(qJ(3));
t183 = Icges(4,5) * t122;
t121 = sin(qJ(4));
t124 = cos(qJ(4));
t125 = cos(qJ(3));
t85 = -Icges(6,6) * t125 + (Icges(6,4) * t124 - Icges(6,2) * t121) * t122;
t86 = -Icges(5,6) * t125 + (Icges(5,4) * t124 - Icges(5,2) * t121) * t122;
t182 = (t85 + t86) * t121;
t181 = t183 / 0.2e1;
t83 = -Icges(6,3) * t125 + (Icges(6,5) * t124 - Icges(6,6) * t121) * t122;
t84 = -Icges(5,3) * t125 + (Icges(5,5) * t124 - Icges(5,6) * t121) * t122;
t180 = t83 + t84;
t118 = qJ(1) + pkin(8);
t115 = sin(t118);
t152 = t115 * t122;
t116 = cos(t118);
t147 = t121 * t125;
t79 = -t115 * t147 - t116 * t124;
t145 = t124 * t125;
t151 = t116 * t121;
t80 = t115 * t145 - t151;
t42 = Icges(6,5) * t80 + Icges(6,6) * t79 + Icges(6,3) * t152;
t46 = Icges(6,4) * t80 + Icges(6,2) * t79 + Icges(6,6) * t152;
t50 = Icges(6,1) * t80 + Icges(6,4) * t79 + Icges(6,5) * t152;
t11 = t42 * t152 + t46 * t79 + t50 * t80;
t150 = t116 * t122;
t81 = t115 * t124 - t116 * t147;
t153 = t115 * t121;
t82 = t116 * t145 + t153;
t43 = Icges(6,5) * t82 + Icges(6,6) * t81 + Icges(6,3) * t150;
t47 = Icges(6,4) * t82 + Icges(6,2) * t81 + Icges(6,6) * t150;
t51 = Icges(6,1) * t82 + Icges(6,4) * t81 + Icges(6,5) * t150;
t12 = t43 * t152 + t47 * t79 + t51 * t80;
t44 = Icges(5,5) * t80 + Icges(5,6) * t79 + Icges(5,3) * t152;
t48 = Icges(5,4) * t80 + Icges(5,2) * t79 + Icges(5,6) * t152;
t52 = Icges(5,1) * t80 + Icges(5,4) * t79 + Icges(5,5) * t152;
t13 = t44 * t152 + t48 * t79 + t52 * t80;
t45 = Icges(5,5) * t82 + Icges(5,6) * t81 + Icges(5,3) * t150;
t49 = Icges(5,4) * t82 + Icges(5,2) * t81 + Icges(5,6) * t150;
t53 = Icges(5,1) * t82 + Icges(5,4) * t81 + Icges(5,5) * t150;
t14 = t45 * t152 + t49 * t79 + t53 * t80;
t87 = -Icges(6,5) * t125 + (Icges(6,1) * t124 - Icges(6,4) * t121) * t122;
t27 = t83 * t152 + t79 * t85 + t80 * t87;
t88 = -Icges(5,5) * t125 + (Icges(5,1) * t124 - Icges(5,4) * t121) * t122;
t28 = t84 * t152 + t79 * t86 + t80 * t88;
t178 = (-t27 - t28) * t125 + ((t12 + t14) * t116 + (t11 + t13) * t115) * t122;
t15 = t42 * t150 + t46 * t81 + t50 * t82;
t16 = t43 * t150 + t47 * t81 + t51 * t82;
t17 = t44 * t150 + t48 * t81 + t52 * t82;
t18 = t45 * t150 + t49 * t81 + t53 * t82;
t29 = t83 * t150 + t81 * t85 + t82 * t87;
t30 = t84 * t150 + t81 * t86 + t82 * t88;
t177 = (-t29 - t30) * t125 + ((t16 + t18) * t116 + (t15 + t17) * t115) * t122;
t19 = -t125 * t42 + (-t121 * t46 + t124 * t50) * t122;
t21 = -t125 * t44 + (-t121 * t48 + t124 * t52) * t122;
t176 = -t19 - t21;
t20 = -t125 * t43 + (-t121 * t47 + t124 * t51) * t122;
t22 = -t125 * t45 + (-t121 * t49 + t124 * t53) * t122;
t175 = t20 + t22;
t174 = (t87 + t88) * t122 * t124;
t112 = pkin(4) * t124 + pkin(3);
t149 = t116 * t125;
t173 = t82 * rSges(6,1) + t81 * rSges(6,2) + pkin(4) * t153 + t112 * t149 + t184 * t150;
t172 = rSges(6,1) * t80 + rSges(6,2) * t79 - pkin(4) * t151;
t113 = t115 ^ 2;
t114 = t116 ^ 2;
t171 = t125 ^ 2;
t98 = rSges(4,1) * t122 + rSges(4,2) * t125;
t170 = m(4) * t98;
t169 = t115 / 0.2e1;
t167 = -t125 / 0.2e1;
t123 = sin(qJ(1));
t166 = pkin(1) * t123;
t165 = pkin(3) * t125;
t164 = -pkin(3) + t112;
t163 = pkin(7) + t120;
t162 = t122 * t182 + t180 * t125 - t174;
t161 = (-t163 * t122 + t164 * t125) * t115 + rSges(6,3) * t152 + t172;
t144 = pkin(3) * t149 + pkin(7) * t150;
t160 = -t144 + t173;
t159 = t113 * (pkin(7) * t122 + t165) + t116 * t144;
t158 = (t163 - rSges(6,3)) * t125 + (rSges(6,1) * t124 - rSges(6,2) * t121 + t164) * t122;
t157 = t116 * rSges(4,3);
t104 = pkin(3) * t122 - pkin(7) * t125;
t90 = -t125 * rSges(5,3) + (rSges(5,1) * t124 - rSges(5,2) * t121) * t122;
t156 = -t104 - t90;
t154 = Icges(4,4) * t125;
t143 = t113 + t114;
t141 = -t104 - t158;
t59 = t82 * rSges(5,1) + t81 * rSges(5,2) + rSges(5,3) * t150;
t126 = cos(qJ(1));
t117 = t126 * pkin(1);
t140 = t116 * pkin(2) + t115 * pkin(6) + t117;
t139 = t116 * pkin(6) - t166;
t138 = -rSges(5,1) * t80 - rSges(5,2) * t79;
t136 = rSges(4,1) * t125 - rSges(4,2) * t122;
t132 = -Icges(4,2) * t122 + t154;
t131 = Icges(4,5) * t125 - Icges(4,6) * t122;
t130 = rSges(4,1) * t149 - rSges(4,2) * t150 + t115 * rSges(4,3);
t128 = t28 / 0.2e1 + t27 / 0.2e1 + t21 / 0.2e1 + t19 / 0.2e1;
t127 = t30 / 0.2e1 + t29 / 0.2e1 + t22 / 0.2e1 + t20 / 0.2e1;
t100 = rSges(2,1) * t126 - t123 * rSges(2,2);
t99 = -t123 * rSges(2,1) - rSges(2,2) * t126;
t95 = Icges(4,6) * t125 + t183;
t92 = rSges(3,1) * t116 - rSges(3,2) * t115 + t117;
t91 = -rSges(3,1) * t115 - rSges(3,2) * t116 - t166;
t65 = Icges(4,3) * t115 + t131 * t116;
t64 = -Icges(4,3) * t116 + t131 * t115;
t63 = t156 * t116;
t62 = t156 * t115;
t61 = t130 + t140;
t60 = t157 + (-pkin(2) - t136) * t115 + t139;
t57 = rSges(5,3) * t152 - t138;
t41 = t116 * t130 + (t136 * t115 - t157) * t115;
t40 = t141 * t116;
t39 = t141 * t115;
t36 = -t125 * t59 - t90 * t150;
t35 = t125 * t57 + t90 * t152;
t34 = t140 + t59 + t144;
t33 = (-t165 - pkin(2) + (-rSges(5,3) - pkin(7)) * t122) * t115 + t138 + t139;
t32 = t140 + t173;
t31 = (-t112 * t125 - t184 * t122 - pkin(2)) * t115 + t139 - t172;
t26 = (-t115 * t59 + t116 * t57) * t122;
t25 = t115 * t57 + t116 * t59 + t159;
t24 = -t160 * t125 - t158 * t150;
t23 = t161 * t125 + t158 * t152;
t10 = (-t160 * t115 + t161 * t116) * t122;
t9 = t161 * t115 + t160 * t116 + t159;
t8 = t115 * t18 - t116 * t17;
t7 = t115 * t16 - t116 * t15;
t6 = t115 * t14 - t116 * t13;
t5 = -t11 * t116 + t115 * t12;
t1 = [Icges(2,3) + Icges(3,3) + (Icges(4,4) * t122 + Icges(4,2) * t125 - t180) * t125 + (Icges(4,1) * t122 + t154 - t182) * t122 + m(5) * (t33 ^ 2 + t34 ^ 2) + m(6) * (t31 ^ 2 + t32 ^ 2) + m(4) * (t60 ^ 2 + t61 ^ 2) + m(3) * (t91 ^ 2 + t92 ^ 2) + m(2) * (t100 ^ 2 + t99 ^ 2) + t174; 0; m(3) + m(4) + m(5) + m(6); m(5) * (t33 * t63 + t34 * t62) + m(6) * (t31 * t40 + t32 * t39) + (t132 * t115 * t167 - t60 * t170 - t128 + (t181 - Icges(4,6) * t167 + t95 / 0.2e1) * t116) * t116 + (t115 * t181 + t125 * (Icges(4,6) * t115 + t132 * t116) / 0.2e1 - t61 * t170 + t95 * t169 + t127) * t115; m(4) * t41 + m(5) * t25 + m(6) * t9; m(6) * (t39 ^ 2 + t40 ^ 2 + t9 ^ 2) + m(5) * (t25 ^ 2 + t62 ^ 2 + t63 ^ 2) + m(4) * (t143 * t98 ^ 2 + t41 ^ 2) + (t113 * t65 + t7 + t8) * t115 + (-t114 * t64 - t5 - t6 + (-t115 * t64 + t116 * t65) * t115) * t116; t162 * t125 + m(5) * (t33 * t35 + t34 * t36) + m(6) * (t23 * t31 + t24 * t32) + (t128 * t115 + t127 * t116) * t122; m(5) * t26 + m(6) * t10; m(6) * (t10 * t9 + t23 * t40 + t24 * t39) + m(5) * (t25 * t26 + t35 * t63 + t36 * t62) + ((t7 / 0.2e1 + t8 / 0.2e1) * t116 + (t6 / 0.2e1 + t5 / 0.2e1) * t115) * t122 + t177 * t169 - t178 * t116 / 0.2e1 + (t175 * t115 + t176 * t116) * t167; m(6) * (t10 ^ 2 + t23 ^ 2 + t24 ^ 2) + m(5) * (t26 ^ 2 + t35 ^ 2 + t36 ^ 2) - t162 * t171 + ((-t175 * t125 + t177) * t116 + (t176 * t125 + t178) * t115) * t122; m(6) * (t115 * t32 + t116 * t31) * t122; -m(6) * t125; m(6) * (-t125 * t9 + (t115 * t39 + t116 * t40) * t122); m(6) * (-t125 * t10 + (t115 * t24 + t116 * t23) * t122); m(6) * (t143 * t122 ^ 2 + t171);];
%% Postprocessing: Reshape Output
% From vec2symmat_5_matlab.m
res = [t1(1), t1(2), t1(4), t1(7), t1(11); t1(2), t1(3), t1(5), t1(8), t1(12); t1(4), t1(5), t1(6), t1(9), t1(13); t1(7), t1(8), t1(9), t1(10), t1(14); t1(11), t1(12), t1(13), t1(14), t1(15);];
Mq = res;
