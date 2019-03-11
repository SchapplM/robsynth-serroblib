% Calculate joint inertia matrix for
% S6RPPRRP2
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
% Datum: 2019-03-09 02:01
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Mq = S6RPPRRP2_inertiaJ_slag_vp1(qJ, ...
  pkin, m, rSges, Icges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(10,1),zeros(7,1),zeros(7,3),zeros(7,6)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPPRRP2_inertiaJ_slag_vp1: qJ has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S6RPPRRP2_inertiaJ_slag_vp1: pkin has to be [10x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6RPPRRP2_inertiaJ_slag_vp1: m has to be [7x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [7,3]), ...
  'S6RPPRRP2_inertiaJ_slag_vp1: rSges has to be [7x3] (double)');
assert(isreal(Icges) && all(size(Icges) == [7 6]), ...
  'S6RPPRRP2_inertiaJ_slag_vp1: Icges has to be [7x6] (double)'); 

%% Symbolic Calculation
% From inertia_joint_joint_floatb_twist_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 02:00:12
% EndTime: 2019-03-09 02:00:16
% DurationCPUTime: 1.61s
% Computational Cost: add. (5496->291), mult. (4724->429), div. (0->0), fcn. (5066->10), ass. (0->139)
t122 = pkin(10) + qJ(4);
t117 = sin(t122);
t183 = Icges(5,5) * t117;
t182 = t183 / 0.2e1;
t178 = rSges(7,3) + qJ(6);
t180 = rSges(7,1) + pkin(5);
t127 = sin(qJ(5));
t123 = qJ(1) + pkin(9);
t120 = cos(t123);
t129 = cos(qJ(5));
t149 = t120 * t129;
t118 = sin(t123);
t119 = cos(t122);
t152 = t119 * t118;
t90 = t127 * t152 + t149;
t150 = t120 * t127;
t91 = t129 * t152 - t150;
t181 = -t178 * t90 - t180 * t91;
t78 = -Icges(6,3) * t119 + (Icges(6,5) * t129 - Icges(6,6) * t127) * t117;
t79 = -Icges(7,2) * t119 + (Icges(7,4) * t129 + Icges(7,6) * t127) * t117;
t179 = -t78 - t79;
t156 = t117 * t118;
t42 = Icges(7,5) * t91 + Icges(7,6) * t156 + Icges(7,3) * t90;
t46 = Icges(7,4) * t91 + Icges(7,2) * t156 + Icges(7,6) * t90;
t50 = Icges(7,1) * t91 + Icges(7,4) * t156 + Icges(7,5) * t90;
t11 = t46 * t156 + t42 * t90 + t50 * t91;
t155 = t117 * t120;
t92 = -t118 * t129 + t119 * t150;
t93 = t118 * t127 + t119 * t149;
t43 = Icges(7,5) * t93 + Icges(7,6) * t155 + Icges(7,3) * t92;
t47 = Icges(7,4) * t93 + Icges(7,2) * t155 + Icges(7,6) * t92;
t51 = Icges(7,1) * t93 + Icges(7,4) * t155 + Icges(7,5) * t92;
t12 = t47 * t156 + t43 * t90 + t51 * t91;
t44 = Icges(6,5) * t91 - Icges(6,6) * t90 + Icges(6,3) * t156;
t48 = Icges(6,4) * t91 - Icges(6,2) * t90 + Icges(6,6) * t156;
t52 = Icges(6,1) * t91 - Icges(6,4) * t90 + Icges(6,5) * t156;
t13 = t44 * t156 - t48 * t90 + t52 * t91;
t45 = Icges(6,5) * t93 - Icges(6,6) * t92 + Icges(6,3) * t155;
t49 = Icges(6,4) * t93 - Icges(6,2) * t92 + Icges(6,6) * t155;
t53 = Icges(6,1) * t93 - Icges(6,4) * t92 + Icges(6,5) * t155;
t14 = t45 * t156 - t49 * t90 + t53 * t91;
t77 = -Icges(7,6) * t119 + (Icges(7,5) * t129 + Icges(7,3) * t127) * t117;
t81 = -Icges(7,4) * t119 + (Icges(7,1) * t129 + Icges(7,5) * t127) * t117;
t29 = t79 * t156 + t77 * t90 + t81 * t91;
t80 = -Icges(6,6) * t119 + (Icges(6,4) * t129 - Icges(6,2) * t127) * t117;
t82 = -Icges(6,5) * t119 + (Icges(6,1) * t129 - Icges(6,4) * t127) * t117;
t30 = t78 * t156 - t80 * t90 + t82 * t91;
t177 = (-t29 - t30) * t119 + ((t12 + t14) * t120 + (t11 + t13) * t118) * t117;
t15 = t46 * t155 + t42 * t92 + t50 * t93;
t16 = t47 * t155 + t43 * t92 + t51 * t93;
t17 = t44 * t155 - t48 * t92 + t52 * t93;
t18 = t45 * t155 - t49 * t92 + t53 * t93;
t31 = t79 * t155 + t77 * t92 + t81 * t93;
t32 = t78 * t155 - t80 * t92 + t82 * t93;
t176 = (-t31 - t32) * t119 + ((t16 + t18) * t120 + (t15 + t17) * t118) * t117;
t19 = -t119 * t46 + (t127 * t42 + t129 * t50) * t117;
t21 = -t119 * t44 + (-t127 * t48 + t129 * t52) * t117;
t175 = -t19 - t21;
t20 = -t119 * t47 + (t127 * t43 + t129 * t51) * t117;
t22 = -t119 * t45 + (-t127 * t49 + t129 * t53) * t117;
t174 = t20 + t22;
t154 = t117 * t127;
t173 = t77 * t154 + (t81 + t82) * t117 * t129;
t115 = t118 ^ 2;
t116 = t120 ^ 2;
t172 = t118 / 0.2e1;
t171 = -t119 / 0.2e1;
t101 = rSges(5,1) * t117 + rSges(5,2) * t119;
t169 = m(5) * t101;
t128 = sin(qJ(1));
t168 = pkin(1) * t128;
t167 = pkin(4) * t119;
t166 = t119 * t179 - t80 * t154 + t173;
t165 = rSges(7,2) * t156 - t181;
t164 = rSges(7,2) * t155 + t178 * t92 + t180 * t93;
t151 = t119 * t120;
t148 = pkin(4) * t151 + pkin(8) * t155;
t162 = t115 * (pkin(8) * t117 + t167) + t120 * t148;
t161 = -rSges(7,2) * t119 + (t178 * t127 + t129 * t180) * t117;
t160 = rSges(4,3) + qJ(3);
t102 = pkin(4) * t117 - pkin(8) * t119;
t84 = -rSges(6,3) * t119 + (rSges(6,1) * t129 - rSges(6,2) * t127) * t117;
t159 = -t102 - t84;
t157 = Icges(5,4) * t119;
t147 = t115 + t116;
t146 = -t102 - t161;
t57 = t93 * rSges(6,1) - t92 * rSges(6,2) + rSges(6,3) * t155;
t125 = cos(pkin(10));
t114 = pkin(3) * t125 + pkin(2);
t145 = -t114 - t167;
t144 = -rSges(6,1) * t91 + rSges(6,2) * t90;
t130 = cos(qJ(1));
t121 = t130 * pkin(1);
t126 = -pkin(7) - qJ(3);
t143 = t120 * t114 - t118 * t126 + t121;
t142 = rSges(5,1) * t119 - rSges(5,2) * t117;
t141 = -t120 * t126 - t168;
t137 = -Icges(5,2) * t117 + t157;
t136 = Icges(5,5) * t119 - Icges(5,6) * t117;
t135 = rSges(5,1) * t151 - rSges(5,2) * t155 + t118 * rSges(5,3);
t124 = sin(pkin(10));
t134 = rSges(4,1) * t125 - rSges(4,2) * t124 + pkin(2);
t133 = t21 / 0.2e1 + t19 / 0.2e1 + t30 / 0.2e1 + t29 / 0.2e1;
t132 = t32 / 0.2e1 + t31 / 0.2e1 + t22 / 0.2e1 + t20 / 0.2e1;
t131 = t143 + t148;
t111 = rSges(2,1) * t130 - t128 * rSges(2,2);
t110 = -t128 * rSges(2,1) - rSges(2,2) * t130;
t98 = Icges(5,6) * t119 + t183;
t96 = rSges(3,1) * t120 - rSges(3,2) * t118 + t121;
t95 = -rSges(3,1) * t118 - rSges(3,2) * t120 - t168;
t70 = Icges(5,3) * t118 + t136 * t120;
t69 = -Icges(5,3) * t120 + t136 * t118;
t65 = t160 * t118 + t134 * t120 + t121;
t64 = -t134 * t118 + t160 * t120 - t168;
t61 = t159 * t120;
t60 = t159 * t118;
t59 = t135 + t143;
t58 = -t168 + (rSges(5,3) - t126) * t120 + (-t114 - t142) * t118;
t55 = rSges(6,3) * t156 - t144;
t41 = t120 * t135 + (-t120 * rSges(5,3) + t142 * t118) * t118;
t40 = t146 * t120;
t39 = t146 * t118;
t36 = -t119 * t57 - t84 * t155;
t35 = t119 * t55 + t84 * t156;
t34 = t131 + t57;
t33 = ((-rSges(6,3) - pkin(8)) * t117 + t145) * t118 + t141 + t144;
t28 = (-t118 * t57 + t120 * t55) * t117;
t27 = t131 + t164;
t26 = ((-rSges(7,2) - pkin(8)) * t117 + t145) * t118 + t141 + t181;
t25 = -t164 * t119 - t161 * t155;
t24 = t165 * t119 + t161 * t156;
t23 = t118 * t55 + t120 * t57 + t162;
t10 = (-t164 * t118 + t165 * t120) * t117;
t9 = t165 * t118 + t164 * t120 + t162;
t8 = t118 * t18 - t120 * t17;
t7 = t118 * t16 - t120 * t15;
t6 = t118 * t14 - t120 * t13;
t5 = -t11 * t120 + t118 * t12;
t1 = [Icges(4,2) * t125 ^ 2 + Icges(2,3) + Icges(3,3) + (Icges(4,1) * t124 + 0.2e1 * Icges(4,4) * t125) * t124 + (Icges(5,1) * t117 - t127 * t80 + t157) * t117 + (Icges(5,4) * t117 + Icges(5,2) * t119 + t179) * t119 + m(7) * (t26 ^ 2 + t27 ^ 2) + m(6) * (t33 ^ 2 + t34 ^ 2) + m(5) * (t58 ^ 2 + t59 ^ 2) + m(4) * (t64 ^ 2 + t65 ^ 2) + m(3) * (t95 ^ 2 + t96 ^ 2) + m(2) * (t110 ^ 2 + t111 ^ 2) + t173; 0; m(3) + m(4) + m(5) + m(6) + m(7); m(7) * (t118 * t26 - t120 * t27) + m(6) * (t118 * t33 - t120 * t34) + m(5) * (t118 * t58 - t120 * t59) + m(4) * (t118 * t64 - t120 * t65); 0; 0.2e1 * (m(4) / 0.2e1 + m(5) / 0.2e1 + m(6) / 0.2e1 + m(7) / 0.2e1) * t147; m(7) * (t26 * t40 + t27 * t39) + m(6) * (t33 * t61 + t34 * t60) + (t137 * t118 * t171 - t58 * t169 - t133 + (t182 - Icges(5,6) * t171 + t98 / 0.2e1) * t120) * t120 + (t118 * t182 + t119 * (Icges(5,6) * t118 + t137 * t120) / 0.2e1 - t59 * t169 + t98 * t172 + t132) * t118; m(5) * t41 + m(6) * t23 + m(7) * t9; m(6) * (t118 * t61 - t120 * t60) + m(7) * (t118 * t40 - t120 * t39); m(7) * (t39 ^ 2 + t40 ^ 2 + t9 ^ 2) + m(6) * (t23 ^ 2 + t60 ^ 2 + t61 ^ 2) + m(5) * (t147 * t101 ^ 2 + t41 ^ 2) + (t115 * t70 + t7 + t8) * t118 + (-t116 * t69 - t5 - t6 + (-t118 * t69 + t120 * t70) * t118) * t120; -t166 * t119 + m(7) * (t24 * t26 + t25 * t27) + m(6) * (t33 * t35 + t34 * t36) + (t133 * t118 + t132 * t120) * t117; m(6) * t28 + m(7) * t10; m(6) * (t118 * t35 - t120 * t36) + m(7) * (t118 * t24 - t120 * t25); m(7) * (t10 * t9 + t24 * t40 + t25 * t39) + m(6) * (t23 * t28 + t35 * t61 + t36 * t60) + ((t7 / 0.2e1 + t8 / 0.2e1) * t120 + (t5 / 0.2e1 + t6 / 0.2e1) * t118) * t117 + t176 * t172 + (t174 * t118 + t175 * t120) * t171 - t177 * t120 / 0.2e1; m(7) * (t10 ^ 2 + t24 ^ 2 + t25 ^ 2) + m(6) * (t28 ^ 2 + t35 ^ 2 + t36 ^ 2) + t166 * t119 ^ 2 + ((-t174 * t119 + t176) * t120 + (t175 * t119 + t177) * t118) * t117; m(7) * (t26 * t92 + t27 * t90); m(7) * t154; m(7) * (t118 * t92 - t120 * t90); m(7) * (t9 * t154 + t39 * t90 + t40 * t92); m(7) * (t10 * t154 + t24 * t92 + t25 * t90); m(7) * (t117 ^ 2 * t127 ^ 2 + t90 ^ 2 + t92 ^ 2);];
%% Postprocessing: Reshape Output
% From vec2symmat_6_matlab.m
res = [t1(1) t1(2) t1(4) t1(7) t1(11) t1(16); t1(2) t1(3) t1(5) t1(8) t1(12) t1(17); t1(4) t1(5) t1(6) t1(9) t1(13) t1(18); t1(7) t1(8) t1(9) t1(10) t1(14) t1(19); t1(11) t1(12) t1(13) t1(14) t1(15) t1(20); t1(16) t1(17) t1(18) t1(19) t1(20) t1(21);];
Mq  = res;
