% Calculate joint inertia matrix for
% S6RPPRRP3
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% pkin [9x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d4,d5,theta2]';
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
% Datum: 2019-03-09 02:04
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Mq = S6RPPRRP3_inertiaJ_slag_vp1(qJ, ...
  pkin, m, rSges, Icges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(9,1),zeros(7,1),zeros(7,3),zeros(7,6)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPPRRP3_inertiaJ_slag_vp1: qJ has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S6RPPRRP3_inertiaJ_slag_vp1: pkin has to be [9x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6RPPRRP3_inertiaJ_slag_vp1: m has to be [7x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [7,3]), ...
  'S6RPPRRP3_inertiaJ_slag_vp1: rSges has to be [7x3] (double)');
assert(isreal(Icges) && all(size(Icges) == [7 6]), ...
  'S6RPPRRP3_inertiaJ_slag_vp1: Icges has to be [7x6] (double)'); 

%% Symbolic Calculation
% From inertia_joint_joint_floatb_twist_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 02:02:43
% EndTime: 2019-03-09 02:02:46
% DurationCPUTime: 1.63s
% Computational Cost: add. (3901->284), mult. (4718->418), div. (0->0), fcn. (5054->8), ass. (0->136)
t130 = cos(qJ(4));
t188 = Icges(5,5) * t130;
t127 = sin(qJ(4));
t187 = Icges(5,6) * t127;
t184 = rSges(7,1) + pkin(5);
t183 = rSges(7,3) + qJ(6);
t125 = qJ(1) + pkin(9);
t122 = sin(t125);
t123 = cos(t125);
t129 = cos(qJ(5));
t126 = sin(qJ(5));
t155 = t126 * t127;
t85 = t122 * t129 + t123 * t155;
t153 = t127 * t129;
t87 = t122 * t126 - t123 * t153;
t186 = -t183 * t85 + t184 * t87;
t185 = t188 / 0.2e1 - t187 / 0.2e1;
t154 = t126 * t130;
t88 = Icges(7,6) * t127 + (Icges(7,5) * t129 + Icges(7,3) * t126) * t130;
t89 = Icges(6,3) * t127 + (Icges(6,5) * t129 - Icges(6,6) * t126) * t130;
t90 = Icges(7,2) * t127 + (Icges(7,4) * t129 + Icges(7,6) * t126) * t130;
t92 = Icges(7,4) * t127 + (Icges(7,1) * t129 + Icges(7,5) * t126) * t130;
t93 = Icges(6,5) * t127 + (Icges(6,1) * t129 - Icges(6,4) * t126) * t130;
t182 = t88 * t154 + (t92 + t93) * t129 * t130 + (t89 + t90) * t127;
t157 = t122 * t130;
t83 = t122 * t155 - t123 * t129;
t84 = t122 * t153 + t123 * t126;
t42 = Icges(7,5) * t84 - Icges(7,6) * t157 + Icges(7,3) * t83;
t46 = Icges(7,4) * t84 - Icges(7,2) * t157 + Icges(7,6) * t83;
t50 = Icges(7,1) * t84 - Icges(7,4) * t157 + Icges(7,5) * t83;
t11 = -t46 * t157 + t42 * t83 + t50 * t84;
t156 = t123 * t130;
t43 = Icges(7,5) * t87 + Icges(7,6) * t156 - Icges(7,3) * t85;
t47 = Icges(7,4) * t87 + Icges(7,2) * t156 - Icges(7,6) * t85;
t51 = Icges(7,1) * t87 + Icges(7,4) * t156 - Icges(7,5) * t85;
t12 = -t47 * t157 + t43 * t83 + t51 * t84;
t44 = Icges(6,5) * t84 - Icges(6,6) * t83 - Icges(6,3) * t157;
t48 = Icges(6,4) * t84 - Icges(6,2) * t83 - Icges(6,6) * t157;
t52 = Icges(6,1) * t84 - Icges(6,4) * t83 - Icges(6,5) * t157;
t13 = -t44 * t157 - t48 * t83 + t52 * t84;
t45 = Icges(6,5) * t87 + Icges(6,6) * t85 + Icges(6,3) * t156;
t49 = Icges(6,4) * t87 + Icges(6,2) * t85 + Icges(6,6) * t156;
t53 = Icges(6,1) * t87 + Icges(6,4) * t85 + Icges(6,5) * t156;
t14 = -t45 * t157 - t49 * t83 + t53 * t84;
t29 = -t90 * t157 + t83 * t88 + t84 * t92;
t91 = Icges(6,6) * t127 + (Icges(6,4) * t129 - Icges(6,2) * t126) * t130;
t30 = -t89 * t157 - t83 * t91 + t84 * t93;
t181 = ((t12 + t14) * t123 + (-t11 - t13) * t122) * t130 + (t29 + t30) * t127;
t15 = t46 * t156 - t42 * t85 + t50 * t87;
t16 = t47 * t156 - t43 * t85 + t51 * t87;
t17 = t44 * t156 + t48 * t85 + t52 * t87;
t18 = t45 * t156 + t49 * t85 + t53 * t87;
t31 = t90 * t156 - t85 * t88 + t87 * t92;
t32 = t89 * t156 + t85 * t91 + t87 * t93;
t180 = ((t16 + t18) * t123 + (-t15 - t17) * t122) * t130 + (t31 + t32) * t127;
t19 = t127 * t46 + (t126 * t42 + t129 * t50) * t130;
t21 = t127 * t44 + (-t126 * t48 + t129 * t52) * t130;
t179 = t19 + t21;
t20 = t127 * t47 + (t126 * t43 + t129 * t51) * t130;
t22 = t127 * t45 + (-t126 * t49 + t129 * t53) * t130;
t178 = t20 + t22;
t177 = t183 * t83 + t184 * t84;
t176 = (rSges(5,1) * t127 + rSges(5,2) * t130) * t123;
t120 = t122 ^ 2;
t121 = t123 ^ 2;
t175 = -pkin(2) - pkin(7);
t174 = t122 / 0.2e1;
t173 = t123 / 0.2e1;
t107 = rSges(5,1) * t130 - rSges(5,2) * t127;
t169 = m(5) * t107;
t128 = sin(qJ(1));
t168 = pkin(1) * t128;
t167 = (-t91 * t154 + t182) * t127;
t166 = -rSges(7,2) * t157 + t177;
t165 = rSges(7,2) * t156 + t186;
t162 = t84 * rSges(6,1) - t83 * rSges(6,2);
t161 = rSges(7,2) * t127 + (t183 * t126 + t129 * t184) * t130;
t158 = t122 * t127;
t151 = t120 + t121;
t149 = rSges(5,1) * t158 + rSges(5,2) * t157 + t123 * rSges(5,3);
t131 = cos(qJ(1));
t124 = t131 * pkin(1);
t148 = t123 * pkin(2) + t122 * qJ(3) + t124;
t147 = (-rSges(7,2) - pkin(8)) * t130;
t146 = (-rSges(6,3) - pkin(8)) * t130;
t145 = t161 * t130;
t144 = t123 * qJ(3) - t168;
t143 = t123 * pkin(7) + t148;
t142 = -rSges(6,1) * t87 - rSges(6,2) * t85;
t112 = pkin(4) * t158;
t138 = t112 + t143;
t135 = Icges(5,5) * t127 + Icges(5,6) * t130;
t134 = -t19 / 0.2e1 - t30 / 0.2e1 - t29 / 0.2e1 - t21 / 0.2e1;
t133 = t22 / 0.2e1 + t20 / 0.2e1 + t32 / 0.2e1 + t31 / 0.2e1;
t113 = t123 * t127 * pkin(4);
t132 = t175 * t122 + t113 + t144;
t111 = pkin(4) * t130 + pkin(8) * t127;
t108 = rSges(2,1) * t131 - t128 * rSges(2,2);
t106 = -t128 * rSges(2,1) - rSges(2,2) * t131;
t103 = -t187 + t188;
t100 = t122 * t111;
t98 = rSges(3,1) * t123 - rSges(3,2) * t122 + t124;
t97 = -rSges(3,1) * t122 - rSges(3,2) * t123 - t168;
t96 = -pkin(8) * t157 + t112;
t95 = rSges(6,3) * t127 + (rSges(6,1) * t129 - rSges(6,2) * t126) * t130;
t75 = t123 * (pkin(8) * t156 - t113);
t67 = Icges(5,3) * t122 - t135 * t123;
t66 = Icges(5,3) * t123 + t135 * t122;
t65 = -rSges(4,2) * t123 + rSges(4,3) * t122 + t148;
t64 = rSges(4,3) * t123 + (rSges(4,2) - pkin(2)) * t122 + t144;
t63 = (-t111 - t95) * t123;
t62 = t122 * t95 + t100;
t59 = t143 + t149;
t58 = t176 + (-rSges(5,3) + t175) * t122 + t144;
t57 = rSges(6,3) * t156 - t142;
t55 = -rSges(6,3) * t157 + t162;
t41 = (-t111 - t161) * t123;
t40 = t161 * t122 + t100;
t39 = -t122 * t149 + (t122 * rSges(5,3) - t176) * t123;
t36 = -t127 * t57 + t95 * t156;
t35 = t127 * t55 + t95 * t157;
t34 = t122 * t146 + t138 + t162;
t33 = t123 * t146 + t132 + t142;
t28 = (-t122 * t57 - t123 * t55) * t130;
t27 = t122 * t147 + t138 + t177;
t26 = t123 * t147 + t132 - t186;
t25 = t123 * t145 - t165 * t127;
t24 = t122 * t145 + t166 * t127;
t23 = t123 * t57 + t75 + (-t55 - t96) * t122;
t10 = (-t165 * t122 - t166 * t123) * t130;
t9 = t75 + t165 * t123 + (-t96 - t166) * t122;
t8 = t122 * t18 + t123 * t17;
t7 = t122 * t16 + t123 * t15;
t6 = t122 * t14 + t123 * t13;
t5 = t11 * t123 + t12 * t122;
t1 = [Icges(4,1) + Icges(2,3) + Icges(3,3) + (Icges(5,1) * t130 - t126 * t91) * t130 + m(7) * (t26 ^ 2 + t27 ^ 2) + m(6) * (t33 ^ 2 + t34 ^ 2) + m(5) * (t58 ^ 2 + t59 ^ 2) + m(4) * (t64 ^ 2 + t65 ^ 2) + m(3) * (t97 ^ 2 + t98 ^ 2) + m(2) * (t106 ^ 2 + t108 ^ 2) + t182 + (-0.2e1 * Icges(5,4) * t130 + Icges(5,2) * t127) * t127; 0; m(3) + m(4) + m(5) + m(6) + m(7); m(7) * (t122 * t26 - t123 * t27) + m(6) * (t122 * t33 - t123 * t34) + m(5) * (t122 * t58 - t123 * t59) + m(4) * (t122 * t64 - t123 * t65); 0; 0.2e1 * (m(4) / 0.2e1 + m(5) / 0.2e1 + m(6) / 0.2e1 + m(7) / 0.2e1) * t151; m(7) * (t26 * t40 + t27 * t41) + m(6) * (t33 * t62 + t34 * t63) + (t103 * t173 + t185 * t123 - t59 * t169 - t134) * t123 + (t103 * t174 + t185 * t122 + t58 * t169 + t133) * t122; m(5) * t39 + m(6) * t23 + m(7) * t9; m(6) * (t122 * t62 - t123 * t63) + m(7) * (t122 * t40 - t123 * t41) + t151 * t169; m(7) * (t40 ^ 2 + t41 ^ 2 + t9 ^ 2) + m(6) * (t23 ^ 2 + t62 ^ 2 + t63 ^ 2) + m(5) * (t151 * t107 ^ 2 + t39 ^ 2) + (t120 * t67 + t7 + t8) * t122 + (t121 * t66 + t5 + t6 + (t122 * t66 + t123 * t67) * t122) * t123; m(7) * (t24 * t27 + t25 * t26) + m(6) * (t33 * t36 + t34 * t35) + (t134 * t122 + t133 * t123) * t130 + t167; m(6) * t28 + m(7) * t10; m(6) * (t122 * t36 - t123 * t35) + m(7) * (t122 * t25 - t123 * t24); m(7) * (t10 * t9 + t24 * t41 + t25 * t40) + m(6) * (t23 * t28 + t35 * t63 + t36 * t62) + ((t7 / 0.2e1 + t8 / 0.2e1) * t123 + (-t5 / 0.2e1 - t6 / 0.2e1) * t122) * t130 + t180 * t174 + t181 * t173 + (t178 * t122 + t179 * t123) * t127 / 0.2e1; t167 * t127 + m(7) * (t10 ^ 2 + t24 ^ 2 + t25 ^ 2) + m(6) * (t28 ^ 2 + t35 ^ 2 + t36 ^ 2) + ((t178 * t127 + t180) * t123 + (-t179 * t127 - t181) * t122) * t130; m(7) * (t26 * t83 - t27 * t85); m(7) * t154; m(7) * (t122 * t83 + t123 * t85); m(7) * (t9 * t154 + t40 * t83 - t41 * t85); m(7) * (t10 * t154 - t24 * t85 + t25 * t83); m(7) * (t126 ^ 2 * t130 ^ 2 + t83 ^ 2 + t85 ^ 2);];
%% Postprocessing: Reshape Output
% From vec2symmat_6_matlab.m
res = [t1(1) t1(2) t1(4) t1(7) t1(11) t1(16); t1(2) t1(3) t1(5) t1(8) t1(12) t1(17); t1(4) t1(5) t1(6) t1(9) t1(13) t1(18); t1(7) t1(8) t1(9) t1(10) t1(14) t1(19); t1(11) t1(12) t1(13) t1(14) t1(15) t1(20); t1(16) t1(17) t1(18) t1(19) t1(20) t1(21);];
Mq  = res;
