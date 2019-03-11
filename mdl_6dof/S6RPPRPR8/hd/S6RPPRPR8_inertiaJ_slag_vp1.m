% Calculate joint inertia matrix for
% S6RPPRPR8
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% pkin [9x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d4,d6,theta3]';
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
% Datum: 2019-03-09 01:56
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Mq = S6RPPRPR8_inertiaJ_slag_vp1(qJ, ...
  pkin, m, rSges, Icges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(9,1),zeros(7,1),zeros(7,3),zeros(7,6)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPPRPR8_inertiaJ_slag_vp1: qJ has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S6RPPRPR8_inertiaJ_slag_vp1: pkin has to be [9x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6RPPRPR8_inertiaJ_slag_vp1: m has to be [7x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [7,3]), ...
  'S6RPPRPR8_inertiaJ_slag_vp1: rSges has to be [7x3] (double)');
assert(isreal(Icges) && all(size(Icges) == [7 6]), ...
  'S6RPPRPR8_inertiaJ_slag_vp1: Icges has to be [7x6] (double)'); 

%% Symbolic Calculation
% From inertia_joint_joint_floatb_twist_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 01:55:20
% EndTime: 2019-03-09 01:55:24
% DurationCPUTime: 1.64s
% Computational Cost: add. (2386->249), mult. (3117->374), div. (0->0), fcn. (3111->8), ass. (0->129)
t187 = -Icges(5,4) - Icges(6,6);
t186 = Icges(5,1) + Icges(6,2);
t185 = Icges(5,2) + Icges(6,3);
t107 = pkin(9) + qJ(4);
t100 = cos(t107);
t184 = t187 * t100;
t99 = sin(t107);
t183 = t187 * t99;
t182 = Icges(6,1) + Icges(5,3);
t181 = t185 * t100 - t183;
t180 = t186 * t99 - t184;
t179 = (-Icges(6,4) + Icges(5,5)) * t99 + (-Icges(6,5) + Icges(5,6)) * t100;
t114 = sin(qJ(1));
t178 = -t114 / 0.2e1;
t167 = t114 / 0.2e1;
t116 = cos(qJ(1));
t177 = -t116 / 0.2e1;
t176 = t116 / 0.2e1;
t175 = t100 / 0.2e1;
t166 = rSges(6,2) * t99;
t118 = -t166 + (-rSges(6,3) - qJ(5)) * t100;
t102 = t116 * qJ(2);
t112 = -pkin(7) - qJ(3);
t110 = sin(pkin(9));
t165 = pkin(3) * t110;
t147 = t114 * t112 + t116 * t165 + t102;
t161 = t116 * t99;
t93 = pkin(4) * t161;
t139 = t93 + t147;
t21 = (-rSges(6,1) - pkin(1)) * t114 + t118 * t116 + t139;
t150 = t116 * pkin(1) + t114 * qJ(2);
t119 = -t116 * t112 + t114 * t165 + t150;
t162 = t114 * t99;
t92 = pkin(4) * t162;
t117 = t119 + t92;
t22 = t116 * rSges(6,1) + t114 * t118 + t117;
t174 = m(6) * (t114 * t21 - t116 * t22);
t115 = cos(qJ(6));
t151 = t115 * t116;
t113 = sin(qJ(6));
t153 = t114 * t113;
t66 = t100 * t151 - t153;
t152 = t114 * t115;
t154 = t113 * t116;
t67 = t100 * t154 + t152;
t138 = -t67 * rSges(7,1) - t66 * rSges(7,2);
t156 = qJ(5) * t100;
t16 = (-pkin(1) - pkin(5)) * t114 + (-t156 + (rSges(7,3) + pkin(8)) * t99) * t116 + t138 + t139;
t155 = t100 * t114;
t142 = qJ(5) * t155;
t105 = t116 * pkin(5);
t68 = -t100 * t152 - t154;
t69 = -t100 * t153 + t151;
t33 = t69 * rSges(7,1) + t68 * rSges(7,2) + rSges(7,3) * t162;
t170 = pkin(8) * t162 + t105 + t33;
t17 = t117 - t142 + t170;
t173 = m(7) * (t114 * t16 - t116 * t17);
t172 = t114 * t182 - t179 * t116;
t171 = t179 * t114 + t116 * t182;
t169 = (rSges(5,1) * t99 + rSges(5,2) * t100) * t116;
t108 = t114 ^ 2;
t109 = t116 ^ 2;
t81 = rSges(5,1) * t100 - rSges(5,2) * t99;
t168 = m(5) * t81;
t94 = t108 + t109;
t84 = m(5) * t94;
t160 = rSges(4,3) + qJ(3);
t149 = m(6) / 0.2e1 + m(7) / 0.2e1;
t43 = Icges(7,3) * t100 + (Icges(7,5) * t113 + Icges(7,6) * t115) * t99;
t44 = Icges(7,6) * t100 + (Icges(7,4) * t113 + Icges(7,2) * t115) * t99;
t45 = Icges(7,5) * t100 + (Icges(7,1) * t113 + Icges(7,4) * t115) * t99;
t148 = t100 * t43 + (t113 * t45 + t115 * t44) * t99;
t146 = rSges(5,1) * t162 + rSges(5,2) * t155 + t116 * rSges(5,3);
t26 = Icges(7,5) * t67 + Icges(7,6) * t66 - Icges(7,3) * t161;
t28 = Icges(7,4) * t67 + Icges(7,2) * t66 - Icges(7,6) * t161;
t30 = Icges(7,1) * t67 + Icges(7,4) * t66 - Icges(7,5) * t161;
t10 = t100 * t26 + (t113 * t30 + t115 * t28) * t99;
t12 = -t161 * t43 + t66 * t44 + t67 * t45;
t145 = -t10 / 0.2e1 - t12 / 0.2e1;
t27 = Icges(7,5) * t69 + Icges(7,6) * t68 + Icges(7,3) * t162;
t29 = Icges(7,4) * t69 + Icges(7,2) * t68 + Icges(7,6) * t162;
t31 = Icges(7,1) * t69 + Icges(7,4) * t68 + Icges(7,5) * t162;
t11 = t100 * t27 + (t113 * t31 + t115 * t29) * t99;
t13 = t162 * t43 + t44 * t68 + t45 * t69;
t144 = t11 / 0.2e1 + t13 / 0.2e1;
t143 = Icges(5,5) * t175 - Icges(6,4) * t100 / 0.2e1 + (-Icges(5,6) / 0.2e1 + Icges(6,5) / 0.2e1) * t99;
t46 = t100 * rSges(7,3) + (rSges(7,1) * t113 + rSges(7,2) * t115) * t99;
t141 = pkin(8) * t100 + t46;
t140 = t84 + (m(4) + m(6) + m(7)) * t94;
t136 = rSges(6,3) * t100 + t166;
t111 = cos(pkin(9));
t131 = rSges(4,1) * t110 + rSges(4,2) * t111;
t19 = t100 * t33 - t162 * t46;
t32 = -rSges(7,3) * t161 - t138;
t20 = -t100 * t32 - t161 * t46;
t129 = t20 * t114 - t116 * t19;
t79 = pkin(4) * t100 + qJ(5) * t99;
t70 = t114 * t79;
t23 = t114 * t141 + t70;
t24 = (-t141 - t79) * t116;
t127 = t23 * t114 - t116 * t24;
t80 = -rSges(6,2) * t100 + rSges(6,3) * t99;
t36 = t114 * t80 + t70;
t37 = (-t79 - t80) * t116;
t126 = t36 * t114 - t116 * t37;
t87 = rSges(2,1) * t116 - t114 * rSges(2,2);
t86 = -t114 * rSges(2,1) - rSges(2,2) * t116;
t65 = -rSges(3,2) * t116 + t114 * rSges(3,3) + t150;
t64 = t116 * rSges(3,3) + t102 + (rSges(3,2) - pkin(1)) * t114;
t61 = t92 - t142;
t60 = t116 * (t116 * t156 - t93);
t41 = t114 * t131 + t116 * t160 + t150;
t40 = t102 + t131 * t116 + (-pkin(1) - t160) * t114;
t35 = t119 + t146;
t34 = t169 + (-rSges(5,3) - pkin(1)) * t114 + t147;
t25 = -t114 * t146 + (t114 * rSges(5,3) - t169) * t116;
t18 = t60 + t136 * t109 + (t114 * t136 - t61) * t114;
t15 = t148 * t100;
t14 = (t114 * t32 + t116 * t33) * t99;
t9 = t60 + (-pkin(8) * t161 + t32) * t116 + (-t61 + t105 - t170) * t114;
t8 = t162 * t27 + t29 * t68 + t31 * t69;
t7 = t162 * t26 + t28 * t68 + t30 * t69;
t6 = -t161 * t27 + t66 * t29 + t67 * t31;
t5 = -t161 * t26 + t66 * t28 + t67 * t30;
t4 = t7 * t114 + t116 * t8;
t3 = t5 * t114 + t116 * t6;
t2 = t13 * t100 + (t114 * t8 - t116 * t7) * t99;
t1 = t12 * t100 + (t114 * t6 - t116 * t5) * t99;
t38 = [Icges(4,1) * t111 ^ 2 + Icges(3,1) + Icges(2,3) + (-0.2e1 * Icges(4,4) * t111 + Icges(4,2) * t110) * t110 + m(7) * (t16 ^ 2 + t17 ^ 2) + m(6) * (t21 ^ 2 + t22 ^ 2) + m(5) * (t34 ^ 2 + t35 ^ 2) + m(4) * (t40 ^ 2 + t41 ^ 2) + m(3) * (t64 ^ 2 + t65 ^ 2) + m(2) * (t86 ^ 2 + t87 ^ 2) + t148 + (t185 * t99 + t184) * t99 + (t186 * t100 + t183) * t100; t173 + t174 + m(5) * (t114 * t34 - t116 * t35) + m(4) * (t114 * t40 - t116 * t41) + m(3) * (t114 * t64 - t116 * t65); m(3) * t94 + t140; m(7) * (t114 * t17 + t116 * t16) + m(6) * (t114 * t22 + t116 * t21) + m(5) * (t114 * t35 + t116 * t34) + m(4) * (t114 * t41 + t116 * t40); 0; t140; m(7) * (t16 * t23 + t17 * t24) + m(6) * (t21 * t36 + t22 * t37) + (-t35 * t168 + t143 * t116 + (Icges(6,5) * t176 + Icges(5,6) * t177 + t178 * t181) * t99 + t144) * t116 + (t34 * t168 + t143 * t114 + (Icges(6,5) * t167 + Icges(5,6) * t178 + t176 * t181) * t99 - t145) * t114 + ((Icges(6,4) * t177 + Icges(5,5) * t176 + t167 * t180) * t116 + (Icges(6,4) * t178 + Icges(5,5) * t167 + t177 * t180) * t114) * t100; m(6) * t126 + m(7) * t127 + t81 * t84; m(6) * (t37 * t114 + t116 * t36) + m(7) * (t24 * t114 + t116 * t23); m(7) * (t23 ^ 2 + t24 ^ 2 + t9 ^ 2) + m(6) * (t18 ^ 2 + t36 ^ 2 + t37 ^ 2) + m(5) * (t81 ^ 2 * t94 + t25 ^ 2) + (t171 * t109 + t4) * t116 + (t3 + t172 * t108 + (t171 * t114 + t172 * t116) * t116) * t114; 0.2e1 * (-t173 / 0.2e1 - t174 / 0.2e1) * t100; -0.2e1 * t149 * t94 * t100; 0; m(7) * (-t100 * t127 + t99 * t9) + m(6) * (-t100 * t126 + t99 * t18); 0.2e1 * t149 * (t100 ^ 2 * t94 + t99 ^ 2); m(7) * (t16 * t20 + t17 * t19) + t15 + (t114 * t144 + t116 * t145) * t99; m(7) * t129; m(7) * (t19 * t114 + t116 * t20); m(7) * (t14 * t9 + t19 * t24 + t20 * t23) + t1 * t167 + (t10 * t114 + t11 * t116) * t175 + t2 * t176 + (t4 * t167 + t3 * t177) * t99; m(7) * (-t100 * t129 + t14 * t99); t100 * t15 + m(7) * (t14 ^ 2 + t19 ^ 2 + t20 ^ 2) + (t114 * t2 - t116 * t1 + t100 * (-t10 * t116 + t11 * t114)) * t99;];
%% Postprocessing: Reshape Output
% From vec2symmat_6_matlab.m
res = [t38(1) t38(2) t38(4) t38(7) t38(11) t38(16); t38(2) t38(3) t38(5) t38(8) t38(12) t38(17); t38(4) t38(5) t38(6) t38(9) t38(13) t38(18); t38(7) t38(8) t38(9) t38(10) t38(14) t38(19); t38(11) t38(12) t38(13) t38(14) t38(15) t38(20); t38(16) t38(17) t38(18) t38(19) t38(20) t38(21);];
Mq  = res;
