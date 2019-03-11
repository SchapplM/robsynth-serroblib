% Calculate joint inertia matrix for
% S6RPPRPR3
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% pkin [10x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d4,d6,theta2,theta5]';
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
% Datum: 2019-03-09 01:45
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Mq = S6RPPRPR3_inertiaJ_slag_vp1(qJ, ...
  pkin, m, rSges, Icges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(10,1),zeros(7,1),zeros(7,3),zeros(7,6)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPPRPR3_inertiaJ_slag_vp1: qJ has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S6RPPRPR3_inertiaJ_slag_vp1: pkin has to be [10x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6RPPRPR3_inertiaJ_slag_vp1: m has to be [7x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [7,3]), ...
  'S6RPPRPR3_inertiaJ_slag_vp1: rSges has to be [7x3] (double)');
assert(isreal(Icges) && all(size(Icges) == [7 6]), ...
  'S6RPPRPR3_inertiaJ_slag_vp1: Icges has to be [7x6] (double)'); 

%% Symbolic Calculation
% From inertia_joint_joint_floatb_twist_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 01:44:00
% EndTime: 2019-03-09 01:44:02
% DurationCPUTime: 1.29s
% Computational Cost: add. (3199->230), mult. (2771->337), div. (0->0), fcn. (2804->10), ass. (0->119)
t104 = qJ(4) + pkin(10);
t99 = sin(t104);
t174 = Icges(6,6) * t99;
t111 = cos(qJ(4));
t173 = Icges(5,5) * t111;
t101 = cos(t104);
t172 = Icges(6,5) * t101;
t108 = sin(qJ(4));
t171 = Icges(5,6) * t108;
t170 = Icges(5,3) + Icges(6,3);
t169 = Icges(5,5) * t108 + Icges(6,5) * t99 + Icges(5,6) * t111 + Icges(6,6) * t101;
t168 = t173 / 0.2e1 + t172 / 0.2e1 - t171 / 0.2e1 - t174 / 0.2e1;
t167 = m(6) + m(7);
t105 = qJ(1) + pkin(9);
t100 = sin(t105);
t102 = cos(t105);
t166 = t169 * t100 + t170 * t102;
t165 = t170 * t100 - t169 * t102;
t164 = (rSges(5,1) * t108 + rSges(5,2) * t111) * t102;
t97 = t100 ^ 2;
t98 = t102 ^ 2;
t162 = pkin(5) * t99;
t160 = t102 / 0.2e1;
t157 = pkin(4) * t111;
t109 = sin(qJ(1));
t156 = t109 * pkin(1);
t110 = cos(qJ(6));
t107 = sin(qJ(6));
t49 = Icges(7,3) * t99 + (Icges(7,5) * t110 - Icges(7,6) * t107) * t101;
t55 = Icges(7,5) * t99 + (Icges(7,1) * t110 - Icges(7,4) * t107) * t101;
t155 = t101 * t110 * t55 + t99 * t49;
t58 = t99 * rSges(7,3) + (rSges(7,1) * t110 - rSges(7,2) * t107) * t101;
t154 = pkin(5) * t101 + pkin(8) * t99 + t58;
t137 = t102 * t110;
t143 = t100 * t107;
t62 = -t143 * t99 + t137;
t138 = t102 * t107;
t141 = t100 * t110;
t63 = t141 * t99 + t138;
t153 = t63 * rSges(7,1) + t62 * rSges(7,2);
t75 = t97 + t98;
t152 = t167 * t75;
t106 = -qJ(5) - pkin(7);
t151 = t102 * t108 * pkin(4) + t100 * t106;
t149 = t100 * t99;
t94 = t102 * rSges(6,3);
t52 = Icges(7,6) * t99 + (Icges(7,4) * t110 - Icges(7,2) * t107) * t101;
t148 = t107 * t52;
t144 = t100 * t101;
t142 = t100 * t108;
t140 = t100 * t111;
t139 = t101 * t102;
t136 = -rSges(6,1) * t149 - rSges(6,2) * t144 - t94;
t135 = rSges(5,1) * t142 + rSges(5,2) * t140 + t102 * rSges(5,3);
t112 = cos(qJ(1));
t103 = t112 * pkin(1);
t134 = t102 * pkin(2) + t100 * qJ(3) + t103;
t26 = Icges(7,5) * t63 + Icges(7,6) * t62 - Icges(7,3) * t144;
t28 = Icges(7,4) * t63 + Icges(7,2) * t62 - Icges(7,6) * t144;
t30 = Icges(7,1) * t63 + Icges(7,4) * t62 - Icges(7,5) * t144;
t10 = t99 * t26 + (-t107 * t28 + t110 * t30) * t101;
t13 = -t144 * t49 + t52 * t62 + t55 * t63;
t133 = -t13 / 0.2e1 - t10 / 0.2e1;
t64 = t138 * t99 + t141;
t65 = -t137 * t99 + t143;
t27 = Icges(7,5) * t65 + Icges(7,6) * t64 + Icges(7,3) * t139;
t29 = Icges(7,4) * t65 + Icges(7,2) * t64 + Icges(7,6) * t139;
t31 = Icges(7,1) * t65 + Icges(7,4) * t64 + Icges(7,5) * t139;
t11 = t99 * t27 + (-t107 * t29 + t110 * t31) * t101;
t14 = t139 * t49 + t52 * t64 + t55 * t65;
t132 = t14 / 0.2e1 + t11 / 0.2e1;
t131 = (-rSges(7,3) - pkin(8)) * t101;
t130 = t102 * qJ(3) - t156;
t128 = -rSges(7,1) * t65 - rSges(7,2) * t64;
t127 = rSges(6,1) * t99 + rSges(6,2) * t101;
t121 = t130 + t151;
t34 = t164 + (-rSges(5,3) - pkin(2) - pkin(7)) * t100 + t130;
t35 = pkin(7) * t102 + t134 + t135;
t114 = m(5) * (t100 * t34 - t102 * t35);
t88 = pkin(4) * t142;
t113 = -t102 * t106 + t134 + t88;
t89 = pkin(4) * t140;
t85 = rSges(2,1) * t112 - t109 * rSges(2,2);
t84 = rSges(5,1) * t111 - rSges(5,2) * t108;
t83 = -t109 * rSges(2,1) - rSges(2,2) * t112;
t78 = pkin(5) * t149;
t73 = rSges(6,1) * t101 - rSges(6,2) * t99;
t67 = rSges(3,1) * t102 - rSges(3,2) * t100 + t103;
t66 = -rSges(3,1) * t100 - rSges(3,2) * t102 - t156;
t59 = t88 + (-pkin(7) - t106) * t102;
t48 = (-t73 - t157) * t102;
t47 = t100 * t73 + t89;
t46 = t102 * (-pkin(7) * t100 - t151);
t39 = -rSges(4,2) * t102 + rSges(4,3) * t100 + t134;
t38 = t102 * rSges(4,3) + (rSges(4,2) - pkin(2)) * t100 + t130;
t33 = rSges(7,3) * t139 - t128;
t32 = -rSges(7,3) * t144 + t153;
t25 = (-t154 - t157) * t102;
t24 = t100 * t154 + t89;
t23 = -t100 * t135 + (t100 * rSges(5,3) - t164) * t102;
t22 = t113 - t136;
t21 = t127 * t102 + (-rSges(6,3) - pkin(2)) * t100 + t121;
t20 = t139 * t58 - t33 * t99;
t19 = t144 * t58 + t32 * t99;
t18 = (-t101 * t148 + t155) * t99;
t17 = t100 * t131 + t113 + t153 + t78;
t16 = -pkin(2) * t100 + (t131 + t162) * t102 + t121 + t128;
t15 = t46 - t127 * t98 + (t136 - t59 + t94) * t100;
t12 = (-t100 * t33 - t102 * t32) * t101;
t9 = t139 * t27 + t29 * t64 + t31 * t65;
t8 = t139 * t26 + t28 * t64 + t30 * t65;
t7 = -t144 * t27 + t29 * t62 + t31 * t63;
t6 = -t144 * t26 + t28 * t62 + t30 * t63;
t5 = t46 + (t33 + (pkin(8) * t101 - t162) * t102) * t102 + (pkin(8) * t144 - t32 - t59 - t78) * t100;
t4 = t100 * t9 + t102 * t8;
t3 = t100 * t7 + t102 * t6;
t2 = t14 * t99 + (-t100 * t8 + t102 * t9) * t101;
t1 = t13 * t99 + (-t100 * t6 + t102 * t7) * t101;
t36 = [Icges(5,1) * t111 ^ 2 + t99 ^ 2 * Icges(6,2) + Icges(4,1) + Icges(2,3) + Icges(3,3) + m(7) * (t16 ^ 2 + t17 ^ 2) + m(6) * (t21 ^ 2 + t22 ^ 2) + m(5) * (t34 ^ 2 + t35 ^ 2) + m(4) * (t38 ^ 2 + t39 ^ 2) + m(3) * (t66 ^ 2 + t67 ^ 2) + m(2) * (t83 ^ 2 + t85 ^ 2) + t155 + (-0.2e1 * Icges(5,4) * t111 + Icges(5,2) * t108) * t108 + (Icges(6,1) * t101 - 0.2e1 * Icges(6,4) * t99 - t148) * t101; 0; m(3) + m(4) + m(5) + t167; m(7) * (t100 * t16 - t102 * t17) + m(6) * (t100 * t21 - t102 * t22) + t114 + m(4) * (t100 * t38 - t102 * t39); 0; 0.2e1 * (m(4) / 0.2e1 + m(5) / 0.2e1) * t75 + t152; m(7) * (t16 * t24 + t17 * t25) + m(6) * (t21 * t47 + t22 * t48) + t84 * t114 + (-t171 + t172 + t173 - t174) * (t97 / 0.2e1 + t98 / 0.2e1) + (t168 * t102 - t133) * t102 + (t168 * t100 + t132) * t100; m(5) * t23 + m(6) * t15 + m(7) * t5; m(6) * (t100 * t47 - t102 * t48) + m(7) * (t100 * t24 - t102 * t25) + m(5) * t75 * t84; m(7) * (t24 ^ 2 + t25 ^ 2 + t5 ^ 2) + m(5) * (t75 * t84 ^ 2 + t23 ^ 2) + m(6) * (t15 ^ 2 + t47 ^ 2 + t48 ^ 2) + (t165 * t97 + t4) * t100 + (t3 + t166 * t98 + (t166 * t100 + t165 * t102) * t100) * t102; m(7) * (t100 * t17 + t102 * t16) + m(6) * (t100 * t22 + t102 * t21); 0; 0; m(7) * (t100 * t25 + t102 * t24) + m(6) * (t100 * t48 + t102 * t47); t152; m(7) * (t16 * t20 + t17 * t19) + t18 + (t100 * t133 + t102 * t132) * t101; m(7) * t12; m(7) * (t100 * t20 - t102 * t19); m(7) * (t12 * t5 + t19 * t25 + t20 * t24) + t99 * (t10 * t102 + t11 * t100) / 0.2e1 + t100 * t2 / 0.2e1 + t1 * t160 + (t4 * t160 - t100 * t3 / 0.2e1) * t101; m(7) * (t100 * t19 + t102 * t20); t99 * t18 + m(7) * (t12 ^ 2 + t19 ^ 2 + t20 ^ 2) + (-t100 * t1 + t102 * t2 + t99 * (-t10 * t100 + t102 * t11)) * t101;];
%% Postprocessing: Reshape Output
% From vec2symmat_6_matlab.m
res = [t36(1) t36(2) t36(4) t36(7) t36(11) t36(16); t36(2) t36(3) t36(5) t36(8) t36(12) t36(17); t36(4) t36(5) t36(6) t36(9) t36(13) t36(18); t36(7) t36(8) t36(9) t36(10) t36(14) t36(19); t36(11) t36(12) t36(13) t36(14) t36(15) t36(20); t36(16) t36(17) t36(18) t36(19) t36(20) t36(21);];
Mq  = res;
