% Calculate joint inertia matrix for
% S4RRPR9
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [4x1]
%   Generalized joint coordinates (joint angles)
% pkin [7x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,d1,d2,d4,theta3]';
% m_mdh [5x1]
%   mass of all robot links (including the base)
% rSges [5x3]
%   center of mass of all robot links (in body frames)
%   rows: links of the robot (starting with base)
%   columns: x-, y-, z-coordinates
% Icges [5x6]
%   inertia of all robot links about their respective center of mass, in body frames
%   rows: links of the robot (starting with base)
%   columns: xx, yy, zz, xy, xz, yz (see inertiavector2matrix.m)
% 
% Output:
% Mq [4x4]
%   inertia matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-31 17:10
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Mq = S4RRPR9_inertiaJ_slag_vp1(qJ, ...
  pkin, m, rSges, Icges)
%% Coder Information
%#codegen
%$cgargs {zeros(4,1),zeros(7,1),zeros(5,1),zeros(5,3),zeros(5,6)}
assert(isreal(qJ) && all(size(qJ) == [4 1]), ...
  'S4RRPR9_inertiaJ_slag_vp1: qJ has to be [4x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [7 1]), ...
  'S4RRPR9_inertiaJ_slag_vp1: pkin has to be [7x1] (double)');
assert(isreal(m) && all(size(m) == [5 1]), ...
  'S4RRPR9_inertiaJ_slag_vp1: m has to be [5x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [5,3]), ...
  'S4RRPR9_inertiaJ_slag_vp1: rSges has to be [5x3] (double)');
assert(isreal(Icges) && all(size(Icges) == [5 6]), ...
  'S4RRPR9_inertiaJ_slag_vp1: Icges has to be [5x6] (double)'); 

%% Symbolic Calculation
% From inertia_joint_joint_floatb_twist_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 17:09:21
% EndTime: 2019-12-31 17:09:24
% DurationCPUTime: 1.15s
% Computational Cost: add. (1814->244), mult. (2801->373), div. (0->0), fcn. (2969->8), ass. (0->118)
t105 = sin(qJ(1));
t144 = t105 / 0.2e1;
t107 = cos(qJ(1));
t142 = t107 / 0.2e1;
t103 = -pkin(6) - qJ(3);
t106 = cos(qJ(2));
t126 = t106 * t107;
t101 = sin(pkin(7));
t128 = t105 * t101;
t104 = sin(qJ(2));
t129 = t104 * t107;
t97 = pkin(7) + qJ(4);
t91 = sin(t97);
t92 = cos(t97);
t58 = t105 * t92 - t91 * t126;
t59 = t105 * t91 + t92 * t126;
t32 = t59 * rSges(5,1) + t58 * rSges(5,2) + rSges(5,3) * t129;
t102 = cos(pkin(7));
t90 = t102 * pkin(3) + pkin(2);
t148 = pkin(3) * t128 - t103 * t129 + t90 * t126 + t32;
t99 = t105 ^ 2;
t147 = t106 ^ 2;
t100 = t107 ^ 2;
t146 = m(4) / 0.2e1;
t145 = m(5) / 0.2e1;
t143 = -t107 / 0.2e1;
t141 = pkin(2) * t106;
t47 = -Icges(5,6) * t106 + (Icges(5,4) * t92 - Icges(5,2) * t91) * t104;
t140 = t91 * t47;
t80 = t104 * pkin(2) - t106 * qJ(3);
t139 = t106 * rSges(4,3) - (rSges(4,1) * t102 - rSges(4,2) * t101) * t104 - t80;
t115 = qJ(3) * t104 + t141;
t137 = pkin(2) * t126 + qJ(3) * t129;
t138 = t107 * t137 + t99 * t115;
t136 = t107 * pkin(1) + t105 * pkin(5);
t135 = t106 * t90;
t134 = t107 * rSges(3,3);
t133 = t100 + t99;
t132 = Icges(3,4) * t104;
t131 = Icges(3,4) * t106;
t130 = t104 * t105;
t127 = t105 * t106;
t125 = t107 * t101;
t124 = t107 * t102;
t49 = -t106 * rSges(5,3) + (rSges(5,1) * t92 - rSges(5,2) * t91) * t104;
t123 = -(qJ(3) + t103) * t106 - (-pkin(2) + t90) * t104 - t49 - t80;
t73 = t105 * t102 - t106 * t125;
t74 = t106 * t124 + t128;
t122 = t74 * rSges(4,1) + t73 * rSges(4,2) + rSges(4,3) * t129;
t121 = pkin(3) * t125;
t56 = -t107 * t92 - t91 * t127;
t57 = -t107 * t91 + t92 * t127;
t25 = Icges(5,5) * t57 + Icges(5,6) * t56 + Icges(5,3) * t130;
t27 = Icges(5,4) * t57 + Icges(5,2) * t56 + Icges(5,6) * t130;
t29 = Icges(5,1) * t57 + Icges(5,4) * t56 + Icges(5,5) * t130;
t10 = -t106 * t25 + (-t27 * t91 + t29 * t92) * t104;
t46 = -Icges(5,3) * t106 + (Icges(5,5) * t92 - Icges(5,6) * t91) * t104;
t48 = -Icges(5,5) * t106 + (Icges(5,1) * t92 - Icges(5,4) * t91) * t104;
t13 = t46 * t130 + t56 * t47 + t57 * t48;
t120 = t10 / 0.2e1 + t13 / 0.2e1;
t26 = Icges(5,5) * t59 + Icges(5,6) * t58 + Icges(5,3) * t129;
t28 = Icges(5,4) * t59 + Icges(5,2) * t58 + Icges(5,6) * t129;
t30 = Icges(5,1) * t59 + Icges(5,4) * t58 + Icges(5,5) * t129;
t11 = -t106 * t26 + (-t28 * t91 + t30 * t92) * t104;
t14 = t46 * t129 + t58 * t47 + t59 * t48;
t119 = t11 / 0.2e1 + t14 / 0.2e1;
t71 = -t101 * t127 - t124;
t72 = t102 * t127 - t125;
t118 = -t72 * rSges(4,1) - t71 * rSges(4,2);
t117 = -t57 * rSges(5,1) - t56 * rSges(5,2);
t116 = rSges(3,1) * t106 - rSges(3,2) * t104;
t112 = Icges(3,1) * t106 - t132;
t111 = -Icges(3,2) * t104 + t131;
t110 = Icges(3,5) * t106 - Icges(3,6) * t104;
t109 = rSges(3,1) * t126 - rSges(3,2) * t129 + t105 * rSges(3,3);
t95 = t107 * pkin(5);
t83 = t107 * rSges(2,1) - t105 * rSges(2,2);
t82 = -t105 * rSges(2,1) - t107 * rSges(2,2);
t81 = t104 * rSges(3,1) + t106 * rSges(3,2);
t77 = Icges(3,5) * t104 + Icges(3,6) * t106;
t61 = Icges(3,3) * t105 + t110 * t107;
t60 = -Icges(3,3) * t107 + t110 * t105;
t54 = -Icges(4,5) * t106 + (Icges(4,1) * t102 - Icges(4,4) * t101) * t104;
t53 = -Icges(4,6) * t106 + (Icges(4,4) * t102 - Icges(4,2) * t101) * t104;
t44 = t109 + t136;
t43 = t134 + t95 + (-pkin(1) - t116) * t105;
t42 = t104 * t92 * t48;
t41 = t139 * t107;
t40 = t139 * t105;
t39 = Icges(4,1) * t74 + Icges(4,4) * t73 + Icges(4,5) * t129;
t38 = Icges(4,1) * t72 + Icges(4,4) * t71 + Icges(4,5) * t130;
t37 = Icges(4,4) * t74 + Icges(4,2) * t73 + Icges(4,6) * t129;
t36 = Icges(4,4) * t72 + Icges(4,2) * t71 + Icges(4,6) * t130;
t35 = Icges(4,5) * t74 + Icges(4,6) * t73 + Icges(4,3) * t129;
t34 = Icges(4,5) * t72 + Icges(4,6) * t71 + Icges(4,3) * t130;
t33 = t107 * t109 + (t116 * t105 - t134) * t105;
t31 = rSges(5,3) * t130 - t117;
t24 = t122 + t136 + t137;
t23 = t95 + (-t141 - pkin(1) + (-rSges(4,3) - qJ(3)) * t104) * t105 + t118;
t22 = t123 * t107;
t21 = t123 * t105;
t20 = -t106 * t32 - t49 * t129;
t19 = t106 * t31 + t49 * t130;
t18 = t136 + t148;
t17 = t121 + t95 + (-t135 - pkin(1) + (-rSges(5,3) + t103) * t104) * t105 + t117;
t16 = -t104 * t140 - t106 * t46 + t42;
t15 = (-t105 * t32 + t107 * t31) * t104;
t12 = t105 * (rSges(4,3) * t130 - t118) + t107 * t122 + t138;
t9 = t26 * t129 + t58 * t28 + t59 * t30;
t8 = t25 * t129 + t58 * t27 + t59 * t29;
t7 = t26 * t130 + t56 * t28 + t57 * t30;
t6 = t25 * t130 + t56 * t27 + t57 * t29;
t5 = (-t137 + t148) * t107 + (-t121 + t31 + (-t103 * t104 - t115 + t135) * t105) * t105 + t138;
t4 = t9 * t105 - t8 * t107;
t3 = t7 * t105 - t6 * t107;
t2 = -t14 * t106 + (t105 * t8 + t107 * t9) * t104;
t1 = -t13 * t106 + (t105 * t6 + t107 * t7) * t104;
t45 = [Icges(2,3) + t42 + (-t46 + t132 - (Icges(4,5) * t102 - Icges(4,6) * t101) * t104 + (Icges(3,2) + Icges(4,3)) * t106) * t106 + (Icges(3,1) * t104 - t101 * t53 + t102 * t54 + t131 - t140) * t104 + m(5) * (t17 ^ 2 + t18 ^ 2) + m(4) * (t23 ^ 2 + t24 ^ 2) + m(2) * (t82 ^ 2 + t83 ^ 2) + m(3) * (t43 ^ 2 + t44 ^ 2); (-t71 * t53 / 0.2e1 - t72 * t54 / 0.2e1 + t77 * t142 - t120) * t107 + (t73 * t53 / 0.2e1 + t74 * t54 / 0.2e1 + t77 * t144 + t119) * t105 + m(5) * (t22 * t17 + t21 * t18) + m(4) * (t41 * t23 + t40 * t24) + m(3) * (-t105 * t44 - t107 * t43) * t81 + ((Icges(3,6) * t142 - t111 * t105 / 0.2e1 + t34 / 0.2e1) * t107 + (Icges(3,6) * t144 + t111 * t142 - t35 / 0.2e1) * t105) * t106 + ((Icges(3,5) * t105 - t101 * t37 + t102 * t39 + t112 * t107) * t144 + (-Icges(3,5) * t107 - t101 * t36 + t102 * t38 + t112 * t105) * t143) * t104; m(5) * (t21 ^ 2 + t22 ^ 2 + t5 ^ 2) + m(4) * (t12 ^ 2 + t40 ^ 2 + t41 ^ 2) + m(3) * (t133 * t81 ^ 2 + t33 ^ 2) + (-t100 * t60 - t3 + (t34 * t130 + t71 * t36 + t72 * t38) * t107) * t107 + (t4 + t99 * t61 + (t35 * t129 + t73 * t37 + t74 * t39) * t105 + (-t105 * t60 + t107 * t61 - t34 * t129 - t35 * t130 - t73 * t36 - t71 * t37 - t74 * t38 - t72 * t39) * t107) * t105; 0.2e1 * ((t105 * t18 + t107 * t17) * t145 + (t105 * t24 + t107 * t23) * t146) * t104; m(5) * (-t106 * t5 + (t105 * t21 + t107 * t22) * t104) + m(4) * (-t106 * t12 + (t105 * t40 + t107 * t41) * t104); 0.2e1 * (t146 + t145) * (t133 * t104 ^ 2 + t147); m(5) * (t19 * t17 + t20 * t18) - t16 * t106 + (t120 * t105 + t119 * t107) * t104; m(5) * (t15 * t5 + t19 * t22 + t20 * t21) + t2 * t144 + t1 * t143 - t106 * (-t10 * t107 + t11 * t105) / 0.2e1 + (t4 * t142 + t3 * t144) * t104; m(5) * (-t15 * t106 + (t105 * t20 + t107 * t19) * t104); m(5) * (t15 ^ 2 + t19 ^ 2 + t20 ^ 2) + t147 * t16 + (t107 * t2 + t105 * t1 - t106 * (t10 * t105 + t107 * t11)) * t104;];
%% Postprocessing: Reshape Output
% From vec2symmat_4_matlab.m
res = [t45(1), t45(2), t45(4), t45(7); t45(2), t45(3), t45(5), t45(8); t45(4), t45(5), t45(6), t45(9); t45(7), t45(8), t45(9), t45(10);];
Mq = res;
