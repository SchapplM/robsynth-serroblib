% Calculate joint inertia matrix for
% S5RRPRR6
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% pkin [9x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d1,d2,d4,d5,theta3]';
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
% Datum: 2019-12-05 18:36
% Revision: 77da58f92bca3eff71542919beafa37024070d86 (2019-12-05)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Mq = S5RRPRR6_inertiaJ_slag_vp1(qJ, ...
  pkin, m, rSges, Icges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(9,1),zeros(6,1),zeros(6,3),zeros(6,6)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RRPRR6_inertiaJ_slag_vp1: qJ has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S5RRPRR6_inertiaJ_slag_vp1: pkin has to be [9x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5RRPRR6_inertiaJ_slag_vp1: m has to be [6x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [6,3]), ...
  'S5RRPRR6_inertiaJ_slag_vp1: rSges has to be [6x3] (double)');
assert(isreal(Icges) && all(size(Icges) == [6 6]), ...
  'S5RRPRR6_inertiaJ_slag_vp1: Icges has to be [6x6] (double)'); 

%% Symbolic Calculation
% From inertia_joint_joint_floatb_twist_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-05 18:35:42
% EndTime: 2019-12-05 18:35:46
% DurationCPUTime: 1.34s
% Computational Cost: add. (4216->256), mult. (4188->375), div. (0->0), fcn. (4406->10), ass. (0->127)
t111 = sin(pkin(9));
t112 = cos(pkin(9));
t113 = sin(qJ(4));
t115 = cos(qJ(4));
t110 = qJ(1) + qJ(2);
t106 = sin(t110);
t138 = t106 * t111;
t108 = cos(t110);
t131 = t112 * t113;
t83 = t106 * t131 + t108 * t115;
t130 = t112 * t115;
t133 = t108 * t113;
t84 = -t106 * t130 + t133;
t48 = Icges(5,5) * t84 + Icges(5,6) * t83 - Icges(5,3) * t138;
t50 = Icges(5,4) * t84 + Icges(5,2) * t83 - Icges(5,6) * t138;
t52 = Icges(5,1) * t84 + Icges(5,4) * t83 - Icges(5,5) * t138;
t79 = -Icges(5,3) * t112 + (Icges(5,5) * t115 - Icges(5,6) * t113) * t111;
t80 = -Icges(5,6) * t112 + (Icges(5,4) * t115 - Icges(5,2) * t113) * t111;
t81 = -Icges(5,5) * t112 + (Icges(5,1) * t115 - Icges(5,4) * t113) * t111;
t155 = -t112 * t48 + (-t113 * t50 + t115 * t52) * t111 - t79 * t138 + t80 * t83 + t81 * t84;
t135 = t108 * t111;
t85 = t106 * t115 - t108 * t131;
t136 = t106 * t113;
t86 = t108 * t130 + t136;
t49 = Icges(5,5) * t86 + Icges(5,6) * t85 + Icges(5,3) * t135;
t51 = Icges(5,4) * t86 + Icges(5,2) * t85 + Icges(5,6) * t135;
t53 = Icges(5,1) * t86 + Icges(5,4) * t85 + Icges(5,5) * t135;
t154 = -t112 * t49 + (-t113 * t51 + t115 * t53) * t111 + t79 * t135 + t80 * t85 + t81 * t86;
t104 = pkin(4) * t115 + pkin(3);
t147 = pkin(3) - t104;
t153 = pkin(7) * t111 + t147 * t112;
t152 = t106 ^ 2;
t151 = t108 ^ 2;
t114 = sin(qJ(1));
t149 = t114 * pkin(1);
t116 = cos(qJ(1));
t148 = t116 * pkin(1);
t109 = qJ(4) + qJ(5);
t105 = sin(t109);
t107 = cos(t109);
t134 = t108 * t112;
t75 = -t105 * t134 + t106 * t107;
t76 = t105 * t106 + t107 * t134;
t124 = -t76 * rSges(6,1) - t75 * rSges(6,2);
t47 = rSges(6,3) * t135 - t124;
t72 = -t112 * rSges(6,3) + (rSges(6,1) * t107 - rSges(6,2) * t105) * t111;
t30 = t112 * t47 + t72 * t135;
t117 = -pkin(8) - pkin(7);
t132 = t111 * t117;
t140 = pkin(4) * t133 + t106 * t132;
t137 = t106 * t112;
t73 = t105 * t137 + t107 * t108;
t74 = t105 * t108 - t107 * t137;
t145 = t74 * rSges(6,1) + t73 * rSges(6,2);
t46 = -rSges(6,3) * t138 + t145;
t146 = -t153 * t106 - t140 - t46;
t144 = t84 * rSges(5,1) + t83 * rSges(5,2);
t70 = -Icges(6,6) * t112 + (Icges(6,4) * t107 - Icges(6,2) * t105) * t111;
t143 = t105 * t70;
t142 = t113 * t80;
t71 = -Icges(6,5) * t112 + (Icges(6,1) * t107 - Icges(6,4) * t105) * t111;
t62 = t111 * t107 * t71;
t69 = -Icges(6,3) * t112 + (Icges(6,5) * t107 - Icges(6,6) * t105) * t111;
t31 = -t111 * t143 - t112 * t69 + t62;
t141 = t31 * t112;
t139 = t106 * t108;
t129 = -t138 / 0.2e1;
t128 = t135 / 0.2e1;
t18 = -t69 * t138 + t70 * t73 + t71 * t74;
t19 = t69 * t135 + t70 * t75 + t71 * t76;
t40 = Icges(6,5) * t74 + Icges(6,6) * t73 - Icges(6,3) * t138;
t42 = Icges(6,4) * t74 + Icges(6,2) * t73 - Icges(6,6) * t138;
t44 = Icges(6,1) * t74 + Icges(6,4) * t73 - Icges(6,5) * t138;
t7 = -t112 * t40 + (-t105 * t42 + t107 * t44) * t111;
t41 = Icges(6,5) * t76 + Icges(6,6) * t75 + Icges(6,3) * t135;
t43 = Icges(6,4) * t76 + Icges(6,2) * t75 + Icges(6,6) * t135;
t45 = Icges(6,1) * t76 + Icges(6,4) * t75 + Icges(6,5) * t135;
t8 = -t112 * t41 + (-t105 * t43 + t107 * t45) * t111;
t127 = (t18 + t7) * t129 + (t19 + t8) * t128;
t126 = -rSges(4,1) * t112 - pkin(2);
t92 = -rSges(3,1) * t108 + t106 * rSges(3,2);
t125 = -t86 * rSges(5,1) - t85 * rSges(5,2);
t91 = -rSges(3,1) * t106 - rSges(3,2) * t108;
t123 = t127 - t141;
t122 = -rSges(6,3) * t111 - t104 * t112 - pkin(2);
t1 = (-t19 * t112 - (t40 * t135 + t42 * t75 + t44 * t76) * t138 + (t41 * t135 + t43 * t75 + t45 * t76) * t135) * t135;
t2 = -t18 * t112 - (-t40 * t138 + t42 * t73 + t44 * t74) * t138 + (-t41 * t138 + t43 * t73 + t45 * t74) * t135;
t3 = -t141 + (-t106 * t7 + t108 * t8) * t111;
t121 = -t112 * t3 - t2 * t138 + t1;
t120 = -pkin(3) * t112 - pkin(2) + (-rSges(5,3) - pkin(7)) * t111;
t101 = t108 * qJ(3);
t60 = rSges(4,2) * t138 + t108 * rSges(4,3) + t126 * t106 + t101;
t61 = rSges(4,2) * t135 + t126 * t108 + (-rSges(4,3) - qJ(3)) * t106;
t65 = t111 * t115 * t81;
t38 = -t111 * t142 - t112 * t79 + t65;
t119 = (-t38 - t31) * t112 + t127 + t155 * t129 + t154 * t128;
t36 = t120 * t106 + t101 + t144;
t27 = t122 * t106 + t101 + t140 + t145;
t118 = Icges(3,3) + t62 + t65 + (Icges(4,2) * t112 - t69 - t79) * t112 + (Icges(4,1) * t111 + 0.2e1 * Icges(4,4) * t112 - t142 - t143) * t111;
t37 = -t106 * qJ(3) + t120 * t108 + t125;
t95 = t108 * t132;
t28 = t95 + (-pkin(4) * t113 - qJ(3)) * t106 + t122 * t108 + t124;
t97 = -rSges(2,1) * t116 + rSges(2,2) * t114;
t96 = -rSges(2,1) * t114 - rSges(2,2) * t116;
t88 = t92 - t148;
t87 = t91 - t149;
t82 = -t112 * rSges(5,3) + (rSges(5,1) * t115 - rSges(5,2) * t113) * t111;
t66 = (pkin(7) + t117) * t112 - t147 * t111;
t63 = t72 * t138;
t59 = t61 - t148;
t58 = t60 - t149;
t57 = rSges(5,3) * t135 - t125;
t56 = -rSges(5,3) * t138 + t144;
t55 = pkin(4) * t136 - t153 * t108 - t95;
t35 = t112 * t57 + t82 * t135;
t34 = -t112 * t56 + t82 * t138;
t33 = t37 - t148;
t32 = t36 - t149;
t29 = -t112 * t46 + t63;
t26 = t28 - t148;
t25 = t27 - t149;
t20 = (-t106 * t57 - t108 * t56) * t111;
t17 = (-t106 * t47 - t108 * t46) * t111;
t10 = t112 * t55 + t66 * t135 + t30;
t9 = t146 * t112 + t66 * t138 + t63;
t4 = (t146 * t108 + (-t47 - t55) * t106) * t111;
t5 = [Icges(2,3) + m(2) * (t96 ^ 2 + t97 ^ 2) + m(3) * (t87 ^ 2 + t88 ^ 2) + m(4) * (t58 ^ 2 + t59 ^ 2) + m(5) * (t32 ^ 2 + t33 ^ 2) + m(6) * (t25 ^ 2 + t26 ^ 2) + t118; m(3) * (t87 * t91 + t88 * t92) + m(4) * (t58 * t60 + t59 * t61) + m(5) * (t32 * t36 + t33 * t37) + m(6) * (t25 * t27 + t26 * t28) + t118; m(6) * (t27 ^ 2 + t28 ^ 2) + m(5) * (t36 ^ 2 + t37 ^ 2) + m(4) * (t60 ^ 2 + t61 ^ 2) + m(3) * (t91 ^ 2 + t92 ^ 2) + t118; m(4) * (t106 * t58 + t108 * t59) + m(5) * (t106 * t32 + t108 * t33) + m(6) * (t106 * t25 + t108 * t26); m(6) * (t106 * t27 + t108 * t28) + m(5) * (t106 * t36 + t108 * t37) + m(4) * (t106 * t60 + t108 * t61); 0.2e1 * (m(4) / 0.2e1 + m(5) / 0.2e1 + m(6) / 0.2e1) * (t151 + t152); m(5) * (t32 * t34 + t33 * t35) + m(6) * (t10 * t26 + t25 * t9) + t119; m(6) * (t10 * t28 + t27 * t9) + m(5) * (t34 * t36 + t35 * t37) + t119; m(5) * (t106 * t34 + t108 * t35) + m(6) * (t10 * t108 + t106 * t9); t1 + (t38 * t112 - t3) * t112 + (-t106 * t2 + (-t106 * (-(t50 * t83 + t52 * t84) * t106 + (t51 * t83 + t53 * t84) * t108) + t108 * (-(t50 * t85 + t52 * t86) * t106 + (t51 * t85 + t53 * t86) * t108) + (-t106 * (-t49 * t139 + t152 * t48) + t108 * (-t48 * t139 + t151 * t49)) * t111) * t111 + (t155 * t106 - t154 * t108) * t112) * t111 + m(5) * (t20 ^ 2 + t34 ^ 2 + t35 ^ 2) + m(6) * (t10 ^ 2 + t4 ^ 2 + t9 ^ 2); m(6) * (t25 * t29 + t26 * t30) + t123; m(6) * (t27 * t29 + t28 * t30) + t123; m(6) * (t106 * t29 + t108 * t30); m(6) * (t10 * t30 + t17 * t4 + t29 * t9) + t121; m(6) * (t17 ^ 2 + t29 ^ 2 + t30 ^ 2) + t121;];
%% Postprocessing: Reshape Output
% From vec2symmat_5_matlab.m
res = [t5(1), t5(2), t5(4), t5(7), t5(11); t5(2), t5(3), t5(5), t5(8), t5(12); t5(4), t5(5), t5(6), t5(9), t5(13); t5(7), t5(8), t5(9), t5(10), t5(14); t5(11), t5(12), t5(13), t5(14), t5(15);];
Mq = res;
