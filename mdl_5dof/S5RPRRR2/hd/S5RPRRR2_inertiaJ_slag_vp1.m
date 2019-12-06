% Calculate joint inertia matrix for
% S5RPRRR2
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% pkin [9x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d1,d3,d4,d5,theta2]';
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
% Datum: 2019-12-05 18:13
% Revision: 77da58f92bca3eff71542919beafa37024070d86 (2019-12-05)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Mq = S5RPRRR2_inertiaJ_slag_vp1(qJ, ...
  pkin, m, rSges, Icges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(9,1),zeros(6,1),zeros(6,3),zeros(6,6)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RPRRR2_inertiaJ_slag_vp1: qJ has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S5RPRRR2_inertiaJ_slag_vp1: pkin has to be [9x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5RPRRR2_inertiaJ_slag_vp1: m has to be [6x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [6,3]), ...
  'S5RPRRR2_inertiaJ_slag_vp1: rSges has to be [6x3] (double)');
assert(isreal(Icges) && all(size(Icges) == [6 6]), ...
  'S5RPRRR2_inertiaJ_slag_vp1: Icges has to be [6x6] (double)'); 

%% Symbolic Calculation
% From inertia_joint_joint_floatb_twist_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-05 18:11:16
% EndTime: 2019-12-05 18:11:21
% DurationCPUTime: 1.10s
% Computational Cost: add. (3083->218), mult. (2260->325), div. (0->0), fcn. (2040->10), ass. (0->116)
t106 = sin(qJ(1));
t101 = t106 ^ 2;
t100 = pkin(9) + qJ(3);
t93 = qJ(4) + t100;
t90 = qJ(5) + t93;
t83 = sin(t90);
t84 = cos(t90);
t131 = rSges(6,1) * t84 - rSges(6,2) * t83;
t87 = sin(t93);
t88 = cos(t93);
t132 = rSges(5,1) * t88 - rSges(5,2) * t87;
t107 = cos(qJ(1));
t102 = t107 ^ 2;
t91 = sin(t100);
t158 = pkin(3) * t91;
t157 = t106 / 0.2e1;
t156 = -t107 / 0.2e1;
t104 = cos(pkin(9));
t89 = t104 * pkin(2) + pkin(1);
t92 = cos(t100);
t155 = rSges(4,1) * t92;
t152 = rSges(4,2) * t91;
t105 = -pkin(6) - qJ(2);
t77 = pkin(3) * t92 + t89;
t71 = t107 * t77;
t149 = t107 * (-t107 * t89 + t71) + (t77 - t89) * t101;
t111 = t106 * rSges(6,3) + t131 * t107;
t16 = t106 * (-t107 * rSges(6,3) + t131 * t106) + t107 * t111;
t112 = t106 * rSges(5,3) + t132 * t107;
t19 = t106 * (-t107 * rSges(5,3) + t132 * t106) + t107 * t112;
t148 = t106 * rSges(4,3) + t107 * t155;
t146 = Icges(4,4) * t91;
t145 = Icges(4,4) * t92;
t144 = Icges(5,4) * t87;
t143 = Icges(5,4) * t88;
t142 = Icges(6,4) * t83;
t141 = Icges(6,4) * t84;
t140 = rSges(3,3) + qJ(2);
t138 = t101 + t102;
t99 = -pkin(7) + t105;
t70 = t87 * rSges(5,1) + t88 * rSges(5,2);
t137 = -t70 - t158;
t65 = t83 * rSges(6,1) + t84 * rSges(6,2);
t136 = -pkin(4) * t87 - t65;
t117 = -Icges(6,2) * t83 + t141;
t42 = Icges(6,6) * t106 + t117 * t107;
t120 = Icges(6,1) * t84 - t142;
t44 = Icges(6,5) * t106 + t120 * t107;
t129 = -t42 * t83 + t44 * t84;
t41 = -Icges(6,6) * t107 + t117 * t106;
t43 = -Icges(6,5) * t107 + t120 * t106;
t130 = t41 * t83 - t43 * t84;
t114 = Icges(6,5) * t84 - Icges(6,6) * t83;
t39 = -Icges(6,3) * t107 + t114 * t106;
t40 = Icges(6,3) * t106 + t114 * t107;
t2 = t106 * (t101 * t40 + (t130 * t107 + (t129 - t39) * t106) * t107);
t3 = t102 * t39 + (t129 * t106 + (t130 - t40) * t107) * t106;
t135 = -t107 * t3 + t2;
t63 = Icges(6,2) * t84 + t142;
t64 = Icges(6,1) * t83 + t141;
t124 = -t63 * t83 + t64 * t84;
t62 = Icges(6,5) * t83 + Icges(6,6) * t84;
t134 = (t106 * t62 + t124 * t107 + t84 * t42 + t83 * t44) * t157 + (t124 * t106 - t107 * t62 + t84 * t41 + t83 * t43) * t156;
t60 = pkin(4) * t88 + t77;
t59 = t107 * t60;
t6 = t107 * (t59 - t71) + t16 + (t60 - t77) * t101;
t133 = -t152 + t155;
t118 = -Icges(5,2) * t87 + t143;
t49 = -Icges(5,6) * t107 + t118 * t106;
t121 = Icges(5,1) * t88 - t144;
t51 = -Icges(5,5) * t107 + t121 * t106;
t128 = t49 * t87 - t51 * t88;
t50 = Icges(5,6) * t106 + t118 * t107;
t52 = Icges(5,5) * t106 + t121 * t107;
t127 = -t50 * t87 + t52 * t88;
t68 = Icges(5,2) * t88 + t144;
t69 = Icges(5,1) * t87 + t143;
t123 = -t68 * t87 + t69 * t88;
t122 = Icges(4,1) * t92 - t146;
t119 = -Icges(4,2) * t91 + t145;
t116 = Icges(4,5) * t92 - Icges(4,6) * t91;
t115 = Icges(5,5) * t88 - Icges(5,6) * t87;
t47 = -Icges(5,3) * t107 + t115 * t106;
t48 = Icges(5,3) * t106 + t115 * t107;
t4 = t106 * (t101 * t48 + (t128 * t107 + (t127 - t47) * t106) * t107);
t5 = t102 * t47 + (t127 * t106 + (t128 - t48) * t107) * t106;
t113 = t2 + t4 + (-t5 - t3) * t107;
t110 = t136 - t158;
t103 = sin(pkin(9));
t109 = rSges(3,1) * t104 - rSges(3,2) * t103 + pkin(1);
t67 = Icges(5,5) * t87 + Icges(5,6) * t88;
t108 = t134 + (t106 * t67 + t123 * t107 + t88 * t50 + t87 * t52) * t157 + (t123 * t106 - t107 * t67 + t88 * t49 + t87 * t51) * t156;
t94 = -pkin(8) + t99;
t81 = t107 * rSges(2,1) - t106 * rSges(2,2);
t80 = -t106 * rSges(2,1) - t107 * rSges(2,2);
t76 = t91 * rSges(4,1) + t92 * rSges(4,2);
t54 = Icges(4,3) * t106 + t116 * t107;
t53 = -Icges(4,3) * t107 + t116 * t106;
t46 = t140 * t106 + t109 * t107;
t45 = -t109 * t106 + t140 * t107;
t36 = t137 * t107;
t35 = t137 * t106;
t32 = t136 * t107;
t31 = t136 * t106;
t30 = -t106 * t105 + (t89 - t152) * t107 + t148;
t29 = (rSges(4,3) - t105) * t107 + (-t133 - t89) * t106;
t28 = t110 * t107;
t27 = t110 * t106;
t24 = -t106 * t99 + t112 + t71;
t23 = (rSges(5,3) - t99) * t107 + (-t132 - t77) * t106;
t22 = t107 * (-t107 * t152 + t148) + (-t107 * rSges(4,3) + t133 * t106) * t106;
t18 = -t106 * t94 + t111 + t59;
t17 = (rSges(6,3) - t94) * t107 + (-t131 - t60) * t106;
t7 = t19 + t149;
t1 = t6 + t149;
t8 = [Icges(3,2) * t104 ^ 2 + t84 * t63 + t83 * t64 + t88 * t68 + t87 * t69 + t92 * (Icges(4,2) * t92 + t146) + t91 * (Icges(4,1) * t91 + t145) + Icges(2,3) + (Icges(3,1) * t103 + 0.2e1 * Icges(3,4) * t104) * t103 + m(6) * (t17 ^ 2 + t18 ^ 2) + m(5) * (t23 ^ 2 + t24 ^ 2) + m(4) * (t29 ^ 2 + t30 ^ 2) + m(3) * (t45 ^ 2 + t46 ^ 2) + m(2) * (t80 ^ 2 + t81 ^ 2); m(6) * (t106 * t17 - t107 * t18) + m(5) * (t106 * t23 - t107 * t24) + m(4) * (t106 * t29 - t107 * t30) + m(3) * (t106 * t45 - t107 * t46); 0.2e1 * (m(3) / 0.2e1 + m(4) / 0.2e1 + m(5) / 0.2e1 + m(6) / 0.2e1) * t138; (t92 * (-Icges(4,6) * t107 + t119 * t106) + t91 * (-Icges(4,5) * t107 + t122 * t106)) * t156 + (t92 * (Icges(4,6) * t106 + t119 * t107) + t91 * (Icges(4,5) * t106 + t122 * t107)) * t157 + m(6) * (t28 * t17 + t27 * t18) + m(5) * (t36 * t23 + t35 * t24) + m(4) * (-t106 * t30 - t107 * t29) * t76 + (t102 / 0.2e1 + t101 / 0.2e1) * (Icges(4,5) * t91 + Icges(4,6) * t92) + t108; m(5) * (t36 * t106 - t35 * t107) + m(6) * (t28 * t106 - t27 * t107); m(6) * (t1 ^ 2 + t27 ^ 2 + t28 ^ 2) + m(5) * (t35 ^ 2 + t36 ^ 2 + t7 ^ 2) + t4 + t106 * t101 * t54 + m(4) * (t138 * t76 ^ 2 + t22 ^ 2) + t135 + (-t102 * t53 - t5 + (-t106 * t53 + t107 * t54) * t106) * t107; m(6) * (t32 * t17 + t31 * t18) + m(5) * (-t106 * t24 - t107 * t23) * t70 + t108; m(6) * (t32 * t106 - t31 * t107); m(6) * (t6 * t1 + t31 * t27 + t32 * t28) + m(5) * (t19 * t7 + (-t106 * t35 - t107 * t36) * t70) + t113; m(5) * (t138 * t70 ^ 2 + t19 ^ 2) + m(6) * (t31 ^ 2 + t32 ^ 2 + t6 ^ 2) + t113; m(6) * (-t106 * t18 - t107 * t17) * t65 + t134; 0; m(6) * (t16 * t1 + (-t106 * t27 - t107 * t28) * t65) + t135; m(6) * (t16 * t6 + (-t106 * t31 - t107 * t32) * t65) + t135; m(6) * (t138 * t65 ^ 2 + t16 ^ 2) + t135;];
%% Postprocessing: Reshape Output
% From vec2symmat_5_matlab.m
res = [t8(1), t8(2), t8(4), t8(7), t8(11); t8(2), t8(3), t8(5), t8(8), t8(12); t8(4), t8(5), t8(6), t8(9), t8(13); t8(7), t8(8), t8(9), t8(10), t8(14); t8(11), t8(12), t8(13), t8(14), t8(15);];
Mq = res;
