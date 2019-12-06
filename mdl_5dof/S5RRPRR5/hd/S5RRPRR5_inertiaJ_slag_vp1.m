% Calculate joint inertia matrix for
% S5RRPRR5
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
% Datum: 2019-12-05 18:34
% Revision: 77da58f92bca3eff71542919beafa37024070d86 (2019-12-05)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Mq = S5RRPRR5_inertiaJ_slag_vp1(qJ, ...
  pkin, m, rSges, Icges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(9,1),zeros(6,1),zeros(6,3),zeros(6,6)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RRPRR5_inertiaJ_slag_vp1: qJ has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S5RRPRR5_inertiaJ_slag_vp1: pkin has to be [9x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5RRPRR5_inertiaJ_slag_vp1: m has to be [6x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [6,3]), ...
  'S5RRPRR5_inertiaJ_slag_vp1: rSges has to be [6x3] (double)');
assert(isreal(Icges) && all(size(Icges) == [6 6]), ...
  'S5RRPRR5_inertiaJ_slag_vp1: Icges has to be [6x6] (double)'); 

%% Symbolic Calculation
% From inertia_joint_joint_floatb_twist_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-05 18:33:37
% EndTime: 2019-12-05 18:33:39
% DurationCPUTime: 0.66s
% Computational Cost: add. (2367->173), mult. (1556->255), div. (0->0), fcn. (1342->10), ass. (0->97)
t92 = qJ(1) + qJ(2);
t87 = sin(t92);
t88 = cos(t92);
t139 = t87 * t88;
t141 = rSges(4,3) + qJ(3);
t93 = sin(pkin(9));
t94 = cos(pkin(9));
t140 = -rSges(4,1) * t94 + rSges(4,2) * t93 - pkin(2);
t82 = t87 ^ 2;
t83 = t88 ^ 2;
t138 = t87 / 0.2e1;
t137 = t88 / 0.2e1;
t91 = pkin(9) + qJ(4);
t86 = qJ(5) + t91;
t79 = sin(t86);
t80 = cos(t86);
t100 = Icges(6,5) * t80 - Icges(6,6) * t79;
t29 = Icges(6,3) * t88 - t100 * t87;
t30 = Icges(6,3) * t87 + t100 * t88;
t136 = t87 * (t29 * t139 + t82 * t30) + t88 * (t30 * t139 + t83 * t29);
t84 = sin(t91);
t85 = cos(t91);
t59 = t84 * rSges(5,1) + t85 * rSges(5,2);
t135 = m(5) * t59;
t52 = t79 * rSges(6,1) + t80 * rSges(6,2);
t134 = m(6) * t52;
t96 = sin(qJ(1));
t133 = t96 * pkin(1);
t97 = cos(qJ(1));
t132 = t97 * pkin(1);
t81 = t94 * pkin(3) + pkin(2);
t131 = rSges(5,1) * t85;
t130 = rSges(6,1) * t80;
t128 = rSges(5,2) * t84;
t127 = rSges(6,2) * t79;
t95 = -pkin(7) - qJ(3);
t62 = pkin(4) * t85 + t81;
t126 = t62 - t81;
t125 = t88 * rSges(6,3) + t87 * t127;
t124 = t88 * rSges(5,3) + t87 * t128;
t123 = t83 + t82;
t122 = Icges(5,4) * t84;
t121 = Icges(5,4) * t85;
t120 = Icges(6,4) * t79;
t119 = Icges(6,4) * t80;
t102 = -Icges(6,2) * t79 + t119;
t104 = Icges(6,1) * t80 - t120;
t50 = Icges(6,2) * t80 + t120;
t51 = Icges(6,1) * t79 + t119;
t107 = t50 * t79 - t51 * t80;
t49 = Icges(6,5) * t79 + Icges(6,6) * t80;
t118 = (-t107 * t88 + t80 * (Icges(6,6) * t87 + t102 * t88) + t79 * (Icges(6,5) * t87 + t104 * t88) + t87 * t49) * t138 + (t107 * t87 + t80 * (Icges(6,6) * t88 - t102 * t87) + t79 * (Icges(6,5) * t88 - t104 * t87) + t88 * t49) * t137;
t116 = pkin(4) * t84 + t52;
t61 = -t88 * rSges(3,1) + t87 * rSges(3,2);
t115 = -t81 - t131;
t114 = -t62 - t130;
t113 = -t87 * rSges(5,3) + t88 * t128;
t112 = -t87 * rSges(6,3) + t88 * t127;
t60 = -t87 * rSges(3,1) - t88 * rSges(3,2);
t55 = Icges(5,2) * t85 + t122;
t56 = Icges(5,1) * t84 + t121;
t106 = t55 * t84 - t56 * t85;
t105 = Icges(5,1) * t85 - t122;
t103 = -Icges(5,2) * t84 + t121;
t101 = Icges(5,5) * t85 - Icges(5,6) * t84;
t54 = Icges(5,5) * t84 + Icges(5,6) * t85;
t99 = t118 + (-t106 * t88 + t85 * (Icges(5,6) * t87 + t103 * t88) + t84 * (Icges(5,5) * t87 + t105 * t88) + t87 * t54) * t138 + (t106 * t87 + t85 * (Icges(5,6) * t88 - t103 * t87) + t84 * (Icges(5,5) * t88 - t105 * t87) + t88 * t54) * t137;
t98 = Icges(4,2) * t94 ^ 2 + t80 * t50 + t79 * t51 + t85 * t55 + t84 * t56 + Icges(3,3) + (Icges(4,1) * t93 + 0.2e1 * Icges(4,4) * t94) * t93;
t24 = t140 * t87 + t141 * t88;
t73 = t87 * t95;
t21 = t115 * t88 + t113 + t73;
t90 = -pkin(8) + t95;
t72 = t87 * t90;
t17 = t114 * t88 + t112 + t72;
t16 = t114 * t87 - t88 * t90 + t125;
t20 = t115 * t87 - t88 * t95 + t124;
t25 = t140 * t88 - t141 * t87;
t69 = -t97 * rSges(2,1) + t96 * rSges(2,2);
t68 = -t96 * rSges(2,1) - t97 * rSges(2,2);
t47 = t61 - t132;
t46 = t60 - t133;
t37 = Icges(5,3) * t87 + t101 * t88;
t36 = Icges(5,3) * t88 - t101 * t87;
t35 = -t130 * t87 + t125;
t28 = t116 * t88;
t27 = t116 * t87;
t26 = t88 * (t130 * t88 - t112);
t23 = t25 - t132;
t22 = t24 - t133;
t19 = t21 - t132;
t18 = t20 - t133;
t15 = t17 - t132;
t14 = t16 - t133;
t13 = t88 * (t131 * t88 - t113) - t87 * (-t131 * t87 + t124);
t12 = -t87 * t35 + t26;
t3 = t26 + (t126 * t88 - t72 + t73) * t88 + (-t35 + t126 * t87 + (t90 - t95) * t88) * t87;
t1 = [Icges(2,3) + m(2) * (t68 ^ 2 + t69 ^ 2) + m(3) * (t46 ^ 2 + t47 ^ 2) + m(4) * (t22 ^ 2 + t23 ^ 2) + m(5) * (t18 ^ 2 + t19 ^ 2) + m(6) * (t14 ^ 2 + t15 ^ 2) + t98; m(3) * (t60 * t46 + t61 * t47) + m(4) * (t24 * t22 + t25 * t23) + m(5) * (t20 * t18 + t21 * t19) + m(6) * (t16 * t14 + t17 * t15) + t98; m(6) * (t16 ^ 2 + t17 ^ 2) + m(5) * (t20 ^ 2 + t21 ^ 2) + m(4) * (t24 ^ 2 + t25 ^ 2) + m(3) * (t60 ^ 2 + t61 ^ 2) + t98; m(4) * (t87 * t22 + t88 * t23) + m(5) * (t87 * t18 + t88 * t19) + m(6) * (t87 * t14 + t88 * t15); m(6) * (t87 * t16 + t88 * t17) + m(5) * (t87 * t20 + t88 * t21) + m(4) * (t87 * t24 + t88 * t25); 0.2e1 * (m(4) / 0.2e1 + m(5) / 0.2e1 + m(6) / 0.2e1) * t123; m(6) * (-t28 * t14 + t27 * t15) + (-t18 * t88 + t19 * t87) * t135 + t99; m(6) * (-t28 * t16 + t27 * t17) + (-t20 * t88 + t21 * t87) * t135 + t99; m(6) * (t27 * t88 - t28 * t87); m(5) * (t123 * t59 ^ 2 + t13 ^ 2) + t88 * (t37 * t139 + t83 * t36) + t87 * (t36 * t139 + t82 * t37) + m(6) * (t27 ^ 2 + t28 ^ 2 + t3 ^ 2) + t136; (-t14 * t88 + t15 * t87) * t134 + t118; (-t16 * t88 + t17 * t87) * t134 + t118; 0; m(6) * (t12 * t3 + (t27 * t87 + t28 * t88) * t52) + t136; m(6) * (t123 * t52 ^ 2 + t12 ^ 2) + t136;];
%% Postprocessing: Reshape Output
% From vec2symmat_5_matlab.m
res = [t1(1), t1(2), t1(4), t1(7), t1(11); t1(2), t1(3), t1(5), t1(8), t1(12); t1(4), t1(5), t1(6), t1(9), t1(13); t1(7), t1(8), t1(9), t1(10), t1(14); t1(11), t1(12), t1(13), t1(14), t1(15);];
Mq = res;
