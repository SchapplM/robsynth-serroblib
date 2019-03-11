% Calculate joint inertia matrix for
% S6RPPPRR2
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% pkin [10x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d5,d6,theta2,theta4]';
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
% Datum: 2019-03-09 01:32
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Mq = S6RPPPRR2_inertiaJ_slag_vp1(qJ, ...
  pkin, m, rSges, Icges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(10,1),zeros(7,1),zeros(7,3),zeros(7,6)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPPPRR2_inertiaJ_slag_vp1: qJ has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S6RPPPRR2_inertiaJ_slag_vp1: pkin has to be [10x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6RPPPRR2_inertiaJ_slag_vp1: m has to be [7x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [7,3]), ...
  'S6RPPPRR2_inertiaJ_slag_vp1: rSges has to be [7x3] (double)');
assert(isreal(Icges) && all(size(Icges) == [7 6]), ...
  'S6RPPPRR2_inertiaJ_slag_vp1: Icges has to be [7x6] (double)'); 

%% Symbolic Calculation
% From inertia_joint_joint_floatb_twist_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 01:31:39
% EndTime: 2019-03-09 01:31:41
% DurationCPUTime: 0.94s
% Computational Cost: add. (2895->198), mult. (2475->304), div. (0->0), fcn. (2534->10), ass. (0->107)
t86 = pkin(10) + qJ(5);
t83 = cos(t86);
t140 = Icges(6,5) * t83;
t81 = sin(t86);
t139 = Icges(6,6) * t81;
t138 = t140 / 0.2e1 - t139 / 0.2e1;
t137 = m(5) + m(7);
t87 = qJ(1) + pkin(9);
t82 = sin(t87);
t79 = t82 ^ 2;
t84 = cos(t87);
t80 = t84 ^ 2;
t66 = t79 + t80;
t136 = (rSges(6,1) * t81 + rSges(6,2) * t83) * t84;
t133 = t84 / 0.2e1;
t58 = m(6) * t66;
t88 = sin(pkin(10));
t132 = pkin(4) * t88;
t131 = pkin(5) * t81;
t92 = sin(qJ(1));
t130 = t92 * pkin(1);
t129 = t81 * t82;
t128 = t82 * t83;
t91 = sin(qJ(6));
t127 = t82 * t91;
t93 = cos(qJ(6));
t126 = t82 * t93;
t125 = t83 * t84;
t124 = t84 * t91;
t123 = t84 * t93;
t46 = Icges(7,6) * t81 + (Icges(7,4) * t93 - Icges(7,2) * t91) * t83;
t122 = t91 * t46;
t45 = Icges(7,3) * t81 + (Icges(7,5) * t93 - Icges(7,6) * t91) * t83;
t47 = Icges(7,5) * t81 + (Icges(7,1) * t93 - Icges(7,4) * t91) * t83;
t121 = t83 * t93 * t47 + t81 * t45;
t48 = t81 * rSges(7,3) + (rSges(7,1) * t93 - rSges(7,2) * t91) * t83;
t120 = t83 * pkin(5) + t81 * pkin(8) + t48;
t51 = -t81 * t127 + t123;
t52 = t81 * t126 + t124;
t119 = t52 * rSges(7,1) + t51 * rSges(7,2);
t116 = Icges(7,5) * t83;
t115 = Icges(7,6) * t83;
t114 = Icges(7,3) * t83;
t113 = rSges(5,3) + qJ(4);
t112 = t137 * t66 + t58;
t111 = rSges(6,1) * t129 + rSges(6,2) * t128 + t84 * rSges(6,3);
t94 = cos(qJ(1));
t85 = t94 * pkin(1);
t110 = t84 * pkin(2) + t82 * qJ(3) + t85;
t13 = -t45 * t128 + t51 * t46 + t52 * t47;
t23 = Icges(7,5) * t52 + Icges(7,6) * t51 - t82 * t114;
t25 = Icges(7,4) * t52 + Icges(7,2) * t51 - t82 * t115;
t27 = Icges(7,1) * t52 + Icges(7,4) * t51 - t82 * t116;
t9 = t81 * t23 + (-t25 * t91 + t27 * t93) * t83;
t109 = -t13 / 0.2e1 - t9 / 0.2e1;
t53 = t81 * t124 + t126;
t54 = -t81 * t123 + t127;
t24 = Icges(7,5) * t54 + Icges(7,6) * t53 + t84 * t114;
t26 = Icges(7,4) * t54 + Icges(7,2) * t53 + t84 * t115;
t28 = Icges(7,1) * t54 + Icges(7,4) * t53 + t84 * t116;
t10 = t81 * t24 + (-t26 * t91 + t28 * t93) * t83;
t14 = t45 * t125 + t53 * t46 + t54 * t47;
t108 = t14 / 0.2e1 + t10 / 0.2e1;
t107 = (-rSges(7,3) - pkin(8)) * t83;
t106 = t84 * qJ(3) - t130;
t89 = cos(pkin(10));
t105 = rSges(5,1) * t88 + rSges(5,2) * t89;
t103 = -t54 * rSges(7,1) - t53 * rSges(7,2);
t90 = -pkin(7) - qJ(4);
t100 = t84 * t132 + t82 * t90 + t106;
t97 = Icges(6,5) * t81 + Icges(6,6) * t83;
t96 = t82 * t132 - t84 * t90 + t110;
t21 = t136 + (-rSges(6,3) - pkin(2)) * t82 + t100;
t22 = t96 + t111;
t95 = m(6) * (t82 * t21 - t84 * t22);
t71 = t94 * rSges(2,1) - t92 * rSges(2,2);
t70 = -t92 * rSges(2,1) - t94 * rSges(2,2);
t69 = pkin(5) * t129;
t64 = t83 * rSges(6,1) - t81 * rSges(6,2);
t56 = t84 * rSges(3,1) - t82 * rSges(3,2) + t85;
t55 = -t82 * rSges(3,1) - t84 * rSges(3,2) - t130;
t39 = Icges(6,3) * t84 + t97 * t82;
t38 = -t84 * rSges(4,2) + t82 * rSges(4,3) + t110;
t37 = t84 * rSges(4,3) + (rSges(4,2) - pkin(2)) * t82 + t106;
t34 = t105 * t82 + t113 * t84 + t110;
t33 = t105 * t84 + (-pkin(2) - t113) * t82 + t106;
t32 = t120 * t84;
t31 = t120 * t82;
t30 = rSges(7,3) * t125 - t103;
t29 = -rSges(7,3) * t128 + t119;
t20 = -t82 * t111 + (t82 * rSges(6,3) - t136) * t84;
t19 = t48 * t125 - t81 * t30;
t18 = t48 * t128 + t81 * t29;
t17 = (-t83 * t122 + t121) * t81;
t16 = t82 * t107 + t119 + t69 + t96;
t15 = -t82 * pkin(2) + (t107 + t131) * t84 + t100 + t103;
t12 = (-t29 * t84 - t30 * t82) * t83;
t11 = (pkin(8) * t128 - t29 - t69) * t82 + (t30 + (pkin(8) * t83 - t131) * t84) * t84;
t8 = t24 * t125 + t53 * t26 + t54 * t28;
t7 = t23 * t125 + t53 * t25 + t54 * t27;
t6 = -t24 * t128 + t51 * t26 + t52 * t28;
t5 = -t23 * t128 + t51 * t25 + t52 * t27;
t4 = t7 * t84 + t8 * t82;
t3 = t5 * t84 + t6 * t82;
t2 = t14 * t81 + (-t7 * t82 + t8 * t84) * t83;
t1 = t13 * t81 + (-t5 * t82 + t6 * t84) * t83;
t35 = [Icges(5,1) * t89 ^ 2 + Icges(4,1) + Icges(2,3) + Icges(3,3) + (-0.2e1 * Icges(5,4) * t89 + Icges(5,2) * t88) * t88 + (Icges(6,1) * t83 - t122) * t83 + m(7) * (t15 ^ 2 + t16 ^ 2) + m(6) * (t21 ^ 2 + t22 ^ 2) + m(4) * (t37 ^ 2 + t38 ^ 2) + m(5) * (t33 ^ 2 + t34 ^ 2) + m(3) * (t55 ^ 2 + t56 ^ 2) + m(2) * (t70 ^ 2 + t71 ^ 2) + t121 + (-0.2e1 * Icges(6,4) * t83 + Icges(6,2) * t81) * t81; 0; m(3) + m(4) + m(6) + t137; m(7) * (t82 * t15 - t84 * t16) + t95 + m(4) * (t82 * t37 - t84 * t38) + m(5) * (t82 * t33 - t84 * t34); 0; m(4) * t66 + t112; m(7) * (t84 * t15 + t82 * t16) + m(6) * (t84 * t21 + t82 * t22) + m(5) * (t84 * t33 + t82 * t34); 0; 0; t112; m(7) * (t31 * t15 - t32 * t16) + t64 * t95 + (t79 / 0.2e1 + t80 / 0.2e1) * (-t139 + t140) + (t138 * t84 - t109) * t84 + (t138 * t82 + t108) * t82; m(6) * t20 + m(7) * t11; m(7) * (t31 * t82 + t32 * t84) + t64 * t58; m(7) * (t31 * t84 - t32 * t82); m(6) * (t66 * t64 ^ 2 + t20 ^ 2) + m(7) * (t11 ^ 2 + t31 ^ 2 + t32 ^ 2) + (t80 * t39 + t3) * t84 + (t82 * t39 * t84 + t4 + t66 * (Icges(6,3) * t82 - t97 * t84)) * t82; m(7) * (t19 * t15 + t18 * t16) + t17 + (t108 * t84 + t109 * t82) * t83; m(7) * t12; m(7) * (-t18 * t84 + t19 * t82); m(7) * (t18 * t82 + t19 * t84); t1 * t133 + t82 * t2 / 0.2e1 + t81 * (t10 * t82 + t9 * t84) / 0.2e1 + m(7) * (t12 * t11 - t18 * t32 + t19 * t31) + (-t82 * t3 / 0.2e1 + t4 * t133) * t83; t81 * t17 + m(7) * (t12 ^ 2 + t18 ^ 2 + t19 ^ 2) + (-t82 * t1 + t84 * t2 + t81 * (t10 * t84 - t82 * t9)) * t83;];
%% Postprocessing: Reshape Output
% From vec2symmat_6_matlab.m
res = [t35(1) t35(2) t35(4) t35(7) t35(11) t35(16); t35(2) t35(3) t35(5) t35(8) t35(12) t35(17); t35(4) t35(5) t35(6) t35(9) t35(13) t35(18); t35(7) t35(8) t35(9) t35(10) t35(14) t35(19); t35(11) t35(12) t35(13) t35(14) t35(15) t35(20); t35(16) t35(17) t35(18) t35(19) t35(20) t35(21);];
Mq  = res;
