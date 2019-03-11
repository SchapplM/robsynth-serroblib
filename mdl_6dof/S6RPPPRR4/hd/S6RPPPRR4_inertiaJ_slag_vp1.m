% Calculate joint inertia matrix for
% S6RPPPRR4
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% pkin [9x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d5,d6,theta3]';
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
% Datum: 2019-03-09 01:36
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Mq = S6RPPPRR4_inertiaJ_slag_vp1(qJ, ...
  pkin, m, rSges, Icges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(9,1),zeros(7,1),zeros(7,3),zeros(7,6)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPPPRR4_inertiaJ_slag_vp1: qJ has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S6RPPPRR4_inertiaJ_slag_vp1: pkin has to be [9x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6RPPPRR4_inertiaJ_slag_vp1: m has to be [7x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [7,3]), ...
  'S6RPPPRR4_inertiaJ_slag_vp1: rSges has to be [7x3] (double)');
assert(isreal(Icges) && all(size(Icges) == [7 6]), ...
  'S6RPPPRR4_inertiaJ_slag_vp1: Icges has to be [7x6] (double)'); 

%% Symbolic Calculation
% From inertia_joint_joint_floatb_twist_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 01:35:23
% EndTime: 2019-03-09 01:35:24
% DurationCPUTime: 1.18s
% Computational Cost: add. (2077->204), mult. (4487->317), div. (0->0), fcn. (5518->8), ass. (0->102)
t92 = cos(qJ(5));
t139 = Icges(6,5) * t92;
t89 = sin(qJ(5));
t134 = -t89 / 0.2e1;
t138 = -t139 / 0.2e1 - Icges(6,6) * t134;
t115 = sin(pkin(9));
t116 = cos(pkin(9));
t90 = sin(qJ(1));
t93 = cos(qJ(1));
t69 = t93 * t115 - t90 * t116;
t136 = t69 ^ 2;
t68 = -t90 * t115 - t93 * t116;
t67 = t68 ^ 2;
t114 = t67 + t136;
t107 = m(5) / 0.2e1 + m(6) / 0.2e1 + m(7) / 0.2e1;
t137 = 0.2e1 * t107;
t135 = -t69 / 0.2e1;
t77 = -rSges(6,1) * t92 + rSges(6,2) * t89;
t132 = m(6) * t77;
t131 = pkin(5) * t89;
t88 = sin(qJ(6));
t91 = cos(qJ(6));
t56 = -Icges(7,5) * t89 + (-Icges(7,1) * t91 + Icges(7,4) * t88) * t92;
t130 = t56 * t91;
t129 = t68 * t92;
t128 = t69 * t89;
t127 = t69 * t92;
t126 = t88 * t89;
t125 = t89 * t91;
t47 = -t69 * t126 - t68 * t91;
t48 = t69 * t125 - t68 * t88;
t124 = t48 * rSges(7,1) + t47 * rSges(7,2);
t57 = -t89 * rSges(7,3) + (-rSges(7,1) * t91 + rSges(7,2) * t88) * t92;
t123 = -pkin(5) * t92 - pkin(8) * t89 + t57;
t122 = t93 * pkin(1) + t90 * qJ(2);
t119 = Icges(7,5) * t92;
t118 = Icges(7,6) * t92;
t117 = Icges(7,3) * t92;
t113 = rSges(6,1) * t128 + rSges(6,2) * t127 - t68 * rSges(6,3);
t112 = t93 * pkin(2) + t122;
t24 = Icges(7,5) * t48 + Icges(7,6) * t47 - t69 * t117;
t26 = Icges(7,4) * t48 + Icges(7,2) * t47 - t69 * t118;
t28 = Icges(7,1) * t48 + Icges(7,4) * t47 - t69 * t119;
t10 = -t89 * t24 + (t26 * t88 - t28 * t91) * t92;
t54 = -Icges(7,3) * t89 + (-Icges(7,5) * t91 + Icges(7,6) * t88) * t92;
t55 = -Icges(7,6) * t89 + (-Icges(7,4) * t91 + Icges(7,2) * t88) * t92;
t16 = -t54 * t127 + t47 * t55 + t48 * t56;
t111 = -t10 / 0.2e1 - t16 / 0.2e1;
t49 = t68 * t126 - t69 * t91;
t50 = -t68 * t125 - t69 * t88;
t25 = Icges(7,5) * t50 + Icges(7,6) * t49 + t68 * t117;
t27 = Icges(7,4) * t50 + Icges(7,2) * t49 + t68 * t118;
t29 = Icges(7,1) * t50 + Icges(7,4) * t49 + t68 * t119;
t11 = -t89 * t25 + (t27 * t88 - t29 * t91) * t92;
t17 = t54 * t129 + t49 * t55 + t50 * t56;
t110 = t11 / 0.2e1 + t17 / 0.2e1;
t109 = (-rSges(7,3) - pkin(8)) * t92;
t108 = -t68 * pkin(3) + t112;
t85 = t93 * qJ(2);
t106 = t85 + (-pkin(1) - pkin(2)) * t90;
t105 = -t50 * rSges(7,1) - t49 * rSges(7,2);
t97 = t68 * qJ(4) + t106;
t94 = t69 * pkin(7) + t97;
t95 = (rSges(6,1) * t89 + rSges(6,2) * t92) * t68;
t22 = (rSges(6,3) + pkin(3)) * t69 + t95 + t94;
t96 = -t68 * pkin(7) + t108;
t23 = t69 * qJ(4) + t113 + t96;
t104 = t22 * t69 - t23 * t68;
t101 = t68 * t93 + t69 * t90;
t98 = Icges(6,5) * t89 + Icges(6,6) * t92;
t78 = rSges(2,1) * t93 - t90 * rSges(2,2);
t76 = -t90 * rSges(2,1) - rSges(2,2) * t93;
t62 = pkin(5) * t128;
t59 = rSges(3,1) * t93 + t90 * rSges(3,3) + t122;
t58 = t93 * rSges(3,3) + t85 + (-rSges(3,1) - pkin(1)) * t90;
t53 = t92 * t88 * t55;
t44 = -rSges(4,1) * t68 - rSges(4,2) * t69 + t112;
t43 = t69 * rSges(4,1) - t68 * rSges(4,2) + t106;
t38 = -Icges(6,3) * t69 - t98 * t68;
t36 = t123 * t68;
t35 = t123 * t69;
t34 = t68 * rSges(5,2) + (rSges(5,3) + qJ(4)) * t69 + t108;
t33 = t68 * rSges(5,3) + (-rSges(5,2) + pkin(3)) * t69 + t97;
t32 = -t92 * t130 - t89 * t54 + t53;
t31 = rSges(7,3) * t129 - t105;
t30 = -rSges(7,3) * t127 + t124;
t20 = t69 * t113 + (t69 * rSges(6,3) + t95) * t68;
t19 = t57 * t129 + t31 * t89;
t18 = t57 * t127 - t30 * t89;
t15 = t62 + (qJ(4) + t109) * t69 + t96 + t124;
t14 = t69 * pkin(3) + (t109 + t131) * t68 + t94 + t105;
t12 = t31 * t127 + t30 * t129;
t9 = (-pkin(8) * t127 + t30 + t62) * t69 + (-t31 + (-pkin(8) * t92 + t131) * t68) * t68;
t8 = t25 * t129 + t27 * t49 + t29 * t50;
t7 = t24 * t129 + t26 * t49 + t28 * t50;
t6 = -t25 * t127 + t27 * t47 + t29 * t48;
t5 = -t24 * t127 + t26 * t47 + t28 * t48;
t4 = -t68 * t7 - t69 * t8;
t3 = -t5 * t68 - t6 * t69;
t2 = -t17 * t89 + (t68 * t8 - t69 * t7) * t92;
t1 = -t16 * t89 + (-t5 * t69 + t6 * t68) * t92;
t13 = [Icges(5,1) + Icges(3,2) + Icges(2,3) + Icges(4,3) + t53 + (Icges(6,1) * t92 - t130) * t92 + m(7) * (t14 ^ 2 + t15 ^ 2) + m(6) * (t22 ^ 2 + t23 ^ 2) + m(5) * (t33 ^ 2 + t34 ^ 2) + m(3) * (t58 ^ 2 + t59 ^ 2) + m(4) * (t43 ^ 2 + t44 ^ 2) + m(2) * (t76 ^ 2 + t78 ^ 2) + (-0.2e1 * Icges(6,4) * t92 + Icges(6,2) * t89 - t54) * t89; m(7) * (t90 * t14 - t15 * t93) + m(6) * (t90 * t22 - t23 * t93) + m(5) * (t90 * t33 - t34 * t93) + m(3) * (t90 * t58 - t59 * t93) + m(4) * (t90 * t43 - t44 * t93); 0.2e1 * (m(3) / 0.2e1 + m(4) / 0.2e1 + t107) * (t90 ^ 2 + t93 ^ 2); 0; 0; m(4) + m(5) + m(6) + m(7); m(7) * (t14 * t69 - t15 * t68) + m(6) * t104 + m(5) * (t33 * t69 - t34 * t68); t101 * t137; 0; t114 * t137; m(7) * (-t14 * t35 + t15 * t36) - t104 * t132 + (t136 / 0.2e1 + t67 / 0.2e1) * (Icges(6,6) * t89 - t139) + (t138 * t69 - t110) * t69 + (t138 * t68 + t111) * t68; m(7) * (-t35 * t90 - t36 * t93) - t101 * t132; -m(6) * t20 - m(7) * t9; m(7) * (-t35 * t69 - t36 * t68) - t114 * t132; m(6) * (t114 * t77 ^ 2 + t20 ^ 2) + m(7) * (t35 ^ 2 + t36 ^ 2 + t9 ^ 2) + (-t136 * t38 - t4) * t69 + (-t69 * t38 * t68 - t3 - t114 * (-Icges(6,3) * t68 + t98 * t69)) * t68; -t32 * t89 + m(7) * (t14 * t19 + t15 * t18) + (t110 * t68 + t111 * t69) * t92; m(7) * (-t18 * t93 + t19 * t90); m(7) * t12; m(7) * (-t18 * t68 + t19 * t69); -t68 * t1 / 0.2e1 + t2 * t135 + (-t10 * t68 - t11 * t69) * t134 + m(7) * (-t12 * t9 + t18 * t36 - t19 * t35) + (t3 * t135 + t68 * t4 / 0.2e1) * t92; t89 ^ 2 * t32 + m(7) * (t12 ^ 2 + t18 ^ 2 + t19 ^ 2) + (-t69 * t1 + t68 * t2 - t89 * (-t10 * t69 + t11 * t68)) * t92;];
%% Postprocessing: Reshape Output
% From vec2symmat_6_matlab.m
res = [t13(1) t13(2) t13(4) t13(7) t13(11) t13(16); t13(2) t13(3) t13(5) t13(8) t13(12) t13(17); t13(4) t13(5) t13(6) t13(9) t13(13) t13(18); t13(7) t13(8) t13(9) t13(10) t13(14) t13(19); t13(11) t13(12) t13(13) t13(14) t13(15) t13(20); t13(16) t13(17) t13(18) t13(19) t13(20) t13(21);];
Mq  = res;
