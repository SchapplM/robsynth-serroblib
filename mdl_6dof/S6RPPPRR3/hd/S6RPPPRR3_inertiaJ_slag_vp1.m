% Calculate joint inertia matrix for
% S6RPPPRR3
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% pkin [10x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d5,d6,theta3,theta4]';
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
% Datum: 2019-03-09 01:34
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Mq = S6RPPPRR3_inertiaJ_slag_vp1(qJ, ...
  pkin, m, rSges, Icges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(10,1),zeros(7,1),zeros(7,3),zeros(7,6)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPPPRR3_inertiaJ_slag_vp1: qJ has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S6RPPPRR3_inertiaJ_slag_vp1: pkin has to be [10x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6RPPPRR3_inertiaJ_slag_vp1: m has to be [7x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [7,3]), ...
  'S6RPPPRR3_inertiaJ_slag_vp1: rSges has to be [7x3] (double)');
assert(isreal(Icges) && all(size(Icges) == [7 6]), ...
  'S6RPPPRR3_inertiaJ_slag_vp1: Icges has to be [7x6] (double)'); 

%% Symbolic Calculation
% From inertia_joint_joint_floatb_twist_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 01:33:28
% EndTime: 2019-03-09 01:33:31
% DurationCPUTime: 1.11s
% Computational Cost: add. (2867->212), mult. (4469->334), div. (0->0), fcn. (5490->10), ass. (0->104)
t90 = pkin(10) + qJ(5);
t84 = sin(t90);
t139 = Icges(6,5) * t84;
t138 = -t139 / 0.2e1;
t113 = sin(pkin(9));
t114 = cos(pkin(9));
t95 = sin(qJ(1));
t97 = cos(qJ(1));
t75 = t97 * t113 - t95 * t114;
t71 = t75 ^ 2;
t74 = -t95 * t113 - t97 * t114;
t72 = t74 ^ 2;
t122 = t71 + t72;
t109 = m(5) / 0.2e1 + m(6) / 0.2e1 + m(7) / 0.2e1;
t137 = 0.2e1 * t109;
t85 = cos(t90);
t129 = t74 * t85;
t130 = t74 * t84;
t94 = sin(qJ(6));
t127 = t85 * t94;
t96 = cos(qJ(6));
t47 = t74 * t127 + t75 * t96;
t126 = t85 * t96;
t48 = -t74 * t126 + t75 * t94;
t30 = t48 * rSges(7,1) + t47 * rSges(7,2) - rSges(7,3) * t130;
t136 = -pkin(5) * t129 - pkin(8) * t130 + t30;
t135 = -t74 / 0.2e1;
t134 = t85 / 0.2e1;
t70 = -rSges(6,1) * t84 - rSges(6,2) * t85;
t133 = m(6) * t70;
t132 = pkin(5) * t85;
t55 = Icges(7,5) * t85 + (-Icges(7,1) * t96 + Icges(7,4) * t94) * t84;
t131 = t55 * t96;
t128 = t75 * t84;
t53 = Icges(7,3) * t85 + (-Icges(7,5) * t96 + Icges(7,6) * t94) * t84;
t54 = Icges(7,6) * t85 + (-Icges(7,4) * t96 + Icges(7,2) * t94) * t84;
t125 = t84 * t94 * t54 + t85 * t53;
t56 = t85 * rSges(7,3) + (-rSges(7,1) * t96 + rSges(7,2) * t94) * t84;
t124 = pkin(5) * t84 - pkin(8) * t85 - t56;
t121 = t97 * pkin(1) + t95 * qJ(2);
t119 = Icges(6,4) * t85;
t118 = Icges(7,5) * t84;
t117 = Icges(7,6) * t84;
t116 = Icges(7,3) * t84;
t115 = rSges(5,3) + qJ(4);
t112 = t97 * pkin(2) + t121;
t45 = t75 * t127 - t74 * t96;
t46 = -t75 * t126 - t74 * t94;
t23 = Icges(7,5) * t46 + Icges(7,6) * t45 - t75 * t116;
t25 = Icges(7,4) * t46 + Icges(7,2) * t45 - t75 * t117;
t27 = Icges(7,1) * t46 + Icges(7,4) * t45 - t75 * t118;
t10 = t85 * t23 + (t25 * t94 - t27 * t96) * t84;
t15 = -t53 * t128 + t45 * t54 + t46 * t55;
t111 = -t10 / 0.2e1 - t15 / 0.2e1;
t24 = Icges(7,5) * t48 + Icges(7,6) * t47 - t74 * t116;
t26 = Icges(7,4) * t48 + Icges(7,2) * t47 - t74 * t117;
t28 = Icges(7,1) * t48 + Icges(7,4) * t47 - t74 * t118;
t11 = t85 * t24 + (t26 * t94 - t28 * t96) * t84;
t16 = -t53 * t130 + t47 * t54 + t48 * t55;
t110 = -t11 / 0.2e1 - t16 / 0.2e1;
t87 = t97 * qJ(2);
t108 = t87 + (-pkin(1) - pkin(2)) * t95;
t92 = cos(pkin(10));
t83 = pkin(4) * t92 + pkin(3);
t93 = -pkin(7) - qJ(4);
t107 = -t74 * t83 - t75 * t93 + t112;
t106 = -rSges(6,1) * t85 + rSges(6,2) * t84;
t105 = -rSges(7,1) * t46 - rSges(7,2) * t45;
t101 = Icges(6,2) * t84 - t119;
t100 = -Icges(6,5) * t85 + Icges(6,6) * t84;
t99 = -rSges(6,1) * t129 + rSges(6,2) * t130 + t75 * rSges(6,3);
t91 = sin(pkin(10));
t98 = rSges(5,1) * t92 - rSges(5,2) * t91 + pkin(3);
t78 = rSges(2,1) * t97 - t95 * rSges(2,2);
t77 = -t95 * rSges(2,1) - rSges(2,2) * t97;
t63 = rSges(3,1) * t97 + t95 * rSges(3,3) + t121;
t62 = rSges(3,3) * t97 + t87 + (-rSges(3,1) - pkin(1)) * t95;
t42 = -rSges(4,1) * t74 - rSges(4,2) * t75 + t112;
t41 = rSges(4,1) * t75 - rSges(4,2) * t74 + t108;
t36 = Icges(6,3) * t75 + t100 * t74;
t34 = t124 * t74;
t33 = t124 * t75;
t32 = t115 * t75 - t98 * t74 + t112;
t31 = t115 * t74 + t98 * t75 + t108;
t29 = -rSges(7,3) * t128 - t105;
t22 = t99 + t107;
t21 = (rSges(6,3) - t93) * t74 + (-t106 + t83) * t75 + t108;
t20 = (-t84 * t131 + t125) * t85;
t19 = t74 * t99 + (-t74 * rSges(6,3) + t106 * t75) * t75;
t18 = t56 * t130 + t30 * t85;
t17 = -t56 * t128 - t29 * t85;
t14 = t107 + t136;
t13 = -t74 * t93 + (t132 + t83 + (rSges(7,3) + pkin(8)) * t84) * t75 + t105 + t108;
t12 = (-t29 * t74 + t30 * t75) * t84;
t9 = t136 * t74 + (t29 + (-pkin(8) * t84 - t132) * t75) * t75;
t8 = -t24 * t130 + t26 * t47 + t28 * t48;
t7 = -t23 * t130 + t25 * t47 + t27 * t48;
t6 = -t24 * t128 + t26 * t45 + t28 * t46;
t5 = -t23 * t128 + t25 * t45 + t27 * t46;
t4 = -t7 * t74 + t75 * t8;
t3 = -t5 * t74 + t6 * t75;
t2 = t16 * t85 + (-t7 * t75 - t74 * t8) * t84;
t1 = t15 * t85 + (-t5 * t75 - t6 * t74) * t84;
t35 = [t92 ^ 2 * Icges(5,2) - t85 * (-Icges(6,4) * t84 - Icges(6,2) * t85) + Icges(3,2) + Icges(2,3) + Icges(4,3) + (Icges(5,1) * t91 + 0.2e1 * Icges(5,4) * t92) * t91 + (Icges(6,1) * t84 + t119 - t131) * t84 + m(7) * (t13 ^ 2 + t14 ^ 2) + m(6) * (t21 ^ 2 + t22 ^ 2) + m(4) * (t41 ^ 2 + t42 ^ 2) + m(5) * (t31 ^ 2 + t32 ^ 2) + m(3) * (t62 ^ 2 + t63 ^ 2) + m(2) * (t77 ^ 2 + t78 ^ 2) + t125; m(7) * (t95 * t13 - t14 * t97) + m(6) * (t95 * t21 - t22 * t97) + m(4) * (t95 * t41 - t42 * t97) + m(5) * (t95 * t31 - t32 * t97) + m(3) * (t95 * t62 - t63 * t97); 0.2e1 * (m(3) / 0.2e1 + m(4) / 0.2e1 + t109) * (t95 ^ 2 + t97 ^ 2); 0; 0; m(4) + m(5) + m(6) + m(7); m(7) * (t13 * t75 - t14 * t74) + m(6) * (t21 * t75 - t22 * t74) + m(5) * (t31 * t75 - t32 * t74); (t74 * t97 + t75 * t95) * t137; 0; t122 * t137; (-t85 * (Icges(6,6) * t75 + t101 * t74) / 0.2e1 + t75 * t138 - t110) * t75 + ((-Icges(6,6) * t74 + t101 * t75) * t134 + t74 * t138 + t111) * t74 + m(7) * (t13 * t34 + t14 * t33) + (-t21 * t74 - t22 * t75) * t133 + (t72 / 0.2e1 + t71 / 0.2e1) * (-Icges(6,6) * t85 - t139); m(7) * (-t33 * t97 + t34 * t95) + (-t74 * t95 + t75 * t97) * t133; -m(6) * t19 - m(7) * t9; m(7) * (-t33 * t74 + t34 * t75); m(6) * (t122 * t70 ^ 2 + t19 ^ 2) + m(7) * (t33 ^ 2 + t34 ^ 2 + t9 ^ 2) + (t71 * t36 + t4) * t75 + (t75 * t36 * t74 - t3 - t122 * (-Icges(6,3) * t74 + t100 * t75)) * t74; m(7) * (t13 * t17 + t14 * t18) + t20 + (t110 * t74 + t111 * t75) * t84; m(7) * (t17 * t95 - t18 * t97); -m(7) * t12; m(7) * (t17 * t75 - t18 * t74); t75 * t2 / 0.2e1 + t1 * t135 + (-t10 * t74 + t11 * t75) * t134 + m(7) * (t12 * t9 + t17 * t34 + t18 * t33) + (t4 * t135 - t75 * t3 / 0.2e1) * t84; t85 * t20 + m(7) * (t12 ^ 2 + t17 ^ 2 + t18 ^ 2) + (-t74 * t2 - t75 * t1 + t85 * (-t10 * t75 - t11 * t74)) * t84;];
%% Postprocessing: Reshape Output
% From vec2symmat_6_matlab.m
res = [t35(1) t35(2) t35(4) t35(7) t35(11) t35(16); t35(2) t35(3) t35(5) t35(8) t35(12) t35(17); t35(4) t35(5) t35(6) t35(9) t35(13) t35(18); t35(7) t35(8) t35(9) t35(10) t35(14) t35(19); t35(11) t35(12) t35(13) t35(14) t35(15) t35(20); t35(16) t35(17) t35(18) t35(19) t35(20) t35(21);];
Mq  = res;
