% Calculate Gravitation load on the joints for
% S6PRRRRP5
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% g [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [12x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,alpha2,alpha3,d2,d3,d4,d5,theta1]';
% m_mdh [7x1]
%   mass of all robot links (including the base)
% rSges [7x3]
%   center of mass of all robot links (in body frames)
%   rows: links of the robot (starting with base)
%   columns: x-, y-, z-coordinates
% 
% Output:
% taug [6x1]
%   joint torques required to compensate gravitation load

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-09 00:27
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function taug = S6PRRRRP5_gravloadJ_floatb_twist_slag_vp1(qJ, g, ...
  pkin, m, rSges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(3,1),zeros(12,1),zeros(7,1),zeros(7,3)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6PRRRRP5_gravloadJ_floatb_twist_slag_vp1: qJ has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6PRRRRP5_gravloadJ_floatb_twist_slag_vp1: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [12 1]), ...
  'S6PRRRRP5_gravloadJ_floatb_twist_slag_vp1: pkin has to be [12x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6PRRRRP5_gravloadJ_floatb_twist_slag_vp1: m has to be [7x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [7,3]), ...
  'S6PRRRRP5_gravloadJ_floatb_twist_slag_vp1: rSges has to be [7x3] (double)');

%% Symbolic Calculation
% From gravload_joint_floatb_twist_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 00:20:24
% EndTime: 2019-03-09 00:20:27
% DurationCPUTime: 1.17s
% Computational Cost: add. (1175->215), mult. (3205->323), div. (0->0), fcn. (4076->14), ass. (0->110)
t114 = cos(pkin(12));
t115 = cos(pkin(7));
t111 = sin(pkin(12));
t125 = sin(qJ(2));
t127 = cos(qJ(2));
t116 = cos(pkin(6));
t97 = t116 * t114;
t80 = t111 * t125 - t127 * t97;
t112 = sin(pkin(7));
t113 = sin(pkin(6));
t94 = t113 * t112;
t143 = t114 * t94 + t80 * t115;
t96 = t116 * t111;
t81 = t114 * t125 + t127 * t96;
t93 = t113 * t111;
t142 = -t112 * t93 + t81 * t115;
t95 = t115 * t113;
t141 = t112 * t116 + t127 * t95;
t126 = cos(qJ(3));
t58 = t111 * t127 + t125 * t97;
t72 = sin(qJ(3));
t30 = t58 * t126 - t143 * t72;
t71 = sin(qJ(4));
t74 = cos(qJ(4));
t75 = t80 * t112 - t114 * t95;
t16 = t30 * t74 + t75 * t71;
t59 = t114 * t127 - t125 * t96;
t32 = t59 * t126 - t142 * t72;
t76 = t81 * t112 + t115 * t93;
t18 = t32 * t74 + t76 * t71;
t99 = t113 * t125;
t47 = t126 * t99 + t141 * t72;
t79 = t116 * t115 - t127 * t94;
t34 = t47 * t74 + t79 * t71;
t140 = g(1) * t18 + g(2) * t16 + g(3) * t34;
t15 = t30 * t71 - t75 * t74;
t17 = t32 * t71 - t76 * t74;
t33 = t47 * t71 - t79 * t74;
t139 = g(1) * t17 + g(2) * t15 + g(3) * t33;
t29 = t143 * t126 + t58 * t72;
t31 = t142 * t126 + t59 * t72;
t46 = -t141 * t126 + t72 * t99;
t138 = t71 * (-g(1) * t31 - g(2) * t29 - g(3) * t46);
t137 = pkin(4) * t74;
t129 = rSges(5,3) + pkin(10);
t128 = rSges(6,3) + pkin(11);
t70 = sin(qJ(5));
t124 = t30 * t70;
t123 = t32 * t70;
t122 = t47 * t70;
t73 = cos(qJ(5));
t68 = pkin(5) * t73 + pkin(4);
t121 = t68 * t74;
t120 = t70 * t74;
t119 = t73 * t74;
t118 = rSges(7,3) + qJ(6) + pkin(11);
t100 = t127 * t113;
t85 = t125 * t94;
t117 = pkin(2) * t100 + pkin(9) * t85;
t86 = t125 * t95;
t54 = t126 * t100 - t72 * t86;
t110 = t54 * pkin(3) + t117;
t109 = pkin(5) * t70 + pkin(10);
t27 = t29 * pkin(3);
t108 = t30 * pkin(10) - t27;
t28 = t31 * pkin(3);
t107 = t32 * pkin(10) - t28;
t45 = t46 * pkin(3);
t106 = t47 * pkin(10) - t45;
t105 = t112 * pkin(9);
t104 = t71 * t112;
t103 = t72 * t115;
t102 = t74 * t112;
t101 = t115 * t126;
t98 = -rSges(5,1) * t74 + rSges(5,2) * t71;
t1 = -t16 * t70 + t29 * t73;
t3 = -t18 * t70 + t31 * t73;
t13 = -t34 * t70 + t46 * t73;
t38 = -t58 * t103 - t80 * t126;
t56 = t80 * pkin(2);
t89 = t38 * pkin(3) + t58 * t105 - t56;
t40 = -t59 * t103 - t81 * t126;
t57 = t81 * pkin(2);
t88 = t40 * pkin(3) + t59 * t105 - t57;
t82 = t112 * rSges(4,3) + t105;
t53 = t72 * t100 + t126 * t86;
t42 = t54 * t74 + t71 * t85;
t41 = t54 * t71 - t74 * t85;
t39 = t59 * t101 - t81 * t72;
t37 = t58 * t101 - t80 * t72;
t26 = t42 * t73 + t53 * t70;
t25 = -t42 * t70 + t53 * t73;
t24 = t59 * t104 + t40 * t74;
t23 = -t59 * t102 + t40 * t71;
t22 = t58 * t104 + t38 * t74;
t21 = -t58 * t102 + t38 * t71;
t20 = -t46 * t119 + t122;
t19 = t46 * t120 + t47 * t73;
t14 = -t34 * t73 - t46 * t70;
t12 = t24 * t73 + t39 * t70;
t11 = -t24 * t70 + t39 * t73;
t10 = t22 * t73 + t37 * t70;
t9 = -t22 * t70 + t37 * t73;
t8 = -t31 * t119 + t123;
t7 = t31 * t120 + t32 * t73;
t6 = -t29 * t119 + t124;
t5 = t29 * t120 + t30 * t73;
t4 = -t18 * t73 - t31 * t70;
t2 = -t16 * t73 - t29 * t70;
t35 = [(-m(2) - m(3) - m(4) - m(5) - m(6) - m(7)) * g(3), -m(3) * (g(1) * (-t81 * rSges(3,1) - t59 * rSges(3,2)) + g(2) * (-t80 * rSges(3,1) - t58 * rSges(3,2)) + g(3) * (rSges(3,1) * t100 - rSges(3,2) * t99)) - m(4) * (g(1) * (t40 * rSges(4,1) - t39 * rSges(4,2) + t82 * t59 - t57) + g(2) * (t38 * rSges(4,1) - t37 * rSges(4,2) + t82 * t58 - t56) + g(3) * (t54 * rSges(4,1) - t53 * rSges(4,2) + rSges(4,3) * t85 + t117)) - m(5) * (g(1) * (t24 * rSges(5,1) - t23 * rSges(5,2) + t129 * t39 + t88) + g(2) * (t22 * rSges(5,1) - t21 * rSges(5,2) + t129 * t37 + t89) + g(3) * (rSges(5,1) * t42 - rSges(5,2) * t41 + t129 * t53 + t110)) - m(6) * (g(1) * (t12 * rSges(6,1) + t11 * rSges(6,2) + t24 * pkin(4) + t39 * pkin(10) + t128 * t23 + t88) + g(2) * (t10 * rSges(6,1) + t9 * rSges(6,2) + t22 * pkin(4) + t37 * pkin(10) + t128 * t21 + t89) + g(3) * (rSges(6,1) * t26 + rSges(6,2) * t25 + pkin(4) * t42 + pkin(10) * t53 + t128 * t41 + t110)) - m(7) * (g(1) * (t12 * rSges(7,1) + t11 * rSges(7,2) + t109 * t39 + t118 * t23 + t24 * t68 + t88) + g(2) * (t10 * rSges(7,1) + t9 * rSges(7,2) + t109 * t37 + t118 * t21 + t22 * t68 + t89) + g(3) * (rSges(7,1) * t26 + rSges(7,2) * t25 + t109 * t53 + t118 * t41 + t42 * t68 + t110)) -m(4) * (g(1) * (-rSges(4,1) * t31 - rSges(4,2) * t32) + g(2) * (-rSges(4,1) * t29 - rSges(4,2) * t30) + g(3) * (-rSges(4,1) * t46 - rSges(4,2) * t47)) - m(5) * (g(1) * (t129 * t32 + t98 * t31 - t28) + g(2) * (t129 * t30 + t98 * t29 - t27) + g(3) * (t129 * t47 + t98 * t46 - t45)) + (-g(1) * (t8 * rSges(7,1) + t7 * rSges(7,2) + pkin(5) * t123 - t31 * t121 + t107) - g(2) * (t6 * rSges(7,1) + t5 * rSges(7,2) + pkin(5) * t124 - t29 * t121 + t108) - g(3) * (t20 * rSges(7,1) + t19 * rSges(7,2) + pkin(5) * t122 - t46 * t121 + t106) - t118 * t138) * m(7) + (-g(1) * (t8 * rSges(6,1) + t7 * rSges(6,2) - t31 * t137 + t107) - g(2) * (t6 * rSges(6,1) + t5 * rSges(6,2) - t29 * t137 + t108) - g(3) * (t20 * rSges(6,1) + t19 * rSges(6,2) - t46 * t137 + t106) - t128 * t138) * m(6), -m(5) * (g(1) * (-rSges(5,1) * t17 - rSges(5,2) * t18) + g(2) * (-rSges(5,1) * t15 - rSges(5,2) * t16) + g(3) * (-rSges(5,1) * t33 - rSges(5,2) * t34)) - m(6) * (t139 * (-rSges(6,1) * t73 + rSges(6,2) * t70 - pkin(4)) + t140 * t128) - m(7) * (t139 * (-rSges(7,1) * t73 + rSges(7,2) * t70 - t68) + t140 * t118) -m(6) * (g(1) * (rSges(6,1) * t3 + rSges(6,2) * t4) + g(2) * (rSges(6,1) * t1 + rSges(6,2) * t2) + g(3) * (rSges(6,1) * t13 + rSges(6,2) * t14)) + (-g(1) * (rSges(7,1) * t3 + rSges(7,2) * t4) - g(2) * (rSges(7,1) * t1 + rSges(7,2) * t2) - g(3) * (rSges(7,1) * t13 + rSges(7,2) * t14) - (g(1) * t3 + g(2) * t1 + g(3) * t13) * pkin(5)) * m(7), -m(7) * t139];
taug  = t35(:);
