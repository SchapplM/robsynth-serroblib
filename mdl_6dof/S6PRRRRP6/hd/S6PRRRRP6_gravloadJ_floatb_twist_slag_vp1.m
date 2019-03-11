% Calculate Gravitation load on the joints for
% S6PRRRRP6
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
% Datum: 2019-03-09 00:36
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function taug = S6PRRRRP6_gravloadJ_floatb_twist_slag_vp1(qJ, g, ...
  pkin, m, rSges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(3,1),zeros(12,1),zeros(7,1),zeros(7,3)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6PRRRRP6_gravloadJ_floatb_twist_slag_vp1: qJ has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6PRRRRP6_gravloadJ_floatb_twist_slag_vp1: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [12 1]), ...
  'S6PRRRRP6_gravloadJ_floatb_twist_slag_vp1: pkin has to be [12x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6PRRRRP6_gravloadJ_floatb_twist_slag_vp1: m has to be [7x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [7,3]), ...
  'S6PRRRRP6_gravloadJ_floatb_twist_slag_vp1: rSges has to be [7x3] (double)');

%% Symbolic Calculation
% From gravload_joint_floatb_twist_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 00:28:44
% EndTime: 2019-03-09 00:28:46
% DurationCPUTime: 1.04s
% Computational Cost: add. (1325->215), mult. (3662->327), div. (0->0), fcn. (4689->14), ass. (0->112)
t146 = rSges(7,2) + pkin(11);
t127 = sin(pkin(6));
t128 = cos(pkin(12));
t109 = t128 * t127;
t129 = cos(pkin(7));
t89 = sin(pkin(7));
t130 = cos(pkin(6));
t112 = t130 * t128;
t126 = sin(pkin(12));
t140 = sin(qJ(2));
t142 = cos(qJ(2));
t97 = -t112 * t142 + t126 * t140;
t152 = t89 * t109 + t97 * t129;
t108 = t127 * t126;
t111 = t130 * t126;
t98 = t111 * t142 + t128 * t140;
t151 = t89 * t108 - t98 * t129;
t110 = t129 * t127;
t150 = t142 * t110 + t130 * t89;
t94 = cos(qJ(4));
t149 = pkin(4) * t94;
t148 = pkin(9) * t89;
t147 = rSges(7,1) + pkin(5);
t144 = rSges(5,3) + pkin(10);
t143 = rSges(6,3) + pkin(11);
t141 = cos(qJ(3));
t79 = t112 * t140 + t126 * t142;
t92 = sin(qJ(3));
t42 = t141 * t152 + t79 * t92;
t91 = sin(qJ(4));
t139 = t42 * t91;
t80 = -t111 * t140 + t128 * t142;
t44 = -t141 * t151 + t80 * t92;
t138 = t44 * t91;
t120 = t127 * t140;
t64 = t120 * t92 - t141 * t150;
t137 = t64 * t91;
t136 = t89 * t91;
t135 = t89 * t94;
t90 = sin(qJ(5));
t134 = t90 * t94;
t93 = cos(qJ(5));
t133 = t93 * t94;
t107 = t89 * t120;
t121 = t142 * t127;
t132 = pkin(2) * t121 + pkin(9) * t107;
t131 = rSges(7,3) + qJ(6);
t102 = t140 * t110;
t75 = -t102 * t92 + t121 * t141;
t125 = t75 * pkin(3) + t132;
t124 = t92 * t129;
t122 = t129 * t141;
t119 = -rSges(5,1) * t94 + rSges(5,2) * t91;
t118 = rSges(6,1) * t93 - rSges(6,2) * t90;
t53 = -t124 * t79 - t141 * t97;
t76 = t97 * pkin(2);
t117 = t53 * pkin(3) + t148 * t79 - t76;
t55 = -t124 * t80 - t141 * t98;
t77 = t98 * pkin(2);
t116 = t55 * pkin(3) + t148 * t80 - t77;
t39 = t42 * pkin(3);
t43 = t141 * t79 - t152 * t92;
t115 = t43 * pkin(10) - pkin(11) * t139 - t42 * t149 - t39;
t40 = t44 * pkin(3);
t45 = t80 * t141 + t151 * t92;
t114 = t45 * pkin(10) - pkin(11) * t138 - t44 * t149 - t40;
t63 = t64 * pkin(3);
t65 = t141 * t120 + t150 * t92;
t113 = t65 * pkin(10) - pkin(11) * t137 - t64 * t149 - t63;
t58 = t107 * t91 + t75 * t94;
t74 = t102 * t141 + t121 * t92;
t106 = t58 * pkin(4) + pkin(10) * t74 + t125;
t26 = t136 * t79 + t53 * t94;
t52 = t122 * t79 - t92 * t97;
t100 = t26 * pkin(4) + pkin(10) * t52 + t117;
t28 = t136 * t80 + t55 * t94;
t54 = t122 * t80 - t92 * t98;
t99 = t28 * pkin(4) + pkin(10) * t54 + t116;
t78 = -t121 * t89 + t129 * t130;
t67 = t108 * t129 + t89 * t98;
t66 = -t109 * t129 + t89 * t97;
t57 = -t107 * t94 + t75 * t91;
t47 = t65 * t94 + t78 * t91;
t46 = -t65 * t91 + t78 * t94;
t41 = t46 * pkin(4);
t30 = t58 * t93 + t74 * t90;
t29 = t58 * t90 - t74 * t93;
t27 = -t135 * t80 + t55 * t91;
t25 = -t135 * t79 + t53 * t91;
t22 = -t133 * t64 + t65 * t90;
t21 = -t134 * t64 - t65 * t93;
t20 = t45 * t94 + t67 * t91;
t19 = -t45 * t91 + t67 * t94;
t18 = t43 * t94 + t66 * t91;
t17 = -t43 * t91 + t66 * t94;
t16 = t47 * t93 + t64 * t90;
t15 = t47 * t90 - t64 * t93;
t14 = t19 * pkin(4);
t13 = t17 * pkin(4);
t12 = t28 * t93 + t54 * t90;
t11 = t28 * t90 - t54 * t93;
t10 = t26 * t93 + t52 * t90;
t9 = t26 * t90 - t52 * t93;
t8 = -t133 * t44 + t45 * t90;
t7 = -t134 * t44 - t45 * t93;
t6 = -t133 * t42 + t43 * t90;
t5 = -t134 * t42 - t43 * t93;
t4 = t20 * t93 + t44 * t90;
t3 = t20 * t90 - t44 * t93;
t2 = t18 * t93 + t42 * t90;
t1 = t18 * t90 - t42 * t93;
t23 = [(-m(2) - m(3) - m(4) - m(5) - m(6) - m(7)) * g(3), -m(3) * (g(1) * (-rSges(3,1) * t98 - t80 * rSges(3,2)) + g(2) * (-rSges(3,1) * t97 - t79 * rSges(3,2)) + g(3) * (rSges(3,1) * t121 - rSges(3,2) * t120)) - m(4) * (g(1) * (rSges(4,1) * t55 - rSges(4,2) * t54 - t77) + g(2) * (rSges(4,1) * t53 - rSges(4,2) * t52 - t76) + g(3) * (t75 * rSges(4,1) - t74 * rSges(4,2) + t132) + (rSges(4,3) * g(3) * t120 + (g(1) * t80 + g(2) * t79) * (rSges(4,3) + pkin(9))) * t89) - m(5) * (g(1) * (rSges(5,1) * t28 - rSges(5,2) * t27 + t144 * t54 + t116) + g(2) * (rSges(5,1) * t26 - rSges(5,2) * t25 + t144 * t52 + t117) + g(3) * (rSges(5,1) * t58 - rSges(5,2) * t57 + t144 * t74 + t125)) - m(6) * (g(1) * (rSges(6,1) * t12 - rSges(6,2) * t11 + t143 * t27 + t99) + g(2) * (rSges(6,1) * t10 - rSges(6,2) * t9 + t143 * t25 + t100) + g(3) * (rSges(6,1) * t30 - rSges(6,2) * t29 + t143 * t57 + t106)) - m(7) * (g(1) * (t11 * t131 + t12 * t147 + t146 * t27 + t99) + g(2) * (t10 * t147 + t131 * t9 + t146 * t25 + t100) + g(3) * (t131 * t29 + t146 * t57 + t147 * t30 + t106)) -m(4) * (g(1) * (-rSges(4,1) * t44 - rSges(4,2) * t45) + g(2) * (-rSges(4,1) * t42 - rSges(4,2) * t43) + g(3) * (-rSges(4,1) * t64 - rSges(4,2) * t65)) - m(5) * (g(1) * (t119 * t44 + t144 * t45 - t40) + g(2) * (t119 * t42 + t144 * t43 - t39) + g(3) * (t119 * t64 + t144 * t65 - t63)) - m(6) * (g(1) * (rSges(6,1) * t8 - rSges(6,2) * t7 - rSges(6,3) * t138 + t114) + g(2) * (rSges(6,1) * t6 - rSges(6,2) * t5 - rSges(6,3) * t139 + t115) + g(3) * (rSges(6,1) * t22 - rSges(6,2) * t21 - rSges(6,3) * t137 + t113)) - m(7) * (g(1) * (-rSges(7,2) * t138 + t131 * t7 + t147 * t8 + t114) + g(2) * (-rSges(7,2) * t139 + t131 * t5 + t147 * t6 + t115) + g(3) * (-rSges(7,2) * t137 + t131 * t21 + t147 * t22 + t113)) -m(5) * (g(1) * (rSges(5,1) * t19 - rSges(5,2) * t20) + g(2) * (rSges(5,1) * t17 - rSges(5,2) * t18) + g(3) * (rSges(5,1) * t46 - rSges(5,2) * t47)) - m(6) * (g(1) * (t118 * t19 + t143 * t20 + t14) + g(2) * (t118 * t17 + t143 * t18 + t13) + g(3) * (t118 * t46 + t143 * t47 + t41)) + (-g(1) * (t146 * t20 + t14) - g(2) * (t146 * t18 + t13) - g(3) * (t146 * t47 + t41) - (g(1) * t19 + g(2) * t17 + g(3) * t46) * (t131 * t90 + t147 * t93)) * m(7), -m(6) * (g(1) * (-rSges(6,1) * t3 - rSges(6,2) * t4) + g(2) * (-rSges(6,1) * t1 - rSges(6,2) * t2) + g(3) * (-rSges(6,1) * t15 - rSges(6,2) * t16)) - m(7) * (g(1) * (t131 * t4 - t147 * t3) + g(2) * (-t1 * t147 + t131 * t2) + g(3) * (t131 * t16 - t147 * t15)) -m(7) * (g(1) * t3 + g(2) * t1 + g(3) * t15)];
taug  = t23(:);
