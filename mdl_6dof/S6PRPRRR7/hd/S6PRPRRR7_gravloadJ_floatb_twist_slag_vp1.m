% Calculate Gravitation load on the joints for
% S6PRPRRR7
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% g [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [14x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,alpha2,alpha3,alpha4,d2,d4,d5,d6,theta1,theta3]';
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
% Datum: 2019-03-08 20:57
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function taug = S6PRPRRR7_gravloadJ_floatb_twist_slag_vp1(qJ, g, ...
  pkin, m, rSges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(3,1),zeros(14,1),zeros(7,1),zeros(7,3)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6PRPRRR7_gravloadJ_floatb_twist_slag_vp1: qJ has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6PRPRRR7_gravloadJ_floatb_twist_slag_vp1: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [14 1]), ...
  'S6PRPRRR7_gravloadJ_floatb_twist_slag_vp1: pkin has to be [14x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6PRPRRR7_gravloadJ_floatb_twist_slag_vp1: m has to be [7x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [7,3]), ...
  'S6PRPRRR7_gravloadJ_floatb_twist_slag_vp1: rSges has to be [7x3] (double)');

%% Symbolic Calculation
% From gravload_joint_floatb_twist_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-08 20:50:39
% EndTime: 2019-03-08 20:50:42
% DurationCPUTime: 1.28s
% Computational Cost: add. (1810->216), mult. (5189->349), div. (0->0), fcn. (6752->18), ass. (0->119)
t144 = rSges(7,3) + pkin(12);
t129 = sin(pkin(8));
t133 = cos(pkin(8));
t130 = sin(pkin(6));
t131 = cos(pkin(14));
t111 = t131 * t130;
t71 = sin(pkin(7));
t104 = t71 * t111;
t127 = sin(pkin(14));
t128 = sin(pkin(13));
t134 = cos(pkin(7));
t135 = cos(pkin(6));
t115 = t135 * t128;
t132 = cos(pkin(13));
t141 = sin(qJ(2));
t143 = cos(qJ(2));
t65 = -t115 * t141 + t132 * t143;
t94 = t115 * t143 + t132 * t141;
t91 = t94 * t131;
t79 = -t104 * t128 + t127 * t65 + t134 * t91;
t113 = t134 * t130;
t87 = t113 * t128 + t71 * t94;
t149 = -t87 * t129 + t79 * t133;
t116 = t135 * t132;
t64 = t116 * t141 + t128 * t143;
t93 = -t116 * t143 + t128 * t141;
t89 = t93 * t131;
t80 = t104 * t132 + t127 * t64 + t134 * t89;
t86 = -t113 * t132 + t71 * t93;
t148 = -t86 * t129 + t80 * t133;
t101 = t134 * t111;
t110 = t130 * t127;
t124 = t71 * t135;
t85 = -t101 * t143 + t110 * t141 - t124 * t131;
t119 = t143 * t130;
t95 = -t119 * t71 + t134 * t135;
t147 = -t95 * t129 + t85 * t133;
t76 = cos(qJ(5));
t146 = t76 * pkin(5);
t145 = rSges(6,3) + pkin(11);
t142 = cos(qJ(4));
t72 = sin(qJ(6));
t140 = t72 * t76;
t75 = cos(qJ(6));
t139 = t75 * t76;
t118 = t130 * t141;
t105 = t71 * t118;
t138 = pkin(2) * t119 + qJ(3) * t105;
t137 = qJ(3) * t71;
t126 = -m(4) - m(5) - m(6) - m(7);
t60 = -t101 * t141 - t110 * t143;
t125 = t60 * t129;
t123 = t71 * t133;
t122 = t71 * t129;
t100 = t134 * t110;
t61 = -t100 * t141 + t111 * t143;
t99 = t133 * t105;
t121 = t61 * pkin(3) + pkin(10) * t99 + t138;
t120 = t133 * t142;
t73 = sin(qJ(5));
t117 = -rSges(6,1) * t76 + rSges(6,2) * t73;
t114 = t134 * t131;
t112 = t134 * t127;
t48 = -t112 * t64 - t89;
t62 = t93 * pkin(2);
t109 = t48 * pkin(3) + t137 * t64 - t62;
t50 = -t112 * t65 - t91;
t63 = t94 * pkin(2);
t108 = t50 * pkin(3) + t137 * t65 - t63;
t107 = rSges(7,1) * t75 - rSges(7,2) * t72 + pkin(5);
t106 = t142 * t122;
t103 = t71 * t110;
t74 = sin(qJ(4));
t98 = t129 * t105;
t36 = -t120 * t60 - t142 * t98 + t61 * t74;
t37 = t61 * t142 + (t133 * t60 + t98) * t74;
t102 = t37 * pkin(4) + t36 * pkin(11) + t121;
t88 = t93 * t127;
t47 = -t114 * t64 + t88;
t19 = -t106 * t64 - t120 * t47 + t48 * t74;
t20 = t48 * t142 + (t122 * t64 + t133 * t47) * t74;
t97 = t20 * pkin(4) + t19 * pkin(11) + t109;
t90 = t94 * t127;
t49 = -t114 * t65 + t90;
t21 = -t106 * t65 - t120 * t49 + t50 * t74;
t22 = t50 * t142 + (t122 * t65 + t133 * t49) * t74;
t96 = t22 * pkin(4) + t21 * pkin(11) + t108;
t34 = t123 * t64 - t129 * t47;
t35 = t123 * t65 - t129 * t49;
t82 = (g(1) * t35 + g(2) * t34 - g(3) * t125) * pkin(10);
t56 = t100 * t143 + t111 * t141 + t124 * t127;
t52 = -t125 + t99;
t42 = t103 * t128 + t131 * t65 - t134 * t90;
t41 = -t103 * t132 + t131 * t64 - t134 * t88;
t40 = t129 * t85 + t133 * t95;
t29 = t129 * t79 + t133 * t87;
t28 = t129 * t80 + t133 * t86;
t27 = t56 * t142 - t147 * t74;
t26 = t142 * t147 + t56 * t74;
t25 = t26 * pkin(4);
t24 = t37 * t76 + t52 * t73;
t23 = t37 * t73 - t52 * t76;
t16 = t42 * t142 - t149 * t74;
t15 = t142 * t149 + t42 * t74;
t14 = t41 * t142 - t148 * t74;
t13 = t142 * t148 + t41 * t74;
t12 = t15 * pkin(4);
t11 = t13 * pkin(4);
t10 = t27 * t76 + t40 * t73;
t9 = -t27 * t73 + t40 * t76;
t8 = t22 * t76 + t35 * t73;
t7 = t22 * t73 - t35 * t76;
t6 = t20 * t76 + t34 * t73;
t5 = t20 * t73 - t34 * t76;
t4 = t16 * t76 + t29 * t73;
t3 = -t16 * t73 + t29 * t76;
t2 = t14 * t76 + t28 * t73;
t1 = -t14 * t73 + t28 * t76;
t17 = [(-m(2) - m(3) + t126) * g(3), -m(3) * (g(1) * (-rSges(3,1) * t94 - t65 * rSges(3,2)) + g(2) * (-rSges(3,1) * t93 - t64 * rSges(3,2)) + g(3) * (rSges(3,1) * t119 - rSges(3,2) * t118)) - m(4) * (g(1) * (rSges(4,1) * t50 + rSges(4,2) * t49 - t63) + g(2) * (rSges(4,1) * t48 + rSges(4,2) * t47 - t62) + g(3) * (t61 * rSges(4,1) + t60 * rSges(4,2) + t138) + (rSges(4,3) * g(3) * t118 + (g(1) * t65 + g(2) * t64) * (rSges(4,3) + qJ(3))) * t71) - m(5) * (g(1) * (t22 * rSges(5,1) - t21 * rSges(5,2) + t35 * rSges(5,3) + t108) + g(2) * (t20 * rSges(5,1) - t19 * rSges(5,2) + t34 * rSges(5,3) + t109) + g(3) * (t37 * rSges(5,1) - t36 * rSges(5,2) + t52 * rSges(5,3) + t121) + t82) - m(6) * (g(1) * (t8 * rSges(6,1) - t7 * rSges(6,2) + t21 * rSges(6,3) + t96) + g(2) * (t6 * rSges(6,1) - t5 * rSges(6,2) + t19 * rSges(6,3) + t97) + g(3) * (t24 * rSges(6,1) - t23 * rSges(6,2) + t36 * rSges(6,3) + t102) + t82) - m(7) * (g(1) * (t8 * pkin(5) + (t21 * t72 + t75 * t8) * rSges(7,1) + (t21 * t75 - t72 * t8) * rSges(7,2) + t96 + t144 * t7) + g(2) * (t6 * pkin(5) + (t19 * t72 + t6 * t75) * rSges(7,1) + (t19 * t75 - t6 * t72) * rSges(7,2) + t97 + t144 * t5) + g(3) * (t24 * pkin(5) + (t24 * t75 + t36 * t72) * rSges(7,1) + (-t24 * t72 + t36 * t75) * rSges(7,2) + t102 + t144 * t23) + t82) t126 * (g(1) * t87 + g(2) * t86 + g(3) * t95) -m(5) * (g(1) * (-rSges(5,1) * t15 - rSges(5,2) * t16) + g(2) * (-rSges(5,1) * t13 - rSges(5,2) * t14) + g(3) * (-rSges(5,1) * t26 - rSges(5,2) * t27)) - m(6) * (g(1) * (t117 * t15 + t145 * t16 - t12) + g(2) * (t117 * t13 + t14 * t145 - t11) + g(3) * (t117 * t26 + t145 * t27 - t25)) + (-g(1) * (-t15 * t146 - t12 + t16 * pkin(11) + (-t139 * t15 + t16 * t72) * rSges(7,1) + (t140 * t15 + t16 * t75) * rSges(7,2)) - g(2) * (-t13 * t146 - t11 + t14 * pkin(11) + (-t13 * t139 + t14 * t72) * rSges(7,1) + (t13 * t140 + t14 * t75) * rSges(7,2)) - g(3) * (-t26 * t146 - t25 + t27 * pkin(11) + (-t139 * t26 + t27 * t72) * rSges(7,1) + (t140 * t26 + t27 * t75) * rSges(7,2)) - (-g(1) * t15 - g(2) * t13 - g(3) * t26) * t73 * t144) * m(7), -m(6) * (g(1) * (rSges(6,1) * t3 - rSges(6,2) * t4) + g(2) * (rSges(6,1) * t1 - rSges(6,2) * t2) + g(3) * (rSges(6,1) * t9 - rSges(6,2) * t10)) - m(7) * (g(3) * (t10 * t144 + t107 * t9) + (t107 * t1 + t144 * t2) * g(2) + (t107 * t3 + t144 * t4) * g(1)) -m(7) * (g(1) * ((t15 * t75 - t4 * t72) * rSges(7,1) + (-t15 * t72 - t4 * t75) * rSges(7,2)) + g(2) * ((t13 * t75 - t2 * t72) * rSges(7,1) + (-t13 * t72 - t2 * t75) * rSges(7,2)) + g(3) * ((-t10 * t72 + t26 * t75) * rSges(7,1) + (-t10 * t75 - t26 * t72) * rSges(7,2)))];
taug  = t17(:);
