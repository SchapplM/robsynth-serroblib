% Calculate Gravitation load on the joints for
% S6RPRRPR9
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% g [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [13x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,alpha2,alpha3,d1,d3,d4,d6,theta2,theta5]';
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
% Datum: 2019-03-09 05:33
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function taug = S6RPRRPR9_gravloadJ_floatb_twist_slag_vp1(qJ, g, ...
  pkin, m, rSges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(3,1),zeros(13,1),zeros(7,1),zeros(7,3)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPRRPR9_gravloadJ_floatb_twist_slag_vp1: qJ has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6RPRRPR9_gravloadJ_floatb_twist_slag_vp1: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [13 1]), ...
  'S6RPRRPR9_gravloadJ_floatb_twist_slag_vp1: pkin has to be [13x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6RPRRPR9_gravloadJ_floatb_twist_slag_vp1: m has to be [7x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [7,3]), ...
  'S6RPRRPR9_gravloadJ_floatb_twist_slag_vp1: rSges has to be [7x3] (double)');

%% Symbolic Calculation
% From gravload_joint_floatb_twist_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 05:27:06
% EndTime: 2019-03-09 05:27:10
% DurationCPUTime: 1.72s
% Computational Cost: add. (1087->180), mult. (2598->266), div. (0->0), fcn. (3289->16), ass. (0->80)
t126 = cos(qJ(3));
t109 = cos(pkin(7));
t71 = cos(qJ(1));
t105 = sin(pkin(12));
t125 = sin(qJ(1));
t108 = cos(pkin(12));
t110 = cos(pkin(6));
t96 = t110 * t108;
t84 = t125 * t105 - t71 * t96;
t106 = sin(pkin(7));
t107 = sin(pkin(6));
t93 = t107 * t106;
t146 = t84 * t109 + t71 * t93;
t95 = t110 * t105;
t49 = t125 * t108 + t71 * t95;
t68 = sin(qJ(3));
t22 = t126 * t146 + t49 * t68;
t66 = sin(qJ(6));
t69 = cos(qJ(6));
t25 = -t49 * t126 + t146 * t68;
t94 = t109 * t107;
t38 = -t84 * t106 + t71 * t94;
t64 = qJ(4) + pkin(13);
t61 = sin(t64);
t62 = cos(t64);
t7 = t25 * t62 + t38 * t61;
t150 = t22 * t69 + t66 * t7;
t149 = -t22 * t66 + t69 * t7;
t67 = sin(qJ(4));
t121 = t38 * t67;
t70 = cos(qJ(4));
t145 = t25 * t70 + t121;
t139 = t25 * t67 - t38 * t70;
t144 = t25 * t61 - t38 * t62;
t78 = t71 * t105 + t125 * t96;
t40 = t78 * t106 + t125 * t94;
t137 = t106 * t110 + t108 * t94;
t92 = t107 * t105;
t37 = t126 * t92 + t137 * t68;
t48 = -t108 * t93 + t110 * t109;
t140 = -t37 * t67 + t48 * t70;
t138 = t78 * t109 - t125 * t93;
t50 = t71 * t108 - t125 * t95;
t27 = t50 * t126 - t138 * t68;
t10 = -t27 * t67 + t40 * t70;
t26 = t126 * t138 + t50 * t68;
t36 = -t126 * t137 + t68 * t92;
t134 = g(1) * t26 + g(2) * t22 + g(3) * t36;
t127 = pkin(11) + rSges(7,3);
t133 = -m(6) - m(7);
t129 = t62 * pkin(5);
t128 = pkin(10) + rSges(5,3);
t119 = t40 * t67;
t116 = t62 * t66;
t115 = t62 * t69;
t60 = pkin(4) * t70 + pkin(3);
t65 = -qJ(5) - pkin(10);
t114 = -t22 * t60 + t25 * t65;
t113 = -t26 * t60 - t27 * t65;
t112 = -t36 * t60 - t37 * t65;
t98 = t107 * t125;
t111 = t71 * pkin(1) + qJ(2) * t98;
t104 = t71 * t107;
t103 = t139 * pkin(4);
t102 = t10 * pkin(4);
t101 = t140 * pkin(4);
t100 = -t125 * pkin(1) + qJ(2) * t104;
t97 = -rSges(6,1) * t62 + rSges(6,2) * t61;
t75 = -t49 * pkin(2) + t38 * pkin(9) + t100;
t74 = t50 * pkin(2) + t40 * pkin(9) + t111;
t73 = pkin(4) * t121 + t22 * t65 + t25 * t60 + t75;
t72 = pkin(4) * t119 - t26 * t65 + t27 * t60 + t74;
t17 = t37 * t62 + t48 * t61;
t16 = -t37 * t61 + t48 * t62;
t11 = t27 * t70 + t119;
t9 = t27 * t62 + t40 * t61;
t8 = t27 * t61 - t40 * t62;
t2 = t26 * t66 + t69 * t9;
t1 = t26 * t69 - t66 * t9;
t3 = [-m(2) * (g(1) * (-t125 * rSges(2,1) - t71 * rSges(2,2)) + g(2) * (t71 * rSges(2,1) - t125 * rSges(2,2))) - m(3) * (g(1) * (-t49 * rSges(3,1) + t84 * rSges(3,2) + rSges(3,3) * t104 + t100) + g(2) * (t50 * rSges(3,1) - t78 * rSges(3,2) + rSges(3,3) * t98 + t111)) - m(4) * (g(1) * (t25 * rSges(4,1) + rSges(4,2) * t22 + t38 * rSges(4,3) + t75) + g(2) * (t27 * rSges(4,1) - t26 * rSges(4,2) + t40 * rSges(4,3) + t74)) - m(5) * (g(1) * (rSges(5,1) * t145 - rSges(5,2) * t139 + t25 * pkin(3) - t128 * t22 + t75) + g(2) * (t11 * rSges(5,1) + t10 * rSges(5,2) + t27 * pkin(3) + t128 * t26 + t74)) - m(6) * (g(1) * (t7 * rSges(6,1) - rSges(6,2) * t144 - rSges(6,3) * t22 + t73) + g(2) * (t9 * rSges(6,1) - t8 * rSges(6,2) + t26 * rSges(6,3) + t72)) - m(7) * (g(1) * (t149 * rSges(7,1) - t150 * rSges(7,2) + t7 * pkin(5) + t127 * t144 + t73) + g(2) * (t2 * rSges(7,1) + t1 * rSges(7,2) + t9 * pkin(5) + t127 * t8 + t72)) (-m(3) - m(4) - m(5) + t133) * (g(1) * t98 - g(2) * t104 + g(3) * t110) -m(4) * (g(1) * (-rSges(4,1) * t26 - rSges(4,2) * t27) + g(2) * (-rSges(4,1) * t22 + rSges(4,2) * t25) + g(3) * (-rSges(4,1) * t36 - rSges(4,2) * t37)) - m(5) * (t134 * (-t70 * rSges(5,1) + t67 * rSges(5,2) - pkin(3)) + (g(1) * t27 - g(2) * t25 + g(3) * t37) * t128) - m(6) * (g(1) * (rSges(6,3) * t27 + t97 * t26 + t113) + g(2) * (-rSges(6,3) * t25 + t97 * t22 + t114) + g(3) * (rSges(6,3) * t37 + t97 * t36 + t112)) + (-g(1) * (-t26 * t129 + (-t26 * t115 + t27 * t66) * rSges(7,1) + (t26 * t116 + t27 * t69) * rSges(7,2) + t113) - g(2) * (-t22 * t129 + (-t22 * t115 - t25 * t66) * rSges(7,1) + (t22 * t116 - t25 * t69) * rSges(7,2) + t114) - g(3) * (-t36 * t129 + (-t36 * t115 + t37 * t66) * rSges(7,1) + (t36 * t116 + t37 * t69) * rSges(7,2) + t112) + t134 * t61 * t127) * m(7), -m(5) * (g(1) * (rSges(5,1) * t10 - rSges(5,2) * t11) + g(2) * (t139 * rSges(5,1) + rSges(5,2) * t145) + g(3) * (t140 * rSges(5,1) + (-t37 * t70 - t48 * t67) * rSges(5,2))) - m(6) * (g(1) * (-rSges(6,1) * t8 - rSges(6,2) * t9 + t102) + g(2) * (rSges(6,1) * t144 + rSges(6,2) * t7 + t103) + g(3) * (rSges(6,1) * t16 - rSges(6,2) * t17 + t101)) + (-g(1) * (t127 * t9 + t102) - g(2) * (-t127 * t7 + t103) - g(3) * (t127 * t17 + t101) - (-g(1) * t8 + g(2) * t144 + g(3) * t16) * (t69 * rSges(7,1) - t66 * rSges(7,2) + pkin(5))) * m(7), t133 * t134, -m(7) * (g(1) * (rSges(7,1) * t1 - rSges(7,2) * t2) + g(2) * (t150 * rSges(7,1) + t149 * rSges(7,2)) + g(3) * ((-t17 * t66 + t36 * t69) * rSges(7,1) + (-t17 * t69 - t36 * t66) * rSges(7,2)))];
taug  = t3(:);
