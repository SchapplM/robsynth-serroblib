% Calculate Gravitation load on the joints for
% S6RRRRPR9
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% g [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [12x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,alpha2,d1,d2,d3,d4,d6,theta5]';
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
% Datum: 2019-03-09 22:59
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function taug = S6RRRRPR9_gravloadJ_floatb_twist_slag_vp1(qJ, g, ...
  pkin, m, rSges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(3,1),zeros(12,1),zeros(7,1),zeros(7,3)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRRRPR9_gravloadJ_floatb_twist_slag_vp1: qJ has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6RRRRPR9_gravloadJ_floatb_twist_slag_vp1: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [12 1]), ...
  'S6RRRRPR9_gravloadJ_floatb_twist_slag_vp1: pkin has to be [12x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6RRRRPR9_gravloadJ_floatb_twist_slag_vp1: m has to be [7x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [7,3]), ...
  'S6RRRRPR9_gravloadJ_floatb_twist_slag_vp1: rSges has to be [7x3] (double)');

%% Symbolic Calculation
% From gravload_joint_floatb_twist_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 22:47:49
% EndTime: 2019-03-09 22:47:52
% DurationCPUTime: 1.44s
% Computational Cost: add. (956->206), mult. (1516->294), div. (0->0), fcn. (1782->14), ass. (0->87)
t118 = pkin(11) + qJ(5) + rSges(7,3);
t117 = qJ(5) + rSges(6,3);
t116 = cos(pkin(6));
t76 = sin(pkin(6));
t80 = sin(qJ(2));
t127 = t76 * t80;
t79 = sin(qJ(3));
t82 = cos(qJ(3));
t150 = t116 * t82 - t79 * t127;
t125 = t76 * t82;
t107 = t80 * t116;
t140 = cos(qJ(2));
t141 = cos(qJ(1));
t81 = sin(qJ(1));
t52 = -t81 * t107 + t141 * t140;
t28 = t81 * t125 - t52 * t79;
t75 = sin(pkin(12));
t77 = cos(pkin(12));
t149 = -rSges(6,1) * t77 + rSges(6,2) * t75 - pkin(4);
t66 = pkin(5) * t77 + pkin(4);
t73 = pkin(12) + qJ(6);
t68 = sin(t73);
t69 = cos(t73);
t148 = -rSges(7,1) * t69 + rSges(7,2) * t68 - t66;
t111 = t76 * t141;
t50 = t141 * t107 + t81 * t140;
t74 = qJ(3) + qJ(4);
t70 = sin(t74);
t71 = cos(t74);
t23 = -t70 * t111 + t50 * t71;
t100 = t116 * t140;
t49 = -t141 * t100 + t80 * t81;
t147 = t23 * t68 - t49 * t69;
t146 = -t23 * t69 - t49 * t68;
t145 = rSges(4,1) * t82 - rSges(4,2) * t79 + pkin(2);
t144 = pkin(5) * t75;
t143 = g(3) * t76;
t142 = -pkin(9) - rSges(4,3);
t135 = rSges(5,2) * t70;
t130 = t49 * t75;
t51 = t81 * t100 + t141 * t80;
t129 = t51 * t75;
t126 = t76 * t81;
t22 = t71 * t111 + t50 * t70;
t124 = -t22 * rSges(5,1) - t23 * rSges(5,2);
t26 = -t71 * t126 + t52 * t70;
t27 = t70 * t126 + t52 * t71;
t123 = -t26 * rSges(5,1) - t27 * rSges(5,2);
t67 = pkin(3) * t82 + pkin(2);
t83 = -pkin(10) - pkin(9);
t122 = -t49 * t67 - t50 * t83;
t121 = -t51 * t67 - t52 * t83;
t43 = -t116 * t71 + t70 * t127;
t44 = t116 * t70 + t71 * t127;
t120 = -t43 * rSges(5,1) - t44 * rSges(5,2);
t119 = t141 * pkin(1) + pkin(8) * t126;
t115 = t79 * t126;
t112 = t71 * t140;
t110 = t76 * t140;
t109 = -t81 * pkin(1) + pkin(8) * t111;
t60 = t79 * t111;
t108 = -t50 * t82 + t60;
t105 = t70 * t110;
t104 = t71 * t110;
t103 = t28 * pkin(3);
t53 = t67 * t110;
t102 = -t83 * t127 + t53;
t101 = pkin(3) * t115 - t51 * t83 + t52 * t67 + t119;
t99 = -rSges(5,1) * t71 + t135;
t98 = rSges(6,1) * t75 + rSges(6,2) * t77;
t97 = t150 * pkin(3);
t95 = pkin(3) * t60 + t49 * t83 - t50 * t67 + t109;
t94 = t118 * t23 + t148 * t22;
t93 = t118 * t27 + t148 * t26;
t92 = t118 * t44 + t148 * t43;
t91 = t82 * t111 + t50 * t79;
t90 = rSges(7,1) * t68 + rSges(7,2) * t69 + t144;
t89 = t91 * pkin(3);
t88 = t117 * t23 + t149 * t22;
t87 = t117 * t27 + t149 * t26;
t86 = t117 * t44 + t149 * t43;
t85 = -t117 * t70 + t149 * t71;
t84 = -t118 * t70 + t148 * t71;
t29 = t52 * t82 + t115;
t3 = t27 * t69 + t51 * t68;
t2 = -t27 * t68 + t51 * t69;
t1 = [-m(2) * (g(1) * (-t81 * rSges(2,1) - t141 * rSges(2,2)) + g(2) * (t141 * rSges(2,1) - t81 * rSges(2,2))) - m(3) * (g(1) * (-t50 * rSges(3,1) + t49 * rSges(3,2) + rSges(3,3) * t111 + t109) + g(2) * (rSges(3,1) * t52 - rSges(3,2) * t51 + rSges(3,3) * t126 + t119)) - m(4) * (g(1) * (t108 * rSges(4,1) + t91 * rSges(4,2) - t50 * pkin(2) + t142 * t49 + t109) + g(2) * (rSges(4,1) * t29 + rSges(4,2) * t28 + pkin(2) * t52 - t142 * t51 + t119)) - m(5) * (g(1) * (-rSges(5,1) * t23 + rSges(5,2) * t22 - rSges(5,3) * t49 + t95) + g(2) * (rSges(5,1) * t27 - rSges(5,2) * t26 + rSges(5,3) * t51 + t101)) - m(6) * (g(1) * (-t23 * pkin(4) + (-t23 * t77 - t130) * rSges(6,1) + (t23 * t75 - t49 * t77) * rSges(6,2) - t117 * t22 + t95) + g(2) * (t27 * pkin(4) + (t27 * t77 + t129) * rSges(6,1) + (-t27 * t75 + t51 * t77) * rSges(6,2) + t117 * t26 + t101)) - m(7) * (g(1) * (t146 * rSges(7,1) + t147 * rSges(7,2) - pkin(5) * t130 - t118 * t22 - t23 * t66 + t95) + g(2) * (rSges(7,1) * t3 + rSges(7,2) * t2 + pkin(5) * t129 + t118 * t26 + t27 * t66 + t101)) -m(3) * (g(1) * (-rSges(3,1) * t51 - rSges(3,2) * t52) + g(2) * (-rSges(3,1) * t49 - rSges(3,2) * t50) + (t140 * rSges(3,1) - rSges(3,2) * t80) * t143) - m(4) * (g(1) * (-t142 * t52 - t145 * t51) + (-t142 * t50 - t145 * t49) * g(2) + (t145 * t140 - t142 * t80) * t143) - m(5) * (g(1) * (rSges(5,3) * t52 + t99 * t51 + t121) + g(2) * (rSges(5,3) * t50 + t99 * t49 + t122) + g(3) * t53 + (rSges(5,1) * t112 - t140 * t135 + (rSges(5,3) - t83) * t80) * t143) - m(6) * (g(1) * (t85 * t51 + t98 * t52 + t121) + g(2) * (t85 * t49 + t98 * t50 + t122) + g(3) * (pkin(4) * t104 + t102 + ((t77 * t112 + t75 * t80) * rSges(6,1) + (-t75 * t112 + t77 * t80) * rSges(6,2)) * t76 + t117 * t105)) - m(7) * (g(1) * (t84 * t51 + t90 * t52 + t121) + g(2) * (t84 * t49 + t90 * t50 + t122) + g(3) * (t66 * t104 + t127 * t144 + t102 + ((t69 * t112 + t68 * t80) * rSges(7,1) + (-t68 * t112 + t69 * t80) * rSges(7,2)) * t76 + t118 * t105)) -m(4) * (g(1) * (rSges(4,1) * t28 - rSges(4,2) * t29) + g(2) * (-t91 * rSges(4,1) + t108 * rSges(4,2)) + g(3) * (t150 * rSges(4,1) + (-t116 * t79 - t80 * t125) * rSges(4,2))) - m(5) * (g(1) * (t103 + t123) + g(2) * (-t89 + t124) + g(3) * (t97 + t120)) - m(6) * (g(1) * (t103 + t87) + g(2) * (-t89 + t88) + g(3) * (t86 + t97)) - m(7) * (g(1) * (t103 + t93) + g(2) * (-t89 + t94) + g(3) * (t92 + t97)) -m(5) * (g(1) * t123 + g(2) * t124 + g(3) * t120) - m(6) * (g(1) * t87 + g(2) * t88 + g(3) * t86) - m(7) * (g(1) * t93 + g(2) * t94 + g(3) * t92) (-m(6) - m(7)) * (g(1) * t26 + g(2) * t22 + g(3) * t43) -m(7) * (g(1) * (rSges(7,1) * t2 - rSges(7,2) * t3) + g(2) * (-t147 * rSges(7,1) + t146 * rSges(7,2)) + g(3) * ((-t69 * t110 - t44 * t68) * rSges(7,1) + (t68 * t110 - t44 * t69) * rSges(7,2)))];
taug  = t1(:);
