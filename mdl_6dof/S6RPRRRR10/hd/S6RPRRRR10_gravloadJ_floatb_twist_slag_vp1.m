% Calculate Gravitation load on the joints for
% S6RPRRRR10
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% g [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [13x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,alpha2,alpha3,d1,d3,d4,d5,d6,theta2]';
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
% Datum: 2019-03-09 07:37
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function taug = S6RPRRRR10_gravloadJ_floatb_twist_slag_vp1(qJ, g, ...
  pkin, m, rSges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(3,1),zeros(13,1),zeros(7,1),zeros(7,3)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPRRRR10_gravloadJ_floatb_twist_slag_vp1: qJ has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6RPRRRR10_gravloadJ_floatb_twist_slag_vp1: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [13 1]), ...
  'S6RPRRRR10_gravloadJ_floatb_twist_slag_vp1: pkin has to be [13x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6RPRRRR10_gravloadJ_floatb_twist_slag_vp1: m has to be [7x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [7,3]), ...
  'S6RPRRRR10_gravloadJ_floatb_twist_slag_vp1: rSges has to be [7x3] (double)');

%% Symbolic Calculation
% From gravload_joint_floatb_twist_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 07:27:46
% EndTime: 2019-03-09 07:27:51
% DurationCPUTime: 1.70s
% Computational Cost: add. (1250->185), mult. (2927->272), div. (0->0), fcn. (3718->16), ass. (0->86)
t145 = cos(qJ(3));
t126 = cos(pkin(13));
t128 = cos(pkin(6));
t114 = t128 * t126;
t123 = sin(pkin(13));
t144 = sin(qJ(1));
t83 = cos(qJ(1));
t103 = -t114 * t83 + t144 * t123;
t124 = sin(pkin(7));
t125 = sin(pkin(6));
t111 = t125 * t124;
t127 = cos(pkin(7));
t165 = t103 * t127 + t83 * t111;
t113 = t128 * t123;
t62 = t113 * t83 + t126 * t144;
t80 = sin(qJ(3));
t39 = -t145 * t62 + t165 * t80;
t112 = t127 * t125;
t52 = -t103 * t124 + t83 * t112;
t77 = qJ(4) + qJ(5);
t74 = sin(t77);
t75 = cos(t77);
t16 = t39 * t75 + t52 * t74;
t36 = t145 * t165 + t62 * t80;
t78 = sin(qJ(6));
t81 = cos(qJ(6));
t169 = t16 * t78 + t36 * t81;
t168 = t16 * t81 - t36 * t78;
t79 = sin(qJ(4));
t140 = t52 * t79;
t82 = cos(qJ(4));
t164 = t39 * t82 + t140;
t158 = t39 * t79 - t52 * t82;
t163 = t39 * t74 - t52 * t75;
t147 = pkin(12) + rSges(7,3);
t94 = t114 * t144 + t123 * t83;
t156 = t144 * t111 - t94 * t127;
t63 = -t113 * t144 + t126 * t83;
t41 = t63 * t145 + t156 * t80;
t87 = t144 * t112 + t94 * t124;
t19 = -t41 * t79 + t87 * t82;
t110 = t125 * t123;
t154 = t126 * t112 + t124 * t128;
t51 = t145 * t110 + t154 * t80;
t61 = -t111 * t126 + t127 * t128;
t159 = -t51 * t79 + t61 * t82;
t157 = t81 * rSges(7,1) - t78 * rSges(7,2) + pkin(5);
t40 = -t156 * t145 + t63 * t80;
t50 = t80 * t110 - t154 * t145;
t153 = g(1) * t40 + g(2) * t36 + g(3) * t50;
t152 = rSges(6,1) * t163 + rSges(6,2) * t16;
t149 = t75 * pkin(5);
t148 = pkin(10) + rSges(5,3);
t17 = t41 * t74 - t75 * t87;
t18 = t41 * t75 + t74 * t87;
t146 = -t17 * rSges(6,1) - t18 * rSges(6,2);
t137 = t75 * t78;
t136 = t75 * t81;
t30 = -t51 * t74 + t61 * t75;
t31 = t51 * t75 + t61 * t74;
t133 = t30 * rSges(6,1) - t31 * rSges(6,2);
t73 = pkin(4) * t82 + pkin(3);
t84 = -pkin(11) - pkin(10);
t132 = -t36 * t73 + t39 * t84;
t131 = -t40 * t73 - t41 * t84;
t130 = -t50 * t73 - t51 * t84;
t116 = t125 * t144;
t129 = t83 * pkin(1) + qJ(2) * t116;
t122 = t83 * t125;
t121 = t158 * pkin(4);
t120 = t19 * pkin(4);
t119 = t159 * pkin(4);
t118 = -pkin(1) * t144 + qJ(2) * t122;
t115 = -rSges(6,1) * t75 + rSges(6,2) * t74;
t101 = -t147 * t16 + t157 * t163;
t100 = t147 * t18 - t157 * t17;
t99 = t147 * t31 + t157 * t30;
t91 = -t62 * pkin(2) + t52 * pkin(9) + t118;
t90 = t63 * pkin(2) + t87 * pkin(9) + t129;
t89 = pkin(4) * t140 + t36 * t84 + t39 * t73 + t91;
t86 = t87 * t79;
t88 = pkin(4) * t86 - t40 * t84 + t41 * t73 + t90;
t20 = t41 * t82 + t86;
t2 = t18 * t81 + t40 * t78;
t1 = -t18 * t78 + t40 * t81;
t3 = [-m(2) * (g(1) * (-rSges(2,1) * t144 - t83 * rSges(2,2)) + g(2) * (t83 * rSges(2,1) - rSges(2,2) * t144)) - m(3) * (g(1) * (-t62 * rSges(3,1) + rSges(3,2) * t103 + rSges(3,3) * t122 + t118) + g(2) * (t63 * rSges(3,1) - rSges(3,2) * t94 + rSges(3,3) * t116 + t129)) - m(4) * (g(1) * (t39 * rSges(4,1) + rSges(4,2) * t36 + t52 * rSges(4,3) + t91) + g(2) * (t41 * rSges(4,1) - t40 * rSges(4,2) + rSges(4,3) * t87 + t90)) - m(5) * (g(1) * (rSges(5,1) * t164 - rSges(5,2) * t158 + t39 * pkin(3) - t148 * t36 + t91) + g(2) * (t20 * rSges(5,1) + t19 * rSges(5,2) + t41 * pkin(3) + t148 * t40 + t90)) - m(6) * (g(1) * (t16 * rSges(6,1) - rSges(6,2) * t163 - rSges(6,3) * t36 + t89) + g(2) * (t18 * rSges(6,1) - t17 * rSges(6,2) + t40 * rSges(6,3) + t88)) - m(7) * (g(1) * (t168 * rSges(7,1) - t169 * rSges(7,2) + t16 * pkin(5) + t147 * t163 + t89) + g(2) * (t2 * rSges(7,1) + t1 * rSges(7,2) + t18 * pkin(5) + t147 * t17 + t88)) (-m(3) - m(4) - m(5) - m(6) - m(7)) * (g(1) * t116 - g(2) * t122 + g(3) * t128) -m(4) * (g(1) * (-rSges(4,1) * t40 - rSges(4,2) * t41) + g(2) * (-rSges(4,1) * t36 + rSges(4,2) * t39) + g(3) * (-rSges(4,1) * t50 - rSges(4,2) * t51)) - m(5) * ((g(1) * t41 - g(2) * t39 + g(3) * t51) * t148 + t153 * (-t82 * rSges(5,1) + t79 * rSges(5,2) - pkin(3))) - m(6) * (g(1) * (rSges(6,3) * t41 + t115 * t40 + t131) + g(2) * (-rSges(6,3) * t39 + t115 * t36 + t132) + g(3) * (rSges(6,3) * t51 + t115 * t50 + t130)) + (-g(1) * (-t40 * t149 + (-t136 * t40 + t41 * t78) * rSges(7,1) + (t137 * t40 + t41 * t81) * rSges(7,2) + t131) - g(2) * (-t36 * t149 + (-t136 * t36 - t39 * t78) * rSges(7,1) + (t137 * t36 - t39 * t81) * rSges(7,2) + t132) - g(3) * (-t50 * t149 + (-t136 * t50 + t51 * t78) * rSges(7,1) + (t137 * t50 + t51 * t81) * rSges(7,2) + t130) + t153 * t74 * t147) * m(7), -m(5) * (g(1) * (rSges(5,1) * t19 - rSges(5,2) * t20) + g(2) * (t158 * rSges(5,1) + rSges(5,2) * t164) + g(3) * (t159 * rSges(5,1) + (-t51 * t82 - t61 * t79) * rSges(5,2))) - m(6) * (g(1) * (t120 + t146) + g(2) * (t121 + t152) + g(3) * (t119 + t133)) - m(7) * (g(1) * (t100 + t120) + g(2) * (t101 + t121) + g(3) * (t119 + t99)) -m(6) * (g(1) * t146 + g(2) * t152 + g(3) * t133) - m(7) * (g(1) * t100 + g(2) * t101 + g(3) * t99) -m(7) * (g(1) * (rSges(7,1) * t1 - rSges(7,2) * t2) + g(2) * (t169 * rSges(7,1) + t168 * rSges(7,2)) + g(3) * ((-t31 * t78 + t50 * t81) * rSges(7,1) + (-t31 * t81 - t50 * t78) * rSges(7,2)))];
taug  = t3(:);
