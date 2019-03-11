% Calculate Gravitation load on the joints for
% S6RRRPPP1
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% g [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [10x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,alpha4,d1,d2,d3,theta4]';
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
% Datum: 2019-03-09 15:20
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function taug = S6RRRPPP1_gravloadJ_floatb_twist_slag_vp1(qJ, g, ...
  pkin, m, rSges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(3,1),zeros(10,1),zeros(7,1),zeros(7,3)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRRPPP1_gravloadJ_floatb_twist_slag_vp1: qJ has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6RRRPPP1_gravloadJ_floatb_twist_slag_vp1: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S6RRRPPP1_gravloadJ_floatb_twist_slag_vp1: pkin has to be [10x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6RRRPPP1_gravloadJ_floatb_twist_slag_vp1: m has to be [7x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [7,3]), ...
  'S6RRRPPP1_gravloadJ_floatb_twist_slag_vp1: rSges has to be [7x3] (double)');

%% Symbolic Calculation
% From gravload_joint_floatb_twist_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 15:14:29
% EndTime: 2019-03-09 15:14:34
% DurationCPUTime: 1.66s
% Computational Cost: add. (637->223), mult. (1689->304), div. (0->0), fcn. (1973->10), ass. (0->100)
t75 = sin(qJ(2));
t77 = cos(qJ(3));
t131 = t75 * t77;
t72 = cos(pkin(10));
t110 = t72 * t131;
t70 = sin(pkin(10));
t74 = sin(qJ(3));
t134 = t74 * t75;
t71 = sin(pkin(6));
t73 = cos(pkin(6));
t78 = cos(qJ(2));
t89 = t73 * t134 + t71 * t78;
t148 = -t70 * t89 + t110;
t141 = rSges(7,1) + pkin(5);
t114 = rSges(7,3) + qJ(6);
t76 = sin(qJ(1));
t79 = cos(qJ(1));
t150 = g(1) * t79 + g(2) * t76;
t126 = t79 * t74;
t128 = t77 * t78;
t41 = t76 * t128 - t126;
t143 = g(2) * t41;
t127 = t78 * t79;
t43 = t77 * t127 + t76 * t74;
t149 = g(1) * t43 + t143;
t147 = -m(6) - m(7);
t145 = g(1) * t76;
t67 = t78 * pkin(2);
t140 = t70 * t73;
t139 = t70 * t77;
t138 = t71 * t75;
t136 = t72 * t73;
t135 = t72 * t74;
t133 = t74 * t78;
t132 = t75 * t76;
t130 = t75 * t79;
t129 = t76 * t78;
t120 = qJ(4) * t73;
t105 = t78 * t120;
t59 = pkin(9) * t129;
t125 = t76 * t105 + t59;
t62 = pkin(9) * t127;
t124 = t79 * t105 + t62;
t123 = t75 * pkin(9) + t67;
t122 = t79 * pkin(1) + t76 * pkin(8);
t121 = qJ(4) * t71;
t119 = qJ(4) * t75;
t117 = rSges(7,2) + qJ(5);
t115 = rSges(6,3) + qJ(5);
t113 = g(3) * t131;
t111 = t71 * t132;
t108 = t72 * t138;
t106 = t74 * t121;
t58 = t73 * t119;
t104 = -pkin(1) - t67;
t102 = pkin(2) * t127 + pkin(9) * t130 + t122;
t101 = t77 * t71 * t119 - pkin(3) * t134;
t40 = t74 * t129 + t77 * t79;
t11 = t41 * t136 - t40 * t70;
t12 = -t41 * t140 - t40 * t72;
t35 = t40 * pkin(3);
t100 = t12 * pkin(4) + qJ(5) * t11 - t35;
t42 = -t78 * t126 + t76 * t77;
t13 = t43 * t136 + t42 * t70;
t14 = -t43 * t140 + t42 * t72;
t37 = t42 * pkin(3);
t99 = t14 * pkin(4) + qJ(5) * t13 + t37;
t98 = pkin(3) * t128 + t78 * t106 + t123 + t58;
t97 = rSges(3,1) * t78 - rSges(3,2) * t75;
t68 = t79 * pkin(8);
t95 = -t41 * pkin(3) - t40 * t121 + t68;
t81 = t70 * t131 + t72 * t89;
t17 = t81 * t76;
t18 = t148 * t76;
t94 = -t18 * pkin(4) - qJ(5) * t17 + t125;
t19 = t81 * t79;
t20 = t148 * t79;
t93 = -t20 * pkin(4) - t19 * qJ(5) + t124;
t23 = t72 * t128 - t133 * t140 + t70 * t138;
t92 = t23 * pkin(4) + t98;
t90 = t73 * t132 + t40 * t71;
t88 = t71 * t134 - t73 * t78;
t6 = -t70 * t111 + t40 * t140 - t41 * t72;
t31 = -t73 * t110 + t70 * t134;
t32 = (t73 * t139 + t135) * t75;
t86 = -t32 * pkin(4) - qJ(5) * t31 + t101;
t5 = t76 * t108 - t40 * t136 - t41 * t70;
t85 = t6 * pkin(4) + qJ(5) * t5 + t95;
t84 = t43 * pkin(3) - t42 * t121 + t79 * t58 + t102;
t8 = t43 * t72 + (t71 * t130 + t42 * t73) * t70;
t83 = t8 * pkin(4) + t84;
t82 = ((-pkin(9) - t120) * t75 + t104) * t145;
t80 = t150 * (-pkin(3) * t77 - pkin(2) - t106) * t75;
t39 = t71 * t133 + t73 * t75;
t34 = t88 * t79;
t33 = t88 * t76;
t25 = t73 * t130 - t42 * t71;
t22 = -t108 + (t73 * t135 + t139) * t78;
t7 = -t79 * t108 - t42 * t136 + t43 * t70;
t1 = [-m(2) * (g(1) * (-t76 * rSges(2,1) - rSges(2,2) * t79) + g(2) * (rSges(2,1) * t79 - t76 * rSges(2,2))) - m(3) * (g(1) * (rSges(3,3) * t79 + t68) + g(2) * (rSges(3,1) * t127 - rSges(3,2) * t130 + t122) + (g(1) * (-pkin(1) - t97) + g(2) * rSges(3,3)) * t76) - m(4) * (g(1) * (-rSges(4,1) * t41 + rSges(4,2) * t40 + t68) + g(2) * (t43 * rSges(4,1) + t42 * rSges(4,2) + rSges(4,3) * t130 + t102) + ((-rSges(4,3) - pkin(9)) * t75 + t104) * t145) - m(5) * (g(1) * (rSges(5,1) * t6 - rSges(5,2) * t5 - rSges(5,3) * t90 + t95) + g(2) * (rSges(5,1) * t8 - rSges(5,2) * t7 + rSges(5,3) * t25 + t84) + t82) - m(6) * (g(1) * (-rSges(6,1) * t90 - rSges(6,2) * t6 + rSges(6,3) * t5 + t85) + g(2) * (rSges(6,1) * t25 - rSges(6,2) * t8 + t115 * t7 + t83) + t82) - m(7) * (g(1) * (rSges(7,2) * t5 + t114 * t6 - t141 * t90 + t85) + g(2) * (t114 * t8 + t117 * t7 + t141 * t25 + t83) + t82) -m(3) * (g(3) * t97 + t150 * (-rSges(3,1) * t75 - rSges(3,2) * t78)) - m(4) * (g(1) * (rSges(4,3) * t127 + t62) + g(2) * (rSges(4,3) * t129 + t59) + g(3) * (rSges(4,1) * t128 - rSges(4,2) * t133 + t123) + (g(3) * rSges(4,3) + t150 * (-rSges(4,1) * t77 + rSges(4,2) * t74 - pkin(2))) * t75) - m(5) * (g(1) * (-t20 * rSges(5,1) + t19 * rSges(5,2) - t34 * rSges(5,3) + t124) + g(2) * (-rSges(5,1) * t18 + rSges(5,2) * t17 - rSges(5,3) * t33 + t125) + g(3) * (rSges(5,1) * t23 - rSges(5,2) * t22 + rSges(5,3) * t39 + t98) + t80) - m(6) * (g(1) * (-t34 * rSges(6,1) + t20 * rSges(6,2) - t19 * rSges(6,3) + t93) + g(2) * (-rSges(6,1) * t33 + rSges(6,2) * t18 - rSges(6,3) * t17 + t94) + g(3) * (rSges(6,1) * t39 - rSges(6,2) * t23 + t115 * t22 + t92) + t80) - m(7) * (g(1) * (-t19 * rSges(7,2) - t114 * t20 - t141 * t34 + t93) + g(2) * (-rSges(7,2) * t17 - t114 * t18 - t141 * t33 + t94) + g(3) * (t114 * t23 + t117 * t22 + t141 * t39 + t92) + t80) -m(4) * (g(1) * (rSges(4,1) * t42 - rSges(4,2) * t43) + g(2) * (-rSges(4,1) * t40 - rSges(4,2) * t41) + g(3) * (-rSges(4,1) * t74 - rSges(4,2) * t77) * t75) - m(5) * (g(1) * (rSges(5,1) * t14 - rSges(5,2) * t13 + t37) + g(2) * (rSges(5,1) * t12 - rSges(5,2) * t11 - t35) + g(3) * (-rSges(5,1) * t32 + rSges(5,2) * t31 + t101) + (rSges(5,3) * t113 + t149 * (rSges(5,3) + qJ(4))) * t71) - m(6) * (g(1) * (-rSges(6,2) * t14 + rSges(6,3) * t13 + t99) + g(2) * (-rSges(6,2) * t12 + rSges(6,3) * t11 + t100) + g(3) * (rSges(6,2) * t32 - rSges(6,3) * t31 + t86) + (rSges(6,1) * t113 + t149 * (rSges(6,1) + qJ(4))) * t71) - m(7) * (g(1) * (rSges(7,2) * t13 + t114 * t14 + t99) + g(2) * (rSges(7,2) * t11 + t114 * t12 + t100) + g(3) * (-rSges(7,2) * t31 - t114 * t32 + t86) + (t141 * t113 + t149 * (qJ(4) + t141)) * t71) (-m(5) + t147) * (g(1) * t25 + g(2) * t90 + g(3) * t88) t147 * (g(1) * t7 + (t113 + t143) * t70 + (g(2) * (t40 * t73 - t111) + g(3) * t89) * t72) -m(7) * (g(1) * t8 - g(2) * t6 + g(3) * t148)];
taug  = t1(:);
