% Calculate Gravitation load on the joints for
% S6RRPRRR4
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% g [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [12x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,alpha2,d1,d2,d4,d5,d6,theta3]';
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
% Datum: 2019-03-09 13:37
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function taug = S6RRPRRR4_gravloadJ_floatb_twist_slag_vp1(qJ, g, ...
  pkin, m, rSges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(3,1),zeros(12,1),zeros(7,1),zeros(7,3)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRPRRR4_gravloadJ_floatb_twist_slag_vp1: qJ has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6RRPRRR4_gravloadJ_floatb_twist_slag_vp1: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [12 1]), ...
  'S6RRPRRR4_gravloadJ_floatb_twist_slag_vp1: pkin has to be [12x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6RRPRRR4_gravloadJ_floatb_twist_slag_vp1: m has to be [7x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [7,3]), ...
  'S6RRPRRR4_gravloadJ_floatb_twist_slag_vp1: rSges has to be [7x3] (double)');

%% Symbolic Calculation
% From gravload_joint_floatb_twist_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 13:27:50
% EndTime: 2019-03-09 13:27:54
% DurationCPUTime: 1.41s
% Computational Cost: add. (942->194), mult. (1908->281), div. (0->0), fcn. (2365->14), ass. (0->88)
t142 = pkin(11) + rSges(7,3);
t117 = cos(pkin(12));
t76 = sin(pkin(12));
t81 = sin(qJ(2));
t85 = cos(qJ(2));
t55 = -t81 * t117 - t85 * t76;
t86 = cos(qJ(1));
t123 = t85 * t86;
t82 = sin(qJ(1));
t128 = t81 * t82;
t78 = cos(pkin(6));
t158 = t78 * t123 - t128;
t77 = sin(pkin(6));
t132 = t77 * t82;
t118 = t55 * t78;
t96 = t85 * t117 - t81 * t76;
t34 = t118 * t82 + t86 * t96;
t80 = sin(qJ(4));
t84 = cos(qJ(4));
t19 = t84 * t132 - t34 * t80;
t121 = t86 * t81;
t125 = t82 * t85;
t52 = -t78 * t125 - t121;
t152 = t52 * pkin(2);
t91 = t78 * t96;
t33 = t55 * t86 - t82 * t91;
t71 = pkin(4) * t84 + pkin(3);
t87 = -pkin(10) - pkin(9);
t157 = t33 * t71 - t34 * t87 + t152;
t48 = t55 * t77;
t156 = t48 * t80 + t78 * t84;
t79 = sin(qJ(6));
t83 = cos(qJ(6));
t155 = rSges(7,1) * t83 - rSges(7,2) * t79 + pkin(5);
t131 = t77 * t86;
t29 = t118 * t86 - t82 * t96;
t75 = qJ(4) + qJ(5);
t73 = sin(t75);
t74 = cos(t75);
t14 = -t73 * t131 - t29 * t74;
t30 = t82 * t55 + t86 * t91;
t154 = t14 * t79 + t30 * t83;
t153 = -t14 * t83 + t30 * t79;
t143 = pkin(9) + rSges(5,3);
t47 = t96 * t77;
t151 = g(1) * t33 + g(2) * t30 + g(3) * t47;
t110 = -t74 * t131 + t29 * t73;
t150 = rSges(6,1) * t110 - t14 * rSges(6,2);
t149 = pkin(2) * t85;
t145 = t74 * pkin(5);
t17 = -t132 * t74 + t34 * t73;
t18 = t132 * t73 + t34 * t74;
t141 = -t17 * rSges(6,1) - t18 * rSges(6,2);
t134 = t74 * t79;
t133 = t74 * t83;
t72 = pkin(1) + t149;
t126 = t82 * t72;
t40 = t48 * t73 + t74 * t78;
t41 = -t48 * t74 + t73 * t78;
t119 = t40 * rSges(6,1) - t41 * rSges(6,2);
t116 = t80 * t132;
t67 = t80 * t131;
t68 = t77 * t149;
t112 = t47 * t71 + t48 * t87 + t68;
t49 = pkin(2) * t78 * t81 + (-pkin(8) - qJ(3)) * t77;
t111 = rSges(4,3) * t77 - t49;
t109 = t29 * t84 + t67;
t66 = t86 * t72;
t108 = -t49 * t82 + t66;
t105 = t158 * pkin(2);
t104 = t19 * pkin(4);
t103 = t156 * pkin(4);
t102 = rSges(6,1) * t74 - rSges(6,2) * t73;
t101 = -t86 * t49 - t126;
t99 = t131 * t84 - t29 * t80;
t98 = t29 * t87 + t30 * t71 + t105;
t97 = pkin(4) * t116 + t33 * t87 + t34 * t71 + t108;
t95 = t99 * pkin(4);
t92 = pkin(4) * t67 + t29 * t71 - t30 * t87 + t101;
t90 = t155 * t110 + t14 * t142;
t89 = t142 * t18 - t155 * t17;
t88 = t142 * t41 + t155 * t40;
t53 = -t128 * t78 + t123;
t51 = -t121 * t78 - t125;
t20 = t34 * t84 + t116;
t2 = t18 * t83 - t33 * t79;
t1 = -t18 * t79 - t33 * t83;
t3 = [-m(2) * (g(1) * (-rSges(2,1) * t82 - rSges(2,2) * t86) + g(2) * (rSges(2,1) * t86 - rSges(2,2) * t82)) - m(3) * (g(1) * (rSges(3,1) * t51 - rSges(3,2) * t158 - pkin(1) * t82) + g(2) * (rSges(3,1) * t53 + rSges(3,2) * t52 + pkin(1) * t86) + (g(1) * t86 + g(2) * t82) * t77 * (rSges(3,3) + pkin(8))) - m(4) * (g(1) * (rSges(4,1) * t29 - rSges(4,2) * t30 + t111 * t86 - t126) + g(2) * (rSges(4,1) * t34 + rSges(4,2) * t33 + t111 * t82 + t66)) - m(5) * (g(1) * (rSges(5,1) * t109 + rSges(5,2) * t99 + pkin(3) * t29 + t143 * t30 + t101) + g(2) * (rSges(5,1) * t20 + rSges(5,2) * t19 + pkin(3) * t34 - t143 * t33 + t108)) - m(6) * (g(1) * (-rSges(6,1) * t14 - rSges(6,2) * t110 + rSges(6,3) * t30 + t92) + g(2) * (rSges(6,1) * t18 - rSges(6,2) * t17 - rSges(6,3) * t33 + t97)) - m(7) * (g(1) * (rSges(7,1) * t153 + rSges(7,2) * t154 - t14 * pkin(5) + t142 * t110 + t92) + g(2) * (rSges(7,1) * t2 + rSges(7,2) * t1 + pkin(5) * t18 + t142 * t17 + t97)) -m(3) * (g(1) * (rSges(3,1) * t52 - rSges(3,2) * t53) + g(2) * (rSges(3,1) * t158 + rSges(3,2) * t51) + g(3) * (rSges(3,1) * t85 - rSges(3,2) * t81) * t77) - m(4) * (g(1) * (rSges(4,1) * t33 - rSges(4,2) * t34 + t152) + g(2) * (rSges(4,1) * t30 + rSges(4,2) * t29 + t105) + g(3) * (rSges(4,1) * t47 + rSges(4,2) * t48 + t68)) - m(5) * (g(1) * (t143 * t34 + t152) + g(2) * (-t143 * t29 + t105) + g(3) * (-t143 * t48 + t68) + t151 * (t84 * rSges(5,1) - t80 * rSges(5,2) + pkin(3))) - m(6) * (g(1) * (rSges(6,3) * t34 + t102 * t33 + t157) + g(2) * (-rSges(6,3) * t29 + t102 * t30 + t98) + g(3) * (-rSges(6,3) * t48 + t102 * t47 + t112)) - m(7) * (g(1) * (t33 * t145 + (t133 * t33 + t34 * t79) * rSges(7,1) + (-t134 * t33 + t34 * t83) * rSges(7,2) + t157) + g(2) * (t30 * t145 + (t133 * t30 - t29 * t79) * rSges(7,1) + (-t134 * t30 - t29 * t83) * rSges(7,2) + t98) + g(3) * (t47 * t145 + (t47 * t133 - t48 * t79) * rSges(7,1) + (-t47 * t134 - t48 * t83) * rSges(7,2) + t112) + t151 * t73 * t142) (-m(4) - m(5) - m(6) - m(7)) * (g(3) * t78 + (g(1) * t82 - g(2) * t86) * t77) -m(5) * (g(1) * (rSges(5,1) * t19 - rSges(5,2) * t20) + g(2) * (-rSges(5,1) * t99 + rSges(5,2) * t109) + g(3) * (t156 * rSges(5,1) + (t48 * t84 - t78 * t80) * rSges(5,2))) - m(6) * (g(1) * (t104 + t141) + g(2) * (-t95 + t150) + g(3) * (t103 + t119)) - m(7) * (g(1) * (t104 + t89) + g(2) * (-t95 + t90) + g(3) * (t103 + t88)) -m(6) * (g(1) * t141 + g(2) * t150 + g(3) * t119) - m(7) * (g(1) * t89 + g(2) * t90 + g(3) * t88) -m(7) * (g(1) * (rSges(7,1) * t1 - rSges(7,2) * t2) + g(2) * (-t154 * rSges(7,1) + rSges(7,2) * t153) + g(3) * ((-t41 * t79 - t47 * t83) * rSges(7,1) + (-t41 * t83 + t47 * t79) * rSges(7,2)))];
taug  = t3(:);
