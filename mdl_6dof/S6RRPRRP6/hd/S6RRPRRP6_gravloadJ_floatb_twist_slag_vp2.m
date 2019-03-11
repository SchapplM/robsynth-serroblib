% Calculate Gravitation load on the joints for
% S6RRPRRP6
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% g [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [11x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,alpha2,d1,d2,d4,d5,theta3]';
% m_mdh [7x1]
%   mass of all robot links (including the base)
% mrSges [7x3]
%  first moment of all robot links (mass times center of mass in body frames)
%  rows: links of the robot (starting with base)
%  columns: x-, y-, z-coordinates
% 
% Output:
% taug [6x1]
%   joint torques required to compensate gravitation load

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-09 12:15
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function taug = S6RRPRRP6_gravloadJ_floatb_twist_slag_vp2(qJ, g, ...
  pkin, m, mrSges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(3,1),zeros(11,1),zeros(7,1),zeros(7,3)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRPRRP6_gravloadJ_floatb_twist_slag_vp2: qJ has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6RRPRRP6_gravloadJ_floatb_twist_slag_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [11 1]), ...
  'S6RRPRRP6_gravloadJ_floatb_twist_slag_vp2: pkin has to be [11x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6RRPRRP6_gravloadJ_floatb_twist_slag_vp2: m has to be [7x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [7,3]), ...
  'S6RRPRRP6_gravloadJ_floatb_twist_slag_vp2: mrSges has to be [7x3] (double)');

%% Symbolic Calculation
% From gravload_joint_floatb_twist_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 12:07:48
% EndTime: 2019-03-09 12:07:51
% DurationCPUTime: 1.17s
% Computational Cost: add. (989->141), mult. (2462->202), div. (0->0), fcn. (3097->12), ass. (0->75)
t102 = m(7) * pkin(5) + mrSges(6,1) + mrSges(7,1);
t99 = -m(7) * qJ(6) + mrSges(6,2) - mrSges(7,3);
t82 = cos(qJ(2));
t83 = cos(qJ(1));
t116 = t82 * t83;
t78 = sin(qJ(2));
t79 = sin(qJ(1));
t119 = t79 * t78;
t75 = cos(pkin(6));
t139 = t75 * t116 - t119;
t74 = sin(pkin(6));
t123 = t74 * t83;
t112 = cos(pkin(11));
t73 = sin(pkin(11));
t92 = -t112 * t78 - t82 * t73;
t56 = t92 * t75;
t91 = t82 * t112 - t78 * t73;
t40 = -t56 * t83 + t79 * t91;
t77 = sin(qJ(4));
t81 = cos(qJ(4));
t18 = -t77 * t123 + t40 * t81;
t55 = t91 * t75;
t39 = t83 * t55 + t79 * t92;
t76 = sin(qJ(5));
t80 = cos(qJ(5));
t1 = t18 * t76 + t39 * t80;
t138 = t18 * t80 - t39 * t76;
t132 = -m(6) - m(7);
t137 = t81 * mrSges(5,1) - mrSges(5,2) * t77 + mrSges(4,1);
t115 = mrSges(5,3) - mrSges(4,2);
t136 = mrSges(6,3) + mrSges(7,2);
t135 = -t102 * t80 + t99 * t76 - mrSges(5,1);
t134 = mrSges(5,2) - t136;
t111 = m(5) - t132;
t133 = pkin(9) * t111 + t115;
t131 = pkin(2) * t82;
t130 = pkin(4) * t81;
t128 = t39 * t77;
t42 = -t55 * t79 + t83 * t92;
t126 = t42 * t77;
t53 = t91 * t74;
t125 = t53 * t77;
t124 = t74 * t79;
t122 = t76 * t81;
t120 = t78 * t83;
t118 = t79 * t82;
t117 = t80 * t81;
t41 = -t79 * t56 - t83 * t91;
t72 = pkin(1) + t131;
t66 = t83 * t72;
t113 = -pkin(3) * t41 + t66;
t70 = t74 * t131;
t108 = m(4) + t111;
t107 = -m(3) * pkin(1) - mrSges(2,1);
t17 = -t81 * t123 - t40 * t77;
t101 = t139 * pkin(2);
t54 = t92 * t74;
t100 = t53 * pkin(3) - pkin(9) * t54 + t70;
t61 = -t118 * t75 - t120;
t93 = pkin(10) * t132 + t134;
t90 = t61 * pkin(2);
t89 = t39 * pkin(3) + pkin(9) * t40 + t101;
t87 = t42 * pkin(3) - t41 * pkin(9) + t90;
t85 = mrSges(2,2) + (-m(3) * pkin(8) - mrSges(3,3) - mrSges(4,3)) * t74 + t108 * (pkin(2) * t75 * t78 + (-pkin(8) - qJ(3)) * t74);
t62 = -t119 * t75 + t116;
t60 = -t120 * t75 - t118;
t46 = -t54 * t81 + t75 * t77;
t45 = t54 * t77 + t75 * t81;
t35 = t40 * pkin(3);
t22 = t124 * t77 - t41 * t81;
t21 = -t124 * t81 - t41 * t77;
t11 = t46 * t76 + t53 * t80;
t6 = t22 * t80 - t42 * t76;
t5 = t22 * t76 + t42 * t80;
t2 = [(-t62 * mrSges(3,1) - t61 * mrSges(3,2) - m(4) * t66 + t41 * mrSges(4,1) - m(5) * t113 - t22 * mrSges(5,1) + t107 * t83 - t102 * t6 + t99 * t5 + t133 * t42 + t85 * t79 + t93 * t21 + t132 * (t22 * pkin(4) + t113)) * g(2) + (-t60 * mrSges(3,1) + t139 * mrSges(3,2) + t40 * mrSges(4,1) + m(5) * t35 + t18 * mrSges(5,1) + t102 * t138 - t133 * t39 - t99 * t1 + (t108 * t72 - t107) * t79 + t85 * t83 + t93 * t17 + t132 * (-pkin(4) * t18 - t35)) * g(1) (-(mrSges(3,1) * t82 - mrSges(3,2) * t78) * t74 - m(4) * t70 - m(5) * t100 + t132 * (pkin(10) * t125 + t53 * t130 + t100) + t115 * t54 - t137 * t53 - t102 * (t117 * t53 - t54 * t76) + t99 * (t122 * t53 + t54 * t80) - t136 * t125) * g(3) + (-m(4) * t101 - m(5) * t89 - mrSges(3,1) * t139 - mrSges(3,2) * t60 + t132 * (pkin(10) * t128 + t39 * t130 + t89) - t102 * (t117 * t39 + t40 * t76) + t99 * (t122 * t39 - t40 * t80) - t137 * t39 - t115 * t40 - t136 * t128) * g(2) + (-m(4) * t90 - m(5) * t87 - mrSges(3,1) * t61 + mrSges(3,2) * t62 + t132 * (pkin(10) * t126 + t42 * t130 + t87) + t99 * (t122 * t42 + t41 * t80) - t137 * t42 + t115 * t41 - t136 * t126 - t102 * (t117 * t42 - t41 * t76)) * g(1) (-g(3) * t75 + (-g(1) * t79 + g(2) * t83) * t74) * t108 (t132 * (t45 * pkin(4) + pkin(10) * t46) + t134 * t46 + t135 * t45) * g(3) + (t132 * (t17 * pkin(4) + pkin(10) * t18) + t134 * t18 + t135 * t17) * g(2) + (t132 * (-t21 * pkin(4) + pkin(10) * t22) + t134 * t22 - t135 * t21) * g(1) (t99 * (t46 * t80 - t53 * t76) + t102 * t11) * g(3) + (t102 * t1 + t99 * t138) * g(2) + (t102 * t5 + t6 * t99) * g(1) (-g(1) * t5 - g(2) * t1 - g(3) * t11) * m(7)];
taug  = t2(:);
