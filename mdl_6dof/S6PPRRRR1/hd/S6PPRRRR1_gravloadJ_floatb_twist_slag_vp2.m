% Calculate Gravitation load on the joints for
% S6PPRRRR1
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% g [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [13x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,alpha2,alpha3,d3,d4,d5,d6,theta1,theta2]';
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

% Quelle: HybrDyn-Toolbox (ehem. IRT-Maple-Toolbox)
% Datum: 2018-11-23 14:52
% Revision: 76f9d5e39f14dc242b53c0d9d3d9db48bd8f37c0
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für mechatronische Systeme, Universität Hannover

function taug = S6PPRRRR1_gravloadJ_floatb_twist_slag_vp2(qJ, g, ...
  pkin, m, mrSges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(3,1),zeros(13,1),zeros(7,1),zeros(7,3)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6PPRRRR1_gravloadJ_floatb_twist_slag_vp2: qJ has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6PPRRRR1_gravloadJ_floatb_twist_slag_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [13 1]), ...
  'S6PPRRRR1_gravloadJ_floatb_twist_slag_vp2: pkin has to be [13x1] (double)');
assert( isreal(m) && all(size(m) == [7 1]), ...
  'S6PPRRRR1_gravloadJ_floatb_twist_slag_vp2: m has to be [7x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [7,3]), ...
  'S6PPRRRR1_gravloadJ_floatb_twist_slag_vp2: mrSges has to be [7x3] (double)');

%% Symbolic Calculation
% From joint_gravload_floatb_twist_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2018-11-23 14:51:50
% EndTime: 2018-11-23 14:51:51
% DurationCPUTime: 0.77s
% Computational Cost: add. (2591->116), mult. (2615->165), div. (0->0), fcn. (2566->24), ass. (0->80)
t143 = mrSges(6,2) - mrSges(7,3);
t68 = sin(qJ(6));
t71 = cos(qJ(6));
t138 = -mrSges(7,1) * t71 + mrSges(7,2) * t68 - mrSges(6,1);
t118 = cos(pkin(6));
t114 = pkin(7) - qJ(3);
t103 = sin(t114);
t113 = pkin(7) + qJ(3);
t96 = sin(t113) / 0.2e1;
t133 = t96 - t103 / 0.2e1;
t111 = pkin(6) + pkin(13);
t102 = cos(t111);
t112 = pkin(6) - pkin(13);
t93 = cos(t112) / 0.2e1;
t56 = t93 - t102 / 0.2e1;
t104 = cos(t114);
t97 = cos(t113) / 0.2e1;
t57 = t97 - t104 / 0.2e1;
t73 = cos(qJ(3));
t101 = sin(t112);
t92 = sin(t111) / 0.2e1;
t79 = t92 + t101 / 0.2e1;
t135 = -t118 * t57 + t133 * t79 + t56 * t73;
t64 = sin(pkin(7));
t67 = cos(pkin(7));
t50 = t118 * t67 - t79 * t64;
t69 = sin(qJ(4));
t72 = cos(qJ(4));
t142 = -t135 * t69 + t50 * t72;
t116 = sin(pkin(12));
t65 = sin(pkin(6));
t105 = t65 * t116;
t117 = cos(pkin(12));
t55 = t92 - t101 / 0.2e1;
t66 = cos(pkin(13));
t49 = -t116 * t55 + t117 * t66;
t115 = sin(pkin(13));
t80 = t93 + t102 / 0.2e1;
t76 = t117 * t115 + t116 * t80;
t136 = -t57 * t105 - t133 * t76 + t49 * t73;
t40 = t67 * t105 + t76 * t64;
t141 = -t136 * t69 + t40 * t72;
t106 = t65 * t117;
t48 = t116 * t66 + t117 * t55;
t75 = t116 * t115 - t117 * t80;
t137 = t57 * t106 - t133 * t75 + t48 * t73;
t39 = -t67 * t106 + t75 * t64;
t140 = -t137 * t69 + t39 * t72;
t134 = -m(6) - m(7);
t63 = qJ(4) + qJ(5);
t61 = sin(t63);
t62 = cos(t63);
t132 = m(5) * pkin(3) + t72 * mrSges(5,1) - t69 * mrSges(5,2) + mrSges(4,1) + (m(7) * pkin(5) - t138) * t62 + (m(7) * pkin(11) - t143) * t61;
t131 = -m(5) * pkin(9) - t68 * mrSges(7,1) - t71 * mrSges(7,2) + mrSges(4,2) - mrSges(5,3) - mrSges(6,3);
t130 = m(3) + m(4) + m(5) - t134;
t11 = -t137 * t61 + t39 * t62;
t12 = t137 * t62 + t39 * t61;
t109 = t11 * pkin(5) + pkin(11) * t12;
t13 = -t136 * t61 + t40 * t62;
t14 = t136 * t62 + t40 * t61;
t108 = t13 * pkin(5) + pkin(11) * t14;
t20 = -t135 * t61 + t50 * t62;
t21 = t135 * t62 + t50 * t61;
t107 = t20 * pkin(5) + pkin(11) * t21;
t100 = t140 * pkin(4);
t99 = t141 * pkin(4);
t98 = t142 * pkin(4);
t90 = t138 * t11 + t143 * t12;
t89 = t138 * t13 + t143 * t14;
t86 = t138 * t20 + t143 * t21;
t83 = t104 / 0.2e1 + t97;
t81 = t96 + t103 / 0.2e1;
t78 = t65 * t81;
t74 = -pkin(10) - pkin(9);
t70 = sin(qJ(3));
t60 = pkin(4) * t72 + pkin(3);
t34 = -t118 * t81 + t56 * t70 - t79 * t83;
t29 = -t116 * t78 + t49 * t70 + t76 * t83;
t26 = t117 * t78 + t48 * t70 + t75 * t83;
t1 = [(-m(2) - t130) * g(3) (-t118 * g(3) + (-g(1) * t116 + g(2) * t117) * t65) * t130 (t134 * (-t135 * t74 - t34 * t60) + t131 * t135 + t132 * t34) * g(3) + (t134 * (-t137 * t74 - t26 * t60) + t131 * t137 + t132 * t26) * g(2) + (t134 * (-t136 * t74 - t29 * t60) + t131 * t136 + t132 * t29) * g(1) (-t142 * mrSges(5,1) - (-t135 * t72 - t50 * t69) * mrSges(5,2) - m(6) * t98 - m(7) * (t107 + t98) + t86) * g(3) + (-t140 * mrSges(5,1) - (-t137 * t72 - t39 * t69) * mrSges(5,2) - m(6) * t100 - m(7) * (t100 + t109) + t90) * g(2) + (-t141 * mrSges(5,1) - (-t136 * t72 - t40 * t69) * mrSges(5,2) - m(6) * t99 - m(7) * (t108 + t99) + t89) * g(1) (-m(7) * t107 + t86) * g(3) + (-m(7) * t109 + t90) * g(2) + (-m(7) * t108 + t89) * g(1), -g(1) * ((-t14 * t68 + t29 * t71) * mrSges(7,1) + (-t14 * t71 - t29 * t68) * mrSges(7,2)) - g(2) * ((-t12 * t68 + t26 * t71) * mrSges(7,1) + (-t12 * t71 - t26 * t68) * mrSges(7,2)) - g(3) * ((-t21 * t68 + t34 * t71) * mrSges(7,1) + (-t21 * t71 - t34 * t68) * mrSges(7,2))];
taug  = t1(:);
