% Calculate Gravitation load on the joints for
% S6PRRRPR8
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% g [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [12x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,alpha2,alpha3,d2,d3,d4,d6,theta1]';
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
% Datum: 2019-03-08 23:54
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function taug = S6PRRRPR8_gravloadJ_floatb_twist_slag_vp2(qJ, g, ...
  pkin, m, mrSges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(3,1),zeros(12,1),zeros(7,1),zeros(7,3)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6PRRRPR8_gravloadJ_floatb_twist_slag_vp2: qJ has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6PRRRPR8_gravloadJ_floatb_twist_slag_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [12 1]), ...
  'S6PRRRPR8_gravloadJ_floatb_twist_slag_vp2: pkin has to be [12x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6PRRRPR8_gravloadJ_floatb_twist_slag_vp2: m has to be [7x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [7,3]), ...
  'S6PRRRPR8_gravloadJ_floatb_twist_slag_vp2: mrSges has to be [7x3] (double)');

%% Symbolic Calculation
% From gravload_joint_floatb_twist_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-08 23:48:14
% EndTime: 2019-03-08 23:48:18
% DurationCPUTime: 1.27s
% Computational Cost: add. (1052->151), mult. (2903->219), div. (0->0), fcn. (3658->14), ass. (0->94)
t145 = m(6) + m(7);
t149 = mrSges(5,1) - mrSges(6,2) + mrSges(7,3);
t64 = sin(qJ(6));
t67 = cos(qJ(6));
t133 = -t64 * mrSges(7,1) - t67 * mrSges(7,2) + mrSges(5,2) - mrSges(6,3);
t135 = m(7) * pkin(11) + t149;
t146 = pkin(4) * t145 + t135;
t65 = sin(qJ(4));
t116 = qJ(5) * t65;
t125 = cos(qJ(3));
t113 = cos(pkin(12));
t114 = cos(pkin(7));
t110 = sin(pkin(12));
t124 = sin(qJ(2));
t126 = cos(qJ(2));
t115 = cos(pkin(6));
t95 = t115 * t113;
t74 = t110 * t124 - t126 * t95;
t111 = sin(pkin(7));
t112 = sin(pkin(6));
t91 = t112 * t111;
t140 = t113 * t91 + t74 * t114;
t54 = t110 * t126 + t124 * t95;
t66 = sin(qJ(3));
t20 = t140 * t125 + t54 * t66;
t68 = cos(qJ(4));
t123 = t20 * t68;
t143 = -pkin(4) * t123 - t20 * t116;
t93 = t115 * t110;
t75 = t113 * t124 + t126 * t93;
t90 = t112 * t110;
t139 = -t111 * t90 + t75 * t114;
t55 = t113 * t126 - t124 * t93;
t22 = t139 * t125 + t55 * t66;
t122 = t22 * t68;
t142 = -pkin(4) * t122 - t22 * t116;
t92 = t114 * t112;
t138 = t115 * t111 + t126 * t92;
t98 = t112 * t124;
t40 = -t138 * t125 + t66 * t98;
t121 = t40 * t68;
t141 = -pkin(4) * t121 - t40 * t116;
t134 = -qJ(5) * t145 + t133;
t132 = -t133 * t65 + t149 * t68 + mrSges(4,1);
t131 = -m(7) * (pkin(5) + pkin(10)) - t67 * mrSges(7,1) + t64 * mrSges(7,2) - mrSges(6,1) + mrSges(4,2) - mrSges(5,3);
t82 = t124 * t92;
t99 = t126 * t112;
t49 = t125 * t82 + t66 * t99;
t129 = pkin(10) * t49;
t100 = t114 * t125;
t28 = t100 * t54 - t66 * t74;
t128 = t28 * pkin(10);
t30 = t100 * t55 - t66 * t75;
t127 = t30 * pkin(10);
t80 = t124 * t91;
t117 = pkin(2) * t99 + pkin(9) * t80;
t50 = t125 * t99 - t66 * t82;
t109 = t50 * pkin(3) + t117;
t17 = t20 * pkin(3);
t21 = t125 * t54 - t140 * t66;
t108 = t21 * pkin(10) - t17;
t18 = t22 * pkin(3);
t23 = t55 * t125 - t139 * t66;
t107 = t23 * pkin(10) - t18;
t39 = t40 * pkin(3);
t41 = t125 * t98 + t138 * t66;
t106 = t41 * pkin(10) - t39;
t105 = t54 * t111;
t104 = t55 * t111;
t103 = t65 * t111;
t102 = t66 * t114;
t101 = t68 * t111;
t89 = -t74 * pkin(2) + pkin(9) * t105;
t88 = -t75 * pkin(2) + pkin(9) * t104;
t33 = t50 * t65 - t68 * t80;
t34 = t50 * t68 + t65 * t80;
t87 = t34 * pkin(4) + qJ(5) * t33 + t109;
t29 = -t102 * t54 - t125 * t74;
t85 = t29 * pkin(3) + t89;
t31 = -t102 * t55 - t125 * t75;
t84 = t31 * pkin(3) + t88;
t10 = t103 * t54 + t29 * t68;
t9 = -t101 * t54 + t29 * t65;
t77 = t10 * pkin(4) + t9 * qJ(5) + t85;
t11 = -t101 * t55 + t31 * t65;
t12 = t103 * t55 + t31 * t68;
t76 = t12 * pkin(4) + t11 * qJ(5) + t84;
t73 = t114 * t115 - t126 * t91;
t70 = t111 * t75 + t114 * t90;
t69 = t111 * t74 - t113 * t92;
t24 = t41 * t65 - t68 * t73;
t5 = t23 * t65 - t68 * t70;
t3 = t21 * t65 - t68 * t69;
t1 = [(-m(2) - m(3) - m(4) - m(5) - t145) * g(3) (-mrSges(3,1) * t99 + mrSges(3,2) * t98 - m(4) * t117 - t50 * mrSges(4,1) - mrSges(4,3) * t80 - m(5) * (t109 + t129) - m(6) * (t87 + t129) - m(7) * t87 + t133 * t33 + t131 * t49 - t135 * t34) * g(3) + (mrSges(3,1) * t74 + t54 * mrSges(3,2) - m(4) * t89 - t29 * mrSges(4,1) - mrSges(4,3) * t105 - m(5) * (t85 + t128) - m(6) * (t77 + t128) - m(7) * t77 + t133 * t9 + t131 * t28 - t135 * t10) * g(2) + (mrSges(3,1) * t75 + t55 * mrSges(3,2) - m(4) * t88 - t31 * mrSges(4,1) - mrSges(4,3) * t104 - m(5) * (t84 + t127) - m(6) * (t76 + t127) - m(7) * t76 + t133 * t11 + t131 * t30 - t135 * t12) * g(1) (-m(5) * t106 - m(6) * (t106 + t141) - m(7) * (-pkin(11) * t121 + t141 - t39) + t131 * t41 + t132 * t40) * g(3) + (-m(5) * t108 - m(6) * (t108 + t143) - m(7) * (-pkin(11) * t123 + t143 - t17) + t131 * t21 + t132 * t20) * g(2) + (-m(5) * t107 - m(6) * (t107 + t142) - m(7) * (-pkin(11) * t122 + t142 - t18) + t131 * t23 + t132 * t22) * g(1) (t134 * (t41 * t68 + t65 * t73) + t146 * t24) * g(3) + (t134 * (t21 * t68 + t65 * t69) + t146 * t3) * g(2) + (t134 * (t23 * t68 + t65 * t70) + t146 * t5) * g(1), t145 * (-g(1) * t5 - g(2) * t3 - g(3) * t24) -g(1) * ((-t22 * t64 + t5 * t67) * mrSges(7,1) + (-t22 * t67 - t5 * t64) * mrSges(7,2)) - g(2) * ((-t20 * t64 + t3 * t67) * mrSges(7,1) + (-t20 * t67 - t3 * t64) * mrSges(7,2)) - g(3) * ((t24 * t67 - t40 * t64) * mrSges(7,1) + (-t24 * t64 - t40 * t67) * mrSges(7,2))];
taug  = t1(:);
