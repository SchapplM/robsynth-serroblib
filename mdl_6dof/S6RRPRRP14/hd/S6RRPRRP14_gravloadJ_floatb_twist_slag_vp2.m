% Calculate Gravitation load on the joints for
% S6RRPRRP14
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% g [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [10x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,alpha2,d1,d2,d4,d5]';
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
% Datum: 2019-03-09 13:11
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function taug = S6RRPRRP14_gravloadJ_floatb_twist_slag_vp2(qJ, g, ...
  pkin, m, mrSges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(3,1),zeros(10,1),zeros(7,1),zeros(7,3)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRPRRP14_gravloadJ_floatb_twist_slag_vp2: qJ has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6RRPRRP14_gravloadJ_floatb_twist_slag_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S6RRPRRP14_gravloadJ_floatb_twist_slag_vp2: pkin has to be [10x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6RRPRRP14_gravloadJ_floatb_twist_slag_vp2: m has to be [7x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [7,3]), ...
  'S6RRPRRP14_gravloadJ_floatb_twist_slag_vp2: mrSges has to be [7x3] (double)');

%% Symbolic Calculation
% From gravload_joint_floatb_twist_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 13:03:47
% EndTime: 2019-03-09 13:03:49
% DurationCPUTime: 1.01s
% Computational Cost: add. (603->127), mult. (1473->179), div. (0->0), fcn. (1735->10), ass. (0->67)
t132 = m(6) + m(7);
t130 = mrSges(6,3) + mrSges(7,2);
t88 = m(7) * pkin(5) + mrSges(6,1) + mrSges(7,1);
t84 = -m(7) * qJ(6) + mrSges(6,2) - mrSges(7,3);
t131 = mrSges(3,2) - mrSges(4,3);
t62 = sin(qJ(4));
t66 = cos(qJ(4));
t133 = -t62 * mrSges(5,1) - t66 * mrSges(5,2) + t131;
t129 = mrSges(3,1) - mrSges(4,2) + mrSges(5,3);
t61 = sin(qJ(5));
t65 = cos(qJ(5));
t128 = t84 * t61 - t88 * t65 - mrSges(5,1);
t127 = mrSges(5,2) - t130;
t105 = m(5) + t132;
t126 = pkin(9) * t105 + t129;
t63 = sin(qJ(2));
t64 = sin(qJ(1));
t67 = cos(qJ(2));
t106 = cos(pkin(6));
t68 = cos(qJ(1));
t90 = t68 * t106;
t43 = t63 * t90 + t64 * t67;
t60 = sin(pkin(6));
t113 = t60 * t68;
t42 = t63 * t64 - t67 * t90;
t76 = t113 * t66 - t42 * t62;
t125 = t43 * t65 + t61 * t76;
t124 = -t43 * t61 + t65 * t76;
t99 = m(4) + t105;
t123 = qJ(3) * t99 - t131;
t122 = t133 - t132 * (-pkin(10) * t66 + qJ(3)) + t130 * t66 + (-m(4) - m(5)) * qJ(3);
t120 = pkin(4) * t62;
t117 = t60 * t63;
t116 = t60 * t64;
t115 = t60 * t66;
t114 = t60 * t67;
t112 = t61 * t62;
t111 = t62 * t65;
t110 = t63 * t65;
t108 = pkin(2) * t114 + qJ(3) * t117;
t107 = t68 * pkin(1) + pkin(8) * t116;
t104 = t62 * t117;
t103 = t63 * t115;
t91 = t64 * t106;
t45 = -t63 * t91 + t67 * t68;
t101 = t45 * pkin(2) + t107;
t100 = pkin(9) * t114 + t108;
t98 = -pkin(1) * t64 + pkin(8) * t113;
t36 = t42 * pkin(2);
t97 = -pkin(9) * t42 - t36;
t44 = t68 * t63 + t67 * t91;
t38 = t44 * pkin(2);
t96 = -pkin(9) * t44 - t38;
t89 = pkin(3) * t116 + t101;
t87 = -t43 * pkin(2) + t98;
t82 = mrSges(2,2) + (-mrSges(4,1) - mrSges(3,3)) * t60;
t77 = pkin(3) * t113 + t87;
t19 = t113 * t62 + t42 * t66;
t74 = -t132 * pkin(10) + t127;
t41 = t106 * t66 - t114 * t62;
t40 = -t106 * t62 - t114 * t66;
t18 = t115 * t64 + t44 * t62;
t17 = t116 * t62 - t44 * t66;
t15 = -t110 * t60 + t41 * t61;
t2 = t18 * t65 + t45 * t61;
t1 = t18 * t61 - t45 * t65;
t3 = [(-t68 * mrSges(2,1) - m(3) * t107 - m(4) * t101 - m(5) * t89 - t18 * mrSges(5,1) - t123 * t44 - t88 * t2 + t84 * t1 + t82 * t64 - t126 * t45 + t74 * t17 - t132 * (t18 * pkin(4) + t89)) * g(2) + (t64 * mrSges(2,1) - m(3) * t98 - m(4) * t87 - m(5) * t77 - t76 * mrSges(5,1) + t123 * t42 - t88 * t124 + t84 * t125 + t82 * t68 + t126 * t43 + t74 * t19 + t132 * (-pkin(4) * t76 - t77)) * g(1) (-m(4) * t108 - m(5) * t100 - t132 * (pkin(4) * t104 - pkin(10) * t103 + t100) + t84 * (t104 * t61 - t114 * t65) + t130 * t103 + (-t88 * (t110 * t62 + t61 * t67) - t129 * t67 + t133 * t63) * t60) * g(3) + (m(4) * t36 - m(5) * t97 - t132 * (t43 * t120 + t97) - t88 * (t111 * t43 - t42 * t61) + t84 * (t112 * t43 + t42 * t65) + t129 * t42 + t122 * t43) * g(2) + (m(4) * t38 - m(5) * t96 - t132 * (t45 * t120 + t96) + t84 * (t112 * t45 + t44 * t65) - t88 * (t111 * t45 - t44 * t61) + t129 * t44 + t122 * t45) * g(1) (-g(1) * t44 - g(2) * t42 + g(3) * t114) * t99 (-t132 * (t40 * pkin(4) + pkin(10) * t41) + t127 * t41 + t128 * t40) * g(3) + (-t132 * (t19 * pkin(4) - pkin(10) * t76) - t127 * t76 + t128 * t19) * g(2) + (-t132 * (-t17 * pkin(4) + pkin(10) * t18) + t127 * t18 - t128 * t17) * g(1) (t84 * (t117 * t61 + t41 * t65) + t88 * t15) * g(3) + (-t84 * t124 - t125 * t88) * g(2) + (t1 * t88 + t2 * t84) * g(1) (-g(1) * t1 + g(2) * t125 - g(3) * t15) * m(7)];
taug  = t3(:);
