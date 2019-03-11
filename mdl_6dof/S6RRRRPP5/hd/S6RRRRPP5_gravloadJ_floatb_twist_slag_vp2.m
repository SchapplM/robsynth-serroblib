% Calculate Gravitation load on the joints for
% S6RRRRPP5
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% g [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [9x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d2,d3,d4]';
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
% Datum: 2019-03-09 21:10
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function taug = S6RRRRPP5_gravloadJ_floatb_twist_slag_vp2(qJ, g, ...
  pkin, m, mrSges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(3,1),zeros(9,1),zeros(7,1),zeros(7,3)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRRRPP5_gravloadJ_floatb_twist_slag_vp2: qJ has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6RRRRPP5_gravloadJ_floatb_twist_slag_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S6RRRRPP5_gravloadJ_floatb_twist_slag_vp2: pkin has to be [9x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6RRRRPP5_gravloadJ_floatb_twist_slag_vp2: m has to be [7x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [7,3]), ...
  'S6RRRRPP5_gravloadJ_floatb_twist_slag_vp2: mrSges has to be [7x3] (double)');

%% Symbolic Calculation
% From gravload_joint_floatb_twist_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 21:05:32
% EndTime: 2019-03-09 21:05:35
% DurationCPUTime: 1.26s
% Computational Cost: add. (567->146), mult. (791->162), div. (0->0), fcn. (775->8), ass. (0->75)
t134 = -mrSges(3,2) + mrSges(5,3) + mrSges(6,2);
t131 = -mrSges(7,2) - mrSges(6,3);
t120 = mrSges(5,2) + t131;
t132 = mrSges(5,1) + mrSges(6,1);
t48 = sin(qJ(2));
t51 = cos(qJ(2));
t119 = t51 * mrSges(3,1) + t134 * t48;
t133 = t48 * mrSges(4,3) + mrSges(2,1) + t119;
t129 = mrSges(7,1) + t132;
t46 = qJ(3) + qJ(4);
t41 = sin(t46);
t42 = cos(t46);
t47 = sin(qJ(3));
t50 = cos(qJ(3));
t128 = m(4) * pkin(2) + t50 * mrSges(4,1) - t47 * mrSges(4,2) - t120 * t41 + t132 * t42;
t49 = sin(qJ(1));
t52 = cos(qJ(1));
t90 = t51 * t52;
t25 = -t47 * t90 + t49 * t50;
t126 = g(1) * t52 + g(2) * t49;
t125 = -m(3) - m(4);
t124 = m(6) + m(7);
t123 = mrSges(2,2) - mrSges(3,3);
t87 = qJ(5) * t41;
t91 = t51 * t42;
t122 = pkin(4) * t91 + t51 * t87;
t101 = t42 * t48;
t28 = qJ(5) * t101;
t121 = -m(7) * t28 + t131 * t101;
t118 = m(7) * pkin(5) + mrSges(7,1);
t89 = t52 * t41;
t21 = -t49 * t42 + t51 * t89;
t22 = t41 * t49 + t42 * t90;
t117 = t120 * t22 + t129 * t21;
t92 = t49 * t51;
t19 = t41 * t92 + t42 * t52;
t20 = t49 * t91 - t89;
t116 = t120 * t20 + t129 * t19;
t114 = t118 + t132;
t110 = pkin(3) * t47;
t107 = g(3) * t48;
t102 = t41 * t48;
t100 = t47 * t49;
t99 = t47 * t52;
t53 = -pkin(9) - pkin(8);
t94 = t48 * t53;
t40 = pkin(3) * t50 + pkin(2);
t30 = t51 * t40;
t88 = t52 * pkin(1) + t49 * pkin(7);
t85 = t52 * t94;
t44 = t52 * pkin(7);
t84 = pkin(3) * t99 + t49 * t94 + t44;
t83 = m(4) * pkin(8) + mrSges(4,3);
t80 = t30 - t94;
t79 = -t19 * pkin(4) + qJ(5) * t20;
t78 = -t21 * pkin(4) + qJ(5) * t22;
t77 = -t40 - t87;
t76 = pkin(3) * t100 + t40 * t90 + t88;
t75 = m(7) * (-pkin(4) - pkin(5)) - mrSges(7,1);
t74 = pkin(2) * t51 + pkin(8) * t48;
t73 = m(7) * (-qJ(6) - t53) - mrSges(7,3);
t69 = -mrSges(5,1) * t41 - mrSges(5,2) * t42;
t67 = t25 * pkin(3);
t23 = t47 * t92 + t50 * t52;
t66 = t75 * t41;
t64 = t73 * t48;
t62 = t23 * pkin(3);
t61 = t22 * pkin(4) + t21 * qJ(5) + t76;
t59 = t67 + t78;
t58 = -t62 + t79;
t26 = t50 * t90 + t100;
t24 = -t50 * t92 + t99;
t16 = t21 * pkin(5);
t13 = t19 * pkin(5);
t1 = [(-t26 * mrSges(4,1) - t25 * mrSges(4,2) - m(5) * (t76 - t85) - m(6) * (t61 - t85) - m(7) * t61 + t125 * t88 + t123 * t49 - t114 * t22 + t120 * t21 + (-m(4) * t74 - t133 - t64) * t52) * g(2) + (-m(5) * t84 - t24 * mrSges(4,1) - t23 * mrSges(4,2) - t124 * (-t20 * pkin(4) - qJ(5) * t19 + t84) + t123 * t52 + t125 * t44 + t114 * t20 - t120 * t19 + (m(3) * pkin(1) - m(4) * (-pkin(1) - t74) + (-m(7) * qJ(6) - mrSges(7,3)) * t48 + (-m(5) - t124) * (-pkin(1) - t30) + t133) * t49) * g(1) (-m(5) * t80 - m(6) * (t80 + t122) - m(7) * (t30 + t122) - t64 - t119) * g(3) + ((-t118 * t42 - t128) * g(3) + t126 * (-t73 - t83 + (m(5) + m(6)) * t53 - t134)) * t51 + (-t83 * g(3) + t126 * (mrSges(3,1) + m(5) * t40 - m(6) * (-pkin(4) * t42 + t77) - m(7) * t77 - t42 * t75 + t128)) * t48 (m(5) * t110 + mrSges(4,1) * t47 + mrSges(4,2) * t50 - t69) * t107 + (-m(6) * t28 + (-m(6) * (-pkin(4) * t41 - t110) + t41 * mrSges(6,1) + m(7) * t110 - t66) * t48 + t121) * g(3) + (t23 * mrSges(4,1) - t24 * mrSges(4,2) + m(5) * t62 - m(6) * t58 - m(7) * (-t13 + t58) + t116) * g(2) + (-t25 * mrSges(4,1) + t26 * mrSges(4,2) - m(5) * t67 - m(6) * t59 - m(7) * (-t16 + t59) + t117) * g(1), -t69 * t107 + (-m(6) * (-pkin(4) * t102 + t28) + mrSges(6,1) * t102 - t48 * t66 + t121) * g(3) + (-m(6) * t79 - m(7) * (-t13 + t79) + t116) * g(2) + (-m(6) * t78 - m(7) * (-t16 + t78) + t117) * g(1), t124 * (-g(1) * t21 - g(2) * t19 - g(3) * t102) (-g(3) * t51 + t126 * t48) * m(7)];
taug  = t1(:);
