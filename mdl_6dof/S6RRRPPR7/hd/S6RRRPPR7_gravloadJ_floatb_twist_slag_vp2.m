% Calculate Gravitation load on the joints for
% S6RRRPPR7
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% g [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [10x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d2,d3,d6,theta5]';
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
% Datum: 2019-03-09 16:02
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function taug = S6RRRPPR7_gravloadJ_floatb_twist_slag_vp2(qJ, g, ...
  pkin, m, mrSges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(3,1),zeros(10,1),zeros(7,1),zeros(7,3)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRRPPR7_gravloadJ_floatb_twist_slag_vp2: qJ has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6RRRPPR7_gravloadJ_floatb_twist_slag_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S6RRRPPR7_gravloadJ_floatb_twist_slag_vp2: pkin has to be [10x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6RRRPPR7_gravloadJ_floatb_twist_slag_vp2: m has to be [7x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [7,3]), ...
  'S6RRRPPR7_gravloadJ_floatb_twist_slag_vp2: mrSges has to be [7x3] (double)');

%% Symbolic Calculation
% From gravload_joint_floatb_twist_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 15:56:34
% EndTime: 2019-03-09 15:56:38
% DurationCPUTime: 1.51s
% Computational Cost: add. (450->153), mult. (938->179), div. (0->0), fcn. (986->10), ass. (0->76)
t35 = cos(pkin(10));
t123 = t35 * mrSges(6,1) + mrSges(5,1);
t122 = t35 * mrSges(6,2) + mrSges(5,3);
t119 = -mrSges(5,2) - mrSges(4,3);
t121 = -mrSges(6,3) - mrSges(7,3);
t36 = -pkin(9) - qJ(5);
t100 = -m(7) * t36 - t121;
t120 = t100 + t119;
t110 = m(6) + m(7);
t104 = m(5) + t110;
t106 = -mrSges(4,2) + mrSges(5,3);
t34 = sin(pkin(10));
t111 = -t106 - m(7) * (pkin(5) * t34 + qJ(4));
t37 = sin(qJ(3));
t40 = cos(qJ(3));
t33 = pkin(10) + qJ(6);
t26 = sin(t33);
t27 = cos(t33);
t54 = t26 * t37 + t27 * t40;
t55 = t26 * t40 - t27 * t37;
t88 = t34 * t40;
t89 = t34 * t37;
t118 = t89 * mrSges(6,1) + t54 * mrSges(7,1) - t88 * mrSges(6,2) - t55 * mrSges(7,2) + (mrSges(4,1) + t123) * t40 + (-mrSges(4,2) + t122) * t37;
t42 = cos(qJ(1));
t39 = sin(qJ(1));
t41 = cos(qJ(2));
t84 = t39 * t41;
t12 = t37 * t84 + t40 * t42;
t81 = t42 * t37;
t83 = t40 * t41;
t13 = t39 * t83 - t81;
t109 = t12 * t27 - t13 * t26;
t58 = t12 * t26 + t13 * t27;
t117 = t109 * mrSges(7,1) - t58 * mrSges(7,2);
t116 = m(7) * pkin(5);
t115 = -m(5) - m(6);
t113 = m(6) * qJ(5);
t38 = sin(qJ(2));
t112 = t119 * t38;
t23 = pkin(5) * t35 + pkin(4);
t96 = m(6) * pkin(4) + m(7) * t23 + mrSges(4,1) + mrSges(5,1);
t107 = mrSges(2,2) - mrSges(3,3);
t105 = g(1) * t42 + g(2) * t39;
t62 = t41 * mrSges(3,1) - t38 * mrSges(3,2);
t103 = t38 * t113 - t62;
t78 = qJ(4) * t37;
t49 = -pkin(3) * t40 - pkin(2) - t78;
t50 = pkin(5) * t89 + t23 * t40;
t95 = pkin(4) * t40;
t101 = (t113 + t120) * t41 + (-m(7) * (t49 - t50) - m(6) * (t49 - t95) - m(5) * t49 + m(4) * pkin(2) + t118) * t38;
t99 = pkin(8) * (-m(4) - t104);
t92 = g(3) * t38;
t28 = t38 * pkin(8);
t30 = t41 * pkin(2);
t91 = t12 * t34;
t85 = t38 * t42;
t82 = t41 * t42;
t80 = t30 + t28;
t79 = t42 * pkin(1) + t39 * pkin(7);
t76 = -pkin(1) - t30;
t74 = -t12 * t35 + t13 * t34;
t14 = -t39 * t40 + t41 * t81;
t15 = t39 * t37 + t40 * t82;
t73 = t14 * t34 + t15 * t35;
t68 = pkin(2) * t82 + pkin(8) * t85 + t79;
t66 = t15 * pkin(3) + t68;
t1 = t14 * t27 - t15 * t26;
t2 = t14 * t26 + t15 * t27;
t65 = mrSges(7,1) * t1 - mrSges(7,2) * t2;
t64 = (-mrSges(7,1) * t55 - mrSges(7,2) * t54) * t38;
t57 = t13 * t35 + t91;
t56 = t14 * t35 - t15 * t34;
t31 = t42 * pkin(7);
t10 = t14 * pkin(3);
t8 = t12 * pkin(3);
t3 = [(-m(3) * t79 - m(4) * t68 - m(7) * t66 - t73 * mrSges(6,1) - t2 * mrSges(7,1) - t56 * mrSges(6,2) - t1 * mrSges(7,2) + t115 * (t14 * qJ(4) + t66) + (-mrSges(2,1) + t103) * t42 + t107 * t39 - t96 * t15 + t111 * t14 + t120 * t85) * g(2) + (t91 * t116 + t57 * mrSges(6,1) + t58 * mrSges(7,1) - t74 * mrSges(6,2) + t109 * mrSges(7,2) - t104 * (-t13 * pkin(3) - qJ(4) * t12 + t31) + t107 * t42 + (-m(3) - m(4)) * t31 + t96 * t13 + t106 * t12 + (m(3) * pkin(1) + mrSges(2,1) + t62 - t110 * t76 + (-m(4) - m(5)) * (t76 - t28) + (-m(6) * (-pkin(8) + qJ(5)) - m(7) * (-pkin(8) - t36) + t121) * t38 - t112) * t39) * g(1), t105 * (mrSges(3,1) * t38 + mrSges(3,2) * t41) + (t101 * t39 + t84 * t99) * g(2) + (t101 * t42 + t82 * t99) * g(1) + (-m(4) * t80 - t104 * (pkin(3) * t83 + t41 * t78 + t80) + t100 * t38 + (-m(6) * t95 - m(7) * t50 - t118) * t41 + t103 + t112) * g(3) -(-mrSges(4,1) * t37 - mrSges(4,2) * t40) * t92 + (t64 + (-t88 * t116 + (m(5) * pkin(3) - m(6) * (-pkin(3) - pkin(4)) - t34 * mrSges(6,2) - m(7) * (-pkin(3) - t23) + t123) * t37 + (-t34 * mrSges(6,1) - t104 * qJ(4) - t122) * t40) * t38) * g(3) + (m(7) * t8 - t74 * mrSges(6,1) - t57 * mrSges(6,2) + t115 * (qJ(4) * t13 - t8) + t111 * t13 + t96 * t12 + t117) * g(2) + (m(7) * t10 + t56 * mrSges(6,1) - t73 * mrSges(6,2) + t65 + t115 * (qJ(4) * t15 - t10) + t111 * t15 + t96 * t14) * g(1), t104 * (-g(1) * t14 - g(2) * t12 - t37 * t92) (-g(3) * t41 + t105 * t38) * t110, -g(1) * t65 - g(2) * t117 - g(3) * t64];
taug  = t3(:);
