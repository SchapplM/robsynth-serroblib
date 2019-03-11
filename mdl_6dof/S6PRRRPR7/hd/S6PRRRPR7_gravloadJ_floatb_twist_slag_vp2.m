% Calculate Gravitation load on the joints for
% S6PRRRPR7
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% g [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [13x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,alpha2,alpha3,d2,d3,d4,d6,theta1,theta5]';
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
% Datum: 2019-03-08 23:47
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function taug = S6PRRRPR7_gravloadJ_floatb_twist_slag_vp2(qJ, g, ...
  pkin, m, mrSges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(3,1),zeros(13,1),zeros(7,1),zeros(7,3)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6PRRRPR7_gravloadJ_floatb_twist_slag_vp2: qJ has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6PRRRPR7_gravloadJ_floatb_twist_slag_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [13 1]), ...
  'S6PRRRPR7_gravloadJ_floatb_twist_slag_vp2: pkin has to be [13x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6PRRRPR7_gravloadJ_floatb_twist_slag_vp2: m has to be [7x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [7,3]), ...
  'S6PRRRPR7_gravloadJ_floatb_twist_slag_vp2: mrSges has to be [7x3] (double)');

%% Symbolic Calculation
% From gravload_joint_floatb_twist_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-08 23:39:02
% EndTime: 2019-03-08 23:39:04
% DurationCPUTime: 1.20s
% Computational Cost: add. (1135->124), mult. (2996->189), div. (0->0), fcn. (3780->16), ass. (0->75)
t130 = -m(5) - m(6);
t53 = pkin(13) + qJ(6);
t51 = sin(t53);
t52 = cos(t53);
t54 = sin(pkin(13));
t55 = cos(pkin(13));
t122 = mrSges(5,1) + m(7) * (pkin(5) * t55 + pkin(4)) + t52 * mrSges(7,1) - t51 * mrSges(7,2) + m(6) * pkin(4) + t55 * mrSges(6,1) - t54 * mrSges(6,2);
t117 = -m(6) * qJ(5) + mrSges(5,2) - mrSges(6,3) + m(7) * (-pkin(11) - qJ(5)) - mrSges(7,3);
t129 = m(6) + m(7);
t123 = m(5) + t129;
t57 = sin(qJ(4));
t59 = cos(qJ(4));
t131 = pkin(3) * t123 - t117 * t57 + t122 * t59 + mrSges(4,1);
t116 = -t54 * mrSges(6,1) - t55 * mrSges(6,2) - m(7) * (pkin(5) * t54 + pkin(10)) - t51 * mrSges(7,1) - t52 * mrSges(7,2) + mrSges(4,2) - mrSges(5,3);
t106 = sin(pkin(7));
t128 = -t106 * mrSges(4,3) + mrSges(3,2);
t108 = cos(pkin(12));
t109 = cos(pkin(7));
t105 = sin(pkin(12));
t112 = sin(qJ(2));
t114 = cos(qJ(2));
t110 = cos(pkin(6));
t90 = t110 * t108;
t67 = t105 * t112 - t114 * t90;
t107 = sin(pkin(6));
t86 = t107 * t106;
t126 = t108 * t86 + t67 * t109;
t88 = t110 * t105;
t68 = t108 * t112 + t114 * t88;
t85 = t107 * t105;
t125 = -t106 * t85 + t68 * t109;
t87 = t109 * t107;
t124 = t110 * t106 + t114 * t87;
t118 = t130 * pkin(10) + t116;
t113 = cos(qJ(3));
t74 = t112 * t86;
t94 = t114 * t107;
t111 = pkin(2) * t94 + pkin(9) * t74;
t58 = sin(qJ(3));
t76 = t112 * t87;
t36 = t113 * t94 - t58 * t76;
t104 = t36 * pkin(3) + t111;
t101 = pkin(9) * t106;
t99 = t57 * t106;
t98 = t58 * t109;
t97 = t59 * t106;
t95 = t109 * t113;
t93 = t107 * t112;
t40 = t105 * t114 + t112 * t90;
t84 = -t67 * pkin(2) + t101 * t40;
t41 = t108 * t114 - t112 * t88;
t83 = -t68 * pkin(2) + t101 * t41;
t20 = -t113 * t67 - t40 * t98;
t81 = t20 * pkin(3) + t84;
t22 = -t113 * t68 - t41 * t98;
t80 = t22 * pkin(3) + t83;
t66 = t110 * t109 - t114 * t86;
t61 = t106 * t68 + t109 * t85;
t60 = t106 * t67 - t108 * t87;
t35 = t113 * t76 + t58 * t94;
t29 = t113 * t93 + t124 * t58;
t28 = -t124 * t113 + t58 * t93;
t21 = t41 * t95 - t58 * t68;
t19 = t40 * t95 - t58 * t67;
t16 = t29 * t59 + t57 * t66;
t15 = t29 * t57 - t59 * t66;
t14 = t41 * t113 - t125 * t58;
t13 = t125 * t113 + t41 * t58;
t12 = t113 * t40 - t126 * t58;
t11 = t126 * t113 + t40 * t58;
t4 = t14 * t59 + t57 * t61;
t3 = t14 * t57 - t59 * t61;
t2 = t12 * t59 + t57 * t60;
t1 = t12 * t57 - t59 * t60;
t5 = [(-m(2) - m(3) - m(4) - t123) * g(3) (-m(4) * t111 - m(7) * t104 - mrSges(3,1) * t94 - t36 * mrSges(4,1) + mrSges(3,2) * t93 - mrSges(4,3) * t74 + t130 * (pkin(10) * t35 + t104) - t122 * (t36 * t59 + t57 * t74) + t116 * t35 + t117 * (t36 * t57 - t59 * t74)) * g(3) + (-m(4) * t84 - m(7) * t81 + mrSges(3,1) * t67 - t20 * mrSges(4,1) + t130 * (t19 * pkin(10) + t81) + t128 * t40 - t122 * (t20 * t59 + t40 * t99) + t116 * t19 + t117 * (t20 * t57 - t40 * t97)) * g(2) + (-m(4) * t83 - m(7) * t80 + mrSges(3,1) * t68 - t22 * mrSges(4,1) + t130 * (t21 * pkin(10) + t80) + t128 * t41 - t122 * (t22 * t59 + t41 * t99) + t116 * t21 + t117 * (t22 * t57 - t41 * t97)) * g(1) (t118 * t29 + t131 * t28) * g(3) + (t131 * t11 + t118 * t12) * g(2) + (t118 * t14 + t131 * t13) * g(1) (t117 * t16 + t122 * t15) * g(3) + (t122 * t1 + t117 * t2) * g(2) + (t117 * t4 + t122 * t3) * g(1), t129 * (-g(1) * t3 - g(2) * t1 - g(3) * t15) -g(1) * ((t13 * t52 - t4 * t51) * mrSges(7,1) + (-t13 * t51 - t4 * t52) * mrSges(7,2)) - g(2) * ((t11 * t52 - t2 * t51) * mrSges(7,1) + (-t11 * t51 - t2 * t52) * mrSges(7,2)) - g(3) * ((-t16 * t51 + t28 * t52) * mrSges(7,1) + (-t16 * t52 - t28 * t51) * mrSges(7,2))];
taug  = t5(:);
