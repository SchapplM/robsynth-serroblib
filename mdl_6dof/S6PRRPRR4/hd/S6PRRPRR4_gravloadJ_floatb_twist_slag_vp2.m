% Calculate Gravitation load on the joints for
% S6PRRPRR4
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% g [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [11x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,alpha2,d2,d3,d5,d6,theta1]';
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
% Datum: 2019-03-08 22:15
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function taug = S6PRRPRR4_gravloadJ_floatb_twist_slag_vp2(qJ, g, ...
  pkin, m, mrSges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(3,1),zeros(11,1),zeros(7,1),zeros(7,3)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6PRRPRR4_gravloadJ_floatb_twist_slag_vp2: qJ has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6PRRPRR4_gravloadJ_floatb_twist_slag_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [11 1]), ...
  'S6PRRPRR4_gravloadJ_floatb_twist_slag_vp2: pkin has to be [11x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6PRRPRR4_gravloadJ_floatb_twist_slag_vp2: m has to be [7x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [7,3]), ...
  'S6PRRPRR4_gravloadJ_floatb_twist_slag_vp2: mrSges has to be [7x3] (double)');

%% Symbolic Calculation
% From gravload_joint_floatb_twist_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-08 22:10:20
% EndTime: 2019-03-08 22:10:23
% DurationCPUTime: 1.00s
% Computational Cost: add. (601->117), mult. (1542->167), div. (0->0), fcn. (1873->12), ass. (0->64)
t117 = -m(6) - m(7);
t116 = mrSges(4,1) + mrSges(5,1);
t114 = mrSges(4,2) - mrSges(5,3);
t113 = -m(7) * pkin(10) + mrSges(6,2) - mrSges(7,3);
t57 = sin(qJ(6));
t61 = cos(qJ(6));
t115 = m(7) * pkin(5) + t61 * mrSges(7,1) - t57 * mrSges(7,2) + mrSges(6,1);
t56 = cos(pkin(11));
t60 = sin(qJ(2));
t64 = cos(qJ(2));
t89 = sin(pkin(11));
t90 = cos(pkin(6));
t68 = t90 * t89;
t42 = t56 * t64 - t60 * t68;
t59 = sin(qJ(3));
t63 = cos(qJ(3));
t55 = sin(pkin(6));
t83 = t55 * t89;
t23 = t42 * t59 - t63 * t83;
t24 = t42 * t63 + t59 * t83;
t58 = sin(qJ(5));
t62 = cos(qJ(5));
t8 = t23 * t58 + t24 * t62;
t121 = t113 * t8 - t115 * (t23 * t62 - t24 * t58);
t82 = t56 * t90;
t40 = t60 * t82 + t89 * t64;
t96 = t55 * t63;
t21 = t40 * t59 + t56 * t96;
t22 = -t56 * t55 * t59 + t40 * t63;
t4 = t21 * t58 + t22 * t62;
t120 = t113 * t4 - t115 * (t21 * t62 - t22 * t58);
t97 = t55 * t60;
t43 = t59 * t97 - t90 * t63;
t44 = t90 * t59 + t60 * t96;
t16 = t43 * t58 + t44 * t62;
t119 = t113 * t16 - t115 * (t43 * t62 - t44 * t58);
t109 = t114 * t59 - t116 * t63 - mrSges(3,1);
t69 = t58 * t59 + t62 * t63;
t118 = -t69 * t115 + (t58 * t63 - t59 * t62) * t113 + t109;
t39 = -t89 * t60 + t64 * t82;
t91 = qJ(4) * t59;
t99 = t39 * t63;
t112 = pkin(3) * t99 + t39 * t91;
t41 = -t56 * t60 - t64 * t68;
t98 = t41 * t63;
t111 = pkin(3) * t98 + t41 * t91;
t110 = m(5) - t117;
t108 = mrSges(3,2) - mrSges(5,2) - mrSges(4,3);
t101 = t57 * mrSges(7,1) + t61 * mrSges(7,2) + mrSges(6,3) + t108 + t117 * (pkin(8) - pkin(9));
t95 = t55 * t64;
t92 = pkin(2) * t95 + pkin(8) * t97;
t88 = t59 * t95;
t87 = t63 * t95;
t35 = t39 * pkin(2);
t85 = pkin(8) * t40 + t35;
t36 = t41 * pkin(2);
t84 = pkin(8) * t42 + t36;
t81 = -t21 * pkin(3) + qJ(4) * t22;
t80 = -t23 * pkin(3) + qJ(4) * t24;
t79 = -t43 * pkin(3) + qJ(4) * t44;
t76 = pkin(3) * t87 + qJ(4) * t88 + t92;
t65 = pkin(4) * t87 - pkin(9) * t97 + t76;
t26 = t69 * t95;
t1 = [(-m(2) - m(3) - m(4) - t110) * g(3) (-m(4) * t92 - m(5) * t76 - m(6) * t65 - t26 * mrSges(6,1) + mrSges(6,3) * t97 - m(7) * (pkin(5) * t26 + t65) - (t26 * t61 - t57 * t97) * mrSges(7,1) - (-t26 * t57 - t61 * t97) * mrSges(7,2) + t113 * (t58 * t87 - t62 * t88) + (t108 * t60 + t109 * t64) * t55) * g(3) + (-m(4) * t85 - m(5) * (t85 + t112) + t117 * (pkin(4) * t99 + t112 + t35) + t101 * t40 + t118 * t39) * g(2) + (-m(4) * t84 - m(5) * (t84 + t111) + t117 * (pkin(4) * t98 + t111 + t36) + t101 * t42 + t118 * t41) * g(1) (-m(5) * t79 + t117 * (-t43 * pkin(4) + t79) + t114 * t44 + t116 * t43 - t119) * g(3) + (-m(5) * t81 + t117 * (-t21 * pkin(4) + t81) + t114 * t22 + t116 * t21 - t120) * g(2) + (-m(5) * t80 + t117 * (-t23 * pkin(4) + t80) + t114 * t24 + t116 * t23 - t121) * g(1), t110 * (-g(1) * t23 - g(2) * t21 - g(3) * t43) t121 * g(1) + t120 * g(2) + t119 * g(3), -g(1) * ((t41 * t61 - t57 * t8) * mrSges(7,1) + (-t41 * t57 - t61 * t8) * mrSges(7,2)) - g(2) * ((t39 * t61 - t4 * t57) * mrSges(7,1) + (-t39 * t57 - t4 * t61) * mrSges(7,2)) - g(3) * ((-t16 * t57 + t61 * t95) * mrSges(7,1) + (-t16 * t61 - t57 * t95) * mrSges(7,2))];
taug  = t1(:);
