% Calculate Gravitation load on the joints for
% S6RRPRPR6
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% g [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [11x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,alpha2,d1,d2,d4,d6,theta3]';
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
% Datum: 2019-03-09 10:44
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function taug = S6RRPRPR6_gravloadJ_floatb_twist_slag_vp2(qJ, g, ...
  pkin, m, mrSges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(3,1),zeros(11,1),zeros(7,1),zeros(7,3)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRPRPR6_gravloadJ_floatb_twist_slag_vp2: qJ has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6RRPRPR6_gravloadJ_floatb_twist_slag_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [11 1]), ...
  'S6RRPRPR6_gravloadJ_floatb_twist_slag_vp2: pkin has to be [11x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6RRPRPR6_gravloadJ_floatb_twist_slag_vp2: m has to be [7x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [7,3]), ...
  'S6RRPRPR6_gravloadJ_floatb_twist_slag_vp2: mrSges has to be [7x3] (double)');

%% Symbolic Calculation
% From gravload_joint_floatb_twist_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 10:38:11
% EndTime: 2019-03-09 10:38:14
% DurationCPUTime: 1.16s
% Computational Cost: add. (798->132), mult. (1960->178), div. (0->0), fcn. (2422->12), ass. (0->78)
t127 = mrSges(5,2) - mrSges(6,3);
t128 = m(6) + m(7);
t74 = -qJ(5) * t128 + t127;
t57 = sin(qJ(6));
t61 = cos(qJ(6));
t79 = -t57 * mrSges(7,1) - t61 * mrSges(7,2);
t118 = t74 + t79;
t133 = mrSges(5,1) - mrSges(6,2) + mrSges(7,3);
t59 = sin(qJ(2));
t63 = cos(qJ(2));
t92 = sin(pkin(11));
t93 = cos(pkin(11));
t69 = t59 * t93 + t63 * t92;
t60 = sin(qJ(1));
t102 = t60 * t59;
t56 = cos(pkin(6));
t64 = cos(qJ(1));
t99 = t63 * t64;
t130 = t56 * t99 - t102;
t119 = m(7) * pkin(10) + t133;
t129 = pkin(4) * t128 + t119;
t42 = t59 * t92 - t63 * t93;
t67 = t56 * t42;
t22 = -t60 * t69 - t64 * t67;
t62 = cos(qJ(4));
t111 = t22 * t62;
t58 = sin(qJ(4));
t94 = qJ(5) * t58;
t126 = pkin(4) * t111 + t22 * t94;
t25 = t60 * t67 - t64 * t69;
t110 = t25 * t62;
t125 = pkin(4) * t110 + t25 * t94;
t55 = sin(pkin(6));
t35 = t42 * t55;
t109 = t35 * t62;
t124 = -pkin(4) * t109 - t35 * t94;
t123 = mrSges(6,1) + mrSges(5,3) - mrSges(4,2);
t91 = m(5) + t128;
t121 = m(7) * pkin(5) + pkin(9) * t91 + t123;
t96 = t69 * t56;
t26 = -t64 * t42 - t60 * t96;
t21 = t60 * t42 - t64 * t96;
t120 = -t61 * mrSges(7,1) + t57 * mrSges(7,2);
t117 = mrSges(4,1) + t133 * t62 + (-t79 - t127) * t58;
t116 = m(7) * (pkin(5) + pkin(9)) - t120 + t123;
t113 = pkin(2) * t63;
t108 = t55 * t60;
t107 = t55 * t64;
t104 = t59 * t64;
t101 = t60 * t63;
t54 = pkin(1) + t113;
t48 = t64 * t54;
t97 = t26 * pkin(3) + t48;
t51 = t55 * t113;
t95 = -t35 * pkin(3) + t51;
t88 = m(4) + t91;
t87 = -m(3) * pkin(1) - mrSges(2,1);
t8 = -t58 * t107 - t21 * t62;
t82 = t130 * pkin(2);
t36 = t69 * t55;
t81 = pkin(9) * t36 + t95;
t75 = t22 * pkin(3) + t82;
t7 = t107 * t62 - t21 * t58;
t40 = -t101 * t56 - t104;
t72 = t40 * pkin(2);
t70 = -pkin(9) * t21 + t75;
t68 = t25 * pkin(3) + t72;
t66 = pkin(9) * t26 + t68;
t65 = mrSges(2,2) + (-m(3) * pkin(8) - mrSges(3,3) - mrSges(4,3)) * t55 + t88 * (pkin(2) * t56 * t59 + (-pkin(8) - qJ(3)) * t55);
t41 = -t102 * t56 + t99;
t39 = -t104 * t56 - t101;
t28 = t36 * t58 - t56 * t62;
t18 = t21 * pkin(3);
t12 = t108 * t58 + t26 * t62;
t11 = -t108 * t62 + t26 * t58;
t2 = t11 * t57 - t25 * t61;
t1 = t11 * t61 + t25 * t57;
t3 = [(-t41 * mrSges(3,1) - t40 * mrSges(3,2) - m(4) * t48 - t26 * mrSges(4,1) - m(5) * t97 - t2 * mrSges(7,1) - t1 * mrSges(7,2) + t87 * t64 + t74 * t11 + t65 * t60 - t119 * t12 + t121 * t25 - t128 * (t12 * pkin(4) + t97)) * g(2) + (-t39 * mrSges(3,1) + t130 * mrSges(3,2) - t21 * mrSges(4,1) - m(5) * t18 + (t54 * t88 - t87) * t60 + t65 * t64 + t119 * t8 - t118 * t7 + (t120 - t121) * t22 + t128 * (pkin(4) * t8 - t18)) * g(1) (-(mrSges(3,1) * t63 - mrSges(3,2) * t59) * t55 - m(4) * t51 - m(5) * t81 - m(6) * (t81 + t124) - m(7) * (-pkin(10) * t109 + t124 + t95) - t116 * t36 + t117 * t35) * g(3) + (-mrSges(3,1) * t130 - mrSges(3,2) * t39 - m(4) * t82 - m(5) * t70 - m(6) * (t70 + t126) - m(7) * (pkin(10) * t111 + t126 + t75) + t116 * t21 - t117 * t22) * g(2) + (-mrSges(3,1) * t40 + mrSges(3,2) * t41 - m(4) * t72 - m(5) * t66 - m(6) * (t66 + t125) - m(7) * (pkin(10) * t110 + t125 + t68) - t116 * t26 - t117 * t25) * g(1) (-t56 * g(3) + (-g(1) * t60 + g(2) * t64) * t55) * t88 (t118 * (t36 * t62 + t56 * t58) + t129 * t28) * g(3) + (t118 * t8 + t129 * t7) * g(2) + (t11 * t129 + t118 * t12) * g(1), t128 * (-g(1) * t11 - g(2) * t7 - g(3) * t28) -g(1) * (mrSges(7,1) * t1 - mrSges(7,2) * t2) - g(2) * ((t22 * t57 + t61 * t7) * mrSges(7,1) + (t22 * t61 - t57 * t7) * mrSges(7,2)) - g(3) * ((t28 * t61 - t35 * t57) * mrSges(7,1) + (-t28 * t57 - t35 * t61) * mrSges(7,2))];
taug  = t3(:);
