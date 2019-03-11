% Calculate Gravitation load on the joints for
% S6RRPRPR4
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% g [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [12x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,alpha2,d1,d2,d4,d6,theta3,theta5]';
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
% Datum: 2019-03-09 10:27
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function taug = S6RRPRPR4_gravloadJ_floatb_twist_slag_vp2(qJ, g, ...
  pkin, m, mrSges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(3,1),zeros(12,1),zeros(7,1),zeros(7,3)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRPRPR4_gravloadJ_floatb_twist_slag_vp2: qJ has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6RRPRPR4_gravloadJ_floatb_twist_slag_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [12 1]), ...
  'S6RRPRPR4_gravloadJ_floatb_twist_slag_vp2: pkin has to be [12x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6RRPRPR4_gravloadJ_floatb_twist_slag_vp2: m has to be [7x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [7,3]), ...
  'S6RRPRPR4_gravloadJ_floatb_twist_slag_vp2: mrSges has to be [7x3] (double)');

%% Symbolic Calculation
% From gravload_joint_floatb_twist_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 10:21:25
% EndTime: 2019-03-09 10:21:28
% DurationCPUTime: 1.22s
% Computational Cost: add. (828->137), mult. (1750->191), div. (0->0), fcn. (2132->14), ass. (0->69)
t63 = sin(qJ(6));
t67 = cos(qJ(6));
t126 = m(7) * pkin(5) + mrSges(7,1) * t67 - mrSges(7,2) * t63 + mrSges(6,1);
t125 = -m(7) * pkin(10) + mrSges(6,2) - mrSges(7,3);
t128 = m(6) + m(7);
t129 = -m(4) - m(5);
t97 = t128 - t129;
t105 = sin(pkin(11));
t106 = cos(pkin(11));
t65 = sin(qJ(2));
t69 = cos(qJ(2));
t74 = t69 * t105 + t65 * t106;
t70 = cos(qJ(1));
t109 = t69 * t70;
t66 = sin(qJ(1));
t111 = t66 * t65;
t61 = cos(pkin(6));
t131 = t61 * t109 - t111;
t60 = sin(pkin(6));
t116 = t60 * t66;
t107 = t74 * t61;
t38 = t65 * t105 - t69 * t106;
t24 = -t107 * t66 - t70 * t38;
t64 = sin(qJ(4));
t68 = cos(qJ(4));
t9 = t68 * t116 - t24 * t64;
t32 = t74 * t60;
t130 = -t32 * t64 + t61 * t68;
t127 = m(5) * pkin(3) + t68 * mrSges(5,1) - t64 * mrSges(5,2) + mrSges(4,1);
t82 = m(5) * pkin(9) - mrSges(4,2) + mrSges(5,3) + mrSges(6,3);
t124 = -mrSges(7,1) * t63 - mrSges(7,2) * t67 - t82;
t110 = t66 * t69;
t113 = t65 * t70;
t36 = -t61 * t110 - t113;
t115 = t60 * t70;
t19 = -t107 * t70 + t66 * t38;
t81 = -t68 * t115 + t19 * t64;
t59 = qJ(4) + pkin(12);
t57 = sin(t59);
t58 = cos(t59);
t123 = t125 * t57 - t126 * t58 - t127;
t120 = pkin(2) * t69;
t52 = t60 * t120;
t104 = t64 * t116;
t51 = t64 * t115;
t96 = -m(3) * pkin(1) - mrSges(2,1);
t4 = -t57 * t115 - t19 * t58;
t3 = -t58 * t115 + t19 * t57;
t72 = t61 * t38;
t23 = t66 * t72 - t70 * t74;
t56 = pkin(1) + t120;
t50 = t70 * t56;
t55 = pkin(4) * t68 + pkin(3);
t62 = -qJ(5) - pkin(9);
t90 = pkin(4) * t104 + t23 * t62 + t24 * t55 + t50;
t89 = -m(3) * pkin(8) - mrSges(3,3) - mrSges(4,3);
t87 = t131 * pkin(2);
t79 = (pkin(2) * t61 * t65 + (-pkin(8) - qJ(3)) * t60) * t97 + mrSges(2,2);
t37 = -t111 * t61 + t109;
t35 = -t113 * t61 - t110;
t31 = t38 * t60;
t26 = t32 * t58 + t57 * t61;
t20 = -t66 * t74 - t70 * t72;
t10 = t24 * t68 + t104;
t8 = t116 * t57 + t24 * t58;
t7 = -t116 * t58 + t24 * t57;
t2 = -t23 * t63 + t67 * t8;
t1 = -t23 * t67 - t63 * t8;
t5 = [(-t37 * mrSges(3,1) - t36 * mrSges(3,2) - m(4) * t50 - t24 * mrSges(4,1) - m(5) * (pkin(3) * t24 + t50) - t10 * mrSges(5,1) - t9 * mrSges(5,2) - m(6) * t90 - t8 * mrSges(6,1) - m(7) * (pkin(5) * t8 + t90) - t2 * mrSges(7,1) - t1 * mrSges(7,2) + t96 * t70 + t125 * t7 + (t60 * t89 + t79) * t66 + t82 * t23) * g(2) + (-t35 * mrSges(3,1) + t131 * mrSges(3,2) - t51 * mrSges(5,1) + t125 * t3 - t127 * t19 + (t56 * t97 - t96) * t66 + ((-mrSges(5,2) * t68 + t89) * t60 + t79) * t70 + t126 * t4 + t124 * t20 + t128 * (-pkin(4) * t51 - t19 * t55 + t20 * t62)) * g(1) (-(mrSges(3,1) * t69 - mrSges(3,2) * t65) * t60 - t128 * (-t31 * t55 - t32 * t62 + t52) + t129 * t52 + t124 * t32 - t123 * t31) * g(3) + (-mrSges(3,1) * t131 - mrSges(3,2) * t35 + t129 * t87 - t128 * (t19 * t62 + t20 * t55 + t87) + t123 * t20 - t124 * t19) * g(2) + (mrSges(3,2) * t37 - t128 * (t23 * t55 - t24 * t62) + t123 * t23 + t124 * t24 + (-t97 * pkin(2) - mrSges(3,1)) * t36) * g(1) (-t61 * g(3) + (-g(1) * t66 + g(2) * t70) * t60) * t97 (-t130 * mrSges(5,1) - (-t32 * t68 - t61 * t64) * mrSges(5,2) + t125 * t26 - t126 * (-t32 * t57 + t58 * t61)) * g(3) + (-t81 * mrSges(5,1) - (t19 * t68 + t51) * mrSges(5,2) - t126 * t3 + t125 * t4) * g(2) + (-t9 * mrSges(5,1) + t10 * mrSges(5,2) + t125 * t8 + t126 * t7) * g(1) + (-g(1) * t9 - g(2) * t81 - g(3) * t130) * t128 * pkin(4), t128 * (g(1) * t23 + g(2) * t20 - g(3) * t31) -g(1) * (mrSges(7,1) * t1 - mrSges(7,2) * t2) - g(2) * ((-t20 * t67 - t4 * t63) * mrSges(7,1) + (t20 * t63 - t4 * t67) * mrSges(7,2)) - g(3) * ((-t26 * t63 + t31 * t67) * mrSges(7,1) + (-t26 * t67 - t31 * t63) * mrSges(7,2))];
taug  = t5(:);
