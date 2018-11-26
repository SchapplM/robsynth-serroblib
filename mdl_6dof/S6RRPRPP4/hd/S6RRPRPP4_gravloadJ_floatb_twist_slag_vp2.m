% Calculate Gravitation load on the joints for
% S6RRPRPP4
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% g [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [9x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d2,d4,theta5]';
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
% Datum: 2018-11-23 16:58
% Revision: 76f9d5e39f14dc242b53c0d9d3d9db48bd8f37c0
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für mechatronische Systeme, Universität Hannover

function taug = S6RRPRPP4_gravloadJ_floatb_twist_slag_vp2(qJ, g, ...
  pkin, m, mrSges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(3,1),zeros(9,1),zeros(7,1),zeros(7,3)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRPRPP4_gravloadJ_floatb_twist_slag_vp2: qJ has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6RRPRPP4_gravloadJ_floatb_twist_slag_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S6RRPRPP4_gravloadJ_floatb_twist_slag_vp2: pkin has to be [9x1] (double)');
assert( isreal(m) && all(size(m) == [7 1]), ...
  'S6RRPRPP4_gravloadJ_floatb_twist_slag_vp2: m has to be [7x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [7,3]), ...
  'S6RRPRPP4_gravloadJ_floatb_twist_slag_vp2: mrSges has to be [7x3] (double)');

%% Symbolic Calculation
% From joint_gravload_floatb_twist_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2018-11-23 16:58:31
% EndTime: 2018-11-23 16:58:32
% DurationCPUTime: 0.88s
% Computational Cost: add. (348->117), mult. (585->133), div. (0->0), fcn. (534->8), ass. (0->58)
t89 = m(7) * pkin(5) + mrSges(6,1) + mrSges(7,1);
t88 = m(7) * qJ(6) - mrSges(6,2) + mrSges(7,3);
t94 = m(6) + m(7);
t97 = mrSges(7,2) + mrSges(6,3);
t34 = qJ(4) + pkin(9);
t27 = sin(t34);
t28 = cos(t34);
t36 = sin(qJ(4));
t39 = cos(qJ(4));
t96 = -t36 * mrSges(5,1) - t39 * mrSges(5,2) - t89 * t27 + t88 * t28;
t95 = -m(4) - m(5);
t37 = sin(qJ(2));
t40 = cos(qJ(2));
t93 = (-mrSges(3,1) + mrSges(4,2)) * t40 + (mrSges(3,2) - mrSges(4,3)) * t37;
t38 = sin(qJ(1));
t41 = cos(qJ(1));
t92 = g(1) * t41 + g(2) * t38;
t57 = m(5) * (-pkin(2) - pkin(8)) - mrSges(5,3);
t91 = (-mrSges(4,3) + t96) * t40 + (-mrSges(4,2) - t57 + (m(4) + t94) * pkin(2) + t97) * t37;
t90 = -t97 * t40 + t93;
t87 = -m(5) * pkin(3) - mrSges(4,1) + mrSges(2,2) - mrSges(3,3);
t85 = pkin(4) * t39;
t82 = g(3) * t40;
t31 = t40 * pkin(2);
t35 = -qJ(5) - pkin(8);
t81 = t35 * t40;
t80 = t36 * t41;
t79 = t38 * t36;
t78 = t38 * t37;
t77 = t38 * t39;
t74 = t40 * t41;
t73 = t41 * t37;
t29 = t37 * qJ(3);
t71 = t31 + t29;
t70 = t41 * pkin(1) + t38 * pkin(7);
t69 = qJ(3) * t40;
t68 = pkin(4) * t79;
t23 = t37 * t36 * pkin(4);
t24 = pkin(4) * t80;
t67 = t36 * t73;
t66 = t39 * t73;
t65 = t37 * t77;
t60 = -pkin(1) - t29;
t59 = pkin(2) * t74 + t41 * t29 + t70;
t49 = t60 - t31;
t32 = t41 * pkin(7);
t26 = pkin(3) + t85;
t21 = t41 * t69;
t20 = t38 * t69;
t8 = -t36 * t78 + t39 * t41;
t7 = t65 + t80;
t6 = t67 + t77;
t5 = t66 - t79;
t4 = -t27 * t78 + t28 * t41;
t3 = t27 * t41 + t28 * t78;
t2 = t27 * t73 + t38 * t28;
t1 = t27 * t38 - t28 * t73;
t9 = [(-m(3) * t70 - t6 * mrSges(5,1) - t5 * mrSges(5,2) + t95 * t59 - t94 * (pkin(4) * t67 + t38 * t26 - t35 * t74 + t59) - t89 * t2 - t88 * t1 + (-m(5) * pkin(8) - mrSges(5,3) - t97) * t74 + (-mrSges(2,1) + t93) * t41 + t87 * t38) * g(2) + (-t8 * mrSges(5,1) + t7 * mrSges(5,2) - t94 * (t41 * t26 + t38 * t81 + t32) - t89 * t4 - t88 * t3 + (-m(3) + t95) * t32 + t87 * t41 + (m(3) * pkin(1) - m(4) * t49 - m(5) * t60 - t57 * t40 + mrSges(2,1) - t94 * (t49 - t23) - t90) * t38) * g(1), t92 * (mrSges(3,1) * t37 + mrSges(3,2) * t40) + (-t94 * (t35 * t78 + t40 * t68 + t20) + t95 * t20 + t91 * t38) * g(2) + (-t94 * (t40 * t24 + t35 * t73 + t21) + t95 * t21 + t91 * t41) * g(1) + (-m(4) * t71 - m(5) * (pkin(8) * t40 + t71) - t40 * mrSges(5,3) - t94 * (t23 + t71 - t81) + t96 * t37 + t90) * g(3) (-t37 * t92 + t82) * (t94 - t95) (mrSges(5,1) * t39 - mrSges(5,2) * t36 + t27 * t88 + t28 * t89 + t85 * t94) * t82 + (-t7 * mrSges(5,1) - t8 * mrSges(5,2) - t94 * (pkin(4) * t65 + t24) + t88 * t4 - t89 * t3) * g(2) + (-t5 * mrSges(5,1) + t6 * mrSges(5,2) - t94 * (pkin(4) * t66 - t68) - t88 * t2 + t89 * t1) * g(1) (-g(3) * t37 - t40 * t92) * t94 (-g(1) * t1 + g(2) * t3 - t28 * t82) * m(7)];
taug  = t9(:);
