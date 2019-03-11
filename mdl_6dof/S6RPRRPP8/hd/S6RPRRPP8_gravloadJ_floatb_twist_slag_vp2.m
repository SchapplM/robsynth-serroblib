% Calculate Gravitation load on the joints for
% S6RPRRPP8
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% g [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [8x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d3,d4]';
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
% Datum: 2019-03-09 04:56
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function taug = S6RPRRPP8_gravloadJ_floatb_twist_slag_vp2(qJ, g, ...
  pkin, m, mrSges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(3,1),zeros(8,1),zeros(7,1),zeros(7,3)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPRRPP8_gravloadJ_floatb_twist_slag_vp2: qJ has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6RPRRPP8_gravloadJ_floatb_twist_slag_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S6RPRRPP8_gravloadJ_floatb_twist_slag_vp2: pkin has to be [8x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6RPRRPP8_gravloadJ_floatb_twist_slag_vp2: m has to be [7x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [7,3]), ...
  'S6RPRRPP8_gravloadJ_floatb_twist_slag_vp2: mrSges has to be [7x3] (double)');

%% Symbolic Calculation
% From gravload_joint_floatb_twist_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 04:53:54
% EndTime: 2019-03-09 04:53:56
% DurationCPUTime: 0.76s
% Computational Cost: add. (249->108), mult. (549->124), div. (0->0), fcn. (522->6), ass. (0->54)
t84 = mrSges(5,1) - mrSges(6,2);
t78 = -mrSges(6,3) - mrSges(7,2);
t75 = mrSges(5,2) + t78;
t79 = mrSges(5,3) + mrSges(6,1);
t27 = cos(qJ(3));
t83 = t79 * t27;
t23 = sin(qJ(4));
t26 = cos(qJ(4));
t82 = -t75 * t23 + t84 * t26;
t81 = m(6) + m(7);
t25 = sin(qJ(1));
t61 = t25 * t26;
t24 = sin(qJ(3));
t28 = cos(qJ(1));
t63 = t24 * t28;
t7 = t23 * t63 + t61;
t55 = t28 * t26;
t62 = t25 * t23;
t8 = t24 * t55 - t62;
t80 = t8 * pkin(4) + t7 * qJ(5);
t77 = -m(5) - t81;
t76 = -mrSges(2,1) + mrSges(3,2) - mrSges(4,3);
t40 = m(7) * (-pkin(4) - qJ(6)) - mrSges(7,3);
t45 = -qJ(5) * t23 - pkin(3);
t74 = -m(7) * t45 - t26 * t40 - m(6) * (-pkin(4) * t26 + t45) + m(5) * pkin(3) + t82;
t73 = -m(7) * pkin(5) - mrSges(7,1);
t72 = m(7) * qJ(6) + mrSges(7,3);
t38 = t24 * mrSges(4,1) + t27 * mrSges(4,2);
t41 = m(7) * (-pkin(5) - pkin(8)) - mrSges(7,1);
t71 = -t27 * t41 + mrSges(2,2) - mrSges(3,3) - t38;
t70 = t72 + t84;
t69 = -pkin(1) - pkin(7);
t68 = g(1) * t25;
t67 = g(2) * t28;
t66 = g(3) * t27;
t64 = t24 * t25;
t60 = t25 * t27;
t59 = t26 * t27;
t56 = t27 * t28;
t53 = pkin(3) * t60 + pkin(8) * t64;
t52 = t28 * pkin(1) + t25 * qJ(2);
t51 = qJ(5) * t27;
t50 = pkin(8) * t60;
t49 = t28 * pkin(7) + t52;
t19 = t28 * qJ(2);
t48 = -t25 * pkin(1) + t19;
t43 = pkin(3) * t64 + t49;
t39 = mrSges(4,1) * t27 - mrSges(4,2) * t24;
t5 = t24 * t62 - t55;
t6 = t23 * t28 + t24 * t61;
t34 = t6 * pkin(4) + qJ(5) * t5 + t43;
t15 = pkin(3) * t63;
t31 = -pkin(8) * t56 + t69 * t25 + t15 + t19;
t1 = [(-m(3) * t52 - m(4) * t49 - m(5) * (t43 - t50) - m(6) * (t34 - t50) - m(7) * t34 - t70 * t6 + t75 * t5 + t76 * t28 + (t71 + t83) * t25) * g(2) + (-m(3) * t48 - m(4) * t19 - m(5) * t31 - m(6) * (t31 + t80) - m(7) * (t15 + t48 + t80) + t79 * t56 - t70 * t8 + t75 * t7 + (-m(4) * t69 + m(7) * pkin(7) - t76) * t25 + t71 * t28) * g(1) (-t68 + t67) * (m(3) + m(4) - t77) -t39 * t68 + (-m(5) * t53 - t81 * (t25 * pkin(4) * t59 + t51 * t62 + t53) + ((-t72 * t26 - t82) * t27 + (t73 - t79) * t24) * t25) * g(1) + (t39 + t74 * t27 + (-t41 + (m(5) + m(6)) * pkin(8) + t79) * t24) * t67 + (t74 * t24 + t38 + (pkin(8) * t77 + t73) * t27 - t83) * g(3) -(-mrSges(5,1) * t23 - mrSges(5,2) * t26) * t66 + ((t78 * t26 + (m(6) * pkin(4) - mrSges(6,2) - t40) * t23) * t27 - t81 * t26 * t51) * g(3) + (-t81 * (t7 * pkin(4) - qJ(5) * t8) - t75 * t8 - t70 * t7) * g(2) + (-t81 * (-t5 * pkin(4) + qJ(5) * t6) + t75 * t6 + t70 * t5) * g(1), t81 * (-g(1) * t5 + g(2) * t7 - t23 * t66) (-g(1) * t6 + g(2) * t8 - g(3) * t59) * m(7)];
taug  = t1(:);
