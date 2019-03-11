% Calculate Gravitation load on the joints for
% S6RPRRPP7
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
% Datum: 2019-03-09 04:52
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function taug = S6RPRRPP7_gravloadJ_floatb_twist_slag_vp2(qJ, g, ...
  pkin, m, mrSges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(3,1),zeros(8,1),zeros(7,1),zeros(7,3)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPRRPP7_gravloadJ_floatb_twist_slag_vp2: qJ has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6RPRRPP7_gravloadJ_floatb_twist_slag_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S6RPRRPP7_gravloadJ_floatb_twist_slag_vp2: pkin has to be [8x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6RPRRPP7_gravloadJ_floatb_twist_slag_vp2: m has to be [7x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [7,3]), ...
  'S6RPRRPP7_gravloadJ_floatb_twist_slag_vp2: mrSges has to be [7x3] (double)');

%% Symbolic Calculation
% From gravload_joint_floatb_twist_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 04:50:17
% EndTime: 2019-03-09 04:50:18
% DurationCPUTime: 0.73s
% Computational Cost: add. (247->107), mult. (544->121), div. (0->0), fcn. (515->6), ass. (0->55)
t85 = mrSges(5,1) + mrSges(6,1);
t79 = -mrSges(6,3) - mrSges(7,2);
t75 = mrSges(5,2) + t79;
t80 = mrSges(5,3) + mrSges(6,2);
t27 = cos(qJ(3));
t84 = t80 * t27;
t23 = sin(qJ(4));
t26 = cos(qJ(4));
t83 = -t75 * t23 + t85 * t26;
t82 = m(6) + m(7);
t25 = sin(qJ(1));
t60 = t25 * t26;
t24 = sin(qJ(3));
t28 = cos(qJ(1));
t62 = t24 * t28;
t7 = t23 * t62 + t60;
t55 = t28 * t26;
t61 = t25 * t23;
t8 = t24 * t55 - t61;
t81 = t8 * pkin(4) + t7 * qJ(5);
t66 = g(2) * t28;
t67 = g(1) * t25;
t78 = t66 - t67;
t77 = -m(5) - t82;
t76 = -mrSges(2,1) + mrSges(3,2) - mrSges(4,3);
t41 = m(7) * (-pkin(4) - pkin(5)) - mrSges(7,1);
t45 = -qJ(5) * t23 - pkin(3);
t68 = pkin(4) * t26;
t74 = -m(7) * t45 - t26 * t41 - m(6) * (t45 - t68) + m(5) * pkin(3) + t83;
t73 = m(7) * pkin(5) + mrSges(7,1);
t72 = m(7) * qJ(6) + mrSges(7,3);
t38 = t24 * mrSges(4,1) + t27 * mrSges(4,2);
t40 = m(7) * (-pkin(8) + qJ(6)) + mrSges(7,3);
t71 = -t40 * t27 + mrSges(2,2) - mrSges(3,3) - t38;
t70 = t73 + t85;
t69 = -pkin(1) - pkin(7);
t65 = g(3) * t27;
t63 = t24 * t25;
t59 = t25 * t27;
t56 = t27 * t28;
t53 = pkin(3) * t59 + pkin(8) * t63;
t52 = t28 * pkin(1) + t25 * qJ(2);
t51 = qJ(5) * t27;
t50 = pkin(8) * t59;
t49 = t28 * pkin(7) + t52;
t19 = t28 * qJ(2);
t48 = -t25 * pkin(1) + t19;
t43 = pkin(3) * t63 + t49;
t39 = mrSges(4,1) * t27 - mrSges(4,2) * t24;
t5 = t24 * t61 - t55;
t6 = t23 * t28 + t24 * t60;
t34 = t6 * pkin(4) + qJ(5) * t5 + t43;
t15 = pkin(3) * t62;
t31 = -pkin(8) * t56 + t69 * t25 + t15 + t19;
t1 = [(-m(3) * t52 - m(4) * t49 - m(5) * (t43 - t50) - m(6) * (t34 - t50) - m(7) * t34 - t70 * t6 + t75 * t5 + t76 * t28 + (t71 + t84) * t25) * g(2) + (-m(3) * t48 - m(4) * t19 - m(5) * t31 - m(6) * (t31 + t81) - m(7) * (t15 + t48 + t81) + t80 * t56 - t70 * t8 + t75 * t7 + (-m(4) * t69 + m(7) * pkin(7) - t76) * t25 + t71 * t28) * g(1), t78 * (m(3) + m(4) - t77) -t39 * t67 + (-m(5) * t53 - t82 * (t51 * t61 + t59 * t68 + t53) + ((-t73 * t26 - t83) * t27 + (t72 - t80) * t24) * t25) * g(1) + (t39 + t74 * t27 + (-t40 + (m(5) + m(6)) * pkin(8) + t80) * t24) * t66 + (t74 * t24 + t38 + (pkin(8) * t77 + t72) * t27 - t84) * g(3) -(-mrSges(5,1) * t23 - mrSges(5,2) * t26) * t65 + ((t79 * t26 + (m(6) * pkin(4) + mrSges(6,1) - t41) * t23) * t27 - t82 * t26 * t51) * g(3) + (-t82 * (t7 * pkin(4) - qJ(5) * t8) - t75 * t8 - t70 * t7) * g(2) + (-t82 * (-t5 * pkin(4) + qJ(5) * t6) + t75 * t6 + t70 * t5) * g(1), t82 * (-g(1) * t5 + g(2) * t7 - t23 * t65) (g(3) * t24 + t78 * t27) * m(7)];
taug  = t1(:);
