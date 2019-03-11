% Calculate Gravitation load on the joints for
% S6PRRPPR3
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% g [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [10x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,alpha2,d2,d3,d6,theta1]';
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
% Datum: 2019-03-08 21:11
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function taug = S6PRRPPR3_gravloadJ_floatb_twist_slag_vp2(qJ, g, ...
  pkin, m, mrSges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(3,1),zeros(10,1),zeros(7,1),zeros(7,3)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6PRRPPR3_gravloadJ_floatb_twist_slag_vp2: qJ has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6PRRPPR3_gravloadJ_floatb_twist_slag_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S6PRRPPR3_gravloadJ_floatb_twist_slag_vp2: pkin has to be [10x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6PRRPPR3_gravloadJ_floatb_twist_slag_vp2: m has to be [7x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [7,3]), ...
  'S6PRRPPR3_gravloadJ_floatb_twist_slag_vp2: mrSges has to be [7x3] (double)');

%% Symbolic Calculation
% From gravload_joint_floatb_twist_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-08 21:08:10
% EndTime: 2019-03-08 21:08:12
% DurationCPUTime: 0.71s
% Computational Cost: add. (392->102), mult. (979->133), div. (0->0), fcn. (1119->10), ass. (0->56)
t85 = m(6) + m(7);
t76 = m(7) * pkin(9) + mrSges(4,1) + mrSges(5,1) - mrSges(6,2) + mrSges(7,3);
t36 = sin(qJ(6));
t39 = cos(qJ(6));
t87 = -t39 * mrSges(7,1) + t36 * mrSges(7,2) - mrSges(6,1) + mrSges(4,2) - mrSges(5,3);
t86 = t36 * mrSges(7,1) + t39 * mrSges(7,2) + mrSges(3,2) - mrSges(5,2) - mrSges(4,3) + mrSges(6,3);
t38 = sin(qJ(2));
t41 = cos(qJ(2));
t35 = cos(pkin(10));
t67 = cos(pkin(6));
t59 = t35 * t67;
t66 = sin(pkin(10));
t19 = -t38 * t66 + t41 * t59;
t37 = sin(qJ(3));
t68 = qJ(4) * t37;
t40 = cos(qJ(3));
t75 = t19 * t40;
t84 = pkin(3) * t75 + t19 * t68;
t47 = t67 * t66;
t21 = -t35 * t38 - t41 * t47;
t74 = t21 * t40;
t83 = pkin(3) * t74 + t21 * t68;
t82 = m(5) + t85;
t79 = -m(7) * (qJ(4) + pkin(5)) + t87;
t78 = -mrSges(3,1) - t76 * t40 + (-m(7) * pkin(5) + t87) * t37;
t77 = t86 - t85 * (pkin(8) - qJ(5));
t34 = sin(pkin(6));
t73 = t34 * t38;
t72 = t34 * t40;
t71 = t34 * t41;
t69 = pkin(2) * t71 + pkin(8) * t73;
t65 = t40 * t71;
t15 = t19 * pkin(2);
t20 = t38 * t59 + t41 * t66;
t64 = pkin(8) * t20 + t15;
t16 = t21 * pkin(2);
t22 = t35 * t41 - t38 * t47;
t63 = pkin(8) * t22 + t16;
t5 = t20 * t37 + t35 * t72;
t2 = t5 * pkin(3);
t6 = -t34 * t35 * t37 + t20 * t40;
t62 = qJ(4) * t6 - t2;
t60 = t34 * t66;
t7 = t22 * t37 - t40 * t60;
t4 = t7 * pkin(3);
t8 = t22 * t40 + t37 * t60;
t61 = qJ(4) * t8 - t4;
t23 = t37 * t73 - t40 * t67;
t18 = t23 * pkin(3);
t24 = t37 * t67 + t38 * t72;
t57 = qJ(4) * t24 - t18;
t55 = pkin(3) * t65 + t68 * t71 + t69;
t17 = t23 * pkin(4);
t3 = t7 * pkin(4);
t1 = t5 * pkin(4);
t9 = [(-m(2) - m(3) - m(4) - t82) * g(3) (-m(4) * t69 - m(5) * t55 - t85 * (pkin(4) * t65 + t55) + (t78 * t41 + (qJ(5) * t85 + t86) * t38) * t34) * g(3) + (-m(4) * t64 - m(5) * (t64 + t84) - t85 * (pkin(4) * t75 + t15 + t84) + t77 * t20 + t78 * t19) * g(2) + (-m(4) * t63 - m(5) * (t63 + t83) - t85 * (pkin(4) * t74 + t16 + t83) + t77 * t22 + t78 * t21) * g(1) (-m(5) * t57 - m(6) * (-t17 + t57) - m(7) * (-t17 - t18) + t79 * t24 + t76 * t23) * g(3) + (-m(5) * t62 - m(6) * (-t1 + t62) - m(7) * (-t1 - t2) + t79 * t6 + t76 * t5) * g(2) + (-m(5) * t61 - m(6) * (-t3 + t61) - m(7) * (-t3 - t4) + t79 * t8 + t76 * t7) * g(1), t82 * (-g(1) * t7 - g(2) * t5 - g(3) * t23) t85 * (-g(1) * t21 - g(2) * t19 - g(3) * t71) -g(1) * ((t21 * t39 - t36 * t7) * mrSges(7,1) + (-t21 * t36 - t39 * t7) * mrSges(7,2)) - g(2) * ((t19 * t39 - t36 * t5) * mrSges(7,1) + (-t19 * t36 - t39 * t5) * mrSges(7,2)) - g(3) * ((-t23 * t36 + t39 * t71) * mrSges(7,1) + (-t23 * t39 - t36 * t71) * mrSges(7,2))];
taug  = t9(:);
