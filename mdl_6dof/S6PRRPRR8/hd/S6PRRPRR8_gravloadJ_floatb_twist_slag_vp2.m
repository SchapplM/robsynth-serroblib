% Calculate Gravitation load on the joints for
% S6PRRPRR8
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% g [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [12x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,alpha2,alpha3,d2,d3,d5,d6,theta1]';
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
% Datum: 2019-03-08 22:41
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function taug = S6PRRPRR8_gravloadJ_floatb_twist_slag_vp2(qJ, g, ...
  pkin, m, mrSges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(3,1),zeros(12,1),zeros(7,1),zeros(7,3)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6PRRPRR8_gravloadJ_floatb_twist_slag_vp2: qJ has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6PRRPRR8_gravloadJ_floatb_twist_slag_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [12 1]), ...
  'S6PRRPRR8_gravloadJ_floatb_twist_slag_vp2: pkin has to be [12x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6PRRPRR8_gravloadJ_floatb_twist_slag_vp2: m has to be [7x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [7,3]), ...
  'S6PRRPRR8_gravloadJ_floatb_twist_slag_vp2: mrSges has to be [7x3] (double)');

%% Symbolic Calculation
% From gravload_joint_floatb_twist_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-08 22:36:08
% EndTime: 2019-03-08 22:36:11
% DurationCPUTime: 1.13s
% Computational Cost: add. (865->114), mult. (2381->182), div. (0->0), fcn. (2973->14), ass. (0->66)
t113 = m(6) + m(7);
t107 = m(5) + t113;
t72 = -qJ(4) * t107 + mrSges(4,2) - mrSges(5,3);
t57 = sin(qJ(6));
t61 = cos(qJ(6));
t66 = -t57 * mrSges(7,1) - t61 * mrSges(7,2) - t113 * pkin(10) - mrSges(4,1) + mrSges(5,2) - mrSges(6,3);
t67 = -m(7) * pkin(5) - t61 * mrSges(7,1) + t57 * mrSges(7,2) - mrSges(6,1);
t112 = m(7) * pkin(11) - mrSges(6,2) + mrSges(7,3);
t110 = pkin(3) * t107 - t66;
t58 = sin(qJ(5));
t62 = cos(qJ(5));
t103 = t112 * t62 + t67 * t58 + t72;
t101 = cos(qJ(3));
t54 = sin(pkin(7));
t100 = t54 * t58;
t99 = t54 * t62;
t55 = sin(pkin(6));
t56 = cos(pkin(12));
t98 = t55 * t56;
t60 = sin(qJ(2));
t97 = t55 * t60;
t63 = cos(qJ(2));
t96 = t55 * t63;
t95 = -mrSges(5,1) - mrSges(4,3);
t87 = t54 * t97;
t92 = pkin(2) * t96 + pkin(9) * t87;
t91 = cos(pkin(6));
t90 = cos(pkin(7));
t89 = sin(pkin(12));
t59 = sin(qJ(3));
t79 = t59 * t90;
t83 = t55 * t101;
t37 = t63 * t83 - t79 * t97;
t86 = t37 * pkin(3) + t92;
t85 = -m(4) - t107;
t82 = t54 * t91;
t81 = t55 * t89;
t80 = t56 * t91;
t76 = t54 * t81;
t75 = t90 * t101;
t73 = t91 * t89;
t65 = mrSges(3,2) + (-pkin(4) * t113 + t85 * pkin(9) + t95) * t54;
t44 = t56 * t63 - t60 * t73;
t43 = -t56 * t60 - t63 * t73;
t42 = t60 * t80 + t89 * t63;
t41 = -t89 * t60 + t63 * t80;
t40 = -t54 * t96 + t91 * t90;
t39 = t43 * pkin(2);
t38 = t41 * pkin(2);
t36 = (t59 * t63 + t60 * t75) * t55;
t29 = -t43 * t54 + t90 * t81;
t28 = -t41 * t54 - t90 * t98;
t27 = t59 * t82 + (t101 * t60 + t63 * t79) * t55;
t26 = -t101 * t82 + t59 * t97 - t75 * t96;
t22 = t43 * t101 - t44 * t79;
t21 = t43 * t59 + t44 * t75;
t20 = t41 * t101 - t42 * t79;
t19 = t41 * t59 + t42 * t75;
t16 = t26 * t58 + t40 * t62;
t14 = t44 * t101 + (t90 * t43 + t76) * t59;
t13 = -t101 * t76 - t43 * t75 + t44 * t59;
t12 = -t54 * t59 * t98 + t42 * t101 + t41 * t79;
t11 = t56 * t54 * t83 - t41 * t75 + t42 * t59;
t4 = t13 * t58 + t29 * t62;
t2 = t11 * t58 + t28 * t62;
t1 = [(-m(2) - m(3) + t85) * g(3) (-m(4) * t92 - m(5) * t86 + t72 * t36 - t112 * (-t36 * t62 + t58 * t87) + t67 * (t36 * t58 + t62 * t87) + t66 * t37 + (-mrSges(3,1) * t63 + (t95 * t54 + mrSges(3,2)) * t60) * t55 - t113 * (pkin(4) * t87 + t86)) * g(3) + (-t41 * mrSges(3,1) - m(4) * t38 + t112 * (-t42 * t100 + t19 * t62) + t72 * t19 + t67 * (t19 * t58 + t42 * t99) + t66 * t20 + t65 * t42 - t107 * (t20 * pkin(3) + t38)) * g(2) + (-t43 * mrSges(3,1) - m(4) * t39 + t112 * (-t44 * t100 + t21 * t62) + t72 * t21 + t67 * (t21 * t58 + t44 * t99) + t66 * t22 + t65 * t44 - t107 * (t22 * pkin(3) + t39)) * g(1) (t103 * t27 + t110 * t26) * g(3) + (t103 * t12 + t110 * t11) * g(2) + (t103 * t14 + t110 * t13) * g(1), t107 * (-g(1) * t13 - g(2) * t11 - g(3) * t26) (-t112 * t16 + t67 * (t26 * t62 - t40 * t58)) * g(3) + (-t112 * t2 + t67 * (t11 * t62 - t28 * t58)) * g(2) + (-t112 * t4 + t67 * (t13 * t62 - t29 * t58)) * g(1), -g(1) * ((t14 * t61 - t4 * t57) * mrSges(7,1) + (-t14 * t57 - t4 * t61) * mrSges(7,2)) - g(2) * ((t12 * t61 - t2 * t57) * mrSges(7,1) + (-t12 * t57 - t2 * t61) * mrSges(7,2)) - g(3) * ((-t16 * t57 + t27 * t61) * mrSges(7,1) + (-t16 * t61 - t27 * t57) * mrSges(7,2))];
taug  = t1(:);
