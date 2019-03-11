% Calculate Gravitation load on the joints for
% S6RPRRPR6
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% g [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [11x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d3,d4,d6,theta2,theta5]';
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
% Datum: 2019-03-09 05:18
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function taug = S6RPRRPR6_gravloadJ_floatb_twist_slag_vp2(qJ, g, ...
  pkin, m, mrSges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(3,1),zeros(11,1),zeros(7,1),zeros(7,3)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPRRPR6_gravloadJ_floatb_twist_slag_vp2: qJ has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6RPRRPR6_gravloadJ_floatb_twist_slag_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [11 1]), ...
  'S6RPRRPR6_gravloadJ_floatb_twist_slag_vp2: pkin has to be [11x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6RPRRPR6_gravloadJ_floatb_twist_slag_vp2: m has to be [7x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [7,3]), ...
  'S6RPRRPR6_gravloadJ_floatb_twist_slag_vp2: mrSges has to be [7x3] (double)');

%% Symbolic Calculation
% From gravload_joint_floatb_twist_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 05:15:25
% EndTime: 2019-03-09 05:15:26
% DurationCPUTime: 0.68s
% Computational Cost: add. (471->116), mult. (470->136), div. (0->0), fcn. (421->12), ass. (0->67)
t95 = mrSges(5,3) + mrSges(6,3) + mrSges(7,3);
t34 = qJ(4) + pkin(11);
t29 = cos(t34);
t41 = cos(qJ(4));
t31 = t41 * pkin(4);
t19 = pkin(5) * t29 + t31;
t17 = pkin(3) + t19;
t30 = qJ(6) + t34;
t22 = sin(t30);
t23 = cos(t30);
t25 = t31 + pkin(3);
t27 = sin(t34);
t39 = sin(qJ(4));
t94 = -m(5) * pkin(3) - m(6) * t25 - m(7) * t17 - mrSges(5,1) * t41 - mrSges(6,1) * t29 - mrSges(7,1) * t23 + mrSges(5,2) * t39 + mrSges(6,2) * t27 + mrSges(7,2) * t22;
t37 = -qJ(5) - pkin(8);
t32 = -pkin(9) + t37;
t93 = -m(5) * pkin(8) + m(6) * t37 + m(7) * t32 - t95;
t40 = sin(qJ(1));
t42 = cos(qJ(1));
t92 = g(1) * t42 + g(2) * t40;
t83 = m(6) + m(7);
t88 = t83 + m(4) + m(5);
t84 = m(6) * pkin(4);
t90 = mrSges(5,1) + t84;
t89 = -m(3) * qJ(2) + mrSges(2,2) - mrSges(3,3) - mrSges(4,3);
t80 = pkin(4) * t39;
t18 = pkin(5) * t27 + t80;
t86 = -m(6) * t80 - m(7) * t18;
t33 = pkin(10) + qJ(3);
t26 = sin(t33);
t36 = cos(pkin(10));
t28 = cos(t33);
t54 = t28 * mrSges(4,1) - t26 * mrSges(4,2);
t85 = m(3) * pkin(1) + t36 * mrSges(3,1) - sin(pkin(10)) * mrSges(3,2) + mrSges(2,1) + t54 + t95 * t26;
t69 = t40 * t22;
t5 = t23 * t42 + t28 * t69;
t68 = t40 * t23;
t6 = t22 * t42 - t28 * t68;
t82 = -t5 * mrSges(7,1) + t6 * mrSges(7,2);
t73 = t28 * t42;
t7 = -t22 * t73 + t68;
t8 = t23 * t73 + t69;
t81 = t7 * mrSges(7,1) - t8 * mrSges(7,2);
t77 = g(3) * t26;
t72 = t29 * t42;
t71 = t39 * t42;
t70 = t40 * t18;
t67 = t40 * t27;
t66 = t40 * t29;
t65 = t40 * t39;
t64 = t40 * t41;
t63 = t41 * t42;
t56 = pkin(3) * t28 + pkin(8) * t26;
t52 = -mrSges(7,1) * t22 - mrSges(7,2) * t23;
t51 = t17 * t28 - t26 * t32;
t50 = t25 * t28 - t26 * t37;
t15 = -t28 * t71 + t64;
t13 = t28 * t65 + t63;
t38 = -pkin(7) - qJ(2);
t24 = pkin(2) * t36 + pkin(1);
t16 = t28 * t63 + t65;
t14 = -t28 * t64 + t71;
t12 = t28 * t72 + t67;
t11 = -t27 * t73 + t66;
t10 = t27 * t42 - t28 * t66;
t9 = t28 * t67 + t72;
t1 = [(-t65 * t84 - m(7) * t70 - t16 * mrSges(5,1) - t12 * mrSges(6,1) - t8 * mrSges(7,1) - t15 * mrSges(5,2) - t11 * mrSges(6,2) - t7 * mrSges(7,2) - t88 * (t42 * t24 - t40 * t38) + t89 * t40 + (-m(5) * t56 - m(6) * t50 - m(7) * t51 - t85) * t42) * g(2) + (-t14 * mrSges(5,1) - t10 * mrSges(6,1) - t6 * mrSges(7,1) - t13 * mrSges(5,2) - t9 * mrSges(6,2) - t5 * mrSges(7,2) + (t88 * t38 + t86 + t89) * t42 + (m(4) * t24 - m(5) * (-t24 - t56) - m(6) * (-t24 - t50) - m(7) * (-t24 - t51) + t85) * t40) * g(1) (-g(1) * t40 + g(2) * t42) * (m(3) + t88) -g(3) * t54 + (t94 * g(3) + t92 * (mrSges(4,2) + t93)) * t28 + (t93 * g(3) + t92 * (mrSges(4,1) - t94)) * t26 (mrSges(5,1) * t39 + mrSges(6,1) * t27 + mrSges(5,2) * t41 + mrSges(6,2) * t29 - t52 - t86) * t77 + (-t14 * mrSges(5,2) + t9 * mrSges(6,1) - t10 * mrSges(6,2) - m(7) * (-t19 * t42 - t28 * t70) - t82 + t90 * t13) * g(2) + (t16 * mrSges(5,2) - t11 * mrSges(6,1) + t12 * mrSges(6,2) - m(7) * (-t18 * t73 + t40 * t19) - t81 - t90 * t15) * g(1) (t28 * g(3) - t26 * t92) * t83, -g(1) * t81 - g(2) * t82 - t52 * t77];
taug  = t1(:);
