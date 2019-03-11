% Calculate Gravitation load on the joints for
% S6RPRRPR1
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
% Datum: 2019-03-09 04:59
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function taug = S6RPRRPR1_gravloadJ_floatb_twist_slag_vp2(qJ, g, ...
  pkin, m, mrSges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(3,1),zeros(11,1),zeros(7,1),zeros(7,3)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPRRPR1_gravloadJ_floatb_twist_slag_vp2: qJ has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6RPRRPR1_gravloadJ_floatb_twist_slag_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [11 1]), ...
  'S6RPRRPR1_gravloadJ_floatb_twist_slag_vp2: pkin has to be [11x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6RPRRPR1_gravloadJ_floatb_twist_slag_vp2: m has to be [7x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [7,3]), ...
  'S6RPRRPR1_gravloadJ_floatb_twist_slag_vp2: mrSges has to be [7x3] (double)');

%% Symbolic Calculation
% From gravload_joint_floatb_twist_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 04:57:32
% EndTime: 2019-03-09 04:57:33
% DurationCPUTime: 0.51s
% Computational Cost: add. (450->94), mult. (339->104), div. (0->0), fcn. (277->12), ass. (0->55)
t33 = qJ(3) + qJ(4);
t26 = pkin(11) + t33;
t20 = sin(t26);
t21 = cos(t26);
t34 = sin(qJ(6));
t67 = t34 * mrSges(7,2);
t95 = t20 * t67 + t21 * (m(7) * pkin(9) + mrSges(7,3));
t27 = sin(t33);
t28 = cos(t33);
t94 = mrSges(5,1) * t27 + mrSges(6,1) * t20 + mrSges(5,2) * t28 + mrSges(6,2) * t21;
t32 = qJ(1) + pkin(10);
t24 = sin(t32);
t25 = cos(t32);
t93 = g(1) * t25 + g(2) * t24;
t92 = -t28 * mrSges(5,1) - t21 * mrSges(6,1) + t27 * mrSges(5,2) + (mrSges(6,2) - mrSges(7,3)) * t20;
t89 = m(6) + m(7);
t56 = t21 * pkin(5) + t20 * pkin(9);
t86 = m(3) + m(4) + m(5);
t37 = cos(qJ(6));
t66 = t37 * mrSges(7,1);
t85 = -(t66 - t67) * t21 + t92;
t38 = cos(qJ(3));
t29 = t38 * pkin(3);
t35 = sin(qJ(3));
t53 = t38 * mrSges(4,1) - t35 * mrSges(4,2);
t83 = mrSges(3,1) + m(5) * (t29 + pkin(2)) + m(4) * pkin(2) + t53 - t92;
t40 = -pkin(8) - pkin(7);
t82 = -m(4) * pkin(7) + m(5) * t40 + mrSges(3,2) - mrSges(4,3) - mrSges(5,3) - mrSges(6,3);
t81 = pkin(4) * t27;
t22 = pkin(4) * t28;
t80 = pkin(5) * t20;
t77 = t35 * pkin(3);
t36 = sin(qJ(1));
t76 = t36 * pkin(1);
t39 = cos(qJ(1));
t30 = t39 * pkin(1);
t72 = t24 * t34;
t71 = t24 * t37;
t70 = t25 * t34;
t69 = t25 * t37;
t65 = t22 + t29;
t63 = t20 * t66;
t62 = t22 + t56;
t59 = t95 * t24;
t58 = t95 * t25;
t13 = -t77 - t81;
t42 = m(7) * (t13 - t80) - t63;
t41 = m(7) * (-t80 - t81) - t63;
t31 = -qJ(5) + t40;
t12 = pkin(2) + t65;
t4 = t21 * t69 + t72;
t3 = -t21 * t70 + t71;
t2 = -t21 * t71 + t70;
t1 = t21 * t72 + t69;
t5 = [(-t39 * mrSges(2,1) - t4 * mrSges(7,1) + t36 * mrSges(2,2) - t3 * mrSges(7,2) - t89 * (t25 * t12 - t24 * t31 + t30) - t86 * t30 + t82 * t24 + (-m(7) * t56 - t83) * t25) * g(2) + (t36 * mrSges(2,1) - t2 * mrSges(7,1) + t39 * mrSges(2,2) - t1 * mrSges(7,2) - t89 * (-t25 * t31 - t76) + t86 * t76 + t82 * t25 + (m(6) * t12 - m(7) * (-t12 - t56) + t83) * t24) * g(1) (-t86 - t89) * g(3), -g(1) * (t42 * t25 + t58) - g(2) * (t42 * t24 + t59) + (-t53 - m(5) * t29 - m(6) * t65 - m(7) * (t29 + t62) + t85) * g(3) + t93 * (m(5) * t77 - m(6) * t13 + mrSges(4,1) * t35 + mrSges(4,2) * t38 + t94) -g(1) * (t41 * t25 + t58) - g(2) * (t41 * t24 + t59) + (-m(6) * t22 - m(7) * t62 + t85) * g(3) + (m(6) * t81 + t94) * t93, t89 * (-g(1) * t24 + g(2) * t25) -g(1) * (t3 * mrSges(7,1) - t4 * mrSges(7,2)) - g(2) * (-t1 * mrSges(7,1) + t2 * mrSges(7,2)) - g(3) * (-mrSges(7,1) * t34 - mrSges(7,2) * t37) * t20];
taug  = t5(:);
