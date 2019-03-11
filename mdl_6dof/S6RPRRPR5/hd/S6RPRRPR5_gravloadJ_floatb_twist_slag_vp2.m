% Calculate Gravitation load on the joints for
% S6RPRRPR5
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% g [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [10x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d3,d4,d6,theta2]';
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
% Datum: 2019-03-09 05:14
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function taug = S6RPRRPR5_gravloadJ_floatb_twist_slag_vp2(qJ, g, ...
  pkin, m, mrSges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(3,1),zeros(10,1),zeros(7,1),zeros(7,3)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPRRPR5_gravloadJ_floatb_twist_slag_vp2: qJ has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6RPRRPR5_gravloadJ_floatb_twist_slag_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S6RPRRPR5_gravloadJ_floatb_twist_slag_vp2: pkin has to be [10x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6RPRRPR5_gravloadJ_floatb_twist_slag_vp2: m has to be [7x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [7,3]), ...
  'S6RPRRPR5_gravloadJ_floatb_twist_slag_vp2: mrSges has to be [7x3] (double)');

%% Symbolic Calculation
% From gravload_joint_floatb_twist_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 05:12:13
% EndTime: 2019-03-09 05:12:15
% DurationCPUTime: 0.61s
% Computational Cost: add. (446->109), mult. (394->119), div. (0->0), fcn. (325->10), ass. (0->58)
t40 = sin(qJ(6));
t42 = cos(qJ(6));
t102 = mrSges(7,1) * t40 + mrSges(7,2) * t42;
t36 = pkin(10) + qJ(3);
t33 = qJ(4) + t36;
t29 = cos(t33);
t101 = t102 * t29;
t28 = sin(t33);
t100 = (-mrSges(5,1) + mrSges(6,2)) * t29 + (mrSges(5,2) - mrSges(6,3)) * t28;
t90 = m(6) + m(7);
t21 = t28 * qJ(5);
t43 = cos(qJ(1));
t78 = t29 * t43;
t99 = pkin(4) * t78 + t43 * t21;
t70 = qJ(5) * t29;
t13 = t43 * t70;
t80 = t28 * t43;
t96 = -m(7) * t13 - mrSges(6,2) * t80 - mrSges(6,3) * t78 - t101 * t43;
t41 = sin(qJ(1));
t11 = t41 * t70;
t81 = t28 * t41;
t95 = -m(7) * t11 - mrSges(6,2) * t81 + (-mrSges(6,3) * t29 - t101) * t41;
t94 = g(1) * t43 + g(2) * t41;
t93 = -t29 * mrSges(7,3) - t102 * t28 + t100;
t38 = cos(pkin(10));
t30 = t38 * pkin(2) + pkin(1);
t31 = sin(t36);
t32 = cos(t36);
t53 = t32 * mrSges(4,1) - t31 * mrSges(4,2);
t92 = -m(4) * t30 - mrSges(2,1) - m(3) * pkin(1) - t38 * mrSges(3,1) + sin(pkin(10)) * mrSges(3,2) - t53 + t100;
t39 = -pkin(7) - qJ(2);
t35 = -pkin(8) + t39;
t91 = -m(3) * qJ(2) + m(4) * t39 - m(7) * (pkin(5) - t35) - mrSges(6,1) + mrSges(2,2) - mrSges(3,3) - mrSges(4,3) - mrSges(5,3);
t89 = pkin(3) * t31;
t27 = pkin(3) * t32;
t86 = g(3) * t29;
t26 = t29 * pkin(4);
t77 = t41 * t40;
t76 = t41 * t42;
t75 = t43 * t40;
t74 = t43 * t42;
t71 = t26 + t21;
t66 = t27 + t71;
t10 = t27 + t30;
t5 = t43 * t10;
t64 = -t41 * t35 + t5;
t60 = -t10 - t21;
t56 = m(7) * (-pkin(4) - pkin(9)) - mrSges(7,3);
t55 = -pkin(4) * t28 - t89;
t50 = mrSges(5,1) * t28 + mrSges(5,2) * t29;
t48 = t56 * t28;
t44 = -m(7) * t89 + t48;
t25 = t29 * pkin(9);
t4 = -t28 * t77 + t74;
t3 = t28 * t76 + t75;
t2 = t28 * t75 + t76;
t1 = t28 * t74 - t77;
t6 = [(-m(5) * t64 - m(6) * (t64 + t99) - m(7) * (pkin(9) * t78 + t5 + t99) - t2 * mrSges(7,1) - t1 * mrSges(7,2) - mrSges(7,3) * t78 + t92 * t43 + t91 * t41) * g(2) + (-t4 * mrSges(7,1) + t3 * mrSges(7,2) + ((m(5) + m(6)) * t35 + t91) * t43 + (m(5) * t10 - m(6) * (t60 - t26) - m(7) * t60 - t56 * t29 - t92) * t41) * g(1) (-g(1) * t41 + g(2) * t43) * (m(3) + m(4) + m(5) + t90) (-m(6) * (t55 * t41 + t11) - t44 * t41 + t95) * g(2) + (-m(6) * (t55 * t43 + t13) - t44 * t43 + t96) * g(1) + (-t53 - m(5) * t27 - m(6) * t66 - m(7) * (t25 + t66) + t93) * g(3) + (m(5) * t89 + mrSges(4,1) * t31 + mrSges(4,2) * t32 + t50) * t94, t94 * t50 + (-m(6) * (-pkin(4) * t81 + t11) - t41 * t48 + t95) * g(2) + (-m(6) * (-pkin(4) * t80 + t13) - t43 * t48 + t96) * g(1) + (-m(6) * t71 - m(7) * (t25 + t71) + t93) * g(3) (-t28 * t94 + t86) * t90, -g(1) * (t1 * mrSges(7,1) - t2 * mrSges(7,2)) - g(2) * (t3 * mrSges(7,1) + t4 * mrSges(7,2)) - (-mrSges(7,1) * t42 + mrSges(7,2) * t40) * t86];
taug  = t6(:);
