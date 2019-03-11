% Calculate Gravitation load on the joints for
% S6RRPRRP4
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% g [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [10x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d2,d4,d5,theta3]';
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
% Datum: 2019-03-09 11:56
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function taug = S6RRPRRP4_gravloadJ_floatb_twist_slag_vp2(qJ, g, ...
  pkin, m, mrSges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(3,1),zeros(10,1),zeros(7,1),zeros(7,3)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRPRRP4_gravloadJ_floatb_twist_slag_vp2: qJ has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6RRPRRP4_gravloadJ_floatb_twist_slag_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S6RRPRRP4_gravloadJ_floatb_twist_slag_vp2: pkin has to be [10x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6RRPRRP4_gravloadJ_floatb_twist_slag_vp2: m has to be [7x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [7,3]), ...
  'S6RRPRRP4_gravloadJ_floatb_twist_slag_vp2: mrSges has to be [7x3] (double)');

%% Symbolic Calculation
% From gravload_joint_floatb_twist_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 11:52:35
% EndTime: 2019-03-09 11:52:37
% DurationCPUTime: 0.94s
% Computational Cost: add. (541->118), mult. (585->136), div. (0->0), fcn. (547->10), ass. (0->65)
t121 = mrSges(6,3) + mrSges(7,2) - mrSges(4,2) + mrSges(5,3);
t37 = qJ(2) + pkin(10);
t32 = sin(t37);
t33 = cos(t37);
t41 = sin(qJ(2));
t44 = cos(qJ(2));
t122 = -t44 * mrSges(3,1) - t33 * mrSges(4,1) + t41 * mrSges(3,2) - t121 * t32;
t118 = mrSges(6,1) + mrSges(7,1);
t116 = mrSges(6,2) - mrSges(7,3);
t120 = m(3) * pkin(1) + mrSges(2,1) - t122;
t119 = m(4) + m(5);
t109 = -m(6) - m(7);
t38 = qJ(4) + qJ(5);
t34 = sin(t38);
t35 = cos(t38);
t40 = sin(qJ(4));
t43 = cos(qJ(4));
t114 = m(5) * pkin(3) + t43 * mrSges(5,1) - t40 * mrSges(5,2) - t116 * t34 + t118 * t35;
t30 = pkin(4) * t43 + pkin(3);
t21 = t33 * t30;
t46 = -pkin(9) - pkin(8);
t83 = t32 * t46;
t112 = t21 - t83;
t42 = sin(qJ(1));
t76 = t42 * t43;
t45 = cos(qJ(1));
t78 = t40 * t45;
t17 = -t33 * t78 + t76;
t111 = -m(7) * qJ(6) - mrSges(7,3);
t74 = t45 * t34;
t77 = t42 * t35;
t13 = t33 * t74 - t77;
t81 = t35 * t45;
t82 = t34 * t42;
t14 = t33 * t81 + t82;
t106 = t116 * t14 + t118 * t13;
t11 = t33 * t82 + t81;
t12 = t33 * t77 - t74;
t105 = t118 * t11 + t116 * t12;
t104 = -m(3) * pkin(7) + mrSges(2,2) - mrSges(3,3) - mrSges(4,3);
t102 = m(7) * pkin(5) + t118;
t101 = -mrSges(6,2) - t111;
t95 = pkin(2) * t41;
t94 = pkin(4) * t40;
t93 = pkin(8) * t32;
t90 = g(3) * t32;
t36 = t44 * pkin(2);
t89 = mrSges(6,2) * t35;
t39 = -qJ(3) - pkin(7);
t80 = t39 * t45;
t79 = t40 * t42;
t75 = t43 * t45;
t70 = t111 * t32 * t35;
t69 = -t11 * pkin(5) + qJ(6) * t12;
t31 = t36 + pkin(1);
t68 = t45 * t31 - t42 * t39;
t66 = -t13 * pkin(5) + qJ(6) * t14;
t63 = pkin(3) * t33 + t93;
t56 = pkin(5) * t35 + qJ(6) * t34;
t55 = t17 * pkin(4);
t15 = t33 * t79 + t75;
t52 = t15 * pkin(4);
t18 = t33 * t75 + t79;
t16 = -t33 * t76 + t78;
t1 = [(-t18 * mrSges(5,1) - t17 * mrSges(5,2) - t119 * t68 + t109 * (pkin(4) * t79 + t68) - t102 * t14 - t101 * t13 + t104 * t42) * g(2) + (m(5) * t80 - t16 * mrSges(5,1) - t15 * mrSges(5,2) + t109 * (pkin(4) * t78 + t42 * t83 - t80) + t102 * t12 + t101 * t11 + (m(4) * t31 - m(5) * (-t31 - t63) + t109 * (-t31 - t21) + t120) * t42) * g(1) + ((-m(5) * t63 + t109 * t112 - t120) * g(2) + (m(4) * t39 + t104) * g(1)) * t45 (-m(4) * t36 - m(5) * (t36 + t93) + t109 * (t36 + t112) + (-m(7) * t56 - t114) * t33 + t122) * g(3) + (g(1) * t45 + g(2) * t42) * (mrSges(3,1) * t41 + mrSges(3,2) * t44 + t119 * t95 + t109 * (-t33 * t46 - t95) + (-m(5) * pkin(8) - t121) * t33 + (mrSges(4,1) + m(6) * t30 - m(7) * (-t30 - t56) + t114) * t32) (-g(1) * t42 + g(2) * t45) * (t119 - t109) -g(3) * ((m(7) * (-pkin(5) * t34 - t94) - t34 * mrSges(7,1)) * t32 - t70) + (m(6) * t94 + mrSges(5,1) * t40 + mrSges(6,1) * t34 + mrSges(5,2) * t43 + t89) * t90 + (t15 * mrSges(5,1) - t16 * mrSges(5,2) + m(6) * t52 - m(7) * (-t52 + t69) + t105) * g(2) + (-t17 * mrSges(5,1) + t18 * mrSges(5,2) - m(6) * t55 - m(7) * (t55 + t66) + t106) * g(1) ((t102 * t34 + t89) * t32 + t70) * g(3) + (-m(7) * t69 + t105) * g(2) + (-m(7) * t66 + t106) * g(1) (-g(1) * t13 - g(2) * t11 - t34 * t90) * m(7)];
taug  = t1(:);
