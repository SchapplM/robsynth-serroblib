% Calculate Gravitation load on the joints for
% S6RPRRPR4
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

% Quelle: HybrDyn-Toolbox (ehem. IRT-Maple-Toolbox)
% Datum: 2018-11-23 16:17
% Revision: 76f9d5e39f14dc242b53c0d9d3d9db48bd8f37c0
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für mechatronische Systeme, Universität Hannover

function taug = S6RPRRPR4_gravloadJ_floatb_twist_slag_vp2(qJ, g, ...
  pkin, m, mrSges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(3,1),zeros(11,1),zeros(7,1),zeros(7,3)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPRRPR4_gravloadJ_floatb_twist_slag_vp2: qJ has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6RPRRPR4_gravloadJ_floatb_twist_slag_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [11 1]), ...
  'S6RPRRPR4_gravloadJ_floatb_twist_slag_vp2: pkin has to be [11x1] (double)');
assert( isreal(m) && all(size(m) == [7 1]), ...
  'S6RPRRPR4_gravloadJ_floatb_twist_slag_vp2: m has to be [7x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [7,3]), ...
  'S6RPRRPR4_gravloadJ_floatb_twist_slag_vp2: mrSges has to be [7x3] (double)');

%% Symbolic Calculation
% From joint_gravload_floatb_twist_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2018-11-23 16:17:23
% EndTime: 2018-11-23 16:17:24
% DurationCPUTime: 0.78s
% Computational Cost: add. (494->108), mult. (422->112), div. (0->0), fcn. (357->12), ass. (0->60)
t117 = m(6) * qJ(5) + mrSges(6,3);
t38 = pkin(10) + qJ(3);
t34 = qJ(4) + t38;
t26 = sin(t34);
t27 = cos(t34);
t39 = sin(pkin(11));
t81 = t39 * mrSges(6,2);
t37 = pkin(11) + qJ(6);
t30 = sin(t37);
t83 = t30 * mrSges(7,2);
t116 = (-t81 - t83) * t26 + (-mrSges(7,3) - t117) * t27;
t41 = cos(pkin(11));
t28 = t41 * pkin(5) + pkin(4);
t43 = -pkin(9) - qJ(5);
t107 = -t26 * t43 + t27 * t28;
t111 = m(7) * t107;
t80 = t41 * mrSges(6,1);
t32 = cos(t37);
t82 = t32 * mrSges(7,1);
t110 = (t80 + t82) * t26;
t109 = t80 - t81;
t108 = -t27 * mrSges(5,1) + (mrSges(5,2) - mrSges(7,3)) * t26;
t97 = m(6) + m(7);
t57 = -t26 * t28 - t27 * t43;
t95 = pkin(4) * t26;
t31 = sin(t38);
t96 = pkin(3) * t31;
t105 = -m(7) * (t57 - t96) - m(6) * (-t95 - t96) + t110;
t45 = sin(qJ(1));
t46 = cos(qJ(1));
t104 = g(1) * t46 + g(2) * t45;
t22 = t26 * mrSges(6,3);
t103 = -t22 + t108 + (-t82 + t83 - t109) * t27;
t102 = t116 * t45;
t101 = t116 * t46;
t100 = m(6) * t95 - m(7) * t57 + t110;
t44 = -pkin(7) - qJ(2);
t36 = -pkin(8) + t44;
t99 = -m(3) * qJ(2) + m(4) * t44 + m(6) * t36 - t39 * mrSges(6,1) - t41 * mrSges(6,2) + mrSges(2,2) - mrSges(3,3) - mrSges(4,3) - mrSges(5,3);
t42 = cos(pkin(10));
t29 = t42 * pkin(2) + pkin(1);
t33 = cos(t38);
t62 = t33 * mrSges(4,1) - t31 * mrSges(4,2);
t98 = -m(4) * t29 - (m(6) * pkin(4) + t109) * t27 - mrSges(2,1) - m(3) * pkin(1) - t42 * mrSges(3,1) + sin(pkin(10)) * mrSges(3,2) - t62 + t108;
t25 = pkin(3) * t33;
t94 = pkin(5) * t39;
t79 = t45 * t30;
t78 = t45 * t32;
t77 = t46 * t30;
t76 = t46 * t32;
t20 = t26 * qJ(5);
t75 = t27 * pkin(4) + t20;
t59 = mrSges(5,1) * t26 + mrSges(5,2) * t27;
t12 = t25 + t29;
t5 = t46 * t12;
t4 = t27 * t76 + t79;
t3 = -t27 * t77 + t78;
t2 = -t27 * t78 + t77;
t1 = t27 * t79 + t76;
t6 = [(-m(6) * t5 - t4 * mrSges(7,1) - t3 * mrSges(7,2) + (-m(5) - m(7)) * (-t45 * t36 + t5) + (-m(7) * t94 + t99) * t45 + (-t117 * t26 - t111 + t98) * t46) * g(2) + (-t2 * mrSges(7,1) - t1 * mrSges(7,2) + (m(5) * t36 - m(7) * (-t36 + t94) + t99) * t46 + (m(5) * t12 - m(6) * (-t12 - t20) + t22 - m(7) * (-t12 - t107) - t98) * t45) * g(1) (-g(1) * t45 + g(2) * t46) * (m(3) + m(4) + m(5) + t97) (t105 * t45 + t102) * g(2) + (t105 * t46 + t101) * g(1) + (-t62 - m(5) * t25 - m(6) * (t25 + t75) - m(7) * (t25 + t107) + t103) * g(3) + (m(5) * t96 + mrSges(4,1) * t31 + mrSges(4,2) * t33 + t59) * t104, t104 * t59 + (t100 * t45 + t102) * g(2) + (t100 * t46 + t101) * g(1) + (-m(6) * t75 + t103 - t111) * g(3) (g(3) * t27 - t104 * t26) * t97, -g(1) * (t3 * mrSges(7,1) - t4 * mrSges(7,2)) - g(2) * (-t1 * mrSges(7,1) + t2 * mrSges(7,2)) - g(3) * (-mrSges(7,1) * t30 - mrSges(7,2) * t32) * t26];
taug  = t6(:);
