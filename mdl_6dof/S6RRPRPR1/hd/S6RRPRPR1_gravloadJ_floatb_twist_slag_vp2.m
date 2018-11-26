% Calculate Gravitation load on the joints for
% S6RRPRPR1
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% g [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [11x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d2,d4,d6,theta3,theta5]';
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
% Datum: 2018-11-23 17:00
% Revision: 76f9d5e39f14dc242b53c0d9d3d9db48bd8f37c0
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für mechatronische Systeme, Universität Hannover

function taug = S6RRPRPR1_gravloadJ_floatb_twist_slag_vp2(qJ, g, ...
  pkin, m, mrSges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(3,1),zeros(11,1),zeros(7,1),zeros(7,3)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRPRPR1_gravloadJ_floatb_twist_slag_vp2: qJ has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6RRPRPR1_gravloadJ_floatb_twist_slag_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [11 1]), ...
  'S6RRPRPR1_gravloadJ_floatb_twist_slag_vp2: pkin has to be [11x1] (double)');
assert( isreal(m) && all(size(m) == [7 1]), ...
  'S6RRPRPR1_gravloadJ_floatb_twist_slag_vp2: m has to be [7x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [7,3]), ...
  'S6RRPRPR1_gravloadJ_floatb_twist_slag_vp2: mrSges has to be [7x3] (double)');

%% Symbolic Calculation
% From joint_gravload_floatb_twist_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2018-11-23 16:59:48
% EndTime: 2018-11-23 16:59:49
% DurationCPUTime: 0.72s
% Computational Cost: add. (510->111), mult. (448->119), div. (0->0), fcn. (379->12), ass. (0->58)
t42 = sin(pkin(11));
t43 = cos(pkin(11));
t112 = mrSges(6,1) * t43 - mrSges(6,2) * t42;
t40 = pkin(11) + qJ(6);
t33 = sin(t40);
t35 = cos(t40);
t122 = -mrSges(7,1) * t35 + mrSges(7,2) * t33 - t112;
t41 = qJ(2) + pkin(10);
t37 = qJ(4) + t41;
t29 = sin(t37);
t30 = cos(t37);
t31 = pkin(5) * t43 + pkin(4);
t45 = -pkin(9) - qJ(5);
t120 = (m(7) * t31 - t122) * t29 + (m(7) * t45 - mrSges(6,3) - mrSges(7,3)) * t30;
t110 = -t29 * t45 + t30 * t31;
t117 = m(7) * t110;
t34 = sin(t41);
t36 = cos(t41);
t46 = sin(qJ(2));
t48 = cos(qJ(2));
t115 = -t48 * mrSges(3,1) - t36 * mrSges(4,1) + t46 * mrSges(3,2) + t34 * mrSges(4,2);
t113 = -pkin(4) * t29 + qJ(5) * t30;
t47 = sin(qJ(1));
t49 = cos(qJ(1));
t106 = g(1) * t49 + g(2) * t47;
t111 = -t30 * mrSges(5,1) + (mrSges(5,2) - mrSges(7,3)) * t29;
t101 = m(6) + m(7);
t108 = t120 * t49;
t107 = t120 * t47;
t25 = t29 * mrSges(6,3);
t105 = t122 * t30 + t111 - t25;
t44 = -qJ(3) - pkin(7);
t39 = -pkin(8) + t44;
t103 = -m(3) * pkin(7) + m(4) * t44 + m(6) * t39 - t42 * mrSges(6,1) - t43 * mrSges(6,2) + mrSges(2,2) - mrSges(3,3) - mrSges(4,3) - mrSges(5,3);
t38 = t48 * pkin(2);
t102 = -(m(6) * pkin(4) + t112) * t30 - mrSges(2,1) - m(4) * (t38 + pkin(1)) - m(3) * pkin(1) + t111 + t115;
t100 = pkin(2) * t46;
t98 = pkin(5) * t42;
t87 = t33 * t49;
t86 = t35 * t49;
t85 = t47 * t33;
t84 = t47 * t35;
t23 = t29 * qJ(5);
t83 = t30 * pkin(4) + t23;
t82 = pkin(3) * t36 + t38;
t71 = t113 * t47;
t70 = t113 * t49;
t64 = mrSges(5,1) * t29 + mrSges(5,2) * t30;
t15 = -pkin(3) * t34 - t100;
t14 = pkin(1) + t82;
t9 = t49 * t15;
t8 = t47 * t15;
t5 = t49 * t14;
t4 = t30 * t86 + t85;
t3 = -t30 * t87 + t84;
t2 = -t30 * t84 + t87;
t1 = t30 * t85 + t86;
t6 = [(-m(6) * t5 - t4 * mrSges(7,1) - t3 * mrSges(7,2) + (-m(5) - m(7)) * (-t47 * t39 + t5) + (-m(7) * t98 + t103) * t47 + (-(m(6) * qJ(5) + mrSges(6,3)) * t29 - t117 + t102) * t49) * g(2) + (-t2 * mrSges(7,1) - t1 * mrSges(7,2) + (m(5) * t39 - m(7) * (-t39 + t98) + t103) * t49 + (m(5) * t14 - m(6) * (-t14 - t23) + t25 - m(7) * (-t14 - t110) - t102) * t47) * g(1) (-m(6) * (t71 + t8) - m(7) * t8 + t107) * g(2) + (-m(6) * (t70 + t9) - m(7) * t9 + t108) * g(1) + (-m(4) * t38 - m(5) * t82 - m(6) * (t82 + t83) - m(7) * (t110 + t82) + t105 + t115) * g(3) + t106 * (m(4) * t100 - m(5) * t15 + mrSges(3,1) * t46 + mrSges(4,1) * t34 + mrSges(3,2) * t48 + mrSges(4,2) * t36 + t64) (-g(1) * t47 + g(2) * t49) * (m(4) + m(5) + t101) t106 * t64 + (-m(6) * t71 + t107) * g(2) + (-m(6) * t70 + t108) * g(1) + (-m(6) * t83 + t105 - t117) * g(3) (g(3) * t30 - t106 * t29) * t101, -g(1) * (mrSges(7,1) * t3 - mrSges(7,2) * t4) - g(2) * (-mrSges(7,1) * t1 + mrSges(7,2) * t2) - g(3) * (-mrSges(7,1) * t33 - mrSges(7,2) * t35) * t29];
taug  = t6(:);
