% Calculate Gravitation load on the joints for
% S6RRPRRP9
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% g [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [11x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,alpha2,d1,d2,d4,d5,theta3]';
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
% Datum: 2018-11-23 17:17
% Revision: 76f9d5e39f14dc242b53c0d9d3d9db48bd8f37c0
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für mechatronische Systeme, Universität Hannover

function taug = S6RRPRRP9_gravloadJ_floatb_twist_slag_vp2(qJ, g, ...
  pkin, m, mrSges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(3,1),zeros(11,1),zeros(7,1),zeros(7,3)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRPRRP9_gravloadJ_floatb_twist_slag_vp2: qJ has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6RRPRRP9_gravloadJ_floatb_twist_slag_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [11 1]), ...
  'S6RRPRRP9_gravloadJ_floatb_twist_slag_vp2: pkin has to be [11x1] (double)');
assert( isreal(m) && all(size(m) == [7 1]), ...
  'S6RRPRRP9_gravloadJ_floatb_twist_slag_vp2: m has to be [7x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [7,3]), ...
  'S6RRPRRP9_gravloadJ_floatb_twist_slag_vp2: mrSges has to be [7x3] (double)');

%% Symbolic Calculation
% From joint_gravload_floatb_twist_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2018-11-23 17:16:43
% EndTime: 2018-11-23 17:16:44
% DurationCPUTime: 1.20s
% Computational Cost: add. (1468->135), mult. (1556->176), div. (0->0), fcn. (1524->16), ass. (0->65)
t68 = cos(qJ(5));
t55 = pkin(5) * t68 + pkin(4);
t129 = -m(6) * pkin(4) - m(7) * t55 - mrSges(5,1);
t74 = -m(6) * pkin(10) + m(7) * (-qJ(6) - pkin(10)) + mrSges(5,2) - mrSges(6,3) - mrSges(7,3);
t128 = mrSges(6,1) + mrSges(7,1);
t123 = mrSges(6,2) + mrSges(7,2);
t60 = sin(pkin(11));
t62 = cos(pkin(11));
t124 = -m(4) * pkin(2) - t62 * mrSges(4,1) + t60 * mrSges(4,2) - mrSges(3,1);
t59 = pkin(11) + qJ(4);
t56 = sin(t59);
t57 = cos(t59);
t127 = t129 * t57 + t74 * t56 + t124;
t115 = cos(qJ(1));
t96 = pkin(6) + qJ(2);
t83 = sin(t96) / 0.2e1;
t97 = pkin(6) - qJ(2);
t88 = sin(t97);
t42 = t83 - t88 / 0.2e1;
t67 = sin(qJ(1));
t69 = cos(qJ(2));
t79 = t115 * t42 + t67 * t69;
t61 = sin(pkin(6));
t93 = t61 * t115;
t14 = -t56 * t93 + t57 * t79;
t66 = sin(qJ(2));
t84 = cos(t96) / 0.2e1;
t89 = cos(t97);
t70 = t89 / 0.2e1 + t84;
t31 = -t115 * t70 + t66 * t67;
t65 = sin(qJ(5));
t126 = t14 * t65 - t31 * t68;
t125 = -t14 * t68 - t31 * t65;
t116 = m(7) * pkin(5);
t122 = -m(5) - m(6) - m(7);
t120 = -t123 * t65 + t128 * t68 - t129;
t119 = -m(4) * qJ(3) + mrSges(3,2) - mrSges(4,3) - mrSges(5,3);
t118 = -t116 - t128;
t112 = t79 * t65;
t78 = t115 * t69 - t67 * t42;
t111 = t78 * t65;
t43 = t84 - t89 / 0.2e1;
t110 = t43 * t65;
t107 = t57 * t65;
t106 = t57 * t68;
t105 = t61 * t67;
t99 = t115 * pkin(1) + pkin(8) * t105;
t98 = cos(pkin(6));
t91 = -pkin(1) * t67 + pkin(8) * t93;
t87 = t60 * t93;
t34 = t115 * t66 + t67 * t70;
t54 = pkin(3) * t62 + pkin(2);
t64 = -pkin(9) - qJ(3);
t86 = t60 * pkin(3) * t105 - t34 * t64 + t54 * t78 + t99;
t18 = t105 * t56 + t57 * t78;
t5 = -t18 * t65 + t34 * t68;
t80 = pkin(3) * t87 + t31 * t64 - t54 * t79 + t91;
t13 = t56 * t79 + t57 * t93;
t73 = t116 * t65 - t119;
t41 = t83 + t88 / 0.2e1;
t26 = -t43 * t57 + t56 * t98;
t25 = -t43 * t56 - t57 * t98;
t17 = -t105 * t57 + t56 * t78;
t6 = t18 * t68 + t34 * t65;
t1 = [(-t115 * mrSges(2,1) - m(5) * t86 - t18 * mrSges(5,1) - m(6) * (pkin(4) * t18 + t86) - m(7) * (t18 * t55 + t86) + t124 * t78 + (mrSges(2,2) + (-mrSges(4,1) * t60 - mrSges(4,2) * t62 - mrSges(3,3)) * t61) * t67 - t128 * t6 - t123 * t5 - t73 * t34 + t74 * t17 + (-m(3) - m(4)) * t99) * g(2) + (t67 * mrSges(2,1) + t115 * mrSges(2,2) - m(3) * t91 + t79 * mrSges(3,1) - mrSges(3,3) * t93 - m(4) * (-pkin(2) * t79 + t91) - (-t62 * t79 + t87) * mrSges(4,1) - (t60 * t79 + t62 * t93) * mrSges(4,2) - m(5) * t80 + t14 * mrSges(5,1) - m(6) * (-pkin(4) * t14 + t80) - m(7) * (-t14 * t55 + t80) - t128 * t125 - t123 * t126 + t73 * t31 - t74 * t13) * g(1) (t110 * t116 - t128 * (t106 * t41 - t110) - t123 * (-t107 * t41 - t43 * t68) + t122 * (t41 * t54 + t43 * t64) - t119 * t43 + t127 * t41) * g(3) + (-t112 * t116 - t128 * (-t106 * t31 + t112) - t123 * (t107 * t31 + t68 * t79) + t122 * (-t31 * t54 - t64 * t79) + t119 * t79 - t127 * t31) * g(2) + (-t111 * t116 - t123 * (t107 * t34 + t68 * t78) + t122 * (-t34 * t54 - t64 * t78) - t128 * (-t106 * t34 + t111) + t119 * t78 - t127 * t34) * g(1) (-g(1) * t34 - g(2) * t31 + g(3) * t41) * (m(4) - t122) (t120 * t25 + t74 * t26) * g(3) + (t120 * t13 + t74 * t14) * g(2) + (t120 * t17 + t74 * t18) * g(1) (-t123 * (-t26 * t68 + t41 * t65) + t118 * (-t26 * t65 - t41 * t68)) * g(3) + (-t118 * t126 - t123 * t125) * g(2) + (t118 * t5 + t123 * t6) * g(1) (-g(1) * t17 - g(2) * t13 - g(3) * t25) * m(7)];
taug  = t1(:);
