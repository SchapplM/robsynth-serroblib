% Calculate Gravitation load on the joints for
% S6RRPRPR5
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% g [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [12x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,alpha2,d1,d2,d4,d6,theta3,theta5]';
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
% Datum: 2018-11-23 17:03
% Revision: 76f9d5e39f14dc242b53c0d9d3d9db48bd8f37c0
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für mechatronische Systeme, Universität Hannover

function taug = S6RRPRPR5_gravloadJ_floatb_twist_slag_vp2(qJ, g, ...
  pkin, m, mrSges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(3,1),zeros(12,1),zeros(7,1),zeros(7,3)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRPRPR5_gravloadJ_floatb_twist_slag_vp2: qJ has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6RRPRPR5_gravloadJ_floatb_twist_slag_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [12 1]), ...
  'S6RRPRPR5_gravloadJ_floatb_twist_slag_vp2: pkin has to be [12x1] (double)');
assert( isreal(m) && all(size(m) == [7 1]), ...
  'S6RRPRPR5_gravloadJ_floatb_twist_slag_vp2: m has to be [7x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [7,3]), ...
  'S6RRPRPR5_gravloadJ_floatb_twist_slag_vp2: mrSges has to be [7x3] (double)');

%% Symbolic Calculation
% From joint_gravload_floatb_twist_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2018-11-23 17:02:46
% EndTime: 2018-11-23 17:02:48
% DurationCPUTime: 1.36s
% Computational Cost: add. (1696->129), mult. (1424->166), div. (0->0), fcn. (1351->22), ass. (0->77)
t55 = pkin(12) + qJ(6);
t47 = sin(t55);
t51 = cos(t55);
t59 = sin(pkin(12));
t61 = cos(pkin(12));
t73 = -mrSges(5,1) - m(6) * pkin(4) - mrSges(6,1) * t61 + mrSges(6,2) * t59 - m(7) * (pkin(5) * t61 + pkin(4));
t137 = -mrSges(7,1) * t51 + mrSges(7,2) * t47 + t73;
t76 = mrSges(5,2) - m(6) * qJ(5) - mrSges(6,3) + m(7) * (-pkin(10) - qJ(5)) - mrSges(7,3);
t128 = m(6) + m(7);
t101 = m(5) + t128;
t120 = -t61 * mrSges(6,2) - (m(7) * pkin(5) + mrSges(6,1)) * t59 + mrSges(4,2) - mrSges(5,3) - pkin(9) * t101;
t136 = t47 * mrSges(7,1) + t51 * mrSges(7,2) - t120;
t63 = sin(qJ(4));
t66 = cos(qJ(4));
t121 = -t137 * t66 - t76 * t63 + mrSges(4,1);
t98 = m(4) + t101;
t57 = pkin(6) + qJ(2);
t43 = cos(t57) / 0.2e1;
t58 = pkin(6) - qJ(2);
t54 = cos(t58);
t32 = t54 / 0.2e1 + t43;
t56 = qJ(2) + pkin(11);
t93 = pkin(6) + t56;
t40 = sin(t93);
t118 = t40 / 0.2e1;
t117 = sin(t57) / 0.2e1;
t50 = sin(t58);
t114 = pkin(2) * t50;
t60 = sin(pkin(6));
t65 = sin(qJ(1));
t111 = t60 * t65;
t68 = cos(qJ(1));
t110 = t60 * t68;
t64 = sin(qJ(2));
t108 = t64 * t68;
t52 = cos(t56);
t107 = t65 * t52;
t106 = t65 * t64;
t105 = t68 * t52;
t38 = pkin(2) * t117;
t103 = t114 / 0.2e1 + t38;
t102 = cos(pkin(6));
t94 = pkin(6) - t56;
t86 = sin(t94);
t100 = t118 - t86 / 0.2e1;
t14 = t100 * t68 + t107;
t4 = -t63 * t110 + t14 * t66;
t30 = t32 * pkin(2);
t92 = -pkin(2) * t106 + t68 * t30;
t87 = cos(t93);
t67 = cos(qJ(2));
t85 = m(3) * pkin(1) + t67 * mrSges(3,1) + mrSges(2,1);
t84 = -pkin(2) * t108 - t65 * t30;
t3 = t110 * t66 + t14 * t63;
t83 = cos(t94) / 0.2e1;
t82 = t86 / 0.2e1;
t75 = -t40 / 0.2e1 + t82;
t72 = t87 / 0.2e1 + t83;
t31 = t117 - t50 / 0.2e1;
t69 = t31 * mrSges(3,1) + mrSges(2,2) + (-m(3) * pkin(8) - mrSges(3,3) - mrSges(4,3)) * t60 + t98 * (-t114 / 0.2e1 + t38 - t60 * (pkin(8) + qJ(3)));
t48 = sin(t56);
t46 = pkin(2) * t67 + pkin(1);
t33 = t68 * t46;
t29 = t83 - t87 / 0.2e1;
t28 = t82 + t118;
t22 = -t65 * t32 - t108;
t21 = -t32 * t68 + t106;
t20 = t102 * t63 + t29 * t66;
t19 = -t102 * t66 + t29 * t63;
t17 = -t100 * t65 + t105;
t16 = t68 * t48 + t65 * t72;
t13 = t48 * t65 - t68 * t72;
t8 = t111 * t63 + t17 * t66;
t7 = -t111 * t66 + t17 * t63;
t2 = t16 * t47 + t51 * t8;
t1 = t16 * t51 - t47 * t8;
t5 = [(-t21 * mrSges(3,2) + (t46 * t98 + t85) * t65 + t69 * t68 - t76 * t3 - t137 * t4 + t136 * t13 + (pkin(3) * t101 + mrSges(4,1)) * t14) * g(1) + (-m(4) * t33 - t17 * mrSges(4,1) - t2 * mrSges(7,1) - t22 * mrSges(3,2) - t1 * mrSges(7,2) + t120 * t16 + t69 * t65 - t85 * t68 + t76 * t7 + t73 * t8 - t101 * (t17 * pkin(3) + t33)) * g(2) (-(t117 + t50 / 0.2e1) * mrSges(3,1) - (t43 - t54 / 0.2e1) * mrSges(3,2) - m(4) * t103 - t101 * (t28 * pkin(3) + t103) - t136 * t29 - t121 * t28) * g(3) + (t21 * mrSges(3,1) - (-t31 * t68 - t65 * t67) * mrSges(3,2) - m(4) * t92 - t101 * (-t13 * pkin(3) + t92) - t136 * (-t68 * t75 + t107) + t121 * t13) * g(2) + (-t22 * mrSges(3,1) - (t65 * t31 - t67 * t68) * mrSges(3,2) - m(4) * t84 - t101 * (-t16 * pkin(3) + t84) - t136 * (t65 * t75 + t105) + t121 * t16) * g(1) ((-g(1) * t65 + g(2) * t68) * t60 - g(3) * t102) * t98 (-t137 * t19 + t76 * t20) * g(3) + (-t137 * t3 + t76 * t4) * g(2) + (-t137 * t7 + t76 * t8) * g(1), t128 * (-g(1) * t7 - g(2) * t3 - g(3) * t19) -g(1) * (mrSges(7,1) * t1 - mrSges(7,2) * t2) - g(2) * ((t13 * t51 - t4 * t47) * mrSges(7,1) + (-t13 * t47 - t4 * t51) * mrSges(7,2)) - g(3) * ((-t20 * t47 - t28 * t51) * mrSges(7,1) + (-t20 * t51 + t28 * t47) * mrSges(7,2))];
taug  = t5(:);
