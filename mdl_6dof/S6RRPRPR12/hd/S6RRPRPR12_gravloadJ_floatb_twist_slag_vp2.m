% Calculate Gravitation load on the joints for
% S6RRPRPR12
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% g [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [11x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,alpha2,d1,d2,d4,d6,theta5]';
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
% Datum: 2018-11-23 17:08
% Revision: 76f9d5e39f14dc242b53c0d9d3d9db48bd8f37c0
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für mechatronische Systeme, Universität Hannover

function taug = S6RRPRPR12_gravloadJ_floatb_twist_slag_vp2(qJ, g, ...
  pkin, m, mrSges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(3,1),zeros(11,1),zeros(7,1),zeros(7,3)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRPRPR12_gravloadJ_floatb_twist_slag_vp2: qJ has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6RRPRPR12_gravloadJ_floatb_twist_slag_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [11 1]), ...
  'S6RRPRPR12_gravloadJ_floatb_twist_slag_vp2: pkin has to be [11x1] (double)');
assert( isreal(m) && all(size(m) == [7 1]), ...
  'S6RRPRPR12_gravloadJ_floatb_twist_slag_vp2: m has to be [7x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [7,3]), ...
  'S6RRPRPR12_gravloadJ_floatb_twist_slag_vp2: mrSges has to be [7x3] (double)');

%% Symbolic Calculation
% From joint_gravload_floatb_twist_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2018-11-23 17:08:18
% EndTime: 2018-11-23 17:08:19
% DurationCPUTime: 1.07s
% Computational Cost: add. (1216->133), mult. (1329->173), div. (0->0), fcn. (1260->16), ass. (0->73)
t123 = m(6) + m(7);
t124 = m(4) + m(5);
t95 = t123 + t124;
t121 = -t95 * qJ(3) + mrSges(3,2) - mrSges(4,3);
t63 = sin(qJ(4));
t67 = cos(qJ(4));
t128 = t63 * mrSges(5,1) + t67 * mrSges(5,2) - t121;
t62 = sin(qJ(6));
t66 = cos(qJ(6));
t122 = m(7) * pkin(5) + t66 * mrSges(7,1) - t62 * mrSges(7,2) + mrSges(6,1);
t120 = m(7) * pkin(10) - mrSges(6,2) + mrSges(7,3);
t59 = sin(pkin(6));
t69 = cos(qJ(1));
t110 = t59 * t69;
t64 = sin(qJ(2));
t65 = sin(qJ(1));
t103 = pkin(6) - qJ(2);
t86 = cos(t103) / 0.2e1;
t102 = pkin(6) + qJ(2);
t93 = cos(t102);
t71 = t86 + t93 / 0.2e1;
t30 = t64 * t65 - t69 * t71;
t126 = t63 * t110 + t30 * t67;
t111 = t59 * t65;
t33 = t69 * t64 + t65 * t71;
t9 = -t63 * t111 + t33 * t67;
t91 = sin(t102);
t84 = t91 / 0.2e1;
t92 = sin(t103);
t85 = t92 / 0.2e1;
t45 = t84 + t85;
t60 = cos(pkin(6));
t125 = t45 * t67 + t60 * t63;
t77 = -m(5) * pkin(9) - mrSges(3,1) + mrSges(4,2) - mrSges(5,3) - mrSges(6,3);
t119 = t62 * mrSges(7,1) + t66 * mrSges(7,2) - t77;
t58 = qJ(4) + pkin(11);
t55 = sin(t58);
t56 = cos(t58);
t118 = t120 * t56 - t122 * t55 - t128;
t117 = pkin(4) * t63;
t116 = t30 * t63;
t114 = t33 * t63;
t68 = cos(qJ(2));
t108 = t65 * t68;
t107 = t69 * t68;
t104 = t69 * pkin(1) + pkin(8) * t111;
t78 = t84 - t92 / 0.2e1;
t34 = -t65 * t78 + t107;
t97 = t34 * pkin(2) + t104;
t94 = -t65 * pkin(1) + pkin(8) * t110;
t90 = -m(5) * pkin(3) - mrSges(4,1) - mrSges(3,3);
t31 = t69 * t78 + t108;
t87 = t31 * pkin(2) - t94;
t54 = pkin(4) * t67 + pkin(3);
t61 = -qJ(5) - pkin(9);
t80 = pkin(4) * t114 + t54 * t111 - t34 * t61 + t97;
t79 = t85 - t91 / 0.2e1;
t7 = t56 * t110 - t30 * t55;
t5 = t55 * t110 + t30 * t56;
t48 = t67 * t110;
t46 = t86 - t93 / 0.2e1;
t42 = t45 * pkin(2);
t35 = t65 * t79 + t107;
t32 = -t69 * t79 + t108;
t28 = t33 * pkin(2);
t26 = t30 * pkin(2);
t15 = -t45 * t55 + t56 * t60;
t10 = t67 * t111 + t114;
t4 = t56 * t111 + t33 * t55;
t3 = t55 * t111 - t33 * t56;
t2 = t34 * t62 + t4 * t66;
t1 = t34 * t66 - t4 * t62;
t6 = [(-t69 * mrSges(2,1) - m(3) * t104 - t10 * mrSges(5,1) - t9 * mrSges(5,2) - m(6) * t80 - t4 * mrSges(6,1) - m(7) * (pkin(5) * t4 + t80) - t2 * mrSges(7,1) - t1 * mrSges(7,2) + t121 * t33 - t120 * t3 + (t90 * t59 + mrSges(2,2)) * t65 + t77 * t34 - t124 * t97) * g(2) + (t65 * mrSges(2,1) - m(3) * t94 - t48 * mrSges(5,1) - t120 * t5 - t122 * t7 + t128 * t30 + (mrSges(2,2) + (mrSges(5,2) * t63 + t90) * t59) * t69 + t119 * t31 + t124 * t87 + t123 * (pkin(4) * t116 - t54 * t110 - t31 * t61 + t87)) * g(1) (-t123 * (t46 * t117 - t45 * t61 + t42) - t124 * t42 + t118 * t46 - t119 * t45) * g(3) + (-t123 * (t32 * t117 + t30 * t61 - t26) + t124 * t26 + t118 * t32 + t119 * t30) * g(2) + (-t123 * (t35 * t117 + t33 * t61 - t28) + t124 * t28 + t118 * t35 + t119 * t33) * g(1) (-g(1) * t33 - g(2) * t30 + g(3) * t45) * t95 (t125 * mrSges(5,1) - (t45 * t63 - t60 * t67) * mrSges(5,2) - t120 * t15 - t122 * (-t45 * t56 - t55 * t60)) * g(3) + (-t126 * mrSges(5,1) - (t48 - t116) * mrSges(5,2) + t120 * t7 - t122 * t5) * g(2) + (-t9 * mrSges(5,1) + t10 * mrSges(5,2) - t120 * t4 + t122 * t3) * g(1) + (-g(1) * t9 - g(2) * t126 + g(3) * t125) * t123 * pkin(4), t123 * (-g(1) * t34 - g(2) * t31 - g(3) * t46) -g(1) * (mrSges(7,1) * t1 - mrSges(7,2) * t2) - g(2) * ((t31 * t66 + t62 * t7) * mrSges(7,1) + (-t31 * t62 + t66 * t7) * mrSges(7,2)) - g(3) * ((-t15 * t62 + t46 * t66) * mrSges(7,1) + (-t15 * t66 - t46 * t62) * mrSges(7,2))];
taug  = t6(:);
