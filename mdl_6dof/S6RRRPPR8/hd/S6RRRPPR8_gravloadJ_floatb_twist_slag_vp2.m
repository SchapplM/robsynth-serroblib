% Calculate Gravitation load on the joints for
% S6RRRPPR8
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% g [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [10x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,alpha2,d1,d2,d3,d6]';
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
% Datum: 2018-11-23 17:38
% Revision: 76f9d5e39f14dc242b53c0d9d3d9db48bd8f37c0
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für mechatronische Systeme, Universität Hannover

function taug = S6RRRPPR8_gravloadJ_floatb_twist_slag_vp2(qJ, g, ...
  pkin, m, mrSges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(3,1),zeros(10,1),zeros(7,1),zeros(7,3)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRRPPR8_gravloadJ_floatb_twist_slag_vp2: qJ has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6RRRPPR8_gravloadJ_floatb_twist_slag_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S6RRRPPR8_gravloadJ_floatb_twist_slag_vp2: pkin has to be [10x1] (double)');
assert( isreal(m) && all(size(m) == [7 1]), ...
  'S6RRRPPR8_gravloadJ_floatb_twist_slag_vp2: m has to be [7x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [7,3]), ...
  'S6RRRPPR8_gravloadJ_floatb_twist_slag_vp2: mrSges has to be [7x3] (double)');

%% Symbolic Calculation
% From joint_gravload_floatb_twist_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2018-11-23 17:38:30
% EndTime: 2018-11-23 17:38:31
% DurationCPUTime: 1.18s
% Computational Cost: add. (1386->151), mult. (1703->178), div. (0->0), fcn. (1677->14), ass. (0->87)
t133 = m(6) + m(7);
t134 = mrSges(3,2) - mrSges(5,2) - mrSges(4,3) + mrSges(6,3) - t133 * (pkin(9) - qJ(5));
t54 = sin(qJ(6));
t58 = cos(qJ(6));
t123 = t54 * mrSges(7,1) + t58 * mrSges(7,2) + t134;
t135 = mrSges(4,2) - mrSges(6,1) - mrSges(5,3);
t120 = m(7) * pkin(10) + mrSges(4,1) + mrSges(5,1) - mrSges(6,2) + mrSges(7,3);
t122 = -m(7) * (qJ(4) + pkin(5)) + t135;
t73 = t58 * mrSges(7,1) - t54 * mrSges(7,2);
t125 = t122 - t73;
t55 = sin(qJ(3));
t105 = qJ(4) * t55;
t116 = cos(qJ(1));
t56 = sin(qJ(2));
t57 = sin(qJ(1));
t101 = pkin(6) + qJ(2);
t78 = cos(t101) / 0.2e1;
t102 = pkin(6) - qJ(2);
t85 = cos(t102);
t63 = t85 / 0.2e1 + t78;
t31 = -t116 * t63 + t56 * t57;
t59 = cos(qJ(3));
t113 = t31 * t59;
t130 = -pkin(3) * t113 - t31 * t105;
t34 = t116 * t56 + t57 * t63;
t112 = t34 * t59;
t129 = -pkin(3) * t112 - t34 * t105;
t83 = sin(t101);
t49 = t83 / 0.2e1;
t84 = sin(t102);
t77 = t84 / 0.2e1;
t43 = t49 + t77;
t111 = t43 * t59;
t128 = pkin(3) * t111 + t43 * t105;
t124 = mrSges(3,1) + t120 * t59 + (m(7) * pkin(5) - t135 + t73) * t55;
t118 = pkin(9) * t31;
t117 = pkin(9) * t34;
t115 = t31 * t54;
t114 = t31 * t58;
t60 = cos(qJ(2));
t110 = t57 * t60;
t107 = t49 - t84 / 0.2e1;
t103 = sin(pkin(6));
t91 = t57 * t103;
t106 = t116 * pkin(1) + pkin(8) * t91;
t104 = cos(pkin(6));
t98 = t116 * t60;
t35 = -t107 * t57 + t98;
t100 = t35 * pkin(2) + t106;
t79 = t116 * t103;
t97 = -pkin(1) * t57 + pkin(8) * t79;
t25 = t31 * pkin(2);
t62 = t77 - t83 / 0.2e1;
t33 = -t116 * t62 + t110;
t96 = pkin(9) * t33 - t25;
t27 = t34 * pkin(2);
t36 = t57 * t62 + t98;
t95 = pkin(9) * t36 - t27;
t42 = t43 * pkin(2);
t44 = t78 - t85 / 0.2e1;
t94 = -pkin(9) * t44 + t42;
t32 = t107 * t116 + t110;
t12 = t32 * t59 - t55 * t79;
t11 = t32 * t55 + t59 * t79;
t4 = t11 * pkin(3);
t93 = qJ(4) * t12 - t4;
t16 = t35 * t59 + t55 * t91;
t15 = t35 * t55 - t59 * t91;
t8 = t15 * pkin(3);
t92 = qJ(4) * t16 - t8;
t29 = -t104 * t59 - t44 * t55;
t24 = t29 * pkin(3);
t30 = t104 * t55 - t44 * t59;
t90 = qJ(4) * t30 - t24;
t89 = t16 * pkin(3) + t100;
t82 = -t32 * pkin(2) + t97;
t71 = -pkin(3) * t12 + t82;
t70 = qJ(4) * t15 + t89;
t66 = -qJ(4) * t11 + t71;
t23 = t29 * pkin(4);
t9 = t16 * pkin(4);
t7 = t15 * pkin(4);
t5 = t12 * pkin(4);
t3 = t11 * pkin(4);
t2 = t15 * t58 - t34 * t54;
t1 = -t15 * t54 - t34 * t58;
t6 = [(-t116 * mrSges(2,1) + t57 * mrSges(2,2) - m(3) * t106 - t35 * mrSges(3,1) - mrSges(3,3) * t91 - m(4) * (t100 + t117) - m(5) * (t70 + t117) - m(6) * (t70 + t9) - m(7) * (t89 + t9) - t2 * mrSges(7,1) - t1 * mrSges(7,2) + t122 * t15 + t134 * t34 - t120 * t16) * g(2) + (t57 * mrSges(2,1) + t116 * mrSges(2,2) - m(3) * t97 + t32 * mrSges(3,1) - mrSges(3,3) * t79 - m(4) * (t82 - t118) - m(5) * (t66 - t118) - m(6) * (-t5 + t66) - m(7) * (-t5 + t71) - t115 * mrSges(7,1) - t114 * mrSges(7,2) - t125 * t11 - t134 * t31 + t120 * t12) * g(1) (-m(4) * t94 - m(5) * (t94 + t128) - t133 * (pkin(4) * t111 + t128 + t42) - t123 * t44 - t124 * t43) * g(3) + (-m(4) * t96 - m(5) * (t96 + t130) - t133 * (-pkin(4) * t113 + t130 - t25) + t123 * t33 + t124 * t31) * g(2) + (-m(4) * t95 - m(5) * (t95 + t129) - t133 * (-pkin(4) * t112 + t129 - t27) + t123 * t36 + t124 * t34) * g(1) (-m(5) * t90 - m(6) * (-t23 + t90) - m(7) * (-t23 - t24) + t125 * t30 + t120 * t29) * g(3) + (-m(5) * t93 - m(6) * (-t3 + t93) - m(7) * (-t3 - t4) + t125 * t12 + t120 * t11) * g(2) + (-m(5) * t92 - m(6) * (-t7 + t92) - m(7) * (-t7 - t8) + t125 * t16 + t120 * t15) * g(1) (m(5) + t133) * (-g(1) * t15 - g(2) * t11 - g(3) * t29) t133 * (g(1) * t34 + g(2) * t31 - g(3) * t43) -g(1) * (mrSges(7,1) * t1 - mrSges(7,2) * t2) - g(2) * ((-t11 * t54 - t114) * mrSges(7,1) + (-t11 * t58 + t115) * mrSges(7,2)) - g(3) * ((-t29 * t54 + t43 * t58) * mrSges(7,1) + (-t29 * t58 - t43 * t54) * mrSges(7,2))];
taug  = t6(:);
