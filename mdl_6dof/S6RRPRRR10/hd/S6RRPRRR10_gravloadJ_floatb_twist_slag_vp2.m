% Calculate Gravitation load on the joints for
% S6RRPRRR10
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% g [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [12x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,alpha2,d1,d2,d4,d5,d6,theta3]';
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
% Datum: 2018-11-23 17:28
% Revision: 76f9d5e39f14dc242b53c0d9d3d9db48bd8f37c0
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für mechatronische Systeme, Universität Hannover

function taug = S6RRPRRR10_gravloadJ_floatb_twist_slag_vp2(qJ, g, ...
  pkin, m, mrSges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(3,1),zeros(12,1),zeros(7,1),zeros(7,3)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRPRRR10_gravloadJ_floatb_twist_slag_vp2: qJ has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6RRPRRR10_gravloadJ_floatb_twist_slag_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [12 1]), ...
  'S6RRPRRR10_gravloadJ_floatb_twist_slag_vp2: pkin has to be [12x1] (double)');
assert( isreal(m) && all(size(m) == [7 1]), ...
  'S6RRPRRR10_gravloadJ_floatb_twist_slag_vp2: m has to be [7x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [7,3]), ...
  'S6RRPRRR10_gravloadJ_floatb_twist_slag_vp2: mrSges has to be [7x3] (double)');

%% Symbolic Calculation
% From joint_gravload_floatb_twist_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2018-11-23 17:27:30
% EndTime: 2018-11-23 17:27:31
% DurationCPUTime: 1.20s
% Computational Cost: add. (1569->135), mult. (1609->173), div. (0->0), fcn. (1583->18), ass. (0->63)
t67 = cos(qJ(5));
t51 = pkin(5) * t67 + pkin(4);
t58 = qJ(5) + qJ(6);
t54 = sin(t58);
t55 = cos(t58);
t64 = sin(qJ(5));
t115 = m(6) * pkin(4) + m(7) * t51 + mrSges(6,1) * t67 + mrSges(7,1) * t55 - mrSges(6,2) * t64 - mrSges(7,2) * t54 + mrSges(5,1);
t74 = mrSges(5,2) + m(7) * (-pkin(11) - pkin(10)) - mrSges(7,3) - m(6) * pkin(10) - mrSges(6,3);
t59 = sin(pkin(12));
t61 = cos(pkin(12));
t121 = -m(4) * pkin(2) - t61 * mrSges(4,1) + t59 * mrSges(4,2) - mrSges(3,1);
t57 = pkin(12) + qJ(4);
t52 = sin(t57);
t53 = cos(t57);
t114 = t115 * t53 - t74 * t52 - t121;
t81 = m(4) * qJ(3) - mrSges(3,2) + mrSges(4,3) + mrSges(5,3);
t122 = t54 * mrSges(7,1) + t67 * mrSges(6,2) + t55 * mrSges(7,2) + t81;
t112 = m(7) * pkin(5);
t118 = mrSges(6,1) + t112;
t117 = m(5) + m(6) + m(7);
t101 = t64 * t112;
t113 = -t64 * mrSges(6,1) - t101 - t122;
t108 = cos(qJ(1));
t99 = pkin(6) + qJ(2);
t87 = sin(t99) / 0.2e1;
t100 = pkin(6) - qJ(2);
t91 = sin(t100);
t38 = t87 - t91 / 0.2e1;
t66 = sin(qJ(1));
t68 = cos(qJ(2));
t79 = t108 * t38 + t66 * t68;
t60 = sin(pkin(6));
t96 = t60 * t108;
t12 = -t52 * t96 + t53 * t79;
t65 = sin(qJ(2));
t88 = cos(t99) / 0.2e1;
t92 = cos(t100);
t72 = t92 / 0.2e1 + t88;
t27 = -t108 * t72 + t65 * t66;
t111 = (-t12 * t54 + t27 * t55) * mrSges(7,1) + (-t12 * t55 - t27 * t54) * mrSges(7,2);
t107 = t60 * t66;
t78 = t108 * t68 - t66 * t38;
t16 = t107 * t52 + t53 * t78;
t30 = t108 * t65 + t66 * t72;
t5 = -t16 * t54 + t30 * t55;
t6 = t16 * t55 + t30 * t54;
t110 = t5 * mrSges(7,1) - t6 * mrSges(7,2);
t39 = t88 - t92 / 0.2e1;
t62 = cos(pkin(6));
t22 = -t39 * t53 + t52 * t62;
t37 = t87 + t91 / 0.2e1;
t109 = (-t22 * t54 - t37 * t55) * mrSges(7,1) + (-t22 * t55 + t37 * t54) * mrSges(7,2);
t102 = t108 * pkin(1) + pkin(8) * t107;
t94 = -t66 * pkin(1) + pkin(8) * t96;
t11 = -t52 * t79 - t53 * t96;
t90 = t59 * t96;
t50 = pkin(3) * t61 + pkin(2);
t63 = -pkin(9) - qJ(3);
t89 = t59 * pkin(3) * t107 - t30 * t63 + t50 * t78 + t102;
t7 = -t16 * t64 + t30 * t67;
t15 = -t107 * t53 + t52 * t78;
t8 = t16 * t67 + t30 * t64;
t1 = [(-t108 * mrSges(2,1) - m(5) * t89 - t16 * mrSges(5,1) - m(6) * (pkin(4) * t16 + t89) - t8 * mrSges(6,1) - t7 * mrSges(6,2) - m(7) * (t16 * t51 + t89) - t6 * mrSges(7,1) - t5 * mrSges(7,2) + t121 * t78 + (mrSges(2,2) + (-mrSges(4,1) * t59 - mrSges(4,2) * t61 - mrSges(3,3)) * t60) * t66 + (-t81 - t101) * t30 + t74 * t15 + (-m(3) - m(4)) * t102) * g(2) + (t66 * mrSges(2,1) + t108 * mrSges(2,2) - m(3) * t94 + t79 * mrSges(3,1) - mrSges(3,3) * t96 - m(4) * (-pkin(2) * t79 + t94) - (-t61 * t79 + t90) * mrSges(4,1) - (t59 * t79 + t61 * t96) * mrSges(4,2) + t115 * t12 + (t118 * t64 + t122) * t27 + t74 * t11 + t117 * (-pkin(3) * t90 - t27 * t63 + t50 * t79 - t94)) * g(1) (-t117 * (t37 * t50 + t39 * t63) - t113 * t39 - t114 * t37) * g(3) + (-t117 * (-t27 * t50 - t63 * t79) + t113 * t79 + t114 * t27) * g(2) + (-t117 * (-t30 * t50 - t63 * t78) + t113 * t78 + t114 * t30) * g(1) (-g(1) * t30 - g(2) * t27 + g(3) * t37) * (m(4) + t117) (t74 * t22 - t115 * (t39 * t52 + t53 * t62)) * g(3) + (-t115 * t11 + t74 * t12) * g(2) + (t115 * t15 + t74 * t16) * g(1) (-(-t22 * t67 + t37 * t64) * mrSges(6,2) - t109 - t118 * (-t22 * t64 - t37 * t67)) * g(3) + (-(-t12 * t67 - t27 * t64) * mrSges(6,2) - t111 - t118 * (-t12 * t64 + t27 * t67)) * g(2) + (mrSges(6,2) * t8 - t118 * t7 - t110) * g(1), -g(1) * t110 - g(2) * t111 - g(3) * t109];
taug  = t1(:);
