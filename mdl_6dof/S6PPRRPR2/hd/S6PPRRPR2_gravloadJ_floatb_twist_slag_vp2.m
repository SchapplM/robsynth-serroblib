% Calculate Gravitation load on the joints for
% S6PPRRPR2
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% g [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [12x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,alpha2,alpha3,d3,d4,d6,theta1,theta2]';
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
% Datum: 2018-11-23 14:50
% Revision: 76f9d5e39f14dc242b53c0d9d3d9db48bd8f37c0
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für mechatronische Systeme, Universität Hannover

function taug = S6PPRRPR2_gravloadJ_floatb_twist_slag_vp2(qJ, g, ...
  pkin, m, mrSges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(3,1),zeros(12,1),zeros(7,1),zeros(7,3)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6PPRRPR2_gravloadJ_floatb_twist_slag_vp2: qJ has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6PPRRPR2_gravloadJ_floatb_twist_slag_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [12 1]), ...
  'S6PPRRPR2_gravloadJ_floatb_twist_slag_vp2: pkin has to be [12x1] (double)');
assert( isreal(m) && all(size(m) == [7 1]), ...
  'S6PPRRPR2_gravloadJ_floatb_twist_slag_vp2: m has to be [7x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [7,3]), ...
  'S6PPRRPR2_gravloadJ_floatb_twist_slag_vp2: mrSges has to be [7x3] (double)');

%% Symbolic Calculation
% From joint_gravload_floatb_twist_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2018-11-23 14:50:17
% EndTime: 2018-11-23 14:50:18
% DurationCPUTime: 0.78s
% Computational Cost: add. (2394->108), mult. (2507->151), div. (0->0), fcn. (2447->22), ass. (0->85)
t120 = m(6) + m(7);
t128 = mrSges(5,1) - mrSges(6,2) + mrSges(7,3);
t49 = sin(qJ(6));
t52 = cos(qJ(6));
t127 = t49 * mrSges(7,1) + t52 * mrSges(7,2) - mrSges(5,2) + mrSges(6,3);
t124 = m(7) * pkin(10) + pkin(4) * t120 + t128;
t91 = pkin(7) + qJ(3);
t77 = sin(t91) / 0.2e1;
t92 = pkin(7) - qJ(3);
t82 = sin(t92);
t116 = t77 - t82 / 0.2e1;
t90 = pkin(6) - pkin(12);
t74 = cos(t90) / 0.2e1;
t89 = pkin(6) + pkin(12);
t81 = cos(t89);
t44 = t74 - t81 / 0.2e1;
t54 = cos(qJ(3));
t73 = sin(t89) / 0.2e1;
t80 = sin(t90);
t62 = t73 + t80 / 0.2e1;
t123 = t116 * t62 + t44 * t54;
t43 = t73 - t80 / 0.2e1;
t48 = cos(pkin(12));
t94 = sin(pkin(11));
t96 = cos(pkin(11));
t38 = -t43 * t94 + t48 * t96;
t63 = t74 + t81 / 0.2e1;
t93 = sin(pkin(12));
t58 = t63 * t94 + t93 * t96;
t122 = -t116 * t58 + t38 * t54;
t37 = t43 * t96 + t48 * t94;
t57 = -t63 * t96 + t93 * t94;
t121 = -t116 * t57 + t37 * t54;
t51 = sin(qJ(3));
t47 = sin(pkin(6));
t64 = t77 + t82 / 0.2e1;
t60 = t47 * t64;
t83 = cos(t91);
t78 = t83 / 0.2e1;
t84 = cos(t92);
t79 = t84 / 0.2e1;
t67 = t79 + t78;
t16 = t37 * t51 + t57 * t67 + t60 * t96;
t53 = cos(qJ(4));
t108 = t16 * t53;
t50 = sin(qJ(4));
t99 = qJ(5) * t50;
t119 = -pkin(4) * t108 - t16 * t99;
t19 = t38 * t51 + t58 * t67 - t60 * t94;
t107 = t19 * t53;
t118 = -pkin(4) * t107 - t19 * t99;
t98 = cos(pkin(6));
t25 = t44 * t51 - t62 * t67 - t64 * t98;
t106 = t25 * t53;
t117 = -pkin(4) * t106 - t25 * t99;
t114 = m(3) + m(4) + m(5) + t120;
t112 = -t120 * qJ(5) - t127;
t111 = t127 * t50 + t128 * t53 + mrSges(4,1);
t110 = -m(7) * (pkin(5) + pkin(9)) - t52 * mrSges(7,1) + t49 * mrSges(7,2) - mrSges(6,1) + mrSges(4,2) - mrSges(5,3);
t97 = cos(pkin(7));
t95 = sin(pkin(7));
t14 = t16 * pkin(3);
t66 = t78 - t84 / 0.2e1;
t61 = t47 * t66;
t18 = t61 * t96 + t121;
t88 = pkin(9) * t18 - t14;
t15 = t19 * pkin(3);
t21 = -t61 * t94 + t122;
t87 = pkin(9) * t21 - t15;
t24 = t25 * pkin(3);
t27 = -t66 * t98 + t123;
t86 = pkin(9) * t27 - t24;
t85 = t47 * t97;
t72 = t79 - t83 / 0.2e1;
t69 = t47 * t72;
t59 = -t62 * t95 + t97 * t98;
t56 = t58 * t95 + t85 * t94;
t55 = t57 * t95 - t85 * t96;
t26 = t72 * t98 + t123;
t20 = t69 * t94 + t122;
t17 = -t69 * t96 + t121;
t8 = t26 * t50 - t53 * t59;
t5 = t20 * t50 - t53 * t56;
t3 = t17 * t50 - t53 * t55;
t1 = [(-m(2) - t114) * g(3) (-t98 * g(3) + (-g(1) * t94 + g(2) * t96) * t47) * t114 (-m(5) * t86 - m(6) * (t86 + t117) - m(7) * (-pkin(10) * t106 + t117 - t24) + t110 * t27 + t111 * t25) * g(3) + (-m(5) * t88 - m(6) * (t88 + t119) - m(7) * (-pkin(10) * t108 + t119 - t14) + t110 * t18 + t111 * t16) * g(2) + (-m(5) * t87 - m(6) * (t87 + t118) - m(7) * (-pkin(10) * t107 + t118 - t15) + t110 * t21 + t111 * t19) * g(1) (t112 * (t26 * t53 + t50 * t59) + t124 * t8) * g(3) + (t112 * (t17 * t53 + t50 * t55) + t124 * t3) * g(2) + (t112 * (t20 * t53 + t50 * t56) + t124 * t5) * g(1), t120 * (-g(1) * t5 - g(2) * t3 - g(3) * t8) -g(1) * ((-t19 * t49 + t5 * t52) * mrSges(7,1) + (-t19 * t52 - t49 * t5) * mrSges(7,2)) - g(2) * ((-t16 * t49 + t3 * t52) * mrSges(7,1) + (-t16 * t52 - t3 * t49) * mrSges(7,2)) - g(3) * ((-t25 * t49 + t52 * t8) * mrSges(7,1) + (-t25 * t52 - t49 * t8) * mrSges(7,2))];
taug  = t1(:);
