% Calculate Gravitation load on the joints for
% S6RRPRRR2
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% g [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [11x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d2,d4,d5,d6,theta3]';
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
% Datum: 2018-11-23 17:21
% Revision: 76f9d5e39f14dc242b53c0d9d3d9db48bd8f37c0
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für mechatronische Systeme, Universität Hannover

function taug = S6RRPRRR2_gravloadJ_floatb_twist_slag_vp2(qJ, g, ...
  pkin, m, mrSges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(3,1),zeros(11,1),zeros(7,1),zeros(7,3)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRPRRR2_gravloadJ_floatb_twist_slag_vp2: qJ has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6RRPRRR2_gravloadJ_floatb_twist_slag_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [11 1]), ...
  'S6RRPRRR2_gravloadJ_floatb_twist_slag_vp2: pkin has to be [11x1] (double)');
assert( isreal(m) && all(size(m) == [7 1]), ...
  'S6RRPRRR2_gravloadJ_floatb_twist_slag_vp2: m has to be [7x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [7,3]), ...
  'S6RRPRRR2_gravloadJ_floatb_twist_slag_vp2: mrSges has to be [7x3] (double)');

%% Symbolic Calculation
% From joint_gravload_floatb_twist_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2018-11-23 17:21:22
% EndTime: 2018-11-23 17:21:22
% DurationCPUTime: 0.75s
% Computational Cost: add. (557->117), mult. (493->133), div. (0->0), fcn. (429->12), ass. (0->71)
t49 = qJ(5) + qJ(6);
t44 = sin(t49);
t45 = cos(t49);
t51 = sin(qJ(5));
t54 = cos(qJ(5));
t137 = -mrSges(6,1) * t54 - mrSges(7,1) * t45 + mrSges(6,2) * t51 + mrSges(7,2) * t44;
t134 = -mrSges(6,3) - mrSges(7,3);
t48 = qJ(2) + pkin(11);
t43 = qJ(4) + t48;
t37 = sin(t43);
t38 = cos(t43);
t39 = pkin(5) * t54 + pkin(4);
t57 = -pkin(10) - pkin(9);
t133 = m(7) * t38 * t57 + (m(7) * t39 - t137) * t37;
t114 = m(7) * pkin(5);
t41 = sin(t48);
t42 = cos(t48);
t52 = sin(qJ(2));
t55 = cos(qJ(2));
t131 = -t55 * mrSges(3,1) - t42 * mrSges(4,1) + t52 * mrSges(3,2) + t41 * mrSges(4,2);
t53 = sin(qJ(1));
t56 = cos(qJ(1));
t129 = g(1) * t56 + g(2) * t53;
t128 = -t38 * mrSges(5,1) + (mrSges(5,2) + t134) * t37;
t121 = -t37 * t57 + t38 * t39;
t125 = t38 * pkin(4) + t37 * pkin(9);
t127 = -m(6) * t125 - m(7) * t121;
t126 = t51 * t114;
t98 = t38 * t56;
t124 = t133 * t56 + t134 * t98;
t99 = t38 * t53;
t123 = t133 * t53 + t134 * t99;
t122 = mrSges(6,1) + t114;
t119 = m(5) + m(6) + m(7);
t118 = t137 * t38 + t128;
t50 = -qJ(3) - pkin(7);
t116 = -m(3) * pkin(7) + m(4) * t50 + mrSges(2,2) - mrSges(3,3) - mrSges(4,3) - mrSges(5,3);
t46 = t55 * pkin(2);
t115 = mrSges(2,1) + m(4) * (t46 + pkin(1)) + m(3) * pkin(1) - t128 - t131;
t94 = t45 * t56;
t97 = t44 * t53;
t5 = t38 * t97 + t94;
t95 = t45 * t53;
t96 = t44 * t56;
t6 = -t38 * t95 + t96;
t113 = -t5 * mrSges(7,1) + t6 * mrSges(7,2);
t7 = -t38 * t96 + t95;
t8 = t38 * t94 + t97;
t112 = t7 * mrSges(7,1) - t8 * mrSges(7,2);
t111 = pkin(2) * t52;
t110 = pkin(4) * t37;
t106 = g(3) * t37;
t93 = t51 * t53;
t92 = t51 * t56;
t91 = t53 * t54;
t90 = t54 * t56;
t88 = pkin(3) * t42 + t46;
t78 = pkin(9) * t99 - t53 * t110;
t77 = pkin(9) * t98 - t56 * t110;
t71 = mrSges(5,1) * t37 + mrSges(5,2) * t38;
t70 = -mrSges(7,1) * t44 - mrSges(7,2) * t45;
t11 = -t38 * t92 + t91;
t9 = t38 * t93 + t90;
t47 = -pkin(8) + t50;
t23 = -pkin(3) * t41 - t111;
t22 = pkin(1) + t88;
t17 = t56 * t23;
t16 = t53 * t23;
t12 = t38 * t90 + t93;
t10 = -t38 * t91 + t92;
t1 = [(-t93 * t114 - t12 * mrSges(6,1) - t8 * mrSges(7,1) - t11 * mrSges(6,2) - t7 * mrSges(7,2) - t119 * (t56 * t22 - t53 * t47) + t116 * t53 + (-t115 + t127) * t56) * g(2) + (-t10 * mrSges(6,1) - t6 * mrSges(7,1) - t9 * mrSges(6,2) - t5 * mrSges(7,2) + (t119 * t47 + t116 - t126) * t56 + (m(5) * t22 - m(6) * (-t22 - t125) - m(7) * (-t22 - t121) + t115) * t53) * g(1) (-m(6) * (t16 + t78) - m(7) * t16 + t123) * g(2) + (-m(6) * (t17 + t77) - m(7) * t17 + t124) * g(1) + (-m(4) * t46 - m(5) * t88 - m(6) * (t88 + t125) - m(7) * (t121 + t88) + t118 + t131) * g(3) + t129 * (m(4) * t111 - m(5) * t23 + mrSges(3,1) * t52 + mrSges(4,1) * t41 + mrSges(3,2) * t55 + mrSges(4,2) * t42 + t71) (-g(1) * t53 + g(2) * t56) * (m(4) + t119) t129 * t71 + (-m(6) * t78 + t123) * g(2) + (-m(6) * t77 + t124) * g(1) + (t118 + t127) * g(3) (mrSges(6,1) * t51 + mrSges(6,2) * t54 + t126 - t70) * t106 + (-t10 * mrSges(6,2) + t122 * t9 - t113) * g(2) + (t12 * mrSges(6,2) - t122 * t11 - t112) * g(1), -g(1) * t112 - g(2) * t113 - t70 * t106];
taug  = t1(:);
