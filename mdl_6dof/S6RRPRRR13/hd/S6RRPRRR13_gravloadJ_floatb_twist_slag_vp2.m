% Calculate Gravitation load on the joints for
% S6RRPRRR13
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% g [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [11x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,alpha2,d1,d2,d4,d5,d6]';
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
% Datum: 2018-11-23 17:30
% Revision: 76f9d5e39f14dc242b53c0d9d3d9db48bd8f37c0
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für mechatronische Systeme, Universität Hannover

function taug = S6RRPRRR13_gravloadJ_floatb_twist_slag_vp2(qJ, g, ...
  pkin, m, mrSges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(3,1),zeros(11,1),zeros(7,1),zeros(7,3)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRPRRR13_gravloadJ_floatb_twist_slag_vp2: qJ has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6RRPRRR13_gravloadJ_floatb_twist_slag_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [11 1]), ...
  'S6RRPRRR13_gravloadJ_floatb_twist_slag_vp2: pkin has to be [11x1] (double)');
assert( isreal(m) && all(size(m) == [7 1]), ...
  'S6RRPRRR13_gravloadJ_floatb_twist_slag_vp2: m has to be [7x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [7,3]), ...
  'S6RRPRRR13_gravloadJ_floatb_twist_slag_vp2: mrSges has to be [7x3] (double)');

%% Symbolic Calculation
% From joint_gravload_floatb_twist_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2018-11-23 17:29:49
% EndTime: 2018-11-23 17:29:50
% DurationCPUTime: 1.01s
% Computational Cost: add. (1387->123), mult. (1620->156), div. (0->0), fcn. (1588->16), ass. (0->67)
t120 = m(5) + m(6);
t94 = m(7) + t120;
t116 = m(4) + t94;
t111 = qJ(3) * t116 - mrSges(3,2) + mrSges(4,3);
t49 = qJ(5) + qJ(6);
t46 = sin(t49);
t47 = cos(t49);
t56 = cos(qJ(5));
t119 = t46 * mrSges(7,1) + t56 * mrSges(6,2) + t47 * mrSges(7,2);
t52 = sin(qJ(5));
t103 = pkin(5) * t52;
t91 = mrSges(4,2) - mrSges(3,1) - mrSges(5,3);
t117 = pkin(2) * t116 + t52 * mrSges(6,1) - t91 + t120 * pkin(9) + m(7) * (pkin(9) + t103) + t119;
t45 = pkin(5) * t56 + pkin(4);
t114 = m(6) * pkin(4) + m(7) * t45 + mrSges(6,1) * t56 + mrSges(7,1) * t47 - mrSges(6,2) * t52 - mrSges(7,2) * t46 + mrSges(5,1);
t118 = m(6) * pkin(10) - m(7) * (-pkin(11) - pkin(10)) - mrSges(5,2) + mrSges(6,3) + mrSges(7,3);
t115 = m(7) * pkin(5) + mrSges(6,1);
t112 = t94 * pkin(9) - t91;
t53 = sin(qJ(4));
t57 = cos(qJ(4));
t107 = -t114 * t53 + t118 * t57 - t111;
t50 = sin(pkin(6));
t55 = sin(qJ(1));
t101 = t50 * t55;
t54 = sin(qJ(2));
t59 = cos(qJ(1));
t93 = pkin(6) - qJ(2);
t79 = cos(t93) / 0.2e1;
t92 = pkin(6) + qJ(2);
t83 = cos(t92);
t64 = t79 + t83 / 0.2e1;
t29 = t59 * t54 + t55 * t64;
t12 = t101 * t57 + t29 * t53;
t81 = sin(t92);
t77 = t81 / 0.2e1;
t82 = sin(t93);
t70 = t77 - t82 / 0.2e1;
t58 = cos(qJ(2));
t97 = t59 * t58;
t30 = -t55 * t70 + t97;
t5 = -t12 * t46 + t30 * t47;
t6 = t12 * t47 + t30 * t46;
t105 = t5 * mrSges(7,1) - t6 * mrSges(7,2);
t100 = t50 * t59;
t26 = t54 * t55 - t59 * t64;
t15 = t100 * t57 - t26 * t53;
t99 = t55 * t58;
t27 = t59 * t70 + t99;
t104 = (t15 * t46 + t27 * t47) * mrSges(7,1) + (t15 * t47 - t27 * t46) * mrSges(7,2);
t78 = t82 / 0.2e1;
t36 = t77 + t78;
t51 = cos(pkin(6));
t25 = -t36 * t53 + t51 * t57;
t37 = t79 - t83 / 0.2e1;
t102 = (-t25 * t46 + t37 * t47) * mrSges(7,1) + (-t25 * t47 - t37 * t46) * mrSges(7,2);
t95 = t59 * pkin(1) + pkin(8) * t101;
t90 = t30 * pkin(2) + t95;
t85 = -t55 * pkin(1) + pkin(8) * t100;
t84 = pkin(3) * t101 + t90;
t80 = -t27 * pkin(2) + t85;
t76 = mrSges(2,2) + (-mrSges(4,1) - mrSges(3,3)) * t50;
t7 = -t12 * t52 + t30 * t56;
t71 = t78 - t81 / 0.2e1;
t13 = t100 * t53 + t26 * t57;
t11 = t101 * t53 - t29 * t57;
t8 = t12 * t56 + t30 * t52;
t1 = [(-t59 * mrSges(2,1) - m(3) * t95 - m(4) * t90 - m(5) * t84 - t12 * mrSges(5,1) - m(6) * (pkin(4) * t12 + t84) - t8 * mrSges(6,1) - t7 * mrSges(6,2) - m(7) * (t12 * t45 + t84) - t6 * mrSges(7,1) - t5 * mrSges(7,2) - t111 * t29 + t76 * t55 - t118 * t11 + (-m(7) * t103 - t112) * t30) * g(2) + (t55 * mrSges(2,1) - m(3) * t85 - m(4) * t80 + t111 * t26 + t76 * t59 - t118 * t13 - t114 * t15 + (t115 * t52 + t112 + t119) * t27 + t94 * (-pkin(3) * t100 - t80)) * g(1) (t107 * t37 - t117 * t36) * g(3) + (t107 * (-t59 * t71 + t99) + t117 * t26) * g(2) + (t107 * (t55 * t71 + t97) + t117 * t29) * g(1) (-g(1) * t29 - g(2) * t26 + g(3) * t36) * t116 (-t118 * t25 - t114 * (-t36 * t57 - t51 * t53)) * g(3) + (-t114 * t13 + t118 * t15) * g(2) + (t114 * t11 - t118 * t12) * g(1) (-(-t25 * t56 - t37 * t52) * mrSges(6,2) - t102 - t115 * (-t25 * t52 + t37 * t56)) * g(3) + (-(t15 * t56 - t27 * t52) * mrSges(6,2) - t104 - t115 * (t15 * t52 + t27 * t56)) * g(2) + (t8 * mrSges(6,2) - t115 * t7 - t105) * g(1), -g(1) * t105 - g(2) * t104 - g(3) * t102];
taug  = t1(:);
