% Calculate Gravitation load on the joints for
% S6PRRPRR4
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% g [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [11x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,alpha2,d2,d3,d5,d6,theta1]';
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
% Datum: 2018-11-23 15:16
% Revision: 76f9d5e39f14dc242b53c0d9d3d9db48bd8f37c0
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für mechatronische Systeme, Universität Hannover

function taug = S6PRRPRR4_gravloadJ_floatb_twist_slag_vp2(qJ, g, ...
  pkin, m, mrSges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(3,1),zeros(11,1),zeros(7,1),zeros(7,3)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6PRRPRR4_gravloadJ_floatb_twist_slag_vp2: qJ has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6PRRPRR4_gravloadJ_floatb_twist_slag_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [11 1]), ...
  'S6PRRPRR4_gravloadJ_floatb_twist_slag_vp2: pkin has to be [11x1] (double)');
assert( isreal(m) && all(size(m) == [7 1]), ...
  'S6PRRPRR4_gravloadJ_floatb_twist_slag_vp2: m has to be [7x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [7,3]), ...
  'S6PRRPRR4_gravloadJ_floatb_twist_slag_vp2: mrSges has to be [7x3] (double)');

%% Symbolic Calculation
% From joint_gravload_floatb_twist_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2018-11-23 15:16:45
% EndTime: 2018-11-23 15:16:46
% DurationCPUTime: 1.12s
% Computational Cost: add. (1429->117), mult. (1818->155), div. (0->0), fcn. (1873->16), ass. (0->71)
t131 = m(6) + m(7);
t64 = sin(qJ(6));
t68 = cos(qJ(6));
t137 = t131 * (pkin(8) - pkin(9)) - t64 * mrSges(7,1) - t68 * mrSges(7,2) - mrSges(3,2) + mrSges(5,2) + mrSges(4,3) - mrSges(6,3);
t130 = mrSges(4,1) + mrSges(5,1);
t128 = mrSges(4,2) - mrSges(5,3);
t127 = -m(7) * pkin(10) + mrSges(6,2) - mrSges(7,3);
t129 = m(7) * pkin(5) + t68 * mrSges(7,1) - t64 * mrSges(7,2) + mrSges(6,1);
t100 = pkin(6) + qJ(2);
t85 = sin(t100);
t58 = t85 / 0.2e1;
t101 = pkin(6) - qJ(2);
t86 = sin(t101);
t106 = t58 - t86 / 0.2e1;
t63 = sin(pkin(11));
t103 = cos(pkin(11));
t71 = cos(qJ(2));
t94 = t103 * t71;
t43 = -t63 * t106 + t94;
t66 = sin(qJ(3));
t70 = cos(qJ(3));
t102 = sin(pkin(6));
t95 = t63 * t102;
t25 = t43 * t66 - t70 * t95;
t26 = t43 * t70 + t66 * t95;
t65 = sin(qJ(5));
t69 = cos(qJ(5));
t8 = t25 * t65 + t26 * t69;
t136 = t127 * t8 - t129 * (t25 * t69 - t26 * t65);
t109 = t63 * t71;
t40 = t103 * t106 + t109;
t75 = t103 * t102;
t23 = t40 * t66 + t70 * t75;
t24 = t40 * t70 - t66 * t75;
t4 = t23 * t65 + t24 * t69;
t135 = t127 * t4 - t129 * (t23 * t69 - t24 * t65);
t104 = cos(pkin(6));
t80 = cos(t100) / 0.2e1;
t87 = cos(t101);
t54 = t80 - t87 / 0.2e1;
t45 = -t104 * t70 - t54 * t66;
t46 = t104 * t66 - t54 * t70;
t16 = t45 * t65 + t46 * t69;
t134 = t127 * t16 - t129 * (t45 * t69 - t46 * t65);
t105 = qJ(4) * t66;
t55 = t87 / 0.2e1 + t80;
t67 = sin(qJ(2));
t39 = t103 * t55 - t63 * t67;
t112 = t39 * t70;
t126 = pkin(3) * t112 + t39 * t105;
t42 = -t103 * t67 - t63 * t55;
t111 = t42 * t70;
t125 = pkin(3) * t111 + t42 * t105;
t79 = t86 / 0.2e1;
t53 = t58 + t79;
t110 = t53 * t70;
t124 = pkin(3) * t110 + t53 * t105;
t123 = m(5) + t131;
t35 = t39 * pkin(2);
t72 = t79 - t85 / 0.2e1;
t41 = -t103 * t72 + t109;
t98 = pkin(8) * t41 + t35;
t36 = t42 * pkin(2);
t44 = t63 * t72 + t94;
t97 = pkin(8) * t44 + t36;
t52 = t53 * pkin(2);
t96 = -pkin(8) * t54 + t52;
t93 = -t23 * pkin(3) + qJ(4) * t24;
t92 = -t25 * pkin(3) + qJ(4) * t26;
t91 = -t45 * pkin(3) + qJ(4) * t46;
t1 = [(-m(2) - m(3) - m(4) - t123) * g(3) (-m(4) * t96 - m(5) * (t96 + t124) - t131 * (pkin(4) * t110 + t124 + t52) + t137 * t54) * g(3) + (-m(4) * t98 - m(5) * (t98 + t126) - t131 * (pkin(4) * t112 + t126 + t35) - t137 * t41) * g(2) + (-m(4) * t97 - m(5) * (t97 + t125) - t131 * (pkin(4) * t111 + t125 + t36) - t137 * t44) * g(1) + (t42 * g(1) + t39 * g(2) + t53 * g(3)) * ((t65 * t70 - t66 * t69) * t127 - mrSges(3,1) - t130 * t70 + t128 * t66 - (t65 * t66 + t69 * t70) * t129) (-m(5) * t91 - t131 * (-t45 * pkin(4) + t91) + t128 * t46 + t130 * t45 - t134) * g(3) + (-m(5) * t93 - t131 * (-t23 * pkin(4) + t93) + t128 * t24 + t130 * t23 - t135) * g(2) + (-m(5) * t92 - t131 * (-t25 * pkin(4) + t92) + t128 * t26 + t130 * t25 - t136) * g(1), t123 * (-g(1) * t25 - g(2) * t23 - g(3) * t45) t136 * g(1) + t135 * g(2) + t134 * g(3), -g(1) * ((t42 * t68 - t64 * t8) * mrSges(7,1) + (-t42 * t64 - t68 * t8) * mrSges(7,2)) - g(2) * ((t39 * t68 - t4 * t64) * mrSges(7,1) + (-t39 * t64 - t4 * t68) * mrSges(7,2)) - g(3) * ((-t16 * t64 + t53 * t68) * mrSges(7,1) + (-t16 * t68 - t53 * t64) * mrSges(7,2))];
taug  = t1(:);
