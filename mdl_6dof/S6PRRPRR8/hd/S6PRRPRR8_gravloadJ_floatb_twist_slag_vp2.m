% Calculate Gravitation load on the joints for
% S6PRRPRR8
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% g [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [12x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,alpha2,alpha3,d2,d3,d5,d6,theta1]';
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
% Datum: 2018-11-23 15:19
% Revision: 76f9d5e39f14dc242b53c0d9d3d9db48bd8f37c0
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für mechatronische Systeme, Universität Hannover

function taug = S6PRRPRR8_gravloadJ_floatb_twist_slag_vp2(qJ, g, ...
  pkin, m, mrSges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(3,1),zeros(12,1),zeros(7,1),zeros(7,3)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6PRRPRR8_gravloadJ_floatb_twist_slag_vp2: qJ has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6PRRPRR8_gravloadJ_floatb_twist_slag_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [12 1]), ...
  'S6PRRPRR8_gravloadJ_floatb_twist_slag_vp2: pkin has to be [12x1] (double)');
assert( isreal(m) && all(size(m) == [7 1]), ...
  'S6PRRPRR8_gravloadJ_floatb_twist_slag_vp2: m has to be [7x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [7,3]), ...
  'S6PRRPRR8_gravloadJ_floatb_twist_slag_vp2: mrSges has to be [7x3] (double)');

%% Symbolic Calculation
% From joint_gravload_floatb_twist_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2018-11-23 15:19:24
% EndTime: 2018-11-23 15:19:25
% DurationCPUTime: 1.20s
% Computational Cost: add. (3025->126), mult. (3101->182), div. (0->0), fcn. (2973->22), ass. (0->82)
t129 = m(6) + m(7);
t134 = m(5) + t129;
t93 = -qJ(4) * t134 + mrSges(4,2) - mrSges(5,3);
t70 = sin(qJ(6));
t74 = cos(qJ(6));
t82 = -t70 * mrSges(7,1) - t74 * mrSges(7,2) - pkin(10) * t129 - mrSges(4,1) + mrSges(5,2) - mrSges(6,3);
t88 = -m(7) * pkin(5) - t74 * mrSges(7,1) + t70 * mrSges(7,2) - mrSges(6,1);
t143 = m(7) * pkin(11) - mrSges(6,2) + mrSges(7,3);
t141 = pkin(3) * t134 - t82;
t114 = pkin(7) + qJ(3);
t100 = sin(t114) / 0.2e1;
t115 = pkin(7) - qJ(3);
t106 = sin(t115);
t135 = t100 - t106 / 0.2e1;
t116 = pkin(6) + qJ(2);
t104 = cos(t116) / 0.2e1;
t117 = pkin(6) - qJ(2);
t110 = cos(t117);
t58 = t104 - t110 / 0.2e1;
t76 = cos(qJ(3));
t101 = sin(t116) / 0.2e1;
t107 = sin(t117);
t95 = t101 + t107 / 0.2e1;
t140 = t135 * t95 - t58 * t76;
t57 = t101 - t107 / 0.2e1;
t66 = sin(pkin(12));
t68 = cos(pkin(12));
t77 = cos(qJ(2));
t51 = t66 * t57 - t68 * t77;
t59 = t110 / 0.2e1 + t104;
t73 = sin(qJ(2));
t97 = t59 * t66 + t68 * t73;
t139 = -t135 * t97 - t51 * t76;
t48 = t68 * t57 + t66 * t77;
t98 = t59 * t68 - t66 * t73;
t138 = t135 * t98 + t48 * t76;
t71 = sin(qJ(5));
t75 = cos(qJ(5));
t130 = t143 * t75 + t88 * t71 + t93;
t67 = sin(pkin(7));
t125 = t67 * t71;
t124 = t67 * t75;
t120 = cos(pkin(6));
t119 = sin(pkin(6));
t113 = m(4) + t134;
t69 = cos(pkin(7));
t111 = t69 * t119;
t109 = cos(t115);
t108 = cos(t114);
t103 = t109 / 0.2e1;
t102 = t108 / 0.2e1;
t96 = t103 - t108 / 0.2e1;
t87 = t96 * t119;
t86 = t103 + t102;
t85 = t102 - t109 / 0.2e1;
t83 = t100 + t106 / 0.2e1;
t81 = t85 * t119;
t80 = t83 * t119;
t79 = -mrSges(3,2) + (pkin(4) * t129 + t113 * pkin(9) + mrSges(5,1) + mrSges(4,3)) * t67;
t72 = sin(qJ(3));
t56 = t95 * pkin(2);
t47 = t120 * t69 - t67 * t95;
t46 = t97 * pkin(2);
t45 = t98 * pkin(2);
t35 = t111 * t66 + t67 * t97;
t34 = -t111 * t68 - t67 * t98;
t33 = t135 * t58 + t76 * t95;
t32 = -t58 * t86 + t72 * t95;
t29 = t120 * t96 + t140;
t28 = -t120 * t83 - t58 * t72 - t86 * t95;
t26 = t135 * t51 - t76 * t97;
t25 = -t51 * t86 - t72 * t97;
t24 = -t135 * t48 + t76 * t98;
t23 = t48 * t86 + t72 * t98;
t17 = t66 * t87 + t139;
t16 = -t51 * t72 - t66 * t80 + t86 * t97;
t14 = -t68 * t87 + t138;
t13 = t48 * t72 + t68 * t80 - t86 * t98;
t10 = t28 * t71 + t47 * t75;
t4 = t16 * t71 + t35 * t75;
t2 = t13 * t71 + t34 * t75;
t1 = [(-m(2) - m(3) - t113) * g(3) (-t95 * mrSges(3,1) - m(4) * t56 + t93 * t32 + t143 * (t125 * t58 + t32 * t75) + t88 * (-t124 * t58 + t32 * t71) + t82 * t33 + t79 * t58) * g(3) + (-t98 * mrSges(3,1) - m(4) * t45 + t143 * (-t125 * t48 + t23 * t75) + t93 * t23 + t88 * (t124 * t48 + t23 * t71) + t82 * t24 - t79 * t48) * g(2) + (t97 * mrSges(3,1) + m(4) * t46 + t143 * (t125 * t51 + t25 * t75) + t93 * t25 + t88 * (-t124 * t51 + t25 * t71) + t82 * t26 + t79 * t51) * g(1) + (-g(2) * (t24 * pkin(3) + t45) - g(1) * (t26 * pkin(3) - t46) - g(3) * (t33 * pkin(3) + t56)) * t134 (t130 * (-t120 * t85 + t140) + t141 * t28) * g(3) + (t130 * (t68 * t81 + t138) + t141 * t13) * g(2) + (t130 * (-t66 * t81 + t139) + t141 * t16) * g(1), t134 * (-g(1) * t16 - g(2) * t13 - g(3) * t28) (t88 * (t28 * t75 - t47 * t71) - t143 * t10) * g(3) + (-t143 * t2 + t88 * (t13 * t75 - t34 * t71)) * g(2) + (-t143 * t4 + t88 * (t16 * t75 - t35 * t71)) * g(1), -g(1) * ((t17 * t74 - t4 * t70) * mrSges(7,1) + (-t17 * t70 - t4 * t74) * mrSges(7,2)) - g(2) * ((t14 * t74 - t2 * t70) * mrSges(7,1) + (-t14 * t70 - t2 * t74) * mrSges(7,2)) - g(3) * ((-t10 * t70 + t29 * t74) * mrSges(7,1) + (-t10 * t74 - t29 * t70) * mrSges(7,2))];
taug  = t1(:);
