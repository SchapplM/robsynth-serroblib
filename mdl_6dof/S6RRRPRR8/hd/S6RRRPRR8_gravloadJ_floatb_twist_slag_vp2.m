% Calculate Gravitation load on the joints for
% S6RRRPRR8
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% g [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [12x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,alpha2,d1,d2,d3,d5,d6,theta4]';
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
% Datum: 2018-11-23 17:56
% Revision: 76f9d5e39f14dc242b53c0d9d3d9db48bd8f37c0
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für mechatronische Systeme, Universität Hannover

function taug = S6RRRPRR8_gravloadJ_floatb_twist_slag_vp2(qJ, g, ...
  pkin, m, mrSges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(3,1),zeros(12,1),zeros(7,1),zeros(7,3)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRRPRR8_gravloadJ_floatb_twist_slag_vp2: qJ has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6RRRPRR8_gravloadJ_floatb_twist_slag_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [12 1]), ...
  'S6RRRPRR8_gravloadJ_floatb_twist_slag_vp2: pkin has to be [12x1] (double)');
assert( isreal(m) && all(size(m) == [7 1]), ...
  'S6RRRPRR8_gravloadJ_floatb_twist_slag_vp2: m has to be [7x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [7,3]), ...
  'S6RRRPRR8_gravloadJ_floatb_twist_slag_vp2: mrSges has to be [7x3] (double)');

%% Symbolic Calculation
% From joint_gravload_floatb_twist_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2018-11-23 17:55:52
% EndTime: 2018-11-23 17:55:53
% DurationCPUTime: 1.44s
% Computational Cost: add. (1640->148), mult. (1713->193), div. (0->0), fcn. (1686->18), ass. (0->68)
t71 = cos(qJ(5));
t55 = pkin(5) * t71 + pkin(4);
t63 = qJ(5) + qJ(6);
t59 = sin(t63);
t60 = cos(t63);
t67 = sin(qJ(5));
t133 = m(6) * pkin(4) + m(7) * t55 + t71 * mrSges(6,1) + t60 * mrSges(7,1) - t67 * mrSges(6,2) - t59 * mrSges(7,2) + mrSges(5,1);
t132 = mrSges(5,2) + m(7) * (-pkin(11) - pkin(10)) - mrSges(7,3) - m(6) * pkin(10) - mrSges(6,3);
t62 = qJ(3) + pkin(12);
t57 = sin(t62);
t58 = cos(t62);
t68 = sin(qJ(3));
t72 = cos(qJ(3));
t131 = m(4) * pkin(2) + t72 * mrSges(4,1) - t68 * mrSges(4,2) - t132 * t57 + t133 * t58 + mrSges(3,1);
t135 = m(5) + m(6) + m(7);
t64 = sin(pkin(6));
t70 = sin(qJ(1));
t119 = t64 * t70;
t124 = cos(qJ(1));
t112 = pkin(6) - qJ(2);
t101 = sin(t112);
t111 = pkin(6) + qJ(2);
t95 = sin(t111) / 0.2e1;
t40 = t95 - t101 / 0.2e1;
t73 = cos(qJ(2));
t86 = t124 * t73 - t70 * t40;
t17 = t72 * t119 - t86 * t68;
t102 = cos(t112);
t96 = cos(t111) / 0.2e1;
t41 = t96 - t102 / 0.2e1;
t65 = cos(pkin(6));
t141 = t41 * t68 + t65 * t72;
t89 = m(4) * pkin(9) - mrSges(3,2) + mrSges(4,3) + mrSges(5,3);
t140 = mrSges(7,1) * t59 + mrSges(6,2) * t71 + mrSges(7,2) * t60 + t89;
t128 = m(7) * pkin(5);
t137 = mrSges(6,1) + t128;
t106 = t64 * t124;
t87 = t124 * t40 + t70 * t73;
t82 = t72 * t106 + t87 * t68;
t113 = t67 * t128;
t130 = -t67 * mrSges(6,1) - t113 - t140;
t12 = -t57 * t106 + t58 * t87;
t69 = sin(qJ(2));
t77 = t102 / 0.2e1 + t96;
t29 = -t124 * t77 + t69 * t70;
t127 = (-t12 * t59 + t29 * t60) * mrSges(7,1) + (-t12 * t60 - t29 * t59) * mrSges(7,2);
t16 = t57 * t119 + t58 * t86;
t32 = t124 * t69 + t70 * t77;
t5 = -t16 * t59 + t32 * t60;
t6 = t16 * t60 + t32 * t59;
t126 = t5 * mrSges(7,1) - t6 * mrSges(7,2);
t24 = -t41 * t58 + t57 * t65;
t39 = t95 + t101 / 0.2e1;
t125 = (-t24 * t59 - t39 * t60) * mrSges(7,1) + (-t24 * t60 + t39 * t59) * mrSges(7,2);
t114 = t124 * pkin(1) + pkin(8) * t119;
t110 = t68 * t119;
t104 = -t70 * pkin(1) + pkin(8) * t106;
t11 = -t58 * t106 - t57 * t87;
t48 = t68 * t106;
t103 = -t72 * t87 + t48;
t56 = pkin(3) * t72 + pkin(2);
t66 = -qJ(4) - pkin(9);
t97 = pkin(3) * t110 - t32 * t66 + t56 * t86 + t114;
t7 = -t16 * t67 + t32 * t71;
t18 = t72 * t86 + t110;
t15 = -t58 * t119 + t57 * t86;
t8 = t16 * t71 + t32 * t67;
t1 = [(-t124 * mrSges(2,1) - m(3) * t114 - t86 * mrSges(3,1) - m(4) * (pkin(2) * t86 + t114) - t18 * mrSges(4,1) - t17 * mrSges(4,2) - m(5) * t97 - t16 * mrSges(5,1) - m(6) * (pkin(4) * t16 + t97) - t8 * mrSges(6,1) - t7 * mrSges(6,2) - m(7) * (t16 * t55 + t97) - t6 * mrSges(7,1) - t5 * mrSges(7,2) + (-mrSges(3,3) * t64 + mrSges(2,2)) * t70 + (-t89 - t113) * t32 + t132 * t15) * g(2) + (t70 * mrSges(2,1) + t124 * mrSges(2,2) - m(3) * t104 + t87 * mrSges(3,1) - mrSges(3,3) * t106 - m(4) * (-pkin(2) * t87 + t104) - t103 * mrSges(4,1) - t82 * mrSges(4,2) + t133 * t12 + (t137 * t67 + t140) * t29 + t132 * t11 + t135 * (-pkin(3) * t48 - t29 * t66 + t56 * t87 - t104)) * g(1) (-t135 * (t39 * t56 + t41 * t66) - t130 * t41 - t131 * t39) * g(3) + (-t135 * (-t29 * t56 - t66 * t87) + t130 * t87 + t131 * t29) * g(2) + (-t135 * (-t32 * t56 - t66 * t86) + t130 * t86 + t131 * t32) * g(1) (-t141 * mrSges(4,1) - (t41 * t72 - t65 * t68) * mrSges(4,2) + t132 * t24 - t133 * (t41 * t57 + t58 * t65)) * g(3) + (t82 * mrSges(4,1) - t103 * mrSges(4,2) - t133 * t11 + t12 * t132) * g(2) + (-mrSges(4,1) * t17 + mrSges(4,2) * t18 + t132 * t16 + t133 * t15) * g(1) + (-g(1) * t17 + g(2) * t82 - g(3) * t141) * t135 * pkin(3), t135 * (-g(1) * t32 - g(2) * t29 + g(3) * t39) (-(-t24 * t71 + t39 * t67) * mrSges(6,2) - t125 - t137 * (-t24 * t67 - t39 * t71)) * g(3) + (-(-t12 * t71 - t29 * t67) * mrSges(6,2) - t127 - t137 * (-t12 * t67 + t29 * t71)) * g(2) + (t8 * mrSges(6,2) - t137 * t7 - t126) * g(1), -g(1) * t126 - g(2) * t127 - g(3) * t125];
taug  = t1(:);
