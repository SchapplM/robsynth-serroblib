% Calculate Gravitation load on the joints for
% S6RRPRRP5
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% g [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [11x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,alpha2,d1,d2,d4,d5,theta3]';
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

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-09 12:06
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function taug = S6RRPRRP5_gravloadJ_floatb_twist_slag_vp2(qJ, g, ...
  pkin, m, mrSges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(3,1),zeros(11,1),zeros(7,1),zeros(7,3)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRPRRP5_gravloadJ_floatb_twist_slag_vp2: qJ has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6RRPRRP5_gravloadJ_floatb_twist_slag_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [11 1]), ...
  'S6RRPRRP5_gravloadJ_floatb_twist_slag_vp2: pkin has to be [11x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6RRPRRP5_gravloadJ_floatb_twist_slag_vp2: m has to be [7x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [7,3]), ...
  'S6RRPRRP5_gravloadJ_floatb_twist_slag_vp2: mrSges has to be [7x3] (double)');

%% Symbolic Calculation
% From gravload_joint_floatb_twist_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 11:57:50
% EndTime: 2019-03-09 11:57:52
% DurationCPUTime: 1.36s
% Computational Cost: add. (888->136), mult. (2176->189), div. (0->0), fcn. (2705->12), ass. (0->73)
t65 = cos(qJ(5));
t55 = pkin(5) * t65 + pkin(4);
t137 = m(6) * pkin(4) + m(7) * t55 + mrSges(5,1);
t136 = m(6) * pkin(10) - m(7) * (-qJ(6) - pkin(10)) - mrSges(5,2) + mrSges(6,3) + mrSges(7,3);
t135 = mrSges(6,1) + mrSges(7,1);
t126 = mrSges(6,2) + mrSges(7,2);
t57 = sin(pkin(11));
t63 = sin(qJ(2));
t67 = cos(qJ(2));
t96 = cos(pkin(11));
t44 = -t57 * t67 - t63 * t96;
t64 = sin(qJ(1));
t68 = cos(qJ(1));
t59 = cos(pkin(6));
t79 = -t57 * t63 + t67 * t96;
t73 = t59 * t79;
t29 = t44 * t68 - t64 * t73;
t97 = t44 * t59;
t30 = t64 * t97 + t68 * t79;
t105 = t64 * t67;
t108 = t63 * t68;
t41 = -t105 * t59 - t108;
t77 = t41 * pkin(2);
t134 = pkin(3) * t29 + pkin(9) * t30 + t77;
t102 = t67 * t68;
t106 = t64 * t63;
t133 = t102 * t59 - t106;
t62 = sin(qJ(4));
t66 = cos(qJ(4));
t132 = t136 * t62 + t137 * t66 + mrSges(4,1);
t58 = sin(pkin(6));
t113 = t58 * t68;
t25 = -t64 * t79 + t68 * t97;
t14 = -t113 * t62 - t25 * t66;
t26 = t44 * t64 + t68 * t73;
t61 = sin(qJ(5));
t131 = t14 * t61 + t26 * t65;
t130 = -t14 * t65 + t26 * t61;
t121 = m(7) * pkin(5);
t129 = -m(5) - m(6);
t127 = mrSges(4,2) - mrSges(5,3);
t95 = m(7) - t129;
t125 = -t126 * t61 + t135 * t65 + t137;
t124 = -t121 - t135;
t123 = pkin(9) * t95 + t121 * t61 - t127;
t120 = pkin(2) * t67;
t119 = t25 * t61;
t116 = t30 * t61;
t37 = t44 * t58;
t115 = t37 * t61;
t114 = t58 * t64;
t112 = t61 * t66;
t104 = t65 * t66;
t56 = pkin(1) + t120;
t49 = t68 * t56;
t98 = pkin(3) * t30 + t49;
t52 = t58 * t120;
t92 = m(4) + t95;
t91 = -m(3) * pkin(1) - mrSges(2,1);
t85 = t133 * pkin(2);
t18 = t114 * t62 + t30 * t66;
t5 = -t18 * t61 - t29 * t65;
t13 = t113 * t66 - t25 * t62;
t69 = mrSges(2,2) + (-m(3) * pkin(8) - mrSges(3,3) - mrSges(4,3)) * t58 + t92 * (pkin(2) * t59 * t63 + (-pkin(8) - qJ(3)) * t58);
t42 = -t106 * t59 + t102;
t40 = -t108 * t59 - t105;
t36 = t79 * t58;
t32 = -t37 * t66 + t59 * t62;
t31 = -t37 * t62 - t59 * t66;
t22 = t25 * pkin(3);
t17 = -t114 * t66 + t30 * t62;
t6 = t18 * t65 - t29 * t61;
t1 = [(-t42 * mrSges(3,1) - t41 * mrSges(3,2) - m(4) * t49 - t30 * mrSges(4,1) - m(5) * t98 - t18 * mrSges(5,1) - m(6) * (pkin(4) * t18 + t98) - m(7) * (t18 * t55 + t98) + t91 * t68 - t135 * t6 - t126 * t5 + t69 * t64 + t123 * t29 - t136 * t17) * g(2) + (-t40 * mrSges(3,1) + t133 * mrSges(3,2) - t25 * mrSges(4,1) - m(5) * t22 + t14 * mrSges(5,1) - m(6) * (-pkin(4) * t14 + t22) - m(7) * (-t14 * t55 + t22) - t135 * t130 - t126 * t131 + (t56 * t92 - t91) * t64 + t69 * t68 - t123 * t26 + t136 * t13) * g(1) (-(mrSges(3,1) * t67 - mrSges(3,2) * t63) * t58 - m(4) * t52 + t115 * t121 - t95 * (pkin(3) * t36 - pkin(9) * t37 + t52) - t127 * t37 - t135 * (t104 * t36 - t115) - t126 * (-t112 * t36 - t37 * t65) - t132 * t36) * g(3) + (t119 * t121 - m(4) * t85 - mrSges(3,1) * t133 - mrSges(3,2) * t40 - t135 * (t104 * t26 - t119) - t95 * (pkin(3) * t26 - pkin(9) * t25 + t85) - t126 * (-t112 * t26 - t25 * t65) - t127 * t25 - t132 * t26) * g(2) + (-mrSges(3,1) * t41 + mrSges(3,2) * t42 - m(4) * t77 - m(7) * (pkin(5) * t116 + t134) - t126 * (-t112 * t29 + t30 * t65) + t129 * t134 + t127 * t30 - t135 * (t104 * t29 + t116) - t132 * t29) * g(1) (-g(3) * t59 + (-g(1) * t64 + g(2) * t68) * t58) * t92 (t125 * t31 - t136 * t32) * g(3) + (t125 * t13 - t136 * t14) * g(2) + (t125 * t17 - t136 * t18) * g(1) (-t126 * (-t32 * t65 + t36 * t61) + t124 * (-t32 * t61 - t36 * t65)) * g(3) + (-t124 * t131 - t126 * t130) * g(2) + (t124 * t5 + t126 * t6) * g(1) (-g(1) * t17 - g(2) * t13 - g(3) * t31) * m(7)];
taug  = t1(:);
