% Calculate Gravitation load on the joints for
% S6RRPRRR5
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
% Datum: 2018-11-23 17:24
% Revision: 76f9d5e39f14dc242b53c0d9d3d9db48bd8f37c0
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für mechatronische Systeme, Universität Hannover

function taug = S6RRPRRR5_gravloadJ_floatb_twist_slag_vp2(qJ, g, ...
  pkin, m, mrSges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(3,1),zeros(12,1),zeros(7,1),zeros(7,3)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRPRRR5_gravloadJ_floatb_twist_slag_vp2: qJ has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6RRPRRR5_gravloadJ_floatb_twist_slag_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [12 1]), ...
  'S6RRPRRR5_gravloadJ_floatb_twist_slag_vp2: pkin has to be [12x1] (double)');
assert( isreal(m) && all(size(m) == [7 1]), ...
  'S6RRPRRR5_gravloadJ_floatb_twist_slag_vp2: m has to be [7x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [7,3]), ...
  'S6RRPRRR5_gravloadJ_floatb_twist_slag_vp2: mrSges has to be [7x3] (double)');

%% Symbolic Calculation
% From joint_gravload_floatb_twist_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2018-11-23 17:23:28
% EndTime: 2018-11-23 17:23:29
% DurationCPUTime: 1.35s
% Computational Cost: add. (1893->149), mult. (1573->191), div. (0->0), fcn. (1510->22), ass. (0->85)
t73 = cos(qJ(5));
t53 = pkin(5) * t73 + pkin(4);
t66 = qJ(5) + qJ(6);
t61 = sin(t66);
t62 = cos(t66);
t69 = sin(qJ(5));
t132 = m(6) * pkin(4) + m(7) * t53 + t73 * mrSges(6,1) + t62 * mrSges(7,1) - t69 * mrSges(6,2) - t61 * mrSges(7,2) + mrSges(5,1);
t83 = mrSges(5,2) + m(7) * (-pkin(11) - pkin(10)) - mrSges(7,3) - m(6) * pkin(10) - mrSges(6,3);
t110 = m(5) + m(6) + m(7);
t131 = t110 * pkin(9) - mrSges(4,2) + mrSges(5,3);
t144 = t61 * mrSges(7,1) + t73 * mrSges(6,2) + t62 * mrSges(7,2) + t131;
t70 = sin(qJ(4));
t74 = cos(qJ(4));
t130 = t132 * t74 - t83 * t70 + mrSges(4,1);
t143 = m(7) * pkin(5);
t141 = t69 * t143;
t134 = mrSges(6,1) + t143;
t64 = pkin(6) + qJ(2);
t52 = cos(t64) / 0.2e1;
t65 = pkin(6) - qJ(2);
t60 = cos(t65);
t40 = t60 / 0.2e1 + t52;
t129 = -t69 * mrSges(6,1) - t141 - t144;
t127 = sin(t64) / 0.2e1;
t67 = sin(pkin(6));
t76 = cos(qJ(1));
t119 = t67 * t76;
t63 = qJ(2) + pkin(12);
t58 = cos(t63);
t72 = sin(qJ(1));
t116 = t72 * t58;
t104 = pkin(6) + t63;
t95 = sin(t104);
t90 = t95 / 0.2e1;
t105 = pkin(6) - t63;
t96 = sin(t105);
t85 = t90 - t96 / 0.2e1;
t22 = t76 * t85 + t116;
t12 = -t70 * t119 + t22 * t74;
t55 = sin(t63);
t92 = cos(t105) / 0.2e1;
t97 = cos(t104);
t81 = t97 / 0.2e1 + t92;
t21 = t55 * t72 - t76 * t81;
t125 = (-t12 * t61 + t21 * t62) * mrSges(7,1) + (-t12 * t62 - t21 * t61) * mrSges(7,2);
t120 = t67 * t72;
t114 = t76 * t58;
t25 = -t72 * t85 + t114;
t16 = t120 * t70 + t25 * t74;
t24 = t76 * t55 + t72 * t81;
t5 = -t16 * t61 + t24 * t62;
t6 = t16 * t62 + t24 * t61;
t124 = t5 * mrSges(7,1) - t6 * mrSges(7,2);
t57 = sin(t65);
t123 = pkin(2) * t57;
t37 = t92 - t97 / 0.2e1;
t68 = cos(pkin(6));
t28 = t37 * t74 + t68 * t70;
t91 = t96 / 0.2e1;
t36 = t91 + t90;
t121 = (-t28 * t61 - t36 * t62) * mrSges(7,1) + (-t28 * t62 + t36 * t61) * mrSges(7,2);
t71 = sin(qJ(2));
t118 = t71 * t72;
t117 = t71 * t76;
t75 = cos(qJ(2));
t54 = pkin(2) * t75 + pkin(1);
t41 = t76 * t54;
t112 = t25 * pkin(3) + t41;
t48 = pkin(2) * t127;
t111 = t123 / 0.2e1 + t48;
t108 = m(4) + t110;
t11 = -t74 * t119 - t22 * t70;
t38 = t40 * pkin(2);
t103 = -pkin(2) * t118 + t76 * t38;
t7 = -t16 * t69 + t24 * t73;
t94 = m(3) * pkin(1) + t75 * mrSges(3,1) + mrSges(2,1);
t93 = -pkin(2) * t117 - t72 * t38;
t86 = -t95 / 0.2e1 + t91;
t39 = t127 - t57 / 0.2e1;
t78 = t39 * mrSges(3,1) + mrSges(2,2) + (-m(3) * pkin(8) - mrSges(3,3) - mrSges(4,3)) * t67 + t108 * (-t123 / 0.2e1 + t48 - t67 * (pkin(8) + qJ(3)));
t30 = -t40 * t72 - t117;
t29 = -t40 * t76 + t118;
t15 = -t120 * t74 + t25 * t70;
t8 = t16 * t73 + t24 * t69;
t1 = [(-t30 * mrSges(3,2) - m(4) * t41 - t25 * mrSges(4,1) - m(5) * t112 - t16 * mrSges(5,1) - m(6) * (pkin(4) * t16 + t112) - t8 * mrSges(6,1) - t7 * mrSges(6,2) - m(7) * (t16 * t53 + t112) - t6 * mrSges(7,1) - t5 * mrSges(7,2) - t94 * t76 + t78 * t72 + (-t131 - t141) * t24 + t83 * t15) * g(2) + (-t29 * mrSges(3,2) + (t108 * t54 + t94) * t72 + t78 * t76 + t132 * t12 + (t134 * t69 + t144) * t21 + t83 * t11 + (pkin(3) * t110 + mrSges(4,1)) * t22) * g(1) (-(t127 + t57 / 0.2e1) * mrSges(3,1) - (t52 - t60 / 0.2e1) * mrSges(3,2) - m(4) * t111 - t110 * (t36 * pkin(3) + t111) + t129 * t37 - t130 * t36) * g(3) + (t29 * mrSges(3,1) - (-t39 * t76 - t72 * t75) * mrSges(3,2) - m(4) * t103 - t110 * (-t21 * pkin(3) + t103) + t129 * (-t76 * t86 + t116) + t130 * t21) * g(2) + (-t30 * mrSges(3,1) - (t39 * t72 - t75 * t76) * mrSges(3,2) - m(4) * t93 - t110 * (-t24 * pkin(3) + t93) + t129 * (t72 * t86 + t114) + t130 * t24) * g(1) (-g(3) * t68 + (-g(1) * t72 + g(2) * t76) * t67) * t108 (t83 * t28 - t132 * (-t37 * t70 + t68 * t74)) * g(3) + (-t132 * t11 + t83 * t12) * g(2) + (t132 * t15 + t83 * t16) * g(1) (-(-t28 * t73 + t36 * t69) * mrSges(6,2) - t121 - t134 * (-t28 * t69 - t36 * t73)) * g(3) + (-(-t12 * t73 - t21 * t69) * mrSges(6,2) - t125 - t134 * (-t12 * t69 + t21 * t73)) * g(2) + (t8 * mrSges(6,2) - t134 * t7 - t124) * g(1), -g(1) * t124 - g(2) * t125 - g(3) * t121];
taug  = t1(:);
