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

% Quelle: HybrDyn-Toolbox (ehem. IRT-Maple-Toolbox)
% Datum: 2018-11-23 17:14
% Revision: 76f9d5e39f14dc242b53c0d9d3d9db48bd8f37c0
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für mechatronische Systeme, Universität Hannover

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
assert( isreal(m) && all(size(m) == [7 1]), ...
  'S6RRPRRP5_gravloadJ_floatb_twist_slag_vp2: m has to be [7x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [7,3]), ...
  'S6RRPRRP5_gravloadJ_floatb_twist_slag_vp2: mrSges has to be [7x3] (double)');

%% Symbolic Calculation
% From joint_gravload_floatb_twist_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2018-11-23 17:13:18
% EndTime: 2018-11-23 17:13:19
% DurationCPUTime: 1.34s
% Computational Cost: add. (1774->152), mult. (1520->198), div. (0->0), fcn. (1451->20), ass. (0->89)
t74 = cos(qJ(5));
t57 = pkin(5) * t74 + pkin(4);
t144 = -m(6) * pkin(4) - m(7) * t57 - mrSges(5,1);
t84 = -m(6) * pkin(10) + m(7) * (-qJ(6) - pkin(10)) + mrSges(5,2) - mrSges(6,3) - mrSges(7,3);
t143 = mrSges(6,1) + mrSges(7,1);
t138 = mrSges(6,2) + mrSges(7,2);
t71 = sin(qJ(4));
t75 = cos(qJ(4));
t142 = t144 * t75 + t84 * t71 - mrSges(4,1);
t68 = sin(pkin(6));
t77 = cos(qJ(1));
t122 = t68 * t77;
t65 = qJ(2) + pkin(11);
t62 = cos(t65);
t73 = sin(qJ(1));
t117 = t73 * t62;
t103 = pkin(6) + t65;
t94 = sin(t103);
t89 = t94 / 0.2e1;
t104 = pkin(6) - t65;
t95 = sin(t104);
t86 = t89 - t95 / 0.2e1;
t26 = t77 * t86 + t117;
t14 = -t71 * t122 + t26 * t75;
t59 = sin(t65);
t91 = cos(t104) / 0.2e1;
t96 = cos(t103);
t79 = t96 / 0.2e1 + t91;
t25 = t59 * t73 - t77 * t79;
t70 = sin(qJ(5));
t141 = t14 * t70 - t25 * t74;
t140 = -t14 * t74 - t25 * t70;
t132 = m(7) * pkin(5);
t108 = m(5) + m(6) + m(7);
t107 = m(4) + t108;
t139 = mrSges(4,2) - mrSges(5,3);
t66 = pkin(6) + qJ(2);
t55 = cos(t66) / 0.2e1;
t67 = pkin(6) - qJ(2);
t64 = cos(t67);
t44 = t64 / 0.2e1 + t55;
t136 = -t138 * t70 + t143 * t74 - t144;
t135 = -t132 - t143;
t134 = pkin(9) * t108 + t132 * t70 - t139;
t131 = sin(t66) / 0.2e1;
t61 = sin(t67);
t129 = pkin(2) * t61;
t90 = t95 / 0.2e1;
t87 = -t94 / 0.2e1 + t90;
t27 = -t77 * t87 + t117;
t126 = t27 * t70;
t114 = t77 * t62;
t30 = t73 * t87 + t114;
t125 = t30 * t70;
t41 = t91 - t96 / 0.2e1;
t124 = t41 * t70;
t123 = t68 * t73;
t121 = t70 * t75;
t72 = sin(qJ(2));
t118 = t72 * t77;
t116 = t73 * t72;
t115 = t74 * t75;
t29 = -t73 * t86 + t114;
t76 = cos(qJ(2));
t58 = pkin(2) * t76 + pkin(1);
t45 = t77 * t58;
t111 = t29 * pkin(3) + t45;
t51 = pkin(2) * t131;
t110 = t129 / 0.2e1 + t51;
t109 = cos(pkin(6));
t42 = t44 * pkin(2);
t102 = -pkin(2) * t116 + t77 * t42;
t18 = t123 * t71 + t29 * t75;
t28 = t77 * t59 + t73 * t79;
t5 = -t18 * t70 + t28 * t74;
t93 = m(3) * pkin(1) + mrSges(3,1) * t76 + mrSges(2,1);
t92 = -pkin(2) * t118 - t73 * t42;
t13 = t122 * t75 + t26 * t71;
t43 = t131 - t61 / 0.2e1;
t78 = mrSges(3,1) * t43 + mrSges(2,2) + (-m(3) * pkin(8) - mrSges(3,3) - mrSges(4,3)) * t68 + t107 * (-t129 / 0.2e1 + t51 - t68 * (pkin(8) + qJ(3)));
t40 = t90 + t89;
t34 = -t73 * t44 - t118;
t33 = -t44 * t77 + t116;
t32 = t109 * t71 + t41 * t75;
t31 = -t109 * t75 + t41 * t71;
t22 = t26 * pkin(3);
t17 = -t123 * t75 + t29 * t71;
t6 = t18 * t74 + t28 * t70;
t1 = [(-t34 * mrSges(3,2) - m(4) * t45 - t29 * mrSges(4,1) - m(5) * t111 - t18 * mrSges(5,1) - m(6) * (pkin(4) * t18 + t111) - m(7) * (t18 * t57 + t111) - t143 * t6 - t138 * t5 - t93 * t77 + t78 * t73 - t134 * t28 + t84 * t17) * g(2) + (-t33 * mrSges(3,2) + t26 * mrSges(4,1) + m(5) * t22 + t14 * mrSges(5,1) - m(6) * (-pkin(4) * t14 - t22) - m(7) * (-t14 * t57 - t22) - t143 * t140 - t138 * t141 + (t107 * t58 + t93) * t73 + t78 * t77 + t134 * t25 - t84 * t13) * g(1) (-(t131 + t61 / 0.2e1) * mrSges(3,1) - (t55 - t64 / 0.2e1) * mrSges(3,2) - m(4) * t110 - t124 * t132 - t108 * (t40 * pkin(3) + pkin(9) * t41 + t110) + t139 * t41 - t143 * (t115 * t40 + t124) - t138 * (-t121 * t40 + t41 * t74) + t142 * t40) * g(3) + (t33 * mrSges(3,1) - (-t43 * t77 - t73 * t76) * mrSges(3,2) - m(4) * t102 - t126 * t132 - t108 * (-t25 * pkin(3) + pkin(9) * t27 + t102) - t143 * (-t115 * t25 + t126) - t138 * (t121 * t25 + t27 * t74) + t139 * t27 - t142 * t25) * g(2) + (-t34 * mrSges(3,1) - (t73 * t43 - t76 * t77) * mrSges(3,2) - m(4) * t92 - t125 * t132 - t138 * (t121 * t28 + t30 * t74) - t108 * (-t28 * pkin(3) + t30 * pkin(9) + t92) + t139 * t30 - t143 * (-t115 * t28 + t125) - t142 * t28) * g(1) ((-g(1) * t73 + g(2) * t77) * t68 - g(3) * t109) * t107 (t136 * t31 + t84 * t32) * g(3) + (t136 * t13 + t84 * t14) * g(2) + (t136 * t17 + t84 * t18) * g(1) (-t138 * (-t32 * t74 + t40 * t70) + t135 * (-t32 * t70 - t40 * t74)) * g(3) + (-t135 * t141 - t138 * t140) * g(2) + (t135 * t5 + t138 * t6) * g(1) (-g(1) * t17 - g(2) * t13 - g(3) * t31) * m(7)];
taug  = t1(:);
