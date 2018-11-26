% Calculate Gravitation load on the joints for
% S6RRRPRP6
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% g [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [11x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,alpha2,d1,d2,d3,d5,theta4]';
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
% Datum: 2018-11-23 17:45
% Revision: 76f9d5e39f14dc242b53c0d9d3d9db48bd8f37c0
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für mechatronische Systeme, Universität Hannover

function taug = S6RRRPRP6_gravloadJ_floatb_twist_slag_vp2(qJ, g, ...
  pkin, m, mrSges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(3,1),zeros(11,1),zeros(7,1),zeros(7,3)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRRPRP6_gravloadJ_floatb_twist_slag_vp2: qJ has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6RRRPRP6_gravloadJ_floatb_twist_slag_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [11 1]), ...
  'S6RRRPRP6_gravloadJ_floatb_twist_slag_vp2: pkin has to be [11x1] (double)');
assert( isreal(m) && all(size(m) == [7 1]), ...
  'S6RRRPRP6_gravloadJ_floatb_twist_slag_vp2: m has to be [7x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [7,3]), ...
  'S6RRRPRP6_gravloadJ_floatb_twist_slag_vp2: mrSges has to be [7x3] (double)');

%% Symbolic Calculation
% From joint_gravload_floatb_twist_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2018-11-23 17:44:26
% EndTime: 2018-11-23 17:44:28
% DurationCPUTime: 1.37s
% Computational Cost: add. (1539->148), mult. (1660->196), div. (0->0), fcn. (1627->16), ass. (0->70)
t72 = cos(qJ(5));
t59 = pkin(5) * t72 + pkin(4);
t144 = -m(6) * pkin(4) - m(7) * t59 - mrSges(5,1);
t80 = -m(6) * pkin(10) + m(7) * (-qJ(6) - pkin(10)) + mrSges(5,2) - mrSges(6,3) - mrSges(7,3);
t143 = mrSges(6,1) + mrSges(7,1);
t137 = mrSges(6,2) + mrSges(7,2);
t136 = m(5) + m(6) + m(7);
t109 = cos(pkin(6));
t107 = pkin(6) + qJ(2);
t92 = cos(t107) / 0.2e1;
t108 = pkin(6) - qJ(2);
t99 = cos(t108);
t45 = t92 - t99 / 0.2e1;
t69 = sin(qJ(3));
t73 = cos(qJ(3));
t142 = t109 * t73 + t45 * t69;
t65 = sin(pkin(6));
t71 = sin(qJ(1));
t116 = t65 * t71;
t129 = cos(qJ(1));
t91 = sin(t107) / 0.2e1;
t98 = sin(t108);
t44 = t91 - t98 / 0.2e1;
t74 = cos(qJ(2));
t86 = t129 * t74 - t71 * t44;
t21 = t73 * t116 - t86 * t69;
t64 = qJ(3) + pkin(11);
t61 = sin(t64);
t62 = cos(t64);
t141 = -m(4) * pkin(2) - t73 * mrSges(4,1) + t69 * mrSges(4,2) + t144 * t62 + t80 * t61 - mrSges(3,1);
t103 = t65 * t129;
t87 = t129 * t44 + t71 * t74;
t14 = -t61 * t103 + t62 * t87;
t70 = sin(qJ(2));
t75 = t99 / 0.2e1 + t92;
t33 = -t129 * t75 + t70 * t71;
t68 = sin(qJ(5));
t140 = t14 * t68 - t33 * t72;
t139 = -t14 * t72 - t33 * t68;
t130 = m(7) * pkin(5);
t82 = t73 * t103 + t87 * t69;
t134 = -t137 * t68 + t143 * t72 - t144;
t133 = -t130 - t143;
t132 = -m(4) * pkin(9) + mrSges(3,2) - mrSges(4,3) - mrSges(5,3);
t125 = t87 * t68;
t123 = t86 * t68;
t122 = t45 * t68;
t118 = t62 * t68;
t117 = t62 * t72;
t110 = t129 * pkin(1) + pkin(8) * t116;
t106 = t69 * t116;
t102 = -pkin(1) * t71 + pkin(8) * t103;
t51 = t69 * t103;
t101 = -t73 * t87 + t51;
t36 = t129 * t70 + t71 * t75;
t60 = pkin(3) * t73 + pkin(2);
t67 = -qJ(4) - pkin(9);
t94 = pkin(3) * t106 - t36 * t67 + t60 * t86 + t110;
t18 = t61 * t116 + t62 * t86;
t5 = -t18 * t68 + t36 * t72;
t88 = pkin(3) * t51 + t33 * t67 - t60 * t87 + t102;
t13 = t62 * t103 + t61 * t87;
t79 = t68 * t130 - t132;
t43 = t91 + t98 / 0.2e1;
t28 = t109 * t61 - t45 * t62;
t27 = -t109 * t62 - t45 * t61;
t22 = t73 * t86 + t106;
t17 = -t62 * t116 + t61 * t86;
t6 = t18 * t72 + t36 * t68;
t1 = [(-t129 * mrSges(2,1) - m(3) * t110 - t86 * mrSges(3,1) - m(4) * (pkin(2) * t86 + t110) - t22 * mrSges(4,1) - t21 * mrSges(4,2) - m(5) * t94 - t18 * mrSges(5,1) - m(6) * (pkin(4) * t18 + t94) - m(7) * (t18 * t59 + t94) + (-mrSges(3,3) * t65 + mrSges(2,2)) * t71 - t143 * t6 - t137 * t5 - t79 * t36 + t80 * t17) * g(2) + (t71 * mrSges(2,1) + t129 * mrSges(2,2) - m(3) * t102 + t87 * mrSges(3,1) - mrSges(3,3) * t103 - m(4) * (-pkin(2) * t87 + t102) - t101 * mrSges(4,1) - t82 * mrSges(4,2) - m(5) * t88 + t14 * mrSges(5,1) - m(6) * (-pkin(4) * t14 + t88) - m(7) * (-t14 * t59 + t88) - t143 * t139 - t137 * t140 + t79 * t33 - t80 * t13) * g(1) (t122 * t130 - t143 * (t43 * t117 - t122) - t137 * (-t43 * t118 - t45 * t72) - t136 * (t43 * t60 + t45 * t67) - t132 * t45 + t141 * t43) * g(3) + (-t125 * t130 - t143 * (-t33 * t117 + t125) - t137 * (t33 * t118 + t72 * t87) - t136 * (-t33 * t60 - t67 * t87) + t132 * t87 - t141 * t33) * g(2) + (-t123 * t130 - t137 * (t36 * t118 + t72 * t86) - t136 * (-t36 * t60 - t67 * t86) - t143 * (-t36 * t117 + t123) + t132 * t86 - t141 * t36) * g(1) (-t142 * mrSges(4,1) - (-t109 * t69 + t45 * t73) * mrSges(4,2) + t80 * t28 + t134 * t27) * g(3) + (t82 * mrSges(4,1) - t101 * mrSges(4,2) + t134 * t13 + t80 * t14) * g(2) + (-mrSges(4,1) * t21 + mrSges(4,2) * t22 + t134 * t17 + t80 * t18) * g(1) + (-g(1) * t21 + g(2) * t82 - g(3) * t142) * t136 * pkin(3), t136 * (-g(1) * t36 - g(2) * t33 + g(3) * t43) (-t137 * (-t28 * t72 + t43 * t68) + t133 * (-t28 * t68 - t43 * t72)) * g(3) + (-t133 * t140 - t137 * t139) * g(2) + (t133 * t5 + t137 * t6) * g(1) (-g(1) * t17 - g(2) * t13 - g(3) * t27) * m(7)];
taug  = t1(:);
