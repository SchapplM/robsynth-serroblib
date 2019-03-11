% Calculate Gravitation load on the joints for
% S6PRRPRR3
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% g [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [13x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,alpha2,alpha3,d2,d3,d5,d6,theta1,theta4]';
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
% Datum: 2019-03-08 22:09
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function taug = S6PRRPRR3_gravloadJ_floatb_twist_slag_vp2(qJ, g, ...
  pkin, m, mrSges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(3,1),zeros(13,1),zeros(7,1),zeros(7,3)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6PRRPRR3_gravloadJ_floatb_twist_slag_vp2: qJ has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6PRRPRR3_gravloadJ_floatb_twist_slag_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [13 1]), ...
  'S6PRRPRR3_gravloadJ_floatb_twist_slag_vp2: pkin has to be [13x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6PRRPRR3_gravloadJ_floatb_twist_slag_vp2: m has to be [7x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [7,3]), ...
  'S6PRRPRR3_gravloadJ_floatb_twist_slag_vp2: mrSges has to be [7x3] (double)');

%% Symbolic Calculation
% From gravload_joint_floatb_twist_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-08 22:02:19
% EndTime: 2019-03-08 22:02:24
% DurationCPUTime: 1.71s
% Computational Cost: add. (1153->159), mult. (3111->251), div. (0->0), fcn. (3967->16), ass. (0->89)
t146 = m(6) + m(7);
t70 = sin(qJ(6));
t74 = cos(qJ(6));
t134 = -m(7) * pkin(5) - t74 * mrSges(7,1) + t70 * mrSges(7,2) - mrSges(6,1);
t139 = -m(7) * pkin(11) + mrSges(6,2) - mrSges(7,3);
t76 = cos(qJ(3));
t77 = cos(qJ(2));
t113 = t76 * t77;
t72 = sin(qJ(3));
t73 = sin(qJ(2));
t116 = t72 * t73;
t69 = cos(pkin(7));
t145 = t69 * t113 - t116;
t136 = t70 * mrSges(7,1) + t74 * mrSges(7,2) - mrSges(5,2) + mrSges(6,3);
t143 = t72 * mrSges(4,2);
t141 = m(5) + t146;
t71 = sin(qJ(5));
t75 = cos(qJ(5));
t140 = t134 * t75 + t139 * t71 - mrSges(5,1);
t110 = cos(pkin(6));
t109 = cos(pkin(13));
t65 = sin(pkin(13));
t56 = -t109 * t72 - t76 * t65;
t66 = sin(pkin(7));
t46 = t56 * t66;
t48 = t56 * t69;
t67 = sin(pkin(6));
t87 = t109 * t76 - t72 * t65;
t22 = t67 * (-t48 * t77 + t73 * t87) - t110 * t46;
t137 = t146 * pkin(10) + t136;
t135 = -m(4) * pkin(2) - t76 * mrSges(4,1) - mrSges(3,1) + t143;
t117 = t69 * t76;
t118 = t69 * t72;
t133 = t118 * mrSges(4,1) + t117 * mrSges(4,2) + mrSges(3,2) + (-m(4) * pkin(9) - mrSges(4,3) - mrSges(5,3)) * t66;
t131 = pkin(3) * t76;
t68 = cos(pkin(12));
t101 = t68 * t110;
t108 = sin(pkin(12));
t52 = t101 * t73 + t108 * t77;
t128 = t52 * t72;
t91 = t110 * t108;
t54 = t68 * t77 - t73 * t91;
t127 = t54 * t72;
t124 = t66 * t71;
t123 = t66 * t73;
t122 = t66 * t75;
t121 = t67 * t68;
t120 = t67 * t73;
t119 = t67 * t77;
t115 = t72 * t77;
t114 = t73 * t76;
t49 = pkin(3) * t118 + (-pkin(9) - qJ(4)) * t66;
t51 = t101 * t77 - t108 * t73;
t64 = pkin(2) + t131;
t112 = -t52 * t49 + t51 * t64;
t53 = -t68 * t73 - t77 * t91;
t111 = -t54 * t49 + t53 * t64;
t107 = pkin(3) * t117;
t106 = t66 * t121;
t105 = t66 * t120;
t102 = t67 * t108;
t100 = t76 * t110;
t98 = t64 * t119 - t120 * t49;
t97 = t66 * t102;
t94 = -pkin(3) * t127 + t53 * t107 + t97 * t131;
t90 = -t51 * t69 + t106;
t88 = (t100 * t66 + t145 * t67) * pkin(3);
t82 = t53 * t69 + t97;
t11 = -t121 * t46 + t48 * t51 - t52 * t87;
t80 = t51 * t107 + (-t106 * t76 - t128) * pkin(3);
t14 = t102 * t46 + t53 * t48 - t54 * t87;
t50 = t110 * t69 - t119 * t66;
t47 = t87 * t69;
t45 = t87 * t66;
t37 = t102 * t69 - t53 * t66;
t36 = -t121 * t69 - t51 * t66;
t33 = (t48 * t73 + t77 * t87) * t67;
t32 = t119 * t56 - t120 * t47;
t28 = t48 * t54 + t53 * t87;
t27 = -t47 * t54 + t53 * t56;
t26 = t48 * t52 + t51 * t87;
t25 = -t47 * t52 + t51 * t56;
t21 = t110 * t45 + (t47 * t77 + t56 * t73) * t67;
t18 = t22 * t75 + t50 * t71;
t15 = t102 * t45 + t53 * t47 + t54 * t56;
t12 = -t121 * t45 + t47 * t51 + t52 * t56;
t4 = -t14 * t75 + t37 * t71;
t2 = -t11 * t75 + t36 * t71;
t1 = [(-m(2) - m(3) - m(4) - t141) * g(3) (-m(5) * t112 - t26 * mrSges(5,1) - t146 * (t26 * pkin(4) - pkin(10) * t25 + t112) + t139 * (-t122 * t52 + t26 * t71) + t134 * (t124 * t52 + t26 * t75) + t136 * t25 + t135 * t51 + t133 * t52) * g(2) + (-m(5) * t111 - t28 * mrSges(5,1) - t146 * (t28 * pkin(4) - pkin(10) * t27 + t111) + t139 * (-t122 * t54 + t28 * t71) + t134 * (t124 * t54 + t28 * t75) + t136 * t27 + t135 * t53 + t133 * t54) * g(1) + ((-mrSges(3,1) * t77 + mrSges(3,2) * t73 - m(4) * (pkin(2) * t77 + pkin(9) * t123) - (-t116 * t69 + t113) * mrSges(4,1) - (-t114 * t69 - t115) * mrSges(4,2) - mrSges(4,3) * t123) * t67 - m(5) * t98 - t33 * mrSges(5,1) - mrSges(5,3) * t105 - t146 * (t33 * pkin(4) - pkin(10) * t32 + t98) + t134 * (t105 * t71 + t33 * t75) + t136 * t32 + t139 * (-t105 * t75 + t33 * t71)) * g(3) (-(mrSges(4,1) * t100 - t110 * t143) * t66 - (t145 * mrSges(4,1) + (-t115 * t69 - t114) * mrSges(4,2)) * t67 - m(5) * t88 - t146 * (t21 * pkin(4) + t88) + t140 * t21 - t137 * t22) * g(3) + (-(-t76 * t90 - t128) * mrSges(4,1) - (-t52 * t76 + t72 * t90) * mrSges(4,2) - m(5) * t80 - t146 * (t12 * pkin(4) + t80) + t140 * t12 + t137 * t11) * g(2) + (-(t76 * t82 - t127) * mrSges(4,1) - (-t54 * t76 - t72 * t82) * mrSges(4,2) - m(5) * t94 - t146 * (t15 * pkin(4) + t94) + t140 * t15 + t137 * t14) * g(1), t141 * (-g(1) * t37 - g(2) * t36 - g(3) * t50) (t139 * t18 + t134 * (-t22 * t71 + t50 * t75)) * g(3) + (t139 * t2 + t134 * (t11 * t71 + t36 * t75)) * g(2) + (t139 * t4 + t134 * (t14 * t71 + t37 * t75)) * g(1), -g(1) * ((-t15 * t74 - t4 * t70) * mrSges(7,1) + (t15 * t70 - t4 * t74) * mrSges(7,2)) - g(2) * ((-t12 * t74 - t2 * t70) * mrSges(7,1) + (t12 * t70 - t2 * t74) * mrSges(7,2)) - g(3) * ((-t18 * t70 - t21 * t74) * mrSges(7,1) + (-t18 * t74 + t21 * t70) * mrSges(7,2))];
taug  = t1(:);
