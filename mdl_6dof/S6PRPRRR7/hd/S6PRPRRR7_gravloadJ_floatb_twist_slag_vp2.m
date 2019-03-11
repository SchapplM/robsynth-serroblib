% Calculate Gravitation load on the joints for
% S6PRPRRR7
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% g [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [14x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,alpha2,alpha3,alpha4,d2,d4,d5,d6,theta1,theta3]';
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
% Datum: 2019-03-08 20:57
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function taug = S6PRPRRR7_gravloadJ_floatb_twist_slag_vp2(qJ, g, ...
  pkin, m, mrSges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(3,1),zeros(14,1),zeros(7,1),zeros(7,3)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6PRPRRR7_gravloadJ_floatb_twist_slag_vp2: qJ has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6PRPRRR7_gravloadJ_floatb_twist_slag_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [14 1]), ...
  'S6PRPRRR7_gravloadJ_floatb_twist_slag_vp2: pkin has to be [14x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6PRPRRR7_gravloadJ_floatb_twist_slag_vp2: m has to be [7x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [7,3]), ...
  'S6PRPRRR7_gravloadJ_floatb_twist_slag_vp2: mrSges has to be [7x3] (double)');

%% Symbolic Calculation
% From gravload_joint_floatb_twist_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-08 20:50:39
% EndTime: 2019-03-08 20:50:43
% DurationCPUTime: 1.50s
% Computational Cost: add. (1809->155), mult. (5207->252), div. (0->0), fcn. (6752->18), ass. (0->101)
t157 = m(6) + m(7);
t71 = sin(qJ(6));
t74 = cos(qJ(6));
t159 = m(7) * pkin(5) + t74 * mrSges(7,1) - t71 * mrSges(7,2) + mrSges(6,1);
t148 = -m(7) * pkin(12) + mrSges(6,2) - mrSges(7,3);
t132 = sin(pkin(6));
t141 = sin(qJ(2));
t119 = t132 * t141;
t70 = sin(pkin(7));
t109 = t70 * t119;
t131 = sin(pkin(8));
t135 = cos(pkin(8));
t133 = cos(pkin(14));
t112 = t133 * t132;
t136 = cos(pkin(7));
t104 = t136 * t112;
t129 = sin(pkin(14));
t111 = t132 * t129;
t143 = cos(qJ(2));
t59 = -t141 * t104 - t143 * t111;
t51 = t135 * t109 - t59 * t131;
t145 = -t71 * mrSges(7,1) - t74 * mrSges(7,2) + mrSges(5,2) - mrSges(6,3);
t72 = sin(qJ(5));
t75 = cos(qJ(5));
t158 = pkin(4) * t157 - t148 * t72 + t159 * t75 + mrSges(5,1);
t108 = t70 * t112;
t130 = sin(pkin(13));
t137 = cos(pkin(6));
t116 = t137 * t130;
t134 = cos(pkin(13));
t64 = -t141 * t116 + t134 * t143;
t97 = t143 * t116 + t134 * t141;
t91 = t97 * t133;
t78 = -t130 * t108 + t64 * t129 + t136 * t91;
t114 = t136 * t132;
t85 = -t130 * t114 - t97 * t70;
t155 = t85 * t131 + t78 * t135;
t117 = t137 * t134;
t63 = t141 * t117 + t130 * t143;
t96 = -t143 * t117 + t130 * t141;
t89 = t96 * t133;
t79 = t134 * t108 + t63 * t129 + t136 * t89;
t84 = t134 * t114 - t96 * t70;
t154 = t84 * t131 + t79 * t135;
t126 = t70 * t137;
t83 = -t143 * t104 + t141 * t111 - t133 * t126;
t120 = t143 * t132;
t98 = t70 * t120 - t137 * t136;
t153 = t98 * t131 + t83 * t135;
t152 = -t70 * mrSges(4,3) + mrSges(3,2);
t149 = m(4) + m(5) + t157;
t146 = -t157 * pkin(11) + t145;
t142 = cos(qJ(4));
t139 = pkin(2) * t120 + qJ(3) * t109;
t138 = qJ(3) * t70;
t125 = t70 * t135;
t124 = t70 * t131;
t123 = -t96 * pkin(2) + t63 * t138;
t122 = -t97 * pkin(2) + t64 * t138;
t121 = t135 * t142;
t115 = t136 * t133;
t113 = t136 * t129;
t110 = t142 * t124;
t107 = t70 * t111;
t103 = t136 * t111;
t101 = t131 * t109;
t60 = -t141 * t103 + t143 * t112;
t100 = t60 * pkin(3) + t51 * pkin(10) + t139;
t88 = t96 * t129;
t46 = -t63 * t115 + t88;
t33 = t63 * t125 - t46 * t131;
t90 = t97 * t129;
t48 = -t64 * t115 + t90;
t34 = t64 * t125 - t48 * t131;
t47 = -t63 * t113 - t89;
t93 = t47 * pkin(3) + t33 * pkin(10) + t123;
t49 = -t64 * t113 - t91;
t92 = t49 * pkin(3) + t34 * pkin(10) + t122;
t73 = sin(qJ(4));
t55 = t143 * t103 + t141 * t112 + t129 * t126;
t41 = t130 * t107 + t64 * t133 - t136 * t90;
t40 = -t134 * t107 + t63 * t133 - t136 * t88;
t39 = t83 * t131 - t98 * t135;
t36 = t60 * t142 + (t135 * t59 + t101) * t73;
t35 = -t142 * t101 - t59 * t121 + t60 * t73;
t29 = t78 * t131 - t85 * t135;
t28 = t79 * t131 - t84 * t135;
t27 = t55 * t142 - t153 * t73;
t26 = t153 * t142 + t55 * t73;
t22 = t49 * t142 + (t64 * t124 + t135 * t48) * t73;
t21 = -t64 * t110 - t48 * t121 + t49 * t73;
t20 = t47 * t142 + (t63 * t124 + t135 * t46) * t73;
t19 = -t63 * t110 - t46 * t121 + t47 * t73;
t16 = t41 * t142 - t155 * t73;
t15 = t155 * t142 + t41 * t73;
t14 = t40 * t142 - t154 * t73;
t13 = t154 * t142 + t40 * t73;
t10 = t27 * t75 + t39 * t72;
t4 = t16 * t75 + t29 * t72;
t2 = t14 * t75 + t28 * t72;
t1 = [(-m(2) - m(3) - t149) * g(3) (-m(4) * t139 - m(5) * t100 - mrSges(3,1) * t120 - t60 * mrSges(4,1) - t36 * mrSges(5,1) + mrSges(3,2) * t119 - t59 * mrSges(4,2) - mrSges(4,3) * t109 - t51 * mrSges(5,3) - t157 * (t36 * pkin(4) + t35 * pkin(11) + t100) - t159 * (t36 * t75 + t51 * t72) + t145 * t35 + t148 * (t36 * t72 - t51 * t75)) * g(3) + (-m(4) * t123 - m(5) * t93 + t96 * mrSges(3,1) - t47 * mrSges(4,1) - t20 * mrSges(5,1) - t46 * mrSges(4,2) - t33 * mrSges(5,3) - t157 * (t20 * pkin(4) + t19 * pkin(11) + t93) + t152 * t63 + t148 * (t20 * t72 - t33 * t75) - t159 * (t20 * t75 + t33 * t72) + t145 * t19) * g(2) + (-m(4) * t122 - m(5) * t92 + t97 * mrSges(3,1) - t49 * mrSges(4,1) - t22 * mrSges(5,1) - t48 * mrSges(4,2) - t34 * mrSges(5,3) - t157 * (t22 * pkin(4) + t21 * pkin(11) + t92) + t148 * (t22 * t72 - t34 * t75) + t152 * t64 - t159 * (t22 * t75 + t34 * t72) + t145 * t21) * g(1) (g(1) * t85 + g(2) * t84 + g(3) * t98) * t149 (t146 * t27 + t158 * t26) * g(3) + (t158 * t13 + t146 * t14) * g(2) + (t146 * t16 + t158 * t15) * g(1) (-t159 * (-t27 * t72 + t39 * t75) + t148 * t10) * g(3) + (t148 * t2 - t159 * (-t14 * t72 + t28 * t75)) * g(2) + (t148 * t4 - t159 * (-t16 * t72 + t29 * t75)) * g(1), -g(1) * ((t15 * t74 - t4 * t71) * mrSges(7,1) + (-t15 * t71 - t4 * t74) * mrSges(7,2)) - g(2) * ((t13 * t74 - t2 * t71) * mrSges(7,1) + (-t13 * t71 - t2 * t74) * mrSges(7,2)) - g(3) * ((-t10 * t71 + t26 * t74) * mrSges(7,1) + (-t10 * t74 - t26 * t71) * mrSges(7,2))];
taug  = t1(:);
