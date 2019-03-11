% Calculate Gravitation load on the joints for
% S6PRRRRR4
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% g [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [13x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,alpha2,alpha3,d2,d3,d4,d5,d6,theta1]';
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
% Datum: 2019-03-09 01:02
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function taug = S6PRRRRR4_gravloadJ_floatb_twist_slag_vp2(qJ, g, ...
  pkin, m, mrSges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(3,1),zeros(13,1),zeros(7,1),zeros(7,3)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6PRRRRR4_gravloadJ_floatb_twist_slag_vp2: qJ has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6PRRRRR4_gravloadJ_floatb_twist_slag_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [13 1]), ...
  'S6PRRRRR4_gravloadJ_floatb_twist_slag_vp2: pkin has to be [13x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6PRRRRR4_gravloadJ_floatb_twist_slag_vp2: m has to be [7x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [7,3]), ...
  'S6PRRRRR4_gravloadJ_floatb_twist_slag_vp2: mrSges has to be [7x3] (double)');

%% Symbolic Calculation
% From gravload_joint_floatb_twist_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 00:54:52
% EndTime: 2019-03-09 00:54:55
% DurationCPUTime: 1.27s
% Computational Cost: add. (1222->152), mult. (2885->231), div. (0->0), fcn. (3624->16), ass. (0->86)
t165 = mrSges(6,2) - mrSges(7,3);
t80 = sin(qJ(6));
t83 = cos(qJ(6));
t160 = -t83 * mrSges(7,1) + t80 * mrSges(7,2) - mrSges(6,1);
t166 = m(7) * pkin(5) - t160;
t121 = -m(7) * pkin(12) + t165;
t131 = sin(pkin(6));
t151 = sin(qJ(2));
t113 = t131 * t151;
t152 = cos(qJ(3));
t133 = cos(pkin(7));
t108 = t133 * t131;
t134 = cos(pkin(6));
t153 = cos(qJ(2));
t79 = sin(pkin(7));
t156 = t153 * t108 + t134 * t79;
t82 = sin(qJ(3));
t49 = t152 * t113 + t156 * t82;
t114 = t153 * t131;
t63 = -t114 * t79 + t133 * t134;
t81 = sin(qJ(4));
t84 = cos(qJ(4));
t163 = -t49 * t81 + t63 * t84;
t130 = sin(pkin(13));
t106 = t131 * t130;
t109 = t134 * t130;
t132 = cos(pkin(13));
t91 = t109 * t153 + t132 * t151;
t157 = t79 * t106 - t91 * t133;
t65 = -t109 * t151 + t132 * t153;
t33 = t65 * t152 + t157 * t82;
t51 = t106 * t133 + t79 * t91;
t162 = -t33 * t81 + t51 * t84;
t107 = t132 * t131;
t110 = t134 * t132;
t90 = -t110 * t153 + t130 * t151;
t158 = t79 * t107 + t90 * t133;
t64 = t110 * t151 + t130 * t153;
t31 = t152 * t64 - t158 * t82;
t50 = -t107 * t133 + t79 * t90;
t161 = -t31 * t81 + t50 * t84;
t154 = -m(6) - m(7);
t94 = -m(5) * pkin(3) - t84 * mrSges(5,1) + t81 * mrSges(5,2) - mrSges(4,1);
t129 = -m(4) - m(5) + t154;
t159 = -pkin(2) * t129 + mrSges(3,1);
t78 = qJ(4) + qJ(5);
t76 = sin(t78);
t77 = cos(t78);
t155 = -t121 * t76 + t166 * t77 - t94;
t92 = -m(5) * pkin(10) - t80 * mrSges(7,1) - t83 * mrSges(7,2) + mrSges(4,2) - mrSges(5,3) - mrSges(6,3);
t144 = t76 * t79;
t143 = t77 * t79;
t105 = t79 * t113;
t135 = pkin(2) * t114 + pkin(9) * t105;
t11 = -t31 * t76 + t50 * t77;
t12 = t31 * t77 + t50 * t76;
t126 = t11 * pkin(5) + pkin(12) * t12;
t13 = -t33 * t76 + t51 * t77;
t14 = t33 * t77 + t51 * t76;
t125 = t13 * pkin(5) + pkin(12) * t14;
t26 = -t49 * t76 + t63 * t77;
t27 = t49 * t77 + t63 * t76;
t124 = t26 * pkin(5) + pkin(12) * t27;
t123 = t82 * t133;
t120 = t161 * pkin(4);
t119 = t162 * pkin(4);
t118 = t163 * pkin(4);
t115 = t133 * t152;
t102 = t160 * t11 + t165 * t12;
t101 = t160 * t13 + t165 * t14;
t100 = t81 * t105;
t97 = t160 * t26 + t165 * t27;
t96 = t151 * t108;
t86 = mrSges(3,2) + (-t84 * mrSges(5,2) - mrSges(4,3) + (pkin(4) * t154 - mrSges(5,1)) * t81 + t129 * pkin(9)) * t79;
t85 = -pkin(11) - pkin(10);
t75 = pkin(4) * t84 + pkin(3);
t59 = t114 * t152 - t82 * t96;
t58 = t114 * t82 + t152 * t96;
t48 = t113 * t82 - t156 * t152;
t41 = -t123 * t65 - t152 * t91;
t40 = t115 * t65 - t82 * t91;
t39 = -t123 * t64 - t152 * t90;
t38 = t115 * t64 - t82 * t90;
t32 = -t157 * t152 + t65 * t82;
t30 = t158 * t152 + t64 * t82;
t1 = [(-m(2) - m(3) + t129) * g(3) (-mrSges(3,1) * t114 + mrSges(3,2) * t113 - m(4) * t135 - t59 * mrSges(4,1) - mrSges(4,3) * t105 - m(5) * (pkin(3) * t59 + t135) - (t59 * t84 + t100) * mrSges(5,1) - (t105 * t84 - t59 * t81) * mrSges(5,2) + t121 * (-t105 * t77 + t59 * t76) - t166 * (t105 * t76 + t59 * t77) + t92 * t58 + t154 * (pkin(4) * t100 - t58 * t85 + t59 * t75 + t135)) * g(3) + (t121 * (-t143 * t64 + t39 * t76) + t94 * t39 - t166 * (t144 * t64 + t39 * t77) + t92 * t38 + t86 * t64 + t159 * t90 + t154 * (-t38 * t85 + t39 * t75)) * g(2) + (t121 * (-t65 * t143 + t41 * t76) + t94 * t41 - t166 * (t65 * t144 + t41 * t77) + t92 * t40 + t86 * t65 + t159 * t91 + t154 * (-t40 * t85 + t41 * t75)) * g(1) (t154 * (-t48 * t75 - t49 * t85) + t92 * t49 + t155 * t48) * g(3) + (t154 * (-t30 * t75 - t31 * t85) + t92 * t31 + t155 * t30) * g(2) + (t154 * (-t32 * t75 - t33 * t85) + t92 * t33 + t155 * t32) * g(1) (-t163 * mrSges(5,1) - (-t49 * t84 - t63 * t81) * mrSges(5,2) - m(6) * t118 - m(7) * (t118 + t124) + t97) * g(3) + (-t161 * mrSges(5,1) - (-t31 * t84 - t50 * t81) * mrSges(5,2) - m(6) * t120 - m(7) * (t120 + t126) + t102) * g(2) + (-t162 * mrSges(5,1) - (-t33 * t84 - t51 * t81) * mrSges(5,2) - m(6) * t119 - m(7) * (t119 + t125) + t101) * g(1) (-m(7) * t124 + t97) * g(3) + (-m(7) * t126 + t102) * g(2) + (-m(7) * t125 + t101) * g(1), -g(1) * ((-t14 * t80 + t32 * t83) * mrSges(7,1) + (-t14 * t83 - t32 * t80) * mrSges(7,2)) - g(2) * ((-t12 * t80 + t30 * t83) * mrSges(7,1) + (-t12 * t83 - t30 * t80) * mrSges(7,2)) - g(3) * ((-t27 * t80 + t48 * t83) * mrSges(7,1) + (-t27 * t83 - t48 * t80) * mrSges(7,2))];
taug  = t1(:);
