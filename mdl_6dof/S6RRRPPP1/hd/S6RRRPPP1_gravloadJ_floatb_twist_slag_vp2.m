% Calculate Gravitation load on the joints for
% S6RRRPPP1
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% g [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [10x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,alpha4,d1,d2,d3,theta4]';
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
% Datum: 2019-03-09 15:20
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function taug = S6RRRPPP1_gravloadJ_floatb_twist_slag_vp2(qJ, g, ...
  pkin, m, mrSges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(3,1),zeros(10,1),zeros(7,1),zeros(7,3)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRRPPP1_gravloadJ_floatb_twist_slag_vp2: qJ has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6RRRPPP1_gravloadJ_floatb_twist_slag_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S6RRRPPP1_gravloadJ_floatb_twist_slag_vp2: pkin has to be [10x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6RRRPPP1_gravloadJ_floatb_twist_slag_vp2: m has to be [7x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [7,3]), ...
  'S6RRRPPP1_gravloadJ_floatb_twist_slag_vp2: mrSges has to be [7x3] (double)');

%% Symbolic Calculation
% From gravload_joint_floatb_twist_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 15:14:34
% EndTime: 2019-03-09 15:14:39
% DurationCPUTime: 1.80s
% Computational Cost: add. (636->171), mult. (1717->227), div. (0->0), fcn. (1973->10), ass. (0->97)
t73 = sin(pkin(6));
t120 = qJ(4) * t73;
t76 = sin(qJ(3));
t109 = t76 * t120;
t75 = cos(pkin(6));
t119 = qJ(4) * t75;
t77 = sin(qJ(2));
t79 = cos(qJ(3));
t80 = cos(qJ(2));
t167 = t80 * t119 + t77 * (-pkin(3) * t79 - pkin(2) - t109);
t162 = -mrSges(5,3) - mrSges(6,1);
t166 = t162 - mrSges(7,1);
t130 = t77 * t79;
t74 = cos(pkin(10));
t112 = t74 * t130;
t72 = sin(pkin(10));
t134 = t76 * t77;
t94 = t75 * t134 + t73 * t80;
t147 = -t72 * t94 + t112;
t156 = m(7) * pkin(5) - t166;
t40 = t134 * t73 - t75 * t80;
t98 = t79 * mrSges(4,1) - t76 * mrSges(4,2);
t164 = t40 * t156 + (mrSges(3,2) - mrSges(4,3)) * t80 + (m(4) * pkin(2) + mrSges(3,1) + t98) * t77;
t154 = m(6) + m(7);
t163 = mrSges(2,2) - mrSges(3,3);
t69 = t80 * pkin(2);
t107 = -pkin(1) - t69;
t78 = sin(qJ(1));
t128 = t78 * t80;
t81 = cos(qJ(1));
t42 = t128 * t76 + t79 * t81;
t125 = t81 * t76;
t127 = t79 * t80;
t43 = t127 * t78 - t125;
t70 = t81 * pkin(8);
t161 = ((-pkin(9) - t119) * t77 + t107) * t78 - t43 * pkin(3) - t120 * t42 + t70;
t126 = t80 * t81;
t63 = pkin(9) * t126;
t159 = t167 * t81 + t63;
t61 = pkin(9) * t128;
t158 = t167 * t78 + t61;
t100 = t80 * mrSges(3,1) - mrSges(3,2) * t77;
t157 = t77 * mrSges(4,3) + t100;
t148 = mrSges(5,2) - mrSges(6,3) - mrSges(7,2);
t145 = m(7) * qJ(6) + mrSges(5,1) - mrSges(6,2) + mrSges(7,3);
t136 = t74 * t75;
t10 = t136 * t43 - t42 * t72;
t143 = t72 * t75;
t11 = -t143 * t43 - t42 * t74;
t153 = t11 * pkin(4) + qJ(5) * t10;
t44 = -t125 * t80 + t78 * t79;
t45 = t126 * t79 + t78 * t76;
t12 = t136 * t45 + t44 * t72;
t13 = -t143 * t45 + t44 * t74;
t149 = t13 * pkin(4) + qJ(5) * t12;
t146 = mrSges(4,2) + (-m(7) * (pkin(5) + qJ(4)) + t166) * t73;
t67 = t77 * pkin(9);
t144 = t43 * t72;
t142 = t72 * t79;
t139 = t73 * t77;
t138 = t73 * t79;
t135 = t74 * t76;
t133 = t76 * t80;
t131 = t77 * t78;
t129 = t77 * t81;
t122 = t69 + t67;
t121 = t81 * pkin(1) + t78 * pkin(8);
t116 = pkin(3) * t134;
t114 = t73 * t131;
t113 = t73 * t130;
t110 = t74 * t139;
t60 = t77 * t119;
t105 = pkin(2) * t126 + pkin(9) * t129 + t121;
t36 = t42 * pkin(3);
t104 = t120 * t43 - t36;
t38 = t44 * pkin(3);
t103 = t120 * t45 + t38;
t32 = -t112 * t75 + t134 * t72;
t33 = (t142 * t75 + t135) * t77;
t55 = qJ(4) * t113;
t102 = -t33 * pkin(4) - qJ(5) * t32 + t55;
t101 = pkin(3) * t127 + t80 * t109 + t122 + t60;
t25 = t131 * t75 + t42 * t73;
t5 = -t72 * t114 + t143 * t42 - t43 * t74;
t85 = t45 * pkin(3) - t120 * t44 + t81 * t60 + t105;
t21 = t130 * t72 + t74 * t94;
t26 = t129 * t75 - t44 * t73;
t23 = t127 * t74 - t133 * t143 + t139 * t72;
t22 = -t110 + (t135 * t75 + t142) * t80;
t19 = t147 * t81;
t18 = t21 * t81;
t17 = t147 * t78;
t16 = t21 * t78;
t7 = t45 * t74 + (t129 * t73 + t44 * t75) * t72;
t6 = -t110 * t81 - t136 * t44 + t45 * t72;
t4 = t110 * t78 - t136 * t42 - t144;
t1 = [(-m(3) * t121 - m(4) * t105 - m(5) * t85 - t45 * mrSges(4,1) - t44 * mrSges(4,2) - mrSges(4,3) * t129 - t154 * (t7 * pkin(4) + qJ(5) * t6 + t85) + (-mrSges(2,1) - t100) * t81 + t163 * t78 - t145 * t7 + t148 * t6 - t156 * t26) * g(2) + (t43 * mrSges(4,1) - t42 * mrSges(4,2) - t161 * m(5) + t163 * t81 + (-m(3) - m(4)) * t70 - t154 * (t5 * pkin(4) + qJ(5) * t4 + t161) + (mrSges(2,1) + m(3) * pkin(1) - m(4) * (t107 - t67) + t157) * t78 - t145 * t5 + t148 * t4 + t156 * t25) * g(1) (-m(4) * t122 - m(5) * t101 - t80 * t98 - t154 * (t23 * pkin(4) + qJ(5) * t22 + t101) - t156 * (t133 * t73 + t75 * t77) - t145 * t23 + t148 * t22 - t157) * g(3) + (-m(4) * t61 - t158 * m(5) - t154 * (-t17 * pkin(4) - qJ(5) * t16 + t158) + t145 * t17 - t148 * t16 + t164 * t78) * g(2) + (-m(4) * t63 - t159 * m(5) - t154 * (-t19 * pkin(4) - t18 * qJ(5) + t159) + t145 * t19 - t148 * t18 + t164 * t81) * g(1) (-m(5) * (t55 - t116) - m(6) * (t102 - t116) - m(7) * t102 + (mrSges(4,1) * t76 + mrSges(4,2) * t79 - m(7) * (-pkin(3) * t76 + pkin(5) * t138) - mrSges(7,1) * t138) * t77 + t162 * t113 + t145 * t33 - t148 * t32) * g(3) + (mrSges(4,1) * t42 - m(5) * t104 - m(6) * (t104 + t153) - m(7) * (-t36 + t153) - t145 * t11 + t148 * t10 + t146 * t43) * g(2) + (-mrSges(4,1) * t44 - m(5) * t103 - m(6) * (t103 + t149) - m(7) * (t38 + t149) - t145 * t13 + t148 * t12 + t146 * t45) * g(1) (m(5) + t154) * (-g(1) * t26 - g(2) * t25 - g(3) * t40) t154 * (-g(1) * t6 - g(2) * (t144 + (t42 * t75 - t114) * t74) - g(3) * t21) (-g(1) * t7 + g(2) * t5 - g(3) * t147) * m(7)];
taug  = t1(:);
