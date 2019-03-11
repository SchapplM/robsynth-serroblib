% Calculate Gravitation load on the joints for
% S6RPRRRR10
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% g [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [13x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,alpha2,alpha3,d1,d3,d4,d5,d6,theta2]';
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
% Datum: 2019-03-09 07:37
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function taug = S6RPRRRR10_gravloadJ_floatb_twist_slag_vp2(qJ, g, ...
  pkin, m, mrSges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(3,1),zeros(13,1),zeros(7,1),zeros(7,3)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPRRRR10_gravloadJ_floatb_twist_slag_vp2: qJ has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6RPRRRR10_gravloadJ_floatb_twist_slag_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [13 1]), ...
  'S6RPRRRR10_gravloadJ_floatb_twist_slag_vp2: pkin has to be [13x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6RPRRRR10_gravloadJ_floatb_twist_slag_vp2: m has to be [7x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [7,3]), ...
  'S6RPRRRR10_gravloadJ_floatb_twist_slag_vp2: mrSges has to be [7x3] (double)');

%% Symbolic Calculation
% From gravload_joint_floatb_twist_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 07:27:46
% EndTime: 2019-03-09 07:27:50
% DurationCPUTime: 1.37s
% Computational Cost: add. (1249->152), mult. (2949->213), div. (0->0), fcn. (3718->16), ass. (0->80)
t162 = mrSges(6,2) - mrSges(7,3);
t78 = sin(qJ(6));
t81 = cos(qJ(6));
t157 = mrSges(7,1) * t81 - mrSges(7,2) * t78 + mrSges(6,1);
t133 = cos(pkin(13));
t135 = cos(pkin(6));
t111 = t135 * t133;
t131 = sin(pkin(13));
t149 = sin(qJ(1));
t83 = cos(qJ(1));
t100 = -t111 * t83 + t149 * t131;
t132 = sin(pkin(7));
t77 = sin(pkin(6));
t123 = t77 * t132;
t134 = cos(pkin(7));
t169 = t100 * t134 + t83 * t123;
t124 = t77 * t134;
t52 = -t100 * t132 + t83 * t124;
t79 = sin(qJ(4));
t143 = t52 * t79;
t150 = cos(qJ(3));
t109 = t135 * t131;
t61 = t109 * t83 + t133 * t149;
t80 = sin(qJ(3));
t39 = -t150 * t61 + t169 * t80;
t82 = cos(qJ(4));
t168 = t39 * t82 + t143;
t159 = t39 * t79 - t52 * t82;
t76 = qJ(4) + qJ(5);
t73 = sin(t76);
t74 = cos(t76);
t13 = t39 * t73 - t52 * t74;
t167 = t39 * t74 + t52 * t73;
t163 = m(7) * pkin(5) + t157;
t122 = -m(7) * pkin(12) + t162;
t95 = t111 * t149 + t131 * t83;
t161 = t149 * t124 + t95 * t132;
t155 = t149 * t123 - t95 * t134;
t62 = -t109 * t149 + t133 * t83;
t41 = t62 * t150 + t155 * t80;
t19 = t161 * t82 - t41 * t79;
t108 = t134 * t133;
t110 = t135 * t132;
t51 = t80 * t110 + (t80 * t108 + t131 * t150) * t77;
t60 = -t123 * t133 + t134 * t135;
t160 = -t51 * t79 + t60 * t82;
t36 = t169 * t150 + t61 * t80;
t156 = -m(6) - m(7);
t112 = -m(5) * pkin(10) + mrSges(4,2) - mrSges(5,3) - mrSges(6,3);
t151 = -t78 * mrSges(7,1) - t81 * mrSges(7,2) + t112;
t152 = m(5) * pkin(3) + t82 * mrSges(5,1) - t79 * mrSges(5,2) - t122 * t73 + t163 * t74 + mrSges(4,1);
t140 = t77 * t83;
t128 = t77 * t149;
t136 = t83 * pkin(1) + qJ(2) * t128;
t127 = t13 * pkin(5) - pkin(12) * t167;
t17 = -t161 * t74 + t41 * t73;
t18 = t161 * t73 + t41 * t74;
t126 = -t17 * pkin(5) + pkin(12) * t18;
t30 = -t51 * t73 + t60 * t74;
t31 = t51 * t74 + t60 * t73;
t125 = t30 * pkin(5) + pkin(12) * t31;
t121 = t159 * pkin(4);
t120 = t19 * pkin(4);
t119 = t160 * pkin(4);
t116 = -pkin(1) * t149 + qJ(2) * t140;
t105 = -t157 * t13 - t162 * t167;
t104 = t157 * t17 + t162 * t18;
t101 = -t157 * t30 + t162 * t31;
t91 = -t61 * pkin(2) + t52 * pkin(9) + t116;
t90 = t62 * pkin(2) + t161 * pkin(9) + t136;
t40 = -t155 * t150 + t62 * t80;
t72 = pkin(4) * t82 + pkin(3);
t84 = -pkin(11) - pkin(10);
t86 = t161 * t79;
t87 = pkin(4) * t86 - t40 * t84 + t41 * t72 + t90;
t50 = t77 * t131 * t80 + (-t108 * t77 - t110) * t150;
t20 = t41 * t82 + t86;
t2 = t18 * t81 + t40 * t78;
t1 = -t18 * t78 + t40 * t81;
t3 = [(-t83 * mrSges(2,1) + t149 * mrSges(2,2) - m(3) * t136 - t62 * mrSges(3,1) + t95 * mrSges(3,2) - mrSges(3,3) * t128 - m(4) * t90 - t41 * mrSges(4,1) - t161 * mrSges(4,3) - m(5) * (t41 * pkin(3) + t90) - t20 * mrSges(5,1) - t19 * mrSges(5,2) - m(6) * t87 - t18 * mrSges(6,1) - m(7) * (t18 * pkin(5) + t87) - t2 * mrSges(7,1) - t1 * mrSges(7,2) + t122 * t17 + t112 * t40) * g(2) + (t149 * mrSges(2,1) + t83 * mrSges(2,2) - m(3) * t116 + t61 * mrSges(3,1) - t100 * mrSges(3,2) - mrSges(3,3) * t140 - m(4) * t91 - t39 * mrSges(4,1) - t52 * mrSges(4,3) - m(5) * (t39 * pkin(3) + t91) - t168 * mrSges(5,1) + t159 * mrSges(5,2) + t122 * t13 - t163 * t167 - t151 * t36 + t156 * (pkin(4) * t143 + t36 * t84 + t39 * t72 + t91)) * g(1) (-t135 * g(3) + (-t149 * g(1) + t83 * g(2)) * t77) * (m(3) + m(4) + m(5) - t156) (t156 * (-t50 * t72 - t51 * t84) + t151 * t51 + t152 * t50) * g(3) + (t156 * (-t36 * t72 + t39 * t84) - t151 * t39 + t152 * t36) * g(2) + (t156 * (-t40 * t72 - t41 * t84) + t151 * t41 + t152 * t40) * g(1) (-t160 * mrSges(5,1) - (-t51 * t82 - t60 * t79) * mrSges(5,2) - m(6) * t119 - m(7) * (t119 + t125) + t101) * g(3) + (-t159 * mrSges(5,1) - t168 * mrSges(5,2) - m(6) * t121 - m(7) * (t121 + t127) + t105) * g(2) + (-t19 * mrSges(5,1) + t20 * mrSges(5,2) - m(6) * t120 - m(7) * (t120 + t126) + t104) * g(1) (-m(7) * t125 + t101) * g(3) + (-m(7) * t127 + t105) * g(2) + (-m(7) * t126 + t104) * g(1), -g(1) * (mrSges(7,1) * t1 - mrSges(7,2) * t2) - g(2) * ((t167 * t78 + t36 * t81) * mrSges(7,1) + (t167 * t81 - t36 * t78) * mrSges(7,2)) - g(3) * ((-t31 * t78 + t50 * t81) * mrSges(7,1) + (-t31 * t81 - t50 * t78) * mrSges(7,2))];
taug  = t3(:);
