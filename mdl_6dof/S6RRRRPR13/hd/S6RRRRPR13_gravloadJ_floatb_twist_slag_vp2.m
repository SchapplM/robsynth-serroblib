% Calculate Gravitation load on the joints for
% S6RRRRPR13
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% g [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [11x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,alpha2,d1,d2,d3,d4,d6]';
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
% Datum: 2018-11-23 18:23
% Revision: 76f9d5e39f14dc242b53c0d9d3d9db48bd8f37c0
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für mechatronische Systeme, Universität Hannover

function taug = S6RRRRPR13_gravloadJ_floatb_twist_slag_vp2(qJ, g, ...
  pkin, m, mrSges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(3,1),zeros(11,1),zeros(7,1),zeros(7,3)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRRRPR13_gravloadJ_floatb_twist_slag_vp2: qJ has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6RRRRPR13_gravloadJ_floatb_twist_slag_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [11 1]), ...
  'S6RRRRPR13_gravloadJ_floatb_twist_slag_vp2: pkin has to be [11x1] (double)');
assert( isreal(m) && all(size(m) == [7 1]), ...
  'S6RRRRPR13_gravloadJ_floatb_twist_slag_vp2: m has to be [7x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [7,3]), ...
  'S6RRRRPR13_gravloadJ_floatb_twist_slag_vp2: mrSges has to be [7x3] (double)');

%% Symbolic Calculation
% From joint_gravload_floatb_twist_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2018-11-23 18:22:38
% EndTime: 2018-11-23 18:22:40
% DurationCPUTime: 1.72s
% Computational Cost: add. (2260->176), mult. (2844->230), div. (0->0), fcn. (2970->16), ass. (0->101)
t173 = -m(7) * pkin(5) - mrSges(5,1) - mrSges(6,1);
t165 = mrSges(5,2) - mrSges(6,3);
t82 = sin(qJ(6));
t86 = cos(qJ(6));
t154 = -t82 * mrSges(7,1) - t86 * mrSges(7,2) + t165;
t153 = -t86 * mrSges(7,1) + t82 * mrSges(7,2) + t173;
t83 = sin(qJ(4));
t87 = cos(qJ(4));
t172 = pkin(4) * t87 + qJ(5) * t83;
t169 = mrSges(7,3) - mrSges(6,2) - mrSges(5,3);
t137 = sin(pkin(6));
t147 = cos(qJ(1));
t115 = t147 * t137;
t146 = cos(qJ(3));
t135 = pkin(6) - qJ(2);
t119 = sin(t135);
t134 = pkin(6) + qJ(2);
t118 = sin(t134);
t77 = t118 / 0.2e1;
t104 = t77 - t119 / 0.2e1;
t145 = sin(qJ(1));
t88 = cos(qJ(2));
t128 = t145 * t88;
t58 = t104 * t147 + t128;
t84 = sin(qJ(3));
t32 = -t84 * t115 + t146 * t58;
t85 = sin(qJ(2));
t113 = cos(t134) / 0.2e1;
t120 = cos(t135);
t92 = t120 / 0.2e1 + t113;
t57 = t145 * t85 - t147 * t92;
t7 = t32 * t83 - t57 * t87;
t8 = t32 * t87 + t57 * t83;
t168 = m(6) + m(7);
t166 = mrSges(3,2) - mrSges(4,3);
t121 = -t146 * t115 - t58 * t84;
t163 = t172 * t121;
t114 = t137 * t145;
t129 = t147 * t88;
t61 = -t104 * t145 + t129;
t35 = -t114 * t146 + t61 * t84;
t162 = t172 * t35;
t71 = t113 - t120 / 0.2e1;
t81 = cos(pkin(6));
t55 = t146 * t81 + t71 * t84;
t161 = t172 * t55;
t160 = mrSges(4,1) * t146 - t84 * mrSges(4,2) + mrSges(3,1);
t159 = t168 * pkin(4) - t153;
t158 = mrSges(4,2) - m(7) * (pkin(10) - pkin(11)) + t169;
t157 = t153 * t87 + t154 * t83 - mrSges(4,1);
t156 = m(7) * pkin(11) + t169;
t149 = pkin(10) * t35;
t148 = t121 * pkin(10);
t143 = t57 * t84;
t60 = t145 * t92 + t147 * t85;
t141 = t60 * t84;
t78 = t119 / 0.2e1;
t70 = t77 + t78;
t140 = t70 * t84;
t139 = t147 * pkin(1) + pkin(8) * t114;
t132 = t146 * pkin(3);
t131 = t83 * t146;
t130 = t87 * t146;
t103 = t78 - t118 / 0.2e1;
t59 = -t103 * t147 + t128;
t127 = -t57 * pkin(2) + pkin(9) * t59;
t62 = t103 * t145 + t129;
t126 = -t60 * pkin(2) + pkin(9) * t62;
t125 = t70 * pkin(2) - pkin(9) * t71;
t26 = t121 * pkin(3);
t124 = pkin(10) * t32 + t26;
t28 = t35 * pkin(3);
t36 = t114 * t84 + t146 * t61;
t123 = pkin(10) * t36 - t28;
t50 = t55 * pkin(3);
t56 = -t146 * t71 + t81 * t84;
t122 = pkin(10) * t56 + t50;
t116 = -pkin(1) * t145 + pkin(8) * t115;
t108 = -pkin(10) * t143 - t57 * t132 + t127;
t107 = -pkin(10) * t141 - t60 * t132 + t126;
t106 = t61 * pkin(2) + pkin(9) * t60 + t139;
t105 = pkin(10) * t140 + t70 * t132 + t125;
t102 = t36 * pkin(3) + t106;
t99 = -t58 * pkin(2) - t57 * pkin(9) + t116;
t95 = -pkin(3) * t32 + t99;
t11 = t36 * t83 - t60 * t87;
t12 = t36 * t87 + t60 * t83;
t93 = t12 * pkin(4) + qJ(5) * t11 + t102;
t91 = -t168 * qJ(5) + t154;
t90 = -pkin(4) * t8 - qJ(5) * t7 + t95;
t38 = t130 * t70 - t71 * t83;
t37 = t131 * t70 + t71 * t87;
t21 = t56 * t87 - t70 * t83;
t20 = t56 * t83 + t70 * t87;
t18 = -t130 * t60 + t62 * t83;
t17 = -t131 * t60 - t62 * t87;
t16 = -t130 * t57 + t59 * t83;
t15 = -t131 * t57 - t59 * t87;
t2 = t11 * t82 + t12 * t86;
t1 = t11 * t86 - t12 * t82;
t3 = [(-t147 * mrSges(2,1) + t145 * mrSges(2,2) - m(3) * t139 - t61 * mrSges(3,1) - mrSges(3,3) * t114 - m(4) * t106 - t36 * mrSges(4,1) - m(5) * (t102 + t149) - m(6) * (t93 + t149) - m(7) * t93 - t2 * mrSges(7,1) - t1 * mrSges(7,2) + t166 * t60 + t173 * t12 + t165 * t11 + t158 * t35) * g(2) + (t145 * mrSges(2,1) + t147 * mrSges(2,2) - m(3) * t116 + t58 * mrSges(3,1) - mrSges(3,3) * t115 - m(4) * t99 + t32 * mrSges(4,1) - m(5) * (t95 + t148) - m(6) * (t90 + t148) - m(7) * t90 - t166 * t57 - t154 * t7 - t153 * t8 + t158 * t121) * g(1) (-m(4) * t125 - m(5) * t105 - t168 * (t38 * pkin(4) + qJ(5) * t37 + t105) - t166 * t71 - t160 * t70 + t153 * t38 + t154 * t37 + t156 * t140) * g(3) + (-m(4) * t127 - m(5) * t108 - t168 * (t16 * pkin(4) + qJ(5) * t15 + t108) + t166 * t59 + t160 * t57 + t153 * t16 + t154 * t15 - t156 * t143) * g(2) + (-m(4) * t126 - m(5) * t107 - t168 * (t18 * pkin(4) + qJ(5) * t17 + t107) + t166 * t62 + t160 * t60 + t153 * t18 + t154 * t17 - t156 * t141) * g(1) (-m(5) * t122 - m(6) * (t122 + t161) - m(7) * (t50 + t161) + t158 * t56 + t157 * t55) * g(3) + (-m(5) * t124 - m(6) * (t124 + t163) - m(7) * (t26 + t163) + t158 * t32 + t157 * t121) * g(2) + (-m(5) * t123 - m(6) * (t123 - t162) - m(7) * (-t28 - t162) + t158 * t36 - t157 * t35) * g(1) (t159 * t20 + t21 * t91) * g(3) + (t159 * t7 + t8 * t91) * g(2) + (t159 * t11 + t12 * t91) * g(1), t168 * (-g(1) * t11 - g(2) * t7 - g(3) * t20) -g(1) * (mrSges(7,1) * t1 - mrSges(7,2) * t2) - g(2) * ((t7 * t86 - t8 * t82) * mrSges(7,1) + (-t7 * t82 - t8 * t86) * mrSges(7,2)) - g(3) * ((t20 * t86 - t21 * t82) * mrSges(7,1) + (-t20 * t82 - t21 * t86) * mrSges(7,2))];
taug  = t3(:);
