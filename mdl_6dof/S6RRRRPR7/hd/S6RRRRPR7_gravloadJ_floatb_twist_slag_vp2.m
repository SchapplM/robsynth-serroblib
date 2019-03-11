% Calculate Gravitation load on the joints for
% S6RRRRPR7
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% g [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [12x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,alpha2,d1,d2,d3,d4,d6,theta5]';
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
% Datum: 2019-03-09 22:37
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function taug = S6RRRRPR7_gravloadJ_floatb_twist_slag_vp2(qJ, g, ...
  pkin, m, mrSges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(3,1),zeros(12,1),zeros(7,1),zeros(7,3)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRRRPR7_gravloadJ_floatb_twist_slag_vp2: qJ has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6RRRRPR7_gravloadJ_floatb_twist_slag_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [12 1]), ...
  'S6RRRRPR7_gravloadJ_floatb_twist_slag_vp2: pkin has to be [12x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6RRRRPR7_gravloadJ_floatb_twist_slag_vp2: m has to be [7x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [7,3]), ...
  'S6RRRRPR7_gravloadJ_floatb_twist_slag_vp2: mrSges has to be [7x3] (double)');

%% Symbolic Calculation
% From gravload_joint_floatb_twist_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 22:26:44
% EndTime: 2019-03-09 22:26:48
% DurationCPUTime: 1.41s
% Computational Cost: add. (950->163), mult. (1327->214), div. (0->0), fcn. (1503->14), ass. (0->78)
t161 = -mrSges(7,3) + mrSges(6,2);
t82 = sin(qJ(6));
t86 = cos(qJ(6));
t160 = -mrSges(7,1) * t86 + mrSges(7,2) * t82 - mrSges(6,1);
t162 = m(7) * pkin(5) - t160;
t111 = -m(7) * pkin(11) + t161;
t90 = -pkin(10) - pkin(9);
t92 = -m(4) * pkin(9) + m(5) * t90 + mrSges(3,2) - mrSges(4,3) - mrSges(5,3) - mrSges(6,3);
t155 = mrSges(7,1) * t82 + mrSges(7,2) * t86 - t92;
t157 = m(6) + m(7);
t125 = cos(pkin(6));
t81 = sin(pkin(6));
t84 = sin(qJ(2));
t139 = t81 * t84;
t80 = qJ(3) + qJ(4);
t75 = sin(t80);
t76 = cos(t80);
t159 = t125 * t76 - t75 * t139;
t85 = sin(qJ(1));
t138 = t81 * t85;
t115 = t85 * t125;
t88 = cos(qJ(2));
t89 = cos(qJ(1));
t54 = -t115 * t84 + t88 * t89;
t23 = t76 * t138 - t54 * t75;
t87 = cos(qJ(3));
t77 = t87 * pkin(3);
t73 = t77 + pkin(2);
t83 = sin(qJ(3));
t154 = m(4) * pkin(2) + m(5) * t73 + t87 * mrSges(4,1) + t76 * mrSges(5,1) - t83 * mrSges(4,2) - t75 * mrSges(5,2) + mrSges(3,1);
t74 = pkin(12) + t80;
t70 = sin(t74);
t71 = cos(t74);
t149 = -t111 * t70 + t162 * t71 + t154;
t156 = -m(5) * pkin(3) - mrSges(4,1);
t38 = t125 * t71 - t139 * t70;
t39 = t125 * t70 + t139 * t71;
t152 = -t159 * mrSges(5,1) - (-t125 * t75 - t139 * t76) * mrSges(5,2) + t161 * t39 + t160 * t38;
t17 = -t138 * t71 + t54 * t70;
t18 = t138 * t70 + t54 * t71;
t24 = t138 * t75 + t54 * t76;
t151 = -t23 * mrSges(5,1) + t24 * mrSges(5,2) - t160 * t17 + t161 * t18;
t135 = t81 * t89;
t114 = t89 * t125;
t52 = t114 * t84 + t85 * t88;
t101 = -t135 * t76 - t52 * t75;
t13 = -t71 * t135 - t52 * t70;
t14 = -t70 * t135 + t52 * t71;
t63 = t75 * t135;
t150 = -t101 * mrSges(5,1) - (-t52 * t76 + t63) * mrSges(5,2) + t161 * t14 + t160 * t13;
t145 = pkin(3) * t83;
t137 = t81 * t87;
t136 = t81 * t88;
t61 = pkin(4) * t75 + t145;
t62 = pkin(4) * t76 + t77;
t130 = t62 * t138 - t54 * t61;
t127 = t125 * t62 - t61 * t139;
t126 = t89 * pkin(1) + pkin(8) * t138;
t119 = t85 * pkin(1) - pkin(8) * t135;
t118 = t13 * pkin(5) + pkin(11) * t14;
t117 = -t17 * pkin(5) + pkin(11) * t18;
t116 = t38 * pkin(5) + pkin(11) * t39;
t112 = -m(5) * t145 - mrSges(3,3);
t110 = t23 * pkin(4);
t109 = -t135 * t62 - t52 * t61;
t53 = t115 * t88 + t89 * t84;
t59 = pkin(2) + t62;
t79 = -qJ(5) + t90;
t108 = t61 * t138 - t53 * t79 + t54 * t59 + t126;
t103 = t159 * pkin(4);
t25 = t137 * t85 - t54 * t83;
t96 = t101 * pkin(4);
t66 = t83 * t135;
t51 = -t114 * t88 + t84 * t85;
t26 = t138 * t83 + t54 * t87;
t2 = t18 * t86 + t53 * t82;
t1 = -t18 * t82 + t53 * t86;
t3 = [(-t89 * mrSges(2,1) - m(3) * t126 - t54 * mrSges(3,1) - m(4) * (pkin(2) * t54 + t126) - t26 * mrSges(4,1) - t25 * mrSges(4,2) - m(5) * (t54 * t73 + t126) - t24 * mrSges(5,1) - t23 * mrSges(5,2) - m(6) * t108 - t18 * mrSges(6,1) - m(7) * (pkin(5) * t18 + t108) - t2 * mrSges(7,1) - t1 * mrSges(7,2) + (t112 * t81 + mrSges(2,2)) * t85 + t111 * t17 + t92 * t53) * g(2) + (t85 * mrSges(2,1) - t66 * mrSges(4,1) - t63 * mrSges(5,1) + t154 * t52 + (mrSges(2,2) + (-mrSges(4,2) * t87 - mrSges(5,2) * t76 + t112) * t81) * t89 + t111 * t13 + t162 * t14 + t155 * t51 + (m(3) + m(4) + m(5)) * t119 + t157 * (-t61 * t135 - t51 * t79 + t52 * t59 + t119)) * g(1) (-t157 * (-t51 * t59 - t52 * t79) - t155 * t52 + t149 * t51) * g(2) + (-t157 * (-t53 * t59 - t54 * t79) - t155 * t54 + t149 * t53) * g(1) + (-t157 * t59 * t136 + (-t149 * t88 + (t157 * t79 - t155) * t84) * t81) * g(3) (-(-t125 * t83 - t137 * t84) * mrSges(4,2) - m(6) * t127 - m(7) * (t116 + t127) + t156 * (t125 * t87 - t139 * t83) + t152) * g(3) + (-(-t52 * t87 + t66) * mrSges(4,2) - m(6) * t109 - m(7) * (t109 + t118) + t156 * (-t135 * t87 - t52 * t83) + t150) * g(2) + (mrSges(4,2) * t26 - m(6) * t130 - m(7) * (t117 + t130) + t156 * t25 + t151) * g(1) (-m(6) * t103 - m(7) * (t103 + t116) + t152) * g(3) + (-m(6) * t96 - m(7) * (t118 + t96) + t150) * g(2) + (-m(6) * t110 - m(7) * (t110 + t117) + t151) * g(1), t157 * (-g(1) * t53 - g(2) * t51 + g(3) * t136) -g(1) * (mrSges(7,1) * t1 - mrSges(7,2) * t2) - g(2) * ((-t14 * t82 + t51 * t86) * mrSges(7,1) + (-t14 * t86 - t51 * t82) * mrSges(7,2)) - g(3) * ((-t136 * t86 - t39 * t82) * mrSges(7,1) + (t136 * t82 - t39 * t86) * mrSges(7,2))];
taug  = t3(:);
