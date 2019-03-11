% Calculate Gravitation load on the joints for
% S6RRRRPR10
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

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-09 23:11
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function taug = S6RRRRPR10_gravloadJ_floatb_twist_slag_vp2(qJ, g, ...
  pkin, m, mrSges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(3,1),zeros(11,1),zeros(7,1),zeros(7,3)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRRRPR10_gravloadJ_floatb_twist_slag_vp2: qJ has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6RRRRPR10_gravloadJ_floatb_twist_slag_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [11 1]), ...
  'S6RRRRPR10_gravloadJ_floatb_twist_slag_vp2: pkin has to be [11x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6RRRRPR10_gravloadJ_floatb_twist_slag_vp2: m has to be [7x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [7,3]), ...
  'S6RRRRPR10_gravloadJ_floatb_twist_slag_vp2: mrSges has to be [7x3] (double)');

%% Symbolic Calculation
% From gravload_joint_floatb_twist_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 23:00:59
% EndTime: 2019-03-09 23:01:04
% DurationCPUTime: 1.69s
% Computational Cost: add. (880->168), mult. (1482->213), div. (0->0), fcn. (1712->12), ass. (0->91)
t174 = mrSges(5,1) - mrSges(6,2);
t172 = mrSges(7,3) + t174;
t76 = sin(qJ(6));
t80 = cos(qJ(6));
t166 = t76 * mrSges(7,1) + t80 * mrSges(7,2);
t168 = mrSges(5,2) - mrSges(6,3);
t177 = t166 - t168;
t173 = -m(4) * pkin(2) - mrSges(3,1);
t169 = m(4) * pkin(9);
t159 = m(6) + m(7);
t126 = cos(pkin(6));
t75 = sin(pkin(6));
t78 = sin(qJ(2));
t141 = t75 * t78;
t77 = sin(qJ(3));
t81 = cos(qJ(3));
t167 = t126 * t81 - t77 * t141;
t139 = t75 * t81;
t109 = t78 * t126;
t148 = cos(qJ(2));
t149 = cos(qJ(1));
t79 = sin(qJ(1));
t54 = -t79 * t109 + t149 * t148;
t27 = t79 * t139 - t54 * t77;
t165 = -m(7) * pkin(11) - t174;
t164 = t80 * mrSges(7,1) - t76 * mrSges(7,2);
t163 = -m(7) * pkin(5) - mrSges(6,1) + mrSges(3,2) - mrSges(5,3);
t162 = t163 - t164;
t158 = t81 * mrSges(4,1) - t77 * mrSges(4,2);
t74 = qJ(3) + qJ(4);
t71 = sin(t74);
t72 = cos(t74);
t161 = t172 * t72 + t177 * t71 + t158 - t173;
t118 = mrSges(4,3) + t169;
t160 = -t118 + t162;
t45 = -t126 * t72 + t71 * t141;
t46 = t126 * t71 + t72 * t141;
t157 = t172 * t45 - t177 * t46;
t140 = t75 * t79;
t25 = -t72 * t140 + t54 * t71;
t26 = t71 * t140 + t54 * t72;
t156 = t172 * t25 - t177 * t26;
t115 = t75 * t149;
t52 = t149 * t109 + t79 * t148;
t21 = t72 * t115 + t52 * t71;
t22 = -t71 * t115 + t52 * t72;
t155 = t172 * t21 - t177 * t22;
t154 = pkin(11) * t25;
t153 = pkin(11) * t45;
t152 = t21 * pkin(11);
t99 = t126 * t148;
t51 = -t149 * t99 + t78 * t79;
t147 = t51 * t72;
t53 = t149 * t78 + t79 * t99;
t146 = t53 * t72;
t70 = pkin(3) * t81 + pkin(2);
t82 = -pkin(10) - pkin(9);
t132 = -t51 * t70 - t52 * t82;
t131 = -t53 * t70 - t54 * t82;
t128 = t149 * pkin(1) + pkin(8) * t140;
t127 = qJ(5) * t71;
t123 = t77 * t140;
t114 = t75 * t148;
t113 = t76 * t148;
t112 = t80 * t148;
t111 = -pkin(1) * t79 + pkin(8) * t115;
t64 = t77 * t115;
t110 = -t52 * t81 + t64;
t107 = -t21 * pkin(4) + t22 * qJ(5);
t106 = -t25 * pkin(4) + qJ(5) * t26;
t105 = -t45 * pkin(4) + qJ(5) * t46;
t104 = -pkin(4) * t147 - t51 * t127 + t132;
t103 = -pkin(4) * t146 - t53 * t127 + t131;
t102 = t72 * t114;
t101 = t27 * pkin(3);
t100 = pkin(3) * t123 - t53 * t82 + t54 * t70 + t128;
t96 = -mrSges(7,3) + t165;
t95 = t167 * pkin(3);
t93 = -t159 * qJ(5) + t168;
t92 = pkin(3) * t64 + t51 * t82 - t52 * t70 + t111;
t88 = t81 * t115 + t52 * t77;
t87 = t101 + t106;
t86 = -t118 + t163;
t85 = t88 * pkin(3);
t84 = t105 + t95;
t83 = t107 - t85;
t56 = t70 * t114;
t28 = t54 * t81 + t123;
t2 = t25 * t76 + t53 * t80;
t1 = t25 * t80 - t53 * t76;
t3 = [(-t149 * mrSges(2,1) - m(3) * t128 - t54 * mrSges(3,1) - m(4) * (pkin(2) * t54 + t128) - t28 * mrSges(4,1) - t27 * mrSges(4,2) - m(5) * t100 - t2 * mrSges(7,1) - t1 * mrSges(7,2) + (-mrSges(3,3) * t75 + mrSges(2,2)) * t79 + t93 * t25 + t96 * t26 + t86 * t53 - t159 * (t26 * pkin(4) + t100)) * g(2) + (t79 * mrSges(2,1) + t149 * mrSges(2,2) - m(3) * t111 + t52 * mrSges(3,1) - mrSges(3,3) * t115 - m(4) * (-pkin(2) * t52 + t111) - t110 * mrSges(4,1) - t88 * mrSges(4,2) - m(5) * t92 - t96 * t22 - (t93 - t166) * t21 + (t164 - t86) * t51 + t159 * (pkin(4) * t22 - t92)) * g(1) (-m(5) * t56 - mrSges(4,3) * t141 - mrSges(7,3) * t102 - t158 * t114 - t159 * (pkin(4) * t102 + t114 * t127 + t56) + ((-t113 * mrSges(7,1) - t112 * mrSges(7,2)) * t71 + (-t169 + (m(5) + t159) * t82 + t162) * t78 + (t165 * t72 + t168 * t71 + t173) * t148) * t75) * g(3) + (-m(5) * t132 - m(6) * t104 - m(7) * (-pkin(11) * t147 + t104) + t160 * t52 + t161 * t51) * g(2) + (-m(5) * t131 - m(6) * t103 - m(7) * (-pkin(11) * t146 + t103) + t160 * t54 + t161 * t53) * g(1) (-t167 * mrSges(4,1) - (-t126 * t77 - t78 * t139) * mrSges(4,2) - m(5) * t95 - m(6) * t84 - m(7) * (t84 - t153) + t157) * g(3) + (t88 * mrSges(4,1) - t110 * mrSges(4,2) + m(5) * t85 - m(6) * t83 - m(7) * (t83 - t152) + t155) * g(2) + (-mrSges(4,1) * t27 + t28 * mrSges(4,2) - m(5) * t101 - m(6) * t87 - m(7) * (t87 - t154) + t156) * g(1) (-m(6) * t105 - m(7) * (t105 - t153) + t157) * g(3) + (-m(6) * t107 - m(7) * (t107 - t152) + t155) * g(2) + (-m(6) * t106 - m(7) * (t106 - t154) + t156) * g(1), t159 * (-g(1) * t25 - g(2) * t21 - g(3) * t45) -g(1) * (mrSges(7,1) * t1 - mrSges(7,2) * t2) - g(2) * ((t21 * t80 - t51 * t76) * mrSges(7,1) + (-t21 * t76 - t51 * t80) * mrSges(7,2)) - g(3) * ((t75 * t113 + t45 * t80) * mrSges(7,1) + (t75 * t112 - t45 * t76) * mrSges(7,2))];
taug  = t3(:);
