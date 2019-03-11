% Calculate Gravitation load on the joints for
% S6RRRPRR11
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% g [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [11x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,alpha2,d1,d2,d3,d5,d6]';
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
% Datum: 2019-03-09 19:37
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function taug = S6RRRPRR11_gravloadJ_floatb_twist_slag_vp2(qJ, g, ...
  pkin, m, mrSges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(3,1),zeros(11,1),zeros(7,1),zeros(7,3)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRRPRR11_gravloadJ_floatb_twist_slag_vp2: qJ has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6RRRPRR11_gravloadJ_floatb_twist_slag_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [11 1]), ...
  'S6RRRPRR11_gravloadJ_floatb_twist_slag_vp2: pkin has to be [11x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6RRRPRR11_gravloadJ_floatb_twist_slag_vp2: m has to be [7x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [7,3]), ...
  'S6RRRPRR11_gravloadJ_floatb_twist_slag_vp2: mrSges has to be [7x3] (double)');

%% Symbolic Calculation
% From gravload_joint_floatb_twist_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 19:24:53
% EndTime: 2019-03-09 19:24:57
% DurationCPUTime: 1.74s
% Computational Cost: add. (845->164), mult. (2114->221), div. (0->0), fcn. (2585->12), ass. (0->80)
t128 = sin(qJ(1));
t129 = cos(qJ(2));
t77 = sin(qJ(2));
t117 = cos(pkin(6));
t130 = cos(qJ(1));
t94 = t117 * t130;
t55 = t128 * t77 - t129 * t94;
t73 = sin(pkin(6));
t114 = t73 * t130;
t56 = t128 * t129 + t77 * t94;
t76 = sin(qJ(3));
t80 = cos(qJ(3));
t31 = t114 * t80 + t56 * t76;
t32 = -t76 * t114 + t56 * t80;
t75 = sin(qJ(5));
t79 = cos(qJ(5));
t6 = t31 * t75 + t32 * t79;
t74 = sin(qJ(6));
t78 = cos(qJ(6));
t159 = t55 * t78 + t6 * t74;
t158 = t55 * t74 - t6 * t78;
t141 = mrSges(3,2) - mrSges(5,2) - mrSges(4,3);
t149 = m(6) + m(7);
t157 = (pkin(9) - pkin(10)) * t149 - mrSges(6,3) - t141;
t148 = mrSges(4,1) + mrSges(5,1);
t146 = mrSges(4,2) - mrSges(5,3);
t143 = -m(7) * pkin(11) + mrSges(6,2) - mrSges(7,3);
t147 = m(7) * pkin(5) + t78 * mrSges(7,1) - t74 * mrSges(7,2) + mrSges(6,1);
t150 = t31 * t79 - t32 * t75;
t156 = t143 * t6 - t147 * t150;
t112 = t73 * t128;
t93 = t117 * t128;
t58 = t129 * t130 - t77 * t93;
t35 = -t112 * t80 + t58 * t76;
t36 = t112 * t76 + t58 * t80;
t12 = t35 * t75 + t36 * t79;
t90 = -t35 * t79 + t36 * t75;
t155 = t12 * t143 + t147 * t90;
t123 = t73 * t77;
t53 = -t117 * t80 + t123 * t76;
t54 = t117 * t76 + t123 * t80;
t20 = t53 * t75 + t54 * t79;
t154 = t143 * t20 - t147 * (t53 * t79 - t54 * t75);
t151 = (-t75 * t80 + t76 * t79) * t143 + mrSges(3,1) + t148 * t80 - t146 * t76 - (-t75 * t76 - t79 * t80) * t147;
t118 = qJ(4) * t76;
t125 = t55 * t80;
t145 = -pkin(3) * t125 - t55 * t118;
t57 = t129 * t93 + t130 * t77;
t124 = t57 * t80;
t144 = -pkin(3) * t124 - t57 * t118;
t134 = t74 * mrSges(7,1) + t78 * mrSges(7,2) - t157;
t132 = pkin(9) * t57;
t131 = t55 * pkin(9);
t113 = t73 * t129;
t120 = pkin(2) * t113 + pkin(9) * t123;
t119 = t130 * pkin(1) + pkin(8) * t112;
t116 = t58 * pkin(2) + t119;
t111 = t76 * t129;
t110 = t80 * t129;
t49 = t55 * pkin(2);
t109 = pkin(9) * t56 - t49;
t51 = t57 * pkin(2);
t108 = pkin(9) * t58 - t51;
t107 = -t31 * pkin(3) + qJ(4) * t32;
t106 = -t35 * pkin(3) + qJ(4) * t36;
t105 = -t53 * pkin(3) + qJ(4) * t54;
t100 = t73 * t110;
t101 = t73 * t111;
t102 = pkin(3) * t100 + qJ(4) * t101 + t120;
t95 = -pkin(1) * t128 + pkin(8) * t114;
t88 = -t56 * pkin(2) + t95;
t87 = t36 * pkin(3) + qJ(4) * t35 + t116;
t84 = t36 * pkin(4) + t87;
t83 = pkin(4) * t100 - pkin(10) * t123 + t102;
t82 = -pkin(3) * t32 - qJ(4) * t31 + t88;
t81 = -pkin(4) * t32 + t82;
t38 = (t110 * t79 + t111 * t75) * t73;
t2 = t12 * t78 - t57 * t74;
t1 = -t12 * t74 - t57 * t78;
t3 = [(-t130 * mrSges(2,1) + t128 * mrSges(2,2) - m(3) * t119 - t58 * mrSges(3,1) - mrSges(3,3) * t112 - m(4) * (t116 + t132) - m(5) * (t87 + t132) - m(6) * t84 - t12 * mrSges(6,1) - m(7) * (pkin(5) * t12 + t84) - t2 * mrSges(7,1) - t1 * mrSges(7,2) - t148 * t36 + t146 * t35 + t143 * t90 - t157 * t57) * g(2) + (t128 * mrSges(2,1) + t130 * mrSges(2,2) - m(3) * t95 + t56 * mrSges(3,1) - mrSges(3,3) * t114 - m(4) * (t88 - t131) - m(5) * (t82 - t131) - m(6) * t81 + t6 * mrSges(6,1) - m(7) * (-pkin(5) * t6 + t81) - t158 * mrSges(7,1) - t159 * mrSges(7,2) + t143 * t150 + t148 * t32 - t146 * t31 + t157 * t55) * g(1) (-m(4) * t120 - m(5) * t102 - m(6) * t83 - t38 * mrSges(6,1) + mrSges(6,3) * t123 - m(7) * (pkin(5) * t38 + t83) - (-t123 * t74 + t38 * t78) * mrSges(7,1) - (-t123 * t78 - t38 * t74) * mrSges(7,2) + t143 * (t100 * t75 - t101 * t79) + (-mrSges(3,1) * t129 - t110 * t148 + t111 * t146 + t141 * t77) * t73) * g(3) + (-m(4) * t109 - m(5) * (t109 + t145) - t149 * (-pkin(4) * t125 + t145 - t49) + t134 * t56 + t151 * t55) * g(2) + (-m(4) * t108 - m(5) * (t108 + t144) - t149 * (-pkin(4) * t124 + t144 - t51) + t134 * t58 + t151 * t57) * g(1) (-m(5) * t105 - t149 * (-t53 * pkin(4) + t105) + t146 * t54 + t148 * t53 - t154) * g(3) + (-m(5) * t107 - t149 * (-t31 * pkin(4) + t107) + t146 * t32 + t148 * t31 - t156) * g(2) + (-m(5) * t106 - t149 * (-t35 * pkin(4) + t106) + t146 * t36 + t148 * t35 - t155) * g(1) (m(5) + t149) * (-g(1) * t35 - g(2) * t31 - g(3) * t53) g(1) * t155 + g(2) * t156 + g(3) * t154, -g(1) * (mrSges(7,1) * t1 - mrSges(7,2) * t2) - g(2) * (-t159 * mrSges(7,1) + t158 * mrSges(7,2)) - g(3) * ((t113 * t78 - t20 * t74) * mrSges(7,1) + (-t113 * t74 - t20 * t78) * mrSges(7,2))];
taug  = t3(:);
