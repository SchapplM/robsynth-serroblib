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

% Quelle: HybrDyn-Toolbox (ehem. IRT-Maple-Toolbox)
% Datum: 2018-11-23 17:32
% Revision: 76f9d5e39f14dc242b53c0d9d3d9db48bd8f37c0
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für mechatronische Systeme, Universität Hannover

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
assert( isreal(m) && all(size(m) == [7 1]), ...
  'S6RRRPPP1_gravloadJ_floatb_twist_slag_vp2: m has to be [7x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [7,3]), ...
  'S6RRRPPP1_gravloadJ_floatb_twist_slag_vp2: mrSges has to be [7x3] (double)');

%% Symbolic Calculation
% From joint_gravload_floatb_twist_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2018-11-23 17:31:55
% EndTime: 2018-11-23 17:31:57
% DurationCPUTime: 1.73s
% Computational Cost: add. (1119->178), mult. (1878->230), div. (0->0), fcn. (1973->14), ass. (0->101)
t81 = sin(pkin(6));
t129 = qJ(4) * t81;
t84 = sin(qJ(3));
t119 = t84 * t129;
t83 = cos(pkin(6));
t128 = qJ(4) * t83;
t85 = sin(qJ(2));
t87 = cos(qJ(3));
t88 = cos(qJ(2));
t166 = t88 * t128 + t85 * (-pkin(3) * t87 - pkin(2) - t119);
t162 = -mrSges(5,3) - mrSges(6,1);
t165 = t162 - mrSges(7,1);
t106 = t87 * mrSges(4,1) - t84 * mrSges(4,2);
t156 = m(7) * pkin(5) - t165;
t143 = t84 * t85;
t48 = t143 * t81 - t83 * t88;
t164 = t48 * t156 + (mrSges(3,2) - mrSges(4,3)) * t88 + (m(4) * pkin(2) + mrSges(3,1) + t106) * t85;
t154 = m(6) + m(7);
t163 = mrSges(2,2) - mrSges(3,3);
t75 = t88 * pkin(2);
t117 = -pkin(1) - t75;
t86 = sin(qJ(1));
t137 = t86 * t88;
t89 = cos(qJ(1));
t51 = t137 * t84 + t87 * t89;
t134 = t89 * t84;
t136 = t87 * t88;
t52 = t136 * t86 - t134;
t76 = t89 * pkin(8);
t161 = -t52 * pkin(3) - t129 * t51 + t76 + ((-pkin(9) - t128) * t85 + t117) * t86;
t135 = t88 * t89;
t65 = pkin(9) * t135;
t159 = t166 * t89 + t65;
t63 = pkin(9) * t137;
t158 = t166 * t86 + t63;
t108 = mrSges(3,1) * t88 - mrSges(3,2) * t85;
t157 = t85 * mrSges(4,3) + t108;
t149 = mrSges(5,2) - mrSges(6,3) - mrSges(7,2);
t147 = m(7) * qJ(6) + mrSges(5,1) - mrSges(6,2) + mrSges(7,3);
t139 = t85 * t87;
t78 = pkin(6) + pkin(10);
t68 = sin(t78) / 0.2e1;
t79 = pkin(6) - pkin(10);
t71 = sin(t79);
t56 = t68 - t71 / 0.2e1;
t82 = cos(pkin(10));
t150 = t82 * t139 - t56 * t143;
t148 = mrSges(4,2) + (-m(7) * (pkin(5) + qJ(4)) + t165) * t81;
t73 = t85 * pkin(9);
t144 = t81 * t87;
t142 = t84 * t88;
t140 = t85 * t86;
t138 = t85 * t89;
t131 = t75 + t73;
t130 = t89 * pkin(1) + t86 * pkin(8);
t69 = cos(t79) / 0.2e1;
t72 = cos(t78);
t127 = t69 - t72 / 0.2e1;
t126 = t69 + t72 / 0.2e1;
t125 = pkin(3) * t143;
t124 = t81 * t139;
t121 = t52 * t129;
t54 = t135 * t87 + t86 * t84;
t120 = t54 * t129;
t62 = t85 * t128;
t116 = pkin(2) * t135 + pkin(9) * t138 + t130;
t115 = t85 * t127;
t114 = t88 * t127;
t113 = t89 * t127;
t80 = sin(pkin(10));
t16 = t126 * t52 - t51 * t80;
t17 = -t51 * t82 - t52 * t56;
t44 = t51 * pkin(3);
t112 = t17 * pkin(4) + qJ(5) * t16 - t44;
t53 = -t134 * t88 + t86 * t87;
t18 = t126 * t54 + t53 * t80;
t19 = t53 * t82 - t54 * t56;
t46 = t53 * pkin(3);
t111 = t19 * pkin(4) + qJ(5) * t18 + t46;
t32 = -t126 * t139 + t143 * t80;
t33 = (t56 * t87 + t82 * t84) * t85;
t58 = qJ(4) * t124;
t110 = -t33 * pkin(4) - qJ(5) * t32 + t58;
t109 = pkin(3) * t136 + t88 * t119 + t131 + t62;
t26 = t140 * t83 + t51 * t81;
t5 = -t86 * t115 + t51 * t56 - t52 * t82;
t99 = t126 * t84 + t87 * t80;
t93 = t54 * pkin(3) - t129 * t53 + t89 * t62 + t116;
t55 = t68 + t71 / 0.2e1;
t4 = -t126 * t51 + t140 * t55 - t52 * t80;
t21 = t88 * t55 + t85 * t99;
t27 = t138 * t83 - t53 * t81;
t23 = t136 * t82 - t142 * t56 + t115;
t22 = -t85 * t55 + t88 * t99;
t13 = -t113 * t88 + t150 * t89;
t12 = t21 * t89;
t11 = (-t114 + t150) * t86;
t10 = t21 * t86;
t7 = t113 * t85 + t53 * t56 + t54 * t82;
t6 = -t126 * t53 - t138 * t55 + t54 * t80;
t1 = [(-m(3) * t130 - m(4) * t116 - m(5) * t93 - t54 * mrSges(4,1) - t53 * mrSges(4,2) - mrSges(4,3) * t138 - t154 * (t7 * pkin(4) + qJ(5) * t6 + t93) + (-mrSges(2,1) - t108) * t89 + t163 * t86 - t147 * t7 + t149 * t6 - t156 * t27) * g(2) + (t52 * mrSges(4,1) - t51 * mrSges(4,2) - t161 * m(5) + t163 * t89 + (-m(3) - m(4)) * t76 - t154 * (t5 * pkin(4) + qJ(5) * t4 + t161) + (mrSges(2,1) + m(3) * pkin(1) - m(4) * (t117 - t73) + t157) * t86 - t147 * t5 + t149 * t4 + t156 * t26) * g(1) (-m(4) * t131 - m(5) * t109 - t106 * t88 - t154 * (t23 * pkin(4) + qJ(5) * t22 + t109) - t156 * (t142 * t81 + t83 * t85) - t147 * t23 + t149 * t22 - t157) * g(3) + (-m(4) * t63 - t158 * m(5) - t154 * (-t11 * pkin(4) - qJ(5) * t10 + t158) + t147 * t11 - t149 * t10 + t164 * t86) * g(2) + (-m(4) * t65 - t159 * m(5) - t154 * (-t13 * pkin(4) - t12 * qJ(5) + t159) + t147 * t13 - t149 * t12 + t164 * t89) * g(1) (-m(5) * (t58 - t125) - m(6) * (t110 - t125) - m(7) * t110 + (mrSges(4,1) * t84 + mrSges(4,2) * t87 - m(7) * (-pkin(3) * t84 + pkin(5) * t144) - mrSges(7,1) * t144) * t85 + t162 * t124 + t147 * t33 - t149 * t32) * g(3) + (mrSges(4,1) * t51 - m(5) * (-t44 + t121) - m(6) * (t112 + t121) - m(7) * t112 - t147 * t17 + t149 * t16 + t148 * t52) * g(2) + (-mrSges(4,1) * t53 - m(5) * (t46 + t120) - m(6) * (t111 + t120) - m(7) * t111 - t147 * t19 + t149 * t18 + t148 * t54) * g(1) (m(5) + t154) * (-g(1) * t27 - g(2) * t26 - g(3) * t48) t154 * (-g(1) * t6 + g(2) * t4 - g(3) * t21) (-g(1) * t7 + g(2) * t5 - g(3) * (-t114 + (-t56 * t84 + t82 * t87) * t85)) * m(7)];
taug  = t1(:);
