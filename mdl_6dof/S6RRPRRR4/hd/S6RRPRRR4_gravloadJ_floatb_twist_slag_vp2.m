% Calculate Gravitation load on the joints for
% S6RRPRRR4
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% g [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [12x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,alpha2,d1,d2,d4,d5,d6,theta3]';
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
% Datum: 2018-11-23 17:23
% Revision: 76f9d5e39f14dc242b53c0d9d3d9db48bd8f37c0
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für mechatronische Systeme, Universität Hannover

function taug = S6RRPRRR4_gravloadJ_floatb_twist_slag_vp2(qJ, g, ...
  pkin, m, mrSges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(3,1),zeros(12,1),zeros(7,1),zeros(7,3)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRPRRR4_gravloadJ_floatb_twist_slag_vp2: qJ has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6RRPRRR4_gravloadJ_floatb_twist_slag_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [12 1]), ...
  'S6RRPRRR4_gravloadJ_floatb_twist_slag_vp2: pkin has to be [12x1] (double)');
assert( isreal(m) && all(size(m) == [7 1]), ...
  'S6RRPRRR4_gravloadJ_floatb_twist_slag_vp2: m has to be [7x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [7,3]), ...
  'S6RRPRRR4_gravloadJ_floatb_twist_slag_vp2: mrSges has to be [7x3] (double)');

%% Symbolic Calculation
% From joint_gravload_floatb_twist_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2018-11-23 17:22:38
% EndTime: 2018-11-23 17:22:39
% DurationCPUTime: 1.21s
% Computational Cost: add. (1715->165), mult. (1374->211), div. (0->0), fcn. (1289->22), ass. (0->93)
t167 = mrSges(6,2) - mrSges(7,3);
t89 = sin(qJ(6));
t93 = cos(qJ(6));
t165 = mrSges(7,1) * t93 - mrSges(7,2) * t89 + mrSges(6,1);
t166 = m(7) * pkin(5) + t165;
t124 = -m(7) * pkin(11) + t167;
t83 = qJ(2) + pkin(12);
t126 = pkin(6) + t83;
t105 = sin(t126) / 0.2e1;
t127 = pkin(6) - t83;
t112 = sin(t127);
t50 = t105 - t112 / 0.2e1;
t78 = cos(t83);
t92 = sin(qJ(1));
t96 = cos(qJ(1));
t117 = -t50 * t92 + t78 * t96;
t87 = sin(pkin(6));
t144 = t87 * t92;
t90 = sin(qJ(4));
t94 = cos(qJ(4));
t19 = -t117 * t90 + t94 * t144;
t106 = cos(t127) / 0.2e1;
t113 = cos(t126);
t52 = t106 - t113 / 0.2e1;
t88 = cos(pkin(6));
t164 = -t52 * t90 + t88 * t94;
t162 = -m(4) - m(5);
t161 = m(6) + m(7);
t160 = m(5) * pkin(3) + t94 * mrSges(5,1) - t90 * mrSges(5,2) + mrSges(4,1);
t84 = pkin(6) + qJ(2);
t72 = cos(t84) / 0.2e1;
t85 = pkin(6) - qJ(2);
t80 = cos(t85);
t55 = t80 / 0.2e1 + t72;
t86 = qJ(4) + qJ(5);
t81 = sin(t86);
t82 = cos(t86);
t34 = -t52 * t81 + t82 * t88;
t35 = t52 * t82 + t81 * t88;
t159 = -t165 * t34 + t167 * t35;
t114 = -m(5) * pkin(9) + mrSges(4,2) - mrSges(5,3) - mrSges(6,3);
t158 = mrSges(7,1) * t89 + mrSges(7,2) * t93 - t114;
t17 = t117 * t81 - t144 * t82;
t18 = t117 * t82 + t144 * t81;
t157 = t165 * t17 + t167 * t18;
t118 = t96 * t50 + t78 * t92;
t143 = t87 * t96;
t13 = -t118 * t81 - t82 * t143;
t14 = t118 * t82 - t81 * t143;
t156 = -t165 * t13 + t167 * t14;
t155 = -t124 * t81 + t166 * t82 + t160;
t153 = sin(t84) / 0.2e1;
t77 = sin(t85);
t150 = pkin(2) * t77;
t91 = sin(qJ(2));
t141 = t91 * t92;
t140 = t91 * t96;
t67 = pkin(2) * t153;
t138 = t150 / 0.2e1 + t67;
t137 = t90 * t144;
t63 = t90 * t143;
t134 = t161 - t162;
t132 = t13 * pkin(5) + pkin(11) * t14;
t131 = -t17 * pkin(5) + pkin(11) * t18;
t130 = t34 * pkin(5) + pkin(11) * t35;
t75 = sin(t83);
t99 = t113 / 0.2e1 + t106;
t39 = t96 * t75 + t92 * t99;
t95 = cos(qJ(2));
t74 = pkin(2) * t95 + pkin(1);
t62 = t96 * t74;
t73 = pkin(4) * t94 + pkin(3);
t97 = -pkin(10) - pkin(9);
t129 = pkin(4) * t137 + t117 * t73 - t39 * t97 + t62;
t125 = -m(3) * pkin(8) - mrSges(3,3) - mrSges(4,3);
t53 = t55 * pkin(2);
t123 = -pkin(2) * t141 + t96 * t53;
t122 = t19 * pkin(4);
t121 = t164 * pkin(4);
t110 = m(3) * pkin(1) + t95 * mrSges(3,1) + mrSges(2,1);
t109 = -pkin(2) * t140 - t53 * t92;
t108 = -t118 * t90 - t143 * t94;
t102 = t108 * pkin(4);
t54 = t153 - t77 / 0.2e1;
t100 = t54 * mrSges(3,1) + t134 * (-t150 / 0.2e1 + t67 - t87 * (pkin(8) + qJ(3))) + mrSges(2,2);
t51 = t112 / 0.2e1 + t105;
t43 = -t55 * t92 - t140;
t42 = -t55 * t96 + t141;
t36 = t75 * t92 - t96 * t99;
t20 = t117 * t94 + t137;
t2 = t18 * t93 + t39 * t89;
t1 = -t18 * t89 + t39 * t93;
t3 = [(-t43 * mrSges(3,2) - m(4) * t62 - t117 * mrSges(4,1) - m(5) * (pkin(3) * t117 + t62) - t20 * mrSges(5,1) - t19 * mrSges(5,2) - m(6) * t129 - t18 * mrSges(6,1) - m(7) * (pkin(5) * t18 + t129) - t2 * mrSges(7,1) - t1 * mrSges(7,2) + t124 * t17 - t110 * t96 + (t125 * t87 + t100) * t92 + t114 * t39) * g(2) + (-t42 * mrSges(3,2) - t63 * mrSges(5,1) + t124 * t13 + t160 * t118 + (t134 * t74 + t110) * t92 + ((-mrSges(5,2) * t94 + t125) * t87 + t100) * t96 + t166 * t14 + t158 * t36 + t161 * (-pkin(4) * t63 + t118 * t73 - t36 * t97)) * g(1) (-(t153 + t77 / 0.2e1) * mrSges(3,1) - (t72 - t80 / 0.2e1) * mrSges(3,2) + t162 * t138 - t161 * (t51 * t73 - t52 * t97 + t138) - t158 * t52 - t155 * t51) * g(3) + (t42 * mrSges(3,1) - (-t54 * t96 - t92 * t95) * mrSges(3,2) + t162 * t123 - t161 * (-t118 * t97 - t36 * t73 + t123) - t158 * t118 + t155 * t36) * g(2) + (-t43 * mrSges(3,1) - (t54 * t92 - t95 * t96) * mrSges(3,2) + t162 * t109 - t161 * (-t117 * t97 - t39 * t73 + t109) - t158 * t117 + t155 * t39) * g(1) (-t88 * g(3) + (-g(1) * t92 + g(2) * t96) * t87) * t134 (-t164 * mrSges(5,1) - (-t52 * t94 - t88 * t90) * mrSges(5,2) - m(6) * t121 - m(7) * (t121 + t130) + t159) * g(3) + (-t108 * mrSges(5,1) - (-t118 * t94 + t63) * mrSges(5,2) - m(6) * t102 - m(7) * (t102 + t132) + t156) * g(2) + (-t19 * mrSges(5,1) + t20 * mrSges(5,2) - m(6) * t122 - m(7) * (t122 + t131) + t157) * g(1) (-m(7) * t130 + t159) * g(3) + (-m(7) * t132 + t156) * g(2) + (-m(7) * t131 + t157) * g(1), -g(1) * (mrSges(7,1) * t1 - mrSges(7,2) * t2) - g(2) * ((-t14 * t89 + t36 * t93) * mrSges(7,1) + (-t14 * t93 - t36 * t89) * mrSges(7,2)) - g(3) * ((-t35 * t89 - t51 * t93) * mrSges(7,1) + (-t35 * t93 + t51 * t89) * mrSges(7,2))];
taug  = t3(:);
