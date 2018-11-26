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

% Quelle: HybrDyn-Toolbox (ehem. IRT-Maple-Toolbox)
% Datum: 2018-11-23 15:07
% Revision: 76f9d5e39f14dc242b53c0d9d3d9db48bd8f37c0
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für mechatronische Systeme, Universität Hannover

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
assert( isreal(m) && all(size(m) == [7 1]), ...
  'S6PRPRRR7_gravloadJ_floatb_twist_slag_vp2: m has to be [7x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [7,3]), ...
  'S6PRPRRR7_gravloadJ_floatb_twist_slag_vp2: mrSges has to be [7x3] (double)');

%% Symbolic Calculation
% From joint_gravload_floatb_twist_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2018-11-23 15:07:06
% EndTime: 2018-11-23 15:07:08
% DurationCPUTime: 1.66s
% Computational Cost: add. (6849->176), mult. (6887->252), div. (0->0), fcn. (6752->30), ass. (0->116)
t178 = m(6) + m(7);
t88 = sin(qJ(6));
t91 = cos(qJ(6));
t183 = -t88 * mrSges(7,1) - t91 * mrSges(7,2) - pkin(11) * t178 + mrSges(5,2) - mrSges(6,3);
t176 = m(5) + t178;
t182 = m(4) + t176;
t181 = m(7) * pkin(5) + t91 * mrSges(7,1) - t88 * mrSges(7,2) + mrSges(6,1);
t144 = -m(7) * pkin(12) + mrSges(6,2) - mrSges(7,3);
t89 = sin(qJ(5));
t92 = cos(qJ(5));
t179 = pkin(4) * t178 - t144 * t89 + t181 * t92 + mrSges(5,1);
t170 = cos(qJ(2));
t169 = sin(qJ(2));
t156 = pkin(6) + qJ(2);
t139 = cos(t156) / 0.2e1;
t157 = pkin(6) - qJ(2);
t149 = cos(t157);
t118 = t149 / 0.2e1 + t139;
t159 = sin(pkin(13));
t161 = cos(pkin(13));
t103 = -t118 * t161 + t159 * t169;
t153 = pkin(7) - pkin(14);
t127 = cos(t153) / 0.2e1;
t152 = pkin(7) + pkin(14);
t143 = cos(t152);
t112 = t127 + t143 / 0.2e1;
t158 = sin(pkin(14));
t147 = sin(t157);
t137 = t147 / 0.2e1;
t146 = sin(t156);
t125 = t137 - t146 / 0.2e1;
t140 = t159 * t170;
t72 = t125 * t161 - t140;
t55 = t103 * t158 + t112 * t72;
t84 = sin(pkin(8));
t168 = t55 * t84;
t105 = t118 * t159 + t161 * t169;
t141 = t161 * t170;
t73 = -t125 * t159 - t141;
t57 = t105 * t158 + t112 * t73;
t167 = t57 * t84;
t136 = t146 / 0.2e1;
t115 = t136 + t137;
t79 = t139 - t149 / 0.2e1;
t66 = t112 * t79 - t115 * t158;
t166 = t66 * t84;
t154 = pkin(8) + qJ(4);
t138 = cos(t154) / 0.2e1;
t155 = pkin(8) - qJ(4);
t148 = cos(t155);
t78 = t138 - t148 / 0.2e1;
t85 = sin(pkin(7));
t165 = t78 * t85;
t87 = cos(pkin(8));
t164 = t85 * t87;
t163 = cos(pkin(6));
t162 = cos(pkin(7));
t160 = sin(pkin(6));
t145 = sin(t155);
t142 = sin(t153);
t135 = sin(t154) / 0.2e1;
t126 = sin(t152) / 0.2e1;
t75 = t126 - t142 / 0.2e1;
t86 = cos(pkin(14));
t56 = -t103 * t86 + t72 * t75;
t70 = t103 * pkin(2);
t132 = t56 * pkin(3) - pkin(10) * t168 - t70;
t58 = -t105 * t86 + t73 * t75;
t71 = t105 * pkin(2);
t131 = t58 * pkin(3) - pkin(10) * t167 - t71;
t67 = t115 * t86 + t79 * t75;
t74 = t115 * pkin(2);
t130 = t67 * pkin(3) - pkin(10) * t166 + t74;
t129 = t161 * t160;
t128 = t160 * t159;
t117 = t148 / 0.2e1 + t138;
t116 = t136 - t147 / 0.2e1;
t114 = t135 + t145 / 0.2e1;
t111 = t126 + t142 / 0.2e1;
t110 = t85 * t114;
t109 = t111 * t160;
t107 = t115 * t85 - t162 * t163;
t106 = -t116 * t159 + t141;
t104 = t116 * t161 + t140;
t102 = -mrSges(3,2) + (m(4) * qJ(3) + mrSges(4,3) + t176 * (pkin(10) * t87 + qJ(3))) * t85;
t101 = -t105 * t85 - t128 * t162;
t100 = -t103 * t85 + t129 * t162;
t99 = t111 * t163 + t112 * t115 + t158 * t79;
t76 = t127 - t143 / 0.2e1;
t59 = t115 * t75 + t163 * t76 - t79 * t86;
t77 = t135 - t145 / 0.2e1;
t93 = cos(qJ(4));
t98 = t107 * t78 + t59 * t93 + t77 * t99;
t97 = t105 * t112 + t106 * t158 - t109 * t159;
t96 = t103 * t112 + t104 * t158 + t109 * t161;
t45 = -t103 * t75 + t104 * t86 - t129 * t76;
t95 = t100 * t78 + t45 * t93 - t77 * t96;
t46 = -t105 * t75 + t106 * t86 + t128 * t76;
t94 = t101 * t78 + t46 * t93 - t77 * t97;
t90 = sin(qJ(4));
t52 = -t164 * t79 - t166;
t44 = -t107 * t87 - t84 * t99;
t39 = -t164 * t73 - t167;
t38 = -t164 * t72 - t168;
t35 = -t101 * t87 + t84 * t97;
t34 = -t100 * t87 + t84 * t96;
t33 = t165 * t79 + t66 * t77 + t67 * t93;
t28 = t107 * t114 - t117 * t99 + t59 * t90;
t26 = t165 * t73 + t57 * t77 + t58 * t93;
t24 = t165 * t72 + t55 * t77 + t56 * t93;
t16 = t101 * t114 + t117 * t97 + t46 * t90;
t13 = t100 * t114 + t117 * t96 + t45 * t90;
t10 = t44 * t89 + t92 * t98;
t4 = t35 * t89 + t92 * t94;
t2 = t34 * t89 + t92 * t95;
t1 = [(-m(2) - m(3) - t182) * g(3) (-t115 * mrSges(3,1) - m(4) * t74 - t67 * mrSges(4,1) - t66 * mrSges(4,2) - m(5) * t130 - t33 * mrSges(5,1) - t52 * mrSges(5,3) - t181 * (t33 * t92 + t52 * t89) + t183 * (t110 * t79 - t117 * t66 + t67 * t90) + t144 * (t33 * t89 - t52 * t92) + t102 * t79 + t178 * (-t33 * pkin(4) - t130)) * g(3) + (t103 * mrSges(3,1) + m(4) * t70 - t56 * mrSges(4,1) - t55 * mrSges(4,2) - m(5) * t132 - t24 * mrSges(5,1) - t38 * mrSges(5,3) + t144 * (t24 * t89 - t38 * t92) - t181 * (t24 * t92 + t38 * t89) + t183 * (t110 * t72 - t117 * t55 + t56 * t90) + t102 * t72 + t178 * (-t24 * pkin(4) - t132)) * g(2) + (t105 * mrSges(3,1) + m(4) * t71 - t58 * mrSges(4,1) - t57 * mrSges(4,2) - m(5) * t131 - t26 * mrSges(5,1) - t39 * mrSges(5,3) + t144 * (t26 * t89 - t39 * t92) - t181 * (t26 * t92 + t39 * t89) + t183 * (t110 * t73 - t117 * t57 + t58 * t90) + t102 * t73 + t178 * (-t26 * pkin(4) - t131)) * g(1) (g(1) * t101 + g(2) * t100 + g(3) * t107) * t182 (t179 * t28 + t183 * t98) * g(3) + (t179 * t13 + t183 * t95) * g(2) + (t179 * t16 + t183 * t94) * g(1) (-t181 * (t44 * t92 - t89 * t98) + t144 * t10) * g(3) + (t144 * t2 - t181 * (t34 * t92 - t89 * t95)) * g(2) + (t144 * t4 - t181 * (t35 * t92 - t89 * t94)) * g(1), -g(1) * ((t16 * t91 - t4 * t88) * mrSges(7,1) + (-t16 * t88 - t4 * t91) * mrSges(7,2)) - g(2) * ((t13 * t91 - t2 * t88) * mrSges(7,1) + (-t13 * t88 - t2 * t91) * mrSges(7,2)) - g(3) * ((-t10 * t88 + t28 * t91) * mrSges(7,1) + (-t10 * t91 - t28 * t88) * mrSges(7,2))];
taug  = t1(:);
