% Calculate Gravitation load on the joints for
% S6RRRPRR9
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% g [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [13x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,alpha2,alpha3,d1,d2,d3,d5,d6,theta4]';
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
% Datum: 2019-03-09 19:14
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function taug = S6RRRPRR9_gravloadJ_floatb_twist_slag_vp2(qJ, g, ...
  pkin, m, mrSges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(3,1),zeros(13,1),zeros(7,1),zeros(7,3)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRRPRR9_gravloadJ_floatb_twist_slag_vp2: qJ has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6RRRPRR9_gravloadJ_floatb_twist_slag_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [13 1]), ...
  'S6RRRPRR9_gravloadJ_floatb_twist_slag_vp2: pkin has to be [13x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6RRRPRR9_gravloadJ_floatb_twist_slag_vp2: m has to be [7x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [7,3]), ...
  'S6RRRPRR9_gravloadJ_floatb_twist_slag_vp2: mrSges has to be [7x3] (double)');

%% Symbolic Calculation
% From gravload_joint_floatb_twist_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 18:56:42
% EndTime: 2019-03-09 18:56:48
% DurationCPUTime: 2.38s
% Computational Cost: add. (1581->211), mult. (4187->310), div. (0->0), fcn. (5341->16), ass. (0->111)
t185 = m(6) + m(7);
t191 = t185 * pkin(11);
t155 = -mrSges(4,3) - mrSges(5,3);
t93 = sin(qJ(6));
t98 = cos(qJ(6));
t182 = m(7) * pkin(5) + t98 * mrSges(7,1) - t93 * mrSges(7,2) + mrSges(6,1);
t181 = -m(7) * pkin(12) + mrSges(6,2) - mrSges(7,3);
t190 = m(4) * pkin(10);
t126 = t93 * mrSges(7,1) + t98 * mrSges(7,2);
t183 = mrSges(6,3) - mrSges(5,2);
t177 = t126 + t183;
t100 = cos(qJ(3));
t101 = cos(qJ(2));
t143 = t100 * t101;
t95 = sin(qJ(3));
t96 = sin(qJ(2));
t156 = t95 * t96;
t92 = cos(pkin(7));
t189 = t92 * t143 - t156;
t150 = t100 * t92;
t102 = cos(qJ(1));
t145 = cos(pkin(6));
t132 = t102 * t145;
t97 = sin(qJ(1));
t71 = t97 * t101 + t132 * t96;
t167 = t71 * t95;
t70 = -t101 * t132 + t96 * t97;
t188 = t70 * t150 + t167;
t91 = sin(pkin(6));
t146 = t102 * t91;
t144 = cos(pkin(13));
t89 = sin(pkin(13));
t115 = t100 * t89 + t144 * t95;
t90 = sin(pkin(7));
t64 = t115 * t90;
t66 = t115 * t92;
t74 = -t100 * t144 + t95 * t89;
t106 = t146 * t64 + t70 * t66 + t71 * t74;
t168 = t70 * t90;
t50 = t146 * t92 - t168;
t94 = sin(qJ(5));
t99 = cos(qJ(5));
t187 = t106 * t99 + t50 * t94;
t186 = t106 * t94 - t50 * t99;
t184 = mrSges(4,2) * t95;
t180 = t181 * t94 - t182 * t99 - mrSges(5,1);
t30 = (t101 * t66 - t74 * t96) * t91 + t145 * t64;
t178 = t177 + t191;
t175 = m(4) * pkin(2) + mrSges(4,1) * t100 + mrSges(3,1) - t184;
t157 = t92 * t95;
t174 = mrSges(4,1) * t157 + mrSges(4,2) * t150 + mrSges(3,2) + (-t190 + t155) * t90;
t173 = pkin(3) * t90;
t169 = pkin(3) * t100;
t135 = t97 * t145;
t72 = -t101 * t135 - t102 * t96;
t166 = t72 * t90;
t73 = t101 * t102 - t135 * t96;
t165 = t73 * t95;
t162 = t90 * t94;
t161 = t90 * t96;
t160 = t90 * t99;
t159 = t91 * t96;
t158 = t91 * t97;
t154 = pkin(10) + qJ(4);
t68 = pkin(3) * t157 - t154 * t90;
t87 = pkin(2) + t169;
t153 = -t71 * t68 - t70 * t87;
t152 = -t73 * t68 + t72 * t87;
t151 = t102 * pkin(1) + pkin(9) * t158;
t149 = t100 * t96;
t148 = t101 * t91;
t147 = t101 * t95;
t142 = t90 * t159;
t141 = t90 * t158;
t140 = t90 * t146;
t136 = -t97 * pkin(1) + pkin(9) * t146;
t133 = t100 * t145;
t131 = -t92 * t190 - mrSges(3,3);
t129 = t87 * t148 - t159 * t68;
t67 = t154 * t92 + t173 * t95;
t128 = t67 * t158 + t72 * t68 + t73 * t87 + t151;
t124 = t141 * t169 + (t150 * t72 - t165) * pkin(3);
t121 = t183 + t191;
t26 = t158 * t64 + t66 * t72 - t73 * t74;
t120 = t26 * pkin(4) + t128;
t119 = t72 * t92 + t141;
t118 = t189 * pkin(3) * t91 + t133 * t173;
t116 = t67 * t146 + t70 * t68 - t71 * t87 + t136;
t113 = -t100 * t71 + t95 * t140 + t157 * t70;
t63 = t74 * t90;
t65 = t74 * t92;
t107 = t115 * t71 - t146 * t63 - t70 * t65;
t105 = (-t100 * t140 - t188) * pkin(3);
t69 = t145 * t92 - t148 * t90;
t52 = t158 * t92 - t166;
t43 = (-t101 * t74 - t66 * t96) * t91;
t42 = -t115 * t148 + t159 * t65;
t40 = t100 * t73 + t119 * t95;
t39 = t100 * t119 - t165;
t36 = -t66 * t73 - t72 * t74;
t35 = -t115 * t72 + t65 * t73;
t34 = -t66 * t71 + t70 * t74;
t33 = t115 * t70 + t65 * t71;
t29 = -t145 * t63 + (-t101 * t65 - t115 * t96) * t91;
t25 = -t115 * t73 - t158 * t63 - t65 * t72;
t14 = t30 * t99 + t69 * t94;
t8 = t26 * t99 + t52 * t94;
t7 = t26 * t94 - t52 * t99;
t2 = -t25 * t93 + t8 * t98;
t1 = -t25 * t98 - t8 * t93;
t3 = [(-t102 * mrSges(2,1) - m(3) * t151 - t73 * mrSges(3,1) - t72 * mrSges(3,2) - m(4) * (pkin(2) * t73 - pkin(10) * t166 + t151) - t40 * mrSges(4,1) - t39 * mrSges(4,2) - m(5) * t128 - t26 * mrSges(5,1) - m(6) * t120 - t8 * mrSges(6,1) - m(7) * (pkin(5) * t8 + t120) - t2 * mrSges(7,1) - t1 * mrSges(7,2) + (t131 * t91 + mrSges(2,2)) * t97 + t181 * t7 + t155 * t52 + t121 * t25) * g(2) + (t97 * mrSges(2,1) - m(3) * t136 + t71 * mrSges(3,1) - t70 * mrSges(3,2) - m(4) * (-t71 * pkin(2) - pkin(10) * t168 + t136) - t113 * mrSges(4,1) - t188 * mrSges(4,2) - m(5) * t116 - t106 * mrSges(5,1) + t155 * t50 + t181 * t186 - t182 * t187 + (t121 + t126) * t107 + (mrSges(2,2) + (-mrSges(4,2) * t100 * t90 + t131) * t91) * t102 + t185 * (-pkin(4) * t106 - t116)) * g(1) (-m(5) * t153 - t34 * mrSges(5,1) - t185 * (t34 * pkin(4) - pkin(11) * t33 + t153) - t182 * (t162 * t71 + t34 * t99) + t177 * t33 + t181 * (-t160 * t71 + t34 * t94) + t175 * t70 + t174 * t71) * g(2) + (-m(5) * t152 - t36 * mrSges(5,1) - t185 * (t36 * pkin(4) - pkin(11) * t35 + t152) - t182 * (t162 * t73 + t36 * t99) + t177 * t35 + t181 * (-t160 * t73 + t36 * t94) - t175 * t72 + t174 * t73) * g(1) + ((-mrSges(3,1) * t101 + mrSges(3,2) * t96 - m(4) * (pkin(2) * t101 + pkin(10) * t161) - (-t156 * t92 + t143) * mrSges(4,1) - (-t149 * t92 - t147) * mrSges(4,2) - mrSges(4,3) * t161) * t91 - m(5) * t129 - t43 * mrSges(5,1) - mrSges(5,3) * t142 - t185 * (t43 * pkin(4) - pkin(11) * t42 + t129) - t182 * (t142 * t94 + t43 * t99) + t177 * t42 + t181 * (-t142 * t99 + t43 * t94)) * g(3) (-(mrSges(4,1) * t133 - t145 * t184) * t90 - (t189 * mrSges(4,1) + (-t147 * t92 - t149) * mrSges(4,2)) * t91 - m(5) * t118 - t185 * (t29 * pkin(4) + t118) + t180 * t29 - t178 * t30) * g(3) + (-(-t167 + (-t70 * t92 - t140) * t100) * mrSges(4,1) - t113 * mrSges(4,2) - m(5) * t105 - t185 * (-pkin(4) * t107 + t105) - t180 * t107 + t178 * t106) * g(2) + (-m(5) * t124 - mrSges(4,1) * t39 + mrSges(4,2) * t40 - t185 * (t25 * pkin(4) + t124) + t180 * t25 - t178 * t26) * g(1) (m(5) + t185) * (-g(1) * t52 + g(2) * t50 - g(3) * t69) (t181 * t14 - t182 * (-t30 * t94 + t69 * t99)) * g(3) + (-t181 * t187 - t182 * t186) * g(2) + (t181 * t8 + t182 * t7) * g(1), -g(1) * (mrSges(7,1) * t1 - mrSges(7,2) * t2) - g(2) * ((t107 * t98 + t187 * t93) * mrSges(7,1) + (-t107 * t93 + t187 * t98) * mrSges(7,2)) - g(3) * ((-t14 * t93 - t29 * t98) * mrSges(7,1) + (-t14 * t98 + t29 * t93) * mrSges(7,2))];
taug  = t3(:);
