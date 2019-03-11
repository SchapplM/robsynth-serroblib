% Calculate Gravitation load on the joints for
% S6RRRRRP11
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% g [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [12x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,alpha2,alpha3,d1,d2,d3,d4,d5]';
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
% Datum: 2019-03-10 02:58
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function taug = S6RRRRRP11_gravloadJ_floatb_twist_slag_vp2(qJ, g, ...
  pkin, m, mrSges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(3,1),zeros(12,1),zeros(7,1),zeros(7,3)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRRRRP11_gravloadJ_floatb_twist_slag_vp2: qJ has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6RRRRRP11_gravloadJ_floatb_twist_slag_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [12 1]), ...
  'S6RRRRRP11_gravloadJ_floatb_twist_slag_vp2: pkin has to be [12x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6RRRRRP11_gravloadJ_floatb_twist_slag_vp2: m has to be [7x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [7,3]), ...
  'S6RRRRRP11_gravloadJ_floatb_twist_slag_vp2: mrSges has to be [7x3] (double)');

%% Symbolic Calculation
% From gravload_joint_floatb_twist_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-10 02:32:40
% EndTime: 2019-03-10 02:32:45
% DurationCPUTime: 2.42s
% Computational Cost: add. (1584->202), mult. (4283->293), div. (0->0), fcn. (5432->14), ass. (0->100)
t157 = cos(pkin(6));
t168 = cos(qJ(2));
t139 = t157 * t168;
t165 = sin(qJ(2));
t166 = sin(qJ(1));
t169 = cos(qJ(1));
t114 = -t139 * t169 + t166 * t165;
t155 = sin(pkin(6));
t156 = cos(pkin(7));
t129 = t156 * t155;
t154 = sin(pkin(7));
t100 = t114 * t154 - t169 * t129;
t167 = cos(qJ(3));
t128 = t155 * t154;
t188 = t114 * t156 + t169 * t128;
t138 = t157 * t165;
t70 = t138 * t169 + t166 * t168;
t90 = sin(qJ(3));
t40 = -t167 * t70 + t188 * t90;
t89 = sin(qJ(4));
t92 = cos(qJ(4));
t18 = t100 * t89 - t40 * t92;
t37 = t188 * t167 + t70 * t90;
t88 = sin(qJ(5));
t91 = cos(qJ(5));
t192 = t18 * t88 - t37 * t91;
t191 = t18 * t91 + t37 * t88;
t187 = t100 * t92 + t40 * t89;
t186 = mrSges(6,1) + mrSges(7,1);
t181 = mrSges(6,2) + mrSges(7,2);
t85 = pkin(5) * t91 + pkin(4);
t185 = m(6) * pkin(4) + m(7) * t85 + mrSges(5,1);
t184 = m(6) * pkin(12) - m(7) * (-qJ(6) - pkin(12)) - mrSges(5,2) + mrSges(6,3) + mrSges(7,3);
t105 = t139 * t166 + t165 * t169;
t93 = t105 * t154 + t166 * t129;
t170 = m(7) * pkin(5);
t182 = mrSges(4,2) - mrSges(5,3);
t180 = -t154 * mrSges(4,3) + mrSges(3,2);
t178 = t105 * t156 - t166 * t128;
t177 = t168 * t129 + t157 * t154;
t176 = -m(5) - m(6) - m(7);
t175 = -t170 - t186;
t174 = -t181 * t88 + t186 * t91 + t185;
t173 = t184 * t89 + t185 * t92 + mrSges(4,1);
t172 = -m(7) * (pkin(5) * t88 + pkin(11)) + t182;
t164 = t40 * t88;
t71 = -t138 * t166 + t168 * t169;
t42 = t71 * t167 - t178 * t90;
t163 = t42 * t88;
t133 = t155 * t165;
t56 = t167 * t133 + t177 * t90;
t162 = t56 * t88;
t161 = t88 * t92;
t160 = t91 * t92;
t117 = t165 * t128;
t135 = t168 * t155;
t159 = pkin(2) * t135 + pkin(10) * t117;
t134 = t155 * t166;
t158 = t169 * pkin(1) + pkin(9) * t134;
t121 = t165 * t129;
t64 = -t121 * t90 + t135 * t167;
t153 = t64 * pkin(3) + t159;
t146 = pkin(10) * t154;
t145 = t89 * t154;
t144 = t90 * t156;
t143 = t92 * t154;
t136 = t169 * t155;
t141 = -pkin(1) * t166 + pkin(9) * t136;
t137 = t156 * t167;
t22 = t42 * t92 + t89 * t93;
t41 = t178 * t167 + t71 * t90;
t5 = -t22 * t88 + t41 * t91;
t63 = t121 * t167 + t135 * t90;
t131 = pkin(11) * t63 + t153;
t127 = -t114 * pkin(2) + t146 * t70;
t126 = -t105 * pkin(2) + t146 * t71;
t46 = -t114 * t167 - t144 * t70;
t125 = t46 * pkin(3) + t127;
t48 = -t105 * t167 - t144 * t71;
t124 = t48 * pkin(3) + t126;
t45 = -t114 * t90 + t137 * t70;
t112 = t45 * pkin(11) + t125;
t47 = -t105 * t90 + t137 * t71;
t111 = t47 * pkin(11) + t124;
t104 = -t128 * t168 + t156 * t157;
t101 = -t70 * pkin(2) - t100 * pkin(10) + t141;
t98 = t40 * pkin(3) + t101;
t97 = t71 * pkin(2) + t93 * pkin(10) + t158;
t96 = t42 * pkin(3) + t97;
t95 = -pkin(11) * t37 + t98;
t94 = t41 * pkin(11) + t96;
t55 = t133 * t90 - t177 * t167;
t50 = t117 * t89 + t64 * t92;
t36 = t104 * t89 + t56 * t92;
t35 = -t104 * t92 + t56 * t89;
t28 = t145 * t71 + t48 * t92;
t26 = t145 * t70 + t46 * t92;
t21 = t42 * t89 - t92 * t93;
t6 = t22 * t91 + t41 * t88;
t1 = [(-t169 * mrSges(2,1) + t166 * mrSges(2,2) - m(3) * t158 - t71 * mrSges(3,1) + t105 * mrSges(3,2) - mrSges(3,3) * t134 - m(4) * t97 - t42 * mrSges(4,1) - t93 * mrSges(4,3) - m(5) * t94 - t22 * mrSges(5,1) - m(6) * (t22 * pkin(4) + t94) - m(7) * (t22 * t85 + t96) - t186 * t6 - t181 * t5 + t172 * t41 - t184 * t21) * g(2) + (t166 * mrSges(2,1) + t169 * mrSges(2,2) - m(3) * t141 + t70 * mrSges(3,1) - t114 * mrSges(3,2) - mrSges(3,3) * t136 - m(4) * t101 - t40 * mrSges(4,1) + t100 * mrSges(4,3) - m(5) * t95 + t18 * mrSges(5,1) - m(6) * (-pkin(4) * t18 + t95) - m(7) * (-t18 * t85 + t98) + t186 * t191 - t172 * t37 - t181 * t192 - t184 * t187) * g(1) (-mrSges(3,1) * t135 + mrSges(3,2) * t133 - m(4) * t159 - t64 * mrSges(4,1) - mrSges(4,3) * t117 - m(5) * t131 - t50 * mrSges(5,1) - m(6) * (pkin(4) * t50 + t131) - m(7) * (t50 * t85 + t153) + t172 * t63 - t186 * (t50 * t91 + t63 * t88) - t181 * (-t50 * t88 + t63 * t91) - t184 * (-t117 * t92 + t64 * t89)) * g(3) + (mrSges(3,1) * t114 - m(4) * t127 - t46 * mrSges(4,1) - m(5) * t112 - t26 * mrSges(5,1) - m(6) * (t26 * pkin(4) + t112) - m(7) * (t26 * t85 + t125) + t180 * t70 + t172 * t45 - t186 * (t26 * t91 + t45 * t88) - t181 * (-t26 * t88 + t45 * t91) - t184 * (-t143 * t70 + t46 * t89)) * g(2) + (mrSges(3,1) * t105 - m(4) * t126 - t48 * mrSges(4,1) - m(5) * t111 - t28 * mrSges(5,1) - m(6) * (t28 * pkin(4) + t111) - m(7) * (t28 * t85 + t124) + t180 * t71 + t172 * t47 - t186 * (t28 * t91 + t47 * t88) - t181 * (-t28 * t88 + t47 * t91) - t184 * (-t143 * t71 + t48 * t89)) * g(1) (-t162 * t170 + t182 * t56 - t186 * (-t160 * t55 + t162) - t181 * (t161 * t55 + t56 * t91) + t176 * (-t55 * pkin(3) + t56 * pkin(11)) + t173 * t55) * g(3) + (t164 * t170 - t186 * (-t160 * t37 - t164) - t181 * (t161 * t37 - t40 * t91) - t182 * t40 + t176 * (-t37 * pkin(3) - pkin(11) * t40) + t173 * t37) * g(2) + (-t163 * t170 - t181 * (t161 * t41 + t42 * t91) + t182 * t42 + t176 * (-t41 * pkin(3) + t42 * pkin(11)) - t186 * (-t160 * t41 + t163) + t173 * t41) * g(1) (t174 * t35 - t184 * t36) * g(3) + (-t174 * t187 - t18 * t184) * g(2) + (t174 * t21 - t184 * t22) * g(1) (-t181 * (-t36 * t91 - t55 * t88) + t175 * (-t36 * t88 + t55 * t91)) * g(3) + (-t175 * t192 + t181 * t191) * g(2) + (t175 * t5 + t181 * t6) * g(1) (-g(1) * t21 + g(2) * t187 - g(3) * t35) * m(7)];
taug  = t1(:);
