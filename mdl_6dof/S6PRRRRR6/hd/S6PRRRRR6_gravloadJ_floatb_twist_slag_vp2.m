% Calculate Gravitation load on the joints for
% S6PRRRRR6
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% g [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [14x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,alpha2,alpha3,alpha4,d2,d3,d4,d5,d6,theta1]';
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
% Datum: 2018-11-23 15:36
% Revision: 76f9d5e39f14dc242b53c0d9d3d9db48bd8f37c0
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für mechatronische Systeme, Universität Hannover

function taug = S6PRRRRR6_gravloadJ_floatb_twist_slag_vp2(qJ, g, ...
  pkin, m, mrSges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(3,1),zeros(14,1),zeros(7,1),zeros(7,3)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6PRRRRR6_gravloadJ_floatb_twist_slag_vp2: qJ has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6PRRRRR6_gravloadJ_floatb_twist_slag_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [14 1]), ...
  'S6PRRRRR6_gravloadJ_floatb_twist_slag_vp2: pkin has to be [14x1] (double)');
assert( isreal(m) && all(size(m) == [7 1]), ...
  'S6PRRRRR6_gravloadJ_floatb_twist_slag_vp2: m has to be [7x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [7,3]), ...
  'S6PRRRRR6_gravloadJ_floatb_twist_slag_vp2: mrSges has to be [7x3] (double)');

%% Symbolic Calculation
% From joint_gravload_floatb_twist_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2018-11-23 15:36:13
% EndTime: 2018-11-23 15:36:15
% DurationCPUTime: 2.31s
% Computational Cost: add. (9277->211), mult. (9224->306), div. (0->0), fcn. (9003->30), ass. (0->117)
t113 = sin(qJ(6));
t116 = cos(qJ(6));
t203 = m(6) + m(7);
t207 = -mrSges(7,1) * t113 - mrSges(7,2) * t116 - pkin(12) * t203 + mrSges(5,2) - mrSges(6,3);
t180 = m(5) + t203;
t206 = m(7) * pkin(5) + mrSges(7,1) * t116 - mrSges(7,2) * t113 + mrSges(6,1);
t165 = -m(7) * pkin(13) + mrSges(6,2) - mrSges(7,3);
t114 = sin(qJ(5));
t117 = cos(qJ(5));
t204 = pkin(4) * t203 - t114 * t165 + t117 * t206 + mrSges(5,1);
t197 = sin(qJ(2));
t196 = sin(qJ(3));
t110 = sin(pkin(8));
t178 = pkin(6) + qJ(2);
t164 = cos(t178) / 0.2e1;
t179 = pkin(6) - qJ(2);
t171 = cos(t179);
t142 = t171 / 0.2e1 + t164;
t185 = sin(pkin(14));
t187 = cos(pkin(14));
t130 = -t142 * t187 + t185 * t197;
t176 = pkin(7) + qJ(3);
t163 = cos(t176) / 0.2e1;
t177 = pkin(7) - qJ(3);
t170 = cos(t177);
t141 = t170 / 0.2e1 + t163;
t161 = sin(t178) / 0.2e1;
t168 = sin(t179);
t101 = t161 - t168 / 0.2e1;
t120 = cos(qJ(2));
t143 = t101 * t187 + t120 * t185;
t77 = t130 * t196 - t141 * t143;
t192 = t110 * t77;
t131 = t142 * t185 + t187 * t197;
t144 = -t101 * t185 + t120 * t187;
t79 = t131 * t196 - t141 * t144;
t191 = t110 * t79;
t104 = t164 - t171 / 0.2e1;
t139 = t161 + t168 / 0.2e1;
t90 = t104 * t141 - t139 * t196;
t190 = t110 * t90;
t189 = cos(pkin(6));
t188 = cos(pkin(7));
t186 = sin(pkin(6));
t174 = pkin(8) + qJ(4);
t162 = cos(t174) / 0.2e1;
t175 = pkin(8) - qJ(4);
t169 = cos(t175);
t102 = t162 - t169 / 0.2e1;
t111 = sin(pkin(7));
t184 = t102 * t111;
t183 = t110 * t114;
t182 = t110 * t117;
t112 = cos(pkin(8));
t181 = t111 * t112;
t167 = sin(t177);
t166 = sin(t175);
t160 = sin(t176) / 0.2e1;
t159 = sin(t174) / 0.2e1;
t158 = t187 * t186;
t157 = t186 * t185;
t100 = t160 - t167 / 0.2e1;
t119 = cos(qJ(3));
t78 = -t100 * t143 - t119 * t130;
t94 = t130 * pkin(2);
t154 = pkin(3) * t78 - pkin(11) * t192 - t94;
t80 = -t100 * t144 - t119 * t131;
t95 = t131 * pkin(2);
t153 = pkin(3) * t80 - pkin(11) * t191 - t95;
t91 = t100 * t104 + t119 * t139;
t98 = t139 * pkin(2);
t152 = pkin(3) * t91 - pkin(11) * t190 + t98;
t145 = -mrSges(4,2) + (pkin(11) * t180 + mrSges(5,3)) * t110;
t140 = t169 / 0.2e1 + t162;
t138 = t160 + t167 / 0.2e1;
t137 = t159 + t166 / 0.2e1;
t135 = t111 * t137;
t134 = t138 * t186;
t132 = t111 * t139 - t188 * t189;
t103 = t163 - t170 / 0.2e1;
t83 = -t100 * t139 + t103 * t189 + t104 * t119;
t129 = -mrSges(3,2) + (m(4) * pkin(10) + mrSges(4,3) + t180 * (pkin(11) * t112 + pkin(10))) * t111;
t128 = -t111 * t131 - t157 * t188;
t127 = -t111 * t130 + t158 * t188;
t126 = t104 * t196 + t138 * t189 + t139 * t141;
t66 = -t100 * t131 - t103 * t157 + t119 * t144;
t64 = -t100 * t130 + t103 * t158 + t119 * t143;
t125 = t131 * t141 - t134 * t185 + t144 * t196;
t124 = t130 * t141 + t134 * t187 + t143 * t196;
t118 = cos(qJ(4));
t99 = t159 - t166 / 0.2e1;
t123 = t102 * t132 - t118 * t83 + t126 * t99;
t122 = t102 * t127 + t118 * t64 - t124 * t99;
t121 = t102 * t128 + t118 * t66 - t125 * t99;
t115 = sin(qJ(4));
t81 = t126 * pkin(3);
t74 = -t104 * t181 - t190;
t63 = t125 * pkin(3);
t62 = t124 * pkin(3);
t61 = -t110 * t126 - t112 * t132;
t54 = t144 * t181 - t191;
t53 = t143 * t181 - t192;
t50 = t110 * t125 - t112 * t128;
t49 = t110 * t124 - t112 * t127;
t48 = t104 * t184 + t118 * t91 + t90 * t99;
t45 = t118 * t126 + t83 * t99;
t40 = -t115 * t83 - t126 * t140 + t132 * t137;
t38 = -t118 * t125 - t66 * t99;
t36 = -t118 * t124 - t64 * t99;
t34 = t118 * t80 - t144 * t184 + t79 * t99;
t32 = t118 * t78 - t143 * t184 + t77 * t99;
t20 = t115 * t66 + t125 * t140 + t128 * t137;
t17 = t115 * t64 + t124 * t140 + t127 * t137;
t14 = t114 * t61 + t117 * t123;
t4 = t114 * t50 + t117 * t121;
t2 = t114 * t49 + t117 * t122;
t1 = [(-m(2) - m(3) - m(4) - t180) * g(3) (-t139 * mrSges(3,1) - m(4) * t98 - t91 * mrSges(4,1) - t90 * mrSges(4,2) - m(5) * t152 - t48 * mrSges(5,1) - t74 * mrSges(5,3) - t206 * (t114 * t74 + t117 * t48) + t207 * (t104 * t135 + t115 * t91 - t140 * t90) + t165 * (t114 * t48 - t117 * t74) + t129 * t104 + t203 * (-pkin(4) * t48 - t152)) * g(3) + (t130 * mrSges(3,1) + m(4) * t94 - t78 * mrSges(4,1) - t77 * mrSges(4,2) - m(5) * t154 - t32 * mrSges(5,1) - t53 * mrSges(5,3) + t165 * (t114 * t32 - t117 * t53) - t206 * (t114 * t53 + t117 * t32) + t207 * (t115 * t78 - t135 * t143 - t140 * t77) - t129 * t143 + t203 * (-pkin(4) * t32 - t154)) * g(2) + (t131 * mrSges(3,1) + m(4) * t95 - t80 * mrSges(4,1) - t79 * mrSges(4,2) - m(5) * t153 - t34 * mrSges(5,1) - t54 * mrSges(5,3) + t165 * (t114 * t34 - t117 * t54) - t206 * (t114 * t54 + t117 * t34) + t207 * (t115 * t80 - t135 * t144 - t140 * t79) - t129 * t144 + t203 * (-pkin(4) * t34 - t153)) * g(1) (-t126 * mrSges(4,1) - m(5) * t81 - t45 * mrSges(5,1) + t145 * t83 - t206 * (t117 * t45 - t183 * t83) + t207 * (t115 * t126 - t140 * t83) + t165 * (t114 * t45 + t182 * t83) - t203 * (pkin(4) * t45 + t81)) * g(3) + (t124 * mrSges(4,1) + m(5) * t62 - t36 * mrSges(5,1) + t165 * (t114 * t36 - t182 * t64) - t145 * t64 - t206 * (t117 * t36 + t183 * t64) + t207 * (-t115 * t124 + t140 * t64) - t203 * (pkin(4) * t36 - t62)) * g(2) + (t125 * mrSges(4,1) + m(5) * t63 - t38 * mrSges(5,1) - t145 * t66 - t206 * (t117 * t38 + t183 * t66) + t207 * (-t115 * t125 + t140 * t66) + t165 * (t114 * t38 - t182 * t66) - t203 * (pkin(4) * t38 - t63)) * g(1) (t123 * t207 + t204 * t40) * g(3) + (t122 * t207 + t17 * t204) * g(2) + (t121 * t207 + t20 * t204) * g(1) (t165 * t14 - t206 * (-t114 * t123 + t117 * t61)) * g(3) + (t165 * t2 - t206 * (-t114 * t122 + t117 * t49)) * g(2) + (t165 * t4 - t206 * (-t114 * t121 + t117 * t50)) * g(1), -g(1) * ((-t113 * t4 + t116 * t20) * mrSges(7,1) + (-t113 * t20 - t116 * t4) * mrSges(7,2)) - g(2) * ((-t113 * t2 + t116 * t17) * mrSges(7,1) + (-t113 * t17 - t116 * t2) * mrSges(7,2)) - g(3) * ((-t113 * t14 + t116 * t40) * mrSges(7,1) + (-t113 * t40 - t116 * t14) * mrSges(7,2))];
taug  = t1(:);
