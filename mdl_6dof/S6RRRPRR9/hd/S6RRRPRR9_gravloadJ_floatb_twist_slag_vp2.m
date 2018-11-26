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

% Quelle: HybrDyn-Toolbox (ehem. IRT-Maple-Toolbox)
% Datum: 2018-11-23 17:58
% Revision: 76f9d5e39f14dc242b53c0d9d3d9db48bd8f37c0
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für mechatronische Systeme, Universität Hannover

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
assert( isreal(m) && all(size(m) == [7 1]), ...
  'S6RRRPRR9_gravloadJ_floatb_twist_slag_vp2: m has to be [7x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [7,3]), ...
  'S6RRRPRR9_gravloadJ_floatb_twist_slag_vp2: mrSges has to be [7x3] (double)');

%% Symbolic Calculation
% From joint_gravload_floatb_twist_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2018-11-23 17:57:01
% EndTime: 2018-11-23 17:57:03
% DurationCPUTime: 2.01s
% Computational Cost: add. (4381->222), mult. (3903->296), div. (0->0), fcn. (3741->28), ass. (0->119)
t126 = sin(qJ(6));
t131 = cos(qJ(6));
t212 = m(6) + m(7);
t160 = -t212 * pkin(11) + mrSges(5,2) - mrSges(6,3);
t215 = -t126 * mrSges(7,1) - t131 * mrSges(7,2) + t160;
t210 = m(7) * pkin(5) + t131 * mrSges(7,1) - t126 * mrSges(7,2) + mrSges(6,1);
t173 = -m(7) * pkin(12) + mrSges(6,2) - mrSges(7,3);
t127 = sin(qJ(5));
t132 = cos(qJ(5));
t119 = qJ(3) + pkin(13);
t116 = cos(t119);
t123 = sin(pkin(6));
t135 = cos(qJ(1));
t189 = t123 * t135;
t129 = sin(qJ(2));
t130 = sin(qJ(1));
t186 = pkin(6) + qJ(2);
t172 = cos(t186) / 0.2e1;
t187 = pkin(6) - qJ(2);
t177 = cos(t187);
t158 = t177 / 0.2e1 + t172;
t73 = t129 * t130 - t135 * t158;
t134 = cos(qJ(2));
t171 = sin(t186) / 0.2e1;
t176 = sin(t187);
t96 = t171 - t176 / 0.2e1;
t74 = t130 * t134 + t135 * t96;
t174 = pkin(7) + t119;
t154 = sin(t174) / 0.2e1;
t175 = pkin(7) - t119;
t165 = sin(t175);
t89 = t154 - t165 / 0.2e1;
t155 = cos(t175) / 0.2e1;
t166 = cos(t174);
t90 = t155 - t166 / 0.2e1;
t23 = -t74 * t116 + t90 * t189 + t73 * t89;
t124 = cos(pkin(7));
t122 = sin(pkin(7));
t196 = t73 * t122;
t58 = t124 * t189 - t196;
t214 = t127 * t58 + t132 * t23;
t213 = t127 * t23 - t58 * t132;
t209 = -t173 * t127 + t210 * t132 + mrSges(5,1);
t207 = m(4) * pkin(10);
t120 = pkin(7) + qJ(3);
t117 = cos(t120);
t206 = t117 / 0.2e1;
t121 = pkin(7) - qJ(3);
t205 = cos(t121);
t204 = sin(t121);
t203 = pkin(3) * t117;
t202 = -mrSges(4,3) - mrSges(5,3);
t133 = cos(qJ(3));
t114 = pkin(3) * t133 + pkin(2);
t181 = sin(t120) / 0.2e1;
t102 = pkin(3) * t181;
t125 = pkin(10) + qJ(4);
t180 = pkin(3) * t204;
t79 = -t180 / 0.2e1 + t102 - t122 * t125;
t201 = -t73 * t114 - t74 * t79;
t76 = -t135 * t129 - t130 * t158;
t78 = t130 * t96 - t134 * t135;
t200 = t76 * t114 + t78 * t79;
t95 = t171 + t176 / 0.2e1;
t99 = t172 - t177 / 0.2e1;
t199 = t95 * t114 + t99 * t79;
t198 = t122 * t76;
t128 = sin(qJ(3));
t197 = t128 * t78;
t195 = t74 * t128;
t194 = t99 * t128;
t193 = cos(pkin(6));
t192 = t122 * t127;
t191 = t122 * t132;
t190 = t123 * t130;
t188 = t135 * pkin(1) + pkin(9) * t190;
t182 = t205 / 0.2e1;
t178 = -t130 * pkin(1) + pkin(9) * t189;
t103 = pkin(3) * t182;
t80 = t103 - t203 / 0.2e1 + t124 * t125;
t170 = -t114 * t78 + t80 * t190 + t76 * t79 + t188;
t98 = t182 + t206;
t169 = t73 * t98 + t195;
t94 = t181 - t204 / 0.2e1;
t168 = -t74 * t133 + t73 * t94;
t167 = -t133 * t78 + t76 * t94;
t91 = t180 / 0.2e1 + t102;
t92 = t103 + t203 / 0.2e1;
t162 = pkin(3) * t197 + t91 * t190 + t76 * t92;
t161 = pkin(3) * t194 + t193 * t91 + t95 * t92;
t144 = -t116 * t78 + t90 * t190 + t76 * t89;
t159 = pkin(4) * t144 + t170;
t153 = -t74 * t114 + t80 * t189 + t73 * t79 + t178;
t97 = t206 - t205 / 0.2e1;
t152 = mrSges(4,1) * t97 - t124 * t207 - mrSges(3,3);
t147 = -pkin(3) * t195 - t91 * t189 - t73 * t92;
t146 = m(4) * pkin(2) + t133 * mrSges(4,1) - t128 * mrSges(4,2) + mrSges(3,1);
t142 = -t99 * t116 + t193 * t90 + t95 * t89;
t140 = -t94 * mrSges(4,1) - t98 * mrSges(4,2) - mrSges(3,2) + (-t202 + t207) * t122;
t139 = t166 / 0.2e1 + t155;
t138 = t165 / 0.2e1 + t154;
t137 = t123 * t138;
t115 = sin(t119);
t19 = t74 * t115 + t135 * t137 + t73 * t139;
t93 = t181 + t204 / 0.2e1;
t72 = -t95 * t122 + t193 * t124;
t60 = t124 * t190 - t198;
t42 = t116 * t95 + t89 * t99;
t37 = -t115 * t99 - t193 * t138 - t95 * t139;
t36 = t116 * t76 + t78 * t89;
t34 = -t116 * t73 - t74 * t89;
t27 = t190 * t93 + t76 * t98 + t197;
t24 = -t115 * t78 - t130 * t137 - t139 * t76;
t14 = t127 * t72 + t132 * t142;
t8 = t127 * t60 + t132 * t144;
t7 = t127 * t144 - t60 * t132;
t2 = t126 * t24 + t131 * t8;
t1 = -t126 * t8 + t131 * t24;
t3 = [(-mrSges(2,1) * t135 - m(3) * t188 + t78 * mrSges(3,1) - t76 * mrSges(3,2) - m(4) * (-pkin(2) * t78 - pkin(10) * t198 + t188) - t167 * mrSges(4,1) - t27 * mrSges(4,2) - m(5) * t170 - t144 * mrSges(5,1) - m(6) * t159 - t8 * mrSges(6,1) - m(7) * (pkin(5) * t8 + t159) - t2 * mrSges(7,1) - t1 * mrSges(7,2) + t173 * t7 + t202 * t60 + t160 * t24 + (t123 * t152 + mrSges(2,2)) * t130) * g(2) + (t130 * mrSges(2,1) - m(3) * t178 + t74 * mrSges(3,1) - t73 * mrSges(3,2) - m(4) * (-t74 * pkin(2) - pkin(10) * t196 + t178) - t168 * mrSges(4,1) - t169 * mrSges(4,2) - m(5) * t153 - t23 * mrSges(5,1) + t202 * t58 + t173 * t213 - t210 * t214 - t215 * t19 + (mrSges(2,2) + (-mrSges(4,2) * t93 + t152) * t123) * t135 + t212 * (-t23 * pkin(4) - t153)) * g(1) (-m(5) * t199 - t42 * mrSges(5,1) - t210 * (t132 * t42 - t192 * t99) + t215 * (t115 * t95 - t139 * t99) + t173 * (t127 * t42 + t191 * t99) - t146 * t95 + t140 * t99 - t212 * (t42 * pkin(4) + t199)) * g(3) + (-m(5) * t201 - t34 * mrSges(5,1) + t173 * (t127 * t34 - t191 * t74) - t210 * (t132 * t34 + t192 * t74) + t215 * (-t115 * t73 + t139 * t74) + t146 * t73 - t140 * t74 - t212 * (t34 * pkin(4) + t201)) * g(2) + (-m(5) * t200 - t36 * mrSges(5,1) - t210 * (t132 * t36 - t192 * t78) + t215 * (t115 * t76 - t139 * t78) + t173 * (t127 * t36 + t191 * t78) - t146 * t76 + t140 * t78 - t212 * (t36 * pkin(4) + t200)) * g(1) (-(t193 * t93 + t95 * t98 + t194) * mrSges(4,1) - (t99 * t133 + t193 * t97 - t95 * t94) * mrSges(4,2) - m(5) * t161 - t212 * (-t37 * pkin(4) + t161) + t215 * t142 + t209 * t37) * g(3) + (-(-t189 * t93 - t169) * mrSges(4,1) - (-t189 * t97 + t168) * mrSges(4,2) - m(5) * t147 - t212 * (-t19 * pkin(4) + t147) - t215 * t23 + t209 * t19) * g(2) + (-t27 * mrSges(4,1) - (t190 * t97 - t167) * mrSges(4,2) - m(5) * t162 - t212 * (-t24 * pkin(4) + t162) + t215 * t144 + t209 * t24) * g(1) (m(5) + t212) * (-g(1) * t60 + g(2) * t58 - g(3) * t72) (t173 * t14 - t210 * (-t127 * t142 + t132 * t72)) * g(3) + (-t173 * t214 - t210 * t213) * g(2) + (t173 * t8 + t210 * t7) * g(1), -g(1) * (mrSges(7,1) * t1 - mrSges(7,2) * t2) - g(2) * ((t126 * t214 + t131 * t19) * mrSges(7,1) + (-t126 * t19 + t131 * t214) * mrSges(7,2)) - g(3) * ((-t126 * t14 + t131 * t37) * mrSges(7,1) + (-t126 * t37 - t131 * t14) * mrSges(7,2))];
taug  = t3(:);
