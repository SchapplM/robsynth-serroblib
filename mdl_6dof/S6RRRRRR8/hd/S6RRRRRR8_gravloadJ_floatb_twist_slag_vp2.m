% Calculate Gravitation load on the joints for
% S6RRRRRR8
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% g [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [13x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,alpha2,alpha3,d1,d2,d3,d4,d5,d6]';
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
% Datum: 2019-03-10 05:15
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function taug = S6RRRRRR8_gravloadJ_floatb_twist_slag_vp2(qJ, g, ...
  pkin, m, mrSges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(3,1),zeros(13,1),zeros(7,1),zeros(7,3)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRRRRR8_gravloadJ_floatb_twist_slag_vp2: qJ has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6RRRRRR8_gravloadJ_floatb_twist_slag_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [13 1]), ...
  'S6RRRRRR8_gravloadJ_floatb_twist_slag_vp2: pkin has to be [13x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6RRRRRR8_gravloadJ_floatb_twist_slag_vp2: m has to be [7x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [7,3]), ...
  'S6RRRRRR8_gravloadJ_floatb_twist_slag_vp2: mrSges has to be [7x3] (double)');

%% Symbolic Calculation
% From gravload_joint_floatb_twist_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-10 04:48:19
% EndTime: 2019-03-10 04:48:26
% DurationCPUTime: 2.34s
% Computational Cost: add. (1608->217), mult. (3795->316), div. (0->0), fcn. (4776->16), ass. (0->103)
t209 = mrSges(6,2) - mrSges(7,3);
t107 = sin(qJ(6));
t110 = cos(qJ(6));
t202 = t110 * mrSges(7,1) - t107 * mrSges(7,2) + mrSges(6,1);
t212 = m(7) * pkin(5) + t202;
t199 = -m(7) * pkin(13) + t209;
t108 = sin(qJ(4));
t111 = cos(qJ(4));
t105 = sin(pkin(7));
t112 = cos(qJ(1));
t106 = sin(pkin(6));
t170 = cos(pkin(7));
t154 = t106 * t170;
t171 = cos(pkin(6));
t191 = cos(qJ(2));
t140 = t171 * t191;
t188 = sin(qJ(2));
t189 = sin(qJ(1));
t83 = -t112 * t140 + t188 * t189;
t206 = -t83 * t105 + t112 * t154;
t109 = sin(qJ(3));
t153 = t109 * t170;
t165 = t106 * t112;
t162 = t105 * t165;
t190 = cos(qJ(3));
t139 = t171 * t188;
t84 = t112 * t139 + t189 * t191;
t41 = -t109 * t162 - t153 * t83 + t190 * t84;
t211 = t108 * t41 + t111 * t206;
t178 = t108 * t206;
t210 = -t111 * t41 + t178;
t123 = t112 * t188 + t140 * t189;
t207 = t123 * t105 + t189 * t154;
t157 = t106 * t189;
t198 = -t105 * t157 + t123 * t170;
t85 = t112 * t191 - t189 * t139;
t45 = -t198 * t109 + t85 * t190;
t19 = -t45 * t108 + t111 * t207;
t155 = t105 * t171;
t65 = t109 * t155 + (t153 * t191 + t190 * t188) * t106;
t158 = t106 * t191;
t82 = -t105 * t158 + t170 * t171;
t205 = -t108 * t65 + t111 * t82;
t104 = qJ(4) + qJ(5);
t101 = sin(t104);
t102 = cos(t104);
t14 = -t101 * t206 + t102 * t41;
t13 = -t101 * t41 - t102 * t206;
t200 = -m(6) - m(7);
t136 = -m(5) * pkin(11) + mrSges(4,2) - mrSges(5,3) - mrSges(6,3);
t197 = -t107 * mrSges(7,1) - t110 * mrSges(7,2) + t136;
t194 = m(5) * pkin(3) + t111 * mrSges(5,1) - t108 * mrSges(5,2) - t199 * t101 + t212 * t102 + mrSges(4,1);
t156 = t106 * t188;
t146 = t105 * t156;
t184 = pkin(2) * t158 + pkin(10) * t146;
t183 = t105 * t84;
t182 = t105 * t85;
t172 = t112 * pkin(1) + pkin(9) * t157;
t169 = t101 * t105;
t168 = t102 * t105;
t167 = t105 * t108;
t166 = t105 * t111;
t161 = t13 * pkin(5) + pkin(13) * t14;
t17 = t101 * t45 - t102 * t207;
t18 = t101 * t207 + t45 * t102;
t160 = -t17 * pkin(5) + pkin(13) * t18;
t34 = -t101 * t65 + t102 * t82;
t35 = t101 * t82 + t102 * t65;
t159 = t34 * pkin(5) + pkin(13) * t35;
t151 = t211 * pkin(4);
t150 = t19 * pkin(4);
t149 = t205 * pkin(4);
t78 = t83 * pkin(2);
t148 = pkin(10) * t183 - t78;
t80 = t123 * pkin(2);
t147 = pkin(10) * t182 - t80;
t144 = -pkin(1) * t189 + pkin(9) * t165;
t141 = t105 * (pkin(4) * t108 + pkin(10));
t138 = t170 * t190;
t137 = t170 * t188;
t133 = t108 * t146;
t131 = -t202 * t13 + t209 * t14;
t130 = t202 * t17 + t209 * t18;
t129 = -t202 * t34 + t209 * t35;
t124 = -t84 * pkin(2) + t206 * pkin(10) + t144;
t40 = t84 * t109 + t138 * t83 + t162 * t190;
t118 = t85 * pkin(2) + t207 * pkin(10) + t172;
t100 = pkin(4) * t111 + pkin(3);
t113 = -pkin(12) - pkin(11);
t115 = t207 * t108;
t44 = t109 * t85 + t198 * t190;
t116 = pkin(4) * t115 + t45 * t100 - t44 * t113 + t118;
t76 = (-t109 * t137 + t190 * t191) * t106;
t75 = (t109 * t191 + t190 * t137) * t106;
t64 = t109 * t156 - t138 * t158 - t155 * t190;
t53 = -t123 * t190 - t153 * t85;
t52 = -t109 * t123 + t138 * t85;
t51 = -t153 * t84 - t190 * t83;
t50 = -t109 * t83 + t138 * t84;
t20 = t45 * t111 + t115;
t2 = t107 * t44 + t110 * t18;
t1 = -t107 * t18 + t110 * t44;
t3 = [(-t112 * mrSges(2,1) + t189 * mrSges(2,2) - m(3) * t172 - t85 * mrSges(3,1) + t123 * mrSges(3,2) - mrSges(3,3) * t157 - m(4) * t118 - t45 * mrSges(4,1) - t207 * mrSges(4,3) - m(5) * (t45 * pkin(3) + t118) - t20 * mrSges(5,1) - t19 * mrSges(5,2) - m(6) * t116 - t18 * mrSges(6,1) - m(7) * (t18 * pkin(5) + t116) - t2 * mrSges(7,1) - t1 * mrSges(7,2) + t199 * t17 + t136 * t44) * g(2) + (t189 * mrSges(2,1) - m(3) * t144 + t84 * mrSges(3,1) - t83 * mrSges(3,2) - m(4) * t124 + t41 * mrSges(4,1) - t206 * mrSges(4,3) - m(5) * (-pkin(3) * t41 + t124) - t210 * mrSges(5,1) - t211 * mrSges(5,2) + t199 * t13 + (-mrSges(3,3) * t106 + mrSges(2,2)) * t112 + t212 * t14 - t197 * t40 + t200 * (pkin(4) * t178 - t100 * t41 + t113 * t40 + t124)) * g(1) (-(mrSges(3,1) * t191 - mrSges(3,2) * t188) * t106 - m(4) * t184 - t76 * mrSges(4,1) - mrSges(4,3) * t146 - m(5) * (pkin(3) * t76 + t184) - (t76 * t111 + t133) * mrSges(5,1) - (-t76 * t108 + t111 * t146) * mrSges(5,2) + t200 * (pkin(4) * t133 + t76 * t100 - t75 * t113 + t184) + t199 * (t101 * t76 - t102 * t146) - t212 * (t101 * t146 + t76 * t102) + t197 * t75) * g(3) + (mrSges(3,1) * t83 + mrSges(3,2) * t84 - m(4) * t148 - t51 * mrSges(4,1) - mrSges(4,3) * t183 - m(5) * (pkin(3) * t51 + t148) - (t111 * t51 + t167 * t84) * mrSges(5,1) - (-t108 * t51 + t166 * t84) * mrSges(5,2) + t200 * (t51 * t100 - t50 * t113 + t141 * t84 - t78) + t199 * (t101 * t51 - t168 * t84) - t212 * (t102 * t51 + t169 * t84) + t197 * t50) * g(2) + (t123 * mrSges(3,1) + t85 * mrSges(3,2) - m(4) * t147 - t53 * mrSges(4,1) - mrSges(4,3) * t182 - m(5) * (pkin(3) * t53 + t147) - (t111 * t53 + t167 * t85) * mrSges(5,1) - (-t108 * t53 + t166 * t85) * mrSges(5,2) + t200 * (t53 * t100 - t52 * t113 + t141 * t85 - t80) + t199 * (t101 * t53 - t168 * t85) - t212 * (t102 * t53 + t169 * t85) + t197 * t52) * g(1) (t200 * (-t64 * t100 - t65 * t113) + t197 * t65 + t194 * t64) * g(3) + (t200 * (-t40 * t100 - t41 * t113) + t197 * t41 + t194 * t40) * g(2) + (t200 * (-t44 * t100 - t45 * t113) + t197 * t45 + t194 * t44) * g(1) (-t205 * mrSges(5,1) - (-t108 * t82 - t111 * t65) * mrSges(5,2) - m(6) * t149 - m(7) * (t149 + t159) + t129) * g(3) + (t211 * mrSges(5,1) - t210 * mrSges(5,2) + m(6) * t151 - m(7) * (-t151 + t161) + t131) * g(2) + (-t19 * mrSges(5,1) + t20 * mrSges(5,2) - m(6) * t150 - m(7) * (t150 + t160) + t130) * g(1) (-m(7) * t159 + t129) * g(3) + (-m(7) * t161 + t131) * g(2) + (-m(7) * t160 + t130) * g(1), -g(1) * (mrSges(7,1) * t1 - mrSges(7,2) * t2) - g(2) * ((-t107 * t14 + t110 * t40) * mrSges(7,1) + (-t107 * t40 - t110 * t14) * mrSges(7,2)) - g(3) * ((-t107 * t35 + t110 * t64) * mrSges(7,1) + (-t107 * t64 - t110 * t35) * mrSges(7,2))];
taug  = t3(:);
