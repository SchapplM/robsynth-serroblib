% Calculate Gravitation load on the joints for
% S6RRPRRR14
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% g [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [14x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,alpha2,alpha3,alpha4,d1,d2,d4,d5,d6,theta3]';
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
% Datum: 2018-12-10 18:39
% Revision: bb42a8b95257d9bc83910d26e849f5825122f662 (2018-12-05)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für mechatronische Systeme, Universität Hannover

function taug = S6RRPRRR14_gravloadJ_floatb_twist_slag_vp2(qJ, g, ...
  pkin, m, mrSges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(3,1),zeros(14,1),zeros(7,1),zeros(7,3)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRPRRR14_gravloadJ_floatb_twist_slag_vp2: qJ has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6RRPRRR14_gravloadJ_floatb_twist_slag_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [14 1]), ...
  'S6RRPRRR14_gravloadJ_floatb_twist_slag_vp2: pkin has to be [14x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6RRPRRR14_gravloadJ_floatb_twist_slag_vp2: m has to be [7x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [7,3]), ...
  'S6RRPRRR14_gravloadJ_floatb_twist_slag_vp2: mrSges has to be [7x3] (double)');

%% Symbolic Calculation
% From gravload_joint_floatb_twist_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2018-12-10 18:10:07
% EndTime: 2018-12-10 18:10:18
% DurationCPUTime: 3.16s
% Computational Cost: add. (9137->223), mult. (9359->306), div. (0->0), fcn. (9298->30), ass. (0->128)
t118 = sin(qJ(5));
t121 = cos(qJ(5));
t198 = pkin(8) + qJ(4);
t180 = cos(t198) / 0.2e1;
t199 = pkin(8) - qJ(4);
t192 = cos(t199);
t100 = t180 - t192 / 0.2e1;
t122 = cos(qJ(4));
t200 = pkin(6) + qJ(2);
t190 = sin(t200);
t178 = t190 / 0.2e1;
t201 = pkin(6) - qJ(2);
t191 = sin(t201);
t155 = t178 - t191 / 0.2e1;
t216 = sin(qJ(1));
t217 = cos(qJ(2));
t184 = t216 * t217;
t218 = cos(qJ(1));
t138 = t155 * t218 + t184;
t196 = pkin(7) + pkin(14);
t174 = sin(t196) / 0.2e1;
t197 = pkin(7) - pkin(14);
t186 = sin(t197);
t151 = t174 + t186 / 0.2e1;
t206 = sin(pkin(6));
t146 = t151 * t206;
t175 = cos(t197) / 0.2e1;
t187 = cos(t196);
t152 = t175 + t187 / 0.2e1;
t205 = sin(pkin(14));
t109 = cos(t200) / 0.2e1;
t193 = cos(t201);
t167 = t193 / 0.2e1 + t109;
t215 = sin(qJ(2));
t93 = -t167 * t218 + t216 * t215;
t220 = t138 * t205 + t146 * t218 + t93 * t152;
t114 = sin(pkin(7));
t207 = cos(pkin(7));
t176 = t207 * t206;
t231 = -t93 * t114 + t218 * t176;
t115 = cos(pkin(14));
t182 = t218 * t206;
t97 = t174 - t186 / 0.2e1;
t98 = t175 - t187 / 0.2e1;
t58 = t115 * t138 - t182 * t98 - t93 * t97;
t177 = sin(t198) / 0.2e1;
t189 = sin(t199);
t99 = t177 - t189 / 0.2e1;
t127 = t100 * t231 + t58 * t122 - t220 * t99;
t113 = sin(pkin(8));
t116 = cos(pkin(8));
t45 = -t113 * t220 + t116 * t231;
t236 = t118 * t45 - t121 * t127;
t235 = -t118 * t127 - t45 * t121;
t117 = sin(qJ(6));
t120 = cos(qJ(6));
t227 = m(6) + m(7);
t168 = -t227 * pkin(12) + mrSges(5,2) - mrSges(6,3);
t234 = -t117 * mrSges(7,1) - t120 * mrSges(7,2) + t168;
t225 = m(7) * pkin(5) + t120 * mrSges(7,1) - t117 * mrSges(7,2) + mrSges(6,1);
t188 = -m(7) * pkin(13) + mrSges(6,2) - mrSges(7,3);
t119 = sin(qJ(4));
t153 = t177 + t189 / 0.2e1;
t156 = t192 / 0.2e1 + t180;
t19 = t119 * t58 + t153 * t231 + t156 * t220;
t185 = t218 * t217;
t139 = -t155 * t216 + t185;
t145 = t216 * t167 + t215 * t218;
t129 = t139 * t205 + t145 * t152 - t146 * t216;
t232 = t145 * t114 + t216 * t176;
t233 = t129 * t113 + t116 * t232;
t230 = pkin(4) * t227 - t118 * t188 + t121 * t225 + mrSges(5,1);
t224 = m(5) + t227;
t179 = t191 / 0.2e1;
t166 = t179 - t190 / 0.2e1;
t94 = t166 * t218 - t184;
t70 = t152 * t94 + t205 * t93;
t213 = t113 * t70;
t95 = -t166 * t216 - t185;
t72 = t145 * t205 + t152 * t95;
t212 = t113 * t72;
t101 = t109 - t193 / 0.2e1;
t154 = t178 + t179;
t82 = t101 * t152 - t154 * t205;
t211 = t113 * t82;
t208 = cos(pkin(6));
t204 = t100 * t114;
t203 = t114 * t116;
t181 = t206 * t216;
t202 = t218 * pkin(1) + pkin(10) * t181;
t183 = -pkin(1) * t216 + pkin(10) * t182;
t71 = -t115 * t93 + t94 * t97;
t89 = t93 * pkin(2);
t171 = t71 * pkin(3) - pkin(11) * t213 - t89;
t73 = -t115 * t145 + t95 * t97;
t91 = t145 * pkin(2);
t170 = t73 * pkin(3) - pkin(11) * t212 - t91;
t83 = t101 * t97 + t115 * t154;
t96 = t154 * pkin(2);
t169 = t83 * pkin(3) - pkin(11) * t211 + t96;
t150 = -t138 * pkin(2) + qJ(3) * t231 + t183;
t147 = t114 * t153;
t143 = -pkin(3) * t58 + pkin(11) * t45 + t150;
t140 = t114 * t154 - t207 * t208;
t137 = t139 * pkin(2) + qJ(3) * t232 + t202;
t136 = -mrSges(3,2) + (m(4) * qJ(3) + mrSges(4,3) + t224 * (pkin(11) * t116 + qJ(3))) * t114;
t132 = t101 * t205 + t151 * t208 + t152 * t154;
t74 = -t101 * t115 + t154 * t97 + t208 * t98;
t130 = t100 * t140 + t74 * t122 + t132 * t99;
t61 = t115 * t139 - t145 * t97 + t181 * t98;
t126 = t61 * pkin(3) + pkin(11) * t233 + t137;
t123 = -t100 * t232 + t61 * t122 - t129 * t99;
t125 = pkin(4) * t123 + t126;
t67 = -t101 * t203 - t211;
t55 = -t113 * t132 - t116 * t140;
t49 = -t203 * t95 - t212;
t48 = -t203 * t94 - t213;
t43 = t101 * t204 + t122 * t83 + t82 * t99;
t36 = t119 * t74 - t132 * t156 + t140 * t153;
t34 = t122 * t73 + t204 * t95 + t72 * t99;
t32 = t122 * t71 + t204 * t94 + t70 * t99;
t24 = t119 * t61 + t129 * t156 - t153 * t232;
t14 = t118 * t55 + t121 * t130;
t8 = t118 * t233 + t121 * t123;
t7 = t118 * t123 - t121 * t233;
t2 = t117 * t24 + t120 * t8;
t1 = -t117 * t8 + t120 * t24;
t3 = [(-t218 * mrSges(2,1) + t216 * mrSges(2,2) - m(3) * t202 - t139 * mrSges(3,1) + t145 * mrSges(3,2) - mrSges(3,3) * t181 - m(4) * t137 - t61 * mrSges(4,1) + t129 * mrSges(4,2) - t232 * mrSges(4,3) - m(5) * t126 - t123 * mrSges(5,1) - t233 * mrSges(5,3) - m(6) * t125 - t8 * mrSges(6,1) - m(7) * (t8 * pkin(5) + t125) - t2 * mrSges(7,1) - t1 * mrSges(7,2) + t188 * t7 + t168 * t24) * g(2) + (t216 * mrSges(2,1) + t218 * mrSges(2,2) - m(3) * t183 + t138 * mrSges(3,1) - t93 * mrSges(3,2) - mrSges(3,3) * t182 - m(4) * t150 + t58 * mrSges(4,1) - t220 * mrSges(4,2) - t231 * mrSges(4,3) - m(5) * t143 + t127 * mrSges(5,1) - t45 * mrSges(5,3) + t188 * t235 - t225 * t236 - t234 * t19 + t227 * (pkin(4) * t127 - t143)) * g(1) (-t154 * mrSges(3,1) - m(4) * t96 - t83 * mrSges(4,1) - t82 * mrSges(4,2) - m(5) * t169 - t43 * mrSges(5,1) - t67 * mrSges(5,3) - t225 * (t118 * t67 + t121 * t43) + t234 * (t101 * t147 + t83 * t119 - t156 * t82) + t188 * (t118 * t43 - t67 * t121) + t136 * t101 + t227 * (-t43 * pkin(4) - t169)) * g(3) + (mrSges(3,1) * t93 + m(4) * t89 - t71 * mrSges(4,1) - t70 * mrSges(4,2) - m(5) * t171 - t32 * mrSges(5,1) - t48 * mrSges(5,3) + t188 * (t118 * t32 - t48 * t121) - t225 * (t118 * t48 + t121 * t32) + t234 * (t71 * t119 + t147 * t94 - t156 * t70) + t136 * t94 + t227 * (-t32 * pkin(4) - t171)) * g(2) + (t145 * mrSges(3,1) + m(4) * t91 - t73 * mrSges(4,1) - t72 * mrSges(4,2) - m(5) * t170 - t34 * mrSges(5,1) - t49 * mrSges(5,3) - t225 * (t118 * t49 + t121 * t34) + t234 * (t73 * t119 + t147 * t95 - t156 * t72) + t188 * (t118 * t34 - t49 * t121) + t136 * t95 + t227 * (-t34 * pkin(4) - t170)) * g(1) (-g(1) * t232 + g(2) * t231 + g(3) * t140) * (m(4) + t224) (t130 * t234 + t230 * t36) * g(3) + (t127 * t234 + t19 * t230) * g(2) + (t123 * t234 + t230 * t24) * g(1) (t188 * t14 - t225 * (-t118 * t130 + t121 * t55)) * g(3) + (-t188 * t236 - t225 * t235) * g(2) + (t188 * t8 + t225 * t7) * g(1), -g(1) * (mrSges(7,1) * t1 - mrSges(7,2) * t2) - g(2) * ((t117 * t236 + t120 * t19) * mrSges(7,1) + (-t117 * t19 + t120 * t236) * mrSges(7,2)) - g(3) * ((-t117 * t14 + t120 * t36) * mrSges(7,1) + (-t117 * t36 - t120 * t14) * mrSges(7,2))];
taug  = t3(:);
