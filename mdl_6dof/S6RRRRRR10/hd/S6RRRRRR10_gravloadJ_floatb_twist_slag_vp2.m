% Calculate Gravitation load on the joints for
% S6RRRRRR10
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% g [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [14x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,alpha2,alpha3,alpha4,d1,d2,d3,d4,d5,d6]';
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
% Datum: 2018-11-23 11:28
% Revision: 76f9d5e39f14dc242b53c0d9d3d9db48bd8f37c0
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für mechatronische Systeme, Universität Hannover

function taug = S6RRRRRR10_gravloadJ_floatb_twist_slag_vp2(qJ, g, ...
  pkin, m, mrSges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(3,1),zeros(14,1),zeros(7,1),zeros(7,3)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRRRRR10_gravloadJ_floatb_twist_slag_vp2: qJ has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6RRRRRR10_gravloadJ_floatb_twist_slag_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [14 1]), ...
  'S6RRRRRR10_gravloadJ_floatb_twist_slag_vp2: pkin has to be [14x1] (double)');
assert( isreal(m) && all(size(m) == [7 1]), ...
  'S6RRRRRR10_gravloadJ_floatb_twist_slag_vp2: m has to be [7x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [7,3]), ...
  'S6RRRRRR10_gravloadJ_floatb_twist_slag_vp2: mrSges has to be [7x3] (double)');

%% Symbolic Calculation
% From joint_gravload_floatb_twist_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2018-11-23 10:29:17
% EndTime: 2018-11-23 10:29:24
% DurationCPUTime: 3.72s
% Computational Cost: add. (11565->257), mult. (11696->360), div. (0->0), fcn. (11549->30), ass. (0->130)
t143 = sin(qJ(5));
t146 = cos(qJ(5));
t217 = pkin(8) + qJ(4);
t199 = sin(t217) / 0.2e1;
t218 = pkin(8) - qJ(4);
t208 = sin(t218);
t122 = t199 - t208 / 0.2e1;
t202 = cos(t217) / 0.2e1;
t211 = cos(t218);
t125 = t202 - t211 / 0.2e1;
t147 = cos(qJ(4));
t150 = cos(qJ(1));
t221 = pkin(6) + qJ(2);
t204 = cos(t221) / 0.2e1;
t222 = pkin(6) - qJ(2);
t213 = cos(t222);
t180 = t213 / 0.2e1 + t204;
t242 = sin(qJ(2));
t243 = sin(qJ(1));
t117 = -t150 * t180 + t242 * t243;
t140 = sin(pkin(7));
t231 = sin(pkin(6));
t232 = cos(pkin(7));
t198 = t232 * t231;
t252 = -t117 * t140 + t150 * t198;
t201 = sin(t221) / 0.2e1;
t210 = sin(t222);
t124 = t201 - t210 / 0.2e1;
t149 = cos(qJ(2));
t118 = t150 * t124 + t149 * t243;
t219 = pkin(7) + qJ(3);
t200 = sin(t219) / 0.2e1;
t220 = pkin(7) - qJ(3);
t209 = sin(t220);
t123 = t200 - t209 / 0.2e1;
t203 = cos(t219) / 0.2e1;
t212 = cos(t220);
t126 = t203 - t212 / 0.2e1;
t148 = cos(qJ(3));
t214 = t150 * t231;
t77 = -t117 * t123 + t118 * t148 + t126 * t214;
t176 = t200 + t209 / 0.2e1;
t170 = t176 * t231;
t179 = t212 / 0.2e1 + t203;
t241 = sin(qJ(3));
t79 = t117 * t179 + t118 * t241 + t150 * t170;
t27 = t122 * t79 - t125 * t252 - t77 * t147;
t139 = sin(pkin(8));
t141 = cos(pkin(8));
t60 = -t139 * t79 + t141 * t252;
t260 = t143 * t60 + t146 * t27;
t259 = t143 * t27 - t60 * t146;
t142 = sin(qJ(6));
t145 = cos(qJ(6));
t250 = m(6) + m(7);
t195 = -pkin(13) * t250 + mrSges(5,2) - mrSges(6,3);
t256 = -t142 * mrSges(7,1) - t145 * mrSges(7,2) + t195;
t255 = m(5) + t250;
t248 = m(7) * pkin(5) + t145 * mrSges(7,1) - t142 * mrSges(7,2) + mrSges(6,1);
t207 = -m(7) * pkin(14) + mrSges(6,2) - mrSges(7,3);
t166 = t150 * t242 + t180 * t243;
t187 = -t124 * t243 + t150 * t149;
t157 = t166 * t179 - t170 * t243 + t187 * t241;
t253 = t166 * t140 + t243 * t198;
t254 = t157 * t139 + t141 * t253;
t251 = pkin(4) * t250 - t207 * t143 + t248 * t146 + mrSges(5,1);
t144 = sin(qJ(4));
t175 = t199 + t208 / 0.2e1;
t178 = t211 / 0.2e1 + t202;
t23 = t144 * t77 + t175 * t252 + t79 * t178;
t92 = t117 * t241 - t118 * t179;
t235 = t139 * t92;
t94 = t166 * t241 - t179 * t187;
t234 = t139 * t94;
t233 = cos(pkin(6));
t127 = t204 - t213 / 0.2e1;
t177 = t201 + t210 / 0.2e1;
t106 = t127 * t179 - t177 * t241;
t230 = t106 * t139;
t227 = t139 * t143;
t226 = t139 * t146;
t225 = t140 * t125;
t224 = t140 * t141;
t205 = t231 * t243;
t223 = t150 * pkin(1) + pkin(10) * t205;
t206 = -pkin(1) * t243 + pkin(10) * t214;
t113 = t117 * pkin(2);
t93 = -t117 * t148 - t118 * t123;
t194 = t93 * pkin(3) - pkin(12) * t235 - t113;
t115 = t166 * pkin(2);
t95 = -t123 * t187 - t148 * t166;
t193 = t95 * pkin(3) - pkin(12) * t234 - t115;
t107 = t127 * t123 + t148 * t177;
t121 = t177 * pkin(2);
t192 = t107 * pkin(3) - pkin(12) * t230 + t121;
t182 = -mrSges(4,2) + (t255 * pkin(12) + mrSges(5,3)) * t139;
t174 = -t118 * pkin(2) + pkin(11) * t252 + t206;
t171 = t140 * t175;
t168 = -t77 * pkin(3) + pkin(12) * t60 + t174;
t165 = t140 * t177 - t232 * t233;
t98 = -t123 * t177 + t126 * t233 + t127 * t148;
t163 = -mrSges(3,2) + (m(4) * pkin(11) + mrSges(4,3) + t255 * (pkin(12) * t141 + pkin(11))) * t140;
t162 = t187 * pkin(2) + pkin(11) * t253 + t223;
t81 = -t123 * t166 - t126 * t205 + t148 * t187;
t159 = t127 * t241 + t176 * t233 + t177 * t179;
t155 = t122 * t159 + t125 * t165 - t147 * t98;
t154 = t81 * pkin(3) + pkin(12) * t254 + t162;
t151 = -t122 * t157 - t125 * t253 + t81 * t147;
t153 = pkin(4) * t151 + t154;
t96 = t159 * pkin(3);
t89 = -t127 * t224 - t230;
t75 = t157 * pkin(3);
t73 = t79 * pkin(3);
t72 = -t139 * t159 - t141 * t165;
t64 = t187 * t224 - t234;
t63 = t118 * t224 - t235;
t58 = t106 * t122 + t107 * t147 + t127 * t225;
t53 = t98 * t122 + t147 * t159;
t48 = -t144 * t98 - t159 * t178 + t165 * t175;
t46 = -t122 * t81 - t147 * t157;
t44 = -t122 * t77 - t147 * t79;
t42 = t122 * t94 + t147 * t95 - t187 * t225;
t40 = -t118 * t225 + t122 * t92 + t147 * t93;
t28 = t144 * t81 + t157 * t178 - t175 * t253;
t18 = t143 * t72 + t146 * t155;
t8 = t143 * t254 + t146 * t151;
t7 = t143 * t151 - t146 * t254;
t2 = t142 * t28 + t145 * t8;
t1 = -t142 * t8 + t145 * t28;
t3 = [(-t150 * mrSges(2,1) + t243 * mrSges(2,2) - m(3) * t223 - t187 * mrSges(3,1) + t166 * mrSges(3,2) - mrSges(3,3) * t205 - m(4) * t162 - t81 * mrSges(4,1) + t157 * mrSges(4,2) - t253 * mrSges(4,3) - m(5) * t154 - t151 * mrSges(5,1) - t254 * mrSges(5,3) - m(6) * t153 - t8 * mrSges(6,1) - m(7) * (t8 * pkin(5) + t153) - t2 * mrSges(7,1) - t1 * mrSges(7,2) + t207 * t7 + t195 * t28) * g(2) + (t243 * mrSges(2,1) - m(3) * t206 + t118 * mrSges(3,1) - t117 * mrSges(3,2) - m(4) * t174 + t77 * mrSges(4,1) - t79 * mrSges(4,2) - t252 * mrSges(4,3) - m(5) * t168 - t27 * mrSges(5,1) - t60 * mrSges(5,3) + t207 * t259 - t248 * t260 - t256 * t23 + (-mrSges(3,3) * t231 + mrSges(2,2)) * t150 + t250 * (-t27 * pkin(4) - t168)) * g(1) (-t177 * mrSges(3,1) - m(4) * t121 - t107 * mrSges(4,1) - t106 * mrSges(4,2) - m(5) * t192 - t58 * mrSges(5,1) - t89 * mrSges(5,3) - t248 * (t143 * t89 + t146 * t58) + t256 * (-t106 * t178 + t107 * t144 + t127 * t171) + t207 * (t143 * t58 - t89 * t146) + t163 * t127 + t250 * (-t58 * pkin(4) - t192)) * g(3) + (mrSges(3,1) * t117 + m(4) * t113 - t93 * mrSges(4,1) - t92 * mrSges(4,2) - m(5) * t194 - t40 * mrSges(5,1) - t63 * mrSges(5,3) + t207 * (t143 * t40 - t63 * t146) - t248 * (t143 * t63 + t146 * t40) + t256 * (-t118 * t171 + t93 * t144 - t178 * t92) - t163 * t118 + t250 * (-t40 * pkin(4) - t194)) * g(2) + (t166 * mrSges(3,1) + m(4) * t115 - t95 * mrSges(4,1) - t94 * mrSges(4,2) - m(5) * t193 - t42 * mrSges(5,1) - t64 * mrSges(5,3) - t248 * (t143 * t64 + t146 * t42) + t256 * (t95 * t144 - t171 * t187 - t178 * t94) + t207 * (t143 * t42 - t64 * t146) - t163 * t187 + t250 * (-t42 * pkin(4) - t193)) * g(1) (-t159 * mrSges(4,1) - m(5) * t96 - t53 * mrSges(5,1) + t182 * t98 - t248 * (t146 * t53 - t227 * t98) + t256 * (t144 * t159 - t178 * t98) + t207 * (t143 * t53 + t226 * t98) - t250 * (t53 * pkin(4) + t96)) * g(3) + (t79 * mrSges(4,1) + m(5) * t73 - t44 * mrSges(5,1) - t182 * t77 - t248 * (t146 * t44 + t227 * t77) + t256 * (-t144 * t79 + t178 * t77) + t207 * (t143 * t44 - t226 * t77) - t250 * (t44 * pkin(4) - t73)) * g(2) + (t157 * mrSges(4,1) + m(5) * t75 - t46 * mrSges(5,1) - t182 * t81 - t248 * (t146 * t46 + t227 * t81) + t256 * (-t144 * t157 + t178 * t81) + t207 * (t143 * t46 - t226 * t81) - t250 * (t46 * pkin(4) - t75)) * g(1) (t155 * t256 + t251 * t48) * g(3) + (t23 * t251 - t256 * t27) * g(2) + (t151 * t256 + t251 * t28) * g(1) (t207 * t18 - t248 * (-t143 * t155 + t146 * t72)) * g(3) + (-t207 * t260 - t248 * t259) * g(2) + (t207 * t8 + t248 * t7) * g(1), -g(1) * (mrSges(7,1) * t1 - mrSges(7,2) * t2) - g(2) * ((t142 * t260 + t145 * t23) * mrSges(7,1) + (-t142 * t23 + t145 * t260) * mrSges(7,2)) - g(3) * ((-t142 * t18 + t145 * t48) * mrSges(7,1) + (-t142 * t48 - t145 * t18) * mrSges(7,2))];
taug  = t3(:);
