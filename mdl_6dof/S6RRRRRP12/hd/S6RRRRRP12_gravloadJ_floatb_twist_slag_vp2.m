% Calculate Gravitation load on the joints for
% S6RRRRRP12
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

% Quelle: HybrDyn-Toolbox (ehem. IRT-Maple-Toolbox)
% Datum: 2018-11-23 18:36
% Revision: 76f9d5e39f14dc242b53c0d9d3d9db48bd8f37c0
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für mechatronische Systeme, Universität Hannover

function taug = S6RRRRRP12_gravloadJ_floatb_twist_slag_vp2(qJ, g, ...
  pkin, m, mrSges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(3,1),zeros(12,1),zeros(7,1),zeros(7,3)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRRRRP12_gravloadJ_floatb_twist_slag_vp2: qJ has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6RRRRRP12_gravloadJ_floatb_twist_slag_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [12 1]), ...
  'S6RRRRRP12_gravloadJ_floatb_twist_slag_vp2: pkin has to be [12x1] (double)');
assert( isreal(m) && all(size(m) == [7 1]), ...
  'S6RRRRRP12_gravloadJ_floatb_twist_slag_vp2: m has to be [7x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [7,3]), ...
  'S6RRRRRP12_gravloadJ_floatb_twist_slag_vp2: mrSges has to be [7x3] (double)');

%% Symbolic Calculation
% From joint_gravload_floatb_twist_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2018-11-23 18:35:56
% EndTime: 2018-11-23 18:35:58
% DurationCPUTime: 2.30s
% Computational Cost: add. (5867->215), mult. (6215->302), div. (0->0), fcn. (6195->22), ass. (0->118)
t124 = sin(qJ(5));
t127 = cos(qJ(5));
t125 = sin(qJ(4));
t128 = cos(qJ(4));
t123 = sin(pkin(7));
t131 = cos(qJ(1));
t197 = pkin(6) + qJ(2);
t175 = cos(t197) / 0.2e1;
t198 = pkin(6) - qJ(2);
t184 = cos(t198);
t154 = t184 / 0.2e1 + t175;
t218 = sin(qJ(2));
t219 = sin(qJ(1));
t145 = -t131 * t154 + t219 * t218;
t208 = sin(pkin(6));
t209 = cos(pkin(7));
t167 = t209 * t208;
t230 = t145 * t123 - t131 * t167;
t196 = pkin(7) - qJ(3);
t183 = cos(t196);
t174 = t183 / 0.2e1;
t195 = pkin(7) + qJ(3);
t182 = cos(t195);
t162 = t174 - t182 / 0.2e1;
t156 = t162 * t208;
t172 = sin(t197) / 0.2e1;
t181 = sin(t198);
t111 = t172 - t181 / 0.2e1;
t130 = cos(qJ(2));
t102 = t131 * t111 + t130 * t219;
t129 = cos(qJ(3));
t171 = sin(t195) / 0.2e1;
t180 = sin(t196);
t223 = t171 - t180 / 0.2e1;
t228 = t102 * t129 - t145 * t223;
t57 = t131 * t156 - t228;
t24 = t125 * t230 - t57 * t128;
t126 = sin(qJ(3));
t149 = t171 + t180 / 0.2e1;
t146 = t149 * t208;
t173 = t182 / 0.2e1;
t153 = t174 + t173;
t53 = t102 * t126 + t131 * t146 + t145 * t153;
t1 = t124 * t24 - t53 * t127;
t234 = t124 * t53 + t127 * t24;
t179 = m(7) * pkin(5) + mrSges(6,1) + mrSges(7,1);
t178 = -m(7) * qJ(6) + mrSges(6,2) - mrSges(7,3);
t23 = t125 * t57 + t128 * t230;
t142 = t131 * t218 + t154 * t219;
t231 = t142 * t123 + t219 * t167;
t104 = -t111 * t219 + t131 * t130;
t229 = t104 * t129 - t142 * t223;
t112 = t175 - t184 / 0.2e1;
t151 = t172 + t181 / 0.2e1;
t227 = -t112 * t129 + t151 * t223;
t220 = -m(6) - m(7);
t226 = mrSges(4,2) - mrSges(5,3);
t225 = mrSges(6,3) + mrSges(7,2);
t224 = t128 * mrSges(5,1) - mrSges(5,2) * t125 + mrSges(4,1);
t222 = mrSges(5,2) - t225;
t221 = t178 * t124 - t179 * t127 - mrSges(5,1);
t217 = pkin(4) * t128;
t69 = -t102 * t223 - t129 * t145;
t97 = t145 * pkin(2);
t216 = t69 * pkin(3) - t97;
t71 = -t104 * t223 - t129 * t142;
t99 = t142 * pkin(2);
t215 = t71 * pkin(3) - t99;
t214 = t125 * t53;
t58 = t104 * t126 + t142 * t153 - t146 * t219;
t213 = t125 * t58;
t210 = cos(pkin(6));
t77 = -t112 * t126 - t149 * t210 - t151 * t153;
t212 = t125 * t77;
t110 = t151 * pkin(2);
t83 = t112 * t223 + t129 * t151;
t211 = t83 * pkin(3) + t110;
t204 = t123 * t125;
t203 = t123 * t128;
t202 = t124 * t128;
t201 = t127 * t128;
t176 = t208 * t219;
t200 = t131 * pkin(1) + pkin(9) * t176;
t199 = -m(5) + t220;
t152 = t173 - t183 / 0.2e1;
t147 = t152 * t208;
t55 = t131 * t147 + t228;
t191 = -t53 * pkin(3) + pkin(11) * t55;
t60 = -t147 * t219 + t229;
t190 = -t58 * pkin(3) + pkin(11) * t60;
t79 = -t152 * t210 + t227;
t189 = -t77 * pkin(3) + pkin(11) * t79;
t185 = t131 * t208;
t177 = -pkin(1) * t219 + pkin(9) * t185;
t160 = pkin(11) * t199 + t226;
t158 = pkin(12) * t220 + t222;
t155 = -mrSges(3,2) + (mrSges(4,3) + (m(4) - t199) * pkin(10)) * t123;
t138 = -t102 * pkin(2) - pkin(10) * t230 + t177;
t137 = t57 * pkin(3) + t138;
t136 = t104 * pkin(2) + pkin(10) * t231 + t200;
t59 = t156 * t219 + t229;
t134 = t59 * pkin(3) + t136;
t101 = -t123 * t151 + t209 * t210;
t82 = -t112 * t153 + t126 * t151;
t78 = t162 * t210 + t227;
t70 = t104 * t153 - t126 * t142;
t68 = t102 * t153 - t126 * t145;
t65 = -t112 * t204 + t128 * t83;
t39 = t101 * t125 + t128 * t78;
t38 = t101 * t128 - t125 * t78;
t36 = t104 * t204 + t128 * t71;
t34 = t102 * t204 + t128 * t69;
t28 = t125 * t231 + t59 * t128;
t27 = t125 * t59 - t128 * t231;
t15 = t124 * t39 - t77 * t127;
t6 = t124 * t58 + t127 * t28;
t5 = t124 * t28 - t58 * t127;
t2 = [(-m(3) * t200 - m(4) * t136 - m(5) * t134 - t131 * mrSges(2,1) - t104 * mrSges(3,1) - t59 * mrSges(4,1) - t28 * mrSges(5,1) + t219 * mrSges(2,2) + t142 * mrSges(3,2) - mrSges(3,3) * t176 - t231 * mrSges(4,3) + t158 * t27 + t160 * t58 + t178 * t5 - t179 * t6 + t220 * (t28 * pkin(4) + t134)) * g(2) + (t219 * mrSges(2,1) + t131 * mrSges(2,2) - m(3) * t177 + t102 * mrSges(3,1) - t145 * mrSges(3,2) - mrSges(3,3) * t185 - m(4) * t138 - t57 * mrSges(4,1) + t230 * mrSges(4,3) - m(5) * t137 + t24 * mrSges(5,1) - t160 * t53 + t179 * t234 - t178 * t1 + t158 * t23 + t220 * (-pkin(4) * t24 + t137)) * g(1) (-t151 * mrSges(3,1) - m(4) * t110 - t83 * mrSges(4,1) - m(5) * t211 - t65 * mrSges(5,1) + t160 * t82 - t179 * (t124 * t82 + t127 * t65) + t178 * (t124 * t65 - t82 * t127) + t155 * t112 + t158 * (t112 * t203 + t125 * t83) + t220 * (t65 * pkin(4) + t211)) * g(3) + (t145 * mrSges(3,1) + m(4) * t97 - t69 * mrSges(4,1) - m(5) * t216 - t34 * mrSges(5,1) + t160 * t68 - t179 * (t124 * t68 + t127 * t34) + t178 * (t124 * t34 - t68 * t127) - t155 * t102 + t158 * (-t102 * t203 + t125 * t69) + t220 * (t34 * pkin(4) + t216)) * g(2) + (t142 * mrSges(3,1) + m(4) * t99 - t71 * mrSges(4,1) - m(5) * t215 - t36 * mrSges(5,1) + t160 * t70 - t179 * (t124 * t70 + t127 * t36) + t178 * (t124 * t36 - t70 * t127) - t155 * t104 + t158 * (-t104 * t203 + t125 * t71) + t220 * (t36 * pkin(4) + t215)) * g(1) (-m(5) * t189 + t226 * t79 + t224 * t77 + t220 * (-pkin(12) * t212 - t77 * t217 + t189) - t179 * (t124 * t79 - t201 * t77) + t225 * t212 + t178 * (-t79 * t127 - t202 * t77)) * g(3) + (-m(5) * t191 + t220 * (-pkin(12) * t214 - t53 * t217 + t191) - t179 * (t124 * t55 - t201 * t53) + t178 * (-t55 * t127 - t202 * t53) + t226 * t55 + t224 * t53 + t225 * t214) * g(2) + (-m(5) * t190 + t220 * (-pkin(12) * t213 - t58 * t217 + t190) + t178 * (-t60 * t127 - t202 * t58) + t226 * t60 + t224 * t58 + t225 * t213 - t179 * (t124 * t60 - t201 * t58)) * g(1) (t220 * (t38 * pkin(4) + pkin(12) * t39) + t222 * t39 + t221 * t38) * g(3) + (t220 * (t23 * pkin(4) + pkin(12) * t24) + t222 * t24 + t221 * t23) * g(2) + (t220 * (-t27 * pkin(4) + pkin(12) * t28) + t222 * t28 - t221 * t27) * g(1) (t178 * (t124 * t77 + t127 * t39) + t179 * t15) * g(3) + (t179 * t1 + t178 * t234) * g(2) + (t178 * t6 + t179 * t5) * g(1) (-g(1) * t5 - g(2) * t1 - g(3) * t15) * m(7)];
taug  = t2(:);
