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

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-10 06:28
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

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
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6RRRRRR10_gravloadJ_floatb_twist_slag_vp2: m has to be [7x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [7,3]), ...
  'S6RRRRRR10_gravloadJ_floatb_twist_slag_vp2: mrSges has to be [7x3] (double)');

%% Symbolic Calculation
% From gravload_joint_floatb_twist_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-10 05:49:45
% EndTime: 2019-03-10 05:49:55
% DurationCPUTime: 3.44s
% Computational Cost: add. (3117->248), mult. (8880->376), div. (0->0), fcn. (11549->18), ass. (0->125)
t123 = sin(qJ(5));
t127 = cos(qJ(5));
t124 = sin(qJ(4));
t218 = cos(qJ(4));
t120 = sin(pkin(8));
t208 = cos(pkin(8));
t210 = cos(pkin(6));
t216 = sin(qJ(2));
t181 = t210 * t216;
t217 = sin(qJ(1));
t219 = cos(qJ(2));
t220 = cos(qJ(1));
t106 = t181 * t220 + t217 * t219;
t125 = sin(qJ(3));
t128 = cos(qJ(3));
t182 = t210 * t219;
t105 = -t182 * t220 + t217 * t216;
t121 = sin(pkin(7));
t207 = sin(pkin(6));
t179 = t220 * t207;
t209 = cos(pkin(7));
t236 = t209 * t105 + t121 * t179;
t226 = t106 * t125 + t128 * t236;
t175 = t209 * t207;
t237 = -t105 * t121 + t220 * t175;
t229 = t120 * t237 + t208 * t226;
t72 = t106 * t128 - t236 * t125;
t24 = -t124 * t229 + t218 * t72;
t51 = -t226 * t120 + t208 * t237;
t242 = t123 * t51 - t127 * t24;
t241 = -t123 * t24 - t51 * t127;
t234 = m(6) + m(7);
t122 = sin(qJ(6));
t126 = cos(qJ(6));
t232 = m(7) * pkin(5) + t126 * mrSges(7,1) - t122 * mrSges(7,2) + mrSges(6,1);
t225 = -m(7) * pkin(14) + mrSges(6,2) - mrSges(7,3);
t238 = t124 * t72;
t107 = -t217 * t181 + t219 * t220;
t153 = t217 * t182 + t216 * t220;
t177 = t207 * t217;
t138 = t121 * t177 - t153 * t209;
t134 = t107 * t125 - t128 * t138;
t139 = t153 * t121 + t217 * t175;
t129 = t134 * t120 + t139 * t208;
t176 = t207 * t216;
t173 = t121 * t176;
t168 = t216 * t175;
t178 = t219 * t207;
t99 = -t125 * t178 - t128 * t168;
t85 = -t99 * t120 + t208 * t173;
t233 = mrSges(5,2) - mrSges(6,3);
t222 = -t122 * mrSges(7,1) - t126 * mrSges(7,2) + t233;
t235 = pkin(4) * t234 - t225 * t123 + t232 * t127 + mrSges(5,1);
t231 = -t120 * t139 + t134 * t208;
t156 = t121 * t210 + t175 * t219;
t140 = t125 * t176 - t128 * t156;
t155 = -t121 * t178 + t209 * t210;
t230 = -t155 * t120 + t140 * t208;
t224 = -t234 * pkin(13) + t222;
t215 = t120 * t72;
t75 = t107 * t128 + t125 * t138;
t214 = t120 * t75;
t94 = t125 * t156 + t128 * t176;
t213 = t120 * t94;
t205 = t106 * t121;
t203 = t107 * t121;
t202 = t120 * t121;
t201 = t120 * t123;
t200 = t120 * t127;
t199 = pkin(2) * t178 + pkin(11) * t173;
t198 = t220 * pkin(1) + pkin(10) * t177;
t196 = t120 * t218;
t194 = t121 * t208;
t193 = t124 * t208;
t192 = t125 * t209;
t191 = t128 * t209;
t189 = -pkin(3) * t226 + pkin(12) * t215;
t188 = -t134 * pkin(3) + pkin(12) * t214;
t187 = -t140 * pkin(3) + pkin(12) * t213;
t186 = t121 * t196;
t185 = -t105 * pkin(2) + pkin(11) * t205;
t184 = -t153 * pkin(2) + pkin(11) * t203;
t183 = -pkin(1) * t217 + pkin(10) * t179;
t180 = t208 * t218;
t171 = t120 * t173;
t100 = -t125 * t168 + t128 * t178;
t162 = t100 * pkin(3) + t85 * pkin(12) + t199;
t80 = t105 * t125 - t106 * t191;
t60 = t106 * t194 - t80 * t120;
t82 = -t107 * t191 + t125 * t153;
t61 = t107 * t194 - t82 * t120;
t160 = -t106 * pkin(2) + t237 * pkin(11) + t183;
t81 = -t105 * t128 - t106 * t192;
t149 = t81 * pkin(3) + pkin(12) * t60 + t185;
t83 = -t107 * t192 - t128 * t153;
t148 = t83 * pkin(3) + pkin(12) * t61 + t184;
t147 = -t72 * pkin(3) + t51 * pkin(12) + t160;
t141 = t107 * pkin(2) + t139 * pkin(11) + t198;
t131 = t75 * pkin(3) + t129 * pkin(12) + t141;
t27 = t124 * t75 + t218 * t231;
t28 = -t124 * t231 + t218 * t75;
t130 = t28 * pkin(4) + t27 * pkin(13) + t131;
t67 = t120 * t140 + t155 * t208;
t59 = t100 * t218 + (t208 * t99 + t171) * t124;
t58 = t100 * t124 - t171 * t218 - t180 * t99;
t54 = -t140 * t218 - t193 * t94;
t53 = -t124 * t140 + t180 * t94;
t47 = -t124 * t230 + t94 * t218;
t46 = t124 * t94 + t218 * t230;
t40 = -t134 * t218 - t193 * t75;
t39 = -t124 * t134 + t180 * t75;
t38 = -t193 * t72 - t218 * t226;
t37 = -t124 * t226 + t180 * t72;
t36 = t83 * t218 + (t107 * t202 + t208 * t82) * t124;
t35 = -t107 * t186 + t83 * t124 - t180 * t82;
t34 = t81 * t218 + (t106 * t202 + t208 * t80) * t124;
t33 = -t106 * t186 + t81 * t124 - t180 * t80;
t25 = -t180 * t226 - t196 * t237 - t238;
t23 = t218 * t229 + t238;
t18 = t123 * t67 + t127 * t47;
t8 = t123 * t129 + t28 * t127;
t7 = t123 * t28 - t127 * t129;
t2 = t122 * t27 + t126 * t8;
t1 = -t122 * t8 + t126 * t27;
t3 = [(-t220 * mrSges(2,1) + t217 * mrSges(2,2) - m(3) * t198 - t107 * mrSges(3,1) + t153 * mrSges(3,2) - mrSges(3,3) * t177 - m(4) * t141 - t75 * mrSges(4,1) + t134 * mrSges(4,2) - t139 * mrSges(4,3) - m(5) * t131 - t28 * mrSges(5,1) - t129 * mrSges(5,3) - m(6) * t130 - t8 * mrSges(6,1) - m(7) * (t8 * pkin(5) + t130) - t2 * mrSges(7,1) - t1 * mrSges(7,2) + t225 * t7 + t233 * t27) * g(2) + (-m(3) * t183 - m(4) * t160 - m(5) * t147 + t217 * mrSges(2,1) + t106 * mrSges(3,1) + t72 * mrSges(4,1) + t24 * mrSges(5,1) + t220 * mrSges(2,2) - t105 * mrSges(3,2) - t226 * mrSges(4,2) - mrSges(3,3) * t179 - t237 * mrSges(4,3) - t51 * mrSges(5,3) - t234 * (-pkin(4) * t24 + t25 * pkin(13) + t147) + t225 * t241 - t232 * t242 + t222 * t25) * g(1) (-m(4) * t199 - m(5) * t162 - mrSges(3,1) * t178 - t100 * mrSges(4,1) - t59 * mrSges(5,1) + mrSges(3,2) * t176 - t99 * mrSges(4,2) - mrSges(4,3) * t173 - t85 * mrSges(5,3) - t234 * (t59 * pkin(4) + pkin(13) * t58 + t162) - t232 * (t123 * t85 + t127 * t59) + t222 * t58 + t225 * (t123 * t59 - t85 * t127)) * g(3) + (-m(4) * t185 - m(5) * t149 + mrSges(3,1) * t105 - t81 * mrSges(4,1) - t34 * mrSges(5,1) + mrSges(3,2) * t106 - t80 * mrSges(4,2) - mrSges(4,3) * t205 - t60 * mrSges(5,3) - t234 * (t34 * pkin(4) + t33 * pkin(13) + t149) + t225 * (t123 * t34 - t60 * t127) - t232 * (t123 * t60 + t127 * t34) + t222 * t33) * g(2) + (-m(4) * t184 - m(5) * t148 + mrSges(3,1) * t153 - t83 * mrSges(4,1) - t36 * mrSges(5,1) + t107 * mrSges(3,2) - t82 * mrSges(4,2) - mrSges(4,3) * t203 - t61 * mrSges(5,3) - t234 * (t36 * pkin(4) + t35 * pkin(13) + t148) - t232 * (t123 * t61 + t127 * t36) + t222 * t35 + t225 * (t123 * t36 - t61 * t127)) * g(1) (-m(5) * t187 + mrSges(4,1) * t140 - t54 * mrSges(5,1) + t94 * mrSges(4,2) - mrSges(5,3) * t213 - t234 * (t54 * pkin(4) + pkin(13) * t53 + t187) - t232 * (t127 * t54 + t201 * t94) + t222 * t53 + t225 * (t123 * t54 - t200 * t94)) * g(3) + (-m(5) * t189 + mrSges(4,1) * t226 - t38 * mrSges(5,1) + t72 * mrSges(4,2) - mrSges(5,3) * t215 - t234 * (t38 * pkin(4) + pkin(13) * t37 + t189) - t232 * (t127 * t38 + t201 * t72) + t222 * t37 + t225 * (t123 * t38 - t200 * t72)) * g(2) + (-m(5) * t188 + mrSges(4,1) * t134 - t40 * mrSges(5,1) + t75 * mrSges(4,2) - mrSges(5,3) * t214 - t234 * (t40 * pkin(4) + pkin(13) * t39 + t188) - t232 * (t127 * t40 + t201 * t75) + t222 * t39 + t225 * (t123 * t40 - t200 * t75)) * g(1) (t224 * t47 + t235 * t46) * g(3) + (t224 * t24 + t235 * t23) * g(2) + (t224 * t28 + t235 * t27) * g(1) (t225 * t18 - t232 * (-t123 * t47 + t127 * t67)) * g(3) + (-t225 * t242 - t232 * t241) * g(2) + (t225 * t8 + t232 * t7) * g(1), -g(1) * (mrSges(7,1) * t1 - mrSges(7,2) * t2) - g(2) * ((t122 * t242 + t126 * t23) * mrSges(7,1) + (-t122 * t23 + t126 * t242) * mrSges(7,2)) - g(3) * ((-t122 * t18 + t126 * t46) * mrSges(7,1) + (-t122 * t46 - t126 * t18) * mrSges(7,2))];
taug  = t3(:);
