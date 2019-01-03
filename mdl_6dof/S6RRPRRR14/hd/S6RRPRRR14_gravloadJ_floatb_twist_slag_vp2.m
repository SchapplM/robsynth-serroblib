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
%   mass of all robot links (leg links until cut joint, platform)
% mrSges [7x3]
%  first moment of all robot links (mass times center of mass in body frames)
%  rows: links of the robot (starting with base)
%  columns: x-, y-, z-coordinates
% 
% Output:
% taug [6x1]
%   joint torques required to compensate gravitation load

% Quelle: HybrDyn-Toolbox (ehem. IRT-Maple-Toolbox)
% Datum: 2019-01-03 10:26
% Revision: 5fdbc45bcf2cc60deefd7ac2d71d743ed41bf7e4 (2018-12-21)
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
% StartTime: 2019-01-03 10:07:52
% EndTime: 2019-01-03 10:08:00
% DurationCPUTime: 2.46s
% Computational Cost: add. (2525->208), mult. (7155->311), div. (0->0), fcn. (9298->18), ass. (0->117)
t101 = sin(qJ(5));
t104 = cos(qJ(5));
t102 = sin(qJ(4));
t188 = cos(qJ(4));
t105 = cos(qJ(1));
t174 = sin(pkin(6));
t177 = cos(pkin(7));
t154 = t177 * t174;
t144 = t105 * t154;
t178 = cos(pkin(6));
t189 = cos(qJ(2));
t161 = t178 * t189;
t186 = sin(qJ(2));
t187 = sin(qJ(1));
t84 = -t105 * t161 + t186 * t187;
t99 = sin(pkin(7));
t135 = t84 * t99 - t144;
t176 = cos(pkin(8));
t175 = cos(pkin(14));
t152 = t175 * t174;
t146 = t99 * t152;
t155 = t177 * t175;
t173 = sin(pkin(14));
t160 = t178 * t186;
t85 = t105 * t160 + t187 * t189;
t54 = t105 * t146 + t155 * t84 + t173 * t85;
t98 = sin(pkin(8));
t197 = -t135 * t98 + t176 * t54;
t151 = t174 * t173;
t145 = t99 * t151;
t153 = t177 * t173;
t53 = -t105 * t145 - t153 * t84 + t175 * t85;
t22 = t197 * t102 - t53 * t188;
t38 = t135 * t176 + t54 * t98;
t4 = t101 * t38 - t104 * t22;
t212 = t101 * t22 + t104 * t38;
t204 = m(6) + m(7);
t211 = t204 * pkin(12);
t100 = sin(qJ(6));
t103 = cos(qJ(6));
t202 = m(7) * pkin(5) + t103 * mrSges(7,1) - t100 * mrSges(7,2) + mrSges(6,1);
t200 = -m(7) * pkin(13) + mrSges(6,2) - mrSges(7,3);
t210 = t102 * t53;
t131 = t105 * t186 + t161 * t187;
t125 = t131 * t175;
t86 = t105 * t189 - t187 * t160;
t111 = t125 * t177 - t146 * t187 + t173 * t86;
t206 = t131 * t99 + t187 * t154;
t207 = t111 * t98 + t176 * t206;
t156 = t174 * t186;
t147 = t99 * t156;
t138 = t177 * t152;
t78 = -t186 * t138 - t151 * t189;
t66 = t176 * t147 - t78 * t98;
t149 = -t100 * mrSges(7,1) - t103 * mrSges(7,2);
t203 = mrSges(5,2) - mrSges(6,3);
t191 = t149 + t203;
t205 = pkin(4) * t204 - t101 * t200 + t104 * t202 + mrSges(5,1);
t201 = -t99 * mrSges(4,3) + mrSges(3,2);
t199 = t111 * t176 - t206 * t98;
t169 = t99 * t178;
t113 = -t138 * t189 + t186 * t151 - t175 * t169;
t158 = t189 * t174;
t130 = t158 * t99 - t177 * t178;
t198 = t113 * t176 + t130 * t98;
t192 = t191 - t211;
t183 = t98 * t99;
t181 = pkin(2) * t158 + qJ(3) * t147;
t157 = t174 * t187;
t180 = t105 * pkin(1) + pkin(10) * t157;
t179 = qJ(3) * t99;
t171 = t98 * t188;
t168 = t99 * t176;
t165 = t99 * t171;
t164 = -t84 * pkin(2) + t179 * t85;
t163 = -t131 * pkin(2) + t179 * t86;
t162 = t105 * t174 * pkin(10) - pkin(1) * t187;
t159 = t176 * t188;
t148 = t203 - t211;
t143 = t98 * t147;
t137 = t177 * t151;
t79 = -t186 * t137 + t152 * t189;
t142 = t79 * pkin(3) + pkin(11) * t66 + t181;
t61 = -t155 * t85 + t173 * t84;
t45 = t168 * t85 - t61 * t98;
t124 = t131 * t173;
t63 = -t155 * t86 + t124;
t46 = t168 * t86 - t63 * t98;
t134 = -t85 * pkin(2) + qJ(3) * t144 - t179 * t84 + t162;
t62 = -t153 * t85 - t175 * t84;
t127 = t62 * pkin(3) + pkin(11) * t45 + t164;
t64 = -t153 * t86 - t125;
t126 = t64 * pkin(3) + pkin(11) * t46 + t163;
t122 = -t53 * pkin(3) - pkin(11) * t38 + t134;
t118 = t86 * pkin(2) + qJ(3) * t206 + t180;
t56 = -t124 * t177 + t145 * t187 + t175 * t86;
t108 = t56 * pkin(3) + pkin(11) * t207 + t118;
t24 = -t199 * t102 + t56 * t188;
t107 = t24 * pkin(4) + t108;
t73 = t137 * t189 + t186 * t152 + t173 * t169;
t50 = t113 * t98 - t130 * t176;
t44 = t79 * t188 + (t176 * t78 + t143) * t102;
t43 = t102 * t79 - t143 * t188 - t159 * t78;
t35 = -t198 * t102 + t73 * t188;
t34 = t102 * t73 + t198 * t188;
t30 = t64 * t188 + (t176 * t63 + t183 * t86) * t102;
t29 = t64 * t102 - t159 * t63 - t165 * t86;
t28 = t62 * t188 + (t176 * t61 + t183 * t85) * t102;
t27 = t62 * t102 - t159 * t61 - t165 * t85;
t23 = t102 * t56 + t199 * t188;
t19 = t197 * t188 + t210;
t14 = t101 * t50 + t104 * t35;
t8 = t101 * t207 + t24 * t104;
t7 = t101 * t24 - t104 * t207;
t2 = t100 * t23 + t103 * t8;
t1 = -t100 * t8 + t103 * t23;
t3 = [(-t105 * mrSges(2,1) + t187 * mrSges(2,2) - m(3) * t180 - t86 * mrSges(3,1) + t131 * mrSges(3,2) - mrSges(3,3) * t157 - m(4) * t118 - t56 * mrSges(4,1) + t111 * mrSges(4,2) - t206 * mrSges(4,3) - m(5) * t108 - t24 * mrSges(5,1) - t207 * mrSges(5,3) - m(6) * t107 - t8 * mrSges(6,1) - m(7) * (t8 * pkin(5) + t107) - t2 * mrSges(7,1) - t1 * mrSges(7,2) + t200 * t7 + t148 * t23) * g(2) + (t187 * mrSges(2,1) - m(3) * t162 + t85 * mrSges(3,1) - t84 * mrSges(3,2) - m(4) * t134 + t53 * mrSges(4,1) - t54 * mrSges(4,2) + t135 * mrSges(4,3) - m(5) * t122 - t22 * mrSges(5,1) + t38 * mrSges(5,3) + t200 * t212 + t202 * t4 + (t148 + t149) * (t135 * t171 - t159 * t54 - t210) + (-mrSges(3,3) * t174 + mrSges(2,2)) * t105 + t204 * (-t22 * pkin(4) - t122)) * g(1) (-m(4) * t181 - m(5) * t142 - mrSges(3,1) * t158 - t79 * mrSges(4,1) - t44 * mrSges(5,1) + mrSges(3,2) * t156 - t78 * mrSges(4,2) - mrSges(4,3) * t147 - t66 * mrSges(5,3) - t204 * (t44 * pkin(4) + pkin(12) * t43 + t142) - t202 * (t101 * t66 + t104 * t44) + t191 * t43 + t200 * (t101 * t44 - t66 * t104)) * g(3) + (-m(4) * t164 - m(5) * t127 + mrSges(3,1) * t84 - t62 * mrSges(4,1) - t28 * mrSges(5,1) - t61 * mrSges(4,2) - t45 * mrSges(5,3) - t204 * (t28 * pkin(4) + t27 * pkin(12) + t127) + t200 * (t101 * t28 - t45 * t104) + t201 * t85 - t202 * (t101 * t45 + t104 * t28) + t191 * t27) * g(2) + (-m(4) * t163 - m(5) * t126 + t131 * mrSges(3,1) - t64 * mrSges(4,1) - t30 * mrSges(5,1) - t63 * mrSges(4,2) - t46 * mrSges(5,3) + t201 * t86 - t204 * (t30 * pkin(4) + t29 * pkin(12) + t126) - t202 * (t101 * t46 + t104 * t30) + t191 * t29 + t200 * (t101 * t30 - t46 * t104)) * g(1) (-g(1) * t206 - g(2) * t135 + g(3) * t130) * (m(4) + m(5) + t204) (t192 * t35 + t205 * t34) * g(3) + (t19 * t205 - t192 * t22) * g(2) + (t192 * t24 + t205 * t23) * g(1) (t200 * t14 - t202 * (-t101 * t35 + t104 * t50)) * g(3) + (t200 * t4 - t202 * t212) * g(2) + (t200 * t8 + t202 * t7) * g(1), -g(1) * (mrSges(7,1) * t1 - mrSges(7,2) * t2) - g(2) * ((-t100 * t4 + t103 * t19) * mrSges(7,1) + (-t100 * t19 - t103 * t4) * mrSges(7,2)) - g(3) * ((-t100 * t14 + t103 * t34) * mrSges(7,1) + (-t100 * t34 - t103 * t14) * mrSges(7,2))];
taug  = t3(:);
