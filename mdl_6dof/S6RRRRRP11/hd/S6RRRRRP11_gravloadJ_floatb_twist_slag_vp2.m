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

% Quelle: HybrDyn-Toolbox (ehem. IRT-Maple-Toolbox)
% Datum: 2018-11-23 18:35
% Revision: 76f9d5e39f14dc242b53c0d9d3d9db48bd8f37c0
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für mechatronische Systeme, Universität Hannover

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
assert( isreal(m) && all(size(m) == [7 1]), ...
  'S6RRRRRP11_gravloadJ_floatb_twist_slag_vp2: m has to be [7x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [7,3]), ...
  'S6RRRRRP11_gravloadJ_floatb_twist_slag_vp2: mrSges has to be [7x3] (double)');

%% Symbolic Calculation
% From joint_gravload_floatb_twist_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2018-11-23 18:34:50
% EndTime: 2018-11-23 18:34:53
% DurationCPUTime: 2.65s
% Computational Cost: add. (5181->218), mult. (5482->298), div. (0->0), fcn. (5432->22), ass. (0->118)
t102 = sin(qJ(5));
t105 = cos(qJ(5));
t103 = sin(qJ(4));
t106 = cos(qJ(4));
t100 = sin(pkin(7));
t168 = pkin(6) + qJ(2);
t151 = cos(t168) / 0.2e1;
t169 = pkin(6) - qJ(2);
t159 = cos(t169);
t133 = t159 / 0.2e1 + t151;
t191 = sin(qJ(2));
t192 = sin(qJ(1));
t193 = cos(qJ(1));
t124 = -t133 * t193 + t192 * t191;
t175 = sin(pkin(6));
t176 = cos(pkin(7));
t146 = t176 * t175;
t116 = t124 * t100 - t193 * t146;
t167 = pkin(7) - qJ(3);
t158 = cos(t167);
t150 = t158 / 0.2e1;
t166 = pkin(7) + qJ(3);
t157 = cos(t166);
t144 = t150 - t157 / 0.2e1;
t137 = t144 * t175;
t107 = cos(qJ(3));
t147 = sin(t166) / 0.2e1;
t155 = sin(t167);
t201 = t147 - t155 / 0.2e1;
t108 = cos(qJ(2));
t148 = sin(t168) / 0.2e1;
t156 = sin(t169);
t87 = t148 - t156 / 0.2e1;
t78 = t108 * t192 + t193 * t87;
t205 = t78 * t107 - t124 * t201;
t41 = t137 * t193 - t205;
t20 = t103 * t116 - t41 * t106;
t104 = sin(qJ(3));
t128 = t147 + t155 / 0.2e1;
t125 = t128 * t175;
t149 = t157 / 0.2e1;
t132 = t150 + t149;
t37 = t78 * t104 + t124 * t132 + t125 * t193;
t216 = t102 * t20 - t105 * t37;
t215 = t102 * t37 + t105 * t20;
t214 = t103 * t41 + t106 * t116;
t213 = mrSges(6,1) + mrSges(7,1);
t202 = mrSges(6,2) + mrSges(7,2);
t98 = pkin(5) * t105 + pkin(4);
t208 = m(6) * pkin(4) + m(7) * t98 + mrSges(5,1);
t207 = m(6) * pkin(12) - m(7) * (-qJ(6) - pkin(12)) - mrSges(5,2) + mrSges(6,3) + mrSges(7,3);
t119 = t133 * t192 + t191 * t193;
t109 = t119 * t100 + t192 * t146;
t80 = t108 * t193 - t192 * t87;
t206 = t80 * t107 - t119 * t201;
t130 = t148 + t156 / 0.2e1;
t88 = t151 - t159 / 0.2e1;
t204 = -t88 * t107 + t130 * t201;
t194 = m(7) * pkin(5);
t203 = mrSges(4,2) - mrSges(5,3);
t170 = -m(5) - m(6) - m(7);
t200 = -t194 - t213;
t199 = -t202 * t102 + t213 * t105 + t208;
t197 = t207 * t103 + t208 * t106 + mrSges(4,1);
t190 = pkin(5) * t102;
t196 = -m(7) * (pkin(11) + t190) + t203;
t50 = -t107 * t124 - t201 * t78;
t74 = t124 * pkin(2);
t187 = t50 * pkin(3) - t74;
t52 = -t107 * t119 - t201 * t80;
t76 = t119 * pkin(2);
t186 = t52 * pkin(3) - t76;
t59 = t107 * t130 + t201 * t88;
t86 = t130 * pkin(2);
t185 = t59 * pkin(3) + t86;
t152 = t175 * t192;
t184 = t193 * pkin(1) + pkin(9) * t152;
t131 = t149 - t158 / 0.2e1;
t126 = t131 * t175;
t39 = t126 * t193 + t205;
t183 = t102 * t39;
t44 = -t126 * t192 + t206;
t182 = t102 * t44;
t177 = cos(pkin(6));
t56 = -t131 * t177 + t204;
t181 = t102 * t56;
t174 = t100 * t103;
t173 = t100 * t106;
t172 = t102 * t106;
t171 = t105 * t106;
t153 = t193 * t175;
t154 = -pkin(1) * t192 + pkin(9) * t153;
t43 = t137 * t192 + t206;
t24 = t103 * t109 + t43 * t106;
t42 = t104 * t80 + t119 * t132 - t125 * t192;
t5 = -t102 * t24 + t105 * t42;
t136 = -mrSges(3,2) + (mrSges(4,3) + (m(4) - t170) * pkin(10)) * t100;
t135 = -m(7) * t190 + pkin(11) * t170 + t203;
t121 = -t100 * t130 + t176 * t177;
t117 = -t78 * pkin(2) - pkin(10) * t116 + t154;
t114 = t41 * pkin(3) + t117;
t113 = t80 * pkin(2) + pkin(10) * t109 + t184;
t112 = t43 * pkin(3) + t113;
t111 = -pkin(11) * t37 + t114;
t110 = t42 * pkin(11) + t112;
t58 = t104 * t130 - t132 * t88;
t55 = t144 * t177 + t204;
t54 = -t104 * t88 - t128 * t177 - t130 * t132;
t51 = -t104 * t119 + t132 * t80;
t49 = -t104 * t124 + t132 * t78;
t46 = t106 * t59 - t174 * t88;
t32 = t103 * t121 + t55 * t106;
t31 = t103 * t55 - t106 * t121;
t30 = t106 * t52 + t174 * t80;
t28 = t106 * t50 + t174 * t78;
t23 = t103 * t43 - t106 * t109;
t6 = t102 * t42 + t105 * t24;
t1 = [(-t193 * mrSges(2,1) + t192 * mrSges(2,2) - m(3) * t184 - t80 * mrSges(3,1) + t119 * mrSges(3,2) - mrSges(3,3) * t152 - m(4) * t113 - t43 * mrSges(4,1) - t109 * mrSges(4,3) - m(5) * t110 - t24 * mrSges(5,1) - m(6) * (t24 * pkin(4) + t110) - m(7) * (t24 * t98 + t112) - t213 * t6 - t202 * t5 + t196 * t42 - t207 * t23) * g(2) + (t192 * mrSges(2,1) + t193 * mrSges(2,2) - m(3) * t154 + t78 * mrSges(3,1) - t124 * mrSges(3,2) - mrSges(3,3) * t153 - m(4) * t117 - t41 * mrSges(4,1) + t116 * mrSges(4,3) - m(5) * t111 + t20 * mrSges(5,1) - m(6) * (-pkin(4) * t20 + t111) - m(7) * (-t20 * t98 + t114) - t196 * t37 + t213 * t215 - t202 * t216 - t207 * t214) * g(1) (-t130 * mrSges(3,1) - m(4) * t86 - t59 * mrSges(4,1) - m(5) * t185 - t46 * mrSges(5,1) - m(6) * (pkin(4) * t46 + t185) - m(7) * (t46 * t98 + t185) + t136 * t88 - t213 * (t102 * t58 + t105 * t46) - t202 * (-t102 * t46 + t105 * t58) + t135 * t58 - t207 * (t103 * t59 + t173 * t88)) * g(3) + (t124 * mrSges(3,1) + m(4) * t74 - t50 * mrSges(4,1) - m(5) * t187 - t28 * mrSges(5,1) - m(6) * (pkin(4) * t28 + t187) - m(7) * (t28 * t98 + t187) - t136 * t78 - t213 * (t102 * t49 + t105 * t28) - t202 * (-t102 * t28 + t105 * t49) + t135 * t49 - t207 * (t103 * t50 - t173 * t78)) * g(2) + (t119 * mrSges(3,1) + m(4) * t76 - t52 * mrSges(4,1) - m(5) * t186 - t30 * mrSges(5,1) - m(6) * (pkin(4) * t30 + t186) - m(7) * (t30 * t98 + t186) - t136 * t80 - t213 * (t102 * t51 + t105 * t30) - t202 * (-t102 * t30 + t105 * t51) + t135 * t51 - t207 * (t103 * t52 - t173 * t80)) * g(1) (-t181 * t194 + t203 * t56 - t213 * (-t171 * t54 + t181) - t202 * (t105 * t56 + t172 * t54) + t170 * (-t54 * pkin(3) + pkin(11) * t56) + t197 * t54) * g(3) + (-t183 * t194 - t213 * (-t171 * t37 + t183) - t202 * (t105 * t39 + t172 * t37) + t203 * t39 + t170 * (-t37 * pkin(3) + pkin(11) * t39) + t197 * t37) * g(2) + (-t182 * t194 - t202 * (t105 * t44 + t172 * t42) + t203 * t44 + t170 * (-t42 * pkin(3) + pkin(11) * t44) - t213 * (-t171 * t42 + t182) + t197 * t42) * g(1) (t199 * t31 - t207 * t32) * g(3) + (-t199 * t214 - t20 * t207) * g(2) + (t199 * t23 - t207 * t24) * g(1) (-t202 * (-t102 * t54 - t105 * t32) + t200 * (-t102 * t32 + t105 * t54)) * g(3) + (-t200 * t216 + t202 * t215) * g(2) + (t200 * t5 + t202 * t6) * g(1) (-g(1) * t23 + g(2) * t214 - g(3) * t31) * m(7)];
taug  = t1(:);
