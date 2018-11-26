% Calculate Gravitation load on the joints for
% S6RRRRRP8
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% g [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [11x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,alpha2,d1,d2,d3,d4,d5]';
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
% Datum: 2018-11-23 18:32
% Revision: 76f9d5e39f14dc242b53c0d9d3d9db48bd8f37c0
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für mechatronische Systeme, Universität Hannover

function taug = S6RRRRRP8_gravloadJ_floatb_twist_slag_vp2(qJ, g, ...
  pkin, m, mrSges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(3,1),zeros(11,1),zeros(7,1),zeros(7,3)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRRRRP8_gravloadJ_floatb_twist_slag_vp2: qJ has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6RRRRRP8_gravloadJ_floatb_twist_slag_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [11 1]), ...
  'S6RRRRRP8_gravloadJ_floatb_twist_slag_vp2: pkin has to be [11x1] (double)');
assert( isreal(m) && all(size(m) == [7 1]), ...
  'S6RRRRRP8_gravloadJ_floatb_twist_slag_vp2: m has to be [7x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [7,3]), ...
  'S6RRRRRP8_gravloadJ_floatb_twist_slag_vp2: mrSges has to be [7x3] (double)');

%% Symbolic Calculation
% From joint_gravload_floatb_twist_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2018-11-23 18:31:50
% EndTime: 2018-11-23 18:31:51
% DurationCPUTime: 1.38s
% Computational Cost: add. (1990->177), mult. (2121->231), div. (0->0), fcn. (2137->16), ass. (0->93)
t111 = sin(qJ(5));
t202 = mrSges(6,2) - mrSges(7,3);
t204 = t111 * t202 - mrSges(5,1);
t196 = mrSges(6,3) + mrSges(7,2);
t201 = mrSges(5,2) - t196;
t203 = mrSges(6,1) + mrSges(7,1);
t112 = sin(qJ(3));
t116 = cos(qJ(3));
t114 = sin(qJ(1));
t117 = cos(qJ(2));
t187 = cos(qJ(1));
t159 = pkin(6) + qJ(2);
t101 = sin(t159) / 0.2e1;
t160 = pkin(6) - qJ(2);
t104 = sin(t160);
t89 = t101 - t104 / 0.2e1;
t127 = -t114 * t89 + t117 * t187;
t109 = sin(pkin(6));
t162 = t109 * t114;
t43 = -t112 * t127 + t116 * t162;
t110 = cos(pkin(6));
t137 = cos(t159) / 0.2e1;
t143 = cos(t160);
t90 = t137 - t143 / 0.2e1;
t199 = t110 * t116 + t112 * t90;
t115 = cos(qJ(5));
t108 = qJ(3) + qJ(4);
t105 = sin(t108);
t106 = cos(t108);
t128 = t114 * t117 + t187 * t89;
t149 = t109 * t187;
t36 = -t105 * t149 + t106 * t128;
t113 = sin(qJ(2));
t120 = t143 / 0.2e1 + t137;
t74 = t113 * t114 - t120 * t187;
t1 = t111 * t36 - t74 * t115;
t198 = t111 * t74 + t115 * t36;
t197 = m(6) + m(7);
t165 = qJ(6) * t111;
t35 = -t105 * t128 - t106 * t149;
t170 = t115 * t35;
t195 = pkin(5) * t170 + t35 * t165;
t39 = t105 * t127 - t106 * t162;
t169 = t115 * t39;
t194 = -pkin(5) * t169 - t39 * t165;
t63 = t105 * t90 + t106 * t110;
t168 = t115 * t63;
t193 = pkin(5) * t168 + t63 * t165;
t192 = -m(4) * pkin(9) + mrSges(3,2) - mrSges(4,3) - mrSges(5,3);
t191 = m(4) * pkin(2) + t116 * mrSges(4,1) + t106 * mrSges(5,1) - t112 * mrSges(4,2) - mrSges(5,2) * t105 + mrSges(3,1);
t64 = t105 * t110 - t106 * t90;
t190 = -t203 * t168 + t201 * t64 + t204 * t63;
t40 = t105 * t162 + t106 * t127;
t189 = t203 * t169 + t201 * t40 - t204 * t39;
t188 = -t203 * t170 + t201 * t36 + t204 * t35;
t142 = m(7) * pkin(5) + t203;
t140 = -m(7) * qJ(6) + t202;
t186 = pkin(4) * t106;
t103 = pkin(3) * t116 + pkin(2);
t118 = -pkin(10) - pkin(9);
t183 = -t74 * t103 - t118 * t128;
t77 = t113 * t187 + t114 * t120;
t182 = -t77 * t103 - t118 * t127;
t88 = t101 + t104 / 0.2e1;
t180 = t88 * t103 + t90 * t118;
t179 = t105 * t74;
t178 = t105 * t77;
t177 = t105 * t88;
t166 = t187 * pkin(1) + pkin(8) * t162;
t164 = t106 * t111;
t163 = t106 * t115;
t154 = t112 * t162;
t152 = t35 * pkin(4) + t36 * pkin(11);
t151 = -t39 * pkin(4) + pkin(11) * t40;
t150 = t63 * pkin(4) + pkin(11) * t64;
t148 = -pkin(1) * t114 + pkin(8) * t149;
t97 = t112 * t149;
t147 = -t116 * t128 + t97;
t141 = t43 * pkin(3);
t139 = t199 * pkin(3);
t138 = pkin(3) * t154 + t103 * t127 - t77 * t118 + t166;
t130 = pkin(3) * t97 - t103 * t128 + t74 * t118 + t148;
t129 = -t197 * pkin(11) + t201;
t124 = t141 + t151;
t123 = t112 * t128 + t116 * t149;
t122 = t139 + t150;
t121 = t123 * pkin(3);
t119 = -t121 + t152;
t44 = t116 * t127 + t154;
t11 = t111 * t64 + t88 * t115;
t6 = t111 * t77 + t115 * t40;
t5 = t111 * t40 - t77 * t115;
t2 = [(-t187 * mrSges(2,1) - m(3) * t166 - t127 * mrSges(3,1) - m(4) * (pkin(2) * t127 + t166) - t44 * mrSges(4,1) - t43 * mrSges(4,2) - m(5) * t138 - t40 * mrSges(5,1) - t142 * t6 + t140 * t5 + (-mrSges(3,3) * t109 + mrSges(2,2)) * t114 + t192 * t77 + t129 * t39 - t197 * (t40 * pkin(4) + t138)) * g(2) + (t114 * mrSges(2,1) + t187 * mrSges(2,2) - m(3) * t148 + t128 * mrSges(3,1) - mrSges(3,3) * t149 - m(4) * (-pkin(2) * t128 + t148) - t147 * mrSges(4,1) - t123 * mrSges(4,2) - m(5) * t130 + t36 * mrSges(5,1) + t142 * t198 - t140 * t1 - t192 * t74 + t129 * t35 + t197 * (pkin(4) * t36 - t130)) * g(1) (-m(5) * t180 - t197 * (pkin(11) * t177 + t88 * t186 + t180) - t142 * (-t111 * t90 + t163 * t88) + t140 * (t90 * t115 + t164 * t88) - t196 * t177 - t192 * t90 - t191 * t88) * g(3) + (-m(5) * t183 - t197 * (-pkin(11) * t179 - t74 * t186 + t183) - t142 * (t111 * t128 - t163 * t74) + t140 * (-t115 * t128 - t164 * t74) + t196 * t179 + t192 * t128 + t191 * t74) * g(2) + (-m(5) * t182 - t197 * (-pkin(11) * t178 - t77 * t186 + t182) + t140 * (-t115 * t127 - t164 * t77) + t196 * t178 - t142 * (t111 * t127 - t163 * t77) + t192 * t127 + t191 * t77) * g(1) (-t199 * mrSges(4,1) - (-t110 * t112 + t116 * t90) * mrSges(4,2) - m(5) * t139 - m(6) * t122 - m(7) * (t122 + t193) + t190) * g(3) + (t123 * mrSges(4,1) - t147 * mrSges(4,2) + m(5) * t121 - m(6) * t119 - m(7) * (t119 + t195) + t188) * g(2) + (-mrSges(4,1) * t43 + mrSges(4,2) * t44 - m(5) * t141 - m(6) * t124 - m(7) * (t124 + t194) + t189) * g(1) (-m(6) * t150 - m(7) * (t150 + t193) + t190) * g(3) + (-m(6) * t152 - m(7) * (t152 + t195) + t188) * g(2) + (-m(6) * t151 - m(7) * (t151 + t194) + t189) * g(1) (t140 * (-t111 * t88 + t115 * t64) + t142 * t11) * g(3) + (t142 * t1 + t140 * t198) * g(2) + (t140 * t6 + t142 * t5) * g(1) (-g(1) * t5 - g(2) * t1 - g(3) * t11) * m(7)];
taug  = t2(:);
