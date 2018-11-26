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

% Quelle: HybrDyn-Toolbox (ehem. IRT-Maple-Toolbox)
% Datum: 2018-11-23 18:44
% Revision: 76f9d5e39f14dc242b53c0d9d3d9db48bd8f37c0
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für mechatronische Systeme, Universität Hannover

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
assert( isreal(m) && all(size(m) == [7 1]), ...
  'S6RRRRRR8_gravloadJ_floatb_twist_slag_vp2: m has to be [7x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [7,3]), ...
  'S6RRRRRR8_gravloadJ_floatb_twist_slag_vp2: mrSges has to be [7x3] (double)');

%% Symbolic Calculation
% From joint_gravload_floatb_twist_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2018-11-23 18:43:44
% EndTime: 2018-11-23 18:43:46
% DurationCPUTime: 2.19s
% Computational Cost: add. (4737->203), mult. (4838->275), div. (0->0), fcn. (4776->24), ass. (0->105)
t124 = cos(qJ(4));
t121 = sin(qJ(4));
t118 = sin(pkin(7));
t119 = cos(pkin(7));
t127 = cos(qJ(1));
t189 = pkin(6) + qJ(2);
t167 = cos(t189) / 0.2e1;
t190 = pkin(6) - qJ(2);
t179 = cos(t190);
t151 = t179 / 0.2e1 + t167;
t210 = sin(qJ(2));
t211 = sin(qJ(1));
t145 = -t127 * t151 + t211 * t210;
t194 = sin(pkin(6));
t180 = t127 * t194;
t74 = t145 * t118 - t119 * t180;
t199 = t121 * t74;
t187 = pkin(7) + qJ(3);
t164 = sin(t187) / 0.2e1;
t188 = pkin(7) - qJ(3);
t176 = sin(t188);
t100 = t164 - t176 / 0.2e1;
t166 = cos(t187) / 0.2e1;
t178 = cos(t188);
t102 = t166 - t178 / 0.2e1;
t125 = cos(qJ(3));
t165 = sin(t189) / 0.2e1;
t177 = sin(t190);
t101 = t165 - t177 / 0.2e1;
t126 = cos(qJ(2));
t91 = t127 * t101 + t211 * t126;
t44 = t145 * t100 - t102 * t180 - t91 * t125;
t228 = -t124 * t44 + t199;
t222 = mrSges(6,2) - mrSges(7,3);
t120 = sin(qJ(6));
t123 = cos(qJ(6));
t217 = mrSges(7,1) * t123 - mrSges(7,2) * t120 + mrSges(6,1);
t117 = qJ(4) + qJ(5);
t114 = sin(t117);
t115 = cos(t117);
t14 = t114 * t74 - t115 * t44;
t13 = t114 * t44 + t115 * t74;
t219 = t121 * t44 + t124 * t74;
t223 = m(7) * pkin(5) + t217;
t175 = -m(7) * pkin(13) + t222;
t142 = t127 * t210 + t211 * t151;
t168 = t194 * t211;
t93 = -t211 * t101 + t127 * t126;
t129 = -t142 * t100 - t102 * t168 + t93 * t125;
t221 = t142 * t118 + t119 * t168;
t19 = -t129 * t121 + t124 * t221;
t103 = t167 - t179 / 0.2e1;
t149 = t165 + t177 / 0.2e1;
t195 = cos(pkin(6));
t138 = t149 * t100 - t195 * t102 - t103 * t125;
t90 = -t149 * t118 + t195 * t119;
t220 = -t121 * t138 + t124 * t90;
t212 = m(6) + m(7);
t186 = m(4) + m(5) + t212;
t216 = pkin(2) * t186 + mrSges(3,1);
t153 = -m(5) * pkin(3) - t124 * mrSges(5,1) + t121 * mrSges(5,2) - mrSges(4,1);
t163 = -m(5) * pkin(11) + mrSges(4,2) - mrSges(5,3) - mrSges(6,3);
t215 = -mrSges(7,1) * t120 - mrSges(7,2) * t123 + t163;
t122 = sin(qJ(3));
t148 = t164 + t176 / 0.2e1;
t146 = t148 * t194;
t150 = t178 / 0.2e1 + t166;
t40 = t91 * t122 + t127 * t146 + t145 * t150;
t213 = -t175 * t114 + t223 * t115 - t153;
t193 = t114 * t118;
t192 = t115 * t118;
t191 = t127 * pkin(1) + pkin(9) * t168;
t183 = t13 * pkin(5) + pkin(13) * t14;
t17 = t114 * t129 - t115 * t221;
t18 = t114 * t221 + t115 * t129;
t182 = -t17 * pkin(5) + pkin(13) * t18;
t30 = -t114 * t138 + t115 * t90;
t31 = t114 * t90 + t115 * t138;
t181 = t30 * pkin(5) + pkin(13) * t31;
t174 = t219 * pkin(4);
t173 = t19 * pkin(4);
t172 = t220 * pkin(4);
t169 = -t211 * pkin(1) + pkin(9) * t180;
t159 = -t217 * t13 + t14 * t222;
t158 = t217 * t17 + t18 * t222;
t157 = -t217 * t30 + t222 * t31;
t140 = -mrSges(3,2) + (t124 * mrSges(5,2) + mrSges(4,3) + (t212 * pkin(4) + mrSges(5,1)) * t121 + t186 * pkin(10)) * t118;
t137 = -t91 * pkin(2) - pkin(10) * t74 + t169;
t135 = t93 * pkin(2) + pkin(10) * t221 + t191;
t113 = pkin(4) * t124 + pkin(3);
t128 = -pkin(12) - pkin(11);
t131 = t221 * t121;
t45 = t122 * t93 + t142 * t150 - t211 * t146;
t132 = pkin(4) * t131 + t113 * t129 - t45 * t128 + t135;
t68 = t103 * t100 + t149 * t125;
t67 = -t103 * t150 + t149 * t122;
t60 = -t103 * t122 - t195 * t148 - t149 * t150;
t58 = -t100 * t93 - t142 * t125;
t57 = -t142 * t122 + t150 * t93;
t56 = -t100 * t91 - t145 * t125;
t55 = -t145 * t122 + t150 * t91;
t20 = t124 * t129 + t131;
t2 = t120 * t45 + t123 * t18;
t1 = -t120 * t18 + t123 * t45;
t3 = [(-t127 * mrSges(2,1) + t211 * mrSges(2,2) - m(3) * t191 - t93 * mrSges(3,1) + t142 * mrSges(3,2) - mrSges(3,3) * t168 - m(4) * t135 - t129 * mrSges(4,1) - t221 * mrSges(4,3) - m(5) * (pkin(3) * t129 + t135) - t20 * mrSges(5,1) - t19 * mrSges(5,2) - m(6) * t132 - t18 * mrSges(6,1) - m(7) * (t18 * pkin(5) + t132) - t2 * mrSges(7,1) - t1 * mrSges(7,2) + t175 * t17 + t163 * t45) * g(2) + (t211 * mrSges(2,1) + t127 * mrSges(2,2) - m(3) * t169 + t91 * mrSges(3,1) - t145 * mrSges(3,2) - mrSges(3,3) * t180 - m(4) * t137 - t44 * mrSges(4,1) + t74 * mrSges(4,3) - m(5) * (t44 * pkin(3) + t137) + t228 * mrSges(5,1) + t219 * mrSges(5,2) + t175 * t13 + t223 * t14 - t215 * t40 - t212 * (-pkin(4) * t199 + t44 * t113 + t128 * t40 + t137)) * g(1) (t175 * (t103 * t192 + t114 * t68) + t153 * t68 - t223 * (-t103 * t193 + t115 * t68) + t215 * t67 + t140 * t103 - t212 * (t68 * t113 - t67 * t128) - t216 * t149) * g(3) + (t175 * (t114 * t56 - t192 * t91) + t153 * t56 - t223 * (t115 * t56 + t193 * t91) + t215 * t55 - t140 * t91 - t212 * (t56 * t113 - t55 * t128) + t216 * t145) * g(2) + (t175 * (t114 * t58 - t192 * t93) + t153 * t58 - t223 * (t115 * t58 + t193 * t93) + t215 * t57 - t140 * t93 - t212 * (t58 * t113 - t57 * t128) + t216 * t142) * g(1) (-t212 * (-t60 * t113 - t128 * t138) + t215 * t138 + t213 * t60) * g(3) + (-t212 * (-t40 * t113 + t128 * t44) - t215 * t44 + t213 * t40) * g(2) + (-t212 * (-t45 * t113 - t128 * t129) + t215 * t129 + t213 * t45) * g(1) (-t220 * mrSges(5,1) - (-t121 * t90 - t124 * t138) * mrSges(5,2) - m(6) * t172 - m(7) * (t172 + t181) + t157) * g(3) + (-t219 * mrSges(5,1) + t228 * mrSges(5,2) - m(6) * t174 - m(7) * (t174 + t183) + t159) * g(2) + (-t19 * mrSges(5,1) + t20 * mrSges(5,2) - m(6) * t173 - m(7) * (t173 + t182) + t158) * g(1) (-m(7) * t181 + t157) * g(3) + (-m(7) * t183 + t159) * g(2) + (-m(7) * t182 + t158) * g(1), -g(1) * (mrSges(7,1) * t1 - mrSges(7,2) * t2) - g(2) * ((-t120 * t14 + t123 * t40) * mrSges(7,1) + (-t120 * t40 - t123 * t14) * mrSges(7,2)) - g(3) * ((-t120 * t31 + t123 * t60) * mrSges(7,1) + (-t120 * t60 - t123 * t31) * mrSges(7,2))];
taug  = t3(:);
