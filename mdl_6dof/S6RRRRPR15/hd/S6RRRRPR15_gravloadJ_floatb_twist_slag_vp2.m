% Calculate Gravitation load on the joints for
% S6RRRRPR15
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% g [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [12x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,alpha2,alpha3,d1,d2,d3,d4,d6]';
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
% Datum: 2018-11-23 18:26
% Revision: 76f9d5e39f14dc242b53c0d9d3d9db48bd8f37c0
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für mechatronische Systeme, Universität Hannover

function taug = S6RRRRPR15_gravloadJ_floatb_twist_slag_vp2(qJ, g, ...
  pkin, m, mrSges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(3,1),zeros(12,1),zeros(7,1),zeros(7,3)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRRRPR15_gravloadJ_floatb_twist_slag_vp2: qJ has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6RRRRPR15_gravloadJ_floatb_twist_slag_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [12 1]), ...
  'S6RRRRPR15_gravloadJ_floatb_twist_slag_vp2: pkin has to be [12x1] (double)');
assert( isreal(m) && all(size(m) == [7 1]), ...
  'S6RRRRPR15_gravloadJ_floatb_twist_slag_vp2: m has to be [7x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [7,3]), ...
  'S6RRRRPR15_gravloadJ_floatb_twist_slag_vp2: mrSges has to be [7x3] (double)');

%% Symbolic Calculation
% From joint_gravload_floatb_twist_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2018-11-23 18:25:06
% EndTime: 2018-11-23 18:25:08
% DurationCPUTime: 2.26s
% Computational Cost: add. (4760->205), mult. (5033->269), div. (0->0), fcn. (4968->22), ass. (0->122)
t102 = cos(qJ(6));
t205 = mrSges(5,2) - mrSges(6,3);
t99 = sin(qJ(6));
t221 = t99 * mrSges(7,1) + t102 * mrSges(7,2) - t205;
t207 = m(6) + m(7);
t220 = -qJ(5) * t207 - t221;
t100 = sin(qJ(4));
t103 = cos(qJ(4));
t106 = cos(qJ(1));
t167 = pkin(6) + qJ(2);
t150 = cos(t167) / 0.2e1;
t168 = pkin(6) - qJ(2);
t157 = cos(t168);
t132 = t157 / 0.2e1 + t150;
t188 = sin(qJ(2));
t189 = sin(qJ(1));
t123 = -t106 * t132 + t189 * t188;
t172 = sin(pkin(6));
t173 = cos(pkin(7));
t143 = t173 * t172;
t98 = sin(pkin(7));
t115 = -t106 * t143 + t123 * t98;
t166 = pkin(7) - qJ(3);
t156 = cos(t166);
t149 = t156 / 0.2e1;
t165 = pkin(7) + qJ(3);
t155 = cos(t165);
t140 = t149 - t155 / 0.2e1;
t134 = t140 * t172;
t104 = cos(qJ(3));
t146 = sin(t165) / 0.2e1;
t153 = sin(t166);
t201 = t146 - t153 / 0.2e1;
t105 = cos(qJ(2));
t147 = sin(t167) / 0.2e1;
t154 = sin(t168);
t86 = t147 - t154 / 0.2e1;
t76 = t189 * t105 + t106 * t86;
t209 = t76 * t104 - t123 * t201;
t34 = t106 * t134 - t209;
t219 = t115 * t100 - t34 * t103;
t9 = t100 * t34 + t115 * t103;
t218 = mrSges(5,1) - mrSges(6,2);
t119 = t106 * t188 + t189 * t132;
t108 = t119 * t98 + t189 * t143;
t198 = m(7) * pkin(12) + mrSges(7,3) + t218;
t211 = pkin(4) * t207 + t198;
t78 = t106 * t105 - t189 * t86;
t210 = t78 * t104 - t119 * t201;
t129 = t147 + t154 / 0.2e1;
t87 = t150 - t157 / 0.2e1;
t208 = -t87 * t104 + t129 * t201;
t206 = m(7) * (pkin(5) + pkin(11));
t171 = qJ(5) * t100;
t101 = sin(qJ(3));
t127 = t146 + t153 / 0.2e1;
t124 = t127 * t172;
t148 = t155 / 0.2e1;
t131 = t149 + t148;
t30 = t76 * t101 + t106 * t124 + t123 * t131;
t181 = t103 * t30;
t204 = -pkin(4) * t181 - t30 * t171;
t35 = t101 * t78 + t119 * t131 - t189 * t124;
t180 = t103 * t35;
t203 = -pkin(4) * t180 - t35 * t171;
t174 = cos(pkin(6));
t50 = -t101 * t87 - t174 * t127 - t129 * t131;
t179 = t103 * t50;
t202 = -pkin(4) * t179 - t50 * t171;
t200 = mrSges(4,2) - mrSges(6,1) - mrSges(5,3);
t196 = -t102 * mrSges(7,1) + t99 * mrSges(7,2) + t200;
t195 = t221 * t100 + t218 * t103 + mrSges(4,1);
t194 = t196 - t206;
t191 = t30 * pkin(11);
t190 = t35 * pkin(11);
t45 = -t123 * t104 - t201 * t76;
t72 = t123 * pkin(2);
t187 = t45 * pkin(3) - t72;
t47 = -t119 * t104 - t201 * t78;
t74 = t119 * pkin(2);
t186 = t47 * pkin(3) - t74;
t55 = t129 * t104 + t201 * t87;
t85 = t129 * pkin(2);
t185 = t55 * pkin(3) + t85;
t151 = t172 * t189;
t184 = t106 * pkin(1) + pkin(9) * t151;
t183 = t100 * t98;
t178 = t103 * t98;
t169 = -m(5) - t207;
t26 = t30 * pkin(3);
t130 = t148 - t156 / 0.2e1;
t125 = t130 * t172;
t32 = t106 * t125 + t209;
t161 = pkin(11) * t32 - t26;
t28 = t35 * pkin(3);
t37 = -t189 * t125 + t210;
t160 = pkin(11) * t37 - t28;
t49 = t50 * pkin(3);
t52 = -t174 * t130 + t208;
t159 = pkin(11) * t52 - t49;
t158 = t106 * t172;
t152 = -t189 * pkin(1) + pkin(9) * t158;
t133 = -mrSges(3,2) + (mrSges(4,3) + (m(4) - t169) * pkin(10)) * t98;
t121 = -m(7) * pkin(5) + t169 * pkin(11) + t196;
t117 = -t129 * t98 + t174 * t173;
t113 = -t76 * pkin(2) - t115 * pkin(10) + t152;
t112 = t34 * pkin(3) + t113;
t111 = t78 * pkin(2) + t108 * pkin(10) + t184;
t36 = t189 * t134 + t210;
t110 = t36 * pkin(3) + t111;
t109 = -pkin(4) * t219 + t9 * qJ(5) + t112;
t11 = t100 * t36 - t108 * t103;
t12 = t108 * t100 + t36 * t103;
t107 = t12 * pkin(4) + t11 * qJ(5) + t110;
t51 = t174 * t140 + t208;
t40 = t103 * t55 - t87 * t183;
t20 = t100 * t51 - t117 * t103;
t18 = t103 * t47 + t183 * t78;
t16 = t103 * t45 + t183 * t76;
t2 = t102 * t35 + t11 * t99;
t1 = t102 * t11 - t35 * t99;
t3 = [(-t106 * mrSges(2,1) + t189 * mrSges(2,2) - m(3) * t184 - t78 * mrSges(3,1) + t119 * mrSges(3,2) - mrSges(3,3) * t151 - m(4) * t111 - t36 * mrSges(4,1) - t108 * mrSges(4,3) - m(5) * (t110 + t190) - m(6) * (t107 + t190) - m(7) * t107 - t2 * mrSges(7,1) - t1 * mrSges(7,2) + t205 * t11 + (t200 - t206) * t35 - t198 * t12) * g(2) + (t189 * mrSges(2,1) + t106 * mrSges(2,2) - m(3) * t152 + t76 * mrSges(3,1) - t123 * mrSges(3,2) - mrSges(3,3) * t158 - m(4) * t113 - t34 * mrSges(4,1) + t115 * mrSges(4,3) - m(5) * (t112 - t191) - m(6) * (t109 - t191) - m(7) * t109 - t221 * t9 - t194 * t30 + t198 * t219) * g(1) (-t129 * mrSges(3,1) - m(4) * t85 - t55 * mrSges(4,1) - m(5) * t185 + t133 * t87 - t198 * t40 + t220 * (t100 * t55 + t87 * t178) + t121 * (t129 * t101 - t87 * t131) - t207 * (t40 * pkin(4) + t185)) * g(3) + (t123 * mrSges(3,1) + m(4) * t72 - t45 * mrSges(4,1) - m(5) * t187 - t133 * t76 - t198 * t16 + t220 * (t100 * t45 - t178 * t76) + t121 * (-t123 * t101 + t131 * t76) - t207 * (t16 * pkin(4) + t187)) * g(2) + (t119 * mrSges(3,1) + m(4) * t74 - t47 * mrSges(4,1) - m(5) * t186 - t133 * t78 - t198 * t18 + t220 * (t100 * t47 - t178 * t78) + t121 * (-t119 * t101 + t131 * t78) - t207 * (t18 * pkin(4) + t186)) * g(1) (-m(5) * t159 - m(6) * (t159 + t202) - m(7) * (-pkin(12) * t179 + t202 - t49) + mrSges(7,3) * t179 + t194 * t52 + t195 * t50) * g(3) + (-m(5) * t161 - m(6) * (t161 + t204) - m(7) * (-pkin(12) * t181 + t204 - t26) + mrSges(7,3) * t181 + t194 * t32 + t195 * t30) * g(2) + (-m(5) * t160 - m(6) * (t160 + t203) - m(7) * (-pkin(12) * t180 + t203 - t28) + mrSges(7,3) * t180 + t194 * t37 + t195 * t35) * g(1) (t220 * (t117 * t100 + t51 * t103) + t211 * t20) * g(3) + (-t211 * t9 + t219 * t220) * g(2) + (t211 * t11 + t12 * t220) * g(1), t207 * (-g(1) * t11 + g(2) * t9 - g(3) * t20) -g(1) * (mrSges(7,1) * t1 - mrSges(7,2) * t2) - g(2) * ((-t102 * t9 - t30 * t99) * mrSges(7,1) + (-t102 * t30 + t9 * t99) * mrSges(7,2)) - g(3) * ((t102 * t20 - t50 * t99) * mrSges(7,1) + (-t102 * t50 - t20 * t99) * mrSges(7,2))];
taug  = t3(:);
