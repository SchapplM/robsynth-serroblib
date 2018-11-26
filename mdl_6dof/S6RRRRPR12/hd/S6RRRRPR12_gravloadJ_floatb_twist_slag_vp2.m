% Calculate Gravitation load on the joints for
% S6RRRRPR12
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% g [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [13x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,alpha2,alpha3,d1,d2,d3,d4,d6,theta5]';
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
% Datum: 2018-11-23 18:22
% Revision: 76f9d5e39f14dc242b53c0d9d3d9db48bd8f37c0
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für mechatronische Systeme, Universität Hannover

function taug = S6RRRRPR12_gravloadJ_floatb_twist_slag_vp2(qJ, g, ...
  pkin, m, mrSges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(3,1),zeros(13,1),zeros(7,1),zeros(7,3)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRRRPR12_gravloadJ_floatb_twist_slag_vp2: qJ has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6RRRRPR12_gravloadJ_floatb_twist_slag_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [13 1]), ...
  'S6RRRRPR12_gravloadJ_floatb_twist_slag_vp2: pkin has to be [13x1] (double)');
assert( isreal(m) && all(size(m) == [7 1]), ...
  'S6RRRRPR12_gravloadJ_floatb_twist_slag_vp2: m has to be [7x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [7,3]), ...
  'S6RRRRPR12_gravloadJ_floatb_twist_slag_vp2: mrSges has to be [7x3] (double)');

%% Symbolic Calculation
% From joint_gravload_floatb_twist_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2018-11-23 18:21:32
% EndTime: 2018-11-23 18:21:34
% DurationCPUTime: 2.11s
% Computational Cost: add. (4331->191), mult. (4431->262), div. (0->0), fcn. (4347->24), ass. (0->94)
t111 = cos(qJ(4));
t108 = sin(qJ(4));
t104 = sin(pkin(7));
t105 = cos(pkin(7));
t114 = cos(qJ(1));
t167 = pkin(6) + qJ(2);
t148 = cos(t167) / 0.2e1;
t168 = pkin(6) - qJ(2);
t160 = cos(t168);
t134 = t160 / 0.2e1 + t148;
t188 = sin(qJ(2));
t189 = sin(qJ(1));
t128 = -t114 * t134 + t189 * t188;
t171 = sin(pkin(6));
t161 = t114 * t171;
t59 = t128 * t104 - t105 * t161;
t178 = t108 * t59;
t112 = cos(qJ(3));
t113 = cos(qJ(2));
t146 = sin(t167) / 0.2e1;
t158 = sin(t168);
t87 = t146 - t158 / 0.2e1;
t77 = t113 * t189 + t114 * t87;
t165 = pkin(7) + qJ(3);
t145 = sin(t165) / 0.2e1;
t166 = pkin(7) - qJ(3);
t157 = sin(t166);
t86 = t145 - t157 / 0.2e1;
t147 = cos(t165) / 0.2e1;
t159 = cos(t166);
t88 = t147 - t159 / 0.2e1;
t29 = -t77 * t112 + t128 * t86 - t161 * t88;
t206 = -t111 * t29 + t178;
t103 = qJ(4) + pkin(13);
t100 = sin(t103);
t101 = cos(t103);
t4 = t100 * t59 - t101 * t29;
t205 = t100 * t29 + t101 * t59;
t198 = t108 * t29 + t111 * t59;
t107 = sin(qJ(6));
t110 = cos(qJ(6));
t196 = m(7) * pkin(5) + mrSges(7,1) * t110 - mrSges(7,2) * t107 + mrSges(6,1);
t195 = -m(7) * pkin(12) + mrSges(6,2) - mrSges(7,3);
t190 = m(6) + m(7);
t200 = pkin(4) * t190;
t125 = t114 * t188 + t134 * t189;
t149 = t171 * t189;
t60 = t125 * t104 + t105 * t149;
t132 = t146 + t158 / 0.2e1;
t172 = cos(pkin(6));
t89 = t148 - t160 / 0.2e1;
t122 = -t89 * t112 + t132 * t86 - t172 * t88;
t76 = -t104 * t132 + t105 * t172;
t199 = -t108 * t122 + t111 * t76;
t79 = t114 * t113 - t189 * t87;
t115 = t79 * t112 - t125 * t86 - t149 * t88;
t9 = -t108 * t115 + t111 * t60;
t164 = m(4) + m(5) + t190;
t197 = pkin(2) * t164 + mrSges(3,1);
t136 = -m(5) * pkin(3) - t111 * mrSges(5,1) + t108 * mrSges(5,2) - mrSges(4,1);
t144 = -m(5) * pkin(11) + mrSges(4,2) - mrSges(5,3) - mrSges(6,3);
t194 = -mrSges(7,1) * t107 - mrSges(7,2) * t110 + t144;
t109 = sin(qJ(3));
t131 = t145 + t157 / 0.2e1;
t129 = t131 * t171;
t133 = t159 / 0.2e1 + t147;
t25 = t77 * t109 + t114 * t129 + t128 * t133;
t191 = -t195 * t100 + t196 * t101 - t136;
t177 = t108 * t60;
t173 = t114 * pkin(1) + pkin(9) * t149;
t170 = t100 * t104;
t169 = t101 * t104;
t151 = -pkin(1) * t189 + pkin(9) * t161;
t123 = -mrSges(3,2) + (t111 * mrSges(5,2) + mrSges(4,3) + (mrSges(5,1) + t200) * t108 + t164 * pkin(10)) * t104;
t120 = -t77 * pkin(2) - pkin(10) * t59 + t151;
t118 = t79 * pkin(2) + pkin(10) * t60 + t173;
t106 = -qJ(5) - pkin(11);
t30 = t109 * t79 + t125 * t133 - t129 * t189;
t99 = pkin(4) * t111 + pkin(3);
t116 = pkin(4) * t177 - t30 * t106 + t115 * t99 + t118;
t53 = t112 * t132 + t89 * t86;
t52 = t109 * t132 - t133 * t89;
t45 = -t109 * t89 - t131 * t172 - t132 * t133;
t43 = -t112 * t125 - t79 * t86;
t42 = -t109 * t125 + t133 * t79;
t41 = -t112 * t128 - t77 * t86;
t40 = -t109 * t128 + t133 * t77;
t16 = t100 * t76 + t101 * t122;
t10 = t111 * t115 + t177;
t8 = t100 * t60 + t101 * t115;
t7 = t100 * t115 - t60 * t101;
t2 = t107 * t30 + t110 * t8;
t1 = -t107 * t8 + t110 * t30;
t3 = [(-t114 * mrSges(2,1) + t189 * mrSges(2,2) - m(3) * t173 - t79 * mrSges(3,1) + t125 * mrSges(3,2) - mrSges(3,3) * t149 - m(4) * t118 - t115 * mrSges(4,1) - t60 * mrSges(4,3) - m(5) * (pkin(3) * t115 + t118) - t10 * mrSges(5,1) - t9 * mrSges(5,2) - m(6) * t116 - t8 * mrSges(6,1) - m(7) * (t8 * pkin(5) + t116) - t2 * mrSges(7,1) - t1 * mrSges(7,2) + t195 * t7 + t144 * t30) * g(2) + (t189 * mrSges(2,1) + t114 * mrSges(2,2) - m(3) * t151 + t77 * mrSges(3,1) - t128 * mrSges(3,2) - mrSges(3,3) * t161 - m(4) * t120 - t29 * mrSges(4,1) + t59 * mrSges(4,3) - m(5) * (t29 * pkin(3) + t120) + t206 * mrSges(5,1) + t198 * mrSges(5,2) + t195 * t205 + t196 * t4 - t194 * t25 - t190 * (-pkin(4) * t178 + t106 * t25 + t29 * t99 + t120)) * g(1) (t195 * (t100 * t53 + t169 * t89) + t136 * t53 - t196 * (t101 * t53 - t170 * t89) + t194 * t52 + t123 * t89 - t190 * (-t52 * t106 + t53 * t99) - t197 * t132) * g(3) + (t195 * (t100 * t41 - t169 * t77) + t136 * t41 - t196 * (t101 * t41 + t170 * t77) + t194 * t40 - t123 * t77 - t190 * (-t40 * t106 + t41 * t99) + t197 * t128) * g(2) + (t195 * (t100 * t43 - t169 * t79) + t136 * t43 - t196 * (t101 * t43 + t170 * t79) + t194 * t42 - t123 * t79 - t190 * (-t42 * t106 + t43 * t99) + t197 * t125) * g(1) (-t190 * (-t106 * t122 - t45 * t99) + t194 * t122 + t191 * t45) * g(3) + (-t190 * (t106 * t29 - t25 * t99) - t194 * t29 + t191 * t25) * g(2) + (-t190 * (-t106 * t115 - t30 * t99) + t194 * t115 + t191 * t30) * g(1) (-t199 * mrSges(5,1) - (-t108 * t76 - t111 * t122) * mrSges(5,2) + t195 * t16 - t196 * (-t100 * t122 + t101 * t76)) * g(3) + (-t198 * mrSges(5,1) + mrSges(5,2) * t206 + t195 * t4 - t196 * t205) * g(2) + (-t9 * mrSges(5,1) + t10 * mrSges(5,2) + t195 * t8 + t196 * t7) * g(1) + (-g(1) * t9 - g(2) * t198 - g(3) * t199) * t200, t190 * (-g(1) * t30 - g(2) * t25 - g(3) * t45) -g(1) * (mrSges(7,1) * t1 - mrSges(7,2) * t2) - g(2) * ((-t107 * t4 + t110 * t25) * mrSges(7,1) + (-t107 * t25 - t110 * t4) * mrSges(7,2)) - g(3) * ((-t107 * t16 + t110 * t45) * mrSges(7,1) + (-t107 * t45 - t110 * t16) * mrSges(7,2))];
taug  = t3(:);
