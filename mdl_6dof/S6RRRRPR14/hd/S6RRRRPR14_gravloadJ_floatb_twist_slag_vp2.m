% Calculate Gravitation load on the joints for
% S6RRRRPR14
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
% Datum: 2018-11-23 18:24
% Revision: 76f9d5e39f14dc242b53c0d9d3d9db48bd8f37c0
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für mechatronische Systeme, Universität Hannover

function taug = S6RRRRPR14_gravloadJ_floatb_twist_slag_vp2(qJ, g, ...
  pkin, m, mrSges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(3,1),zeros(13,1),zeros(7,1),zeros(7,3)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRRRPR14_gravloadJ_floatb_twist_slag_vp2: qJ has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6RRRRPR14_gravloadJ_floatb_twist_slag_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [13 1]), ...
  'S6RRRRPR14_gravloadJ_floatb_twist_slag_vp2: pkin has to be [13x1] (double)');
assert( isreal(m) && all(size(m) == [7 1]), ...
  'S6RRRRPR14_gravloadJ_floatb_twist_slag_vp2: m has to be [7x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [7,3]), ...
  'S6RRRRPR14_gravloadJ_floatb_twist_slag_vp2: mrSges has to be [7x3] (double)');

%% Symbolic Calculation
% From joint_gravload_floatb_twist_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2018-11-23 18:23:40
% EndTime: 2018-11-23 18:23:42
% DurationCPUTime: 2.15s
% Computational Cost: add. (4952->175), mult. (5189->237), div. (0->0), fcn. (5136->24), ass. (0->99)
t85 = sin(pkin(13));
t87 = cos(pkin(13));
t174 = -m(7) * (pkin(5) * t87 + pkin(4)) - m(6) * pkin(4) - mrSges(6,1) * t87 + mrSges(6,2) * t85 - mrSges(5,1);
t84 = pkin(13) + qJ(6);
t81 = sin(t84);
t82 = cos(t84);
t203 = -mrSges(7,1) * t82 + mrSges(7,2) * t81 + t174;
t155 = pkin(6) + qJ(2);
t141 = cos(t155) / 0.2e1;
t156 = pkin(6) - qJ(2);
t149 = cos(t156);
t122 = t149 / 0.2e1 + t141;
t171 = sin(qJ(2));
t172 = sin(qJ(1));
t173 = cos(qJ(1));
t113 = -t122 * t173 + t172 * t171;
t158 = sin(pkin(6));
t159 = cos(pkin(7));
t134 = t159 * t158;
t86 = sin(pkin(7));
t101 = t113 * t86 - t173 * t134;
t154 = pkin(7) - qJ(3);
t148 = cos(t154);
t140 = t148 / 0.2e1;
t153 = pkin(7) + qJ(3);
t147 = cos(t153);
t133 = t140 - t147 / 0.2e1;
t124 = t133 * t158;
t137 = sin(t153) / 0.2e1;
t145 = sin(t154);
t183 = t137 - t145 / 0.2e1;
t138 = sin(t155) / 0.2e1;
t146 = sin(t156);
t69 = t138 - t146 / 0.2e1;
t93 = cos(qJ(2));
t60 = t172 * t93 + t173 * t69;
t92 = cos(qJ(3));
t193 = -t113 * t183 + t60 * t92;
t23 = t124 * t173 - t193;
t89 = sin(qJ(4));
t91 = cos(qJ(4));
t4 = t101 * t89 - t23 * t91;
t202 = t101 * t91 + t23 * t89;
t191 = -m(5) - m(6);
t180 = -m(6) * qJ(5) + mrSges(5,2) - mrSges(6,3) + m(7) * (-pkin(12) - qJ(5)) - mrSges(7,3);
t105 = t172 * t122 + t171 * t173;
t94 = t105 * t86 + t172 * t134;
t135 = t81 * mrSges(7,1) + t82 * mrSges(7,2);
t196 = -t87 * mrSges(6,2) + mrSges(4,2) - mrSges(5,3);
t175 = -t85 * mrSges(6,1) - m(7) * (pkin(5) * t85 + pkin(11)) + t196;
t198 = -t135 + t175;
t190 = m(6) + m(7);
t182 = m(5) + t190;
t197 = pkin(3) * t182 - t180 * t89 - t203 * t91 + mrSges(4,1);
t62 = -t172 * t69 + t173 * t93;
t194 = -t105 * t183 + t62 * t92;
t119 = t138 + t146 / 0.2e1;
t70 = t141 - t149 / 0.2e1;
t192 = t119 * t183 - t70 * t92;
t117 = t137 + t145 / 0.2e1;
t114 = t117 * t158;
t139 = t147 / 0.2e1;
t121 = t140 + t139;
t90 = sin(qJ(3));
t19 = t113 * t121 + t114 * t173 + t60 * t90;
t177 = t191 * pkin(11) + t198;
t167 = t86 * t89;
t166 = t86 * t91;
t142 = t158 * t172;
t161 = t173 * pkin(1) + pkin(9) * t142;
t160 = cos(pkin(6));
t143 = t173 * t158;
t144 = -pkin(1) * t172 + pkin(9) * t143;
t123 = -mrSges(3,2) + (mrSges(4,3) + (m(4) + t182) * pkin(10)) * t86;
t120 = t139 - t148 / 0.2e1;
t115 = t120 * t158;
t107 = -t119 * t86 + t159 * t160;
t104 = (-m(7) * pkin(5) - mrSges(6,1)) * t85 - t182 * pkin(11) - t135 + t196;
t102 = -t60 * pkin(2) - t101 * pkin(10) + t144;
t99 = t23 * pkin(3) + t102;
t98 = t62 * pkin(2) + t94 * pkin(10) + t161;
t25 = t124 * t172 + t194;
t97 = t25 * pkin(3) + t98;
t68 = t119 * pkin(2);
t58 = t105 * pkin(2);
t56 = t113 * pkin(2);
t41 = t119 * t92 + t183 * t70;
t37 = t133 * t160 + t192;
t36 = -t117 * t160 - t119 * t121 - t70 * t90;
t34 = -t105 * t92 - t183 * t62;
t32 = -t113 * t92 - t183 * t60;
t24 = t105 * t121 - t114 * t172 + t62 * t90;
t14 = t107 * t89 + t37 * t91;
t13 = -t107 * t91 + t37 * t89;
t8 = t25 * t91 + t89 * t94;
t7 = t25 * t89 - t91 * t94;
t2 = t24 * t81 + t8 * t82;
t1 = t24 * t82 - t8 * t81;
t3 = [(-m(3) * t161 - m(4) * t98 - m(7) * t97 - mrSges(2,1) * t173 - t62 * mrSges(3,1) - t25 * mrSges(4,1) - t2 * mrSges(7,1) + t172 * mrSges(2,2) + mrSges(3,2) * t105 - t1 * mrSges(7,2) - mrSges(3,3) * t142 - mrSges(4,3) * t94 + t191 * (t24 * pkin(11) + t97) + t174 * t8 + t175 * t24 + t180 * t7) * g(2) + (-m(3) * t144 - m(4) * t102 - m(7) * t99 + t172 * mrSges(2,1) + t60 * mrSges(3,1) - t23 * mrSges(4,1) + mrSges(2,2) * t173 - mrSges(3,2) * t113 - mrSges(3,3) * t143 + mrSges(4,3) * t101 + t191 * (-pkin(11) * t19 + t99) - t203 * t4 - t198 * t19 + t180 * t202) * g(1) (-t119 * mrSges(3,1) - m(4) * t68 - t41 * mrSges(4,1) + t123 * t70 + t203 * (-t167 * t70 + t41 * t91) + t104 * (t119 * t90 - t121 * t70) + t180 * (t166 * t70 + t41 * t89)) * g(3) + (t113 * mrSges(3,1) + m(4) * t56 - t32 * mrSges(4,1) - t123 * t60 + t180 * (-t166 * t60 + t32 * t89) + t203 * (t167 * t60 + t32 * t91) + t104 * (-t113 * t90 + t121 * t60)) * g(2) + (t105 * mrSges(3,1) + m(4) * t58 - t34 * mrSges(4,1) - t123 * t62 + t203 * (t167 * t62 + t34 * t91) + t104 * (-t105 * t90 + t121 * t62) + t180 * (-t166 * t62 + t34 * t89)) * g(1) + (-g(2) * (t32 * pkin(3) - t56) - g(1) * (t34 * pkin(3) - t58) - g(3) * (t41 * pkin(3) + t68)) * t182 (t177 * (-t120 * t160 + t192) + t197 * t36) * g(3) + (t177 * (t115 * t173 + t193) + t197 * t19) * g(2) + (t177 * (-t115 * t172 + t194) + t197 * t24) * g(1) (-t13 * t203 + t180 * t14) * g(3) + (t180 * t4 + t202 * t203) * g(2) + (t180 * t8 - t203 * t7) * g(1), t190 * (-g(1) * t7 + g(2) * t202 - g(3) * t13) -g(1) * (mrSges(7,1) * t1 - mrSges(7,2) * t2) - g(2) * ((t19 * t82 - t4 * t81) * mrSges(7,1) + (-t19 * t81 - t4 * t82) * mrSges(7,2)) - g(3) * ((-t14 * t81 + t36 * t82) * mrSges(7,1) + (-t14 * t82 - t36 * t81) * mrSges(7,2))];
taug  = t3(:);
