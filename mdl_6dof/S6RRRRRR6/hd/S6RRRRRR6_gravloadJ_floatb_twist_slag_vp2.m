% Calculate Gravitation load on the joints for
% S6RRRRRR6
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% g [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [12x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,alpha2,d1,d2,d3,d4,d5,d6]';
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
% Datum: 2018-11-23 18:42
% Revision: 76f9d5e39f14dc242b53c0d9d3d9db48bd8f37c0
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für mechatronische Systeme, Universität Hannover

function taug = S6RRRRRR6_gravloadJ_floatb_twist_slag_vp2(qJ, g, ...
  pkin, m, mrSges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(3,1),zeros(12,1),zeros(7,1),zeros(7,3)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRRRRR6_gravloadJ_floatb_twist_slag_vp2: qJ has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6RRRRRR6_gravloadJ_floatb_twist_slag_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [12 1]), ...
  'S6RRRRRR6_gravloadJ_floatb_twist_slag_vp2: pkin has to be [12x1] (double)');
assert( isreal(m) && all(size(m) == [7 1]), ...
  'S6RRRRRR6_gravloadJ_floatb_twist_slag_vp2: m has to be [7x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [7,3]), ...
  'S6RRRRRR6_gravloadJ_floatb_twist_slag_vp2: mrSges has to be [7x3] (double)');

%% Symbolic Calculation
% From joint_gravload_floatb_twist_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2018-11-23 18:41:29
% EndTime: 2018-11-23 18:41:31
% DurationCPUTime: 1.72s
% Computational Cost: add. (1886->171), mult. (1947->217), div. (0->0), fcn. (1941->18), ass. (0->84)
t89 = qJ(5) + qJ(6);
t84 = sin(t89);
t86 = cos(t89);
t93 = sin(qJ(5));
t97 = cos(qJ(5));
t197 = -mrSges(6,1) * t97 - mrSges(7,1) * t86 + mrSges(6,2) * t93 + mrSges(7,2) * t84 - mrSges(5,1);
t82 = pkin(5) * t97 + pkin(4);
t200 = m(6) * pkin(4) + m(7) * t82 - t197;
t192 = mrSges(5,2) - mrSges(6,3) - mrSges(7,3);
t100 = -pkin(12) - pkin(11);
t107 = -m(6) * pkin(11) + m(7) * t100 + t192;
t90 = qJ(3) + qJ(4);
t85 = sin(t90);
t87 = cos(t90);
t94 = sin(qJ(3));
t98 = cos(qJ(3));
t173 = m(4) * pkin(2) + t98 * mrSges(4,1) - t94 * mrSges(4,2) - t107 * t85 + t200 * t87 + mrSges(3,1);
t168 = cos(qJ(1));
t140 = pkin(6) + qJ(2);
t122 = sin(t140) / 0.2e1;
t141 = pkin(6) - qJ(2);
t127 = sin(t141);
t67 = t122 - t127 / 0.2e1;
t96 = sin(qJ(1));
t99 = cos(qJ(2));
t112 = t168 * t99 - t67 * t96;
t91 = sin(pkin(6));
t163 = t91 * t96;
t35 = -t112 * t94 + t163 * t98;
t123 = cos(t140) / 0.2e1;
t128 = cos(t141);
t68 = t123 - t128 / 0.2e1;
t92 = cos(pkin(6));
t190 = t68 * t94 + t92 * t98;
t113 = t168 * t67 + t96 * t99;
t135 = t91 * t168;
t29 = -t113 * t85 - t135 * t87;
t30 = t113 * t87 - t135 * t85;
t189 = t192 * t30 + t197 * t29;
t33 = t112 * t85 - t163 * t87;
t34 = t112 * t87 + t163 * t85;
t188 = t192 * t34 - t197 * t33;
t50 = t68 * t85 + t87 * t92;
t51 = -t68 * t87 + t85 * t92;
t187 = t192 * t51 + t197 * t50;
t115 = m(4) * pkin(9) - mrSges(3,2) + mrSges(4,3) + mrSges(5,3);
t186 = t84 * mrSges(7,1) + t97 * mrSges(6,2) + t86 * mrSges(7,2) + t115;
t172 = m(7) * pkin(5);
t182 = mrSges(6,1) + t172;
t179 = m(5) + m(6) + m(7);
t142 = t93 * t172;
t174 = -mrSges(6,1) * t93 - t142 - t186;
t104 = t128 / 0.2e1 + t123;
t95 = sin(qJ(2));
t56 = -t104 * t168 + t95 * t96;
t171 = (-t30 * t84 + t56 * t86) * mrSges(7,1) + (-t30 * t86 - t56 * t84) * mrSges(7,2);
t59 = t104 * t96 + t168 * t95;
t5 = -t34 * t84 + t59 * t86;
t6 = t34 * t86 + t59 * t84;
t170 = mrSges(7,1) * t5 - mrSges(7,2) * t6;
t66 = t122 + t127 / 0.2e1;
t169 = (-t51 * t84 - t66 * t86) * mrSges(7,1) + (-t51 * t86 + t66 * t84) * mrSges(7,2);
t158 = -t100 * t30 + t29 * t82;
t157 = -t100 * t34 - t33 * t82;
t150 = -t100 * t51 + t50 * t82;
t143 = pkin(1) * t168 + pkin(8) * t163;
t139 = t94 * t163;
t134 = -pkin(1) * t96 + pkin(8) * t135;
t133 = pkin(4) * t29 + pkin(11) * t30;
t132 = -pkin(4) * t33 + pkin(11) * t34;
t131 = pkin(4) * t50 + pkin(11) * t51;
t75 = t94 * t135;
t129 = -t113 * t98 + t75;
t126 = t35 * pkin(3);
t125 = t190 * pkin(3);
t101 = -pkin(10) - pkin(9);
t83 = pkin(3) * t98 + pkin(2);
t124 = pkin(3) * t139 - t101 * t59 + t112 * t83 + t143;
t7 = -t34 * t93 + t59 * t97;
t109 = t113 * t94 + t135 * t98;
t105 = t109 * pkin(3);
t36 = t112 * t98 + t139;
t8 = t34 * t97 + t59 * t93;
t1 = [(-t168 * mrSges(2,1) - m(3) * t143 - t112 * mrSges(3,1) - m(4) * (pkin(2) * t112 + t143) - t36 * mrSges(4,1) - t35 * mrSges(4,2) - m(5) * t124 - t34 * mrSges(5,1) - m(6) * (pkin(4) * t34 + t124) - t8 * mrSges(6,1) - t7 * mrSges(6,2) - m(7) * (t34 * t82 + t124) - t6 * mrSges(7,1) - t5 * mrSges(7,2) + (-mrSges(3,3) * t91 + mrSges(2,2)) * t96 + (-t115 - t142) * t59 + t107 * t33) * g(2) + (t96 * mrSges(2,1) + t168 * mrSges(2,2) - m(3) * t134 + t113 * mrSges(3,1) - mrSges(3,3) * t135 - m(4) * (-pkin(2) * t113 + t134) - t129 * mrSges(4,1) - t109 * mrSges(4,2) + t200 * t30 + (t182 * t93 + t186) * t56 + t107 * t29 + t179 * (-pkin(3) * t75 - t101 * t56 + t113 * t83 - t134)) * g(1) (-t179 * (t101 * t68 + t66 * t83) - t174 * t68 - t173 * t66) * g(3) + (-t179 * (-t101 * t113 - t56 * t83) + t174 * t113 + t173 * t56) * g(2) + (-t179 * (-t101 * t112 - t59 * t83) + t174 * t112 + t173 * t59) * g(1) (-t190 * mrSges(4,1) - (t68 * t98 - t92 * t94) * mrSges(4,2) - m(5) * t125 - m(6) * (t125 + t131) - m(7) * (t125 + t150) + t187) * g(3) + (mrSges(4,1) * t109 - mrSges(4,2) * t129 + m(5) * t105 - m(6) * (-t105 + t133) - m(7) * (-t105 + t158) + t189) * g(2) + (-mrSges(4,1) * t35 + mrSges(4,2) * t36 - m(5) * t126 - m(6) * (t126 + t132) - m(7) * (t126 + t157) + t188) * g(1) (-m(6) * t131 - m(7) * t150 + t187) * g(3) + (-m(6) * t133 - m(7) * t158 + t189) * g(2) + (-m(6) * t132 - m(7) * t157 + t188) * g(1) (-(-t51 * t97 + t66 * t93) * mrSges(6,2) - t169 - t182 * (-t51 * t93 - t66 * t97)) * g(3) + (-(-t30 * t97 - t56 * t93) * mrSges(6,2) - t171 - t182 * (-t30 * t93 + t56 * t97)) * g(2) + (mrSges(6,2) * t8 - t182 * t7 - t170) * g(1), -g(1) * t170 - g(2) * t171 - g(3) * t169];
taug  = t1(:);
