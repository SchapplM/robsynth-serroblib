% Calculate Gravitation load on the joints for
% S6RRRPRR15
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% g [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [12x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,alpha2,alpha3,d1,d2,d3,d5,d6]';
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
% Datum: 2018-11-23 18:03
% Revision: 76f9d5e39f14dc242b53c0d9d3d9db48bd8f37c0
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für mechatronische Systeme, Universität Hannover

function taug = S6RRRPRR15_gravloadJ_floatb_twist_slag_vp2(qJ, g, ...
  pkin, m, mrSges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(3,1),zeros(12,1),zeros(7,1),zeros(7,3)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRRPRR15_gravloadJ_floatb_twist_slag_vp2: qJ has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6RRRPRR15_gravloadJ_floatb_twist_slag_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [12 1]), ...
  'S6RRRPRR15_gravloadJ_floatb_twist_slag_vp2: pkin has to be [12x1] (double)');
assert( isreal(m) && all(size(m) == [7 1]), ...
  'S6RRRPRR15_gravloadJ_floatb_twist_slag_vp2: m has to be [7x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [7,3]), ...
  'S6RRRPRR15_gravloadJ_floatb_twist_slag_vp2: mrSges has to be [7x3] (double)');

%% Symbolic Calculation
% From joint_gravload_floatb_twist_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2018-11-23 18:03:04
% EndTime: 2018-11-23 18:03:06
% DurationCPUTime: 1.66s
% Computational Cost: add. (4033->164), mult. (4229->227), div. (0->0), fcn. (4125->22), ass. (0->96)
t165 = m(6) + m(7);
t170 = m(5) + t165;
t123 = -qJ(4) * t170 + mrSges(4,2) - mrSges(5,3);
t124 = -pkin(11) * t165 - mrSges(4,1) + mrSges(5,2) - mrSges(6,3);
t92 = sin(qJ(6));
t97 = cos(qJ(6));
t106 = -t92 * mrSges(7,1) - t97 * mrSges(7,2) + t124;
t173 = m(7) * pkin(5) + mrSges(7,1) * t97 - mrSges(7,2) * t92 + mrSges(6,1);
t180 = m(7) * pkin(12) - mrSges(6,2) + mrSges(7,3);
t101 = cos(qJ(1));
t90 = sin(pkin(6));
t151 = t101 * t90;
t147 = pkin(6) + qJ(2);
t133 = cos(t147) / 0.2e1;
t148 = pkin(6) - qJ(2);
t139 = cos(t148);
t77 = t139 / 0.2e1 + t133;
t95 = sin(qJ(2));
t96 = sin(qJ(1));
t64 = -t101 * t77 + t96 * t95;
t89 = sin(pkin(7));
t91 = cos(pkin(7));
t46 = t91 * t151 - t64 * t89;
t157 = t90 * t96;
t67 = -t101 * t95 - t96 * t77;
t48 = t91 * t157 - t67 * t89;
t179 = pkin(3) * t170 - t106;
t145 = pkin(7) + qJ(3);
t129 = sin(t145) / 0.2e1;
t146 = pkin(7) - qJ(3);
t135 = sin(t146);
t171 = t129 - t135 / 0.2e1;
t100 = cos(qJ(2));
t130 = sin(t147) / 0.2e1;
t136 = sin(t148);
t75 = t130 - t136 / 0.2e1;
t68 = t100 * t101 - t96 * t75;
t99 = cos(qJ(3));
t178 = t171 * t67 + t68 * t99;
t65 = t96 * t100 + t101 * t75;
t177 = -t171 * t64 + t65 * t99;
t126 = t130 + t136 / 0.2e1;
t76 = t133 - t139 / 0.2e1;
t176 = t126 * t171 - t76 * t99;
t108 = t129 + t135 / 0.2e1;
t104 = t90 * t108;
t137 = cos(t145);
t131 = t137 / 0.2e1;
t138 = cos(t146);
t132 = t138 / 0.2e1;
t111 = t132 + t131;
t94 = sin(qJ(3));
t21 = t101 * t104 + t111 * t64 + t65 * t94;
t93 = sin(qJ(5));
t98 = cos(qJ(5));
t175 = t21 * t98 + t46 * t93;
t6 = t21 * t93 - t46 * t98;
t166 = -t173 * t93 + t180 * t98 + t123;
t159 = t89 * t93;
t158 = t89 * t98;
t156 = -mrSges(4,3) - mrSges(5,1);
t152 = t101 * pkin(1) + pkin(9) * t157;
t150 = cos(pkin(6));
t141 = -pkin(1) * t96 + pkin(9) * t151;
t140 = -mrSges(3,3) * t90 + mrSges(2,2);
t127 = t132 - t137 / 0.2e1;
t121 = t68 * pkin(2) + t48 * pkin(10) + t152;
t118 = t90 * t127;
t27 = t118 * t96 + t178;
t116 = t27 * pkin(3) + t121;
t114 = -t65 * pkin(2) + t46 * pkin(10) + t141;
t113 = t48 * pkin(4) + t116;
t22 = -t101 * t118 + t177;
t112 = -pkin(3) * t22 + t114;
t110 = t131 - t138 / 0.2e1;
t105 = t90 * t110;
t103 = -mrSges(3,2) + (t165 * pkin(4) + (m(4) + t170) * pkin(10) - t156) * t89;
t74 = t126 * pkin(2);
t63 = -t126 * t89 + t150 * t91;
t61 = t67 * pkin(2);
t59 = t64 * pkin(2);
t43 = t126 * t99 + t171 * t76;
t42 = -t111 * t76 + t126 * t94;
t39 = t127 * t150 + t176;
t38 = -t108 * t150 - t111 * t126 - t76 * t94;
t36 = -t171 * t68 + t67 * t99;
t35 = t111 * t68 + t67 * t94;
t34 = -t171 * t65 - t64 * t99;
t33 = t111 * t65 - t64 * t94;
t26 = -t104 * t96 - t111 * t67 + t68 * t94;
t14 = t38 * t93 + t63 * t98;
t8 = t26 * t93 + t48 * t98;
t7 = -t26 * t98 + t48 * t93;
t2 = t27 * t92 + t8 * t97;
t1 = t27 * t97 - t8 * t92;
t3 = [(-t101 * mrSges(2,1) - m(3) * t152 - t68 * mrSges(3,1) - t67 * mrSges(3,2) - m(4) * t121 - m(5) * t116 - m(6) * t113 - t8 * mrSges(6,1) - m(7) * (pkin(5) * t8 + t113) - t2 * mrSges(7,1) - t1 * mrSges(7,2) + t140 * t96 - t180 * t7 + t156 * t48 + t123 * t26 + t124 * t27) * g(2) + (t96 * mrSges(2,1) - m(3) * t141 + t65 * mrSges(3,1) - t64 * mrSges(3,2) - m(4) * t114 - m(5) * t112 + t156 * t46 - t180 * t175 - t123 * t21 + t140 * t101 + t173 * t6 - t106 * t22 + t165 * (-t46 * pkin(4) - t112)) * g(1) (-t126 * mrSges(3,1) - m(4) * t74 + t123 * t42 + t180 * (t159 * t76 + t42 * t98) - t173 * (-t158 * t76 + t42 * t93) + t106 * t43 + t103 * t76) * g(3) + (mrSges(3,1) * t64 + m(4) * t59 + t180 * (-t159 * t65 + t33 * t98) + t123 * t33 - t173 * (t158 * t65 + t33 * t93) + t106 * t34 - t103 * t65) * g(2) + (-mrSges(3,1) * t67 - m(4) * t61 + t123 * t35 + t180 * (-t159 * t68 + t35 * t98) - t173 * (t158 * t68 + t35 * t93) + t106 * t36 - t103 * t68) * g(1) + (-g(2) * (t34 * pkin(3) - t59) - g(1) * (t36 * pkin(3) + t61) - g(3) * (t43 * pkin(3) + t74)) * t170 (t166 * (-t110 * t150 + t176) + t179 * t38) * g(3) + (t166 * (t101 * t105 + t177) + t179 * t21) * g(2) + (t166 * (-t105 * t96 + t178) + t179 * t26) * g(1), t170 * (-g(1) * t26 - g(2) * t21 - g(3) * t38) (-t180 * t14 - t173 * (t38 * t98 - t63 * t93)) * g(3) + (-t173 * t175 - t180 * t6) * g(2) + (t173 * t7 - t180 * t8) * g(1), -g(1) * (mrSges(7,1) * t1 - mrSges(7,2) * t2) - g(2) * ((t22 * t97 - t6 * t92) * mrSges(7,1) + (-t22 * t92 - t6 * t97) * mrSges(7,2)) - g(3) * ((-t14 * t92 + t39 * t97) * mrSges(7,1) + (-t14 * t97 - t39 * t92) * mrSges(7,2))];
taug  = t3(:);
