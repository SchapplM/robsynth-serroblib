% Calculate Gravitation load on the joints for
% S6RRRRRR5
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
% Datum: 2018-11-23 18:41
% Revision: 76f9d5e39f14dc242b53c0d9d3d9db48bd8f37c0
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für mechatronische Systeme, Universität Hannover

function taug = S6RRRRRR5_gravloadJ_floatb_twist_slag_vp2(qJ, g, ...
  pkin, m, mrSges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(3,1),zeros(12,1),zeros(7,1),zeros(7,3)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRRRRR5_gravloadJ_floatb_twist_slag_vp2: qJ has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6RRRRRR5_gravloadJ_floatb_twist_slag_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [12 1]), ...
  'S6RRRRRR5_gravloadJ_floatb_twist_slag_vp2: pkin has to be [12x1] (double)');
assert( isreal(m) && all(size(m) == [7 1]), ...
  'S6RRRRRR5_gravloadJ_floatb_twist_slag_vp2: m has to be [7x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [7,3]), ...
  'S6RRRRRR5_gravloadJ_floatb_twist_slag_vp2: mrSges has to be [7x3] (double)');

%% Symbolic Calculation
% From joint_gravload_floatb_twist_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2018-11-23 18:40:22
% EndTime: 2018-11-23 18:40:23
% DurationCPUTime: 1.30s
% Computational Cost: add. (1794->178), mult. (1680->224), div. (0->0), fcn. (1646->18), ass. (0->93)
t178 = mrSges(6,2) - mrSges(7,3);
t91 = sin(qJ(6));
t95 = cos(qJ(6));
t176 = mrSges(7,1) * t95 - mrSges(7,2) * t91 + mrSges(6,1);
t177 = m(7) * pkin(5) + t176;
t122 = -m(7) * pkin(12) + t178;
t89 = sin(pkin(6));
t94 = sin(qJ(1));
t151 = t89 * t94;
t157 = cos(qJ(1));
t97 = cos(qJ(2));
t132 = t157 * t97;
t139 = pkin(6) - qJ(2);
t123 = sin(t139);
t138 = pkin(6) + qJ(2);
t81 = sin(t138);
t161 = t81 / 0.2e1;
t137 = t161 - t123 / 0.2e1;
t52 = -t137 * t94 + t132;
t88 = qJ(3) + qJ(4);
t82 = sin(t88);
t83 = cos(t88);
t23 = t83 * t151 - t52 * t82;
t118 = cos(t138) / 0.2e1;
t124 = cos(t139);
t66 = t118 - t124 / 0.2e1;
t90 = cos(pkin(6));
t175 = t66 * t82 + t83 * t90;
t173 = m(6) + m(7);
t172 = m(5) * pkin(3) + mrSges(4,1);
t98 = -pkin(10) - pkin(9);
t102 = -m(4) * pkin(9) + m(5) * t98 + mrSges(3,2) - mrSges(4,3) - mrSges(5,3) - mrSges(6,3);
t171 = mrSges(7,1) * t91 + mrSges(7,2) * t95 - t102;
t84 = qJ(5) + t88;
t78 = sin(t84);
t79 = cos(t84);
t38 = t66 * t78 + t79 * t90;
t39 = -t66 * t79 + t78 * t90;
t170 = -t176 * t38 + t178 * t39;
t17 = -t151 * t79 + t52 * t78;
t18 = t151 * t78 + t52 * t79;
t169 = t176 * t17 + t178 * t18;
t133 = t89 * t157;
t150 = t94 * t97;
t49 = t137 * t157 + t150;
t13 = -t79 * t133 - t49 * t78;
t14 = -t78 * t133 + t49 * t79;
t168 = -t176 * t13 + t178 * t14;
t167 = -t175 * mrSges(5,1) - (t66 * t83 - t82 * t90) * mrSges(5,2) + t170;
t24 = t151 * t82 + t52 * t83;
t166 = -t23 * mrSges(5,1) + t24 * mrSges(5,2) + t169;
t105 = t133 * t83 + t49 * t82;
t126 = t82 * t133 - t49 * t83;
t165 = t105 * mrSges(5,1) - t126 * mrSges(5,2) + t168;
t96 = cos(qJ(3));
t85 = t96 * pkin(3);
t80 = t85 + pkin(2);
t92 = sin(qJ(3));
t164 = m(4) * pkin(2) + m(5) * t80 + t96 * mrSges(4,1) + t83 * mrSges(5,1) - t92 * mrSges(4,2) - t82 * mrSges(5,2) - t122 * t78 + t177 * t79 + mrSges(3,1);
t159 = pkin(3) * t92;
t68 = pkin(4) * t82 + t159;
t69 = pkin(4) * t83 + t85;
t145 = t69 * t151 - t52 * t68;
t141 = t66 * t68 + t90 * t69;
t140 = t157 * pkin(1) + pkin(8) * t151;
t130 = -t94 * pkin(1) + pkin(8) * t133;
t129 = t13 * pkin(5) + t14 * pkin(12);
t128 = -t17 * pkin(5) + pkin(12) * t18;
t127 = t38 * pkin(5) + pkin(12) * t39;
t72 = t92 * t133;
t125 = -t49 * t96 + t72;
t121 = t23 * pkin(4);
t120 = t175 * pkin(4);
t100 = t124 / 0.2e1 + t118;
t93 = sin(qJ(2));
t51 = t100 * t94 + t157 * t93;
t64 = pkin(2) + t69;
t87 = -pkin(11) + t98;
t119 = t68 * t151 - t51 * t87 + t52 * t64 + t140;
t117 = t123 / 0.2e1;
t110 = -t133 * t69 - t49 * t68;
t25 = t151 * t96 - t52 * t92;
t108 = t117 - t81 / 0.2e1;
t104 = t133 * t96 + t49 * t92;
t101 = t105 * pkin(4);
t65 = t161 + t117;
t53 = t108 * t94 + t132;
t50 = -t108 * t157 + t150;
t48 = -t100 * t157 + t93 * t94;
t26 = t151 * t92 + t52 * t96;
t2 = t18 * t95 + t51 * t91;
t1 = -t18 * t91 + t51 * t95;
t3 = [(-t157 * mrSges(2,1) - m(3) * t140 - t52 * mrSges(3,1) - m(4) * (pkin(2) * t52 + t140) - t26 * mrSges(4,1) - t25 * mrSges(4,2) - m(5) * (t52 * t80 + t140) - t24 * mrSges(5,1) - t23 * mrSges(5,2) - m(6) * t119 - t18 * mrSges(6,1) - m(7) * (pkin(5) * t18 + t119) - t2 * mrSges(7,1) - t1 * mrSges(7,2) + (mrSges(2,2) + (-m(5) * t159 - mrSges(3,3)) * t89) * t94 + t122 * t17 + t102 * t51) * g(2) + (t94 * mrSges(2,1) + t157 * mrSges(2,2) - m(3) * t130 + t49 * mrSges(3,1) - mrSges(3,3) * t133 - m(4) * (-pkin(2) * t49 + t130) - t125 * mrSges(4,1) - t104 * mrSges(4,2) - m(5) * (pkin(3) * t72 - t49 * t80 + t130) - t126 * mrSges(5,1) - t105 * mrSges(5,2) + t122 * t13 + t177 * t14 + t171 * t48 + t173 * (-t68 * t133 - t48 * t87 + t49 * t64 - t130)) * g(1) (-t173 * (t65 * t64 + t66 * t87) + t171 * t66 - t164 * t65) * g(3) + (-t173 * (-t48 * t64 - t50 * t87) - t171 * t50 + t164 * t48) * g(2) + (-t173 * (-t51 * t64 - t53 * t87) - t171 * t53 + t164 * t51) * g(1) (-(t66 * t96 - t90 * t92) * mrSges(4,2) - m(6) * t141 - m(7) * (t127 + t141) - t172 * (t66 * t92 + t90 * t96) + t167) * g(3) + (-t125 * mrSges(4,2) - m(6) * t110 - m(7) * (t110 + t129) + t172 * t104 + t165) * g(2) + (mrSges(4,2) * t26 - m(6) * t145 - m(7) * (t128 + t145) - t172 * t25 + t166) * g(1) (-m(6) * t120 - m(7) * (t120 + t127) + t167) * g(3) + (m(6) * t101 - m(7) * (-t101 + t129) + t165) * g(2) + (-m(6) * t121 - m(7) * (t121 + t128) + t166) * g(1) (-m(7) * t127 + t170) * g(3) + (-m(7) * t129 + t168) * g(2) + (-m(7) * t128 + t169) * g(1), -g(1) * (mrSges(7,1) * t1 - mrSges(7,2) * t2) - g(2) * ((-t14 * t91 + t48 * t95) * mrSges(7,1) + (-t14 * t95 - t48 * t91) * mrSges(7,2)) - g(3) * ((-t39 * t91 - t65 * t95) * mrSges(7,1) + (-t39 * t95 + t65 * t91) * mrSges(7,2))];
taug  = t3(:);
