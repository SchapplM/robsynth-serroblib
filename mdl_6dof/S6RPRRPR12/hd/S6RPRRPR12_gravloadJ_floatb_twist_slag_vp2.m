% Calculate Gravitation load on the joints for
% S6RPRRPR12
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% g [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [12x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,alpha2,alpha3,d1,d3,d4,d6,theta2]';
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
% Datum: 2018-11-23 16:23
% Revision: 76f9d5e39f14dc242b53c0d9d3d9db48bd8f37c0
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für mechatronische Systeme, Universität Hannover

function taug = S6RPRRPR12_gravloadJ_floatb_twist_slag_vp2(qJ, g, ...
  pkin, m, mrSges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(3,1),zeros(12,1),zeros(7,1),zeros(7,3)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPRRPR12_gravloadJ_floatb_twist_slag_vp2: qJ has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6RPRRPR12_gravloadJ_floatb_twist_slag_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [12 1]), ...
  'S6RPRRPR12_gravloadJ_floatb_twist_slag_vp2: pkin has to be [12x1] (double)');
assert( isreal(m) && all(size(m) == [7 1]), ...
  'S6RPRRPR12_gravloadJ_floatb_twist_slag_vp2: m has to be [7x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [7,3]), ...
  'S6RPRRPR12_gravloadJ_floatb_twist_slag_vp2: mrSges has to be [7x3] (double)');

%% Symbolic Calculation
% From joint_gravload_floatb_twist_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2018-11-23 16:23:10
% EndTime: 2018-11-23 16:23:11
% DurationCPUTime: 1.56s
% Computational Cost: add. (3528->156), mult. (3779->202), div. (0->0), fcn. (3757->22), ass. (0->103)
t161 = mrSges(5,2) - mrSges(6,3);
t70 = sin(qJ(6));
t73 = cos(qJ(6));
t175 = t70 * mrSges(7,1) + t73 * mrSges(7,2) - t161;
t129 = pkin(7) - qJ(3);
t119 = cos(t129);
t113 = t119 / 0.2e1;
t128 = pkin(7) + qJ(3);
t118 = cos(t128);
t104 = t113 - t118 / 0.2e1;
t68 = sin(pkin(6));
t101 = t68 * t104;
t145 = cos(qJ(1));
t111 = sin(t128) / 0.2e1;
t117 = sin(t129);
t157 = t111 - t117 / 0.2e1;
t144 = sin(qJ(1));
t126 = pkin(6) + pkin(12);
t107 = sin(t126) / 0.2e1;
t127 = pkin(6) - pkin(12);
t115 = sin(t127);
t57 = t107 - t115 / 0.2e1;
t69 = cos(pkin(12));
t51 = t144 * t69 + t145 * t57;
t75 = cos(qJ(3));
t130 = sin(pkin(12));
t108 = cos(t127) / 0.2e1;
t116 = cos(t126);
t95 = t108 + t116 / 0.2e1;
t90 = t144 * t130 - t145 * t95;
t163 = t157 * t90 - t51 * t75;
t28 = t101 * t145 + t163;
t71 = sin(qJ(4));
t74 = cos(qJ(4));
t131 = cos(pkin(7));
t120 = t68 * t131;
t67 = sin(pkin(7));
t84 = -t145 * t120 + t90 * t67;
t174 = -t28 * t74 + t71 * t84;
t9 = t28 * t71 + t74 * t84;
t162 = m(6) + m(7);
t170 = mrSges(5,1) - mrSges(6,2) + mrSges(7,3);
t86 = t130 * t145 + t144 * t95;
t78 = t144 * t120 + t86 * t67;
t153 = m(7) * pkin(11) + t170;
t166 = pkin(4) * t162 + t153;
t58 = t108 - t116 / 0.2e1;
t94 = t107 + t115 / 0.2e1;
t165 = t157 * t94 + t58 * t75;
t52 = -t144 * t57 + t145 * t69;
t164 = -t157 * t86 + t52 * t75;
t133 = qJ(5) * t71;
t72 = sin(qJ(3));
t96 = t111 + t117 / 0.2e1;
t92 = t68 * t96;
t112 = t118 / 0.2e1;
t99 = t113 + t112;
t24 = t145 * t92 + t51 * t72 + t90 * t99;
t143 = t24 * t74;
t160 = -pkin(4) * t143 - t24 * t133;
t29 = -t144 * t92 + t52 * t72 + t86 * t99;
t142 = t29 * t74;
t159 = -pkin(4) * t142 - t29 * t133;
t132 = cos(pkin(6));
t35 = -t132 * t96 + t58 * t72 - t94 * t99;
t141 = t35 * t74;
t158 = -pkin(4) * t141 - t35 * t133;
t152 = -m(7) * (pkin(5) + pkin(10)) - mrSges(6,1) + mrSges(4,2) - mrSges(5,3);
t151 = -t162 * qJ(5) - t175;
t150 = t170 * t74 + t175 * t71 + mrSges(4,1);
t149 = -t73 * mrSges(7,1) + t70 * mrSges(7,2) + t152;
t147 = t24 * pkin(10);
t146 = t29 * pkin(10);
t124 = t68 * t144;
t134 = t145 * pkin(1) + qJ(2) * t124;
t125 = t68 * t145;
t20 = t24 * pkin(3);
t98 = t112 - t119 / 0.2e1;
t93 = t68 * t98;
t26 = t145 * t93 - t163;
t123 = pkin(10) * t26 - t20;
t22 = t29 * pkin(3);
t31 = -t144 * t93 + t164;
t122 = pkin(10) * t31 - t22;
t34 = t35 * pkin(3);
t37 = -t132 * t98 + t165;
t121 = pkin(10) * t37 - t34;
t114 = -pkin(1) * t144 + qJ(2) * t125;
t87 = t131 * t132 - t67 * t94;
t82 = -t51 * pkin(2) - pkin(9) * t84 + t114;
t81 = t28 * pkin(3) + t82;
t80 = t52 * pkin(2) + pkin(9) * t78 + t134;
t30 = t101 * t144 + t164;
t79 = t30 * pkin(3) + t80;
t77 = -pkin(4) * t174 + t9 * qJ(5) + t81;
t11 = t30 * t71 - t74 * t78;
t12 = t30 * t74 + t71 * t78;
t76 = t12 * pkin(4) + t11 * qJ(5) + t79;
t36 = t104 * t132 + t165;
t14 = t36 * t71 - t74 * t87;
t2 = t11 * t70 + t29 * t73;
t1 = t11 * t73 - t29 * t70;
t3 = [(-t145 * mrSges(2,1) + t144 * mrSges(2,2) - m(3) * t134 - t52 * mrSges(3,1) + t86 * mrSges(3,2) - mrSges(3,3) * t124 - m(4) * t80 - t30 * mrSges(4,1) - t78 * mrSges(4,3) - m(5) * (t79 + t146) - m(6) * (t76 + t146) - m(7) * t76 - t2 * mrSges(7,1) - t1 * mrSges(7,2) + t161 * t11 + t152 * t29 - t153 * t12) * g(2) + (t144 * mrSges(2,1) + t145 * mrSges(2,2) - m(3) * t114 + t51 * mrSges(3,1) - t90 * mrSges(3,2) - mrSges(3,3) * t125 - m(4) * t82 - t28 * mrSges(4,1) + t84 * mrSges(4,3) - m(5) * (t81 - t147) - m(6) * (t77 - t147) - m(7) * t77 - t175 * t9 - t149 * t24 + t153 * t174) * g(1) (-t132 * g(3) + (-t144 * g(1) + t145 * g(2)) * t68) * (m(3) + m(4) + m(5) + t162) (-m(5) * t121 - m(6) * (t121 + t158) - m(7) * (-pkin(11) * t141 + t158 - t34) + t149 * t37 + t150 * t35) * g(3) + (-m(5) * t123 - m(6) * (t123 + t160) - m(7) * (-pkin(11) * t143 + t160 - t20) + t149 * t26 + t150 * t24) * g(2) + (-m(5) * t122 - m(6) * (t122 + t159) - m(7) * (-pkin(11) * t142 + t159 - t22) + t149 * t31 + t150 * t29) * g(1) (t151 * (t36 * t74 + t71 * t87) + t166 * t14) * g(3) + (t151 * t174 - t166 * t9) * g(2) + (t11 * t166 + t151 * t12) * g(1), t162 * (-g(1) * t11 + g(2) * t9 - g(3) * t14) -g(1) * (mrSges(7,1) * t1 - mrSges(7,2) * t2) - g(2) * ((-t24 * t70 - t73 * t9) * mrSges(7,1) + (-t24 * t73 + t70 * t9) * mrSges(7,2)) - g(3) * ((t14 * t73 - t35 * t70) * mrSges(7,1) + (-t14 * t70 - t35 * t73) * mrSges(7,2))];
taug  = t3(:);
