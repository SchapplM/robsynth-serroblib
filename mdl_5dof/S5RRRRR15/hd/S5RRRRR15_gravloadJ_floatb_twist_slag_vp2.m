% Calculate Gravitation load on the joints for
% S5RRRRR15
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% g [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [11x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,alpha2,alpha5,d1,d2,d3,d4,d5]';
% m [6x1]
%   mass of all robot links (including the base)
% mrSges [6x3]
%  first moment of all robot links (mass times center of mass in body frames)
%  rows: links of the robot (starting with base)
%  columns: x-, y-, z-coordinates
% 
% Output:
% taug [5x1]
%   joint torques required to compensate gravitation load

% Quelle: HybrDyn-Toolbox
% Datum: 2024-09-27 22:28
% Revision: 4128c0dfad36ceb08eaa7a9fbc8f797c888f59be (2024-07-29)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function taug = S5RRRRR15_gravloadJ_floatb_twist_slag_vp2(qJ, g, ...
  pkin, m, mrSges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(3,1),zeros(11,1),zeros(6,1),zeros(6,3)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RRRRR15_gravloadJ_floatb_twist_slag_vp2: qJ has to be [5x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S5RRRRR15_gravloadJ_floatb_twist_slag_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [11 1]), ...
  'S5RRRRR15_gravloadJ_floatb_twist_slag_vp2: pkin has to be [11x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5RRRRR15_gravloadJ_floatb_twist_slag_vp2: m has to be [6x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [6,3]), ...
  'S5RRRRR15_gravloadJ_floatb_twist_slag_vp2: mrSges has to be [6x3] (double)');

%% Symbolic Calculation
% From gravload_joint_floatb_twist_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2024-09-27 22:24:02
% EndTime: 2024-09-27 22:24:04
% DurationCPUTime: 1.08s
% Computational Cost: add. (907->182), mult. (1119->269), div. (0->0), fcn. (1151->26), ass. (0->132)
t165 = m(4) * pkin(2);
t172 = -mrSges(3,1) - t165;
t102 = cos(qJ(2));
t101 = cos(qJ(3));
t100 = cos(qJ(4));
t90 = sin(pkin(6));
t161 = pkin(11) * t90;
t95 = sin(qJ(4));
t66 = -t95 * pkin(4) + t100 * t161;
t59 = t66 * t101;
t65 = -pkin(4) * t100 - t95 * t161;
t62 = pkin(3) - t65;
t96 = sin(qJ(3));
t19 = -t96 * t62 + t59;
t16 = t19 * t102;
t155 = t66 * t96;
t18 = -t62 * t101 - t155;
t17 = pkin(2) - t18;
t97 = sin(qJ(2));
t6 = -t17 * t97 + t16;
t89 = qJ(2) + qJ(3);
t84 = pkin(5) - t89;
t115 = -qJ(4) + t84;
t109 = sin(t115);
t83 = pkin(5) + t89;
t77 = qJ(4) + t83;
t122 = sin(t77) / 0.2e1;
t158 = mrSges(6,3) * t90;
t87 = qJ(4) + t89;
t79 = sin(t87);
t126 = t79 * t158;
t92 = cos(pkin(6));
t99 = cos(qJ(5));
t151 = t92 * t99;
t94 = sin(qJ(5));
t152 = t92 * t94;
t71 = cos(t77) / 0.2e1;
t72 = cos(t115);
t80 = cos(t87);
t91 = sin(pkin(5));
t171 = -((-t79 * t152 + t80 * t99) * mrSges(6,1) + (-t79 * t151 - t80 * t94) * mrSges(6,2) + t126) * t91 - (t122 + t109 / 0.2e1) * mrSges(5,1) - (t71 - t72 / 0.2e1) * mrSges(5,2);
t103 = cos(qJ(1));
t136 = t103 * t79;
t98 = sin(qJ(1));
t148 = t98 * t80;
t149 = t98 * t79;
t52 = t72 / 0.2e1 + t71;
t26 = -t103 * t52 + t149;
t133 = t103 * t94;
t93 = cos(pkin(5));
t116 = t93 * t133;
t144 = t99 * t98;
t43 = t92 * t116 + t144;
t131 = t103 * t99;
t117 = t92 * t131;
t146 = t98 * t94;
t44 = t93 * t117 - t146;
t45 = -t93 * t131 + t92 * t146;
t46 = t92 * t144 + t116;
t51 = t122 - t109 / 0.2e1;
t170 = -(-t43 * t79 - t45 * t80) * mrSges(6,1) - (-t44 * t79 - t46 * t80) * mrSges(6,2) - (t93 * t136 + t148) * t158 + t26 * mrSges(5,1) - (-t103 * t51 - t148) * mrSges(5,2);
t135 = t103 * t80;
t27 = -t98 * t52 - t136;
t124 = t93 * t146;
t41 = t92 * t124 - t131;
t123 = t93 * t144;
t42 = -t92 * t123 - t133;
t47 = t92 * t133 + t123;
t48 = -t117 + t124;
t169 = -(t41 * t79 - t47 * t80) * mrSges(6,1) - (-t42 * t79 + t48 * t80) * mrSges(6,2) - (-t93 * t149 + t135) * t158 - t27 * mrSges(5,1) - (t98 * t51 - t135) * mrSges(5,2);
t73 = sin(t83);
t74 = sin(t84);
t75 = cos(t83);
t76 = cos(t84);
t168 = -(t73 / 0.2e1 + t74 / 0.2e1) * mrSges(4,1) - (t75 / 0.2e1 - t76 / 0.2e1) * mrSges(4,2) + t171;
t159 = -t103 / 0.2e1;
t85 = sin(t89);
t147 = t98 * t85;
t64 = t76 + t75;
t37 = t64 * t159 + t147;
t63 = t73 - t74;
t86 = cos(t89);
t167 = t37 * mrSges(4,1) - (t63 * t159 - t98 * t86) * mrSges(4,2) + t170;
t164 = t63 / 0.2e1;
t134 = t103 * t85;
t38 = -t134 - t98 * t64 / 0.2e1;
t166 = -t38 * mrSges(4,1) - (-t103 * t86 + t98 * t164) * mrSges(4,2) + t169;
t163 = pkin(8) + pkin(9);
t162 = m(5) * t91;
t160 = t97 * pkin(2);
t82 = t102 * pkin(2) + pkin(1);
t156 = t19 * t97;
t154 = t90 * t91;
t153 = t90 * t93;
t150 = t96 * t97;
t145 = t98 * t97;
t137 = t102 * t91;
t132 = t103 * t97;
t130 = t17 * t102;
t129 = t98 * t102;
t128 = t103 * t102;
t127 = pkin(10) + t163;
t125 = t98 * t154;
t118 = t103 * t154;
t114 = -t41 * t80 - t47 * t79;
t113 = -t44 * t80 + t46 * t79;
t112 = t130 + t156;
t24 = t65 * t101 - t155;
t25 = t65 * t96 + t59;
t111 = t102 * t24 - t25 * t97;
t81 = t101 * pkin(3) + pkin(2);
t110 = -pkin(3) * t150 + t102 * t81;
t57 = -t93 * t129 - t132;
t108 = t93 * t128 - t145;
t107 = pkin(3) * (t101 * t102 - t150);
t106 = -m(3) * pkin(8) - t92 * mrSges(6,3) - mrSges(3,3) - mrSges(4,3) - mrSges(5,3);
t105 = -t80 * mrSges(6,3) * t153 + m(4) * (t93 * t160 - t91 * t163) + m(5) * (-t91 * t127 + (t102 * t96 * pkin(3) + t97 * t81) * t93) - m(6) * (t91 * (t92 * pkin(11) + t127) + t6 * t93) + mrSges(4,1) * t164 + t51 * mrSges(5,1) + mrSges(2,2);
t104 = m(3) * pkin(1) + m(4) * t82 + m(5) * (pkin(3) * t86 + t82) + m(6) * (pkin(1) + t112) + t86 * mrSges(4,1) + t80 * mrSges(5,1) + mrSges(2,1) + t126;
t68 = -pkin(3) * t85 - t160;
t58 = -t93 * t145 + t128;
t56 = -t93 * t132 - t129;
t40 = t93 * t107;
t39 = t110 * t93;
t15 = t91 * t156;
t10 = t94 * t118 - t43 * t80 + t45 * t79;
t9 = t99 * t125 + t42 * t80 + t48 * t79;
t8 = t25 * t102 + t24 * t97;
t7 = t18 * t97 + t16;
t4 = t111 * t93;
t3 = (t102 * t18 - t156) * t93;
t2 = t112 * t93;
t1 = [(-t58 * mrSges(3,1) - t57 * mrSges(3,2) - t38 * mrSges(4,2) - t27 * mrSges(5,2) - t114 * mrSges(6,1) - t9 * mrSges(6,2) - t104 * t103 + ((-t90 * t94 * mrSges(6,1) + t106) * t91 + t105) * t98) * g(2) + (-t56 * mrSges(3,1) + t108 * mrSges(3,2) - t37 * mrSges(4,2) - t26 * mrSges(5,2) - t10 * mrSges(6,1) - t113 * mrSges(6,2) + t104 * t98 + ((-t90 * t99 * mrSges(6,2) + t106) * t91 + t105) * t103) * g(1), (-(mrSges(3,1) * t102 - mrSges(3,2) * t97) * t91 - t137 * t165 - t110 * t162 - m(6) * (t91 * t130 + t15) + t168) * g(3) + (-t56 * mrSges(3,2) - m(5) * (t103 * t39 + t98 * t68) - m(6) * (t2 * t103 + t6 * t98) + t167 + t172 * t108) * g(2) + (t58 * mrSges(3,2) - m(5) * (t103 * t68 - t98 * t39) - m(6) * (t6 * t103 - t2 * t98) + t172 * t57 + t166) * g(1), (-t107 * t162 - m(6) * (-t18 * t137 + t15) + t168) * g(3) + (-m(5) * (-pkin(3) * t147 + t103 * t40) - m(6) * (-t3 * t103 + t7 * t98) + t167) * g(2) + (-m(5) * (-pkin(3) * t134 - t98 * t40) - m(6) * (t7 * t103 + t3 * t98) + t166) * g(1), (m(6) * t111 * t91 + t171) * g(3) + (-m(6) * (-t4 * t103 + t8 * t98) + t170) * g(2) + (-m(6) * (t8 * t103 + t4 * t98) + t169) * g(1), -g(1) * (t9 * mrSges(6,1) + (-t94 * t125 - t114) * mrSges(6,2)) - g(2) * ((-t99 * t118 - t113) * mrSges(6,1) + t10 * mrSges(6,2)) - g(3) * ((mrSges(6,1) * t99 - mrSges(6,2) * t94) * t153 + ((t80 * t151 - t79 * t94) * mrSges(6,1) + (-t80 * t152 - t79 * t99) * mrSges(6,2)) * t91)];
taug = t1(:);
