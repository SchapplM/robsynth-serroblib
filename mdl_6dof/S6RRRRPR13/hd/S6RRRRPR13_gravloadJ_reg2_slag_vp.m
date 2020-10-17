% Calculate inertial parameters regressor of gravitation load for
% S6RRRRPR13
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% g_base [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [11x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,alpha2,d1,d2,d3,d4,d6]';
% 
% Output:
% taug_reg [6x(6*10)]
%   inertial parameter regressor of gravitation joint torque vector

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-10 00:07
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function taug_reg = S6RRRRPR13_gravloadJ_reg2_slag_vp(qJ, g, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(3,1),zeros(11,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRRRPR13_gravloadJ_reg2_slag_vp: qJ has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6RRRRPR13_gravloadJ_reg2_slag_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [11 1]), ...
  'S6RRRRPR13_gravloadJ_reg2_slag_vp: pkin has to be [11x1] (double)');

%% Symbolic Calculation
% From gravload_joint_fixb_regressor_matlab.m
% OptimizationMode: 2
% StartTime: 2019-05-08 01:35:38
% EndTime: 2019-05-08 01:35:43
% DurationCPUTime: 1.36s
% Computational Cost: add. (898->182), mult. (2353->262), div. (0->0), fcn. (2970->12), ass. (0->107)
t154 = cos(qJ(1));
t86 = sin(pkin(6));
t134 = t86 * t154;
t137 = cos(pkin(6));
t119 = t137 * t154;
t153 = sin(qJ(1));
t90 = sin(qJ(2));
t94 = cos(qJ(2));
t69 = t119 * t90 + t153 * t94;
t89 = sin(qJ(3));
t93 = cos(qJ(3));
t42 = -t89 * t134 + t69 * t93;
t68 = -t119 * t94 + t153 * t90;
t88 = sin(qJ(4));
t92 = cos(qJ(4));
t17 = t42 * t88 - t68 * t92;
t18 = t42 * t92 + t68 * t88;
t87 = sin(qJ(6));
t91 = cos(qJ(6));
t171 = t17 * t91 - t18 * t87;
t170 = t17 * t87 + t18 * t91;
t127 = -t93 * t134 - t69 * t89;
t138 = qJ(5) * t88;
t152 = t127 * t92;
t169 = pkin(4) * t152 + t127 * t138;
t133 = t86 * t153;
t118 = t137 * t153;
t71 = -t118 * t90 + t154 * t94;
t45 = -t133 * t93 + t71 * t89;
t151 = t45 * t92;
t168 = -pkin(4) * t151 - t45 * t138;
t145 = t86 * t90;
t66 = t137 * t93 - t145 * t89;
t150 = t66 * t92;
t167 = pkin(4) * t150 + t66 * t138;
t166 = -g(1) * t45 + g(2) * t127 + g(3) * t66;
t144 = t86 * t94;
t67 = t137 * t89 + t145 * t93;
t39 = t144 * t92 + t67 * t88;
t46 = t133 * t89 + t71 * t93;
t70 = t118 * t94 + t154 * t90;
t21 = t46 * t88 - t70 * t92;
t22 = t46 * t92 + t70 * t88;
t4 = t21 * t87 + t22 * t91;
t135 = t88 * t144;
t40 = t67 * t92 - t135;
t165 = g(2) * t170 + g(3) * (t39 * t87 + t40 * t91) + g(1) * t4;
t3 = t21 * t91 - t22 * t87;
t164 = g(2) * t171 + g(3) * (t39 * t91 - t40 * t87) + g(1) * t3;
t161 = pkin(10) - pkin(11);
t160 = pkin(3) * t93;
t156 = t127 * pkin(10);
t155 = t45 * pkin(10);
t148 = t68 * t89;
t146 = t70 * t89;
t143 = t88 * t93;
t142 = t92 * t93;
t141 = t93 * t94;
t140 = pkin(2) * t144 + pkin(9) * t145;
t139 = t154 * pkin(1) + pkin(8) * t133;
t136 = t89 * t144;
t132 = -t68 * pkin(2) + t69 * pkin(9);
t131 = -t70 * pkin(2) + t71 * pkin(9);
t35 = t127 * pkin(3);
t130 = t42 * pkin(10) + t35;
t37 = t45 * pkin(3);
t129 = t46 * pkin(10) - t37;
t61 = t66 * pkin(3);
t128 = t67 * pkin(10) + t61;
t126 = -t17 * pkin(4) + t18 * qJ(5);
t125 = -t21 * pkin(4) + t22 * qJ(5);
t124 = -t39 * pkin(4) + t40 * qJ(5);
t123 = t86 * pkin(3) * t141 + pkin(10) * t136 + t140;
t122 = -pkin(1) * t153 + pkin(8) * t134;
t121 = -g(1) * t17 + g(2) * t21;
t12 = g(1) * t127 + g(2) * t45;
t120 = g(1) * t68 - g(2) * t70;
t111 = -pkin(10) * t148 - t68 * t160 + t132;
t110 = -pkin(10) * t146 - t70 * t160 + t131;
t109 = t71 * pkin(2) + t70 * pkin(9) + t139;
t108 = t46 * pkin(3) + t109;
t2 = g(1) * t21 + g(2) * t17 + g(3) * t39;
t107 = g(1) * t22 + g(2) * t18 + g(3) * t40;
t26 = -t143 * t68 - t69 * t92;
t28 = -t143 * t70 - t71 * t92;
t48 = t135 * t93 - t145 * t92;
t106 = g(1) * t28 + g(2) * t26 + g(3) * t48;
t10 = g(1) * t46 + g(2) * t42 + g(3) * t67;
t49 = (t141 * t92 + t88 * t90) * t86;
t104 = t49 * pkin(4) + t48 * qJ(5) + t123;
t103 = g(1) * t154 + g(2) * t153;
t102 = -t69 * pkin(2) - t68 * pkin(9) + t122;
t101 = -g(1) * t70 - g(2) * t68 + g(3) * t144;
t100 = g(1) * t71 + g(2) * t69 + g(3) * t145;
t27 = -t142 * t68 + t69 * t88;
t99 = t27 * pkin(4) + t26 * qJ(5) + t111;
t29 = -t142 * t70 + t71 * t88;
t98 = t29 * pkin(4) + t28 * qJ(5) + t110;
t97 = -pkin(3) * t42 + t102;
t96 = t22 * pkin(4) + t21 * qJ(5) + t108;
t95 = -pkin(4) * t18 - qJ(5) * t17 + t97;
t23 = t101 * t89;
t8 = t166 * t92;
t7 = t166 * t88;
t6 = g(1) * t18 - g(2) * t22;
t5 = -g(1) * t29 - g(2) * t27 - g(3) * t49;
t1 = [0, 0, 0, 0, 0, 0, g(1) * t153 - g(2) * t154, t103, 0, 0, 0, 0, 0, 0, 0, 0, g(1) * t69 - g(2) * t71, -t120, -t103 * t86, -g(1) * t122 - g(2) * t139, 0, 0, 0, 0, 0, 0, g(1) * t42 - g(2) * t46, t12, t120, -g(1) * t102 - g(2) * t109, 0, 0, 0, 0, 0, 0, t6, t121, -t12, -g(1) * (t97 + t156) - g(2) * (t108 + t155) 0, 0, 0, 0, 0, 0, t6, -t12, -t121, -g(1) * (t95 + t156) - g(2) * (t96 + t155) 0, 0, 0, 0, 0, 0, g(1) * t170 - g(2) * t4, g(1) * t171 - g(2) * t3, t12, -g(1) * (-pkin(5) * t18 + t127 * t161 + t95) - g(2) * (t22 * pkin(5) + t161 * t45 + t96); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t101, t100, 0, 0, 0, 0, 0, 0, 0, 0, -t101 * t93, t23, -t100, -g(1) * t131 - g(2) * t132 - g(3) * t140, 0, 0, 0, 0, 0, 0, t5, t106, -t23, -g(1) * t110 - g(2) * t111 - g(3) * t123, 0, 0, 0, 0, 0, 0, t5, -t23, -t106, -g(1) * t98 - g(2) * t99 - g(3) * t104, 0, 0, 0, 0, 0, 0, -g(1) * (t28 * t87 + t29 * t91) - g(2) * (t26 * t87 + t27 * t91) - g(3) * (t48 * t87 + t49 * t91) -g(1) * (t28 * t91 - t29 * t87) - g(2) * (t26 * t91 - t27 * t87) - g(3) * (t48 * t91 - t49 * t87) t23, -g(1) * (t29 * pkin(5) + pkin(11) * t146 + t98) - g(2) * (t27 * pkin(5) + pkin(11) * t148 + t99) - g(3) * (t49 * pkin(5) - pkin(11) * t136 + t104); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t166, t10, 0, 0, 0, 0, 0, 0, 0, 0, -t8, t7, -t10, -g(1) * t129 - g(2) * t130 - g(3) * t128, 0, 0, 0, 0, 0, 0, -t8, -t10, -t7, -g(1) * (t129 + t168) - g(2) * (t130 + t169) - g(3) * (t128 + t167) 0, 0, 0, 0, 0, 0, -t166 * (t87 * t88 + t91 * t92) t166 * (t87 * t92 - t88 * t91) t10, -g(1) * (-pkin(5) * t151 + t161 * t46 + t168 - t37) - g(2) * (pkin(5) * t152 + t161 * t42 + t169 + t35) - g(3) * (pkin(5) * t150 + t161 * t67 + t167 + t61); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t2, t107, 0, 0, 0, 0, 0, 0, 0, 0, t2, 0, -t107, -g(1) * t125 - g(2) * t126 - g(3) * t124, 0, 0, 0, 0, 0, 0, t164, -t165, 0, -g(1) * (-t21 * pkin(5) + t125) - g(2) * (-t17 * pkin(5) + t126) - g(3) * (-t39 * pkin(5) + t124); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t2, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t2; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t164, t165, 0, 0;];
taug_reg  = t1;
