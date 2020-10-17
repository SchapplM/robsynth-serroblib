% Calculate inertial parameters regressor of gravitation load for
% S6PRRRRP6
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% g_base [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [12x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,alpha2,alpha3,d2,d3,d4,d5,theta1]';
% 
% Output:
% taug_reg [6x(6*10)]
%   inertial parameter regressor of gravitation joint torque vector

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-09 00:36
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function taug_reg = S6PRRRRP6_gravloadJ_reg2_slag_vp(qJ, g, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(3,1),zeros(12,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6PRRRRP6_gravloadJ_reg2_slag_vp: qJ has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6PRRRRP6_gravloadJ_reg2_slag_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [12 1]), ...
  'S6PRRRRP6_gravloadJ_reg2_slag_vp: pkin has to be [12x1] (double)');

%% Symbolic Calculation
% From gravload_joint_fixb_regressor_matlab.m
% OptimizationMode: 2
% StartTime: 2019-05-05 10:16:51
% EndTime: 2019-05-05 10:16:53
% DurationCPUTime: 0.92s
% Computational Cost: add. (1282->172), mult. (3629->265), div. (0->0), fcn. (4689->14), ass. (0->119)
t101 = sin(qJ(4));
t105 = cos(qJ(4));
t165 = -pkin(4) * t105 - pkin(11) * t101;
t97 = sin(pkin(7));
t99 = cos(pkin(12));
t164 = t97 * t99;
t103 = sin(qJ(2));
t106 = cos(qJ(2));
t151 = sin(pkin(12));
t153 = cos(pkin(6));
t127 = t153 * t151;
t112 = t99 * t103 + t106 * t127;
t98 = sin(pkin(6));
t139 = t98 * t151;
t152 = cos(pkin(7));
t163 = t112 * t152 - t97 * t139;
t162 = pkin(9) * t97;
t161 = cos(qJ(3));
t156 = t103 * t98;
t148 = t97 * t156;
t154 = t106 * t98;
t158 = pkin(2) * t154 + pkin(9) * t148;
t157 = t101 * t97;
t155 = t105 * t97;
t100 = sin(qJ(5));
t150 = t100 * t105;
t104 = cos(qJ(5));
t149 = t104 * t105;
t147 = t98 * t161;
t102 = sin(qJ(3));
t138 = t99 * t153;
t111 = t151 * t103 - t106 * t138;
t107 = t111 * t152;
t87 = t103 * t138 + t151 * t106;
t50 = t87 * t102 + t161 * t107 + t147 * t164;
t51 = t87 * t161 + (-t98 * t164 - t107) * t102;
t146 = -t50 * pkin(3) + t51 * pkin(10);
t88 = -t103 * t127 + t99 * t106;
t52 = t88 * t102 + t163 * t161;
t53 = -t163 * t102 + t88 * t161;
t145 = -t52 * pkin(3) + t53 * pkin(10);
t132 = t152 * t161;
t137 = t153 * t97;
t72 = t102 * t156 - t132 * t154 - t161 * t137;
t140 = t98 * t152;
t73 = t103 * t147 + (t106 * t140 + t137) * t102;
t144 = -t72 * pkin(3) + t73 * pkin(10);
t74 = t111 * t97 - t99 * t140;
t25 = -t51 * t101 + t74 * t105;
t26 = t74 * t101 + t51 * t105;
t143 = t25 * pkin(4) + t26 * pkin(11);
t75 = t112 * t97 + t152 * t139;
t27 = -t53 * t101 + t75 * t105;
t28 = t75 * t101 + t53 * t105;
t142 = t27 * pkin(4) + t28 * pkin(11);
t86 = t153 * t152 - t97 * t154;
t54 = -t73 * t101 + t86 * t105;
t55 = t86 * t101 + t73 * t105;
t141 = t54 * pkin(4) + t55 * pkin(11);
t136 = t102 * t152;
t135 = -t111 * pkin(2) + t87 * t162;
t134 = -t112 * pkin(2) + t88 * t162;
t131 = t165 * t50 + t146;
t130 = t165 * t52 + t145;
t129 = t165 * t72 + t144;
t82 = (t102 * t106 + t103 * t132) * t98;
t83 = (-t103 * t136 + t161 * t106) * t98;
t128 = t83 * pkin(3) + t82 * pkin(10) + t158;
t126 = pkin(5) * t104 + qJ(6) * t100;
t10 = t28 * t100 - t52 * t104;
t23 = t55 * t100 - t72 * t104;
t8 = t26 * t100 - t50 * t104;
t1 = g(1) * t10 + g(2) * t8 + g(3) * t23;
t11 = t52 * t100 + t28 * t104;
t24 = t72 * t100 + t55 * t104;
t9 = t50 * t100 + t26 * t104;
t125 = g(1) * t11 + g(2) * t9 + g(3) * t24;
t13 = -t51 * t104 - t50 * t150;
t15 = -t53 * t104 - t52 * t150;
t29 = -t73 * t104 - t72 * t150;
t124 = g(1) * t15 + g(2) * t13 + g(3) * t29;
t61 = -t111 * t161 - t87 * t136;
t34 = t61 * t105 + t87 * t157;
t60 = -t111 * t102 + t87 * t132;
t17 = t34 * t100 - t60 * t104;
t63 = -t112 * t161 - t88 * t136;
t36 = t63 * t105 + t88 * t157;
t62 = -t112 * t102 + t88 * t132;
t19 = t36 * t100 - t62 * t104;
t66 = t101 * t148 + t83 * t105;
t37 = t66 * t100 - t82 * t104;
t123 = g(1) * t19 + g(2) * t17 + g(3) * t37;
t122 = g(1) * t27 + g(2) * t25 + g(3) * t54;
t121 = g(1) * t28 + g(2) * t26 + g(3) * t55;
t33 = t61 * t101 - t87 * t155;
t35 = t63 * t101 - t88 * t155;
t65 = t83 * t101 - t105 * t148;
t120 = g(1) * t35 + g(2) * t33 + g(3) * t65;
t119 = g(1) * t52 + g(2) * t50 + g(3) * t72;
t118 = g(1) * t53 + g(2) * t51 + g(3) * t73;
t117 = g(1) * t62 + g(2) * t60 + g(3) * t82;
t116 = t61 * pkin(3) + t60 * pkin(10) + t135;
t115 = t63 * pkin(3) + t62 * pkin(10) + t134;
t114 = g(1) * t88 + g(2) * t87 + g(3) * t156;
t113 = t66 * pkin(4) + t65 * pkin(11) + t128;
t110 = t34 * pkin(4) + t33 * pkin(11) + t116;
t109 = t36 * pkin(4) + t35 * pkin(11) + t115;
t38 = t82 * t100 + t66 * t104;
t30 = t73 * t100 - t72 * t149;
t20 = t62 * t100 + t36 * t104;
t18 = t60 * t100 + t34 * t104;
t16 = t53 * t100 - t52 * t149;
t14 = t51 * t100 - t50 * t149;
t12 = t119 * t101;
t5 = t122 * t104;
t4 = t122 * t100;
t3 = -g(1) * t20 - g(2) * t18 - g(3) * t38;
t2 = -g(1) * t16 - g(2) * t14 - g(3) * t30;
t6 = [0, 0, 0, 0, 0, 0, 0, 0, 0, -g(3), 0, 0, 0, 0, 0, 0, 0, 0, 0, -g(3), 0, 0, 0, 0, 0, 0, 0, 0, 0, -g(3), 0, 0, 0, 0, 0, 0, 0, 0, 0, -g(3), 0, 0, 0, 0, 0, 0, 0, 0, 0, -g(3), 0, 0, 0, 0, 0, 0, 0, 0, 0, -g(3); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, g(1) * t112 + g(2) * t111 - g(3) * t154, t114, 0, 0, 0, 0, 0, 0, 0, 0, -g(1) * t63 - g(2) * t61 - g(3) * t83, t117, -t114 * t97, -g(1) * t134 - g(2) * t135 - g(3) * t158, 0, 0, 0, 0, 0, 0, -g(1) * t36 - g(2) * t34 - g(3) * t66, t120, -t117, -g(1) * t115 - g(2) * t116 - g(3) * t128, 0, 0, 0, 0, 0, 0, t3, t123, -t120, -g(1) * t109 - g(2) * t110 - g(3) * t113, 0, 0, 0, 0, 0, 0, t3, -t120, -t123, -g(1) * (t20 * pkin(5) + t19 * qJ(6) + t109) - g(2) * (t18 * pkin(5) + t17 * qJ(6) + t110) - g(3) * (t38 * pkin(5) + t37 * qJ(6) + t113); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t119, t118, 0, 0, 0, 0, 0, 0, 0, 0, t119 * t105, -t12, -t118, -g(1) * t145 - g(2) * t146 - g(3) * t144, 0, 0, 0, 0, 0, 0, t2, t124, t12, -g(1) * t130 - g(2) * t131 - g(3) * t129, 0, 0, 0, 0, 0, 0, t2, t12, -t124, -g(1) * (t16 * pkin(5) + t15 * qJ(6) + t130) - g(2) * (t14 * pkin(5) + t13 * qJ(6) + t131) - g(3) * (t30 * pkin(5) + t29 * qJ(6) + t129); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t122, t121, 0, 0, 0, 0, 0, 0, 0, 0, -t5, t4, -t121, -g(1) * t142 - g(2) * t143 - g(3) * t141, 0, 0, 0, 0, 0, 0, -t5, -t121, -t4, -g(1) * (t126 * t27 + t142) - g(2) * (t126 * t25 + t143) - g(3) * (t126 * t54 + t141); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t1, t125, 0, 0, 0, 0, 0, 0, 0, 0, t1, 0, -t125, -g(1) * (-t10 * pkin(5) + t11 * qJ(6)) - g(2) * (-t8 * pkin(5) + t9 * qJ(6)) - g(3) * (-t23 * pkin(5) + t24 * qJ(6)); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t1;];
taug_reg  = t6;
