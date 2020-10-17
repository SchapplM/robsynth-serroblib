% Calculate inertial parameters regressor of gravitation load for
% S6RPRRRP12
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% g_base [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [12x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,alpha2,alpha3,d1,d3,d4,d5,theta2]';
% 
% Output:
% taug_reg [6x(6*10)]
%   inertial parameter regressor of gravitation joint torque vector

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-09 06:53
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function taug_reg = S6RPRRRP12_gravloadJ_reg2_slag_vp(qJ, g, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(3,1),zeros(12,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPRRRP12_gravloadJ_reg2_slag_vp: qJ has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6RPRRRP12_gravloadJ_reg2_slag_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [12 1]), ...
  'S6RPRRRP12_gravloadJ_reg2_slag_vp: pkin has to be [12x1] (double)');

%% Symbolic Calculation
% From gravload_joint_fixb_regressor_matlab.m
% OptimizationMode: 2
% StartTime: 2019-05-06 02:22:04
% EndTime: 2019-05-06 02:22:08
% DurationCPUTime: 1.30s
% Computational Cost: add. (1321->134), mult. (3685->214), div. (0->0), fcn. (4786->14), ass. (0->93)
t143 = cos(pkin(12));
t145 = cos(pkin(6));
t122 = t145 * t143;
t140 = sin(pkin(12));
t149 = sin(qJ(1));
t89 = cos(qJ(1));
t106 = -t89 * t122 + t149 * t140;
t142 = sin(pkin(6));
t144 = cos(pkin(7));
t120 = t144 * t142;
t141 = sin(pkin(7));
t157 = t106 * t141 - t89 * t120;
t150 = cos(qJ(3));
t119 = t142 * t141;
t162 = t106 * t144 + t89 * t119;
t121 = t145 * t140;
t72 = t89 * t121 + t149 * t143;
t86 = sin(qJ(3));
t52 = -t72 * t150 + t162 * t86;
t85 = sin(qJ(4));
t88 = cos(qJ(4));
t28 = -t157 * t85 + t52 * t88;
t49 = t162 * t150 + t72 * t86;
t84 = sin(qJ(5));
t87 = cos(qJ(5));
t11 = t28 * t84 + t49 * t87;
t12 = t28 * t87 - t49 * t84;
t161 = t157 * t88 + t52 * t85;
t100 = t149 * t122 + t89 * t140;
t158 = t100 * t141 + t149 * t120;
t156 = -pkin(4) * t88 - pkin(11) * t85;
t155 = -t100 * t144 + t149 * t119;
t153 = t143 * t120 + t141 * t145;
t148 = t84 * t88;
t147 = t87 * t88;
t127 = t142 * t149;
t146 = t89 * pkin(1) + qJ(2) * t127;
t139 = -t49 * pkin(3) - pkin(10) * t52;
t73 = -t149 * t121 + t89 * t143;
t53 = -t155 * t150 + t73 * t86;
t54 = t73 * t150 + t155 * t86;
t138 = -t53 * pkin(3) + pkin(10) * t54;
t118 = t142 * t140;
t62 = t86 * t118 - t153 * t150;
t63 = t150 * t118 + t153 * t86;
t137 = -t62 * pkin(3) + pkin(10) * t63;
t136 = pkin(4) * t161 - pkin(11) * t28;
t29 = -t158 * t88 + t54 * t85;
t30 = t158 * t85 + t54 * t88;
t135 = -t29 * pkin(4) + t30 * pkin(11);
t71 = -t143 * t119 + t145 * t144;
t47 = -t63 * t85 + t71 * t88;
t48 = t63 * t88 + t71 * t85;
t134 = t47 * pkin(4) + t48 * pkin(11);
t133 = t89 * t142;
t132 = -t149 * pkin(1) + qJ(2) * t133;
t13 = t30 * t84 - t53 * t87;
t131 = g(1) * t11 + g(2) * t13;
t130 = g(1) * t161 + g(2) * t29;
t129 = -g(1) * t49 + g(2) * t53;
t126 = pkin(5) * t87 + qJ(6) * t84;
t125 = t156 * t49 + t139;
t124 = t156 * t53 + t138;
t123 = t156 * t62 + t137;
t19 = t48 * t84 - t62 * t87;
t1 = g(1) * t13 - g(2) * t11 + g(3) * t19;
t14 = t30 * t87 + t53 * t84;
t20 = t48 * t87 + t62 * t84;
t114 = g(1) * t14 - g(2) * t12 + g(3) * t20;
t15 = -t49 * t148 + t52 * t87;
t17 = -t53 * t148 - t54 * t87;
t31 = -t62 * t148 - t63 * t87;
t113 = g(1) * t17 + g(2) * t15 + g(3) * t31;
t112 = g(1) * t29 - g(2) * t161 - g(3) * t47;
t111 = g(1) * t30 - g(2) * t28 + g(3) * t48;
t110 = g(1) * t53 + g(2) * t49 + g(3) * t62;
t109 = g(1) * t54 - g(2) * t52 + g(3) * t63;
t96 = -t72 * pkin(2) - t157 * pkin(9) + t132;
t95 = t73 * pkin(2) + t158 * pkin(9) + t146;
t94 = t52 * pkin(3) - pkin(10) * t49 + t96;
t93 = t54 * pkin(3) + t53 * pkin(10) + t95;
t91 = t28 * pkin(4) + pkin(11) * t161 + t94;
t90 = t30 * pkin(4) + t29 * pkin(11) + t93;
t68 = -g(1) * t127 + g(2) * t133 - g(3) * t145;
t32 = -t62 * t147 + t63 * t84;
t18 = -t53 * t147 + t54 * t84;
t16 = -t49 * t147 - t52 * t84;
t8 = t110 * t85;
t5 = t112 * t87;
t4 = t112 * t84;
t3 = -g(1) * t12 - g(2) * t14;
t2 = -g(1) * t18 - g(2) * t16 - g(3) * t32;
t6 = [0, 0, 0, 0, 0, 0, g(1) * t149 - g(2) * t89, g(1) * t89 + g(2) * t149, 0, 0, 0, 0, 0, 0, 0, 0, g(1) * t72 - g(2) * t73, -g(1) * t106 + g(2) * t100, -g(1) * t133 - g(2) * t127, -g(1) * t132 - g(2) * t146, 0, 0, 0, 0, 0, 0, -g(1) * t52 - g(2) * t54, t129, g(1) * t157 - g(2) * t158, -g(1) * t96 - g(2) * t95, 0, 0, 0, 0, 0, 0, -g(1) * t28 - g(2) * t30, t130, -t129, -g(1) * t94 - g(2) * t93, 0, 0, 0, 0, 0, 0, t3, t131, -t130, -g(1) * t91 - g(2) * t90, 0, 0, 0, 0, 0, 0, t3, -t130, -t131, -g(1) * (t12 * pkin(5) + t11 * qJ(6) + t91) - g(2) * (t14 * pkin(5) + t13 * qJ(6) + t90); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t68, 0, 0, 0, 0, 0, 0, 0, 0, 0, t68, 0, 0, 0, 0, 0, 0, 0, 0, 0, t68, 0, 0, 0, 0, 0, 0, 0, 0, 0, t68, 0, 0, 0, 0, 0, 0, 0, 0, 0, t68; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t110, t109, 0, 0, 0, 0, 0, 0, 0, 0, t110 * t88, -t8, -t109, -g(1) * t138 - g(2) * t139 - g(3) * t137, 0, 0, 0, 0, 0, 0, t2, t113, t8, -g(1) * t124 - g(2) * t125 - g(3) * t123, 0, 0, 0, 0, 0, 0, t2, t8, -t113, -g(1) * (pkin(5) * t18 + qJ(6) * t17 + t124) - g(2) * (pkin(5) * t16 + qJ(6) * t15 + t125) - g(3) * (pkin(5) * t32 + qJ(6) * t31 + t123); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t112, t111, 0, 0, 0, 0, 0, 0, 0, 0, t5, -t4, -t111, -g(1) * t135 - g(2) * t136 - g(3) * t134, 0, 0, 0, 0, 0, 0, t5, -t111, t4, -g(1) * (-t126 * t29 + t135) - g(2) * (t126 * t161 + t136) - g(3) * (t126 * t47 + t134); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t1, t114, 0, 0, 0, 0, 0, 0, 0, 0, t1, 0, -t114, -g(1) * (-pkin(5) * t13 + qJ(6) * t14) - g(2) * (pkin(5) * t11 - qJ(6) * t12) - g(3) * (-t19 * pkin(5) + t20 * qJ(6)); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t1;];
taug_reg  = t6;
