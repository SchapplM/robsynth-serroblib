% Calculate inertial parameters regressor of gravitation load for
% S6RRRPRR15
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% g_base [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [12x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,alpha2,alpha3,d1,d2,d3,d5,d6]';
% 
% Output:
% taug_reg [6x(6*10)]
%   inertial parameter regressor of gravitation joint torque vector

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-09 20:42
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function taug_reg = S6RRRPRR15_gravloadJ_reg2_slag_vp(qJ, g, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(3,1),zeros(12,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRRPRR15_gravloadJ_reg2_slag_vp: qJ has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6RRRPRR15_gravloadJ_reg2_slag_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [12 1]), ...
  'S6RRRPRR15_gravloadJ_reg2_slag_vp: pkin has to be [12x1] (double)');

%% Symbolic Calculation
% From gravload_joint_fixb_regressor_matlab.m
t136 = cos(pkin(7));
t84 = sin(pkin(6));
t129 = t84 * t136;
t137 = cos(pkin(6));
t93 = cos(qJ(1));
t126 = t93 * t137;
t88 = sin(qJ(2));
t89 = sin(qJ(1));
t92 = cos(qJ(2));
t64 = -t92 * t126 + t88 * t89;
t83 = sin(pkin(7));
t163 = t93 * t129 - t64 * t83;
t150 = cos(qJ(3));
t113 = t136 * t150;
t132 = t84 * t150;
t122 = t83 * t132;
t65 = t88 * t126 + t89 * t92;
t87 = sin(qJ(3));
t27 = t64 * t113 + t93 * t122 + t65 * t87;
t86 = sin(qJ(5));
t91 = cos(qJ(5));
t11 = -t163 * t91 + t27 * t86;
t128 = t87 * t136;
t142 = t84 * t93;
t28 = -t83 * t87 * t142 - t64 * t128 + t65 * t150;
t85 = sin(qJ(6));
t90 = cos(qJ(6));
t167 = t11 * t85 - t28 * t90;
t166 = t11 * t90 + t28 * t85;
t127 = t89 * t137;
t66 = -t92 * t127 - t93 * t88;
t49 = t89 * t129 - t66 * t83;
t162 = t163 * t86 + t27 * t91;
t67 = -t88 * t127 + t92 * t93;
t161 = -g(1) * t67 - g(2) * t65;
t144 = t84 * t89;
t32 = t67 * t150 + (t136 * t66 + t83 * t144) * t87;
t130 = t83 * t137;
t44 = t87 * t130 + (t92 * t128 + t150 * t88) * t84;
t105 = g(1) * t32 + g(2) * t28 + g(3) * t44;
t159 = pkin(10) * t83;
t153 = t27 * pkin(11);
t31 = -t66 * t113 - t89 * t122 + t67 * t87;
t152 = t31 * pkin(11);
t143 = t84 * t92;
t145 = t84 * t88;
t43 = -t113 * t143 - t150 * t130 + t87 * t145;
t151 = t43 * pkin(11);
t147 = t83 * t86;
t146 = t83 * t91;
t141 = t85 * t86;
t140 = t86 * t90;
t133 = t83 * t145;
t139 = pkin(2) * t143 + pkin(10) * t133;
t138 = t93 * pkin(1) + pkin(9) * t144;
t135 = t65 * t159;
t134 = t67 * t159;
t131 = -pkin(1) * t89 + pkin(9) * t142;
t21 = t27 * pkin(3);
t125 = qJ(4) * t28 - t21;
t23 = t31 * pkin(3);
t124 = qJ(4) * t32 - t23;
t42 = t43 * pkin(3);
t123 = qJ(4) * t44 - t42;
t12 = -t31 * t91 + t49 * t86;
t119 = g(1) * t162 + g(2) * t12;
t36 = t65 * t113 - t64 * t87;
t37 = -t65 * t128 - t64 * t150;
t59 = t64 * pkin(2);
t118 = t37 * pkin(3) + t36 * qJ(4) - t59;
t38 = t67 * t113 + t66 * t87;
t39 = -t67 * t128 + t66 * t150;
t61 = t66 * pkin(2);
t117 = t39 * pkin(3) + t38 * qJ(4) + t61;
t116 = -g(1) * t27 + g(2) * t31;
t115 = -g(1) * t28 + g(2) * t32;
t114 = g(1) * t93 + g(2) * t89;
t57 = (t88 * t113 + t87 * t92) * t84;
t58 = -t128 * t145 + t92 * t132;
t112 = t58 * pkin(3) + t57 * qJ(4) + t139;
t110 = t37 * pkin(11) + t118;
t109 = t39 * pkin(11) + t117;
t63 = t137 * t136 - t83 * t143;
t25 = t43 * t91 - t63 * t86;
t108 = g(1) * t12 - g(2) * t162 - g(3) * t25;
t13 = t31 * t86 + t49 * t91;
t26 = t43 * t86 + t63 * t91;
t107 = g(1) * t13 + g(2) * t11 + g(3) * t26;
t14 = -t65 * t147 + t36 * t91;
t16 = -t67 * t147 + t38 * t91;
t40 = t86 * t133 - t57 * t91;
t106 = g(1) * t16 + g(2) * t14 - g(3) * t40;
t5 = g(1) * t31 + g(2) * t27 + g(3) * t43;
t104 = g(1) * t38 + g(2) * t36 + g(3) * t57;
t103 = g(1) * t39 + g(2) * t37 + g(3) * t58;
t102 = t67 * pkin(2) + t49 * pkin(10) + t138;
t101 = g(3) * t145 - t161;
t100 = -t65 * pkin(2) + t163 * pkin(10) + t131;
t99 = pkin(4) * t133 + t58 * pkin(11) + t112;
t98 = t161 * (pkin(4) + pkin(10)) * t83;
t97 = t32 * pkin(3) + qJ(4) * t31 + t102;
t96 = -pkin(3) * t28 - qJ(4) * t27 + t100;
t95 = t49 * pkin(4) + pkin(11) * t32 + t97;
t94 = pkin(4) * t163 - pkin(11) * t28 + t96;
t41 = t91 * t133 + t57 * t86;
t33 = t101 * t83;
t18 = -g(1) * t163 - g(2) * t49;
t17 = t67 * t146 + t38 * t86;
t15 = t65 * t146 + t36 * t86;
t3 = t13 * t90 + t32 * t85;
t2 = -t13 * t85 + t32 * t90;
t1 = t105 * t91;
t4 = [0, 0, 0, 0, 0, 0, g(1) * t89 - g(2) * t93, t114, 0, 0, 0, 0, 0, 0, 0, 0, g(1) * t65 - g(2) * t67, -g(1) * t64 - g(2) * t66, -t114 * t84, -g(1) * t131 - g(2) * t138, 0, 0, 0, 0, 0, 0, -t115, t116, t18, -g(1) * t100 - g(2) * t102, 0, 0, 0, 0, 0, 0, t18, t115, -t116, -g(1) * t96 - g(2) * t97, 0, 0, 0, 0, 0, 0, g(1) * t11 - g(2) * t13, t119, -t115, -g(1) * t94 - g(2) * t95, 0, 0, 0, 0, 0, 0, g(1) * t166 - g(2) * t3, -g(1) * t167 - g(2) * t2, -t119, -g(1) * (-pkin(5) * t11 + pkin(12) * t162 + t94) - g(2) * (pkin(5) * t13 + pkin(12) * t12 + t95); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -g(1) * t66 + g(2) * t64 - g(3) * t143, t101, 0, 0, 0, 0, 0, 0, 0, 0, -t103, t104, -t33, -g(1) * (t61 + t134) - g(2) * (-t59 + t135) - g(3) * t139, 0, 0, 0, 0, 0, 0, -t33, t103, -t104, -g(1) * (t117 + t134) - g(2) * (t118 + t135) - g(3) * t112, 0, 0, 0, 0, 0, 0, -g(1) * t17 - g(2) * t15 - g(3) * t41, -t106, -t103, -g(1) * t109 - g(2) * t110 - g(3) * t99 + t98, 0, 0, 0, 0, 0, 0, -g(1) * (t17 * t90 + t39 * t85) - g(2) * (t15 * t90 + t37 * t85) - g(3) * (t41 * t90 + t58 * t85) -g(1) * (-t17 * t85 + t39 * t90) - g(2) * (-t15 * t85 + t37 * t90) - g(3) * (-t41 * t85 + t58 * t90) t106, -g(1) * (t17 * pkin(5) - t16 * pkin(12) + t109) - g(2) * (t15 * pkin(5) - t14 * pkin(12) + t110) - g(3) * (t41 * pkin(5) + t40 * pkin(12) + t99) + t98; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t5, t105, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t5, -t105, -g(1) * t124 - g(2) * t125 - g(3) * t123, 0, 0, 0, 0, 0, 0, -t105 * t86, -t1, t5, -g(1) * (t124 - t152) - g(2) * (t125 - t153) - g(3) * (t123 - t151) 0, 0, 0, 0, 0, 0, -g(1) * (t32 * t140 - t31 * t85) - g(2) * (t28 * t140 - t27 * t85) - g(3) * (t44 * t140 - t43 * t85) -g(1) * (-t32 * t141 - t31 * t90) - g(2) * (-t28 * t141 - t27 * t90) - g(3) * (-t44 * t141 - t43 * t90) t1, -g(1) * (-t23 - t152) - g(2) * (-t21 - t153) - g(3) * (-t42 - t151) - t105 * (pkin(5) * t86 - pkin(12) * t91 + qJ(4)); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t5, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t5, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t5; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t108, t107, 0, 0, 0, 0, 0, 0, 0, 0, t108 * t90, -t108 * t85, -t107, -g(1) * (-pkin(5) * t12 + pkin(12) * t13) - g(2) * (pkin(5) * t162 + pkin(12) * t11) - g(3) * (pkin(5) * t25 + pkin(12) * t26); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -g(1) * t2 + g(2) * t167 - g(3) * (-t26 * t85 + t44 * t90) g(1) * t3 + g(2) * t166 - g(3) * (-t26 * t90 - t44 * t85) 0, 0;];
taug_reg  = t4;
