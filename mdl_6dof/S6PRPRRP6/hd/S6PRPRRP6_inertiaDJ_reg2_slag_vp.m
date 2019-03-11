% Calculate inertial parameters regressor of joint inertia matrix time derivative for
% S6PRPRRP6
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% qJD [6x1]
%   Generalized joint velocities
% pkin [10x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,alpha2,d2,d4,d5,theta1]';
% 
% Output:
% MMD_reg [((6+1)*6/2)x(6*10)]
%   inertial parameter regressor of inerta matrix time derivative
%   (only lower left triangular matrix (including diagonal) due to symmetry

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-08 20:21
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function MMD_reg = S6PRPRRP6_inertiaDJ_reg2_slag_vp(qJ, qJD, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(10,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6PRPRRP6_inertiaDJ_reg2_slag_vp: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6PRPRRP6_inertiaDJ_reg2_slag_vp: qJD has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S6PRPRRP6_inertiaDJ_reg2_slag_vp: pkin has to be [10x1] (double)');

%% Symbolic Calculation
% From inertiaD_joint_joint_fixb_regressor_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-08 20:20:56
% EndTime: 2019-03-08 20:21:03
% DurationCPUTime: 2.25s
% Computational Cost: add. (1271->212), mult. (3229->359), div. (0->0), fcn. (2791->8), ass. (0->142)
t82 = sin(qJ(5));
t76 = t82 ^ 2;
t84 = cos(qJ(5));
t78 = t84 ^ 2;
t156 = t76 + t78;
t85 = cos(qJ(4));
t73 = t85 * qJD(4);
t44 = t156 * t73;
t167 = sin(qJ(2));
t80 = sin(pkin(6));
t125 = t80 * t167;
t168 = cos(qJ(2));
t126 = t80 * t168;
t83 = sin(qJ(4));
t116 = t83 * t126;
t81 = cos(pkin(6));
t43 = t81 * t85 - t116;
t100 = t84 * t125 - t43 * t82;
t24 = t82 * t125 + t43 * t84;
t122 = qJD(2) * t168;
t115 = t80 * t122;
t42 = t85 * t126 + t81 * t83;
t121 = qJD(2) * t167;
t69 = t80 * t121;
t21 = qJD(4) * t42 - t83 * t69;
t8 = t100 * qJD(5) + t82 * t115 - t21 * t84;
t6 = t8 * t84;
t151 = qJD(5) * t24;
t7 = -t84 * t115 - t21 * t82 + t151;
t95 = t6 + t7 * t82 + (-t100 * t84 - t24 * t82) * qJD(5);
t157 = t76 - t78;
t147 = t83 * qJD(4);
t22 = -qJD(4) * t116 - t85 * t69 + t81 * t73;
t166 = t22 * t85;
t176 = t42 * t147 - t166;
t77 = t83 ^ 2;
t79 = t85 ^ 2;
t119 = (t77 - t79) * qJD(4);
t58 = t157 * qJD(5);
t169 = pkin(9) * t85;
t113 = pkin(4) * t83 - t169;
t104 = qJ(3) + t113;
t86 = -pkin(2) - pkin(8);
t163 = t83 * t86;
t62 = t84 * t163;
t159 = t82 * t104 + t62;
t170 = pkin(9) * t83;
t114 = pkin(4) * t85 + t170;
t98 = t114 * qJD(4) + qJD(3);
t89 = -qJD(5) * t159 + t84 * t98;
t110 = pkin(5) * t82 - qJ(6) * t84;
t39 = t110 * qJD(5) - t82 * qJD(6);
t162 = t85 * t39;
t103 = t110 - t86;
t38 = t103 * t85;
t111 = pkin(5) * t84 + qJ(6) * t82;
t53 = -pkin(4) - t111;
t175 = qJD(4) * (t53 * t83 + t169) - qJD(5) * t38 - t162;
t174 = t111 * qJD(5) - t84 * qJD(6);
t173 = 0.2e1 * qJD(3);
t172 = 0.2e1 * qJD(6);
t171 = pkin(4) * t84;
t16 = t42 * t22;
t165 = t42 * t83;
t164 = t82 * t86;
t161 = t85 * t86;
t160 = qJ(3) * t115 + qJD(3) * t125;
t158 = pkin(9) * t44;
t154 = t77 + t79;
t153 = qJD(4) * t38;
t152 = qJD(4) * t84;
t150 = qJD(5) * t82;
t74 = qJD(5) * t84;
t149 = qJD(5) * t85;
t148 = qJD(5) * t86;
t145 = qJ(3) * qJD(4);
t144 = qJ(6) * qJD(4);
t143 = -0.2e1 * pkin(4) * qJD(5);
t142 = t82 * t163;
t101 = t84 * t104;
t127 = t86 * t73;
t141 = qJD(5) * t101 + t84 * t127 + t82 * t98;
t139 = pkin(9) * t150;
t138 = pkin(9) * t74;
t137 = t84 * t147;
t136 = t82 * t149;
t135 = t82 * t148;
t134 = t83 * t74;
t133 = t84 * t149;
t132 = t82 * t147;
t131 = t82 * t74;
t130 = t86 * t147;
t129 = t84 * t73;
t128 = t83 * t73;
t124 = t85 * t144;
t123 = -pkin(5) + t164;
t120 = qJD(4) * t167;
t67 = 0.2e1 * t128;
t118 = t84 * t132;
t117 = t79 * t131;
t27 = qJ(6) * t83 + t159;
t28 = -t84 * (qJ(3) - t169) + (t123 - t171) * t83;
t109 = t27 * t84 + t28 * t82;
t108 = t27 * t82 - t28 * t84;
t35 = t101 - t142;
t107 = t159 * t82 + t35 * t84;
t106 = -t159 * t84 + t35 * t82;
t12 = t42 * t150 - t22 * t84;
t102 = -0.2e1 * t100 * t7 + 0.2e1 * t24 * t8 + 0.2e1 * t16;
t15 = -t103 * t147 + t174 * t85;
t99 = -t15 + (t53 * t85 - t170) * qJD(5);
t10 = t123 * t73 - t89;
t9 = t124 + (qJD(6) - t135) * t83 + t141;
t93 = -t108 * qJD(5) + t10 * t82 + t9 * t84;
t13 = t83 * t135 - t141;
t14 = -t82 * t127 + t89;
t92 = -t107 * qJD(5) - t13 * t84 - t14 * t82;
t2 = -t21 * t83 - t166 + (t43 * t85 + t165) * qJD(4);
t91 = t95 * pkin(9);
t90 = -t7 * t83 + t82 * t166 + t42 * t133 + (t100 * t85 - t82 * t165) * qJD(4);
t88 = t100 * t137 + t24 * t132 + (t7 * t84 - t8 * t82 + (t100 * t82 - t24 * t84) * qJD(5)) * t85;
t87 = -t100 * t134 + t24 * t129 + t83 * t6 + (-t100 * t73 + (t7 - t151) * t83) * t82 + t176;
t75 = qJ(3) * t173;
t66 = -0.2e1 * t131;
t65 = 0.2e1 * t131;
t50 = -t132 + t133;
t49 = t82 * t73 + t134;
t48 = t154 * t74;
t47 = -t136 - t137;
t46 = t83 * t150 - t129;
t45 = t154 * t150;
t34 = -0.2e1 * t78 * t128 - 0.2e1 * t117;
t33 = -0.2e1 * t76 * t128 + 0.2e1 * t117;
t31 = t157 * t149 + t118;
t30 = -t82 * t119 + t83 * t133;
t29 = 0.4e1 * t85 * t131 - t157 * t147;
t26 = -0.2e1 * t84 * t119 - 0.2e1 * t83 * t136;
t25 = 0.2e1 * t85 * t118 + t79 * t58;
t20 = (-0.1e1 + t156) * t67;
t11 = t22 * t82 + t42 * t74;
t1 = (t42 * t152 + t8) * t83 + (qJD(4) * t24 + t12) * t85;
t3 = [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0.2e1 * t121 * t168 * t80 ^ 2 - 0.2e1 * t43 * t21 + 0.2e1 * t16, 0, 0, 0, 0, 0, 0, 0, 0, 0, t102, 0, 0, 0, 0, 0, 0, 0, 0, 0, t102; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t69, -t115, 0, 0, 0, 0, 0, 0, 0, 0, 0, t69, t115, -pkin(2) * t69 + t160, 0, 0, 0, 0, 0, 0 (t120 * t85 + t122 * t83) * t80 (-t120 * t83 + t122 * t85) * t80, -t2, t2 * t86 + t160, 0, 0, 0, 0, 0, 0, t90, -t1, t88, t100 * t14 - t24 * t13 + t8 * t159 + t176 * t86 - t7 * t35, 0, 0, 0, 0, 0, 0, t90, t88, t1, -t10 * t100 + t15 * t42 + t22 * t38 + t24 * t9 + t27 * t8 + t28 * t7; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t173, t75, -0.2e1 * t128, 0.2e1 * t119, 0, t67, 0, 0, 0.2e1 * qJD(3) * t83 + 0.2e1 * t145 * t85, 0.2e1 * qJD(3) * t85 - 0.2e1 * t145 * t83, 0, t75, t34, 0.2e1 * t25, t26, t33, -0.2e1 * t30, t67, -0.2e1 * t79 * t84 * t148 + 0.2e1 * t14 * t83 + 0.2e1 * (t35 + 0.2e1 * t142) * t73, 0.2e1 * t79 * t135 + 0.2e1 * t13 * t83 + 0.2e1 * (-t159 + 0.2e1 * t62) * t73, 0.2e1 * t107 * t147 + 0.2e1 * (qJD(5) * t106 + t13 * t82 - t14 * t84) * t85, -0.2e1 * t128 * t86 ^ 2 - 0.2e1 * t13 * t159 + 0.2e1 * t35 * t14, t34, t26, -0.2e1 * t25, t67, 0.2e1 * t30, t33, 0.2e1 * (-t153 * t82 - t10) * t83 + 0.2e1 * (-qJD(4) * t28 + t15 * t82 + t38 * t74) * t85, 0.2e1 * t108 * t147 + 0.2e1 * (-qJD(5) * t109 + t10 * t84 - t82 * t9) * t85, 0.2e1 * (t152 * t38 + t9) * t83 + 0.2e1 * (qJD(4) * t27 - t15 * t84 + t150 * t38) * t85, 0.2e1 * t10 * t28 + 0.2e1 * t15 * t38 + 0.2e1 * t27 * t9; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t69, 0, 0, 0, 0, 0, 0, 0, 0, 0, t2, 0, 0, 0, 0, 0, 0, 0, 0, 0, t87, 0, 0, 0, 0, 0, 0, 0, 0, 0, t87; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t48, t45, 0, -t106 * t73 + (t92 - 0.2e1 * t127) * t83, 0, 0, 0, 0, 0, 0, -t48, 0, -t45 (qJD(4) * t109 - t15) * t85 + (t93 + t153) * t83; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t20, 0, 0, 0, 0, 0, 0, 0, 0, 0, t20; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t22, t21, 0, 0, 0, 0, 0, 0, 0, 0, t12, t11, t95, -t22 * pkin(4) + t91, 0, 0, 0, 0, 0, 0, t12, t95, -t11, t22 * t53 + t42 * t39 + t91; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t147, 0, -t73, 0, -t130, -t127, 0, 0, -t31, -t29, t49, t31, -t46, 0 (-t114 * t84 - t161 * t82) * qJD(5) + (t113 * t82 - t62) * qJD(4) (t114 * t82 - t161 * t84) * qJD(5) + (-t84 * t169 + (t164 + t171) * t83) * qJD(4), t92, -pkin(4) * t130 + pkin(9) * t92, -t31, t49, t29, 0, t46, t31, -t175 * t82 + t99 * t84, t93, t175 * t84 + t99 * t82, pkin(9) * t93 + t15 * t53 + t38 * t39; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t147, -t73, 0, 0, 0, 0, 0, 0, 0, 0, t47, -t50, t44, -pkin(4) * t147 + t158, 0, 0, 0, 0, 0, 0, t47, t44, t50, t147 * t53 + t158 - t162; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t65, -0.2e1 * t58, 0, t66, 0, 0, t82 * t143, t84 * t143, 0, 0, t65, 0, 0.2e1 * t58, 0, 0, t66, 0.2e1 * t150 * t53 - 0.2e1 * t39 * t84, 0, -0.2e1 * t39 * t82 - 0.2e1 * t53 * t74, 0.2e1 * t53 * t39; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t7, -t8, 0, 0, 0, 0, 0, 0, 0, 0, -t7, 0, t8, -pkin(5) * t7 + qJ(6) * t8 + qJD(6) * t24; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t47, 0, -t50, t73, t14, t13, 0, 0, 0, t47, 0, t73, t50, 0 (0.2e1 * pkin(5) - t164) * t73 + t89 (pkin(5) * t147 - qJ(6) * t149) * t84 + (t83 * t144 + (pkin(5) * qJD(5) - qJD(6)) * t85) * t82, 0.2e1 * t124 + (t172 - t135) * t83 + t141, -pkin(5) * t10 + qJ(6) * t9 + qJD(6) * t27; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t49, t46, 0, 0, 0, 0, 0, 0, 0, 0, -t49, 0, -t46, -t110 * t73 - t174 * t83; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t74, 0, -t150, 0, -t138, t139, 0, 0, 0, t74, 0, 0, t150, 0, -t138, -t174, -t139, -t174 * pkin(9); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t172, qJ(6) * t172; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t7; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t73, t47, 0, t10; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t49; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t74, 0, t138; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0;];
MMD_reg  = t3;
