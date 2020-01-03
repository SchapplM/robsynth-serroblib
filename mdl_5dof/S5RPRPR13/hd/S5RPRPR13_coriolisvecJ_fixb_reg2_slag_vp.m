% Calculate inertial parameters regressor of coriolis joint torque vector for
% S5RPRPR13
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% qJD [5x1]
%   Generalized joint velocities
% pkin [8x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d1,d3,d5,theta2]';
% 
% Output:
% tauc_reg [5x(5*10)]
%   inertial parameter regressor of coriolis joint torque vector

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-31 18:33
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauc_reg = S5RPRPR13_coriolisvecJ_fixb_reg2_slag_vp(qJ, qJD, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(8,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RPRPR13_coriolisvecJ_fixb_reg2_slag_vp: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RPRPR13_coriolisvecJ_fixb_reg2_slag_vp: qJD has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S5RPRPR13_coriolisvecJ_fixb_reg2_slag_vp: pkin has to be [8x1] (double)');

%% Symbolic Calculation
% From coriolisvec_joint_fixb_regressor_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 18:32:53
% EndTime: 2019-12-31 18:32:59
% DurationCPUTime: 1.75s
% Computational Cost: add. (2776->256), mult. (7190->331), div. (0->0), fcn. (5192->6), ass. (0->139)
t165 = cos(qJ(3));
t95 = sin(pkin(8));
t96 = cos(pkin(8));
t98 = sin(qJ(3));
t77 = t165 * t95 + t98 * t96;
t172 = t77 * qJD(1);
t177 = qJD(5) + t172;
t132 = t165 * t96;
t122 = qJD(1) * t132;
t153 = t95 * t98;
t134 = qJD(1) * t153;
t68 = -t122 + t134;
t97 = sin(qJ(5));
t99 = cos(qJ(5));
t48 = qJD(3) * t97 - t99 * t68;
t125 = t177 * t48;
t140 = qJD(5) * t99;
t141 = qJD(5) * t97;
t73 = t77 * qJD(3);
t57 = qJD(1) * t73;
t31 = qJD(3) * t141 - t68 * t140 - t97 * t57;
t179 = t31 - t125;
t142 = qJD(3) * t98;
t133 = t95 * t142;
t87 = qJD(3) * t122;
t56 = qJD(1) * t133 - t87;
t149 = qJ(4) * t56;
t115 = -qJD(4) * t172 + t149;
t168 = pkin(3) + pkin(7);
t11 = t168 * t57 + t115;
t130 = qJD(2) * t165;
t120 = qJD(1) * t130;
t129 = qJD(3) * t165;
t135 = qJD(1) * qJD(2);
t131 = t98 * t135;
t152 = pkin(6) + qJ(2);
t81 = t152 * t95;
t78 = qJD(1) * t81;
t82 = t152 * t96;
t79 = qJD(1) * t82;
t23 = t95 * t120 + t79 * t129 + t96 * t131 - t78 * t142;
t14 = -pkin(4) * t56 + t23;
t92 = -pkin(2) * t96 - pkin(1);
t80 = t92 * qJD(1) + qJD(2);
t104 = -qJ(4) * t172 + t80;
t16 = t168 * t68 + t104;
t43 = t165 * t78 + t79 * t98;
t139 = -qJD(4) - t43;
t138 = pkin(4) * t172 - t139;
t20 = -t168 * qJD(3) + t138;
t6 = t16 * t99 + t20 * t97;
t2 = -qJD(5) * t6 - t11 * t97 + t99 * t14;
t171 = t177 * t6 + t2;
t117 = t16 * t97 - t20 * t99;
t1 = -t117 * qJD(5) + t11 * t99 + t14 * t97;
t119 = t117 * t177 + t1;
t127 = t97 * t177;
t52 = t99 * t56;
t110 = -t127 * t177 - t52;
t178 = t172 * qJD(3);
t169 = t68 ^ 2;
t66 = t172 ^ 2;
t175 = -t169 - t66;
t174 = -t169 + t66;
t45 = t56 * t77;
t72 = -t96 * t129 + t133;
t173 = -t172 * t72 - t45;
t34 = (qJD(2) * t95 + qJD(3) * t82) * t98 + t81 * t129 - t96 * t130;
t166 = t68 * pkin(4);
t44 = t165 * t79 - t98 * t78;
t40 = -qJD(3) * qJ(4) - t44;
t24 = -t40 - t166;
t170 = t168 * t56 + t177 * t24;
t167 = pkin(3) * t57;
t123 = -t96 * t120 + t78 * t129 + t95 * t131 + t79 * t142;
t21 = -qJD(3) * qJD(4) + t123;
t12 = -pkin(4) * t57 - t21;
t164 = t12 * t99;
t46 = t165 * t81 + t82 * t98;
t163 = t23 * t46;
t162 = t31 * t99;
t50 = qJD(3) * t99 + t68 * t97;
t53 = t99 * t57;
t32 = qJD(5) * t50 - t53;
t161 = t32 * t97;
t33 = pkin(3) * t68 + t104;
t160 = t33 * t172;
t159 = t48 * t68;
t158 = t50 * t48;
t157 = t50 * t68;
t156 = t68 * t172;
t76 = -t132 + t153;
t154 = t76 * t97;
t26 = t99 * t32;
t150 = t95 ^ 2 + t96 ^ 2;
t148 = qJ(4) * t68;
t146 = qJD(3) * t34;
t47 = t165 * t82 - t98 * t81;
t35 = qJD(2) * t77 + qJD(3) * t47;
t145 = qJD(3) * t35;
t144 = qJD(3) * t72;
t143 = qJD(3) * t73;
t136 = qJD(5) * t168;
t128 = t150 * qJD(1) ^ 2;
t126 = t99 * t177;
t124 = t177 * t50;
t121 = t117 * t97 + t6 * t99;
t111 = -qJ(4) * t77 + t92;
t28 = t168 * t76 + t111;
t36 = pkin(4) * t77 + t46;
t10 = t28 * t99 + t36 * t97;
t9 = -t28 * t97 + t36 * t99;
t116 = t57 * t76 + t68 * t73;
t114 = qJ(4) * t72 - qJD(4) * t77;
t113 = 0.2e1 * t150 * t135;
t109 = t76 * t140 + t73 * t97;
t108 = t76 * t141 - t73 * t99;
t107 = qJD(3) * t44 - t23;
t106 = -t126 * t177 + t56 * t97;
t103 = -t172 * t73 + t56 * t76 - t57 * t77 + t68 * t72;
t102 = t172 * t35 + t23 * t77 + t34 * t68 - t46 * t56 - t47 * t57;
t58 = qJD(3) * t68;
t42 = pkin(3) * t76 + t111;
t41 = pkin(3) * t172 + t148;
t39 = t56 - t58;
t38 = -qJD(3) * pkin(3) - t139;
t37 = -t76 * pkin(4) + t47;
t30 = t44 - t166;
t27 = pkin(3) * t73 + t114;
t25 = t168 * t172 + t148;
t19 = t115 + t167;
t18 = -t72 * pkin(4) + t35;
t17 = -pkin(4) * t73 - t34;
t15 = t168 * t73 + t114;
t8 = t25 * t99 + t30 * t97;
t7 = -t25 * t97 + t30 * t99;
t4 = -qJD(5) * t10 - t15 * t97 + t18 * t99;
t3 = qJD(5) * t9 + t15 * t99 + t18 * t97;
t5 = [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t113, qJ(2) * t113, t173, t103, -t144, t116, -t143, 0, t57 * t92 + t73 * t80 - t145, -t56 * t92 - t72 * t80 + t146, t123 * t76 - t43 * t72 - t44 * t73 + t102, -t123 * t47 - t34 * t44 + t35 * t43 + t163, 0, t144, t143, t173, t103, t116, t21 * t76 - t38 * t72 + t40 * t73 + t102, -t19 * t76 - t27 * t68 - t33 * t73 - t42 * t57 + t145, -t172 * t27 - t19 * t77 + t33 * t72 + t42 * t56 - t146, t19 * t42 - t21 * t47 + t27 * t33 + t34 * t40 + t35 * t38 + t163, t109 * t50 - t31 * t154, (-t48 * t97 + t50 * t99) * t73 + (-t162 - t161 + (-t48 * t99 - t50 * t97) * qJD(5)) * t76, t109 * t177 - t56 * t154 - t31 * t77 - t50 * t72, t108 * t48 - t76 * t26, -t108 * t177 - t32 * t77 + t48 * t72 - t76 * t52, -t177 * t72 - t45, t108 * t24 + t117 * t72 - t76 * t164 + t17 * t48 + t177 * t4 + t2 * t77 + t32 * t37 - t56 * t9, -t1 * t77 + t10 * t56 + t109 * t24 + t12 * t154 + t17 * t50 - t177 * t3 - t31 * t37 + t6 * t72, -t10 * t32 - t3 * t48 + t31 * t9 - t4 * t50 + t121 * t73 + (t1 * t99 - t2 * t97 + (t117 * t99 - t6 * t97) * qJD(5)) * t76, t1 * t10 - t117 * t4 + t12 * t37 + t17 * t24 + t2 * t9 + t3 * t6; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t128, -qJ(2) * t128, 0, 0, 0, 0, 0, 0, 0.2e1 * t178, t87 + (-t68 - t134) * qJD(3), t175, -t172 * t43 + t44 * t68, 0, 0, 0, 0, 0, 0, t175, -0.2e1 * t178, t56 + t58, t167 + t149 - t40 * t68 + (-qJD(4) - t38) * t172, 0, 0, 0, 0, 0, 0, t106 + t159, t157 - t110, t99 * t124 - t179 * t97 - t26, t119 * t99 - t171 * t97 + t24 * t68; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t156, t174, t87 + (t68 - t134) * qJD(3), -t156, 0, 0, -t172 * t80 + t107, -qJD(3) * t43 + t68 * t80 + t123, 0, 0, 0, t39, 0, t156, t174, -t156, pkin(3) * t56 - qJ(4) * t57 + (-t40 - t44) * t172 + (t38 + t139) * t68, t41 * t68 - t107 + t160, -t33 * t68 + t41 * t172 + (0.2e1 * qJD(4) + t43) * qJD(3) - t123, -pkin(3) * t23 - qJ(4) * t21 + t139 * t40 - t33 * t41 - t38 * t44, -t50 * t127 - t162, -t26 - t50 * t126 + (t31 + t125) * t97, t110 + t157, t99 * t125 + t161, t106 - t159, t177 * t68, qJ(4) * t32 + t12 * t97 - t117 * t68 + (t97 * t136 - t7) * t177 + t138 * t48 + t170 * t99, -qJ(4) * t31 + t164 - t6 * t68 + (t136 * t99 + t8) * t177 + t138 * t50 - t170 * t97, t48 * t8 + t50 * t7 + (-t168 * t31 - t6 * t172 - t2 + (t168 * t48 - t6) * qJD(5)) * t99 + (t168 * t32 - t117 * t172 - t1 + (-t168 * t50 - t117) * qJD(5)) * t97, qJ(4) * t12 + t117 * t7 - t6 * t8 + t138 * t24 - (qJD(5) * t121 + t1 * t97 + t2 * t99) * t168; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t39, -t156, -qJD(3) ^ 2 - t66, qJD(3) * t40 + t160 + t23, 0, 0, 0, 0, 0, 0, -qJD(3) * t48 + t110, -qJD(3) * t50 + t106, t179 * t99 + (-t32 + t124) * t97, -qJD(3) * t24 + t119 * t97 + t171 * t99; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t158, -t48 ^ 2 + t50 ^ 2, -t179, -t158, t53 + (-qJD(5) + t177) * t50, -t56, -t24 * t50 + t171, t24 * t48 - t119, 0, 0;];
tauc_reg = t5;
