% Calculate minimal parameter regressor of joint inertia matrix time derivative for
% S5RRRRP8
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% qJD [5x1]
%   Generalized joint velocities
% pkin [8x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d1,d2,d3,d4]';
% 
% Output:
% MMD_reg [((5+1)*5/2)x28]
%   minimal parameter regressor of inertia matrix time derivative
%   (only lower left triangular matrix (including diagonal) due to symmetry

% Quelle: HybrDyn-Toolbox
% Datum: 2021-01-16 00:22
% Revision: 24b2e7d74a0c1a3b64fa2f8f5ad758691ad61af3 (2021-01-15)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function MMD_reg = S5RRRRP8_inertiaDJ_regmin_slag_vp(qJ, qJD, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(8,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RRRRP8_inertiaDJ_regmin_slag_vp: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RRRRP8_inertiaDJ_regmin_slag_vp: qJD has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S5RRRRP8_inertiaDJ_regmin_slag_vp: pkin has to be [8x1] (double)');

%% Symbolic Calculation
% From inertiaD_joint_joint_fixb_regressor_minpar_matlab.m
% OptimizationMode: 2
% StartTime: 2021-01-16 00:21:20
% EndTime: 2021-01-16 00:21:27
% DurationCPUTime: 1.59s
% Computational Cost: add. (1890->202), mult. (4772->373), div. (0->0), fcn. (4013->6), ass. (0->112)
t144 = pkin(7) + pkin(8);
t85 = cos(qJ(3));
t148 = t144 * t85;
t82 = sin(qJ(4));
t152 = t148 * t82;
t86 = cos(qJ(2));
t126 = t86 * qJD(2);
t83 = sin(qJ(3));
t113 = t83 * t126;
t129 = qJD(3) * t85;
t84 = sin(qJ(2));
t151 = t84 * t129 + t113;
t150 = -0.4e1 * t84;
t130 = qJD(3) * t83;
t117 = t84 * t130;
t139 = t83 * t84;
t123 = t82 * t139;
t141 = cos(qJ(4));
t103 = t141 * qJD(4);
t147 = t141 * qJD(3) + t103;
t105 = qJD(2) * t141;
t96 = t86 * t105;
t16 = t83 * t96 - t82 * t117 - qJD(4) * t123 + (t82 * t126 + t147 * t84) * t85;
t49 = t141 * t83 + t82 * t85;
t37 = t49 * t84;
t149 = t16 * qJ(5) + t37 * qJD(5);
t143 = pkin(6) * t83;
t109 = -pkin(3) - t143;
t138 = t84 * t85;
t95 = -t86 * pkin(2) - t84 * pkin(7);
t54 = -pkin(1) + t95;
t47 = t85 * t54;
t28 = -pkin(8) * t138 + t109 * t86 + t47;
t137 = t85 * t86;
t68 = pkin(6) * t137;
t133 = t83 * t54 + t68;
t33 = -pkin(8) * t139 + t133;
t136 = t141 * t33 + t82 * t28;
t79 = t84 ^ 2;
t101 = (-t86 ^ 2 + t79) * qJD(2);
t80 = t85 ^ 2;
t132 = t83 ^ 2 - t80;
t102 = t132 * qJD(3);
t146 = qJD(3) + qJD(4);
t128 = qJD(3) * t86;
t116 = t83 * t128;
t76 = t84 * qJD(2);
t91 = t85 * t76 + t116;
t94 = pkin(2) * t84 - pkin(7) * t86;
t92 = t94 * t83;
t19 = t91 * pkin(6) - qJD(2) * t92 - t54 * t129;
t140 = t82 * t83;
t145 = t146 * t140;
t142 = t85 * pkin(2);
t134 = -t144 * t140 + t141 * t148;
t53 = pkin(3) * t139 + t84 * pkin(6);
t127 = qJD(4) * t82;
t125 = -0.2e1 * pkin(1) * qJD(2);
t124 = -0.2e1 * pkin(2) * qJD(3);
t13 = -pkin(8) * t151 - t19;
t87 = (-t68 + (pkin(8) * t84 - t54) * t83) * qJD(3) + (-t86 * t148 + (-t109 + t142) * t84) * qJD(2);
t122 = -t28 * t103 - t141 * t13 - t82 * t87;
t74 = pkin(6) * t126;
t34 = t151 * pkin(3) + t74;
t121 = pkin(3) * t130;
t120 = pkin(3) * t127;
t118 = t141 * pkin(3);
t114 = t85 * t128;
t112 = t83 * t129;
t111 = t84 * t126;
t110 = t85 * t126;
t73 = -t85 * pkin(3) - pkin(2);
t108 = t141 * t85;
t107 = t141 * t28 - t82 * t33;
t97 = t144 * t141;
t51 = t83 * t97;
t106 = -t51 - t152;
t100 = 0.2e1 * t111;
t99 = t83 * t110;
t98 = pkin(3) * t103;
t93 = qJD(3) * t97;
t3 = t33 * t127 + t122;
t17 = qJD(3) * t152 + qJD(4) * t51 + t127 * t148 + t83 * t93;
t4 = -t136 * qJD(4) - t82 * t13 + t141 * t87;
t89 = t86 * t98 + (-pkin(3) * t76 + qJD(4) * t33) * t82 + t122;
t18 = -t103 * t148 + t144 * t145 - t85 * t93;
t32 = t146 * t49;
t15 = t82 * t113 + t32 * t84 - t85 * t96;
t38 = t84 * t108 - t123;
t88 = t15 * qJ(5) - t38 * qJD(5) + t4;
t75 = pkin(4) * t76;
t1 = t75 + t88;
t72 = t118 + pkin(4);
t70 = -0.2e1 * t98;
t69 = -0.2e1 * t120;
t64 = -0.2e1 * t111;
t61 = t86 * t120;
t48 = -t108 + t140;
t35 = t48 * pkin(4) + t73;
t31 = -t147 * t85 + t145;
t29 = t37 * pkin(4) + t53;
t23 = t32 * pkin(4) + t121;
t22 = -t48 * qJ(5) + t134;
t21 = -t49 * qJ(5) + t106;
t20 = -t133 * qJD(3) + (pkin(6) * t139 + t85 * t94) * qJD(2);
t9 = t16 * pkin(4) + t34;
t8 = -t37 * qJ(5) + t136;
t7 = -t86 * pkin(4) - t38 * qJ(5) + t107;
t6 = t31 * qJ(5) - t49 * qJD(5) + t18;
t5 = t32 * qJ(5) + t48 * qJD(5) + t17;
t2 = t3 + t149;
t10 = [0, 0, 0, t100, -0.2e1 * t101, 0, 0, 0, t84 * t125, t86 * t125, 0.2e1 * t80 * t111 - 0.2e1 * t79 * t112, 0.2e1 * t79 * t102 + t99 * t150, 0.2e1 * t85 * t101 + 0.2e1 * t84 * t116, -0.2e1 * t83 * t101 + 0.2e1 * t84 * t114, t64, 0.2e1 * t47 * t76 - 0.2e1 * t20 * t86 + 0.2e1 * (t83 * t111 + t79 * t129) * pkin(6), -0.2e1 * t19 * t86 - 0.2e1 * t133 * t76 + 0.2e1 * (t85 * t100 - t79 * t130) * pkin(6), -0.2e1 * t38 * t15, 0.2e1 * t15 * t37 - 0.2e1 * t38 * t16, 0.2e1 * t15 * t86 + 0.2e1 * t38 * t76, 0.2e1 * t16 * t86 - 0.2e1 * t37 * t76, t64, 0.2e1 * t107 * t76 + 0.2e1 * t53 * t16 + 0.2e1 * t34 * t37 - 0.2e1 * t4 * t86, -0.2e1 * t136 * t76 - 0.2e1 * t53 * t15 - 0.2e1 * t3 * t86 + 0.2e1 * t34 * t38, -0.2e1 * t1 * t86 + 0.2e1 * t29 * t16 + 0.2e1 * t9 * t37 + 0.2e1 * t7 * t76, -0.2e1 * t29 * t15 - 0.2e1 * t2 * t86 + 0.2e1 * t9 * t38 - 0.2e1 * t8 * t76, -0.2e1 * t1 * t38 + 0.2e1 * t7 * t15 - 0.2e1 * t8 * t16 + 0.2e1 * t2 * t37, 0.2e1 * t7 * t1 - 0.2e1 * t8 * t2 + 0.2e1 * t29 * t9; 0, 0, 0, 0, 0, t126, -t76, 0, -t74, pkin(6) * t76, -t84 * t102 + t99, t112 * t150 - t132 * t126, t83 * t76 - t114, t91, 0, (pkin(7) * t137 + (-t142 + t143) * t84) * qJD(3) + (t95 * t83 - t68) * qJD(2), (pkin(6) * t138 + t92) * qJD(3) + (t86 * t143 + t95 * t85) * qJD(2), -t15 * t49 - t38 * t31, t15 * t48 - t49 * t16 + t31 * t37 - t38 * t32, t31 * t86 + t49 * t76, t32 * t86 - t48 * t76, 0, t106 * t76 + t37 * t121 + t73 * t16 - t18 * t86 + t53 * t32 + t34 * t48, t38 * t121 - t134 * t76 - t73 * t15 - t17 * t86 - t53 * t31 + t34 * t49, t35 * t16 + t21 * t76 + t23 * t37 + t29 * t32 + t9 * t48 - t6 * t86, -t35 * t15 - t22 * t76 + t23 * t38 - t29 * t31 + t9 * t49 - t5 * t86, -t1 * t49 + t21 * t15 - t22 * t16 + t2 * t48 + t7 * t31 - t8 * t32 + t5 * t37 - t6 * t38, t1 * t21 - t2 * t22 + t29 * t23 + t9 * t35 - t8 * t5 + t7 * t6; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0.2e1 * t112, -0.2e1 * t102, 0, 0, 0, t83 * t124, t85 * t124, -0.2e1 * t49 * t31, 0.2e1 * t31 * t48 - 0.2e1 * t49 * t32, 0, 0, 0, 0.2e1 * t48 * t121 + 0.2e1 * t73 * t32, 0.2e1 * t49 * t121 - 0.2e1 * t73 * t31, 0.2e1 * t23 * t48 + 0.2e1 * t35 * t32, 0.2e1 * t23 * t49 - 0.2e1 * t35 * t31, 0.2e1 * t21 * t31 - 0.2e1 * t22 * t32 + 0.2e1 * t5 * t48 - 0.2e1 * t6 * t49, 0.2e1 * t21 * t6 - 0.2e1 * t22 * t5 + 0.2e1 * t35 * t23; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t110 - t117, -t151, t76, t20, t19, 0, 0, -t15, -t16, t76, pkin(3) * t105 * t84 + t4 + t61, t89, t72 * t76 + t1 + t61, t89 + t149, t72 * t15 + (-t16 * t82 + (-t141 * t37 + t38 * t82) * qJD(4)) * pkin(3), t1 * t72 + (-t2 * t82 + (t141 * t8 - t7 * t82) * qJD(4)) * pkin(3); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t129, -t130, 0, -pkin(7) * t129, pkin(7) * t130, 0, 0, -t31, -t32, 0, t18, t17, t6, t5, t72 * t31 + (-t32 * t82 + (-t141 * t48 + t49 * t82) * qJD(4)) * pkin(3), t6 * t72 + (-t5 * t82 + (t141 * t22 - t21 * t82) * qJD(4)) * pkin(3); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t69, t70, t69, t70, 0, 0.2e1 * (t118 - t72) * t120; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t15, -t16, t76, t4, t3, 0.2e1 * t75 + t88, t2, t15 * pkin(4), t1 * pkin(4); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t31, -t32, 0, t18, t17, t6, t5, t31 * pkin(4), t6 * pkin(4); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t120, -t98, -t120, -t98, 0, -pkin(4) * t120; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t16, -t15, 0, t9; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t32, -t31, 0, t23; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0;];
MMD_reg = t10;
