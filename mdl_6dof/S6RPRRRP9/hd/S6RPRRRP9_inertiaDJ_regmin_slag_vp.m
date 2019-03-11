% Calculate minimal parameter regressor of joint inertia matrix time derivative for
% S6RPRRRP9
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% qJD [6x1]
%   Generalized joint velocities
% pkin [9x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d3,d4,d5]';
% 
% Output:
% MMD_reg [((6+1)*6/2)x29]
%   minimal parameter regressor of inerta matrix time derivative
%   (only lower left triangular matrix (including diagonal) due to symmetry

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-09 06:29
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function MMD_reg = S6RPRRRP9_inertiaDJ_regmin_slag_vp(qJ, qJD, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(9,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPRRRP9_inertiaDJ_regmin_slag_vp: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RPRRRP9_inertiaDJ_regmin_slag_vp: qJD has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S6RPRRRP9_inertiaDJ_regmin_slag_vp: pkin has to be [9x1] (double)');

%% Symbolic Calculation
% From inertiaD_joint_joint_fixb_regressor_minpar_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 06:28:31
% EndTime: 2019-03-09 06:28:36
% DurationCPUTime: 1.61s
% Computational Cost: add. (1963->220), mult. (4521->394), div. (0->0), fcn. (3894->6), ass. (0->120)
t91 = sin(qJ(3));
t126 = t91 * qJD(3);
t90 = sin(qJ(4));
t113 = t90 * t126;
t94 = cos(qJ(3));
t130 = qJD(4) * t94;
t93 = cos(qJ(4));
t115 = t93 * t130;
t46 = t113 - t115;
t89 = sin(qJ(5));
t144 = t89 * t90;
t92 = cos(qJ(5));
t56 = -t92 * t93 + t144;
t39 = t56 * t91;
t95 = -pkin(1) - pkin(7);
t142 = t90 * t95;
t106 = pkin(4) - t142;
t140 = t93 * t94;
t96 = t91 * pkin(3) - t94 * pkin(8);
t63 = qJ(2) + t96;
t53 = t93 * t63;
t29 = -pkin(9) * t140 + t106 * t91 + t53;
t143 = t90 * t94;
t52 = t90 * t63;
t141 = t91 * t93;
t74 = t95 * t141;
t149 = -t74 - t52;
t34 = -pkin(9) * t143 - t149;
t30 = t92 * t34;
t137 = t89 * t29 + t30;
t147 = pkin(8) + pkin(9);
t71 = t147 * t90;
t72 = t147 * t93;
t136 = -t89 * t71 + t92 * t72;
t110 = t93 * t126;
t117 = t90 * t130;
t44 = -t110 - t117;
t86 = t91 ^ 2;
t88 = t94 ^ 2;
t100 = (t86 - t88) * qJD(3);
t87 = t93 ^ 2;
t135 = t90 ^ 2 - t87;
t101 = t135 * qJD(4);
t124 = qJD(4) + qJD(5);
t148 = 2 * qJD(2);
t146 = t92 * pkin(4);
t57 = t89 * t93 + t92 * t90;
t33 = t124 * t57;
t145 = t33 * t91;
t139 = t94 * t33;
t138 = t94 * t95;
t133 = t86 + t88;
t132 = qJD(4) * t90;
t131 = qJD(4) * t93;
t129 = qJD(4) * t95;
t128 = qJD(5) * t89;
t127 = qJD(5) * t92;
t84 = t94 * qJD(3);
t125 = qJ(2) * qJD(3);
t123 = -0.2e1 * pkin(3) * qJD(4);
t122 = t89 * t143;
t121 = t92 * t140;
t120 = pkin(4) * t132;
t119 = pkin(4) * t128;
t118 = pkin(4) * t127;
t116 = t90 * t129;
t114 = t57 * t84;
t112 = t90 * t84;
t111 = t90 * t131;
t109 = t93 * t84;
t108 = t91 * t84;
t107 = t95 * t84;
t83 = -t93 * pkin(4) - pkin(3);
t105 = qJD(4) * t147;
t97 = pkin(3) * t94 + pkin(8) * t91;
t54 = t97 * qJD(3) + qJD(2);
t42 = t93 * t54;
t12 = t42 + (-t74 + (pkin(9) * t94 - t63) * t90) * qJD(4) + (pkin(9) * t141 + t106 * t94) * qJD(3);
t21 = -t93 * t107 + t91 * t116 - t63 * t131 - t90 * t54;
t14 = t46 * pkin(9) - t21;
t104 = t92 * t12 - t89 * t14;
t103 = t92 * t29 - t89 * t34;
t102 = -t92 * t71 - t89 * t72;
t55 = pkin(4) * t143 - t138;
t99 = t90 * t107;
t98 = t90 * t110;
t79 = t95 * t126;
t35 = -t46 * pkin(4) + t79;
t3 = -t89 * t12 - t29 * t127 + t34 * t128 - t92 * t14;
t61 = t90 * t105;
t62 = t93 * t105;
t19 = t71 * t127 + t72 * t128 + t92 * t61 + t89 * t62;
t4 = -t137 * qJD(5) + t104;
t20 = -t136 * qJD(5) + t89 * t61 - t92 * t62;
t82 = pkin(5) + t146;
t76 = 0.2e1 * t108;
t45 = t91 * t131 + t112;
t43 = t91 * t132 - t109;
t40 = t121 - t122;
t38 = t57 * t94;
t37 = t57 * t91;
t36 = t56 * pkin(5) + t83;
t32 = t124 * t144 - t93 * t127 - t92 * t131;
t31 = t38 * pkin(5) + t55;
t26 = t33 * pkin(5) + t120;
t24 = -t56 * qJ(6) + t136;
t23 = -t57 * qJ(6) + t102;
t22 = t149 * qJD(4) + t42 - t99;
t18 = -qJD(5) * t122 - t92 * t113 + t124 * t121 + t44 * t89;
t17 = t124 * t39 - t114;
t16 = t92 * t110 - t89 * t113 + t139;
t15 = -t92 * t109 + t89 * t112 + t145;
t9 = t18 * pkin(5) + t35;
t8 = -t38 * qJ(6) + t137;
t7 = t91 * pkin(5) - t40 * qJ(6) + t103;
t6 = t32 * qJ(6) - t57 * qJD(6) + t20;
t5 = -t33 * qJ(6) - t56 * qJD(6) - t19;
t2 = -t18 * qJ(6) - t38 * qJD(6) - t3;
t1 = pkin(5) * t84 + t16 * qJ(6) - t40 * qJD(6) + t4;
t10 = [0, 0, 0, 0, t148, qJ(2) * t148, -0.2e1 * t108, 0.2e1 * t100, 0, 0, 0, 0.2e1 * qJD(2) * t91 + 0.2e1 * t94 * t125, 0.2e1 * qJD(2) * t94 - 0.2e1 * t91 * t125, -0.2e1 * t87 * t108 - 0.2e1 * t88 * t111, 0.2e1 * t88 * t101 + 0.4e1 * t94 * t98, -0.2e1 * t100 * t93 - 0.2e1 * t91 * t117, 0.2e1 * t90 * t100 - 0.2e1 * t91 * t115, t76, -0.2e1 * t88 * t93 * t129 + 0.2e1 * t53 * t84 + 0.2e1 * (t22 + t99) * t91, 0.2e1 * t88 * t116 + 0.2e1 * t21 * t91 + 0.2e1 * (-t52 + t74) * t84, -0.2e1 * t40 * t16, 0.2e1 * t16 * t38 - 0.2e1 * t40 * t18, -0.2e1 * t16 * t91 + 0.2e1 * t40 * t84, -0.2e1 * t18 * t91 - 0.2e1 * t38 * t84, t76, 0.2e1 * t103 * t84 + 0.2e1 * t55 * t18 + 0.2e1 * t35 * t38 + 0.2e1 * t4 * t91, -0.2e1 * t137 * t84 - 0.2e1 * t55 * t16 + 0.2e1 * t3 * t91 + 0.2e1 * t35 * t40, -0.2e1 * t1 * t40 + 0.2e1 * t7 * t16 - 0.2e1 * t8 * t18 - 0.2e1 * t2 * t38, 0.2e1 * t7 * t1 + 0.2e1 * t8 * t2 + 0.2e1 * t31 * t9; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t133 * t131, t133 * t132, 0, 0, 0, 0, 0, t17 * t91 - t94 * t18 + (-t37 * t94 + t38 * t91) * qJD(3), t15 * t91 + t94 * t16 + (t39 * t94 + t40 * t91) * qJD(3), t15 * t38 - t37 * t16 - t17 * t40 + t39 * t18, -t1 * t37 + t126 * t31 - t8 * t15 + t7 * t17 - t2 * t39 - t9 * t94; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0.2e1 * t39 * t15 - 0.2e1 * t37 * t17 - 0.2e1 * t108; 0, 0, 0, 0, 0, 0, 0, 0, -t126, -t84, 0, -t79, -t107, -t94 * t101 - t98, -0.4e1 * t94 * t111 + t135 * t126, t45, -t43, 0 (-t90 * t138 - t93 * t97) * qJD(4) + (t90 * t96 - t74) * qJD(3) (-t93 * t138 + t90 * t97) * qJD(4) + (-pkin(8) * t140 + (pkin(3) * t93 + t142) * t91) * qJD(3), -t16 * t57 - t40 * t32, t16 * t56 - t57 * t18 + t32 * t38 - t40 * t33, -t32 * t91 + t114, -t56 * t84 - t145, 0, t102 * t84 + t120 * t38 + t83 * t18 + t20 * t91 + t55 * t33 + t35 * t56, t120 * t40 - t136 * t84 - t83 * t16 + t19 * t91 - t55 * t32 + t35 * t57, -t1 * t57 + t23 * t16 - t24 * t18 - t2 * t56 + t7 * t32 - t8 * t33 - t5 * t38 - t6 * t40, t1 * t23 + t2 * t24 + t31 * t26 + t9 * t36 + t8 * t5 + t7 * t6; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t126, -t84, 0, 0, 0, 0, 0, t44, t46, 0, 0, 0, 0, 0, t126 * t56 - t139, t126 * t57 + t94 * t32, t15 * t56 - t17 * t57 - t37 * t32 + t39 * t33, t126 * t36 - t15 * t24 + t17 * t23 - t94 * t26 - t37 * t6 - t39 * t5; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0.2e1 * t111, -0.2e1 * t101, 0, 0, 0, t90 * t123, t93 * t123, -0.2e1 * t57 * t32, 0.2e1 * t32 * t56 - 0.2e1 * t57 * t33, 0, 0, 0, 0.2e1 * t120 * t56 + 0.2e1 * t83 * t33, 0.2e1 * t120 * t57 - 0.2e1 * t83 * t32, 0.2e1 * t23 * t32 - 0.2e1 * t24 * t33 - 0.2e1 * t5 * t56 - 0.2e1 * t6 * t57, 0.2e1 * t23 * t6 + 0.2e1 * t24 * t5 + 0.2e1 * t36 * t26; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t44, t46, t84, t22, t21, 0, 0, -t16, -t18, t84, t84 * t146 + (-t30 + (-t91 * pkin(4) - t29) * t89) * qJD(5) + t104 (-t127 * t91 - t84 * t89) * pkin(4) + t3, t82 * t16 + (-t18 * t89 + (-t38 * t92 + t40 * t89) * qJD(5)) * pkin(4), t1 * t82 + (t2 * t89 + (-t7 * t89 + t8 * t92) * qJD(5)) * pkin(4); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t45, t43, 0, 0, 0, 0, 0, t17, t15, 0, t17 * t82 + (-t15 * t89 + (t37 * t89 - t39 * t92) * qJD(5)) * pkin(4); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t131, -t132, 0, -pkin(8) * t131, pkin(8) * t132, 0, 0, -t32, -t33, 0, t20, t19, t82 * t32 + (-t33 * t89 + (-t56 * t92 + t57 * t89) * qJD(5)) * pkin(4), t6 * t82 + (t5 * t89 + (-t23 * t89 + t24 * t92) * qJD(5)) * pkin(4); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -0.2e1 * t119, -0.2e1 * t118, 0, 0.2e1 * (-t82 + t146) * t119; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t16, -t18, t84, t4, t3, pkin(5) * t16, t1 * pkin(5); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t17, t15, 0, t17 * pkin(5); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t32, -t33, 0, t20, t19, pkin(5) * t32, t6 * pkin(5); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t119, -t118, 0, -pkin(5) * t119; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t9; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t126; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t26; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0;];
MMD_reg  = t10;
