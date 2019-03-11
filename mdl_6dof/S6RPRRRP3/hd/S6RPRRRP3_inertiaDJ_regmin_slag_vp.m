% Calculate minimal parameter regressor of joint inertia matrix time derivative for
% S6RPRRRP3
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% qJD [6x1]
%   Generalized joint velocities
% pkin [10x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d3,d4,d5,theta2]';
% 
% Output:
% MMD_reg [((6+1)*6/2)x29]
%   minimal parameter regressor of inerta matrix time derivative
%   (only lower left triangular matrix (including diagonal) due to symmetry

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-09 06:06
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function MMD_reg = S6RPRRRP3_inertiaDJ_regmin_slag_vp(qJ, qJD, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(10,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPRRRP3_inertiaDJ_regmin_slag_vp: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RPRRRP3_inertiaDJ_regmin_slag_vp: qJD has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S6RPRRRP3_inertiaDJ_regmin_slag_vp: pkin has to be [10x1] (double)');

%% Symbolic Calculation
% From inertiaD_joint_joint_fixb_regressor_minpar_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 06:05:16
% EndTime: 2019-03-09 06:05:21
% DurationCPUTime: 1.88s
% Computational Cost: add. (2595->250), mult. (5987->399), div. (0->0), fcn. (5111->8), ass. (0->121)
t96 = cos(qJ(3));
t136 = t96 * qJD(3);
t93 = sin(qJ(4));
t129 = t93 * t136;
t95 = cos(qJ(4));
t140 = qJD(4) * t95;
t94 = sin(qJ(3));
t164 = t94 * t140 + t129;
t163 = -0.4e1 * t94;
t156 = cos(qJ(5));
t150 = t94 * t95;
t78 = sin(pkin(10)) * pkin(1) + pkin(7);
t154 = t78 * t93;
t109 = -pkin(3) * t96 - pkin(8) * t94;
t80 = -cos(pkin(10)) * pkin(1) - pkin(2);
t55 = t109 + t80;
t44 = t95 * t55;
t29 = -pkin(9) * t150 + t44 + (-pkin(4) - t154) * t96;
t151 = t93 * t94;
t43 = t93 * t55;
t149 = t95 * t96;
t60 = t78 * t149;
t161 = -t43 - t60;
t31 = -pkin(9) * t151 - t161;
t92 = sin(qJ(5));
t162 = t156 * t31 + t92 * t29;
t115 = t156 * qJD(5);
t160 = t156 * qJD(4) + t115;
t88 = t94 ^ 2;
t113 = (-t96 ^ 2 + t88) * qJD(3);
t89 = t95 ^ 2;
t143 = t93 ^ 2 - t89;
t114 = t143 * qJD(4);
t159 = qJD(4) + qJD(5);
t139 = qJD(4) * t96;
t126 = t93 * t139;
t86 = t94 * qJD(3);
t47 = t95 * t86 + t126;
t108 = pkin(3) * t94 - pkin(8) * t96;
t61 = t108 * qJD(3);
t17 = -t55 * t140 + t47 * t78 - t93 * t61;
t121 = t78 * t86;
t144 = t93 * t121 + t95 * t61;
t157 = pkin(4) * t94;
t11 = (-pkin(9) * t149 + t157) * qJD(3) + (-t60 + (pkin(9) * t94 - t55) * t93) * qJD(4) + t144;
t15 = -pkin(9) * t164 - t17;
t4 = -qJD(5) * t162 + t156 * t11 - t92 * t15;
t97 = 2 * qJD(6);
t158 = -pkin(9) - pkin(8);
t117 = qJD(3) * t156;
t110 = t96 * t117;
t152 = t92 * t95;
t59 = t156 * t93 + t152;
t34 = t159 * t59;
t19 = -t95 * t110 + t92 * t129 + t34 * t94;
t120 = t156 * t95;
t134 = t92 * t151;
t42 = t94 * t120 - t134;
t155 = t42 * t19;
t153 = t92 * t93;
t58 = -t120 + t153;
t147 = t19 * t58 - t42 * t34;
t45 = pkin(4) * t151 + t94 * t78;
t141 = qJD(4) * t93;
t138 = qJD(5) * t92;
t137 = qJD(6) * t96;
t135 = -0.2e1 * pkin(3) * qJD(4);
t133 = 0.2e1 * qJD(3) * t80;
t62 = t78 * t136;
t35 = t164 * pkin(4) + t62;
t132 = pkin(4) * t141;
t131 = pkin(5) * t86;
t130 = pkin(4) * t138;
t128 = t95 * t136;
t127 = t94 * t141;
t124 = t95 * t139;
t123 = t93 * t140;
t122 = t94 * t136;
t84 = -pkin(4) * t95 - pkin(3);
t118 = t158 * qJD(4);
t112 = t93 * t128;
t111 = t158 * t156;
t20 = t93 * t110 - t92 * t127 - qJD(5) * t134 + (t92 * t136 + t160 * t94) * t95;
t33 = t159 * t153 - t160 * t95;
t41 = t59 * t94;
t107 = -t20 * t59 + t33 * t41;
t106 = t93 * t111;
t105 = qJD(4) * t111;
t67 = t158 * t95;
t21 = -qJD(5) * t106 - t93 * t105 - t118 * t152 - t67 * t138;
t37 = t158 * t153 - t156 * t67;
t104 = t21 * t96 + t37 * t86;
t22 = -t67 * t115 - t95 * t105 + (qJD(5) * t158 + t118) * t153;
t36 = -t92 * t67 - t106;
t103 = t22 * t96 - t36 * t86;
t25 = t33 * t96 + t59 * t86;
t24 = -t34 * t96 + t58 * t86;
t102 = t156 * t29 - t92 * t31;
t3 = -t92 * t11 - t29 * t115 + t31 * t138 - t156 * t15;
t81 = qJ(6) * t86;
t100 = -t3 + t81;
t99 = -pkin(5) * t20 - qJ(6) * t19 + qJD(6) * t42;
t98 = t96 * t130 + t4;
t85 = pkin(4) * t115;
t83 = -t156 * pkin(4) - pkin(5);
t79 = pkin(4) * t92 + qJ(6);
t77 = -0.2e1 * t130;
t73 = t85 + qJD(6);
t72 = -0.2e1 * t122;
t49 = t93 * t86 - t124;
t46 = t127 - t128;
t32 = pkin(5) * t58 - qJ(6) * t59 + t84;
t23 = pkin(5) * t41 - qJ(6) * t42 + t45;
t18 = t161 * qJD(4) + t144;
t10 = pkin(5) * t34 + qJ(6) * t33 - qJD(6) * t59 + t132;
t7 = t96 * pkin(5) - t102;
t6 = -qJ(6) * t96 + t162;
t5 = -t99 + t35;
t2 = -t131 - t4;
t1 = t100 - t137;
t8 = [0, 0, 0, 0, 0.2e1 * t122, -0.2e1 * t113, 0, 0, 0, t94 * t133, t96 * t133, 0.2e1 * t89 * t122 - 0.2e1 * t88 * t123, t112 * t163 + 0.2e1 * t88 * t114, 0.2e1 * t95 * t113 + 0.2e1 * t94 * t126, -0.2e1 * t93 * t113 + 0.2e1 * t94 * t124, t72, 0.2e1 * t44 * t86 - 0.2e1 * t18 * t96 + 0.2e1 * (t93 * t122 + t88 * t140) * t78, -0.2e1 * t88 * t78 * t141 - 0.2e1 * t17 * t96 + 0.2e1 * (-t43 + t60) * t86, -0.2e1 * t155, 0.2e1 * t41 * t19 - 0.2e1 * t20 * t42, 0.2e1 * t19 * t96 + 0.2e1 * t42 * t86, 0.2e1 * t20 * t96 - 0.2e1 * t41 * t86, t72, 0.2e1 * t102 * t86 + 0.2e1 * t45 * t20 + 0.2e1 * t35 * t41 - 0.2e1 * t4 * t96, -0.2e1 * t162 * t86 - 0.2e1 * t45 * t19 - 0.2e1 * t3 * t96 + 0.2e1 * t35 * t42, 0.2e1 * t2 * t96 + 0.2e1 * t20 * t23 + 0.2e1 * t41 * t5 - 0.2e1 * t7 * t86, -0.2e1 * t1 * t41 - 0.2e1 * t19 * t7 + 0.2e1 * t2 * t42 - 0.2e1 * t20 * t6, -0.2e1 * t1 * t96 + 0.2e1 * t19 * t23 - 0.2e1 * t42 * t5 + 0.2e1 * t6 * t86, 0.2e1 * t1 * t6 + 0.2e1 * t2 * t7 + 0.2e1 * t23 * t5; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t1 * t42 - t19 * t6 + t2 * t41 + t20 * t7 + t23 * t86 - t5 * t96; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0.2e1 * t20 * t41 - 0.2e1 * t122 - 0.2e1 * t155; 0, 0, 0, 0, 0, 0, t136, -t86, 0, -t62, t121, -t94 * t114 + t112, t123 * t163 - t143 * t136, t49, t47, 0 (pkin(8) * t149 + (-pkin(3) * t95 + t154) * t94) * qJD(4) + (t109 * t93 - t60) * qJD(3) (t108 * t93 + t78 * t150) * qJD(4) + (t109 * t95 + t96 * t154) * qJD(3), -t19 * t59 - t33 * t42, t107 + t147, t25, -t24, 0, t41 * t132 + t20 * t84 + t34 * t45 + t35 * t58 + t103, t42 * t132 - t19 * t84 - t33 * t45 + t35 * t59 - t104, t10 * t41 + t20 * t32 + t23 * t34 + t5 * t58 + t103, -t1 * t58 - t19 * t36 + t2 * t59 - t20 * t37 + t21 * t41 + t22 * t42 - t33 * t7 - t34 * t6, -t10 * t42 + t19 * t32 + t23 * t33 - t5 * t59 + t104, t1 * t37 + t10 * t23 + t2 * t36 - t21 * t6 + t22 * t7 + t32 * t5; 0, 0, 0, 0, 0, 0, 0, 0, 0, -t86, -t136, 0, 0, 0, 0, 0, -t47, t49, 0, 0, 0, 0, 0, t24, t25, t24, -t107 + t147, -t25, -t10 * t96 - t19 * t37 + t20 * t36 - t21 * t42 + t22 * t41 + t32 * t86; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0.2e1 * t123, -0.2e1 * t114, 0, 0, 0, t93 * t135, t95 * t135, -0.2e1 * t59 * t33, 0.2e1 * t33 * t58 - 0.2e1 * t34 * t59, 0, 0, 0, 0.2e1 * t58 * t132 + 0.2e1 * t34 * t84, 0.2e1 * t59 * t132 - 0.2e1 * t33 * t84, 0.2e1 * t10 * t58 + 0.2e1 * t32 * t34, 0.2e1 * t21 * t58 + 0.2e1 * t22 * t59 - 0.2e1 * t33 * t36 - 0.2e1 * t34 * t37, -0.2e1 * t10 * t59 + 0.2e1 * t32 * t33, 0.2e1 * t10 * t32 - 0.2e1 * t21 * t37 + 0.2e1 * t22 * t36; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t46, -t164, t86, t18, t17, 0, 0, -t19, -t20, t86, t117 * t157 + t98 (t115 * t96 - t92 * t86) * pkin(4) + t3 (pkin(5) - t83) * t86 + t98, t130 * t42 - t19 * t83 - t20 * t79 - t41 * t73, t79 * t86 + (-qJD(6) - t73) * t96 + t100, t1 * t79 + t130 * t7 + t2 * t83 + t6 * t73; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t164, t46, 0, 0, 0, 0, 0, -t20, t19, -t20, 0, -t19, t130 * t41 - t19 * t79 + t20 * t83 + t42 * t73; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t140, -t141, 0, -pkin(8) * t140, pkin(8) * t141, 0, 0, -t33, -t34, 0, -t22, t21, -t22, t130 * t59 - t33 * t83 - t34 * t79 - t58 * t73, -t21, t130 * t36 - t21 * t79 + t22 * t83 + t37 * t73; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t77, -0.2e1 * t85, t77, 0, 0.2e1 * t73, 0.2e1 * t130 * t83 + 0.2e1 * t73 * t79; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t19, -t20, t86, t4, t3, t4 + 0.2e1 * t131, pkin(5) * t19 - qJ(6) * t20 - qJD(6) * t41, -t3 + 0.2e1 * t81 - 0.2e1 * t137, -pkin(5) * t2 + qJ(6) * t1 + qJD(6) * t6; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t20, t19, -t20, 0, -t19, t99; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t33, -t34, 0, -t22, t21, -t22, pkin(5) * t33 - qJ(6) * t34 - qJD(6) * t58, -t21, -pkin(5) * t22 - qJ(6) * t21 + qJD(6) * t37; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t130, -t85, -t130, 0, t97 + t85, -pkin(5) * t130 + qJ(6) * t73 + qJD(6) * t79; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t97, qJ(6) * t97; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t86, -t19, 0, t2; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t20; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t33, 0, t22; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t130; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0;];
MMD_reg  = t8;
