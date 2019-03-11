% Calculate minimal parameter regressor of joint inertia matrix time derivative for
% S6RPRRRP10
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
% MMD_reg [((6+1)*6/2)x31]
%   minimal parameter regressor of inerta matrix time derivative
%   (only lower left triangular matrix (including diagonal) due to symmetry

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-09 06:33
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function MMD_reg = S6RPRRRP10_inertiaDJ_regmin_slag_vp(qJ, qJD, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(9,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPRRRP10_inertiaDJ_regmin_slag_vp: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RPRRRP10_inertiaDJ_regmin_slag_vp: qJD has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S6RPRRRP10_inertiaDJ_regmin_slag_vp: pkin has to be [9x1] (double)');

%% Symbolic Calculation
% From inertiaD_joint_joint_fixb_regressor_minpar_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 06:32:42
% EndTime: 2019-03-09 06:32:49
% DurationCPUTime: 2.54s
% Computational Cost: add. (2539->274), mult. (5718->442), div. (0->0), fcn. (4837->6), ass. (0->131)
t162 = sin(qJ(5));
t128 = qJD(5) * t162;
t90 = sin(qJ(4));
t167 = (t162 * qJD(4) + t128) * t90;
t92 = cos(qJ(4));
t163 = cos(qJ(5));
t126 = t163 * qJD(5);
t98 = t163 * qJD(4) + t126;
t29 = -t98 * t92 + t167;
t91 = sin(qJ(3));
t172 = t29 * t91;
t93 = cos(qJ(3));
t149 = qJD(4) * t93;
t137 = t92 * t149;
t147 = t91 * qJD(3);
t170 = -t90 * t147 + t137;
t109 = pkin(3) * t91 - pkin(8) * t93;
t57 = qJ(2) + t109;
t51 = t90 * t57;
t157 = t91 * t92;
t94 = -pkin(1) - pkin(7);
t69 = t94 * t157;
t168 = -t69 - t51;
t55 = t162 * t90 - t163 * t92;
t87 = t91 ^ 2;
t89 = t93 ^ 2;
t123 = (t87 - t89) * qJD(3);
t88 = t92 ^ 2;
t154 = t90 ^ 2 - t88;
t124 = t154 * qJD(4);
t129 = qJD(3) * t162;
t113 = t93 * t129;
t130 = qJD(3) * t163;
t115 = t93 * t130;
t56 = t162 * t92 + t163 * t90;
t30 = (qJD(4) + qJD(5)) * t56;
t161 = t30 * t91;
t14 = t90 * t113 - t92 * t115 + t161;
t116 = t91 * t129;
t117 = t91 * t130;
t160 = t30 * t93;
t15 = -t90 * t116 + t92 * t117 + t160;
t38 = t55 * t91;
t39 = t55 * t93;
t166 = (-t38 * t93 + t39 * t91) * qJD(3) - t14 * t91 - t15 * t93;
t158 = t90 * t94;
t131 = pkin(4) - t158;
t156 = t92 * t93;
t52 = t92 * t57;
t27 = -pkin(9) * t156 + t131 * t91 + t52;
t159 = t90 * t93;
t31 = -pkin(9) * t159 - t168;
t100 = t162 * t27 + t163 * t31;
t110 = pkin(3) * t93 + pkin(8) * t91;
t53 = t110 * qJD(3) + qJD(2);
t41 = t92 * t53;
t11 = t41 + (-t69 + (pkin(9) * t93 - t57) * t90) * qJD(4) + (pkin(9) * t157 + t131 * t93) * qJD(3);
t85 = t93 * qJD(3);
t134 = t94 * t85;
t148 = qJD(4) * t94;
t138 = t90 * t148;
t150 = qJD(4) * t92;
t21 = -t92 * t134 + t91 * t138 - t57 * t150 - t90 * t53;
t13 = -pkin(9) * t170 - t21;
t111 = -t163 * t11 + t162 * t13;
t4 = -t100 * qJD(5) - t111;
t165 = 0.2e1 * qJD(2);
t95 = 2 * qJD(6);
t164 = -pkin(9) - pkin(8);
t155 = t93 * t94;
t152 = t87 + t89;
t151 = qJD(4) * t90;
t146 = qJ(2) * qJD(3);
t145 = -0.2e1 * pkin(3) * qJD(4);
t144 = pkin(4) * t151;
t143 = pkin(5) * t85;
t142 = t162 * pkin(4);
t140 = t92 * t147;
t139 = t90 * t149;
t136 = t90 * t150;
t135 = t91 * t85;
t82 = -pkin(4) * t92 - pkin(3);
t54 = pkin(4) * t159 - t155;
t122 = t90 * t140;
t121 = t90 * t134;
t120 = pkin(4) * t128;
t119 = t164 * t163;
t118 = t164 * t162;
t106 = t90 * t119;
t105 = qJD(4) * t119;
t104 = qJD(4) * t118;
t67 = t164 * t92;
t18 = -qJD(5) * t106 - t92 * t104 - t90 * t105 - t67 * t128;
t33 = t90 * t118 - t163 * t67;
t103 = -t18 * t91 + t33 * t85;
t19 = -t67 * t126 - t92 * t105 + (qJD(5) * t118 + t104) * t90;
t32 = -t162 * t67 - t106;
t102 = -t19 * t91 - t32 * t85;
t101 = t56 * t147 + t29 * t93;
t75 = t94 * t147;
t34 = t170 * pkin(4) + t75;
t99 = -t162 * t31 + t163 * t27;
t3 = -t162 * t11 - t27 * t126 + t31 * t128 - t163 * t13;
t16 = t113 * t92 + t90 * t115 - t172;
t17 = -t92 * t116 - t90 * t117 + t98 * t156 - t167 * t93;
t36 = t56 * t91;
t37 = t56 * t93;
t97 = -t16 * t91 + t37 * t147 + (-qJD(3) * t36 - t17) * t93;
t79 = qJ(6) * t85;
t84 = t91 * qJD(6);
t1 = t79 + t84 - t3;
t96 = (-t91 * t142 - t100) * qJD(5) - t111;
t83 = pkin(4) * t126;
t81 = -t163 * pkin(4) - pkin(5);
t78 = t142 + qJ(6);
t77 = -0.2e1 * t120;
t72 = t83 + qJD(6);
t71 = 0.2e1 * t135;
t44 = t91 * t150 + t90 * t85;
t43 = -t139 - t140;
t42 = t91 * t151 - t92 * t85;
t26 = pkin(5) * t55 - qJ(6) * t56 + t82;
t23 = t55 * t147 - t160;
t22 = t168 * qJD(4) - t121 + t41;
t20 = pkin(5) * t37 + qJ(6) * t39 + t54;
t10 = -t91 * pkin(5) - t99;
t9 = qJ(6) * t91 + t100;
t6 = pkin(5) * t30 + qJ(6) * t29 - qJD(6) * t56 + t144;
t5 = pkin(5) * t17 + qJ(6) * t15 + qJD(6) * t39 + t34;
t2 = -t143 - t4;
t7 = [0, 0, 0, 0, t165, qJ(2) * t165, -0.2e1 * t135, 0.2e1 * t123, 0, 0, 0, 0.2e1 * qJD(2) * t91 + 0.2e1 * t93 * t146, 0.2e1 * qJD(2) * t93 - 0.2e1 * t91 * t146, -0.2e1 * t88 * t135 - 0.2e1 * t89 * t136, 0.4e1 * t93 * t122 + 0.2e1 * t89 * t124, -0.2e1 * t92 * t123 - 0.2e1 * t91 * t139, 0.2e1 * t90 * t123 - 0.2e1 * t91 * t137, t71, -0.2e1 * t89 * t92 * t148 + 0.2e1 * t52 * t85 + 0.2e1 * (t22 + t121) * t91, 0.2e1 * t89 * t138 + 0.2e1 * t21 * t91 + 0.2e1 * (-t51 + t69) * t85, 0.2e1 * t39 * t15, 0.2e1 * t15 * t37 + 0.2e1 * t17 * t39, -0.2e1 * t15 * t91 - 0.2e1 * t39 * t85, -0.2e1 * t17 * t91 - 0.2e1 * t37 * t85, t71, 0.2e1 * t54 * t17 + 0.2e1 * t34 * t37 + 0.2e1 * t4 * t91 + 0.2e1 * t99 * t85, -0.2e1 * t100 * t85 - 0.2e1 * t54 * t15 + 0.2e1 * t3 * t91 - 0.2e1 * t34 * t39, -0.2e1 * t10 * t85 + 0.2e1 * t17 * t20 - 0.2e1 * t2 * t91 + 0.2e1 * t37 * t5, -0.2e1 * t1 * t37 - 0.2e1 * t10 * t15 - 0.2e1 * t17 * t9 - 0.2e1 * t2 * t39, 0.2e1 * t1 * t91 + 0.2e1 * t15 * t20 + 0.2e1 * t39 * t5 + 0.2e1 * t9 * t85, 0.2e1 * t1 * t9 + 0.2e1 * t10 * t2 + 0.2e1 * t20 * t5; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t152 * t150, t152 * t151, 0, 0, 0, 0, 0, t97, -t166, t97, t14 * t37 - t15 * t36 - t16 * t39 + t17 * t38, t166, -t1 * t38 + t10 * t16 - t14 * t9 + t20 * t147 + t2 * t36 - t5 * t93; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0.2e1 * t14 * t38 + 0.2e1 * t16 * t36 - 0.2e1 * t135; 0, 0, 0, 0, 0, 0, 0, 0, -t147, -t85, 0, -t75, -t134, -t93 * t124 - t122, -0.4e1 * t93 * t136 + t154 * t147, t44, -t42, 0 (-t110 * t92 - t90 * t155) * qJD(4) + (t109 * t90 - t69) * qJD(3) (t110 * t90 - t92 * t155) * qJD(4) + (-pkin(8) * t156 + (pkin(3) * t92 + t158) * t91) * qJD(3), -t15 * t56 + t29 * t39, t15 * t55 - t17 * t56 + t29 * t37 + t30 * t39, t56 * t85 - t172, -t55 * t85 - t161, 0, t144 * t37 + t17 * t82 + t30 * t54 + t34 * t55 + t102, -t144 * t39 - t15 * t82 - t29 * t54 + t34 * t56 - t103, t17 * t26 + t20 * t30 + t37 * t6 + t5 * t55 + t102, -t1 * t55 - t10 * t29 - t15 * t32 - t17 * t33 + t18 * t37 - t19 * t39 + t2 * t56 - t30 * t9, t15 * t26 + t20 * t29 + t39 * t6 - t5 * t56 + t103, t1 * t33 + t10 * t19 - t18 * t9 + t2 * t32 + t20 * t6 + t26 * t5; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t147, -t85, 0, 0, 0, 0, 0, t43, -t170, 0, 0, 0, 0, 0, t23, t101, t23, t14 * t55 + t16 * t56 - t29 * t36 + t30 * t38, -t101, -t14 * t33 + t26 * t147 + t16 * t32 + t18 * t38 + t19 * t36 - t6 * t93; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0.2e1 * t136, -0.2e1 * t124, 0, 0, 0, t90 * t145, t92 * t145, -0.2e1 * t56 * t29, 0.2e1 * t29 * t55 - 0.2e1 * t30 * t56, 0, 0, 0, 0.2e1 * t144 * t55 + 0.2e1 * t30 * t82, 0.2e1 * t144 * t56 - 0.2e1 * t29 * t82, 0.2e1 * t26 * t30 + 0.2e1 * t55 * t6, 0.2e1 * t18 * t55 + 0.2e1 * t19 * t56 - 0.2e1 * t29 * t32 - 0.2e1 * t30 * t33, 0.2e1 * t26 * t29 - 0.2e1 * t56 * t6, -0.2e1 * t18 * t33 + 0.2e1 * t19 * t32 + 0.2e1 * t26 * t6; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t43, -t170, t85, t22, t21, 0, 0, -t15, -t17, t85, pkin(4) * t115 + t96 (-t91 * t126 - t113) * pkin(4) + t3 (pkin(5) - t81) * t85 + t96, -t120 * t39 - t81 * t15 - t78 * t17 - t72 * t37, t72 * t91 + t78 * t85 + t1, t1 * t78 + t10 * t120 + t2 * t81 + t9 * t72; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t44, t42, 0, 0, 0, 0, 0, -t16, t14, -t16, 0, -t14, t120 * t36 - t14 * t78 + t16 * t81 - t38 * t72; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t150, -t151, 0, -pkin(8) * t150, pkin(8) * t151, 0, 0, -t29, -t30, 0, -t19, t18, -t19, t120 * t56 - t81 * t29 - t78 * t30 - t72 * t55, -t18, t120 * t32 - t18 * t78 + t19 * t81 + t33 * t72; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t77, -0.2e1 * t83, t77, 0, 0.2e1 * t72, 0.2e1 * t120 * t81 + 0.2e1 * t78 * t72; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t15, -t17, t85, t4, t3, t4 + 0.2e1 * t143, pkin(5) * t15 - qJ(6) * t17 - qJD(6) * t37, 0.2e1 * t79 + 0.2e1 * t84 - t3, -pkin(5) * t2 + qJ(6) * t1 + qJD(6) * t9; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t16, t14, -t16, 0, -t14, -pkin(5) * t16 - qJ(6) * t14 - qJD(6) * t38; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t29, -t30, 0, -t19, t18, -t19, pkin(5) * t29 - qJ(6) * t30 - qJD(6) * t55, -t18, -pkin(5) * t19 - qJ(6) * t18 + qJD(6) * t33; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t120, -t83, -t120, 0, t95 + t83, -pkin(5) * t120 + t72 * qJ(6) + t78 * qJD(6); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t95, qJ(6) * t95; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t85, -t15, 0, t2; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t16; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t29, 0, t19; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t120; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0;];
MMD_reg  = t7;
