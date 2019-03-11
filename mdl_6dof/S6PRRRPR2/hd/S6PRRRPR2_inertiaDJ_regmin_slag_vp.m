% Calculate minimal parameter regressor of joint inertia matrix time derivative for
% S6PRRRPR2
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% qJD [6x1]
%   Generalized joint velocities
% pkin [12x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,alpha2,d2,d3,d4,d6,theta1,theta5]';
% 
% Output:
% MMD_reg [((6+1)*6/2)x29]
%   minimal parameter regressor of inerta matrix time derivative
%   (only lower left triangular matrix (including diagonal) due to symmetry

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-08 23:09
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function MMD_reg = S6PRRRPR2_inertiaDJ_regmin_slag_vp(qJ, qJD, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(12,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6PRRRPR2_inertiaDJ_regmin_slag_vp: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6PRRRPR2_inertiaDJ_regmin_slag_vp: qJD has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [12 1]), ...
  'S6PRRRPR2_inertiaDJ_regmin_slag_vp: pkin has to be [12x1] (double)');

%% Symbolic Calculation
% From inertiaD_joint_joint_fixb_regressor_minpar_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-08 23:08:58
% EndTime: 2019-03-08 23:09:04
% DurationCPUTime: 1.88s
% Computational Cost: add. (2374->216), mult. (6071->389), div. (0->0), fcn. (5958->12), ass. (0->137)
t105 = sin(pkin(12));
t107 = cos(pkin(12));
t148 = t105 ^ 2 + t107 ^ 2;
t171 = t148 * qJD(5);
t172 = 0.2e1 * t171;
t114 = cos(qJ(4));
t145 = qJD(4) * t114;
t139 = pkin(3) * t145;
t93 = qJD(5) + t139;
t170 = t148 * t93;
t169 = qJD(3) + qJD(4);
t168 = -pkin(9) - pkin(8);
t167 = pkin(3) * t114;
t102 = t107 * pkin(10);
t113 = cos(qJ(6));
t109 = sin(qJ(6));
t153 = t105 * t109;
t81 = -t113 * t107 + t153;
t71 = t81 * qJD(6);
t96 = -pkin(5) * t107 - pkin(4);
t166 = t96 * t71;
t152 = t105 * t113;
t82 = t107 * t109 + t152;
t72 = t82 * qJD(6);
t165 = t96 * t72;
t110 = sin(qJ(4));
t111 = sin(qJ(3));
t115 = cos(qJ(3));
t83 = t110 * t111 - t114 * t115;
t59 = t169 * t83;
t159 = t105 * t59;
t138 = qJD(3) * t168;
t130 = t115 * t138;
t131 = t111 * t138;
t91 = t168 * t111;
t92 = t168 * t115;
t66 = t110 * t91 - t114 * t92;
t42 = qJD(4) * t66 + t110 * t131 - t114 * t130;
t25 = -pkin(5) * t159 + t42;
t84 = t110 * t115 + t111 * t114;
t158 = t105 * t84;
t65 = -t110 * t92 - t114 * t91;
t47 = pkin(5) * t158 + t65;
t164 = t25 * t81 + t47 * t72;
t163 = t25 * t82 - t47 * t71;
t144 = t111 * qJD(3);
t141 = pkin(3) * t144;
t60 = t169 * t84;
t29 = pkin(4) * t60 + qJ(5) * t59 - qJD(5) * t84 + t141;
t146 = qJD(4) * t110;
t41 = -t110 * t130 - t114 * t131 - t145 * t91 - t146 * t92;
t14 = t105 * t29 - t107 * t41;
t99 = -pkin(3) * t115 - pkin(2);
t54 = pkin(4) * t83 - qJ(5) * t84 + t99;
t35 = t105 * t54 + t107 * t66;
t140 = pkin(3) * t146;
t88 = t96 - t167;
t162 = t81 * t140 + t88 * t72;
t161 = t82 * t140 - t88 * t71;
t157 = t107 * t42;
t156 = t107 * t59;
t106 = sin(pkin(6));
t116 = cos(qJ(2));
t150 = t106 * t116;
t136 = qJD(2) * t150;
t108 = cos(pkin(6));
t112 = sin(qJ(2));
t151 = t106 * t112;
t74 = t108 * t111 + t115 * t151;
t117 = qJD(3) * t74 + t111 * t136;
t73 = t108 * t115 - t111 * t151;
t49 = t110 * t73 + t114 * t74;
t62 = qJD(3) * t73 + t115 * t136;
t23 = qJD(4) * t49 + t110 * t62 + t114 * t117;
t155 = t23 * t107;
t154 = qJD(6) * t84;
t147 = qJD(2) * t112;
t143 = t115 * qJD(3);
t142 = -0.2e1 * pkin(2) * qJD(3);
t137 = t106 * t147;
t13 = t105 * t41 + t107 * t29;
t34 = -t105 * t66 + t107 * t54;
t133 = t105 * t140;
t132 = t107 * t140;
t48 = t110 * t74 - t114 * t73;
t129 = t23 * t84 - t48 * t59;
t128 = t42 * t84 - t59 * t65;
t6 = -t105 * t13 + t107 * t14;
t22 = t110 * t117 - t114 * t62 - t145 * t73 + t146 * t74;
t15 = t105 * t22 + t107 * t137;
t16 = t105 * t137 - t107 * t22;
t7 = -t105 * t15 + t107 * t16;
t24 = pkin(5) * t83 - t84 * t102 + t34;
t28 = -pkin(10) * t158 + t35;
t127 = t109 * t28 - t113 * t24;
t126 = t109 * t24 + t113 * t28;
t43 = -t105 * t49 - t107 * t150;
t44 = -t105 * t150 + t107 * t49;
t125 = t109 * t44 - t113 * t43;
t124 = t109 * t43 + t113 * t44;
t95 = pkin(3) * t110 + qJ(5);
t78 = (-pkin(10) - t95) * t105;
t79 = t107 * t95 + t102;
t123 = t109 * t79 - t113 * t78;
t122 = t109 * t78 + t113 * t79;
t89 = (-pkin(10) - qJ(5)) * t105;
t90 = qJ(5) * t107 + t102;
t121 = t109 * t90 - t113 * t89;
t120 = t109 * t89 + t113 * t90;
t119 = pkin(4) * t59 - qJ(5) * t60 - qJD(5) * t83;
t98 = -pkin(4) - t167;
t118 = t140 * t84 - t59 * t98 - t60 * t95 - t83 * t93;
t58 = -0.2e1 * t82 * t71;
t51 = t81 * t84;
t50 = t82 * t84;
t46 = -qJD(5) * t82 - qJD(6) * t120;
t45 = qJD(5) * t81 + qJD(6) * t121;
t38 = t42 * t105;
t36 = 0.2e1 * t71 * t81 - 0.2e1 * t72 * t82;
t33 = -qJD(6) * t122 - t82 * t93;
t32 = qJD(6) * t123 + t81 * t93;
t31 = -t60 * t81 - t72 * t83;
t30 = t60 * t82 - t71 * t83;
t19 = t23 * t105;
t18 = -t59 * t152 - t153 * t154 + (-t109 * t59 + t113 * t154) * t107;
t17 = -t154 * t82 + t59 * t81;
t12 = t23 * t82 - t48 * t71;
t11 = t23 * t81 + t48 * t72;
t10 = t17 * t82 + t51 * t71;
t9 = pkin(10) * t159 + t14;
t8 = pkin(5) * t60 + pkin(10) * t156 + t13;
t5 = -qJD(6) * t124 - t109 * t16 + t113 * t15;
t4 = qJD(6) * t125 - t109 * t15 - t113 * t16;
t3 = -t17 * t81 - t18 * t82 + t50 * t71 + t51 * t72;
t2 = -qJD(6) * t126 - t109 * t9 + t113 * t8;
t1 = qJD(6) * t127 - t109 * t8 - t113 * t9;
t20 = [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0.2e1 * t15 * t43 + 0.2e1 * t16 * t44 + 0.2e1 * t23 * t48, 0, 0, 0, 0, 0, 0, 0; 0, 0, -t137, -t136, 0, 0, 0, 0, 0 (-t115 * t147 - t116 * t144) * t106 (t111 * t147 - t116 * t143) * t106, 0, 0, 0, 0, 0 (-t116 * t60 + t147 * t83) * t106 (t116 * t59 + t147 * t84) * t106, t105 * t129 + t15 * t83 + t43 * t60, t129 * t107 - t16 * t83 - t44 * t60 (-t15 * t84 + t43 * t59) * t107 + (-t16 * t84 + t44 * t59) * t105, t13 * t43 + t14 * t44 + t15 * t34 + t16 * t35 + t23 * t65 + t42 * t48, 0, 0, 0, 0, 0, -t125 * t60 + t48 * t18 + t23 * t50 + t5 * t83, -t124 * t60 + t48 * t17 - t23 * t51 + t4 * t83; 0, 0, 0, 0, 0.2e1 * t111 * t143, 0.2e1 * (-t111 ^ 2 + t115 ^ 2) * qJD(3), 0, 0, 0, t111 * t142, t115 * t142, -0.2e1 * t84 * t59, 0.2e1 * t59 * t83 - 0.2e1 * t60 * t84, 0, 0, 0, 0.2e1 * t141 * t83 + 0.2e1 * t60 * t99, 0.2e1 * t141 * t84 - 0.2e1 * t59 * t99, 0.2e1 * t105 * t128 + 0.2e1 * t13 * t83 + 0.2e1 * t34 * t60, 0.2e1 * t128 * t107 - 0.2e1 * t14 * t83 - 0.2e1 * t35 * t60, 0.2e1 * (-t13 * t84 + t34 * t59) * t107 + 0.2e1 * (-t14 * t84 + t35 * t59) * t105, 0.2e1 * t13 * t34 + 0.2e1 * t14 * t35 + 0.2e1 * t42 * t65, -0.2e1 * t51 * t17, -0.2e1 * t17 * t50 + 0.2e1 * t18 * t51, 0.2e1 * t17 * t83 - 0.2e1 * t51 * t60, -0.2e1 * t18 * t83 - 0.2e1 * t50 * t60, 0.2e1 * t83 * t60, -0.2e1 * t127 * t60 + 0.2e1 * t47 * t18 + 0.2e1 * t2 * t83 + 0.2e1 * t25 * t50, 0.2e1 * t1 * t83 - 0.2e1 * t126 * t60 + 0.2e1 * t47 * t17 - 0.2e1 * t25 * t51; 0, 0, 0, 0, 0, 0, 0, 0, 0, -t117, -t62, 0, 0, 0, 0, 0, -t23, t22, -t155, t19, t7, t48 * t140 + t23 * t98 + (t16 * t95 + t44 * t93) * t107 + (-t15 * t95 - t43 * t93) * t105, 0, 0, 0, 0, 0, t11, t12; 0, 0, 0, 0, 0, 0, t143, -t144, 0, -pkin(8) * t143, pkin(8) * t144, 0, 0, -t59, -t60, 0, -t42, t41, t105 * t118 - t157, t118 * t107 + t38, t6, t65 * t140 + t42 * t98 + (t14 * t95 + t35 * t93) * t107 + (-t13 * t95 - t34 * t93) * t105, t10, t3, t30, t31, 0, -t123 * t60 + t140 * t50 + t88 * t18 + t33 * t83 + t164, -t122 * t60 - t140 * t51 + t88 * t17 + t32 * t83 + t163; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -0.2e1 * t140, -0.2e1 * t139, -0.2e1 * t132, 0.2e1 * t133, 0.2e1 * t170, 0.2e1 * t140 * t98 + 0.2e1 * t170 * t95, t58, t36, 0, 0, 0, 0.2e1 * t162, 0.2e1 * t161; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t23, t22, -t155, t19, t7, -pkin(4) * t23 + (-t105 * t43 + t107 * t44) * qJD(5) + t7 * qJ(5), 0, 0, 0, 0, 0, t11, t12; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t59, -t60, 0, -t42, t41, t105 * t119 - t157, t119 * t107 + t38, t6, -pkin(4) * t42 + (-t105 * t34 + t107 * t35) * qJD(5) + t6 * qJ(5), t10, t3, t30, t31, 0, -t121 * t60 + t96 * t18 + t46 * t83 + t164, -t120 * t60 + t96 * t17 + t45 * t83 + t163; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t140, -t139, -t132, t133, t171 + t170, -pkin(4) * t140 + qJ(5) * t170 + t171 * t95, t58, t36, 0, 0, 0, t162 + t165, t161 - t166; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t172, qJ(5) * t172, t58, t36, 0, 0, 0, 0.2e1 * t165, -0.2e1 * t166; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t23, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t159, -t156, 0, t42, 0, 0, 0, 0, 0, t18, t17; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t140, 0, 0, 0, 0, 0, t72, -t71; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t72, -t71; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t5, t4; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t17, -t18, t60, t2, t1; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t71, -t72, 0, t33, t32; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t71, -t72, 0, t46, t45; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0;];
MMD_reg  = t20;
