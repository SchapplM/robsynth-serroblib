% Calculate minimal parameter regressor of joint inertia matrix time derivative for
% S6RRPPRR7
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% qJD [6x1]
%   Generalized joint velocities
% pkin [10x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,alpha2,d1,d2,d5,d6]';
% 
% Output:
% MMD_reg [((6+1)*6/2)x32]
%   minimal parameter regressor of inerta matrix time derivative
%   (only lower left triangular matrix (including diagonal) due to symmetry

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-09 09:21
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function MMD_reg = S6RRPPRR7_inertiaDJ_regmin_slag_vp(qJ, qJD, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(10,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRPPRR7_inertiaDJ_regmin_slag_vp: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RRPPRR7_inertiaDJ_regmin_slag_vp: qJD has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S6RRPPRR7_inertiaDJ_regmin_slag_vp: pkin has to be [10x1] (double)');

%% Symbolic Calculation
% From inertiaD_joint_joint_fixb_regressor_minpar_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 09:20:59
% EndTime: 2019-03-09 09:21:05
% DurationCPUTime: 1.97s
% Computational Cost: add. (1520->236), mult. (4070->442), div. (0->0), fcn. (3634->8), ass. (0->135)
t93 = sin(qJ(2));
t152 = qJD(2) * t93;
t89 = sin(pkin(6));
t130 = t89 * t152;
t153 = cos(pkin(6));
t121 = qJD(2) * t153;
t96 = cos(qJ(2));
t114 = t96 * t121;
t64 = pkin(1) * t114;
t101 = -(pkin(8) * t152 + qJD(4) * t96) * t89 + qJ(4) * t130 + t64;
t77 = t153 * qJD(3);
t23 = -t101 - t77;
t164 = t89 * t96;
t155 = qJ(3) * t93;
t44 = -pkin(2) * t164 + (-pkin(1) - t155) * t89;
t36 = pkin(3) * t164 - t44;
t22 = (pkin(4) * t93 + pkin(9) * t96) * t89 + t36;
t154 = qJ(4) * t89;
t126 = pkin(1) * t153;
t165 = t89 * t93;
t45 = -t153 * pkin(2) + pkin(8) * t165 - t96 * t126;
t30 = -t153 * pkin(3) - t93 * t154 + t45;
t26 = -t153 * pkin(9) + t30;
t92 = sin(qJ(5));
t95 = cos(qJ(5));
t171 = t92 * t22 + t95 * t26;
t74 = t93 * t126;
t170 = pkin(8) * t164 + t74;
t50 = -t153 * t92 - t95 * t164;
t94 = cos(qJ(6));
t87 = t94 ^ 2;
t91 = sin(qJ(6));
t158 = t91 ^ 2 - t87;
t123 = qJD(6) * t158;
t151 = qJD(2) * t96;
t66 = t89 * t151;
t159 = qJ(3) * t66 + qJD(3) * t165;
t97 = -pkin(2) - pkin(3);
t84 = -pkin(9) + t97;
t19 = (pkin(4) * t96 + t84 * t93) * t89 * qJD(2) + t159;
t28 = -qJD(4) * t165 + (t74 + (pkin(8) - qJ(4)) * t164) * qJD(2);
t6 = -qJD(5) * t171 + t95 * t19 - t92 * t28;
t169 = 0.2e1 * t89;
t98 = 0.2e1 * qJD(3);
t109 = t94 * t165 - t50 * t91;
t49 = -t153 * t95 + t92 * t164;
t147 = qJD(5) * t49;
t32 = t95 * t130 + t147;
t13 = t109 * qJD(6) + t32 * t94 + t91 * t66;
t168 = t13 * t91;
t167 = t84 * t92;
t166 = t84 * t95;
t163 = t92 * t94;
t162 = t94 * t95;
t90 = qJ(3) + pkin(4);
t86 = t92 ^ 2;
t88 = t95 ^ 2;
t157 = t86 - t88;
t156 = t86 + t88;
t150 = qJD(3) * t96;
t149 = qJD(5) * t109;
t34 = t91 * t165 + t50 * t94;
t148 = qJD(5) * t34;
t79 = qJD(5) * t92;
t146 = qJD(5) * t94;
t145 = qJD(5) * t95;
t144 = qJD(6) * t91;
t143 = qJD(6) * t94;
t142 = qJD(6) * t95;
t141 = -0.2e1 * pkin(5) * qJD(6);
t140 = t92 * t165;
t138 = t91 * t166;
t137 = t84 * t162;
t43 = t153 * qJ(3) + t170;
t83 = t89 ^ 2;
t136 = t83 * t151;
t135 = t91 * t147;
t134 = t49 * t146;
t133 = qJD(6) * t84 * t86;
t132 = t91 * t142;
t131 = t94 * t142;
t129 = t91 * t143;
t128 = t92 * t145;
t127 = t94 * t145;
t47 = t170 * qJD(2);
t125 = t47 * t153;
t122 = t157 * qJD(5);
t120 = t95 * t66;
t119 = t91 * t127;
t118 = t95 * pkin(5) + t92 * pkin(10);
t117 = -pkin(5) * t92 + pkin(10) * t95;
t12 = t34 * qJD(6) + t32 * t91 - t94 * t66;
t116 = -t12 + t135;
t115 = -t13 + t134;
t11 = pkin(10) * t165 + t171;
t35 = t96 * t154 - t43;
t29 = t153 * pkin(4) - t35;
t15 = -t49 * pkin(5) - t50 * pkin(10) + t29;
t8 = t94 * t11 + t91 * t15;
t112 = t95 * t22 - t92 * t26;
t110 = -t109 * t94 + t34 * t91;
t46 = pkin(8) * t130 - t64;
t56 = t118 + t90;
t39 = t91 * t56 + t137;
t10 = -pkin(5) * t165 - t112;
t4 = -pkin(5) * t66 - t6;
t108 = t10 * t143 + t4 * t91;
t107 = t10 * t144 - t4 * t94;
t31 = t50 * qJD(5) + t92 * t130;
t106 = t49 * t143 - t31 * t91;
t105 = t49 * t144 + t31 * t94;
t5 = -t22 * t145 - t92 * t19 + t26 * t79 - t95 * t28;
t52 = t92 * t144 - t127;
t53 = t92 * t146 + t132;
t103 = t106 - t149;
t102 = -t105 + t148;
t100 = pkin(10) * t66 - t5;
t99 = t31 * pkin(5) - t32 * pkin(10) - t23;
t82 = qJ(3) * t98;
t76 = 0.2e1 * t77;
t57 = 0.2e1 * t93 * t136;
t55 = -t91 * t79 + t131;
t54 = t92 * t143 + t91 * t145;
t42 = -qJD(5) * t140 + t120;
t41 = (-t93 * t145 - t92 * t151) * t89;
t40 = -t46 + t77;
t38 = t94 * t56 - t138;
t37 = pkin(2) * t130 - t159;
t27 = t97 * t130 + t159;
t18 = t94 * qJD(3) - t39 * qJD(6) + (t94 * t117 + t91 * t167) * qJD(5);
t17 = -t91 * (t117 * qJD(5) + qJD(3)) - t56 * t143 + t53 * t84;
t7 = -t91 * t11 + t94 * t15;
t2 = -t8 * qJD(6) - t91 * t100 + t94 * t99;
t1 = -t100 * t94 + t11 * t144 - t15 * t143 - t91 * t99;
t3 = [0, 0, 0, t57, 0.2e1 * (-t93 ^ 2 + t96 ^ 2) * t83 * qJD(2), t114 * t169, -0.2e1 * t121 * t165, 0, -0.2e1 * t83 * pkin(1) * t152 - 0.2e1 * t125, -0.2e1 * pkin(1) * t136 + 0.2e1 * t46 * t153, -0.2e1 * t125 + 0.2e1 * (t44 * t152 - t37 * t96) * t89 (t40 * t96 + t47 * t93 + (-t43 * t93 + t45 * t96) * qJD(2)) * t169, 0.2e1 * t40 * t153 + 0.2e1 * (-t44 * t151 - t37 * t93) * t89, 0.2e1 * t44 * t37 + 0.2e1 * t43 * t40 + 0.2e1 * t45 * t47, -0.2e1 * t23 * t153 + 0.2e1 * (t36 * t151 + t27 * t93) * t89, 0.2e1 * t28 * t153 + 0.2e1 * (t36 * t152 - t27 * t96) * t89 (t23 * t96 - t28 * t93 + (-t30 * t96 - t35 * t93) * qJD(2)) * t169, 0.2e1 * t35 * t23 + 0.2e1 * t36 * t27 + 0.2e1 * t30 * t28, 0.2e1 * t50 * t32, -0.2e1 * t50 * t31 + 0.2e1 * t32 * t49 (t50 * t151 + t32 * t93) * t169 (t49 * t151 - t31 * t93) * t169, t57, 0.2e1 * t23 * t49 + 0.2e1 * t29 * t31 + 0.2e1 * (t112 * t151 + t6 * t93) * t89, -0.2e1 * t23 * t50 + 0.2e1 * t29 * t32 + 0.2e1 * (-t151 * t171 + t5 * t93) * t89, 0.2e1 * t34 * t13, 0.2e1 * t109 * t13 - 0.2e1 * t12 * t34, -0.2e1 * t13 * t49 + 0.2e1 * t31 * t34, 0.2e1 * t109 * t31 + 0.2e1 * t12 * t49, -0.2e1 * t49 * t31, 0.2e1 * t10 * t12 - 0.2e1 * t109 * t4 - 0.2e1 * t2 * t49 + 0.2e1 * t31 * t7, -0.2e1 * t1 * t49 + 0.2e1 * t10 * t13 - 0.2e1 * t31 * t8 + 0.2e1 * t34 * t4; 0, 0, 0, 0, 0, t66, -t130, 0, -t47, t46, -t47 (t150 + (-pkin(2) * t96 - t155) * qJD(2)) * t89, -t46 + t76, -pkin(2) * t47 + qJ(3) * t40 + qJD(3) * t43, t76 + t101, t28 (-t150 + (-t96 * t97 + t155) * qJD(2)) * t89, -qJ(3) * t23 - qJD(3) * t35 + t28 * t97, -t50 * t145 - t32 * t92, t92 * t31 - t32 * t95 + (-t49 * t95 + t50 * t92) * qJD(5), t41, -t42, 0, -t66 * t167 - qJD(3) * t49 - t23 * t95 + t90 * t31 + (-t165 * t166 - t29 * t92) * qJD(5), -t84 * t120 + qJD(3) * t50 + t23 * t92 + t90 * t32 + (t84 * t140 - t29 * t95) * qJD(5), -t13 * t163 + t52 * t34, t110 * t145 + (t12 * t94 + t168 + (t109 * t91 + t34 * t94) * qJD(6)) * t92 (t13 + t134) * t95 + (-t105 - t148) * t92 (-t12 - t135) * t95 + (-t106 - t149) * t92, t31 * t95 + t49 * t79, -t18 * t49 + t38 * t31 + (t2 + (-t10 * t91 - t109 * t84) * qJD(5)) * t95 + (-qJD(5) * t7 + t12 * t84 - t108) * t92, -t17 * t49 - t39 * t31 + (t1 + (-t10 * t94 + t34 * t84) * qJD(5)) * t95 + (qJD(5) * t8 + t13 * t84 + t107) * t92; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t98, t82, t98, 0, 0, t82, 0.2e1 * t128, -0.2e1 * t122, 0, 0, 0, 0.2e1 * qJD(3) * t95 - 0.2e1 * t90 * t79, -0.2e1 * qJD(3) * t92 - 0.2e1 * t90 * t145, 0.2e1 * t87 * t128 - 0.2e1 * t86 * t129, -0.4e1 * t92 * t119 + 0.2e1 * t86 * t123, 0.2e1 * t92 * t132 + 0.2e1 * t157 * t146, -0.2e1 * t122 * t91 + 0.2e1 * t92 * t131, -0.2e1 * t128, -0.2e1 * t94 * t133 + 0.2e1 * t18 * t95 + 0.2e1 * (-t38 - 0.2e1 * t138) * t79, 0.2e1 * t91 * t133 + 0.2e1 * t17 * t95 + 0.2e1 * (t39 - 0.2e1 * t137) * t79; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t66, 0, t47, 0, 0, -t66, t28, 0, 0, 0, 0, 0, t41, -t42, 0, 0, 0, 0, 0, t103 * t95 - t116 * t92, t102 * t95 - t115 * t92; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t156 * t143, t156 * t144; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t66, t130, 0, t27, 0, 0, 0, 0, 0, t42, t41, 0, 0, 0, 0, 0, t103 * t92 + t116 * t95, t102 * t92 + t115 * t95; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t32, -t31, t66, t6, t5, t34 * t143 + t168, -t110 * qJD(6) - t91 * t12 + t13 * t94, -t106, t105, 0, -pkin(5) * t12 + pkin(10) * t106 + t107, -pkin(5) * t13 - pkin(10) * t105 + t108; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t145, t79, 0, -t84 * t145, t84 * t79, t92 * t123 - t119, 0.4e1 * t92 * t129 + t158 * t145, t55, -t53, 0 (-pkin(10) * t162 + (pkin(5) * t94 + t84 * t91) * t92) * qJD(6) + (t118 * t91 - t137) * qJD(5) (t117 * t91 + t84 * t163) * qJD(6) + (t118 * t94 + t138) * qJD(5); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t145, t79, 0, 0, 0, 0, 0, t52, t54; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t79, -t145, 0, 0, 0, 0, 0, -t53, -t55; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0.2e1 * t129, -0.2e1 * t123, 0, 0, 0, t91 * t141, t94 * t141; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t13, -t12, t31, t2, t1; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t52, t54, -t79, t18, t17; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t55, t53; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t54, t52; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t143, -t144, 0, -pkin(10) * t143, pkin(10) * t144; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0;];
MMD_reg  = t3;
