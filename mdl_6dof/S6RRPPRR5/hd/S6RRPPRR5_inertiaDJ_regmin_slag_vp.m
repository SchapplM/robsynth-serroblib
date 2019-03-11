% Calculate minimal parameter regressor of joint inertia matrix time derivative for
% S6RRPPRR5
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
% Datum: 2019-03-09 09:12
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function MMD_reg = S6RRPPRR5_inertiaDJ_regmin_slag_vp(qJ, qJD, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(10,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRPPRR5_inertiaDJ_regmin_slag_vp: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RRPPRR5_inertiaDJ_regmin_slag_vp: qJD has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S6RRPPRR5_inertiaDJ_regmin_slag_vp: pkin has to be [10x1] (double)');

%% Symbolic Calculation
% From inertiaD_joint_joint_fixb_regressor_minpar_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 09:10:59
% EndTime: 2019-03-09 09:11:04
% DurationCPUTime: 1.67s
% Computational Cost: add. (1466->237), mult. (4002->439), div. (0->0), fcn. (3562->8), ass. (0->132)
t100 = cos(qJ(2));
t92 = sin(pkin(6));
t153 = t100 * t92;
t97 = sin(qJ(2));
t154 = qJ(3) * t97;
t45 = -pkin(2) * t153 + (-pkin(1) - t154) * t92;
t37 = pkin(3) * t153 - t45;
t24 = (pkin(4) * t100 - pkin(9) * t97) * t92 + t37;
t93 = cos(pkin(6));
t166 = pkin(1) * t93;
t106 = pkin(8) * t153 + t97 * t166;
t44 = t93 * qJ(3) + t106;
t36 = -qJ(4) * t153 + t44;
t30 = -pkin(9) * t93 + t36;
t96 = sin(qJ(5));
t99 = cos(qJ(5));
t168 = t96 * t24 + t99 * t30;
t162 = t92 * t97;
t48 = t106 * qJD(2);
t141 = qJD(2) * t100;
t65 = t92 * t141;
t29 = -qJ(4) * t65 - qJD(4) * t162 + t48;
t98 = cos(qJ(6));
t90 = t98 ^ 2;
t95 = sin(qJ(6));
t156 = t95 ^ 2 - t90;
t122 = qJD(6) * t156;
t157 = qJ(3) * t65 + qJD(3) * t162;
t101 = -pkin(2) - pkin(3);
t87 = pkin(4) - t101;
t20 = (-pkin(9) * t100 - t87 * t97) * t92 * qJD(2) + t157;
t152 = qJD(2) * t97;
t126 = t92 * t152;
t75 = t100 * t166;
t62 = qJD(2) * t75;
t105 = qJ(4) * t126 + t62 + (-pkin(8) * t152 - qJD(4) * t100) * t92;
t77 = t93 * qJD(3);
t25 = t77 + t105;
t6 = -qJD(5) * t168 + t99 * t20 - t96 * t25;
t167 = 0.2e1 * t92;
t102 = 0.2e1 * qJD(3);
t51 = t99 * t162 - t93 * t96;
t110 = t98 * t153 - t51 * t95;
t50 = t96 * t162 + t93 * t99;
t33 = -qJD(5) * t50 + t99 * t65;
t14 = t110 * qJD(6) - t95 * t126 + t33 * t98;
t165 = t14 * t95;
t164 = t14 * t98;
t163 = t48 * t93;
t94 = qJ(3) - pkin(9);
t161 = t94 * t95;
t160 = t94 * t99;
t89 = t96 ^ 2;
t155 = -t99 ^ 2 + t89;
t151 = qJD(3) * t98;
t150 = qJD(3) * t99;
t149 = qJD(5) * t110;
t35 = t95 * t153 + t51 * t98;
t148 = qJD(5) * t35;
t147 = qJD(5) * t95;
t80 = qJD(5) * t96;
t146 = qJD(5) * t98;
t145 = qJD(5) * t99;
t144 = qJD(6) * t95;
t143 = qJD(6) * t96;
t79 = qJD(6) * t98;
t142 = qJD(6) * t99;
t140 = qJD(3) * t100;
t139 = qJD(5) * t100;
t138 = -0.2e1 * pkin(5) * qJD(6);
t137 = t95 * t160;
t136 = t98 * t160;
t135 = -0.2e1 * qJD(5) * t87;
t46 = -t93 * pkin(2) + pkin(8) * t162 - t75;
t134 = pkin(10) * t79;
t133 = t96 * t152;
t132 = t99 * t152;
t131 = t50 * t147;
t130 = t50 * t146;
t129 = t98 * t145;
t128 = t95 * t142;
t127 = t98 * t142;
t125 = t95 * t79;
t124 = t96 * t145;
t86 = t92 ^ 2;
t123 = t86 * t141;
t121 = t155 * qJD(5);
t120 = t95 * t129;
t119 = t97 * t123;
t31 = -t93 * pkin(3) - qJ(4) * t162 + t46;
t118 = t99 * pkin(5) + t96 * pkin(10);
t117 = -pkin(5) * t96 + pkin(10) * t99;
t11 = pkin(10) * t153 + t168;
t26 = t93 * pkin(4) - t31;
t15 = pkin(5) * t50 - pkin(10) * t51 + t26;
t8 = t11 * t98 + t15 * t95;
t115 = t24 * t99 - t30 * t96;
t113 = -t110 * t98 + t35 * t95;
t47 = pkin(8) * t126 - t62;
t57 = t118 + t87;
t40 = t57 * t95 + t136;
t10 = -pkin(5) * t153 - t115;
t4 = pkin(5) * t126 - t6;
t112 = t10 * t79 + t4 * t95;
t111 = t10 * t144 - t4 * t98;
t32 = -t93 * t80 + (t96 * t141 + t97 * t145) * t92;
t109 = t32 * t95 + t50 * t79;
t19 = t50 * t144 - t32 * t98;
t108 = -t96 * qJD(3) - t94 * t145;
t107 = t94 * t80 - t150;
t5 = -t24 * t145 - t96 * t20 - t99 * t25 + t30 * t80;
t52 = t95 * t143 - t129;
t53 = t96 * t146 + t128;
t104 = pkin(10) * t126 + t5;
t103 = t32 * pkin(5) - t33 * pkin(10) - t29;
t85 = qJ(3) * t102;
t76 = 0.2e1 * t77;
t55 = -t95 * t80 + t127;
t54 = t95 * t145 + t96 * t79;
t43 = (t96 * t139 + t132) * t92;
t42 = (-t99 * t139 + t133) * t92;
t41 = -t47 + t77;
t39 = t57 * t98 - t137;
t38 = pkin(2) * t126 - t157;
t28 = t101 * t126 + t157;
t17 = -t95 * t150 - t40 * qJD(6) + (t98 * t117 + t96 * t161) * qJD(5);
t16 = -t117 * t147 - t98 * t150 + t53 * t94 - t57 * t79;
t13 = t35 * qJD(6) + t98 * t126 + t33 * t95;
t7 = -t11 * t95 + t15 * t98;
t2 = -t8 * qJD(6) + t98 * t103 + t95 * t104;
t1 = -t95 * t103 + t98 * t104 + t11 * t144 - t15 * t79;
t3 = [0, 0, 0, 0.2e1 * t119, 0.2e1 * (t100 ^ 2 - t97 ^ 2) * t86 * qJD(2), 0.2e1 * t93 * t65, -0.2e1 * t93 * t126, 0, -0.2e1 * pkin(1) * t86 * t152 - 0.2e1 * t163, -0.2e1 * pkin(1) * t123 + 0.2e1 * t47 * t93, -0.2e1 * t163 + 0.2e1 * (-t100 * t38 + t45 * t152) * t92 (t100 * t41 + t48 * t97 + (t100 * t46 - t44 * t97) * qJD(2)) * t167, 0.2e1 * t41 * t93 + 0.2e1 * (-t45 * t141 - t38 * t97) * t92, 0.2e1 * t38 * t45 + 0.2e1 * t41 * t44 + 0.2e1 * t46 * t48, -0.2e1 * t29 * t93 + 0.2e1 * (t100 * t28 - t37 * t152) * t92, 0.2e1 * t25 * t93 + 0.2e1 * (t37 * t141 + t28 * t97) * t92 (-t100 * t25 - t29 * t97 + (-t100 * t31 + t36 * t97) * qJD(2)) * t167, 0.2e1 * t25 * t36 + 0.2e1 * t28 * t37 + 0.2e1 * t29 * t31, 0.2e1 * t51 * t33, -0.2e1 * t32 * t51 - 0.2e1 * t33 * t50 (t100 * t33 - t51 * t152) * t167 (-t100 * t32 + t50 * t152) * t167, -0.2e1 * t119, 0.2e1 * t26 * t32 - 0.2e1 * t29 * t50 + 0.2e1 * (t6 * t100 - t115 * t152) * t92, 0.2e1 * t26 * t33 - 0.2e1 * t29 * t51 + 0.2e1 * (t5 * t100 + t152 * t168) * t92, 0.2e1 * t35 * t14, 0.2e1 * t110 * t14 - 0.2e1 * t13 * t35, 0.2e1 * t14 * t50 + 0.2e1 * t32 * t35, 0.2e1 * t110 * t32 - 0.2e1 * t13 * t50, 0.2e1 * t50 * t32, 0.2e1 * t10 * t13 - 0.2e1 * t110 * t4 + 0.2e1 * t2 * t50 + 0.2e1 * t32 * t7, 0.2e1 * t1 * t50 + 0.2e1 * t10 * t14 - 0.2e1 * t32 * t8 + 0.2e1 * t35 * t4; 0, 0, 0, 0, 0, t65, -t126, 0, -t48, t47, -t48 (t140 + (-pkin(2) * t100 - t154) * qJD(2)) * t92, -t47 + t76, -pkin(2) * t48 + qJ(3) * t41 + qJD(3) * t44, -t29, t76 + t105 (-t140 + (-t100 * t101 + t154) * qJD(2)) * t92, qJ(3) * t25 + qJD(3) * t36 + t101 * t29, -t51 * t145 - t33 * t96, t96 * t32 - t33 * t99 + (t50 * t99 + t51 * t96) * qJD(5), t42, t43, 0, -t26 * t80 - t29 * t99 + t87 * t32 + (t108 * t100 + t94 * t133) * t92, -t26 * t145 + t29 * t96 + t87 * t33 + (t107 * t100 + t94 * t132) * t92, -t96 * t164 + t35 * t52, t113 * t145 + (t13 * t98 + t165 + (t110 * t95 + t35 * t98) * qJD(6)) * t96 (t14 - t130) * t99 + (t19 - t148) * t96 (-t13 + t131) * t99 + (t109 - t149) * t96, t32 * t99 - t50 * t80, t17 * t50 + t39 * t32 + (t2 + (-t10 * t95 - t110 * t94) * qJD(5)) * t99 + (-qJD(3) * t110 - qJD(5) * t7 + t13 * t94 - t112) * t96, t16 * t50 - t40 * t32 + (t1 + (-t10 * t98 + t35 * t94) * qJD(5)) * t99 + (qJD(3) * t35 + qJD(5) * t8 + t14 * t94 + t111) * t96; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t102, t85, 0, t102, 0, t85, 0.2e1 * t124, -0.2e1 * t121, 0, 0, 0, t96 * t135, t99 * t135, 0.2e1 * t124 * t90 - 0.2e1 * t125 * t89, -0.4e1 * t120 * t96 + 0.2e1 * t122 * t89, 0.2e1 * t128 * t96 + 0.2e1 * t146 * t155, -0.2e1 * t121 * t95 + 0.2e1 * t127 * t96, -0.2e1 * t124, 0.2e1 * t17 * t99 + 0.2e1 * (-qJD(3) * t95 - t94 * t79) * t89 + 0.2e1 * (-t39 - 0.2e1 * t137) * t80, 0.2e1 * t16 * t99 + 0.2e1 * (t94 * t144 - t151) * t89 + 0.2e1 * (t40 - 0.2e1 * t136) * t80; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t65, 0, t48, 0, 0, -t65, t29, 0, 0, 0, 0, 0, -t32, -t33, 0, 0, 0, 0, 0, t19, t109; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t80, t145, 0, 0, 0, 0, 0, t53, t55; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t126, t65, 0, t28, 0, 0, 0, 0, 0, -t43, t42, 0, 0, 0, 0, 0 (-t13 - t131) * t99 + (-t109 - t149) * t96 (-t14 - t130) * t99 + (t19 + t148) * t96; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t33, -t32, -t126, t6, t5, t35 * t79 + t165, -qJD(6) * t113 - t95 * t13 + t164, t109, -t19, 0, -pkin(5) * t13 - pkin(10) * t109 + t111, -pkin(5) * t14 + pkin(10) * t19 + t112; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t145, t80, 0, t108, t107, t122 * t96 - t120, 0.4e1 * t125 * t96 + t145 * t156, t55, -t53, 0 (-t134 + (pkin(5) * t95 - t94 * t98) * qJD(5)) * t99 + (pkin(10) * t147 - t151 + (pkin(5) * t98 + t161) * qJD(6)) * t96 (qJD(5) * t118 + t94 * t143) * t98 + (qJD(6) * t117 - t108) * t95; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t80, -t145, 0, 0, 0, 0, 0, -t53, -t55; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0.2e1 * t125, -0.2e1 * t122, 0, 0, 0, t95 * t138, t98 * t138; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t14, -t13, t32, t2, t1; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t52, t54, -t80, t17, t16; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t144, t79; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t54, t52; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t79, -t144, 0, -t134, pkin(10) * t144; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0;];
MMD_reg  = t3;
