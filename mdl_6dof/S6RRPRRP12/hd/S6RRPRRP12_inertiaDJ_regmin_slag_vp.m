% Calculate minimal parameter regressor of joint inertia matrix time derivative for
% S6RRPRRP12
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% qJD [6x1]
%   Generalized joint velocities
% pkin [9x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d2,d4,d5]';
% 
% Output:
% MMD_reg [((6+1)*6/2)x32]
%   minimal parameter regressor of inerta matrix time derivative
%   (only lower left triangular matrix (including diagonal) due to symmetry

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-09 12:54
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function MMD_reg = S6RRPRRP12_inertiaDJ_regmin_slag_vp(qJ, qJD, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(9,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRPRRP12_inertiaDJ_regmin_slag_vp: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RRPRRP12_inertiaDJ_regmin_slag_vp: qJD has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S6RRPRRP12_inertiaDJ_regmin_slag_vp: pkin has to be [9x1] (double)');

%% Symbolic Calculation
% From inertiaD_joint_joint_fixb_regressor_minpar_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 12:53:32
% EndTime: 2019-03-09 12:53:37
% DurationCPUTime: 1.92s
% Computational Cost: add. (2773->239), mult. (5836->409), div. (0->0), fcn. (4800->6), ass. (0->137)
t168 = pkin(3) + pkin(7);
t165 = sin(qJ(5));
t166 = cos(qJ(5));
t94 = sin(qJ(4));
t96 = cos(qJ(4));
t52 = t165 * t96 + t166 * t94;
t97 = cos(qJ(2));
t37 = t52 * t97;
t95 = sin(qJ(2));
t159 = qJ(3) * t95;
t169 = pkin(2) + pkin(8);
t109 = -t169 * t97 - t159;
t49 = -pkin(1) + t109;
t62 = t168 * t95;
t53 = t94 * t62;
t162 = t96 * t49 + t53;
t129 = t165 * qJD(5);
t174 = t165 * qJD(4) + t129;
t131 = t166 * qJD(5);
t173 = t166 * qJD(4) + t131;
t93 = t97 ^ 2;
t128 = qJD(2) * (t95 ^ 2 - t93);
t90 = t94 ^ 2;
t161 = -t96 ^ 2 + t90;
t127 = t161 * qJD(4);
t172 = qJD(4) + qJD(5);
t63 = t168 * t97;
t153 = t63 * qJD(4);
t157 = qJD(3) * t97;
t171 = t109 * qJD(2) + t153 + t157;
t137 = pkin(9) * t97 - t49;
t54 = t96 * t62;
t25 = pkin(4) * t95 + t137 * t94 + t54;
t163 = t96 * t97;
t27 = -pkin(9) * t163 + t162;
t107 = t165 * t25 + t166 * t27;
t151 = t96 * qJD(2);
t141 = t95 * t151;
t154 = qJD(4) * t97;
t144 = t94 * t154;
t108 = t141 + t144;
t155 = qJD(4) * t96;
t156 = qJD(4) * t94;
t152 = t95 * qJD(2);
t126 = pkin(2) * t152 - qJD(3) * t95;
t158 = qJ(3) * t97;
t35 = (pkin(8) * t95 - t158) * qJD(2) + t126;
t86 = t97 * qJD(2);
t81 = pkin(7) * t86;
t56 = pkin(3) * t86 + t81;
t14 = -t62 * t155 + t49 * t156 - t96 * t35 - t94 * t56;
t11 = t108 * pkin(9) - t14;
t135 = -t35 * t94 + t96 * t56;
t167 = pkin(4) * t97;
t8 = (-pkin(9) * t94 * t95 + t167) * qJD(2) + (t137 * t96 - t53) * qJD(4) + t135;
t119 = t165 * t11 - t166 * t8;
t4 = -t107 * qJD(5) - t119;
t170 = 0.2e1 * qJD(3);
t98 = 2 * qJD(6);
t29 = -t173 * t94 - t174 * t96;
t138 = t165 * t94;
t139 = t166 * t96;
t51 = t138 - t139;
t164 = t51 * t29;
t77 = t94 * pkin(4) + qJ(3);
t67 = pkin(4) * t155 + qJD(3);
t150 = qJ(3) * qJD(4);
t149 = pkin(9) + t169;
t148 = -0.2e1 * pkin(1) * qJD(2);
t43 = pkin(4) * t163 + t63;
t147 = pkin(5) * t86;
t146 = pkin(7) * t152;
t145 = t165 * pkin(4);
t143 = t96 * t154;
t142 = t95 * t86;
t140 = t94 * t155;
t136 = qJD(4) * t169;
t134 = qJD(2) * t166;
t133 = qJD(2) * t165;
t125 = t149 * t96;
t124 = t94 * t141;
t123 = t97 * t139;
t122 = pkin(4) * t129;
t121 = t95 * t133;
t120 = t95 * t134;
t118 = -pkin(2) * t97 - t159;
t18 = t96 * t121 + (t174 * t97 + t120) * t94 - t172 * t123;
t117 = -t18 * t51 - t29 * t37;
t30 = t173 * t96 - t174 * t94;
t116 = -t30 * t52 + t164;
t115 = t149 * t156;
t113 = t166 * t125;
t112 = t165 * t125;
t57 = t149 * t94;
t16 = t172 * t113 - t165 * t115 - t57 * t129;
t32 = -t166 * t57 - t112;
t111 = -t16 * t95 + t32 * t86;
t17 = -t172 * t112 - t166 * t115 - t57 * t131;
t31 = -t165 * t57 + t113;
t110 = -t17 * t95 - t31 * t86;
t21 = t29 * t95 - t51 * t86;
t22 = t30 * t95 + t52 * t86;
t3 = -t166 * t11 + t27 * t129 - t25 * t131 - t165 * t8;
t106 = -t165 * t27 + t166 * t25;
t105 = pkin(5) * t29 + qJ(6) * t30 + qJD(6) * t52;
t78 = qJ(6) * t86;
t85 = t95 * qJD(6);
t1 = t78 + t85 - t3;
t10 = qJ(6) * t95 + t107;
t12 = -t95 * pkin(5) - t106;
t2 = -t147 - t4;
t104 = t1 * t52 + t10 * t30 - t12 * t29 + t2 * t51;
t55 = t168 * t152;
t103 = -t55 + (t169 * t95 - t158) * qJD(4);
t102 = t16 * t52 - t17 * t51 + t29 * t31 - t30 * t32;
t101 = t118 * qJD(2) + t157;
t33 = -pkin(4) * t144 + (-pkin(4) * t96 - t168) * t152;
t82 = pkin(4) * t131;
t66 = t82 + qJD(6);
t76 = t145 + qJ(6);
t80 = -t166 * pkin(4) - pkin(5);
t100 = -t51 * t122 + t80 * t29 - t30 * t76 - t52 * t66;
t99 = (-t95 * t145 - t107) * qJD(5) - t119;
t75 = -0.2e1 * t122;
t65 = 0.2e1 * t142;
t58 = -pkin(1) + t118;
t45 = -t95 * t155 - t94 * t86;
t44 = -t95 * t156 + t96 * t86;
t39 = -qJ(3) * t86 + t126;
t36 = t97 * t138 - t123;
t28 = pkin(5) * t52 + qJ(6) * t51 + t77;
t20 = -pkin(5) * t36 + qJ(6) * t37 + t43;
t19 = t96 * t120 - t94 * t121 + t172 * t37;
t15 = -t162 * qJD(4) + t135;
t13 = pkin(5) * t30 - qJ(6) * t29 + qJD(6) * t51 + t67;
t5 = -pkin(5) * t19 - qJ(6) * t18 + qJD(6) * t37 + t33;
t6 = [0, 0, 0, t65, -0.2e1 * t128, 0, 0, 0, t95 * t148, t97 * t148, 0, -0.2e1 * t58 * t152 + 0.2e1 * t39 * t97, -0.2e1 * t39 * t95 - 0.2e1 * t58 * t86, 0.2e1 * t58 * t39, 0.2e1 * t93 * t140 - 0.2e1 * t90 * t142, -0.4e1 * t97 * t124 - 0.2e1 * t93 * t127, 0.2e1 * t94 * t128 - 0.2e1 * t95 * t143, 0.2e1 * t96 * t128 + 0.2e1 * t95 * t144, t65, 0.2e1 * (-t63 * t151 + t15) * t95 + 0.2e1 * ((-t49 * t94 + t54) * qJD(2) - t55 * t96 - t94 * t153) * t97, 0.2e1 * (qJD(2) * t63 * t94 + t14) * t95 + 0.2e1 * (-t162 * qJD(2) - t96 * t153 + t55 * t94) * t97, -0.2e1 * t37 * t18, 0.2e1 * t18 * t36 - 0.2e1 * t19 * t37, 0.2e1 * t18 * t95 - 0.2e1 * t37 * t86, 0.2e1 * t19 * t95 + 0.2e1 * t36 * t86, t65, 0.2e1 * t106 * t86 - 0.2e1 * t43 * t19 - 0.2e1 * t33 * t36 + 0.2e1 * t4 * t95, -0.2e1 * t107 * t86 + 0.2e1 * t43 * t18 + 0.2e1 * t3 * t95 - 0.2e1 * t33 * t37, -0.2e1 * t12 * t86 - 0.2e1 * t19 * t20 - 0.2e1 * t2 * t95 - 0.2e1 * t36 * t5, 0.2e1 * t1 * t36 + 0.2e1 * t10 * t19 + 0.2e1 * t12 * t18 - 0.2e1 * t2 * t37, 0.2e1 * t1 * t95 + 0.2e1 * t10 * t86 - 0.2e1 * t18 * t20 + 0.2e1 * t37 * t5, 0.2e1 * t1 * t10 + 0.2e1 * t12 * t2 + 0.2e1 * t20 * t5; 0, 0, 0, 0, 0, t86, -t152, 0, -t81, t146, t101, t81, -t146, t101 * pkin(7), t97 * t127 + t124, 0.4e1 * t97 * t140 - t161 * t152, t44, t45, 0, t103 * t94 + t171 * t96, t103 * t96 - t171 * t94, t117, -t18 * t52 - t19 * t51 + t29 * t36 + t30 * t37, t21, -t22, 0, -t19 * t77 + t30 * t43 + t33 * t52 - t36 * t67 + t110, t18 * t77 + t29 * t43 - t33 * t51 - t37 * t67 - t111, -t13 * t36 - t19 * t28 + t20 * t30 + t5 * t52 + t110, -t16 * t36 - t17 * t37 + t18 * t31 + t19 * t32 - t104, t13 * t37 - t18 * t28 - t20 * t29 + t5 * t51 + t111, t1 * t32 - t10 * t16 + t12 * t17 + t13 * t20 + t2 * t31 + t28 * t5; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t170, qJ(3) * t170, -0.2e1 * t140, 0.2e1 * t127, 0, 0, 0, 0.2e1 * qJD(3) * t94 + 0.2e1 * t96 * t150, 0.2e1 * qJD(3) * t96 - 0.2e1 * t94 * t150, -0.2e1 * t164, -0.2e1 * t29 * t52 + 0.2e1 * t30 * t51, 0, 0, 0, 0.2e1 * t30 * t77 + 0.2e1 * t52 * t67, 0.2e1 * t29 * t77 - 0.2e1 * t51 * t67, 0.2e1 * t13 * t52 + 0.2e1 * t28 * t30, 0.2e1 * t102, 0.2e1 * t13 * t51 - 0.2e1 * t28 * t29, 0.2e1 * t13 * t28 - 0.2e1 * t16 * t32 + 0.2e1 * t17 * t31; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t86, 0, 0, t81, 0, 0, 0, 0, 0, t44, t45, 0, 0, 0, 0, 0, t21, -t22, t21, t19 * t52 + t30 * t36 - t117, t22, t104; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0.2e1 * t116, 0, -t102; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -0.2e1 * t116; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t94 * t152 - t143, t108, t86, t15, t14, 0, 0, t18, t19, t86, t134 * t167 + t99 (-t95 * t131 - t97 * t133) * pkin(4) + t3 (pkin(5) - t80) * t86 + t99, -t37 * t122 + t80 * t18 + t76 * t19 + t66 * t36, t66 * t95 + t76 * t86 + t1, t1 * t76 + t10 * t66 + t12 * t122 + t2 * t80; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t156, -t155, 0, t94 * t136, t96 * t136, 0, 0, t29, -t30, 0, -t17, t16, -t17, t100, -t16, t122 * t31 - t16 * t76 + t17 * t80 + t32 * t66; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t156, -t155, 0, 0, 0, 0, 0, t29, -t30, t29, 0, t30, -t100; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t75, -0.2e1 * t82, t75, 0, 0.2e1 * t66, 0.2e1 * t122 * t80 + 0.2e1 * t76 * t66; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t18, t19, t86, t4, t3, t4 + 0.2e1 * t147, -pkin(5) * t18 + qJ(6) * t19 + qJD(6) * t36, 0.2e1 * t78 + 0.2e1 * t85 - t3, -pkin(5) * t2 + qJ(6) * t1 + qJD(6) * t10; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t29, -t30, 0, -t17, t16, -t17, -t105, -t16, -pkin(5) * t17 - qJ(6) * t16 + qJD(6) * t32; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t29, -t30, t29, 0, t30, t105; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t122, -t82, -t122, 0, t98 + t82, -pkin(5) * t122 + t66 * qJ(6) + t76 * qJD(6); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t98, qJ(6) * t98; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t86, t18, 0, t2; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t29, 0, t17; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t29; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t122; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0;];
MMD_reg  = t6;
