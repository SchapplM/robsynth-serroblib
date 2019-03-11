% Calculate minimal parameter regressor of joint inertia matrix time derivative for
% S6PRRPRR2
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% qJD [6x1]
%   Generalized joint velocities
% pkin [12x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,alpha2,d2,d3,d5,d6,theta1,theta4]';
% 
% Output:
% MMD_reg [((6+1)*6/2)x27]
%   minimal parameter regressor of inerta matrix time derivative
%   (only lower left triangular matrix (including diagonal) due to symmetry

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-08 22:01
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function MMD_reg = S6PRRPRR2_inertiaDJ_regmin_slag_vp(qJ, qJD, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(12,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6PRRPRR2_inertiaDJ_regmin_slag_vp: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6PRRPRR2_inertiaDJ_regmin_slag_vp: qJD has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [12 1]), ...
  'S6PRRPRR2_inertiaDJ_regmin_slag_vp: pkin has to be [12x1] (double)');

%% Symbolic Calculation
% From inertiaD_joint_joint_fixb_regressor_minpar_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-08 22:00:16
% EndTime: 2019-03-08 22:00:20
% DurationCPUTime: 1.29s
% Computational Cost: add. (1817->198), mult. (4617->376), div. (0->0), fcn. (4645->12), ass. (0->134)
t138 = cos(pkin(12));
t86 = sin(qJ(3));
t112 = t138 * t86;
t81 = sin(pkin(12));
t90 = cos(qJ(3));
t63 = t81 * t90 + t112;
t84 = sin(qJ(6));
t85 = sin(qJ(5));
t88 = cos(qJ(6));
t89 = cos(qJ(5));
t66 = t84 * t89 + t88 * t85;
t35 = t66 * t63;
t111 = t138 * t90;
t150 = t81 * t86;
t62 = -t111 + t150;
t77 = -t90 * pkin(3) - pkin(2);
t40 = t62 * pkin(4) - t63 * pkin(9) + t77;
t141 = -qJ(4) - pkin(8);
t69 = t141 * t90;
t48 = -t138 * t69 + t141 * t150;
t43 = t89 * t48;
t140 = t85 * t40 + t43;
t135 = qJD(5) * t85;
t131 = t86 * qJD(3);
t57 = qJD(3) * t111 - t81 * t131;
t142 = t89 * t57;
t94 = -t63 * t135 + t142;
t129 = qJD(5) + qJD(6);
t157 = 0.2e1 * qJD(5);
t56 = t63 * qJD(3);
t156 = t56 * pkin(5);
t155 = t62 * pkin(5);
t134 = qJD(5) * t89;
t110 = qJD(3) * t141;
t55 = t90 * qJD(4) + t86 * t110;
t93 = -t86 * qJD(4) + t90 * t110;
t30 = t138 * t55 + t81 * t93;
t78 = pkin(3) * t131;
t31 = t56 * pkin(4) - t57 * pkin(9) + t78;
t8 = -t40 * t134 + t48 * t135 - t89 * t30 - t85 * t31;
t145 = t85 * t57;
t95 = t63 * t134 + t145;
t7 = -t95 * pkin(10) - t8;
t154 = t88 * t7;
t74 = t81 * pkin(3) + pkin(9);
t153 = pkin(10) + t74;
t152 = t63 * t85;
t151 = t63 * t89;
t82 = sin(pkin(6));
t87 = sin(qJ(2));
t149 = t82 * t87;
t91 = cos(qJ(2));
t148 = t82 * t91;
t15 = -pkin(10) * t152 + t140;
t147 = t84 * t15;
t146 = t84 * t85;
t144 = t88 * t15;
t143 = t89 * t56;
t80 = t89 ^ 2;
t139 = t85 ^ 2 - t80;
t137 = qJD(2) * t87;
t136 = qJD(2) * t91;
t133 = qJD(6) * t84;
t132 = qJD(6) * t88;
t130 = t90 * qJD(3);
t128 = -0.2e1 * pkin(2) * qJD(3);
t75 = -t138 * pkin(3) - pkin(4);
t127 = t75 * t157;
t126 = pkin(5) * t135;
t125 = pkin(5) * t133;
t124 = pkin(5) * t132;
t123 = t91 * t131;
t121 = t82 * t137;
t120 = t82 * t136;
t119 = t85 * t134;
t115 = -t85 * t30 + t89 * t31;
t6 = -pkin(10) * t142 + t156 + (-t43 + (pkin(10) * t63 - t40) * t85) * qJD(5) + t115;
t118 = t88 * t6 - t84 * t7;
t114 = t89 * t40 - t85 * t48;
t14 = -pkin(10) * t151 + t114 + t155;
t117 = -t14 - t155;
t116 = -0.4e1 * t85 * t151;
t29 = -t138 * t93 + t81 * t55;
t47 = -t141 * t112 - t81 * t69;
t113 = qJD(5) * t153;
t109 = t139 * qJD(5);
t108 = t88 * t14 - t147;
t107 = t84 * t14 + t144;
t83 = cos(pkin(6));
t58 = t90 * t149 + t83 * t86;
t97 = -t86 * t149 + t83 * t90;
t34 = t138 * t58 + t81 * t97;
t98 = t85 * t148 - t89 * t34;
t99 = t89 * t148 + t85 * t34;
t106 = t84 * t98 - t88 * t99;
t105 = -t84 * t99 - t88 * t98;
t44 = t129 * t146 - t89 * t132 - t88 * t134;
t104 = t44 * t62 - t56 * t66;
t103 = -t56 * t74 + t57 * t75;
t59 = t153 * t85;
t60 = t153 * t89;
t102 = -t88 * t59 - t84 * t60;
t101 = -t84 * t59 + t88 * t60;
t100 = t62 * t74 - t63 * t75;
t65 = -t88 * t89 + t146;
t96 = t62 * t134 + t56 * t85;
t92 = -t58 * qJD(3) - t86 * t120;
t67 = -t89 * pkin(5) + t75;
t61 = t63 ^ 2;
t54 = t89 * t113;
t53 = t85 * t113;
t46 = t97 * qJD(3) + t90 * t120;
t45 = t129 * t66;
t42 = 0.2e1 * t62 * t56;
t39 = -t62 * t135 + t143;
t36 = t65 * t63;
t33 = -t138 * t97 + t81 * t58;
t28 = pkin(5) * t152 + t47;
t22 = t138 * t46 + t81 * t92;
t21 = -t138 * t92 + t81 * t46;
t19 = -t45 * t62 - t65 * t56;
t18 = -t101 * qJD(6) + t84 * t53 - t88 * t54;
t17 = -t102 * qJD(6) + t88 * t53 + t84 * t54;
t16 = t95 * pkin(5) + t29;
t13 = -t133 * t152 + (t129 * t151 + t145) * t88 + t94 * t84;
t12 = -t129 * t35 - t65 * t57;
t11 = t98 * qJD(5) + t89 * t121 - t85 * t22;
t10 = t99 * qJD(5) - t85 * t121 - t89 * t22;
t9 = -t140 * qJD(5) + t115;
t4 = -t105 * qJD(6) + t84 * t10 + t88 * t11;
t3 = -t106 * qJD(6) + t88 * t10 - t84 * t11;
t2 = -t107 * qJD(6) + t118;
t1 = -t108 * qJD(6) - t84 * t6 - t154;
t5 = [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -0.2e1 * t82 ^ 2 * t87 * t136 + 0.2e1 * t33 * t21 + 0.2e1 * t34 * t22, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, -t121, -t120, 0, 0, 0, 0, 0 (-t90 * t137 - t123) * t82 (-t91 * t130 + t86 * t137) * t82, t21 * t63 - t22 * t62 + t33 * t57 - t34 * t56, t21 * t47 + t22 * t48 + t33 * t29 + t34 * t30 + (-pkin(3) * t123 + t77 * t137) * t82, 0, 0, 0, 0, 0, t11 * t62 + t21 * t152 + t95 * t33 - t56 * t99, t10 * t62 + t21 * t151 + t94 * t33 + t56 * t98, 0, 0, 0, 0, 0, t106 * t56 + t33 * t13 + t21 * t35 + t4 * t62, -t105 * t56 + t33 * t12 - t21 * t36 + t3 * t62; 0, 0, 0, 0, 0.2e1 * t86 * t130, 0.2e1 * (-t86 ^ 2 + t90 ^ 2) * qJD(3), 0, 0, 0, t86 * t128, t90 * t128, 0.2e1 * t29 * t63 - 0.2e1 * t30 * t62 + 0.2e1 * t47 * t57 - 0.2e1 * t48 * t56, 0.2e1 * t29 * t47 + 0.2e1 * t30 * t48 + 0.2e1 * t77 * t78, 0.2e1 * t57 * t63 * t80 - 0.2e1 * t61 * t119, t139 * t61 * t157 + t57 * t116, 0.2e1 * t63 * t143 + 0.2e1 * t94 * t62, -0.2e1 * t56 * t152 - 0.2e1 * t95 * t62, t42, 0.2e1 * t114 * t56 + 0.2e1 * t29 * t152 + 0.2e1 * t95 * t47 + 0.2e1 * t9 * t62, -0.2e1 * t140 * t56 + 0.2e1 * t29 * t151 + 0.2e1 * t94 * t47 + 0.2e1 * t8 * t62, -0.2e1 * t36 * t12, -0.2e1 * t12 * t35 + 0.2e1 * t13 * t36, 0.2e1 * t12 * t62 - 0.2e1 * t36 * t56, -0.2e1 * t13 * t62 - 0.2e1 * t35 * t56, t42, 0.2e1 * t108 * t56 + 0.2e1 * t28 * t13 + 0.2e1 * t16 * t35 + 0.2e1 * t2 * t62, 0.2e1 * t1 * t62 - 0.2e1 * t107 * t56 + 0.2e1 * t28 * t12 - 0.2e1 * t16 * t36; 0, 0, 0, 0, 0, 0, 0, 0, 0, t92, -t46, 0 (-t138 * t21 + t22 * t81) * pkin(3), 0, 0, 0, 0, 0, t33 * t135 - t21 * t89, t33 * t134 + t21 * t85, 0, 0, 0, 0, 0, t21 * t65 + t33 * t45, t21 * t66 - t33 * t44; 0, 0, 0, 0, 0, 0, t130, -t131, 0, -pkin(8) * t130, pkin(8) * t131 (-t138 * t57 - t56 * t81) * pkin(3) (-t138 * t29 + t30 * t81) * pkin(3), -t63 * t109 + t85 * t142, qJD(5) * t116 - t139 * t57, t96, t39, 0, -t29 * t89 + t103 * t85 + (-t100 * t89 + t47 * t85) * qJD(5), t29 * t85 + t103 * t89 + (t100 * t85 + t47 * t89) * qJD(5), t12 * t66 + t36 * t44, -t12 * t65 - t13 * t66 + t35 * t44 + t36 * t45, -t104, t19, 0, t102 * t56 + t35 * t126 + t67 * t13 + t16 * t65 + t18 * t62 + t28 * t45, -t101 * t56 + t67 * t12 - t36 * t126 + t16 * t66 + t17 * t62 - t28 * t44; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0.2e1 * t119, -0.2e1 * t109, 0, 0, 0, t85 * t127, t89 * t127, -0.2e1 * t66 * t44, 0.2e1 * t44 * t65 - 0.2e1 * t45 * t66, 0, 0, 0, 0.2e1 * t65 * t126 + 0.2e1 * t45 * t67, 0.2e1 * t66 * t126 - 0.2e1 * t44 * t67; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t121, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t78, 0, 0, 0, 0, 0, t39, -t96, 0, 0, 0, 0, 0, t19, t104; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t11, t10, 0, 0, 0, 0, 0, t4, t3; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t94, -t95, t56, t9, t8, 0, 0, t12, -t13, t56, t88 * t156 + (t117 * t84 - t144) * qJD(6) + t118, -t154 + (-t6 - t156) * t84 + (t117 * t88 + t147) * qJD(6); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t134, -t135, 0, -t74 * t134, t74 * t135, 0, 0, -t44, -t45, 0, t18, t17; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t135, -t134, 0, 0, 0, 0, 0, -t45, t44; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -0.2e1 * t125, -0.2e1 * t124; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t4, t3; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t12, -t13, t56, t2, t1; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t44, -t45, 0, t18, t17; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t45, t44; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t125, -t124; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0;];
MMD_reg  = t5;
