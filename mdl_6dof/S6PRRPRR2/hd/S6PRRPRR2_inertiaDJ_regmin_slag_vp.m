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
% MMD_reg [((6+1)*6/2)x29]
%   minimal parameter regressor of inertia matrix time derivative
%   (only lower left triangular matrix (including diagonal) due to symmetry

% Quelle: HybrDyn-Toolbox
% Datum: 2021-01-16 03:49
% Revision: 24b2e7d74a0c1a3b64fa2f8f5ad758691ad61af3 (2021-01-15)
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
% StartTime: 2021-01-16 03:46:37
% EndTime: 2021-01-16 03:46:45
% DurationCPUTime: 1.55s
% Computational Cost: add. (1853->206), mult. (4733->391), div. (0->0), fcn. (4749->12), ass. (0->136)
t138 = cos(pkin(12));
t85 = sin(qJ(3));
t111 = t138 * t85;
t81 = sin(pkin(12));
t89 = cos(qJ(3));
t63 = t81 * t89 + t111;
t83 = sin(qJ(6));
t84 = sin(qJ(5));
t87 = cos(qJ(6));
t88 = cos(qJ(5));
t66 = t83 * t88 + t87 * t84;
t35 = t66 * t63;
t110 = t138 * t89;
t152 = t81 * t85;
t62 = -t110 + t152;
t77 = -t89 * pkin(3) - pkin(2);
t40 = t62 * pkin(4) - t63 * pkin(9) + t77;
t142 = qJ(4) + pkin(8);
t69 = t142 * t89;
t48 = t138 * t69 - t142 * t152;
t43 = t88 * t48;
t141 = t84 * t40 + t43;
t132 = qJD(5) * t84;
t135 = qJD(3) * t85;
t57 = qJD(3) * t110 - t81 * t135;
t143 = t88 * t57;
t94 = -t63 * t132 + t143;
t128 = qJD(5) + qJD(6);
t159 = 0.2e1 * qJD(5);
t56 = t63 * qJD(3);
t158 = t56 * pkin(5);
t157 = t62 * pkin(5);
t131 = qJD(5) * t88;
t109 = qJD(3) * t142;
t55 = t89 * qJD(4) - t85 * t109;
t92 = -t85 * qJD(4) - t89 * t109;
t30 = t138 * t55 + t81 * t92;
t78 = pkin(3) * t135;
t31 = t56 * pkin(4) - t57 * pkin(9) + t78;
t8 = -t40 * t131 + t48 * t132 - t88 * t30 - t84 * t31;
t146 = t84 * t57;
t95 = t63 * t131 + t146;
t7 = -t95 * pkin(10) - t8;
t156 = t87 * t7;
t74 = t81 * pkin(3) + pkin(9);
t155 = pkin(10) + t74;
t154 = t63 * t84;
t153 = t63 * t88;
t82 = sin(pkin(6));
t86 = sin(qJ(2));
t151 = t82 * t86;
t90 = cos(qJ(2));
t150 = t82 * t90;
t15 = -pkin(10) * t154 + t141;
t149 = t83 * t15;
t148 = t83 * t84;
t147 = t84 * t56;
t145 = t87 * t15;
t144 = t88 * t56;
t80 = t88 ^ 2;
t140 = t84 ^ 2 - t80;
t139 = cos(pkin(6));
t137 = qJD(2) * t86;
t136 = qJD(2) * t90;
t134 = qJD(3) * t89;
t133 = qJD(3) * t90;
t130 = qJD(6) * t83;
t129 = qJD(6) * t87;
t127 = -0.2e1 * pkin(2) * qJD(3);
t75 = -t138 * pkin(3) - pkin(4);
t126 = t75 * t159;
t125 = pkin(5) * t132;
t124 = pkin(5) * t130;
t123 = pkin(5) * t129;
t122 = t85 * t133;
t120 = t82 * t137;
t119 = t82 * t136;
t118 = t84 * t131;
t114 = -t84 * t30 + t88 * t31;
t6 = -pkin(10) * t143 + t158 + (-t43 + (pkin(10) * t63 - t40) * t84) * qJD(5) + t114;
t117 = t87 * t6 - t83 * t7;
t113 = t88 * t40 - t84 * t48;
t14 = -pkin(10) * t153 + t113 + t157;
t116 = -t14 - t157;
t115 = -0.4e1 * t84 * t153;
t29 = -t138 * t92 + t81 * t55;
t47 = t142 * t111 + t81 * t69;
t112 = qJD(5) * t155;
t108 = t140 * qJD(5);
t107 = t87 * t14 - t149;
t106 = t83 * t14 + t145;
t58 = t139 * t85 + t89 * t151;
t93 = t139 * t89 - t85 * t151;
t34 = t138 * t58 + t81 * t93;
t97 = t84 * t150 - t88 * t34;
t98 = t88 * t150 + t84 * t34;
t105 = t83 * t97 - t87 * t98;
t104 = -t83 * t98 - t87 * t97;
t44 = t128 * t148 - t88 * t129 - t87 * t131;
t103 = t44 * t62 - t66 * t56;
t102 = -t56 * t74 + t57 * t75;
t59 = t155 * t84;
t60 = t155 * t88;
t101 = -t87 * t59 - t83 * t60;
t100 = -t83 * t59 + t87 * t60;
t99 = t62 * t74 - t63 * t75;
t65 = -t87 * t88 + t148;
t96 = t62 * t131 + t147;
t91 = -t58 * qJD(3) - t85 * t119;
t67 = -t88 * pkin(5) + t75;
t61 = t63 ^ 2;
t54 = t88 * t112;
t53 = t84 * t112;
t46 = t93 * qJD(3) + t89 * t119;
t45 = t128 * t66;
t42 = 0.2e1 * t62 * t56;
t39 = -t62 * t132 + t144;
t36 = t65 * t63;
t33 = -t138 * t93 + t81 * t58;
t28 = pkin(5) * t154 + t47;
t22 = t138 * t46 + t81 * t91;
t21 = -t138 * t91 + t81 * t46;
t19 = -t45 * t62 - t65 * t56;
t18 = -t100 * qJD(6) + t83 * t53 - t87 * t54;
t17 = -t101 * qJD(6) + t87 * t53 + t83 * t54;
t16 = t95 * pkin(5) + t29;
t13 = -t130 * t154 + (t128 * t153 + t146) * t87 + t94 * t83;
t12 = -t128 * t35 - t65 * t57;
t11 = t97 * qJD(5) + t88 * t120 - t84 * t22;
t10 = t98 * qJD(5) - t84 * t120 - t88 * t22;
t9 = -t141 * qJD(5) + t114;
t4 = -t104 * qJD(6) + t83 * t10 + t87 * t11;
t3 = -t105 * qJD(6) + t87 * t10 - t83 * t11;
t2 = -t106 * qJD(6) + t117;
t1 = -t107 * qJD(6) - t83 * t6 - t156;
t5 = [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -0.2e1 * t82 ^ 2 * t86 * t136 + 0.2e1 * t33 * t21 + 0.2e1 * t34 * t22, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, -t120, -t119, 0, 0, 0, 0, 0, (-t89 * t137 - t122) * t82, (-t89 * t133 + t85 * t137) * t82, (t62 * t137 - t56 * t90) * t82, (t63 * t137 - t57 * t90) * t82, t21 * t63 - t22 * t62 + t33 * t57 - t34 * t56, t21 * t47 + t22 * t48 + t33 * t29 + t34 * t30 + (-pkin(3) * t122 + t77 * t137) * t82, 0, 0, 0, 0, 0, t11 * t62 + t21 * t154 + t33 * t95 - t56 * t98, t10 * t62 + t21 * t153 + t33 * t94 + t56 * t97, 0, 0, 0, 0, 0, t105 * t56 + t33 * t13 + t21 * t35 + t4 * t62, -t104 * t56 + t33 * t12 - t21 * t36 + t3 * t62; 0, 0, 0, 0, 0.2e1 * t85 * t134, 0.2e1 * (-t85 ^ 2 + t89 ^ 2) * qJD(3), 0, 0, 0, t85 * t127, t89 * t127, 0.2e1 * t77 * t56 + 0.2e1 * t62 * t78, 0.2e1 * t77 * t57 + 0.2e1 * t63 * t78, 0.2e1 * t29 * t63 - 0.2e1 * t30 * t62 + 0.2e1 * t47 * t57 - 0.2e1 * t48 * t56, 0.2e1 * t47 * t29 + 0.2e1 * t48 * t30 + 0.2e1 * t77 * t78, 0.2e1 * t80 * t63 * t57 - 0.2e1 * t118 * t61, t140 * t61 * t159 + t57 * t115, 0.2e1 * t63 * t144 + 0.2e1 * t62 * t94, -0.2e1 * t63 * t147 - 0.2e1 * t62 * t95, t42, 0.2e1 * t113 * t56 + 0.2e1 * t29 * t154 + 0.2e1 * t47 * t95 + 0.2e1 * t9 * t62, -0.2e1 * t141 * t56 + 0.2e1 * t29 * t153 + 0.2e1 * t94 * t47 + 0.2e1 * t8 * t62, -0.2e1 * t36 * t12, -0.2e1 * t12 * t35 + 0.2e1 * t36 * t13, 0.2e1 * t12 * t62 - 0.2e1 * t36 * t56, -0.2e1 * t13 * t62 - 0.2e1 * t35 * t56, t42, 0.2e1 * t107 * t56 + 0.2e1 * t28 * t13 + 0.2e1 * t16 * t35 + 0.2e1 * t2 * t62, 0.2e1 * t1 * t62 - 0.2e1 * t106 * t56 + 0.2e1 * t28 * t12 - 0.2e1 * t16 * t36; 0, 0, 0, 0, 0, 0, 0, 0, 0, t91, -t46, -t21, -t22, 0, (-t138 * t21 + t22 * t81) * pkin(3), 0, 0, 0, 0, 0, t33 * t132 - t21 * t88, t33 * t131 + t21 * t84, 0, 0, 0, 0, 0, t21 * t65 + t33 * t45, t21 * t66 - t33 * t44; 0, 0, 0, 0, 0, 0, t134, -t135, 0, -pkin(8) * t134, pkin(8) * t135, -t29, -t30, (-t138 * t57 - t56 * t81) * pkin(3), (-t138 * t29 + t30 * t81) * pkin(3), -t108 * t63 + t84 * t143, qJD(5) * t115 - t140 * t57, t96, t39, 0, -t29 * t88 + t102 * t84 + (t47 * t84 - t88 * t99) * qJD(5), t29 * t84 + t102 * t88 + (t47 * t88 + t84 * t99) * qJD(5), t12 * t66 + t36 * t44, -t12 * t65 - t66 * t13 + t44 * t35 + t36 * t45, -t103, t19, 0, t101 * t56 + t125 * t35 + t67 * t13 + t16 * t65 + t18 * t62 + t28 * t45, -t100 * t56 + t67 * t12 - t125 * t36 + t16 * t66 + t17 * t62 - t28 * t44; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0.2e1 * t118, -0.2e1 * t108, 0, 0, 0, t84 * t126, t88 * t126, -0.2e1 * t66 * t44, 0.2e1 * t44 * t65 - 0.2e1 * t45 * t66, 0, 0, 0, 0.2e1 * t125 * t65 + 0.2e1 * t45 * t67, 0.2e1 * t125 * t66 - 0.2e1 * t44 * t67; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t120, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t56, t57, 0, t78, 0, 0, 0, 0, 0, t39, -t96, 0, 0, 0, 0, 0, t19, t103; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t11, t10, 0, 0, 0, 0, 0, t4, t3; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t94, -t95, t56, t9, t8, 0, 0, t12, -t13, t56, t87 * t158 + (t116 * t83 - t145) * qJD(6) + t117, -t156 + (-t6 - t158) * t83 + (t116 * t87 + t149) * qJD(6); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t131, -t132, 0, -t74 * t131, t74 * t132, 0, 0, -t44, -t45, 0, t18, t17; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t132, -t131, 0, 0, 0, 0, 0, -t45, t44; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -0.2e1 * t124, -0.2e1 * t123; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t4, t3; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t12, -t13, t56, t2, t1; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t44, -t45, 0, t18, t17; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t45, t44; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t124, -t123; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0;];
MMD_reg = t5;
