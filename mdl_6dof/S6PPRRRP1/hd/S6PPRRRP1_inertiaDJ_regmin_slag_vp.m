% Calculate minimal parameter regressor of joint inertia matrix time derivative for
% S6PPRRRP1
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% qJD [6x1]
%   Generalized joint velocities
% pkin [12x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,alpha2,alpha3,d3,d4,d5,theta1,theta2]';
% 
% Output:
% MMD_reg [((6+1)*6/2)x23]
%   minimal parameter regressor of inertia matrix time derivative
%   (only lower left triangular matrix (including diagonal) due to symmetry

% Quelle: HybrDyn-Toolbox
% Datum: 2021-01-16 00:51
% Revision: 24b2e7d74a0c1a3b64fa2f8f5ad758691ad61af3 (2021-01-15)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function MMD_reg = S6PPRRRP1_inertiaDJ_regmin_slag_vp(qJ, qJD, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(12,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6PPRRRP1_inertiaDJ_regmin_slag_vp: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6PPRRRP1_inertiaDJ_regmin_slag_vp: qJD has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [12 1]), ...
  'S6PPRRRP1_inertiaDJ_regmin_slag_vp: pkin has to be [12x1] (double)');

%% Symbolic Calculation
% From inertiaD_joint_joint_fixb_regressor_minpar_matlab.m
% OptimizationMode: 2
% StartTime: 2021-01-16 00:50:11
% EndTime: 2021-01-16 00:50:20
% DurationCPUTime: 1.90s
% Computational Cost: add. (1425->202), mult. (4382->378), div. (0->0), fcn. (4490->12), ass. (0->128)
t145 = sin(qJ(3));
t62 = cos(qJ(3));
t128 = cos(pkin(7));
t129 = cos(pkin(6));
t57 = sin(pkin(7));
t127 = sin(pkin(6));
t84 = cos(pkin(12)) * t127;
t74 = t128 * t84 + t129 * t57;
t83 = t127 * sin(pkin(12));
t150 = -t145 * t83 + t74 * t62;
t68 = qJD(3) * t150;
t59 = sin(qJ(4));
t152 = -0.4e1 * t59;
t61 = cos(qJ(4));
t70 = -t74 * t145 - t62 * t83;
t73 = t129 * t128 - t57 * t84;
t17 = -t59 * t70 - t61 * t73;
t63 = -qJD(4) * t17 + t61 * t68;
t151 = qJD(5) * t150 - t63;
t140 = t57 * t62;
t102 = qJD(3) * t140;
t99 = t57 * t145;
t34 = -t61 * t128 + t59 * t99;
t72 = -qJD(4) * t34 + t61 * t102;
t64 = t70 * qJD(3);
t149 = t59 * qJD(6) + (pkin(9) * qJD(5) + qJ(6) * qJD(4)) * t61;
t148 = 0.2e1 * qJD(5);
t60 = cos(qJ(5));
t147 = pkin(5) * t60;
t58 = sin(qJ(5));
t146 = pkin(9) * t58;
t121 = qJD(4) * t61;
t111 = pkin(9) * t121;
t52 = qJD(5) * t60;
t37 = t58 * t121 + t59 * t52;
t31 = t37 * pkin(5) + t111;
t144 = t31 * t58;
t143 = t31 * t60;
t136 = -qJ(6) - pkin(10);
t45 = t136 * t58;
t142 = t45 * t59;
t46 = t136 * t60;
t141 = t46 * t59;
t139 = t58 * t61;
t138 = t59 * t60;
t137 = t60 * t61;
t89 = -pkin(4) * t61 - pkin(10) * t59;
t82 = -pkin(3) + t89;
t78 = qJD(5) * t82;
t88 = pkin(4) * t59 - pkin(10) * t61;
t79 = t88 * qJD(4);
t135 = -t58 * t79 - t60 * t78;
t123 = qJD(4) * t59;
t110 = t58 * t123;
t134 = pkin(9) * t110 + t60 * t79;
t49 = pkin(9) * t137;
t133 = t58 * t82 + t49;
t53 = t58 ^ 2;
t55 = t60 ^ 2;
t132 = t53 - t55;
t54 = t59 ^ 2;
t131 = -t61 ^ 2 + t54;
t130 = qJ(6) * t59;
t124 = qJD(4) * t58;
t122 = qJD(4) * t60;
t120 = qJD(4) * t62;
t119 = qJD(5) * t58;
t118 = qJD(5) * t61;
t116 = -0.2e1 * pkin(3) * qJD(4);
t115 = -0.2e1 * pkin(4) * qJD(5);
t114 = t58 * t140;
t113 = pkin(5) * t123;
t112 = pkin(5) * t119;
t109 = t60 * t121;
t108 = t17 * t119;
t107 = t34 * t119;
t43 = (pkin(5) * t58 + pkin(9)) * t59;
t106 = t43 * t119;
t105 = t59 * t119;
t104 = t58 * t118;
t103 = t60 * t118;
t101 = t58 * t52;
t100 = t59 * t121;
t50 = -pkin(4) - t147;
t98 = -t50 + t147;
t97 = qJD(3) * t145;
t95 = t132 * qJD(5);
t94 = t131 * qJD(4);
t93 = 0.2e1 * t100;
t91 = t58 * t109;
t90 = t57 * t97;
t87 = pkin(5) * t53 + t50 * t60;
t18 = t59 * t73 - t61 * t70;
t11 = -t150 * t60 - t58 * t18;
t12 = -t150 * t58 + t60 * t18;
t86 = -t11 * t60 - t12 * t58;
t35 = t59 * t128 + t61 * t99;
t29 = t60 * t35 - t114;
t80 = t60 * t140 + t58 * t35;
t85 = -t29 * t58 + t60 * t80;
t10 = t18 * qJD(4) + t59 * t68;
t7 = t10 * t58 + t17 * t52;
t8 = -t10 * t60 + t108;
t27 = t35 * qJD(4) + t59 * t102;
t19 = t27 * t58 + t34 * t52;
t20 = -t27 * t60 + t107;
t36 = -t105 + t109;
t77 = t59 * t122 + t104;
t69 = qJ(6) * t105 - t149 * t60 - t58 * t78 + t134;
t67 = qJD(4) * t150;
t42 = t60 * t82;
t33 = -t58 * qJD(6) + t136 * t52;
t32 = -t60 * qJD(6) - t136 * t119;
t30 = -t58 * t130 + t133;
t26 = -t60 * t130 + t42 + (-pkin(5) - t146) * t61;
t22 = -t133 * qJD(5) + t134;
t21 = t77 * pkin(9) + t135;
t16 = (pkin(9) * qJD(4) + qJ(6) * qJD(5)) * t138 + t149 * t58 + t135;
t15 = qJD(5) * t114 - t35 * t52 - t58 * t72 + t60 * t90;
t14 = t80 * qJD(5) - t58 * t90 - t60 * t72;
t13 = t69 + t113;
t6 = (t34 * t124 - t15) * t61 + (-qJD(4) * t80 + t19) * t59;
t5 = (t34 * t122 - t14) * t61 + (-qJD(4) * t29 - t20) * t59;
t4 = t151 * t58 - t18 * t52 - t60 * t64;
t3 = t18 * t119 + t151 * t60 + t58 * t64;
t2 = (t17 * t124 - t4) * t61 + (qJD(4) * t11 + t7) * t59;
t1 = (t17 * t122 - t3) * t61 + (-qJD(4) * t12 - t8) * t59;
t9 = [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0.2e1 * t10 * t17 + 0.2e1 * t11 * t4 - 0.2e1 * t12 * t3; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t10 * t34 + t11 * t15 - t12 * t14 + t17 * t27 - t29 * t3 - t4 * t80; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -0.2e1 * t14 * t29 - 0.2e1 * t15 * t80 + 0.2e1 * t27 * t34; 0, 0, 0, t64, -t68, 0, 0, 0, 0, 0, -t59 * t67 + t61 * t64, -t59 * t64 - t61 * t67, 0, 0, 0, 0, 0, t2, t1, t2, t1, t86 * t121 + (t3 * t58 - t4 * t60 + (t11 * t58 - t12 * t60) * qJD(5)) * t59, t10 * t43 + t11 * t13 - t12 * t16 + t17 * t31 + t26 * t4 - t3 * t30; 0, 0, 0, -t90, -t102, 0, 0, 0, 0, 0, (-t59 * t120 - t61 * t97) * t57, (-t61 * t120 + t59 * t97) * t57, 0, 0, 0, 0, 0, t6, t5, t6, t5, t85 * t121 + (t14 * t58 - t15 * t60 + (-t29 * t60 - t58 * t80) * qJD(5)) * t59, -t13 * t80 - t14 * t30 + t15 * t26 - t16 * t29 + t27 * t43 + t31 * t34; 0, 0, 0, 0, 0, t93, -0.2e1 * t94, 0, 0, 0, t59 * t116, t61 * t116, 0.2e1 * t55 * t100 - 0.2e1 * t54 * t101, t132 * t54 * t148 + t91 * t152, 0.2e1 * t59 * t104 + 0.2e1 * t131 * t122, 0.2e1 * t59 * t103 - 0.2e1 * t58 * t94, -0.2e1 * t100, 0.2e1 * t42 * t123 - 0.2e1 * t22 * t61 + 0.2e1 * (t100 * t58 + t54 * t52) * pkin(9), -0.2e1 * t21 * t61 - 0.2e1 * t133 * t123 + 0.2e1 * (-t54 * t119 + t60 * t93) * pkin(9), 0.2e1 * (t43 * t124 - t13) * t61 + 0.2e1 * (qJD(4) * t26 + t43 * t52 + t144) * t59, 0.2e1 * (t43 * t122 - t16) * t61 + 0.2e1 * (-qJD(4) * t30 - t106 + t143) * t59, 0.2e1 * (-t26 * t60 - t30 * t58) * t121 + 0.2e1 * (-t13 * t60 + t16 * t58 + (t26 * t58 - t30 * t60) * qJD(5)) * t59, 0.2e1 * t13 * t26 - 0.2e1 * t16 * t30 + 0.2e1 * t31 * t43; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t10, -t63, 0, 0, 0, 0, 0, t8, t7, t8, t7, qJD(5) * t86 - t3 * t60 - t4 * t58, pkin(5) * t108 + t10 * t50 + t11 * t33 - t12 * t32 + t3 * t46 + t4 * t45; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t27, -t72, 0, 0, 0, 0, 0, t20, t19, t20, t19, qJD(5) * t85 - t14 * t60 - t15 * t58, pkin(5) * t107 + t14 * t46 + t15 * t45 + t27 * t50 - t29 * t32 - t33 * t80; 0, 0, 0, 0, 0, 0, 0, t121, -t123, 0, -t111, pkin(9) * t123, -t59 * t95 + t91, t101 * t152 - t132 * t121, -t103 + t110, t77, 0, (pkin(10) * t137 + (-pkin(4) * t60 + t146) * t59) * qJD(5) + (t58 * t89 - t49) * qJD(4), (pkin(9) * t138 + t58 * t88) * qJD(5) + (pkin(9) * t139 + t60 * t89) * qJD(4), -t143 - t33 * t61 + (t50 * t139 + t142) * qJD(4) + (t43 * t58 + t59 * t87) * qJD(5), t144 - t32 * t61 + (t50 * t137 + t141) * qJD(4) + (t58 * t59 * t98 + t43 * t60) * qJD(5), (-t45 * t121 - t33 * t59 - t16 + (-t26 + t141) * qJD(5)) * t60 + (t46 * t121 + t32 * t59 - t13 + (-t30 + t142) * qJD(5)) * t58, pkin(5) * t106 + t13 * t45 + t16 * t46 + t26 * t33 - t30 * t32 + t31 * t50; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0.2e1 * t101, -0.2e1 * t95, 0, 0, 0, t58 * t115, t60 * t115, -0.2e1 * t98 * t119, t87 * t148, -0.2e1 * t32 * t60 - 0.2e1 * t33 * t58 + 0.2e1 * (-t45 * t60 + t46 * t58) * qJD(5), 0.2e1 * t112 * t50 + 0.2e1 * t32 * t46 + 0.2e1 * t33 * t45; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t4, t3, t4, t3, 0, t4 * pkin(5); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t15, t14, t15, t14, 0, t15 * pkin(5); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t36, -t37, t123, t22, t21, t69 + 0.2e1 * t113, t16, -t36 * pkin(5), t13 * pkin(5); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t52, -t119, 0, -pkin(10) * t52, pkin(10) * t119, t33, t32, -pkin(5) * t52, t33 * pkin(5); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t10; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t27; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t37, t36, 0, t31; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t119, t52, 0, t112; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0;];
MMD_reg = t9;
