% Calculate minimal parameter regressor of joint inertia matrix time derivative for
% S6RRRRPR1
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% qJD [6x1]
%   Generalized joint velocities
% pkin [11x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d2,d3,d4,d6,theta5]';
% 
% Output:
% MMD_reg [((6+1)*6/2)x33]
%   minimal parameter regressor of inerta matrix time derivative
%   (only lower left triangular matrix (including diagonal) due to symmetry

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-09 21:56
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function MMD_reg = S6RRRRPR1_inertiaDJ_regmin_slag_vp(qJ, qJD, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(11,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRRRPR1_inertiaDJ_regmin_slag_vp: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RRRRPR1_inertiaDJ_regmin_slag_vp: qJD has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [11 1]), ...
  'S6RRRRPR1_inertiaDJ_regmin_slag_vp: pkin has to be [11x1] (double)');

%% Symbolic Calculation
% From inertiaD_joint_joint_fixb_regressor_minpar_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 21:54:55
% EndTime: 2019-03-09 21:55:01
% DurationCPUTime: 1.82s
% Computational Cost: add. (6141->199), mult. (13555->354), div. (0->0), fcn. (13584->10), ass. (0->139)
t90 = sin(qJ(3));
t91 = sin(qJ(2));
t94 = cos(qJ(3));
t95 = cos(qJ(2));
t64 = t90 * t91 - t94 * t95;
t65 = t90 * t95 + t91 * t94;
t89 = sin(qJ(4));
t93 = cos(qJ(4));
t117 = t64 * t89 - t65 * t93;
t160 = pkin(7) + pkin(8);
t72 = t160 * t91;
t73 = t160 * t95;
t115 = -t72 * t94 - t73 * t90;
t108 = pkin(9) * t65 - t115;
t114 = t72 * t90 - t73 * t94;
t109 = pkin(9) * t64 + t114;
t98 = -t108 * t93 + t109 * t89;
t97 = qJ(5) * t117 + t98;
t144 = qJD(4) * t93;
t154 = t93 * t90;
t103 = (-t90 * t144 + (-t89 * t94 - t154) * qJD(3)) * pkin(2);
t145 = qJD(4) * t89;
t81 = pkin(2) * t94 + pkin(3);
t101 = -t145 * t81 + t103;
t92 = cos(qJ(6));
t86 = t92 ^ 2;
t88 = sin(qJ(6));
t151 = t88 ^ 2 - t86;
t128 = t151 * qJD(6);
t162 = qJD(2) + qJD(3);
t134 = qJD(2) * t160;
t66 = t91 * t134;
t67 = t95 * t134;
t116 = t94 * t66 + t90 * t67;
t47 = t162 * t65;
t110 = t47 * pkin(9) + t116;
t130 = t90 * t66 - t67 * t94;
t46 = t162 * t64;
t113 = t46 * pkin(9) + t130;
t161 = (t114 * t89 + t115 * t93) * qJD(3) - t93 * t110 + t89 * t113;
t119 = t89 * t46 - t93 * t47;
t105 = -qJD(4) * t117 - t119;
t148 = cos(pkin(11));
t118 = t64 * t93 + t65 * t89;
t99 = t108 * t89 + t109 * t93;
t24 = -qJ(5) * t118 - t99;
t87 = sin(pkin(11));
t17 = -t148 * t97 + t24 * t87;
t9 = t119 * qJ(5) + qJD(4) * t97 - t118 * qJD(5) + t161;
t14 = t89 * t110 + t93 * t113 + t99 * qJD(4) + (t114 * t93 - t115 * t89) * qJD(3);
t26 = -qJD(4) * t118 - t93 * t46 - t89 * t47;
t96 = -t26 * qJ(5) + qJD(5) * t117 + t14;
t4 = -t148 * t96 + t87 * t9;
t84 = qJD(6) * t92;
t159 = t17 * t84 + t4 * t88;
t32 = -t117 * t148 - t118 * t87;
t158 = t32 * t88;
t157 = t32 * t92;
t156 = t87 * t89;
t20 = -t105 * t87 + t148 * t26;
t155 = t92 * t20;
t150 = pkin(2) * qJD(3);
t139 = t94 * t150;
t142 = pkin(2) * t89 * t90;
t44 = -t81 * t144 - t93 * t139 + (qJD(3) + qJD(4)) * t142;
t28 = -t101 * t148 - t44 * t87;
t58 = t81 * t93 + pkin(4) - t142;
t61 = pkin(2) * t154 + t81 * t89;
t39 = t148 * t58 - t61 * t87;
t37 = -pkin(5) - t39;
t153 = t28 * t88 + t37 * t84;
t80 = pkin(3) * t93 + pkin(4);
t59 = -pkin(3) * t156 + t148 * t80;
t54 = -pkin(5) - t59;
t129 = t148 * t89;
t149 = pkin(3) * qJD(4);
t56 = (t87 * t93 + t129) * t149;
t152 = t54 * t84 + t56 * t88;
t40 = t148 * t61 + t58 * t87;
t60 = pkin(3) * t129 + t80 * t87;
t147 = qJD(2) * t91;
t146 = qJD(2) * t95;
t143 = qJD(6) * t88;
t141 = -0.2e1 * pkin(1) * qJD(2);
t83 = pkin(2) * t147;
t140 = t90 * t150;
t138 = pkin(3) * t145;
t137 = pkin(3) * t144;
t135 = t88 * t84;
t82 = -t95 * pkin(2) - pkin(1);
t42 = t47 * pkin(3) + t83;
t133 = -0.4e1 * t88 * t157;
t35 = t37 * t143;
t132 = -t28 * t92 + t35;
t48 = t54 * t143;
t131 = -t56 * t92 + t48;
t18 = t148 * t24 + t87 * t97;
t53 = t64 * pkin(3) + t82;
t104 = pkin(4) * t118 + t53;
t31 = -t117 * t87 + t118 * t148;
t21 = t31 * pkin(5) - t32 * pkin(10) + t104;
t127 = t18 * t92 + t21 * t88;
t126 = t18 * t88 - t21 * t92;
t19 = t105 * t148 + t26 * t87;
t77 = pkin(4) * t87 + pkin(10);
t78 = -pkin(4) * t148 - pkin(5);
t125 = -t19 * t77 + t20 * t78;
t29 = t101 * t87 - t148 * t44;
t124 = t28 * t32 - t29 * t31;
t38 = pkin(10) + t40;
t123 = t31 * t38 - t32 * t37;
t55 = pkin(10) + t60;
t122 = t31 * t55 - t32 * t54;
t57 = (t148 * t93 - t156) * t149;
t121 = -t31 * t57 + t32 * t56;
t120 = t31 * t77 - t32 * t78;
t12 = t19 * t88 + t31 * t84;
t112 = t20 * t88 + t32 * t84;
t111 = -t143 * t32 + t155;
t107 = -t19 * t38 + t20 * t37 + t124;
t106 = -t19 * t55 + t20 * t54 + t121;
t22 = pkin(4) * t105 + t42;
t75 = 0.2e1 * t135;
t69 = t78 * t84;
t68 = t78 * t143;
t63 = -0.2e1 * t128;
t34 = qJD(3) * t114 + t130;
t33 = -qJD(3) * t115 + t116;
t30 = t32 ^ 2;
t15 = t17 * t143;
t13 = -qJD(4) * t98 - t161;
t11 = -t143 * t31 + t19 * t92;
t10 = -t128 * t32 + t155 * t88;
t7 = qJD(6) * t133 - t151 * t20;
t6 = t19 * pkin(5) - t20 * pkin(10) + t22;
t5 = t148 * t9 + t87 * t96;
t2 = -qJD(6) * t127 - t5 * t88 + t6 * t92;
t1 = qJD(6) * t126 - t5 * t92 - t6 * t88;
t3 = [0, 0, 0, 0.2e1 * t91 * t146, 0.2e1 * (-t91 ^ 2 + t95 ^ 2) * qJD(2), 0, 0, 0, t91 * t141, t95 * t141, -0.2e1 * t65 * t46, 0.2e1 * t46 * t64 - 0.2e1 * t47 * t65, 0, 0, 0, 0.2e1 * t47 * t82 + 0.2e1 * t64 * t83, -0.2e1 * t46 * t82 + 0.2e1 * t65 * t83, -0.2e1 * t117 * t26, 0.2e1 * t105 * t117 - 0.2e1 * t118 * t26, 0, 0, 0, 0.2e1 * t105 * t53 + 0.2e1 * t118 * t42, -0.2e1 * t117 * t42 + 0.2e1 * t26 * t53, 0.2e1 * t17 * t20 - 0.2e1 * t18 * t19 - 0.2e1 * t31 * t5 + 0.2e1 * t32 * t4, 0.2e1 * t104 * t22 + 0.2e1 * t17 * t4 + 0.2e1 * t18 * t5, 0.2e1 * t20 * t32 * t86 - 0.2e1 * t135 * t30, 0.2e1 * t128 * t30 + t133 * t20, 0.2e1 * t111 * t31 + 0.2e1 * t157 * t19, -0.2e1 * t112 * t31 - 0.2e1 * t158 * t19, 0.2e1 * t31 * t19, 0.2e1 * t112 * t17 - 0.2e1 * t126 * t19 + 0.2e1 * t158 * t4 + 0.2e1 * t2 * t31, 0.2e1 * t1 * t31 + 0.2e1 * t111 * t17 - 0.2e1 * t127 * t19 + 0.2e1 * t157 * t4; 0, 0, 0, 0, 0, t146, -t147, 0, -pkin(7) * t146, pkin(7) * t147, 0, 0, -t46, -t47, 0, t34, t33, 0, 0, t26, -t105, 0, t14, t13, -t19 * t40 - t20 * t39 + t124, t17 * t28 + t18 * t29 - t39 * t4 + t40 * t5, t10, t7, t12, t11, 0, t15 + (-qJD(6) * t123 - t4) * t92 + t107 * t88, t107 * t92 + t123 * t143 + t159; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -0.2e1 * t140, -0.2e1 * t139, 0, 0, 0, 0, 0, 0.2e1 * t101, 0.2e1 * t44, 0, -0.2e1 * t28 * t39 + 0.2e1 * t29 * t40, t75, t63, 0, 0, 0, 0.2e1 * t132, 0.2e1 * t153; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t46, -t47, 0, t34, t33, 0, 0, t26, -t105, 0, t14, t13, -t19 * t60 - t20 * t59 + t121, t17 * t56 + t18 * t57 - t4 * t59 + t5 * t60, t10, t7, t12, t11, 0, t15 + (-qJD(6) * t122 - t4) * t92 + t106 * t88, t106 * t92 + t122 * t143 + t159; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t140, -t139, 0, 0, 0, 0, 0 (-pkin(3) - t81) * t145 + t103, t44 - t137, 0, -t28 * t59 + t29 * t60 - t39 * t56 + t40 * t57, t75, t63, 0, 0, 0, t35 + t48 + (-t28 - t56) * t92, t152 + t153; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -0.2e1 * t138, -0.2e1 * t137, 0, -0.2e1 * t56 * t59 + 0.2e1 * t57 * t60, t75, t63, 0, 0, 0, 0.2e1 * t131, 0.2e1 * t152; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t26, -t105, 0, t14, t13 (-t148 * t20 - t19 * t87) * pkin(4) (-t148 * t4 + t5 * t87) * pkin(4), t10, t7, t12, t11, 0, t15 + t125 * t88 + (-qJD(6) * t120 - t4) * t92, t120 * t143 + t125 * t92 + t159; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t101, t44, 0 (-t148 * t28 + t29 * t87) * pkin(4), t75, t63, 0, 0, 0, t132 + t68, t69 + t153; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t138, -t137, 0 (-t148 * t56 + t57 * t87) * pkin(4), t75, t63, 0, 0, 0, t131 + t68, t69 + t152; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t75, t63, 0, 0, 0, 0.2e1 * t68, 0.2e1 * t69; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t22, 0, 0, 0, 0, 0, t11, -t12; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t111, -t112, t19, t2, t1; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t84, -t143, 0, -t29 * t88 - t38 * t84, t143 * t38 - t29 * t92; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t84, -t143, 0, -t55 * t84 - t57 * t88, t143 * t55 - t57 * t92; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t84, -t143, 0, -t77 * t84, t77 * t143; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t143, -t84; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0;];
MMD_reg  = t3;
