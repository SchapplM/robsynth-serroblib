% Calculate minimal parameter regressor of joint inertia matrix time derivative for
% S6RPRRRP8
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
% Datum: 2019-03-09 06:25
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function MMD_reg = S6RPRRRP8_inertiaDJ_regmin_slag_vp(qJ, qJD, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(9,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPRRRP8_inertiaDJ_regmin_slag_vp: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RPRRRP8_inertiaDJ_regmin_slag_vp: qJD has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S6RPRRRP8_inertiaDJ_regmin_slag_vp: pkin has to be [9x1] (double)');

%% Symbolic Calculation
% From inertiaD_joint_joint_fixb_regressor_minpar_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 06:24:34
% EndTime: 2019-03-09 06:24:39
% DurationCPUTime: 1.56s
% Computational Cost: add. (2600->198), mult. (5174->336), div. (0->0), fcn. (4648->6), ass. (0->127)
t156 = cos(qJ(4));
t113 = t156 * qJD(3);
t155 = sin(qJ(4));
t115 = qJD(3) * t155;
t85 = sin(qJ(3));
t87 = cos(qJ(3));
t169 = -t113 * t85 - t115 * t87;
t84 = sin(qJ(5));
t136 = qJD(5) * t84;
t107 = t85 * t115;
t114 = qJD(4) * t155;
t112 = t156 * qJD(4);
t163 = t113 + t112;
t36 = -t114 * t85 + t163 * t87 - t107;
t58 = t155 * t87 + t156 * t85;
t147 = t58 * t36;
t35 = -t112 * t85 - t114 * t87 + t169;
t57 = t155 * t85 - t156 * t87;
t150 = t57 * t35;
t54 = t58 ^ 2;
t55 = t57 ^ 2;
t86 = cos(qJ(5));
t168 = (t54 + t55) * t136 + 0.2e1 * (-t147 + t150) * t86;
t74 = pkin(3) * t85 + qJ(2);
t30 = pkin(4) * t58 + pkin(9) * t57 + t74;
t88 = -pkin(1) - pkin(7);
t157 = pkin(8) - t88;
t60 = t157 * t85;
t61 = t157 * t87;
t39 = -t155 * t61 - t156 * t60;
t167 = t30 * t84 + t39 * t86;
t82 = t84 ^ 2;
t83 = t86 ^ 2;
t165 = t82 + t83;
t140 = t82 - t83;
t111 = qJD(5) * t140;
t132 = t87 * qJD(3);
t67 = pkin(3) * t132 + qJD(2);
t13 = pkin(4) * t36 - pkin(9) * t35 + t67;
t18 = -t157 * t107 - t60 * t114 + t163 * t61;
t5 = -qJD(5) * t167 + t13 * t86 + t18 * t84;
t162 = 0.2e1 * qJD(2);
t161 = 2 * qJD(6);
t160 = pkin(5) * t36;
t159 = pkin(9) * t36;
t158 = pkin(9) * t58;
t154 = t35 * t84;
t153 = t35 * t86;
t76 = pkin(3) * t155 + pkin(9);
t152 = t36 * t76;
t151 = t36 * t84;
t149 = t57 * t84;
t148 = t57 * t86;
t146 = t58 * t76;
t145 = t84 * t86;
t19 = t39 * qJD(4) + t157 * t169;
t38 = -t155 * t60 + t156 * t61;
t80 = qJD(5) * t86;
t144 = t19 * t84 + t38 * t80;
t109 = pkin(3) * t114;
t134 = t84 * qJD(6);
t46 = pkin(5) * t136 - qJ(6) * t80 - t134;
t40 = t109 + t46;
t142 = -t40 - t46;
t122 = t156 * pkin(3);
t77 = -t122 - pkin(4);
t141 = t109 * t84 + t77 * t80;
t139 = pkin(3) * qJD(4);
t138 = qJ(6) * t36;
t102 = pkin(5) * t84 - qJ(6) * t86;
t20 = -t102 * t57 + t38;
t137 = qJD(5) * t20;
t135 = qJD(6) * t58;
t133 = t85 * qJD(3);
t131 = qJ(2) * qJD(3);
t130 = -0.2e1 * t151;
t127 = 0.2e1 * t149 * t35 - t55 * t80;
t126 = pkin(4) * t136;
t125 = pkin(4) * t80;
t124 = pkin(9) * t136;
t123 = pkin(9) * t80;
t121 = t58 * t80;
t120 = t84 * t80;
t119 = t84 * t156;
t118 = t86 * t156;
t16 = t165 * t36;
t117 = 0.4e1 * t57 * t145;
t110 = pkin(3) * t112;
t8 = qJ(6) * t58 + t167;
t101 = t30 * t86 - t39 * t84;
t9 = -pkin(5) * t58 - t101;
t105 = t8 * t86 + t84 * t9;
t104 = -t8 * t84 + t86 * t9;
t103 = -pkin(5) * t86 - qJ(6) * t84;
t64 = -pkin(4) + t103;
t52 = -t122 + t64;
t98 = -t35 * t52 + t57 * t40;
t97 = -t35 * t64 + t46 * t57;
t96 = t57 * t77 + t146;
t24 = t57 * t80 - t154;
t95 = t136 * t57 + t153;
t23 = t121 + t151;
t21 = t136 * t58 - t36 * t86;
t4 = -t13 * t84 + t136 * t39 + t18 * t86 - t30 * t80;
t94 = -t109 * t86 + t136 * t77;
t93 = -t97 - t159;
t45 = qJD(5) * t103 + qJD(6) * t86;
t6 = t102 * t35 + t45 * t57 + t19;
t92 = -t6 + (-t57 * t64 - t158) * qJD(5);
t91 = -t6 + (-t52 * t57 - t146) * qJD(5);
t44 = t165 * t156 * t139;
t2 = t135 - t4 + t138;
t3 = -t160 - t5;
t1 = qJD(5) * t104 + t2 * t86 + t3 * t84;
t90 = -t110 * t58 - t152 - t98;
t89 = t35 * t77 - t152 + (-t155 * t57 - t156 * t58) * t139;
t66 = 0.2e1 * t120;
t56 = -0.2e1 * t111;
t51 = t64 * t136;
t43 = t52 * t136;
t42 = t110 * t84 + t76 * t80;
t41 = -t110 * t86 + t136 * t76;
t31 = t38 * t136;
t17 = t20 * t136;
t12 = t111 * t57 + t145 * t35;
t7 = qJD(5) * t117 - t140 * t35;
t10 = [0, 0, 0, 0, t162, qJ(2) * t162, -0.2e1 * t85 * t132, 0.2e1 * (t85 ^ 2 - t87 ^ 2) * qJD(3), 0, 0, 0, 0.2e1 * qJD(2) * t85 + 0.2e1 * t131 * t87, 0.2e1 * qJD(2) * t87 - 0.2e1 * t131 * t85, -0.2e1 * t150, -0.2e1 * t35 * t58 + 0.2e1 * t36 * t57, 0, 0, 0, 0.2e1 * t36 * t74 + 0.2e1 * t58 * t67, 0.2e1 * t35 * t74 - 0.2e1 * t57 * t67, -0.2e1 * t120 * t55 - 0.2e1 * t150 * t83, 0.2e1 * t111 * t55 + t117 * t35, 0.2e1 * t153 * t58 + 0.2e1 * t21 * t57, -0.2e1 * t154 * t58 + 0.2e1 * t23 * t57, 0.2e1 * t147, 0.2e1 * t101 * t36 - 0.2e1 * t149 * t19 - 0.2e1 * t24 * t38 + 0.2e1 * t5 * t58, -0.2e1 * t148 * t19 - 0.2e1 * t167 * t36 + 0.2e1 * t38 * t95 + 0.2e1 * t4 * t58, -0.2e1 * t149 * t6 - 0.2e1 * t20 * t24 - 0.2e1 * t3 * t58 - 0.2e1 * t36 * t9, 0.2e1 * t104 * t35 + 0.2e1 * (qJD(5) * t105 + t2 * t84 - t3 * t86) * t57, 0.2e1 * t148 * t6 + 0.2e1 * t2 * t58 - 0.2e1 * t20 * t95 + 0.2e1 * t36 * t8, 0.2e1 * t2 * t8 + 0.2e1 * t20 * t6 + 0.2e1 * t3 * t9; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t130 * t58 - t54 * t80 + t127, t168 (-t121 + t130) * t58 + t127, 0, -t168, t1 * t58 + t105 * t36 - t20 * t35 + t57 * t6; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0.2e1 * t16 * t58 - 0.2e1 * t150; 0, 0, 0, 0, 0, 0, 0, 0, -t133, -t132, 0, -t88 * t133, -t88 * t132, 0, 0, t35, -t36, 0, -t19, t18, t12, t7, t23, -t21, 0, t31 + (-qJD(5) * t96 - t19) * t86 + t89 * t84, t136 * t96 + t86 * t89 + t144, t84 * t90 + t86 * t91 + t17, t1, t91 * t84 + (-t90 - t137) * t86, t20 * t40 + t6 * t52 + (t118 * t8 + t119 * t9) * t139 + t1 * t76; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t133, -t132, 0, 0, 0, 0, 0, t35, -t36, 0, 0, 0, 0, 0, t95, t24, t95, t16, -t24, t16 * t76 + t44 * t58 + t98; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -0.2e1 * t109, -0.2e1 * t110, t66, t56, 0, 0, 0, 0.2e1 * t94, 0.2e1 * t141, -0.2e1 * t40 * t86 + 0.2e1 * t43, 0.2e1 * t44, -0.2e1 * t40 * t84 - 0.2e1 * t52 * t80, 0.2e1 * t52 * t40 + 0.2e1 * t44 * t76; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t35, -t36, 0, -t19, t18, t12, t7, t23, -t21, 0, t31 + (-pkin(4) * t35 - t159) * t84 + (-t19 + (pkin(4) * t57 - t158) * qJD(5)) * t86, -pkin(4) * t95 + pkin(9) * t21 + t144, t84 * t93 + t86 * t92 + t17, t1, t92 * t84 + (-t93 - t137) * t86, pkin(9) * t1 + t20 * t46 + t6 * t64; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t35, -t36, 0, 0, 0, 0, 0, t95, t24, t95, t16, -t24, pkin(9) * t16 + t97; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t109, -t110, t66, t56, 0, 0, 0, t94 - t126, -t125 + t141, t142 * t86 + t43 + t51, t44, t142 * t84 + (-t52 - t64) * t80, pkin(9) * t44 + t40 * t64 + t52 * t46; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t66, t56, 0, 0, 0, -0.2e1 * t126, -0.2e1 * t125, -0.2e1 * t46 * t86 + 0.2e1 * t51, 0, -0.2e1 * t46 * t84 - 0.2e1 * t64 * t80, 0.2e1 * t64 * t46; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t95, t24, t36, t5, t4, t5 + 0.2e1 * t160, t103 * t35 + (-qJD(5) * t102 + t134) * t57, 0.2e1 * t135 - t4 + 0.2e1 * t138, -pkin(5) * t3 + qJ(6) * t2 + qJD(6) * t8; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t23, t21, -t23, 0, -t21, -t102 * t36 + t45 * t58; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t80, -t136, 0, -t42, t41, -t42, t45, -t41 (-pkin(5) * t119 + qJ(6) * t118) * t139 + t45 * t76; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t80, -t136, 0, -t123, t124, -t123, t45, -t124, t45 * pkin(9); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t161, qJ(6) * t161; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t36, t95, 0, t3; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t23; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t80, 0, t42; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t80, 0, t123; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0;];
MMD_reg  = t10;
