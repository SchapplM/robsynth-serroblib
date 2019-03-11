% Calculate minimal parameter regressor of joint inertia matrix time derivative for
% S6RRPRPP4
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% qJD [6x1]
%   Generalized joint velocities
% pkin [9x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d2,d4,theta5]';
% 
% Output:
% MMD_reg [((6+1)*6/2)x27]
%   minimal parameter regressor of inerta matrix time derivative
%   (only lower left triangular matrix (including diagonal) due to symmetry

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-09 10:02
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function MMD_reg = S6RRPRPP4_inertiaDJ_regmin_slag_vp(qJ, qJD, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(9,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRPRPP4_inertiaDJ_regmin_slag_vp: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RRPRPP4_inertiaDJ_regmin_slag_vp: qJD has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S6RRPRPP4_inertiaDJ_regmin_slag_vp: pkin has to be [9x1] (double)');

%% Symbolic Calculation
% From inertiaD_joint_joint_fixb_regressor_minpar_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 10:01:40
% EndTime: 2019-03-09 10:01:44
% DurationCPUTime: 1.43s
% Computational Cost: add. (2209->202), mult. (4699->346), div. (0->0), fcn. (3760->6), ass. (0->117)
t101 = cos(qJ(4));
t158 = sin(pkin(9));
t128 = qJD(4) * t158;
t159 = cos(pkin(9));
t129 = qJD(4) * t159;
t99 = sin(qJ(4));
t50 = -t101 * t129 + t99 * t128;
t51 = -t101 * t128 - t99 * t129;
t60 = t159 * t101 - t158 * t99;
t61 = -t158 * t101 - t159 * t99;
t173 = 0.2e1 * t50 * t61 + 0.2e1 * t60 * t51;
t166 = pkin(3) + pkin(7);
t172 = (-t158 * t50 + t159 * t51) * pkin(4);
t102 = cos(qJ(2));
t103 = -pkin(2) - pkin(8);
t100 = sin(qJ(2));
t154 = t100 * qJ(3);
t118 = -t102 * t103 + t154;
t58 = -pkin(1) - t118;
t72 = t166 * t100;
t162 = t101 * t58 + t99 * t72;
t98 = t102 ^ 2;
t133 = qJD(2) * (t100 ^ 2 - t98);
t95 = t99 ^ 2;
t161 = -t101 ^ 2 + t95;
t132 = t161 * qJD(4);
t145 = t101 * qJD(5);
t151 = qJ(5) - t103;
t155 = qJD(4) * t99;
t109 = t151 * t155 - t145;
t127 = t151 * t101;
t45 = -qJD(4) * t127 - t99 * qJD(5);
t29 = -t159 * t109 + t158 * t45;
t30 = t158 * t109 + t159 * t45;
t66 = t151 * t99;
t38 = t159 * t127 - t158 * t66;
t39 = -t158 * t127 - t159 * t66;
t170 = t29 * t60 + t30 * t61 + t38 * t51 + t39 * t50;
t144 = t102 * qJD(3);
t73 = t166 * t102;
t152 = t73 * qJD(4);
t169 = qJD(2) * t118 - t144 - t152;
t168 = 0.2e1 * qJD(3);
t167 = 2 * qJD(6);
t147 = t100 * qJD(2);
t136 = t101 * t147;
t149 = qJD(4) * t102;
t140 = t99 * t149;
t111 = t136 + t140;
t150 = qJD(4) * t101;
t126 = pkin(2) * t147 - t100 * qJD(3);
t157 = qJ(3) * t102;
t43 = (pkin(8) * t100 - t157) * qJD(2) + t126;
t91 = t102 * qJD(2);
t88 = pkin(7) * t91;
t65 = pkin(3) * t91 + t88;
t17 = -t101 * t43 - t72 * t150 + t58 * t155 - t99 * t65;
t11 = t111 * qJ(5) - t102 * t145 - t17;
t130 = qJ(5) * t102 - t58;
t56 = t101 * t65;
t9 = pkin(4) * t91 + t56 + t130 * t150 + (-qJ(5) * t147 - qJD(4) * t72 + qJD(5) * t102 - t43) * t99;
t4 = t159 * t11 + t158 * t9;
t63 = t101 * t72;
t31 = t100 * pkin(4) + t130 * t99 + t63;
t153 = t101 * t102;
t35 = -qJ(5) * t153 + t162;
t15 = t158 * t31 + t159 * t35;
t86 = t99 * pkin(4) + qJ(3);
t156 = qJD(2) * t73;
t76 = pkin(4) * t150 + qJD(3);
t148 = qJD(4) * t103;
t146 = t100 * qJD(6);
t143 = -0.2e1 * pkin(1) * qJD(2);
t142 = qJ(6) * t91 + t4;
t52 = pkin(4) * t153 + t73;
t141 = pkin(7) * t147;
t139 = t99 * t150;
t138 = t101 * t149;
t137 = t100 * t91;
t135 = t38 * t29 + t39 * t30;
t125 = t99 * t136;
t123 = t159 * t147;
t122 = t158 * t147;
t3 = -t158 * t11 + t159 * t9;
t119 = -t102 * pkin(2) - t154;
t47 = t61 * t102;
t32 = -qJD(4) * t47 + t101 * t123 - t99 * t122;
t33 = -t101 * t122 - t50 * t102 - t99 * t123;
t46 = t60 * t102;
t115 = t29 * t47 - t30 * t46 + t39 * t32 - t38 * t33;
t113 = -t61 * t32 + t60 * t33 + t50 * t46 - t51 * t47;
t81 = t158 * pkin(4) + qJ(6);
t85 = -t159 * pkin(4) - pkin(5);
t112 = qJD(6) * t61 + t81 * t50 + t85 * t51;
t14 = -t158 * t35 + t159 * t31;
t64 = t166 * t147;
t108 = -t64 + (-t100 * t103 - t157) * qJD(4);
t107 = 0.2e1 * t170;
t1 = t142 + t146;
t12 = t100 * qJ(6) + t15;
t13 = -t100 * pkin(5) - t14;
t2 = -pkin(5) * t91 - t3;
t106 = -t1 * t61 - t12 * t50 - t13 * t51 - t2 * t60;
t105 = t14 * t51 - t15 * t50 + t3 * t60 - t4 * t61;
t104 = t119 * qJD(2) + t144;
t41 = -pkin(4) * t140 + (-pkin(4) * t101 - t166) * t147;
t75 = 0.2e1 * t137;
t67 = -pkin(1) + t119;
t54 = -t100 * t150 - t99 * t91;
t53 = -t100 * t155 + t101 * t91;
t49 = -qJ(3) * t91 + t126;
t36 = -t61 * pkin(5) - t60 * qJ(6) + t86;
t24 = t46 * pkin(5) - t47 * qJ(6) + t52;
t23 = -t50 * pkin(5) - t51 * qJ(6) - t60 * qJD(6) + t76;
t18 = -t162 * qJD(4) - t99 * t43 + t56;
t5 = -t32 * pkin(5) + t33 * qJ(6) - t47 * qJD(6) + t41;
t6 = [0, 0, 0, t75, -0.2e1 * t133, 0, 0, 0, t100 * t143, t102 * t143, 0, 0.2e1 * t49 * t102 - 0.2e1 * t67 * t147, -0.2e1 * t49 * t100 - 0.2e1 * t67 * t91, 0.2e1 * t67 * t49, -0.2e1 * t95 * t137 + 0.2e1 * t98 * t139, -0.4e1 * t102 * t125 - 0.2e1 * t98 * t132, -0.2e1 * t100 * t138 + 0.2e1 * t99 * t133, 0.2e1 * t100 * t140 + 0.2e1 * t101 * t133, t75, 0.2e1 * (-t101 * t156 + t18) * t100 + 0.2e1 * ((-t99 * t58 + t63) * qJD(2) - t64 * t101 - t99 * t152) * t102, 0.2e1 * (t99 * t156 + t17) * t100 + 0.2e1 * (-t162 * qJD(2) - t73 * t150 + t64 * t99) * t102, 0.2e1 * t14 * t33 + 0.2e1 * t15 * t32 - 0.2e1 * t3 * t47 - 0.2e1 * t4 * t46, 0.2e1 * t14 * t3 + 0.2e1 * t15 * t4 + 0.2e1 * t52 * t41, -0.2e1 * t2 * t100 - 0.2e1 * t13 * t91 - 0.2e1 * t24 * t32 + 0.2e1 * t5 * t46, -0.2e1 * t1 * t46 + 0.2e1 * t12 * t32 - 0.2e1 * t13 * t33 + 0.2e1 * t2 * t47, 0.2e1 * t1 * t100 + 0.2e1 * t12 * t91 + 0.2e1 * t24 * t33 - 0.2e1 * t5 * t47, 0.2e1 * t12 * t1 + 0.2e1 * t13 * t2 + 0.2e1 * t24 * t5; 0, 0, 0, 0, 0, t91, -t147, 0, -t88, t141, t104, t88, -t141, t104 * pkin(7), t102 * t132 + t125, 0.4e1 * t99 * t138 - t161 * t147, t53, t54, 0, -t169 * t101 + t108 * t99, t108 * t101 + t169 * t99, -t105 + t115, -t14 * t29 + t15 * t30 - t3 * t38 + t4 * t39 + t41 * t86 + t52 * t76, -t29 * t100 + t23 * t46 - t24 * t50 - t36 * t32 - t38 * t91 - t5 * t61, -t106 + t115, t30 * t100 - t23 * t47 - t24 * t51 + t36 * t33 + t39 * t91 - t5 * t60, t1 * t39 + t12 * t30 + t13 * t29 + t2 * t38 + t24 * t23 + t5 * t36; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t168, qJ(3) * t168, -0.2e1 * t139, 0.2e1 * t132, 0, 0, 0, 0.2e1 * qJ(3) * t150 + 0.2e1 * qJD(3) * t99, -0.2e1 * qJ(3) * t155 + 0.2e1 * qJD(3) * t101, t107, 0.2e1 * t86 * t76 + 0.2e1 * t135, -0.2e1 * t23 * t61 - 0.2e1 * t36 * t50, t107, -0.2e1 * t23 * t60 - 0.2e1 * t36 * t51, 0.2e1 * t36 * t23 + 0.2e1 * t135; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t91, 0, 0, t88, 0, 0, 0, 0, 0, t53, t54, t113, t105, t51 * t100 + t60 * t91, t113, -t50 * t100 - t61 * t91, t106; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t173, -t170, 0, -t173, 0, -t170; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t173, 0, 0, 0, t173; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t99 * t147 - t138, t111, t91, t18, t17 (t158 * t32 + t159 * t33) * pkin(4) (t158 * t4 + t159 * t3) * pkin(4) (pkin(5) - t85) * t91 + t3, -qJD(6) * t46 + t81 * t32 - t85 * t33, t81 * t91 + t142 + 0.2e1 * t146, t12 * qJD(6) + t1 * t81 + t2 * t85; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t155, -t150, 0, -t99 * t148, -t101 * t148, -t172 (t158 * t30 - t159 * t29) * pkin(4), -t29, t112, t30, t39 * qJD(6) + t29 * t85 + t30 * t81; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t155, -t150, 0, t172, t51, 0, -t50, -t112; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t167, t81 * t167; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t41, -t32, 0, t33, t5; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t76, -t50, 0, -t51, t23; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t91, -t33, 0, t2; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t51, 0, t29; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t51; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0;];
MMD_reg  = t6;
