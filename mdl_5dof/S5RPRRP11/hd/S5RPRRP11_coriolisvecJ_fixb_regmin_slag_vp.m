% Calculate minimal parameter regressor of coriolis joint torque vector for
% S5RPRRP11
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% qJD [5x1]
%   Generalized joint velocities
% pkin [8x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d1,d3,d4,theta2]';
% 
% Output:
% tauc_reg [5x25]
%   minimal parameter regressor of coriolis joint torque vector

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-31 18:55
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauc_reg = S5RPRRP11_coriolisvecJ_fixb_regmin_slag_vp(qJ, qJD, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(8,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RPRRP11_coriolisvecJ_fixb_regmin_slag_vp: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RPRRP11_coriolisvecJ_fixb_regmin_slag_vp: qJD has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S5RPRRP11_coriolisvecJ_fixb_regmin_slag_vp: pkin has to be [8x1] (double)');

%% Symbolic Calculation
% From coriolisvec_joint_fixb_regressor_minpar_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 18:54:32
% EndTime: 2019-12-31 18:54:38
% DurationCPUTime: 1.75s
% Computational Cost: add. (2506->247), mult. (6547->326), div. (0->0), fcn. (4752->6), ass. (0->111)
t92 = sin(pkin(8));
t93 = cos(pkin(8));
t95 = sin(qJ(3));
t97 = cos(qJ(3));
t164 = -t92 * t95 + t97 * t93;
t166 = t164 * qJD(1);
t136 = pkin(6) + qJ(2);
t81 = t136 * t92;
t77 = qJD(1) * t81;
t82 = t136 * t93;
t78 = qJD(1) * t82;
t49 = -t95 * t77 + t97 * t78;
t165 = qJD(3) * t49;
t96 = cos(qJ(4));
t127 = qJD(4) * t96;
t94 = sin(qJ(4));
t129 = qJD(3) * t94;
t76 = t92 * t97 + t93 * t95;
t158 = t76 * qJD(1);
t28 = t158 * t127 + (qJD(4) + t166) * t129;
t65 = qJD(4) - t166;
t163 = t158 * qJD(3);
t87 = -pkin(2) * t93 - pkin(1);
t79 = t87 * qJD(1) + qJD(2);
t30 = -pkin(3) * t166 - pkin(7) * t158 + t79;
t44 = qJD(3) * pkin(7) + t49;
t13 = t30 * t94 + t44 * t96;
t8 = qJ(5) * t65 + t13;
t151 = t65 * t8;
t128 = qJD(4) * t94;
t101 = t164 * qJD(2);
t160 = -t77 * t97 - t95 * t78;
t21 = qJD(1) * t101 + qJD(3) * t160;
t73 = t76 * qJD(3);
t62 = qJD(1) * t73;
t72 = t164 * qJD(3);
t99 = qJD(1) * t72;
t35 = t62 * pkin(3) - pkin(7) * t99;
t118 = t44 * t127 + t30 * t128 + t94 * t21 - t96 * t35;
t154 = pkin(4) * t62;
t2 = t118 - t154;
t161 = -t2 + t151;
t47 = -pkin(3) * t164 - pkin(7) * t76 + t87;
t53 = -t81 * t95 + t82 * t97;
t132 = t94 * t47 + t96 * t53;
t52 = t97 * t81 + t95 * t82;
t43 = -qJD(3) * pkin(3) - t160;
t126 = t96 * qJD(3);
t54 = t158 * t94 - t126;
t56 = t158 * t96 + t129;
t14 = pkin(4) * t54 - qJ(5) * t56 + t43;
t153 = pkin(7) * t62;
t157 = t65 * t14 - t153;
t156 = t56 ^ 2;
t155 = t65 ^ 2;
t124 = qJD(1) * qJD(2);
t22 = t76 * t124 + t165;
t27 = -qJD(4) * t126 + t128 * t158 - t96 * t99;
t4 = pkin(4) * t28 + qJ(5) * t27 - qJD(5) * t56 + t22;
t152 = t4 * t94;
t150 = t13 * t65;
t149 = t27 * t94;
t148 = t54 * t166;
t147 = t56 * t54;
t120 = t56 * t65;
t146 = t56 * t158;
t145 = t65 * t94;
t144 = t158 * t54;
t143 = t76 * t96;
t59 = t94 * t62;
t61 = t96 * t62;
t114 = pkin(4) * t94 - qJ(5) * t96;
t135 = t94 * qJD(5) - t65 * t114 + t49;
t134 = -t54 * t127 - t94 * t28;
t45 = pkin(3) * t158 - pkin(7) * t166;
t133 = t160 * t96 + t94 * t45;
t131 = t92 ^ 2 + t93 ^ 2;
t130 = qJ(5) * t62;
t12 = t30 * t96 - t44 * t94;
t125 = qJD(5) - t12;
t123 = pkin(7) * t128;
t121 = t131 * qJD(1) ^ 2;
t7 = -pkin(4) * t65 + t125;
t117 = t7 * t96 - t8 * t94;
t106 = t30 * t127 - t44 * t128 + t96 * t21 + t94 * t35;
t1 = qJD(5) * t65 + t106 + t130;
t116 = t65 * t7 + t1;
t115 = pkin(4) * t96 + qJ(5) * t94;
t112 = t59 + (-t166 * t96 + t127) * t65;
t111 = 0.2e1 * t131 * t124;
t110 = -t65 * t128 + t145 * t166 + t61;
t109 = t76 * t127 + t72 * t94;
t108 = t76 * t128 - t72 * t96;
t107 = t14 * t56 + t118;
t31 = -t52 * qJD(3) + t101;
t46 = pkin(3) * t73 - pkin(7) * t72;
t105 = t47 * t127 - t53 * t128 + t96 * t31 + t94 * t46;
t104 = t65 * t43 - t153;
t32 = t76 * qJD(2) + t53 * qJD(3);
t80 = -pkin(3) - t115;
t23 = pkin(4) * t56 + qJ(5) * t54;
t18 = t114 * t76 + t52;
t16 = pkin(4) * t164 - t47 * t96 + t53 * t94;
t15 = -qJ(5) * t164 + t132;
t11 = t54 * t65 - t27;
t10 = -pkin(4) * t158 + t160 * t94 - t45 * t96;
t9 = qJ(5) * t158 + t133;
t6 = t114 * t72 + (t115 * qJD(4) - qJD(5) * t96) * t76 + t32;
t5 = -t73 * pkin(4) + t132 * qJD(4) + t94 * t31 - t96 * t46;
t3 = qJ(5) * t73 - qJD(5) * t164 + t105;
t17 = [0, 0, 0, 0, 0, t111, qJ(2) * t111, t158 * t72 + t76 * t99, -t158 * t73 + t164 * t99 + t166 * t72 - t76 * t62, t72 * qJD(3), -t73 * qJD(3), 0, -qJD(3) * t32 + t62 * t87 + t73 * t79, t79 * t72 + (t166 * t87 - t31) * qJD(3), -t108 * t56 - t27 * t143, (-t54 * t96 - t56 * t94) * t72 + (t149 - t28 * t96 + (t54 * t94 - t56 * t96) * qJD(4)) * t76, -t108 * t65 + t164 * t27 + t56 * t73 + t76 * t61, -t109 * t65 + t164 * t28 - t54 * t73 - t76 * t59, -t164 * t62 + t65 * t73, t118 * t164 + t12 * t73 + t32 * t54 + t52 * t28 + ((-qJD(4) * t53 + t46) * t65 + t47 * t62 + t43 * qJD(4) * t76) * t96 + ((-qJD(4) * t47 - t31) * t65 - t53 * t62 + t22 * t76 + t43 * t72) * t94, -t105 * t65 + t106 * t164 - t108 * t43 - t13 * t73 - t132 * t62 + t22 * t143 - t52 * t27 + t32 * t56, t109 * t14 + t76 * t152 - t16 * t62 + t164 * t2 + t18 * t28 - t5 * t65 + t54 * t6 - t7 * t73, -t15 * t28 - t16 * t27 - t3 * t54 + t5 * t56 + t117 * t72 + (-t1 * t94 + t2 * t96 + (-t7 * t94 - t8 * t96) * qJD(4)) * t76, -t1 * t164 + t108 * t14 - t4 * t143 + t15 * t62 + t18 * t27 + t3 * t65 - t56 * t6 + t73 * t8, t1 * t15 + t14 * t6 + t16 * t2 + t18 * t4 + t3 * t8 + t5 * t7; 0, 0, 0, 0, 0, -t121, -qJ(2) * t121, 0, 0, 0, 0, 0, 0.2e1 * t163, 0.2e1 * t166 * qJD(3), 0, 0, 0, 0, 0, t110 - t144, -t155 * t96 - t146 - t59, -t145 * t65 - t144 + t61, (t27 + t148) * t96 + t94 * t120 + t134, t112 + t146, t116 * t94 - t14 * t158 + t161 * t96; 0, 0, 0, 0, 0, 0, 0, -t158 * t166, t158 ^ 2 - t166 ^ 2, 0, 0, 0, -t158 * t79 + t165 - t22, (-qJD(2) - t79) * t166, t96 * t120 - t149, (-t27 + t148) * t96 - t56 * t145 + t134, t112 - t146, t110 + t144, -t65 * t158, -pkin(3) * t28 - t12 * t158 - t49 * t54 + (-t22 + (-pkin(7) * qJD(4) - t45) * t65) * t96 + (t160 * t65 + t104) * t94, pkin(3) * t27 + t13 * t158 + t22 * t94 - t49 * t56 + (t123 + t133) * t65 + t104 * t96, t28 * t80 - t4 * t96 + t7 * t158 + (-pkin(7) * t127 + t10) * t65 - t135 * t54 + t157 * t94, -t10 * t56 + t54 * t9 + ((qJD(4) * t56 - t28) * pkin(7) + t116) * t96 + ((qJD(4) * t54 - t27) * pkin(7) - t161) * t94, t27 * t80 - t152 - t158 * t8 + (-t9 - t123) * t65 + t135 * t56 - t157 * t96, -t7 * t10 + t4 * t80 - t8 * t9 - t135 * t14 + (t117 * qJD(4) + t1 * t96 + t2 * t94) * pkin(7); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t147, -t54 ^ 2 + t156, t11, t120 - t28, t62, -t43 * t56 - t118 + t150, t12 * t65 + t43 * t54 - t106, -t23 * t54 - t107 + t150 + 0.2e1 * t154, pkin(4) * t27 - t28 * qJ(5) + (-t13 + t8) * t56 + (t7 - t125) * t54, 0.2e1 * t130 - t14 * t54 + t23 * t56 + (0.2e1 * qJD(5) - t12) * t65 + t106, -t2 * pkin(4) + t1 * qJ(5) + t125 * t8 - t7 * t13 - t14 * t23; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t147 - t163, t11, -t155 - t156, t107 - t151 - t154;];
tauc_reg = t17;
