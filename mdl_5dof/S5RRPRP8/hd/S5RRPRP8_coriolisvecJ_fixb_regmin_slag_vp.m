% Calculate minimal parameter regressor of coriolis joint torque vector for
% S5RRPRP8
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% qJD [5x1]
%   Generalized joint velocities
% pkin [7x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d1,d2,d4]';
% 
% Output:
% tauc_reg [5x23]
%   minimal parameter regressor of coriolis joint torque vector

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-31 20:04
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauc_reg = S5RRPRP8_coriolisvecJ_fixb_regmin_slag_vp(qJ, qJD, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(7,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RRPRP8_coriolisvecJ_fixb_regmin_slag_vp: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RRPRP8_coriolisvecJ_fixb_regmin_slag_vp: qJD has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [7 1]), ...
  'S5RRPRP8_coriolisvecJ_fixb_regmin_slag_vp: pkin has to be [7x1] (double)');

%% Symbolic Calculation
% From coriolisvec_joint_fixb_regressor_minpar_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 20:04:27
% EndTime: 2019-12-31 20:04:30
% DurationCPUTime: 0.93s
% Computational Cost: add. (1014->173), mult. (2415->240), div. (0->0), fcn. (1427->4), ass. (0->114)
t91 = sin(qJ(2));
t130 = qJD(1) * t91;
t76 = pkin(6) * t130;
t160 = -pkin(7) * t130 + qJD(3) + t76;
t92 = cos(qJ(4));
t125 = qJD(4) * t92;
t90 = sin(qJ(4));
t126 = qJD(4) * t90;
t93 = cos(qJ(2));
t127 = qJD(2) * t93;
t159 = t91 * t125 - t93 * t126 + t90 * t127;
t48 = t91 * t90 + t93 * t92;
t103 = t48 * qJD(4);
t123 = qJD(1) * qJD(2);
t118 = t93 * t123;
t119 = t91 * t123;
t13 = qJD(1) * t103 - t92 * t118 - t90 * t119;
t33 = t48 * qJD(1);
t84 = qJD(2) - qJD(4);
t155 = t33 * t84 + t13;
t129 = qJD(1) * t93;
t124 = qJ(3) * qJD(1);
t38 = -qJD(1) * pkin(1) - pkin(2) * t129 - t91 * t124;
t23 = pkin(3) * t129 - t38;
t35 = -t90 * t129 + t92 * t130;
t94 = -pkin(2) - pkin(3);
t120 = t94 * qJD(2);
t24 = t120 + t160;
t77 = pkin(6) * t129;
t54 = -pkin(7) * t129 + t77;
t86 = qJD(2) * qJ(3);
t37 = t54 + t86;
t109 = t90 * t24 + t92 * t37;
t128 = qJD(2) * t91;
t150 = pkin(6) - pkin(7);
t53 = t150 * t128;
t85 = qJD(2) * qJD(3);
t27 = -qJD(1) * t53 + t85;
t71 = pkin(6) * t118;
t45 = -pkin(7) * t118 + t71;
t99 = -t109 * qJD(4) - t90 * t27 + t92 * t45;
t158 = -t23 * t35 + t99;
t157 = -0.2e1 * t123;
t57 = t92 * qJ(3) + t90 * t94;
t156 = t57 * qJD(4) + t160 * t90 + t92 * t54;
t14 = t159 * qJD(1) - t92 * t119;
t154 = t35 * t84 + t14;
t107 = -t90 * qJ(3) + t92 * t94;
t153 = t107 * qJD(4) + t160 * t92 - t90 * t54;
t151 = t35 ^ 2;
t32 = t33 ^ 2;
t152 = -t32 + t151;
t117 = t92 * t24 - t90 * t37;
t131 = t35 * qJ(5);
t6 = t117 - t131;
t5 = -t84 * pkin(4) + t6;
t149 = t5 - t6;
t148 = -t131 + t153;
t132 = t33 * qJ(5);
t147 = t132 - t156;
t144 = t35 * t33;
t96 = qJD(1) ^ 2;
t142 = t93 * t96;
t95 = qJD(2) ^ 2;
t141 = t95 * t91;
t140 = t95 * t93;
t80 = t91 * qJD(3);
t136 = qJ(3) * t118 + qJD(1) * t80;
t135 = t93 * t86 + t80;
t87 = t91 ^ 2;
t134 = -t93 ^ 2 + t87;
t133 = qJD(2) * pkin(2);
t59 = -t93 * pkin(2) - t91 * qJ(3) - pkin(1);
t122 = t91 * t142;
t62 = t150 * t93;
t46 = t93 * pkin(3) - t59;
t114 = qJD(1) * t59 + t38;
t113 = pkin(1) * t157;
t112 = qJD(3) - t133;
t111 = t84 ^ 2;
t110 = t91 * t120;
t61 = t150 * t91;
t108 = -t90 * t61 - t92 * t62;
t74 = t93 * t124;
t30 = t94 * t130 + t74;
t22 = pkin(2) * t119 - t136;
t31 = pkin(2) * t128 - t135;
t106 = -pkin(6) * t95 - qJD(1) * t31 - t22;
t105 = -t24 * t125 + t37 * t126 - t92 * t27 - t90 * t45;
t55 = qJD(2) * t62;
t104 = t61 * t125 - t62 * t126 - t92 * t53 + t90 * t55;
t19 = t110 + t135;
t17 = qJD(1) * t110 + t136;
t101 = t23 * t33 + t105;
t100 = t14 * pkin(4) + t17;
t98 = t108 * qJD(4) + t90 * t53 + t92 * t55;
t56 = -pkin(6) * t119 + t85;
t58 = t112 + t76;
t60 = t77 + t86;
t97 = t56 * t93 + (t58 * t93 + (-t60 + t77) * t91) * qJD(2);
t51 = -pkin(4) + t107;
t50 = pkin(2) * t130 - t74;
t49 = -t93 * t90 + t91 * t92;
t16 = t48 * qJD(2) - t103;
t15 = -t92 * t128 + t159;
t12 = t33 * pkin(4) + qJD(5) + t23;
t11 = -t48 * qJ(5) - t108;
t10 = -t49 * qJ(5) + t92 * t61 - t90 * t62;
t7 = t109 - t132;
t4 = -t16 * qJ(5) - t49 * qJD(5) + t98;
t3 = -t15 * qJ(5) - t48 * qJD(5) + t104;
t2 = t13 * qJ(5) - t35 * qJD(5) + t99;
t1 = -t14 * qJ(5) - t33 * qJD(5) - t105;
t8 = [0, 0, 0, 0.2e1 * t91 * t118, t134 * t157, t140, -t141, 0, -pkin(6) * t140 + t91 * t113, pkin(6) * t141 + t93 * t113, t106 * t93 + t114 * t128, t97, t106 * t91 - t114 * t127, t97 * pkin(6) + t22 * t59 + t38 * t31, -t13 * t49 + t35 * t16, t13 * t48 - t49 * t14 - t35 * t15 - t16 * t33, -t16 * t84, t15 * t84, 0, t46 * t14 + t23 * t15 + t17 * t48 + t19 * t33 - t98 * t84, t104 * t84 - t46 * t13 + t23 * t16 + t17 * t49 + t19 * t35, -t1 * t48 + t10 * t13 - t11 * t14 - t7 * t15 - t5 * t16 - t2 * t49 - t3 * t33 - t4 * t35, t1 * t11 + t7 * t3 + t2 * t10 + t5 * t4 + t100 * (t48 * pkin(4) + t46) + t12 * (t15 * pkin(4) + t19); 0, 0, 0, -t122, t134 * t96, 0, 0, 0, t96 * pkin(1) * t91, pkin(1) * t142, (-t38 * t91 + t50 * t93) * qJD(1), ((t60 - t86) * t91 + (t112 - t58) * t93) * qJD(1), 0.2e1 * t85 + (t38 * t93 + t50 * t91) * qJD(1), t56 * qJ(3) + t60 * qJD(3) - t38 * t50 + (t60 * t91 + (-t58 - t133) * t93) * qJD(1) * pkin(6), -t144, -t152, t155, t154, 0, t156 * t84 - t30 * t33 - t158, t153 * t84 - t30 * t35 - t101, t51 * t13 - t57 * t14 + (-t7 - t147) * t35 + (t5 - t148) * t33, t1 * t57 + t2 * t51 - t12 * (-t35 * pkin(4) + t30) + t148 * t7 + t147 * t5; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t122, 0, -t87 * t96 - t95, -t60 * qJD(2) + t38 * t130 + t71, 0, 0, 0, 0, 0, -t90 * t111 - t33 * t130, -t92 * t111 - t35 * t130, -t154 * t90 + t155 * t92, -t12 * t130 + (-t84 * t7 + t2) * t92 + (t84 * t5 + t1) * t90; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t144, t152, -t155, -t154, 0, -t109 * t84 + t158, -t117 * t84 + t101, pkin(4) * t13 - t149 * t33, t149 * t7 + (-t12 * t35 + t2) * pkin(4); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t32 - t151, t7 * t33 + t5 * t35 + t100;];
tauc_reg = t8;
