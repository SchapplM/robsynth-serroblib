% Calculate minimal parameter regressor of coriolis joint torque vector for
% S5RPRRP9
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
% Datum: 2019-12-31 18:50
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauc_reg = S5RPRRP9_coriolisvecJ_fixb_regmin_slag_vp(qJ, qJD, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(8,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RPRRP9_coriolisvecJ_fixb_regmin_slag_vp: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RPRRP9_coriolisvecJ_fixb_regmin_slag_vp: qJD has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S5RPRRP9_coriolisvecJ_fixb_regmin_slag_vp: pkin has to be [8x1] (double)');

%% Symbolic Calculation
% From coriolisvec_joint_fixb_regressor_minpar_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 18:49:11
% EndTime: 2019-12-31 18:49:15
% DurationCPUTime: 0.98s
% Computational Cost: add. (2255->188), mult. (6090->244), div. (0->0), fcn. (4635->6), ass. (0->112)
t140 = cos(qJ(4));
t87 = cos(pkin(8));
t90 = cos(qJ(3));
t128 = t90 * t87;
t86 = sin(pkin(8));
t89 = sin(qJ(3));
t130 = t89 * t86;
t105 = -t128 + t130;
t64 = t105 * qJD(1);
t68 = t90 * t86 + t89 * t87;
t65 = t68 * qJD(1);
t88 = sin(qJ(4));
t42 = t140 * t64 + t88 * t65;
t85 = qJD(3) + qJD(4);
t135 = t42 * t85;
t110 = qJD(4) * t140;
t121 = qJD(4) * t88;
t117 = qJD(1) * qJD(3);
t112 = t90 * t117;
t113 = t89 * t117;
t124 = t86 * t112 + t87 * t113;
t76 = t87 * t112;
t60 = -t86 * t113 + t76;
t19 = t64 * t110 + t65 * t121 + t88 * t124 - t140 * t60;
t7 = -t19 + t135;
t136 = t42 ^ 2;
t100 = t140 * t65 - t88 * t64;
t149 = t100 ^ 2;
t153 = -t136 + t149;
t80 = -t87 * pkin(2) - pkin(1);
t71 = t80 * qJD(1) + qJD(2);
t49 = t64 * pkin(3) + t71;
t12 = t42 * pkin(4) - qJ(5) * t100 + t49;
t152 = t12 * t42;
t151 = t49 * t42;
t134 = t100 * t42;
t137 = t100 * t85;
t20 = t100 * qJD(4) + t140 * t124 + t88 * t60;
t150 = -t20 + t137;
t21 = pkin(4) * t100 + t42 * qJ(5);
t127 = pkin(6) + qJ(2);
t72 = t127 * t86;
t69 = qJD(1) * t72;
t73 = t127 * t87;
t70 = qJD(1) * t73;
t147 = -t90 * t69 - t89 * t70;
t146 = -t68 * pkin(7) - t89 * t73;
t28 = -t124 * pkin(7) - qJD(2) * t64 + t147 * qJD(3);
t107 = t89 * t69 - t90 * t70;
t99 = t68 * qJD(2);
t96 = qJD(1) * t99;
t29 = -t60 * pkin(7) + t107 * qJD(3) - t96;
t37 = -t65 * pkin(7) + t147;
t36 = qJD(3) * pkin(3) + t37;
t38 = -t64 * pkin(7) - t107;
t2 = t38 * t110 + t36 * t121 - t140 * t29 + t88 * t28;
t102 = t100 * t12 + t2;
t145 = -t49 * t100 - t2;
t129 = t90 * t72;
t39 = -t129 + t146;
t106 = t89 * t72 - t90 * t73;
t40 = -t105 * pkin(7) - t106;
t101 = t140 * t39 - t88 * t40;
t122 = qJD(2) * t86;
t125 = qJD(2) * t128 - qJD(3) * t129;
t30 = t146 * qJD(3) - t89 * t122 + t125;
t66 = t105 * qJD(3);
t92 = t106 * qJD(3) - t99;
t31 = t66 * pkin(7) + t92;
t4 = t101 * qJD(4) + t140 * t30 + t88 * t31;
t143 = t4 * t85;
t18 = t140 * t40 + t88 * t39;
t5 = t18 * qJD(4) - t140 * t31 + t88 * t30;
t142 = t5 * t85;
t115 = t140 * t38;
t11 = t88 * t36 + t115;
t139 = t11 * t85;
t133 = t88 * t38;
t14 = t140 * t37 - t133;
t126 = -pkin(3) * t110 - qJD(5) + t14;
t123 = t86 ^ 2 + t87 ^ 2;
t120 = t65 * qJD(3);
t10 = t140 * t36 - t133;
t119 = qJD(5) - t10;
t118 = qJD(1) * qJD(2);
t116 = qJD(1) * t130;
t114 = t123 * qJD(1) ^ 2;
t109 = -t36 * t110 + t38 * t121 - t140 * t28 - t88 * t29;
t13 = t88 * t37 + t115;
t108 = pkin(3) * t121 - t13;
t82 = t85 * qJD(5);
t1 = t82 - t109;
t104 = 0.2e1 * t123 * t118;
t103 = t10 * t85 + t109;
t98 = t68 * qJD(3);
t97 = t140 * t105;
t95 = t19 + t135;
t53 = t105 * pkin(3) + t80;
t48 = -t88 * t105 + t140 * t68;
t3 = t124 * pkin(3) + t20 * pkin(4) + t19 * qJ(5) - qJD(5) * t100;
t93 = t20 + t137;
t81 = -t140 * pkin(3) - pkin(4);
t79 = t88 * pkin(3) + qJ(5);
t47 = t88 * t68 + t97;
t27 = t48 * qJD(4) + t140 * t98 - t88 * t66;
t26 = qJD(4) * t97 + t68 * t121 + t140 * t66 + t88 * t98;
t16 = t47 * pkin(4) - t48 * qJ(5) + t53;
t15 = t65 * pkin(3) + t21;
t9 = t85 * qJ(5) + t11;
t8 = -t85 * pkin(4) + t119;
t6 = pkin(3) * t98 + t27 * pkin(4) + t26 * qJ(5) - t48 * qJD(5);
t17 = [0, 0, 0, 0, 0, t104, qJ(2) * t104, t60 * t68 - t65 * t66, -t60 * t105 - t68 * t124 + t66 * t64 - t65 * t98, -t66 * qJD(3), -t68 * qJD(3) ^ 2, 0, t80 * t124 + (t71 * t68 + t92) * qJD(3), t80 * t60 - t71 * t66 - ((-qJD(3) * t73 - t122) * t89 + t125) * qJD(3), -t100 * t26 - t19 * t48, -t100 * t27 + t19 * t47 - t48 * t20 + t26 * t42, -t26 * t85, -t27 * t85, 0, t53 * t20 + t49 * t27 - t142 + (t124 * t47 + t42 * t98) * pkin(3), -t53 * t19 - t49 * t26 - t143 + (t100 * t98 + t124 * t48) * pkin(3), t12 * t27 + t16 * t20 + t3 * t47 + t6 * t42 - t142, -t1 * t47 + t100 * t5 + t101 * t19 - t18 * t20 + t2 * t48 - t8 * t26 - t9 * t27 - t4 * t42, -t100 * t6 + t12 * t26 + t16 * t19 - t3 * t48 + t143, t1 * t18 - t101 * t2 + t12 * t6 + t3 * t16 + t9 * t4 + t8 * t5; 0, 0, 0, 0, 0, -t114, -qJ(2) * t114, 0, 0, 0, 0, 0, t120 + t124, t76 + (-t64 - t116) * qJD(3), 0, 0, 0, 0, 0, t93, -t95, t93, -t136 - t149, t95, -t100 * t8 + t42 * t9 + t3; 0, 0, 0, 0, 0, 0, 0, t65 * t64, -t64 ^ 2 + t65 ^ 2, t76 + (t64 - t116) * qJD(3), t120 - t124, 0, -t71 * t65 - t96, t105 * t118 + t71 * t64, t134, t153, t7, t150, 0, t13 * t85 + (-t85 * t121 - t42 * t65) * pkin(3) + t145, t14 * t85 + t151 + (-t100 * t65 - t85 * t110) * pkin(3) + t109, -t108 * t85 - t15 * t42 - t102, -t81 * t19 - t79 * t20 + (t108 + t9) * t100 + (t126 + t8) * t42, t100 * t15 - t126 * t85 + t1 - t152, t1 * t79 + t108 * t8 - t12 * t15 - t126 * t9 + t2 * t81; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t134, t153, t7, t150, 0, t139 + t145, t103 + t151, -t21 * t42 - t102 + t139, pkin(4) * t19 - t20 * qJ(5) + (-t11 + t9) * t100 + (t8 - t119) * t42, t100 * t21 - t103 - t152 + 0.2e1 * t82, -t2 * pkin(4) + t1 * qJ(5) - t8 * t11 + t119 * t9 - t12 * t21; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t134, t7, -t85 ^ 2 - t149, -t9 * t85 + t102;];
tauc_reg = t17;
