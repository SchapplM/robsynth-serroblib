% Calculate minimal parameter regressor of coriolis joint torque vector for
% S5RPRRP6
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
% tauc_reg [5x22]
%   minimal parameter regressor of coriolis joint torque vector

% Quelle: HybrDyn-Toolbox
% Datum: 2021-01-15 18:09
% Revision: 24b2e7d74a0c1a3b64fa2f8f5ad758691ad61af3 (2021-01-15)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauc_reg = S5RPRRP6_coriolisvecJ_fixb_regmin_slag_vp(qJ, qJD, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(8,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RPRRP6_coriolisvecJ_fixb_regmin_slag_vp: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RPRRP6_coriolisvecJ_fixb_regmin_slag_vp: qJD has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S5RPRRP6_coriolisvecJ_fixb_regmin_slag_vp: pkin has to be [8x1] (double)');

%% Symbolic Calculation
% From coriolisvec_joint_fixb_regressor_minpar_matlab.m
% OptimizationMode: 2
% StartTime: 2021-01-15 18:08:27
% EndTime: 2021-01-15 18:08:33
% DurationCPUTime: 1.11s
% Computational Cost: add. (1470->234), mult. (3454->324), div. (0->0), fcn. (2061->6), ass. (0->130)
t81 = sin(qJ(3));
t130 = qJD(1) * t81;
t80 = sin(qJ(4));
t111 = t80 * t130;
t82 = cos(qJ(4));
t122 = t82 * qJD(3);
t56 = t111 - t122;
t83 = cos(qJ(3));
t120 = t83 * qJD(1);
t68 = -qJD(4) + t120;
t152 = t56 * t68;
t110 = t83 * t122;
t117 = qJD(3) * qJD(4);
t24 = -qJD(1) * t110 + qJD(4) * t111 - t82 * t117;
t167 = -t24 + t152;
t123 = t80 * qJD(3);
t58 = t82 * t130 + t123;
t150 = t58 * t68;
t125 = qJD(4) * t82;
t112 = t81 * t125;
t88 = t83 * t123 + t112;
t25 = t88 * qJD(1) + t80 * t117;
t166 = t25 - t150;
t69 = sin(pkin(8)) * pkin(1) + pkin(6);
t63 = t69 * qJD(1);
t165 = t83 * qJD(2) - t81 * t63;
t126 = qJD(4) * t80;
t113 = t81 * t126;
t76 = t81 ^ 2;
t94 = qJD(1) * t76 - t68 * t83;
t164 = -t68 * t113 - t94 * t122;
t163 = t58 ^ 2;
t74 = t81 * qJD(2);
t36 = t83 * t63 + t74;
t29 = qJD(3) * pkin(7) + t36;
t70 = -cos(pkin(8)) * pkin(1) - pkin(2);
t50 = -t83 * pkin(3) - t81 * pkin(7) + t70;
t32 = t50 * qJD(1);
t11 = -t80 * t29 + t82 * t32;
t8 = -t58 * qJ(5) + t11;
t7 = -t68 * pkin(4) + t8;
t162 = t7 - t8;
t161 = pkin(4) * t80;
t160 = t56 * pkin(4);
t148 = t80 * t32;
t12 = t82 * t29 + t148;
t9 = -t56 * qJ(5) + t12;
t159 = t9 * t68;
t128 = qJD(3) * t83;
t31 = qJD(3) * t74 + t63 * t128;
t13 = t25 * pkin(4) + t31;
t158 = t13 * t80;
t157 = t13 * t82;
t28 = -qJD(3) * pkin(3) - t165;
t156 = t28 * t80;
t155 = t28 * t82;
t154 = t31 * t80;
t153 = t31 * t82;
t151 = t56 * t81;
t149 = t68 * t82;
t147 = t80 * t83;
t146 = t81 * t82;
t145 = t82 * t83;
t144 = t83 * t25;
t84 = qJD(3) ^ 2;
t143 = t84 * t81;
t142 = t84 * t83;
t141 = -qJ(5) - pkin(7);
t103 = qJD(4) * t141;
t97 = pkin(3) * t81 - pkin(7) * t83;
t61 = t97 * qJD(1);
t105 = -t165 * t80 + t82 * t61;
t131 = qJ(5) * t82;
t91 = pkin(4) * t81 - t83 * t131;
t140 = t91 * qJD(1) + t80 * qJD(5) - t82 * t103 + t105;
t116 = t80 * t120;
t121 = t82 * qJD(5);
t137 = t165 * t82 + t80 * t61;
t139 = -qJ(5) * t116 - t80 * t103 - t121 + t137;
t138 = -t56 * t110 - t25 * t146;
t62 = t97 * qJD(3);
t136 = t50 * t125 + t80 * t62;
t135 = t81 * t69 * t123 + t82 * t62;
t60 = t69 * t145;
t134 = t80 * t50 + t60;
t133 = -t83 ^ 2 + t76;
t132 = qJ(5) * t81;
t64 = qJD(1) * t70;
t129 = qJD(3) * t81;
t127 = qJD(4) * t56;
t124 = t28 * qJD(4);
t118 = qJD(1) * qJD(3);
t115 = t58 * t128;
t114 = t68 * t126;
t108 = t81 * t118;
t107 = t58 * t129 + t24 * t83;
t30 = t165 * qJD(3);
t49 = qJD(1) * t62;
t106 = t80 * t30 - t82 * t49;
t104 = -qJD(5) - t160;
t102 = -t32 * t125 + t29 * t126 - t82 * t30 - t80 * t49;
t100 = t68 * t112;
t99 = t58 * t112;
t98 = pkin(4) * t108;
t96 = -t7 * t82 - t80 * t9;
t95 = t7 * t80 - t82 * t9;
t92 = 0.2e1 * qJD(3) * t64;
t90 = t25 * qJ(5) + t102;
t89 = t94 * t80;
t87 = -t12 * qJD(4) - t106;
t86 = t24 * qJ(5) + t87;
t85 = qJD(1) ^ 2;
t73 = -t82 * pkin(4) - pkin(3);
t66 = t141 * t82;
t65 = t141 * t80;
t53 = t56 ^ 2;
t44 = (t69 + t161) * t81;
t40 = t82 * t50;
t20 = pkin(4) * t116 + t36;
t19 = t88 * pkin(4) + t69 * t128;
t17 = -t104 + t28;
t16 = -t80 * t132 + t134;
t15 = -t81 * t131 + t40 + (-t69 * t80 - pkin(4)) * t83;
t6 = t100 - t144 + (-t89 + t151) * qJD(3);
t5 = t107 + t164;
t4 = (-qJ(5) * qJD(4) - qJD(3) * t69) * t146 + (-qJD(5) * t81 + (-qJ(5) * qJD(3) - qJD(4) * t69) * t83) * t80 + t136;
t3 = -t81 * t121 + t91 * qJD(3) + (-t60 + (-t50 + t132) * t80) * qJD(4) + t135;
t2 = -t56 * qJD(5) - t90;
t1 = -t58 * qJD(5) + t86 + t98;
t10 = [0, 0, 0, 0, 0.2e1 * t83 * t108, -0.2e1 * t133 * t118, t142, -t143, 0, -t69 * t142 + t81 * t92, t69 * t143 + t83 * t92, -t24 * t146 + (t110 - t113) * t58, -t99 + (-t115 + (t24 + t127) * t81) * t80 + t138, t107 - t164, t100 + t144 + (-t89 - t151) * qJD(3), (-t68 - t120) * t129, -(-t50 * t126 + t135) * t68 + ((t56 * t69 + t156) * qJD(3) + (t148 + (t68 * t69 + t29) * t82) * qJD(4) + t106) * t83 + (t82 * t124 + t69 * t25 + t154 + ((-t69 * t147 + t40) * qJD(1) + t11) * qJD(3)) * t81, t136 * t68 + (-t69 * t114 + (t58 * t69 + t155) * qJD(3) - t102) * t83 + (-t80 * t124 - t69 * t24 + t153 + (-t134 * qJD(1) - t69 * t149 - t12) * qJD(3)) * t81, t19 * t56 + t44 * t25 - t3 * t68 + (t17 * t123 - t1) * t83 + (t17 * t125 + t158 + (qJD(1) * t15 + t7) * qJD(3)) * t81, t19 * t58 - t44 * t24 + t4 * t68 + (t17 * t122 + t2) * t83 + (-t17 * t126 + t157 + (-qJD(1) * t16 - t9) * qJD(3)) * t81, t15 * t24 - t16 * t25 - t3 * t58 - t4 * t56 + t96 * t128 + (t95 * qJD(4) - t1 * t82 - t2 * t80) * t81, t1 * t15 + t13 * t44 + t2 * t16 + t17 * t19 + t7 * t3 + t9 * t4; 0, 0, 0, 0, 0, 0, 0, 0, 0, -t143, -t142, 0, 0, 0, 0, 0, t6, t5, t6, t5, t99 + (t115 + (-t24 + t127) * t81) * t80 + t138, (-t95 * qJD(3) - t13) * t83 + (qJD(3) * t17 + t96 * qJD(4) - t1 * t80 + t2 * t82) * t81; 0, 0, 0, 0, -t81 * t85 * t83, t133 * t85, 0, 0, 0, t36 * qJD(3) - t64 * t130 - t31, -t64 * t120, -t58 * t149 - t24 * t80, -t166 * t80 + t167 * t82, -t68 * t125 + (t68 * t145 + (-t58 + t123) * t81) * qJD(1), t114 + (-t68 * t147 + (t56 + t122) * t81) * qJD(1), t68 * t130, -pkin(3) * t25 - t153 + t105 * t68 - t36 * t56 + (pkin(7) * t149 + t156) * qJD(4) + (-t11 * t81 + (-pkin(7) * t129 - t28 * t83) * t80) * qJD(1), pkin(3) * t24 + t154 - t137 * t68 - t36 * t58 + (-t80 * pkin(7) * t68 + t155) * qJD(4) + (-t28 * t145 + (-pkin(7) * t122 + t12) * t81) * qJD(1), -t157 - t20 * t56 + t73 * t25 + t140 * t68 + (t17 + t160) * t126 + (-t17 * t147 + (qJD(3) * t65 - t7) * t81) * qJD(1), t158 - t20 * t58 - t73 * t24 - t139 * t68 + (t58 * t161 + t17 * t82) * qJD(4) + (-t17 * t145 + (qJD(3) * t66 + t9) * t81) * qJD(1), t65 * t24 + t66 * t25 + t140 * t58 + t139 * t56 + (t68 * t7 + t2) * t82 + (-t1 + t159) * t80, t1 * t65 + t13 * t73 - t2 * t66 - t139 * t9 - t140 * t7 + (pkin(4) * t126 - t20) * t17; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t58 * t56, -t53 + t163, -t24 - t152, -t150 - t25, t108, -t12 * t68 - t28 * t58 + t87, -t11 * t68 + t28 * t56 + t102, 0.2e1 * t98 - t159 + (t104 - t17) * t58 + t86, -t163 * pkin(4) - t8 * t68 + (qJD(5) + t17) * t56 + t90, t24 * pkin(4) - t162 * t56, t162 * t9 + (-t17 * t58 + t1) * pkin(4); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t166, t167, -t53 - t163, t9 * t56 + t7 * t58 + t13;];
tauc_reg = t10;
