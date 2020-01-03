% Calculate minimal parameter regressor of coriolis joint torque vector for
% S5RRPRP6
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% qJD [5x1]
%   Generalized joint velocities
% pkin [8x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d1,d2,d4,theta3]';
% 
% Output:
% tauc_reg [5x21]
%   minimal parameter regressor of coriolis joint torque vector

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-31 19:59
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauc_reg = S5RRPRP6_coriolisvecJ_fixb_regmin_slag_vp(qJ, qJD, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(8,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RRPRP6_coriolisvecJ_fixb_regmin_slag_vp: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RRPRP6_coriolisvecJ_fixb_regmin_slag_vp: qJD has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S5RRPRP6_coriolisvecJ_fixb_regmin_slag_vp: pkin has to be [8x1] (double)');

%% Symbolic Calculation
% From coriolisvec_joint_fixb_regressor_minpar_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 19:58:27
% EndTime: 2019-12-31 19:58:32
% DurationCPUTime: 1.26s
% Computational Cost: add. (2062->230), mult. (5323->328), div. (0->0), fcn. (3693->6), ass. (0->130)
t132 = (qJD(1) * qJD(2));
t164 = -2 * t132;
t102 = sin(qJ(4));
t100 = sin(pkin(8));
t103 = sin(qJ(2));
t136 = qJD(1) * t103;
t101 = cos(pkin(8));
t105 = cos(qJ(2));
t140 = t101 * t105;
t71 = qJD(1) * t140 - t100 * t136;
t68 = qJD(4) - t71;
t121 = t102 * t68;
t104 = cos(qJ(4));
t83 = t100 * t105 + t101 * t103;
t73 = t83 * qJD(1);
t56 = t102 * qJD(2) + t104 * t73;
t163 = t56 * t121;
t127 = -t105 * pkin(2) - pkin(1);
t115 = t127 * qJD(1);
t87 = qJD(3) + t115;
t27 = -t71 * pkin(3) - t73 * pkin(7) + t87;
t153 = -qJ(3) - pkin(6);
t89 = t153 * t105;
t86 = qJD(1) * t89;
t146 = t101 * t86;
t147 = qJD(2) * pkin(2);
t88 = t153 * t103;
t85 = qJD(1) * t88;
t79 = t85 + t147;
t46 = t100 * t79 - t146;
t41 = qJD(2) * pkin(7) + t46;
t13 = t102 * t27 + t104 * t41;
t125 = t105 * t132;
t126 = t103 * t132;
t112 = -t100 * t126 + t101 * t125;
t72 = t83 * qJD(2);
t66 = qJD(1) * t72;
t92 = pkin(2) * t126;
t26 = t66 * pkin(3) - t112 * pkin(7) + t92;
t19 = t104 * t26;
t124 = qJD(2) * t153;
t69 = t105 * qJD(3) + t103 * t124;
t63 = t69 * qJD(1);
t70 = -t103 * qJD(3) + t105 * t124;
t64 = t70 * qJD(1);
t22 = t100 * t64 + t101 * t63;
t108 = -t13 * qJD(4) - t102 * t22 + t19;
t133 = t104 * qJD(2);
t135 = qJD(4) * t102;
t23 = -qJD(4) * t133 - t104 * t112 + t73 * t135;
t1 = t66 * pkin(4) + t23 * qJ(5) - t56 * qJD(5) + t108;
t54 = t102 * t73 - t133;
t8 = -t54 * qJ(5) + t13;
t162 = t68 * t8 + t1;
t161 = t56 ^ 2;
t12 = -t102 * t41 + t104 * t27;
t7 = -t56 * qJ(5) + t12;
t5 = t68 * pkin(4) + t7;
t160 = t5 - t7;
t94 = t100 * pkin(2) + pkin(7);
t143 = qJ(5) + t94;
t123 = qJD(4) * t143;
t141 = qJ(5) * t104;
t35 = pkin(2) * t136 + t73 * pkin(3) - t71 * pkin(7);
t31 = t104 * t35;
t76 = t100 * t86;
t48 = t101 * t85 + t76;
t159 = -t73 * pkin(4) - t104 * t123 + t71 * t141 - t31 + (-qJD(5) + t48) * t102;
t21 = t100 * t63 - t101 * t64;
t158 = t21 * t83;
t45 = t101 * t79 + t76;
t40 = -qJD(2) * pkin(3) - t45;
t82 = t100 * t103 - t140;
t75 = t82 * qJD(2);
t157 = t40 * t75;
t156 = t54 * t71;
t155 = t56 * t73;
t154 = t73 * t54;
t134 = qJD(4) * t104;
t109 = t102 * t112;
t24 = t56 * qJD(4) + t109;
t152 = -t102 * t24 - t54 * t134;
t151 = t102 * t35 + t104 * t48;
t44 = t82 * pkin(3) - t83 * pkin(7) + t127;
t52 = t100 * t88 - t101 * t89;
t49 = t104 * t52;
t150 = t102 * t44 + t49;
t142 = qJ(5) * t102;
t149 = t104 * qJD(5) - t102 * t123 + t71 * t142 - t151;
t148 = t103 ^ 2 - t105 ^ 2;
t145 = t102 * t66;
t144 = t23 * t102;
t107 = qJD(1) ^ 2;
t139 = t105 * t107;
t106 = qJD(2) ^ 2;
t138 = t106 * t103;
t137 = t106 * t105;
t34 = t100 * t70 + t101 * t69;
t130 = t103 * t147;
t36 = t72 * pkin(3) + t75 * pkin(7) + t130;
t131 = t102 * t36 + t104 * t34 + t44 * t134;
t129 = t83 * t135;
t128 = t83 * t134;
t95 = -t101 * pkin(2) - pkin(3);
t33 = t100 * t69 - t101 * t70;
t47 = t100 * t85 - t146;
t51 = -t100 * t89 - t101 * t88;
t122 = t104 * t68;
t120 = pkin(1) * t164;
t119 = qJD(4) * t94 * t68 + t21;
t111 = t102 * t26 + t104 * t22 + t27 * t134 - t41 * t135;
t2 = -t24 * qJ(5) - t54 * qJD(5) + t111;
t118 = -t68 * t5 + t2;
t117 = -t52 * t66 + t158;
t116 = t66 * t83 - t68 * t75;
t114 = qJ(5) * t75 - qJD(5) * t83;
t113 = t104 * t66 + t71 * t121 - t68 * t135;
t11 = t24 * pkin(4) + t21;
t110 = t68 * t40 - t94 * t66;
t81 = t143 * t104;
t80 = t143 * t102;
t53 = t54 ^ 2;
t39 = t104 * t44;
t32 = t104 * t36;
t15 = t54 * pkin(4) + qJD(5) + t40;
t14 = -t83 * t142 + t150;
t9 = t82 * pkin(4) - t102 * t52 - t83 * t141 + t39;
t4 = -qJ(5) * t128 + (-qJD(4) * t52 + t114) * t102 + t131;
t3 = t72 * pkin(4) - t102 * t34 + t32 + t114 * t104 + (-t49 + (qJ(5) * t83 - t44) * t102) * qJD(4);
t6 = [0, 0, 0, 0.2e1 * t103 * t125, t148 * t164, t137, -t138, 0, -pkin(6) * t137 + t103 * t120, pkin(6) * t138 + t105 * t120, t51 * t112 - t22 * t82 + t33 * t73 + t34 * t71 + t45 * t75 - t46 * t72 + t117, t21 * t51 + t22 * t52 - t45 * t33 + t46 * t34 + (t87 + t115) * t130, -t56 * t129 + (-t23 * t83 - t56 * t75) * t104, -(-t102 * t56 - t104 * t54) * t75 + (t144 - t104 * t24 + (t102 * t54 - t104 * t56) * qJD(4)) * t83, t116 * t104 - t68 * t129 - t23 * t82 + t56 * t72, -t116 * t102 - t68 * t128 - t24 * t82 - t54 * t72, t66 * t82 + t68 * t72, (-t134 * t52 + t32) * t68 + t39 * t66 + (-t134 * t41 + t19) * t82 + t12 * t72 + t33 * t54 + t51 * t24 + t40 * t128 + ((-qJD(4) * t44 - t34) * t68 + (-qJD(4) * t27 - t22) * t82 - t157 + t117) * t102, -(-t135 * t52 + t131) * t68 - t150 * t66 - t111 * t82 - t13 * t72 + t33 * t56 - t51 * t23 - t40 * t129 + (-t157 + t158) * t104, -t14 * t24 + t9 * t23 - t3 * t56 - t4 * t54 - (-t102 * t8 - t104 * t5) * t75 + (-t1 * t104 - t2 * t102 + (t102 * t5 - t104 * t8) * qJD(4)) * t83, t2 * t14 + t8 * t4 + t1 * t9 + t5 * t3 + t11 * (t102 * t83 * pkin(4) + t51) + t15 * ((-t102 * t75 + t128) * pkin(4) + t33); 0, 0, 0, -t103 * t139, t148 * t107, 0, 0, 0, t107 * pkin(1) * t103, pkin(1) * t139, (t46 - t47) * t73 + (-t48 + t45) * t71 + (-t100 * t66 - t101 * t112) * pkin(2), t45 * t47 - t46 * t48 + (t100 * t22 - t101 * t21 - t87 * t136) * pkin(2), t56 * t122 - t144, (-t23 + t156) * t104 - t163 + t152, t68 * t122 + t145 - t155, t113 + t154, -t68 * t73, -t12 * t73 + t95 * t24 - t31 * t68 - t47 * t54 - t119 * t104 + (t48 * t68 + t110) * t102, t119 * t102 + t110 * t104 + t13 * t73 + t151 * t68 - t95 * t23 - t47 * t56, -t162 * t102 + t118 * t104 - t149 * t54 - t159 * t56 - t80 * t23 - t81 * t24, t2 * t81 - t1 * t80 + t11 * (-t104 * pkin(4) + t95) + t149 * t8 + t159 * t5 + (pkin(4) * t121 - t47) * t15; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t71 ^ 2 - t73 ^ 2, t45 * t73 - t46 * t71 + t92, 0, 0, 0, 0, 0, t113 - t154, -t104 * t68 ^ 2 - t145 - t155, (t23 + t156) * t104 + t163 + t152, t118 * t102 + t162 * t104 - t15 * t73; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t56 * t54, -t53 + t161, t54 * t68 - t23, -t109 + (-qJD(4) + t68) * t56, t66, t13 * t68 - t40 * t56 + t108, t12 * t68 + t40 * t54 - t111, pkin(4) * t23 - t160 * t54, t160 * t8 + (-t15 * t56 + t1) * pkin(4); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t53 - t161, t5 * t56 + t8 * t54 + t11;];
tauc_reg = t6;
