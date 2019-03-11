% Calculate minimal parameter regressor of coriolis joint torque vector for
% S6PRPRPR1
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% qJD [6x1]
%   Generalized joint velocities
% pkin [12x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,alpha2,d2,d4,d6,theta1,theta3,theta5]';
% 
% Output:
% tauc_reg [6x21]
%   minimal parameter regressor of coriolis joint torque vector

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-08 19:28
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauc_reg = S6PRPRPR1_coriolisvecJ_fixb_regmin_slag_vp(qJ, qJD, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(12,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6PRPRPR1_coriolisvecJ_fixb_regmin_slag_vp: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6PRPRPR1_coriolisvecJ_fixb_regmin_slag_vp: qJD has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [12 1]), ...
  'S6PRPRPR1_coriolisvecJ_fixb_regmin_slag_vp: pkin has to be [12x1] (double)');

%% Symbolic Calculation
% From coriolisvec_joint_fixb_regressor_minpar_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-08 19:28:07
% EndTime: 2019-03-08 19:28:08
% DurationCPUTime: 1.31s
% Computational Cost: add. (2023->221), mult. (5254->331), div. (0->0), fcn. (4219->12), ass. (0->136)
t109 = cos(qJ(2));
t100 = sin(pkin(6));
t154 = qJD(1) * t100;
t144 = t109 * t154;
t102 = cos(pkin(11));
t106 = sin(qJ(2));
t145 = t106 * t154;
t83 = t102 * t145;
t99 = sin(pkin(11));
t60 = t99 * t144 + t83;
t52 = qJD(2) * t60;
t105 = sin(qJ(4));
t153 = qJD(2) * t105;
t101 = cos(pkin(12));
t108 = cos(qJ(4));
t157 = t101 * t108;
t98 = sin(pkin(12));
t70 = qJD(2) * t157 - t98 * t153;
t69 = qJD(6) - t70;
t183 = qJD(6) - t69;
t152 = qJD(4) * t105;
t182 = pkin(4) * t152 - t60;
t81 = qJD(2) * pkin(2) + t144;
t50 = t99 * t81 + t83;
t135 = t50 + (pkin(8) + qJ(5)) * qJD(2);
t103 = cos(pkin(6));
t87 = qJD(1) * t103 + qJD(3);
t34 = -t135 * t105 + t108 * t87;
t64 = (t102 * t109 - t106 * t99) * t100;
t62 = qJD(2) * t64;
t53 = qJD(1) * t62;
t134 = qJD(2) * qJD(5) + t53;
t35 = t105 * t87 + t135 * t108;
t181 = -t35 * qJD(4) - t134 * t105;
t14 = qJD(4) * t34 + t134 * t108;
t3 = -t101 * t181 + t98 * t14;
t78 = t101 * t105 + t108 * t98;
t72 = t78 * qJD(2);
t90 = pkin(4) * t98 + pkin(9);
t180 = (pkin(4) * t153 + pkin(5) * t72 - pkin(9) * t70 + qJD(6) * t90) * t69 + t3;
t125 = -t105 * t98 + t157;
t91 = pkin(2) * t99 + pkin(8);
t160 = qJ(5) + t91;
t140 = t160 * t105;
t76 = t160 * t108;
t41 = t101 * t76 - t98 * t140;
t71 = t78 * qJD(4);
t66 = qJD(2) * t71;
t131 = t3 * t78 - t41 * t66;
t137 = qJD(4) * t160;
t116 = -t105 * qJD(5) - t108 * t137;
t57 = t108 * qJD(5) - t105 * t137;
t82 = t99 * t145;
t63 = t102 * t144 - t82;
t167 = t101 * t57 + t98 * t116 - t125 * t63;
t146 = -pkin(4) * t108 - pkin(3);
t49 = t102 * t81 - t82;
t43 = t146 * qJD(2) + qJD(5) - t49;
t17 = -t70 * pkin(5) - t72 * pkin(9) + t43;
t177 = pkin(2) * t102;
t120 = t146 - t177;
t39 = -pkin(5) * t125 - pkin(9) * t78 + t120;
t4 = t101 * t14 + t181 * t98;
t176 = t35 * t98;
t29 = qJD(4) * pkin(4) + t34;
t9 = t101 * t29 - t176;
t7 = -qJD(4) * pkin(5) - t9;
t74 = t125 * qJD(4);
t179 = (qJD(6) * t17 + t4) * t125 + t7 * t74 + (-qJD(6) * t39 - t167) * t69 + t131;
t104 = sin(qJ(6));
t107 = cos(qJ(6));
t56 = qJD(4) * t104 + t107 * t72;
t149 = qJD(2) * qJD(4);
t142 = t108 * t149;
t143 = t105 * t149;
t67 = t101 * t142 - t98 * t143;
t33 = t56 * qJD(6) + t104 * t67;
t170 = t66 * t78;
t129 = t69 * t74 + t170;
t151 = qJD(6) * t104;
t147 = t78 * t151;
t178 = -t129 * t107 + t69 * t147;
t175 = t39 * t66;
t150 = t107 * qJD(4);
t54 = t104 * t72 - t150;
t174 = t54 * t69;
t173 = t54 * t72;
t172 = t56 * t69;
t171 = t56 * t72;
t32 = qJD(6) * t150 + t107 * t67 - t72 * t151;
t169 = -t125 * t32 + t56 * t71;
t27 = t101 * t35;
t10 = t98 * t29 + t27;
t168 = -t101 * t116 + t98 * t57 - t78 * t63;
t166 = pkin(5) * t71 - pkin(9) * t74 + t182;
t165 = t105 ^ 2 - t108 ^ 2;
t164 = t104 * t66;
t136 = t107 * t69;
t161 = t32 * t104;
t159 = qJD(6) * t78;
t111 = qJD(2) ^ 2;
t158 = t100 * t111;
t110 = qJD(4) ^ 2;
t156 = t110 * t105;
t155 = t110 * t108;
t45 = pkin(4) * t143 + t52;
t47 = -qJD(2) * pkin(3) - t49;
t139 = -qJD(2) * t47 - t53;
t130 = t125 * t33 - t71 * t54;
t8 = qJD(4) * pkin(9) + t10;
t128 = t104 * t8 - t107 * t17;
t2 = t104 * t17 + t107 * t8;
t65 = (t102 * t106 + t109 * t99) * t100;
t123 = t103 * t108 - t105 * t65;
t46 = t103 * t105 + t108 * t65;
t19 = t101 * t46 + t98 * t123;
t127 = -t104 * t64 + t107 * t19;
t126 = -t104 * t19 - t107 * t64;
t119 = t107 * t66 + (t104 * t70 - t151) * t69;
t118 = t110 * t91;
t117 = qJD(4) * (qJD(2) * (-pkin(3) - t177) + t47 + t63);
t12 = t101 * t34 - t176;
t115 = -t90 * t66 + (t12 + t7) * t69;
t113 = -t129 * t104 - t159 * t136;
t92 = -pkin(4) * t101 - pkin(5);
t61 = qJD(2) * t65;
t40 = t101 * t140 + t98 * t76;
t23 = t123 * qJD(4) + t62 * t108;
t22 = -t46 * qJD(4) - t62 * t105;
t18 = -t101 * t123 + t98 * t46;
t16 = pkin(5) * t66 - pkin(9) * t67 + t45;
t15 = t107 * t16;
t11 = t34 * t98 + t27;
t6 = t101 * t23 + t22 * t98;
t5 = -t101 * t22 + t23 * t98;
t1 = [0, 0, -t106 * t158, -t109 * t158, -t49 * t61 + t50 * t62 - t52 * t64 + t53 * t65, 0, 0, 0, 0, 0, t22 * qJD(4) + (-t108 * t61 - t64 * t152) * qJD(2), -t23 * qJD(4) + (-qJD(4) * t108 * t64 + t105 * t61) * qJD(2), t18 * t67 - t19 * t66 + t5 * t72 + t6 * t70, t10 * t6 + t18 * t3 + t19 * t4 + t43 * t61 - t45 * t64 - t5 * t9, 0, 0, 0, 0, 0 (-t127 * qJD(6) - t104 * t6 + t61 * t107) * t69 + t126 * t66 + t5 * t54 + t18 * t33 -(qJD(6) * t126 + t61 * t104 + t107 * t6) * t69 - t127 * t66 + t5 * t56 + t18 * t32; 0, 0, 0, 0, t49 * t60 - t50 * t63 + (-t102 * t52 + t53 * t99) * pkin(2), 0.2e1 * t105 * t142, -0.2e1 * t165 * t149, t155, -t156, 0, t105 * t117 - t118 * t108, t118 * t105 + t108 * t117, -t10 * t71 + t125 * t4 + t167 * t70 + t168 * t72 + t40 * t67 - t9 * t74 + t131, t167 * t10 + t45 * t120 - t168 * t9 + t182 * t43 + t3 * t40 + t4 * t41, -t56 * t147 + (t32 * t78 + t56 * t74) * t107 (-t104 * t56 - t107 * t54) * t74 + (-t161 - t107 * t33 + (t104 * t54 - t107 * t56) * qJD(6)) * t78, t169 - t178, t113 + t130, -t125 * t66 + t69 * t71, -t128 * t71 - t15 * t125 + t40 * t33 + t168 * t54 + (t175 + t166 * t69 + (t125 * t8 - t41 * t69 + t7 * t78) * qJD(6)) * t107 + t179 * t104, -t2 * t71 + t40 * t32 + t168 * t56 + (-t175 + (-qJD(6) * t8 + t16) * t125 - t7 * t159 + (qJD(6) * t41 - t166) * t69) * t104 + t179 * t107; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t156, -t155, -t125 * t67 + t70 * t74 + t71 * t72 - t170, t10 * t74 - t125 * t3 + t4 * t78 - t71 * t9, 0, 0, 0, 0, 0, t113 - t130, t169 + t178; 0, 0, 0, 0, 0, -t105 * t111 * t108, t165 * t111, 0, 0, 0, t139 * t105, t139 * t108 (t10 - t11) * t72 + (-t12 + t9) * t70 + (-t101 * t67 - t66 * t98) * pkin(4), -t10 * t12 + t9 * t11 + (-t101 * t3 - t43 * t153 + t4 * t98) * pkin(4), t56 * t136 + t161 (t32 - t174) * t107 + (-t33 - t172) * t104, t69 * t136 + t164 - t171, t119 + t173, -t69 * t72, t115 * t104 - t107 * t180 - t11 * t54 + t128 * t72 + t92 * t33, t104 * t180 + t115 * t107 - t11 * t56 + t2 * t72 + t92 * t32; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t70 ^ 2 - t72 ^ 2, -t10 * t70 + t72 * t9 + t45, 0, 0, 0, 0, 0, t119 - t173, -t69 ^ 2 * t107 - t164 - t171; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t56 * t54, -t54 ^ 2 + t56 ^ 2, t32 + t174, t172 - t33, t66, -t104 * t4 - t183 * t2 - t7 * t56 + t15, -t104 * t16 - t107 * t4 + t183 * t128 + t7 * t54;];
tauc_reg  = t1;
