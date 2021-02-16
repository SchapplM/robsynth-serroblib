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
% tauc_reg [6x23]
%   minimal parameter regressor of coriolis joint torque vector

% Quelle: HybrDyn-Toolbox
% Datum: 2021-01-16 01:08
% Revision: 24b2e7d74a0c1a3b64fa2f8f5ad758691ad61af3 (2021-01-15)
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
% StartTime: 2021-01-16 01:06:19
% EndTime: 2021-01-16 01:06:27
% DurationCPUTime: 1.72s
% Computational Cost: add. (2208->250), mult. (5786->368), div. (0->0), fcn. (4641->12), ass. (0->145)
t100 = sin(pkin(11));
t109 = cos(qJ(2));
t101 = sin(pkin(6));
t157 = qJD(1) * t101;
t145 = t109 * t157;
t102 = cos(pkin(11));
t106 = sin(qJ(2));
t146 = t106 * t157;
t83 = t102 * t146;
t59 = t100 * t145 + t83;
t51 = qJD(2) * t59;
t105 = sin(qJ(4));
t156 = qJD(2) * t105;
t108 = cos(qJ(4));
t162 = cos(pkin(12));
t140 = t162 * t108;
t87 = qJD(2) * t140;
t99 = sin(pkin(12));
t69 = t99 * t156 - t87;
t68 = qJD(6) + t69;
t189 = qJD(6) - t68;
t169 = qJD(4) * pkin(4);
t188 = t105 * t169 - t59;
t81 = qJD(2) * pkin(2) + t145;
t49 = t100 * t81 + t83;
t47 = qJD(2) * pkin(8) + t49;
t135 = qJ(5) * qJD(2) + t47;
t103 = cos(pkin(6));
t88 = qJD(1) * t103 + qJD(3);
t33 = -t135 * t105 + t108 * t88;
t63 = (t100 * t106 - t102 * t109) * t101;
t62 = qJD(2) * t63;
t141 = t162 * t105;
t78 = t99 * t108 + t141;
t72 = t78 * qJD(2);
t52 = qJD(1) * t62;
t134 = qJD(2) * qJD(5) - t52;
t166 = t105 * t88;
t34 = t135 * t108 + t166;
t187 = -t34 * qJD(4) - t134 * t105;
t150 = pkin(4) * t156;
t14 = t33 * qJD(4) + t134 * t108;
t3 = t14 * t99 - t162 * t187;
t91 = pkin(4) * t99 + pkin(9);
t186 = (pkin(5) * t72 + pkin(9) * t69 + qJD(6) * t91 + t150) * t68 + t3;
t164 = t99 * t105;
t118 = t140 - t164;
t92 = pkin(2) * t100 + pkin(8);
t163 = qJ(5) + t92;
t76 = t163 * t108;
t40 = t162 * t76 - t163 * t164;
t71 = t78 * qJD(4);
t65 = qJD(2) * t71;
t131 = t3 * t78 - t40 * t65;
t147 = -pkin(4) * t108 - pkin(3);
t82 = t100 * t146;
t48 = t102 * t81 - t82;
t42 = t147 * qJD(2) + qJD(5) - t48;
t17 = pkin(5) * t69 - pkin(9) * t72 + t42;
t137 = qJD(4) * t163;
t117 = -t105 * qJD(5) - t108 * t137;
t56 = t108 * qJD(5) - t105 * t137;
t61 = t102 * t145 - t82;
t173 = -t99 * t117 + t118 * t61 - t162 * t56;
t183 = pkin(2) * t102;
t85 = t147 - t183;
t38 = -pkin(5) * t118 - pkin(9) * t78 + t85;
t144 = t162 * t14;
t4 = t187 * t99 + t144;
t175 = t99 * t34;
t28 = t33 + t169;
t9 = t162 * t28 - t175;
t7 = -qJD(4) * pkin(5) - t9;
t74 = t118 * qJD(4);
t185 = (qJD(6) * t17 + t4) * t118 + t7 * t74 + (-qJD(6) * t38 + t173) * t68 + t131;
t104 = sin(qJ(6));
t107 = cos(qJ(6));
t55 = qJD(4) * t104 + t107 * t72;
t151 = qJD(2) * qJD(4);
t143 = t105 * t151;
t86 = t99 * t143;
t66 = qJD(4) * t87 - t86;
t32 = qJD(6) * t55 + t104 * t66;
t176 = t65 * t78;
t129 = t68 * t74 + t176;
t153 = qJD(6) * t104;
t148 = t78 * t153;
t184 = -t107 * t129 + t68 * t148;
t182 = pkin(4) * t105;
t181 = t38 * t65;
t152 = t107 * qJD(4);
t53 = t104 * t72 - t152;
t180 = t53 * t68;
t179 = t53 * t72;
t178 = t55 * t68;
t177 = t55 * t72;
t31 = qJD(6) * t152 + t107 * t66 - t72 * t153;
t174 = -t118 * t31 + t55 * t71;
t26 = t162 * t34;
t10 = t99 * t28 + t26;
t172 = -t162 * t117 + t56 * t99 - t78 * t61;
t171 = pkin(5) * t71 - pkin(9) * t74 + t188;
t170 = t105 ^ 2 - t108 ^ 2;
t168 = t104 * t65;
t136 = t107 * t68;
t165 = t31 * t104;
t161 = qJD(6) * t78;
t111 = qJD(2) ^ 2;
t160 = t101 * t111;
t110 = qJD(4) ^ 2;
t159 = t110 * t105;
t158 = t110 * t108;
t155 = qJD(2) * t108;
t154 = qJD(4) * t108;
t44 = pkin(4) * t143 + t51;
t46 = -qJD(2) * pkin(3) - t48;
t139 = -qJD(2) * t46 + t52;
t130 = t118 * t32 - t71 * t53;
t8 = qJD(4) * pkin(9) + t10;
t128 = t104 * t8 - t107 * t17;
t2 = t104 * t17 + t107 * t8;
t64 = (t100 * t109 + t102 * t106) * t101;
t125 = t103 * t108 - t105 * t64;
t45 = t103 * t105 + t108 * t64;
t19 = t99 * t125 + t162 * t45;
t127 = t104 * t63 + t107 * t19;
t126 = -t104 * t19 + t107 * t63;
t121 = t107 * t65 + (-t104 * t69 - t153) * t68;
t120 = t110 * t92;
t119 = qJD(4) * (qJD(2) * (-pkin(3) - t183) + t46 + t61);
t12 = t162 * t33 - t175;
t116 = -t91 * t65 + (t12 + t7) * t68;
t115 = -qJD(4) * t45 + t105 * t62;
t113 = -t104 * t129 - t161 * t136;
t93 = -t162 * pkin(4) - pkin(5);
t60 = qJD(2) * t64;
t39 = t163 * t141 + t76 * t99;
t22 = qJD(4) * t125 - t108 * t62;
t18 = -t162 * t125 + t45 * t99;
t16 = pkin(5) * t65 - pkin(9) * t66 + t44;
t15 = t107 * t16;
t11 = t33 * t99 + t26;
t6 = t99 * t115 + t162 * t22;
t5 = -t162 * t115 + t22 * t99;
t1 = [0, 0, -t106 * t160, -t109 * t160, -t48 * t60 - t49 * t62 + t51 * t63 - t52 * t64, 0, 0, 0, 0, 0, -t60 * t155 + (-t64 * t154 + (-qJD(4) * t103 + 0.2e1 * t62) * t105) * qJD(4), -t22 * qJD(4) + (t105 * t60 + t63 * t154) * qJD(2), -qJD(4) * t5 + t60 * t69 + t63 * t65, -qJD(4) * t6 + t60 * t72 + t63 * t66, t18 * t66 - t19 * t65 + t5 * t72 - t6 * t69, t10 * t6 + t18 * t3 + t19 * t4 + t42 * t60 + t44 * t63 - t5 * t9, 0, 0, 0, 0, 0, (-qJD(6) * t127 - t104 * t6 + t107 * t60) * t68 + t126 * t65 + t5 * t53 + t18 * t32, -(qJD(6) * t126 + t104 * t60 + t107 * t6) * t68 - t127 * t65 + t5 * t55 + t18 * t31; 0, 0, 0, 0, t48 * t59 - t49 * t61 + (-t100 * t52 - t102 * t51) * pkin(2), 0.2e1 * t108 * t143, -0.2e1 * t170 * t151, t158, -t159, 0, t105 * t119 - t108 * t120, t105 * t120 + t108 * t119, t42 * t71 - t44 * t118 - t59 * t69 + t65 * t85 + (t69 * t182 - t172) * qJD(4), t42 * t74 + t44 * t78 - t59 * t72 + t66 * t85 + (t72 * t182 + t173) * qJD(4), -t10 * t71 + t118 * t4 + t172 * t72 + t173 * t69 + t39 * t66 - t74 * t9 + t131, -t173 * t10 - t172 * t9 + t188 * t42 + t3 * t39 + t4 * t40 + t44 * t85, -t55 * t148 + (t31 * t78 + t55 * t74) * t107, (-t104 * t55 - t107 * t53) * t74 + (-t165 - t107 * t32 + (t104 * t53 - t107 * t55) * qJD(6)) * t78, t174 - t184, t113 + t130, -t118 * t65 + t68 * t71, -t128 * t71 - t15 * t118 + t39 * t32 + t172 * t53 + (t181 + t171 * t68 + (t118 * t8 - t40 * t68 + t7 * t78) * qJD(6)) * t107 + t185 * t104, -t2 * t71 + t39 * t31 + t172 * t55 + (-t181 + (-qJD(6) * t8 + t16) * t118 - t7 * t161 + (qJD(6) * t40 - t171) * t68) * t104 + t185 * t107; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t159, -t158, -t71 * qJD(4), -t74 * qJD(4), -t118 * t66 - t69 * t74 + t71 * t72 - t176, t10 * t74 - t118 * t3 + t4 * t78 - t71 * t9, 0, 0, 0, 0, 0, t113 - t130, t174 + t184; 0, 0, 0, 0, 0, -t105 * t111 * t108, t170 * t111, 0, 0, 0, t139 * t105, t139 * t108, qJD(4) * t11 - t69 * t150 - t42 * t72 - t3, -t144 + t42 * t69 + (-qJD(2) * pkin(4) * t72 + t99 * t134) * t105 + (-t99 * (-qJ(5) * t155 - t108 * t47 - t166) + t12) * qJD(4), (t10 - t11) * t72 + (t12 - t9) * t69 + (-t162 * t66 - t65 * t99) * pkin(4), -t10 * t12 + t9 * t11 + (-t156 * t42 - t162 * t3 + t4 * t99) * pkin(4), t136 * t55 + t165, (t31 - t180) * t107 + (-t32 - t178) * t104, t136 * t68 + t168 - t177, t121 + t179, -t68 * t72, t116 * t104 - t186 * t107 - t11 * t53 + t128 * t72 + t93 * t32, t186 * t104 + t116 * t107 - t11 * t55 + t2 * t72 + t93 * t31; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0.2e1 * t72 * qJD(4), -t86 + (t87 - t69) * qJD(4), -t69 ^ 2 - t72 ^ 2, t10 * t69 + t72 * t9 + t44, 0, 0, 0, 0, 0, t121 - t179, -t107 * t68 ^ 2 - t168 - t177; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t55 * t53, -t53 ^ 2 + t55 ^ 2, t31 + t180, t178 - t32, t65, -t104 * t4 - t189 * t2 - t55 * t7 + t15, -t104 * t16 - t107 * t4 + t189 * t128 + t53 * t7;];
tauc_reg = t1;
