% Calculate minimal parameter regressor of coriolis joint torque vector for
% S6RRPPRP3
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% qJD [6x1]
%   Generalized joint velocities
% pkin [8x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d2,d5]';
% 
% Output:
% tauc_reg [6x27]
%   minimal parameter regressor of coriolis joint torque vector

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-09 08:36
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauc_reg = S6RRPPRP3_coriolisvecJ_fixb_regmin_slag_vp(qJ, qJD, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(8,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRPPRP3_coriolisvecJ_fixb_regmin_slag_vp: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RRPPRP3_coriolisvecJ_fixb_regmin_slag_vp: qJD has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S6RRPPRP3_coriolisvecJ_fixb_regmin_slag_vp: pkin has to be [8x1] (double)');

%% Symbolic Calculation
% From coriolisvec_joint_fixb_regressor_minpar_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 08:35:33
% EndTime: 2019-03-09 08:35:39
% DurationCPUTime: 1.95s
% Computational Cost: add. (1888->294), mult. (4046->395), div. (0->0), fcn. (2106->4), ass. (0->159)
t106 = -pkin(2) - pkin(3);
t141 = qJD(2) * t106;
t103 = sin(qJ(2));
t161 = qJD(1) * t103;
t83 = pkin(7) * t161;
t59 = -qJ(4) * t161 + t83;
t145 = qJD(3) + t59;
t43 = t141 + t145;
t102 = sin(qJ(5));
t104 = cos(qJ(5));
t152 = t104 * qJD(2);
t105 = cos(qJ(2));
t160 = qJD(1) * t105;
t56 = t102 * t160 - t152;
t75 = qJD(5) + t161;
t187 = t56 * t75;
t154 = qJD(5) * t105;
t139 = t102 * t154;
t34 = -qJD(5) * t152 + (t103 * t152 + t139) * qJD(1);
t193 = t34 - t187;
t159 = qJD(2) * t102;
t57 = t104 * t160 + t159;
t186 = t57 * t75;
t150 = qJD(1) * qJD(2);
t137 = t103 * t150;
t35 = qJD(5) * t57 - t102 * t137;
t113 = t35 - t186;
t95 = -pkin(8) + t106;
t114 = pkin(4) * t105 + t103 * t95;
t112 = t114 * qJD(2);
t136 = t105 * t150;
t87 = t103 * qJD(3);
t179 = qJ(3) * t136 + qJD(1) * t87;
t17 = qJD(1) * t112 + t179;
t16 = t104 * t17;
t157 = qJD(2) * t105;
t140 = qJ(4) * t157;
t153 = t103 * qJD(4);
t74 = pkin(7) * t136;
t40 = t74 + (-t140 - t153) * qJD(1);
t124 = pkin(4) * t103 + pkin(8) * t105;
t53 = -qJD(1) * pkin(1) - pkin(2) * t160 - qJ(3) * t161;
t39 = pkin(3) * t160 + qJD(4) - t53;
t24 = qJD(1) * t124 + t39;
t37 = qJD(2) * t95 + t145;
t9 = t102 * t24 + t104 * t37;
t111 = -qJD(5) * t9 - t102 * t40 + t16;
t1 = pkin(5) * t136 - t34 * qJ(6) + t57 * qJD(6) + t111;
t155 = qJD(5) * t104;
t156 = qJD(5) * t102;
t115 = -t102 * t17 - t104 * t40 - t24 * t155 + t156 * t37;
t2 = qJ(6) * t35 + qJD(6) * t56 - t115;
t8 = -t102 * t37 + t104 * t24;
t6 = qJ(6) * t57 + t8;
t5 = pkin(5) * t75 + t6;
t7 = qJ(6) * t56 + t9;
t195 = -(t7 * t75 + t1) * t102 - (t5 * t75 - t2) * t104;
t194 = -0.2e1 * t150;
t84 = pkin(7) * t160;
t61 = -qJ(4) * t160 + t84;
t96 = qJD(2) * qJ(3);
t51 = -t61 - t96;
t158 = qJD(2) * t103;
t147 = pkin(7) * t158;
t116 = t105 * qJD(4) + t147;
t72 = qJ(4) * t137;
t94 = qJD(2) * qJD(3);
t32 = qJD(1) * t116 - t72 - t94;
t190 = t57 ^ 2;
t189 = t5 - t6;
t188 = pkin(5) * t102;
t185 = pkin(7) - qJ(4);
t101 = qJ(3) + pkin(4);
t169 = qJ(6) - t95;
t130 = qJD(5) * t169;
t151 = t104 * qJD(6);
t80 = qJ(3) * t160;
t27 = qJD(1) * t114 + t80;
t182 = t102 * t27 + t104 * t61;
t184 = -t151 - t182 + (qJ(6) * t161 + t130) * t102;
t167 = t103 * t104;
t118 = pkin(5) * t105 - qJ(6) * t167;
t135 = -t102 * t61 + t104 * t27;
t183 = -qJD(1) * t118 + t102 * qJD(6) + t104 * t130 - t135;
t66 = -t105 * pkin(2) - t103 * qJ(3) - pkin(1);
t54 = t105 * pkin(3) - t66;
t36 = t124 + t54;
t68 = t185 * t103;
t58 = t104 * t68;
t181 = t102 * t36 + t58;
t178 = qJ(3) * t157 + t87;
t97 = t103 ^ 2;
t98 = t105 ^ 2;
t177 = t97 - t98;
t176 = qJD(2) * pkin(2);
t175 = t102 * t75;
t174 = t103 * t75;
t173 = t104 * t75;
t172 = t32 * t102;
t171 = t32 * t104;
t170 = t34 * t102;
t168 = qJ(6) * t105;
t108 = qJD(1) ^ 2;
t166 = t105 * t108;
t107 = qJD(2) ^ 2;
t165 = t107 * t103;
t164 = t107 * t105;
t162 = -qJD(4) - t39;
t23 = t112 + t178;
t49 = t157 * t185 - t153;
t149 = t102 * t23 + t104 * t49 + t36 * t155;
t148 = t102 * t174;
t146 = t75 * t167;
t144 = t75 * t156;
t143 = t75 * t155;
t142 = t103 * t166;
t138 = t104 * t154;
t134 = -t102 * t68 + t104 * t36;
t125 = t103 * t141;
t28 = qJD(1) * t125 + t179;
t38 = t125 + t178;
t133 = qJD(1) * t38 + t28;
t132 = qJD(1) * t54 + t39;
t131 = qJD(1) * t66 + t53;
t129 = t162 * t103;
t128 = pkin(1) * t194;
t127 = qJD(3) - t176;
t126 = 0.2e1 * t136;
t123 = -t102 * t5 + t104 * t7;
t122 = qJD(1) * t98 - t174;
t41 = pkin(2) * t137 - t179;
t50 = pkin(2) * t158 - t178;
t119 = -pkin(7) * t107 - qJD(1) * t50 - t41;
t46 = qJD(2) * pkin(4) - t51;
t117 = -t103 * t46 - t157 * t95;
t110 = qJD(5) * t123 + t1 * t104 + t2 * t102;
t13 = -t35 * pkin(5) - t32;
t64 = -pkin(7) * t137 + t94;
t65 = t127 + t83;
t67 = t84 + t96;
t109 = t64 * t105 + (t105 * t65 + (-t67 + t84) * t103) * qJD(2);
t91 = 0.2e1 * t94;
t90 = t105 * qJ(4);
t81 = qJ(4) * t158;
t71 = -t108 * t97 - t107;
t69 = pkin(7) * t105 - t90;
t63 = t169 * t104;
t62 = t169 * t102;
t60 = pkin(2) * t161 - t80;
t55 = t56 ^ 2;
t48 = t116 - t81;
t47 = t106 * t161 + t80;
t22 = -pkin(5) * t56 + qJD(6) + t46;
t20 = t104 * t23;
t14 = t102 * t168 + t181;
t12 = pkin(5) * t103 + t104 * t168 + t134;
t4 = qJ(6) * t138 + (-qJ(6) * t158 - qJD(5) * t68 + qJD(6) * t105) * t102 + t149;
t3 = t105 * t151 - t102 * t49 + t20 + t118 * qJD(2) + (-t58 + (-t36 - t168) * t102) * qJD(5);
t10 = [0, 0, 0, t103 * t126, t177 * t194, t164, -t165, 0, -pkin(7) * t164 + t103 * t128, pkin(7) * t165 + t105 * t128, t105 * t119 + t131 * t158, t109, t103 * t119 - t131 * t157, pkin(7) * t109 + t41 * t66 + t53 * t50, t133 * t103 + (t105 * t132 - t48) * qJD(2), -t133 * t105 + (t103 * t132 + t49) * qJD(2), -t40 * t103 + t32 * t105 + (-t103 * t51 - t105 * t43) * qJD(2) + (-t103 * t49 + t105 * t48 + (t103 * t69 - t105 * t68) * qJD(2)) * qJD(1), t28 * t54 - t32 * t69 + t38 * t39 + t40 * t68 + t43 * t49 + t48 * t51, -t57 * t139 + (-t105 * t34 - t158 * t57) * t104 (t102 * t57 + t104 * t56) * t158 + (t170 - t104 * t35 + (t102 * t56 - t104 * t57) * qJD(5)) * t105, t75 * t139 + t34 * t103 + (-t104 * t122 - t105 * t57) * qJD(2), t75 * t138 + t35 * t103 + (t102 * t122 + t105 * t56) * qJD(2) (t75 + t161) * t157 (-t155 * t68 + t20) * t75 + (-t155 * t37 + t16) * t103 + t48 * t56 - t69 * t35 + ((-qJD(5) * t36 - t49) * t75 + (qJD(2) * t46 - qJD(5) * t24 - t40) * t103) * t102 + (-t46 * t155 + t172 + (qJD(1) * t134 + t8) * qJD(2)) * t105 -(-t156 * t68 + t149) * t75 + t48 * t57 + t69 * t34 + (t152 * t46 + t115) * t103 + (t46 * t156 + t171 + (-qJD(1) * t181 - t9) * qJD(2)) * t105, -t12 * t34 + t14 * t35 + t3 * t57 + t4 * t56 + (-t102 * t7 - t104 * t5) * t158 + t110 * t105, t2 * t14 + t7 * t4 + t1 * t12 + t5 * t3 - t13 * t90 + t22 * (t158 * t188 - t147 + t81) + (t13 * (pkin(7) - t188) + t22 * (-pkin(5) * t155 - qJD(4))) * t105; 0, 0, 0, -t142, t177 * t108, 0, 0, 0, t108 * pkin(1) * t103, pkin(1) * t166 (-t103 * t53 + t105 * t60) * qJD(1) ((t67 - t96) * t103 + (t127 - t65) * t105) * qJD(1), t91 + (t103 * t60 + t105 * t53) * qJD(1), t64 * qJ(3) + t67 * qJD(3) - t53 * t60 + (t103 * t67 + (-t65 - t176) * t105) * qJD(1) * pkin(7), t59 * qJD(2) + t72 + t91 + (t162 * t105 + (-pkin(7) * qJD(2) - t47) * t103) * qJD(1), -t61 * qJD(2) + t74 + ((-qJ(4) * qJD(2) + t47) * t105 + t129) * qJD(1), 0, -t32 * qJ(3) + t40 * t106 - t145 * t51 - t39 * t47 - t43 * t61, t173 * t57 - t170 (-t34 - t187) * t104 + (-t35 - t186) * t102, -t143 + (-t146 + (t57 - t159) * t105) * qJD(1), t144 + (t148 + (-t56 - t152) * t105) * qJD(1), -t75 * t160, -t101 * t35 - t171 - t135 * t75 - t145 * t56 + (-t102 * t46 - t173 * t95) * qJD(5) + (t102 * t117 - t8 * t105) * qJD(1), t101 * t34 + t172 + t182 * t75 - t145 * t57 + (-t104 * t46 + t175 * t95) * qJD(5) + (t104 * t117 + t9 * t105) * qJD(1), t183 * t57 + t184 * t56 - t62 * t34 - t63 * t35 - t195, -t2 * t63 + t1 * t62 + t13 * (pkin(5) * t104 + t101) + t184 * t7 + t183 * t5 + (-pkin(5) * t175 + t145) * t22; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t142, 0, t71, -qJD(2) * t67 + t161 * t53 + t74, t71, t142, 0, t51 * qJD(2) + t74 + (t129 - t140) * qJD(1), 0, 0, 0, 0, 0, -t143 + qJD(2) * t56 + (-t102 * t157 - t146) * qJD(1), t144 + qJD(2) * t57 + (-t105 * t152 + t148) * qJD(1), t102 * t193 + t113 * t104, -t22 * qJD(2) + t195; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t126, 0.2e1 * t137 (-t97 - t98) * t108 (-t105 * t51 + (t43 + t141) * t103) * qJD(1) + t179, 0, 0, 0, 0, 0, -t144 + (-t148 + (-t56 + t152) * t105) * qJD(1), -t143 + (-t146 + (-t57 - t159) * t105) * qJD(1), t113 * t102 - t193 * t104 (t103 * t123 + t105 * t22) * qJD(1) + t110; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t57 * t56, -t55 + t190, t193, t113, t136, t46 * t57 + t9 * t75 + t111, -t46 * t56 + t75 * t8 + t115, -pkin(5) * t34 + t189 * t56, t189 * t7 + (t22 * t57 + t1) * pkin(5); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t55 - t190, -t5 * t57 - t7 * t56 + t13;];
tauc_reg  = t10;
