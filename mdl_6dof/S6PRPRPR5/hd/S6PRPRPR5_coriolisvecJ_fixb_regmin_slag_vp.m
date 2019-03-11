% Calculate minimal parameter regressor of coriolis joint torque vector for
% S6PRPRPR5
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% qJD [6x1]
%   Generalized joint velocities
% pkin [11x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,alpha2,d2,d4,d6,theta1,theta3]';
% 
% Output:
% tauc_reg [6x26]
%   minimal parameter regressor of coriolis joint torque vector

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-08 19:45
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauc_reg = S6PRPRPR5_coriolisvecJ_fixb_regmin_slag_vp(qJ, qJD, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(11,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6PRPRPR5_coriolisvecJ_fixb_regmin_slag_vp: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6PRPRPR5_coriolisvecJ_fixb_regmin_slag_vp: qJD has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [11 1]), ...
  'S6PRPRPR5_coriolisvecJ_fixb_regmin_slag_vp: pkin has to be [11x1] (double)');

%% Symbolic Calculation
% From coriolisvec_joint_fixb_regressor_minpar_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-08 19:45:03
% EndTime: 2019-03-08 19:45:09
% DurationCPUTime: 1.95s
% Computational Cost: add. (1884->264), mult. (5035->359), div. (0->0), fcn. (4005->10), ass. (0->147)
t104 = sin(pkin(11));
t106 = cos(pkin(11));
t109 = sin(qJ(4));
t192 = cos(qJ(4));
t86 = t192 * t104 + t109 * t106;
t117 = qJD(2) * t86;
t207 = qJD(6) + t117;
t108 = sin(qJ(6));
t111 = cos(qJ(6));
t151 = t192 * t106;
t140 = qJD(2) * t151;
t173 = t104 * t109;
t149 = qJD(2) * t173;
t75 = -t140 + t149;
t53 = qJD(4) * t108 - t111 * t75;
t211 = t207 * t53;
t210 = qJD(6) - t207;
t141 = t108 * t207;
t158 = qJD(4) * t109;
t147 = t104 * t158;
t94 = qJD(4) * t140;
t68 = qJD(2) * t147 - t94;
t128 = -t111 * t68 - t141 * t207;
t209 = t151 - t173;
t208 = t117 * qJD(4);
t177 = pkin(8) * qJD(2);
t110 = sin(qJ(2));
t105 = sin(pkin(6));
t162 = qJD(1) * t105;
t150 = t110 * t162;
t90 = qJD(2) * qJ(3) + t150;
t107 = cos(pkin(6));
t161 = qJD(1) * t107;
t96 = t106 * t161;
t47 = t96 + (-t90 - t177) * t104;
t57 = t104 * t161 + t106 * t90;
t48 = t106 * t177 + t57;
t179 = -t109 * t48 + t192 * t47;
t201 = qJD(5) - t179;
t112 = cos(qJ(2));
t171 = t105 * t112;
t116 = t86 * t171;
t184 = pkin(8) + qJ(3);
t91 = t184 * t104;
t92 = t184 * t106;
t126 = -t109 * t91 + t192 * t92;
t180 = -qJD(1) * t116 + t86 * qJD(3) + t126 * qJD(4);
t145 = qJD(4) * t192;
t79 = -t106 * t145 + t147;
t204 = -t79 * qJ(5) + t86 * qJD(5) + t150;
t203 = t180 * qJD(4);
t131 = t151 * t171;
t146 = t112 * t162;
t181 = qJD(1) * t131 + (qJD(3) * t104 + qJD(4) * t92) * t109 - t146 * t173 - qJD(3) * t151 + t91 * t145;
t202 = t181 * qJD(4);
t134 = qJD(3) - t146;
t21 = t109 * t47 + t192 * t48;
t15 = -qJD(4) * qJ(5) - t21;
t194 = pkin(5) * t75;
t10 = -t15 - t194;
t195 = pkin(4) + pkin(9);
t200 = t195 * t68 + (t10 - t21 + t194) * t207;
t159 = qJD(2) * t110;
t172 = t105 * t110;
t73 = -t104 * t172 + t106 * t107;
t74 = t104 * t107 + t106 * t172;
t35 = t109 * t73 + t192 * t74;
t17 = qJD(2) * t116 + t35 * qJD(4);
t80 = t86 * qJD(4);
t69 = qJD(2) * t80;
t199 = t105 * (-t112 * t69 + t75 * t159) - t17 * qJD(4);
t160 = qJD(2) * t105;
t16 = -t73 * t145 - qJD(2) * t131 + (t104 * t112 * t160 + qJD(4) * t74) * t109;
t198 = t105 * (-t112 * t68 - t117 * t159) - t16 * qJD(4);
t197 = t75 ^ 2;
t196 = t117 ^ 2;
t193 = t69 * pkin(4);
t191 = t10 * t209;
t101 = -t106 * pkin(3) - pkin(2);
t70 = t101 * qJD(2) + t134;
t115 = -qJ(5) * t117 + t70;
t24 = t75 * pkin(4) + t115;
t190 = t24 * t117;
t129 = -t86 * qJ(5) + t101;
t27 = -t195 * t209 + t129;
t189 = t27 * t68;
t55 = qJD(4) * t111 + t108 * t75;
t188 = t55 * t75;
t187 = t68 * t209;
t186 = t75 * t53;
t185 = t117 * t75;
t183 = -t80 * pkin(5) - t181;
t182 = -t79 * pkin(5) + t180;
t176 = qJD(2) * pkin(2);
t156 = qJD(6) * t111;
t157 = qJD(6) * t108;
t28 = -qJD(4) * t157 + t108 * t69 + t75 * t156;
t175 = t28 * t111;
t174 = t75 * qJ(5);
t114 = qJD(2) ^ 2;
t170 = t105 * t114;
t165 = pkin(5) * t117 + t201;
t163 = t104 ^ 2 + t106 ^ 2;
t155 = t209 * t156;
t154 = t110 * t170;
t153 = t112 * t170;
t148 = t105 * t159;
t93 = qJD(1) * t148;
t144 = t68 * qJ(5) + t93;
t84 = (qJD(3) + t146) * qJD(2);
t143 = t163 * t84;
t142 = t47 * t145 - t48 * t158 + t209 * t84;
t7 = t48 * t145 + t47 * t158 + t86 * t84;
t49 = t109 * t92 + t192 * t91;
t139 = -t195 * t80 + t204;
t138 = -pkin(4) * t80 + t204;
t137 = t207 * t80 + t187;
t13 = t195 * t75 + t115;
t8 = -t195 * qJD(4) + t165;
t2 = t108 * t8 + t111 * t13;
t136 = t108 * t13 - t111 * t8;
t133 = t104 * (-t104 * t90 + t96) - t106 * t57;
t127 = -qJD(5) * t117 + t144;
t34 = t109 * t74 - t192 * t73;
t123 = -t108 * t34 + t111 * t171;
t122 = t108 * t171 + t111 * t34;
t121 = qJD(4) * t21 - t7;
t6 = -qJD(4) * qJD(5) - t142;
t3 = -pkin(5) * t69 - t6;
t120 = t3 + (t195 * t207 + t174) * t207;
t32 = t86 * pkin(5) + t49;
t119 = -t10 * t80 + t209 * t3 - t32 * t68;
t118 = -t111 * t207 ^ 2 + t108 * t68;
t89 = t134 - t176;
t71 = qJD(4) * t75;
t60 = t111 * t69;
t46 = t68 * t86;
t37 = -pkin(4) * t209 + t129;
t36 = pkin(4) * t117 + t174;
t33 = pkin(5) * t209 + t126;
t29 = t55 * qJD(6) - t60;
t19 = t127 + t193;
t14 = -qJD(4) * pkin(4) + t201;
t9 = t195 * t69 + t127;
t5 = -pkin(5) * t68 + t7;
t4 = t111 * t5;
t1 = [0, 0, -t154, -t153, -t106 * t154, t104 * t154, t163 * t153 (-t104 * t73 + t106 * t74) * t84 + (t110 * t89 + (-t133 - t150) * t112) * t160, 0, 0, 0, 0, 0, t199, -t198, t117 * t17 + t16 * t75 - t34 * t68 - t35 * t69, -t199, t198, t14 * t17 + t15 * t16 + t34 * t7 - t35 * t6 + (-t112 * t19 + t159 * t24) * t105, 0, 0, 0, 0, 0 (qJD(6) * t123 - t108 * t148 + t111 * t17) * t207 - t122 * t68 - t16 * t53 + t35 * t29 -(qJD(6) * t122 + t108 * t17 + t111 * t148) * t207 - t123 * t68 - t16 * t55 + t35 * t28; 0, 0, 0, 0, 0, 0, qJD(2) * t134 * t163 + t143, -t133 * qJD(3) + qJ(3) * t143 + (t133 * t112 + (-t89 - t176) * t110) * t162, -t117 * t79 - t46, -t117 * t80 - t69 * t86 + t75 * t79 - t187, -t79 * qJD(4), -t80 * qJD(4), 0, t101 * t69 + t70 * t80 - t203 + (-qJD(2) * t209 - t75) * t150, -t101 * t68 - t70 * t79 + t202, t117 * t180 - t126 * t69 - t14 * t79 + t15 * t80 + t181 * t75 - t209 * t6 - t49 * t68 + t7 * t86, t138 * t75 + t19 * t209 - t24 * t80 - t37 * t69 + t203, t117 * t138 - t19 * t86 + t24 * t79 + t37 * t68 - t202, -t126 * t6 - t138 * t24 + t180 * t14 + t181 * t15 + t19 * t37 + t49 * t7, -t55 * t155 + (-t209 * t28 + t55 * t80) * t108 (-t108 * t53 + t111 * t55) * t80 - (-t108 * t29 + t175 + (-t108 * t55 - t111 * t53) * qJD(6)) * t209, t108 * t137 - t155 * t207 + t28 * t86 - t55 * t79, t157 * t207 * t209 + t111 * t137 - t29 * t86 + t53 * t79, -t207 * t79 - t46, t136 * t79 + t33 * t29 + t4 * t86 + t183 * t53 + (t139 * t207 - t9 * t86 + t189) * t108 + (t182 * t207 + t119) * t111 + ((-t108 * t32 - t111 * t27) * t207 - t2 * t86 - t108 * t191) * qJD(6), t2 * t79 + t33 * t28 + t183 * t55 + (t189 - (qJD(6) * t8 + t9) * t86 - qJD(6) * t191 + (-qJD(6) * t32 + t139) * t207) * t111 + (-(-qJD(6) * t13 + t5) * t86 + (qJD(6) * t27 - t182) * t207 - t119) * t108; 0, 0, 0, 0, 0, 0, -t163 * t114, t133 * qJD(2) + t93, 0, 0, 0, 0, 0, 0.2e1 * t208, t94 + (-t75 - t149) * qJD(4), -t196 - t197, -0.2e1 * t208, t68 + t71, t193 - t15 * t75 + (-qJD(5) - t14) * t117 + t144, 0, 0, 0, 0, 0, t118 + t186, t188 - t128; 0, 0, 0, 0, 0, 0, 0, 0, t185, t196 - t197, t94 + (t75 - t149) * qJD(4), 0, 0, -t117 * t70 + t121, qJD(4) * t179 + t70 * t75 - t142, pkin(4) * t68 - qJ(5) * t69 + (-t15 - t21) * t117 + (t14 - t201) * t75, t36 * t75 - t121 + t190, -t24 * t75 + t36 * t117 + (0.2e1 * qJD(5) - t179) * qJD(4) + t142, -t7 * pkin(4) - t6 * qJ(5) - t14 * t21 - t15 * t201 - t24 * t36, -t141 * t55 + t175 (-t207 * t55 - t29) * t111 + (-t28 + t211) * t108, t128 + t188, t118 - t186, t207 * t75, qJ(5) * t29 + t120 * t108 + t111 * t200 - t136 * t75 + t165 * t53, qJ(5) * t28 - t108 * t200 + t120 * t111 + t165 * t55 - t2 * t75; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t68 + t71, -t185, -qJD(4) ^ 2 - t196, qJD(4) * t15 + t190 + t7, 0, 0, 0, 0, 0, -qJD(4) * t53 + t128, -qJD(4) * t55 + t118; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t55 * t53, -t53 ^ 2 + t55 ^ 2, t28 + t211, -t210 * t55 + t60, -t68, -t10 * t55 - t108 * t9 - t2 * t210 + t4, t10 * t53 - t108 * t5 - t111 * t9 + t136 * t210;];
tauc_reg  = t1;
