% Calculate minimal parameter regressor of coriolis joint torque vector for
% S5RRPRR14
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% qJD [5x1]
%   Generalized joint velocities
% pkin [10x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,alpha2,d1,d2,d4,d5,theta3]';
% 
% Output:
% tauc_reg [5x28]
%   minimal parameter regressor of coriolis joint torque vector

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-31 20:40
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauc_reg = S5RRPRR14_coriolisvecJ_fixb_regmin_slag_vp(qJ, qJD, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(10,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RRPRR14_coriolisvecJ_fixb_regmin_slag_vp: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RRPRR14_coriolisvecJ_fixb_regmin_slag_vp: qJD has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S5RRPRR14_coriolisvecJ_fixb_regmin_slag_vp: pkin has to be [10x1] (double)');

%% Symbolic Calculation
% From coriolisvec_joint_fixb_regressor_minpar_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 20:38:38
% EndTime: 2019-12-31 20:38:49
% DurationCPUTime: 3.85s
% Computational Cost: add. (4767->361), mult. (13117->536), div. (0->0), fcn. (10630->10), ass. (0->186)
t171 = sin(qJ(5));
t174 = cos(qJ(5));
t170 = cos(pkin(5));
t238 = qJD(1) * t170;
t159 = qJD(2) + t238;
t167 = sin(pkin(10));
t169 = cos(pkin(10));
t173 = sin(qJ(2));
t168 = sin(pkin(5));
t239 = qJD(1) * t168;
t225 = t173 * t239;
t114 = t159 * t169 - t167 * t225;
t115 = t159 * t167 + t169 * t225;
t172 = sin(qJ(4));
t175 = cos(qJ(4));
t69 = -t175 * t114 + t115 * t172;
t66 = qJD(5) + t69;
t271 = t66 ^ 2;
t136 = t167 * t175 + t169 * t172;
t176 = cos(qJ(2));
t245 = t168 * t176;
t185 = t136 * t245;
t181 = qJD(2) * t185;
t195 = t114 * t172 + t115 * t175;
t264 = qJD(4) * t195;
t37 = qJD(1) * t181 + t264;
t272 = -t271 * t171 + t174 * t37;
t237 = qJD(1) * t176;
t224 = t168 * t237;
t152 = -qJD(4) + t224;
t270 = t69 * t152;
t243 = t175 * t169;
t135 = t167 * t172 - t243;
t184 = t135 * t245;
t103 = qJD(1) * t184;
t131 = t135 * qJD(4);
t242 = t103 - t131;
t241 = -qJD(1) * t185 + t136 * qJD(4);
t49 = t174 * t152 + t171 * t195;
t269 = t195 * t49;
t196 = t152 * t171 - t174 * t195;
t268 = t195 * t196;
t256 = pkin(8) + qJ(3);
t147 = t256 * t167;
t148 = t256 * t169;
t101 = -t147 * t172 + t148 * t175;
t244 = t169 * t176;
t188 = (pkin(3) * t173 - pkin(8) * t244) * t168;
t204 = pkin(2) * t173 - qJ(3) * t176;
t125 = t204 * t239;
t229 = pkin(1) * t238;
t126 = -pkin(7) * t225 + t176 * t229;
t77 = t169 * t125 - t167 * t126;
t57 = qJD(1) * t188 + t77;
t211 = t167 * t224;
t78 = t167 * t125 + t169 * t126;
t63 = -pkin(8) * t211 + t78;
t267 = qJD(3) * t136 + qJD(4) * t101 - t172 * t63 + t175 * t57;
t194 = -t147 * t175 - t148 * t172;
t266 = -qJD(3) * t135 + qJD(4) * t194 - t172 * t57 - t175 * t63;
t265 = t152 * t195;
t127 = pkin(7) * t224 + t173 * t229;
t104 = qJ(3) * t159 + t127;
t122 = (-pkin(2) * t176 - qJ(3) * t173 - pkin(1)) * t168;
t109 = qJD(1) * t122;
t58 = -t104 * t167 + t169 * t109;
t30 = -pkin(3) * t224 - pkin(8) * t115 + t58;
t59 = t169 * t104 + t167 * t109;
t38 = pkin(8) * t114 + t59;
t14 = t172 * t30 + t175 * t38;
t183 = qJD(2) * t188;
t107 = (qJD(2) * t204 - qJD(3) * t173) * t168;
t91 = qJD(1) * t107;
t230 = qJD(1) * qJD(2);
t220 = t168 * t230;
t208 = t173 * t220;
t228 = pkin(1) * qJD(2) * t170;
t212 = qJD(1) * t228;
t189 = -pkin(7) * t208 + t176 * t212;
t92 = qJD(3) * t159 + t189;
t53 = -t167 * t92 + t169 * t91;
t31 = qJD(1) * t183 + t53;
t209 = t176 * t220;
t193 = t167 * t209;
t54 = t167 * t91 + t169 * t92;
t39 = -pkin(8) * t193 + t54;
t179 = -qJD(4) * t14 - t172 * t39 + t175 * t31;
t4 = -pkin(4) * t208 - t179;
t263 = (pkin(4) * t195 + t66 * pkin(9)) * t66 + t4;
t163 = -pkin(3) * t169 - pkin(2);
t85 = pkin(4) * t135 - pkin(9) * t136 + t163;
t95 = pkin(3) * t211 + t127;
t262 = (-t241 * pkin(4) + t242 * pkin(9) + qJD(5) * t101 + t95) * t66 - t85 * t37;
t246 = t168 * t173;
t130 = t167 * t170 + t169 * t246;
t260 = pkin(1) * t173;
t121 = pkin(7) * t245 + (qJ(3) + t260) * t170;
t73 = -t121 * t167 + t169 * t122;
t44 = -pkin(3) * t245 - pkin(8) * t130 + t73;
t129 = t167 * t246 - t170 * t169;
t74 = t169 * t121 + t167 * t122;
t56 = -pkin(8) * t129 + t74;
t198 = t172 * t44 + t175 * t56;
t236 = qJD(2) * t173;
t223 = t168 * t236;
t192 = -pkin(7) * t223 + t176 * t228;
t113 = qJD(3) * t170 + t192;
t61 = t169 * t107 - t167 * t113;
t42 = t183 + t61;
t222 = qJD(2) * t245;
t210 = t167 * t222;
t62 = t167 * t107 + t169 * t113;
t52 = -pkin(8) * t210 + t62;
t261 = -qJD(4) * t198 - t172 * t52 + t175 * t42;
t234 = qJD(4) * t175;
t36 = t114 * t234 + t209 * t243 + (-qJD(4) * t115 - t193) * t172;
t17 = -qJD(5) * t196 + t171 * t36 - t174 * t208;
t120 = pkin(7) * t209 + t173 * t212;
t88 = pkin(3) * t193 + t120;
t10 = pkin(4) * t37 - pkin(9) * t36 + t88;
t12 = -pkin(9) * t152 + t14;
t97 = -pkin(2) * t159 + qJD(3) - t126;
t67 = -pkin(3) * t114 + t97;
t18 = pkin(4) * t69 - pkin(9) * t195 + t67;
t203 = t12 * t171 - t174 * t18;
t235 = qJD(4) * t172;
t191 = t172 * t31 + t175 * t39 + t30 * t234 - t235 * t38;
t3 = pkin(9) * t208 + t191;
t1 = -qJD(5) * t203 + t171 * t10 + t174 * t3;
t259 = t49 * t66;
t258 = t196 * t66;
t254 = pkin(4) * t225 + t267;
t232 = qJD(5) * t174;
t233 = qJD(5) * t171;
t16 = -t152 * t232 + t171 * t208 + t174 * t36 - t195 * t233;
t253 = t16 * t171;
t252 = t171 * t37;
t251 = t173 * t58;
t250 = t173 * t59;
t249 = t136 * t174;
t164 = t168 ^ 2;
t248 = t164 * qJD(1) ^ 2;
t247 = t167 * t176;
t128 = pkin(7) * t222 + t173 * t228;
t240 = t173 ^ 2 - t176 ^ 2;
t231 = qJD(2) - t159;
t227 = t176 * t248;
t226 = t171 * t245;
t96 = pkin(3) * t210 + t128;
t221 = t164 * t230;
t215 = t174 * t66;
t214 = t159 + t238;
t213 = 0.2e1 * t221;
t206 = -0.2e1 * pkin(1) * t221;
t6 = t12 * t174 + t171 * t18;
t20 = -pkin(9) * t245 + t198;
t81 = t175 * t129 + t130 * t172;
t82 = -t129 * t172 + t130 * t175;
t124 = pkin(7) * t246 + (-pkin(1) * t176 - pkin(2)) * t170;
t83 = t129 * pkin(3) + t124;
t25 = t81 * pkin(4) - t82 * pkin(9) + t83;
t202 = t171 * t25 + t174 * t20;
t201 = -t171 * t20 + t174 * t25;
t13 = -t172 * t38 + t175 * t30;
t199 = -t172 * t56 + t175 * t44;
t64 = t171 * t82 + t174 * t245;
t190 = t172 * t42 + t175 * t52 + t44 * t234 - t235 * t56;
t79 = -t103 * t171 - t174 * t225;
t187 = -t131 * t171 + t136 * t232 - t79;
t80 = -t103 * t174 + t171 * t225;
t186 = -t131 * t174 - t136 * t233 - t80;
t11 = pkin(4) * t152 - t13;
t182 = -pkin(9) * t37 + (t11 + t13) * t66;
t2 = -qJD(5) * t6 + t174 * t10 - t171 * t3;
t180 = -qJ(3) * t236 + (-pkin(2) * qJD(2) + qJD(3) - t97) * t176;
t178 = -t101 * t37 + t4 * t136 + (pkin(9) * t225 - qJD(5) * t85 - t266) * t66;
t65 = t174 * t82 - t226;
t47 = qJD(4) * t82 + t181;
t46 = -qJD(2) * t184 - qJD(4) * t81;
t22 = -qJD(5) * t226 + t171 * t46 - t174 * t223 + t232 * t82;
t21 = -qJD(5) * t64 + t171 * t223 + t174 * t46;
t19 = pkin(4) * t245 - t199;
t15 = pkin(4) * t47 - pkin(9) * t46 + t96;
t8 = -pkin(4) * t223 - t261;
t7 = pkin(9) * t223 + t190;
t5 = [0, 0, 0, t173 * t176 * t213, -t240 * t213, t214 * t222, -t214 * t223, 0, -t120 * t170 - t128 * t159 + t173 * t206, -t159 * t192 - t170 * t189 + t176 * t206, -t114 * t128 + t120 * t129 + ((-qJD(1) * t61 - t53) * t176 + (t97 * t247 + t251 + (t124 * t247 + t173 * t73) * qJD(1)) * qJD(2)) * t168, t115 * t128 + t120 * t130 + ((qJD(1) * t62 + t54) * t176 + (t97 * t244 - t250 + (t124 * t244 - t173 * t74) * qJD(1)) * qJD(2)) * t168, t114 * t62 - t115 * t61 - t129 * t54 - t130 * t53 + (-t167 * t59 - t169 * t58 + (-t167 * t74 - t169 * t73) * qJD(1)) * t222, t120 * t124 + t128 * t97 + t53 * t73 + t54 * t74 + t58 * t61 + t59 * t62, t195 * t46 + t36 * t82, -t195 * t47 - t36 * t81 - t37 * t82 - t46 * t69, -t46 * t152 + (-t176 * t36 + (qJD(1) * t82 + t195) * t236) * t168, t47 * t152 + (t176 * t37 + (-qJD(1) * t81 - t69) * t236) * t168, (-t152 * t168 - t164 * t237) * t236, -t261 * t152 + t96 * t69 + t83 * t37 + t88 * t81 + t67 * t47 + (-t179 * t176 + (qJD(1) * t199 + t13) * t236) * t168, t190 * t152 + t96 * t195 + t83 * t36 + t88 * t82 + t67 * t46 + (t191 * t176 + (-t198 * qJD(1) - t14) * t236) * t168, t16 * t65 - t196 * t21, -t16 * t64 - t17 * t65 + t196 * t22 - t21 * t49, t16 * t81 - t196 * t47 + t21 * t66 + t37 * t65, -t17 * t81 - t22 * t66 - t37 * t64 - t47 * t49, t37 * t81 + t47 * t66, (-qJD(5) * t202 + t174 * t15 - t171 * t7) * t66 + t201 * t37 + t2 * t81 - t203 * t47 + t8 * t49 + t19 * t17 + t4 * t64 + t11 * t22, -(qJD(5) * t201 + t171 * t15 + t174 * t7) * t66 - t202 * t37 - t1 * t81 - t6 * t47 - t8 * t196 + t19 * t16 + t4 * t65 + t11 * t21; 0, 0, 0, -t173 * t227, t240 * t248, t231 * t224, -t231 * t225, 0, t127 * t159 + t248 * t260 - t120, pkin(1) * t227 + t126 * t159 - t189, t127 * t114 - t120 * t169 + (t167 * t180 + t176 * t77 - t251) * t239, -t115 * t127 + t120 * t167 + (t169 * t180 - t176 * t78 + t250) * t239, -t78 * t114 + t77 * t115 + (qJD(3) * t114 + t224 * t58 + t54) * t169 + (qJD(3) * t115 + t224 * t59 - t53) * t167, -t120 * pkin(2) - t97 * t127 - t58 * t77 - t59 * t78 + (-t167 * t58 + t169 * t59) * qJD(3) + (-t53 * t167 + t54 * t169) * qJ(3), t36 * t136 + t195 * t242, -t36 * t135 - t136 * t37 - t195 * t241 - t242 * t69, -t242 * t152 + (qJD(2) * t136 - t195) * t225, t241 * t152 + (-qJD(2) * t135 + t69) * t225, t152 * t225, t88 * t135 + t163 * t37 - t95 * t69 + t241 * t67 + t267 * t152 + (qJD(2) * t194 - t13) * t225, t88 * t136 + t163 * t36 - t95 * t195 + t242 * t67 + t266 * t152 + (-qJD(2) * t101 + t14) * t225, t16 * t249 - t186 * t196, t80 * t49 - t196 * t79 - (t171 * t196 - t174 * t49) * t131 + (-t253 - t17 * t174 + (t171 * t49 + t174 * t196) * qJD(5)) * t136, t16 * t135 + t186 * t66 - t196 * t241 + t249 * t37, -t17 * t135 - t136 * t252 - t187 * t66 - t241 * t49, t37 * t135 + t241 * t66, t187 * t11 + t2 * t135 - t17 * t194 + t178 * t171 - t262 * t174 - t203 * t241 + t254 * t49, -t1 * t135 + t186 * t11 - t16 * t194 + t262 * t171 + t178 * t174 - t196 * t254 - t241 * t6; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, (qJD(2) * t167 - t115) * t224, (qJD(2) * t169 - t114) * t224, -t114 ^ 2 - t115 ^ 2, -t114 * t59 + t115 * t58 + t120, 0, 0, 0, 0, 0, t37 - t265, t36 + t270, 0, 0, 0, 0, 0, -t269 + t272, -t174 * t271 - t252 + t268; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t195 * t69, t195 ^ 2 - t69 ^ 2, t36 - t270, -t136 * t209 - t264 - t265, t208, -t14 * t152 - t195 * t67 + t179, -t13 * t152 + t67 * t69 - t191, -t196 * t215 + t253, (t16 - t259) * t174 + (-t17 + t258) * t171, t215 * t66 + t252 + t268, t269 + t272, -t66 * t195, -pkin(4) * t17 - t14 * t49 + t182 * t171 - t263 * t174 + t195 * t203, -pkin(4) * t16 + t14 * t196 + t263 * t171 + t182 * t174 + t6 * t195; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t196 * t49, t196 ^ 2 - t49 ^ 2, t16 + t259, -t17 - t258, t37, t11 * t196 + t6 * t66 + t2, t11 * t49 - t203 * t66 - t1;];
tauc_reg = t5;
