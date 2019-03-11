% Calculate minimal parameter regressor of coriolis joint torque vector for
% S6RRPRPR2
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% qJD [6x1]
%   Generalized joint velocities
% pkin [10x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d2,d4,d6,theta3]';
% 
% Output:
% tauc_reg [6x30]
%   minimal parameter regressor of coriolis joint torque vector

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-09 10:14
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauc_reg = S6RRPRPR2_coriolisvecJ_fixb_regmin_slag_vp(qJ, qJD, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(10,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRPRPR2_coriolisvecJ_fixb_regmin_slag_vp: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RRPRPR2_coriolisvecJ_fixb_regmin_slag_vp: qJD has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S6RRPRPR2_coriolisvecJ_fixb_regmin_slag_vp: pkin has to be [10x1] (double)');

%% Symbolic Calculation
% From coriolisvec_joint_fixb_regressor_minpar_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 10:13:54
% EndTime: 2019-03-09 10:14:05
% DurationCPUTime: 3.52s
% Computational Cost: add. (5434->326), mult. (14141->415), div. (0->0), fcn. (10896->8), ass. (0->198)
t163 = qJD(2) + qJD(4);
t167 = sin(qJ(6));
t170 = cos(qJ(6));
t169 = sin(qJ(2));
t171 = cos(qJ(2));
t237 = sin(pkin(10));
t238 = cos(pkin(10));
t189 = t237 * t169 - t238 * t171;
t138 = t189 * qJD(1);
t188 = -t238 * t169 - t237 * t171;
t139 = t188 * qJD(1);
t168 = sin(qJ(4));
t263 = cos(qJ(4));
t94 = -t263 * t138 + t139 * t168;
t79 = t163 * t170 - t167 * t94;
t255 = t79 * t94;
t179 = qJD(2) * t139;
t224 = qJD(1) * qJD(2);
t192 = t189 * t224;
t219 = qJD(4) * t263;
t227 = qJD(4) * t168;
t194 = -t138 * t219 + t139 * t227 + t168 * t179 - t263 * t192;
t196 = t168 * t138 + t263 * t139;
t277 = qJD(6) - t196;
t293 = t167 * t277;
t281 = -t170 * t194 + t277 * t293;
t300 = -t255 - t281;
t76 = t163 * t167 + t170 * t94;
t256 = t76 * t94;
t292 = t170 * t277;
t296 = -t167 * t194 - t277 * t292;
t299 = t296 + t256;
t182 = qJD(2) * t188;
t177 = t263 * t182;
t175 = qJD(1) * t177 + t196 * qJD(4) + t168 * t192;
t225 = qJD(6) * t170;
t226 = qJD(6) * t167;
t29 = -t163 * t226 - t167 * t175 - t225 * t94;
t28 = t29 * t170;
t298 = -t293 * t79 + t28;
t294 = t277 * t76;
t49 = t170 * t175;
t30 = t79 * qJD(6) + t49;
t297 = -t292 * t79 + (-t29 + t294) * t167 - t170 * t30;
t265 = t94 * pkin(5);
t261 = t139 * pkin(8);
t251 = -qJ(3) - pkin(7);
t155 = t251 * t171;
t151 = qJD(1) * t155;
t141 = t237 * t151;
t154 = t251 * t169;
t150 = qJD(1) * t154;
t250 = qJD(2) * pkin(2);
t145 = t150 + t250;
t97 = t238 * t145 + t141;
t73 = qJD(2) * pkin(3) + t261 + t97;
t262 = t138 * pkin(8);
t216 = t238 * t151;
t98 = t237 * t145 - t216;
t74 = t98 - t262;
t39 = t168 * t73 + t263 * t74;
t36 = -qJ(5) * t163 - t39;
t20 = -t36 + t265;
t295 = t20 * t277;
t244 = t196 * t163;
t291 = t175 - t244;
t290 = t175 + t244;
t246 = t163 * t94;
t186 = t194 - t246;
t253 = t94 ^ 2;
t286 = t196 ^ 2;
t289 = -t253 + t286;
t288 = qJD(6) - t277;
t285 = pkin(4) * t196;
t284 = pkin(5) * t196;
t283 = t277 * t94;
t243 = qJ(5) * t94;
t282 = t196 * t20;
t252 = t196 * t94;
t222 = -t171 * pkin(2) - pkin(1);
t205 = t222 * qJD(1);
t153 = qJD(3) + t205;
t103 = t138 * pkin(3) + t153;
t178 = qJ(5) * t196 + t103;
t40 = -pkin(4) * t94 + t178;
t258 = t40 * t196;
t221 = pkin(2) * t237;
t157 = t168 * t221;
t220 = t238 * pkin(2);
t159 = t220 + pkin(3);
t101 = -t237 * t150 + t216;
t80 = t101 + t262;
t102 = t238 * t150 + t141;
t81 = t102 + t261;
t280 = qJD(4) * t157 - t159 * t219 + t168 * t80 + t263 * t81;
t279 = t103 * t196;
t266 = pkin(4) + pkin(9);
t278 = t196 * t266;
t38 = t168 * t74 - t263 * t73;
t232 = qJD(5) + t38;
t217 = qJD(2) * t251;
t135 = t171 * qJD(3) + t169 * t217;
t116 = t135 * qJD(1);
t136 = -t169 * qJD(3) + t171 * t217;
t117 = t136 * qJD(1);
t75 = -t237 * t116 + t238 * t117;
t59 = t192 * pkin(8) + t75;
t78 = t238 * t116 + t237 * t117;
t60 = pkin(8) * t179 + t78;
t213 = -t168 * t59 - t73 * t219 + t74 * t227 - t263 * t60;
t12 = -t163 * qJD(5) + t213;
t3 = pkin(5) * t175 - t12;
t233 = t232 - t284;
t18 = -t266 * t163 + t233;
t23 = -t266 * t94 + t178;
t8 = t167 * t18 + t170 * t23;
t276 = t3 * t170 + t8 * t94;
t275 = t40 * t94 - t12;
t274 = -t103 * t94 + t213;
t203 = t167 * t23 - t170 * t18;
t273 = t3 * t167 + t20 * t225 + t203 * t94;
t272 = -0.2e1 * t224;
t241 = -qJD(5) + t280;
t187 = t168 * t159 + t263 * t221;
t239 = t187 * qJD(4) - t168 * t81 + t263 * t80;
t13 = t168 * t60 + t74 * t219 + t73 * t227 - t263 * t59;
t269 = -t239 * t163 - t13;
t133 = -t263 * t159 - pkin(4) + t157;
t129 = -pkin(9) + t133;
t268 = -(t239 - t265) * t277 - t129 * t194;
t180 = t263 * t189;
t99 = -t168 * t188 + t180;
t260 = t20 * t99;
t184 = t168 * t189;
t100 = -t188 * t263 - t184;
t119 = t189 * pkin(3) + t222;
t176 = -t100 * qJ(5) + t119;
t31 = t266 * t99 + t176;
t259 = t31 * t194;
t257 = t194 * t99;
t183 = qJD(2) * t189;
t86 = -t237 * t135 + t238 * t136;
t65 = pkin(8) * t183 + t86;
t87 = t238 * t135 + t237 * t136;
t66 = pkin(8) * t182 + t87;
t105 = t238 * t154 + t237 * t155;
t88 = pkin(8) * t188 + t105;
t106 = t237 * t154 - t238 * t155;
t89 = -t189 * pkin(8) + t106;
t15 = -t168 * t65 - t88 * t219 + t89 * t227 - t263 * t66;
t248 = t15 * t163;
t197 = t168 * t88 + t263 * t89;
t16 = t197 * qJD(4) + t168 * t66 - t263 * t65;
t247 = t16 * t163;
t242 = -t241 - t284;
t174 = qJD(1) ^ 2;
t236 = t171 * t174;
t173 = qJD(2) ^ 2;
t235 = t173 * t169;
t234 = t173 * t171;
t229 = t169 ^ 2 - t171 ^ 2;
t228 = qJD(1) * t169;
t161 = t169 * t250;
t223 = t99 * t225;
t218 = t169 * t224;
t110 = pkin(2) * t228 - pkin(3) * t139;
t208 = pkin(1) * t272;
t45 = t168 * t89 - t263 * t88;
t54 = -qJD(4) * t184 - t168 * t183 - t188 * t219 - t177;
t206 = t277 * t54 + t257;
t204 = t266 * t194 + (t39 + t265) * t277;
t200 = t110 - t243;
t202 = (-qJD(6) * t129 + t200 - t278) * t277;
t201 = (qJD(6) * t266 - t243 - t278) * t277;
t198 = t163 * t39 - t13;
t32 = t100 * pkin(5) + t45;
t195 = -t194 * t32 + t20 * t54 + t3 * t99;
t185 = t194 + t246;
t111 = -pkin(3) * t182 + t161;
t158 = pkin(2) * t218;
t104 = -pkin(3) * t179 + t158;
t53 = t163 * t180 - t168 * t182 - t188 * t227;
t17 = t54 * pkin(4) + t53 * qJ(5) - t100 * qJD(5) + t111;
t14 = -pkin(4) * t175 - qJ(5) * t194 + qJD(5) * t196 + t104;
t132 = qJ(5) + t187;
t52 = -t243 - t285;
t44 = t99 * pkin(4) + t176;
t43 = t200 - t285;
t37 = t194 * t100;
t35 = -pkin(4) * t163 + t232;
t33 = -t99 * pkin(5) + t197;
t11 = t54 * pkin(9) + t17;
t10 = -t53 * pkin(5) + t16;
t9 = -pkin(5) * t54 - t15;
t6 = -pkin(9) * t175 + t14;
t5 = pkin(5) * t194 + t13;
t4 = t170 * t5;
t1 = [0, 0, 0, 0.2e1 * t171 * t218, t229 * t272, t234, -t235, 0, -pkin(7) * t234 + t169 * t208, pkin(7) * t235 + t171 * t208, t105 * t192 + t106 * t179 - t87 * t138 + t86 * t139 + t75 * t188 + t98 * t182 + (qJD(2) * t97 - t78) * t189, t75 * t105 + t78 * t106 + t97 * t86 + t98 * t87 + (t153 + t205) * t161, t196 * t53 + t37, t100 * t175 + t196 * t54 - t53 * t94 - t257, -t53 * t163, -t54 * t163, 0, t103 * t54 + t104 * t99 - t111 * t94 - t119 * t175 - t247, t100 * t104 - t103 * t53 - t111 * t196 + t119 * t194 + t248, t100 * t13 + t12 * t99 - t15 * t94 - t16 * t196 + t175 * t197 + t194 * t45 - t35 * t53 + t36 * t54, -t14 * t99 + t17 * t94 + t175 * t44 - t40 * t54 + t247, -t100 * t14 + t17 * t196 - t194 * t44 + t40 * t53 - t248, -t12 * t197 + t13 * t45 + t14 * t44 + t15 * t36 + t16 * t35 + t17 * t40, t79 * t223 + (t29 * t99 + t54 * t79) * t167 (-t167 * t76 + t170 * t79) * t54 + (-t167 * t30 + t28 + (-t167 * t79 - t170 * t76) * qJD(6)) * t99, t100 * t29 + t206 * t167 + t223 * t277 - t53 * t79, -t226 * t277 * t99 - t100 * t30 + t206 * t170 + t53 * t76, -t277 * t53 + t37, t4 * t100 + t33 * t30 + t203 * t53 + t9 * t76 + (-t6 * t100 - t11 * t277 - t259) * t167 + (t10 * t277 - t195) * t170 + ((-t167 * t32 - t170 * t31) * t277 - t8 * t100 + t167 * t260) * qJD(6), t33 * t29 + t8 * t53 + t9 * t79 + (-(qJD(6) * t32 + t11) * t277 - t259 - (qJD(6) * t18 + t6) * t100 + qJD(6) * t260) * t170 + (-(-qJD(6) * t31 + t10) * t277 - (-qJD(6) * t23 + t5) * t100 + t195) * t167; 0, 0, 0, -t169 * t236, t229 * t174, 0, 0, 0, t174 * pkin(1) * t169, pkin(1) * t236, t179 * t221 + t192 * t220 + (-t101 - t98) * t139 - (-t102 + t97) * t138, -t97 * t101 - t98 * t102 + (-t153 * t228 + t237 * t78 + t238 * t75) * pkin(2), t252, t289, t186, t291, 0, t110 * t94 + t269 + t279, t110 * t196 + t280 * t163 + t274, t132 * t175 + t133 * t194 + (-t241 - t35) * t94 + (-t239 + t36) * t196, -t43 * t94 - t258 - t269, -t241 * t163 - t196 * t43 + t275, -t12 * t132 + t13 * t133 + t239 * t35 + t241 * t36 - t40 * t43, t298, t297, t300, t299, -t283, t132 * t30 + t242 * t76 + t167 * t202 + (-t268 - t282) * t170 + t273, t132 * t29 + t242 * t79 + t170 * t202 + (t268 - t295) * t167 + t276; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t138 ^ 2 - t139 ^ 2, t138 * t98 - t139 * t97 + t158, 0, 0, 0, 0, 0, -t290, t185, -t253 - t286, t290, -t185, t196 * t35 + t36 * t94 + t14, 0, 0, 0, 0, 0, t296 - t256, -t255 + t281; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t252, t289, t186, t291, 0, t198 + t279, -t163 * t38 + t274, -pkin(4) * t194 + qJ(5) * t175 - (-t36 - t39) * t196 - (t35 - t232) * t94, -t52 * t94 - t198 - t258, t232 * t163 - t196 * t52 + t275, -pkin(4) * t13 - qJ(5) * t12 - t232 * t36 - t35 * t39 - t40 * t52, t298, t297, t300, t299, -t283, qJ(5) * t30 + t233 * t76 + t167 * t201 + (-t204 - t282) * t170 + t273, qJ(5) * t29 + t233 * t79 + t170 * t201 + (t204 - t295) * t167 + t276; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t186, -t252, -t163 ^ 2 - t286, t163 * t36 + t13 - t258, 0, 0, 0, 0, 0, -t163 * t76 - t281, -t163 * t79 + t296; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t79 * t76, -t76 ^ 2 + t79 ^ 2, t29 + t294, -t288 * t79 - t49, t194, -t167 * t6 - t20 * t79 - t288 * t8 + t4, -t167 * t5 - t170 * t6 + t20 * t76 + t288 * t203;];
tauc_reg  = t1;
