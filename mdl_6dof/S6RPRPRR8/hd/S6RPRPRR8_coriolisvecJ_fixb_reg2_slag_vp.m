% Calculate inertial parameters regressor of coriolis joint torque vector for
% S6RPRPRR8
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% qJD [6x1]
%   Generalized joint velocities
% pkin [10x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d3,d5,d6,theta4]';
% 
% Output:
% tauc_reg [6x(6*10)]
%   inertial parameter regressor of coriolis joint torque vector

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-09 04:00
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauc_reg = S6RPRPRR8_coriolisvecJ_fixb_reg2_slag_vp(qJ, qJD, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(10,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPRPRR8_coriolisvecJ_fixb_reg2_slag_vp: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RPRPRR8_coriolisvecJ_fixb_reg2_slag_vp: qJD has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S6RPRPRR8_coriolisvecJ_fixb_reg2_slag_vp: pkin has to be [10x1] (double)');

%% Symbolic Calculation
% From coriolisvec_joint_fixb_regressor_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 03:59:41
% EndTime: 2019-03-09 03:59:48
% DurationCPUTime: 4.45s
% Computational Cost: add. (9848->460), mult. (21284->616), div. (0->0), fcn. (14947->8), ass. (0->225)
t189 = sin(qJ(5));
t192 = cos(qJ(5));
t249 = qJD(5) * t192;
t250 = qJD(5) * t189;
t187 = sin(pkin(10));
t194 = -pkin(1) - pkin(7);
t164 = t194 * qJD(1) + qJD(2);
t193 = cos(qJ(3));
t244 = t193 * qJD(4);
t190 = sin(qJ(3));
t252 = qJD(3) * t190;
t205 = qJ(4) * t252 - t244;
t197 = qJD(1) * t205 - t164 * t252;
t245 = t190 * qJD(4);
t251 = qJD(3) * t193;
t112 = t164 * t251 + (-qJ(4) * t251 - t245) * qJD(1);
t282 = cos(pkin(10));
t227 = t282 * t112;
t52 = t187 * t197 + t227;
t224 = qJD(3) * t282;
t215 = qJD(1) * t224;
t163 = t193 * t215;
t255 = qJD(1) * t190;
t234 = t187 * t255;
t130 = qJD(3) * t234 - t163;
t242 = qJD(1) * qJD(3);
t232 = t193 * t242;
t131 = t187 * t232 + t190 * t215;
t183 = qJD(1) * qJD(2);
t156 = pkin(3) * t232 + t183;
t64 = -pkin(4) * t130 + pkin(8) * t131 + t156;
t254 = qJD(1) * t193;
t135 = -qJ(4) * t254 + t193 * t164;
t125 = qJD(3) * pkin(3) + t135;
t134 = (-qJ(4) * qJD(1) + t164) * t190;
t226 = t282 * t134;
t75 = t187 * t125 + t226;
t67 = qJD(3) * pkin(8) + t75;
t204 = -t187 * t193 - t282 * t190;
t140 = t204 * qJD(1);
t225 = t282 * t193;
t143 = qJD(1) * t225 - t234;
t157 = pkin(3) * t255 + qJD(1) * qJ(2) + qJD(4);
t76 = -pkin(4) * t140 - pkin(8) * t143 + t157;
t15 = t189 * t64 + t192 * t52 + t76 * t249 - t67 * t250;
t315 = qJD(5) - t140;
t37 = -t189 * t67 + t192 * t76;
t322 = -t315 * t37 + t15;
t38 = t189 * t76 + t192 * t67;
t16 = -t38 * qJD(5) - t189 * t52 + t192 * t64;
t321 = t315 * t38 + t16;
t113 = -t192 * qJD(3) + t143 * t189;
t115 = qJD(3) * t189 + t143 * t192;
t188 = sin(qJ(6));
t191 = cos(qJ(6));
t212 = t113 * t188 - t191 * t115;
t54 = t191 * t113 + t115 * t188;
t301 = t54 * t212;
t240 = qJD(5) + qJD(6);
t247 = qJD(6) * t191;
t265 = t188 * t189;
t306 = t191 * t192 - t265;
t287 = t140 * t306 - t191 * t249 - t192 * t247 + t240 * t265;
t264 = t188 * t192;
t154 = t189 * t191 + t264;
t105 = t240 * t154;
t286 = -t154 * t140 + t105;
t172 = pkin(3) * t187 + pkin(8);
t300 = pkin(9) + t172;
t228 = qJD(5) * t300;
t273 = t140 * t189;
t122 = t187 * t134;
t86 = t282 * t135 - t122;
t238 = pkin(3) * t254;
t88 = pkin(4) * t143 - pkin(8) * t140 + t238;
t40 = t189 * t88 + t192 * t86;
t320 = -pkin(9) * t273 + t189 * t228 + t40;
t39 = -t189 * t86 + t192 * t88;
t319 = pkin(5) * t143 + t39 + (-pkin(9) * t140 + t228) * t192;
t318 = -qJD(3) * qJD(5) + t131;
t317 = t250 - t273;
t222 = t192 * t315;
t263 = t189 * t130;
t316 = t222 * t315 - t263;
t314 = t212 ^ 2 - t54 ^ 2;
t68 = t143 * t250 + t318 * t192;
t6 = -pkin(5) * t130 + pkin(9) * t68 + t16;
t28 = -pkin(9) * t115 + t37;
t25 = pkin(5) * t315 + t28;
t29 = -pkin(9) * t113 + t38;
t291 = t191 * t29;
t8 = t188 * t25 + t291;
t236 = t143 * t249 - t318 * t189;
t9 = -t236 * pkin(9) + t15;
t2 = -t8 * qJD(6) - t188 * t9 + t191 * t6;
t74 = t282 * t125 - t122;
t66 = -qJD(3) * pkin(4) - t74;
t46 = t113 * pkin(5) + t66;
t313 = t46 * t212 + t2;
t129 = qJD(6) + t315;
t248 = qJD(6) * t188;
t23 = t113 * t247 + t115 * t248 + t188 * t236 + t191 * t68;
t312 = t129 * t54 - t23;
t1 = (qJD(6) * t25 + t9) * t191 + t188 * t6 - t29 * t248;
t311 = t46 * t54 - t1;
t199 = t212 * qJD(6) + t188 * t68 - t191 * t236;
t310 = -t129 * t212 + t199;
t258 = qJ(4) - t194;
t158 = t258 * t190;
t159 = t258 * t193;
t107 = -t282 * t158 - t187 * t159;
t100 = t192 * t107;
t266 = t187 * t190;
t152 = t225 - t266;
t174 = t190 * pkin(3) + qJ(2);
t99 = -pkin(4) * t204 - pkin(8) * t152 + t174;
t48 = t189 * t99 + t100;
t308 = t306 * t204;
t307 = t154 * t204;
t305 = -t212 * t286 + t23 * t306;
t304 = -t287 * t129 - t154 * t130;
t303 = t143 ^ 2;
t302 = 0.2e1 * t183;
t149 = t300 * t189;
t150 = t300 * t192;
t98 = -t149 * t188 + t150 * t191;
t299 = qJD(6) * t98 - t320 * t188 + t319 * t191;
t97 = -t149 * t191 - t150 * t188;
t298 = -qJD(6) * t97 + t319 * t188 + t320 * t191;
t106 = -t158 * t187 + t282 * t159;
t51 = t187 * t112 - t282 * t197;
t297 = t106 * t51;
t296 = t143 * t54;
t295 = t152 * t51;
t293 = t188 * t29;
t292 = t189 * t51;
t290 = t51 * t192;
t289 = t212 * t143;
t288 = t68 * t189;
t142 = t187 * t252 - t193 * t224;
t285 = qJD(1) * t306 - t154 * t142 - t240 * t308;
t284 = t154 * qJD(1) + t142 * t306 - t240 * t307;
t61 = t189 * t236;
t283 = -t113 * t249 - t61;
t281 = t113 * t140;
t280 = t113 * t143;
t279 = t113 * t189;
t278 = t113 * t192;
t277 = t115 * t113;
t276 = t115 * t143;
t275 = t115 * t189;
t274 = t115 * t192;
t272 = t143 * t140;
t145 = t204 * qJD(3);
t271 = t145 * t189;
t270 = t145 * t192;
t101 = t204 * t130;
t269 = t152 * t189;
t268 = t152 * t192;
t119 = t192 * t130;
t195 = qJD(3) ^ 2;
t261 = t195 * t190;
t260 = t195 * t193;
t196 = qJD(1) ^ 2;
t259 = t196 * qJ(2);
t257 = t190 ^ 2 - t193 ^ 2;
t256 = -t195 - t196;
t253 = qJD(3) * t140;
t246 = t142 * qJD(3);
t165 = pkin(3) * t251 + qJD(2);
t239 = 0.2e1 * qJD(1);
t237 = t193 * t196 * t190;
t235 = t152 * t250;
t132 = -qJD(3) * t159 - t245;
t202 = t258 * t252 - t244;
t80 = t282 * t132 + t187 * t202;
t87 = -pkin(4) * t142 - pkin(8) * t145 + t165;
t229 = -t189 * t80 + t192 * t87;
t47 = -t107 * t189 + t192 * t99;
t79 = t187 * t132 - t282 * t202;
t85 = t135 * t187 + t226;
t223 = t189 * t315;
t221 = -qJD(5) * t204 + qJD(1);
t220 = t154 * t199 + t287 * t54;
t219 = t190 * t232;
t218 = t317 * pkin(5) - t85;
t217 = -t286 * t129 - t130 * t306;
t173 = -t282 * pkin(3) - pkin(4);
t216 = t236 * t192;
t36 = -pkin(5) * t204 - pkin(9) * t268 + t47;
t41 = -pkin(9) * t269 + t48;
t18 = -t188 * t41 + t191 * t36;
t19 = t188 * t36 + t191 * t41;
t214 = t189 * t38 + t192 * t37;
t213 = t189 * t37 - t192 * t38;
t211 = t274 + t279;
t210 = t140 * t142 + t101;
t209 = -t131 * t152 + t143 * t145;
t208 = -t315 * t317 - t119;
t207 = t152 * t249 + t271;
t206 = -t235 + t270;
t21 = -t107 * t250 + t189 * t87 + t192 * t80 + t99 * t249;
t203 = t172 * t130 + t315 * t66;
t200 = -t142 * t75 + t145 * t74 - t204 * t52 - t295;
t198 = -t214 * qJD(5) + t15 * t192 - t16 * t189;
t180 = qJ(2) * t302;
t160 = -t192 * pkin(5) + t173;
t139 = t140 ^ 2;
t133 = t145 * qJD(3);
t93 = t306 * t152;
t91 = t154 * t152;
t72 = pkin(5) * t269 + t106;
t43 = pkin(5) * t207 + t79;
t34 = t236 * pkin(5) + t51;
t33 = t145 * t264 - t188 * t235 - t248 * t269 + (t240 * t268 + t271) * t191;
t31 = t105 * t152 - t306 * t145;
t22 = -t48 * qJD(5) + t229;
t17 = -pkin(9) * t207 + t21;
t14 = -pkin(9) * t270 - pkin(5) * t142 + (-t100 + (pkin(9) * t152 - t99) * t189) * qJD(5) + t229;
t11 = t191 * t28 - t293;
t10 = -t188 * t28 - t291;
t7 = t191 * t25 - t293;
t4 = -t19 * qJD(6) + t191 * t14 - t188 * t17;
t3 = t18 * qJD(6) + t188 * t14 + t191 * t17;
t5 = [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t302, t180, -0.2e1 * t219, 0.2e1 * t257 * t242, -t261, 0.2e1 * t219, -t260, 0, -t194 * t261 + (qJ(2) * t251 + qJD(2) * t190) * t239, -t194 * t260 + (-qJ(2) * t252 + qJD(2) * t193) * t239, 0, t180, t209, t130 * t152 - t131 * t204 + t140 * t145 + t142 * t143, t133, t210, t246, 0, -qJD(3) * t79 - t130 * t174 - t140 * t165 - t142 * t157 - t156 * t204, -qJD(3) * t80 - t131 * t174 + t143 * t165 + t145 * t157 + t152 * t156, -t106 * t131 + t107 * t130 + t140 * t80 + t143 * t79 - t200, t107 * t52 + t156 * t174 + t157 * t165 - t74 * t79 + t75 * t80 + t297, t115 * t206 - t68 * t268 -(t275 + t278) * t145 + (-t216 + t288 + (-t274 + t279) * qJD(5)) * t152, -t115 * t142 - t119 * t152 + t204 * t68 + t206 * t315, t113 * t207 + t152 * t61, t113 * t142 + t152 * t263 + t204 * t236 - t207 * t315, -t142 * t315 + t101, t22 * t315 - t47 * t130 - t16 * t204 - t37 * t142 + t79 * t113 + t106 * t236 + t66 * t271 + (t249 * t66 + t292) * t152, t66 * t270 - t106 * t68 + t115 * t79 + t130 * t48 - t315 * t21 + t142 * t38 + t15 * t204 + (-t250 * t66 + t290) * t152, -t21 * t113 - t48 * t236 - t22 * t115 + t47 * t68 - t214 * t145 + (qJD(5) * t213 - t15 * t189 - t16 * t192) * t152, t15 * t48 + t16 * t47 + t21 * t38 + t22 * t37 + t66 * t79 + t297, t212 * t31 - t23 * t93, t199 * t93 + t212 * t33 + t23 * t91 + t31 * t54, -t129 * t31 - t130 * t93 + t142 * t212 + t204 * t23, -t199 * t91 + t33 * t54, -t129 * t33 + t130 * t91 + t142 * t54 - t199 * t204, -t129 * t142 + t101, t129 * t4 - t130 * t18 - t142 * t7 - t199 * t72 - t2 * t204 + t33 * t46 + t34 * t91 + t43 * t54, t1 * t204 - t129 * t3 + t130 * t19 + t142 * t8 - t212 * t43 - t23 * t72 - t31 * t46 + t34 * t93, -t1 * t91 + t18 * t23 + t19 * t199 - t2 * t93 + t212 * t4 - t3 * t54 + t31 * t7 - t33 * t8, t1 * t19 + t18 * t2 + t3 * t8 + t34 * t72 + t4 * t7 + t43 * t46; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t196, -t259, 0, 0, 0, 0, 0, 0, t256 * t190, t256 * t193, 0, -t259, 0, 0, 0, 0, 0, 0, qJD(1) * t140 + t133, -qJD(1) * t143 + t246, -t209 - t210, -qJD(1) * t157 + t200, 0, 0, 0, 0, 0, 0, -t204 * t263 - t145 * t113 - t152 * t236 + (t189 * t142 - t192 * t221) * t315, -t204 * t119 - t115 * t145 + t152 * t68 + (t142 * t192 + t189 * t221) * t315 (-t275 + t278) * t142 + t211 * qJD(1) - (qJD(5) * t211 - t216 - t288) * t204, -qJD(1) * t214 + t142 * t213 - t145 * t66 - t198 * t204 - t295, 0, 0, 0, 0, 0, 0, -t285 * t129 - t130 * t307 - t145 * t54 + t152 * t199, t284 * t129 - t130 * t308 + t145 * t212 + t152 * t23, -t199 * t308 - t212 * t285 + t23 * t307 + t284 * t54, -t1 * t308 - t145 * t46 - t152 * t34 + t2 * t307 - t284 * t8 - t285 * t7; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t237, -t257 * t196, 0, -t237, 0, 0, -t193 * t259, t190 * t259, 0, 0, -t272, -t139 + t303, -t253 - t131, t272, -t163 + (t143 + t234) * qJD(3), 0, qJD(3) * t85 + t140 * t238 - t143 * t157 - t51, -t227 - t157 * t140 + (t164 * t266 + t86) * qJD(3) + (-t193 * pkin(3) * t143 - t187 * t205) * qJD(1) (t75 - t85) * t143 - (-t74 + t86) * t140 + (t130 * t187 + t282 * t131) * pkin(3), t74 * t85 - t75 * t86 + (-t157 * t254 + t187 * t52 - t282 * t51) * pkin(3), t115 * t222 - t288 (-t68 + t281) * t192 - t315 * t275 + t283, -t276 + t316, t113 * t223 - t216, t208 + t280, -t315 * t143, t173 * t236 - t290 - t37 * t143 - t85 * t113 + (-t172 * t249 - t39) * t315 + t203 * t189, -t115 * t85 + t143 * t38 - t173 * t68 + t292 + (t172 * t250 + t40) * t315 + t203 * t192, t40 * t113 + t39 * t115 + (-t172 * t236 + t15 + t37 * t140 + (t115 * t172 - t37) * qJD(5)) * t192 + (t38 * t140 - t172 * t68 - t16 + (t113 * t172 - t38) * qJD(5)) * t189, t172 * t198 + t173 * t51 - t37 * t39 - t38 * t40 - t66 * t85, -t23 * t154 + t212 * t287, t220 - t305, t289 + t304, t199 * t306 + t286 * t54, t217 + t296, -t129 * t143, -t299 * t129 - t130 * t97 - t143 * t7 - t160 * t199 + t218 * t54 + t286 * t46 - t306 * t34, t298 * t129 + t130 * t98 + t143 * t8 + t154 * t34 - t160 * t23 - t212 * t218 - t287 * t46, t1 * t306 - t154 * t2 + t199 * t98 - t212 * t299 + t23 * t97 - t286 * t8 + t287 * t7 + t298 * t54, t1 * t98 + t160 * t34 + t2 * t97 + t218 * t46 - t298 * t8 - t299 * t7; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t163 + (t143 - t234) * qJD(3), t253 - t131, -t139 - t303, -t140 * t75 + t143 * t74 + t156, 0, 0, 0, 0, 0, 0, t208 - t280, -t276 - t316 (t68 + t281) * t192 + t115 * t223 + t283, -t143 * t66 + t322 * t189 + t321 * t192, 0, 0, 0, 0, 0, 0, t217 - t296, t289 - t304, t220 + t305, t1 * t154 - t143 * t46 + t2 * t306 - t286 * t7 - t287 * t8; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t277, -t113 ^ 2 + t115 ^ 2, t113 * t315 - t68, -t277, t115 * t315 - t236, -t130, -t115 * t66 + t321, t113 * t66 - t322, 0, 0, -t301, t314, t312, t301, t310, -t130, -t10 * t129 + (-t115 * t54 - t129 * t248 - t130 * t191) * pkin(5) + t313, t11 * t129 + (t115 * t212 - t129 * t247 + t130 * t188) * pkin(5) + t311, -t10 * t212 + t11 * t54 - t212 * t8 - t54 * t7 + (t188 * t199 + t191 * t23 + (-t188 * t212 - t191 * t54) * qJD(6)) * pkin(5), -t10 * t7 - t11 * t8 + (t1 * t188 - t115 * t46 + t191 * t2 + (-t188 * t7 + t191 * t8) * qJD(6)) * pkin(5); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t301, t314, t312, t301, t310, -t130, t129 * t8 + t313, t129 * t7 + t311, 0, 0;];
tauc_reg  = t5;
