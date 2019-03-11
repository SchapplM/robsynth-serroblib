% Calculate minimal parameter regressor of coriolis joint torque vector for
% S6PRRPRR5
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% qJD [6x1]
%   Generalized joint velocities
% pkin [12x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,alpha2,d2,d3,d5,d6,theta1,theta4]';
% 
% Output:
% tauc_reg [6x29]
%   minimal parameter regressor of coriolis joint torque vector

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-08 22:21
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauc_reg = S6PRRPRR5_coriolisvecJ_fixb_regmin_slag_vp(qJ, qJD, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(12,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6PRRPRR5_coriolisvecJ_fixb_regmin_slag_vp: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6PRRPRR5_coriolisvecJ_fixb_regmin_slag_vp: qJD has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [12 1]), ...
  'S6PRRPRR5_coriolisvecJ_fixb_regmin_slag_vp: pkin has to be [12x1] (double)');

%% Symbolic Calculation
% From coriolisvec_joint_fixb_regressor_minpar_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-08 22:20:57
% EndTime: 2019-03-08 22:21:06
% DurationCPUTime: 4.42s
% Computational Cost: add. (4371->408), mult. (11320->602), div. (0->0), fcn. (9030->12), ass. (0->211)
t202 = cos(qJ(3));
t265 = qJD(2) * t202;
t183 = -qJD(5) + t265;
t178 = -qJD(6) + t183;
t196 = sin(qJ(6));
t200 = cos(qJ(6));
t192 = sin(pkin(12));
t194 = cos(pkin(12));
t255 = t194 * qJD(3);
t198 = sin(qJ(3));
t266 = qJD(2) * t198;
t154 = t192 * t266 - t255;
t264 = qJD(3) * t192;
t156 = t194 * t266 + t264;
t197 = sin(qJ(5));
t201 = cos(qJ(5));
t93 = t154 * t197 - t156 * t201;
t94 = t201 * t154 + t156 * t197;
t305 = t196 * t93 - t200 * t94;
t319 = t178 * t305;
t226 = pkin(3) * t198 - qJ(4) * t202;
t138 = qJD(3) * t226 - t198 * qJD(4);
t263 = qJD(3) * t198;
t253 = pkin(8) * t263;
t176 = t192 * t253;
t199 = sin(qJ(2));
t193 = sin(pkin(6));
t269 = qJD(1) * t193;
t203 = cos(qJ(2));
t278 = t202 * t203;
t275 = t194 * t138 - (-t192 * t278 + t194 * t199) * t269 + t176;
t318 = t192 * t138 - (t192 * t199 + t194 * t278) * t269;
t280 = t194 * t202;
t219 = pkin(4) * t198 - pkin(9) * t280;
t213 = t219 * qJD(3);
t317 = -t213 - t275;
t281 = t194 * t198;
t285 = t192 * t202;
t316 = (-pkin(8) * t281 - pkin(9) * t285) * qJD(3) + t318;
t249 = t192 * t265;
t163 = t226 * qJD(2);
t246 = t199 * t269;
t167 = qJD(2) * pkin(8) + t246;
t195 = cos(pkin(6));
t268 = qJD(1) * t195;
t301 = -t198 * t167 + t202 * t268;
t77 = t192 * t163 + t194 * t301;
t63 = -pkin(9) * t249 + t77;
t315 = qJD(4) * t194 - t63;
t236 = -qJD(3) * pkin(3) + qJD(4);
t115 = -t301 + t236;
t78 = pkin(4) * t154 + t115;
t34 = pkin(5) * t94 + t78;
t314 = t34 * t305;
t222 = t196 * t94 + t200 * t93;
t313 = t222 * t305;
t312 = t93 * t183;
t311 = t94 * t183;
t310 = t178 * t222;
t279 = t201 * t194;
t287 = t192 * t197;
t161 = -t279 + t287;
t214 = t161 * t202;
t273 = qJD(2) * t214 - t161 * qJD(5);
t162 = t192 * t201 + t194 * t197;
t215 = t162 * t202;
t272 = -qJD(2) * t215 + t162 * qJD(5);
t309 = t222 ^ 2 - t305 ^ 2;
t256 = qJD(6) * t200;
t257 = qJD(6) * t196;
t254 = qJD(2) * qJD(3);
t242 = t202 * t254;
t230 = t192 * t242;
t258 = qJD(5) * t201;
t53 = -t154 * t258 + t242 * t279 + (-qJD(5) * t156 - t230) * t197;
t212 = qJD(3) * t215;
t302 = qJD(5) * t93;
t54 = qJD(2) * t212 - t302;
t7 = -t196 * t54 + t200 * t53 - t94 * t256 + t257 * t93;
t308 = t7 + t319;
t186 = t198 * t254;
t179 = t198 * t268;
t124 = t202 * t167 + t179;
t119 = qJD(3) * qJ(4) + t124;
t169 = -pkin(3) * t202 - qJ(4) * t198 - pkin(2);
t245 = t203 * t269;
t125 = qJD(2) * t169 - t245;
t57 = -t119 * t192 + t194 * t125;
t35 = -pkin(4) * t265 - pkin(9) * t156 + t57;
t58 = t194 * t119 + t192 * t125;
t43 = -pkin(9) * t154 + t58;
t17 = t197 * t35 + t201 * t43;
t106 = (t138 + t246) * qJD(2);
t231 = qJD(2) * t245;
t84 = t202 * t231 + (qJD(4) + t301) * qJD(3);
t32 = t194 * t106 - t192 * t84;
t25 = qJD(2) * t213 + t32;
t33 = t192 * t106 + t194 * t84;
t29 = -pkin(9) * t230 + t33;
t239 = -t197 * t29 + t201 * t25;
t207 = -qJD(5) * t17 + t239;
t2 = pkin(5) * t186 - t53 * pkin(10) + t207;
t260 = qJD(5) * t197;
t218 = t197 * t25 + t201 * t29 + t35 * t258 - t260 * t43;
t3 = -pkin(10) * t54 + t218;
t251 = -t196 * t3 + t200 * t2;
t16 = -t197 * t43 + t201 * t35;
t12 = pkin(10) * t93 + t16;
t10 = -pkin(5) * t183 + t12;
t13 = -pkin(10) * t94 + t17;
t291 = t13 * t200;
t5 = t10 * t196 + t291;
t307 = -qJD(6) * t5 + t34 * t222 + t251;
t206 = qJD(6) * t222 - t196 * t53 - t200 * t54;
t306 = t206 + t310;
t153 = t194 * t169;
t102 = -pkin(9) * t281 + t153 + (-pkin(8) * t192 - pkin(4)) * t202;
t127 = pkin(8) * t280 + t192 * t169;
t286 = t192 * t198;
t114 = -pkin(9) * t286 + t127;
t304 = t102 * t258 - t114 * t260 - t197 * t317 + t316 * t201;
t303 = t316 * t197 + t201 * t317;
t296 = pkin(9) + qJ(4);
t171 = t296 * t192;
t172 = t296 * t194;
t271 = -t197 * t171 + t201 * t172;
t220 = qJD(4) * t192 + qJD(5) * t172;
t76 = t194 * t163 - t192 * t301;
t56 = qJD(2) * t219 + t76;
t300 = -t171 * t258 + t315 * t201 + (-t220 - t56) * t197;
t11 = t13 * t257;
t241 = qJD(6) * t10 + t3;
t299 = t196 * t2 + t200 * t241 - t11;
t288 = t197 * t102 + t201 * t114;
t259 = qJD(5) * t198;
t81 = -qJD(3) * t214 - t162 * t259;
t298 = -pkin(5) * t263 + t81 * pkin(10) + qJD(5) * t288 + t303;
t82 = t258 * t281 - t259 * t287 + t212;
t297 = -pkin(10) * t82 + t304;
t98 = t200 * t161 + t162 * t196;
t295 = -qJD(6) * t98 - t196 * t272 + t200 * t273;
t99 = -t161 * t196 + t162 * t200;
t294 = qJD(6) * t99 + t196 * t273 + t200 * t272;
t292 = qJD(2) * pkin(2);
t262 = qJD(3) * t202;
t88 = qJD(3) * t179 + t167 * t262 + t198 * t231;
t290 = t192 * t88;
t289 = t194 * t88;
t284 = t193 * t199;
t283 = t193 * t203;
t205 = qJD(2) ^ 2;
t282 = t193 * t205;
t204 = qJD(3) ^ 2;
t277 = t204 * t198;
t276 = t204 * t202;
t234 = t194 * t253;
t274 = -t234 + t318;
t148 = (pkin(4) * t192 + pkin(8)) * t262;
t164 = pkin(4) * t286 + t198 * pkin(8);
t270 = t198 ^ 2 - t202 ^ 2;
t267 = qJD(2) * t193;
t252 = t199 * t282;
t105 = pkin(4) * t249 + t124;
t185 = -pkin(4) * t194 - pkin(3);
t248 = t199 * t267;
t247 = t203 * t267;
t244 = pkin(5) * t272 - t105;
t237 = t201 * t102 - t114 * t197;
t235 = -t201 * t171 - t172 * t197;
t66 = pkin(4) * t230 + t88;
t233 = t198 * t247;
t232 = t202 * t247;
t52 = t201 * t56;
t80 = -pkin(10) * t161 + t271;
t229 = pkin(5) * t266 + pkin(10) * t273 + t162 * qJD(4) + t271 * qJD(5) + qJD(6) * t80 - t197 * t63 + t52;
t79 = -pkin(10) * t162 + t235;
t228 = -pkin(10) * t272 + qJD(6) * t79 + t300;
t168 = -t245 - t292;
t227 = -t168 - t245;
t136 = t161 * t198;
t26 = -pkin(5) * t202 + pkin(10) * t136 + t237;
t135 = t162 * t198;
t28 = -pkin(10) * t135 + t288;
t225 = t196 * t26 + t200 * t28;
t145 = t195 * t198 + t202 * t284;
t110 = -t145 * t192 - t194 * t283;
t111 = t145 * t194 - t192 * t283;
t44 = t110 * t201 - t111 * t197;
t45 = t110 * t197 + t111 * t201;
t224 = -t196 * t45 + t200 * t44;
t223 = t196 * t44 + t200 * t45;
t71 = t200 * t135 - t136 * t196;
t72 = -t135 * t196 - t136 * t200;
t144 = -t195 * t202 + t198 * t284;
t211 = qJD(3) * (-t227 - t292);
t208 = -qJ(4) * t263 + (-t115 + t236) * t202;
t131 = pkin(5) * t161 + t185;
t126 = -pkin(8) * t285 + t153;
t113 = qJD(3) * t145 + t233;
t112 = -qJD(3) * t144 + t232;
t104 = pkin(5) * t135 + t164;
t75 = t112 * t194 + t192 * t248;
t74 = -t112 * t192 + t194 * t248;
t55 = pkin(5) * t82 + t148;
t22 = pkin(5) * t54 + t66;
t21 = qJD(6) * t72 + t196 * t81 + t200 * t82;
t20 = -qJD(6) * t71 - t196 * t82 + t200 * t81;
t15 = -qJD(5) * t45 - t197 * t75 + t201 * t74;
t14 = qJD(5) * t44 + t197 * t74 + t201 * t75;
t4 = t10 * t200 - t13 * t196;
t1 = [0, 0, -t252, -t203 * t282, 0, 0, 0, 0, 0, -t202 * t252 + (-t113 - t233) * qJD(3), t198 * t252 + (-t112 - t232) * qJD(3), t113 * t154 + (-t202 * t74 + (t110 * t198 + t144 * t285) * qJD(3)) * qJD(2), t113 * t156 + (t202 * t75 + (-t111 * t198 + t144 * t280) * qJD(3)) * qJD(2), -t75 * t154 - t74 * t156 + (-t110 * t194 - t111 * t192) * t242, t110 * t32 + t111 * t33 + t113 * t115 + t144 * t88 + t57 * t74 + t58 * t75, 0, 0, 0, 0, 0, t113 * t94 + t144 * t54 - t15 * t183 + t186 * t44, -t113 * t93 + t14 * t183 + t144 * t53 - t186 * t45, 0, 0, 0, 0, 0 -(-qJD(6) * t223 - t196 * t14 + t200 * t15) * t178 + t224 * t186 - t113 * t305 - t144 * t206 (qJD(6) * t224 + t200 * t14 + t196 * t15) * t178 - t223 * t186 - t113 * t222 + t144 * t7; 0, 0, 0, 0, 0.2e1 * t202 * t186, -0.2e1 * t270 * t254, t276, -t277, 0, -pkin(8) * t276 + t198 * t211, pkin(8) * t277 + t202 * t211 (-t154 * t245 + t290 + (qJD(2) * t126 + t57) * qJD(3)) * t198 + (-t32 + (pkin(8) * t154 + t115 * t192) * qJD(3) + (t176 - t275) * qJD(2)) * t202 (-t156 * t245 + t289 + (-qJD(2) * t127 - t58) * qJD(3)) * t198 + (t33 + (pkin(8) * t156 + t115 * t194) * qJD(3) + (t234 + t274) * qJD(2)) * t202 (-t192 * t33 - t194 * t32) * t198 - t275 * t156 - t274 * t154 + (-t192 * t58 - t194 * t57 + (-t126 * t194 - t127 * t192) * qJD(2)) * t262, -t115 * t198 * t245 + t126 * t32 + t127 * t33 + t274 * t58 + t275 * t57 + (t115 * t262 + t198 * t88) * pkin(8), -t136 * t53 - t81 * t93, -t135 * t53 + t136 * t54 - t81 * t94 + t82 * t93, -t81 * t183 - t53 * t202 + (-qJD(2) * t136 - t93) * t263, t82 * t183 + t54 * t202 + (-qJD(2) * t135 - t94) * t263 (-t183 - t265) * t263, -t239 * t202 + t148 * t94 + t164 * t54 + t66 * t135 + t78 * t82 + t303 * t183 + (t17 * t202 + t183 * t288) * qJD(5) + (-t94 * t245 + (qJD(2) * t237 + t16) * qJD(3)) * t198, t218 * t202 - t148 * t93 + t164 * t53 - t66 * t136 + t78 * t81 + t304 * t183 + (t93 * t245 + (-t288 * qJD(2) - t17) * qJD(3)) * t198, -t20 * t222 + t7 * t72, t20 * t305 + t206 * t72 + t21 * t222 - t7 * t71, -t20 * t178 - t7 * t202 + (qJD(2) * t72 - t222) * t263, t21 * t178 - t206 * t202 + (-qJD(2) * t71 + t305) * t263 (-t178 - t265) * t263, -t251 * t202 - t55 * t305 - t104 * t206 + t22 * t71 + t34 * t21 + (t297 * t196 + t298 * t200) * t178 + (t178 * t225 + t202 * t5) * qJD(6) + (t305 * t245 + ((-t196 * t28 + t200 * t26) * qJD(2) + t4) * qJD(3)) * t198, t299 * t202 - t55 * t222 + t104 * t7 + t22 * t72 + t34 * t20 + ((qJD(6) * t26 + t297) * t200 + (-qJD(6) * t28 - t298) * t196) * t178 + (t222 * t245 + (-qJD(2) * t225 - t5) * qJD(3)) * t198; 0, 0, 0, 0, -t198 * t205 * t202, t270 * t205, 0, 0, 0, qJD(3) * t124 - t168 * t266 - t88, t227 * t265, -t124 * t154 - t289 + (t192 * t208 - t198 * t57 + t202 * t76) * qJD(2), -t124 * t156 + t290 + (t194 * t208 + t198 * t58 - t202 * t77) * qJD(2), t154 * t77 + t156 * t76 + (-qJD(4) * t154 + t265 * t57 + t33) * t194 + (qJD(4) * t156 + t265 * t58 - t32) * t192, -pkin(3) * t88 - t115 * t124 - t57 * t76 - t58 * t77 + (-t192 * t57 + t194 * t58) * qJD(4) + (-t192 * t32 + t194 * t33) * qJ(4), t53 * t162 - t273 * t93, -t53 * t161 - t162 * t54 + t272 * t93 - t273 * t94, -t273 * t183 + (qJD(3) * t162 + t93) * t266, t272 * t183 + (-qJD(3) * t161 + t94) * t266, t183 * t266, -t105 * t94 + t66 * t161 + t185 * t54 + t272 * t78 + (t52 + t220 * t201 + (-qJD(5) * t171 + t315) * t197) * t183 + (qJD(3) * t235 - t16) * t266, t105 * t93 + t66 * t162 + t185 * t53 + t273 * t78 + t300 * t183 + (-qJD(3) * t271 + t17) * t266, -t222 * t295 + t7 * t99, t206 * t99 + t222 * t294 + t295 * t305 - t7 * t98, -t295 * t178 + (qJD(3) * t99 + t222) * t266, t294 * t178 + (-qJD(3) * t98 - t305) * t266, t178 * t266, -t131 * t206 + t22 * t98 - t244 * t305 + t294 * t34 + (t196 * t228 + t200 * t229) * t178 + ((-t196 * t80 + t200 * t79) * qJD(3) - t4) * t266, t131 * t7 + t22 * t99 - t244 * t222 + t295 * t34 + (-t196 * t229 + t200 * t228) * t178 + (-(t196 * t79 + t200 * t80) * qJD(3) + t5) * t266; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0 (-t156 + t264) * t265 (t154 + t255) * t265, -t154 ^ 2 - t156 ^ 2, t154 * t58 + t156 * t57 + t88, 0, 0, 0, 0, 0, t54 + t312, t53 + t311, 0, 0, 0, 0, 0, -t206 + t310, t7 - t319; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t93 * t94, t93 ^ 2 - t94 ^ 2, t53 - t311, -t162 * t242 + t302 + t312, t186, -t17 * t183 + t78 * t93 + t207, -t16 * t183 + t78 * t94 - t218, t313, t309, t308, t306, t186 (-t12 * t196 - t291) * t178 + (t178 * t257 + t186 * t200 - t305 * t93) * pkin(5) + t307, -t314 + t11 + (t13 * t178 - t2) * t196 + (-t12 * t178 - t241) * t200 + (t178 * t256 - t186 * t196 - t222 * t93) * pkin(5); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t313, t309, t308, t306, t186, -t5 * t178 + t307, -t4 * t178 - t299 - t314;];
tauc_reg  = t1;
