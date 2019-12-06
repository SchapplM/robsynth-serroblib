% Calculate inertial parameters regressor of inverse dynamics joint torque vector for
% S5RRRRR2
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% qJD [5x1]
%   Generalized joint velocities
% qJDD [5x1]
%   Generalized joint accelerations
% g [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [2x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a4]';
% 
% Output:
% tau_reg [5x(5*10)]
%   inertial parameter regressor of inverse dynamics joint torque vector

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-05 18:54
% Revision: 77da58f92bca3eff71542919beafa37024070d86 (2019-12-05)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tau_reg = S5RRRRR2_invdynJ_fixb_reg2_slag_vp(qJ, qJD, qJDD, g, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(5,1),zeros(3,1),zeros(2,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RRRRR2_invdynJ_fixb_reg2_slag_vp: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RRRRR2_invdynJ_fixb_reg2_slag_vp: qJD has to be [5x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [5 1]), ...
  'S5RRRRR2_invdynJ_fixb_reg2_slag_vp: qJDD has to be [5x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S5RRRRR2_invdynJ_fixb_reg2_slag_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [2 1]), ...
  'S5RRRRR2_invdynJ_fixb_reg2_slag_vp: pkin has to be [2x1] (double)');

%% Symbolic Calculation
% From invdyn_joint_fixb_regressor_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-05 18:53:49
% EndTime: 2019-12-05 18:53:58
% DurationCPUTime: 4.29s
% Computational Cost: add. (4014->383), mult. (8590->544), div. (0->0), fcn. (6181->14), ass. (0->249)
t164 = sin(qJ(3));
t157 = t164 ^ 2;
t168 = cos(qJ(3));
t159 = t168 ^ 2;
t274 = t157 + t159;
t348 = pkin(1) * t274;
t163 = sin(qJ(4));
t165 = sin(qJ(2));
t241 = qJD(4) * t164 * t165;
t265 = qJDD(1) * t165;
t169 = cos(qJ(2));
t271 = qJD(2) * t169;
t178 = ((t168 * t271 - t241) * qJD(1) + t168 * t265) * pkin(1);
t316 = qJD(3) * pkin(2);
t347 = (qJD(4) * t316 + t178) * t163;
t332 = cos(qJ(4));
t238 = t332 * qJD(4);
t346 = t332 * qJD(3) + t238;
t273 = qJD(1) * t165;
t256 = pkin(1) * t273;
t195 = -t164 * t256 + t316;
t227 = t168 * t256;
t341 = -t163 * t227 + t332 * t195;
t345 = t341 * qJD(4);
t156 = qJD(1) + qJD(2);
t230 = qJD(1) * (-qJD(2) + t156);
t161 = qJ(1) + qJ(2);
t151 = sin(t161);
t153 = cos(t161);
t276 = g(1) * t153 + g(2) * t151;
t179 = (t169 * t230 - t265) * pkin(1) + t276;
t155 = qJDD(1) + qJDD(2);
t243 = t332 * t168;
t281 = t164 * t155;
t244 = t332 * t164;
t101 = t163 * t168 + t244;
t261 = qJD(3) + qJD(4);
t342 = t261 * t101;
t43 = -t155 * t243 + t156 * t342 + t163 * t281;
t162 = sin(qJ(5));
t167 = cos(qJ(5));
t266 = qJD(5) * t167;
t202 = t332 * t227;
t268 = qJD(3) * t165;
t226 = pkin(1) * t164 * t268;
t208 = qJD(1) * t226;
t66 = qJDD(3) * pkin(2) + (-t164 * t265 + (-t164 * t271 - t168 * t268) * qJD(1)) * pkin(1);
t185 = qJD(4) * t202 - t163 * t208 - t332 * t66;
t30 = t185 + t347;
t343 = t30 * t162 - t266 * t341;
t340 = t346 * t168;
t166 = sin(qJ(1));
t170 = cos(qJ(1));
t338 = g(1) * t166 - g(2) * t170;
t160 = qJ(3) + qJ(4);
t152 = cos(t160);
t150 = sin(t160);
t295 = t150 * t153;
t296 = t150 * t151;
t337 = g(1) * t295 + g(2) * t296 - g(3) * t152;
t336 = t165 * (-0.1e1 + t274);
t144 = g(1) * t151;
t264 = qJDD(1) * t169;
t148 = pkin(1) * t264;
t335 = t144 + t148;
t259 = qJDD(3) + qJDD(4);
t282 = t163 * t164;
t214 = t261 * t282;
t280 = t168 * t155;
t223 = -t155 * t244 - t340 * t156 - t163 * t280;
t42 = t156 * t214 + t223;
t91 = t101 * t156;
t72 = t162 * t261 + t167 * t91;
t27 = qJD(5) * t72 - t162 * t42 - t167 * t259;
t220 = t332 * t271;
t334 = (qJD(1) * t220 + t332 * t265) * pkin(1) * t168 - t332 * t208 + t345;
t331 = pkin(1) * t165;
t329 = pkin(2) * t168;
t327 = g(2) * t153;
t142 = g(3) * t150;
t324 = g(3) * t168;
t229 = t167 * t261;
t70 = t162 * t91 - t229;
t251 = t156 * t282;
t89 = -t156 * t243 + t251;
t82 = qJD(5) + t89;
t323 = t70 * t82;
t322 = t72 * t70;
t321 = t72 * t82;
t320 = t341 * t70;
t319 = t341 * t72;
t318 = t82 * t91;
t317 = t91 * t89;
t307 = t163 * t66;
t29 = t307 + t334;
t77 = t163 * t195 + t202;
t272 = qJD(1) * t169;
t255 = pkin(1) * t272;
t95 = -t156 * t329 - t255;
t51 = -t162 * t77 + t167 * t95;
t298 = qJD(5) * t51;
t254 = qJD(2) * t331;
t225 = qJD(1) * t254;
t269 = qJD(3) * t164;
t242 = t156 * t269;
t69 = t225 - t148 + (t242 - t280) * pkin(2);
t15 = t162 * t69 + t167 * t29 + t298;
t315 = t15 * t162;
t12 = t15 * t167;
t314 = t162 * t27;
t41 = qJDD(5) + t43;
t313 = t162 * t41;
t312 = t162 * t70;
t311 = t162 * t72;
t73 = t214 - t340;
t310 = t162 * t73;
t309 = t162 * t82;
t308 = t162 * t89;
t306 = t167 * t41;
t305 = t167 * t70;
t304 = t167 * t72;
t303 = t167 * t73;
t232 = t167 * t82;
t302 = t167 * t89;
t267 = qJD(5) * t162;
t26 = -qJD(5) * t229 - t162 * t259 + t167 * t42 + t91 * t267;
t301 = t26 * t162;
t300 = t26 * t167;
t299 = t27 * t167;
t25 = t30 * t101;
t52 = t162 * t95 + t167 * t77;
t297 = qJD(5) * t52;
t294 = t151 * t152;
t293 = t151 * t164;
t292 = t151 * t167;
t291 = t151 * t168;
t290 = t152 * t153;
t289 = t152 * t162;
t288 = t153 * t162;
t287 = t153 * t164;
t286 = t153 * t167;
t285 = t153 * t168;
t284 = t155 * t169;
t283 = t156 * t164;
t279 = -g(1) * t296 + g(2) * t295;
t278 = g(2) * t287 + t164 * t225;
t132 = g(1) * t291;
t277 = t168 * t148 + t132;
t275 = t157 - t159;
t270 = qJD(3) * t156;
t173 = pkin(1) ^ 2;
t263 = qJDD(1) * t173;
t262 = qJDD(3) * t168;
t260 = -qJDD(1) - t155;
t257 = pkin(2) * t283;
t252 = qJD(5) * t168 * t82;
t154 = t156 ^ 2;
t250 = t164 * t154 * t168;
t248 = g(1) * t290 + g(2) * t294 + t142;
t247 = t162 * t332;
t246 = t167 * t332;
t245 = (-t82 + t89) * t341;
t240 = t101 * t267;
t237 = qJD(1) * t271;
t193 = pkin(1) * t101;
t79 = t193 * t272;
t235 = pkin(2) * t132 + t341 * t79;
t65 = t167 * t69;
t16 = -t162 * t29 - t297 + t65;
t234 = -t16 - t297;
t231 = t69 * t101 - t95 * t73 + t279;
t228 = (-qJD(1) - t156) * qJD(2);
t222 = t82 * t238;
t221 = t168 * t242;
t196 = t243 - t282;
t78 = t196 * t256;
t217 = t341 * t78 + (g(1) * t287 + g(2) * t293) * pkin(2);
t212 = t162 * t52 + t167 * t51;
t211 = t162 * t51 - t167 * t52;
t210 = t304 + t312;
t209 = -0.2e1 * t270;
t129 = -pkin(1) * t169 - t329;
t88 = t196 * t331;
t60 = t129 * t167 - t162 * t88;
t61 = t129 * t162 + t167 * t88;
t207 = t169 ^ 2 * t263 + t338 * pkin(1);
t206 = t196 * t29 + t341 * t73 - t342 * t77 + t25 - t276;
t205 = -t51 * t302 - t52 * t308 + t12 - t248;
t204 = -qJD(5) * t95 + t142 - t29;
t203 = -t332 * t30 - t324;
t201 = -t276 + (t237 + t265) * t348;
t200 = -t168 * t41 + t82 * t269;
t199 = -t82 * t266 - t313;
t198 = -t82 * t267 + t306;
t197 = g(1) * t294 - g(2) * t290 - t196 * t69 + t342 * t95;
t192 = t169 * t196;
t191 = t230 * t331 - t327;
t84 = -t152 * t292 + t288;
t86 = t151 * t162 + t152 * t286;
t190 = -g(1) * t84 - g(2) * t86 + t343 * t101 - t16 * t196 + t310 * t341 + t342 * t51;
t39 = -t163 * t226 + (-t163 * t241 + t164 * t220 + (t163 * t271 + t346 * t165) * t168) * pkin(1);
t87 = t165 * t193;
t189 = -g(1) * (-pkin(1) * t166 - pkin(2) * t291) - g(2) * (pkin(1) * t170 + pkin(2) * t285) + t30 * t87 - t39 * t341;
t186 = -t212 * qJD(5) - t16 * t162;
t184 = -t341 * t267 - t51 * t91 + (-t30 + t337) * t167;
t83 = t151 * t289 + t286;
t85 = -t152 * t288 + t292;
t182 = -g(1) * t83 - g(2) * t85 + t15 * t196 + t167 * t25 - t52 * t342 - (-t240 - t303) * t341;
t181 = -t150 * t162 * t276 + g(3) * t289 + t52 * t91 + t343;
t177 = -t95 * t91 - t185 + t337;
t175 = t51 * t240 + (t234 * t167 - t315) * t101 + t212 * t73 - t279;
t174 = t95 * t89 + t248 - t334;
t171 = qJD(3) ^ 2;
t158 = t165 ^ 2;
t118 = -t164 * t171 + t262;
t117 = qJDD(3) * t164 + t168 * t171;
t105 = pkin(2) * t269 + t254;
t93 = t155 * t159 - 0.2e1 * t221;
t92 = t155 * t157 + 0.2e1 * t221;
t81 = qJD(1) * pkin(1) * t192;
t80 = t101 * t256;
t75 = 0.2e1 * t164 * t280 - 0.2e1 * t275 * t270;
t64 = t162 * t256 + t167 * t81;
t63 = -t162 * t81 + t167 * t256;
t59 = t162 * t257 - t167 * t80;
t58 = t162 * t80 + t167 * t257;
t50 = t196 * t259 - t261 * t342;
t49 = t101 * t259 - t73 * t261;
t44 = -t89 ^ 2 + t91 ^ 2;
t38 = (qJD(2) * t192 - t165 * t342) * pkin(1);
t32 = t91 * t261 - t43;
t31 = -t223 + (-t251 + t89) * t261;
t21 = -qJD(5) * t61 + t167 * t105 - t162 * t38;
t20 = qJD(5) * t60 + t162 * t105 + t167 * t38;
t19 = -t196 * t43 + t342 * t89;
t18 = -t101 * t42 - t73 * t91;
t17 = -t196 * t41 + t342 * t82;
t14 = t82 * t232 - t72 * t91 + t313;
t13 = -t82 ^ 2 * t162 + t70 * t91 + t306;
t11 = t70 * t309 - t299;
t10 = t72 * t232 - t301;
t7 = -t70 * t310 + (t70 * t266 + t314) * t101;
t6 = -t72 * t303 + (-t72 * t267 - t300) * t101;
t5 = -t101 * t43 - t196 * t42 - t342 * t91 + t73 * t89;
t4 = t101 * t199 + t196 * t27 + t73 * t309 - t342 * t70;
t3 = t101 * t198 + t196 * t26 - t73 * t232 + t342 * t72;
t2 = (-t26 - t323) * t167 + (-t27 - t321) * t162;
t1 = (t305 + t311) * t73 + (t301 - t299 + (-t304 + t312) * qJD(5)) * t101;
t8 = [0, 0, 0, 0, 0, qJDD(1), t338, g(1) * t170 + g(2) * t166, 0, 0, 0, 0, 0, 0, 0, t155, -t327 + (t165 * t228 + t284) * pkin(1) + t335, (t260 * t165 + t169 * t228) * pkin(1) + t276, 0, t158 * t263 + t207, t92, t75, t117, t93, t118, 0, -g(2) * t285 + ((t284 + (-t171 + t228) * t165) * t168 + (-qJDD(3) * t165 + t169 * t209) * t164) * pkin(1) + t277, -g(1) * t293 + ((-t262 + (qJD(2) * t156 + t171) * t164) * t165 + (t164 * t260 + t168 * t209) * t169) * pkin(1) + t278, t201 + (t155 * t165 + t156 * t271) * t348, (qJDD(1) * t158 * t274 + 0.2e1 * t237 * t336) * t173 + t207, t18, t5, t49, t19, t50, 0, t105 * t89 + t129 * t43 - t259 * t87 - t261 * t39 + t197, t105 * t91 - t129 * t42 - t259 * t88 - t261 * t38 + t231, -t38 * t89 + t39 * t91 - t42 * t87 - t43 * t88 + t206, t105 * t95 + t129 * t69 + t29 * t88 + t38 * t77 + t189, t6, t1, t3, t7, t4, t17, t21 * t82 + t27 * t87 + t39 * t70 + t41 * t60 + t190, -t20 * t82 - t26 * t87 + t39 * t72 - t41 * t61 + t182, -t20 * t70 - t21 * t72 + t26 * t60 - t27 * t61 + t175, t15 * t61 + t16 * t60 + t20 * t52 + t21 * t51 + t189; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t155, t191 + t335, t179, 0, 0, t92, t75, t117, t93, t118, 0, t168 * t191 + t277, (-t144 + (-t156 * t273 - t264) * pkin(1)) * t164 + t278, -t156 * t255 * t274 + t201, -t173 * qJD(1) ^ 2 * t169 * t336, t18, t5, t49, t19, t50, 0, -t89 * t256 + t79 * t261 + (-t168 * t43 + t269 * t89) * pkin(2) + t197, -t91 * t256 + t81 * t261 + (t168 * t42 + t269 * t91) * pkin(2) + t231, -t79 * t91 + t81 * t89 + t206, -t95 * t256 - t77 * t81 + (t95 * t269 + (-t69 - t327) * t168) * pkin(2) + t235, t6, t1, t3, t7, t4, t17, -t63 * t82 - t70 * t79 + (t162 * t252 + t167 * t200) * pkin(2) + t190, t64 * t82 - t72 * t79 + (-t162 * t200 + t167 * t252) * pkin(2) + t182, t63 * t72 + t64 * t70 + (-t210 * t269 + (t314 - t300 + (t305 - t311) * qJD(5)) * t168) * pkin(2) + t175, -t51 * t63 - t52 * t64 + (t212 * t269 + (qJD(5) * t211 - t16 * t167 - t315 - t327) * t168) * pkin(2) + t235; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t250, t275 * t154, t281, t250, t280, qJDD(3), t164 * t179 - t324, g(3) * t164 + t168 * t179, 0, 0, t317, t44, t31, -t317, t32, t259, t78 * t261 + (t332 * t259 - t89 * t283) * pkin(2) + ((-0.2e1 * qJD(3) - qJD(4)) * qJD(4) * pkin(2) - t178) * t163 + t177, -t91 * t257 + (-pkin(2) * t259 - t66) * t163 + t174 + (-pkin(2) * t238 - t80) * t261, (t77 - t78) * t91 + (-t341 - t80) * t89 + (t332 * t42 - t163 * t43 + (t163 * t91 - t332 * t89) * qJD(4)) * pkin(2), t77 * t80 + (-t95 * t283 + t163 * t29 + (-t163 * t341 + t332 * t77) * qJD(4) + t203) * pkin(2) + t217, t10, t2, t14, t11, t13, -t318, -t341 * t308 - t58 * t82 - t78 * t70 + (-t162 * t222 - t332 * t27 + (qJD(4) * t70 + t199) * t163) * pkin(2) + t184, -t341 * t302 + t59 * t82 - t78 * t72 + (-t167 * t222 + t332 * t26 + (qJD(4) * t72 - t198) * t163) * pkin(2) + t181, t58 * t72 + t59 * t70 + ((-t246 * t70 + t247 * t72) * qJD(4) + (qJD(5) * t210 - t299 - t301) * t163) * pkin(2) + t186 + t205, -t51 * t58 - t52 * t59 + ((t246 * t52 - t247 * t51) * qJD(4) + (t12 + t186 - t345) * t163 + t203) * pkin(2) + t217; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t317, t44, t31, -t317, t32, t259, t77 * t261 + t177 - t347, t261 * t341 + t174 - t307, 0, 0, t10, t2, t14, t11, t13, -t318, -t162 * t245 - t77 * t70 + t184, -t167 * t245 - t77 * t72 + t181, (-t298 + t320) * t167 + (t234 - t319) * t162 + t205, -(-t211 - t77) * t341; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t322, -t70 ^ 2 + t72 ^ 2, -t26 + t323, -t322, -t27 + t321, t41, -g(1) * t85 + g(2) * t83 + t162 * t204 - t266 * t77 + t52 * t82 + t319 + t65, g(1) * t86 - g(2) * t84 + t51 * t82 - t320 + (qJD(5) * t77 - t69) * t162 + t204 * t167, 0, 0;];
tau_reg = t8;
