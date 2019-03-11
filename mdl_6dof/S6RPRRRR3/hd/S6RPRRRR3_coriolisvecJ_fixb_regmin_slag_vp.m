% Calculate minimal parameter regressor of coriolis joint torque vector for
% S6RPRRRR3
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% qJD [6x1]
%   Generalized joint velocities
% pkin [11x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d3,d4,d5,d6,theta2]';
% 
% Output:
% tauc_reg [6x32]
%   minimal parameter regressor of coriolis joint torque vector

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-09 07:03
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauc_reg = S6RPRRRR3_coriolisvecJ_fixb_regmin_slag_vp(qJ, qJD, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(11,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPRRRR3_coriolisvecJ_fixb_regmin_slag_vp: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RPRRRR3_coriolisvecJ_fixb_regmin_slag_vp: qJD has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [11 1]), ...
  'S6RPRRRR3_coriolisvecJ_fixb_regmin_slag_vp: pkin has to be [11x1] (double)');

%% Symbolic Calculation
% From coriolisvec_joint_fixb_regressor_minpar_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 07:02:30
% EndTime: 2019-03-09 07:02:44
% DurationCPUTime: 5.62s
% Computational Cost: add. (5723->410), mult. (13595->570), div. (0->0), fcn. (9797->10), ass. (0->221)
t204 = sin(qJ(5));
t210 = cos(qJ(3));
t284 = qJD(1) * t210;
t187 = -qJD(4) + t284;
t209 = cos(qJ(4));
t205 = sin(qJ(4));
t282 = qJD(3) * t205;
t206 = sin(qJ(3));
t285 = qJD(1) * t206;
t162 = t209 * t285 + t282;
t189 = sin(pkin(11)) * pkin(1) + pkin(7);
t177 = t189 * qJD(1);
t197 = t206 * qJD(2);
t133 = t210 * t177 + t197;
t123 = qJD(3) * pkin(8) + t133;
t190 = -cos(pkin(11)) * pkin(1) - pkin(2);
t154 = -pkin(3) * t210 - pkin(8) * t206 + t190;
t127 = t154 * qJD(1);
t71 = -t123 * t205 + t209 * t127;
t63 = -pkin(9) * t162 + t71;
t56 = -pkin(4) * t187 + t63;
t208 = cos(qJ(5));
t273 = t209 * qJD(3);
t160 = t205 * t285 - t273;
t300 = t205 * t127;
t72 = t123 * t209 + t300;
t64 = -pkin(9) * t160 + t72;
t62 = t208 * t64;
t20 = t204 * t56 + t62;
t101 = t208 * t160 + t162 * t204;
t345 = pkin(10) * t101;
t14 = t20 - t345;
t203 = sin(qJ(6));
t275 = qJD(6) * t203;
t12 = t14 * t275;
t207 = cos(qJ(6));
t234 = t160 * t204 - t208 * t162;
t51 = t207 * t101 - t203 * t234;
t329 = qJD(2) * t210 - t206 * t177;
t122 = -qJD(3) * pkin(3) - t329;
t98 = pkin(4) * t160 + t122;
t57 = pkin(5) * t101 + t98;
t349 = t51 * t57 + t12;
t272 = qJD(1) * qJD(3);
t193 = t206 * t272;
t271 = qJD(3) * qJD(4);
t279 = qJD(4) * t205;
t262 = t206 * t279;
t327 = t210 * t273 - t262;
t116 = t327 * qJD(1) + t209 * t271;
t124 = t329 * qJD(3);
t238 = pkin(3) * t206 - pkin(8) * t210;
t172 = t238 * qJD(3);
t153 = qJD(1) * t172;
t246 = t205 * t124 - t209 * t153;
t214 = -qJD(4) * t72 - t246;
t25 = pkin(4) * t193 - pkin(9) * t116 + t214;
t278 = qJD(4) * t209;
t261 = t206 * t278;
t280 = qJD(3) * t210;
t264 = t205 * t280;
t218 = t261 + t264;
t117 = qJD(1) * t218 + t205 * t271;
t268 = t209 * t124 + t127 * t278 + t205 * t153;
t223 = -t123 * t279 + t268;
t28 = -pkin(9) * t117 + t223;
t255 = -t204 * t28 + t208 * t25;
t216 = -qJD(5) * t20 + t255;
t276 = qJD(5) * t208;
t277 = qJD(5) * t204;
t36 = t208 * t116 - t204 * t117 - t160 * t276 - t162 * t277;
t2 = pkin(5) * t193 - pkin(10) * t36 + t216;
t348 = -t203 * t2 + t349;
t235 = t101 * t203 + t207 * t234;
t342 = t235 * t51;
t339 = t235 ^ 2 - t51 ^ 2;
t183 = -qJD(5) + t187;
t176 = -qJD(6) + t183;
t213 = qJD(5) * t234 - t116 * t204 - t208 * t117;
t274 = qJD(6) * t207;
t8 = -t101 * t274 + t203 * t213 + t207 * t36 + t234 * t275;
t337 = -t176 * t51 + t8;
t250 = -t204 * t25 - t208 * t28 - t56 * t276 + t64 * t277;
t3 = pkin(10) * t213 - t250;
t267 = t207 * t2 - t203 * t3;
t347 = t57 * t235 + t267;
t215 = qJD(6) * t235 - t203 * t36 + t207 * t213;
t332 = t176 * t235 + t215;
t265 = t205 * t284;
t324 = pkin(8) + pkin(9);
t266 = qJD(4) * t324;
t169 = t238 * qJD(1);
t291 = t205 * t169 + t209 * t329;
t341 = pkin(9) * t265 - t205 * t266 - t291;
t296 = t209 * t210;
t230 = pkin(4) * t206 - pkin(9) * t296;
t245 = t209 * t169 - t205 * t329;
t346 = qJD(1) * t230 + t209 * t266 + t245;
t344 = pkin(10) * t234;
t340 = t234 * t101;
t163 = t204 * t205 - t208 * t209;
t224 = t163 * t210;
t326 = qJD(4) + qJD(5);
t293 = qJD(1) * t224 - t326 * t163;
t164 = t204 * t209 + t205 * t208;
t292 = (-t284 + t326) * t164;
t338 = -t101 ^ 2 + t234 ^ 2;
t336 = -t101 * t183 + t36;
t335 = t101 * t98 + t250;
t334 = t98 * t234 + t216;
t60 = t204 * t64;
t19 = t208 * t56 - t60;
t13 = t19 + t344;
t10 = -pkin(5) * t183 + t13;
t315 = t207 * t14;
t5 = t203 * t10 + t315;
t333 = -qJD(6) * t5 + t347;
t331 = t183 * t234 + t213;
t330 = t346 * t208;
t136 = t164 * t206;
t239 = -t133 + (-t265 + t279) * pkin(4);
t180 = t324 * t205;
t181 = t324 * t209;
t288 = -t204 * t180 + t208 * t181;
t328 = -t180 * t276 - t181 * t277 - t346 * t204 + t208 * t341;
t199 = t206 ^ 2;
t233 = qJD(1) * t199 - t187 * t210;
t325 = -t187 * t262 - t233 * t273;
t66 = -qJD(3) * t224 - t326 * t136;
t297 = t206 * t209;
t299 = t205 * t206;
t67 = -t277 * t299 + (t326 * t297 + t264) * t208 + t327 * t204;
t137 = t163 * t206;
t85 = -t136 * t203 - t137 * t207;
t16 = qJD(6) * t85 + t203 * t66 + t207 * t67;
t84 = t207 * t136 - t137 * t203;
t323 = t16 * t176 - t84 * t193;
t105 = t207 * t163 + t164 * t203;
t322 = -qJD(6) * t105 - t292 * t203 + t293 * t207;
t106 = -t163 * t203 + t164 * t207;
t321 = qJD(6) * t106 + t293 * t203 + t292 * t207;
t320 = t208 * t63 - t60;
t139 = t209 * t154;
t303 = t189 * t205;
t86 = -pkin(9) * t297 + t139 + (-pkin(4) - t303) * t210;
t168 = t189 * t296;
t287 = t205 * t154 + t168;
t93 = -pkin(9) * t299 + t287;
t318 = t204 * t86 + t208 * t93;
t317 = t292 * pkin(5) + t239;
t316 = t207 * t10;
t314 = -t136 * t193 + t67 * t183;
t313 = t116 * t205;
t312 = t117 * t210;
t311 = t122 * t205;
t310 = t122 * t209;
t125 = qJD(3) * t197 + t177 * t280;
t309 = t125 * t205;
t308 = t125 * t209;
t307 = t160 * t187;
t306 = t160 * t206;
t305 = t162 * t187;
t304 = t187 * t209;
t302 = t203 * t204;
t301 = t204 * t207;
t298 = t205 * t210;
t211 = qJD(3) ^ 2;
t295 = t211 * t206;
t294 = t211 * t210;
t290 = t154 * t278 + t205 * t172;
t281 = qJD(3) * t206;
t289 = t209 * t172 + t281 * t303;
t142 = pkin(4) * t299 + t206 * t189;
t286 = -t210 ^ 2 + t199;
t178 = qJD(1) * t190;
t270 = pkin(4) * qJD(5) * t176;
t109 = t218 * pkin(4) + t189 * t280;
t196 = -pkin(4) * t209 - pkin(3);
t259 = -t210 * t8 - t235 * t281;
t257 = qJD(6) * t10 + t3;
t43 = t230 * qJD(3) + (-t168 + (pkin(9) * t206 - t154) * t205) * qJD(4) + t289;
t45 = (-t206 * t273 - t210 * t279) * t189 - t218 * pkin(9) + t290;
t254 = -t204 * t45 + t208 * t43;
t253 = -t204 * t63 - t62;
t252 = -t204 * t93 + t208 * t86;
t251 = -t210 * t36 - t234 * t281;
t248 = -t116 * t210 + t162 * t281;
t247 = t187 * t189 + t123;
t244 = -t208 * t180 - t181 * t204;
t242 = t187 * t261;
t79 = pkin(4) * t117 + t125;
t89 = -pkin(10) * t163 + t288;
t241 = pkin(5) * t285 + t293 * pkin(10) + t288 * qJD(5) + qJD(6) * t89 + t341 * t204 + t330;
t88 = -pkin(10) * t164 + t244;
t240 = -t292 * pkin(10) + qJD(6) * t88 + t328;
t29 = -pkin(5) * t210 + pkin(10) * t137 + t252;
t30 = -pkin(10) * t136 + t318;
t237 = t203 * t29 + t207 * t30;
t231 = 0.2e1 * qJD(3) * t178;
t195 = pkin(4) * t208 + pkin(5);
t229 = pkin(4) * t301 + t195 * t203;
t228 = -pkin(4) * t302 + t195 * t207;
t227 = -t210 * t215 - t51 * t281;
t226 = -t101 * t281 - t210 * t213;
t225 = t204 * t43 + t208 * t45 + t86 * t276 - t93 * t277;
t221 = t233 * t205;
t15 = -qJD(6) * t84 - t203 * t67 + t207 * t66;
t220 = t15 * t176 - t85 * t193;
t219 = t137 * t193 + t183 * t66;
t212 = qJD(1) ^ 2;
t130 = pkin(5) * t163 + t196;
t94 = pkin(5) * t136 + t142;
t73 = pkin(4) * t162 - pkin(5) * t234;
t40 = pkin(5) * t67 + t109;
t21 = -pkin(5) * t213 + t79;
t18 = t320 + t344;
t17 = t253 + t345;
t7 = -pkin(10) * t67 + t225;
t6 = pkin(5) * t281 - pkin(10) * t66 - qJD(5) * t318 + t254;
t4 = -t14 * t203 + t316;
t1 = [0, 0, 0, 0, 0.2e1 * t210 * t193, -0.2e1 * t286 * t272, t294, -t295, 0, -t189 * t294 + t206 * t231, t189 * t295 + t210 * t231, t116 * t297 + t327 * t162 (-t160 * t209 - t162 * t205) * t280 + (-t313 - t117 * t209 + (t160 * t205 - t162 * t209) * qJD(4)) * t206, t248 - t325, t242 + t312 + (-t221 - t306) * qJD(3) (-t187 - t284) * t281 -(-t154 * t279 + t289) * t187 + ((t160 * t189 + t311) * qJD(3) + (t209 * t247 + t300) * qJD(4) + t246) * t210 + (t122 * t278 + t189 * t117 + t309 + ((-t189 * t298 + t139) * qJD(1) + t71) * qJD(3)) * t206, t290 * t187 + (-t247 * t279 + (t162 * t189 + t310) * qJD(3) + t268) * t210 + (-t122 * t279 + t189 * t116 + t308 + (-t287 * qJD(1) - t189 * t304 - t72) * qJD(3)) * t206, -t137 * t36 - t234 * t66, -t101 * t66 - t136 * t36 - t137 * t213 + t234 * t67, -t219 + t251, t226 + t314 (-t183 - t284) * t281, -t254 * t183 - t255 * t210 + t109 * t101 - t142 * t213 + t79 * t136 + t98 * t67 + (t183 * t318 + t20 * t210) * qJD(5) + (qJD(1) * t252 + t19) * t281, t225 * t183 - t250 * t210 - t109 * t234 + t142 * t36 - t79 * t137 + t98 * t66 + (-t318 * qJD(1) - t20) * t281, -t15 * t235 + t8 * t85, -t15 * t51 + t16 * t235 + t215 * t85 - t8 * t84, -t220 + t259, t227 + t323 (-t176 - t284) * t281 -(-t203 * t7 + t207 * t6) * t176 - t267 * t210 + t40 * t51 - t94 * t215 + t21 * t84 + t57 * t16 + (t176 * t237 + t210 * t5) * qJD(6) + ((-t203 * t30 + t207 * t29) * qJD(1) + t4) * t281, -t12 * t210 + t57 * t15 + t21 * t85 - t40 * t235 + t94 * t8 + ((-qJD(6) * t30 + t6) * t176 + t2 * t210) * t203 + ((qJD(6) * t29 + t7) * t176 + t257 * t210) * t207 + (-qJD(1) * t237 - t5) * t281; 0, 0, 0, 0, 0, 0, 0, 0, 0, -t295, -t294, 0, 0, 0, 0, 0, t242 - t312 + (-t221 + t306) * qJD(3), t248 + t325, 0, 0, 0, 0, 0, -t226 + t314, t219 + t251, 0, 0, 0, 0, 0, -t227 + t323, t220 + t259; 0, 0, 0, 0, -t206 * t212 * t210, t286 * t212, 0, 0, 0, qJD(3) * t133 - t178 * t285 - t125, -t178 * t284, -t162 * t304 + t313 (t116 + t307) * t209 + (-t117 + t305) * t205, -t187 * t278 + (t187 * t296 + (-t162 + t282) * t206) * qJD(1), t187 * t279 + (-t187 * t298 + (t160 + t273) * t206) * qJD(1), t187 * t285, -pkin(3) * t117 - t308 + t245 * t187 - t133 * t160 + (pkin(8) * t304 + t311) * qJD(4) + (-t71 * t206 + (-pkin(8) * t281 - t122 * t210) * t205) * qJD(1), -pkin(3) * t116 + t309 - t291 * t187 - t133 * t162 + (-pkin(8) * t187 * t205 + t310) * qJD(4) + (-t122 * t296 + (-pkin(8) * t273 + t72) * t206) * qJD(1), t164 * t36 - t234 * t293, -t293 * t101 - t163 * t36 + t164 * t213 + t234 * t292, -t293 * t183 + (qJD(3) * t164 + t234) * t285, t292 * t183 + (-qJD(3) * t163 + t101) * t285, t183 * t285, t79 * t163 - t196 * t213 + t292 * t98 + (t181 * t276 + (-qJD(5) * t180 + t341) * t204 + t330) * t183 + t239 * t101 + (qJD(3) * t244 - t19) * t285, t79 * t164 + t196 * t36 + t293 * t98 + t328 * t183 - t239 * t234 + (-t288 * qJD(3) + t20) * t285, t106 * t8 - t235 * t322, -t105 * t8 + t106 * t215 + t235 * t321 - t322 * t51, -t322 * t176 + (qJD(3) * t106 + t235) * t285, t321 * t176 + (-qJD(3) * t105 + t51) * t285, t176 * t285, t21 * t105 - t130 * t215 + t321 * t57 + t317 * t51 + (t203 * t240 + t207 * t241) * t176 + ((-t203 * t89 + t207 * t88) * qJD(3) - t4) * t285, t21 * t106 + t130 * t8 + t322 * t57 - t317 * t235 + (-t203 * t241 + t207 * t240) * t176 + (-(t203 * t88 + t207 * t89) * qJD(3) + t5) * t285; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t162 * t160, -t160 ^ 2 + t162 ^ 2, t116 - t307, -t117 - t305, t193, -t122 * t162 - t187 * t72 + t214, t122 * t160 - t187 * t71 - t223, -t340, t338, t336, t331, t193, t253 * t183 + (-t101 * t162 + t183 * t277 + t193 * t208) * pkin(4) + t334, -t320 * t183 + (t162 * t234 + t183 * t276 - t193 * t204) * pkin(4) + t335, -t342, t339, t337, t332, t193, t228 * t193 + (t17 * t207 - t18 * t203) * t176 - t73 * t51 - (-t203 * t208 - t301) * t270 + (t176 * t229 - t5) * qJD(6) + t347, -t229 * t193 - t207 * t3 - (t17 * t203 + t18 * t207) * t176 + t73 * t235 + (t207 * t208 - t302) * t270 + (t176 * t228 - t316) * qJD(6) + t348; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t340, t338, t336, t331, t193, -t183 * t20 + t334, -t183 * t19 + t335, -t342, t339, t337, t332, t193 (-t13 * t203 - t315) * t176 + (t176 * t275 + t193 * t207 + t234 * t51) * pkin(5) + t333 (t14 * t176 - t2) * t203 + (-t13 * t176 - t257) * t207 + (t176 * t274 - t193 * t203 - t234 * t235) * pkin(5) + t349; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t342, t339, t337, t332, t193, -t176 * t5 + t333, -t176 * t4 - t207 * t257 + t348;];
tauc_reg  = t1;
