% Calculate minimal parameter regressor of inverse dynamics joint torque vector for
% S5RRRRP8
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
% pkin [8x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d1,d2,d3,d4]';
% 
% Output:
% tau_reg [5x28]
%   minimal parameter regressor of inverse dynamics joint torque vector

% Quelle: HybrDyn-Toolbox
% Datum: 2021-01-16 00:22
% Revision: 24b2e7d74a0c1a3b64fa2f8f5ad758691ad61af3 (2021-01-15)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tau_reg = S5RRRRP8_invdynJ_fixb_regmin_slag_vp(qJ, qJD, qJDD, g, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(5,1),zeros(3,1),zeros(8,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RRRRP8_invdynJ_fixb_regmin_slag_vp: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RRRRP8_invdynJ_fixb_regmin_slag_vp: qJD has to be [5x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [5 1]), ...
  'S5RRRRP8_invdynJ_fixb_regmin_slag_vp: qJDD has to be [5x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S5RRRRP8_invdynJ_fixb_regmin_slag_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S5RRRRP8_invdynJ_fixb_regmin_slag_vp: pkin has to be [8x1] (double)');

%% Symbolic Calculation
% From invdyn_joint_fixb_regressor_minpar_matlab.m
% OptimizationMode: 2
% StartTime: 2021-01-16 00:21:20
% EndTime: 2021-01-16 00:21:35
% DurationCPUTime: 4.33s
% Computational Cost: add. (5123->453), mult. (11504->588), div. (0->0), fcn. (7987->10), ass. (0->234)
t189 = sin(qJ(2));
t192 = cos(qJ(2));
t225 = pkin(2) * t189 - pkin(7) * t192;
t129 = t225 * qJD(1);
t188 = sin(qJ(3));
t106 = t188 * t129;
t318 = pkin(7) + pkin(8);
t251 = qJD(3) * t318;
t191 = cos(qJ(3));
t283 = t189 * t191;
t284 = t188 * t192;
t342 = -t188 * t251 - t106 - (-pkin(6) * t283 - pkin(8) * t284) * qJD(1);
t280 = t191 * t192;
t217 = pkin(3) * t189 - pkin(8) * t280;
t268 = qJD(1) * t189;
t244 = t188 * t268;
t271 = pkin(6) * t244 + t191 * t129;
t341 = qJD(1) * t217 + t191 * t251 + t271;
t258 = t191 * qJD(2);
t121 = -t244 + t258;
t259 = t188 * qJD(2);
t122 = t191 * t268 + t259;
t187 = sin(qJ(4));
t262 = qJD(3) * t189;
t239 = qJD(1) * t262;
t256 = qJD(1) * qJD(2);
t240 = t192 * t256;
t254 = t189 * qJDD(1);
t334 = qJD(2) * qJD(3) + t240 + t254;
t229 = t334 * t188 + t191 * t239;
t205 = t191 * qJDD(2) - t229;
t314 = cos(qJ(4));
t241 = t314 * qJD(4);
t260 = qJD(4) * t187;
t55 = (qJDD(2) - t239) * t188 + t334 * t191;
t16 = -t121 * t241 + t122 * t260 - t187 * t205 - t314 * t55;
t257 = t192 * qJD(1);
t164 = -qJD(3) + t257;
t152 = -qJD(4) + t164;
t63 = -t314 * t121 + t187 * t122;
t296 = t152 * t63;
t340 = -t16 - t296;
t212 = t187 * t121 + t314 * t122;
t319 = t212 ^ 2;
t61 = t63 ^ 2;
t339 = -t61 + t319;
t286 = t187 * t188;
t211 = t314 * t191 - t286;
t325 = qJD(3) + qJD(4);
t326 = t314 * qJD(3) + t241;
t298 = -t326 * t191 + t211 * t257 + t325 * t286;
t124 = t187 * t191 + t314 * t188;
t77 = t325 * t124;
t297 = -t124 * t257 + t77;
t338 = t212 * t63;
t337 = t63 * qJ(5);
t17 = qJD(4) * t212 + t187 * t55 - t314 * t205;
t292 = t212 * t152;
t336 = -t17 - t292;
t245 = t192 * t259;
t261 = qJD(3) * t191;
t335 = t189 * t261 + t245;
t186 = qJ(3) + qJ(4);
t180 = cos(t186);
t178 = t192 * qJDD(1);
t238 = t189 * t256;
t327 = -t238 + t178;
t118 = qJDD(3) - t327;
t132 = t225 * qJD(2);
t138 = -t192 * pkin(2) - t189 * pkin(7) - pkin(1);
t78 = qJD(1) * t132 + qJDD(1) * t138;
t69 = t191 * t78;
t113 = t138 * qJD(1);
t175 = pkin(6) * t257;
t143 = qJD(2) * pkin(7) + t175;
t73 = t188 * t113 + t191 * t143;
t96 = t327 * pkin(6) + qJDD(2) * pkin(7);
t11 = t118 * pkin(3) - t55 * pkin(8) - qJD(3) * t73 - t188 * t96 + t69;
t263 = qJD(3) * t188;
t208 = t113 * t261 - t143 * t263 + t188 * t78 + t191 * t96;
t15 = pkin(8) * t205 + t208;
t72 = t191 * t113 - t188 * t143;
t44 = -t122 * pkin(8) + t72;
t38 = -t164 * pkin(3) + t44;
t45 = t121 * pkin(8) + t73;
t237 = -t187 * t11 - t314 * t15 - t38 * t241 + t45 * t260;
t310 = g(3) * t189;
t179 = sin(t186);
t193 = cos(qJ(1));
t279 = t193 * t179;
t190 = sin(qJ(1));
t281 = t190 * t192;
t88 = -t180 * t281 + t279;
t278 = t193 * t180;
t90 = t190 * t179 + t192 * t278;
t204 = g(1) * t90 - g(2) * t88 + t180 * t310 + t237;
t142 = -qJD(2) * pkin(2) + pkin(6) * t268;
t80 = -t121 * pkin(3) + t142;
t333 = t63 * t80 + t204;
t112 = qJDD(4) + t118;
t313 = pkin(3) * t152;
t332 = -t187 * pkin(3) * t112 + t241 * t313;
t331 = qJ(5) * t212;
t307 = t188 * pkin(3);
t330 = pkin(3) * t263 - t257 * t307 - t175;
t144 = t318 * t188;
t145 = t318 * t191;
t272 = -t187 * t144 + t314 * t145;
t329 = t272 * qJD(4) + t342 * t187 + t341 * t314;
t328 = -t144 * t241 - t145 * t260 - t341 * t187 + t342 * t314;
t224 = g(1) * t193 + g(2) * t190;
t309 = g(3) * t192;
t203 = -t189 * t224 + t309;
t173 = pkin(6) * t254;
t290 = qJDD(2) * pkin(2);
t97 = pkin(6) * t240 + t173 - t290;
t324 = -qJD(3) * pkin(7) * t164 + t203 + t97;
t295 = t16 * qJ(5);
t98 = t112 * pkin(4);
t323 = -t212 * qJD(5) + t295 + t98;
t87 = t179 * t281 + t278;
t282 = t190 * t180;
t89 = -t192 * t279 + t282;
t322 = -g(1) * t89 + g(2) * t87 + t179 * t310;
t43 = t314 * t45;
t19 = t187 * t38 + t43;
t243 = t314 * t11 - t187 * t15;
t200 = -qJD(4) * t19 + t243;
t196 = t200 + t322;
t321 = -t80 * t212 + t196;
t236 = -t63 * pkin(4) - qJD(5);
t34 = -t236 + t80;
t320 = -t34 * t212 + t322;
t41 = t187 * t45;
t18 = t314 * t38 - t41;
t9 = t18 - t331;
t6 = -pkin(4) * t152 + t9;
t317 = -t9 + t6;
t308 = t122 * pkin(3);
t181 = t189 * pkin(6);
t182 = t191 * pkin(3);
t306 = pkin(2) + t182;
t136 = pkin(4) * t179 + t307;
t305 = pkin(6) + t136;
t304 = -t297 * qJ(5) + qJD(5) * t211 + t328;
t303 = -pkin(4) * t268 + t298 * qJ(5) - t124 * qJD(5) - t329;
t302 = t314 * t44 - t41;
t301 = t297 * pkin(4) + t330;
t120 = t191 * t138;
t71 = -pkin(8) * t283 + t120 + (-pkin(6) * t188 - pkin(3)) * t192;
t166 = pkin(6) * t280;
t270 = t188 * t138 + t166;
t285 = t188 * t189;
t79 = -pkin(8) * t285 + t270;
t299 = t187 * t71 + t314 * t79;
t294 = t17 * qJ(5);
t293 = t55 * t188;
t291 = qJD(5) * t63;
t289 = t121 * t164;
t288 = t122 * t164;
t287 = t122 * t191;
t277 = t193 * t188;
t276 = t193 * t191;
t274 = t188 * t132 + t138 * t261;
t273 = t191 * t132 + t259 * t181;
t137 = pkin(4) * t180 + t182;
t133 = pkin(3) * t285 + t181;
t184 = t189 ^ 2;
t269 = -t192 ^ 2 + t184;
t267 = qJD(2) * t121;
t266 = qJD(2) * t122;
t265 = qJD(2) * t189;
t264 = qJD(2) * t192;
t81 = t335 * pkin(3) + pkin(6) * t264;
t250 = t164 * t258;
t249 = t164 * t263;
t248 = t188 * t262;
t247 = t164 * t261;
t235 = -t187 * t44 - t43;
t232 = -t187 * t79 + t314 * t71;
t231 = -qJD(3) * t113 - t96;
t230 = -t314 * t144 - t187 * t145;
t228 = t314 * t264;
t227 = -g(1) * t87 - g(2) * t89;
t226 = -g(1) * t88 - g(2) * t90;
t223 = g(1) * t190 - g(2) * t193;
t222 = t143 * t261 - t69;
t221 = -pkin(7) * t118 + qJD(3) * t142;
t128 = pkin(2) + t137;
t183 = -qJ(5) - t318;
t220 = t128 * t192 - t183 * t189;
t218 = -t180 * t309 + (g(1) * t278 + g(2) * t282) * t189;
t216 = pkin(1) + t220;
t215 = t224 * t179;
t214 = -0.2e1 * pkin(1) * t256 - pkin(6) * qJDD(2);
t29 = t217 * qJD(2) + (-t166 + (pkin(8) * t189 - t138) * t188) * qJD(3) + t273;
t31 = -t335 * pkin(8) + (-t189 * t258 - t192 * t263) * pkin(6) + t274;
t213 = t187 * t29 + t71 * t241 - t79 * t260 + t314 * t31;
t210 = t188 * t118 - t247;
t209 = t191 * t118 + t249;
t195 = qJD(1) ^ 2;
t206 = pkin(1) * t195 + t224;
t194 = qJD(2) ^ 2;
t202 = 0.2e1 * qJDD(1) * pkin(1) - pkin(6) * t194 + t223;
t201 = t204 + t294;
t199 = -t299 * qJD(4) - t187 * t31 + t314 * t29;
t37 = -pkin(3) * t205 + t97;
t5 = t17 * pkin(4) + qJDD(5) + t37;
t171 = t314 * pkin(3) + pkin(4);
t151 = t179 * t309;
t104 = t190 * t188 + t192 * t276;
t103 = t190 * t191 - t192 * t277;
t102 = -t190 * t280 + t277;
t101 = t188 * t281 + t276;
t94 = t211 * t189;
t93 = t124 * t189;
t84 = -pkin(4) * t211 - t306;
t74 = t93 * pkin(4) + t133;
t49 = qJ(5) * t211 + t272;
t48 = -t124 * qJ(5) + t230;
t40 = pkin(4) * t212 + t308;
t33 = t188 * t228 - t187 * t248 - t260 * t285 + (t187 * t264 + t326 * t189) * t191;
t32 = t187 * t245 + t77 * t189 - t191 * t228;
t26 = pkin(4) * t33 + t81;
t25 = -qJ(5) * t93 + t299;
t24 = -pkin(4) * t192 - qJ(5) * t94 + t232;
t13 = t302 - t331;
t12 = t235 + t337;
t10 = t19 - t337;
t4 = -qJ(5) * t33 - qJD(5) * t93 + t213;
t3 = pkin(4) * t265 + t32 * qJ(5) - t94 * qJD(5) + t199;
t2 = -t237 - t291 - t294;
t1 = t200 + t323;
t7 = [qJDD(1), t223, t224, qJDD(1) * t184 + 0.2e1 * t192 * t238, 0.2e1 * t189 * t178 - 0.2e1 * t269 * t256, qJDD(2) * t189 + t192 * t194, qJDD(2) * t192 - t189 * t194, 0, t189 * t214 + t192 * t202, -t189 * t202 + t192 * t214, t55 * t283 + (t192 * t258 - t248) * t122, (t121 * t191 - t122 * t188) * t264 + (t191 * t205 - t293 + (-t121 * t188 - t287) * qJD(3)) * t189, (-t55 - t250) * t192 + (t209 + t266) * t189, (t164 * t259 - t205) * t192 + (-t210 + t267) * t189, -t118 * t192 - t164 * t265, -(-t138 * t263 + t273) * t164 + t120 * t118 - g(1) * t102 - g(2) * t104 + ((t247 - t267) * pkin(6) + (-pkin(6) * t118 + qJD(2) * t142 - t231) * t188 + t222) * t192 + (-pkin(6) * t205 + t72 * qJD(2) + t142 * t261 + t97 * t188) * t189, t274 * t164 - t270 * t118 - g(1) * t101 - g(2) * t103 + (t142 * t258 + (-t249 + t266) * pkin(6) + t208) * t192 + (-t142 * t263 - t73 * qJD(2) + t97 * t191 + (t55 - t250) * pkin(6)) * t189, -t16 * t94 - t212 * t32, t16 * t93 - t17 * t94 - t212 * t33 + t32 * t63, t112 * t94 + t152 * t32 + t16 * t192 + t212 * t265, -t112 * t93 + t152 * t33 + t17 * t192 - t63 * t265, -t112 * t192 - t152 * t265, t112 * t232 + t133 * t17 - t152 * t199 + t18 * t265 - t192 * t200 + t80 * t33 + t37 * t93 + t81 * t63 + t226, -t299 * t112 - t133 * t16 + t213 * t152 - t19 * t265 - t237 * t192 + t212 * t81 - t80 * t32 + t37 * t94 + t227, -t1 * t192 + t112 * t24 - t152 * t3 + t17 * t74 + t26 * t63 + t6 * t265 + t33 * t34 + t5 * t93 + t226, -t10 * t265 - t112 * t25 + t152 * t4 - t16 * t74 + t192 * t2 + t212 * t26 - t32 * t34 + t5 * t94 + t227, -t1 * t94 - t10 * t33 + t24 * t16 - t25 * t17 + t189 * t223 - t2 * t93 - t212 * t3 + t6 * t32 - t4 * t63, t1 * t24 + t10 * t4 + t2 * t25 + t34 * t26 + t6 * t3 + t5 * t74 + (-g(1) * t305 - g(2) * t216) * t193 + (g(1) * t216 - g(2) * t305) * t190; 0, 0, 0, -t189 * t195 * t192, t269 * t195, t254, t178, qJDD(2), t189 * t206 - t173 - t309, t310 + (-pkin(6) * qJDD(1) + t206) * t192, -t164 * t287 + t293, (t55 - t289) * t191 + (t205 + t288) * t188, (-t122 * t189 + t164 * t280) * qJD(1) + t210, (-t121 * t189 - t164 * t284) * qJD(1) + t209, t164 * t268, -pkin(2) * t229 + t271 * t164 + t221 * t188 + (-t72 * t189 + (pkin(6) * t121 - t142 * t188) * t192) * qJD(1) + (t290 - t324) * t191, -pkin(2) * t55 - t106 * t164 + t221 * t191 + (-t142 * t280 + t73 * t189 + (-t122 * t192 + t164 * t283) * pkin(6)) * qJD(1) + t324 * t188, -t16 * t124 - t212 * t298, -t124 * t17 - t16 * t211 - t212 * t297 + t298 * t63, t124 * t112 + t298 * t152 - t212 * t268, t112 * t211 + t297 * t152 + t63 * t268, t152 * t268, t112 * t230 + t329 * t152 - t17 * t306 - t18 * t268 - t211 * t37 + t297 * t80 + t330 * t63 + t218, -t272 * t112 + t306 * t16 + t37 * t124 + t151 - t298 * t80 + t330 * t212 + t328 * t152 + (t19 * qJD(1) - t215) * t189, t48 * t112 - t303 * t152 + t84 * t17 - t211 * t5 - t6 * t268 + t297 * t34 + t301 * t63 + t218, -t49 * t112 + t5 * t124 - t84 * t16 + t151 + t301 * t212 - t298 * t34 + t304 * t152 + (qJD(1) * t10 - t215) * t189, -t1 * t124 - t297 * t10 + t48 * t16 - t49 * t17 - t224 * t192 + t2 * t211 - t212 * t303 + t298 * t6 - t304 * t63 - t310, -g(3) * t220 + t1 * t48 + t304 * t10 + t2 * t49 + t301 * t34 + t303 * t6 + t5 * t84 + t224 * (t128 * t189 + t183 * t192); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t122 * t121, -t121 ^ 2 + t122 ^ 2, t55 + t289, t205 - t288, t118, -g(1) * t103 + g(2) * t101 - t142 * t122 - t73 * t164 + (t231 + t310) * t188 - t222, g(1) * t104 - g(2) * t102 + g(3) * t283 - t121 * t142 - t164 * t72 - t208, t338, t339, t340, t336, t112, t235 * t152 + (t314 * t112 - t122 * t63 + t152 * t260) * pkin(3) + t321, -t302 * t152 - t212 * t308 + t332 + t333, t171 * t112 + t12 * t152 - t40 * t63 + (-t43 + (-t38 + t313) * t187) * qJD(4) + t243 + t320 + t323, -t13 * t152 - t212 * t40 + t34 * t63 + t201 + t291 + t332, t10 * t212 + t12 * t212 + t13 * t63 + t171 * t16 - t6 * t63 + (-t17 * t187 + (t187 * t212 - t314 * t63) * qJD(4)) * pkin(3), t1 * t171 - t10 * t13 - t6 * t12 - t34 * t40 - g(1) * (-t136 * t192 * t193 + t137 * t190) - g(2) * (-t136 * t281 - t137 * t193) + t136 * t310 + (t2 * t187 + (t314 * t10 - t187 * t6) * qJD(4)) * pkin(3); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t338, t339, t340, t336, t112, -t19 * t152 + t321, -t152 * t18 + t333, t295 - t10 * t152 + 0.2e1 * t98 + (t236 - t34) * t212 + t196, -t319 * pkin(4) - t9 * t152 + (qJD(5) + t34) * t63 + t201, t16 * pkin(4) - t317 * t63, t317 * t10 + (t1 + t320) * pkin(4); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t17 - t292, -t16 + t296, -t61 - t319, t10 * t63 + t212 * t6 + t203 + t5;];
tau_reg = t7;
