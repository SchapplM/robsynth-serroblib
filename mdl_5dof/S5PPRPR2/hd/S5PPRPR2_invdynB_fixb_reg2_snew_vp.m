% Calculate inertial parameters regressor of inverse dynamics base forces vector with Newton-Euler for
% S5PPRPR2
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
%   pkin=[a2,a3,a4,a5,d3,d5,theta1,theta2]';
% 
% Output:
% tauB_reg [6x(6*10)]
%   inertial parameter regressor of inverse dynamics base forces vector

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-05 15:03
% Revision: 77da58f92bca3eff71542919beafa37024070d86 (2019-12-05)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauB_reg = S5PPRPR2_invdynB_fixb_reg2_snew_vp(qJ, qJD, qJDD, g, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(5,1),zeros(3,1),zeros(8,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5PPRPR2_invdynB_fixb_reg2_snew_vp: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5PPRPR2_invdynB_fixb_reg2_snew_vp: qJD has to be [5x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [5 1]), ...
  'S5PPRPR2_invdynB_fixb_reg2_snew_vp: qJDD has to be [5x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S5PPRPR2_invdynB_fixb_reg2_snew_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S5PPRPR2_invdynB_fixb_reg2_snew_vp: pkin has to be [8x1] (double)');

%% Symbolic Calculation
% From invdyn_fixb_NewtonEuler_tauB_regressor_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-05 15:03:23
% EndTime: 2019-12-05 15:03:30
% DurationCPUTime: 3.95s
% Computational Cost: add. (5879->298), mult. (9956->436), div. (0->0), fcn. (7134->8), ass. (0->201)
t259 = sin(pkin(7));
t261 = cos(pkin(7));
t241 = t259 * g(1) - t261 * g(2);
t234 = -qJDD(2) + t241;
t265 = sin(qJ(3));
t269 = qJD(3) ^ 2;
t267 = cos(qJ(3));
t291 = t267 * qJDD(3);
t237 = -t265 * t269 + t291;
t195 = pkin(5) * t237 + t265 * t234;
t258 = sin(pkin(8));
t260 = cos(pkin(8));
t293 = t265 * qJDD(3);
t238 = t267 * t269 + t293;
t276 = t258 * t237 + t260 * t238;
t287 = -pkin(5) * t238 + t267 * t234;
t132 = -qJ(2) * t276 - t258 * t195 + t260 * t287;
t348 = t259 * t132;
t347 = t261 * t132;
t278 = t260 * t237 - t258 * t238;
t325 = qJ(2) * t278 + t260 * t195 + t258 * t287;
t346 = t259 * t325;
t345 = t261 * t325;
t343 = t259 * t276;
t342 = t259 * t278;
t341 = t261 * t276;
t340 = t261 * t278;
t242 = t261 * g(1) + t259 * g(2);
t299 = g(3) - qJDD(1);
t217 = -t258 * t242 + t260 * t299;
t218 = -t260 * t242 - t258 * t299;
t157 = t267 * t217 + t265 * t218;
t158 = -t265 * t217 + t267 * t218;
t279 = t265 * t157 + t267 * t158;
t121 = t267 * t157 - t265 * t158;
t313 = t260 * t121;
t99 = -t258 * t279 + t313;
t319 = t258 * t121;
t100 = t260 * t279 + t319;
t281 = -pkin(1) * t278 - pkin(2) * t237 - qJ(1) * t276 + t157;
t339 = pkin(1) * t276 + pkin(2) * t238 - qJ(1) * t278;
t270 = (2 * qJD(4) * qJD(3)) + t158;
t296 = qJDD(3) * qJ(4);
t150 = -t269 * pkin(3) + t270 + t296;
t257 = qJDD(3) * pkin(3);
t151 = -t269 * qJ(4) + qJDD(4) + t157 - t257;
t115 = t265 * t150 - t267 * t151;
t280 = t267 * t150 + t265 * t151;
t91 = -t258 * t115 + t260 * t280;
t89 = t260 * t115 + t258 * t280;
t332 = t259 * t299;
t329 = t261 * t299;
t324 = pkin(3) + pkin(6);
t264 = sin(qJ(5));
t255 = t264 ^ 2;
t321 = t255 * t269;
t266 = cos(qJ(5));
t256 = t266 ^ 2;
t320 = t256 * t269;
t314 = t259 * t234;
t310 = t260 * t234;
t219 = t261 * t234;
t148 = -t269 * pkin(6) + t150;
t307 = t264 * t148;
t290 = t266 * t269 * t264;
t243 = qJDD(5) + t290;
t306 = t264 * t243;
t298 = t255 + t256;
t236 = t298 * qJDD(3);
t305 = t265 * t236;
t304 = t266 * t148;
t303 = t266 * t243;
t244 = qJDD(5) - t290;
t302 = t266 * t244;
t300 = t267 * t236;
t297 = qJD(3) * qJD(5);
t295 = t261 * qJDD(3);
t294 = t264 * qJDD(3);
t292 = t266 * qJDD(3);
t289 = t264 * t297;
t288 = t266 * t297;
t149 = -qJDD(3) * pkin(6) + t151;
t141 = -t266 * t149 - t264 * t234;
t142 = t264 * t149 - t266 * t234;
t103 = -t266 * t141 + t264 * t142;
t239 = t298 * t269;
t191 = -t265 * t239 + t300;
t192 = -t267 * t239 - t305;
t145 = t260 * t191 + t258 * t192;
t146 = -t258 * t191 + t260 * t192;
t285 = -pkin(1) * t145 - pkin(2) * t191 + qJ(1) * t146 + qJ(4) * t239 - t324 * t236 + t103;
t284 = -t270 - 0.2e1 * t296 - t339;
t283 = -qJDD(4) + 0.2e1 * t257 - t281;
t282 = t158 + t339;
t156 = t258 * t217 + t260 * t218;
t194 = -t259 * t241 - t261 * t242;
t268 = qJD(5) ^ 2;
t274 = -t268 - t320;
t273 = t265 * t290;
t272 = t267 * t290;
t104 = t264 * t141 + t266 * t142;
t155 = t260 * t217 - t258 * t218;
t193 = t261 * t241 - t259 * t242;
t250 = t259 * qJDD(3);
t247 = t268 - t320;
t246 = -t268 - t321;
t245 = -t268 + t321;
t240 = (-t255 + t256) * t269;
t233 = -0.2e1 * t289 + t292;
t232 = -t289 + t292;
t231 = -t288 - t294;
t230 = 0.2e1 * t288 + t294;
t229 = pkin(1) * t234;
t228 = t264 * t244;
t227 = t298 * t297;
t216 = t267 * qJDD(5) - t265 * t227;
t215 = t265 * qJDD(5) + t267 * t227;
t214 = -t264 * t232 - t256 * t297;
t213 = -t266 * t231 - t255 * t297;
t211 = -t264 * t274 - t303;
t210 = t264 * t247 - t302;
t209 = (-t232 + t289) * t266;
t208 = t266 * t246 - t228;
t207 = -t266 * t245 + t306;
t206 = t266 * t274 - t306;
t205 = -t266 * t247 - t228;
t204 = t264 * t246 + t302;
t203 = -t264 * t245 - t303;
t202 = (t231 - t288) * t264;
t182 = t266 * t230 + t264 * t233;
t181 = t264 * t230 - t266 * t233;
t172 = -t265 * t213 - t272;
t171 = -t265 * t214 + t272;
t170 = t267 * t213 - t273;
t169 = t267 * t214 + t273;
t168 = -t265 * t205 + t266 * t291;
t167 = -t265 * t203 - t264 * t291;
t166 = t267 * t205 + t265 * t292;
t165 = t267 * t203 - t264 * t293;
t164 = t265 * t206 + t267 * t233;
t163 = t265 * t204 + t267 * t230;
t162 = -t267 * t206 + t265 * t233;
t161 = -t267 * t204 + t265 * t230;
t160 = -t265 * t181 + t267 * t240;
t159 = t267 * t181 + t265 * t240;
t153 = -t258 * t215 + t260 * t216;
t144 = t261 * t156 - t314;
t143 = t259 * t156 + t219;
t140 = -t258 * t170 + t260 * t172;
t139 = -t258 * t169 + t260 * t171;
t138 = -t258 * t166 + t260 * t168;
t137 = -t258 * t165 + t260 * t167;
t129 = -t258 * t162 + t260 * t164;
t128 = -t258 * t161 + t260 * t163;
t127 = t260 * t162 + t258 * t164;
t126 = t260 * t161 + t258 * t163;
t123 = -t258 * t159 + t260 * t160;
t118 = pkin(2) * t234 + pkin(5) * t279;
t113 = t261 * t129 + t259 * t211;
t112 = t261 * t128 + t259 * t208;
t111 = t259 * t129 - t261 * t211;
t110 = t259 * t128 - t261 * t208;
t109 = pkin(4) * t206 - qJ(4) * t211 - t142;
t108 = pkin(4) * t230 - t324 * t208 + t304;
t107 = pkin(4) * t204 - qJ(4) * t208 - t141;
t106 = pkin(4) * t233 - t324 * t211 - t307;
t105 = -pkin(5) * t115 + (-pkin(3) * t265 + qJ(4) * t267) * t234;
t102 = pkin(5) * t280 + (pkin(3) * t267 + qJ(4) * t265 + pkin(2)) * t234;
t101 = -pkin(4) * t239 - t104;
t97 = -pkin(4) * t300 - pkin(5) * t191 - t265 * t101;
t96 = -pkin(4) * t305 + pkin(5) * t192 + t267 * t101;
t95 = t265 * t103 + t267 * t148;
t94 = -t267 * t103 + t265 * t148;
t93 = t261 * t100 - t314;
t92 = t259 * t100 + t219;
t88 = -pkin(1) * t127 - pkin(2) * t162 - qJ(4) * t233 + t324 * t206 - t304;
t87 = -pkin(1) * t126 - pkin(2) * t161 - qJ(4) * t230 + t324 * t204 - t307;
t86 = t261 * t91 - t314;
t85 = t259 * t91 + t219;
t83 = -pkin(5) * t162 - t265 * t106 + t267 * t109;
t82 = -pkin(5) * t161 + t267 * t107 - t265 * t108;
t81 = pkin(1) * t99 + pkin(2) * t121;
t80 = pkin(4) * t103 - qJ(4) * t104;
t79 = -pkin(2) * t211 + pkin(5) * t164 + t267 * t106 + t265 * t109;
t78 = -pkin(2) * t208 + pkin(5) * t163 + t265 * t107 + t267 * t108;
t77 = pkin(4) * t148 - t324 * t104;
t76 = pkin(5) * t313 + qJ(2) * t99 - t258 * t118;
t75 = -t258 * t94 + t260 * t95;
t74 = t258 * t95 + t260 * t94;
t73 = -pkin(1) * t89 - pkin(2) * t115 + pkin(3) * t151 - qJ(4) * t150;
t72 = -qJ(2) * t145 - t258 * t96 + t260 * t97;
t71 = -qJ(2) * t89 - t258 * t102 + t260 * t105;
t70 = t259 * t104 + t261 * t75;
t69 = -t261 * t104 + t259 * t75;
t68 = -qJ(2) * t127 - t258 * t79 + t260 * t83;
t67 = -qJ(2) * t126 - t258 * t78 + t260 * t82;
t66 = -pkin(5) * t94 - t265 * t77 + t267 * t80;
t65 = -pkin(1) * t74 - pkin(2) * t94 - qJ(4) * t148 + t324 * t103;
t64 = -pkin(2) * t104 + pkin(5) * t95 + t265 * t80 + t267 * t77;
t63 = -qJ(2) * t74 - t258 * t64 + t260 * t66;
t1 = [0, 0, 0, 0, 0, 0, 0, 0, 0, -g(1), 0, 0, 0, 0, 0, 0, 0, 0, 0, t194, 0, 0, 0, 0, 0, 0, 0, 0, 0, t144, 0, 0, 0, 0, 0, 0, -t341, -t340, 0, t93, 0, 0, 0, 0, 0, 0, 0, t341, t340, t86, 0, 0, 0, 0, 0, 0, t112, t113, t261 * t146, t70; 0, 0, 0, 0, 0, 0, 0, 0, 0, -g(2), 0, 0, 0, 0, 0, 0, 0, 0, 0, t193, 0, 0, 0, 0, 0, 0, 0, 0, 0, t143, 0, 0, 0, 0, 0, 0, -t343, -t342, 0, t92, 0, 0, 0, 0, 0, 0, 0, t343, t342, t85, 0, 0, 0, 0, 0, 0, t110, t111, t259 * t146, t69; 0, 0, 0, 0, 0, 0, 0, 0, 0, -g(3), 0, 0, 0, 0, 0, 0, 0, 0, 0, -t299, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t155, 0, 0, 0, 0, 0, 0, t278, -t276, 0, -t99, 0, 0, 0, 0, 0, 0, 0, -t278, t276, t89, 0, 0, 0, 0, 0, 0, t126, t127, t145, t74; 0, 0, 0, 0, 0, 0, 0, -g(3), g(2), 0, 0, 0, 0, 0, 0, 0, -t332, -t329, -t193, -qJ(1) * t193, 0, 0, 0, 0, 0, 0, -t259 * t217 - t258 * t219, -t259 * t218 - t260 * t219, t261 * t155, -qJ(1) * t143 - (pkin(1) * t259 - qJ(2) * t261) * t155, 0, 0, t340, 0, -t341, t250, -t281 * t259 - t345, -t282 * t259 - t347, t261 * t99, -qJ(1) * t92 - t259 * t81 + t261 * t76, t250, -t340, t341, 0, 0, 0, -t261 * t89, -t259 * t283 + t345, -t259 * t284 + t347, -qJ(1) * t85 - t259 * t73 + t261 * t71, t261 * t139 - t259 * t209, t261 * t123 - t259 * t182, t261 * t138 - t259 * t210, t261 * t140 - t259 * t202, t261 * t137 - t259 * t207, t261 * t153, -qJ(1) * t110 - t259 * t87 + t261 * t67, -qJ(1) * t111 - t259 * t88 + t261 * t68, -t259 * t285 + t261 * t72, -qJ(1) * t69 - t259 * t65 + t261 * t63; 0, 0, 0, 0, 0, 0, g(3), 0, -g(1), 0, 0, 0, 0, 0, 0, 0, t329, -t332, t194, qJ(1) * t194, 0, 0, 0, 0, 0, 0, t261 * t217 - t258 * t314, t261 * t218 - t259 * t310, t259 * t155, qJ(1) * t144 - (-pkin(1) * t261 - qJ(2) * t259) * t155, 0, 0, t342, 0, -t343, -t295, t281 * t261 - t346, t282 * t261 - t348, t259 * t99, qJ(1) * t93 + t259 * t76 + t261 * t81, -t295, -t342, t343, 0, 0, 0, -t259 * t89, t261 * t283 + t346, t261 * t284 + t348, qJ(1) * t86 + t259 * t71 + t261 * t73, t259 * t139 + t261 * t209, t259 * t123 + t261 * t182, t259 * t138 + t261 * t210, t259 * t140 + t261 * t202, t259 * t137 + t261 * t207, t259 * t153, qJ(1) * t112 + t259 * t67 + t261 * t87, qJ(1) * t113 + t259 * t68 + t261 * t88, t259 * t72 + t261 * t285, qJ(1) * t70 + t259 * t63 + t261 * t65; 0, 0, 0, 0, 0, 0, -g(2), g(1), 0, 0, 0, 0, 0, 0, 0, 0, t241, t242, 0, 0, 0, 0, 0, 0, 0, 0, t310, -t258 * t234, t156, qJ(2) * t156 + t229, 0, 0, t276, 0, t278, 0, t132, -t325, t100, pkin(5) * t319 + qJ(2) * t100 + t260 * t118 + t229, 0, -t276, -t278, 0, 0, 0, t91, -t132, t325, qJ(2) * t91 + t260 * t102 + t258 * t105 + t229, t260 * t169 + t258 * t171, t260 * t159 + t258 * t160, t260 * t166 + t258 * t168, t260 * t170 + t258 * t172, t260 * t165 + t258 * t167, t260 * t215 + t258 * t216, -pkin(1) * t208 + qJ(2) * t128 + t258 * t82 + t260 * t78, -pkin(1) * t211 + qJ(2) * t129 + t258 * t83 + t260 * t79, qJ(2) * t146 + t258 * t97 + t260 * t96, -pkin(1) * t104 + qJ(2) * t75 + t258 * t66 + t260 * t64;];
tauB_reg = t1;
