% Calculate inertial parameters regressor of inverse dynamics base forces vector with Newton-Euler for
% S4RRRR1
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [4x1]
%   Generalized joint coordinates (joint angles)
% qJD [4x1]
%   Generalized joint velocities
% qJDD [4x1]
%   Generalized joint accelerations
% g [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [7x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,d1,d2,d3,d4]';
% 
% Output:
% tauB_reg [6x(5*10)]
%   inertial parameter regressor of inverse dynamics base forces vector

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-31 17:22
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauB_reg = S4RRRR1_invdynB_fixb_reg2_snew_vp(qJ, qJD, qJDD, g, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(4,1),zeros(4,1),zeros(4,1),zeros(3,1),zeros(7,1)}
assert(isreal(qJ) && all(size(qJ) == [4 1]), ...
  'S4RRRR1_invdynB_fixb_reg2_snew_vp: qJ has to be [4x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [4 1]), ...
  'S4RRRR1_invdynB_fixb_reg2_snew_vp: qJD has to be [4x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [4 1]), ...
  'S4RRRR1_invdynB_fixb_reg2_snew_vp: qJDD has to be [4x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S4RRRR1_invdynB_fixb_reg2_snew_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [7 1]), ...
  'S4RRRR1_invdynB_fixb_reg2_snew_vp: pkin has to be [7x1] (double)');

%% Symbolic Calculation
% From invdyn_fixb_NewtonEuler_tauB_regressor_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 17:22:15
% EndTime: 2019-12-31 17:22:20
% DurationCPUTime: 2.99s
% Computational Cost: add. (9381->278), mult. (12125->396), div. (0->0), fcn. (7328->8), ass. (0->190)
t289 = qJD(1) + qJD(2);
t284 = qJD(3) + t289;
t282 = t284 ^ 2;
t297 = cos(qJ(3));
t288 = qJDD(1) + qJDD(2);
t283 = qJDD(3) + t288;
t293 = sin(qJ(3));
t323 = t293 * t283;
t249 = t297 * t282 + t323;
t318 = t297 * t283;
t252 = t293 * t282 - t318;
t294 = sin(qJ(2));
t298 = cos(qJ(2));
t203 = t298 * t249 - t294 * t252;
t236 = pkin(6) * t249 - t297 * g(3);
t340 = pkin(6) * t252 - t293 * g(3);
t167 = pkin(5) * t203 + t298 * t236 - t294 * t340;
t208 = t294 * t249 + t298 * t252;
t295 = sin(qJ(1));
t299 = cos(qJ(1));
t175 = t295 * t203 + t299 * t208;
t348 = pkin(5) * t208 + t294 * t236 + t298 * t340;
t357 = pkin(4) * t175 + t295 * t167 + t299 * t348;
t339 = t299 * t203 - t295 * t208;
t356 = pkin(4) * t339 + t299 * t167 - t295 * t348;
t276 = t295 * g(1) - t299 * g(2);
t266 = qJDD(1) * pkin(1) + t276;
t277 = t299 * g(1) + t295 * g(2);
t302 = qJD(1) ^ 2;
t267 = -t302 * pkin(1) - t277;
t222 = t294 * t266 + t298 * t267;
t287 = t289 ^ 2;
t213 = -t287 * pkin(2) + t222;
t304 = t298 * t266 - t294 * t267;
t303 = t288 * pkin(2) + t304;
t178 = t293 * t213 - t297 * t303;
t179 = t297 * t213 + t293 * t303;
t311 = t293 * t178 + t297 * t179;
t139 = t297 * t178 - t293 * t179;
t317 = t298 * t139;
t119 = -t294 * t311 + t317;
t322 = t294 * t139;
t342 = t298 * t311 + t322;
t109 = t295 * t119 + t299 * t342;
t353 = t299 * t119 - t295 * t342;
t259 = t298 * t287 + t294 * t288;
t262 = t294 * t287 - t298 * t288;
t216 = t295 * t259 + t299 * t262;
t240 = pkin(5) * t259 - t298 * g(3);
t341 = pkin(5) * t262 - t294 * g(3);
t350 = pkin(4) * t216 + t295 * t240 + t299 * t341;
t305 = t299 * t259 - t295 * t262;
t349 = pkin(4) * t305 + t299 * t240 - t295 * t341;
t310 = t298 * t222 - t294 * t304;
t183 = -t294 * t222 - t298 * t304;
t316 = t299 * t183;
t343 = -t295 * t310 + t316;
t321 = t295 * t183;
t142 = t299 * t310 + t321;
t292 = sin(qJ(4));
t290 = t292 ^ 2;
t328 = t290 * t282;
t169 = -t283 * pkin(3) - t282 * pkin(7) + t178;
t327 = t292 * t169;
t296 = cos(qJ(4));
t272 = t296 * t282 * t292;
t263 = qJDD(4) + t272;
t326 = t292 * t263;
t264 = qJDD(4) - t272;
t325 = t292 * t264;
t324 = t292 * t283;
t320 = t296 * t169;
t319 = t296 * t264;
t278 = t296 * t283;
t291 = t296 ^ 2;
t315 = t290 + t291;
t314 = qJD(4) * t284;
t313 = t292 * t314;
t312 = t296 * t314;
t170 = -t282 * pkin(3) + t283 * pkin(7) + t179;
t159 = -t292 * g(3) + t296 * t170;
t158 = t296 * g(3) + t292 * t170;
t133 = t292 * t158 + t296 * t159;
t232 = -t295 * t276 - t299 * t277;
t308 = t293 * t272;
t307 = t297 * t272;
t274 = t299 * qJDD(1) - t295 * t302;
t306 = -pkin(4) * t274 - t295 * g(3);
t132 = t296 * t158 - t292 * t159;
t231 = t299 * t276 - t295 * t277;
t301 = qJD(4) ^ 2;
t300 = pkin(1) * g(3);
t279 = t291 * t282;
t273 = t295 * qJDD(1) + t299 * t302;
t271 = -t279 - t301;
t270 = t279 - t301;
t269 = -t301 - t328;
t268 = t301 - t328;
t256 = t296 * t263;
t255 = -pkin(4) * t273 + t299 * g(3);
t254 = t279 - t328;
t253 = t279 + t328;
t247 = t315 * t283;
t246 = t278 - 0.2e1 * t313;
t245 = t278 - t313;
t244 = t312 + t324;
t243 = 0.2e1 * t312 + t324;
t242 = t315 * t314;
t230 = t293 * qJDD(4) + t297 * t242;
t229 = -t297 * qJDD(4) + t293 * t242;
t228 = -t292 * t269 - t319;
t227 = -t292 * t268 + t256;
t226 = t296 * t271 - t326;
t225 = t296 * t270 - t325;
t224 = t296 * t269 - t325;
t223 = t292 * t271 + t256;
t220 = t296 * t244 - t290 * t314;
t219 = -t292 * t245 - t291 * t314;
t206 = t297 * t247 - t293 * t253;
t202 = t293 * t247 + t297 * t253;
t201 = -t292 * t243 + t296 * t246;
t200 = t297 * t227 + t292 * t323;
t199 = t297 * t225 + t278 * t293;
t198 = t293 * t227 - t292 * t318;
t197 = t293 * t225 - t296 * t318;
t196 = t297 * t220 - t308;
t195 = t297 * t219 + t308;
t194 = t293 * t220 + t307;
t193 = t293 * t219 - t307;
t192 = t297 * t228 + t293 * t243;
t191 = t297 * t226 - t293 * t246;
t190 = t293 * t228 - t297 * t243;
t189 = t293 * t226 + t297 * t246;
t188 = -t294 * t229 + t298 * t230;
t187 = t298 * t229 + t294 * t230;
t186 = t297 * t201 - t293 * t254;
t185 = t293 * t201 + t297 * t254;
t180 = pkin(5) * t310 + t300;
t174 = -t294 * t202 + t298 * t206;
t171 = t298 * t202 + t294 * t206;
t163 = -t294 * t198 + t298 * t200;
t162 = -t294 * t197 + t298 * t199;
t161 = t298 * t198 + t294 * t200;
t160 = t298 * t197 + t294 * t199;
t156 = -t294 * t194 + t298 * t196;
t155 = -t294 * t193 + t298 * t195;
t154 = t298 * t194 + t294 * t196;
t153 = t298 * t193 + t294 * t195;
t152 = -t294 * t190 + t298 * t192;
t151 = -t294 * t189 + t298 * t191;
t150 = t298 * t190 + t294 * t192;
t149 = t298 * t189 + t294 * t191;
t148 = -pkin(7) * t224 + t320;
t147 = -pkin(7) * t223 + t327;
t146 = -pkin(3) * t224 + t159;
t145 = -pkin(3) * t223 + t158;
t144 = -t294 * t185 + t298 * t186;
t143 = t298 * t185 + t294 * t186;
t136 = pkin(2) * g(3) + pkin(6) * t311;
t135 = -t295 * t171 + t299 * t174;
t134 = t299 * t171 + t295 * t174;
t130 = -t295 * t150 + t299 * t152;
t129 = -t295 * t149 + t299 * t151;
t128 = t299 * t150 + t295 * t152;
t127 = t299 * t149 + t295 * t151;
t126 = -pkin(6) * t202 + t297 * t132;
t125 = pkin(6) * t206 + t293 * t132;
t124 = t297 * t133 + t293 * t169;
t123 = t293 * t133 - t297 * t169;
t122 = -pkin(6) * t190 - t293 * t146 + t297 * t148;
t121 = -pkin(6) * t189 - t293 * t145 + t297 * t147;
t116 = -pkin(2) * t224 + pkin(6) * t192 + t297 * t146 + t293 * t148;
t115 = -pkin(2) * t223 + pkin(6) * t191 + t297 * t145 + t293 * t147;
t114 = -pkin(5) * t171 - t294 * t125 + t298 * t126;
t113 = pkin(5) * t174 + t298 * t125 + t294 * t126;
t112 = -t294 * t123 + t298 * t124;
t111 = t298 * t123 + t294 * t124;
t110 = pkin(5) * t119 + pkin(6) * t317 - t294 * t136;
t107 = pkin(5) * t342 + pkin(6) * t322 + t298 * t136 + t300;
t106 = -pkin(6) * t123 - (pkin(3) * t293 - pkin(7) * t297) * t132;
t105 = -pkin(5) * t150 - t294 * t116 + t298 * t122;
t104 = -pkin(5) * t149 - t294 * t115 + t298 * t121;
t103 = -pkin(1) * t224 + pkin(5) * t152 + t298 * t116 + t294 * t122;
t102 = -pkin(1) * t223 + pkin(5) * t151 + t298 * t115 + t294 * t121;
t101 = pkin(6) * t124 - (-pkin(3) * t297 - pkin(7) * t293 - pkin(2)) * t132;
t100 = -t295 * t111 + t299 * t112;
t99 = t299 * t111 + t295 * t112;
t98 = -pkin(5) * t111 - t294 * t101 + t298 * t106;
t97 = pkin(1) * t132 + pkin(5) * t112 + t298 * t101 + t294 * t106;
t1 = [0, 0, 0, 0, 0, 0, 0, 0, 0, -g(1), 0, 0, 0, 0, 0, 0, -t273, -t274, 0, t232, 0, 0, 0, 0, 0, 0, -t305, t216, 0, t142, 0, 0, 0, 0, 0, 0, -t339, t175, 0, t109, 0, 0, 0, 0, 0, 0, t129, t130, t135, t100; 0, 0, 0, 0, 0, 0, 0, 0, 0, -g(2), 0, 0, 0, 0, 0, 0, t274, -t273, 0, t231, 0, 0, 0, 0, 0, 0, -t216, -t305, 0, -t343, 0, 0, 0, 0, 0, 0, -t175, -t339, 0, -t353, 0, 0, 0, 0, 0, 0, t127, t128, t134, t99; 0, 0, 0, 0, 0, 0, 0, 0, 0, -g(3), 0, 0, 0, 0, 0, 0, 0, 0, 0, -g(3), 0, 0, 0, 0, 0, 0, 0, 0, 0, -g(3), 0, 0, 0, 0, 0, 0, 0, 0, 0, -g(3), 0, 0, 0, 0, 0, 0, t223, t224, 0, -t132; 0, 0, 0, 0, 0, 0, 0, -g(3), g(2), 0, 0, 0, t274, 0, -t273, 0, t306, -t255, -t231, -pkin(4) * t231, 0, 0, -t216, 0, -t305, 0, t350, t349, t343, pkin(4) * t343 + pkin(5) * t316 - t295 * t180, 0, 0, -t175, 0, -t339, 0, t357, t356, t353, pkin(4) * t353 - t295 * t107 + t299 * t110, -t295 * t154 + t299 * t156, -t295 * t143 + t299 * t144, -t295 * t161 + t299 * t163, -t295 * t153 + t299 * t155, -t295 * t160 + t299 * t162, -t295 * t187 + t299 * t188, -pkin(4) * t127 - t295 * t102 + t299 * t104, -pkin(4) * t128 - t295 * t103 + t299 * t105, -pkin(4) * t134 - t295 * t113 + t299 * t114, -pkin(4) * t99 - t295 * t97 + t299 * t98; 0, 0, 0, 0, 0, 0, g(3), 0, -g(1), 0, 0, 0, t273, 0, t274, 0, t255, t306, t232, pkin(4) * t232, 0, 0, t305, 0, -t216, 0, -t349, t350, t142, pkin(4) * t142 + pkin(5) * t321 + t299 * t180, 0, 0, t339, 0, -t175, 0, -t356, t357, t109, pkin(4) * t109 + t299 * t107 + t295 * t110, t299 * t154 + t295 * t156, t299 * t143 + t295 * t144, t299 * t161 + t295 * t163, t299 * t153 + t295 * t155, t299 * t160 + t295 * t162, t299 * t187 + t295 * t188, pkin(4) * t129 + t299 * t102 + t295 * t104, pkin(4) * t130 + t299 * t103 + t295 * t105, pkin(4) * t135 + t299 * t113 + t295 * t114, pkin(4) * t100 + t295 * t98 + t299 * t97; 0, 0, 0, 0, 0, 0, -g(2), g(1), 0, 0, 0, 0, 0, 0, 0, qJDD(1), t276, t277, 0, 0, 0, 0, 0, 0, 0, t288, -pkin(1) * t262 + t304, -pkin(1) * t259 - t222, 0, -pkin(1) * t183, 0, 0, 0, 0, 0, t283, -pkin(1) * t208 - pkin(2) * t252 - t178, -pkin(1) * t203 - pkin(2) * t249 - t179, 0, -pkin(1) * t119 - pkin(2) * t139, (t244 + t312) * t292, t296 * t243 + t292 * t246, t296 * t268 + t326, (t245 - t313) * t296, t292 * t270 + t319, 0, pkin(1) * t149 + pkin(2) * t189 + pkin(3) * t246 + pkin(7) * t226 - t320, pkin(1) * t150 + pkin(2) * t190 - pkin(3) * t243 + pkin(7) * t228 + t327, pkin(1) * t171 + pkin(2) * t202 + pkin(3) * t253 + pkin(7) * t247 + t133, pkin(1) * t111 + pkin(2) * t123 - pkin(3) * t169 + pkin(7) * t133;];
tauB_reg = t1;
