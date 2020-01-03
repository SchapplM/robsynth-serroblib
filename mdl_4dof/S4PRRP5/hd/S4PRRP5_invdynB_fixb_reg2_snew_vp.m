% Calculate inertial parameters regressor of inverse dynamics base forces vector with Newton-Euler for
% S4PRRP5
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
% pkin [6x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,d2,d3,theta1]';
% 
% Output:
% tauB_reg [6x(5*10)]
%   inertial parameter regressor of inverse dynamics base forces vector

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-31 16:29
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauB_reg = S4PRRP5_invdynB_fixb_reg2_snew_vp(qJ, qJD, qJDD, g, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(4,1),zeros(4,1),zeros(4,1),zeros(3,1),zeros(6,1)}
assert(isreal(qJ) && all(size(qJ) == [4 1]), ...
  'S4PRRP5_invdynB_fixb_reg2_snew_vp: qJ has to be [4x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [4 1]), ...
  'S4PRRP5_invdynB_fixb_reg2_snew_vp: qJD has to be [4x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [4 1]), ...
  'S4PRRP5_invdynB_fixb_reg2_snew_vp: qJDD has to be [4x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S4PRRP5_invdynB_fixb_reg2_snew_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [6 1]), ...
  'S4PRRP5_invdynB_fixb_reg2_snew_vp: pkin has to be [6x1] (double)');

%% Symbolic Calculation
% From invdyn_fixb_NewtonEuler_tauB_regressor_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 16:29:23
% EndTime: 2019-12-31 16:29:26
% DurationCPUTime: 1.42s
% Computational Cost: add. (3008->256), mult. (6419->355), div. (0->0), fcn. (3750->6), ass. (0->201)
t268 = sin(qJ(3));
t270 = cos(qJ(3));
t273 = qJD(2) ^ 2;
t251 = t268 * t273 * t270;
t289 = qJDD(3) + t251;
t324 = t289 * pkin(3);
t265 = sin(pkin(6));
t266 = cos(pkin(6));
t242 = t265 * g(1) - t266 * g(2);
t226 = t266 * t242;
t243 = t266 * g(1) + t265 * g(2);
t194 = -t265 * t243 + t226;
t263 = g(3) - qJDD(1);
t269 = sin(qJ(2));
t271 = cos(qJ(2));
t221 = -t271 * t243 - t269 * t263;
t197 = -t273 * pkin(2) + qJDD(2) * pkin(5) + t221;
t168 = t268 * t197 + t270 * t242;
t294 = qJD(2) * qJD(3);
t284 = t270 * t294;
t293 = t268 * qJDD(2);
t234 = t284 + t293;
t224 = t234 * qJ(4);
t323 = -t224 - t168 + t324;
t272 = qJD(3) ^ 2;
t261 = t268 ^ 2;
t313 = t261 * t273;
t248 = -t272 - t313;
t246 = qJDD(3) - t251;
t303 = t268 * t246;
t210 = t270 * t248 - t303;
t322 = pkin(2) * t210;
t262 = t270 ^ 2;
t259 = t262 * t273;
t250 = -t259 - t272;
t304 = t268 * t289;
t212 = t270 * t250 - t304;
t286 = t268 * t294;
t291 = t270 * qJDD(2);
t236 = -0.2e1 * t286 + t291;
t169 = t269 * t212 + t271 * t236;
t321 = pkin(4) * t169;
t299 = t270 * t246;
t214 = -t268 * t248 - t299;
t233 = 0.2e1 * t284 + t293;
t170 = t269 * t214 - t271 * t233;
t320 = pkin(4) * t170;
t296 = t261 + t262;
t237 = t296 * qJDD(2);
t240 = t259 + t313;
t192 = t269 * t237 + t271 * t240;
t319 = pkin(4) * t192;
t300 = t270 * t289;
t208 = t268 * t250 + t300;
t318 = pkin(5) * t208;
t317 = pkin(5) * t210;
t172 = t271 * t212 - t269 * t236;
t141 = t265 * t172 - t266 * t208;
t316 = qJ(1) * t141;
t173 = t271 * t214 + t269 * t233;
t142 = t265 * t173 - t266 * t210;
t315 = qJ(1) * t142;
t193 = t271 * t237 - t269 * t240;
t314 = qJ(1) * t193;
t290 = t271 * qJDD(2);
t239 = -t269 * t273 + t290;
t312 = t265 * t239;
t311 = t265 * t242;
t309 = t265 * t263;
t186 = t266 * t193;
t308 = t266 * t239;
t307 = t266 * t263;
t137 = (qJ(4) * qJD(3) * t270 - 0.2e1 * qJD(4) * t268) * qJD(2) + t323;
t306 = t268 * t137;
t220 = -t269 * t243 + t271 * t263;
t196 = -qJDD(2) * pkin(2) - t273 * pkin(5) + t220;
t305 = t268 * t196;
t302 = t270 * t137;
t301 = t270 * t196;
t298 = -pkin(1) * t208 + pkin(4) * t172;
t297 = -pkin(1) * t210 + pkin(4) * t173;
t171 = t270 * t197 - t268 * t242;
t295 = qJD(2) * t268;
t292 = t269 * qJDD(2);
t288 = 0.2e1 * qJD(2) * qJD(4);
t285 = t268 * t292;
t283 = t268 * t290;
t238 = t271 * t273 + t292;
t282 = pkin(1) * t239 + qJ(1) * t238 - t220;
t281 = -pkin(1) * t238 + qJ(1) * t239 - t221;
t165 = t269 * t220 + t271 * t221;
t195 = -t266 * t243 - t311;
t280 = t269 * t251;
t279 = t271 * t251;
t146 = -pkin(2) * t208 + t168;
t201 = pkin(4) * t238 - t271 * t242;
t200 = -pkin(4) * t239 - t269 * t242;
t131 = t270 * t168 - t268 * t171;
t132 = t268 * t168 + t270 * t171;
t164 = t271 * t220 - t269 * t221;
t235 = -t286 + t291;
t244 = qJD(3) * pkin(3) - qJ(4) * t295;
t277 = t235 * qJ(4) - qJD(3) * t244 + t270 * t288 + t171;
t276 = -pkin(1) * t169 - pkin(2) * t236 - pkin(5) * t212;
t275 = -pkin(1) * t170 + pkin(2) * t233 - pkin(5) * t214;
t274 = -pkin(1) * t192 - pkin(2) * t240 - pkin(5) * t237;
t158 = -t235 * pkin(3) - qJ(4) * t259 + t244 * t295 + qJDD(4) + t196;
t253 = t268 * t288;
t249 = t259 - t272;
t247 = t272 - t313;
t241 = t259 - t313;
t232 = t296 * t294;
t223 = t266 * t238;
t222 = t265 * t238;
t219 = t269 * qJDD(3) + t271 * t232;
t218 = t270 * t234 - t261 * t294;
t217 = -t271 * qJDD(3) + t269 * t232;
t216 = -t268 * t235 - t262 * t294;
t213 = -t268 * t247 + t300;
t211 = t270 * t249 - t303;
t209 = -t270 * t247 - t304;
t207 = -t268 * t249 - t299;
t206 = (-t234 - t284) * t268;
t205 = (-t235 + t286) * t270;
t204 = t266 * t219;
t203 = t265 * t219;
t202 = -pkin(3) * t233 - qJ(4) * t246;
t191 = pkin(4) * t193;
t188 = -t268 * t233 + t270 * t236;
t187 = -t270 * t233 - t268 * t236;
t185 = t265 * t193;
t182 = qJ(1) * t186;
t181 = t271 * t218 - t280;
t180 = t271 * t216 + t280;
t179 = t269 * t218 + t279;
t178 = t269 * t216 - t279;
t177 = t271 * t213 + t285;
t176 = t271 * t211 + t269 * t291;
t175 = t269 * t213 - t283;
t174 = t269 * t211 - t270 * t290;
t162 = t271 * t188 - t269 * t241;
t161 = t269 * t188 + t271 * t241;
t160 = t301 - t317;
t159 = t305 - t318;
t157 = t266 * t165 - t311;
t156 = t265 * t165 + t226;
t155 = t266 * t181 - t265 * t206;
t154 = t266 * t180 - t265 * t205;
t153 = t265 * t181 + t266 * t206;
t152 = t265 * t180 + t266 * t205;
t151 = t266 * t177 - t265 * t209;
t150 = t266 * t176 - t265 * t207;
t149 = t265 * t177 + t266 * t209;
t148 = t265 * t176 + t266 * t207;
t147 = t171 - t322;
t145 = -qJ(4) * t248 + t158;
t144 = t266 * t173 + t265 * t210;
t143 = t266 * t172 + t265 * t208;
t140 = qJ(1) * t144;
t139 = qJ(1) * t143;
t138 = -pkin(3) * t259 + t277;
t136 = pkin(3) * t236 + qJ(4) * t250 - t158;
t135 = t266 * t162 - t265 * t187;
t134 = t265 * t162 + t266 * t187;
t133 = t253 + (-t284 + t293) * qJ(4) - t323;
t129 = qJ(4) * t291 + (t240 - t259) * pkin(3) + t277;
t128 = t276 + t301;
t127 = t275 - t305;
t126 = -t322 + (-t248 - t259) * pkin(3) + t277;
t125 = -qJ(4) * t284 + t146 + t224 + t253 - 0.2e1 * t324;
t124 = t271 * t131 - t319;
t123 = t271 * t132 + t269 * t196;
t122 = t269 * t132 - t271 * t196;
t121 = -qJ(4) * t300 - t268 * t136 - t318;
t120 = t270 * t145 - t268 * t202 - t317;
t119 = -pkin(3) * t158 + qJ(4) * t138;
t118 = -t132 + t274;
t117 = t270 * t138 - t306;
t116 = t268 * t138 + t302;
t115 = -t269 * t147 + t271 * t160 - t320;
t114 = -t269 * t146 + t271 * t159 - t321;
t113 = qJ(4) * t304 - t270 * t136 + t276;
t112 = -t268 * t145 - t270 * t202 + t275;
t111 = -t268 * t129 + t270 * t133;
t110 = t266 * t123 - t131 * t265;
t109 = t265 * t123 + t131 * t266;
t108 = -pkin(3) * t285 + t271 * t111 - t319;
t107 = t271 * t117 + t269 * t158;
t106 = t269 * t117 - t271 * t158;
t105 = -t270 * t129 - t268 * t133 + t274;
t104 = -pkin(1) * t122 + pkin(2) * t196 - pkin(5) * t132;
t103 = -pkin(2) * t116 - pkin(3) * t137;
t102 = t271 * t120 - t269 * t126 - t320;
t101 = t271 * t121 - t269 * t125 - t321;
t100 = -pkin(4) * t122 - (pkin(2) * t269 - pkin(5) * t271) * t131;
t99 = -pkin(5) * t116 - qJ(4) * t302 - t268 * t119;
t98 = t266 * t107 + t265 * t116;
t97 = t265 * t107 - t266 * t116;
t96 = -pkin(1) * t106 + pkin(2) * t158 - pkin(5) * t117 + qJ(4) * t306 - t270 * t119;
t95 = -pkin(4) * t106 - t269 * t103 + t271 * t99;
t1 = [0, 0, 0, 0, 0, 0, 0, 0, 0, -g(1), 0, 0, 0, 0, 0, 0, 0, 0, 0, t195, 0, 0, 0, 0, 0, 0, -t223, -t308, 0, t157, 0, 0, 0, 0, 0, 0, t143, t144, t186, t110, 0, 0, 0, 0, 0, 0, t143, t144, t186, t98; 0, 0, 0, 0, 0, 0, 0, 0, 0, -g(2), 0, 0, 0, 0, 0, 0, 0, 0, 0, t194, 0, 0, 0, 0, 0, 0, -t222, -t312, 0, t156, 0, 0, 0, 0, 0, 0, t141, t142, t185, t109, 0, 0, 0, 0, 0, 0, t141, t142, t185, t97; 0, 0, 0, 0, 0, 0, 0, 0, 0, -g(3), 0, 0, 0, 0, 0, 0, 0, 0, 0, -t263, 0, 0, 0, 0, 0, 0, t239, -t238, 0, -t164, 0, 0, 0, 0, 0, 0, t169, t170, t192, t122, 0, 0, 0, 0, 0, 0, t169, t170, t192, t106; 0, 0, 0, 0, 0, 0, 0, -g(3), g(2), 0, 0, 0, 0, 0, 0, 0, -t309, -t307, -t194, -qJ(1) * t194, 0, 0, t308, 0, -t223, t265 * qJDD(2), t266 * t200 + t282 * t265, t266 * t201 + t265 * t281, t266 * t164, -qJ(1) * t156 - (pkin(1) * t265 - pkin(4) * t266) * t164, t155, t135, t151, t154, t150, t204, t266 * t114 - t265 * t128 - t316, t266 * t115 - t265 * t127 - t315, t266 * t124 + (-t118 - t314) * t265, -qJ(1) * t109 + t266 * t100 - t265 * t104, t155, t135, t151, t154, t150, t204, t266 * t101 - t265 * t113 - t316, t266 * t102 - t265 * t112 - t315, t266 * t108 + (-t105 - t314) * t265, -qJ(1) * t97 - t265 * t96 + t266 * t95; 0, 0, 0, 0, 0, 0, g(3), 0, -g(1), 0, 0, 0, 0, 0, 0, 0, t307, -t309, t195, qJ(1) * t195, 0, 0, t312, 0, -t222, -t266 * qJDD(2), t265 * t200 - t282 * t266, t265 * t201 - t266 * t281, t265 * t164, qJ(1) * t157 - (-pkin(1) * t266 - pkin(4) * t265) * t164, t153, t134, t149, t152, t148, t203, t265 * t114 + t266 * t128 + t139, t265 * t115 + t266 * t127 + t140, t266 * t118 + t265 * t124 + t182, qJ(1) * t110 + t265 * t100 + t266 * t104, t153, t134, t149, t152, t148, t203, t265 * t101 + t266 * t113 + t139, t265 * t102 + t266 * t112 + t140, t266 * t105 + t265 * t108 + t182, qJ(1) * t98 + t265 * t95 + t266 * t96; 0, 0, 0, 0, 0, 0, -g(2), g(1), 0, 0, 0, 0, 0, 0, 0, 0, t242, t243, 0, 0, 0, 0, t238, 0, t239, 0, -t201, t200, t165, pkin(1) * t242 + pkin(4) * t165, t179, t161, t175, t178, t174, t217, t271 * t146 + t269 * t159 + t298, t271 * t147 + t269 * t160 + t297, t269 * t131 + t191, pkin(4) * t123 - (-pkin(2) * t271 - pkin(5) * t269 - pkin(1)) * t131, t179, t161, t175, t178, t174, t217, t269 * t121 + t271 * t125 + t298, t269 * t120 + t271 * t126 + t297, pkin(3) * t283 + t269 * t111 + t191, -pkin(1) * t116 + pkin(4) * t107 + t271 * t103 + t269 * t99;];
tauB_reg = t1;
