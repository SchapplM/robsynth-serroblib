% Calculate inertial parameters regressor of inverse dynamics base forces vector with Newton-Euler for
% S4RRRP2
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
%   pkin=[a2,a3,a4,d1,d2,d3]';
% 
% Output:
% tauB_reg [6x(5*10)]
%   inertial parameter regressor of inverse dynamics base forces vector

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-31 17:13
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauB_reg = S4RRRP2_invdynB_fixb_reg2_snew_vp(qJ, qJD, qJDD, g, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(4,1),zeros(4,1),zeros(4,1),zeros(3,1),zeros(6,1)}
assert(isreal(qJ) && all(size(qJ) == [4 1]), ...
  'S4RRRP2_invdynB_fixb_reg2_snew_vp: qJ has to be [4x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [4 1]), ...
  'S4RRRP2_invdynB_fixb_reg2_snew_vp: qJD has to be [4x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [4 1]), ...
  'S4RRRP2_invdynB_fixb_reg2_snew_vp: qJDD has to be [4x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S4RRRP2_invdynB_fixb_reg2_snew_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [6 1]), ...
  'S4RRRP2_invdynB_fixb_reg2_snew_vp: pkin has to be [6x1] (double)');

%% Symbolic Calculation
% From invdyn_fixb_NewtonEuler_tauB_regressor_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 17:13:11
% EndTime: 2019-12-31 17:13:14
% DurationCPUTime: 2.21s
% Computational Cost: add. (5384->270), mult. (8010->361), div. (0->0), fcn. (4360->6), ass. (0->207)
t307 = qJD(1) + qJD(2);
t305 = t307 ^ 2;
t314 = cos(qJ(2));
t306 = qJDD(1) + qJDD(2);
t311 = sin(qJ(2));
t348 = t311 * t306;
t275 = t314 * t305 + t348;
t342 = t314 * t306;
t278 = t311 * t305 - t342;
t312 = sin(qJ(1));
t315 = cos(qJ(1));
t228 = t312 * t275 + t315 * t278;
t255 = pkin(5) * t275 - t314 * g(3);
t371 = pkin(5) * t278 - t311 * g(3);
t376 = pkin(4) * t228 + t312 * t255 + t315 * t371;
t320 = t315 * t275 - t312 * t278;
t375 = pkin(4) * t320 + t315 * t255 - t312 * t371;
t294 = t315 * g(1) + t312 * g(2);
t317 = qJD(1) ^ 2;
t282 = -t317 * pkin(1) - t294;
t293 = t312 * g(1) - t315 * g(2);
t319 = qJDD(1) * pkin(1) + t293;
t232 = t311 * t282 - t314 * t319;
t233 = t314 * t282 + t311 * t319;
t325 = t311 * t232 + t314 * t233;
t191 = t314 * t232 - t311 * t233;
t341 = t315 * t191;
t372 = -t312 * t325 + t341;
t347 = t312 * t191;
t152 = t315 * t325 + t347;
t310 = sin(qJ(3));
t313 = cos(qJ(3));
t292 = t313 * t305 * t310;
t337 = qJDD(3) + t292;
t370 = t337 * pkin(3);
t219 = -t305 * pkin(2) + t306 * pkin(6) + t233;
t206 = t313 * g(3) + t310 * t219;
t338 = qJD(3) * t313;
t328 = t307 * t338;
t349 = t310 * t306;
t266 = t328 + t349;
t257 = t266 * qJ(4);
t367 = -t257 - t206 + t370;
t316 = qJD(3) ^ 2;
t308 = t310 ^ 2;
t354 = t308 * t305;
t289 = -t316 - t354;
t285 = qJDD(3) - t292;
t350 = t310 * t285;
t241 = t313 * t289 - t350;
t366 = pkin(2) * t241;
t309 = t313 ^ 2;
t299 = t309 * t305;
t291 = -t299 - t316;
t351 = t310 * t337;
t243 = t313 * t291 - t351;
t339 = qJD(3) * t307;
t329 = t310 * t339;
t343 = t313 * t306;
t268 = -0.2e1 * t329 + t343;
t202 = t311 * t243 + t314 * t268;
t204 = t314 * t243 - t311 * t268;
t161 = t315 * t202 + t312 * t204;
t365 = pkin(4) * t161;
t344 = t313 * t285;
t245 = -t310 * t289 - t344;
t265 = 0.2e1 * t328 + t349;
t203 = t311 * t245 - t314 * t265;
t205 = t314 * t245 + t311 * t265;
t162 = t315 * t203 + t312 * t205;
t364 = pkin(4) * t162;
t340 = t308 + t309;
t273 = t340 * t306;
t279 = t299 + t354;
t224 = t311 * t273 + t314 * t279;
t227 = t314 * t273 - t311 * t279;
t182 = t315 * t224 + t312 * t227;
t363 = pkin(4) * t182;
t362 = pkin(5) * t202;
t361 = pkin(5) * t203;
t360 = pkin(5) * t224;
t274 = t313 * t337;
t239 = t310 * t291 + t274;
t359 = pkin(6) * t239;
t358 = pkin(6) * t241;
t355 = t307 * t310;
t170 = (qJ(4) * t338 - 0.2e1 * qJD(4) * t310) * t307 + t367;
t353 = t310 * t170;
t218 = -t306 * pkin(2) - t305 * pkin(6) + t232;
t352 = t310 * t218;
t346 = t313 * t170;
t345 = t313 * t218;
t336 = 0.2e1 * qJD(4) * t307;
t335 = t310 * t348;
t334 = t310 * t342;
t333 = pkin(1) * t202 + pkin(2) * t268 + pkin(6) * t243;
t332 = pkin(1) * t203 - pkin(2) * t265 + pkin(6) * t245;
t330 = pkin(1) * t224 + pkin(2) * t279 + pkin(6) * t273;
t327 = -pkin(1) * t239 + pkin(5) * t204;
t326 = -pkin(1) * t241 + pkin(5) * t205;
t207 = -t310 * g(3) + t313 * t219;
t167 = t310 * t206 + t313 * t207;
t251 = -t312 * t293 - t315 * t294;
t324 = t311 * t292;
t323 = t314 * t292;
t184 = -pkin(2) * t239 + t206;
t287 = t315 * qJDD(1) - t312 * t317;
t322 = -pkin(4) * t287 - t312 * g(3);
t166 = t313 * t206 - t310 * t207;
t250 = t315 * t293 - t312 * t294;
t267 = -t329 + t343;
t283 = qJD(3) * pkin(3) - qJ(4) * t355;
t318 = t267 * qJ(4) - qJD(3) * t283 + t313 * t336 + t207;
t180 = -t267 * pkin(3) - qJ(4) * t299 + t283 * t355 + qJDD(4) + t218;
t296 = t310 * t336;
t290 = t299 - t316;
t288 = t316 - t354;
t286 = t312 * qJDD(1) + t315 * t317;
t280 = t299 - t354;
t263 = -pkin(4) * t286 + t315 * g(3);
t262 = t340 * t339;
t249 = t311 * qJDD(3) + t314 * t262;
t248 = -t314 * qJDD(3) + t311 * t262;
t247 = t313 * t266 - t308 * t339;
t246 = -t310 * t267 - t309 * t339;
t244 = -t310 * t288 + t274;
t242 = t313 * t290 - t350;
t240 = t313 * t288 + t351;
t238 = t310 * t290 + t344;
t235 = (t266 + t328) * t310;
t234 = (t267 - t329) * t313;
t231 = -pkin(3) * t265 - qJ(4) * t285;
t222 = pkin(5) * t227;
t221 = -t310 * t265 + t313 * t268;
t220 = t313 * t265 + t310 * t268;
t217 = t314 * t244 + t335;
t216 = t314 * t242 + t311 * t343;
t215 = t311 * t244 - t334;
t214 = t311 * t242 - t313 * t342;
t213 = t314 * t247 - t324;
t212 = t314 * t246 + t324;
t211 = t311 * t247 + t323;
t210 = t311 * t246 - t323;
t196 = -t312 * t248 + t315 * t249;
t195 = t315 * t248 + t312 * t249;
t194 = t314 * t221 - t311 * t280;
t193 = t311 * t221 + t314 * t280;
t188 = pkin(1) * g(3) + pkin(5) * t325;
t187 = t345 - t358;
t186 = t352 - t359;
t185 = t207 - t366;
t183 = -t312 * t224 + t315 * t227;
t181 = pkin(4) * t183;
t179 = -t312 * t215 + t315 * t217;
t178 = -t312 * t214 + t315 * t216;
t177 = t315 * t215 + t312 * t217;
t176 = t315 * t214 + t312 * t216;
t175 = -t312 * t211 + t315 * t213;
t174 = -t312 * t210 + t315 * t212;
t173 = t315 * t211 + t312 * t213;
t172 = t315 * t210 + t312 * t212;
t171 = -pkin(3) * t299 + t318;
t169 = -qJ(4) * t289 + t180;
t168 = t296 + (-t328 + t349) * qJ(4) - t367;
t164 = -t312 * t203 + t315 * t205;
t163 = -t312 * t202 + t315 * t204;
t160 = pkin(4) * t164;
t159 = pkin(4) * t163;
t158 = pkin(3) * t268 + qJ(4) * t291 - t180;
t157 = qJ(4) * t343 + (t279 - t299) * pkin(3) + t318;
t156 = -t312 * t193 + t315 * t194;
t155 = t315 * t193 + t312 * t194;
t154 = -t366 + (-t289 - t299) * pkin(3) + t318;
t153 = -qJ(4) * t328 + t184 + t257 + t296 - 0.2e1 * t370;
t150 = t314 * t166 - t360;
t149 = t311 * t166 + t222;
t148 = t314 * t167 + t311 * t218;
t147 = t311 * t167 - t314 * t218;
t146 = -qJ(4) * t274 - t310 * t158 - t359;
t145 = t313 * t169 - t310 * t231 - t358;
t144 = -pkin(3) * t180 + qJ(4) * t171;
t143 = t313 * t171 - t353;
t142 = t310 * t171 + t346;
t141 = -t311 * t185 + t314 * t187 - t361;
t140 = -t311 * t184 + t314 * t186 - t362;
t139 = -t310 * t157 + t313 * t168;
t138 = t314 * t185 + t311 * t187 + t326;
t137 = t314 * t184 + t311 * t186 + t327;
t136 = -pkin(3) * t335 + t314 * t139 - t360;
t135 = pkin(3) * t334 + t311 * t139 + t222;
t134 = t314 * t143 + t311 * t180;
t133 = t311 * t143 - t314 * t180;
t132 = -pkin(2) * t142 - pkin(3) * t170;
t131 = -t312 * t147 + t315 * t148;
t130 = t315 * t147 + t312 * t148;
t129 = t314 * t145 - t311 * t154 - t361;
t128 = t314 * t146 - t311 * t153 - t362;
t127 = -pkin(5) * t147 - (pkin(2) * t311 - pkin(6) * t314) * t166;
t126 = t311 * t145 + t314 * t154 + t326;
t125 = t311 * t146 + t314 * t153 + t327;
t124 = pkin(5) * t148 - (-pkin(2) * t314 - pkin(6) * t311 - pkin(1)) * t166;
t123 = -pkin(6) * t142 - qJ(4) * t346 - t310 * t144;
t122 = -t312 * t133 + t315 * t134;
t121 = t315 * t133 + t312 * t134;
t120 = -pkin(5) * t133 + t314 * t123 - t311 * t132;
t119 = -pkin(1) * t142 + pkin(5) * t134 + t311 * t123 + t314 * t132;
t1 = [0, 0, 0, 0, 0, 0, 0, 0, 0, -g(1), 0, 0, 0, 0, 0, 0, -t286, -t287, 0, t251, 0, 0, 0, 0, 0, 0, -t320, t228, 0, t152, 0, 0, 0, 0, 0, 0, t163, t164, t183, t131, 0, 0, 0, 0, 0, 0, t163, t164, t183, t122; 0, 0, 0, 0, 0, 0, 0, 0, 0, -g(2), 0, 0, 0, 0, 0, 0, t287, -t286, 0, t250, 0, 0, 0, 0, 0, 0, -t228, -t320, 0, -t372, 0, 0, 0, 0, 0, 0, t161, t162, t182, t130, 0, 0, 0, 0, 0, 0, t161, t162, t182, t121; 0, 0, 0, 0, 0, 0, 0, 0, 0, -g(3), 0, 0, 0, 0, 0, 0, 0, 0, 0, -g(3), 0, 0, 0, 0, 0, 0, 0, 0, 0, -g(3), 0, 0, 0, 0, 0, 0, t239, t241, 0, -t166, 0, 0, 0, 0, 0, 0, t239, t241, 0, t142; 0, 0, 0, 0, 0, 0, 0, -g(3), g(2), 0, 0, 0, t287, 0, -t286, 0, t322, -t263, -t250, -pkin(4) * t250, 0, 0, -t228, 0, -t320, 0, t376, t375, t372, pkin(4) * t372 + pkin(5) * t341 - t312 * t188, t175, t156, t179, t174, t178, t196, -t312 * t137 + t315 * t140 - t365, -t312 * t138 + t315 * t141 - t364, -t312 * t149 + t315 * t150 - t363, -pkin(4) * t130 - t312 * t124 + t315 * t127, t175, t156, t179, t174, t178, t196, -t312 * t125 + t315 * t128 - t365, -t312 * t126 + t315 * t129 - t364, -t312 * t135 + t315 * t136 - t363, -pkin(4) * t121 - t312 * t119 + t315 * t120; 0, 0, 0, 0, 0, 0, g(3), 0, -g(1), 0, 0, 0, t286, 0, t287, 0, t263, t322, t251, pkin(4) * t251, 0, 0, t320, 0, -t228, 0, -t375, t376, t152, pkin(4) * t152 + pkin(5) * t347 + t315 * t188, t173, t155, t177, t172, t176, t195, t315 * t137 + t312 * t140 + t159, t315 * t138 + t312 * t141 + t160, t315 * t149 + t312 * t150 + t181, pkin(4) * t131 + t315 * t124 + t312 * t127, t173, t155, t177, t172, t176, t195, t315 * t125 + t312 * t128 + t159, t315 * t126 + t312 * t129 + t160, t315 * t135 + t312 * t136 + t181, pkin(4) * t122 + t315 * t119 + t312 * t120; 0, 0, 0, 0, 0, 0, -g(2), g(1), 0, 0, 0, 0, 0, 0, 0, qJDD(1), t293, t294, 0, 0, 0, 0, 0, 0, 0, t306, -pkin(1) * t278 - t232, -pkin(1) * t275 - t233, 0, -pkin(1) * t191, t235, t220, t240, t234, t238, 0, t333 - t345, t332 + t352, t167 + t330, pkin(1) * t147 - pkin(2) * t218 + pkin(6) * t167, t235, t220, t240, t234, t238, 0, -qJ(4) * t351 + t313 * t158 + t333, t310 * t169 + t313 * t231 + t332, t313 * t157 + t310 * t168 + t330, pkin(1) * t133 - pkin(2) * t180 + pkin(6) * t143 - qJ(4) * t353 + t313 * t144;];
tauB_reg = t1;
