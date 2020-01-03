% Calculate inertial parameters regressor of inverse dynamics base forces vector with Newton-Euler for
% S5RPPRP4
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
% pkin [7x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d1,d4,theta3]';
% 
% Output:
% tauB_reg [6x(6*10)]
%   inertial parameter regressor of inverse dynamics base forces vector

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-31 17:52
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauB_reg = S5RPPRP4_invdynB_fixb_reg2_snew_vp(qJ, qJD, qJDD, g, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(5,1),zeros(3,1),zeros(7,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RPPRP4_invdynB_fixb_reg2_snew_vp: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RPPRP4_invdynB_fixb_reg2_snew_vp: qJD has to be [5x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [5 1]), ...
  'S5RPPRP4_invdynB_fixb_reg2_snew_vp: qJDD has to be [5x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S5RPPRP4_invdynB_fixb_reg2_snew_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [7 1]), ...
  'S5RPPRP4_invdynB_fixb_reg2_snew_vp: pkin has to be [7x1] (double)');

%% Symbolic Calculation
% From invdyn_fixb_NewtonEuler_tauB_regressor_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 17:52:22
% EndTime: 2019-12-31 17:52:27
% DurationCPUTime: 2.44s
% Computational Cost: add. (6054->303), mult. (11376->389), div. (0->0), fcn. (4972->6), ass. (0->211)
t304 = sin(pkin(7));
t305 = cos(pkin(7));
t312 = qJD(1) ^ 2;
t271 = -t304 * qJDD(1) + t305 * t312;
t272 = t305 * qJDD(1) + t304 * t312;
t308 = sin(qJ(1));
t310 = cos(qJ(1));
t223 = t310 * t271 + t308 * t272;
t303 = g(3) + qJDD(3);
t252 = qJ(3) * t272 + t304 * t303;
t324 = qJ(3) * t271 + t305 * t303;
t375 = -pkin(5) * t223 + t308 * t252 + t310 * t324;
t374 = -2 * qJD(1) * qJD(5);
t367 = pkin(1) + pkin(2);
t309 = cos(qJ(4));
t342 = qJD(1) * qJD(4);
t333 = t309 * t342;
t307 = sin(qJ(4));
t340 = t307 * qJDD(1);
t268 = -t333 - t340;
t300 = qJDD(1) * qJ(2);
t282 = t310 * g(1) + t308 * g(2);
t323 = (2 * qJD(2) * qJD(1)) - t282;
t319 = t300 + t323;
t248 = -t367 * t312 + t319;
t281 = t308 * g(1) - t310 * g(2);
t322 = -qJDD(2) + t281;
t317 = -t312 * qJ(2) - t322;
t313 = -t367 * qJDD(1) + t317;
t198 = t305 * t248 + t304 * t313;
t194 = -t312 * pkin(3) - qJDD(1) * pkin(6) + t198;
t345 = -t307 * t194 + t309 * t303;
t327 = t307 * t374 - t345 - (-t268 - t333) * qJ(5);
t288 = t309 * t312 * t307;
t338 = qJDD(4) + t288;
t370 = t338 * pkin(4);
t162 = -t327 + t370;
t332 = -t308 * t271 + t310 * t272;
t372 = -pkin(5) * t332 + t310 * t252 - t308 * t324;
t197 = t304 * t248 - t305 * t313;
t160 = t305 * t197 - t304 * t198;
t161 = t304 * t197 + t305 * t198;
t137 = t310 * t160 + t308 * t161;
t371 = t308 * t160 - t310 * t161;
t302 = t309 ^ 2;
t298 = t302 * t312;
t311 = qJD(4) ^ 2;
t286 = -t298 - t311;
t347 = t309 * t338;
t236 = t307 * t286 + t347;
t366 = pkin(3) * t236;
t301 = t307 ^ 2;
t356 = t301 * t312;
t284 = -t311 - t356;
t280 = qJDD(4) - t288;
t352 = t307 * t280;
t238 = t309 * t284 - t352;
t365 = pkin(3) * t238;
t353 = t307 * t338;
t240 = t309 * t286 - t353;
t334 = t307 * t342;
t339 = t309 * qJDD(1);
t269 = -0.2e1 * t334 + t339;
t205 = t304 * t240 - t305 * t269;
t207 = t305 * t240 + t304 * t269;
t169 = -t310 * t205 + t308 * t207;
t364 = pkin(5) * t169;
t346 = t309 * t280;
t242 = -t307 * t284 - t346;
t267 = 0.2e1 * t333 + t340;
t206 = t304 * t242 + t305 * t267;
t208 = t305 * t242 - t304 * t267;
t170 = -t310 * t206 + t308 * t208;
t363 = pkin(5) * t170;
t344 = -t301 - t302;
t273 = t344 * qJDD(1);
t276 = t298 + t356;
t227 = t304 * t273 + t305 * t276;
t228 = t305 * t273 - t304 * t276;
t189 = -t310 * t227 + t308 * t228;
t362 = pkin(5) * t189;
t361 = pkin(6) * t236;
t360 = pkin(6) * t238;
t359 = qJ(3) * t227;
t358 = qJ(3) * t228;
t357 = qJDD(1) * pkin(1);
t355 = t307 * t162;
t193 = qJDD(1) * pkin(3) - t312 * pkin(6) + t197;
t354 = t307 * t193;
t349 = t309 * t162;
t348 = t309 * t193;
t187 = t309 * t194 + t307 * t303;
t343 = qJD(1) * t307;
t336 = t304 * t340;
t335 = t305 * t340;
t255 = -t312 * pkin(1) + t319;
t256 = -t317 + t357;
t210 = t310 * t255 - t308 * t256;
t231 = -t308 * t281 - t310 * t282;
t331 = t304 * t288;
t330 = t305 * t288;
t274 = t308 * qJDD(1) + t310 * t312;
t258 = -pkin(5) * t274 + t310 * g(3);
t275 = t310 * qJDD(1) - t308 * t312;
t257 = pkin(5) * t275 + t308 * g(3);
t326 = qJ(2) * t236 - qJ(3) * t205;
t325 = qJ(2) * t238 - qJ(3) * t206;
t150 = -t307 * t187 - t309 * t345;
t151 = t309 * t187 - t307 * t345;
t209 = t308 * t255 + t310 * t256;
t230 = t310 * t281 - t308 * t282;
t321 = -qJ(3) * t207 + t367 * t236;
t320 = -qJ(3) * t208 + t367 * t238;
t270 = t334 - t339;
t278 = qJD(4) * pkin(4) + qJ(5) * t343;
t318 = t270 * qJ(5) - qJD(4) * t278 + t309 * t374 + t187;
t316 = pkin(3) * t269 - pkin(6) * t240 + qJ(2) * t207 - t367 * t205;
t315 = -pkin(3) * t267 - pkin(6) * t242 + qJ(2) * t208 - t367 * t206;
t314 = -pkin(3) * t276 - pkin(6) * t273 + qJ(2) * t228 - t367 * t227;
t175 = -t270 * pkin(4) - qJ(5) * t298 - t278 * t343 + qJDD(5) + t193;
t285 = t298 - t311;
t283 = t311 - t356;
t277 = t298 - t356;
t265 = t344 * t342;
t250 = t309 * t268 + t301 * t342;
t249 = -t307 * t270 + t302 * t342;
t246 = t304 * qJDD(4) + t305 * t265;
t245 = t305 * qJDD(4) - t304 * t265;
t241 = -t307 * t283 + t347;
t239 = t309 * t285 - t352;
t237 = -t309 * t283 - t353;
t235 = -t307 * t285 - t346;
t234 = (-t268 + t333) * t307;
t233 = (-t270 - t334) * t309;
t232 = pkin(4) * t267 - qJ(5) * t280;
t220 = t307 * t267 - t309 * t269;
t219 = t309 * t267 + t307 * t269;
t218 = t305 * t250 - t331;
t217 = t305 * t249 + t331;
t216 = -t304 * t250 - t330;
t215 = -t304 * t249 + t330;
t214 = t305 * t241 - t336;
t213 = t305 * t239 - t304 * t339;
t212 = -t304 * t241 - t335;
t211 = -t304 * t239 - t305 * t339;
t200 = t305 * t220 - t304 * t277;
t199 = -t304 * t220 - t305 * t277;
t196 = -t308 * t245 + t310 * t246;
t195 = t310 * t245 + t308 * t246;
t190 = t308 * t227 + t310 * t228;
t188 = pkin(5) * t190;
t185 = -t308 * t216 + t310 * t218;
t184 = -t308 * t215 + t310 * t217;
t183 = t310 * t216 + t308 * t218;
t182 = t310 * t215 + t308 * t217;
t181 = -t308 * t212 + t310 * t214;
t180 = -t308 * t211 + t310 * t213;
t179 = t310 * t212 + t308 * t214;
t178 = t310 * t211 + t308 * t213;
t177 = t348 - t360;
t176 = t354 - t361;
t174 = t187 - t365;
t173 = -t345 - t366;
t172 = t308 * t206 + t310 * t208;
t171 = t308 * t205 + t310 * t207;
t168 = pkin(5) * t172;
t167 = pkin(5) * t171;
t166 = -qJ(5) * t284 + t175;
t165 = -pkin(4) * t298 + t318;
t164 = -t308 * t199 + t310 * t200;
t163 = t310 * t199 + t308 * t200;
t157 = -qJ(5) * t340 - t162;
t156 = -pkin(4) * t269 + qJ(5) * t286 - t175;
t155 = qJ(2) * t303 + qJ(3) * t160;
t154 = -qJ(3) * t161 + t367 * t303;
t153 = -qJ(5) * t339 + (t276 - t298) * pkin(4) + t318;
t152 = -t365 + (-t284 - t298) * pkin(4) + t318;
t148 = t327 - t366 - 0.2e1 * t370;
t147 = -qJ(5) * t347 - t307 * t156 - t361;
t146 = t309 * t166 - t307 * t232 - t360;
t145 = t305 * t150 - t359;
t144 = -t304 * t150 - t358;
t143 = -pkin(4) * t175 + qJ(5) * t165;
t142 = t305 * t151 + t304 * t193;
t141 = t304 * t151 - t305 * t193;
t140 = t309 * t165 - t355;
t139 = t307 * t165 + t349;
t136 = -t307 * t153 + t309 * t157;
t135 = -t304 * t174 + t305 * t177 + t325;
t134 = -t304 * t173 + t305 * t176 + t326;
t133 = pkin(4) * t336 + t305 * t136 - t359;
t132 = pkin(4) * t335 - t304 * t136 - t358;
t131 = -t305 * t174 - t304 * t177 + t320;
t130 = -t305 * t173 - t304 * t176 + t321;
t129 = t305 * t140 + t304 * t175;
t128 = t304 * t140 - t305 * t175;
t127 = -pkin(3) * t139 - pkin(4) * t162;
t126 = t305 * t146 - t304 * t152 + t325;
t125 = t305 * t147 - t304 * t148 + t326;
t124 = -t304 * t146 - t305 * t152 + t320;
t123 = -t304 * t147 - t305 * t148 + t321;
t122 = t308 * t141 + t310 * t142;
t121 = -t310 * t141 + t308 * t142;
t120 = -pkin(6) * t139 - qJ(5) * t349 - t307 * t143;
t119 = t308 * t128 + t310 * t129;
t118 = -t310 * t128 + t308 * t129;
t117 = -qJ(3) * t141 - (pkin(3) * t304 - pkin(6) * t305 + qJ(2)) * t150;
t116 = -qJ(3) * t142 - (pkin(3) * t305 + pkin(6) * t304 + t367) * t150;
t115 = qJ(2) * t139 - qJ(3) * t128 + t305 * t120 - t304 * t127;
t114 = -qJ(3) * t129 - t304 * t120 - t305 * t127 + t367 * t139;
t1 = [0, 0, 0, 0, 0, 0, 0, 0, 0, -g(1), 0, 0, 0, 0, 0, 0, -t274, -t275, 0, t231, 0, 0, 0, 0, 0, 0, -t274, 0, t275, t210, 0, 0, 0, 0, 0, 0, -t223, t332, 0, -t371, 0, 0, 0, 0, 0, 0, t171, t172, t190, t122, 0, 0, 0, 0, 0, 0, t171, t172, t190, t119; 0, 0, 0, 0, 0, 0, 0, 0, 0, -g(2), 0, 0, 0, 0, 0, 0, t275, -t274, 0, t230, 0, 0, 0, 0, 0, 0, t275, 0, t274, t209, 0, 0, 0, 0, 0, 0, t332, t223, 0, t137, 0, 0, 0, 0, 0, 0, t169, t170, t189, t121, 0, 0, 0, 0, 0, 0, t169, t170, t189, t118; 0, 0, 0, 0, 0, 0, 0, 0, 0, -g(3), 0, 0, 0, 0, 0, 0, 0, 0, 0, -g(3), 0, 0, 0, 0, 0, 0, 0, 0, 0, -g(3), 0, 0, 0, 0, 0, 0, 0, 0, 0, -t303, 0, 0, 0, 0, 0, 0, -t236, -t238, 0, t150, 0, 0, 0, 0, 0, 0, -t236, -t238, 0, -t139; 0, 0, 0, 0, 0, 0, 0, -g(3), g(2), 0, 0, 0, t275, 0, -t274, 0, -t257, -t258, -t230, -pkin(5) * t230, 0, t275, 0, 0, t274, 0, -t257, -t209, t258, -pkin(5) * t209 + (-pkin(1) * t308 + qJ(2) * t310) * g(3), 0, 0, -t332, 0, -t223, 0, t372, t375, t137, -pkin(5) * t137 - t308 * t154 + t310 * t155, t185, t164, t181, t184, t180, t196, -t308 * t130 + t310 * t134 - t364, -t308 * t131 + t310 * t135 - t363, -t308 * t144 + t310 * t145 - t362, -pkin(5) * t121 - t308 * t116 + t310 * t117, t185, t164, t181, t184, t180, t196, -t308 * t123 + t310 * t125 - t364, -t308 * t124 + t310 * t126 - t363, -t308 * t132 + t310 * t133 - t362, -pkin(5) * t118 - t308 * t114 + t310 * t115; 0, 0, 0, 0, 0, 0, g(3), 0, -g(1), 0, 0, 0, t274, 0, t275, 0, t258, -t257, t231, pkin(5) * t231, 0, t274, 0, 0, -t275, 0, t258, t210, t257, pkin(5) * t210 + (pkin(1) * t310 + qJ(2) * t308) * g(3), 0, 0, -t223, 0, t332, 0, t375, -t372, t371, -pkin(5) * t371 + t310 * t154 + t308 * t155, t183, t163, t179, t182, t178, t195, t310 * t130 + t308 * t134 + t167, t310 * t131 + t308 * t135 + t168, t310 * t144 + t308 * t145 + t188, pkin(5) * t122 + t310 * t116 + t308 * t117, t183, t163, t179, t182, t178, t195, t310 * t123 + t308 * t125 + t167, t310 * t124 + t308 * t126 + t168, t310 * t132 + t308 * t133 + t188, pkin(5) * t119 + t310 * t114 + t308 * t115; 0, 0, 0, 0, 0, 0, -g(2), g(1), 0, 0, 0, 0, 0, 0, 0, qJDD(1), t281, t282, 0, 0, 0, 0, 0, qJDD(1), 0, 0, t322 + 0.2e1 * t357, 0, 0.2e1 * t300 + t323, pkin(1) * t256 + qJ(2) * t255, 0, 0, 0, 0, 0, qJDD(1), -qJ(2) * t271 + t367 * t272 + t197, qJ(2) * t272 + t367 * t271 + t198, 0, qJ(2) * t161 + t160 * t367, t234, t219, t237, t233, t235, 0, t316 + t348, t315 - t354, -t151 + t314, pkin(3) * t193 - pkin(6) * t151 + qJ(2) * t142 - t367 * t141, t234, t219, t237, t233, t235, 0, qJ(5) * t353 - t309 * t156 + t316, -t307 * t166 - t309 * t232 + t315, -t309 * t153 - t307 * t157 + t314, pkin(3) * t175 - pkin(6) * t140 + qJ(2) * t129 + qJ(5) * t355 - t367 * t128 - t309 * t143;];
tauB_reg = t1;
