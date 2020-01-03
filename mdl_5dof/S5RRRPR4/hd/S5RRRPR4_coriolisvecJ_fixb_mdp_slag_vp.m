% Calculate minimal parameter regressor of Coriolis joint torque vector for
% S5RRRPR4
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% qJD [5x1]
%   Generalized joint velocities
% pkin [8x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d1,d2,d3,d5]';
% MDP [24x1]
%   Minimal dynamic parameter vector (fixed base model)
%   see S5RRRPR4_convert_par2_MPV_fixb.m
% 
% Output:
% tauc [5x1]
%   joint torques required to compensate Coriolis and centrifugal load

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-31 21:12
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauc = S5RRRPR4_coriolisvecJ_fixb_mdp_slag_vp(qJ, qJD, pkin, MDP)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(8,1),zeros(24,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RRRPR4_coriolisvecJ_fixb_mdp_slag_vp: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RRRPR4_coriolisvecJ_fixb_mdp_slag_vp: qJD has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S5RRRPR4_coriolisvecJ_fixb_mdp_slag_vp: pkin has to be [8x1] (double)');
assert(isreal(MDP) && all(size(MDP) == [24 1]), ...
  'S5RRRPR4_coriolisvecJ_fixb_mdp_slag_vp: MDP has to be [24x1] (double)'); 

%% Symbolic Calculation
% From coriolisvec_joint_fixb_mdp_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 21:11:44
% EndTime: 2019-12-31 21:11:49
% DurationCPUTime: 2.19s
% Computational Cost: add. (1123->238), mult. (1907->324), div. (0->0), fcn. (1039->6), ass. (0->124)
t266 = sin(qJ(5));
t267 = sin(qJ(3));
t269 = cos(qJ(5));
t270 = cos(qJ(3));
t229 = -t266 * t270 + t267 * t269;
t310 = qJD(3) * t270;
t311 = qJD(3) * t267;
t193 = t229 * qJD(5) + t266 * t310 - t269 * t311;
t261 = qJD(1) + qJD(2);
t228 = t266 * t267 + t269 * t270;
t278 = t228 * qJD(5);
t295 = t261 * t310;
t296 = t261 * t311;
t181 = -t261 * t278 + t266 * t296 + t269 * t295;
t182 = t193 * t261;
t211 = t228 * t261;
t213 = t229 * t261;
t260 = qJD(3) - qJD(5);
t349 = (-t211 * t260 + t181) * MDP(20) - (t213 * t260 + t182) * MDP(21) - (t211 ^ 2 - t213 ^ 2) * MDP(19) + t213 * t211 * MDP(18);
t345 = MDP(7) * t267;
t264 = t267 ^ 2;
t265 = t270 ^ 2;
t344 = (t264 - t265) * MDP(8);
t268 = sin(qJ(2));
t332 = pkin(1) * qJD(1);
t302 = t268 * t332;
t233 = pkin(7) * t261 + t302;
t271 = cos(qJ(2));
t331 = pkin(1) * qJD(2);
t298 = qJD(1) * t331;
t285 = t271 * t298;
t242 = t270 * t285;
t262 = qJD(3) * qJD(4);
t315 = t242 + t262;
t195 = -t233 * t311 + t315;
t241 = t267 * t285;
t200 = t233 * t310 + t241;
t341 = t195 * t270 + t200 * t267;
t256 = t267 * qJ(4);
t258 = t270 * pkin(3);
t340 = t256 + t258;
t221 = t267 * t233;
t339 = qJD(4) + t221;
t338 = qJD(5) + t260;
t336 = pkin(3) + pkin(4);
t335 = pkin(7) - pkin(8);
t273 = qJD(3) ^ 2;
t334 = pkin(7) * t273;
t251 = pkin(1) * t268 + pkin(7);
t333 = -pkin(8) + t251;
t330 = qJ(4) * t270;
t259 = t261 ^ 2;
t327 = t259 * t270;
t326 = t260 * t271;
t325 = t261 * t267;
t324 = t261 * t270;
t222 = t270 * t233;
t255 = t267 * qJD(4);
t279 = -t336 * t267 + t330;
t286 = t268 * t298;
t183 = -t286 + (t279 * qJD(3) + t255) * t261;
t290 = pkin(2) + t256;
t301 = t271 * t332;
t188 = t301 + (t336 * t270 + t290) * t261;
t321 = t183 * t228 + t188 * t193;
t194 = t228 * qJD(3) - t278;
t320 = t183 * t229 + t188 * t194;
t234 = -pkin(2) * t261 - t301;
t319 = t234 * t310 + t267 * t286;
t318 = t301 * t311 + t302 * t324;
t313 = t264 + t265;
t312 = MDP(17) * t251;
t263 = qJD(3) * qJ(4);
t309 = t261 * MDP(13);
t308 = t273 * MDP(10);
t306 = -pkin(8) * t325 + t339;
t305 = -MDP(12) - MDP(14);
t304 = MDP(13) - MDP(16);
t303 = pkin(2) + t340;
t300 = t268 * t331;
t299 = t271 * t331;
t246 = t335 * t270;
t288 = -qJD(3) * pkin(3) + qJD(4);
t210 = t221 + t288;
t297 = t210 * t310 + t341;
t216 = pkin(3) * t311 - qJ(4) * t310 - t255;
t252 = -pkin(1) * t271 - pkin(2);
t226 = t333 * t270;
t283 = pkin(3) * t267 - t330;
t189 = t286 + (t283 * qJD(3) - t255) * t261;
t289 = -t189 - t334;
t208 = t216 + t300;
t287 = -t208 * t261 - t189;
t206 = -pkin(8) * t324 + t222;
t223 = t252 - t340;
t192 = -t336 * qJD(3) + t306;
t199 = t206 + t263;
t282 = t269 * t192 - t266 * t199;
t281 = -t266 * t192 - t269 * t199;
t215 = t222 + t263;
t280 = t210 * t267 + t215 * t270;
t207 = -pkin(4) * t311 - t216;
t277 = (-t181 * t228 - t182 * t229 - t193 * t213 - t194 * t211) * MDP(19) + (t181 * t229 + t194 * t213) * MDP(18) - 0.2e1 * t261 * qJD(3) * t344 + 0.2e1 * t295 * t345 + t273 * t270 * MDP(9) + (-t194 * MDP(20) + t193 * MDP(21)) * t260;
t184 = (pkin(8) * t261 - t233) * t311 + t315;
t187 = -pkin(8) * t295 + t200;
t275 = -t269 * t184 - t266 * t187 + t188 * t211;
t274 = -t266 * t184 + t269 * t187 - t188 * t213;
t257 = t270 * pkin(4);
t253 = pkin(8) * t311;
t245 = t335 * t267;
t231 = qJD(3) * t246;
t230 = -pkin(7) * t311 + t253;
t225 = t333 * t267;
t224 = t257 + t303;
t219 = t234 * t311;
t217 = t283 * t261;
t214 = -t223 + t257;
t204 = t279 * t261;
t202 = qJD(3) * t226 + t267 * t299;
t201 = -t251 * t311 + t270 * t299 + t253;
t198 = -t301 + (-t290 - t258) * t261;
t197 = t207 - t300;
t191 = t198 * t311;
t1 = [(t189 * t223 + t198 * t208) * MDP(17) + (-t308 + t287 * MDP(16) + (t200 * MDP(17) + t304 * t273) * t251 + ((MDP(12) * t252 + MDP(14) * t223) * t261 + (-MDP(15) - t312) * t215) * qJD(3)) * t267 + (t197 * t211 + t214 * t182 - (-t201 * t266 + t202 * t269 + (-t225 * t266 - t226 * t269) * qJD(5)) * t260 + t321) * MDP(23) + (t287 * MDP(14) + (t195 * MDP(17) + t305 * t273) * t251 + (t252 * t309 + (-t223 * t261 - t198) * MDP(16) + t210 * t312) * qJD(3)) * t270 + (t197 * t213 + t214 * t181 + (t201 * t269 + t202 * t266 + (t225 * t269 - t226 * t266) * qJD(5)) * t260 + t320) * MDP(24) + ((t267 * t309 + (MDP(12) * t270 + MDP(5)) * (-qJD(1) - t261)) * t268 + (-qJD(1) * MDP(6) + t280 * MDP(17) + (t313 * MDP(15) - MDP(6)) * t261 + (t305 * t267 - t304 * t270) * qJD(3)) * t271) * t331 + t277 + t191 * MDP(14) + t219 * MDP(12) + t319 * MDP(13) + t297 * MDP(15); (-t303 * t296 + t191 + (-t216 * t261 + t289) * t270 + t318) * MDP(14) + (-t189 * t303 + t198 * t216 + (-t198 * t268 - t280 * t271) * t332 + ((t210 * t270 - t215 * t267) * qJD(3) + t341) * pkin(7)) * MDP(17) + (t207 * t211 + t224 * t182 - (-t230 * t266 + t231 * t269 + (-t245 * t266 - t246 * t269) * qJD(5)) * t260 + (t268 * t211 + t229 * t326) * t332 + t321) * MDP(23) + (t207 * t213 + t224 * t181 + (t230 * t269 + t231 * t266 + (t245 * t269 - t246 * t266) * qJD(5)) * t260 + (t268 * t213 - t228 * t326) * t332 + t320) * MDP(24) + (-t313 * t261 * t301 - t215 * t311 + t297) * MDP(15) + t277 + ((t261 * t303 - t198 - t301) * t310 + ((-t216 + t302) * t261 + t289) * t267) * MDP(16) - t267 * t308 + (-pkin(2) * t295 + t267 * t334 + (-t268 * t325 + t271 * t310) * t332 + t319) * MDP(13) + (-pkin(2) * t296 + t219 + (-t286 - t334) * t270 + t318) * MDP(12) + (t268 * MDP(5) + t271 * MDP(6)) * (-qJD(2) + t261) * t332; -t327 * t345 + t259 * t344 + (-t234 * t325 - t241) * MDP(12) + (-t234 * t324 - t242) * MDP(13) + (-t241 + (-t198 * t267 + t217 * t270) * t261) * MDP(14) + ((t215 - t263) * t267 + (-t210 + t288) * t270) * t261 * MDP(15) + (t242 + 0.2e1 * t262 + (t198 * t270 + t217 * t267) * t261) * MDP(16) + (-pkin(3) * t200 + qJ(4) * t195 - t198 * t217 - t210 * t222 + t215 * t339) * MDP(17) + (-t204 * t211 + (t269 * t206 + t306 * t266) * t260 + (-(-qJ(4) * t269 + t266 * t336) * t260 - t281) * qJD(5) - t274) * MDP(23) + (-t204 * t213 + (-t266 * t206 + t306 * t269) * t260 + ((-qJ(4) * t266 - t269 * t336) * t260 + t282) * qJD(5) - t275) * MDP(24) - t349; (-t259 * t264 - t273) * MDP(16) + (-qJD(3) * t215 + t200) * MDP(17) + (-MDP(14) * t327 + (t198 * MDP(17) - t211 * MDP(23) - t213 * MDP(24)) * t261) * t267 + (-MDP(23) * t266 - MDP(24) * t269) * t260 ^ 2; (t281 * t338 + t274) * MDP(23) + (-t282 * t338 + t275) * MDP(24) + t349;];
tauc = t1;
