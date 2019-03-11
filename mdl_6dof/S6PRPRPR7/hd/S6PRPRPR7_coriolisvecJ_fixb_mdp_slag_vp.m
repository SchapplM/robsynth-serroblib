% Calculate minimal parameter regressor of Coriolis joint torque vector for
% S6PRPRPR7
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% qJD [6x1]
%   Generalized joint velocities
% pkin [10x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,alpha2,d2,d4,d6,theta1]';
% MDP [25x1]
%   Minimal dynamic parameter vector (fixed base model)
%   see S6PRPRPR7_convert_par2_MPV_fixb.m
% 
% Output:
% tauc [6x1]
%   joint torques required to compensate Coriolis and centrifugal load

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-08 19:54
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauc = S6PRPRPR7_coriolisvecJ_fixb_mdp_slag_vp(qJ, qJD, pkin, MDP)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(10,1),zeros(25,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6PRPRPR7_coriolisvecJ_fixb_mdp_slag_vp: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6PRPRPR7_coriolisvecJ_fixb_mdp_slag_vp: qJD has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S6PRPRPR7_coriolisvecJ_fixb_mdp_slag_vp: pkin has to be [10x1] (double)');
assert(isreal(MDP) && all(size(MDP) == [25 1]), ...
  'S6PRPRPR7_coriolisvecJ_fixb_mdp_slag_vp: MDP has to be [25x1] (double)'); 

%% Symbolic Calculation
% From coriolisvec_joint_fixb_mdp_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-08 19:53:52
% EndTime: 2019-03-08 19:53:57
% DurationCPUTime: 2.42s
% Computational Cost: add. (1042->289), mult. (2318->392), div. (0->0), fcn. (1451->8), ass. (0->139)
t263 = -pkin(2) - pkin(8);
t261 = cos(qJ(2));
t254 = sin(pkin(6));
t329 = qJD(1) * t254;
t302 = t261 * t329;
t282 = qJD(3) - t302;
t215 = qJD(2) * t263 + t282;
t257 = sin(qJ(4));
t260 = cos(qJ(4));
t255 = cos(pkin(6));
t328 = qJD(1) * t255;
t201 = -t215 * t260 + t257 * t328;
t355 = qJD(5) + t201;
t252 = t257 ^ 2;
t253 = t260 ^ 2;
t354 = MDP(9) * (t252 - t253);
t320 = qJD(4) * t260;
t322 = qJD(4) * t257;
t353 = pkin(4) * t320 + qJ(5) * t322;
t195 = -qJD(4) * pkin(4) + t355;
t256 = sin(qJ(6));
t259 = cos(qJ(6));
t352 = MDP(24) * t259 - MDP(25) * t256;
t202 = t215 * t257 + t260 * t328;
t324 = qJD(4) * qJ(5);
t196 = -t324 - t202;
t262 = -pkin(4) - pkin(9);
t319 = qJD(4) * t262;
t325 = qJD(2) * t260;
t306 = pkin(5) * t325;
t311 = t306 + t355;
t258 = sin(qJ(2));
t303 = t258 * t329;
t326 = qJD(2) * t257;
t284 = pkin(4) * t326 + t303;
t290 = -qJ(5) * t260 + qJ(3);
t206 = qJD(2) * t290 + t284;
t245 = qJD(6) + t325;
t312 = t259 * MDP(25);
t273 = MDP(24) * t256 + t312;
t351 = MDP(18) * t206 - t245 * t273;
t301 = qJD(2) * t254 * t258;
t286 = qJD(1) * t301;
t296 = qJD(4) * t328;
t304 = t215 * t320 + (t286 - t296) * t257;
t187 = -qJD(4) * qJD(5) - t304;
t334 = t215 * t322 + t260 * t296;
t188 = -t260 * t286 + t334;
t280 = -t187 * t257 - t188 * t260;
t350 = (t195 * t257 - t196 * t260) * qJD(4) + t280;
t307 = MDP(14) - MDP(17);
t308 = MDP(13) - MDP(16);
t349 = -t257 * t308 - t260 * t307;
t348 = pkin(5) - t263;
t347 = qJD(2) * pkin(2);
t184 = (qJD(5) - t306) * qJD(4) + t304;
t346 = t184 * t256;
t345 = t184 * t259;
t317 = qJD(6) * t256;
t309 = qJD(2) * qJD(4);
t294 = t260 * t309;
t297 = t259 * t326;
t333 = qJD(6) * t297 + t256 * t294;
t203 = -qJD(4) * t317 + t333;
t344 = t203 * t257;
t343 = t203 * t259;
t323 = qJD(4) * t256;
t221 = -t297 + t323;
t342 = t221 * t245;
t341 = t221 * t257;
t321 = qJD(4) * t259;
t223 = t256 * t326 + t321;
t340 = t223 * t245;
t339 = t245 * t260;
t338 = t245 * t262;
t337 = t254 * t261;
t336 = t259 * t260;
t264 = qJD(4) ^ 2;
t335 = t263 * t264;
t295 = t257 * t309;
t332 = pkin(4) * t294 + qJ(5) * t295;
t224 = pkin(4) * t325 + qJ(5) * t326;
t327 = qJD(2) * qJ(3);
t318 = qJD(5) * t260;
t316 = qJD(6) * t259;
t265 = qJD(2) ^ 2;
t305 = t257 * t260 * t265;
t300 = t245 * t317;
t299 = t245 * t316;
t298 = t257 * t316;
t293 = MDP(23) * t326;
t292 = qJD(4) * t348;
t287 = t260 * t303;
t186 = (-pkin(5) * t322 - t287) * qJD(2) + t334;
t270 = qJD(3) + (qJD(4) * pkin(9) - qJD(5)) * t260;
t189 = (t270 + t302) * qJD(2) + t332;
t289 = t186 * t259 - t189 * t256;
t288 = t245 + t325;
t225 = t303 + t327;
t285 = -t225 + t303;
t283 = qJD(4) * t202 - t334;
t281 = qJD(3) + t302;
t190 = t311 + t319;
t271 = pkin(9) * t257 + t290;
t199 = qJD(2) * t271 + t284;
t182 = t190 * t259 - t199 * t256;
t183 = t190 * t256 + t199 * t259;
t251 = t257 * pkin(4);
t216 = t251 + t271;
t229 = t348 * t260;
t278 = t216 * t259 + t229 * t256;
t277 = -qJD(2) * t252 + t339;
t276 = t245 * t256;
t275 = MDP(18) * t196 - MDP(25) * t223;
t194 = -pkin(5) * t326 + t202;
t213 = t255 * t257 + t260 * t337;
t214 = t255 * t260 - t257 * t337;
t269 = t285 - t327;
t227 = t251 + t290;
t268 = -qJD(2) * t227 - t206 + t303;
t192 = (t281 - t318) * qJD(2) + t332;
t211 = qJD(3) - t318 + t353;
t267 = t335 - t192 + (-t211 + t302) * qJD(2);
t219 = t281 * qJD(2);
t266 = qJD(2) * t282 + t219 - t335;
t238 = t259 * t294;
t236 = t256 * t295;
t228 = t348 * t257;
t220 = t282 - t347;
t218 = t260 * t292;
t217 = t257 * t292;
t212 = pkin(9) * t325 + t224;
t205 = t270 + t353;
t204 = qJD(6) * t223 - t238;
t200 = t206 * t325;
t198 = qJD(4) * t214 - t260 * t301;
t197 = qJD(4) * t213 - t257 * t301;
t191 = t194 + t324;
t1 = [(-t187 * t214 + t188 * t213 + t195 * t198 + t196 * t197) * MDP(18) + ((t198 * t259 - t213 * t317) * t245 - t197 * t221 + t214 * t204) * MDP(24) + (-(t198 * t256 + t213 * t316) * t245 - t197 * t223 + t214 * t203) * MDP(25) + (t197 * t257 + t198 * t260) * MDP(15) * qJD(2) + (-t308 * t198 + t307 * t197 + (-t214 * t260 * MDP(15) + (-MDP(15) - t352) * t257 * t213) * qJD(2)) * qJD(4) + (((MDP(7) * t225 + t351) * qJD(2) + (-MDP(4) + MDP(6) - t349) * t265) * t261 + (t192 * MDP(18) + t219 * MDP(7) + (-MDP(3) + MDP(5)) * t265 - t352 * t245 * qJD(6) + ((t220 - t302) * MDP(7) + (t308 * t260 + (t273 - t307) * t257) * qJD(4)) * qJD(2)) * t258) * t254; 0.2e1 * qJD(2) * qJD(3) * MDP(6) + (qJ(3) * t219 + qJD(3) * t225 + (-t225 * t261 + (-t220 - t347) * t258) * t329) * MDP(7) - 0.2e1 * t257 * MDP(8) * t294 + 0.2e1 * t309 * t354 + (t257 * t266 - t269 * t320) * MDP(13) + (t260 * t266 + t269 * t322) * MDP(14) + ((t252 + t253) * t286 - t350) * MDP(15) + (t257 * t267 + t268 * t320) * MDP(16) + (t260 * t267 - t268 * t322) * MDP(17) + (t192 * t227 + t206 * t211 + (-t206 * t261 + (t195 * t260 + t196 * t257) * t258) * t329 + t350 * t263) * MDP(18) + (t256 * t344 + (t256 * t320 + t298) * t223) * MDP(19) + ((-t221 * t256 + t223 * t259) * t320 + (t343 - t204 * t256 + (-t221 * t259 - t223 * t256) * qJD(6)) * t257) * MDP(20) + (t245 * t298 + t203 * t260 + (-t223 * t257 + t256 * t277) * qJD(4)) * MDP(21) + (-t257 * t300 - t204 * t260 + (t259 * t277 + t341) * qJD(4)) * MDP(22) - t288 * MDP(23) * t322 + (-t218 * t221 - t228 * t204 + (-t221 * t303 + t191 * t317 - t345 + (-(-t216 * t256 + t229 * t259) * qJD(2) - t182) * qJD(4)) * t257 + (-qJD(6) * t183 - t191 * t321 + t289) * t260 + (-t205 * t256 - t217 * t259 - t278 * qJD(6) - (-t256 * t261 - t258 * t336) * t329) * t245) * MDP(24) + (-t228 * t203 - t218 * t223 + (-(qJD(6) * t190 + t189) * t260 + (-qJD(6) * t229 - t205 + t302) * t245) * t259 + (-(-qJD(6) * t216 - t217) * t245 + (qJD(4) * t191 + qJD(6) * t199 - t245 * t303 - t186) * t260) * t256 + (-t223 * t303 + t191 * t316 + t346 + (qJD(2) * t278 + t183) * qJD(4)) * t257) * MDP(25) + (-MDP(10) * t257 - MDP(11) * t260) * t264; -t265 * MDP(6) + t280 * MDP(18) + (t204 * t257 + t260 * t300) * MDP(24) + (t260 * t299 + t344) * MDP(25) + (MDP(7) * t285 - t351) * qJD(2) + ((MDP(24) * t221 - t275) * t260 + (t195 * MDP(18) + t288 * t352) * t257) * qJD(4) + t349 * (t264 + t265); MDP(8) * t305 - t265 * t354 + (t285 * t325 + t283) * MDP(13) + (-qJD(4) * t201 + t225 * t326 - t304) * MDP(14) + (t200 + (t224 * t257 - t287) * qJD(2) - t283) * MDP(16) + ((0.2e1 * qJD(5) + t201) * qJD(4) + (-t206 * t257 + t224 * t260) * qJD(2) + t304) * MDP(17) + (-pkin(4) * t188 - qJ(5) * t187 - t195 * t202 - t196 * t355 - t206 * t224) * MDP(18) + (-t223 * t276 + t343) * MDP(19) + ((-t204 - t340) * t259 + (-t203 + t342) * t256) * MDP(20) + (-t300 + (-t256 * t339 + (t223 - t321) * t257) * qJD(2)) * MDP(21) + (-t299 + t236 + (-t245 * t336 - t341) * qJD(2)) * MDP(22) + t245 * t293 + (qJ(5) * t204 + t346 - (t194 * t259 - t212 * t256) * t245 + t311 * t221 + (t191 * t259 - t256 * t338) * qJD(6) + (t191 * t336 + (-t259 * t319 + t182) * t257) * qJD(2)) * MDP(24) + (qJ(5) * t203 + t345 + (t194 * t256 + t212 * t259) * t245 + t311 * t223 + (-t191 * t256 - t259 * t338) * qJD(6) + (-t183 * t257 + (-t191 * t260 + t257 * t319) * t256) * qJD(2)) * MDP(25); -MDP(16) * t305 + (-t253 * t265 - t264) * MDP(17) + (t188 + t200) * MDP(18) + t236 * MDP(25) + ((-t221 - t297) * MDP(24) + t275) * qJD(4) + (-MDP(24) * t276 - t245 * t312) * t245; t223 * t221 * MDP(19) + (-t221 ^ 2 + t223 ^ 2) * MDP(20) + (t333 + t342) * MDP(21) + (t238 + t340) * MDP(22) - qJD(4) * t293 + (t183 * t245 - t191 * t223 + t289) * MDP(24) + (t182 * t245 - t186 * t256 - t189 * t259 + t191 * t221) * MDP(25) + (-MDP(21) * t323 - MDP(22) * t223 - MDP(24) * t183 - MDP(25) * t182) * qJD(6);];
tauc  = t1;
