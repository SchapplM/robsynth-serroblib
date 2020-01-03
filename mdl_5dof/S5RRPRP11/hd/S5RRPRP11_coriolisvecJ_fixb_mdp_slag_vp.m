% Calculate minimal parameter regressor of Coriolis joint torque vector for
% S5RRPRP11
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% qJD [5x1]
%   Generalized joint velocities
% pkin [7x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d1,d2,d4]';
% MDP [25x1]
%   Minimal dynamic parameter vector (fixed base model)
%   see S5RRPRP11_convert_par2_MPV_fixb.m
% 
% Output:
% tauc [5x1]
%   joint torques required to compensate Coriolis and centrifugal load

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-31 20:14
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauc = S5RRPRP11_coriolisvecJ_fixb_mdp_slag_vp(qJ, qJD, pkin, MDP)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(7,1),zeros(25,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RRPRP11_coriolisvecJ_fixb_mdp_slag_vp: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RRPRP11_coriolisvecJ_fixb_mdp_slag_vp: qJD has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [7 1]), ...
  'S5RRPRP11_coriolisvecJ_fixb_mdp_slag_vp: pkin has to be [7x1] (double)');
assert(isreal(MDP) && all(size(MDP) == [25 1]), ...
  'S5RRPRP11_coriolisvecJ_fixb_mdp_slag_vp: MDP has to be [25x1] (double)'); 

%% Symbolic Calculation
% From coriolisvec_joint_fixb_mdp_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 20:14:09
% EndTime: 2019-12-31 20:14:15
% DurationCPUTime: 2.49s
% Computational Cost: add. (1649->324), mult. (3660->431), div. (0->0), fcn. (1939->4), ass. (0->142)
t341 = pkin(3) + pkin(6);
t263 = sin(qJ(2));
t322 = qJD(1) * t263;
t250 = qJD(4) + t322;
t347 = t250 ^ 2;
t309 = qJD(1) * qJD(2);
t346 = -0.2e1 * t309;
t260 = t263 ^ 2;
t265 = cos(qJ(2));
t261 = t265 ^ 2;
t345 = (t260 - t261) * MDP(5);
t266 = -pkin(2) - pkin(7);
t295 = -qJ(3) * t263 - pkin(1);
t223 = t266 * t265 + t295;
t241 = t341 * t263;
t262 = sin(qJ(4));
t264 = cos(qJ(4));
t344 = t264 * t223 + t262 * t241;
t298 = t263 * t309;
t320 = qJD(2) * t262;
t321 = qJD(1) * t265;
t225 = t264 * t321 + t320;
t315 = qJD(4) * t225;
t200 = -t262 * t298 + t315;
t251 = pkin(6) * t322;
t343 = qJD(3) + t251;
t308 = MDP(20) + MDP(22);
t307 = -MDP(21) + MDP(24);
t319 = qJD(2) * t263;
t254 = pkin(2) * t319;
t286 = pkin(7) * t263 - qJ(3) * t265;
t316 = qJD(3) * t263;
t270 = qJD(2) * t286 - t316;
t205 = t254 + t270;
t317 = qJD(2) * t265;
t234 = t341 * t317;
t342 = -qJD(4) * t344 - t205 * t262 + t234 * t264;
t340 = qJD(2) * pkin(2);
t299 = t262 * t321;
t313 = qJD(4) * t264;
t201 = qJD(2) * t313 - qJD(4) * t299 - t264 * t298;
t232 = t341 * t319;
t258 = qJD(2) * qJD(3);
t212 = -qJD(1) * t232 + t258;
t318 = qJD(2) * t264;
t227 = -t299 + t318;
t180 = pkin(4) * t201 + qJ(5) * t200 - qJD(5) * t227 + t212;
t339 = t180 * t262;
t338 = t180 * t264;
t337 = t200 * t264;
t336 = t212 * t262;
t335 = t212 * t264;
t334 = t225 * t250;
t333 = t227 * t265;
t332 = t250 * t263;
t331 = t250 * t266;
t330 = t263 * t264;
t267 = qJD(2) ^ 2;
t329 = t263 * t267;
t328 = t265 * t267;
t268 = qJD(1) ^ 2;
t327 = t265 * t268;
t288 = pkin(4) * t264 + qJ(5) * t262;
t277 = -pkin(3) - t288;
t326 = -qJD(4) * t288 + qJD(5) * t264 + t277 * t322 - t343;
t255 = pkin(2) * t322;
t215 = qJD(1) * t286 + t255;
t252 = pkin(6) * t321;
t233 = pkin(3) * t321 + t252;
t325 = t264 * t215 + t262 * t233;
t242 = t341 * t265;
t239 = -pkin(2) * t265 + t295;
t219 = qJD(1) * t239;
t259 = qJD(2) * qJ(3);
t314 = qJD(4) * t262;
t312 = qJD(5) * t250;
t311 = pkin(3) * t322 + t343;
t208 = t223 * qJD(1);
t210 = t266 * qJD(2) + t311;
t185 = -t208 * t262 + t210 * t264;
t310 = qJD(5) - t185;
t306 = t262 * t331;
t305 = t264 * t331;
t304 = t263 * t327;
t249 = pkin(2) * t298;
t197 = qJD(1) * t270 + t249;
t297 = t265 * t309;
t248 = pkin(6) * t297;
t222 = pkin(3) * t297 + t248;
t303 = -t264 * t197 - t210 * t313 - t262 * t222;
t218 = t259 + t233;
t302 = t266 * t317;
t301 = t250 * t314;
t300 = t265 * t313;
t296 = MDP(19) * t317;
t294 = pkin(1) * t346;
t293 = qJD(3) - t340;
t292 = pkin(4) * t297;
t291 = -t262 * t197 - t208 * t313 - t210 * t314 + t264 * t222;
t290 = qJ(5) * t297;
t244 = t264 * t297;
t178 = -t291 - t292;
t186 = t208 * t264 + t210 * t262;
t183 = qJ(5) * t250 + t186;
t289 = -t183 * t322 + t178;
t287 = -pkin(4) * t262 + qJ(5) * t264;
t182 = -pkin(4) * t250 + t310;
t285 = t182 * t262 + t183 * t264;
t283 = -t215 * t262 + t233 * t264;
t281 = -t223 * t262 + t241 * t264;
t280 = -qJD(1) * t261 + t332;
t279 = -0.2e1 * qJD(2) * t219;
t278 = t227 * t250;
t275 = -qJ(3) * t317 - t316;
t206 = qJD(1) * t275 + t249;
t217 = t254 + t275;
t276 = pkin(6) * t267 + qJD(1) * t217 + t206;
t274 = t186 * t250 + t291;
t273 = t208 * t314 + t303;
t272 = t264 * t205 - t223 * t314 + t262 * t234 + t241 * t313;
t271 = t185 * t250 + t273;
t235 = pkin(6) * t298 - t258;
t237 = t251 + t293;
t240 = -t252 - t259;
t269 = -t235 * t265 + (t237 * t265 + (t240 + t252) * t263) * qJD(2);
t238 = qJ(3) - t287;
t236 = t266 * t244;
t230 = -qJ(3) * t321 + t255;
t209 = t219 * t322;
t204 = t265 * t288 + t242;
t195 = pkin(4) * t227 + qJ(5) * t225;
t192 = -pkin(4) * t263 - t281;
t191 = qJ(5) * t263 + t344;
t190 = -pkin(4) * t321 - t283;
t189 = qJ(5) * t321 + t325;
t188 = pkin(4) * t225 - qJ(5) * t227 + t218;
t187 = t334 - t200;
t184 = (qJD(4) * t287 + qJD(5) * t262) * t265 + (-pkin(6) + t277) * t319;
t181 = -pkin(4) * t317 - t342;
t179 = qJ(5) * t317 + qJD(5) * t263 + t272;
t177 = -t273 + t290 + t312;
t1 = [0.2e1 * t263 * MDP(4) * t297 + t345 * t346 + MDP(6) * t328 - MDP(7) * t329 + (-pkin(6) * t328 + t263 * t294) * MDP(9) + (pkin(6) * t329 + t265 * t294) * MDP(10) + t269 * MDP(11) + (t263 * t279 + t265 * t276) * MDP(12) + (-t263 * t276 + t265 * t279) * MDP(13) + (pkin(6) * t269 + t206 * t239 + t217 * t219) * MDP(14) + (t200 * t262 * t265 + (t262 * t319 - t300) * t227) * MDP(15) + ((-t225 * t262 + t227 * t264) * t319 + (t337 + t201 * t262 + (t225 * t264 + t227 * t262) * qJD(4)) * t265) * MDP(16) + (-t250 * t300 - t200 * t263 + (t262 * t280 + t333) * qJD(2)) * MDP(17) + (t265 * t301 - t201 * t263 + (-t225 * t265 + t264 * t280) * qJD(2)) * MDP(18) + (t250 + t322) * t296 + (t342 * t250 - t232 * t225 + t242 * t201 + (-t218 * t318 + t291) * t263 + (-t218 * t314 + t335 + (qJD(1) * t281 + t185) * qJD(2)) * t265) * MDP(20) + (-t272 * t250 - t232 * t227 - t242 * t200 + ((qJD(2) * t218 + qJD(4) * t208) * t262 + t303) * t263 + (-t218 * t313 - t336 + (-qJD(1) * t344 - t186) * qJD(2)) * t265) * MDP(21) + (-t181 * t250 + t184 * t225 + t201 * t204 + (-t188 * t318 - t178) * t263 + (-t188 * t314 + t338 + (-qJD(1) * t192 - t182) * qJD(2)) * t265) * MDP(22) + (-t179 * t225 + t181 * t227 - t191 * t201 - t192 * t200 + t285 * t319 + (-t177 * t264 - t178 * t262 + (-t182 * t264 + t183 * t262) * qJD(4)) * t265) * MDP(23) + (t179 * t250 - t184 * t227 + t200 * t204 + (-t188 * t320 + t177) * t263 + (t188 * t313 + t339 + (qJD(1) * t191 + t183) * qJD(2)) * t265) * MDP(24) + (t177 * t191 + t178 * t192 + t179 * t183 + t180 * t204 + t181 * t182 + t184 * t188) * MDP(25); -MDP(4) * t304 + t268 * t345 + ((-t240 - t259) * t263 + (-t237 + t293) * t265) * qJD(1) * MDP(11) + (-t230 * t321 + t209) * MDP(12) + (0.2e1 * t258 + (t219 * t265 + t230 * t263) * qJD(1)) * MDP(13) + (-qJ(3) * t235 - qJD(3) * t240 - t219 * t230 + (-t240 * t263 + (-t237 - t340) * t265) * qJD(1) * pkin(6)) * MDP(14) + (-t262 * t278 - t337) * MDP(15) + ((-t201 - t278) * t264 + (t200 + t334) * t262) * MDP(16) + (-t301 + t244 + (-t262 * t332 - t333) * qJD(1)) * MDP(17) + (-t250 * t313 + (-t250 * t330 + (t225 - t320) * t265) * qJD(1)) * MDP(18) - t250 * MDP(19) * t321 + (t236 + qJ(3) * t201 + t336 - t283 * t250 + t311 * t225 + (t218 * t264 - t306) * qJD(4) + (-t185 * t265 + t218 * t330) * qJD(1)) * MDP(20) + (-qJ(3) * t200 + t335 + t325 * t250 + t311 * t227 + (-t218 * t262 - t305) * qJD(4) + (t186 * t265 + (-t218 * t263 - t302) * t262) * qJD(1)) * MDP(21) + (t339 + t190 * t250 + t201 * t238 + t236 - t326 * t225 + (t188 * t264 - t306) * qJD(4) + (t182 * t265 + t188 * t330) * qJD(1)) * MDP(22) + (t189 * t225 - t190 * t227 + (t200 * t266 + (-t225 * t266 - t183) * qJD(4) + t289) * t264 + (-t182 * t322 - t201 * t266 - t177 + (t227 * t266 - t182) * qJD(4)) * t262) * MDP(23) + (-t338 - t189 * t250 + t200 * t238 + t326 * t227 + (t188 * t262 + t305) * qJD(4) + (-t183 * t265 + (t188 * t263 + t302) * t262) * qJD(1)) * MDP(24) + (t180 * t238 - t182 * t190 - t183 * t189 - t326 * t188 + (qJD(4) * t285 + t177 * t262 - t178 * t264) * t266) * MDP(25) + (t268 * t263 * MDP(9) + MDP(10) * t327) * pkin(1); MDP(12) * t304 + (-t260 * t268 - t267) * MDP(13) + (t209 + t248) * MDP(14) + t308 * t244 + (t240 * MDP(14) - t188 * MDP(25) - t225 * t308 + t227 * t307) * qJD(2) + ((-t225 * t322 + t200 - t315) * MDP(23) + (qJD(4) * t183 - t289) * MDP(25) + t307 * t347) * t264 + ((qJD(4) * t227 - t201) * MDP(23) + (qJD(4) * t182 + t177) * MDP(25) + ((t227 * MDP(23) + MDP(25) * t182) * t263 + t307 * t317) * qJD(1) - t308 * t347) * t262; t187 * MDP(17) - t201 * MDP(18) + qJD(1) * t296 + t274 * MDP(20) + t271 * MDP(21) + (t274 + 0.2e1 * t292) * MDP(22) + (pkin(4) * t200 - qJ(5) * t201) * MDP(23) + (-t271 + 0.2e1 * t290 + 0.2e1 * t312) * MDP(24) + (-pkin(4) * t178 + qJ(5) * t177 - t182 * t186 + t183 * t310 - t188 * t195) * MDP(25) + (t250 * MDP(18) - t218 * MDP(20) - t188 * MDP(22) + (t183 - t186) * MDP(23) + t195 * MDP(24) + MDP(16) * t227) * t227 + (t227 * MDP(15) + t218 * MDP(21) - t195 * MDP(22) + (t182 - t310) * MDP(23) - t188 * MDP(24) - MDP(16) * t225) * t225; (t225 * t227 - t297) * MDP(22) + t187 * MDP(23) + (-t227 ^ 2 - t347) * MDP(24) + (-t183 * t250 + t188 * t227 + t178) * MDP(25);];
tauc = t1;
