% Calculate minimal parameter regressor of Coriolis joint torque vector for
% S5PRRRP8
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% qJD [5x1]
%   Generalized joint velocities
% pkin [9x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,alpha2,d2,d3,d4,theta1]';
% MDP [22x1]
%   Minimal dynamic parameter vector (fixed base model)
%   see S5PRRRP8_convert_par2_MPV_fixb.m
% 
% Output:
% tauc [5x1]
%   joint torques required to compensate Coriolis and centrifugal load

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-05 17:01
% Revision: 77da58f92bca3eff71542919beafa37024070d86 (2019-12-05)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauc = S5PRRRP8_coriolisvecJ_fixb_mdp_slag_vp(qJ, qJD, pkin, MDP)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(9,1),zeros(22,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5PRRRP8_coriolisvecJ_fixb_mdp_slag_vp: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5PRRRP8_coriolisvecJ_fixb_mdp_slag_vp: qJD has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S5PRRRP8_coriolisvecJ_fixb_mdp_slag_vp: pkin has to be [9x1] (double)');
assert(isreal(MDP) && all(size(MDP) == [22 1]), ...
  'S5PRRRP8_coriolisvecJ_fixb_mdp_slag_vp: MDP has to be [22x1] (double)'); 

%% Symbolic Calculation
% From coriolisvec_joint_fixb_mdp_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-05 17:00:42
% EndTime: 2019-12-05 17:00:49
% DurationCPUTime: 2.76s
% Computational Cost: add. (1628->313), mult. (4098->439), div. (0->0), fcn. (2794->8), ass. (0->132)
t263 = cos(qJ(3));
t313 = qJD(2) * t263;
t250 = -qJD(4) + t313;
t260 = sin(qJ(3));
t350 = MDP(5) * t260;
t255 = t260 ^ 2;
t349 = (-t263 ^ 2 + t255) * MDP(6);
t261 = sin(qJ(2));
t257 = sin(pkin(5));
t317 = qJD(1) * t257;
t297 = t261 * t317;
t243 = qJD(2) * pkin(7) + t297;
t258 = cos(pkin(5));
t264 = cos(qJ(2));
t315 = qJD(2) * t257;
t293 = t264 * t315;
t311 = qJD(3) * t260;
t316 = qJD(1) * t263;
t200 = -t243 * t311 + (qJD(3) * t258 + t293) * t316;
t330 = t258 * t260;
t249 = qJD(1) * t330;
t223 = t263 * t243 + t249;
t214 = qJD(3) * pkin(8) + t223;
t281 = pkin(3) * t260 - pkin(8) * t263;
t241 = t281 * qJD(3);
t221 = (t241 + t297) * qJD(2);
t246 = -pkin(3) * t263 - pkin(8) * t260 - pkin(2);
t296 = t264 * t317;
t224 = t246 * qJD(2) - t296;
t259 = sin(qJ(4));
t262 = cos(qJ(4));
t308 = qJD(4) * t262;
t309 = qJD(4) * t259;
t285 = -t259 * t200 - t214 * t308 + t262 * t221 - t224 * t309;
t303 = qJD(2) * qJD(3);
t288 = t260 * t303;
t180 = -pkin(4) * t288 - t285;
t188 = t214 * t262 + t224 * t259;
t185 = -qJ(5) * t250 + t188;
t348 = t185 * t250 + t180;
t326 = t263 * t264;
t347 = -(t259 * t326 - t261 * t262) * t317 - t241 * t262 + t246 * t309;
t346 = -(t259 * t261 + t262 * t326) * t317 + t259 * t241 + t246 * t308;
t222 = -t260 * t243 + t258 * t316;
t345 = MDP(17) + MDP(19);
t344 = MDP(18) - MDP(21);
t343 = qJD(2) * pkin(2);
t284 = t260 * t293;
t310 = qJD(3) * t263;
t201 = qJD(1) * t284 + qJD(3) * t249 + t243 * t310;
t287 = t263 * t303;
t290 = t260 * t309;
t302 = qJD(3) * qJD(4);
t217 = qJD(2) * t290 + (-t287 - t302) * t262;
t286 = t259 * t302;
t289 = t260 * t308;
t218 = t286 + (t259 * t310 + t289) * qJD(2);
t312 = qJD(3) * t259;
t314 = qJD(2) * t260;
t238 = t262 * t314 + t312;
t181 = pkin(4) * t218 + qJ(5) * t217 - qJD(5) * t238 + t201;
t342 = t181 * t259;
t341 = t181 * t262;
t339 = t201 * t259;
t338 = t201 * t262;
t337 = t217 * t259;
t305 = t262 * qJD(3);
t236 = t259 * t314 - t305;
t336 = t236 * t250;
t335 = t246 * t262;
t334 = t250 * t259;
t333 = t250 * t262;
t332 = t257 * t261;
t331 = t257 * t264;
t329 = t259 * t263;
t265 = qJD(3) ^ 2;
t328 = t260 * t265;
t327 = t262 * t263;
t325 = t263 * t265;
t307 = qJD(4) * t263;
t324 = qJ(5) * t311 - qJD(5) * t263 + (-t259 * t307 - t260 * t305) * pkin(7) + t346;
t323 = -pkin(4) * t311 + (-t259 * t311 + t262 * t307) * pkin(7) + t347;
t277 = pkin(4) * t259 - qJ(5) * t262;
t322 = qJD(5) * t259 + t250 * t277 + t223;
t240 = t281 * qJD(2);
t321 = t262 * t222 + t259 * t240;
t319 = pkin(7) * t327 + t259 * t246;
t306 = qJD(5) * t250;
t187 = -t214 * t259 + t224 * t262;
t304 = qJD(5) - t187;
t301 = pkin(8) * t334;
t300 = pkin(8) * t333;
t299 = pkin(8) * t311;
t298 = pkin(8) * t305;
t294 = t261 * t315;
t292 = t250 * t309;
t291 = t250 * t308;
t283 = t236 * t296;
t282 = t238 * t296;
t244 = -t296 - t343;
t279 = -t244 - t296;
t278 = pkin(4) * t262 + qJ(5) * t259;
t184 = pkin(4) * t250 + t304;
t276 = t184 * t262 - t185 * t259;
t275 = -t222 * t259 + t240 * t262;
t274 = qJD(2) * t255 - t250 * t263;
t273 = pkin(7) + t277;
t229 = t263 * t332 + t330;
t206 = t229 * t259 + t262 * t331;
t207 = t229 * t262 - t259 * t331;
t228 = -t258 * t263 + t260 * t332;
t270 = -t262 * t200 + t214 * t309 - t259 * t221 - t224 * t308;
t213 = -qJD(3) * pkin(3) - t222;
t268 = qJD(3) * (-t279 - t343);
t267 = -t187 * t250 + t270;
t266 = qJD(2) ^ 2;
t245 = -pkin(3) - t278;
t225 = t273 * t260;
t216 = -t335 + (pkin(7) * t259 + pkin(4)) * t263;
t215 = -qJ(5) * t263 + t319;
t205 = t229 * qJD(3) + t284;
t204 = -t228 * qJD(3) + t263 * t293;
t203 = pkin(4) * t238 + qJ(5) * t236;
t194 = -t217 - t336;
t193 = (t278 * qJD(4) - qJD(5) * t262) * t260 + t273 * t310;
t192 = -pkin(4) * t314 - t275;
t191 = qJ(5) * t314 + t321;
t189 = pkin(4) * t236 - qJ(5) * t238 + t213;
t183 = -t206 * qJD(4) + t204 * t262 + t259 * t294;
t182 = t207 * qJD(4) + t204 * t259 - t262 * t294;
t179 = qJ(5) * t288 - t270 - t306;
t1 = [(t182 * t238 - t183 * t236 - t206 * t217 - t207 * t218) * MDP(20) + (t179 * t207 + t180 * t206 + t181 * t228 + t182 * t184 + t183 * t185 + t189 * t205) * MDP(22) + (-t205 * MDP(10) - t204 * MDP(11) + (-t206 * t345 - t207 * t344) * t314) * qJD(3) + ((-MDP(10) * t260 - MDP(11) * t263) * t264 * t303 + (-t264 * MDP(4) + (-MDP(10) * t263 + MDP(11) * t260 - MDP(3)) * t261) * t266) * t257 + t345 * (t182 * t250 + t205 * t236 + t228 * t218) + t344 * (t183 * t250 + t205 * t238 - t217 * t228); 0.2e1 * t287 * t350 - 0.2e1 * t303 * t349 + MDP(7) * t325 - MDP(8) * t328 + (-pkin(7) * t325 + t260 * t268) * MDP(10) + (pkin(7) * t328 + t263 * t268) * MDP(11) + (-t217 * t260 * t262 + (t263 * t305 - t290) * t238) * MDP(12) + ((-t236 * t262 - t238 * t259) * t310 + (t337 - t218 * t262 + (t236 * t259 - t238 * t262) * qJD(4)) * t260) * MDP(13) + (t250 * t290 + t217 * t263 + (t238 * t260 + t274 * t262) * qJD(3)) * MDP(14) + (t250 * t289 + t218 * t263 + (-t236 * t260 - t274 * t259) * qJD(3)) * MDP(15) + (-t250 - t313) * MDP(16) * t311 + (t347 * t250 + (t213 * t312 + (qJD(3) * t236 + t291) * pkin(7) - t285) * t263 + (-t283 + t213 * t308 + pkin(7) * t218 + t339 + (-pkin(7) * t334 + (-pkin(7) * t329 + t335) * qJD(2) + t187) * qJD(3)) * t260) * MDP(17) + (t346 * t250 + (t213 * t305 + (qJD(3) * t238 - t292) * pkin(7) - t270) * t263 + (-t282 - t213 * t309 - pkin(7) * t217 + t338 + (-pkin(7) * t333 - t319 * qJD(2) - t188) * qJD(3)) * t260) * MDP(18) + (t193 * t236 + t218 * t225 + (t189 * t312 + t180) * t263 + t323 * t250 + (-t283 + t189 * t308 + t342 + (-qJD(2) * t216 - t184) * qJD(3)) * t260) * MDP(19) + (-t215 * t218 - t216 * t217 + t323 * t238 - t324 * t236 + t276 * t310 + (-t179 * t259 + t180 * t262 + (-t184 * t259 - t185 * t262) * qJD(4)) * t260) * MDP(20) + (-t193 * t238 + t217 * t225 + (-t189 * t305 - t179) * t263 - t324 * t250 + (t282 + t189 * t309 - t341 + (qJD(2) * t215 + t185) * qJD(3)) * t260) * MDP(21) + (t179 * t215 + t180 * t216 + t181 * t225 + (-t260 * t296 + t193) * t189 + t324 * t185 + t323 * t184) * MDP(22); (qJD(3) * t223 - t244 * t314 - t201) * MDP(10) + t279 * t313 * MDP(11) + (-t238 * t333 - t337) * MDP(12) + ((-t217 + t336) * t262 + (t238 * t250 - t218) * t259) * MDP(13) + (-t291 + (t250 * t327 + (-t238 + t312) * t260) * qJD(2)) * MDP(14) + (t292 + (-t250 * t329 + (t236 + t305) * t260) * qJD(2)) * MDP(15) + t250 * MDP(16) * t314 + (-pkin(3) * t218 - t338 + t275 * t250 - t223 * t236 + (t213 * t259 + t300) * qJD(4) + (-t187 * t260 + (-t213 * t263 - t299) * t259) * qJD(2)) * MDP(17) + (pkin(3) * t217 + t339 - t321 * t250 - t223 * t238 + (t213 * t262 - t301) * qJD(4) + (-t213 * t327 + (t188 - t298) * t260) * qJD(2)) * MDP(18) + (-t341 - t192 * t250 + t218 * t245 - t322 * t236 + (t189 * t259 + t300) * qJD(4) + (t184 * t260 + (-t189 * t263 - t299) * t259) * qJD(2)) * MDP(19) + (t191 * t236 - t192 * t238 + (t179 - t250 * t184 + (qJD(4) * t238 - t218) * pkin(8)) * t262 + ((qJD(4) * t236 - t217) * pkin(8) + t348) * t259) * MDP(20) + (-t342 + t191 * t250 + t217 * t245 + t322 * t238 + (-t189 * t262 + t301) * qJD(4) + (t189 * t327 + (-t185 + t298) * t260) * qJD(2)) * MDP(21) + (t181 * t245 - t184 * t192 - t185 * t191 - t322 * t189 + (t276 * qJD(4) + t179 * t262 + t180 * t259) * pkin(8)) * MDP(22) + (-t263 * t350 + t349) * t266; t194 * MDP(14) - MDP(15) * t286 + t267 * MDP(18) + (pkin(4) * t217 - qJ(5) * t218) * MDP(20) + (-t267 - 0.2e1 * t306) * MDP(21) + (-pkin(4) * t180 + qJ(5) * t179 - t184 * t188 + t304 * t185 - t189 * t203) * MDP(22) + (-t250 * MDP(15) - t213 * MDP(17) - t189 * MDP(19) + (t185 - t188) * MDP(20) + t203 * MDP(21) + MDP(13) * t238) * t238 + (-MDP(15) * t289 + (-MDP(15) * t329 + (0.2e1 * pkin(4) * MDP(19) + 0.2e1 * qJ(5) * MDP(21) + MDP(16)) * t260) * qJD(3)) * qJD(2) + (t238 * MDP(12) + t213 * MDP(18) - t203 * MDP(19) + (t184 - t304) * MDP(20) - t189 * MDP(21) - MDP(13) * t236) * t236 + t345 * (-t188 * t250 + t285); (t236 * t238 - t288) * MDP(19) + t194 * MDP(20) + (-t238 ^ 2 - t250 ^ 2) * MDP(21) + (t189 * t238 + t348) * MDP(22);];
tauc = t1;
