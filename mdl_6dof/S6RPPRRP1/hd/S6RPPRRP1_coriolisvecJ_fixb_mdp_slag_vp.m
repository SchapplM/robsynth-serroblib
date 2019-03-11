% Calculate minimal parameter regressor of Coriolis joint torque vector for
% S6RPPRRP1
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% qJD [6x1]
%   Generalized joint velocities
% pkin [10x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d4,d5,theta2,theta3]';
% MDP [24x1]
%   Minimal dynamic parameter vector (fixed base model)
%   see S6RPPRRP1_convert_par2_MPV_fixb.m
% 
% Output:
% tauc [6x1]
%   joint torques required to compensate Coriolis and centrifugal load

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-09 01:59
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauc = S6RPPRRP1_coriolisvecJ_fixb_mdp_slag_vp(qJ, qJD, pkin, MDP)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(10,1),zeros(24,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPPRRP1_coriolisvecJ_fixb_mdp_slag_vp: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RPPRRP1_coriolisvecJ_fixb_mdp_slag_vp: qJD has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S6RPPRRP1_coriolisvecJ_fixb_mdp_slag_vp: pkin has to be [10x1] (double)');
assert(isreal(MDP) && all(size(MDP) == [24 1]), ...
  'S6RPPRRP1_coriolisvecJ_fixb_mdp_slag_vp: MDP has to be [24x1] (double)'); 

%% Symbolic Calculation
% From coriolisvec_joint_fixb_mdp_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 01:58:47
% EndTime: 2019-03-09 01:58:52
% DurationCPUTime: 2.25s
% Computational Cost: add. (2725->283), mult. (6503->380), div. (0->0), fcn. (4672->8), ass. (0->121)
t270 = sin(pkin(9)) * pkin(1) + qJ(3);
t263 = t270 * qJD(1);
t280 = cos(pkin(10));
t274 = t280 * qJD(2);
t278 = sin(pkin(10));
t339 = pkin(7) * qJD(1);
t237 = t274 + (-t263 - t339) * t278;
t246 = qJD(2) * t278 + t263 * t280;
t238 = t280 * t339 + t246;
t283 = sin(qJ(4));
t285 = cos(qJ(4));
t202 = t237 * t283 + t238 * t285;
t352 = qJD(4) * t202;
t301 = qJD(1) * (t278 ^ 2 + t280 ^ 2);
t351 = MDP(7) * t301;
t261 = t278 * t285 + t280 * t283;
t254 = t261 * qJD(1);
t284 = cos(qJ(5));
t313 = qJD(1) * t285;
t269 = t280 * t313;
t268 = qJD(4) * t269;
t314 = qJD(1) * t283;
t306 = t278 * t314;
t294 = qJD(4) * t306 - t268;
t310 = t284 * qJD(4);
t282 = sin(qJ(5));
t312 = qJD(5) * t282;
t209 = -qJD(5) * t310 + t254 * t312 + t284 * t294;
t241 = qJD(4) * t282 + t254 * t284;
t256 = t261 * qJD(4);
t247 = qJD(1) * t256;
t199 = qJD(4) * pkin(8) + t202;
t262 = -cos(pkin(9)) * pkin(1) - pkin(3) * t280 - pkin(2);
t250 = qJD(1) * t262 + qJD(3);
t253 = t269 - t306;
t208 = -pkin(4) * t253 - pkin(8) * t254 + t250;
t188 = t199 * t284 + t208 * t282;
t260 = t278 * t283 - t280 * t285;
t289 = t260 * qJD(3);
t343 = qJD(1) * t289;
t344 = t237 * t285 - t238 * t283;
t193 = qJD(4) * t344 - t343;
t217 = pkin(4) * t247 + pkin(8) * t294;
t213 = t284 * t217;
t287 = -qJD(5) * t188 - t193 * t282 + t213;
t177 = pkin(5) * t247 + qJ(6) * t209 - qJD(6) * t241 + t287;
t239 = t254 * t282 - t310;
t183 = -qJ(6) * t239 + t188;
t249 = qJD(5) - t253;
t350 = t183 * t249 + t177;
t210 = qJD(5) * t241 - t282 * t294;
t311 = qJD(5) * t284;
t290 = t193 * t284 - t199 * t312 + t208 * t311 + t217 * t282;
t178 = -qJ(6) * t210 - qJD(6) * t239 + t290;
t187 = -t199 * t282 + t208 * t284;
t182 = -qJ(6) * t241 + t187;
t181 = pkin(5) * t249 + t182;
t349 = -t181 * t249 + t178;
t340 = -qJ(6) - pkin(8);
t348 = qJ(6) * t253 + qJD(5) * t340;
t347 = t278 * t313 + t280 * t314;
t341 = pkin(7) + t270;
t257 = t341 * t278;
t258 = t341 * t280;
t225 = t257 * t285 + t258 * t283;
t342 = t241 ^ 2;
t337 = t209 * t282;
t335 = t239 * t253;
t334 = t239 * t282;
t333 = t241 * t249;
t332 = t241 * t284;
t255 = t260 * qJD(4);
t331 = t255 * t282;
t330 = t255 * t284;
t328 = t261 * t282;
t327 = t261 * t284;
t326 = t282 * t247;
t325 = t282 * t249;
t226 = -t257 * t283 + t258 * t285;
t219 = t284 * t226;
t244 = t284 * t247;
t324 = t181 - t182;
t323 = -t210 * t327 + t239 * t330;
t228 = pkin(4) * t254 - pkin(8) * t253;
t322 = t228 * t282 + t284 * t344;
t321 = -t210 * t282 - t239 * t311;
t320 = -t209 * t260 + t241 * t256;
t220 = pkin(4) * t260 - pkin(8) * t261 + t262;
t319 = t220 * t282 + t219;
t318 = t253 * t325 + t244;
t317 = qJD(6) * t284 + t282 * t348 - t322;
t223 = t284 * t228;
t316 = -pkin(5) * t254 - t223 + t348 * t284 + (-qJD(6) + t344) * t282;
t308 = t241 * t331;
t204 = -qJD(4) * t225 - t289;
t229 = pkin(4) * t256 + pkin(8) * t255;
t307 = t204 * t284 + t220 * t311 + t229 * t282;
t303 = t261 * t311;
t300 = t249 * t284;
t194 = qJD(3) * t347 + t352;
t299 = -t181 * t284 - t183 * t282;
t298 = -t210 * t260 - t239 * t256;
t296 = (-t263 * t278 + t274) * t278 - t246 * t280;
t295 = qJ(6) * t255 - qJD(6) * t261;
t198 = -qJD(4) * pkin(4) - t344;
t293 = t303 - t331;
t292 = -t261 * t312 - t330;
t185 = pkin(5) * t210 + t194;
t291 = -pkin(8) * t247 + t198 * t249;
t205 = qJD(3) * t261 + qJD(4) * t226;
t265 = t340 * t284;
t264 = t340 * t282;
t236 = t239 ^ 2;
t224 = t284 * t229;
t216 = t284 * t220;
t191 = pkin(5) * t239 + qJD(6) + t198;
t190 = -qJ(6) * t328 + t319;
t189 = pkin(5) * t260 - qJ(6) * t327 - t226 * t282 + t216;
t180 = -qJ(6) * t303 + (-qJD(5) * t226 + t295) * t282 + t307;
t179 = pkin(5) * t256 - t204 * t282 + t224 + t295 * t284 + (-t219 + (qJ(6) * t261 - t220) * t282) * qJD(5);
t1 = [(-t254 * t255 - t261 * t294) * MDP(9) + (-t261 * t247 - t255 * t253 - t254 * t256 + t260 * t294) * MDP(10) + (t247 * t262 + t250 * t256) * MDP(14) + (-t250 * t255 - t262 * t294) * MDP(15) + (-t209 * t327 + t241 * t292) * MDP(16) + (t308 + (t337 + (-t332 + t334) * qJD(5)) * t261 + t323) * MDP(17) + (t244 * t261 + t249 * t292 + t320) * MDP(18) + (-t249 * t293 - t261 * t326 + t298) * MDP(19) + (t247 * t260 + t249 * t256) * MDP(20) + ((-t226 * t311 + t224) * t249 + t216 * t247 + (-t199 * t311 + t213) * t260 + t187 * t256 + t205 * t239 + t225 * t210 + t198 * t303 + ((-qJD(5) * t220 - t204) * t249 - t226 * t247 + (-qJD(5) * t208 - t193) * t260 + t194 * t261 - t198 * t255) * t282) * MDP(21) + (-(-t226 * t312 + t307) * t249 - t319 * t247 - t290 * t260 - t188 * t256 + t205 * t241 - t225 * t209 + t194 * t327 + t292 * t198) * MDP(22) + (-t179 * t241 - t180 * t239 + t189 * t209 - t190 * t210 - t299 * t255 + (-t177 * t284 - t178 * t282 + (t181 * t282 - t183 * t284) * qJD(5)) * t261) * MDP(23) + (t178 * t190 + t183 * t180 + t177 * t189 + t181 * t179 + t185 * (pkin(5) * t328 + t225) + t191 * (pkin(5) * t293 + t205)) * MDP(24) + (0.2e1 * t351 + (t270 * t301 - t296) * MDP(8)) * qJD(3) + (-MDP(11) * t255 - MDP(12) * t256 - MDP(14) * t205 - MDP(15) * t204) * qJD(4); (t255 * t325 - t298) * MDP(21) + (t249 * t330 + t320) * MDP(22) + (-t308 + t323) * MDP(23) + (t181 * t331 - t183 * t330 + t185 * t260 + t191 * t256) * MDP(24) + (-MDP(14) * t256 + MDP(15) * t255) * qJD(4) + (-MDP(23) * t337 + (-t177 * t282 + t178 * t284) * MDP(24) + (-MDP(21) * t282 - MDP(22) * t284) * t247 + ((t332 + t334) * MDP(23) + t299 * MDP(24) + (-MDP(21) * t284 + MDP(22) * t282) * t249) * qJD(5)) * t261; t268 * MDP(15) + t318 * MDP(21) + t321 * MDP(23) + (-MDP(21) * t239 - MDP(22) * t241 - MDP(24) * t191) * t254 + (-qJD(5) * t249 * MDP(21) - t247 * MDP(22) + MDP(23) * t333 + MDP(24) * t349) * t282 + ((t209 + t335) * MDP(23) + t350 * MDP(24) - t249 ^ 2 * MDP(22)) * t284 + ((t254 + t347) * MDP(14) + (t253 - t306) * MDP(15)) * qJD(4) + (t296 * MDP(8) - t351) * qJD(1); -t253 ^ 2 * MDP(10) + (t268 + (-t253 - t306) * qJD(4)) * MDP(11) + (-t194 + t352) * MDP(14) + (-t250 * t253 + t343) * MDP(15) + (t241 * t300 - t337) * MDP(16) + ((-t209 + t335) * t284 - t241 * t325 + t321) * MDP(17) + (t249 * t300 + t326) * MDP(18) + (-t249 * t312 + t318) * MDP(19) + (-pkin(4) * t210 - t194 * t284 - t202 * t239 + (-pkin(8) * t311 - t223) * t249 + (t249 * t344 + t291) * t282) * MDP(21) + (pkin(4) * t209 + t194 * t282 - t202 * t241 + (pkin(8) * t312 + t322) * t249 + t291 * t284) * MDP(22) + (t209 * t264 + t210 * t265 - t317 * t239 - t316 * t241 - t282 * t350 + t284 * t349) * MDP(23) + (-t178 * t265 + t177 * t264 + t185 * (-pkin(5) * t284 - pkin(4)) + (pkin(5) * t325 - t202) * t191 + t317 * t183 + t316 * t181) * MDP(24) + (MDP(10) * t254 - MDP(14) * t250 - MDP(18) * t241 + MDP(19) * t239 - MDP(20) * t249 - MDP(21) * t187 + MDP(22) * t188 - MDP(9) * t253) * t254; t241 * t239 * MDP(16) + (-t236 + t342) * MDP(17) + (t239 * t249 - t209) * MDP(18) + (-t210 + t333) * MDP(19) + t247 * MDP(20) + (t188 * t249 - t198 * t241 + t287) * MDP(21) + (t187 * t249 + t198 * t239 - t290) * MDP(22) + (pkin(5) * t209 - t239 * t324) * MDP(23) + (t324 * t183 + (-t191 * t241 + t177) * pkin(5)) * MDP(24); (-t236 - t342) * MDP(23) + (t181 * t241 + t183 * t239 + t185) * MDP(24);];
tauc  = t1;
