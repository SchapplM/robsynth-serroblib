% Calculate minimal parameter regressor of Coriolis joint torque vector for
% S5RRPPR5
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% qJD [5x1]
%   Generalized joint velocities
% pkin [8x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d1,d2,d5,theta3]';
% MDP [23x1]
%   Minimal dynamic parameter vector (fixed base model)
%   see S5RRPPR5_convert_par2_MPV_fixb.m
% 
% Output:
% tauc [5x1]
%   joint torques required to compensate Coriolis and centrifugal load

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-31 19:30
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauc = S5RRPPR5_coriolisvecJ_fixb_mdp_slag_vp(qJ, qJD, pkin, MDP)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(8,1),zeros(23,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RRPPR5_coriolisvecJ_fixb_mdp_slag_vp: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RRPPR5_coriolisvecJ_fixb_mdp_slag_vp: qJD has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S5RRPPR5_coriolisvecJ_fixb_mdp_slag_vp: pkin has to be [8x1] (double)');
assert(isreal(MDP) && all(size(MDP) == [23 1]), ...
  'S5RRPPR5_coriolisvecJ_fixb_mdp_slag_vp: MDP has to be [23x1] (double)'); 

%% Symbolic Calculation
% From coriolisvec_joint_fixb_mdp_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 19:30:01
% EndTime: 2019-12-31 19:30:06
% DurationCPUTime: 2.09s
% Computational Cost: add. (1241->235), mult. (3234->308), div. (0->0), fcn. (2221->6), ass. (0->120)
t279 = sin(pkin(8));
t280 = cos(pkin(8));
t282 = sin(qJ(2));
t284 = cos(qJ(2));
t257 = t279 * t284 + t280 * t282;
t246 = t257 * qJD(2);
t237 = qJD(1) * t246;
t311 = qJD(1) * qJD(2);
t307 = t282 * t311;
t266 = t279 * t307;
t306 = t284 * t311;
t238 = t280 * t306 - t266;
t281 = sin(qJ(5));
t283 = cos(qJ(5));
t317 = qJD(1) * t284;
t308 = t280 * t317;
t314 = t282 * qJD(1);
t244 = t279 * t314 - t308;
t247 = t257 * qJD(1);
t349 = t283 * t244 - t247 * t281;
t176 = t349 * qJD(5) + t281 * t237 + t283 * t238;
t340 = t244 * t281 + t283 * t247;
t177 = qJD(5) * t340 - t283 * t237 + t238 * t281;
t274 = qJD(2) - qJD(5);
t329 = t349 * t274;
t330 = t340 * t274;
t354 = -(t177 + t330) * MDP(20) - t349 * t340 * MDP(17) + (t340 ^ 2 - t349 ^ 2) * MDP(18) + (t176 + t329) * MDP(19);
t309 = pkin(2) * t284 + pkin(1);
t301 = t309 * qJD(1);
t263 = qJD(3) - t301;
t352 = -qJ(4) * t247 + t263;
t333 = -qJ(3) - pkin(6);
t305 = qJD(2) * t333;
t240 = qJD(3) * t284 + t282 * t305;
t227 = t240 * qJD(1);
t241 = -qJD(3) * t282 + t284 * t305;
t228 = t241 * qJD(1);
t195 = t280 * t227 + t279 * t228;
t275 = qJD(2) * qJD(4);
t192 = t275 + t195;
t181 = pkin(7) * t237 + t192;
t194 = t227 * t279 - t280 * t228;
t182 = -pkin(7) * t238 + t194;
t337 = -pkin(3) - pkin(4);
t184 = t244 * t337 - t352;
t347 = t281 * t181 - t283 * t182 + t184 * t340;
t345 = -0.2e1 * t311;
t344 = MDP(4) * t282;
t343 = MDP(5) * (t282 ^ 2 - t284 ^ 2);
t264 = t333 * t282;
t261 = qJD(1) * t264;
t265 = t333 * t284;
t262 = qJD(1) * t265;
t326 = t279 * t262;
t218 = t280 * t261 + t326;
t312 = qJD(4) - t218;
t339 = qJD(5) + t274;
t338 = -t283 * t181 - t281 * t182 - t184 * t349;
t243 = t247 ^ 2;
t336 = pkin(3) * t237;
t335 = pkin(7) * t244;
t334 = pkin(7) * t247;
t219 = -t280 * t264 - t265 * t279;
t331 = t194 * t219;
t325 = t280 * t262;
t324 = t280 * t284;
t285 = qJD(2) ^ 2;
t323 = t282 * t285;
t320 = t284 * t285;
t286 = qJD(1) ^ 2;
t319 = t284 * t286;
t202 = t280 * t240 + t279 * t241;
t255 = qJD(2) * pkin(2) + t261;
t214 = t279 * t255 - t325;
t220 = t279 * t264 - t280 * t265;
t316 = qJD(2) * t282;
t313 = -t334 + t312;
t310 = pkin(2) * t316;
t211 = qJD(2) * qJ(4) + t214;
t273 = -pkin(2) * t280 - pkin(3);
t304 = pkin(1) * t345;
t270 = pkin(2) * t307;
t303 = qJ(4) * t238 - t270;
t201 = t240 * t279 - t280 * t241;
t213 = t255 * t280 + t326;
t217 = t261 * t279 - t325;
t300 = qJD(4) - t213;
t185 = qJD(2) * t337 + t300 - t334;
t191 = t211 + t335;
t297 = t283 * t185 - t281 * t191;
t296 = -t281 * t185 - t283 * t191;
t256 = t279 * t282 - t324;
t293 = t256 * t283 - t257 * t281;
t216 = t256 * t281 + t257 * t283;
t292 = qJ(4) * t257 + t309;
t291 = -pkin(2) * t314 - qJ(4) * t244;
t290 = qJD(4) * t247 + t303;
t249 = qJD(2) * t324 - t279 * t316;
t289 = qJ(4) * t249 + qJD(4) * t257 - t310;
t288 = t194 * t257 + t201 * t247 - t202 * t244 + t219 * t238 - t220 * t237;
t271 = pkin(2) * t279 + qJ(4);
t269 = -pkin(4) + t273;
t212 = pkin(3) * t256 - t292;
t205 = -qJD(2) * pkin(3) + t300;
t204 = pkin(7) * t256 + t220;
t203 = -pkin(7) * t257 + t219;
t200 = pkin(3) * t247 - t291;
t199 = pkin(3) * t244 + t352;
t196 = t217 + t335;
t193 = t256 * t337 + t292;
t190 = pkin(3) * t246 - t289;
t188 = pkin(7) * t246 + t202;
t187 = -pkin(7) * t249 + t201;
t186 = t247 * t337 + t291;
t183 = -t290 + t336;
t180 = t246 * t337 + t289;
t179 = qJD(5) * t216 - t283 * t246 + t249 * t281;
t178 = qJD(5) * t293 + t246 * t281 + t249 * t283;
t175 = t237 * t337 + t290;
t1 = [0.2e1 * t306 * t344 + t343 * t345 + MDP(6) * t320 - MDP(7) * t323 + (-pkin(6) * t320 + t282 * t304) * MDP(9) + (pkin(6) * t323 + t284 * t304) * MDP(10) + (-t195 * t256 - t213 * t249 - t214 * t246 + t288) * MDP(11) + (t331 + t195 * t220 - t213 * t201 + t214 * t202 + (t263 - t301) * t310) * MDP(12) + (-qJD(2) * t201 + t183 * t256 + t190 * t244 + t199 * t246 + t212 * t237) * MDP(13) + (-t192 * t256 + t205 * t249 - t211 * t246 + t288) * MDP(14) + (qJD(2) * t202 - t183 * t257 - t190 * t247 - t199 * t249 - t212 * t238) * MDP(15) + (t183 * t212 + t190 * t199 + t192 * t220 + t201 * t205 + t202 * t211 + t331) * MDP(16) + (t176 * t216 + t178 * t340) * MDP(17) + (t176 * t293 - t177 * t216 + t178 * t349 - t179 * t340) * MDP(18) + (-t175 * t293 + t193 * t177 + t184 * t179 - t180 * t349) * MDP(22) + (t175 * t216 + t193 * t176 + t184 * t178 + t180 * t340) * MDP(23) + (-t178 * MDP(19) + t179 * MDP(20) + (-t187 * t283 + t188 * t281 - (-t203 * t281 - t204 * t283) * qJD(5)) * MDP(22) + (t187 * t281 + t188 * t283 + (t203 * t283 - t204 * t281) * qJD(5)) * MDP(23)) * t274; -t319 * t344 + t286 * t343 + ((t214 - t217) * t247 + (-t213 + t218) * t244 + (-t237 * t279 - t238 * t280) * pkin(2)) * MDP(11) + (t213 * t217 - t214 * t218 + (-t194 * t280 + t195 * t279 - t263 * t314) * pkin(2)) * MDP(12) + (qJD(2) * t217 - t199 * t247 - t200 * t244 - t194) * MDP(13) + (-t237 * t271 + t238 * t273 + (t211 - t217) * t247 + (t205 - t312) * t244) * MDP(14) + (-qJD(2) * t218 - t199 * t244 + t200 * t247 + t195 + 0.2e1 * t275) * MDP(15) + (t192 * t271 + t194 * t273 - t199 * t200 - t205 * t217 + t211 * t312) * MDP(16) + (t186 * t349 + (t283 * t196 + t313 * t281) * t274 + (-(-t269 * t281 - t271 * t283) * t274 - t296) * qJD(5) + t347) * MDP(22) + (-t186 * t340 + (-t281 * t196 + t313 * t283) * t274 + ((t269 * t283 - t271 * t281) * t274 + t297) * qJD(5) - t338) * MDP(23) + (MDP(9) * t282 * t286 + MDP(10) * t319) * pkin(1) - t354; (t213 * t247 + t214 * t244 + t270) * MDP(12) + t266 * MDP(15) + (t336 + t211 * t244 + (-qJD(4) - t205) * t247 - t303) * MDP(16) + (-t177 + t330) * MDP(22) + (-t176 + t329) * MDP(23) + ((t279 * t317 + t280 * t314 + t247) * MDP(13) + (t244 - t308) * MDP(15)) * qJD(2) + (MDP(11) + MDP(14)) * (-t244 ^ 2 - t243); (-t266 + (t244 + t308) * qJD(2)) * MDP(14) + (-t243 - t285) * MDP(15) + (-qJD(2) * t211 + t194) * MDP(16) + (t244 * MDP(13) + t199 * MDP(16) + MDP(22) * t349 - MDP(23) * t340) * t247 + (-t281 * MDP(22) - t283 * MDP(23)) * t274 ^ 2; (t339 * t296 - t347) * MDP(22) + (-t339 * t297 + t338) * MDP(23) + t354;];
tauc = t1;
