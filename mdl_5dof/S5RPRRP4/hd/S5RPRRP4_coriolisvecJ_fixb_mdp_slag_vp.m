% Calculate minimal parameter regressor of Coriolis joint torque vector for
% S5RPRRP4
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% qJD [5x1]
%   Generalized joint velocities
% pkin [8x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d1,d3,d4,theta2]';
% MDP [23x1]
%   Minimal dynamic parameter vector (fixed base model)
%   see S5RPRRP4_convert_par2_MPV_fixb.m
% 
% Output:
% tauc [5x1]
%   joint torques required to compensate Coriolis and centrifugal load

% Quelle: HybrDyn-Toolbox
% Datum: 2020-01-03 11:51
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauc = S5RPRRP4_coriolisvecJ_fixb_mdp_slag_vp(qJ, qJD, pkin, MDP)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(8,1),zeros(23,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RPRRP4_coriolisvecJ_fixb_mdp_slag_vp: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RPRRP4_coriolisvecJ_fixb_mdp_slag_vp: qJD has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S5RPRRP4_coriolisvecJ_fixb_mdp_slag_vp: pkin has to be [8x1] (double)');
assert(isreal(MDP) && all(size(MDP) == [23 1]), ...
  'S5RPRRP4_coriolisvecJ_fixb_mdp_slag_vp: MDP has to be [23x1] (double)'); 

%% Symbolic Calculation
% From coriolisvec_joint_fixb_mdp_matlab.m
% OptimizationMode: 2
% StartTime: 2020-01-03 11:50:29
% EndTime: 2020-01-03 11:50:34
% DurationCPUTime: 1.70s
% Computational Cost: add. (1604->210), mult. (4205->298), div. (0->0), fcn. (2876->6), ass. (0->122)
t257 = sin(pkin(8));
t296 = qJD(3) + qJD(4);
t342 = t257 * t296;
t260 = sin(qJ(3));
t262 = cos(qJ(3));
t341 = (t260 * MDP(10) + t262 * MDP(11)) * t257;
t253 = t257 ^ 2;
t331 = 0.2e1 * t253;
t259 = sin(qJ(4));
t261 = cos(qJ(4));
t232 = t259 * t262 + t260 * t261;
t223 = t232 * t257;
t340 = (t260 ^ 2 - t262 ^ 2) * MDP(9);
t319 = t261 * t262;
t320 = t259 * t260;
t231 = t319 - t320;
t258 = cos(pkin(8));
t306 = qJD(1) * t258;
t316 = (t296 - t306) * t231;
t211 = t296 * t232;
t267 = qJD(1) * t232;
t339 = t258 * t267 - t211;
t307 = qJD(1) * t257;
t292 = t260 * t307;
t279 = t259 * t292;
t291 = t262 * t307;
t280 = t261 * t291;
t219 = -t279 + t280;
t213 = t219 * qJ(5);
t235 = -pkin(2) * t258 - pkin(6) * t257 - pkin(1);
t226 = t235 * qJD(1) + qJD(2);
t222 = t262 * t226;
t321 = t258 * t260;
t330 = pkin(7) * t257;
t266 = -qJ(2) * t321 - t262 * t330;
t207 = t266 * qJD(1) + t222;
t246 = -qJD(3) + t306;
t193 = -pkin(3) * t246 + t207;
t327 = qJ(2) * t262;
t294 = t258 * t327;
t208 = -pkin(7) * t292 + qJD(1) * t294 + t226 * t260;
t200 = t259 * t208;
t286 = t261 * t193 - t200;
t338 = t213 - t286;
t254 = t258 ^ 2;
t337 = t331 + t254;
t336 = -MDP(13) * t260 - MDP(14) * t262;
t288 = -t235 + t330;
t333 = t288 * t260 - t294;
t332 = t219 ^ 2;
t329 = MDP(7) * qJ(2);
t328 = qJ(2) * t260;
t216 = t257 * t267;
t326 = qJ(5) * t216;
t227 = pkin(3) * t292 + qJ(2) * t307;
t205 = pkin(4) * t216 + qJD(5) + t227;
t325 = t205 * t219;
t324 = t227 * t219;
t263 = qJD(1) ^ 2;
t323 = t253 * t263;
t322 = t257 * t260;
t202 = t261 * t208;
t241 = -qJD(4) + t246;
t174 = -pkin(4) * t241 - t338;
t318 = t174 + t338;
t317 = t261 * t207 - t200;
t298 = qJD(1) * qJD(2);
t289 = t258 * t298;
t302 = qJD(3) * t262;
t314 = t226 * t302 + t262 * t289;
t304 = qJD(2) * t262;
t313 = t235 * t302 + t258 * t304;
t312 = t296 * t279;
t244 = t257 * pkin(3) * t302;
t225 = qJD(1) * t244 + t257 * t298;
t229 = t257 * qJD(2) + t244;
t230 = pkin(3) * t322 + t257 * qJ(2);
t311 = t253 + t254;
t305 = qJD(2) * t260;
t303 = qJD(3) * t260;
t301 = qJD(4) * t259;
t300 = qJD(4) * t261;
t299 = qJD(3) + t246;
t293 = qJ(2) * t303;
t290 = t258 * t305;
t265 = t266 * qJD(3);
t188 = qJD(1) * t265 + t314;
t189 = -t226 * t303 + (-t290 + (pkin(7) * t322 - t294) * qJD(3)) * qJD(1);
t287 = -t259 * t188 + t261 * t189;
t203 = t265 + t313;
t204 = t333 * qJD(3) - t290;
t285 = -t203 * t259 + t261 * t204;
t284 = -t207 * t259 - t202;
t283 = qJD(1) * t299;
t282 = t262 * t253 * t260 * MDP(8);
t281 = -t261 * t188 - t259 * t189 - t193 * t300 + t208 * t301;
t273 = t319 * t342;
t192 = qJD(1) * t273 - t312;
t278 = pkin(4) * t192 + t225;
t276 = -t193 * t259 - t202;
t209 = -t288 * t262 + (-pkin(3) - t328) * t258;
t275 = -t209 * t259 + t261 * t333;
t270 = t227 * t216 + t281;
t214 = t216 ^ 2;
t269 = t219 * t216 * MDP(15) + (-t211 * t307 - t216 * t241) * MDP(17) + (-t219 * t241 - t296 * t280 + t312) * MDP(18) + (-t214 + t332) * MDP(16);
t268 = t261 * t203 + t259 * t204 + t209 * t300 + t301 * t333;
t264 = t276 * qJD(4) + t287;
t198 = t296 * t223;
t250 = pkin(3) * t261 + pkin(4);
t224 = t231 * t257;
t199 = -t320 * t342 + t273;
t191 = qJD(1) * t198;
t180 = -qJ(5) * t223 - t275;
t179 = -pkin(4) * t258 - qJ(5) * t224 + t209 * t261 + t259 * t333;
t178 = -t213 + t317;
t177 = t284 + t326;
t176 = -t276 - t326;
t173 = qJ(5) * t198 + t275 * qJD(4) - qJD(5) * t224 + t285;
t172 = -qJ(5) * t199 - qJD(5) * t223 + t268;
t171 = qJ(5) * t191 - qJD(5) * t219 + t264;
t170 = -qJ(5) * t192 - qJD(5) * t216 - t281;
t1 = [(t246 * t290 + t337 * qJD(1) * (qJ(2) * t302 + t305)) * MDP(13) + ((-t258 * t293 + t313) * t246 + t314 * t258 + (-t337 * t293 + t304 * t331) * qJD(1)) * MDP(14) + (-t191 * t224 - t198 * t219) * MDP(15) + (t191 * t223 - t192 * t224 + t198 * t216 - t199 * t219) * MDP(16) + (t191 * t258 + t198 * t241) * MDP(17) + (t192 * t258 + t199 * t241) * MDP(18) + (-t285 * t241 - t287 * t258 + t229 * t216 + t230 * t192 + t225 * t223 + t227 * t199 + (-t275 * t241 - t276 * t258) * qJD(4)) * MDP(20) + (-t230 * t191 - t227 * t198 + t229 * t219 + t225 * t224 + t268 * t241 - t281 * t258) * MDP(21) + (-t170 * t223 - t171 * t224 - t172 * t216 - t173 * t219 + t174 * t198 - t176 * t199 + t179 * t191 - t180 * t192) * MDP(22) + (t170 * t180 + t176 * t172 + t171 * t179 + t174 * t173 + t278 * (pkin(4) * t223 + t230) + t205 * (pkin(4) * t199 + t229)) * MDP(23) + 0.2e1 * (MDP(6) + t329) * t311 * t298 + ((-(-t235 * t260 - t294) * t246 + t226 * t321) * MDP(13) + (t331 * t340 - 0.2e1 * t282) * qJD(1) + (t246 + t306) * t341) * qJD(3); (t191 * t231 - t192 * t232 - t316 * t216 - t219 * t339) * MDP(22) + (t170 * t232 + t171 * t231 + t339 * t174 + t316 * t176) * MDP(23) + (-MDP(20) * t216 - MDP(21) * t219 - MDP(23) * t205) * t307 + (-MDP(20) * t339 + t316 * MDP(21)) * t241 + (-t254 * MDP(6) + (-MDP(6) + t336) * t253 - t311 * t329) * t263 + t336 * t246 ^ 2; t263 * t282 - t323 * t340 + ((-t299 * t226 - t289) * t260 + (-t258 * t283 - t323) * t327) * MDP(13) + (-t222 * t246 + (t299 * t306 + t323) * t328 - t314) * MDP(14) + (t284 * t241 - pkin(3) * t216 * t291 - t324 + (-t202 + (pkin(3) * t241 - t193) * t259) * qJD(4) + t287) * MDP(20) + (-t317 * t241 + (-t219 * t291 + t241 * t300) * pkin(3) + t270) * MDP(21) + (t191 * t250 + (t176 + t177) * t219 + (-t174 + t178) * t216 + (-t192 * t259 + (-t216 * t261 + t219 * t259) * qJD(4)) * pkin(3)) * MDP(22) + (-pkin(4) * t325 + t171 * t250 - t174 * t177 - t176 * t178 + (-t205 * t291 + t170 * t259 + (-t174 * t259 + t176 * t261) * qJD(4)) * pkin(3)) * MDP(23) + t269 - t283 * t341; (t276 * t241 + t264 - t324) * MDP(20) + (-t286 * t241 + t270) * MDP(21) + (pkin(4) * t191 - t318 * t216) * MDP(22) + (t318 * t176 + (t171 - t325) * pkin(4)) * MDP(23) + t269; (-t214 - t332) * MDP(22) + (t174 * t219 + t176 * t216 + t278) * MDP(23);];
tauc = t1;
