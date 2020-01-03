% Calculate minimal parameter regressor of Coriolis joint torque vector for
% S5RRPRP7
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% qJD [5x1]
%   Generalized joint velocities
% pkin [8x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d1,d2,d4,theta3]';
% MDP [23x1]
%   Minimal dynamic parameter vector (fixed base model)
%   see S5RRPRP7_convert_par2_MPV_fixb.m
% 
% Output:
% tauc [5x1]
%   joint torques required to compensate Coriolis and centrifugal load

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-31 20:02
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauc = S5RRPRP7_coriolisvecJ_fixb_mdp_slag_vp(qJ, qJD, pkin, MDP)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(8,1),zeros(23,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RRPRP7_coriolisvecJ_fixb_mdp_slag_vp: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RRPRP7_coriolisvecJ_fixb_mdp_slag_vp: qJD has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S5RRPRP7_coriolisvecJ_fixb_mdp_slag_vp: pkin has to be [8x1] (double)');
assert(isreal(MDP) && all(size(MDP) == [23 1]), ...
  'S5RRPRP7_coriolisvecJ_fixb_mdp_slag_vp: MDP has to be [23x1] (double)'); 

%% Symbolic Calculation
% From coriolisvec_joint_fixb_mdp_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 20:01:37
% EndTime: 2019-12-31 20:01:43
% DurationCPUTime: 2.52s
% Computational Cost: add. (2777->326), mult. (7015->435), div. (0->0), fcn. (4811->6), ass. (0->130)
t273 = sin(pkin(8));
t276 = sin(qJ(2));
t311 = qJD(1) * t276;
t274 = cos(pkin(8));
t278 = cos(qJ(2));
t322 = t274 * t278;
t245 = qJD(1) * t322 - t273 * t311;
t241 = qJD(4) - t245;
t306 = qJD(1) * qJD(2);
t340 = -0.2e1 * t306;
t339 = t276 * MDP(4);
t338 = (t276 ^ 2 - t278 ^ 2) * MDP(5);
t255 = t273 * t276 - t322;
t256 = t273 * t278 + t274 * t276;
t304 = -pkin(2) * t278 - pkin(1);
t219 = pkin(3) * t255 - pkin(7) * t256 + t304;
t334 = -qJ(3) - pkin(6);
t261 = t334 * t276;
t262 = t334 * t278;
t227 = t261 * t273 - t262 * t274;
t275 = sin(qJ(4));
t277 = cos(qJ(4));
t314 = t275 * t219 + t277 * t227;
t259 = qJD(1) * t262;
t250 = t273 * t259;
t258 = qJD(1) * t261;
t333 = qJD(2) * pkin(2);
t253 = t258 + t333;
t220 = t253 * t274 + t250;
t215 = -qJD(2) * pkin(3) - t220;
t247 = t256 * qJD(1);
t308 = t277 * qJD(2);
t228 = t247 * t275 - t308;
t230 = qJD(2) * t275 + t247 * t277;
t186 = pkin(4) * t228 - qJ(5) * t230 + t215;
t246 = t256 * qJD(2);
t239 = qJD(1) * t246;
t267 = pkin(2) * t273 + pkin(7);
t324 = t267 * t239;
t337 = t241 * t186 - t324;
t301 = t278 * t306;
t302 = t276 * t306;
t289 = -t273 * t302 + t274 * t301;
t202 = qJD(4) * t230 + t275 * t289;
t336 = t230 ^ 2;
t335 = pkin(4) * t239;
t332 = qJ(5) * t239;
t300 = qJD(2) * t334;
t242 = qJD(3) * t278 + t276 * t300;
t237 = t242 * qJD(1);
t243 = -t276 * qJD(3) + t278 * t300;
t281 = qJD(1) * t243;
t199 = t237 * t273 - t274 * t281;
t310 = qJD(4) * t275;
t201 = -qJD(4) * t308 + t247 * t310 - t277 * t289;
t177 = pkin(4) * t202 + qJ(5) * t201 - qJD(5) * t230 + t199;
t331 = t177 * t275;
t294 = t304 * qJD(1);
t260 = qJD(3) + t294;
t204 = -pkin(3) * t245 - pkin(7) * t247 + t260;
t323 = t274 * t259;
t221 = t273 * t253 - t323;
t216 = qJD(2) * pkin(7) + t221;
t185 = t204 * t275 + t216 * t277;
t330 = t185 * t241;
t329 = t201 * t275;
t328 = t228 * t230;
t327 = t228 * t245;
t296 = t230 * t241;
t326 = t241 * t275;
t325 = t256 * t277;
t234 = t275 * t239;
t279 = qJD(2) ^ 2;
t321 = t276 * t279;
t235 = t277 * t239;
t320 = t278 * t279;
t280 = qJD(1) ^ 2;
t319 = t278 * t280;
t222 = t258 * t273 - t323;
t292 = pkin(4) * t275 - qJ(5) * t277;
t318 = qJD(5) * t275 - t241 * t292 + t222;
t309 = qJD(4) * t277;
t317 = -t275 * t202 - t228 * t309;
t316 = t245 * t326 + t235;
t211 = pkin(2) * t311 + pkin(3) * t247 - pkin(7) * t245;
t223 = t258 * t274 + t250;
t315 = t275 * t211 + t277 * t223;
t313 = t241 * t309 + t234;
t184 = t204 * t277 - t216 * t275;
t307 = qJD(5) - t184;
t305 = t276 * t333;
t268 = -pkin(2) * t274 - pkin(3);
t303 = t267 * t310;
t299 = pkin(1) * t340;
t200 = t274 * t237 + t273 * t281;
t265 = pkin(2) * t302;
t203 = t239 * pkin(3) - pkin(7) * t289 + t265;
t285 = t277 * t200 + t275 * t203 + t204 * t309 - t216 * t310;
t173 = qJD(5) * t241 + t285 + t332;
t179 = -pkin(4) * t241 + t307;
t298 = -t179 * t245 + t173;
t295 = t275 * t200 - t277 * t203 + t204 * t310 + t216 * t309;
t174 = t295 - t335;
t180 = qJ(5) * t241 + t185;
t297 = t180 * t245 + t174;
t209 = t242 * t273 - t274 * t243;
t226 = -t274 * t261 - t262 * t273;
t293 = pkin(4) * t277 + qJ(5) * t275;
t291 = t179 * t277 - t180 * t275;
t290 = t199 * t256 - t227 * t239;
t249 = t255 * qJD(2);
t288 = -t249 * t275 + t256 * t309;
t287 = t249 * t277 + t256 * t310;
t286 = t186 * t230 + t295;
t210 = t242 * t274 + t243 * t273;
t212 = pkin(3) * t246 + pkin(7) * t249 + t305;
t284 = t277 * t210 + t275 * t212 + t219 * t309 - t227 * t310;
t283 = t215 * t241 - t324;
t254 = t268 - t293;
t191 = pkin(4) * t230 + qJ(5) * t228;
t190 = t256 * t292 + t226;
t188 = -pkin(4) * t255 - t219 * t277 + t227 * t275;
t187 = qJ(5) * t255 + t314;
t183 = t228 * t241 - t201;
t182 = -pkin(4) * t247 - t211 * t277 + t223 * t275;
t181 = qJ(5) * t247 + t315;
t178 = -t292 * t249 + (qJD(4) * t293 - qJD(5) * t277) * t256 + t209;
t176 = -pkin(4) * t246 + qJD(4) * t314 + t210 * t275 - t212 * t277;
t175 = qJ(5) * t246 + qJD(5) * t255 + t284;
t1 = [0.2e1 * t301 * t339 + t338 * t340 + MDP(6) * t320 - MDP(7) * t321 + (-pkin(6) * t320 + t276 * t299) * MDP(9) + (pkin(6) * t321 + t278 * t299) * MDP(10) + (-t200 * t255 + t209 * t247 + t210 * t245 + t220 * t249 - t221 * t246 + t226 * t289 + t290) * MDP(11) + (t199 * t226 + t200 * t227 - t220 * t209 + t221 * t210 + (t260 + t294) * t305) * MDP(12) + (-t201 * t325 - t230 * t287) * MDP(13) + (-(-t228 * t277 - t230 * t275) * t249 + (t329 - t202 * t277 + (t228 * t275 - t230 * t277) * qJD(4)) * t256) * MDP(14) + (-t201 * t255 + t230 * t246 + t256 * t235 - t241 * t287) * MDP(15) + (-t202 * t255 - t228 * t246 - t256 * t234 - t241 * t288) * MDP(16) + (t239 * t255 + t241 * t246) * MDP(17) + (-t295 * t255 + t184 * t246 + t209 * t228 + t226 * t202 + ((-qJD(4) * t227 + t212) * t241 + t219 * t239 + t215 * qJD(4) * t256) * t277 + ((-qJD(4) * t219 - t210) * t241 - t215 * t249 + t290) * t275) * MDP(18) + (-t185 * t246 + t199 * t325 - t226 * t201 + t209 * t230 - t287 * t215 - t314 * t239 - t284 * t241 - t285 * t255) * MDP(19) + (-t174 * t255 - t176 * t241 + t178 * t228 - t179 * t246 + t186 * t288 - t188 * t239 + t190 * t202 + t256 * t331) * MDP(20) + (-t175 * t228 + t176 * t230 - t187 * t202 - t188 * t201 - t291 * t249 + (-t173 * t275 + t174 * t277 + (-t179 * t275 - t180 * t277) * qJD(4)) * t256) * MDP(21) + (t173 * t255 + t175 * t241 - t177 * t325 - t178 * t230 + t180 * t246 + t186 * t287 + t187 * t239 + t190 * t201) * MDP(22) + (t173 * t187 + t174 * t188 + t175 * t180 + t176 * t179 + t177 * t190 + t178 * t186) * MDP(23); -t319 * t339 + t280 * t338 + ((t221 - t222) * t247 + (-t223 + t220) * t245 + (-t273 * t239 - t274 * t289) * pkin(2)) * MDP(11) + (t220 * t222 - t221 * t223 + (-t199 * t274 + t200 * t273 - t260 * t311) * pkin(2)) * MDP(12) + (t277 * t296 - t329) * MDP(13) + ((-t201 + t327) * t277 - t230 * t326 + t317) * MDP(14) + (-t241 * t245 * t277 - t230 * t247 + t313) * MDP(15) + (t228 * t247 - t241 * t310 + t316) * MDP(16) - t241 * t247 * MDP(17) + (-t184 * t247 + t268 * t202 - t222 * t228 + (-t199 + (-qJD(4) * t267 - t211) * t241) * t277 + (t223 * t241 + t283) * t275) * MDP(18) + (t185 * t247 + t199 * t275 - t268 * t201 - t222 * t230 + (t303 + t315) * t241 + t283 * t277) * MDP(19) + (-t177 * t277 + t179 * t247 + t202 * t254 + (-t267 * t309 + t182) * t241 - t318 * t228 + t337 * t275) * MDP(20) + (t181 * t228 - t182 * t230 + (-t202 * t267 + (t230 * t267 + t179) * qJD(4) + t298) * t277 + (-t201 * t267 + (t228 * t267 - t180) * qJD(4) + t297) * t275) * MDP(21) + (-t331 - t180 * t247 + t201 * t254 + (-t181 - t303) * t241 + t318 * t230 - t337 * t277) * MDP(22) + (t177 * t254 - t179 * t182 - t180 * t181 - t318 * t186 + (qJD(4) * t291 + t173 * t277 + t174 * t275) * t267) * MDP(23) + (t280 * t276 * MDP(9) + MDP(10) * t319) * pkin(1); -t245 ^ 2 * MDP(11) + (-t221 * t245 + t265) * MDP(12) + t316 * MDP(18) + t317 * MDP(21) + t313 * MDP(22) + (-MDP(11) * t247 + t220 * MDP(12) - t186 * MDP(23) + (-MDP(19) + MDP(22)) * t230 + (-MDP(18) - MDP(20)) * t228) * t247 + (t239 * MDP(20) + (t201 + t327) * MDP(21) + (qJD(4) * t180 - t297) * MDP(23) + (-t241 * MDP(19) - t245 * MDP(22)) * t241) * t277 + (-t239 * MDP(19) + (qJD(4) * t179 + t298) * MDP(23) + MDP(21) * t296 + (-qJD(4) * MDP(18) - t241 * MDP(20)) * t241) * t275; MDP(13) * t328 + (-t228 ^ 2 + t336) * MDP(14) + t183 * MDP(15) + (-t202 + t296) * MDP(16) + t239 * MDP(17) + (-t215 * t230 - t295 + t330) * MDP(18) + (t184 * t241 + t215 * t228 - t285) * MDP(19) + (-t191 * t228 - t286 + t330 + 0.2e1 * t335) * MDP(20) + (pkin(4) * t201 - qJ(5) * t202 + (t180 - t185) * t230 + (t179 - t307) * t228) * MDP(21) + (0.2e1 * t332 - t186 * t228 + t191 * t230 + (0.2e1 * qJD(5) - t184) * t241 + t285) * MDP(22) + (-pkin(4) * t174 + qJ(5) * t173 - t179 * t185 + t180 * t307 - t186 * t191) * MDP(23); (-qJD(2) * t247 + t328) * MDP(20) + t183 * MDP(21) + (-t241 ^ 2 - t336) * MDP(22) + (-t180 * t241 + t286 - t335) * MDP(23);];
tauc = t1;
