% Calculate Coriolis joint torque vector for
% S5PRRPR4
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% qJD [5x1]
%   Generalized joint velocities
% pkin [9x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d2,d3,d5,theta1,theta4]';
% MDP [22x1]
%   Minimal dynamic parameter vector (fixed base model)
%   see S5PRRPR4_convert_par2_MPV_fixb.m
% 
% Output:
% tauc [5x1]
%   joint torques required to compensate Coriolis and centrifugal load

% Quelle: HybrDyn-Toolbox
% Datum: 2021-01-15 15:53
% Revision: 24b2e7d74a0c1a3b64fa2f8f5ad758691ad61af3 (2021-01-15)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauc = S5PRRPR4_coriolisvecJ_fixb_mdp_slag_vp(qJ, qJD, pkin, MDP)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(9,1),zeros(22,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5PRRPR4_coriolisvecJ_fixb_mdp_slag_vp: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5PRRPR4_coriolisvecJ_fixb_mdp_slag_vp: qJD has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S5PRRPR4_coriolisvecJ_fixb_mdp_slag_vp: pkin has to be [9x1] (double)');
assert(isreal(MDP) && all(size(MDP) == [22 1]), ...
  'S5PRRPR4_coriolisvecJ_fixb_mdp_slag_vp: MDP has to be [22x1] (double)'); 

%% Symbolic Calculation
% From coriolisvec_joint_fixb_mdp_matlab.m
% OptimizationMode: 2
% StartTime: 2021-01-15 15:52:48
% EndTime: 2021-01-15 15:52:54
% DurationCPUTime: 2.05s
% Computational Cost: add. (1159->225), mult. (3031->325), div. (0->0), fcn. (2206->8), ass. (0->110)
t264 = sin(pkin(9));
t265 = cos(pkin(9));
t270 = cos(qJ(3));
t302 = qJD(2) * t270;
t291 = t265 * t302;
t267 = sin(qJ(3));
t303 = qJD(2) * t267;
t237 = t264 * t303 - t291;
t269 = cos(qJ(5));
t226 = t269 * t237;
t246 = t264 * t270 + t265 * t267;
t239 = t246 * qJD(3);
t229 = qJD(2) * t239;
t294 = qJD(2) * qJD(3);
t290 = t267 * t294;
t252 = t264 * t290;
t289 = t270 * t294;
t230 = t265 * t289 - t252;
t240 = t246 * qJD(2);
t266 = sin(qJ(5));
t299 = qJD(5) * t266;
t168 = -qJD(5) * t226 - t229 * t266 + t230 * t269 - t240 * t299;
t197 = -t240 * t266 - t226;
t278 = t237 * t266 - t240 * t269;
t274 = qJD(5) * t278 - t229 * t269 - t230 * t266;
t261 = qJD(3) + qJD(5);
t312 = t197 * t261;
t313 = t278 * t261;
t329 = t197 * t278 * MDP(16) + (-t197 ^ 2 + t278 ^ 2) * MDP(17) + (t168 - t312) * MDP(18) + (t274 - t313) * MDP(19);
t328 = -0.2e1 * t294;
t268 = sin(qJ(2));
t296 = t268 * qJD(1);
t253 = qJD(2) * pkin(6) + t296;
t285 = qJ(4) * qJD(2) + t253;
t233 = t285 * t267;
t218 = qJD(3) * pkin(3) - t233;
t234 = t285 * t270;
t311 = t265 * t234;
t182 = t218 * t264 + t311;
t317 = pkin(7) * t237;
t176 = t182 - t317;
t260 = -pkin(3) * t270 - pkin(2);
t271 = cos(qJ(2));
t304 = qJD(1) * t271;
t243 = qJD(2) * t260 + qJD(4) - t304;
t204 = pkin(4) * t237 + t243;
t327 = t176 * t299 - t204 * t197;
t282 = qJD(4) + t304;
t301 = qJD(3) * t267;
t199 = -t253 * t301 + (-qJ(4) * t301 + t270 * t282) * qJD(2);
t300 = qJD(3) * t270;
t200 = -t253 * t300 + (-qJ(4) * t300 - t267 * t282) * qJD(2);
t170 = -t199 * t264 + t200 * t265;
t166 = -pkin(7) * t230 + t170;
t171 = t199 * t265 + t200 * t264;
t167 = -pkin(7) * t229 + t171;
t326 = t166 * t269 - t266 * t167 + t204 * t278;
t324 = MDP(5) * t267;
t323 = (t267 ^ 2 - t270 ^ 2) * MDP(6);
t315 = qJ(4) + pkin(6);
t288 = qJD(3) * t315;
t235 = qJD(4) * t270 - t267 * t288;
t236 = -qJD(4) * t267 - t270 * t288;
t308 = -t235 * t264 + t236 * t265 + t246 * t304;
t245 = t264 * t267 - t265 * t270;
t276 = t245 * t271;
t307 = qJD(1) * t276 + t235 * t265 + t236 * t264;
t322 = pkin(3) * t301 - t296;
t321 = t267 * MDP(10) + t270 * MDP(11);
t242 = t245 * qJD(3);
t320 = qJD(5) - t261;
t319 = pkin(3) * t264;
t318 = pkin(3) * t267;
t316 = pkin(7) * t240;
t314 = qJD(2) * pkin(2);
t214 = t264 * t234;
t272 = qJD(3) ^ 2;
t310 = t267 * t272;
t309 = t270 * t272;
t184 = -t233 * t265 - t214;
t250 = t315 * t267;
t251 = t315 * t270;
t206 = -t250 * t264 + t251 * t265;
t244 = pkin(3) * t290 + qJD(2) * t296;
t293 = pkin(3) * t303;
t181 = t218 * t265 - t214;
t183 = t233 * t264 - t311;
t205 = -t250 * t265 - t251 * t264;
t284 = pkin(4) * t239 + t322;
t281 = -pkin(7) * t242 + qJD(5) * (-pkin(7) * t245 + t206) - t308;
t280 = pkin(7) * t239 - qJD(5) * (-pkin(7) * t246 + t205) - t307;
t175 = qJD(3) * pkin(4) + t181 - t316;
t279 = -t266 * t175 - t269 * t176;
t202 = t245 * t269 + t246 * t266;
t203 = -t245 * t266 + t246 * t269;
t231 = t246 * t268;
t275 = -0.2e1 * qJD(3) * t314;
t273 = qJD(2) ^ 2;
t258 = pkin(3) * t265 + pkin(4);
t232 = t245 * t268;
t222 = pkin(4) * t245 + t260;
t209 = pkin(4) * t240 + t293;
t201 = pkin(4) * t229 + t244;
t186 = -qJD(2) * t276 - qJD(3) * t231;
t185 = -t240 * t271 + t242 * t268;
t178 = t184 - t316;
t177 = t183 + t317;
t173 = qJD(5) * t203 + t239 * t269 - t242 * t266;
t172 = -qJD(5) * t202 - t239 * t266 - t242 * t269;
t1 = [(-t185 * t240 - t186 * t237 + t229 * t232 + t230 * t231) * MDP(14) + (-t170 * t231 - t171 * t232 + t181 * t185 + t182 * t186) * MDP(15) + ((t185 * t269 - t186 * t266) * MDP(21) - (t185 * t266 + t186 * t269) * MDP(22) + ((t231 * t266 + t232 * t269) * MDP(21) - (-t231 * t269 + t232 * t266) * MDP(22)) * qJD(5)) * t261 + (MDP(12) * t185 - MDP(13) * t186) * qJD(3) + (-t229 * MDP(12) - t230 * MDP(13) - t244 * MDP(15) + MDP(21) * t274 - t168 * MDP(22) - t273 * MDP(4) + t321 * t328) * t271 + (-t273 * MDP(3) + (MDP(12) * t237 + MDP(13) * t240 + MDP(15) * t243 - MDP(21) * t197 - MDP(22) * t278) * qJD(2) + (-MDP(10) * t270 + MDP(11) * t267) * (t272 + t273)) * t268; 0.2e1 * t289 * t324 + t323 * t328 + MDP(7) * t309 - MDP(8) * t310 + (-pkin(6) * t309 + t267 * t275) * MDP(10) + (pkin(6) * t310 + t270 * t275) * MDP(11) + (-t237 * t296 + t229 * t260 + t239 * t243 + t244 * t245 + (t237 * t318 + t308) * qJD(3)) * MDP(12) + (-t240 * t296 + t230 * t260 - t242 * t243 + t244 * t246 + (t240 * t318 - t307) * qJD(3)) * MDP(13) + (-t170 * t246 - t171 * t245 + t181 * t242 - t182 * t239 - t205 * t230 - t206 * t229 - t237 * t307 - t240 * t308) * MDP(14) + (t170 * t205 + t171 * t206 + t308 * t181 + t307 * t182 + t243 * t322 + t244 * t260) * MDP(15) + (t168 * t203 - t172 * t278) * MDP(16) + (-t168 * t202 + t172 * t197 + t173 * t278 + t203 * t274) * MDP(17) + (t204 * t173 - t197 * t284 + t201 * t202 - t222 * t274) * MDP(21) + (t222 * t168 + t204 * t172 + t201 * t203 - t278 * t284) * MDP(22) + (t172 * MDP(18) - t173 * MDP(19) + (t266 * t280 - t269 * t281) * MDP(21) + (t266 * t281 + t269 * t280) * MDP(22)) * t261; (-qJD(3) * t183 - t237 * t293 - t240 * t243 + t170) * MDP(12) + (qJD(3) * t184 + t237 * t243 - t240 * t293 - t171) * MDP(13) + ((t182 + t183) * t240 + (-t181 + t184) * t237 + (-t229 * t264 - t230 * t265) * pkin(3)) * MDP(14) + (-t181 * t183 - t182 * t184 + (t170 * t265 + t171 * t264 - t243 * t303) * pkin(3)) * MDP(15) + (-(t177 * t269 - t178 * t266) * t261 + t209 * t197 + ((-t258 * t266 - t269 * t319) * t261 + t279) * qJD(5) + t326) * MDP(21) + (-t269 * t167 - t266 * t166 + (t177 * t266 + t178 * t269) * t261 + t209 * t278 + (-(t258 * t269 - t266 * t319) * t261 - t269 * t175) * qJD(5) + t327) * MDP(22) + t321 * qJD(2) * t314 + (-t270 * t324 + t323) * t273 + t329; -t252 * MDP(13) + (-t237 ^ 2 - t240 ^ 2) * MDP(14) + (t181 * t240 + t182 * t237 + t244) * MDP(15) + (-t274 - t313) * MDP(21) + (t168 + t312) * MDP(22) + ((t264 * t302 + t265 * t303 + t240) * MDP(12) + (-t237 + t291) * MDP(13)) * qJD(3); (t279 * t320 + t326) * MDP(21) + ((-t176 * t261 - t166) * t266 + (-t175 * t320 - t167) * t269 + t327) * MDP(22) + t329;];
tauc = t1;
