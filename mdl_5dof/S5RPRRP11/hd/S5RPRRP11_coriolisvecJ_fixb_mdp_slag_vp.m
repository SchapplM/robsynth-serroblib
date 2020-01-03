% Calculate minimal parameter regressor of Coriolis joint torque vector for
% S5RPRRP11
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
% MDP [25x1]
%   Minimal dynamic parameter vector (fixed base model)
%   see S5RPRRP11_convert_par2_MPV_fixb.m
% 
% Output:
% tauc [5x1]
%   joint torques required to compensate Coriolis and centrifugal load

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-31 18:55
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauc = S5RPRRP11_coriolisvecJ_fixb_mdp_slag_vp(qJ, qJD, pkin, MDP)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(8,1),zeros(25,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RPRRP11_coriolisvecJ_fixb_mdp_slag_vp: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RPRRP11_coriolisvecJ_fixb_mdp_slag_vp: qJD has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S5RPRRP11_coriolisvecJ_fixb_mdp_slag_vp: pkin has to be [8x1] (double)');
assert(isreal(MDP) && all(size(MDP) == [25 1]), ...
  'S5RPRRP11_coriolisvecJ_fixb_mdp_slag_vp: MDP has to be [25x1] (double)'); 

%% Symbolic Calculation
% From coriolisvec_joint_fixb_mdp_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 18:54:38
% EndTime: 2019-12-31 18:54:45
% DurationCPUTime: 3.70s
% Computational Cost: add. (2559->299), mult. (6608->390), div. (0->0), fcn. (4752->6), ass. (0->105)
t258 = sin(pkin(8));
t259 = cos(pkin(8));
t261 = sin(qJ(3));
t263 = cos(qJ(3));
t327 = -t258 * t261 + t263 * t259;
t332 = t327 * qJD(1);
t316 = pkin(6) + qJ(2);
t247 = t316 * t258;
t243 = qJD(1) * t247;
t248 = t316 * t259;
t244 = qJD(1) * t248;
t215 = -t261 * t243 + t263 * t244;
t331 = qJD(3) * t215;
t262 = cos(qJ(4));
t290 = qJD(4) * t262;
t260 = sin(qJ(4));
t292 = qJD(3) * t260;
t242 = t258 * t263 + t259 * t261;
t321 = t242 * qJD(1);
t194 = t321 * t290 + (qJD(4) + t332) * t292;
t330 = (MDP(7) * qJ(2) + MDP(6)) * (t258 ^ 2 + t259 ^ 2);
t231 = qJD(4) - t332;
t267 = t327 * qJD(2);
t323 = -t243 * t263 - t261 * t244;
t187 = qJD(1) * t267 + qJD(3) * t323;
t253 = -pkin(2) * t259 - pkin(1);
t245 = qJD(1) * t253 + qJD(2);
t196 = -pkin(3) * t332 - pkin(7) * t321 + t245;
t239 = t242 * qJD(3);
t228 = qJD(1) * t239;
t238 = t327 * qJD(3);
t265 = qJD(1) * t238;
t201 = t228 * pkin(3) - pkin(7) * t265;
t210 = qJD(3) * pkin(7) + t215;
t291 = qJD(4) * t260;
t271 = t262 * t187 + t196 * t290 + t260 * t201 - t210 * t291;
t314 = qJ(5) * t228;
t167 = qJD(5) * t231 + t271 + t314;
t178 = t196 * t262 - t210 * t260;
t288 = qJD(5) - t178;
t173 = -pkin(4) * t231 + t288;
t329 = t231 * t173 + t167;
t281 = t260 * t187 + t196 * t291 - t262 * t201 + t210 * t290;
t318 = pkin(4) * t228;
t168 = t281 - t318;
t179 = t196 * t260 + t210 * t262;
t174 = qJ(5) * t231 + t179;
t312 = t174 * t231;
t328 = -t168 + t312;
t213 = -pkin(3) * t327 - pkin(7) * t242 + t253;
t219 = -t247 * t261 + t248 * t263;
t295 = t260 * t213 + t262 * t219;
t218 = t247 * t263 + t261 * t248;
t209 = -qJD(3) * pkin(3) - t323;
t289 = t262 * qJD(3);
t220 = t260 * t321 - t289;
t222 = t262 * t321 + t292;
t180 = pkin(4) * t220 - qJ(5) * t222 + t209;
t317 = pkin(7) * t228;
t320 = t180 * t231 - t317;
t319 = t222 ^ 2;
t286 = qJD(1) * qJD(2);
t188 = t242 * t286 + t331;
t193 = -qJD(4) * t289 - t262 * t265 + t291 * t321;
t170 = pkin(4) * t194 + qJ(5) * t193 - qJD(5) * t222 + t188;
t313 = t170 * t260;
t311 = t179 * t231;
t310 = t193 * t260;
t309 = t220 * t222;
t308 = t220 * t332;
t282 = t222 * t231;
t307 = t231 * t260;
t306 = t242 * t262;
t225 = t260 * t228;
t227 = t262 * t228;
t279 = pkin(4) * t260 - qJ(5) * t262;
t299 = qJD(5) * t260 - t231 * t279 + t215;
t298 = -t260 * t194 - t220 * t290;
t297 = t307 * t332 + t227;
t211 = pkin(3) * t321 - pkin(7) * t332;
t296 = t260 * t211 + t262 * t323;
t294 = t231 * t290 + t225;
t284 = pkin(7) * t291;
t280 = pkin(4) * t262 + qJ(5) * t260;
t278 = t173 * t262 - t174 * t260;
t275 = t238 * t260 + t242 * t290;
t274 = -t238 * t262 + t242 * t291;
t273 = t209 * t231 - t317;
t272 = t180 * t222 + t281;
t197 = -t218 * qJD(3) + t267;
t212 = pkin(3) * t239 - pkin(7) * t238;
t270 = t262 * t197 + t260 * t212 + t213 * t290 - t219 * t291;
t198 = qJD(2) * t242 + qJD(3) * t219;
t246 = -pkin(3) - t280;
t189 = pkin(4) * t222 + qJ(5) * t220;
t184 = t242 * t279 + t218;
t182 = pkin(4) * t327 - t213 * t262 + t219 * t260;
t181 = -qJ(5) * t327 + t295;
t177 = t220 * t231 - t193;
t176 = -pkin(4) * t321 - t211 * t262 + t260 * t323;
t175 = qJ(5) * t321 + t296;
t172 = t279 * t238 + (qJD(4) * t280 - qJD(5) * t262) * t242 + t198;
t171 = -pkin(4) * t239 + t295 * qJD(4) + t197 * t260 - t212 * t262;
t169 = qJ(5) * t239 - qJD(5) * t327 + t270;
t1 = [(t238 * t321 + t242 * t265) * MDP(8) + (-t242 * t228 + t238 * t332 - t239 * t321 + t265 * t327) * MDP(9) + (t228 * t253 + t239 * t245) * MDP(13) + t245 * t238 * MDP(14) + (-t193 * t306 - t222 * t274) * MDP(15) + ((-t220 * t262 - t222 * t260) * t238 + (t310 - t194 * t262 + (t220 * t260 - t222 * t262) * qJD(4)) * t242) * MDP(16) + (t193 * t327 + t222 * t239 + t227 * t242 - t231 * t274) * MDP(17) + (t194 * t327 - t220 * t239 - t225 * t242 - t231 * t275) * MDP(18) + (-t228 * t327 + t231 * t239) * MDP(19) + (t281 * t327 + t178 * t239 + t198 * t220 + t218 * t194 + ((-qJD(4) * t219 + t212) * t231 + t213 * t228 + t209 * qJD(4) * t242) * t262 + ((-qJD(4) * t213 - t197) * t231 - t219 * t228 + t188 * t242 + t209 * t238) * t260) * MDP(20) + (-t179 * t239 + t188 * t306 - t218 * t193 + t198 * t222 - t209 * t274 - t228 * t295 - t231 * t270 + t271 * t327) * MDP(21) + (t168 * t327 - t171 * t231 + t172 * t220 - t173 * t239 + t180 * t275 - t182 * t228 + t184 * t194 + t242 * t313) * MDP(22) + (-t169 * t220 + t171 * t222 - t181 * t194 - t182 * t193 + t278 * t238 + (-t167 * t260 + t168 * t262 + (-t173 * t260 - t174 * t262) * qJD(4)) * t242) * MDP(23) + (-t167 * t327 + t169 * t231 - t170 * t306 - t172 * t222 + t174 * t239 + t180 * t274 + t181 * t228 + t184 * t193) * MDP(24) + (t167 * t181 + t168 * t182 + t169 * t174 + t170 * t184 + t171 * t173 + t172 * t180) * MDP(25) + 0.2e1 * t286 * t330 + (t238 * MDP(10) - t239 * MDP(11) - t198 * MDP(13) + (t253 * t332 - t197) * MDP(14)) * qJD(3); t297 * MDP(20) + t298 * MDP(23) + t294 * MDP(24) + (-t180 * MDP(25) + (-MDP(21) + MDP(24)) * t222 + (-MDP(20) - MDP(22)) * t220) * t321 + (t228 * MDP(22) + (t193 + t308) * MDP(23) + t328 * MDP(25) + (-t231 * MDP(21) - MDP(24) * t332) * t231) * t262 + (-t228 * MDP(21) + t329 * MDP(25) + MDP(23) * t282 + (-qJD(4) * MDP(20) - t231 * MDP(22)) * t231) * t260 + (t321 * MDP(13) + t332 * MDP(14) + (MDP(13) * t242 + MDP(14) * t327) * qJD(1)) * qJD(3) - qJD(1) ^ 2 * t330; (-t188 + t331) * MDP(13) + (t262 * t282 - t310) * MDP(15) + ((-t193 + t308) * t262 - t222 * t307 + t298) * MDP(16) + t294 * MDP(17) + (-t231 * t291 + t297) * MDP(18) + (-pkin(3) * t194 - t215 * t220 + (-t188 + (-pkin(7) * qJD(4) - t211) * t231) * t262 + (t231 * t323 + t273) * t260) * MDP(20) + (pkin(3) * t193 + t188 * t260 - t215 * t222 + (t284 + t296) * t231 + t273 * t262) * MDP(21) + (-t170 * t262 + t194 * t246 + (-pkin(7) * t290 + t176) * t231 - t299 * t220 + t320 * t260) * MDP(22) + (t175 * t220 - t176 * t222 + ((qJD(4) * t222 - t194) * pkin(7) + t329) * t262 + ((qJD(4) * t220 - t193) * pkin(7) - t328) * t260) * MDP(23) + (-t313 + t193 * t246 + (-t175 - t284) * t231 + t299 * t222 - t320 * t262) * MDP(24) + (t170 * t246 - t173 * t176 - t174 * t175 - t299 * t180 + (qJD(4) * t278 + t167 * t262 + t168 * t260) * pkin(7)) * MDP(25) + (-MDP(13) * t245 - MDP(17) * t222 + MDP(18) * t220 - MDP(19) * t231 - MDP(20) * t178 + MDP(21) * t179 + MDP(22) * t173 - MDP(24) * t174 + MDP(9) * t321) * t321 + ((-qJD(2) - t245) * MDP(14) - t231 * t262 * MDP(17) - MDP(8) * t321 - MDP(9) * t332) * t332; MDP(15) * t309 + (-t220 ^ 2 + t319) * MDP(16) + t177 * MDP(17) + (-t194 + t282) * MDP(18) + t228 * MDP(19) + (-t209 * t222 - t281 + t311) * MDP(20) + (t178 * t231 + t209 * t220 - t271) * MDP(21) + (-t189 * t220 - t272 + t311 + 0.2e1 * t318) * MDP(22) + (pkin(4) * t193 - qJ(5) * t194 + (t174 - t179) * t222 + (t173 - t288) * t220) * MDP(23) + (0.2e1 * t314 - t180 * t220 + t189 * t222 + (0.2e1 * qJD(5) - t178) * t231 + t271) * MDP(24) + (-pkin(4) * t168 + qJ(5) * t167 - t173 * t179 + t174 * t288 - t180 * t189) * MDP(25); (-qJD(3) * t321 + t309) * MDP(22) + t177 * MDP(23) + (-t231 ^ 2 - t319) * MDP(24) + (t272 - t312 - t318) * MDP(25);];
tauc = t1;
