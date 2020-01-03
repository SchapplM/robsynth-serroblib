% Calculate minimal parameter regressor of Coriolis joint torque vector for
% S5RRPPR7
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
%   see S5RRPPR7_convert_par2_MPV_fixb.m
% 
% Output:
% tauc [5x1]
%   joint torques required to compensate Coriolis and centrifugal load

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-31 19:37
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauc = S5RRPPR7_coriolisvecJ_fixb_mdp_slag_vp(qJ, qJD, pkin, MDP)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(8,1),zeros(23,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RRPPR7_coriolisvecJ_fixb_mdp_slag_vp: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RRPPR7_coriolisvecJ_fixb_mdp_slag_vp: qJD has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S5RRPPR7_coriolisvecJ_fixb_mdp_slag_vp: pkin has to be [8x1] (double)');
assert(isreal(MDP) && all(size(MDP) == [23 1]), ...
  'S5RRPPR7_coriolisvecJ_fixb_mdp_slag_vp: MDP has to be [23x1] (double)'); 

%% Symbolic Calculation
% From coriolisvec_joint_fixb_mdp_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 19:36:26
% EndTime: 2019-12-31 19:36:31
% DurationCPUTime: 2.30s
% Computational Cost: add. (1424->270), mult. (3679->356), div. (0->0), fcn. (2505->6), ass. (0->130)
t246 = sin(pkin(8));
t249 = sin(qJ(2));
t286 = qJD(1) * t249;
t247 = cos(pkin(8));
t251 = cos(qJ(2));
t292 = t247 * t251;
t216 = -qJD(1) * t292 + t246 * t286;
t250 = cos(qJ(5));
t248 = sin(qJ(5));
t285 = qJD(2) * t248;
t194 = -t250 * t216 + t285;
t227 = t246 * t251 + t247 * t249;
t219 = t227 * qJD(1);
t313 = qJD(5) + t219;
t314 = t194 * t313;
t196 = qJD(2) * t250 + t216 * t248;
t268 = t313 * t196;
t277 = qJD(1) * qJD(2);
t312 = -0.2e1 * t277;
t311 = MDP(4) * t249;
t310 = MDP(5) * (t249 ^ 2 - t251 ^ 2);
t309 = t248 * t313;
t301 = -qJ(3) - pkin(6);
t233 = t301 * t251;
t230 = qJD(1) * t233;
t222 = t246 * t230;
t232 = t301 * t249;
t229 = qJD(1) * t232;
t191 = t229 * t247 + t222;
t280 = -qJD(4) + t191;
t225 = qJD(2) * pkin(2) + t229;
t293 = t247 * t230;
t189 = t246 * t225 - t293;
t186 = -qJD(2) * qJ(4) - t189;
t303 = pkin(4) * t216;
t170 = -t186 - t303;
t190 = t229 * t246 - t293;
t275 = t249 * t277;
t234 = t246 * t275;
t274 = t251 * t277;
t209 = t247 * t274 - t234;
t241 = -pkin(2) * t247 - pkin(3);
t237 = -pkin(7) + t241;
t307 = t237 * t209 + (t170 - t190 + t303) * t313;
t215 = t219 ^ 2;
t305 = pkin(3) + pkin(7);
t218 = t227 * qJD(2);
t208 = qJD(1) * t218;
t304 = pkin(3) * t208;
t302 = pkin(4) * t219;
t226 = t246 * t249 - t292;
t276 = -pkin(2) * t251 - pkin(1);
t262 = -qJ(4) * t227 + t276;
t172 = t226 * t305 + t262;
t300 = t172 * t209;
t273 = qJD(2) * t301;
t213 = qJD(3) * t251 + t249 * t273;
t203 = t213 * qJD(1);
t214 = -qJD(3) * t249 + t251 * t273;
t204 = t214 * qJD(1);
t173 = t203 * t246 - t247 * t204;
t192 = -t247 * t232 - t233 * t246;
t299 = t173 * t192;
t283 = qJD(5) * t248;
t282 = qJD(5) * t250;
t288 = t248 * t208 + t216 * t282;
t177 = -qJD(2) * t283 + t288;
t298 = t177 * t250;
t297 = t194 * t216;
t296 = t196 * t216;
t295 = t209 * t248;
t294 = t226 * t248;
t252 = qJD(2) ^ 2;
t291 = t249 * t252;
t202 = t250 * t209;
t290 = t251 * t252;
t253 = qJD(1) ^ 2;
t289 = t251 * t253;
t174 = t247 * t203 + t246 * t204;
t284 = qJD(2) * t249;
t281 = t219 * MDP(14);
t279 = t302 - t280;
t243 = pkin(2) * t284;
t272 = pkin(1) * t312;
t238 = pkin(2) * t275;
t271 = -qJ(4) * t209 + t238;
t270 = pkin(2) * t286 + qJ(4) * t216;
t181 = t213 * t246 - t247 * t214;
t188 = t225 * t247 + t222;
t269 = t250 * t313;
t267 = MDP(22) * t313;
t266 = t276 * qJD(1);
t265 = qJD(4) - t188;
t171 = -qJD(2) * qJD(4) - t174;
t231 = qJD(3) + t266;
t255 = -qJ(4) * t219 + t231;
t164 = t216 * t305 + t255;
t165 = -qJD(2) * t305 + t265 + t302;
t157 = t164 * t250 + t165 * t248;
t264 = t164 * t248 - t165 * t250;
t182 = t213 * t247 + t214 * t246;
t193 = t232 * t246 - t233 * t247;
t179 = pkin(3) * t216 + t255;
t261 = t179 * t219 + t173;
t260 = t218 * t248 + t226 * t282;
t259 = -qJD(4) * t219 + t271;
t221 = qJD(2) * t292 - t246 * t284;
t258 = -qJ(4) * t221 - qJD(4) * t227 + t243;
t160 = -pkin(4) * t208 - t171;
t257 = t160 + (-qJD(5) * t237 + t305 * t219 + t270) * t313;
t183 = pkin(4) * t227 + t192;
t256 = t160 * t226 + t170 * t218 - t183 * t209;
t254 = t173 * t227 + t181 * t219 - t182 * t216 + t192 * t209 - t193 * t208;
t239 = pkin(2) * t246 + qJ(4);
t211 = qJD(2) * t216;
t201 = t250 * t208;
t187 = pkin(3) * t226 + t262;
t185 = -qJD(2) * pkin(3) + t265;
t184 = -pkin(4) * t226 + t193;
t180 = pkin(3) * t219 + t270;
t178 = qJD(5) * t196 - t201;
t169 = pkin(3) * t218 + t258;
t168 = -pkin(4) * t218 + t182;
t167 = pkin(4) * t221 + t181;
t163 = t259 + t304;
t162 = pkin(4) * t209 + t173;
t161 = t250 * t162;
t159 = t218 * t305 + t258;
t158 = t208 * t305 + t259;
t1 = [0.2e1 * t274 * t311 + t310 * t312 + MDP(6) * t290 - MDP(7) * t291 + (-pkin(6) * t290 + t249 * t272) * MDP(9) + (pkin(6) * t291 + t251 * t272) * MDP(10) + (-t174 * t226 - t188 * t221 - t189 * t218 + t254) * MDP(11) + (t299 + t174 * t193 - t188 * t181 + t189 * t182 + (t231 + t266) * t243) * MDP(12) + (t171 * t226 + t185 * t221 + t186 * t218 + t254) * MDP(13) + (qJD(2) * t181 - t163 * t226 - t169 * t216 - t179 * t218 - t187 * t208) * MDP(14) + (qJD(2) * t182 - t163 * t227 - t169 * t219 - t179 * t221 - t187 * t209) * MDP(15) + (t163 * t187 + t169 * t179 - t171 * t193 + t181 * t185 - t182 * t186 + t299) * MDP(16) + (t177 * t294 + t196 * t260) * MDP(17) + ((-t194 * t248 + t196 * t250) * t218 + (t298 - t178 * t248 + (-t194 * t250 - t196 * t248) * qJD(5)) * t226) * MDP(18) + (t177 * t227 + t196 * t221 + t209 * t294 + t260 * t313) * MDP(19) + (t226 * t202 - t178 * t227 - t194 * t221 + (t218 * t250 - t226 * t283) * t313) * MDP(20) + (t209 * t227 + t221 * t313) * MDP(21) + (-t264 * t221 + t161 * t227 + t168 * t194 + t184 * t178 + (-t158 * t227 - t159 * t313 - t300) * t248 + (t167 * t313 - t256) * t250 + ((-t172 * t250 - t183 * t248) * t313 - t157 * t227 + t170 * t294) * qJD(5)) * MDP(22) + (-t157 * t221 + t168 * t196 + t184 * t177 + (-(qJD(5) * t183 + t159) * t313 - t300 - (qJD(5) * t165 + t158) * t227 + t170 * qJD(5) * t226) * t250 + (-(-qJD(5) * t172 + t167) * t313 - (-qJD(5) * t164 + t162) * t227 + t256) * t248) * MDP(23); -t289 * t311 + t253 * t310 + ((t189 - t190) * t219 + (-t188 + t191) * t216 + (-t208 * t246 - t209 * t247) * pkin(2)) * MDP(11) + (t188 * t190 - t189 * t191 + (-t173 * t247 + t174 * t246 - t231 * t286) * pkin(2)) * MDP(12) + (-t208 * t239 + t209 * t241 + (-t186 - t190) * t219 + (t185 + t280) * t216) * MDP(13) + (-qJD(2) * t190 + t180 * t216 + t261) * MDP(14) + (-t179 * t216 + t180 * t219 + (0.2e1 * qJD(4) - t191) * qJD(2) + t174) * MDP(15) + (-t171 * t239 + t173 * t241 - t179 * t180 - t185 * t190 + t186 * t280) * MDP(16) + (-t248 * t268 + t298) * MDP(17) + ((-t178 - t268) * t250 + (-t177 + t314) * t248) * MDP(18) + (-t309 * t313 + t202 + t296) * MDP(19) + (-t269 * t313 - t295 - t297) * MDP(20) + t313 * t216 * MDP(21) + (t239 * t178 + t279 * t194 - t216 * t264 + t257 * t248 + t250 * t307) * MDP(22) + (-t157 * t216 + t239 * t177 + t279 * t196 - t248 * t307 + t257 * t250) * MDP(23) + (MDP(9) * t249 * t253 + MDP(10) * t289) * pkin(1); (t188 * t219 + t189 * t216 + t238) * MDP(12) + (t211 + t234) * MDP(15) + (t304 - t186 * t216 + (-qJD(4) - t185) * t219 + t271) * MDP(16) + (-t295 + t297) * MDP(22) + (-t202 + t296) * MDP(23) + (MDP(23) * t309 - t250 * t267) * t313 + (-t281 + (-MDP(14) * t227 - MDP(15) * t292) * qJD(1)) * qJD(2) + (MDP(11) + MDP(13)) * (-t216 ^ 2 - t215); (t209 + t211) * MDP(13) - t216 * t281 + (-t215 - t252) * MDP(15) + (qJD(2) * t186 + t261) * MDP(16) + (-qJD(2) * t194 + t202) * MDP(22) + (-qJD(2) * t196 - t295) * MDP(23) + (-MDP(23) * t269 - t248 * t267) * t313; t196 * t194 * MDP(17) + (-t194 ^ 2 + t196 ^ 2) * MDP(18) + (t288 + t314) * MDP(19) + (t201 + t268) * MDP(20) + t209 * MDP(21) + (t157 * t313 - t158 * t248 - t170 * t196 + t161) * MDP(22) + (-t158 * t250 - t162 * t248 + t170 * t194 - t264 * t313) * MDP(23) + (-MDP(19) * t285 - MDP(20) * t196 - MDP(22) * t157 + MDP(23) * t264) * qJD(5);];
tauc = t1;
