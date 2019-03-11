% Calculate minimal parameter regressor of Coriolis joint torque vector for
% S6PPPRRR1
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% qJD [6x1]
%   Generalized joint velocities
% pkin [14x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,alpha2,alpha3,alpha4,d4,d5,d6,theta1,theta2,theta3]';
% MDP [20x1]
%   Minimal dynamic parameter vector (fixed base model)
%   see S6PPPRRR1_convert_par2_MPV_fixb.m
% 
% Output:
% tauc [6x1]
%   joint torques required to compensate Coriolis and centrifugal load

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-08 18:41
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauc = S6PPPRRR1_coriolisvecJ_fixb_mdp_slag_vp(qJ, qJD, pkin, MDP)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(14,1),zeros(20,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6PPPRRR1_coriolisvecJ_fixb_mdp_slag_vp: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6PPPRRR1_coriolisvecJ_fixb_mdp_slag_vp: qJD has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [14 1]), ...
  'S6PPPRRR1_coriolisvecJ_fixb_mdp_slag_vp: pkin has to be [14x1] (double)');
assert(isreal(MDP) && all(size(MDP) == [20 1]), ...
  'S6PPPRRR1_coriolisvecJ_fixb_mdp_slag_vp: MDP has to be [20x1] (double)'); 

%% Symbolic Calculation
% From coriolisvec_joint_fixb_mdp_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-08 18:40:38
% EndTime: 2019-03-08 18:40:44
% DurationCPUTime: 2.31s
% Computational Cost: add. (1996->245), mult. (5409->407), div. (0->0), fcn. (5240->16), ass. (0->132)
t232 = cos(pkin(6));
t217 = qJD(1) * t232 + qJD(2);
t223 = sin(pkin(14));
t226 = sin(pkin(7));
t227 = sin(pkin(6));
t229 = cos(pkin(13));
t231 = cos(pkin(7));
t290 = t229 * t231;
t224 = sin(pkin(13));
t228 = cos(pkin(14));
t296 = t224 * t228;
t193 = t217 * t223 * t226 + (t223 * t290 + t296) * t227 * qJD(1);
t235 = sin(qJ(4));
t238 = cos(qJ(4));
t243 = (-t223 * t224 + t228 * t290) * t227;
t293 = t226 * t228;
t192 = qJD(1) * t243 + t217 * t293;
t267 = t226 * t227 * t229;
t205 = -qJD(1) * t267 + t217 * t231 + qJD(3);
t225 = sin(pkin(8));
t230 = cos(pkin(8));
t252 = t192 * t230 + t205 * t225;
t242 = t235 * t193 - t252 * t238;
t237 = cos(qJ(5));
t284 = MDP(13) * t237;
t234 = sin(qJ(5));
t313 = MDP(7) * t234;
t221 = t234 ^ 2;
t312 = MDP(8) * (-t237 ^ 2 + t221);
t278 = qJD(4) * t238;
t281 = qJD(4) * t235;
t174 = t193 * t278 + t252 * t281;
t178 = t238 * t193 + t235 * t252;
t310 = qJD(4) * t178 - t174;
t291 = t228 * t230;
t294 = t225 * t238;
t309 = t226 * (-t223 * t235 + t238 * t291) + t231 * t294;
t292 = t226 * t232;
t197 = t227 * t296 + (t227 * t290 + t292) * t223;
t196 = t228 * t292 + t243;
t207 = t231 * t232 - t267;
t251 = t196 * t230 + t207 * t225;
t308 = -t197 * t235 + t238 * t251;
t307 = qJD(4) * pkin(4);
t176 = qJD(4) * pkin(10) + t178;
t185 = -t192 * t225 + t205 * t230;
t167 = t176 * t237 + t185 * t234;
t173 = t242 * qJD(4);
t161 = t167 * qJD(5) - t234 * t173;
t233 = sin(qJ(6));
t306 = t161 * t233;
t236 = cos(qJ(6));
t305 = t161 * t236;
t274 = qJD(6) * t233;
t264 = t234 * t274;
t270 = qJD(4) * qJD(5);
t262 = t237 * t270;
t269 = qJD(5) * qJD(6);
t286 = (t262 + t269) * t236;
t203 = -qJD(4) * t264 + t286;
t302 = t203 * t233;
t282 = qJD(4) * t234;
t265 = t233 * t282;
t271 = t236 * qJD(5);
t210 = t265 - t271;
t279 = qJD(4) * t237;
t218 = -qJD(6) + t279;
t300 = t210 * t218;
t277 = qJD(5) * t233;
t280 = qJD(4) * t236;
t212 = t234 * t280 + t277;
t299 = t212 * t218;
t298 = t218 * t233;
t297 = t218 * t236;
t295 = t225 * t235;
t289 = t233 * t237;
t288 = t236 * t237;
t259 = pkin(5) * t234 - pkin(11) * t237;
t214 = t259 * qJD(5);
t287 = t178 - t214;
t276 = qJD(5) * t234;
t275 = qJD(5) * t237;
t273 = qJD(6) * t236;
t272 = qJD(6) * t238;
t266 = t225 * t278;
t263 = t218 * t273;
t261 = MDP(18) * t276;
t165 = qJD(5) * pkin(11) + t167;
t215 = -pkin(5) * t237 - pkin(11) * t234 - pkin(4);
t172 = t215 * qJD(4) + t242;
t159 = t165 * t236 + t172 * t233;
t258 = t165 * t233 - t172 * t236;
t182 = t197 * t238 + t235 * t251;
t186 = -t196 * t225 + t207 * t230;
t169 = t182 * t237 + t186 * t234;
t257 = t169 * t236 - t233 * t308;
t256 = -t169 * t233 - t236 * t308;
t255 = t176 * t234 - t185 * t237;
t168 = t182 * t234 - t186 * t237;
t200 = t231 * t295 + (t223 * t238 + t235 * t291) * t226;
t206 = -t225 * t293 + t230 * t231;
t188 = t200 * t237 + t206 * t234;
t254 = t188 * t236 - t233 * t309;
t253 = -t188 * t233 - t236 * t309;
t187 = t200 * t234 - t206 * t237;
t250 = qJD(4) * t221 - t218 * t237;
t239 = qJD(5) ^ 2;
t249 = pkin(10) * t239 - t310;
t175 = t242 - t307;
t248 = qJD(5) * (t175 - t242 - t307);
t209 = t230 * t234 + t237 * t295;
t208 = -t230 * t237 + t234 * t295;
t246 = -MDP(12) * t237 + MDP(13) * t234 - MDP(5);
t160 = -qJD(5) * t255 - t237 * t173;
t164 = -qJD(5) * pkin(5) + t255;
t241 = qJD(5) * t164 + qJD(6) * t172 + t218 * t242 + t160;
t240 = qJD(4) ^ 2;
t213 = t259 * qJD(4);
t204 = t233 * t269 + (t233 * t275 + t234 * t273) * qJD(4);
t202 = qJD(5) * t209 + t234 * t266;
t201 = -qJD(5) * t208 + t237 * t266;
t195 = t200 * qJD(4);
t194 = t309 * qJD(4);
t184 = qJD(5) * t188 + t234 * t194;
t183 = -qJD(5) * t187 + t237 * t194;
t180 = t182 * qJD(4);
t179 = t308 * qJD(4);
t171 = qJD(4) * t214 + t174;
t170 = t236 * t171;
t163 = -qJD(5) * t168 + t179 * t237;
t162 = qJD(5) * t169 + t179 * t234;
t1 = [(-(-qJD(6) * t257 - t163 * t233 + t180 * t236) * t218 + t162 * t210 + t168 * t204) * MDP(19) + ((qJD(6) * t256 + t163 * t236 + t180 * t233) * t218 + t162 * t212 + t168 * t203) * MDP(20) + (-MDP(12) * t162 - MDP(13) * t163) * qJD(5) + (-t179 * MDP(6) + t246 * t180 + (-t308 * t284 + (-MDP(12) * t308 + MDP(19) * t256 - MDP(20) * t257) * t234) * qJD(5)) * qJD(4); (-(-qJD(6) * t254 - t233 * t183 + t236 * t195) * t218 + t184 * t210 + t187 * t204) * MDP(19) + ((qJD(6) * t253 + t236 * t183 + t233 * t195) * t218 + t184 * t212 + t187 * t203) * MDP(20) + (-MDP(12) * t184 - MDP(13) * t183) * qJD(5) + (-t194 * MDP(6) + t246 * t195 + (-t309 * t284 + (-MDP(12) * t309 + MDP(19) * t253 - MDP(20) * t254) * t234) * qJD(5)) * qJD(4); (-(-t233 * t201 - t209 * t273) * t218 + t202 * t210 + t208 * t204) * MDP(19) + ((t236 * t201 - t209 * t274) * t218 + t202 * t212 + t208 * t203) * MDP(20) + ((-(t233 * t272 + t235 * t280) * MDP(19) + (t233 * t281 - t236 * t272) * MDP(20)) * t218 + (-t238 * MDP(6) + t235 * t246) * t240) * t225 + (-t202 * MDP(12) - t201 * MDP(13) + (-t284 * t294 + (-MDP(12) * t294 + (-t209 * t233 - t236 * t294) * MDP(19) - (t209 * t236 - t233 * t294) * MDP(20)) * t234) * qJD(4)) * qJD(5); t310 * MDP(5) + 0.2e1 * t262 * t313 - 0.2e1 * t270 * t312 + (t234 * t248 - t237 * t249) * MDP(12) + (t234 * t249 + t237 * t248) * MDP(13) + (t203 * t236 * t234 + (t237 * t271 - t264) * t212) * MDP(14) + ((-t210 * t236 - t212 * t233) * t275 + (-t302 - t204 * t236 + (t210 * t233 - t212 * t236) * qJD(6)) * t234) * MDP(15) + (t218 * t264 - t237 * t203 + (t234 * t212 + t236 * t250) * qJD(5)) * MDP(16) + (t234 * t263 + t204 * t237 + (-t234 * t210 - t233 * t250) * qJD(5)) * MDP(17) + (-t218 - t279) * t261 + ((t215 * t274 + t236 * t287) * t218 + (t165 * t273 - t170 + (qJD(5) * t210 + t263) * pkin(10) + t241 * t233) * t237 + (t164 * t273 + pkin(10) * t204 + t306 + t242 * t210 + (-pkin(10) * t298 + (-pkin(10) * t289 + t215 * t236) * qJD(4) - t258) * qJD(5)) * t234) * MDP(19) + ((t215 * t273 - t233 * t287) * t218 + (qJD(5) * pkin(10) * t212 + (t171 + (-pkin(10) * t218 - t165) * qJD(6)) * t233 + t241 * t236) * t237 + (-t164 * t274 + pkin(10) * t203 + t305 + t242 * t212 + (-pkin(10) * t297 - (pkin(10) * t288 + t215 * t233) * qJD(4) - t159) * qJD(5)) * t234) * MDP(20) + (-MDP(10) * t234 + MDP(9) * t237) * t239; (-t212 * t297 + t302) * MDP(14) + ((t203 + t300) * t236 + (-t204 + t299) * t233) * MDP(15) + (-t263 + (t218 * t288 + (-t212 + t277) * t234) * qJD(4)) * MDP(16) + (t218 * t274 + (-t218 * t289 + (t210 + t271) * t234) * qJD(4)) * MDP(17) + t218 * MDP(18) * t282 + (-pkin(5) * t204 - t305 + (t213 * t236 + t233 * t255) * t218 - t167 * t210 + (pkin(11) * t297 + t164 * t233) * qJD(6) + (t258 * t234 + (-pkin(11) * t276 - t164 * t237) * t233) * qJD(4)) * MDP(19) + (-pkin(5) * t203 + t306 - (t213 * t233 - t236 * t255) * t218 - t167 * t212 + (-pkin(11) * t298 + t164 * t236) * qJD(6) + (-t164 * t288 + (-pkin(11) * t271 + t159) * t234) * qJD(4)) * MDP(20) + (-t237 * t313 + t312) * t240 + (MDP(12) * t234 + t284) * (-qJD(4) * t175 + t173); t210 * t212 * MDP(14) + (-t210 ^ 2 + t212 ^ 2) * MDP(15) + (t286 - t300) * MDP(16) + (-t233 * t262 - t299) * MDP(17) + qJD(4) * t261 + (-t159 * t218 - t233 * t160 - t164 * t212 + t170) * MDP(19) + (-t236 * t160 + t164 * t210 - t233 * t171 + t218 * t258) * MDP(20) + (-MDP(16) * t265 - MDP(17) * t212 - MDP(19) * t159 + MDP(20) * t258) * qJD(6);];
tauc  = t1;
