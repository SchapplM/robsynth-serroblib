% Calculate minimal parameter regressor of Coriolis joint torque vector for
% S6PPRPRR1
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% qJD [6x1]
%   Generalized joint velocities
% pkin [13x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,alpha2,alpha3,d3,d5,d6,theta1,theta2,theta4]';
% MDP [20x1]
%   Minimal dynamic parameter vector (fixed base model)
%   see S6PPRPRR1_convert_par2_MPV_fixb.m
% 
% Output:
% tauc [6x1]
%   joint torques required to compensate Coriolis and centrifugal load

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-08 18:44
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauc = S6PPRPRR1_coriolisvecJ_fixb_mdp_slag_vp(qJ, qJD, pkin, MDP)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(13,1),zeros(20,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6PPRPRR1_coriolisvecJ_fixb_mdp_slag_vp: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6PPRPRR1_coriolisvecJ_fixb_mdp_slag_vp: qJD has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [13 1]), ...
  'S6PPRPRR1_coriolisvecJ_fixb_mdp_slag_vp: pkin has to be [13x1] (double)');
assert(isreal(MDP) && all(size(MDP) == [20 1]), ...
  'S6PPRPRR1_coriolisvecJ_fixb_mdp_slag_vp: MDP has to be [20x1] (double)'); 

%% Symbolic Calculation
% From coriolisvec_joint_fixb_mdp_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-08 18:43:41
% EndTime: 2019-03-08 18:43:46
% DurationCPUTime: 2.24s
% Computational Cost: add. (1550->244), mult. (4350->384), div. (0->0), fcn. (3876->14), ass. (0->137)
t235 = sin(qJ(5));
t238 = cos(qJ(5));
t290 = MDP(13) * t238;
t314 = MDP(12) * t235 + t290;
t233 = cos(pkin(6));
t228 = sin(pkin(7));
t236 = sin(qJ(3));
t301 = t228 * t236;
t227 = sin(pkin(12));
t229 = sin(pkin(6));
t239 = cos(qJ(3));
t231 = cos(pkin(12));
t232 = cos(pkin(7));
t299 = t231 * t232;
t312 = (t227 * t239 + t236 * t299) * t229;
t197 = t233 * t301 + t312;
t218 = qJD(1) * t233 + qJD(2);
t193 = qJD(1) * t312 + t218 * t301;
t224 = t235 ^ 2;
t313 = MDP(8) * (-t238 ^ 2 + t224);
t226 = sin(pkin(13));
t230 = cos(pkin(13));
t205 = (t226 * t236 - t230 * t239) * t228;
t237 = cos(qJ(6));
t234 = sin(qJ(6));
t285 = qJD(5) * t234;
t288 = qJD(3) * t235;
t214 = t237 * t288 + t285;
t311 = MDP(14) * t214;
t310 = pkin(3) * t230;
t190 = t230 * t193;
t300 = t228 * t239;
t211 = t218 * t300;
t289 = qJD(1) * t229;
t274 = t231 * t289;
t302 = t227 * t236;
t192 = t232 * t239 * t274 - t289 * t302 + t211;
t191 = qJD(3) * pkin(3) + t192;
t177 = t226 * t191 + t190;
t175 = qJD(3) * pkin(9) + t177;
t202 = t218 * t232 - t228 * t274 + qJD(4);
t166 = t175 * t238 + t202 * t235;
t252 = t239 * t299 - t302;
t247 = t252 * t229;
t188 = (qJD(1) * t247 + t211) * qJD(3);
t242 = qJD(3) * t193;
t173 = t230 * t188 - t226 * t242;
t160 = qJD(5) * t166 + t173 * t235;
t309 = t160 * t234;
t308 = t160 * t237;
t282 = qJD(6) * t234;
t272 = t235 * t282;
t278 = qJD(3) * qJD(5);
t269 = t238 * t278;
t277 = qJD(5) * qJD(6);
t293 = (t269 + t277) * t237;
t200 = -qJD(3) * t272 + t293;
t307 = t200 * t234;
t281 = qJD(6) * t237;
t283 = qJD(5) * t238;
t201 = t234 * t277 + (t234 * t283 + t235 * t281) * qJD(3);
t306 = t201 * t238;
t273 = t234 * t288;
t279 = t237 * qJD(5);
t212 = t273 - t279;
t287 = qJD(3) * t238;
t219 = -qJD(6) + t287;
t305 = t212 * t219;
t304 = t212 * t235;
t303 = t214 * t219;
t189 = t226 * t193;
t298 = t234 * t219;
t297 = t234 * t238;
t296 = t237 * t219;
t295 = t237 * t238;
t178 = t192 * t226 + t190;
t263 = pkin(5) * t235 - pkin(10) * t238;
t216 = t263 * qJD(5);
t294 = t178 - t216;
t220 = pkin(3) * t226 + pkin(9);
t286 = qJD(5) * t220;
t284 = qJD(5) * t235;
t259 = t175 * t235 - t202 * t238;
t163 = -qJD(5) * pkin(5) + t259;
t280 = t163 * qJD(6);
t271 = t219 * t281;
t270 = MDP(18) * t284;
t164 = qJD(5) * pkin(10) + t166;
t268 = t219 * t220 + t164;
t172 = t188 * t226 + t230 * t242;
t176 = t191 * t230 - t189;
t267 = -t200 * t238 + t214 * t284;
t265 = t219 * t272;
t264 = t235 * t271;
t254 = -pkin(5) * t238 - pkin(10) * t235 - pkin(4);
t169 = qJD(3) * t254 - t176;
t158 = t164 * t237 + t169 * t234;
t262 = t164 * t234 - t169 * t237;
t196 = t233 * t300 + t247;
t183 = t196 * t226 + t197 * t230;
t208 = -t228 * t229 * t231 + t232 * t233;
t171 = t183 * t238 + t208 * t235;
t182 = -t196 * t230 + t197 * t226;
t261 = t171 * t237 + t182 * t234;
t260 = -t171 * t234 + t182 * t237;
t170 = t183 * t235 - t208 * t238;
t206 = (t226 * t239 + t230 * t236) * t228;
t199 = t206 * t238 + t232 * t235;
t258 = t199 * t237 + t205 * t234;
t257 = -t199 * t234 + t205 * t237;
t198 = t206 * t235 - t232 * t238;
t255 = qJD(3) * t224 - t219 * t238;
t253 = -MDP(12) * t238 + MDP(13) * t235;
t240 = qJD(5) ^ 2;
t250 = qJD(3) * t178 - t220 * t240 - t172;
t174 = -qJD(3) * pkin(4) - t176;
t179 = t192 * t230 - t189;
t249 = qJD(5) * (qJD(3) * (-pkin(4) - t310) + t174 + t179);
t248 = t255 * t234;
t159 = -qJD(5) * t259 + t173 * t238;
t244 = qJD(5) * t163 + qJD(6) * t169 - t179 * t219 + t159;
t241 = qJD(3) ^ 2;
t215 = t263 * qJD(3);
t209 = t254 - t310;
t204 = qJD(3) * t205;
t203 = qJD(3) * t206;
t195 = t197 * qJD(3);
t194 = t196 * qJD(3);
t185 = qJD(5) * t199 - t204 * t235;
t184 = -qJD(5) * t198 - t204 * t238;
t181 = t194 * t230 - t195 * t226;
t180 = t194 * t226 + t195 * t230;
t168 = qJD(3) * t216 + t172;
t167 = t237 * t168;
t162 = qJD(5) * t171 + t181 * t235;
t161 = -qJD(5) * t170 + t181 * t238;
t1 = [(t172 * t182 + t173 * t183 - t176 * t180 + t177 * t181) * MDP(6) + (-(-qJD(6) * t261 - t161 * t234 + t180 * t237) * t219 + t162 * t212 + t170 * t201) * MDP(19) + ((qJD(6) * t260 + t161 * t237 + t180 * t234) * t219 + t162 * t214 + t170 * t200) * MDP(20) + (-MDP(12) * t162 - MDP(13) * t161) * qJD(5) + (-t195 * MDP(4) - t194 * MDP(5) + t253 * t180 + (t182 * t290 + (t182 * MDP(12) + MDP(19) * t260 - MDP(20) * t261) * t235) * qJD(5)) * qJD(3); (t172 * t205 + t173 * t206 - t176 * t203 - t177 * t204) * MDP(6) + (-(-qJD(6) * t258 - t184 * t234 + t203 * t237) * t219 + t185 * t212 + t198 * t201) * MDP(19) + ((qJD(6) * t257 + t184 * t237 + t203 * t234) * t219 + t185 * t214 + t198 * t200) * MDP(20) + (-MDP(4) * t236 - MDP(5) * t239) * t241 * t228 + (-MDP(12) * t185 - MDP(13) * t184) * qJD(5) + (t253 * t203 + (t205 * t290 + (t205 * MDP(12) + MDP(19) * t257 - MDP(20) * t258) * t235) * qJD(5)) * qJD(3); (t176 * t178 - t177 * t179 + (-t172 * t230 + t173 * t226) * pkin(3)) * MDP(6) - 0.2e1 * t278 * t313 - t272 * t311 + (-t212 * t237 - t214 * t234) * t283 * MDP(15) + (t255 * t279 + t265 + t267) * MDP(16) + (t264 + t306 + (-t248 - t304) * qJD(5)) * MDP(17) - t287 * t270 + (-t270 + (t209 * t282 + t237 * t294) * MDP(19) + (t209 * t281 - t234 * t294) * MDP(20)) * t219 + (t240 * MDP(9) + t250 * MDP(12) + t249 * MDP(13) + t279 * t311 + (t212 * t286 + t234 * t244 + t268 * t281 - t167) * MDP(19) + (t214 * t286 + (-qJD(6) * t268 + t168) * t234 + t244 * t237) * MDP(20)) * t238 + (-t252 * t289 + t192 - t211) * qJD(3) * MDP(5) + (0.2e1 * MDP(7) * t269 - t240 * MDP(10) + t249 * MDP(12) - t250 * MDP(13) + t200 * t237 * MDP(14) + (-t307 - t201 * t237 + (t212 * t234 - t214 * t237) * qJD(6)) * MDP(15) + (t237 * t280 + t309 - t179 * t212 + t220 * t201 + (-t220 * t298 + (t209 * t237 - t220 * t297) * qJD(3) - t262) * qJD(5)) * MDP(19) + (-t234 * t280 + t308 - t179 * t214 + t220 * t200 + (-t220 * t296 - (t209 * t234 + t220 * t295) * qJD(3) - t158) * qJD(5)) * MDP(20)) * t235; (t264 - t306) * MDP(19) + (-t265 + t267) * MDP(20) - t314 * t240 + ((-t248 + t304) * MDP(19) - t255 * MDP(20) * t237) * qJD(5); (-t214 * t296 + t307) * MDP(14) + ((t200 + t305) * t237 + (-t201 + t303) * t234) * MDP(15) + (-t271 + (t219 * t295 + (-t214 + t285) * t235) * qJD(3)) * MDP(16) + (t219 * t282 + (-t219 * t297 + (t212 + t279) * t235) * qJD(3)) * MDP(17) + t219 * MDP(18) * t288 + (-pkin(5) * t201 - t308 + (t215 * t237 + t234 * t259) * t219 - t166 * t212 + (pkin(10) * t296 + t163 * t234) * qJD(6) + (t262 * t235 + (-pkin(10) * t284 - t163 * t238) * t234) * qJD(3)) * MDP(19) + (-pkin(5) * t200 + t309 - (t215 * t234 - t237 * t259) * t219 - t166 * t214 + (-pkin(10) * t298 + t163 * t237) * qJD(6) + (-t163 * t295 + (-pkin(10) * t279 + t158) * t235) * qJD(3)) * MDP(20) + (-t235 * t238 * MDP(7) + t313) * t241 + t314 * (-qJD(3) * t174 - t173); t212 * t311 + (-t212 ^ 2 + t214 ^ 2) * MDP(15) + (t293 - t305) * MDP(16) + (-t234 * t269 - t303) * MDP(17) + qJD(3) * t270 + (-t158 * t219 - t159 * t234 - t163 * t214 + t167) * MDP(19) + (-t159 * t237 + t163 * t212 - t168 * t234 + t219 * t262) * MDP(20) + (-MDP(16) * t273 - MDP(17) * t214 - MDP(19) * t158 + MDP(20) * t262) * qJD(6);];
tauc  = t1;
