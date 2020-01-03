% Calculate minimal parameter regressor of Coriolis joint torque vector for
% S5RRPRP6
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
% MDP [21x1]
%   Minimal dynamic parameter vector (fixed base model)
%   see S5RRPRP6_convert_par2_MPV_fixb.m
% 
% Output:
% tauc [5x1]
%   joint torques required to compensate Coriolis and centrifugal load

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-31 19:59
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauc = S5RRPRP6_coriolisvecJ_fixb_mdp_slag_vp(qJ, qJD, pkin, MDP)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(8,1),zeros(21,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RRPRP6_coriolisvecJ_fixb_mdp_slag_vp: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RRPRP6_coriolisvecJ_fixb_mdp_slag_vp: qJD has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S5RRPRP6_coriolisvecJ_fixb_mdp_slag_vp: pkin has to be [8x1] (double)');
assert(isreal(MDP) && all(size(MDP) == [21 1]), ...
  'S5RRPRP6_coriolisvecJ_fixb_mdp_slag_vp: MDP has to be [21x1] (double)'); 

%% Symbolic Calculation
% From coriolisvec_joint_fixb_mdp_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 19:58:32
% EndTime: 2019-12-31 19:58:37
% DurationCPUTime: 1.97s
% Computational Cost: add. (2102->270), mult. (5372->378), div. (0->0), fcn. (3693->6), ass. (0->126)
t264 = sin(pkin(8));
t265 = cos(pkin(8));
t267 = sin(qJ(2));
t269 = cos(qJ(2));
t247 = t264 * t269 + t265 * t267;
t237 = t247 * qJD(1);
t268 = cos(qJ(4));
t292 = qJD(1) * qJD(2);
t286 = t269 * t292;
t287 = t267 * t292;
t278 = -t264 * t287 + t265 * t286;
t293 = t268 * qJD(2);
t266 = sin(qJ(4));
t296 = qJD(4) * t266;
t187 = -qJD(4) * t293 + t237 * t296 - t268 * t278;
t220 = qJD(2) * t266 + t237 * t268;
t236 = t247 * qJD(2);
t230 = qJD(1) * t236;
t294 = t267 * qJD(1);
t309 = t265 * t269;
t235 = qJD(1) * t309 - t264 * t294;
t289 = -pkin(2) * t269 - pkin(1);
t281 = t289 * qJD(1);
t251 = qJD(3) + t281;
t191 = -pkin(3) * t235 - pkin(7) * t237 + t251;
t319 = -qJ(3) - pkin(6);
t252 = t319 * t267;
t249 = qJD(1) * t252;
t318 = qJD(2) * pkin(2);
t243 = t249 + t318;
t253 = t319 * t269;
t250 = qJD(1) * t253;
t310 = t265 * t250;
t210 = t264 * t243 - t310;
t205 = qJD(2) * pkin(7) + t210;
t177 = t191 * t266 + t205 * t268;
t256 = pkin(2) * t287;
t190 = t230 * pkin(3) - pkin(7) * t278 + t256;
t183 = t268 * t190;
t285 = qJD(2) * t319;
t233 = qJD(3) * t269 + t267 * t285;
t227 = t233 * qJD(1);
t234 = -qJD(3) * t267 + t269 * t285;
t228 = t234 * qJD(1);
t186 = t227 * t265 + t228 * t264;
t272 = -qJD(4) * t177 - t186 * t266 + t183;
t165 = pkin(4) * t230 + qJ(5) * t187 - qJD(5) * t220 + t272;
t218 = t237 * t266 - t293;
t172 = -qJ(5) * t218 + t177;
t232 = qJD(4) - t235;
t327 = t172 * t232 + t165;
t188 = qJD(4) * t220 + t266 * t278;
t295 = qJD(4) * t268;
t275 = t268 * t186 + t266 * t190 + t191 * t295 - t205 * t296;
t166 = -qJ(5) * t188 - qJD(5) * t218 + t275;
t176 = t268 * t191 - t205 * t266;
t171 = -qJ(5) * t220 + t176;
t169 = pkin(4) * t232 + t171;
t326 = -t169 * t232 + t166;
t325 = -0.2e1 * t292;
t324 = MDP(4) * t267;
t323 = MDP(5) * (t267 ^ 2 - t269 ^ 2);
t320 = t220 ^ 2;
t317 = t187 * t266;
t316 = t218 * t235;
t315 = t220 * t232;
t314 = t232 * t266;
t313 = t235 * t266;
t312 = t247 * t266;
t311 = t247 * t268;
t240 = t264 * t250;
t270 = qJD(2) ^ 2;
t308 = t267 * t270;
t216 = t252 * t264 - t253 * t265;
t213 = t268 * t216;
t225 = t268 * t230;
t307 = t269 * t270;
t271 = qJD(1) ^ 2;
t306 = t269 * t271;
t258 = pkin(2) * t264 + pkin(7);
t305 = qJ(5) + t258;
t304 = t169 - t171;
t303 = -t266 * t188 - t218 * t295;
t302 = t232 * t313 + t225;
t199 = pkin(2) * t294 + pkin(3) * t237 - pkin(7) * t235;
t212 = t249 * t265 + t240;
t301 = t266 * t199 + t268 * t212;
t246 = t264 * t267 - t309;
t208 = pkin(3) * t246 - pkin(7) * t247 + t289;
t300 = t266 * t208 + t213;
t283 = qJD(4) * t305;
t299 = qJ(5) * t313 + qJD(5) * t268 - t266 * t283 - t301;
t195 = t268 * t199;
t298 = -pkin(4) * t237 - t195 + (qJ(5) * t235 - t283) * t268 + (-qJD(5) + t212) * t266;
t291 = t267 * t318;
t198 = t233 * t265 + t234 * t264;
t239 = t246 * qJD(2);
t200 = pkin(3) * t236 + pkin(7) * t239 + t291;
t290 = t268 * t198 + t266 * t200 + t208 * t295;
t259 = -pkin(2) * t265 - pkin(3);
t288 = t247 * t295;
t284 = pkin(1) * t325;
t185 = t227 * t264 - t265 * t228;
t197 = t233 * t264 - t265 * t234;
t209 = t243 * t265 + t240;
t211 = t249 * t264 - t310;
t215 = -t265 * t252 - t253 * t264;
t282 = t232 * t268;
t280 = t185 * t247 - t216 * t230;
t279 = qJ(5) * t239 - qJD(5) * t247;
t175 = pkin(4) * t188 + t185;
t204 = -qJD(2) * pkin(3) - t209;
t277 = -t239 * t266 + t288;
t276 = -t239 * t268 - t247 * t296;
t274 = t204 * t232 - t258 * t230;
t245 = t305 * t268;
t244 = t305 * t266;
t217 = t218 ^ 2;
t203 = t268 * t208;
t196 = t268 * t200;
t179 = pkin(4) * t218 + qJD(5) + t204;
t178 = -qJ(5) * t312 + t300;
t173 = pkin(4) * t246 - qJ(5) * t311 - t216 * t266 + t203;
t168 = -qJ(5) * t288 + (-qJD(4) * t216 + t279) * t266 + t290;
t167 = pkin(4) * t236 - t198 * t266 + t196 + t279 * t268 + (-t213 + (qJ(5) * t247 - t208) * t266) * qJD(4);
t1 = [0.2e1 * t286 * t324 + t323 * t325 + MDP(6) * t307 - MDP(7) * t308 + (-pkin(6) * t307 + t267 * t284) * MDP(9) + (pkin(6) * t308 + t269 * t284) * MDP(10) + (-t186 * t246 + t197 * t237 + t198 * t235 + t209 * t239 - t210 * t236 + t215 * t278 + t280) * MDP(11) + (t185 * t215 + t186 * t216 - t209 * t197 + t210 * t198 + (t251 + t281) * t291) * MDP(12) + (-t187 * t311 + t220 * t276) * MDP(13) + (-(-t218 * t268 - t220 * t266) * t239 + (t317 - t188 * t268 + (t218 * t266 - t220 * t268) * qJD(4)) * t247) * MDP(14) + (-t187 * t246 + t220 * t236 + t225 * t247 + t232 * t276) * MDP(15) + (-t188 * t246 - t218 * t236 - t230 * t312 - t232 * t277) * MDP(16) + (t230 * t246 + t232 * t236) * MDP(17) + ((-t216 * t295 + t196) * t232 + t203 * t230 + (-t205 * t295 + t183) * t246 + t176 * t236 + t197 * t218 + t215 * t188 + t204 * t288 + ((-qJD(4) * t208 - t198) * t232 + (-qJD(4) * t191 - t186) * t246 - t204 * t239 + t280) * t266) * MDP(18) + (-(-t216 * t296 + t290) * t232 - t300 * t230 - t275 * t246 - t177 * t236 + t197 * t220 - t215 * t187 + t185 * t311 + t276 * t204) * MDP(19) + (-t167 * t220 - t168 * t218 + t173 * t187 - t178 * t188 - (-t169 * t268 - t172 * t266) * t239 + (-t165 * t268 - t166 * t266 + (t169 * t266 - t172 * t268) * qJD(4)) * t247) * MDP(20) + (t166 * t178 + t172 * t168 + t165 * t173 + t169 * t167 + t175 * (pkin(4) * t312 + t215) + t179 * (pkin(4) * t277 + t197)) * MDP(21); -t306 * t324 + t271 * t323 + ((t210 - t211) * t237 + (-t212 + t209) * t235 + (-t264 * t230 - t265 * t278) * pkin(2)) * MDP(11) + (t209 * t211 - t210 * t212 + (-t185 * t265 + t186 * t264 - t251 * t294) * pkin(2)) * MDP(12) + (t220 * t282 - t317) * MDP(13) + ((-t187 + t316) * t268 - t220 * t314 + t303) * MDP(14) + (-t220 * t237 + t230 * t266 + t232 * t282) * MDP(15) + (t218 * t237 - t232 * t296 + t302) * MDP(16) - t232 * t237 * MDP(17) + (-t176 * t237 - t185 * t268 + t259 * t188 - t211 * t218 + (-t258 * t295 - t195) * t232 + (t212 * t232 + t274) * t266) * MDP(18) + (t177 * t237 + t185 * t266 - t259 * t187 - t211 * t220 + (t258 * t296 + t301) * t232 + t274 * t268) * MDP(19) + (-t187 * t244 - t188 * t245 - t299 * t218 - t298 * t220 - t327 * t266 + t326 * t268) * MDP(20) + (t166 * t245 - t165 * t244 + t175 * (-pkin(4) * t268 + t259) + (pkin(4) * t314 - t211) * t179 + t299 * t172 + t298 * t169) * MDP(21) + (MDP(9) * t267 * t271 + MDP(10) * t306) * pkin(1); -t235 ^ 2 * MDP(11) + (-t210 * t235 + t256) * MDP(12) + t302 * MDP(18) + t303 * MDP(20) + (-MDP(11) * t237 + MDP(12) * t209 - MDP(18) * t218 - MDP(19) * t220 - MDP(21) * t179) * t237 + (-qJD(4) * t232 * MDP(18) - t230 * MDP(19) + MDP(20) * t315 + t326 * MDP(21)) * t266 + ((t187 + t316) * MDP(20) + t327 * MDP(21) - t232 ^ 2 * MDP(19)) * t268; t220 * t218 * MDP(13) + (-t217 + t320) * MDP(14) + (t218 * t232 - t187) * MDP(15) + (-t188 + t315) * MDP(16) + t230 * MDP(17) + (t177 * t232 - t204 * t220 + t272) * MDP(18) + (t176 * t232 + t204 * t218 - t275) * MDP(19) + (pkin(4) * t187 - t218 * t304) * MDP(20) + (t304 * t172 + (-t179 * t220 + t165) * pkin(4)) * MDP(21); (-t217 - t320) * MDP(20) + (t169 * t220 + t172 * t218 + t175) * MDP(21);];
tauc = t1;
