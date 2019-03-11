% Calculate minimal parameter regressor of Coriolis joint torque vector for
% S6RPRPPR8
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% qJD [6x1]
%   Generalized joint velocities
% pkin [8x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d3,d6]';
% MDP [28x1]
%   Minimal dynamic parameter vector (fixed base model)
%   see S6RPRPPR8_convert_par2_MPV_fixb.m
% 
% Output:
% tauc [6x1]
%   joint torques required to compensate Coriolis and centrifugal load

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-09 03:00
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauc = S6RPRPPR8_coriolisvecJ_fixb_mdp_slag_vp(qJ, qJD, pkin, MDP)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(8,1),zeros(28,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPRPPR8_coriolisvecJ_fixb_mdp_slag_vp: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RPRPPR8_coriolisvecJ_fixb_mdp_slag_vp: qJD has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S6RPRPPR8_coriolisvecJ_fixb_mdp_slag_vp: pkin has to be [8x1] (double)');
assert(isreal(MDP) && all(size(MDP) == [28 1]), ...
  'S6RPRPPR8_coriolisvecJ_fixb_mdp_slag_vp: MDP has to be [28x1] (double)'); 

%% Symbolic Calculation
% From coriolisvec_joint_fixb_mdp_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 03:00:03
% EndTime: 2019-03-09 03:00:08
% DurationCPUTime: 2.40s
% Computational Cost: add. (1068->297), mult. (2054->392), div. (0->0), fcn. (964->4), ass. (0->137)
t228 = -pkin(3) - pkin(4);
t229 = -pkin(1) - pkin(7);
t207 = qJD(1) * t229 + qJD(2);
t227 = cos(qJ(3));
t186 = (qJ(5) * qJD(1) + t207) * t227;
t320 = qJD(4) - t186;
t173 = qJD(3) * t228 + t320;
t285 = qJD(1) * t227;
t210 = qJD(6) + t285;
t224 = sin(qJ(6));
t226 = cos(qJ(6));
t321 = (-MDP(27) * t226 + MDP(28) * t224) * t210;
t299 = t207 * t227;
t188 = (qJD(4) + t299) * qJD(3);
t311 = qJD(3) * pkin(3);
t255 = -qJD(4) + t311;
t189 = -t255 - t299;
t225 = sin(qJ(3));
t198 = t225 * t207;
t220 = qJD(3) * qJ(4);
t191 = t198 + t220;
t304 = t191 * t227;
t319 = (t225 * (-t189 + t299) - t304) * qJD(3) - t188 * t225;
t267 = 0.2e1 * qJD(1);
t318 = MDP(7) * t225;
t221 = t225 ^ 2;
t222 = t227 ^ 2;
t317 = MDP(8) * (t221 - t222);
t316 = qJ(2) * MDP(6) + MDP(5);
t219 = -pkin(8) + t228;
t237 = pkin(5) * t227 + t219 * t225 - qJ(2);
t214 = qJ(4) * t285;
t272 = qJD(5) + t214;
t166 = qJD(1) * t237 + t272;
t270 = qJD(1) * qJD(5);
t283 = qJD(3) * t225;
t196 = t207 * t283;
t271 = qJD(1) * qJD(3);
t259 = t225 * t271;
t291 = qJ(5) * t259 + t196;
t172 = -t227 * t270 + t291;
t286 = qJD(1) * t225;
t213 = qJ(5) * t286;
t181 = -t213 - t191;
t175 = qJD(3) * pkin(5) - t181;
t217 = t227 * qJ(4);
t179 = t217 + t237;
t292 = qJ(5) + t229;
t253 = qJD(3) * t292;
t182 = -qJD(5) * t227 + t225 * t253;
t313 = -(qJD(6) * t179 + t182) * t210 + (qJD(3) * t175 - qJD(6) * t166 - t172) * t227;
t223 = qJ(4) + pkin(5);
t309 = qJ(4) * t225;
t168 = t225 * t270 + (qJD(4) + t186) * qJD(3);
t308 = t168 * t224;
t307 = t168 * t226;
t261 = t224 * t286;
t282 = qJD(3) * t226;
t193 = t261 + t282;
t258 = t227 * t271;
t206 = t226 * t258;
t176 = -qJD(6) * t193 + t206;
t306 = t176 * t224;
t303 = t193 * t210;
t302 = t193 * t225;
t260 = t226 * t286;
t275 = t224 * qJD(3);
t195 = t260 - t275;
t301 = t195 * t210;
t300 = t195 * t225;
t298 = t210 * t226;
t297 = t224 * t210;
t231 = qJD(1) ^ 2;
t296 = t225 * t231;
t295 = t226 * t227;
t294 = t227 * t231;
t230 = qJD(3) ^ 2;
t293 = t229 * t230;
t289 = t230 + t231;
t284 = qJD(3) * t181;
t281 = qJD(3) * t227;
t280 = qJD(6) * t224;
t279 = qJD(6) * t226;
t278 = qJD(6) * t227;
t277 = t175 * qJD(6);
t256 = pkin(3) * t225 + qJ(2);
t190 = qJD(1) * t256 - t214;
t276 = t190 * MDP(17);
t216 = t227 * qJD(4);
t247 = t225 * t228 - qJ(2);
t180 = qJD(1) * t247 + t272;
t273 = qJD(5) + t180;
t269 = MDP(14) - MDP(19);
t268 = MDP(16) + MDP(18);
t265 = t227 * t297;
t264 = t210 * t295;
t263 = t225 * t280;
t262 = t225 * t279;
t257 = MDP(26) * t283;
t211 = qJD(1) * t216;
t240 = t227 * t228 - t309;
t233 = qJD(3) * t240 - qJD(2);
t167 = qJD(1) * t233 + t211;
t178 = t216 + t233;
t251 = qJD(1) * t178 + t167;
t192 = t217 + t247;
t250 = qJD(1) * t192 + t180;
t246 = -t210 * t279 + t224 * t259;
t245 = t210 * t280 + t226 * t259;
t244 = pkin(3) * t227 + t309;
t171 = qJD(3) * t219 + t320;
t161 = t166 * t226 - t171 * t224;
t162 = t166 * t224 + t171 * t226;
t200 = t292 * t227;
t235 = t219 * t227 - t223 * t225;
t232 = qJD(3) * t235 - qJD(2);
t243 = (qJD(6) * t200 + t216 + t232) * t210;
t242 = qJD(1) * t221 - t210 * t227;
t201 = -t217 + t256;
t241 = qJD(3) * (qJD(1) * t201 + t190);
t234 = qJD(3) * t244 + qJD(2);
t174 = qJD(1) * t234 - t211;
t184 = -t216 + t234;
t239 = qJD(1) * t184 + t174 - t293;
t218 = 0.2e1 * qJD(3) * qJD(4);
t212 = qJD(6) * t275;
t199 = t292 * t225;
t197 = t244 * qJD(1);
t187 = t240 * qJD(1);
t185 = t198 + t213;
t183 = qJD(5) * t225 + t227 * t253;
t177 = -t212 + (t227 * t275 + t262) * qJD(1);
t170 = t235 * qJD(1);
t169 = t173 * t283;
t164 = qJD(1) * t232 + t211;
t163 = t226 * t164;
t1 = [-0.2e1 * t258 * t318 + 0.2e1 * t271 * t317 + (-t225 * t293 + (qJ(2) * t281 + qJD(2) * t225) * t267) * MDP(12) + (-t227 * t293 + (-qJ(2) * t283 + qJD(2) * t227) * t267) * MDP(13) + (t225 * t239 + t227 * t241) * MDP(14) + t319 * MDP(15) + (t225 * t241 - t227 * t239) * MDP(16) + (t174 * t201 + t184 * t190 - t229 * t319) * MDP(17) + (t251 * t227 + (-t225 * t250 + t183) * qJD(3)) * MDP(18) + (t251 * t225 + (t227 * t250 + t182) * qJD(3)) * MDP(19) + (t168 * t225 + t169 + (-t172 - t284) * t227 + (-t182 * t227 + t183 * t225 + (t199 * t227 - t200 * t225) * qJD(3)) * qJD(1)) * MDP(20) + (t167 * t192 + t168 * t199 - t172 * t200 + t173 * t182 + t178 * t180 - t181 * t183) * MDP(21) + (t176 * t225 * t226 + (t226 * t281 - t263) * t195) * MDP(22) + ((-t193 * t226 - t195 * t224) * t281 + (-t306 - t177 * t226 + (t193 * t224 - t195 * t226) * qJD(6)) * t225) * MDP(23) + (-t210 * t263 + t176 * t227 + (-t226 * t242 - t300) * qJD(3)) * MDP(24) + (-t210 * t262 - t177 * t227 + (t224 * t242 + t302) * qJD(3)) * MDP(25) + (-t210 - t285) * t257 + (t163 * t227 + t199 * t177 + t183 * t193 + (-t171 * t278 + t243) * t226 + t313 * t224 + (t226 * t277 + t308 + (-(t179 * t226 + t200 * t224) * qJD(1) - t161) * qJD(3)) * t225) * MDP(27) + (t199 * t176 + t183 * t195 + (-t243 - (-qJD(6) * t171 + t164) * t227) * t224 + t313 * t226 + (-t224 * t277 + t307 + ((t179 * t224 - t200 * t226) * qJD(1) + t162) * qJD(3)) * t225) * MDP(28) + t316 * qJD(2) * t267 + (-MDP(10) * t227 - MDP(9) * t225) * t230; -qJD(1) * t276 + (qJD(1) * t180 + t169) * MDP(21) - t316 * t231 - (qJD(1) + t278) * t321 + (-t172 * MDP(21) + (-MDP(13) + t268) * t289 + ((t191 - t198) * MDP(17) - t181 * MDP(21) + (t193 - t261) * MDP(27) + (t195 - t260) * MDP(28)) * qJD(3)) * t227 + (t188 * MDP(17) + t168 * MDP(21) + t177 * MDP(27) + t176 * MDP(28) + (-MDP(12) - t269) * t289 + (MDP(17) * t189 + (-t224 * MDP(27) - t226 * MDP(28)) * t210) * qJD(3)) * t225; t294 * t318 - t231 * t317 + t218 * MDP(16) + (qJ(4) * t188 + qJD(4) * t191 - t190 * t197 + (-t304 + (-t189 - t311) * t225) * t207) * MDP(17) + (t218 + (-t186 + t299) * qJD(3)) * MDP(18) + (-qJD(3) * t185 + t291) * MDP(19) + (qJ(4) * t168 + t172 * t228 - t173 * t185 - t180 * t187 - t181 * t320) * MDP(21) + (-t195 * t298 - t306) * MDP(22) + ((-t176 + t303) * t226 + (t177 + t301) * t224) * MDP(23) + t246 * MDP(24) + t245 * MDP(25) + t210 * MDP(26) * t286 + (t223 * t177 + t307 - (t170 * t226 - t185 * t224) * t210 + t320 * t193 + (-t175 * t224 - t219 * t298) * qJD(6)) * MDP(27) + (t223 * t176 - t308 + (t170 * t224 + t185 * t226) * t210 + t320 * t195 + (-t175 * t226 + t219 * t297) * qJD(6)) * MDP(28) + (-MDP(12) * t294 + MDP(13) * t296) * qJ(2) + ((-t190 * t227 - t197 * t225) * MDP(14) + ((t191 - t220) * t227 + (t189 + t255) * t225) * MDP(15) + (-t190 * t225 + t197 * t227) * MDP(16) + ((qJ(5) * qJD(3) - t187) * t227 + t273 * t225) * MDP(18) + (-t187 * t225 - t227 * t273) * MDP(19) + (t181 + t185 + t220) * t227 * MDP(20) + (-t264 + t300) * MDP(24) + (t265 - t302) * MDP(25) + (t161 * t225 + (-t175 * t227 + t219 * t283) * t224) * MDP(27) + (-t175 * t295 + (t219 * t282 - t162) * t225) * MDP(28)) * qJD(1); (-qJD(3) * t191 + t196) * MDP(17) + (t284 + t291) * MDP(21) + (-qJD(3) * t193 + t246) * MDP(27) + (-qJD(3) * t195 + t245) * MDP(28) + t268 * (-t222 * t231 - t230) + (t269 * t296 + (-t273 * MDP(21) + t276 + t321) * qJD(1)) * t227; t211 * MDP(21) - t245 * MDP(27) + t246 * MDP(28) + (-t221 - t222) * MDP(20) * t231 + ((t173 * t227 + t181 * t225 - qJD(2)) * MDP(21) + (-t265 - t302) * MDP(27) + (-t264 - t300) * MDP(28) + ((-MDP(21) * qJ(4) - 0.2e1 * MDP(18)) * t225 + (MDP(21) * t228 + 0.2e1 * MDP(19)) * t227) * qJD(3)) * qJD(1); t195 * t193 * MDP(22) + (-t193 ^ 2 + t195 ^ 2) * MDP(23) + (t206 + t303) * MDP(24) + (-t224 * t258 + t212 + t301) * MDP(25) - qJD(1) * t257 + (t162 * t210 - t172 * t224 - t175 * t195 + t163) * MDP(27) + (t161 * t210 - t164 * t224 - t172 * t226 + t175 * t193) * MDP(28) + (-MDP(24) * t193 - MDP(25) * t260 - MDP(27) * t162 - MDP(28) * t161) * qJD(6);];
tauc  = t1;
