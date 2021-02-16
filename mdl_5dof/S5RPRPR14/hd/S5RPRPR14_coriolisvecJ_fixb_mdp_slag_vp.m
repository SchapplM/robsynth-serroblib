% Calculate Coriolis joint torque vector for
% S5RPRPR14
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% qJD [5x1]
%   Generalized joint velocities
% pkin [8x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d1,d3,d5,theta4]';
% MDP [24x1]
%   Minimal dynamic parameter vector (fixed base model)
%   see S5RPRPR14_convert_par2_MPV_fixb.m
% 
% Output:
% tauc [5x1]
%   joint torques required to compensate Coriolis and centrifugal load

% Quelle: HybrDyn-Toolbox
% Datum: 2021-01-15 12:17
% Revision: d12c3222fdeb2c5f3b3c8fa5751e113be2fc3aae (2021-01-15)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauc = S5RPRPR14_coriolisvecJ_fixb_mdp_slag_vp(qJ, qJD, pkin, MDP)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(8,1),zeros(24,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RPRPR14_coriolisvecJ_fixb_mdp_slag_vp: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RPRPR14_coriolisvecJ_fixb_mdp_slag_vp: qJD has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S5RPRPR14_coriolisvecJ_fixb_mdp_slag_vp: pkin has to be [8x1] (double)');
assert(isreal(MDP) && all(size(MDP) == [24 1]), ...
  'S5RPRPR14_coriolisvecJ_fixb_mdp_slag_vp: MDP has to be [24x1] (double)'); 

%% Symbolic Calculation
% From coriolisvec_joint_fixb_mdp_matlab.m
% OptimizationMode: 2
% StartTime: 2021-01-15 12:17:10
% EndTime: 2021-01-15 12:17:17
% DurationCPUTime: 2.23s
% Computational Cost: add. (1310->245), mult. (2866->344), div. (0->0), fcn. (1868->6), ass. (0->111)
t230 = cos(qJ(5));
t262 = t230 * qJD(3);
t231 = cos(qJ(3));
t292 = cos(pkin(8));
t253 = t292 * t231;
t246 = qJD(1) * t253;
t227 = sin(pkin(8));
t229 = sin(qJ(3));
t270 = qJD(1) * t229;
t256 = t227 * t270;
t197 = t246 - t256;
t228 = sin(qJ(5));
t280 = t197 * t228;
t177 = -t262 + t280;
t241 = -t227 * t231 - t229 * t292;
t295 = t241 * qJD(1);
t299 = qJD(5) - t295;
t302 = t177 * t299;
t179 = qJD(3) * t228 + t197 * t230;
t301 = t179 * t299;
t250 = t230 * t299;
t206 = qJD(3) * t246;
t187 = qJD(3) * t256 - t206;
t277 = t228 * t187;
t300 = t250 * t299 - t277;
t188 = qJD(3) * t295;
t298 = MDP(8) * (t229 ^ 2 - t231 ^ 2);
t296 = qJ(2) * MDP(6) + MDP(5);
t232 = -pkin(1) - pkin(6);
t207 = qJD(1) * t232 + qJD(2);
t266 = qJD(4) * t229;
t267 = qJD(3) * t231;
t176 = t207 * t267 + (-qJ(4) * t267 - t266) * qJD(1);
t265 = qJD(4) * t231;
t268 = qJD(3) * t229;
t242 = qJ(4) * t268 - t265;
t235 = qJD(1) * t242 - t207 * t268;
t154 = t176 * t227 - t235 * t292;
t215 = pkin(3) * t227 + pkin(7);
t269 = qJD(1) * t231;
t258 = pkin(3) * t269;
t294 = t299 * (pkin(4) * t197 - pkin(7) * t295 + qJD(5) * t215 + t258) + t154;
t254 = t292 * t176;
t155 = t227 * t235 + t254;
t191 = -qJ(4) * t269 + t207 * t231;
t186 = qJD(3) * pkin(3) + t191;
t190 = (-qJ(4) * qJD(1) + t207) * t229;
t279 = t227 * t190;
t162 = t186 * t292 - t279;
t158 = -qJD(3) * pkin(4) - t162;
t204 = pkin(3) * t270 + qJD(1) * qJ(2) + qJD(4);
t164 = -pkin(4) * t295 - pkin(7) * t197 + t204;
t274 = qJ(4) - t232;
t251 = t274 * t231;
t189 = -qJD(3) * t251 - t266;
t239 = t268 * t274 - t265;
t166 = t189 * t292 + t227 * t239;
t278 = t227 * t229;
t200 = t253 - t278;
t216 = pkin(3) * t229 + qJ(2);
t172 = -pkin(4) * t241 - pkin(7) * t200 + t216;
t252 = qJD(3) * t292;
t196 = -t227 * t267 - t229 * t252;
t205 = t274 * t229;
t175 = -t205 * t292 - t227 * t251;
t285 = t175 * t187;
t290 = t154 * t200;
t293 = t158 * t196 - (qJD(5) * t172 + t166) * t299 + (qJD(5) * t164 + t155) * t241 + t285 + t290;
t289 = t158 * t200;
t264 = qJD(5) * t228;
t273 = qJD(5) * t262 + t188 * t230;
t160 = -t197 * t264 + t273;
t288 = t160 * t200;
t287 = t160 * t228;
t286 = t172 * t187;
t284 = t177 * t197;
t283 = t179 * t197;
t282 = t187 * t241;
t281 = t188 * t228;
t182 = t230 * t187;
t234 = qJD(1) ^ 2;
t276 = t231 * t234;
t233 = qJD(3) ^ 2;
t275 = t232 * t233;
t184 = t292 * t190;
t163 = t186 * t227 + t184;
t223 = qJD(1) * qJD(2);
t260 = qJD(1) * qJD(3);
t255 = t231 * t260;
t203 = pkin(3) * t255 + t223;
t263 = qJD(5) * t230;
t208 = pkin(3) * t267 + qJD(2);
t259 = 0.2e1 * qJD(1);
t159 = qJD(3) * pkin(7) + t163;
t153 = t159 * t230 + t164 * t228;
t245 = t159 * t228 - t164 * t230;
t244 = -t182 + (t228 * t295 - t264) * t299;
t243 = t196 * t230 - t200 * t264;
t168 = t191 * t292 - t279;
t237 = t215 * t187 + (t158 + t168) * t299;
t195 = t227 * t268 - t231 * t252;
t236 = -t155 * t241 + t162 * t196 - t163 * t195 - t290;
t217 = -pkin(3) * t292 - pkin(4);
t174 = -t205 * t227 + t251 * t292;
t169 = -pkin(4) * t195 - pkin(7) * t196 + t208;
t167 = t191 * t227 + t184;
t165 = t189 * t227 - t239 * t292;
t161 = qJD(5) * t179 + t281;
t157 = -pkin(4) * t187 - pkin(7) * t188 + t203;
t156 = t230 * t157;
t1 = [0.2e1 * t260 * t298 - t233 * t231 * MDP(10) + qJ(2) * t267 * t259 * MDP(12) + (-t231 * t275 + (-qJ(2) * t268 + qJD(2) * t231) * t259) * MDP(13) + (-qJD(3) * t165 - t187 * t216 - t195 * t204 - t203 * t241 - t208 * t295) * MDP(14) + (-qJD(3) * t166 + t188 * t216 + t196 * t204 + t197 * t208 + t200 * t203) * MDP(15) + (t165 * t197 + t166 * t295 + t174 * t188 - t236 + t285) * MDP(16) + (t154 * t174 + t155 * t175 - t162 * t165 + t163 * t166 + t203 * t216 + t204 * t208) * MDP(17) + (t179 * t243 + t230 * t288) * MDP(18) + ((-t177 * t230 - t179 * t228) * t196 + (-t287 - t161 * t230 + (t177 * t228 - t179 * t230) * qJD(5)) * t200) * MDP(19) + (-t160 * t241 - t179 * t195 - t182 * t200 + t243 * t299) * MDP(20) + (t200 * t277 + t161 * t241 + t177 * t195 + (-t196 * t228 - t200 * t263) * t299) * MDP(21) + (-t195 * t299 + t282) * MDP(22) + (t245 * t195 - t156 * t241 + t174 * t161 + t165 * t177 + (t169 * t299 - t286 + (t159 * t241 - t175 * t299 + t289) * qJD(5)) * t230 + t293 * t228) * MDP(23) + (t153 * t195 + t174 * t160 + t165 * t179 + (-(-qJD(5) * t175 + t169) * t299 + t286 + (-qJD(5) * t159 + t157) * t241 - qJD(5) * t289) * t228 + t293 * t230) * MDP(24) + 0.2e1 * t296 * t223 + (-0.2e1 * MDP(7) * t255 - t233 * MDP(9) + (qJD(2) * t259 - t275) * MDP(12)) * t229; (qJD(1) * t295 + qJD(3) * t196) * MDP(14) + (-qJD(1) * t197 + qJD(3) * t195) * MDP(15) + (-t188 * t200 - t195 * t295 - t196 * t197 - t282) * MDP(16) + (-qJD(1) * t204 + t236) * MDP(17) + (-t161 * t200 - t177 * t196 - t241 * t277) * MDP(23) + (-t179 * t196 - t182 * t241 - t288) * MDP(24) - t296 * t234 + ((-qJD(1) * t230 + t195 * t228 + t241 * t263) * MDP(23) + (qJD(1) * t228 + t195 * t230 - t241 * t264) * MDP(24)) * t299 + (MDP(12) * t229 + MDP(13) * t231) * (-t233 - t234); t229 * MDP(7) * t276 - t234 * t298 + (qJD(3) * t167 - t197 * t204 + t258 * t295 - t154) * MDP(14) + (-t254 - t204 * t295 + (t207 * t278 + t168) * qJD(3) + (-t231 * pkin(3) * t197 - t227 * t242) * qJD(1)) * MDP(15) + ((t163 - t167) * t197 - (-t162 + t168) * t295 + (t187 * t227 - t188 * t292) * pkin(3)) * MDP(16) + (t162 * t167 - t163 * t168 + (-t154 * t292 + t155 * t227 - t204 * t269) * pkin(3)) * MDP(17) + (t179 * t250 + t287) * MDP(18) + ((t160 - t302) * t230 + (-t161 - t301) * t228) * MDP(19) + (-t283 + t300) * MDP(20) + (t244 + t284) * MDP(21) - t299 * t197 * MDP(22) + (t217 * t161 - t167 * t177 + t197 * t245 + t237 * t228 - t230 * t294) * MDP(23) + (t153 * t197 + t217 * t160 - t167 * t179 + t228 * t294 + t237 * t230) * MDP(24) + (MDP(13) * t229 * t234 - MDP(12) * t276) * qJ(2); (t206 + (t197 - t256) * qJD(3)) * MDP(14) + 0.2e1 * t188 * MDP(15) + (-t197 ^ 2 - t295 ^ 2) * MDP(16) + (t162 * t197 - t163 * t295 + t203) * MDP(17) + (t244 - t284) * MDP(23) + (-t283 - t300) * MDP(24); t179 * t177 * MDP(18) + (-t177 ^ 2 + t179 ^ 2) * MDP(19) + (t273 + t302) * MDP(20) + (-t281 + t301) * MDP(21) - t187 * MDP(22) + (t153 * t299 - t155 * t228 - t158 * t179 + t156) * MDP(23) + (-t155 * t230 - t157 * t228 + t158 * t177 - t245 * t299) * MDP(24) + (-MDP(20) * t280 - MDP(21) * t179 - MDP(23) * t153 + MDP(24) * t245) * qJD(5);];
tauc = t1;
