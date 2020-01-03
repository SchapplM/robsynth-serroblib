% Calculate minimal parameter regressor of Coriolis joint torque vector for
% S5RPRPR11
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% qJD [5x1]
%   Generalized joint velocities
% pkin [8x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d1,d3,d5,theta2]';
% MDP [25x1]
%   Minimal dynamic parameter vector (fixed base model)
%   see S5RPRPR11_convert_par2_MPV_fixb.m
% 
% Output:
% tauc [5x1]
%   joint torques required to compensate Coriolis and centrifugal load

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-31 18:28
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauc = S5RPRPR11_coriolisvecJ_fixb_mdp_slag_vp(qJ, qJD, pkin, MDP)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(8,1),zeros(25,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RPRPR11_coriolisvecJ_fixb_mdp_slag_vp: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RPRPR11_coriolisvecJ_fixb_mdp_slag_vp: qJD has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S5RPRPR11_coriolisvecJ_fixb_mdp_slag_vp: pkin has to be [8x1] (double)');
assert(isreal(MDP) && all(size(MDP) == [25 1]), ...
  'S5RPRPR11_coriolisvecJ_fixb_mdp_slag_vp: MDP has to be [25x1] (double)'); 

%% Symbolic Calculation
% From coriolisvec_joint_fixb_mdp_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 18:27:58
% EndTime: 2019-12-31 18:28:02
% DurationCPUTime: 1.81s
% Computational Cost: add. (1120->216), mult. (2960->284), div. (0->0), fcn. (2162->6), ass. (0->102)
t257 = cos(pkin(8));
t311 = cos(qJ(3));
t289 = t311 * t257;
t279 = qJD(1) * t289;
t245 = qJD(3) * t279;
t256 = sin(pkin(8));
t259 = sin(qJ(3));
t294 = qJD(3) * t259;
t287 = t256 * t294;
t215 = qJD(1) * t287 - t245;
t236 = t311 * t256 + t259 * t257;
t230 = t236 * qJD(3);
t216 = qJD(1) * t230;
t258 = sin(qJ(5));
t260 = cos(qJ(5));
t300 = t256 * t259;
t288 = qJD(1) * t300;
t225 = -t279 + t288;
t227 = t236 * qJD(1);
t319 = t225 * t258 + t260 * t227;
t165 = qJD(5) * t319 - t215 * t258 - t260 * t216;
t272 = -t260 * t225 + t227 * t258;
t280 = qJD(5) * t272 + t260 * t215 - t258 * t216;
t253 = qJD(3) - qJD(5);
t303 = t272 * t253;
t304 = t319 * t253;
t328 = -(t165 + t304) * MDP(22) + t272 * MDP(19) * t319 + (-t272 ^ 2 + t319 ^ 2) * MDP(20) - (t280 + t303) * MDP(21);
t290 = pkin(2) * t257 + pkin(1);
t239 = -t290 * qJD(1) + qJD(2);
t325 = -qJ(4) * t227 + t239;
t309 = pkin(6) + qJ(2);
t241 = t309 * t257;
t238 = qJD(1) * t241;
t222 = t259 * t238;
t240 = t309 * t256;
t237 = qJD(1) * t240;
t203 = -t311 * t237 - t222;
t292 = qJD(4) - t203;
t254 = qJD(3) * qJD(4);
t291 = qJD(1) * qJD(2);
t286 = qJD(2) * t311;
t277 = qJD(1) * t286;
t285 = qJD(3) * t311;
t297 = -t237 * t285 + t257 * t277;
t177 = t254 + (-qJD(3) * t238 - t256 * t291) * t259 + t297;
t167 = pkin(7) * t216 + t177;
t284 = t259 * t291;
t178 = -t237 * t294 + t238 * t285 + t256 * t277 + t257 * t284;
t170 = pkin(7) * t215 + t178;
t261 = -pkin(3) - pkin(4);
t172 = t261 * t225 - t325;
t322 = t258 * t167 - t260 * t170 + t172 * t319;
t206 = -t259 * t240 + t311 * t241;
t293 = -pkin(7) * t227 + t292;
t318 = qJD(5) + t253;
t317 = MDP(13) + MDP(15);
t316 = -MDP(14) + MDP(17);
t315 = (MDP(7) * qJ(2) + MDP(6)) * (t256 ^ 2 + t257 ^ 2);
t314 = -t260 * t167 - t258 * t170 + t172 * t272;
t313 = (t203 + t222) * qJD(3) + t256 * t284 - t297;
t312 = t227 ^ 2;
t310 = pkin(3) * t216;
t307 = qJ(4) * t215;
t306 = qJ(4) * t225;
t204 = -t259 * t237 + t311 * t238;
t184 = pkin(7) * t225 + t204;
t205 = t311 * t240 + t259 * t241;
t176 = t261 * qJD(3) + t293;
t255 = qJD(3) * qJ(4);
t179 = t184 + t255;
t274 = t260 * t176 - t258 * t179;
t273 = -t258 * t176 - t260 * t179;
t235 = -t289 + t300;
t270 = t235 * t260 - t236 * t258;
t202 = t235 * t258 + t236 * t260;
t269 = qJD(4) * t227 - t307;
t229 = -t257 * t285 + t287;
t268 = -qJ(4) * t229 + qJD(4) * t236;
t267 = qJ(4) * t236 + t290;
t266 = qJD(3) * t204 - t178;
t186 = -t240 * t285 + t257 * t286 + (-qJD(2) * t256 - qJD(3) * t241) * t259;
t187 = t236 * qJD(2) + t206 * qJD(3);
t220 = t225 ^ 2;
t200 = pkin(3) * t235 - t267;
t199 = pkin(3) * t227 + t306;
t198 = t255 + t204;
t197 = t245 + (t225 - t288) * qJD(3);
t195 = -qJD(3) * pkin(3) + t292;
t189 = pkin(7) * t235 + t206;
t188 = -t236 * pkin(7) + t205;
t185 = pkin(3) * t225 + t325;
t182 = t261 * t235 + t267;
t181 = pkin(3) * t230 - t268;
t180 = t261 * t227 - t306;
t175 = -t269 + t310;
t174 = t229 * pkin(7) + t187;
t173 = pkin(7) * t230 + t186;
t171 = t261 * t230 + t268;
t169 = t202 * qJD(5) - t229 * t258 - t260 * t230;
t168 = t270 * qJD(5) - t229 * t260 + t230 * t258;
t166 = t261 * t216 + t269;
t1 = [(-t215 * t236 - t227 * t229) * MDP(8) + (t215 * t235 - t216 * t236 + t225 * t229 - t227 * t230) * MDP(9) + (-t216 * t290 + t230 * t239) * MDP(13) + (t215 * t290 - t229 * t239) * MDP(14) + (t175 * t235 + t181 * t225 + t185 * t230 + t200 * t216) * MDP(15) + (-t177 * t235 + t178 * t236 - t186 * t225 + t187 * t227 - t195 * t229 - t198 * t230 - t205 * t215 - t206 * t216) * MDP(16) + (-t175 * t236 - t181 * t227 + t185 * t229 + t200 * t215) * MDP(17) + (t175 * t200 + t177 * t206 + t178 * t205 + t181 * t185 + t186 * t198 + t187 * t195) * MDP(18) + (t168 * t319 - t202 * t280) * MDP(19) + (-t165 * t202 - t168 * t272 - t169 * t319 - t270 * t280) * MDP(20) + (t182 * t165 - t166 * t270 + t172 * t169 + t171 * t272) * MDP(24) + (t166 * t202 + t172 * t168 + t171 * t319 - t182 * t280) * MDP(25) + (-t168 * MDP(21) + t169 * MDP(22) + (t173 * t258 - t174 * t260) * MDP(24) + (t173 * t260 + t174 * t258) * MDP(25) + ((t188 * t258 + t189 * t260) * MDP(24) + (t188 * t260 - t189 * t258) * MDP(25)) * qJD(5)) * t253 + (-t229 * MDP(10) - t230 * MDP(11) + t316 * t186 - t317 * t187) * qJD(3) + 0.2e1 * t291 * t315; (-t220 - t312) * MDP(16) + (t310 + t307 + t198 * t225 + (-qJD(4) - t195) * t227) * MDP(18) + (-t165 + t304) * MDP(24) + (t280 - t303) * MDP(25) - qJD(1) ^ 2 * t315 + 0.2e1 * t317 * qJD(3) * t227 + t316 * (-t245 + (t225 + t288) * qJD(3)); t227 * t225 * MDP(8) + (-t220 + t312) * MDP(9) + t197 * MDP(10) + (-t227 * t239 + t266) * MDP(13) + (t225 * t239 + t313) * MDP(14) + (-t185 * t227 - t199 * t225 + t266) * MDP(15) + (pkin(3) * t215 - qJ(4) * t216 + (t198 - t204) * t227 + (t195 - t292) * t225) * MDP(16) + (-t185 * t225 + t199 * t227 + 0.2e1 * t254 - t313) * MDP(17) + (-pkin(3) * t178 + qJ(4) * t177 - t185 * t199 - t195 * t204 + t292 * t198) * MDP(18) + (-t180 * t272 + (t260 * t184 + t293 * t258) * t253 + (-(-qJ(4) * t260 - t258 * t261) * t253 - t273) * qJD(5) + t322) * MDP(24) + (-t180 * t319 + (-t258 * t184 + t293 * t260) * t253 + ((-qJ(4) * t258 + t260 * t261) * t253 + t274) * qJD(5) - t314) * MDP(25) - t328; t197 * MDP(16) - qJD(3) ^ 2 * MDP(17) + (-qJD(3) * t198 + t178) * MDP(18) + (MDP(15) * t225 - MDP(17) * t227 + MDP(18) * t185 - MDP(24) * t272 - MDP(25) * t319) * t227 + (-t258 * MDP(24) - t260 * MDP(25)) * t253 ^ 2; (t318 * t273 - t322) * MDP(24) + (-t318 * t274 + t314) * MDP(25) + t328;];
tauc = t1;
