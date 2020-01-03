% Calculate minimal parameter regressor of Coriolis joint torque vector for
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
% MDP [22x1]
%   Minimal dynamic parameter vector (fixed base model)
%   see S5RPRPR14_convert_par2_MPV_fixb.m
% 
% Output:
% tauc [5x1]
%   joint torques required to compensate Coriolis and centrifugal load

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-31 18:35
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauc = S5RPRPR14_coriolisvecJ_fixb_mdp_slag_vp(qJ, qJD, pkin, MDP)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(8,1),zeros(22,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RPRPR14_coriolisvecJ_fixb_mdp_slag_vp: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RPRPR14_coriolisvecJ_fixb_mdp_slag_vp: qJD has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S5RPRPR14_coriolisvecJ_fixb_mdp_slag_vp: pkin has to be [8x1] (double)');
assert(isreal(MDP) && all(size(MDP) == [22 1]), ...
  'S5RPRPR14_coriolisvecJ_fixb_mdp_slag_vp: MDP has to be [22x1] (double)'); 

%% Symbolic Calculation
% From coriolisvec_joint_fixb_mdp_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 18:35:19
% EndTime: 2019-12-31 18:35:24
% DurationCPUTime: 1.68s
% Computational Cost: add. (1204->216), mult. (2634->309), div. (0->0), fcn. (1732->6), ass. (0->104)
t224 = cos(qJ(5));
t254 = t224 * qJD(3);
t221 = sin(pkin(8));
t225 = cos(qJ(3));
t285 = cos(pkin(8));
t246 = t285 * t225;
t223 = sin(qJ(3));
t262 = qJD(1) * t223;
t192 = qJD(1) * t246 - t221 * t262;
t222 = sin(qJ(5));
t273 = t192 * t222;
t173 = -t254 + t273;
t235 = t221 * t225 + t223 * t285;
t288 = t235 * qJD(1);
t292 = qJD(5) + t288;
t295 = t173 * t292;
t175 = qJD(3) * t222 + t192 * t224;
t294 = t175 * t292;
t243 = t292 * t224;
t183 = t192 * qJD(3);
t271 = t222 * t183;
t293 = -t243 * t292 - t271;
t291 = MDP(8) * (t223 ^ 2 - t225 ^ 2);
t289 = qJ(2) * MDP(6) + MDP(5);
t226 = -pkin(1) - pkin(6);
t201 = qJD(1) * t226 + qJD(2);
t258 = qJD(4) * t223;
t259 = qJD(3) * t225;
t172 = t201 * t259 + (-qJ(4) * t259 - t258) * qJD(1);
t257 = qJD(4) * t225;
t260 = qJD(3) * t223;
t229 = -t201 * t260 + (qJ(4) * t260 - t257) * qJD(1);
t150 = t172 * t221 - t229 * t285;
t210 = pkin(3) * t221 + pkin(7);
t261 = qJD(1) * t225;
t287 = t292 * (pkin(3) * t261 + pkin(4) * t192 + pkin(7) * t288 + qJD(5) * t210) + t150;
t151 = t172 * t285 + t221 * t229;
t187 = -qJ(4) * t261 + t225 * t201;
t182 = qJD(3) * pkin(3) + t187;
t186 = (-qJ(4) * qJD(1) + t201) * t223;
t272 = t221 * t186;
t158 = t182 * t285 - t272;
t154 = -qJD(3) * pkin(4) - t158;
t198 = pkin(3) * t262 + qJD(1) * qJ(2) + qJD(4);
t160 = pkin(4) * t288 - pkin(7) * t192 + t198;
t267 = qJ(4) - t226;
t244 = t267 * t225;
t185 = -qJD(3) * t244 - t258;
t233 = t260 * t267 - t257;
t162 = t185 * t285 + t221 * t233;
t195 = -t221 * t223 + t246;
t268 = t223 * pkin(3) + qJ(2);
t168 = pkin(4) * t235 - pkin(7) * t195 + t268;
t245 = qJD(3) * t285;
t191 = -t221 * t259 - t223 * t245;
t199 = t267 * t223;
t171 = -t199 * t285 - t221 * t244;
t278 = t171 * t183;
t283 = t150 * t195;
t286 = t154 * t191 - (qJD(5) * t168 + t162) * t292 - (qJD(5) * t160 + t151) * t235 - t278 + t283;
t282 = t154 * t195;
t256 = qJD(5) * t222;
t184 = qJD(3) * t288;
t265 = qJD(5) * t254 - t224 * t184;
t156 = -t192 * t256 + t265;
t281 = t156 * t195;
t280 = t156 * t222;
t279 = t168 * t183;
t277 = t173 * t192;
t276 = t175 * t192;
t275 = t183 * t235;
t274 = t184 * t222;
t178 = t224 * t183;
t228 = qJD(1) ^ 2;
t270 = t225 * t228;
t227 = qJD(3) ^ 2;
t269 = t226 * t227;
t180 = t285 * t186;
t159 = t221 * t182 + t180;
t217 = qJD(1) * qJD(2);
t251 = qJD(1) * qJD(3);
t247 = t225 * t251;
t266 = pkin(3) * t247 + t217;
t255 = qJD(5) * t224;
t253 = pkin(3) * t259 + qJD(2);
t250 = 0.2e1 * qJD(1);
t155 = qJD(3) * pkin(7) + t159;
t149 = t155 * t224 + t160 * t222;
t238 = t155 * t222 - t160 * t224;
t237 = t178 + (-t222 * t288 - t256) * t292;
t236 = t191 * t224 - t195 * t256;
t164 = t187 * t285 - t272;
t231 = -t210 * t183 + (t154 + t164) * t292;
t190 = t221 * t260 - t225 * t245;
t230 = t151 * t235 + t158 * t191 - t159 * t190 - t283;
t211 = -pkin(3) * t285 - pkin(4);
t170 = -t199 * t221 + t244 * t285;
t165 = -pkin(4) * t190 - pkin(7) * t191 + t253;
t163 = t187 * t221 + t180;
t161 = t185 * t221 - t233 * t285;
t157 = qJD(5) * t175 - t274;
t153 = pkin(4) * t183 + pkin(7) * t184 + t266;
t152 = t224 * t153;
t1 = [0.2e1 * t251 * t291 - t227 * t225 * MDP(10) + qJ(2) * t259 * t250 * MDP(12) + (-t225 * t269 + (-qJ(2) * t260 + qJD(2) * t225) * t250) * MDP(13) + (t161 * t192 - t162 * t288 - t170 * t184 - t230 - t278) * MDP(14) + (t150 * t170 + t151 * t171 - t158 * t161 + t159 * t162 + t198 * t253 + t266 * t268) * MDP(15) + (t175 * t236 + t224 * t281) * MDP(16) + ((-t173 * t224 - t175 * t222) * t191 + (-t280 - t157 * t224 + (t173 * t222 - t175 * t224) * qJD(5)) * t195) * MDP(17) + (t156 * t235 - t175 * t190 + t178 * t195 + t236 * t292) * MDP(18) + (-t195 * t271 - t157 * t235 + t173 * t190 + (-t191 * t222 - t195 * t255) * t292) * MDP(19) + (-t190 * t292 + t275) * MDP(20) + (t238 * t190 + t152 * t235 + t170 * t157 + t161 * t173 + (t165 * t292 + t279 + (-t155 * t235 - t171 * t292 + t282) * qJD(5)) * t224 + t286 * t222) * MDP(21) + (t149 * t190 + t170 * t156 + t161 * t175 + (-(-qJD(5) * t171 + t165) * t292 - t279 - (-qJD(5) * t155 + t153) * t235 - qJD(5) * t282) * t222 + t286 * t224) * MDP(22) + 0.2e1 * t289 * t217 + (-0.2e1 * MDP(7) * t247 - t227 * MDP(9) + (qJD(2) * t250 - t269) * MDP(12)) * t223; (t184 * t195 + t190 * t288 - t191 * t192 - t275) * MDP(14) + (-qJD(1) * t198 + t230) * MDP(15) + (-t157 * t195 - t173 * t191 - t235 * t271) * MDP(21) + (-t175 * t191 - t178 * t235 - t281) * MDP(22) - t289 * t228 + ((-qJD(1) * t224 + t190 * t222 - t235 * t255) * MDP(21) + (qJD(1) * t222 + t190 * t224 + t235 * t256) * MDP(22)) * t292 + (MDP(12) * t223 + MDP(13) * t225) * (-t227 - t228); t223 * MDP(7) * t270 - t228 * t291 + ((t159 - t163) * t192 - (t158 - t164) * t288 + (-t183 * t221 + t184 * t285) * pkin(3)) * MDP(14) + (t158 * t163 - t159 * t164 + (-t150 * t285 + t151 * t221 - t198 * t261) * pkin(3)) * MDP(15) + (t175 * t243 + t280) * MDP(16) + ((t156 - t295) * t224 + (-t157 - t294) * t222) * MDP(17) + (-t276 - t293) * MDP(18) + (t237 + t277) * MDP(19) - t292 * t192 * MDP(20) + (t211 * t157 - t163 * t173 + t192 * t238 + t231 * t222 - t287 * t224) * MDP(21) + (t149 * t192 + t211 * t156 - t163 * t175 + t287 * t222 + t231 * t224) * MDP(22) + (MDP(13) * t223 * t228 - MDP(12) * t270) * qJ(2); (-t192 ^ 2 - t288 ^ 2) * MDP(14) + (t158 * t192 + t159 * t288 + t266) * MDP(15) + (t237 - t277) * MDP(21) + (-t276 + t293) * MDP(22); t175 * t173 * MDP(16) + (-t173 ^ 2 + t175 ^ 2) * MDP(17) + (t265 + t295) * MDP(18) + (t274 + t294) * MDP(19) + t183 * MDP(20) + (t149 * t292 - t151 * t222 - t154 * t175 + t152) * MDP(21) + (-t151 * t224 - t153 * t222 + t154 * t173 - t238 * t292) * MDP(22) + (-MDP(18) * t273 - MDP(19) * t175 - MDP(21) * t149 + MDP(22) * t238) * qJD(5);];
tauc = t1;
