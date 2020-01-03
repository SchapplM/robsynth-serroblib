% Calculate minimal parameter regressor of Coriolis joint torque vector for
% S5RRPRP8
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% qJD [5x1]
%   Generalized joint velocities
% pkin [7x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d1,d2,d4]';
% MDP [23x1]
%   Minimal dynamic parameter vector (fixed base model)
%   see S5RRPRP8_convert_par2_MPV_fixb.m
% 
% Output:
% tauc [5x1]
%   joint torques required to compensate Coriolis and centrifugal load

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-31 20:04
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauc = S5RRPRP8_coriolisvecJ_fixb_mdp_slag_vp(qJ, qJD, pkin, MDP)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(7,1),zeros(23,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RRPRP8_coriolisvecJ_fixb_mdp_slag_vp: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RRPRP8_coriolisvecJ_fixb_mdp_slag_vp: qJD has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [7 1]), ...
  'S5RRPRP8_coriolisvecJ_fixb_mdp_slag_vp: pkin has to be [7x1] (double)');
assert(isreal(MDP) && all(size(MDP) == [23 1]), ...
  'S5RRPRP8_coriolisvecJ_fixb_mdp_slag_vp: MDP has to be [23x1] (double)'); 

%% Symbolic Calculation
% From coriolisvec_joint_fixb_mdp_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 20:04:30
% EndTime: 2019-12-31 20:04:34
% DurationCPUTime: 1.44s
% Computational Cost: add. (1055->216), mult. (2466->285), div. (0->0), fcn. (1427->4), ass. (0->111)
t250 = sin(qJ(4));
t251 = sin(qJ(2));
t253 = cos(qJ(2));
t252 = cos(qJ(4));
t283 = qJD(4) * t252;
t284 = qJD(4) * t250;
t285 = qJD(2) * t253;
t322 = t250 * t285 + t251 * t283 - t253 * t284;
t288 = qJD(1) * t251;
t236 = pkin(6) * t288;
t321 = -pkin(7) * t288 + qJD(3) + t236;
t208 = t251 * t250 + t253 * t252;
t193 = t208 * qJD(1);
t192 = t193 ^ 2;
t287 = qJD(1) * t253;
t195 = -t250 * t287 + t252 * t288;
t309 = t195 ^ 2;
t282 = qJD(1) * qJD(2);
t277 = t251 * t282;
t174 = t322 * qJD(1) - t252 * t277;
t244 = qJD(2) - qJD(4);
t314 = t195 * t244 + t174;
t260 = t208 * qJD(4);
t276 = t253 * t282;
t173 = qJD(1) * t260 - t250 * t277 - t252 * t276;
t319 = t193 * t244 + t173;
t320 = -t319 * MDP(17) + t195 * t193 * MDP(15) - (-t309 + t192) * MDP(16) - t314 * MDP(18);
t198 = -qJD(1) * pkin(1) - pkin(2) * t287 - qJ(3) * t288;
t183 = pkin(3) * t287 - t198;
t286 = qJD(2) * t251;
t308 = pkin(6) - pkin(7);
t213 = t308 * t286;
t245 = qJD(2) * qJD(3);
t187 = -qJD(1) * t213 + t245;
t231 = pkin(6) * t276;
t205 = -pkin(7) * t276 + t231;
t254 = -pkin(2) - pkin(3);
t280 = t254 * qJD(2);
t184 = t280 + t321;
t237 = pkin(6) * t287;
t214 = -pkin(7) * t287 + t237;
t246 = qJD(2) * qJ(3);
t197 = t214 + t246;
t265 = t250 * t184 + t252 * t197;
t257 = -t265 * qJD(4) - t250 * t187 + t252 * t205;
t318 = -t183 * t195 + t257;
t316 = pkin(6) * MDP(14);
t247 = t251 ^ 2;
t315 = (-t253 ^ 2 + t247) * MDP(5);
t266 = -t250 * qJ(3) + t252 * t254;
t313 = t266 * qJD(4) - t250 * t214 + t321 * t252;
t217 = t252 * qJ(3) + t250 * t254;
t312 = t217 * qJD(4) + t252 * t214 + t321 * t250;
t307 = pkin(1) * MDP(9);
t306 = t174 * pkin(4);
t305 = pkin(1) * MDP(10);
t304 = qJD(2) * pkin(2);
t302 = t193 * qJ(5);
t300 = t195 * qJ(5);
t256 = qJD(1) ^ 2;
t298 = t253 * t256;
t273 = t252 * t184 - t250 * t197;
t166 = t273 - t300;
t165 = -t244 * pkin(4) + t166;
t297 = t165 - t166;
t296 = -t300 + t313;
t295 = t302 - t312;
t240 = t251 * qJD(3);
t291 = qJ(3) * t276 + qJD(1) * t240;
t290 = qJ(3) * t285 + t240;
t219 = -t253 * pkin(2) - t251 * qJ(3) - pkin(1);
t222 = t308 * t253;
t281 = t254 * MDP(23);
t275 = MDP(12) + t316;
t274 = qJD(3) - t304;
t215 = qJD(2) * t222;
t270 = t250 * t213 + t252 * t215;
t206 = t253 * pkin(3) - t219;
t182 = pkin(2) * t277 - t291;
t191 = pkin(2) * t286 - t290;
t269 = -qJD(1) * t191 - t182;
t268 = t244 ^ 2;
t267 = t251 * t280;
t221 = t308 * t251;
t264 = -t250 * t221 - t252 * t222;
t234 = qJ(3) * t287;
t190 = t254 * t288 + t234;
t263 = t208 * pkin(4) + t206;
t262 = -t184 * t283 - t252 * t187 + t197 * t284 - t250 * t205;
t261 = -t252 * t213 + t250 * t215 + t221 * t283 - t222 * t284;
t177 = qJD(1) * t267 + t291;
t258 = t183 * t193 + t262;
t255 = qJD(2) ^ 2;
t220 = t237 + t246;
t218 = t236 + t274;
t216 = -pkin(6) * t277 + t245;
t211 = -pkin(4) + t266;
t210 = pkin(2) * t288 - t234;
t209 = -t253 * t250 + t251 * t252;
t179 = t267 + t290;
t176 = t208 * qJD(2) - t260;
t175 = -t252 * t286 + t322;
t172 = t193 * pkin(4) + qJD(5) + t183;
t171 = -t208 * qJ(5) - t264;
t170 = -t209 * qJ(5) + t252 * t221 - t250 * t222;
t167 = t265 - t302;
t164 = -t176 * qJ(5) + t264 * qJD(4) - t209 * qJD(5) + t270;
t163 = -t175 * qJ(5) - t208 * qJD(5) + t261;
t162 = t173 * qJ(5) - t195 * qJD(5) + t257;
t161 = -t174 * qJ(5) - t193 * qJD(5) - t262;
t1 = [(t182 * t219 + t198 * t191) * MDP(14) + (-t173 * t209 + t195 * t176) * MDP(15) + (t173 * t208 - t209 * t174 - t195 * t175 - t176 * t193) * MDP(16) + (t206 * t174 + t183 * t175 + t177 * t208 + t179 * t193) * MDP(20) + (-t206 * t173 + t183 * t176 + t177 * t209 + t179 * t195) * MDP(21) + (-t161 * t208 - t162 * t209 - t163 * t193 - t164 * t195 - t165 * t176 - t167 * t175 + t170 * t173 - t171 * t174) * MDP(22) + (t161 * t171 + t167 * t163 + t162 * t170 + t165 * t164 + (t291 + t306) * t263 + t172 * (t175 * pkin(4) + t290)) * MDP(23) + (-t176 * MDP(17) + t175 * MDP(18) + (t221 * t284 + t222 * t283 - t270) * MDP(20) + t261 * MDP(21)) * t244 - 0.2e1 * t282 * t315 + (t269 * MDP(13) + (-MDP(7) + (MDP(10) - MDP(13)) * pkin(6)) * t255 + (t198 * MDP(11) + t172 * t281 - t275 * t220 + (t219 * MDP(11) + t263 * t281 - 0.2e1 * t307) * qJD(1)) * qJD(2)) * t251 + (t255 * MDP(6) + t269 * MDP(11) + t216 * MDP(12) + (t216 * MDP(14) + (-MDP(11) - MDP(9)) * t255) * pkin(6) + (-t198 * MDP(13) + t275 * t218 + (-0.2e1 * t305 - t219 * MDP(13) + (t275 * pkin(6) + (2 * MDP(4))) * t251) * qJD(1)) * qJD(2)) * t253; t256 * t315 + t298 * t305 + 0.2e1 * t245 * MDP(13) + (t216 * qJ(3) + t220 * qJD(3) - t198 * t210) * MDP(14) + (-t190 * t193 + t312 * t244 - t318) * MDP(20) + (-t190 * t195 + t313 * t244 - t258) * MDP(21) + (t211 * t173 - t217 * t174 + (-t167 - t295) * t195 + (t165 - t296) * t193) * MDP(22) + (t161 * t217 + t162 * t211 - t172 * (-t195 * pkin(4) + t190) + t296 * t167 + t295 * t165) * MDP(23) + (-MDP(4) * t298 + t256 * t307) * t251 + ((-t198 * t251 + t210 * t253) * MDP(11) + ((t220 - t246) * t251 + (-t218 + t274) * t253) * MDP(12) + (t198 * t253 + t210 * t251) * MDP(13) + (t220 * t251 + (-t218 - t304) * t253) * t316) * qJD(1) - t320; (-t247 * t256 - t255) * MDP(13) + (-t220 * qJD(2) + t231) * MDP(14) + (-MDP(11) * t298 + (MDP(14) * t198 - MDP(20) * t193 - MDP(21) * t195 - MDP(23) * t172) * qJD(1)) * t251 + (t319 * MDP(22) + (-t244 * t167 + t162) * MDP(23) - MDP(21) * t268) * t252 + (-t314 * MDP(22) + (t244 * t165 + t161) * MDP(23) - MDP(20) * t268) * t250; (-t265 * t244 + t318) * MDP(20) + (-t273 * t244 + t258) * MDP(21) + (pkin(4) * t173 - t297 * t193) * MDP(22) + (t297 * t167 + (-t172 * t195 + t162) * pkin(4)) * MDP(23) + t320; (-t192 - t309) * MDP(22) + (t165 * t195 + t167 * t193 + t177 + t306) * MDP(23);];
tauc = t1;
