% Calculate minimal parameter regressor of Coriolis joint torque vector for
% S5RPRRP7
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
% MDP [22x1]
%   Minimal dynamic parameter vector (fixed base model)
%   see S5RPRRP7_convert_par2_MPV_fixb.m
% 
% Output:
% tauc [5x1]
%   joint torques required to compensate Coriolis and centrifugal load

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-31 18:46
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauc = S5RPRRP7_coriolisvecJ_fixb_mdp_slag_vp(qJ, qJD, pkin, MDP)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(8,1),zeros(22,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RPRRP7_coriolisvecJ_fixb_mdp_slag_vp: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RPRRP7_coriolisvecJ_fixb_mdp_slag_vp: qJD has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S5RPRRP7_coriolisvecJ_fixb_mdp_slag_vp: pkin has to be [8x1] (double)');
assert(isreal(MDP) && all(size(MDP) == [22 1]), ...
  'S5RPRRP7_coriolisvecJ_fixb_mdp_slag_vp: MDP has to be [22x1] (double)'); 

%% Symbolic Calculation
% From coriolisvec_joint_fixb_mdp_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 18:45:41
% EndTime: 2019-12-31 18:45:45
% DurationCPUTime: 2.00s
% Computational Cost: add. (1512->286), mult. (3450->394), div. (0->0), fcn. (1986->6), ass. (0->118)
t245 = cos(qJ(3));
t292 = qJD(1) * t245;
t230 = -qJD(4) + t292;
t243 = sin(qJ(3));
t323 = MDP(5) * t243;
t238 = t243 ^ 2;
t322 = (-t245 ^ 2 + t238) * MDP(6);
t232 = sin(pkin(8)) * pkin(1) + pkin(6);
t226 = t232 * qJD(1);
t236 = t243 * qJD(2);
t203 = t245 * t226 + t236;
t196 = qJD(3) * pkin(7) + t203;
t320 = qJD(2) * t245 - t243 * t226;
t197 = t320 * qJD(3);
t233 = -cos(pkin(8)) * pkin(1) - pkin(2);
t213 = -pkin(3) * t245 - pkin(7) * t243 + t233;
t199 = t213 * qJD(1);
t261 = pkin(3) * t243 - pkin(7) * t245;
t224 = t261 * qJD(3);
t212 = qJD(1) * t224;
t242 = sin(qJ(4));
t244 = cos(qJ(4));
t285 = qJD(4) * t244;
t286 = qJD(4) * t242;
t262 = t196 * t285 + t242 * t197 + t199 * t286 - t244 * t212;
t280 = qJD(1) * qJD(3);
t266 = t243 * t280;
t168 = -pkin(4) * t266 + t262;
t177 = t196 * t244 + t199 * t242;
t173 = -qJ(5) * t230 + t177;
t321 = t173 * t230 + t168;
t303 = t244 * t245;
t296 = t242 * t213 + t232 * t303;
t278 = MDP(17) + MDP(19);
t265 = t245 * t280;
t269 = t243 * t286;
t279 = qJD(3) * qJD(4);
t191 = qJD(1) * t269 + (-t265 - t279) * t244;
t264 = t242 * t279;
t268 = t243 * t285;
t288 = qJD(3) * t245;
t192 = t264 + (t288 * t242 + t268) * qJD(1);
t198 = qJD(3) * t236 + t226 * t288;
t290 = qJD(3) * t242;
t293 = qJD(1) * t243;
t220 = t244 * t293 + t290;
t169 = pkin(4) * t192 + qJ(5) * t191 - qJD(5) * t220 + t198;
t319 = t169 * t242;
t318 = t169 * t244;
t195 = -qJD(3) * pkin(3) - t320;
t316 = t195 * t242;
t315 = t195 * t244;
t314 = t198 * t242;
t313 = t198 * t244;
t312 = t213 * t244;
t282 = t244 * qJD(3);
t218 = t242 * t293 - t282;
t311 = t218 * t230;
t310 = t224 * t244;
t309 = t230 * t244;
t308 = t232 * t242;
t307 = t232 * t244;
t306 = t242 * t245;
t305 = t243 * t244;
t246 = qJD(3) ^ 2;
t304 = t243 * t246;
t302 = t245 * t246;
t259 = pkin(4) * t242 - qJ(5) * t244;
t301 = qJD(5) * t242 + t230 * t259 + t203;
t271 = t245 * t282;
t300 = -t192 * t305 - t218 * t271;
t223 = t261 * qJD(1);
t299 = t242 * t223 + t244 * t320;
t298 = t213 * t285 + t242 * t224;
t267 = t238 * t280;
t228 = t244 * t267;
t297 = -t230 * t271 + t228;
t294 = MDP(20) * t242;
t227 = qJD(1) * t233;
t289 = qJD(3) * t243;
t287 = qJD(4) * t218;
t284 = qJD(5) * t230;
t283 = t195 * qJD(4);
t176 = -t196 * t242 + t199 * t244;
t281 = qJD(5) - t176;
t277 = MDP(18) - MDP(21);
t276 = pkin(7) * t230 * t242;
t275 = pkin(7) * t309;
t274 = pkin(7) * t289;
t273 = pkin(7) * t282;
t272 = t244 * t197 + t199 * t285 + t242 * t212;
t270 = t230 * t285;
t263 = pkin(4) + t308;
t260 = pkin(4) * t244 + qJ(5) * t242;
t249 = t196 * t286 - t272;
t167 = qJ(5) * t266 - t249 - t284;
t258 = t167 * t244 + t168 * t242;
t172 = pkin(4) * t230 + t281;
t257 = t172 * t244 - t173 * t242;
t256 = t172 * t242 + t173 * t244;
t255 = t223 * t244 - t242 * t320;
t253 = 0.2e1 * qJD(3) * t227;
t251 = t232 + t259;
t248 = -t176 * t230 + t249;
t225 = -pkin(3) - t260;
t207 = t220 * t289;
t193 = t251 * t243;
t185 = pkin(4) * t220 + qJ(5) * t218;
t183 = t245 * t263 - t312;
t182 = -qJ(5) * t245 + t296;
t180 = -pkin(4) * t293 - t255;
t179 = qJ(5) * t293 + t299;
t178 = -t191 - t311;
t175 = pkin(4) * t218 - qJ(5) * t220 + t195;
t174 = (qJD(4) * t260 - qJD(5) * t244) * t243 + t251 * t288;
t171 = t296 * qJD(4) - t263 * t289 - t310;
t170 = (-t232 * t286 - qJD(5)) * t245 + (qJ(5) - t307) * t289 + t298;
t1 = [0.2e1 * t265 * t323 - 0.2e1 * t280 * t322 + MDP(7) * t302 - MDP(8) * t304 + (-t232 * t302 + t243 * t253) * MDP(10) + (t232 * t304 + t245 * t253) * MDP(11) + (-t191 * t305 + (-t269 + t271) * t220) * MDP(12) + (-t220 * t268 + (-t220 * t288 + (t191 + t287) * t243) * t242 + t300) * MDP(13) + (t191 * t245 + t230 * t269 + t207 + t297) * MDP(14) + (t230 * t268 + t192 * t245 + (-t218 * t243 + (-qJD(1) * t238 + t230 * t245) * t242) * qJD(3)) * MDP(15) + (-t230 - t292) * MDP(16) * t289 + (-(-t213 * t286 + t310) * t230 + (t232 * t270 + (t218 * t232 + t316) * qJD(3) + t262) * t245 + (t244 * t283 + t232 * t192 + t314 + (-t230 * t308 + (-t232 * t306 + t312) * qJD(1) + t176) * qJD(3)) * t243) * MDP(17) + (t298 * t230 + ((-t230 * t232 - t196) * t286 + (t220 * t232 + t315) * qJD(3) + t272) * t245 + (-t242 * t283 - t232 * t191 + t313 + (-qJD(1) * t296 - t230 * t307 - t177) * qJD(3)) * t243) * MDP(18) + (t171 * t230 + t174 * t218 + t192 * t193 + (t175 * t290 + t168) * t245 + (t175 * t285 + t319 + (-qJD(1) * t183 - t172) * qJD(3)) * t243) * MDP(19) + (-t170 * t218 + t171 * t220 - t182 * t192 - t183 * t191 + t257 * t288 + (-qJD(4) * t256 - t167 * t242 + t168 * t244) * t243) * MDP(20) + (-t170 * t230 - t174 * t220 + t191 * t193 + (-t175 * t282 - t167) * t245 + (t175 * t286 - t318 + (qJD(1) * t182 + t173) * qJD(3)) * t243) * MDP(21) + (t167 * t182 + t168 * t183 + t169 * t193 + t170 * t173 + t171 * t172 + t174 * t175) * MDP(22); (t207 - t228) * MDP(18) + t300 * MDP(20) + t297 * MDP(21) + (-t246 * MDP(11) - t169 * MDP(22) - t278 * t192 + t277 * t191 + (t220 * t294 + t256 * MDP(22) + (t244 * MDP(18) + t242 * t278) * t230) * qJD(3)) * t245 + (-t246 * MDP(10) - t191 * t294 - qJD(3) * t220 * MDP(21) + (qJD(3) * t175 + t258) * MDP(22) + ((t218 * t242 + t220 * t244) * MDP(20) + t257 * MDP(22) + (-t242 * t277 + t244 * t278) * t230) * qJD(4)) * t243 + t278 * (t218 * t289 - t242 * t267); (qJD(3) * t203 - t227 * t293 - t198) * MDP(10) - t227 * t292 * MDP(11) + (-t191 * t242 - t220 * t309) * MDP(12) + ((-t191 + t311) * t244 + (t220 * t230 - t192) * t242) * MDP(13) + (-t270 + (t230 * t303 + (-t220 + t290) * t243) * qJD(1)) * MDP(14) + (t230 * t286 + (-t230 * t306 + (t218 + t282) * t243) * qJD(1)) * MDP(15) + t230 * MDP(16) * t293 + (-pkin(3) * t192 - t313 + t255 * t230 - t203 * t218 + (t275 + t316) * qJD(4) + (-t176 * t243 + (-t195 * t245 - t274) * t242) * qJD(1)) * MDP(17) + (pkin(3) * t191 + t314 - t299 * t230 - t203 * t220 + (-t276 + t315) * qJD(4) + (-t195 * t303 + (t177 - t273) * t243) * qJD(1)) * MDP(18) + (-t318 - t180 * t230 + t192 * t225 - t301 * t218 + (t175 * t242 + t275) * qJD(4) + (t172 * t243 + (-t175 * t245 - t274) * t242) * qJD(1)) * MDP(19) + (t179 * t218 - t180 * t220 + (t167 - t230 * t172 + (qJD(4) * t220 - t192) * pkin(7)) * t244 + ((-t191 + t287) * pkin(7) + t321) * t242) * MDP(20) + (-t319 + t179 * t230 + t191 * t225 + t301 * t220 + (-t175 * t244 + t276) * qJD(4) + (t175 * t303 + (-t173 + t273) * t243) * qJD(1)) * MDP(21) + (t169 * t225 - t172 * t180 - t173 * t179 - t301 * t175 + (t257 * qJD(4) + t258) * pkin(7)) * MDP(22) + (-t245 * t323 + t322) * qJD(1) ^ 2; t178 * MDP(14) - MDP(15) * t264 + t248 * MDP(18) + (pkin(4) * t191 - qJ(5) * t192) * MDP(20) + (-t248 - 0.2e1 * t284) * MDP(21) + (-pkin(4) * t168 + qJ(5) * t167 - t172 * t177 + t173 * t281 - t175 * t185) * MDP(22) + (-t230 * MDP(15) - t195 * MDP(17) - t175 * MDP(19) + (t173 - t177) * MDP(20) + t185 * MDP(21) + MDP(13) * t220) * t220 + (-MDP(15) * t268 + (-MDP(15) * t306 + (0.2e1 * pkin(4) * MDP(19) + 0.2e1 * qJ(5) * MDP(21) + MDP(16)) * t243) * qJD(3)) * qJD(1) + (t220 * MDP(12) + t195 * MDP(18) - t185 * MDP(19) + (t172 - t281) * MDP(20) - t175 * MDP(21) - MDP(13) * t218) * t218 + t278 * (-t177 * t230 - t262); (t218 * t220 - t266) * MDP(19) + t178 * MDP(20) + (-t220 ^ 2 - t230 ^ 2) * MDP(21) + (t175 * t220 + t321) * MDP(22);];
tauc = t1;
