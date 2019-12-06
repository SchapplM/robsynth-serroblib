% Calculate minimal parameter regressor of Coriolis joint torque vector for
% S5PRRPR5
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% qJD [5x1]
%   Generalized joint velocities
% pkin [10x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,alpha2,d2,d3,d5,theta1,theta4]';
% MDP [20x1]
%   Minimal dynamic parameter vector (fixed base model)
%   see S5PRRPR5_convert_par2_MPV_fixb.m
% 
% Output:
% tauc [5x1]
%   joint torques required to compensate Coriolis and centrifugal load

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-05 16:28
% Revision: 77da58f92bca3eff71542919beafa37024070d86 (2019-12-05)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauc = S5PRRPR5_coriolisvecJ_fixb_mdp_slag_vp(qJ, qJD, pkin, MDP)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(10,1),zeros(20,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5PRRPR5_coriolisvecJ_fixb_mdp_slag_vp: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5PRRPR5_coriolisvecJ_fixb_mdp_slag_vp: qJD has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S5PRRPR5_coriolisvecJ_fixb_mdp_slag_vp: pkin has to be [10x1] (double)');
assert(isreal(MDP) && all(size(MDP) == [20 1]), ...
  'S5PRRPR5_coriolisvecJ_fixb_mdp_slag_vp: MDP has to be [20x1] (double)'); 

%% Symbolic Calculation
% From coriolisvec_joint_fixb_mdp_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-05 16:28:00
% EndTime: 2019-12-05 16:28:09
% DurationCPUTime: 1.82s
% Computational Cost: add. (1330->235), mult. (3572->354), div. (0->0), fcn. (2714->10), ass. (0->120)
t241 = sin(qJ(3));
t316 = MDP(5) * t241;
t244 = cos(qJ(3));
t315 = (t241 ^ 2 - t244 ^ 2) * MDP(6);
t242 = sin(qJ(2));
t237 = sin(pkin(5));
t288 = qJD(1) * t237;
t274 = t242 * t288;
t221 = qJD(2) * pkin(7) + t274;
t245 = cos(qJ(2));
t273 = t245 * t288;
t260 = qJD(4) + t273;
t239 = cos(pkin(5));
t287 = qJD(1) * t239;
t272 = t241 * t287;
t314 = (-t221 * t244 - t272) * qJD(3) + (-qJD(3) * t244 * qJ(4) - t260 * t241) * qJD(2);
t284 = qJD(3) * t241;
t313 = pkin(3) * t284 - t274;
t312 = t241 * MDP(10) + t244 * MDP(11);
t227 = t244 * t287;
t175 = (-t221 * t241 + t227) * qJD(3) + (-qJ(4) * t284 + t260 * t244) * qJD(2);
t236 = sin(pkin(10));
t238 = cos(pkin(10));
t161 = t175 * t236 - t238 * t314;
t286 = qJD(2) * t241;
t296 = t238 * t244;
t212 = qJD(2) * t296 - t236 * t286;
t209 = qJD(5) - t212;
t219 = t236 * t244 + t238 * t241;
t214 = t219 * qJD(2);
t230 = pkin(3) * t236 + pkin(8);
t311 = (pkin(3) * t286 + pkin(4) * t214 - pkin(8) * t212 + qJD(5) * t230) * t209 + t161;
t162 = t238 * t175 + t314 * t236;
t265 = qJ(4) * qJD(2) + t221;
t193 = -t265 * t241 + t227;
t190 = qJD(3) * pkin(3) + t193;
t194 = t265 * t244 + t272;
t304 = t194 * t236;
t165 = t190 * t238 - t304;
t163 = -qJD(3) * pkin(4) - t165;
t275 = -pkin(3) * t244 - pkin(2);
t208 = t275 * qJD(2) + qJD(4) - t273;
t174 = -pkin(4) * t212 - pkin(8) * t214 + t208;
t218 = t236 * t241 - t296;
t186 = pkin(4) * t218 - pkin(8) * t219 + t275;
t216 = t218 * qJD(3);
t309 = -qJ(4) - pkin(7);
t223 = t309 * t244;
t270 = t309 * t241;
t198 = -t238 * t223 + t236 * t270;
t213 = t219 * qJD(3);
t206 = qJD(2) * t213;
t259 = t161 * t219 - t198 * t206;
t267 = qJD(3) * t309;
t210 = qJD(4) * t244 + t241 * t267;
t253 = -qJD(4) * t241 + t244 * t267;
t291 = t238 * t210 + t218 * t273 + t236 * t253;
t310 = -(qJD(5) * t174 + t162) * t218 - t163 * t216 + (-qJD(5) * t186 - t291) * t209 + t259;
t308 = qJD(2) * pkin(2);
t307 = t163 * t219;
t240 = sin(qJ(5));
t283 = qJD(5) * t240;
t277 = qJD(2) * qJD(3);
t268 = t244 * t277;
t269 = t241 * t277;
t207 = -t236 * t269 + t238 * t268;
t243 = cos(qJ(5));
t279 = t243 * qJD(3);
t290 = qJD(5) * t279 + t243 * t207;
t176 = -t214 * t283 + t290;
t306 = t176 * t240;
t305 = t186 * t206;
t298 = t214 * t240;
t199 = -t279 + t298;
t303 = t199 * t209;
t302 = t199 * t214;
t201 = qJD(3) * t240 + t214 * t243;
t301 = t201 * t209;
t300 = t201 * t214;
t299 = t207 * t240;
t297 = t237 * t242;
t188 = t238 * t194;
t295 = t240 * t206;
t246 = qJD(3) ^ 2;
t294 = t241 * t246;
t202 = t243 * t206;
t293 = t244 * t246;
t292 = t210 * t236 - t219 * t273 - t238 * t253;
t166 = t236 * t190 + t188;
t211 = pkin(3) * t269 + qJD(2) * t274;
t285 = qJD(2) * t242;
t282 = qJD(5) * t243;
t281 = qJD(5) * t245;
t271 = qJD(2) * t237 * t245;
t266 = t209 * t243;
t262 = pkin(4) * t213 + pkin(8) * t216 + t313;
t164 = qJD(3) * pkin(8) + t166;
t160 = t164 * t243 + t174 * t240;
t258 = t164 * t240 - t174 * t243;
t257 = t202 + (t212 * t240 - t283) * t209;
t217 = t239 * t241 + t244 * t297;
t256 = t239 * t244 - t241 * t297;
t255 = -t243 * t216 - t219 * t283;
t251 = -0.2e1 * qJD(3) * t308;
t170 = t193 * t238 - t304;
t250 = -t230 * t206 + (t163 + t170) * t209;
t247 = qJD(2) ^ 2;
t231 = -pkin(3) * t238 - pkin(4);
t197 = -t223 * t236 - t238 * t270;
t192 = -t217 * qJD(3) - t241 * t271;
t191 = t256 * qJD(3) + t244 * t271;
t183 = t238 * t217 + t236 * t256;
t182 = t217 * t236 - t238 * t256;
t177 = t201 * qJD(5) + t299;
t172 = pkin(4) * t206 - pkin(8) * t207 + t211;
t171 = t243 * t172;
t169 = t191 * t238 + t192 * t236;
t168 = t193 * t236 + t188;
t167 = t191 * t236 - t238 * t192;
t1 = [(t167 * t214 + t169 * t212 + t182 * t207 - t183 * t206) * MDP(12) + (t161 * t182 + t162 * t183 - t165 * t167 + t166 * t169) * MDP(13) + ((-t169 * t240 - t183 * t282) * t209 - t183 * t295 + t167 * t199 + t182 * t177) * MDP(19) + (-(t169 * t243 - t183 * t283) * t209 - t183 * t202 + t167 * t201 + t182 * t176) * MDP(20) + (MDP(10) * t192 - MDP(11) * t191) * qJD(3) + ((t208 * t285 - t211 * t245) * MDP(13) + ((t240 * t281 + t243 * t285) * t209 - t245 * t202) * MDP(19) + (-(t240 * t285 - t243 * t281) * t209 + t245 * t295) * MDP(20) - t312 * t245 * t277 + (-t245 * MDP(4) + (-MDP(10) * t244 + MDP(11) * t241 - MDP(3)) * t242) * t247) * t237; 0.2e1 * t268 * t316 - 0.2e1 * t277 * t315 + MDP(7) * t293 - MDP(8) * t294 + (-pkin(7) * t293 + t241 * t251) * MDP(10) + (pkin(7) * t294 + t244 * t251) * MDP(11) + (-t162 * t218 + t165 * t216 - t166 * t213 + t197 * t207 + t291 * t212 + t292 * t214 + t259) * MDP(12) + (t161 * t197 + t162 * t198 - t292 * t165 + t291 * t166 + t313 * t208 + t211 * t275) * MDP(13) + (t176 * t219 * t243 + t255 * t201) * MDP(14) + (-(-t199 * t243 - t201 * t240) * t216 + (-t306 - t177 * t243 + (t199 * t240 - t201 * t243) * qJD(5)) * t219) * MDP(15) + (t176 * t218 + t201 * t213 + t219 * t202 + t255 * t209) * MDP(16) + (-t219 * t295 - t177 * t218 - t199 * t213 + (t240 * t216 - t219 * t282) * t209) * MDP(17) + (t206 * t218 + t209 * t213) * MDP(18) + (-t258 * t213 + t171 * t218 + t197 * t177 + t292 * t199 + (t305 + t262 * t209 + (-t164 * t218 - t198 * t209 + t307) * qJD(5)) * t243 + t310 * t240) * MDP(19) + (-t160 * t213 + t197 * t176 + t292 * t201 + (-t305 - (-qJD(5) * t164 + t172) * t218 - qJD(5) * t307 + (qJD(5) * t198 - t262) * t209) * t240 + t310 * t243) * MDP(20); ((t166 - t168) * t214 + (t165 - t170) * t212 + (-t206 * t236 - t207 * t238) * pkin(3)) * MDP(12) + (t165 * t168 - t166 * t170 + (-t161 * t238 + t162 * t236 - t208 * t286) * pkin(3)) * MDP(13) + (t201 * t266 + t306) * MDP(14) + ((t176 - t303) * t243 + (-t177 - t301) * t240) * MDP(15) + (t209 * t266 + t295 - t300) * MDP(16) + (t257 + t302) * MDP(17) - t209 * t214 * MDP(18) + (-t168 * t199 + t231 * t177 + t214 * t258 + t250 * t240 - t311 * t243) * MDP(19) + (t160 * t214 - t168 * t201 + t231 * t176 + t311 * t240 + t250 * t243) * MDP(20) + t312 * t308 * qJD(2) + (-t244 * t316 + t315) * t247; (-t212 ^ 2 - t214 ^ 2) * MDP(12) + (t165 * t214 - t166 * t212 + t211) * MDP(13) + (t257 - t302) * MDP(19) + (-t209 ^ 2 * t243 - t295 - t300) * MDP(20); t201 * t199 * MDP(14) + (-t199 ^ 2 + t201 ^ 2) * MDP(15) + (t290 + t303) * MDP(16) + (-t299 + t301) * MDP(17) + t206 * MDP(18) + (t160 * t209 - t162 * t240 - t163 * t201 + t171) * MDP(19) + (-t162 * t243 + t163 * t199 - t240 * t172 - t209 * t258) * MDP(20) + (-MDP(16) * t298 - t201 * MDP(17) - t160 * MDP(19) + t258 * MDP(20)) * qJD(5);];
tauc = t1;
