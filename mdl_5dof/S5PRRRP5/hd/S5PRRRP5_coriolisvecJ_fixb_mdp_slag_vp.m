% Calculate Coriolis joint torque vector for
% S5PRRRP5
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% qJD [5x1]
%   Generalized joint velocities
% pkin [8x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d2,d3,d4,theta1]';
% MDP [22x1]
%   Minimal dynamic parameter vector (fixed base model)
%   see S5PRRRP5_convert_par2_MPV_fixb.m
% 
% Output:
% tauc [5x1]
%   joint torques required to compensate Coriolis and centrifugal load

% Quelle: HybrDyn-Toolbox
% Datum: 2021-01-15 16:34
% Revision: 24b2e7d74a0c1a3b64fa2f8f5ad758691ad61af3 (2021-01-15)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauc = S5PRRRP5_coriolisvecJ_fixb_mdp_slag_vp(qJ, qJD, pkin, MDP)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(8,1),zeros(22,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5PRRRP5_coriolisvecJ_fixb_mdp_slag_vp: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5PRRRP5_coriolisvecJ_fixb_mdp_slag_vp: qJD has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S5PRRRP5_coriolisvecJ_fixb_mdp_slag_vp: pkin has to be [8x1] (double)');
assert(isreal(MDP) && all(size(MDP) == [22 1]), ...
  'S5PRRRP5_coriolisvecJ_fixb_mdp_slag_vp: MDP has to be [22x1] (double)'); 

%% Symbolic Calculation
% From coriolisvec_joint_fixb_mdp_matlab.m
% OptimizationMode: 2
% StartTime: 2021-01-15 16:33:48
% EndTime: 2021-01-15 16:33:53
% DurationCPUTime: 1.54s
% Computational Cost: add. (1296->212), mult. (3237->289), div. (0->0), fcn. (2190->6), ass. (0->116)
t237 = cos(qJ(3));
t235 = sin(qJ(2));
t285 = qJD(1) * t235;
t220 = qJD(2) * pkin(6) + t285;
t267 = pkin(7) * qJD(2) + t220;
t196 = t267 * t237;
t233 = sin(qJ(4));
t188 = t233 * t196;
t234 = sin(qJ(3));
t195 = t267 * t234;
t191 = qJD(3) * pkin(3) - t195;
t236 = cos(qJ(4));
t264 = t191 * t236 - t188;
t296 = t234 * t236;
t212 = t233 * t237 + t296;
t203 = t212 * qJD(2);
t300 = t203 * qJ(5);
t164 = -t300 + t264;
t230 = qJD(3) + qJD(4);
t315 = t230 * t212;
t178 = t315 * qJD(2);
t274 = qJD(2) * qJD(3);
t316 = -0.2e1 * t274;
t281 = qJD(3) * t234;
t238 = cos(qJ(2));
t284 = qJD(1) * t238;
t183 = -t220 * t281 + (-pkin(7) * t281 + t237 * t284) * qJD(2);
t280 = qJD(3) * t237;
t184 = -t220 * t280 + (-pkin(7) * t280 - t234 * t284) * qJD(2);
t279 = qJD(4) * t233;
t248 = -(qJD(4) * t191 + t183) * t236 - t233 * t184 + t196 * t279;
t314 = -qJ(5) * t178 - t248;
t313 = MDP(5) * t234;
t312 = MDP(6) * (t234 ^ 2 - t237 ^ 2);
t294 = t236 * t237;
t297 = t233 * t234;
t211 = -t294 + t297;
t198 = t211 * t235;
t307 = pkin(6) + pkin(7);
t271 = qJD(3) * t307;
t213 = t234 * t271;
t214 = t237 * t271;
t215 = t307 * t234;
t216 = t307 * t237;
t255 = t215 * t233 - t216 * t236;
t311 = qJD(4) * t255 + t212 * t284 + t213 * t233 - t214 * t236;
t252 = t211 * t238;
t278 = qJD(4) * t236;
t310 = -qJD(1) * t252 + t213 * t236 + t214 * t233 + t215 * t278 + t216 * t279;
t272 = pkin(3) * t281;
t253 = t272 - t285;
t309 = t234 * MDP(10) + t237 * MDP(11);
t308 = t203 ^ 2;
t306 = pkin(3) * t230;
t305 = qJD(2) * pkin(2);
t257 = t230 * t297;
t269 = t237 * t274;
t270 = qJD(2) * t294;
t288 = qJD(4) * t270 + t236 * t269;
t177 = qJD(2) * t257 - t288;
t304 = qJ(5) * t177;
t283 = qJD(2) * t234;
t201 = t233 * t283 - t270;
t302 = qJ(5) * t201;
t301 = t201 * t230;
t299 = t203 * t230;
t229 = -pkin(3) * t237 - pkin(2);
t208 = qJD(2) * t229 - t284;
t298 = t208 * t203;
t239 = qJD(3) ^ 2;
t295 = t234 * t239;
t190 = t236 * t196;
t293 = t237 * t239;
t292 = -qJ(5) * t315 - qJD(5) * t211 - t310;
t181 = -t236 * t280 - t237 * t278 + t257;
t291 = qJ(5) * t181 - qJD(5) * t212 + t311;
t163 = pkin(4) * t230 + t164;
t290 = t163 - t164;
t289 = -t195 * t236 - t188;
t282 = qJD(2) * t235;
t209 = qJD(1) * t282 + qJD(2) * t272;
t268 = pkin(4) * t201 + qJD(5);
t180 = t208 + t268;
t275 = qJD(5) + t180;
t273 = pkin(3) * t283;
t266 = -t233 * t183 + t184 * t236;
t263 = t195 * t233 - t190;
t262 = pkin(4) * t315 + t253;
t259 = t230 * t234;
t258 = MDP(19) * t230;
t172 = pkin(4) * t178 + t209;
t256 = -t191 * t233 - t190;
t200 = t201 ^ 2;
t251 = t203 * t201 * MDP(12) + (-qJD(2) * t233 * t259 + t288 + t301) * MDP(14) + (t299 - t178) * MDP(15) + (-t200 + t308) * MDP(13);
t249 = -0.2e1 * qJD(3) * t305;
t247 = qJD(4) * t256 + t266;
t245 = t208 * t201 + t248;
t244 = t247 + t304;
t243 = (-t190 + (-t191 - t306) * t233) * qJD(4) + t266;
t241 = t201 * t275 - t314;
t240 = qJD(2) ^ 2;
t228 = pkin(3) * t236 + pkin(4);
t217 = t278 * t306;
t197 = t212 * t235;
t194 = pkin(4) * t211 + t229;
t186 = pkin(4) * t203 + t273;
t175 = -qJ(5) * t211 - t255;
t174 = -qJ(5) * t212 - t215 * t236 - t216 * t233;
t169 = t198 * t230 - t238 * t203;
t168 = -qJD(2) * t252 - t235 * t315;
t167 = -t300 + t289;
t166 = t263 + t302;
t165 = -t256 - t302;
t158 = -qJD(5) * t203 + t244;
t157 = -qJD(5) * t201 + t314;
t1 = [(-t168 * t201 - t169 * t203 - t177 * t197 + t178 * t198) * MDP(21) + (-t157 * t198 - t158 * t197 + t163 * t169 + t165 * t168) * MDP(22) + (MDP(17) + MDP(19)) * (t169 * t230 - t178 * t238 + t201 * t282) + (MDP(18) + MDP(20)) * (-t168 * t230 + t177 * t238 + t203 * t282) + (-t172 * MDP(22) - t240 * MDP(4) + t309 * t316) * t238 + (qJD(2) * t180 * MDP(22) - t240 * MDP(3) + (-MDP(10) * t237 + MDP(11) * t234) * (t239 + t240)) * t235; 0.2e1 * t269 * t313 + t312 * t316 + MDP(7) * t293 - MDP(8) * t295 + (-pkin(6) * t293 + t234 * t249) * MDP(10) + (pkin(6) * t295 + t237 * t249) * MDP(11) + (-t177 * t212 - t181 * t203) * MDP(12) + (t177 * t211 - t178 * t212 + t181 * t201 - t203 * t315) * MDP(13) + (t229 * t178 + t253 * t201 + t208 * t315 + t209 * t211) * MDP(17) + (-t229 * t177 - t208 * t181 + t253 * t203 + t209 * t212) * MDP(18) + (t172 * t211 + t178 * t194 + t180 * t315 + t201 * t262) * MDP(19) + (t172 * t212 - t177 * t194 - t180 * t181 + t203 * t262) * MDP(20) + (-t157 * t211 - t158 * t212 + t163 * t181 - t165 * t315 + t174 * t177 - t175 * t178 - t201 * t292 - t203 * t291) * MDP(21) + (t157 * t175 + t158 * t174 + t163 * t291 + t165 * t292 + t172 * t194 + t180 * t262) * MDP(22) + (-t181 * MDP(14) - MDP(15) * t315 + MDP(17) * t311 + MDP(18) * t310 + t291 * MDP(19) - t292 * MDP(20)) * t230; (-t201 * t273 - t230 * t263 + t243 - t298) * MDP(17) + (-t203 * t273 + t230 * t289 - t217 + t245) * MDP(18) + (-t166 * t230 - t186 * t201 - t203 * t275 + t243 + t304) * MDP(19) + (t167 * t230 - t186 * t203 - t217 + t241) * MDP(20) + (t177 * t228 + (t165 + t166) * t203 + (-t163 + t167) * t201 + (-t178 * t233 + (-t201 * t236 + t203 * t233) * qJD(4)) * pkin(3)) * MDP(21) + (t158 * t228 - t163 * t166 - t165 * t167 - t180 * t186 + (t157 * t233 + (-t163 * t233 + t165 * t236) * qJD(4)) * pkin(3)) * MDP(22) + t251 + t309 * qJD(2) * t305 + (-t237 * t313 + t312) * t240; (-t230 * t256 + t247 - t298) * MDP(17) + (t230 * t264 + t245) * MDP(18) + (t165 * t230 + (-t180 - t268) * t203 + t244) * MDP(19) + (-pkin(4) * t308 + t164 * t230 + t241) * MDP(20) + (pkin(4) * t177 - t201 * t290) * MDP(21) + (t290 * t165 + (-t180 * t203 + t158) * pkin(4)) * MDP(22) + t251; MDP(19) * t299 + (t288 - t301) * MDP(20) + (-t200 - t308) * MDP(21) + (t163 * t203 + t165 * t201 + t172) * MDP(22) + (t258 * t296 + (-MDP(20) * t259 + t237 * t258) * t233) * qJD(2);];
tauc = t1;
