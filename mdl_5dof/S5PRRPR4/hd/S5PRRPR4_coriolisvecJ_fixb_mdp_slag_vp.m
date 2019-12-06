% Calculate minimal parameter regressor of Coriolis joint torque vector for
% S5PRRPR4
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% qJD [5x1]
%   Generalized joint velocities
% pkin [9x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d2,d3,d5,theta1,theta4]';
% MDP [20x1]
%   Minimal dynamic parameter vector (fixed base model)
%   see S5PRRPR4_convert_par2_MPV_fixb.m
% 
% Output:
% tauc [5x1]
%   joint torques required to compensate Coriolis and centrifugal load

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-05 16:24
% Revision: 77da58f92bca3eff71542919beafa37024070d86 (2019-12-05)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauc = S5PRRPR4_coriolisvecJ_fixb_mdp_slag_vp(qJ, qJD, pkin, MDP)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(9,1),zeros(20,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5PRRPR4_coriolisvecJ_fixb_mdp_slag_vp: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5PRRPR4_coriolisvecJ_fixb_mdp_slag_vp: qJD has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S5PRRPR4_coriolisvecJ_fixb_mdp_slag_vp: pkin has to be [9x1] (double)');
assert(isreal(MDP) && all(size(MDP) == [20 1]), ...
  'S5PRRPR4_coriolisvecJ_fixb_mdp_slag_vp: MDP has to be [20x1] (double)'); 

%% Symbolic Calculation
% From coriolisvec_joint_fixb_mdp_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-05 16:23:27
% EndTime: 2019-12-05 16:23:33
% DurationCPUTime: 1.46s
% Computational Cost: add. (1043->192), mult. (2723->286), div. (0->0), fcn. (1998->8), ass. (0->104)
t253 = sin(pkin(9));
t254 = cos(pkin(9));
t256 = sin(qJ(3));
t259 = cos(qJ(3));
t266 = t253 * t256 - t254 * t259;
t230 = t266 * qJD(2);
t258 = cos(qJ(5));
t219 = t258 * t230;
t238 = t253 * t259 + t254 * t256;
t231 = t238 * qJD(3);
t222 = qJD(2) * t231;
t283 = qJD(2) * qJD(3);
t279 = t259 * t283;
t280 = t256 * t283;
t223 = -t253 * t280 + t254 * t279;
t232 = t238 * qJD(2);
t255 = sin(qJ(5));
t288 = qJD(5) * t255;
t161 = -qJD(5) * t219 - t255 * t222 + t258 * t223 - t232 * t288;
t268 = -t230 * t255 + t258 * t232;
t162 = qJD(5) * t268 + t258 * t222 + t223 * t255;
t190 = -t232 * t255 - t219;
t250 = qJD(3) + qJD(5);
t300 = t190 * t250;
t301 = t268 * t250;
t317 = (-t162 + t301) * MDP(17) + (-t190 ^ 2 + t268 ^ 2) * MDP(15) - t190 * t268 * MDP(14) + (t161 - t300) * MDP(16);
t316 = -0.2e1 * t283;
t257 = sin(qJ(2));
t285 = t257 * qJD(1);
t245 = qJD(2) * pkin(6) + t285;
t275 = qJ(4) * qJD(2) + t245;
t226 = t275 * t256;
t211 = qJD(3) * pkin(3) - t226;
t227 = t275 * t259;
t299 = t254 * t227;
t175 = t253 * t211 + t299;
t305 = pkin(7) * t230;
t169 = t175 - t305;
t281 = -pkin(3) * t259 - pkin(2);
t260 = cos(qJ(2));
t292 = qJD(1) * t260;
t235 = qJD(2) * t281 + qJD(4) - t292;
t197 = pkin(4) * t230 + t235;
t315 = t169 * t288 - t197 * t190;
t312 = MDP(5) * t256;
t311 = MDP(6) * (t256 ^ 2 - t259 ^ 2);
t303 = -qJ(4) - pkin(6);
t278 = qJD(3) * t303;
t228 = qJD(4) * t259 + t256 * t278;
t229 = -qJD(4) * t256 + t259 * t278;
t296 = -t228 * t253 + t254 * t229 + t238 * t292;
t295 = t254 * t228 + t253 * t229 + t266 * t292;
t290 = qJD(3) * t256;
t310 = pkin(3) * t290 - t285;
t309 = t256 * MDP(10) + t259 * MDP(11);
t234 = t266 * qJD(3);
t308 = qJD(5) - t250;
t272 = qJD(4) + t292;
t192 = -t245 * t290 + (-qJ(4) * t290 + t259 * t272) * qJD(2);
t289 = qJD(3) * t259;
t193 = -t245 * t289 + (-qJ(4) * t289 - t256 * t272) * qJD(2);
t163 = -t192 * t253 + t254 * t193;
t159 = -pkin(7) * t223 + t163;
t164 = t254 * t192 + t253 * t193;
t160 = -pkin(7) * t222 + t164;
t307 = t258 * t159 - t255 * t160 - t197 * t268;
t306 = pkin(3) * t253;
t304 = pkin(7) * t232;
t302 = qJD(2) * pkin(2);
t207 = t253 * t227;
t261 = qJD(3) ^ 2;
t298 = t256 * t261;
t297 = t259 * t261;
t177 = -t254 * t226 - t207;
t242 = t303 * t256;
t243 = t303 * t259;
t199 = t253 * t242 - t254 * t243;
t236 = pkin(3) * t280 + qJD(2) * t285;
t291 = qJD(2) * t256;
t174 = t254 * t211 - t207;
t176 = t226 * t253 - t299;
t198 = t254 * t242 + t243 * t253;
t274 = pkin(4) * t231 + t310;
t271 = -pkin(7) * t234 + qJD(5) * (-pkin(7) * t266 + t199) - t296;
t270 = pkin(7) * t231 - qJD(5) * (-pkin(7) * t238 + t198) - t295;
t168 = qJD(3) * pkin(4) + t174 - t304;
t269 = -t255 * t168 - t258 * t169;
t267 = -t238 * t255 - t258 * t266;
t196 = t238 * t258 - t255 * t266;
t264 = -0.2e1 * qJD(3) * t302;
t262 = qJD(2) ^ 2;
t248 = pkin(3) * t254 + pkin(4);
t225 = t266 * t257;
t224 = t238 * t257;
t215 = pkin(4) * t266 + t281;
t202 = pkin(3) * t291 + pkin(4) * t232;
t194 = pkin(4) * t222 + t236;
t179 = -t230 * t260 - t231 * t257;
t178 = -t232 * t260 + t257 * t234;
t171 = t177 - t304;
t170 = t176 + t305;
t166 = qJD(5) * t196 + t258 * t231 - t234 * t255;
t165 = qJD(5) * t267 - t231 * t255 - t234 * t258;
t1 = [(-t178 * t232 - t179 * t230 + t222 * t225 + t223 * t224) * MDP(12) + (-t163 * t224 - t164 * t225 + t174 * t178 + t175 * t179) * MDP(13) + ((t178 * t258 - t179 * t255) * MDP(19) - (t178 * t255 + t179 * t258) * MDP(20) + ((t224 * t255 + t225 * t258) * MDP(19) - (-t224 * t258 + t225 * t255) * MDP(20)) * qJD(5)) * t250 + (-t236 * MDP(13) - t162 * MDP(19) - t161 * MDP(20) - t262 * MDP(4) + t309 * t316) * t260 + (-t262 * MDP(3) + (MDP(13) * t235 - MDP(19) * t190 + MDP(20) * t268) * qJD(2) + (-MDP(10) * t259 + MDP(11) * t256) * (t261 + t262)) * t257; 0.2e1 * t279 * t312 + t311 * t316 + MDP(7) * t297 - MDP(8) * t298 + (-pkin(6) * t297 + t256 * t264) * MDP(10) + (pkin(6) * t298 + t259 * t264) * MDP(11) + (-t163 * t238 - t164 * t266 + t174 * t234 - t175 * t231 - t198 * t223 - t199 * t222 - t230 * t295 - t232 * t296) * MDP(12) + (t163 * t198 + t164 * t199 + t296 * t174 + t295 * t175 + t310 * t235 + t236 * t281) * MDP(13) + (t161 * t196 + t165 * t268) * MDP(14) + (t161 * t267 - t162 * t196 + t165 * t190 - t166 * t268) * MDP(15) + (t215 * t162 + t197 * t166 - t190 * t274 - t194 * t267) * MDP(19) + (t215 * t161 + t197 * t165 + t194 * t196 + t268 * t274) * MDP(20) + (t165 * MDP(16) - t166 * MDP(17) + (t255 * t270 - t258 * t271) * MDP(19) + (t255 * t271 + t258 * t270) * MDP(20)) * t250; ((t175 + t176) * t232 - (t174 - t177) * t230 + (-t222 * t253 - t223 * t254) * pkin(3)) * MDP(12) + (-t174 * t176 - t175 * t177 + (t163 * t254 + t164 * t253 - t235 * t291) * pkin(3)) * MDP(13) + (-(t170 * t258 - t171 * t255) * t250 + t202 * t190 + ((-t248 * t255 - t258 * t306) * t250 + t269) * qJD(5) + t307) * MDP(19) + (-t258 * t160 - t255 * t159 + (t170 * t255 + t171 * t258) * t250 - t202 * t268 + (-(t248 * t258 - t255 * t306) * t250 - t258 * t168) * qJD(5) + t315) * MDP(20) + t309 * qJD(2) * t302 + (-t259 * t312 + t311) * t262 + t317; (-t230 ^ 2 - t232 ^ 2) * MDP(12) + (t174 * t232 + t175 * t230 + t236) * MDP(13) + (t162 + t301) * MDP(19) + (t161 + t300) * MDP(20); (t308 * t269 + t307) * MDP(19) + ((-t169 * t250 - t159) * t255 + (-t308 * t168 - t160) * t258 + t315) * MDP(20) + t317;];
tauc = t1;
