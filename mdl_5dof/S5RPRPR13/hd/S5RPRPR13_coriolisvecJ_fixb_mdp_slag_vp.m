% Calculate minimal parameter regressor of Coriolis joint torque vector for
% S5RPRPR13
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
%   see S5RPRPR13_convert_par2_MPV_fixb.m
% 
% Output:
% tauc [5x1]
%   joint torques required to compensate Coriolis and centrifugal load

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-31 18:33
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauc = S5RPRPR13_coriolisvecJ_fixb_mdp_slag_vp(qJ, qJD, pkin, MDP)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(8,1),zeros(25,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RPRPR13_coriolisvecJ_fixb_mdp_slag_vp: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RPRPR13_coriolisvecJ_fixb_mdp_slag_vp: qJD has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S5RPRPR13_coriolisvecJ_fixb_mdp_slag_vp: pkin has to be [8x1] (double)');
assert(isreal(MDP) && all(size(MDP) == [25 1]), ...
  'S5RPRPR13_coriolisvecJ_fixb_mdp_slag_vp: MDP has to be [25x1] (double)'); 

%% Symbolic Calculation
% From coriolisvec_joint_fixb_mdp_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 18:32:57
% EndTime: 2019-12-31 18:33:01
% DurationCPUTime: 2.35s
% Computational Cost: add. (1291->254), mult. (3389->327), div. (0->0), fcn. (2442->6), ass. (0->113)
t239 = cos(pkin(8));
t302 = cos(qJ(3));
t272 = t302 * t239;
t261 = qJD(1) * t272;
t238 = sin(pkin(8));
t241 = sin(qJ(3));
t287 = t238 * t241;
t271 = qJD(1) * t287;
t211 = -t261 + t271;
t242 = cos(qJ(5));
t240 = sin(qJ(5));
t280 = qJD(3) * t240;
t191 = -t242 * t211 + t280;
t220 = t238 * t302 + t241 * t239;
t309 = t220 * qJD(1);
t311 = qJD(5) + t309;
t313 = t191 * t311;
t193 = qJD(3) * t242 + t211 * t240;
t263 = t311 * t193;
t230 = qJD(3) * t261;
t279 = qJD(3) * t241;
t270 = t238 * t279;
t199 = qJD(1) * t270 - t230;
t195 = t242 * t199;
t265 = t240 * t311;
t312 = -t265 * t311 - t195;
t299 = pkin(6) + qJ(2);
t224 = t299 * t238;
t221 = qJD(1) * t224;
t225 = t299 * t239;
t222 = qJD(1) * t225;
t285 = -t302 * t221 - t241 * t222;
t308 = qJD(4) - t285;
t307 = (t238 ^ 2 + t239 ^ 2) * (MDP(7) * qJ(2) + MDP(6));
t187 = -t241 * t221 + t302 * t222;
t183 = -qJD(3) * qJ(4) - t187;
t300 = pkin(4) * t211;
t169 = -t183 - t300;
t303 = pkin(3) + pkin(7);
t306 = t303 * t199 + (t169 - t187 + t300) * t311;
t268 = qJD(3) * t302;
t269 = qJD(2) * t302;
t178 = (qJD(2) * t238 + qJD(3) * t225) * t241 + t224 * t268 - t239 * t269;
t305 = t211 ^ 2;
t304 = t309 ^ 2;
t216 = t220 * qJD(3);
t200 = qJD(1) * t216;
t301 = pkin(3) * t200;
t297 = qJ(4) * t199;
t296 = qJ(4) * t211;
t219 = -t272 + t287;
t235 = -pkin(2) * t239 - pkin(1);
t253 = -qJ(4) * t220 + t235;
t172 = t219 * t303 + t253;
t295 = t172 * t199;
t278 = qJD(5) * t240;
t277 = qJD(5) * t242;
t286 = t240 * t200 + t211 * t277;
t175 = -qJD(3) * t278 + t286;
t294 = t175 * t242;
t223 = t235 * qJD(1) + qJD(2);
t245 = -qJ(4) * t309 + t223;
t177 = pkin(3) * t211 + t245;
t293 = t177 * t309;
t292 = t191 * t211;
t291 = t193 * t211;
t290 = t199 * t240;
t289 = t211 * t309;
t288 = t219 * t240;
t282 = qJD(3) * t178;
t252 = -t241 * t224 + t225 * t302;
t179 = t220 * qJD(2) + t252 * qJD(3);
t281 = qJD(3) * t179;
t275 = pkin(4) * t309 + t308;
t273 = qJD(1) * qJD(2);
t267 = t241 * t273;
t264 = t242 * t311;
t260 = qJD(1) * t269;
t262 = t221 * t268 + t222 * t279 + t238 * t267 - t239 * t260;
t168 = -t221 * t279 + t222 * t268 + t238 * t260 + t239 * t267;
t189 = t224 * t302 + t241 * t225;
t162 = t211 * t303 + t245;
t166 = -qJD(3) * t303 + t275;
t156 = t162 * t242 + t166 * t240;
t259 = t162 * t240 - t166 * t242;
t257 = -qJD(4) * t309 + t297;
t215 = -t239 * t268 + t270;
t256 = qJ(4) * t215 - qJD(4) * t220;
t251 = t216 * t240 + t219 * t277;
t250 = qJD(3) * t187 - t168;
t167 = -qJD(3) * qJD(4) + t262;
t158 = -pkin(4) * t200 - t167;
t249 = t158 + (t303 * t311 + t296) * t311;
t180 = t220 * pkin(4) + t189;
t248 = t158 * t219 + t169 * t216 + t180 * t199;
t247 = -t264 * t311 + t290;
t201 = qJD(3) * t211;
t196 = t242 * t200;
t188 = t199 * t220;
t185 = pkin(3) * t219 + t253;
t184 = pkin(3) * t309 + t296;
t182 = -qJD(3) * pkin(3) + t308;
t181 = -t219 * pkin(4) + t252;
t176 = t193 * qJD(5) - t196;
t171 = pkin(3) * t216 + t256;
t165 = t257 + t301;
t164 = -t215 * pkin(4) + t179;
t163 = -pkin(4) * t216 - t178;
t161 = t216 * t303 + t256;
t160 = -pkin(4) * t199 + t168;
t159 = t242 * t160;
t157 = t200 * t303 + t257;
t1 = [(-t215 * t309 - t188) * MDP(8) + (t199 * t219 - t200 * t220 + t211 * t215 - t216 * t309) * MDP(9) + (t200 * t235 + t216 * t223 - t281) * MDP(13) + (-t199 * t235 - t215 * t223 + t282) * MDP(14) + (t167 * t219 + t168 * t220 + t178 * t211 + t179 * t309 - t182 * t215 + t183 * t216 - t189 * t199 - t200 * t252) * MDP(15) + (-t165 * t219 - t171 * t211 - t177 * t216 - t185 * t200 + t281) * MDP(16) + (-t165 * t220 - t171 * t309 + t177 * t215 + t185 * t199 - t282) * MDP(17) + (t165 * t185 - t167 * t252 + t168 * t189 + t171 * t177 + t178 * t183 + t179 * t182) * MDP(18) + (t175 * t288 + t251 * t193) * MDP(19) + ((-t191 * t240 + t193 * t242) * t216 + (t294 - t176 * t240 + (-t191 * t242 - t193 * t240) * qJD(5)) * t219) * MDP(20) + (t175 * t220 - t193 * t215 - t199 * t288 + t251 * t311) * MDP(21) + (-t219 * t195 - t176 * t220 + t191 * t215 + (t216 * t242 - t219 * t278) * t311) * MDP(22) + (-t215 * t311 - t188) * MDP(23) + (t259 * t215 + t159 * t220 + t163 * t191 + t181 * t176 + (-t157 * t220 - t161 * t311 + t295) * t240 + (t164 * t311 - t248) * t242 + ((-t172 * t242 - t180 * t240) * t311 - t156 * t220 + t169 * t288) * qJD(5)) * MDP(24) + (t156 * t215 + t163 * t193 + t181 * t175 + (-(qJD(5) * t180 + t161) * t311 + t295 - (qJD(5) * t166 + t157) * t220 + t169 * qJD(5) * t219) * t242 + (-(-qJD(5) * t172 + t164) * t311 - (-qJD(5) * t162 + t160) * t220 + t248) * t240) * MDP(25) + 0.2e1 * t273 * t307 + (-t215 * MDP(10) - t216 * MDP(11)) * qJD(3); t230 * MDP(14) + (-t304 - t305) * MDP(15) + (t199 + t201) * MDP(17) + (t301 + t297 - t183 * t211 + (-qJD(4) - t182) * t309) * MDP(18) + (t247 + t292) * MDP(24) + (t291 - t312) * MDP(25) - qJD(1) ^ 2 * t307 + ((-t211 - t271) * MDP(14) + (0.2e1 * MDP(13) - 0.2e1 * MDP(16)) * t309) * qJD(3); MDP(8) * t289 + (t304 - t305) * MDP(9) + (t230 + (t211 - t271) * qJD(3)) * MDP(10) + (-t223 * t309 + t250) * MDP(13) + (qJD(3) * t285 + t211 * t223 + t262) * MDP(14) + (pkin(3) * t199 - qJ(4) * t200 + (-t183 - t187) * t309 + (t182 - t308) * t211) * MDP(15) + (t184 * t211 - t250 + t293) * MDP(16) + (-t177 * t211 + t184 * t309 + (0.2e1 * qJD(4) - t285) * qJD(3) - t262) * MDP(17) + (-pkin(3) * t168 - qJ(4) * t167 - t177 * t184 - t182 * t187 - t183 * t308) * MDP(18) + (-t240 * t263 + t294) * MDP(19) + ((-t176 - t263) * t242 + (-t175 + t313) * t240) * MDP(20) + (t291 + t312) * MDP(21) + (t247 - t292) * MDP(22) + t311 * t211 * MDP(23) + (qJ(4) * t176 + t275 * t191 - t211 * t259 + t249 * t240 + t306 * t242) * MDP(24) + (qJ(4) * t175 - t156 * t211 + t275 * t193 - t306 * t240 + t249 * t242) * MDP(25); (-t199 + t201) * MDP(15) - MDP(16) * t289 + (-qJD(3) ^ 2 - t304) * MDP(17) + (qJD(3) * t183 + t168 + t293) * MDP(18) + (-qJD(3) * t191 - t195) * MDP(24) + (-qJD(3) * t193 + t290) * MDP(25) + (-MDP(24) * t265 - MDP(25) * t264) * t311; t193 * t191 * MDP(19) + (-t191 ^ 2 + t193 ^ 2) * MDP(20) + (t286 + t313) * MDP(21) + (t196 + t263) * MDP(22) - t199 * MDP(23) + (t156 * t311 - t157 * t240 - t169 * t193 + t159) * MDP(24) + (-t157 * t242 - t160 * t240 + t169 * t191 - t259 * t311) * MDP(25) + (-MDP(21) * t280 - t193 * MDP(22) - t156 * MDP(24) + t259 * MDP(25)) * qJD(5);];
tauc = t1;
