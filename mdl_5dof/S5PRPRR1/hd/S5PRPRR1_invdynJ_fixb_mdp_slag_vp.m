% Calculate vector of inverse dynamics joint torques for
% S5PRPRR1
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% qJD [5x1]
%   Generalized joint velocities
% qJDD [5x1]
%   Generalized joint accelerations
% g [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [9x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d2,d4,d5,theta1,theta3]';
% MDP [22x1]
%   Minimal dynamic parameter vector (fixed base model)
%   see S5PRPRR1_convert_par2_MPV_fixb.m
% 
% Output:
% tau [5x1]
%   joint torques of inverse dynamics (contains inertial, gravitational Coriolis and centrifugal forces)

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-05 15:43
% Revision: 77da58f92bca3eff71542919beafa37024070d86 (2019-12-05)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tau = S5PRPRR1_invdynJ_fixb_mdp_slag_vp(qJ, qJD, qJDD, g, pkin, MDP)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(5,1),zeros(3,1),zeros(9,1),zeros(22,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5PRPRR1_invdynJ_fixb_mdp_slag_vp: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5PRPRR1_invdynJ_fixb_mdp_slag_vp: qJD has to be [5x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [5 1]), ...
  'S5PRPRR1_invdynJ_fixb_mdp_slag_vp: qJDD has to be [5x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S5PRPRR1_invdynJ_fixb_mdp_slag_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S5PRPRR1_invdynJ_fixb_mdp_slag_vp: pkin has to be [9x1] (double)');
assert(isreal(MDP) && all(size(MDP) == [22 1]), ...
  'S5PRPRR1_invdynJ_fixb_mdp_slag_vp: MDP has to be [22x1] (double)'); 

%% Symbolic Calculation
% From invdyn_joint_fixb_mdp_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-05 15:43:11
% EndTime: 2019-12-05 15:43:15
% DurationCPUTime: 1.50s
% Computational Cost: add. (1113->223), mult. (2416->291), div. (0->0), fcn. (1914->12), ass. (0->112)
t256 = sin(pkin(9));
t259 = sin(qJ(4));
t257 = cos(pkin(9));
t261 = cos(qJ(4));
t300 = t261 * t257;
t214 = t256 * t259 - t300;
t206 = t214 * qJD(2);
t260 = cos(qJ(5));
t215 = t256 * t261 + t257 * t259;
t207 = t215 * qJD(2);
t258 = sin(qJ(5));
t302 = t207 * t258;
t178 = t260 * t206 + t302;
t255 = qJD(4) + qJD(5);
t304 = t178 * t255;
t254 = pkin(8) + qJ(2);
t246 = sin(t254);
t248 = cos(t254);
t278 = g(1) * t248 + g(2) * t246;
t309 = pkin(6) + qJ(3);
t224 = t309 * t256;
t244 = t257 * qJD(1);
t203 = -qJD(2) * t224 + t244;
t295 = qJD(2) * t257;
t221 = qJ(3) * t295 + t256 * qJD(1);
t204 = pkin(6) * t295 + t221;
t270 = -t203 * t259 - t204 * t261;
t173 = -pkin(7) * t206 - t270;
t238 = -pkin(3) * t257 - pkin(2);
t223 = t238 * qJD(2) + qJD(3);
t191 = pkin(4) * t206 + t223;
t253 = pkin(9) + qJ(4);
t249 = qJ(5) + t253;
t236 = sin(t249);
t237 = cos(t249);
t292 = qJD(5) * t258;
t319 = g(3) * t236 + t173 * t292 + t191 * t178 + t278 * t237;
t250 = qJDD(4) + qJDD(5);
t269 = -t206 * t258 + t260 * t207;
t318 = t250 * MDP(20) + t178 * MDP(16) * t269 + (-t178 ^ 2 + t269 ^ 2) * MDP(17);
t305 = t269 * t255;
t316 = t261 * t203 - t204 * t259;
t225 = t309 * t257;
t299 = -t259 * t224 + t261 * t225;
t289 = qJD(2) * qJD(3);
t315 = qJ(3) * qJDD(2) + t289;
t277 = g(1) * t246 - g(2) * t248;
t314 = t277 - qJDD(3);
t294 = qJD(2) * t259;
t283 = t256 * t294;
t286 = qJDD(2) * t261;
t287 = qJDD(2) * t259;
t293 = qJD(4) * t261;
t284 = t256 * t286 + t257 * t287 + t293 * t295;
t185 = -qJD(4) * t283 + t284;
t242 = t257 * qJDD(1);
t193 = t242 + (-t309 * qJDD(2) - t289) * t256;
t288 = qJDD(2) * t257;
t202 = qJ(3) * t288 + t256 * qJDD(1) + t257 * t289;
t194 = pkin(6) * t288 + t202;
t281 = t261 * t193 - t259 * t194;
t161 = qJDD(4) * pkin(4) - pkin(7) * t185 + t270 * qJD(4) + t281;
t209 = t215 * qJD(4);
t229 = t257 * t286;
t276 = -t256 * t287 + t229;
t186 = qJD(2) * t209 - t276;
t271 = t259 * t193 + t261 * t194;
t162 = -pkin(7) * t186 + t316 * qJD(4) + t271;
t313 = -g(3) * t237 + t260 * t161 - t258 * t162 - t191 * t269 + t278 * t236;
t282 = t185 * t258 + t260 * t186;
t165 = t269 * qJD(5) + t282;
t312 = pkin(4) * t209;
t308 = qJDD(2) * pkin(2);
t172 = -pkin(7) * t207 + t316;
t169 = qJD(4) * pkin(4) + t172;
t307 = t169 * t260;
t306 = t173 * t260;
t301 = t257 * MDP(5);
t298 = t256 ^ 2 + t257 ^ 2;
t296 = qJD(2) * t256;
t291 = qJD(5) * t260;
t285 = t260 * t185 - t258 * t186 - t206 * t291;
t280 = -t261 * t224 - t225 * t259;
t189 = t260 * t214 + t215 * t258;
t208 = t214 * qJD(4);
t166 = -t189 * qJD(5) - t208 * t260 - t209 * t258;
t190 = -t214 * t258 + t215 * t260;
t275 = t166 * t255 + t190 * t250;
t274 = -t169 * t258 - t306;
t175 = -pkin(7) * t215 + t280;
t176 = -pkin(7) * t214 + t299;
t273 = t175 * t260 - t176 * t258;
t272 = t175 * t258 + t176 * t260;
t268 = t308 + t314;
t222 = t238 * qJDD(2) + qJDD(3);
t164 = -t207 * t292 + t285;
t265 = -t224 * t293 + qJD(3) * t300 + (-qJD(3) * t256 - qJD(4) * t225) * t259;
t201 = -t315 * t256 + t242;
t264 = -t201 * t256 + t202 * t257 - t278;
t263 = -t215 * qJD(3) - t299 * qJD(4);
t247 = cos(t253);
t245 = sin(t253);
t220 = -qJ(3) * t296 + t244;
t196 = pkin(4) * t214 + t238;
t188 = -qJD(4) * t209 - qJDD(4) * t214;
t187 = -qJD(4) * t208 + qJDD(4) * t215;
t174 = pkin(4) * t186 + t222;
t171 = pkin(7) * t208 + t263;
t170 = -pkin(7) * t209 + t265;
t167 = t190 * qJD(5) - t208 * t258 + t260 * t209;
t163 = -t167 * t255 - t189 * t250;
t1 = [(qJDD(1) - g(3)) * MDP(1) + (t201 * t257 + t202 * t256 - g(3)) * MDP(8) + t188 * MDP(14) - t187 * MDP(15) + t163 * MDP(21) - t275 * MDP(22); qJDD(2) * MDP(2) + t277 * MDP(3) + t278 * MDP(4) + (t315 * t298 + t264) * MDP(7) + ((-t220 * t256 + t221 * t257) * qJD(3) + t268 * pkin(2) + t264 * qJ(3)) * MDP(8) + (t185 * t215 - t207 * t208) * MDP(9) + (-t185 * t214 - t186 * t215 + t206 * t208 - t207 * t209) * MDP(10) + t187 * MDP(11) + t188 * MDP(12) + (t263 * qJD(4) + t280 * qJDD(4) + t238 * t186 + t223 * t209 + t222 * t214 + t277 * t247) * MDP(14) + (-t265 * qJD(4) - t299 * qJDD(4) + t238 * t185 - t223 * t208 + t222 * t215 - t277 * t245) * MDP(15) + (t164 * t190 + t166 * t269) * MDP(16) + (-t164 * t189 - t165 * t190 - t166 * t178 - t167 * t269) * MDP(17) + t275 * MDP(18) + t163 * MDP(19) + (t178 * t312 + t196 * t165 + t174 * t189 + t191 * t167 + (-t272 * qJD(5) - t170 * t258 + t171 * t260) * t255 + t273 * t250 + t277 * t237) * MDP(21) + (t269 * t312 + t196 * t164 + t174 * t190 + t191 * t166 - (t273 * qJD(5) + t170 * t260 + t171 * t258) * t255 - t272 * t250 - t277 * t236) * MDP(22) + (-MDP(6) * t256 + t301) * (t268 + t308); (t220 * t296 - t221 * t295 - t314) * MDP(8) - t229 * MDP(14) + t284 * MDP(15) + (t165 + t305) * MDP(21) + (t164 - t304) * MDP(22) - t298 * MDP(7) * qJD(2) ^ 2 + (-t301 - pkin(2) * MDP(8) + (MDP(14) * t259 + MDP(6)) * t256) * qJDD(2) + ((t257 * t294 + t261 * t296 + t207) * MDP(14) + (-t206 - t283) * MDP(15)) * qJD(4); t207 * t206 * MDP(9) + (-t206 ^ 2 + t207 ^ 2) * MDP(10) + (t284 + (t206 - t283) * qJD(4)) * MDP(11) + t276 * MDP(12) + qJDD(4) * MDP(13) + (-g(3) * t247 - t223 * t207 + t278 * t245 + t281) * MDP(14) + (g(3) * t245 + t223 * t206 + t278 * t247 - t271) * MDP(15) + (t164 + t304) * MDP(18) + (-t165 + t305) * MDP(19) + (-(-t172 * t258 - t306) * t255 + t274 * qJD(5) + (-t178 * t207 + t250 * t260 - t255 * t292) * pkin(4) + t313) * MDP(21) + ((-t173 * t255 - t161) * t258 + (-qJD(5) * t169 + t172 * t255 - t162) * t260 + (-t207 * t269 - t250 * t258 - t255 * t291) * pkin(4) + t319) * MDP(22) + t318; (t285 + t304) * MDP(18) + (-t282 + t305) * MDP(19) + (-t274 * t255 + t313) * MDP(21) + (-t260 * t162 - t258 * t161 + (-t173 * t258 + t307) * t255 + t319) * MDP(22) + (-MDP(18) * t302 - t269 * MDP(19) + t274 * MDP(21) - MDP(22) * t307) * qJD(5) + t318;];
tau = t1;
