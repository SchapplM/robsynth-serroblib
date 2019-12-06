% Calculate vector of inverse dynamics joint torques for
% S5RPRPR3
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
%   pkin=[a2,a3,a4,a5,d1,d3,d5,theta2,theta4]';
% MDP [18x1]
%   Minimal dynamic parameter vector (fixed base model)
%   see S5RPRPR3_convert_par2_MPV_fixb.m
% 
% Output:
% tau [5x1]
%   joint torques of inverse dynamics (contains inertial, gravitational Coriolis and centrifugal forces)

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-05 17:52
% Revision: 77da58f92bca3eff71542919beafa37024070d86 (2019-12-05)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tau = S5RPRPR3_invdynJ_fixb_mdp_slag_vp(qJ, qJD, qJDD, g, pkin, MDP)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(5,1),zeros(3,1),zeros(9,1),zeros(18,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RPRPR3_invdynJ_fixb_mdp_slag_vp: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RPRPR3_invdynJ_fixb_mdp_slag_vp: qJD has to be [5x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [5 1]), ...
  'S5RPRPR3_invdynJ_fixb_mdp_slag_vp: qJDD has to be [5x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S5RPRPR3_invdynJ_fixb_mdp_slag_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S5RPRPR3_invdynJ_fixb_mdp_slag_vp: pkin has to be [9x1] (double)');
assert(isreal(MDP) && all(size(MDP) == [18 1]), ...
  'S5RPRPR3_invdynJ_fixb_mdp_slag_vp: MDP has to be [18x1] (double)'); 

%% Symbolic Calculation
% From invdyn_joint_fixb_mdp_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-05 17:51:39
% EndTime: 2019-12-05 17:51:43
% DurationCPUTime: 1.62s
% Computational Cost: add. (1208->215), mult. (2001->298), div. (0->0), fcn. (1186->14), ass. (0->113)
t231 = cos(pkin(8));
t215 = pkin(1) * t231 + pkin(2);
t229 = sin(pkin(8));
t304 = pkin(1) * t229;
t273 = qJD(3) * t304;
t309 = -qJD(1) * t273 + t215 * qJDD(1);
t225 = qJ(1) + pkin(8);
t219 = qJ(3) + t225;
t213 = sin(t219);
t214 = cos(t219);
t308 = g(2) * t214 + g(3) * t213;
t200 = t215 * qJD(1);
t307 = qJD(3) * t200 + qJDD(1) * t304;
t224 = qJD(1) + qJD(3);
t221 = qJDD(1) + qJDD(3);
t302 = qJ(4) * t221;
t233 = sin(qJ(3));
t236 = cos(qJ(3));
t305 = -t309 * t233 - t307 * t236;
t167 = qJD(4) * t224 + t302 - t305;
t228 = sin(pkin(9));
t230 = cos(pkin(9));
t164 = qJDD(2) * t228 + t167 * t230;
t163 = -t230 * qJDD(2) + t167 * t228;
t301 = t163 * t228;
t306 = t164 * t230 + t301;
t284 = t233 * t215 + t236 * t304;
t232 = sin(qJ(5));
t235 = cos(qJ(5));
t249 = MDP(17) * t232 + MDP(18) * t235;
t274 = qJD(1) * t304;
t183 = t200 * t236 - t233 * t274;
t258 = qJD(4) - t183;
t303 = pkin(3) * t221;
t184 = t200 * t233 + t236 * t274;
t176 = qJ(4) * t224 + t184;
t173 = -t230 * qJD(2) + t176 * t228;
t300 = t173 * t228;
t193 = -pkin(4) * t230 - pkin(7) * t228 - pkin(3);
t294 = t215 * t236;
t264 = -t233 * t304 + t294;
t177 = t193 - t264;
t299 = t177 * t235;
t298 = t184 * t224;
t248 = qJD(3) * t294 - t233 * t273;
t185 = qJD(4) + t248;
t297 = t185 * t224;
t296 = t185 * t232;
t186 = t284 * qJD(3);
t295 = t186 * t224;
t293 = t221 * t230;
t292 = t221 * t232;
t291 = t224 * t230;
t290 = t228 * t232;
t289 = t230 * t232;
t288 = t230 * t235;
t287 = t232 * t235;
t286 = t308 * t230;
t285 = g(2) * t213 - g(3) * t214;
t222 = t228 ^ 2;
t223 = t230 ^ 2;
t283 = t222 + t223;
t227 = t235 ^ 2;
t282 = t232 ^ 2 - t227;
t278 = qJD(5) * t232;
t277 = qJD(5) * t235;
t199 = -qJDD(5) + t293;
t276 = t199 * MDP(16);
t203 = -qJD(5) + t291;
t275 = -qJD(5) - t203;
t271 = t224 * t277;
t270 = -t213 * pkin(3) + t214 * qJ(4);
t269 = t283 * t221;
t261 = -t307 * t233 + t309 * t236;
t250 = qJDD(4) - t261;
t162 = t193 * t221 + t250;
t268 = t235 * t162 - t232 * t164;
t266 = t199 - t293;
t265 = t199 + t293;
t263 = t258 * t235;
t234 = sin(qJ(1));
t237 = cos(qJ(1));
t259 = g(2) * t237 + g(3) * t234;
t257 = -t214 * pkin(3) - t213 * qJ(4);
t256 = t298 + t303;
t255 = t232 * t162 + t235 * t164;
t170 = t193 * t224 + t258;
t174 = qJD(2) * t228 + t176 * t230;
t254 = t170 * t235 - t174 * t232;
t253 = -t170 * t232 - t174 * t235;
t252 = t174 * t230 + t300;
t188 = -pkin(3) - t264;
t251 = -t188 * t221 - t295;
t247 = t285 + t306;
t246 = -t261 - t308;
t245 = -qJ(4) * t289 + t193 * t235;
t179 = t213 * t289 + t214 * t235;
t181 = -t213 * t235 + t214 * t289;
t243 = -g(2) * t181 - g(3) * t179 + (t254 * qJD(5) + t255) * t230 + t235 * t301;
t180 = t213 * t288 - t214 * t232;
t182 = -t213 * t232 - t214 * t288;
t242 = -g(2) * t182 + g(3) * t180 + t163 * t290 + t277 * t300;
t241 = qJ(4) * t277 + t258 * t232;
t190 = t228 * t278 * t291;
t240 = t230 * t276 + (t190 + (t203 * t278 - t265 * t235) * t228) * MDP(14) + (t265 * t232 + (t203 + t291) * t277) * t228 * MDP(15) + t221 * MDP(5) + (0.2e1 * (t282 * t224 * qJD(5) - t221 * t287) * MDP(13) + (t221 * t227 - 0.2e1 * t232 * t271) * MDP(12)) * t222;
t239 = -t285 + t305;
t220 = t224 ^ 2;
t187 = qJ(4) + t284;
t175 = -pkin(3) * t224 + t258;
t169 = t250 - t303;
t168 = t169 * t228;
t157 = t253 * qJD(5) + t268;
t1 = [((-t169 + t251) * t230 + t286) * MDP(8) + (t259 + (t229 ^ 2 + t231 ^ 2) * qJDD(1) * pkin(1)) * pkin(1) * MDP(4) + (-(-t177 * t278 + t186 * t235) * t203 - t199 * t299 + (-(-t187 * t277 - t296) * t203 + t187 * t232 * t199 - t157) * t230 + (t224 * t296 + (t271 + t292) * t187) * t222 + t242) * MDP(17) + (t187 * t269 + t283 * t297 + t247) * MDP(10) + qJDD(1) * MDP(1) + (t264 * t221 - t246 - t295) * MDP(6) + t259 * MDP(2) + (-g(2) * t234 + g(3) * t237) * MDP(3) + (-t284 * t221 - t248 * t224 + t239) * MDP(7) + ((t185 * t288 + t186 * t232) * t203 + (t177 * t232 + t187 * t288) * t199 + (t187 * t221 + t297) * t235 * t222 + (t203 * t299 + (-t300 + (-t203 * t230 - t222 * t224) * t187) * t232) * qJD(5) + t243) * MDP(18) + (t168 + (-t251 - t308) * t228) * MDP(9) + t240 + (t169 * t188 + t175 * t186 - g(2) * (-pkin(2) * cos(t225) - t237 * pkin(1) + t257) - g(3) * (-pkin(2) * sin(t225) - t234 * pkin(1) + t270) + t306 * t187 + t252 * t185) * MDP(11); (qJDD(2) - g(1)) * MDP(4) + (-t163 * t230 - g(1)) * MDP(11) + t190 * MDP(18) + (t164 * MDP(11) + (-qJD(5) * t203 * MDP(18) + t266 * MDP(17)) * t232 + (t266 * MDP(18) + (t203 - t291) * MDP(17) * qJD(5)) * t235) * t228; (-t246 + t298) * MDP(6) + (t183 * t224 + t239) * MDP(7) + ((-t169 + t256) * t230 + t286) * MDP(8) + (t168 + (-t256 - t308) * t228) * MDP(9) + (t258 * t224 * t283 + qJ(4) * t269 + t247) * MDP(10) + (-t169 * pkin(3) - t175 * t184 - g(2) * t257 - g(3) * t270 + (t164 * qJ(4) + t258 * t174) * t230 + (t163 * qJ(4) + t258 * t173) * t228) * MDP(11) + (-t245 * t199 - t157 * t230 + (t235 * t184 + t193 * t278 + t241 * t230) * t203 + (qJ(4) * t292 + t241 * t224) * t222 + t242) * MDP(17) + ((qJ(4) * t288 + t193 * t232) * t199 + (-t232 * t184 + t230 * t263) * t203 + (-t173 * t290 + t245 * t203) * qJD(5) + (t235 * t302 + (-qJ(4) * t278 + t263) * t224) * t222 + t243) * MDP(18) + t240; (-t252 * t224 + qJDD(4) + t246) * MDP(11) + (-pkin(3) * MDP(11) - MDP(8) * t230 + MDP(9) * t228) * t221 + (-MDP(17) * t235 + MDP(18) * t232) * t199 + (-t223 * MDP(10) + (-MDP(10) - t249) * t222) * t220 - t249 * t203 ^ 2; -t276 + (-g(2) * t179 + g(3) * t181 + t253 * t203 + t268) * MDP(17) + (-g(2) * t180 - g(3) * t182 - t254 * t203 - t255) * MDP(18) + (t253 * MDP(17) - t254 * MDP(18)) * qJD(5) + (MDP(12) * t287 - t282 * MDP(13)) * t222 * t220 + ((MDP(14) * t235 - MDP(15) * t232) * t221 + t249 * g(1) + ((t275 * MDP(15) - t173 * MDP(17)) * t235 + (t275 * MDP(14) + t173 * MDP(18)) * t232) * t224) * t228;];
tau = t1;
