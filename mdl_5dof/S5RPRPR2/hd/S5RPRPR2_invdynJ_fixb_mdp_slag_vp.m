% Calculate vector of inverse dynamics joint torques for
% S5RPRPR2
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
%   see S5RPRPR2_convert_par2_MPV_fixb.m
% 
% Output:
% tau [5x1]
%   joint torques of inverse dynamics (contains inertial, gravitational Coriolis and centrifugal forces)

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-05 17:50
% Revision: 77da58f92bca3eff71542919beafa37024070d86 (2019-12-05)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tau = S5RPRPR2_invdynJ_fixb_mdp_slag_vp(qJ, qJD, qJDD, g, pkin, MDP)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(5,1),zeros(3,1),zeros(9,1),zeros(18,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RPRPR2_invdynJ_fixb_mdp_slag_vp: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RPRPR2_invdynJ_fixb_mdp_slag_vp: qJD has to be [5x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [5 1]), ...
  'S5RPRPR2_invdynJ_fixb_mdp_slag_vp: qJDD has to be [5x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S5RPRPR2_invdynJ_fixb_mdp_slag_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S5RPRPR2_invdynJ_fixb_mdp_slag_vp: pkin has to be [9x1] (double)');
assert(isreal(MDP) && all(size(MDP) == [18 1]), ...
  'S5RPRPR2_invdynJ_fixb_mdp_slag_vp: MDP has to be [18x1] (double)'); 

%% Symbolic Calculation
% From invdyn_joint_fixb_mdp_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-05 17:49:44
% EndTime: 2019-12-05 17:49:48
% DurationCPUTime: 1.41s
% Computational Cost: add. (1090->192), mult. (1784->244), div. (0->0), fcn. (1138->16), ass. (0->108)
t259 = cos(pkin(8));
t238 = pkin(1) * t259 + pkin(2);
t257 = sin(pkin(8));
t319 = pkin(1) * t257;
t298 = qJD(3) * t319;
t326 = -qJD(1) * t298 + t238 * qJDD(1);
t255 = qJ(1) + pkin(8);
t247 = qJ(3) + t255;
t236 = sin(t247);
t237 = cos(t247);
t325 = g(2) * t237 + g(3) * t236;
t224 = t238 * qJD(1);
t324 = qJD(3) * t224 + qJDD(1) * t319;
t256 = sin(pkin(9));
t258 = cos(pkin(9));
t260 = sin(qJ(5));
t263 = cos(qJ(5));
t207 = t256 * t263 + t258 * t260;
t254 = qJD(1) + qJD(3);
t199 = t207 * t254;
t301 = t256 ^ 2 + t258 ^ 2;
t323 = t301 * t254;
t250 = qJDD(1) + qJDD(3);
t261 = sin(qJ(3));
t264 = cos(qJ(3));
t320 = -t326 * t261 - t324 * t264;
t172 = qJ(4) * t250 + qJD(4) * t254 - t320;
t242 = t258 * qJDD(2);
t168 = -t172 * t256 + t242;
t169 = t256 * qJDD(2) + t258 * t172;
t322 = -t168 * t256 + t169 * t258;
t303 = g(2) * t236 - g(3) * t237;
t302 = t261 * t238 + t264 * t319;
t299 = qJD(1) * t319;
t195 = t224 * t261 + t264 * t299;
t190 = qJ(4) * t254 + t195;
t181 = t258 * qJD(2) - t190 * t256;
t182 = t256 * qJD(2) + t258 * t190;
t321 = -t181 * t256 + t182 * t258;
t194 = t264 * t224 - t261 * t299;
t283 = qJD(4) - t194;
t318 = pkin(3) * t250;
t317 = pkin(4) * t258;
t248 = t258 * pkin(7);
t313 = t195 * t254;
t200 = t302 * qJD(3);
t312 = t200 * t254;
t253 = pkin(9) + qJ(5);
t246 = cos(t253);
t311 = t236 * t246;
t310 = t237 * t246;
t309 = t238 * t264;
t308 = t256 * t260;
t305 = t263 * t258;
t304 = t325 * t258;
t296 = t254 * t308;
t295 = t254 * t305;
t294 = qJD(5) * t295 + t207 * t250;
t239 = -pkin(3) - t317;
t293 = -t236 * pkin(3) + t237 * qJ(4);
t292 = t301 * t250;
t290 = -t261 * t319 + t309;
t287 = -t324 * t261 + t326 * t264;
t274 = qJDD(4) - t287;
t170 = t239 * t250 + t274;
t183 = t239 * t254 + t283;
t204 = t207 * qJD(5);
t206 = -t305 + t308;
t288 = g(2) * t310 + g(3) * t311 + t170 * t206 + t183 * t204;
t202 = -pkin(3) - t290;
t262 = sin(qJ(1));
t265 = cos(qJ(1));
t285 = g(2) * t265 + g(3) * t262;
t219 = t250 * t305;
t284 = -t250 * t308 + t219;
t282 = -t237 * pkin(3) - t236 * qJ(4);
t281 = t313 + t318;
t179 = -qJD(5) * t296 + t294;
t180 = t254 * t204 - t284;
t203 = t206 * qJD(5);
t186 = -qJD(5) * t203 + qJDD(5) * t207;
t187 = -qJD(5) * t204 - qJDD(5) * t206;
t197 = -t295 + t296;
t280 = (-t179 * t206 - t180 * t207 + t197 * t203 - t199 * t204) * MDP(13) + (t179 * t207 - t199 * t203) * MDP(12) + t186 * MDP(14) + t187 * MDP(15) + t250 * MDP(5);
t201 = qJ(4) + t302;
t191 = (-pkin(7) - t201) * t256;
t192 = t201 * t258 + t248;
t279 = t191 * t263 - t192 * t260;
t278 = t191 * t260 + t192 * t263;
t277 = -t202 * t250 - t312;
t215 = (-pkin(7) - qJ(4)) * t256;
t216 = qJ(4) * t258 + t248;
t276 = t215 * t263 - t216 * t260;
t275 = t215 * t260 + t216 * t263;
t273 = qJD(3) * t309 - t261 * t298;
t272 = t303 + t322;
t271 = t287 + t325;
t245 = sin(t253);
t268 = t170 * t207 - t183 * t203 - t245 * t325;
t267 = -t303 + t320;
t196 = qJD(4) + t273;
t193 = t202 - t317;
t189 = -pkin(3) * t254 + t283;
t174 = t274 - t318;
t173 = t174 * t256;
t164 = t250 * t248 + t169;
t163 = t242 + (-pkin(7) * t250 - t172) * t256;
t1 = [qJDD(1) * MDP(1) + t285 * MDP(2) + (-g(2) * t262 + g(3) * t265) * MDP(3) + (t285 + (t257 ^ 2 + t259 ^ 2) * qJDD(1) * pkin(1)) * pkin(1) * MDP(4) + (t290 * t250 + t271 - t312) * MDP(6) + (-t302 * t250 - t273 * t254 + t267) * MDP(7) + ((-t174 + t277) * t258 + t304) * MDP(8) + (t173 + (-t277 - t325) * t256) * MDP(9) + (t196 * t323 + t201 * t292 + t272) * MDP(10) + (t174 * t202 + t189 * t200 - g(2) * (-pkin(2) * cos(t255) - t265 * pkin(1) + t282) - g(3) * (-pkin(2) * sin(t255) - t262 * pkin(1) + t293) + t322 * t201 + t321 * t196) * MDP(11) + (t200 * t197 + t193 * t180 + t279 * qJDD(5) + (-t278 * qJD(5) - t207 * t196) * qJD(5) + t288) * MDP(17) + (t200 * t199 + t193 * t179 - t278 * qJDD(5) + (-t279 * qJD(5) + t206 * t196) * qJD(5) + t268) * MDP(18) + t280; (qJDD(2) - g(1)) * MDP(4) + (t168 * t258 + t169 * t256 - g(1)) * MDP(11) + t187 * MDP(17) - t186 * MDP(18); (t271 + t313) * MDP(6) + (t194 * t254 + t267) * MDP(7) + ((-t174 + t281) * t258 + t304) * MDP(8) + (t173 + (-t281 - t325) * t256) * MDP(9) + (qJ(4) * t292 + t283 * t323 + t272) * MDP(10) + (-t174 * pkin(3) - t189 * t195 - g(2) * t282 - g(3) * t293 + (t169 * qJ(4) + t283 * t182) * t258 + (-t168 * qJ(4) - t283 * t181) * t256) * MDP(11) + (t239 * t180 + t276 * qJDD(5) - t195 * t197 + (-t275 * qJD(5) - t283 * t207) * qJD(5) + t288) * MDP(17) + (t239 * t179 - t275 * qJDD(5) - t195 * t199 + (-t276 * qJD(5) + t283 * t206) * qJD(5) + t268) * MDP(18) + t280; (qJDD(4) - t271) * MDP(11) - t219 * MDP(17) + t294 * MDP(18) + (-pkin(3) * MDP(11) - t258 * MDP(8) + (t260 * MDP(17) + MDP(9)) * t256) * t250 + (0.2e1 * t199 * MDP(17) + (-t197 - t296) * MDP(18)) * qJD(5) + (-MDP(10) * t323 - t321 * MDP(11)) * t254; t199 * t197 * MDP(12) + (-t197 ^ 2 + t199 ^ 2) * MDP(13) + t284 * MDP(15) + qJDD(5) * MDP(16) + (-g(1) * t246 + t263 * t163 - t260 * t164 - t183 * t199 - t303 * t245) * MDP(17) + (g(1) * t245 - g(2) * t311 + g(3) * t310 - t260 * t163 - t263 * t164 + t183 * t197) * MDP(18) + (t294 + (t197 - t296) * qJD(5)) * MDP(14);];
tau = t1;
