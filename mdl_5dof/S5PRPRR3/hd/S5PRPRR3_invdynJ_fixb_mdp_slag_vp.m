% Calculate vector of inverse dynamics joint torques for
% S5PRPRR3
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
% MDP [19x1]
%   Minimal dynamic parameter vector (fixed base model)
%   see S5PRPRR3_convert_par2_MPV_fixb.m
% 
% Output:
% tau [5x1]
%   joint torques of inverse dynamics (contains inertial, gravitational Coriolis and centrifugal forces)

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-05 15:48
% Revision: 77da58f92bca3eff71542919beafa37024070d86 (2019-12-05)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tau = S5PRPRR3_invdynJ_fixb_mdp_slag_vp(qJ, qJD, qJDD, g, pkin, MDP)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(5,1),zeros(3,1),zeros(9,1),zeros(19,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5PRPRR3_invdynJ_fixb_mdp_slag_vp: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5PRPRR3_invdynJ_fixb_mdp_slag_vp: qJD has to be [5x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [5 1]), ...
  'S5PRPRR3_invdynJ_fixb_mdp_slag_vp: qJDD has to be [5x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S5PRPRR3_invdynJ_fixb_mdp_slag_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S5PRPRR3_invdynJ_fixb_mdp_slag_vp: pkin has to be [9x1] (double)');
assert(isreal(MDP) && all(size(MDP) == [19 1]), ...
  'S5PRPRR3_invdynJ_fixb_mdp_slag_vp: MDP has to be [19x1] (double)'); 

%% Symbolic Calculation
% From invdyn_joint_fixb_mdp_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-05 15:47:43
% EndTime: 2019-12-05 15:47:48
% DurationCPUTime: 1.75s
% Computational Cost: add. (830->205), mult. (1684->287), div. (0->0), fcn. (1272->14), ass. (0->117)
t256 = cos(qJ(2));
t304 = qJD(1) * t256;
t219 = qJD(2) * pkin(2) + t304;
t247 = sin(pkin(9));
t253 = sin(qJ(2));
t305 = qJD(1) * t253;
t222 = t247 * t305;
t249 = cos(pkin(9));
t190 = t219 * t249 - t222;
t255 = cos(qJ(4));
t289 = -pkin(4) * t255 - pkin(3);
t183 = t289 * qJD(2) - t190;
t198 = t249 * t304 - t222;
t325 = t183 + t198;
t243 = qJ(2) + pkin(9);
t234 = sin(t243);
t235 = cos(t243);
t248 = sin(pkin(8));
t250 = cos(pkin(8));
t280 = g(1) * t250 + g(2) * t248;
t322 = -g(3) * t235 + t280 * t234;
t242 = qJD(4) + qJD(5);
t294 = qJD(1) * qJD(2);
t324 = qJDD(1) * t253 + t256 * t294;
t237 = t256 * qJDD(1);
t202 = qJDD(2) * pkin(2) - t253 * t294 + t237;
t179 = t202 * t249 - t324 * t247;
t175 = -qJDD(2) * pkin(3) - t179;
t223 = t249 * t305;
t196 = t247 * t304 + t223;
t230 = pkin(2) * t247 + pkin(6);
t318 = pkin(2) * t249;
t231 = -pkin(3) - t318;
t257 = qJD(4) ^ 2;
t323 = qJD(2) * t196 - qJDD(2) * t231 - t230 * t257 - t175 + t322;
t251 = sin(qJ(5));
t254 = cos(qJ(5));
t290 = qJDD(2) * t255;
t252 = sin(qJ(4));
t291 = qJDD(2) * t252;
t208 = t251 * t255 + t252 * t254;
t321 = t242 * t208;
t171 = qJD(2) * t321 + t251 * t291 - t254 * t290;
t207 = t251 * t252 - t254 * t255;
t268 = t242 * t207;
t191 = t247 * t219 + t223;
t283 = t191 + (pkin(6) + pkin(7)) * qJD(2);
t177 = t255 * qJD(3) - t283 * t252;
t319 = t242 * t255;
t178 = qJD(3) * t252 + t283 * t255;
t317 = g(3) * t234;
t315 = pkin(7) + t230;
t314 = t178 * t254;
t241 = qJDD(4) + qJDD(5);
t313 = t207 * t241;
t312 = t208 * t241;
t246 = qJ(4) + qJ(5);
t239 = sin(t246);
t311 = t239 * t248;
t310 = t239 * t250;
t240 = cos(t246);
t309 = t240 * t248;
t308 = t240 * t250;
t307 = qJDD(1) - g(3);
t244 = t252 ^ 2;
t306 = -t255 ^ 2 + t244;
t205 = t247 * t253 - t249 * t256;
t303 = qJD(2) * t205;
t302 = qJD(2) * t252;
t301 = qJD(2) * t255;
t299 = qJD(4) * t252;
t297 = qJD(5) * t251;
t296 = qJD(5) * t254;
t293 = qJD(2) * qJD(4);
t180 = t247 * t202 + t324 * t249;
t288 = t251 * t302;
t287 = t254 * t301;
t285 = t255 * t293;
t284 = qJD(4) * t315;
t176 = qJDD(2) * pkin(6) + t180;
t282 = pkin(7) * qJDD(2) + t176;
t281 = pkin(4) * t299 - t196;
t279 = -g(1) * t248 + g(2) * t250;
t174 = qJD(4) * pkin(4) + t177;
t276 = -t174 * t251 - t314;
t275 = -t242 * t268 + t312;
t203 = t315 * t252;
t204 = t315 * t255;
t274 = -t203 * t254 - t204 * t251;
t273 = -t203 * t251 + t204 * t254;
t206 = t247 * t256 + t249 * t253;
t272 = -qJDD(3) - t279;
t170 = qJD(5) * t287 - t242 * t288 + t251 * t290 + (t285 + t291) * t254;
t199 = -t287 + t288;
t201 = -t251 * t301 - t254 * t302;
t267 = -t201 * t199 * MDP(13) + (t199 * t242 + t170) * MDP(15) + (-t201 * t242 - t171) * MDP(16) + (-t199 ^ 2 + t201 ^ 2) * MDP(14) + t241 * MDP(17);
t195 = t206 * qJD(2);
t266 = qJD(2) * t195 + qJDD(2) * t205 + t206 * t257;
t265 = 0.2e1 * t303 * qJD(4) - qJDD(4) * t206;
t264 = -g(3) * t256 + t280 * t253;
t184 = -qJD(2) * pkin(3) - t190;
t263 = -qJDD(4) * t230 + (qJD(2) * t231 + t184 + t198) * qJD(4);
t261 = -t184 * qJD(2) + t280 * t235 - t176 + t317;
t236 = t255 * qJDD(3);
t164 = qJDD(4) * pkin(4) - t178 * qJD(4) - t282 * t252 + t236;
t260 = -g(1) * (-t235 * t308 - t311) - g(2) * (-t235 * t309 + t310) + t183 * t199 + t178 * t297 + t240 * t317 + (-t178 * t242 - t164) * t251;
t165 = t177 * qJD(4) + t252 * qJDD(3) + t282 * t255;
t259 = -g(1) * (-t235 * t310 + t309) - g(2) * (-t235 * t311 - t308) + t276 * qJD(5) + t254 * t164 - t251 * t165 + t183 * t201 + t239 * t317;
t258 = qJD(2) ^ 2;
t214 = qJDD(4) * t255 - t252 * t257;
t213 = qJDD(4) * t252 + t255 * t257;
t212 = t289 - t318;
t194 = t255 * t284;
t193 = t252 * t284;
t169 = -t242 * t321 - t313;
t168 = (t252 * t293 - t290) * pkin(4) + t175;
t1 = [t307 * MDP(1) + (qJDD(2) * t256 - t253 * t258) * MDP(3) + (-qJDD(2) * t253 - t256 * t258) * MDP(4) + (-t179 * t205 + t180 * t206 - t190 * t195 - t191 * t303 - g(3)) * MDP(5) + (t265 * t252 - t266 * t255) * MDP(11) + (t266 * t252 + t265 * t255) * MDP(12) + (t205 * t171 + t195 * t199 + t303 * t321 + ((t251 * t299 + t252 * t297 - t319 * t254) * t242 - t312) * t206) * MDP(18) + (t205 * t170 - t195 * t201 - t303 * t268 + (-(-t319 * t251 - t252 * t296 - t254 * t299) * t242 + t313) * t206) * MDP(19); qJDD(2) * MDP(2) + (t237 + t264) * MDP(3) + (-t307 * t253 + t280 * t256) * MDP(4) + (t190 * t196 - t191 * t198 + (t179 * t249 + t180 * t247 + t264) * pkin(2)) * MDP(5) + (qJDD(2) * t244 + 0.2e1 * t252 * t285) * MDP(6) + 0.2e1 * (t252 * t290 - t306 * t293) * MDP(7) + t213 * MDP(8) + t214 * MDP(9) + (t263 * t252 + t323 * t255) * MDP(11) + (-t323 * t252 + t263 * t255) * MDP(12) + (t170 * t208 + t201 * t268) * MDP(13) + (-t170 * t207 - t171 * t208 + t199 * t268 + t201 * t321) * MDP(14) + t275 * MDP(15) + t169 * MDP(16) + ((-t273 * qJD(5) + t193 * t251 - t194 * t254) * t242 + t274 * t241 + t212 * t171 + t168 * t207 + t281 * t199 + t322 * t240 + t325 * t321) * MDP(18) + (-(t274 * qJD(5) - t193 * t254 - t194 * t251) * t242 - t273 * t241 + t212 * t170 + t168 * t208 - t281 * t201 - t322 * t239 - t325 * t268) * MDP(19); t214 * MDP(11) - t213 * MDP(12) + t169 * MDP(18) - t275 * MDP(19) - t272 * MDP(5); MDP(8) * t291 + MDP(9) * t290 + qJDD(4) * MDP(10) + (t261 * t252 + t279 * t255 + t236) * MDP(11) + (t272 * t252 + t261 * t255) * MDP(12) + (-(-t177 * t251 - t314) * t242 + (-t199 * t302 + t254 * t241 - t242 * t297) * pkin(4) + t259) * MDP(18) + ((-qJD(5) * t174 + t177 * t242 - t165) * t254 + (t201 * t302 - t251 * t241 - t242 * t296) * pkin(4) + t260) * MDP(19) + t267 + (-t252 * t255 * MDP(6) + t306 * MDP(7)) * t258; (-t276 * t242 + t259) * MDP(18) + ((-t165 + (-qJD(5) + t242) * t174) * t254 + t260) * MDP(19) + t267;];
tau = t1;
