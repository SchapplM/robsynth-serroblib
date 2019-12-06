% Calculate vector of inverse dynamics joint torques for
% S5RPPRR3
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
%   pkin=[a2,a3,a4,a5,d1,d4,d5,theta2,theta3]';
% MDP [22x1]
%   Minimal dynamic parameter vector (fixed base model)
%   see S5RPPRR3_convert_par2_MPV_fixb.m
% 
% Output:
% tau [5x1]
%   joint torques of inverse dynamics (contains inertial, gravitational Coriolis and centrifugal forces)

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-05 17:42
% Revision: 77da58f92bca3eff71542919beafa37024070d86 (2019-12-05)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tau = S5RPPRR3_invdynJ_fixb_mdp_slag_vp(qJ, qJD, qJDD, g, pkin, MDP)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(5,1),zeros(3,1),zeros(9,1),zeros(22,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RPPRR3_invdynJ_fixb_mdp_slag_vp: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RPPRR3_invdynJ_fixb_mdp_slag_vp: qJD has to be [5x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [5 1]), ...
  'S5RPPRR3_invdynJ_fixb_mdp_slag_vp: qJDD has to be [5x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S5RPPRR3_invdynJ_fixb_mdp_slag_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S5RPPRR3_invdynJ_fixb_mdp_slag_vp: pkin has to be [9x1] (double)');
assert(isreal(MDP) && all(size(MDP) == [22 1]), ...
  'S5RPPRR3_invdynJ_fixb_mdp_slag_vp: MDP has to be [22x1] (double)'); 

%% Symbolic Calculation
% From invdyn_joint_fixb_mdp_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-05 17:42:12
% EndTime: 2019-12-05 17:42:17
% DurationCPUTime: 1.76s
% Computational Cost: add. (1262->237), mult. (2656->311), div. (0->0), fcn. (2041->16), ass. (0->118)
t260 = sin(pkin(9));
t265 = sin(qJ(4));
t262 = cos(pkin(9));
t268 = cos(qJ(4));
t308 = t268 * t262;
t226 = t260 * t265 - t308;
t218 = t226 * qJD(1);
t267 = cos(qJ(5));
t227 = t260 * t268 + t262 * t265;
t219 = t227 * qJD(1);
t264 = sin(qJ(5));
t310 = t219 * t264;
t181 = t267 * t218 + t310;
t258 = qJD(4) + qJD(5);
t312 = t181 * t258;
t259 = qJ(1) + pkin(8);
t250 = sin(t259);
t252 = cos(t259);
t287 = -g(2) * t250 + g(3) * t252;
t261 = sin(pkin(8));
t239 = pkin(1) * t261 + qJ(3);
t232 = t239 * qJD(1);
t248 = t262 * qJD(2);
t200 = t248 + (-pkin(6) * qJD(1) - t232) * t260;
t207 = t260 * qJD(2) + t262 * t232;
t304 = qJD(1) * t262;
t201 = pkin(6) * t304 + t207;
t279 = -t200 * t265 - t201 * t268;
t172 = -pkin(7) * t218 - t279;
t263 = cos(pkin(8));
t243 = -pkin(1) * t263 - pkin(2);
t231 = -pkin(3) * t262 + t243;
t216 = t231 * qJD(1) + qJD(3);
t188 = pkin(4) * t218 + t216;
t257 = pkin(9) + qJ(4);
t253 = qJ(5) + t257;
t241 = sin(t253);
t242 = cos(t253);
t301 = qJD(5) * t264;
t325 = g(1) * t241 + t172 * t301 + t188 * t181 + t287 * t242;
t254 = qJDD(4) + qJDD(5);
t277 = -t218 * t264 + t267 * t219;
t324 = t254 * MDP(20) + t181 * t277 * MDP(16) + (-t181 ^ 2 + t277 ^ 2) * MDP(17);
t313 = t277 * t258;
t322 = t268 * t200 - t201 * t265;
t317 = pkin(6) + t239;
t222 = t317 * t260;
t223 = t317 * t262;
t307 = -t265 * t222 + t268 * t223;
t303 = qJD(1) * t265;
t292 = t260 * t303;
t295 = qJDD(1) * t268;
t296 = qJDD(1) * t265;
t302 = qJD(4) * t268;
t293 = t260 * t295 + t262 * t296 + t302 * t304;
t189 = -qJD(4) * t292 + t293;
t224 = qJD(1) * qJD(3) + t239 * qJDD(1);
t246 = t262 * qJDD(2);
t316 = pkin(6) * qJDD(1);
t196 = t246 + (-t224 - t316) * t260;
t203 = t260 * qJDD(2) + t262 * t224;
t197 = t262 * t316 + t203;
t290 = t268 * t196 - t265 * t197;
t164 = qJDD(4) * pkin(4) - pkin(7) * t189 + t279 * qJD(4) + t290;
t221 = t227 * qJD(4);
t236 = t262 * t295;
t285 = -t260 * t296 + t236;
t190 = qJD(1) * t221 - t285;
t280 = t265 * t196 + t268 * t197;
t165 = -pkin(7) * t190 + qJD(4) * t322 + t280;
t321 = -g(1) * t242 + t267 * t164 - t264 * t165 - t188 * t277 + t287 * t241;
t291 = t189 * t264 + t267 * t190;
t168 = t277 * qJD(5) + t291;
t318 = t221 * pkin(4);
t171 = -pkin(7) * t219 + t322;
t170 = qJD(4) * pkin(4) + t171;
t315 = t170 * t267;
t314 = t172 * t267;
t309 = t262 * MDP(5);
t306 = t260 ^ 2 + t262 ^ 2;
t305 = qJD(1) * t260;
t300 = qJD(5) * t267;
t297 = qJDD(1) * t243;
t294 = t267 * t189 - t264 * t190 - t218 * t300;
t289 = -t268 * t222 - t223 * t265;
t288 = g(2) * t252 + g(3) * t250;
t266 = sin(qJ(1));
t269 = cos(qJ(1));
t286 = g(2) * t269 + g(3) * t266;
t284 = -t170 * t264 - t314;
t193 = t267 * t226 + t227 * t264;
t220 = t226 * qJD(4);
t173 = -t193 * qJD(5) - t220 * t267 - t221 * t264;
t194 = -t226 * t264 + t227 * t267;
t283 = t173 * t258 + t194 * t254;
t178 = -pkin(7) * t227 + t289;
t179 = -pkin(7) * t226 + t307;
t282 = t178 * t267 - t179 * t264;
t281 = t178 * t264 + t179 * t267;
t202 = -t224 * t260 + t246;
t278 = -t202 * t260 + t203 * t262;
t167 = -t219 * t301 + t294;
t214 = t231 * qJDD(1) + qJDD(3);
t274 = -t222 * t302 + qJD(3) * t308 + (-qJD(3) * t260 - qJD(4) * t223) * t265;
t272 = -t227 * qJD(3) - qJD(4) * t307;
t251 = cos(t257);
t249 = sin(t257);
t230 = qJDD(3) + t297;
t206 = -t232 * t260 + t248;
t199 = pkin(4) * t226 + t231;
t192 = -qJD(4) * t221 - qJDD(4) * t226;
t191 = -qJD(4) * t220 + qJDD(4) * t227;
t177 = pkin(4) * t190 + t214;
t176 = pkin(7) * t220 + t272;
t175 = -pkin(7) * t221 + t274;
t174 = t194 * qJD(5) - t220 * t264 + t267 * t221;
t166 = -t174 * t258 - t193 * t254;
t1 = [qJDD(1) * MDP(1) + t286 * MDP(2) + (-g(2) * t266 + g(3) * t269) * MDP(3) + (t286 + (t261 ^ 2 + t263 ^ 2) * qJDD(1) * pkin(1)) * pkin(1) * MDP(4) + (t224 * t306 + t278 - t287) * MDP(7) + (t230 * t243 - g(2) * (-pkin(1) * t269 - pkin(2) * t252 - qJ(3) * t250) - g(3) * (-pkin(1) * t266 - pkin(2) * t250 + qJ(3) * t252) + t278 * t239 + (-t206 * t260 + t207 * t262) * qJD(3)) * MDP(8) + (t189 * t227 - t219 * t220) * MDP(9) + (-t189 * t226 - t190 * t227 + t218 * t220 - t219 * t221) * MDP(10) + t191 * MDP(11) + t192 * MDP(12) + (t272 * qJD(4) + t289 * qJDD(4) + t231 * t190 + t214 * t226 + t216 * t221 + t288 * t251) * MDP(14) + (-t274 * qJD(4) - t307 * qJDD(4) + t231 * t189 + t214 * t227 - t216 * t220 - t288 * t249) * MDP(15) + (t167 * t194 + t173 * t277) * MDP(16) + (-t167 * t193 - t168 * t194 - t173 * t181 - t174 * t277) * MDP(17) + t283 * MDP(18) + t166 * MDP(19) + (t181 * t318 + t199 * t168 + t177 * t193 + t188 * t174 + (-t281 * qJD(5) - t175 * t264 + t176 * t267) * t258 + t282 * t254 + t288 * t242) * MDP(21) + (t277 * t318 + t199 * t167 + t177 * t194 + t188 * t173 - (t282 * qJD(5) + t175 * t267 + t176 * t264) * t258 - t281 * t254 - t288 * t241) * MDP(22) + (-t260 * MDP(6) + t309) * (-t230 + t288 - t297); (qJDD(2) - g(1)) * MDP(4) + (t202 * t262 + t203 * t260 - g(1)) * MDP(8) + t192 * MDP(14) - t191 * MDP(15) + t166 * MDP(21) - t283 * MDP(22); (t206 * t305 - t207 * t304 + qJDD(3) - t288) * MDP(8) - t236 * MDP(14) + t293 * MDP(15) + (t168 + t313) * MDP(21) + (t167 - t312) * MDP(22) - t306 * MDP(7) * qJD(1) ^ 2 + (-t309 + t243 * MDP(8) + (MDP(14) * t265 + MDP(6)) * t260) * qJDD(1) + ((t262 * t303 + t268 * t305 + t219) * MDP(14) + (-t218 - t292) * MDP(15)) * qJD(4); t219 * t218 * MDP(9) + (-t218 ^ 2 + t219 ^ 2) * MDP(10) + (t293 + (t218 - t292) * qJD(4)) * MDP(11) + t285 * MDP(12) + qJDD(4) * MDP(13) + (-g(1) * t251 - t216 * t219 + t287 * t249 + t290) * MDP(14) + (g(1) * t249 + t216 * t218 + t287 * t251 - t280) * MDP(15) + (t167 + t312) * MDP(18) + (-t168 + t313) * MDP(19) + (-(-t171 * t264 - t314) * t258 + t284 * qJD(5) + (-t181 * t219 + t267 * t254 - t258 * t301) * pkin(4) + t321) * MDP(21) + ((-t172 * t258 - t164) * t264 + (-qJD(5) * t170 + t171 * t258 - t165) * t267 + (-t219 * t277 - t264 * t254 - t258 * t300) * pkin(4) + t325) * MDP(22) + t324; (t294 + t312) * MDP(18) + (-t291 + t313) * MDP(19) + (-t284 * t258 + t321) * MDP(21) + (-t267 * t165 - t264 * t164 + (-t172 * t264 + t315) * t258 + t325) * MDP(22) + (-MDP(18) * t310 - t277 * MDP(19) + t284 * MDP(21) - MDP(22) * t315) * qJD(5) + t324;];
tau = t1;
