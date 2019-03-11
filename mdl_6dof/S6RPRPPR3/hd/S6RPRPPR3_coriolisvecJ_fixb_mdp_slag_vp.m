% Calculate minimal parameter regressor of Coriolis joint torque vector for
% S6RPRPPR3
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% qJD [6x1]
%   Generalized joint velocities
% pkin [9x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d3,d6,theta2]';
% MDP [26x1]
%   Minimal dynamic parameter vector (fixed base model)
%   see S6RPRPPR3_convert_par2_MPV_fixb.m
% 
% Output:
% tauc [6x1]
%   joint torques required to compensate Coriolis and centrifugal load

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-09 02:45
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauc = S6RPRPPR3_coriolisvecJ_fixb_mdp_slag_vp(qJ, qJD, pkin, MDP)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(9,1),zeros(26,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPRPPR3_coriolisvecJ_fixb_mdp_slag_vp: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RPRPPR3_coriolisvecJ_fixb_mdp_slag_vp: qJD has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S6RPRPPR3_coriolisvecJ_fixb_mdp_slag_vp: pkin has to be [9x1] (double)');
assert(isreal(MDP) && all(size(MDP) == [26 1]), ...
  'S6RPRPPR3_coriolisvecJ_fixb_mdp_slag_vp: MDP has to be [26x1] (double)'); 

%% Symbolic Calculation
% From coriolisvec_joint_fixb_mdp_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 02:45:28
% EndTime: 2019-03-09 02:45:35
% DurationCPUTime: 2.58s
% Computational Cost: add. (1097->289), mult. (2359->380), div. (0->0), fcn. (1219->6), ass. (0->130)
t229 = -cos(pkin(9)) * pkin(1) - pkin(2);
t219 = qJD(1) * t229;
t228 = sin(pkin(9)) * pkin(1) + pkin(7);
t218 = t228 * qJD(1);
t261 = sin(qJ(3));
t213 = t261 * t218;
t263 = cos(qJ(3));
t323 = -t263 * qJD(2) + t213;
t343 = qJD(4) + t323;
t194 = -qJD(3) * pkin(3) + t343;
t264 = -pkin(3) - pkin(4);
t294 = qJD(3) * t264;
t317 = qJD(1) * t261;
t192 = -qJ(5) * t317 + t323;
t303 = -qJD(4) - t192;
t181 = t294 - t303;
t204 = t261 * qJD(2) + t263 * t218;
t253 = qJD(3) * qJ(4);
t198 = t253 + t204;
t300 = qJD(2) * qJD(3);
t310 = qJD(3) * t263;
t201 = t218 * t310 + t261 * t300;
t342 = (-qJD(3) * t198 + t201) * MDP(15);
t260 = sin(qJ(6));
t262 = cos(qJ(6));
t341 = MDP(25) * t260 + MDP(26) * t262;
t254 = t261 ^ 2;
t255 = t263 ^ 2;
t340 = (t254 - t255) * MDP(6);
t326 = -qJ(5) + t228;
t316 = qJD(1) * t263;
t193 = -qJ(5) * t316 + t204;
t186 = -t193 - t253;
t338 = (-t323 + t213) * qJD(3);
t199 = -pkin(3) * t316 - qJ(4) * t317 + t219;
t184 = pkin(4) * t316 + qJD(5) - t199;
t278 = pkin(5) * t261 + pkin(8) * t263;
t176 = qJD(1) * t278 + t184;
t309 = qJD(5) * t261;
t179 = (-qJ(5) * t310 - t309) * qJD(1) + t201;
t182 = qJD(3) * pkin(5) - t186;
t208 = -t263 * pkin(3) - t261 * qJ(4) + t229;
t202 = t263 * pkin(4) - t208;
t188 = t202 + t278;
t197 = t310 * t326 - t309;
t227 = qJD(6) + t317;
t337 = -(qJD(6) * t188 + t197) * t227 + (qJD(3) * t182 - qJD(6) * t176 - t179) * t261;
t232 = t263 * t300;
t251 = qJD(3) * qJD(4);
t311 = qJD(3) * t261;
t189 = -t218 * t311 + t232 + t251;
t301 = qJD(1) * qJD(3);
t290 = t261 * t301;
t225 = qJ(5) * t290;
t308 = qJD(5) * t263;
t177 = qJD(1) * t308 - t189 - t225;
t336 = t177 * t260;
t335 = t177 * t262;
t304 = t262 * qJD(3);
t289 = t262 * t301;
t306 = qJD(6) * t260;
t291 = t263 * t306;
t322 = qJD(1) * t291 + t261 * t289;
t190 = -qJD(6) * t304 + t322;
t334 = t190 * t260;
t215 = t260 * t316 - t304;
t333 = t215 * t227;
t332 = t215 * t263;
t312 = qJD(3) * t260;
t216 = t262 * t316 + t312;
t331 = t216 * t227;
t330 = t227 * t261;
t329 = t227 * t262;
t265 = qJD(3) ^ 2;
t328 = t261 * t265;
t327 = t263 * t265;
t325 = t190 * t261 - t216 * t310;
t305 = qJD(6) * t262;
t292 = t227 * t305;
t324 = t260 * t255 * t301 + t263 * t292;
t241 = t261 * qJD(4);
t288 = t263 * t301;
t321 = qJ(4) * t288 + qJD(1) * t241;
t320 = t232 + 0.2e1 * t251;
t319 = qJ(4) * t310 + t241;
t315 = qJD(3) * t186;
t313 = qJD(3) * t215;
t252 = -pkin(8) + t264;
t178 = qJD(3) * t252 - t303;
t307 = qJD(6) * t178;
t302 = -qJD(5) - t184;
t299 = -MDP(12) + MDP(17);
t298 = MDP(14) + MDP(16);
t297 = t260 * t330;
t296 = t261 * t329;
t266 = qJD(1) ^ 2;
t295 = t261 * t263 * t266;
t293 = t227 * t306;
t287 = MDP(24) * t316;
t284 = qJD(1) * t202 + t184;
t280 = t261 * t294;
t185 = qJD(1) * t280 + t321;
t195 = t280 + t319;
t283 = qJD(1) * t195 + t185;
t282 = qJD(1) * t208 + t199;
t277 = qJD(3) * t204 - t201;
t172 = t176 * t262 - t178 * t260;
t173 = t176 * t260 + t178 * t262;
t209 = t326 * t261;
t270 = pkin(5) * t263 + t252 * t261;
t268 = t270 * qJD(3);
t276 = (-qJD(6) * t209 + t268 + t319) * t227;
t275 = 0.2e1 * qJD(3) * t219;
t274 = -MDP(25) * t262 + MDP(26) * t260;
t200 = pkin(3) * t290 - t321;
t207 = pkin(3) * t311 - t319;
t271 = -qJD(1) * t207 - t228 * t265 - t200;
t267 = t189 * t263 + t201 * t261 + (t194 * t263 - t198 * t261) * qJD(3);
t259 = qJ(4) + pkin(5);
t236 = qJ(4) * t316;
t221 = t260 * t290;
t217 = pkin(3) * t317 - t236;
t210 = t326 * t263;
t205 = t264 * t317 + t236;
t196 = t326 * t311 + t308;
t191 = t216 * qJD(6) - t221;
t183 = qJD(1) * t270 + t236;
t175 = qJD(1) * t268 + t321;
t174 = t262 * t175;
t1 = [0.2e1 * t261 * MDP(5) * t288 - 0.2e1 * t301 * t340 + MDP(7) * t327 - MDP(8) * t328 + (-t228 * t327 + t261 * t275) * MDP(10) + (t228 * t328 + t263 * t275) * MDP(11) + (t263 * t271 + t282 * t311) * MDP(12) + t267 * MDP(13) + (t261 * t271 - t282 * t310) * MDP(14) + (t199 * t207 + t200 * t208 + t228 * t267) * MDP(15) + (t283 * t261 + (t263 * t284 - t196) * qJD(3)) * MDP(16) + (-t283 * t263 + (t261 * t284 + t197) * qJD(3)) * MDP(17) + (t177 * t263 - t179 * t261 + (-t181 * t263 - t186 * t261) * qJD(3) + (t196 * t263 - t197 * t261 + (-t209 * t263 + t210 * t261) * qJD(3)) * qJD(1)) * MDP(18) + (-t177 * t210 + t179 * t209 + t181 * t197 + t184 * t195 + t185 * t202 + t186 * t196) * MDP(19) + (-t190 * t262 * t263 + (-t261 * t304 - t291) * t216) * MDP(20) + ((t215 * t262 + t216 * t260) * t311 + (t334 - t191 * t262 + (t215 * t260 - t216 * t262) * qJD(6)) * t263) * MDP(21) + (t227 * t291 + (-qJD(1) * t255 + t330) * t304 + t325) * MDP(22) + (t191 * t261 + (-t297 + t332) * qJD(3) + t324) * MDP(23) + (t227 + t317) * MDP(24) * t310 + (t174 * t261 - t210 * t191 + t196 * t215 + (-t261 * t307 + t276) * t262 + t337 * t260 + (-t182 * t305 + t336 + ((t188 * t262 - t209 * t260) * qJD(1) + t172) * qJD(3)) * t263) * MDP(25) + (t210 * t190 + t196 * t216 + (-t276 - (t175 - t307) * t261) * t260 + t337 * t262 + (t182 * t306 + t335 + (-(t188 * t260 + t209 * t262) * qJD(1) - t173) * qJD(3)) * t263) * MDP(26); t324 * MDP(25) + (t255 * t289 + t325) * MDP(26) + (-t342 + (-t179 - t315) * MDP(19) - MDP(25) * t313 - MDP(26) * t293 + (-MDP(11) + t298) * t265) * t263 + ((qJD(3) * t194 + t189) * MDP(15) + (qJD(3) * t181 - t177) * MDP(19) + (-t227 * t312 - t191) * MDP(25) - t227 * MDP(26) * t304 + (-MDP(10) + t299) * t265) * t261; -MDP(5) * t295 + t266 * t340 + (-t219 * t317 + t277) * MDP(10) + (-t219 * t316 - t232 + t338) * MDP(11) + t277 * MDP(12) + (t320 - t338) * MDP(14) + (-pkin(3) * t201 + qJ(4) * t189 - t194 * t204 + t343 * t198 - t199 * t217) * MDP(15) + (t225 + (t192 - t213) * qJD(3) + t320) * MDP(16) + (-qJD(3) * t193 + t201) * MDP(17) + (-qJ(4) * t177 + t179 * t264 - t181 * t193 - t184 * t205 + t186 * t303) * MDP(19) + (t216 * t329 - t334) * MDP(20) + ((-t190 - t333) * t262 + (-t191 - t331) * t260) * MDP(21) - t292 * MDP(22) + t293 * MDP(23) - t227 * t287 + (-t259 * t191 - t335 - (t183 * t262 - t193 * t260) * t227 + t303 * t215 + (-t182 * t260 - t252 * t329) * qJD(6)) * MDP(25) + (t259 * t190 + t336 + (t183 * t260 + t193 * t262) * t227 + t303 * t216 + (t227 * t252 * t260 - t182 * t262) * qJD(6)) * MDP(26) + (-t296 * MDP(22) + t297 * MDP(23) + t341 * (-t182 * t261 - t252 * t310) + (-t199 * MDP(12) + t217 * MDP(14) - t205 * MDP(16) + t302 * MDP(17)) * t261 + (t217 * MDP(12) + t199 * MDP(14) + t302 * MDP(16) + (-qJ(5) * qJD(3) + t205) * MDP(17) + (t216 - t312) * MDP(22) + (-t215 - t304) * MDP(23) - t172 * MDP(25) + t173 * MDP(26)) * t263) * qJD(1); t342 + (t315 + t201) * MDP(19) + (-t292 + t313) * MDP(25) + (qJD(3) * t216 + t293) * MDP(26) + t299 * t295 + t298 * (-t254 * t266 - t265) + ((-MDP(19) * qJ(5) - t341) * t310 + (t199 * MDP(15) + MDP(19) * t302 + t227 * t274) * t261) * qJD(1); t321 * MDP(19) + (-t254 - t255) * MDP(18) * t266 - t341 * t227 * qJD(6) + ((t181 * t261 - t186 * t263) * MDP(19) + (-t297 - t332) * MDP(25) + (-t216 * t263 - t296) * MDP(26) + ((0.2e1 * MDP(16) - t274) * t263 + (MDP(19) * t264 + 0.2e1 * MDP(17)) * t261) * qJD(3)) * qJD(1); t216 * t215 * MDP(20) + (-t215 ^ 2 + t216 ^ 2) * MDP(21) + (t322 - t333) * MDP(22) + (-t221 - t331) * MDP(23) + qJD(3) * t287 + (t173 * t227 - t179 * t260 + t182 * t216 + t174) * MDP(25) + (t172 * t227 - t175 * t260 - t179 * t262 - t182 * t215) * MDP(26) + (-MDP(22) * t304 + MDP(23) * t216 - MDP(25) * t173 - MDP(26) * t172) * qJD(6);];
tauc  = t1;
