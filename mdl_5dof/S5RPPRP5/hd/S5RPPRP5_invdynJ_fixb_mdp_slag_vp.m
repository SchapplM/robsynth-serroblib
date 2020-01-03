% Calculate vector of inverse dynamics joint torques for
% S5RPPRP5
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
% pkin [7x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d1,d4,theta2]';
% MDP [22x1]
%   Minimal dynamic parameter vector (fixed base model)
%   see S5RPPRP5_convert_par2_MPV_fixb.m
% 
% Output:
% tau [5x1]
%   joint torques of inverse dynamics (contains inertial, gravitational Coriolis and centrifugal forces)

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-31 17:54
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tau = S5RPPRP5_invdynJ_fixb_mdp_slag_vp(qJ, qJD, qJDD, g, pkin, MDP)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(5,1),zeros(3,1),zeros(7,1),zeros(22,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RPPRP5_invdynJ_fixb_mdp_slag_vp: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RPPRP5_invdynJ_fixb_mdp_slag_vp: qJD has to be [5x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [5 1]), ...
  'S5RPPRP5_invdynJ_fixb_mdp_slag_vp: qJDD has to be [5x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S5RPPRP5_invdynJ_fixb_mdp_slag_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [7 1]), ...
  'S5RPPRP5_invdynJ_fixb_mdp_slag_vp: pkin has to be [7x1] (double)');
assert(isreal(MDP) && all(size(MDP) == [22 1]), ...
  'S5RPPRP5_invdynJ_fixb_mdp_slag_vp: MDP has to be [22x1] (double)'); 

%% Symbolic Calculation
% From invdyn_joint_fixb_mdp_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 17:53:48
% EndTime: 2019-12-31 17:53:51
% DurationCPUTime: 1.89s
% Computational Cost: add. (1002->271), mult. (2233->328), div. (0->0), fcn. (1547->6), ass. (0->117)
t291 = MDP(18) - MDP(21);
t253 = sin(pkin(7));
t254 = cos(pkin(7));
t210 = -t254 * pkin(2) - t253 * qJ(3) - pkin(1);
t204 = t254 * pkin(3) - t210;
t257 = sin(qJ(1));
t259 = cos(qJ(1));
t307 = g(1) * t259 + g(2) * t257;
t256 = sin(qJ(4));
t258 = cos(qJ(4));
t207 = t253 * t256 + t254 * t258;
t200 = t207 * qJD(1);
t292 = MDP(17) + MDP(19);
t297 = qJD(1) * qJD(3);
t229 = t253 * t297;
t295 = qJDD(1) * t254;
t296 = qJDD(1) * t253;
t322 = -pkin(2) * t295 - qJ(3) * t296 - t229;
t305 = qJD(1) * t254;
t286 = t256 * t305;
t214 = qJD(4) * t286;
t321 = -t207 * qJDD(1) + t214;
t320 = t200 ^ 2;
t306 = qJD(1) * t253;
t288 = t258 * t306;
t202 = -t286 + t288;
t319 = t202 ^ 2;
t318 = g(1) * t257;
t248 = g(2) * t259;
t317 = -pkin(6) + qJ(2);
t252 = qJDD(1) * pkin(1);
t316 = qJDD(4) * pkin(4);
t315 = t200 * t202;
t213 = t317 * t254;
t209 = qJD(1) * t213;
t314 = t209 * t256;
t313 = t253 * t258;
t312 = t253 * t259;
t311 = t254 * t257;
t310 = t254 * t259;
t251 = t254 ^ 2;
t298 = qJD(1) * qJD(2);
t290 = 0.2e1 * t298;
t224 = t251 * t290;
t294 = qJDD(1) * qJ(2) ^ 2;
t309 = qJ(2) * t224 + t251 * t294;
t308 = t259 * pkin(1) + t257 * qJ(2);
t304 = qJD(3) * t253;
t217 = qJ(2) * t306 + qJD(3);
t206 = -pkin(6) * t306 + t217;
t177 = t206 * t256 + t209 * t258;
t303 = qJD(4) * t177;
t302 = qJD(4) * t256;
t301 = qJD(4) * t258;
t176 = t206 * t258 - t314;
t300 = qJD(5) - t176;
t299 = qJ(2) * qJDD(1);
t293 = qJDD(4) * qJ(5);
t239 = qJDD(2) - t252;
t205 = qJ(2) * t296 + t253 * t298 + qJDD(3);
t188 = -pkin(6) * t296 + t205;
t189 = (t317 * qJDD(1) + t298) * t254;
t289 = t256 * t188 + t258 * t189 + t206 * t301;
t287 = t253 * t301;
t244 = t259 * qJ(2);
t285 = -pkin(1) * t257 + t244;
t284 = -t248 + t318;
t250 = t253 ^ 2;
t261 = qJD(1) ^ 2;
t283 = (-t250 - t251) * t261;
t282 = t176 + t314;
t281 = -t258 * t188 + t256 * t189 + t206 * t302 + t209 * t301;
t280 = 0.2e1 * t251 * t299 + t224 - t307;
t279 = pkin(2) * t310 + qJ(3) * t312 + t308;
t197 = -qJD(1) * pkin(1) - pkin(2) * t305 - qJ(3) * t306 + qJD(2);
t208 = -t254 * t256 + t313;
t276 = -pkin(4) * t207 + qJ(5) * t208;
t275 = t256 * t295 - t258 * t296;
t185 = pkin(3) * t305 - t197;
t212 = t317 * t253;
t274 = t212 * t258 - t213 * t256;
t181 = t212 * t256 + t213 * t258;
t184 = t239 + t322;
t272 = -t284 + t239;
t271 = -t239 + t252 - t248;
t270 = -qJDD(1) * t210 - t184 - t248;
t179 = pkin(3) * t295 - t184;
t198 = t207 * qJD(4);
t269 = (t298 + t299) * t250;
t193 = t207 * t257;
t195 = t207 * t259;
t268 = g(1) * t195 + g(2) * t193 + g(3) * t208 - t289;
t168 = t207 * qJD(2) + t274 * qJD(4);
t192 = t256 * t311 - t257 * t313;
t194 = t256 * t310 - t258 * t312;
t267 = g(1) * t192 - g(2) * t194 + qJD(4) * t168 + qJDD(4) * t181;
t169 = -t208 * qJD(2) + t181 * qJD(4);
t266 = g(1) * t193 - g(2) * t195 - qJD(4) * t169 + qJDD(4) * t274;
t265 = g(1) * t194 + g(2) * t192 + g(3) * t207 - t281;
t174 = qJD(1) * t198 + t275;
t175 = qJD(1) * t287 - t321;
t264 = -pkin(4) * t175 - qJ(5) * t174 - t179;
t164 = pkin(4) * t200 - qJ(5) * t202 + t185;
t263 = t164 * t202 + qJDD(5) - t265;
t260 = qJD(4) ^ 2;
t245 = g(3) * t254;
t226 = g(1) * t311;
t199 = -t254 * t302 + t287;
t173 = pkin(4) * t202 + qJ(5) * t200;
t172 = qJD(4) * qJ(5) + t177;
t171 = -qJD(4) * pkin(4) + t300;
t170 = -t276 + t204;
t165 = pkin(4) * t199 + qJ(5) * t198 - qJD(5) * t208 + t304;
t163 = qJDD(5) + t281 - t316;
t162 = t293 + (qJD(5) - t314) * qJD(4) + t289;
t161 = -qJD(5) * t202 - t264;
t1 = [qJDD(1) * MDP(1) + t284 * MDP(2) + t307 * MDP(3) + (t271 * t254 + t226) * MDP(4) + (-t271 - t318) * t253 * MDP(5) + (0.2e1 * t269 + t280) * MDP(6) + (-t239 * pkin(1) - g(1) * t285 - g(2) * t308 + (qJ(2) * t290 + t294) * t250 + t309) * MDP(7) + (t226 + (t270 + t229) * t254) * MDP(8) + (t205 * t253 + t269 + t280) * MDP(9) + (t250 * t297 + (t270 + t318) * t253) * MDP(10) + (t184 * t210 - g(1) * (-pkin(2) * t311 + t285) - g(2) * t279 + (qJ(2) * t205 + qJ(3) * t318 + qJD(2) * t217 - qJD(3) * t197) * t253 + t309) * MDP(11) + (-t174 * t208 - t198 * t202) * MDP(12) + (t174 * t207 - t175 * t208 + t198 * t200 - t199 * t202) * MDP(13) + (-qJD(4) * t198 + qJDD(4) * t208) * MDP(14) + (-qJD(4) * t199 - qJDD(4) * t207) * MDP(15) + (t175 * t204 + t179 * t207 + t185 * t199 + t200 * t304 + t266) * MDP(17) + (-t174 * t204 + t179 * t208 - t185 * t198 + t202 * t304 - t267) * MDP(18) + (t161 * t207 + t164 * t199 + t165 * t200 + t170 * t175 + t266) * MDP(19) + (-t162 * t207 + t163 * t208 - t168 * t200 + t169 * t202 - t171 * t198 - t172 * t199 + t174 * t274 - t175 * t181 + t307) * MDP(20) + (-t161 * t208 + t164 * t198 - t165 * t202 + t170 * t174 + t267) * MDP(21) + (t162 * t181 + t172 * t168 + t161 * t170 + t164 * t165 - t163 * t274 + t171 * t169 - g(1) * (-pkin(4) * t193 - pkin(6) * t259 - qJ(5) * t192 + t244) - g(2) * (pkin(3) * t310 + pkin(4) * t195 + qJ(5) * t194 + t279) + (g(2) * pkin(6) + g(1) * t204) * t257) * MDP(22); (qJ(2) * t283 + t272) * MDP(7) + (-qJ(2) * t251 * t261 - t217 * t306 + t272 + t322) * MDP(11) + (t319 + t320) * MDP(20) + (-t172 * t200 + (qJD(5) + t171) * t202 + t264 - t284) * MDP(22) + t292 * t214 + (MDP(6) + MDP(9)) * t283 + t291 * (0.2e1 * qJD(4) * t200 + t275) + ((-t292 * t258 - MDP(4) - MDP(8)) * t254 + (-t292 * t256 - MDP(10) + MDP(5)) * t253) * qJDD(1) + t292 * qJD(4) * (-t202 - t288); -t250 * t261 * MDP(10) + (t245 + t205) * MDP(11) + (t174 * t258 - t175 * t256) * MDP(20) + (t162 * t256 - t163 * t258 + t245) * MDP(22) + t292 * (qJDD(4) * t258 - t200 * t306 - t256 * t260) + ((-t200 * t258 + t202 * t256) * MDP(20) + (t171 * t256 + t172 * t258) * MDP(22)) * qJD(4) + (-t261 * t254 * MDP(8) + qJDD(1) * MDP(9) + (t197 * MDP(11) - t164 * MDP(22) - t291 * t202) * qJD(1) - (MDP(11) + MDP(22)) * t307) * t253 - t291 * (qJDD(4) * t256 + t258 * t260); MDP(12) * t315 + (t319 - t320) * MDP(13) - t275 * MDP(14) + ((t202 - t288) * qJD(4) + t321) * MDP(15) + qJDD(4) * MDP(16) + (-t185 * t202 + t265 + t303) * MDP(17) + (t282 * qJD(4) + t185 * t200 + t268) * MDP(18) + (-t173 * t200 - t263 + t303 + 0.2e1 * t316) * MDP(19) + (pkin(4) * t174 - qJ(5) * t175 + (t172 - t177) * t202 + (t171 - t300) * t200) * MDP(20) + (0.2e1 * t293 - t164 * t200 + t173 * t202 + (0.2e1 * qJD(5) - t282) * qJD(4) - t268) * MDP(21) + (t162 * qJ(5) - t163 * pkin(4) - t164 * t173 - t171 * t177 - g(1) * (-pkin(4) * t194 + qJ(5) * t195) - g(2) * (-pkin(4) * t192 + qJ(5) * t193) - g(3) * t276 + t300 * t172) * MDP(22); (-qJDD(4) + t315) * MDP(19) - t275 * MDP(20) + (-t260 - t319) * MDP(21) + (-qJD(4) * t172 + t263 - t316) * MDP(22);];
tau = t1;
