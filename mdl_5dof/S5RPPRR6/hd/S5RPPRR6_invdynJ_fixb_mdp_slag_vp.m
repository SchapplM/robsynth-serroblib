% Calculate vector of inverse dynamics joint torques for
% S5RPPRR6
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
%   see S5RPPRR6_convert_par2_MPV_fixb.m
% 
% Output:
% tau [5x1]
%   joint torques of inverse dynamics (contains inertial, gravitational Coriolis and centrifugal forces)

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-31 17:58
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tau = S5RPPRR6_invdynJ_fixb_mdp_slag_vp(qJ, qJD, qJDD, g, pkin, MDP)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(5,1),zeros(3,1),zeros(9,1),zeros(22,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RPPRR6_invdynJ_fixb_mdp_slag_vp: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RPPRR6_invdynJ_fixb_mdp_slag_vp: qJD has to be [5x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [5 1]), ...
  'S5RPPRR6_invdynJ_fixb_mdp_slag_vp: qJDD has to be [5x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S5RPPRR6_invdynJ_fixb_mdp_slag_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S5RPPRR6_invdynJ_fixb_mdp_slag_vp: pkin has to be [9x1] (double)');
assert(isreal(MDP) && all(size(MDP) == [22 1]), ...
  'S5RPPRR6_invdynJ_fixb_mdp_slag_vp: MDP has to be [22x1] (double)'); 

%% Symbolic Calculation
% From invdyn_joint_fixb_mdp_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 17:58:12
% EndTime: 2019-12-31 17:58:16
% DurationCPUTime: 2.20s
% Computational Cost: add. (1403->266), mult. (2997->350), div. (0->0), fcn. (2231->14), ass. (0->127)
t243 = sin(pkin(9));
t248 = sin(qJ(4));
t289 = qJD(1) * t248;
t276 = t243 * t289;
t245 = cos(pkin(9));
t251 = cos(qJ(4));
t288 = qJD(1) * t251;
t223 = t245 * t288;
t281 = qJDD(1) * t251;
t282 = qJDD(1) * t248;
t277 = qJD(4) * t223 + t243 * t281 + t245 * t282;
t183 = -qJD(4) * t276 + t277;
t313 = qJD(4) * qJD(5) + t183;
t244 = sin(pkin(8));
t224 = pkin(1) * t244 + qJ(3);
t218 = t224 * qJD(1);
t231 = t245 * qJD(2);
t189 = t231 + (-pkin(6) * qJD(1) - t218) * t243;
t201 = qJD(2) * t243 + t218 * t245;
t290 = qJD(1) * t245;
t190 = pkin(6) * t290 + t201;
t169 = t189 * t248 + t190 * t251;
t212 = qJD(1) * qJD(3) + qJDD(1) * t224;
t229 = t245 * qJDD(2);
t306 = pkin(6) * qJDD(1);
t187 = t229 + (-t212 - t306) * t243;
t195 = qJDD(2) * t243 + t212 * t245;
t188 = t245 * t306 + t195;
t265 = -t187 * t251 + t188 * t248;
t158 = -qJDD(4) * pkin(4) + qJD(4) * t169 + t265;
t206 = t223 - t276;
t203 = qJD(5) - t206;
t215 = t243 * t251 + t245 * t248;
t207 = t215 * qJD(1);
t241 = pkin(9) + qJ(4);
t232 = sin(t241);
t234 = cos(t241);
t242 = qJ(1) + pkin(8);
t233 = sin(t242);
t235 = cos(t242);
t270 = g(1) * t235 + g(2) * t233;
t256 = -g(3) * t234 + t232 * t270;
t312 = t256 - (pkin(4) * t207 + pkin(7) * t203) * t203 - t158;
t168 = t189 * t251 - t190 * t248;
t166 = -qJD(4) * pkin(4) - t168;
t214 = t243 * t248 - t245 * t251;
t307 = pkin(6) + t224;
t210 = t307 * t243;
t211 = t307 * t245;
t262 = -t210 * t251 - t211 * t248;
t170 = -qJD(3) * t214 + qJD(4) * t262;
t246 = cos(pkin(8));
t226 = -pkin(1) * t246 - pkin(2);
t217 = -pkin(3) * t245 + t226;
t174 = pkin(4) * t214 - pkin(7) * t215 + t217;
t177 = -t210 * t248 + t211 * t251;
t209 = t215 * qJD(4);
t222 = t245 * t281;
t267 = -t243 * t282 + t222;
t184 = qJD(1) * t209 - t267;
t180 = qJDD(5) + t184;
t208 = t214 * qJD(4);
t264 = t187 * t248 + t188 * t251;
t157 = qJDD(4) * pkin(7) + qJD(4) * t168 + t264;
t205 = qJD(1) * t217 + qJD(3);
t172 = -pkin(4) * t206 - pkin(7) * t207 + t205;
t273 = qJD(5) * t172 + t157;
t310 = t158 * t215 - t166 * t208 - t177 * t180 - (qJD(5) * t174 + t170) * t203 - t273 * t214;
t309 = g(3) * t232;
t247 = sin(qJ(5));
t250 = cos(qJ(5));
t278 = t247 * qJDD(4) + t250 * t313;
t286 = qJD(5) * t247;
t164 = -t207 * t286 + t278;
t305 = t164 * t247;
t304 = t174 * t180;
t191 = -qJD(4) * t250 + t207 * t247;
t303 = t191 * t203;
t302 = t191 * t207;
t193 = qJD(4) * t247 + t207 * t250;
t301 = t193 * t203;
t300 = t193 * t207;
t299 = (-t218 * t243 + t231) * t243;
t298 = t233 * t247;
t297 = t233 * t250;
t296 = t235 * t247;
t295 = t235 * t250;
t294 = t245 * MDP(5);
t293 = t247 * t180;
t175 = t250 * t180;
t292 = t164 * t214 + t193 * t209;
t291 = t243 ^ 2 + t245 ^ 2;
t287 = qJD(5) * t215;
t283 = qJDD(1) * t226;
t280 = t215 * t293;
t279 = t215 * t175;
t274 = t203 * t250;
t202 = qJDD(1) * t217 + qJDD(3);
t162 = pkin(4) * t184 - pkin(7) * t183 + t202;
t167 = qJD(4) * pkin(7) + t169;
t272 = qJD(5) * t167 - t162;
t269 = g(1) * t233 - g(2) * t235;
t249 = sin(qJ(1));
t252 = cos(qJ(1));
t268 = g(1) * t249 - g(2) * t252;
t237 = t250 * qJDD(4);
t165 = qJD(5) * t193 + t183 * t247 - t237;
t266 = -t165 * t214 - t191 * t209;
t194 = -t212 * t243 + t229;
t263 = -t194 * t243 + t195 * t245;
t261 = t175 + (t206 * t247 - t286) * t203;
t260 = t208 * t247 - t250 * t287;
t259 = t208 * t250 + t215 * t286;
t255 = -pkin(7) * t180 + (t166 + t168) * t203;
t216 = qJDD(3) + t283;
t199 = t234 * t295 + t298;
t198 = -t234 * t296 + t297;
t197 = -t234 * t297 + t296;
t196 = t234 * t298 + t295;
t186 = -qJD(4) * t209 - qJDD(4) * t214;
t185 = -qJD(4) * t208 + qJDD(4) * t215;
t182 = pkin(4) * t209 + pkin(7) * t208;
t171 = qJD(3) * t215 + qJD(4) * t177;
t161 = t250 * t162;
t160 = t167 * t250 + t172 * t247;
t159 = -t167 * t247 + t172 * t250;
t1 = [qJDD(1) * MDP(1) + t268 * MDP(2) + (g(1) * t252 + g(2) * t249) * MDP(3) + (t268 + (t244 ^ 2 + t246 ^ 2) * qJDD(1) * pkin(1)) * pkin(1) * MDP(4) + (t212 * t291 + t263 - t270) * MDP(7) + (t216 * t226 - g(1) * (-pkin(1) * t249 - pkin(2) * t233 + qJ(3) * t235) - g(2) * (pkin(1) * t252 + pkin(2) * t235 + qJ(3) * t233) + t263 * t224 + (t201 * t245 - t299) * qJD(3)) * MDP(8) + (t183 * t215 - t207 * t208) * MDP(9) + (-t183 * t214 - t184 * t215 - t206 * t208 - t207 * t209) * MDP(10) + t185 * MDP(11) + t186 * MDP(12) + (-qJD(4) * t171 + qJDD(4) * t262 + t184 * t217 + t202 * t214 + t205 * t209 + t234 * t269) * MDP(14) + (-qJD(4) * t170 - qJDD(4) * t177 + t183 * t217 + t202 * t215 - t205 * t208 - t232 * t269) * MDP(15) + (t164 * t215 * t250 - t193 * t259) * MDP(16) + (-(-t191 * t250 - t193 * t247) * t208 + (-t305 - t165 * t250 + (t191 * t247 - t193 * t250) * qJD(5)) * t215) * MDP(17) + (-t203 * t259 + t279 + t292) * MDP(18) + (t203 * t260 + t266 - t280) * MDP(19) + (t180 * t214 + t203 * t209) * MDP(20) + (-g(1) * t197 - g(2) * t199 + t159 * t209 + t161 * t214 - t262 * t165 + t171 * t191 + (t304 + t182 * t203 + (t166 * t215 - t167 * t214 - t177 * t203) * qJD(5)) * t250 + t310 * t247) * MDP(21) + (-g(1) * t196 - g(2) * t198 - t160 * t209 - t262 * t164 + t171 * t193 + (-(-qJD(5) * t177 + t182) * t203 - t304 + t272 * t214 - t166 * t287) * t247 + t310 * t250) * MDP(22) + (-t243 * MDP(6) + t294) * (-t216 + t269 - t283); (qJDD(2) - g(3)) * MDP(4) + (t194 * t245 + t195 * t243 - g(3)) * MDP(8) + t186 * MDP(14) - t185 * MDP(15) + (-t266 - t280) * MDP(21) + (-t279 + t292) * MDP(22) + (MDP(21) * t260 + MDP(22) * t259) * t203; (qJD(1) * t299 - t201 * t290 + qJDD(3) - t269) * MDP(8) - t222 * MDP(14) + t277 * MDP(15) + (t261 - t302) * MDP(21) + (-t203 ^ 2 * t250 - t293 - t300) * MDP(22) - t291 * MDP(7) * qJD(1) ^ 2 + (-t294 + t226 * MDP(8) + (MDP(14) * t248 + MDP(6)) * t243) * qJDD(1) + ((t243 * t288 + t245 * t289 + t207) * MDP(14) + (t206 - t276) * MDP(15)) * qJD(4); -t206 ^ 2 * MDP(10) + ((-t206 - t276) * qJD(4) + t277) * MDP(11) + t267 * MDP(12) + qJDD(4) * MDP(13) + (t256 - t265) * MDP(14) + (-t205 * t206 + t234 * t270 - t264 + t309) * MDP(15) + (t193 * t274 + t305) * MDP(16) + ((t164 - t303) * t250 + (-t165 - t301) * t247) * MDP(17) + (t203 * t274 + t293 - t300) * MDP(18) + (t261 + t302) * MDP(19) + (-pkin(4) * t165 - t169 * t191 + t255 * t247 + t250 * t312) * MDP(21) + (-pkin(4) * t164 - t169 * t193 - t247 * t312 + t255 * t250) * MDP(22) + (MDP(10) * t207 - MDP(14) * t205 - MDP(20) * t203 - MDP(21) * t159 + MDP(22) * t160 - MDP(9) * t206) * t207; t193 * t191 * MDP(16) + (-t191 ^ 2 + t193 ^ 2) * MDP(17) + (t278 + t303) * MDP(18) + (t237 + t301) * MDP(19) + t180 * MDP(20) + (-g(1) * t198 + g(2) * t196 + t160 * t203 - t166 * t193 + t161) * MDP(21) + (g(1) * t199 - g(2) * t197 + t159 * t203 + t166 * t191) * MDP(22) + ((-t157 + t309) * MDP(22) + (-MDP(19) * t207 - MDP(21) * t167 - MDP(22) * t172) * qJD(5)) * t250 + (-qJD(5) * t207 * MDP(18) - t313 * MDP(19) + (-t273 + t309) * MDP(21) + t272 * MDP(22)) * t247;];
tau = t1;
