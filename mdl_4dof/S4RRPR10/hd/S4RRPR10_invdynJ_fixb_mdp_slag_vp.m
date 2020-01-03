% Calculate vector of inverse dynamics joint torques for
% S4RRPR10
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [4x1]
%   Generalized joint coordinates (joint angles)
% qJD [4x1]
%   Generalized joint velocities
% qJDD [4x1]
%   Generalized joint accelerations
% g [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [6x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,d1,d2,d4]';
% MDP [21x1]
%   Minimal dynamic parameter vector (fixed base model)
%   see S4RRPR10_convert_par2_MPV_fixb.m
% 
% Output:
% tau [4x1]
%   joint torques of inverse dynamics (contains inertial, gravitational Coriolis and centrifugal forces)

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-31 17:12
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tau = S4RRPR10_invdynJ_fixb_mdp_slag_vp(qJ, qJD, qJDD, g, pkin, MDP)
%% Coder Information
%#codegen
%$cgargs {zeros(4,1),zeros(4,1),zeros(4,1),zeros(3,1),zeros(6,1),zeros(21,1)}
assert(isreal(qJ) && all(size(qJ) == [4 1]), ...
  'S4RRPR10_invdynJ_fixb_mdp_slag_vp: qJ has to be [4x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [4 1]), ...
  'S4RRPR10_invdynJ_fixb_mdp_slag_vp: qJD has to be [4x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [4 1]), ...
  'S4RRPR10_invdynJ_fixb_mdp_slag_vp: qJDD has to be [4x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S4RRPR10_invdynJ_fixb_mdp_slag_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [6 1]), ...
  'S4RRPR10_invdynJ_fixb_mdp_slag_vp: pkin has to be [6x1] (double)');
assert(isreal(MDP) && all(size(MDP) == [21 1]), ...
  'S4RRPR10_invdynJ_fixb_mdp_slag_vp: MDP has to be [21x1] (double)'); 

%% Symbolic Calculation
% From invdyn_joint_fixb_mdp_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 17:12:00
% EndTime: 2019-12-31 17:12:02
% DurationCPUTime: 2.03s
% Computational Cost: add. (697->255), mult. (1511->335), div. (0->0), fcn. (864->6), ass. (0->129)
t298 = 2 * qJD(2);
t211 = cos(qJ(2));
t208 = sin(qJ(2));
t285 = qJ(3) * t208;
t236 = pkin(2) * t211 + t285;
t229 = pkin(1) + t236;
t169 = t229 * qJD(1);
t297 = t229 * qJDD(1);
t292 = pkin(3) + pkin(5);
t296 = t292 * t208;
t204 = g(3) * t211;
t266 = qJD(1) * t208;
t295 = -t169 * t266 + t204;
t252 = qJD(1) * qJD(2);
t244 = t211 * t252;
t251 = qJDD(1) * t208;
t224 = t244 + t251;
t175 = qJDD(4) + t224;
t210 = cos(qJ(4));
t167 = t210 * t175;
t207 = sin(qJ(4));
t195 = qJD(4) + t266;
t255 = qJD(4) * t195;
t294 = -t207 * t255 + t167;
t209 = sin(qJ(1));
t212 = cos(qJ(1));
t237 = g(1) * t209 - g(2) * t212;
t238 = g(1) * t212 + g(2) * t209;
t265 = qJD(1) * t211;
t199 = pkin(5) * t265;
t264 = qJD(2) * qJ(3);
t186 = -t199 - t264;
t200 = pkin(3) * t265;
t168 = -t186 + t200;
t213 = -pkin(2) - pkin(6);
t293 = t213 * t175 + (t168 - t199 - t200) * t195;
t287 = g(3) * t208;
t286 = pkin(5) * qJDD(2);
t284 = qJ(3) * t211;
t283 = qJDD(2) * pkin(2);
t261 = qJD(2) * t207;
t176 = t210 * t265 + t261;
t250 = qJDD(1) * t211;
t245 = t208 * t252;
t269 = t210 * qJDD(2) + t207 * t245;
t154 = -t176 * qJD(4) - t207 * t250 + t269;
t282 = t154 * t210;
t197 = pkin(5) * t250;
t249 = qJDD(2) * qJ(3);
t239 = -t197 - t249;
t158 = pkin(3) * t250 + (-qJD(1) * t296 + qJD(3)) * qJD(2) - t239;
t281 = t158 * t211;
t280 = t176 * t195;
t259 = qJD(2) * t210;
t178 = -t207 * t265 + t259;
t279 = t178 * t195;
t278 = t195 * t207;
t277 = t207 * t175;
t276 = t207 * t209;
t275 = t207 * t211;
t274 = t207 * t212;
t273 = t208 * t210;
t272 = t209 * t210;
t271 = t210 * t212;
t215 = qJD(1) ^ 2;
t270 = t211 * t215;
t254 = qJD(4) * t211;
t243 = qJD(1) * t254;
t268 = -t207 * t243 - t210 * t245;
t205 = t208 ^ 2;
t206 = t211 ^ 2;
t267 = t205 - t206;
t263 = qJD(2) * t176;
t262 = qJD(2) * t178;
t260 = qJD(2) * t208;
t258 = qJD(2) * t211;
t257 = qJD(3) * t208;
t174 = t213 * t211 - pkin(1) - t285;
t160 = t174 * qJD(1);
t256 = qJD(4) * t160;
t198 = pkin(5) * t266;
t253 = pkin(3) * t266 + qJD(3) + t198;
t196 = pkin(5) * t251;
t248 = qJDD(3) + t196;
t188 = t292 * t211;
t242 = pkin(5) * t244 + t248;
t241 = -qJD(2) * pkin(2) + qJD(3);
t194 = pkin(2) * t245;
t234 = pkin(6) * t208 - t284;
t219 = t234 * qJD(2) - t257;
t150 = t219 * qJD(1) + t174 * qJDD(1) + t194;
t163 = t213 * qJD(2) + t253;
t240 = qJD(4) * t163 + t150;
t235 = pkin(2) * t208 - t284;
t233 = -t256 + t204;
t152 = t160 * t210 + t163 * t207;
t184 = t198 + t241;
t231 = t184 * t211 + t186 * t208;
t228 = qJD(2) * qJD(4) + t250;
t164 = t242 - t283;
t227 = -0.2e1 * pkin(1) * t252 - t286;
t226 = -t210 * t255 - t277;
t225 = -qJ(3) * t258 - t257;
t223 = pkin(1) * t215 + t238;
t214 = qJD(2) ^ 2;
t222 = pkin(5) * t214 - t237;
t221 = t169 * t298 + t286;
t220 = 0.2e1 * qJDD(1) * pkin(1) - t222;
t156 = t225 * qJD(1) + t194 - t297;
t201 = pkin(2) * t260;
t166 = t201 + t225;
t218 = qJD(1) * t166 + t156 + t222 - t297;
t202 = pkin(2) * t266;
t217 = -t287 + t158 - t238 * t211 + (t234 * qJD(1) - qJD(4) * t213 + t202) * t195;
t161 = (-qJD(3) + t198) * qJD(2) + t239;
t216 = t231 * qJD(2) - t161 * t211 + t164 * t208 - t238;
t183 = qJD(2) * t188;
t181 = qJD(2) * t296;
t179 = -qJ(3) * t265 + t202;
t173 = -t208 * t276 + t271;
t172 = t208 * t272 + t274;
t171 = t208 * t274 + t272;
t170 = t208 * t271 - t276;
t159 = t201 + t219;
t157 = t224 * pkin(3) + t213 * qJDD(2) + t242;
t155 = qJDD(2) * t207 + t228 * t210 + t268;
t153 = t210 * t157;
t151 = -t160 * t207 + t163 * t210;
t1 = [qJDD(1) * MDP(1) + t237 * MDP(2) + t238 * MDP(3) + (qJDD(1) * t205 + 0.2e1 * t208 * t244) * MDP(4) + 0.2e1 * (t208 * t250 - t267 * t252) * MDP(5) + (qJDD(2) * t208 + t211 * t214) * MDP(6) + (qJDD(2) * t211 - t208 * t214) * MDP(7) + (t227 * t208 + t220 * t211) * MDP(9) + (-t220 * t208 + t227 * t211) * MDP(10) + ((t205 + t206) * qJDD(1) * pkin(5) + t216) * MDP(11) + (t221 * t208 + t218 * t211) * MDP(12) + (-t218 * t208 + t221 * t211) * MDP(13) + (t216 * pkin(5) - t169 * t166 + (-t156 + t237) * t229) * MDP(14) + (-t154 * t275 + (t207 * t260 - t210 * t254) * t178) * MDP(15) + ((-t176 * t207 + t178 * t210) * t260 + (-t282 + t155 * t207 + (t176 * t210 + t178 * t207) * qJD(4)) * t211) * MDP(16) + ((t195 * t261 + t154) * t208 + (t226 + t262) * t211) * MDP(17) + ((t195 * t259 - t155) * t208 + (-t263 - t294) * t211) * MDP(18) + (t175 * t208 + t195 * t258) * MDP(19) + ((-t159 * t207 + t183 * t210) * t195 + (-t174 * t207 + t210 * t296) * t175 + (-t150 * t207 + t153) * t208 - t181 * t176 + t188 * t155 + t210 * t281 - g(1) * t173 - g(2) * t171 + (t151 * t211 - t168 * t273) * qJD(2) + ((-t174 * t210 - t207 * t296) * t195 - t152 * t208 - t168 * t275) * qJD(4)) * MDP(20) + (-t152 * t258 + g(1) * t172 - g(2) * t170 + t188 * t154 - t181 * t178 + (-(qJD(4) * t296 + t159) * t195 - t174 * t175 - t240 * t208 - t168 * t254) * t210 + (-(-qJD(4) * t174 + t183) * t195 - t296 * t175 - t281 + (qJD(2) * t168 - t157 + t256) * t208) * t207) * MDP(21); -t208 * MDP(4) * t270 + t267 * MDP(5) * t215 + MDP(6) * t251 + MDP(7) * t250 + qJDD(2) * MDP(8) + (t223 * t208 - t196 - t204) * MDP(9) + (t223 * t211 - t197 + t287) * MDP(10) + (-t235 * qJDD(1) + ((-t186 - t264) * t208 + (-t184 + t241) * t211) * qJD(1)) * MDP(11) + (-t179 * t265 - t238 * t208 + t248 - 0.2e1 * t283 + t295) * MDP(12) + (0.2e1 * t249 + qJD(3) * t298 + t197 + (qJD(1) * t179 - g(3)) * t208 + (-qJD(1) * t169 - t238) * t211) * MDP(13) + (-t231 * qJD(1) * pkin(5) - t164 * pkin(2) - g(3) * t236 - t161 * qJ(3) - t186 * qJD(3) + t169 * t179 + t238 * t235) * MDP(14) + (-t178 * t278 + t282) * MDP(15) + ((-t155 - t279) * t210 + (-t154 + t280) * t207) * MDP(16) + ((-t178 * t211 - t208 * t278) * qJD(1) + t294) * MDP(17) + ((t176 * t211 - t195 * t273) * qJD(1) + t226) * MDP(18) - t195 * MDP(19) * t265 + (qJ(3) * t155 - t151 * t265 + t253 * t176 + t217 * t207 + t293 * t210) * MDP(20) + (qJ(3) * t154 + t152 * t265 + t253 * t178 - t293 * t207 + t217 * t210) * MDP(21); qJDD(2) * MDP(12) + (-t205 * t215 - t214) * MDP(13) + (qJD(2) * t186 + t164 + t295) * MDP(14) + (t167 - t263) * MDP(20) + (-t262 - t277) * MDP(21) + (qJDD(1) * MDP(11) + MDP(12) * t270 - t238 * MDP(14)) * t208 + (-t195 * MDP(21) * t210 - MDP(20) * t278) * t195; t178 * t176 * MDP(15) + (-t176 ^ 2 + t178 ^ 2) * MDP(16) + (t269 + t280) * MDP(17) + (-t268 + t279) * MDP(18) + t175 * MDP(19) + (-g(1) * t170 - g(2) * t172 + t152 * t195 - t168 * t178 + t153) * MDP(20) + (g(1) * t171 - g(2) * t173 + t151 * t195 + t168 * t176) * MDP(21) + (-MDP(17) * t243 - t228 * MDP(18) + t233 * MDP(20) - t240 * MDP(21)) * t210 + (-t228 * MDP(17) - qJDD(2) * MDP(18) - t240 * MDP(20) + (-t157 - t233) * MDP(21)) * t207;];
tau = t1;
