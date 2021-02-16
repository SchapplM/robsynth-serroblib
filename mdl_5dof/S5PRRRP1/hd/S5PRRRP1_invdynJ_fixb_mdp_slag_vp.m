% Calculate vector of inverse dynamics joint torques for
% S5PRRRP1
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
% pkin [8x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d2,d3,d4,theta1]';
% MDP [18x1]
%   Minimal dynamic parameter vector (fixed base model)
%   see S5PRRRP1_convert_par2_MPV_fixb.m
% 
% Output:
% tau [5x1]
%   joint torques of inverse dynamics (contains inertial, gravitational Coriolis and centrifugal forces)

% Quelle: HybrDyn-Toolbox
% Datum: 2021-01-15 16:15
% Revision: 24b2e7d74a0c1a3b64fa2f8f5ad758691ad61af3 (2021-01-15)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tau = S5PRRRP1_invdynJ_fixb_mdp_slag_vp(qJ, qJD, qJDD, g, pkin, MDP)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(5,1),zeros(3,1),zeros(8,1),zeros(18,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5PRRRP1_invdynJ_fixb_mdp_slag_vp: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5PRRRP1_invdynJ_fixb_mdp_slag_vp: qJD has to be [5x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [5 1]), ...
  'S5PRRRP1_invdynJ_fixb_mdp_slag_vp: qJDD has to be [5x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S5PRRRP1_invdynJ_fixb_mdp_slag_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S5PRRRP1_invdynJ_fixb_mdp_slag_vp: pkin has to be [8x1] (double)');
assert(isreal(MDP) && all(size(MDP) == [18 1]), ...
  'S5PRRRP1_invdynJ_fixb_mdp_slag_vp: MDP has to be [18x1] (double)'); 

%% Symbolic Calculation
% From invdyn_joint_fixb_mdp_matlab.m
% OptimizationMode: 2
% StartTime: 2021-01-15 16:14:44
% EndTime: 2021-01-15 16:14:48
% DurationCPUTime: 1.13s
% Computational Cost: add. (994->235), mult. (1380->279), div. (0->0), fcn. (695->8), ass. (0->128)
t291 = qJDD(1) - g(3);
t218 = cos(qJ(4));
t212 = qJD(2) + qJD(3);
t217 = sin(qJ(3));
t283 = pkin(2) * qJD(2);
t257 = t217 * t283;
t173 = pkin(7) * t212 + t257;
t244 = qJ(5) * t212 + t173;
t235 = t244 * t218;
t211 = pkin(8) + qJ(2);
t204 = qJ(3) + t211;
t194 = sin(t204);
t190 = g(2) * t194;
t195 = cos(t204);
t192 = g(1) * t195;
t269 = t190 + t192;
t216 = sin(qJ(4));
t261 = qJD(4) * t216;
t248 = t212 * t261;
t287 = pkin(4) * t218;
t197 = pkin(3) + t287;
t210 = qJDD(2) + qJDD(3);
t276 = t197 * t210;
t290 = -pkin(4) * t248 + t276;
t219 = cos(qJ(3));
t289 = pkin(2) * t219;
t213 = t216 ^ 2;
t288 = pkin(4) * t213;
t191 = g(1) * t194;
t286 = g(2) * t195;
t285 = g(3) * t218;
t284 = t210 * pkin(3);
t215 = qJ(5) + pkin(7);
t282 = qJD(4) * pkin(4);
t281 = qJDD(4) * pkin(4);
t265 = qJD(2) * t219;
t254 = pkin(2) * t265;
t174 = -pkin(3) * t212 - t254;
t280 = t174 * t212;
t183 = -t197 - t289;
t279 = t183 * t212;
t278 = t194 * t218;
t277 = t195 * t216;
t275 = t210 * t216;
t274 = t210 * t218;
t273 = t212 * t216;
t272 = t216 * t218;
t196 = pkin(2) * t217 + pkin(7);
t271 = -qJ(5) - t196;
t206 = t218 * qJD(1);
t152 = -t244 * t216 + t206;
t148 = t152 + t282;
t270 = t148 - t152;
t264 = qJD(3) * t217;
t256 = pkin(2) * t264;
t268 = qJD(2) * t256 - qJDD(2) * t289;
t214 = t218 ^ 2;
t267 = -t213 - t214;
t266 = t213 - t214;
t263 = qJD(3) * t219;
t262 = qJD(4) * t212;
t260 = qJD(4) * t218;
t259 = qJDD(2) * t217;
t258 = qJDD(4) * t196;
t255 = pkin(2) * t263;
t165 = t173 * t261;
t227 = qJD(2) * t263 + t259;
t160 = t227 * pkin(2) + t210 * pkin(7);
t240 = -qJD(4) * qJD(1) - t160;
t228 = -qJ(5) * t210 + t240;
t224 = qJD(5) * t212 - t228;
t145 = -t165 + (-qJ(5) * t262 + qJDD(1)) * t216 + t224 * t218;
t253 = t145 * t218 - t269;
t147 = qJDD(5) + t268 - t290;
t158 = -t197 * t212 + qJD(5) - t254;
t179 = g(2) * t277;
t252 = t147 * t216 + t158 * t260 + t179;
t159 = t268 - t284;
t251 = t159 * t216 + t174 * t260 + t179;
t180 = g(1) * t278;
t239 = qJD(4) * t254;
t241 = t212 * t257;
t250 = t216 * t239 + t218 * t241 + t180;
t249 = t212 * t264;
t154 = t158 * t261;
t247 = -t147 - t286;
t246 = -t159 - t286;
t245 = qJD(4) * t215;
t243 = qJD(4) * t271;
t242 = 0.2e1 * t212 * t260;
t220 = qJD(4) ^ 2;
t238 = -pkin(7) * t220 + t284;
t201 = sin(t211);
t202 = cos(t211);
t237 = g(1) * t201 - g(2) * t202;
t175 = qJDD(4) * t216 + t218 * t220;
t176 = qJDD(4) * t218 - t216 * t220;
t236 = 0.2e1 * (t210 * t272 - t266 * t262) * MDP(9) + (t210 * t213 + t216 * t242) * MDP(8) + t175 * MDP(10) + t176 * MDP(11) + t210 * MDP(5);
t153 = t216 * qJD(1) + t235;
t234 = t148 * t216 - t153 * t218;
t171 = pkin(4) * t261 + t256;
t233 = t171 * t212 + t183 * t210;
t232 = -t191 + t268 + t286;
t203 = t218 * qJDD(1);
t231 = g(1) * t277 + t216 * t190 + t203 - t285;
t198 = -pkin(3) - t289;
t230 = t198 * t212 - t255;
t229 = -pkin(3) * t262 - pkin(7) * qJDD(4);
t226 = -t241 - t191;
t225 = g(2) * t278 + t218 * t192 - t291 * t216 + t165;
t223 = pkin(2) * t249 + t196 * t220 + t198 * t210;
t222 = (-qJD(5) - t158) * t212 + t228;
t221 = -g(1) * (-t194 * t197 + t215 * t195) - g(2) * (t194 * t215 + t195 * t197);
t209 = t212 ^ 2;
t207 = t218 * qJ(5);
t205 = t218 * qJD(5);
t188 = pkin(7) * t218 + t207;
t187 = t215 * t216;
t185 = t218 * t239;
t167 = t196 * t218 + t207;
t166 = t271 * t216;
t163 = t174 * t261;
t162 = -t216 * qJD(5) - t218 * t245;
t161 = -t216 * t245 + t205;
t151 = (-qJD(5) - t255) * t216 + t218 * t243;
t150 = t216 * t243 + t218 * t255 + t205;
t144 = -qJD(4) * t235 - t224 * t216 + t203 + t281;
t1 = [t291 * MDP(1) + (-t234 * qJD(4) + t144 * t218 + t145 * t216 - g(3)) * MDP(18) + (MDP(13) + MDP(15)) * t176 + (-MDP(14) - MDP(16)) * t175; -t232 * MDP(6) + ((-t223 + t246) * MDP(13) - MDP(14) * t258 + (-t233 + t247) * MDP(15) + (t150 * t212 + t167 * t210) * MDP(17) + (t230 * MDP(14) + MDP(16) * t279 + (-t166 * t212 - t148) * MDP(17)) * qJD(4)) * t218 + (-MDP(13) * t258 + (t223 - t191) * MDP(14) + (t233 - t191) * MDP(16) + (-t151 * t212 - t166 * t210 - t144) * MDP(17) + (t230 * MDP(13) + MDP(15) * t279 + (-t167 * t212 - t153) * MDP(17)) * qJD(4)) * t216 + ((t210 * t219 - t249) * MDP(6) + (-t210 * t217 - t212 * t263 - t227) * MDP(7) + t237 * MDP(18)) * pkin(2) + (t144 * t166 + t145 * t167 + t147 * t183 + t148 * t151 + t153 * t150 + t158 * t171 + t221) * MDP(18) + t251 * MDP(14) + t253 * MDP(17) + t269 * MDP(7) + t236 + (t151 * qJD(4) + t166 * qJDD(4) + t154 + t180) * MDP(15) + t237 * MDP(3) + (g(1) * t202 + g(2) * t201) * MDP(4) + (t163 + t180) * MDP(13) + (-t150 * qJD(4) - t167 * qJDD(4) + t252) * MDP(16) + qJDD(2) * MDP(2); (-t232 + t241) * MDP(6) + ((-t259 + (-qJD(3) + t212) * t265) * pkin(2) + t269) * MDP(7) + (t163 + t229 * t216 + (t238 + t246) * t218 + t250) * MDP(13) + (t185 + t229 * t218 + (t226 - t238) * t216 + t251) * MDP(14) + (-qJDD(4) * t187 + t154 + (-t197 * t273 + t162) * qJD(4) + (t247 + t290) * t218 + t250) * MDP(15) + (-qJDD(4) * t188 + t185 + (t226 - t276) * t216 + (-t161 + (-t197 * t218 + t288) * t212) * qJD(4) + t252) * MDP(16) + ((-qJD(4) * t148 + t188 * t210) * t218 + (-qJD(4) * t153 + t187 * t210 - t144) * t216 + (t161 * t218 - t162 * t216 + (t187 * t218 - t188 * t216) * qJD(4) + t267 * t254) * t212 + t253) * MDP(17) + (t145 * t188 + t153 * t161 - t144 * t187 + t148 * t162 - t147 * t197 + pkin(4) * t154 + (-t158 * t217 + t234 * t219) * t283 + t221) * MDP(18) + t236; -t209 * MDP(8) * t272 + t266 * t209 * MDP(9) + MDP(10) * t275 + MDP(11) * t274 + qJDD(4) * MDP(12) + ((-t160 - t280) * t216 + t231) * MDP(13) + ((-t173 * t216 + t206) * qJD(4) + (t240 - t280) * t218 + t225) * MDP(14) + (0.2e1 * t281 + (t153 - t235) * qJD(4) + (t209 * t287 + t222) * t216 + t231) * MDP(15) + (-t209 * t288 + (qJ(5) * t273 + t152) * qJD(4) + t222 * t218 + t225) * MDP(16) + (-pkin(4) * t275 + (t270 - t282) * t218 * t212) * MDP(17) + (t270 * t153 + (-t285 + t144 + (-t158 * t212 + t269) * t216) * pkin(4)) * MDP(18); (0.2e1 * t248 - t274) * MDP(15) + (t242 + t275) * MDP(16) + (t234 * t212 - t191 - t247) * MDP(18) + t267 * MDP(17) * t209;];
tau = t1;
