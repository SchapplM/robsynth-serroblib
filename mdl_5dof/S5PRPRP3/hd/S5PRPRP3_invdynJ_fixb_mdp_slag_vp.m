% Calculate vector of inverse dynamics joint torques for
% S5PRPRP3
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
%   pkin=[a2,a3,a4,a5,d2,d4,theta1,theta3]';
% MDP [16x1]
%   Minimal dynamic parameter vector (fixed base model)
%   see S5PRPRP3_convert_par2_MPV_fixb.m
% 
% Output:
% tau [5x1]
%   joint torques of inverse dynamics (contains inertial, gravitational Coriolis and centrifugal forces)

% Quelle: HybrDyn-Toolbox
% Datum: 2021-01-15 15:14
% Revision: 24b2e7d74a0c1a3b64fa2f8f5ad758691ad61af3 (2021-01-15)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tau = S5PRPRP3_invdynJ_fixb_mdp_slag_vp(qJ, qJD, qJDD, g, pkin, MDP)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(5,1),zeros(3,1),zeros(8,1),zeros(16,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5PRPRP3_invdynJ_fixb_mdp_slag_vp: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5PRPRP3_invdynJ_fixb_mdp_slag_vp: qJD has to be [5x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [5 1]), ...
  'S5PRPRP3_invdynJ_fixb_mdp_slag_vp: qJDD has to be [5x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S5PRPRP3_invdynJ_fixb_mdp_slag_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S5PRPRP3_invdynJ_fixb_mdp_slag_vp: pkin has to be [8x1] (double)');
assert(isreal(MDP) && all(size(MDP) == [16 1]), ...
  'S5PRPRP3_invdynJ_fixb_mdp_slag_vp: MDP has to be [16x1] (double)'); 

%% Symbolic Calculation
% From invdyn_joint_fixb_mdp_matlab.m
% OptimizationMode: 2
% StartTime: 2021-01-15 15:14:04
% EndTime: 2021-01-15 15:14:07
% DurationCPUTime: 1.28s
% Computational Cost: add. (828->217), mult. (1651->274), div. (0->0), fcn. (1131->10), ass. (0->113)
t212 = sin(qJ(2));
t214 = cos(qJ(2));
t249 = qJD(1) * qJD(2);
t281 = qJDD(1) * t212 + t214 * t249;
t213 = cos(qJ(4));
t255 = qJD(1) * t214;
t188 = qJD(2) * pkin(2) + t255;
t208 = cos(pkin(8));
t256 = qJD(1) * t212;
t192 = t208 * t256;
t206 = sin(pkin(8));
t162 = t206 * t188 + t192;
t155 = qJD(2) * pkin(6) + t162;
t257 = qJ(5) * qJD(2);
t237 = t155 + t257;
t230 = t237 * t213;
t207 = sin(pkin(7));
t209 = cos(pkin(7));
t240 = -g(1) * t207 + g(2) * t209;
t232 = g(1) * t209 + g(2) * t207;
t211 = sin(qJ(4));
t204 = t211 ^ 2;
t205 = t213 ^ 2;
t258 = t204 + t205;
t239 = MDP(15) * t258;
t203 = qJ(2) + pkin(8);
t198 = cos(t203);
t197 = sin(t203);
t273 = g(3) * t197;
t280 = -t232 * t198 - t273;
t279 = pkin(2) * t208;
t278 = pkin(4) * t204;
t277 = pkin(4) * t213;
t272 = g(3) * t198;
t271 = g(3) * t211;
t270 = qJD(4) * pkin(4);
t269 = qJDD(4) * pkin(4);
t268 = t207 * t211;
t267 = t207 * t213;
t266 = t209 * t211;
t265 = t209 * t213;
t216 = qJD(2) ^ 2;
t264 = t213 * t216;
t194 = pkin(2) * t206 + pkin(6);
t263 = qJ(5) + t194;
t262 = qJDD(1) - g(3);
t201 = t213 * qJD(3);
t148 = -t237 * t211 + t201;
t145 = t148 + t270;
t261 = -t148 + t145;
t191 = t206 * t256;
t171 = t208 * t255 - t191;
t250 = qJD(4) * t213;
t260 = t171 * t250 + t198 * t271;
t259 = t204 - t205;
t161 = t188 * t208 - t191;
t154 = -qJD(2) * pkin(3) - t161;
t254 = qJD(2) * t154;
t175 = t206 * t212 - t208 * t214;
t253 = qJD(2) * t175;
t252 = qJD(2) * t213;
t251 = qJD(4) * t211;
t248 = qJD(2) * qJD(4);
t246 = qJDD(2) * t211;
t245 = qJDD(2) * t213;
t244 = MDP(11) + MDP(13);
t243 = MDP(12) + MDP(14);
t200 = t214 * qJDD(1);
t172 = qJDD(2) * pkin(2) - t212 * t249 + t200;
t151 = t206 * t172 + t281 * t208;
t196 = pkin(3) + t277;
t241 = t211 * t248;
t238 = qJD(4) * t263;
t236 = 0.2e1 * t213 * t248;
t147 = qJDD(2) * pkin(6) + t151;
t235 = -qJD(4) * qJD(3) - t147;
t150 = t172 * t208 - t281 * t206;
t144 = pkin(4) * t241 - t196 * qJDD(2) + qJDD(5) - t150;
t180 = -t196 - t279;
t234 = qJDD(2) * t180 + t144;
t169 = t206 * t255 + t192;
t233 = t169 * t252 + t171 * t251 + (g(1) * t265 + g(2) * t267) * t197;
t149 = qJD(3) * t211 + t230;
t231 = t145 * t211 - t149 * t213;
t176 = t206 * t214 + t208 * t212;
t229 = t232 * t197;
t195 = -pkin(3) - t279;
t215 = qJD(4) ^ 2;
t228 = t194 * t215 - t150 + (-pkin(3) + t195) * qJDD(2);
t199 = t213 * qJDD(3);
t227 = -g(1) * (-t198 * t266 + t267) - g(2) * (-t198 * t268 - t265) + t197 * t271 + t199;
t226 = -qJ(5) * qJDD(2) + t235;
t168 = t176 * qJD(2);
t224 = qJD(2) * t168 + qJDD(2) * t175 + t176 * t215;
t223 = 0.2e1 * t253 * qJD(4) - qJDD(4) * t176;
t222 = -qJDD(4) * t194 + (qJD(2) * t195 + t154) * qJD(4);
t221 = -g(3) * t214 + t232 * t212;
t220 = -qJD(2) * t169 - t229;
t153 = t155 * t251;
t219 = -g(1) * (-t198 * t265 - t268) - g(2) * (-t198 * t267 + t266) - t211 * qJDD(3) + t153 + t213 * t273;
t218 = qJD(2) * qJD(5) - t226;
t152 = -t196 * qJD(2) + qJD(5) - t161;
t217 = (-qJD(5) - t152) * qJD(2) + t226;
t210 = -qJ(5) - pkin(6);
t182 = qJDD(4) * t213 - t211 * t215;
t181 = qJDD(4) * t211 + t213 * t215;
t174 = t263 * t213;
t173 = t263 * t211;
t157 = -qJD(5) * t211 - t213 * t238;
t156 = qJD(5) * t213 - t211 * t238;
t143 = -t153 + (-qJ(5) * t248 + qJDD(3)) * t211 + t218 * t213;
t142 = -qJD(4) * t230 - t218 * t211 + t199 + t269;
t1 = [t262 * MDP(1) + (qJDD(2) * t214 - t212 * t216) * MDP(3) + (-qJDD(2) * t212 - t214 * t216) * MDP(4) + (-t150 * t175 - t161 * t168 - g(3)) * MDP(5) + (t144 * t175 + t152 * t168 - g(3)) * MDP(16) + t244 * (t223 * t211 - t224 * t213) + t243 * (t224 * t211 + t223 * t213) - (-t231 * MDP(16) + t162 * MDP(5) + qJD(2) * t239) * t253 + (t151 * MDP(5) + (-t142 * t211 + t143 * t213 - t145 * t250 - t149 * t251) * MDP(16) + qJDD(2) * t239) * t176; qJDD(2) * MDP(2) + (t200 + t221) * MDP(3) + (-t262 * t212 + t232 * t214) * MDP(4) + (t161 * t169 - t162 * t171 + (t150 * t208 + t151 * t206 + t221) * pkin(2)) * MDP(5) + (qJDD(2) * t204 + t211 * t236) * MDP(6) + 0.2e1 * (t211 * t245 - t259 * t248) * MDP(7) + t181 * MDP(8) + t182 * MDP(9) + (t222 * t211 + (-t228 - t272) * t213 + t233) * MDP(11) + (t222 * t213 + (t220 + t228) * t211 + t260) * MDP(12) + (-qJDD(4) * t173 + (-t234 - t272) * t213 + (t157 + (t152 + (t180 - t277) * qJD(2)) * t211) * qJD(4) + t233) * MDP(13) + (-qJDD(4) * t174 + (t152 * t213 - t156 + (t180 * t213 + t278) * qJD(2)) * qJD(4) + (t220 + t234) * t211 + t260) * MDP(14) + ((-qJD(4) * t145 + qJDD(2) * t174 + t143) * t213 + (-qJD(4) * t149 + qJDD(2) * t173 - t142) * t211 + (t156 * t213 - t157 * t211 - t258 * t171 + (t173 * t213 - t174 * t211) * qJD(4)) * qJD(2) + t280) * MDP(15) + (t143 * t174 - t142 * t173 + t144 * t180 - g(3) * (pkin(2) * t214 + t196 * t198 - t197 * t210) + (pkin(4) * t251 - t169) * t152 + (-t171 * t213 + t156) * t149 + (t171 * t211 + t157) * t145 + t232 * (pkin(2) * t212 + t196 * t197 + t198 * t210)) * MDP(16); (qJDD(3) + t240) * MDP(5) + (-t231 * qJD(4) + t142 * t213 + t143 * t211 + t240) * MDP(16) + t244 * t182 - t243 * t181; -t211 * MDP(6) * t264 + t259 * t216 * MDP(7) + MDP(8) * t246 + MDP(9) * t245 + qJDD(4) * MDP(10) + ((-t147 - t254) * t211 + t227) * MDP(11) + ((-t155 * t211 + t201) * qJD(4) + (t235 - t254) * t213 + t219) * MDP(12) + (0.2e1 * t269 + (t149 - t230) * qJD(4) + (pkin(4) * t264 + t217) * t211 + t227) * MDP(13) + (-t216 * t278 + (t211 * t257 + t148) * qJD(4) + t217 * t213 + t219) * MDP(14) + (-pkin(4) * t246 + (t261 - t270) * t252) * MDP(15) + (t261 * t149 + (t142 + t240 * t213 + (-t152 * qJD(2) - t280) * t211) * pkin(4)) * MDP(16); (0.2e1 * t241 - t245) * MDP(13) + (t236 + t246) * MDP(14) + (t231 * qJD(2) + t144 - t229 + t272) * MDP(16) - t216 * t239;];
tau = t1;
