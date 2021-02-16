% Calculate vector of inverse dynamics joint torques for
% S5PRRPP1
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
%   pkin=[a2,a3,a4,a5,d2,d3,theta1,theta4]';
% MDP [19x1]
%   Minimal dynamic parameter vector (fixed base model)
%   see S5PRRPP1_convert_par2_MPV_fixb.m
% 
% Output:
% tau [5x1]
%   joint torques of inverse dynamics (contains inertial, gravitational Coriolis and centrifugal forces)

% Quelle: HybrDyn-Toolbox
% Datum: 2021-01-15 15:23
% Revision: 24b2e7d74a0c1a3b64fa2f8f5ad758691ad61af3 (2021-01-15)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tau = S5PRRPP1_invdynJ_fixb_mdp_slag_vp(qJ, qJD, qJDD, g, pkin, MDP)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(5,1),zeros(3,1),zeros(8,1),zeros(19,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5PRRPP1_invdynJ_fixb_mdp_slag_vp: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5PRRPP1_invdynJ_fixb_mdp_slag_vp: qJD has to be [5x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [5 1]), ...
  'S5PRRPP1_invdynJ_fixb_mdp_slag_vp: qJDD has to be [5x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S5PRRPP1_invdynJ_fixb_mdp_slag_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S5PRRPP1_invdynJ_fixb_mdp_slag_vp: pkin has to be [8x1] (double)');
assert(isreal(MDP) && all(size(MDP) == [19 1]), ...
  'S5PRRPP1_invdynJ_fixb_mdp_slag_vp: MDP has to be [19x1] (double)'); 

%% Symbolic Calculation
% From invdyn_joint_fixb_mdp_matlab.m
% OptimizationMode: 2
% StartTime: 2021-01-15 15:22:44
% EndTime: 2021-01-15 15:22:49
% DurationCPUTime: 1.89s
% Computational Cost: add. (1325->260), mult. (2734->314), div. (0->0), fcn. (1753->8), ass. (0->106)
t220 = pkin(7) + qJ(2);
t214 = sin(t220);
t216 = cos(t220);
t286 = g(1) * t214 - g(2) * t216;
t249 = g(1) * t216 + g(2) * t214;
t229 = cos(qJ(3));
t281 = pkin(3) * t229;
t213 = pkin(2) + t281;
t247 = t213 * qJDD(2);
t288 = MDP(14) + MDP(17);
t225 = sin(pkin(8));
t226 = cos(pkin(8));
t228 = sin(qJ(3));
t194 = t225 * t229 + t226 * t228;
t190 = t194 * qJD(2);
t261 = MDP(12) + MDP(16);
t276 = qJ(4) + pkin(6);
t203 = t276 * t229;
t256 = t276 * t228;
t175 = t226 * t203 - t225 * t256;
t221 = qJ(3) + pkin(8);
t215 = sin(t221);
t285 = -qJDD(3) * t175 - t215 * t286;
t218 = t229 * qJDD(1);
t252 = t276 * qJD(3);
t250 = qJD(2) * t252;
t283 = qJD(1) * qJD(3) + qJD(2) * qJD(4) + qJDD(2) * t276;
t159 = qJDD(3) * pkin(3) - t228 * t283 - t229 * t250 + t218;
t163 = (qJDD(1) - t250) * t228 + t283 * t229;
t147 = t159 * t225 + t163 * t226;
t217 = cos(t221);
t284 = g(3) * t215 + t249 * t217 - t147;
t186 = t190 ^ 2;
t282 = pkin(3) * t228;
t277 = g(3) * t229;
t275 = qJDD(3) * pkin(4);
t184 = qJD(1) * t228 + qJD(2) * t203;
t274 = t184 * t225;
t273 = t225 * t228;
t177 = t226 * t184;
t272 = t226 * t229;
t271 = qJDD(1) - g(3);
t146 = t159 * t226 - t163 * t225;
t182 = qJD(1) * t229 - qJD(2) * t256;
t180 = qJD(3) * pkin(3) + t182;
t162 = t180 * t225 + t177;
t223 = t228 ^ 2;
t270 = -t229 ^ 2 + t223;
t269 = qJD(2) * t228;
t165 = t182 * t225 + t177;
t268 = qJD(3) * t165;
t267 = qJD(3) * t228;
t167 = t182 * t226 - t274;
t266 = qJD(5) - t167;
t265 = qJD(2) * qJD(3);
t264 = qJDD(2) * t228;
t263 = qJDD(2) * t229;
t260 = MDP(13) - MDP(18);
t255 = t228 * t265;
t259 = pkin(3) * t255 + qJDD(4);
t258 = pkin(3) * t267;
t257 = qJD(2) * t272;
t254 = t229 * t265;
t246 = pkin(4) * t217 + qJ(5) * t215;
t161 = t180 * t226 - t274;
t243 = -0.2e1 * pkin(2) * t265 - pkin(6) * qJDD(3);
t242 = -g(3) * t217 + t215 * t249 + t146;
t200 = -qJD(2) * t213 + qJD(4);
t189 = t194 * qJD(3);
t174 = t203 * t225 + t226 * t256;
t240 = -qJDD(3) * t174 + t217 * t286;
t238 = -qJD(4) * t228 - t229 * t252;
t237 = qJDD(2) * t194 - t225 * t255;
t230 = qJD(3) ^ 2;
t236 = 0.2e1 * qJDD(2) * pkin(2) - pkin(6) * t230 + t286;
t231 = qJD(2) ^ 2;
t235 = pkin(2) * t231 - pkin(6) * qJDD(2) + t249;
t205 = t226 * t263;
t170 = qJD(2) * t189 + t225 * t264 - t205;
t171 = t226 * t254 + t237;
t234 = pkin(4) * t170 - qJ(5) * t171 - qJD(5) * t190 + t259;
t187 = t225 * t269 - t257;
t160 = pkin(4) * t187 - qJ(5) * t190 + t200;
t233 = -t160 * t190 - qJDD(5) + t242;
t183 = qJD(4) * t229 - t228 * t252;
t166 = t183 * t225 - t226 * t238;
t168 = t226 * t183 + t225 * t238;
t232 = t166 * t190 - t168 * t187 - t170 * t175 + t171 * t174 - t249;
t222 = qJDD(3) * qJ(5);
t212 = -pkin(3) * t226 - pkin(4);
t209 = pkin(3) * t225 + qJ(5);
t202 = qJDD(3) * t229 - t228 * t230;
t201 = qJDD(3) * t228 + t229 * t230;
t199 = t216 * t213;
t193 = -t272 + t273;
t192 = qJD(3) * t272 - t225 * t267;
t181 = -t247 + t259;
t169 = pkin(4) * t193 - qJ(5) * t194 - t213;
t164 = pkin(3) * t269 + pkin(4) * t190 + qJ(5) * t187;
t156 = qJD(3) * qJ(5) + t162;
t154 = -qJD(3) * pkin(4) + qJD(5) - t161;
t150 = pkin(4) * t189 - qJ(5) * t192 - qJD(5) * t194 + t258;
t145 = qJDD(5) - t146 - t275;
t144 = qJD(3) * qJD(5) + t147 + t222;
t143 = -t247 + t234;
t1 = [t271 * MDP(1) + t202 * MDP(10) - t201 * MDP(11) + (-t146 * t193 + t147 * t194 - t161 * t189 + t162 * t192 - g(3)) * MDP(15) + (t144 * t194 + t145 * t193 + t154 * t189 + t156 * t192 - g(3)) * MDP(19) - t260 * (qJD(3) * t192 + qJDD(3) * t194) + t261 * (-qJD(3) * t189 - qJDD(3) * t193) + t288 * (-t170 * t194 + t171 * t193 - t187 * t192 + t189 * t190); qJDD(2) * MDP(2) + t286 * MDP(3) + t249 * MDP(4) + (qJDD(2) * t223 + 0.2e1 * t228 * t254) * MDP(5) + 0.2e1 * (t228 * t263 - t265 * t270) * MDP(6) + t201 * MDP(7) + t202 * MDP(8) + (t228 * t243 + t229 * t236) * MDP(10) + (-t228 * t236 + t229 * t243) * MDP(11) + (-t170 * t213 + t181 * t193 + t189 * t200 + (t187 * t282 - t166) * qJD(3) + t240) * MDP(12) + (-t171 * t213 + t181 * t194 + t192 * t200 + (t190 * t282 - t168) * qJD(3) + t285) * MDP(13) + (-t146 * t194 - t147 * t193 - t161 * t192 - t162 * t189 + t232) * MDP(14) + (t147 * t175 + t162 * t168 - t146 * t174 - t161 * t166 - t181 * t213 + t200 * t258 - g(1) * (-t213 * t214 + t216 * t276) - g(2) * (t214 * t276 + t199)) * MDP(15) + (-qJD(3) * t166 + t143 * t193 + t150 * t187 + t160 * t189 + t169 * t170 + t240) * MDP(16) + (-t144 * t193 + t145 * t194 + t154 * t192 - t156 * t189 + t232) * MDP(17) + (qJD(3) * t168 - t143 * t194 - t150 * t190 - t160 * t192 - t169 * t171 - t285) * MDP(18) + (-g(2) * t199 + t143 * t169 + t144 * t175 + t145 * t174 + t160 * t150 + t154 * t166 + t156 * t168 + (-g(1) * t276 - g(2) * t246) * t216 + (-g(1) * (-t213 - t246) - g(2) * t276) * t214) * MDP(19); MDP(7) * t264 + MDP(8) * t263 + qJDD(3) * MDP(9) + (t228 * t235 + t218 - t277) * MDP(10) + (-t228 * t271 + t229 * t235) * MDP(11) + (t268 - t190 * t200 + (qJDD(3) * t226 - t187 * t269) * pkin(3) + t242) * MDP(12) + (qJD(3) * t167 + t187 * t200 + (-qJDD(3) * t225 - t190 * t269) * pkin(3) + t284) * MDP(13) + ((t162 - t165) * t190 + (-t161 + t167) * t187 + (-t170 * t225 - t171 * t226) * pkin(3)) * MDP(14) + (t161 * t165 - t162 * t167 + (-t277 + t146 * t226 + t147 * t225 + (-qJD(2) * t200 + t249) * t228) * pkin(3)) * MDP(15) + (t268 - t164 * t187 + (pkin(4) - t212) * qJDD(3) + t233) * MDP(16) + (-t170 * t209 + t171 * t212 + (t156 - t165) * t190 + (t154 - t266) * t187) * MDP(17) + (qJDD(3) * t209 - t160 * t187 + t164 * t190 + t222 + (0.2e1 * qJD(5) - t167) * qJD(3) - t284) * MDP(18) + (t144 * t209 + t145 * t212 - t160 * t164 - t154 * t165 - g(3) * (t246 + t281) + t266 * t156 + t249 * (pkin(4) * t215 - qJ(5) * t217 + t282)) * MDP(19) + (-MDP(5) * t228 * t229 + MDP(6) * t270) * t231; (t161 * t190 + t162 * t187 + t259 - t286) * MDP(15) + (-t154 * t190 + t156 * t187 + t234 - t286) * MDP(19) + t260 * ((-t187 + t257) * qJD(3) + t237) + (-MDP(15) - MDP(19)) * t247 + t288 * (-t187 ^ 2 - t186) + (0.2e1 * qJD(3) * t190 + qJDD(2) * t273 - t205) * t261; (t187 * t190 - qJDD(3)) * MDP(16) + ((t187 + t257) * qJD(3) + t237) * MDP(17) + (-t186 - t230) * MDP(18) + (-qJD(3) * t156 - t233 - t275) * MDP(19);];
tau = t1;
