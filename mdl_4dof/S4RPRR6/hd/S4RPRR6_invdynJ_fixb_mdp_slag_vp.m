% Calculate vector of inverse dynamics joint torques for
% S4RPRR6
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
% pkin [7x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,d1,d3,d4,theta2]';
% MDP [21x1]
%   Minimal dynamic parameter vector (fixed base model)
%   see S4RPRR6_convert_par2_MPV_fixb.m
% 
% Output:
% tau [4x1]
%   joint torques of inverse dynamics (contains inertial, gravitational Coriolis and centrifugal forces)

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-31 16:52
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tau = S4RPRR6_invdynJ_fixb_mdp_slag_vp(qJ, qJD, qJDD, g, pkin, MDP)
%% Coder Information
%#codegen
%$cgargs {zeros(4,1),zeros(4,1),zeros(4,1),zeros(3,1),zeros(7,1),zeros(21,1)}
assert(isreal(qJ) && all(size(qJ) == [4 1]), ...
  'S4RPRR6_invdynJ_fixb_mdp_slag_vp: qJ has to be [4x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [4 1]), ...
  'S4RPRR6_invdynJ_fixb_mdp_slag_vp: qJD has to be [4x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [4 1]), ...
  'S4RPRR6_invdynJ_fixb_mdp_slag_vp: qJDD has to be [4x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S4RPRR6_invdynJ_fixb_mdp_slag_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [7 1]), ...
  'S4RPRR6_invdynJ_fixb_mdp_slag_vp: pkin has to be [7x1] (double)');
assert(isreal(MDP) && all(size(MDP) == [21 1]), ...
  'S4RPRR6_invdynJ_fixb_mdp_slag_vp: MDP has to be [21x1] (double)'); 

%% Symbolic Calculation
% From invdyn_joint_fixb_mdp_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 16:52:45
% EndTime: 2019-12-31 16:52:47
% DurationCPUTime: 1.37s
% Computational Cost: add. (943->201), mult. (2238->268), div. (0->0), fcn. (1732->12), ass. (0->101)
t221 = sin(pkin(7));
t224 = sin(qJ(3));
t222 = cos(pkin(7));
t227 = cos(qJ(3));
t264 = t227 * t222;
t189 = t221 * t224 - t264;
t181 = t189 * qJD(1);
t226 = cos(qJ(4));
t266 = t221 * t227;
t190 = t222 * t224 + t266;
t182 = t190 * qJD(1);
t223 = sin(qJ(4));
t268 = t182 * t223;
t159 = t226 * t181 + t268;
t220 = qJD(3) + qJD(4);
t269 = t159 * t220;
t225 = sin(qJ(1));
t228 = cos(qJ(1));
t245 = g(1) * t228 + g(2) * t225;
t274 = pkin(5) + qJ(2);
t199 = t274 * t221;
t191 = qJD(1) * t199;
t200 = t274 * t222;
t192 = qJD(1) * t200;
t237 = t191 * t224 - t192 * t227;
t154 = -pkin(6) * t181 - t237;
t210 = -pkin(2) * t222 - pkin(1);
t194 = t210 * qJD(1) + qJD(2);
t170 = pkin(3) * t181 + t194;
t219 = pkin(7) + qJ(3);
t215 = qJ(4) + t219;
t208 = sin(t215);
t209 = cos(t215);
t258 = qJD(4) * t223;
t285 = g(3) * t208 + t154 * t258 + t170 * t159 + t245 * t209;
t216 = qJDD(3) + qJDD(4);
t238 = -t181 * t223 + t226 * t182;
t284 = t216 * MDP(19) + t159 * MDP(15) * t238 + (-t159 ^ 2 + t238 ^ 2) * MDP(16);
t270 = t238 * t220;
t282 = -t227 * t191 - t192 * t224;
t263 = -t224 * t199 + t227 * t200;
t281 = qJ(2) * qJDD(1);
t244 = g(1) * t225 - g(2) * t228;
t280 = qJDD(2) - t244;
t260 = qJD(1) * t224;
t251 = t221 * t260;
t254 = qJDD(1) * t227;
t255 = qJDD(1) * t224;
t259 = qJD(3) * t227;
t252 = t221 * t254 + (qJD(1) * t259 + t255) * t222;
t166 = -qJD(3) * t251 + t252;
t256 = qJD(1) * qJD(2);
t278 = t274 * qJDD(1) + t256;
t174 = t278 * t221;
t175 = t278 * t222;
t249 = -t227 * t174 - t224 * t175;
t143 = qJDD(3) * pkin(3) - pkin(6) * t166 + t237 * qJD(3) + t249;
t184 = t190 * qJD(3);
t204 = t222 * t254;
t243 = -t221 * t255 + t204;
t167 = qJD(1) * t184 - t243;
t239 = -t224 * t174 + t227 * t175;
t144 = -pkin(6) * t167 + qJD(3) * t282 + t239;
t279 = -g(3) * t209 + t226 * t143 - t223 * t144 - t170 * t238 + t245 * t208;
t250 = t166 * t223 + t226 * t167;
t146 = t238 * qJD(4) + t250;
t277 = pkin(3) * t184;
t273 = qJDD(1) * pkin(1);
t153 = -pkin(6) * t182 + t282;
t152 = qJD(3) * pkin(3) + t153;
t272 = t152 * t226;
t271 = t154 * t226;
t265 = t222 * MDP(4);
t262 = t221 ^ 2 + t222 ^ 2;
t257 = qJD(4) * t226;
t253 = t226 * t166 - t223 * t167 - t181 * t257;
t248 = -t227 * t199 - t200 * t224;
t246 = 0.2e1 * t262;
t242 = -t152 * t223 - t271;
t156 = -pkin(6) * t190 + t248;
t157 = -pkin(6) * t189 + t263;
t241 = t156 * t226 - t157 * t223;
t240 = t156 * t223 + t157 * t226;
t168 = t226 * t189 + t190 * t223;
t169 = -t189 * t223 + t190 * t226;
t236 = t273 - t280;
t193 = t210 * qJDD(1) + qJDD(2);
t145 = -t182 * t258 + t253;
t233 = -t199 * t259 + qJD(2) * t264 + (-qJD(2) * t221 - qJD(3) * t200) * t224;
t232 = t246 * t256 - t245;
t231 = -t190 * qJD(2) - qJD(3) * t263;
t214 = cos(t219);
t213 = sin(t219);
t183 = t189 * qJD(3);
t172 = pkin(3) * t189 + t210;
t155 = pkin(3) * t167 + t193;
t150 = pkin(6) * t183 + t231;
t149 = -pkin(6) * t184 + t233;
t148 = t169 * qJD(4) - t183 * t223 + t226 * t184;
t147 = -t168 * qJD(4) - t183 * t226 - t184 * t223;
t1 = [qJDD(1) * MDP(1) + t244 * MDP(2) + t245 * MDP(3) + (t246 * t281 + t232) * MDP(6) + (t236 * pkin(1) + (t262 * t281 + t232) * qJ(2)) * MDP(7) + (t166 * t190 - t182 * t183) * MDP(8) + (-t166 * t189 - t167 * t190 + t181 * t183 - t182 * t184) * MDP(9) + (-qJD(3) * t183 + qJDD(3) * t190) * MDP(10) + (-qJD(3) * t184 - qJDD(3) * t189) * MDP(11) + (t231 * qJD(3) + t248 * qJDD(3) + t210 * t167 + t194 * t184 + t193 * t189 + t244 * t214) * MDP(13) + (-t233 * qJD(3) - t263 * qJDD(3) + t210 * t166 - t194 * t183 + t193 * t190 - t244 * t213) * MDP(14) + (t145 * t169 + t147 * t238) * MDP(15) + (-t145 * t168 - t146 * t169 - t147 * t159 - t148 * t238) * MDP(16) + (t147 * t220 + t169 * t216) * MDP(17) + (-t148 * t220 - t168 * t216) * MDP(18) + (t159 * t277 + t172 * t146 + t155 * t168 + t170 * t148 + (-t240 * qJD(4) - t149 * t223 + t150 * t226) * t220 + t241 * t216 + t244 * t209) * MDP(20) + (t238 * t277 + t172 * t145 + t155 * t169 + t170 * t147 - (t241 * qJD(4) + t149 * t226 + t150 * t223) * t220 - t240 * t216 - t244 * t208) * MDP(21) + (-MDP(5) * t221 + t265) * (t236 + t273); t280 * MDP(7) - t204 * MDP(13) + t252 * MDP(14) + (t146 + t270) * MDP(20) + (t145 - t269) * MDP(21) + (-t265 - pkin(1) * MDP(7) + (MDP(13) * t224 + MDP(5)) * t221) * qJDD(1) + ((qJD(1) * t266 + t222 * t260 + t182) * MDP(13) + (-t181 - t251) * MDP(14)) * qJD(3) + (-MDP(7) * qJ(2) - MDP(6)) * qJD(1) ^ 2 * t262; t182 * t181 * MDP(8) + (-t181 ^ 2 + t182 ^ 2) * MDP(9) + (t252 + (t181 - t251) * qJD(3)) * MDP(10) + t243 * MDP(11) + qJDD(3) * MDP(12) + (-g(3) * t214 - t194 * t182 + t245 * t213 + t249) * MDP(13) + (g(3) * t213 + t194 * t181 + t245 * t214 - t239) * MDP(14) + (t145 + t269) * MDP(17) + (-t146 + t270) * MDP(18) + (-(-t153 * t223 - t271) * t220 + t242 * qJD(4) + (-t159 * t182 + t216 * t226 - t220 * t258) * pkin(3) + t279) * MDP(20) + ((-t154 * t220 - t143) * t223 + (-qJD(4) * t152 + t153 * t220 - t144) * t226 + (-t182 * t238 - t216 * t223 - t220 * t257) * pkin(3) + t285) * MDP(21) + t284; (t253 + t269) * MDP(17) + (-t250 + t270) * MDP(18) + (-t242 * t220 + t279) * MDP(20) + (-t226 * t144 - t223 * t143 + (-t154 * t223 + t272) * t220 + t285) * MDP(21) + (-MDP(17) * t268 - t238 * MDP(18) + t242 * MDP(20) - MDP(21) * t272) * qJD(4) + t284;];
tau = t1;
