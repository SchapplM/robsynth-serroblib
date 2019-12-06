% Calculate vector of inverse dynamics joint torques for
% S5PRRPR1
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
%   pkin=[a2,a3,a4,a5,d2,d3,d5,theta1,theta4]';
% MDP [18x1]
%   Minimal dynamic parameter vector (fixed base model)
%   see S5PRRPR1_convert_par2_MPV_fixb.m
% 
% Output:
% tau [5x1]
%   joint torques of inverse dynamics (contains inertial, gravitational Coriolis and centrifugal forces)

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-05 16:16
% Revision: 77da58f92bca3eff71542919beafa37024070d86 (2019-12-05)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tau = S5PRRPR1_invdynJ_fixb_mdp_slag_vp(qJ, qJD, qJDD, g, pkin, MDP)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(5,1),zeros(3,1),zeros(9,1),zeros(18,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5PRRPR1_invdynJ_fixb_mdp_slag_vp: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5PRRPR1_invdynJ_fixb_mdp_slag_vp: qJD has to be [5x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [5 1]), ...
  'S5PRRPR1_invdynJ_fixb_mdp_slag_vp: qJDD has to be [5x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S5PRRPR1_invdynJ_fixb_mdp_slag_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S5PRRPR1_invdynJ_fixb_mdp_slag_vp: pkin has to be [9x1] (double)');
assert(isreal(MDP) && all(size(MDP) == [18 1]), ...
  'S5PRRPR1_invdynJ_fixb_mdp_slag_vp: MDP has to be [18x1] (double)'); 

%% Symbolic Calculation
% From invdyn_joint_fixb_mdp_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-05 16:15:56
% EndTime: 2019-12-05 16:15:59
% DurationCPUTime: 1.10s
% Computational Cost: add. (826->185), mult. (1138->232), div. (0->0), fcn. (739->12), ass. (0->98)
t221 = sin(pkin(9));
t222 = cos(pkin(9));
t258 = t221 ^ 2 + t222 ^ 2;
t224 = sin(qJ(3));
t256 = qJD(3) * t224;
t251 = pkin(2) * t256;
t226 = cos(qJ(3));
t277 = pkin(2) * t226;
t259 = -qJD(2) * t251 + qJDD(2) * t277;
t246 = qJDD(4) - t259;
t215 = qJDD(2) + qJDD(3);
t276 = pkin(3) * t215;
t163 = t246 - t276;
t219 = pkin(8) + qJ(2);
t212 = qJ(3) + t219;
t198 = cos(t212);
t275 = g(2) * t198;
t280 = t163 + t275;
t223 = sin(qJ(5));
t225 = cos(qJ(5));
t173 = t221 * t225 + t222 * t223;
t220 = qJD(2) + qJD(3);
t166 = t173 * t220;
t279 = t258 * t220;
t254 = qJDD(2) * t224;
t255 = qJD(3) * t226;
t158 = t215 * qJ(4) + t220 * qJD(4) + (qJD(2) * t255 + t254) * pkin(2);
t205 = t222 * qJDD(1);
t150 = -t158 * t221 + t205;
t151 = t221 * qJDD(1) + t222 * t158;
t239 = -t150 * t221 + t151 * t222;
t197 = sin(t212);
t193 = g(1) * t197;
t278 = t193 + t259;
t274 = pkin(2) * qJD(2);
t253 = t224 * t274;
t238 = t258 * (qJ(4) * t220 + t253);
t168 = t173 * qJD(5);
t257 = qJD(2) * t226;
t242 = -pkin(2) * t257 + qJD(4);
t213 = t222 * pkin(7);
t218 = pkin(9) + qJ(5);
t208 = sin(t218);
t270 = t197 * t208;
t210 = cos(t218);
t269 = t197 * t210;
t268 = t198 * t208;
t267 = t198 * t210;
t266 = t221 * t223;
t263 = t225 * t222;
t262 = t280 * t221;
t261 = t198 * pkin(3) + t197 * qJ(4);
t260 = g(1) * t198 + g(2) * t197;
t250 = t220 * t266;
t249 = t220 * t263;
t248 = qJD(5) * t249 + t173 * t215;
t200 = -t222 * pkin(4) - pkin(3);
t247 = t220 * t256;
t245 = -pkin(3) * t197 + t198 * qJ(4);
t244 = t258 * t215;
t243 = t220 * t253;
t186 = t215 * t263;
t241 = -t215 * t266 + t186;
t139 = -qJD(5) * t250 + t248;
t140 = t168 * t220 - t241;
t172 = -t263 + t266;
t167 = t172 * qJD(5);
t152 = -qJD(5) * t167 + qJDD(5) * t173;
t153 = -qJD(5) * t168 - qJDD(5) * t172;
t164 = -t249 + t250;
t240 = (-t139 * t172 - t140 * t173 + t164 * t167 - t166 * t168) * MDP(13) + (t139 * t173 - t166 * t167) * MDP(12) + t152 * MDP(14) + t153 * MDP(15) + t215 * MDP(5);
t199 = pkin(2) * t224 + qJ(4);
t170 = (-pkin(7) - t199) * t221;
t171 = t199 * t222 + t213;
t237 = t170 * t225 - t171 * t223;
t236 = t170 * t223 + t171 * t225;
t182 = (-pkin(7) - qJ(4)) * t221;
t183 = qJ(4) * t222 + t213;
t235 = t182 * t225 - t183 * t223;
t234 = t182 * t223 + t183 * t225;
t233 = -t275 + t278;
t232 = -t260 + t239;
t156 = t200 * t215 + t246;
t162 = t200 * t220 + t242;
t231 = -g(1) * t270 + g(2) * t268 + t156 * t173 - t162 * t167;
t230 = g(1) * t269 - g(2) * t267 + t156 * t172 + t162 * t168;
t229 = t243 - t275;
t201 = -pkin(3) - t277;
t228 = pkin(2) * t247 + t201 * t215;
t211 = cos(t219);
t209 = sin(t219);
t190 = pkin(2) * t255 + qJD(4);
t181 = t200 - t277;
t180 = t222 * t193;
t174 = -t220 * pkin(3) + t242;
t142 = t213 * t215 + t151;
t141 = t205 + (-pkin(7) * t215 - t158) * t221;
t1 = [(qJDD(1) - g(3)) * MDP(1) + (t150 * t222 + t151 * t221 - g(3)) * MDP(11) + t153 * MDP(17) - t152 * MDP(18); qJDD(2) * MDP(2) + (g(1) * t209 - g(2) * t211) * MDP(3) + (g(1) * t211 + g(2) * t209) * MDP(4) + ((t215 * t226 - t247) * pkin(2) + t233) * MDP(6) + (((-qJDD(2) - t215) * t224 + (-qJD(2) - t220) * t255) * pkin(2) + t260) * MDP(7) + (t180 + (-t228 - t280) * t222) * MDP(8) + ((t228 - t193) * t221 + t262) * MDP(9) + (t190 * t279 + t199 * t244 + t232) * MDP(10) + (t163 * t201 + t174 * t251 - g(1) * (-pkin(2) * t209 + t245) - g(2) * (pkin(2) * t211 + t261) + t239 * t199 + t238 * t190) * MDP(11) + (t164 * t251 + t181 * t140 + t237 * qJDD(5) + (-qJD(5) * t236 - t173 * t190) * qJD(5) + t230) * MDP(17) + (t166 * t251 + t181 * t139 - t236 * qJDD(5) + (-qJD(5) * t237 + t172 * t190) * qJD(5) + t231) * MDP(18) + t240; (t229 + t278) * MDP(6) + ((-t254 + (-qJD(3) + t220) * t257) * pkin(2) + t260) * MDP(7) + (t180 + (-t163 + t229 + t276) * t222) * MDP(8) + ((-t193 - t243 - t276) * t221 + t262) * MDP(9) + (qJ(4) * t244 + t242 * t279 + t232) * MDP(10) + (-t163 * pkin(3) - g(1) * t245 - g(2) * t261 + t238 * qJD(4) + t239 * qJ(4) + (-t174 * t224 - t226 * t238) * t274) * MDP(11) + (t200 * t140 + t235 * qJDD(5) + (-qJD(4) * t173 - qJD(5) * t234) * qJD(5) + (-t224 * t164 + t226 * t168) * t274 + t230) * MDP(17) + (t200 * t139 - t234 * qJDD(5) + (qJD(4) * t172 - qJD(5) * t235) * qJD(5) + (-t224 * t166 - t167 * t226) * t274 + t231) * MDP(18) + t240; (qJDD(4) - t233) * MDP(11) - t186 * MDP(17) + t248 * MDP(18) + (-pkin(3) * MDP(11) - t222 * MDP(8) + (t223 * MDP(17) + MDP(9)) * t221) * t215 + (0.2e1 * t166 * MDP(17) + (-t164 - t250) * MDP(18)) * qJD(5) + (-MDP(10) * t279 - t238 * MDP(11)) * t220; t166 * t164 * MDP(12) + (-t164 ^ 2 + t166 ^ 2) * MDP(13) + t241 * MDP(15) + qJDD(5) * MDP(16) + (g(1) * t268 + g(2) * t270 - g(3) * t210 + t225 * t141 - t223 * t142 - t162 * t166) * MDP(17) + (g(1) * t267 + g(2) * t269 + g(3) * t208 - t223 * t141 - t225 * t142 + t162 * t164) * MDP(18) + (t248 + (t164 - t250) * qJD(5)) * MDP(14);];
tau = t1;
