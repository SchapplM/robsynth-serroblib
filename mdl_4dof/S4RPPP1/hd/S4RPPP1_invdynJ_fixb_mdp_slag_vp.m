% Calculate vector of inverse dynamics joint torques for
% S4RPPP1
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
%   pkin=[a2,a3,a4,alpha2,d1,theta2]';
% MDP [15x1]
%   Minimal dynamic parameter vector (fixed base model)
%   see S4RPPP1_convert_par2_MPV_fixb.m
% 
% Output:
% tau [4x1]
%   joint torques of inverse dynamics (contains inertial, gravitational Coriolis and centrifugal forces)

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-08 18:26
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tau = S4RPPP1_invdynJ_fixb_mdp_slag_vp(qJ, qJD, qJDD, g, pkin, MDP)
%% Coder Information
%#codegen
%$cgargs {zeros(4,1),zeros(4,1),zeros(4,1),zeros(3,1),zeros(6,1),zeros(15,1)}
assert(isreal(qJ) && all(size(qJ) == [4 1]), ...
  'S4RPPP1_invdynJ_fixb_mdp_slag_vp: qJ has to be [4x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [4 1]), ...
  'S4RPPP1_invdynJ_fixb_mdp_slag_vp: qJD has to be [4x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [4 1]), ...
  'S4RPPP1_invdynJ_fixb_mdp_slag_vp: qJDD has to be [4x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S4RPPP1_invdynJ_fixb_mdp_slag_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [6 1]), ...
  'S4RPPP1_invdynJ_fixb_mdp_slag_vp: pkin has to be [6x1] (double)');
assert(isreal(MDP) && all(size(MDP) == [15 1]), ...
  'S4RPPP1_invdynJ_fixb_mdp_slag_vp: MDP has to be [15x1] (double)'); 

%% Symbolic Calculation
% From invdyn_joint_fixb_mdp_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-08 18:26:31
% EndTime: 2019-03-08 18:26:32
% DurationCPUTime: 0.82s
% Computational Cost: add. (475->194), mult. (1291->254), div. (0->0), fcn. (939->6), ass. (0->115)
t216 = sin(pkin(4));
t212 = t216 ^ 2;
t217 = cos(pkin(6));
t285 = t212 * t217;
t284 = MDP(7) + MDP(11) + MDP(15);
t218 = cos(pkin(4));
t283 = pkin(1) * t218;
t219 = sin(qJ(1));
t282 = g(1) * t219;
t281 = qJ(3) * t218;
t215 = sin(pkin(6));
t211 = t215 ^ 2;
t280 = t211 * t212;
t279 = t215 * t216;
t278 = t216 * t217;
t277 = t216 * t219;
t220 = cos(qJ(1));
t276 = t216 * t220;
t275 = t218 * t219;
t274 = t218 * t220;
t221 = qJD(1) ^ 2;
t273 = t218 * t221;
t264 = qJD(1) * qJD(2);
t250 = t216 * t264;
t189 = t215 * t250;
t261 = qJDD(1) * t216;
t248 = t215 * t261;
t272 = qJ(2) * t248 + t189;
t202 = qJ(2) * t278;
t269 = qJD(1) * t218;
t252 = t215 * t269;
t171 = pkin(1) * t252 + qJD(1) * t202;
t182 = t215 * t283 + t202;
t271 = t220 * pkin(1) + qJ(2) * t277;
t268 = qJD(2) * t216;
t185 = qJD(3) * t218 + t217 * t268;
t270 = qJD(1) * t185;
t267 = MDP(8) + MDP(12);
t266 = MDP(9) - MDP(14);
t253 = qJD(1) * t279;
t192 = qJ(2) * t253;
t265 = qJD(3) + t192;
t263 = qJD(1) * qJD(3);
t262 = qJDD(1) * t215;
t260 = qJDD(1) * t218;
t259 = MDP(10) + MDP(13);
t258 = t217 * t283;
t257 = pkin(1) * t260;
t256 = t215 * t285;
t255 = t217 * t273;
t191 = t217 * t250;
t247 = t217 * t261;
t158 = qJ(2) * t247 + t215 * t257 + t191;
t254 = -pkin(1) * t217 - pkin(2);
t251 = t212 * t264;
t249 = t215 * t263;
t246 = qJDD(3) + t272;
t245 = -qJ(3) * t215 - pkin(1);
t244 = -pkin(1) * t219 + qJ(2) * t276;
t240 = t254 * t218;
t152 = qJDD(1) * t240 + t246;
t201 = qJ(2) * t279;
t165 = t201 + t240;
t243 = qJDD(1) * t165 + t152;
t234 = -pkin(2) * t217 + t245;
t153 = qJDD(2) + (t234 * qJDD(1) - t249) * t216;
t166 = t234 * t216;
t242 = qJDD(1) * t166 + t153;
t241 = -qJ(4) + t254;
t177 = t215 * t219 - t217 * t274;
t179 = t215 * t220 + t217 * t275;
t239 = g(1) * t177 - g(2) * t179;
t178 = t215 * t274 + t217 * t219;
t180 = -t215 * t275 + t217 * t220;
t238 = g(1) * t178 - g(2) * t180;
t237 = g(1) * t220 + g(2) * t219;
t150 = -qJ(3) * t260 - t218 * t263 - t158;
t235 = -qJD(3) * t215 - qJD(4) * t217;
t233 = t180 * pkin(2) + qJ(3) * t179 + t271;
t232 = pkin(3) * t278 + t281;
t200 = -pkin(1) * t261 + qJDD(2);
t231 = pkin(1) * qJDD(1) * t212 - t200 * t216;
t196 = pkin(3) * t248;
t147 = t196 + (-qJD(1) * qJD(4) + qJDD(1) * t241) * t218 + t246;
t222 = pkin(3) * t279 + t218 * t241;
t156 = t201 + t222;
t184 = -qJD(4) * t218 + t215 * t268;
t230 = qJD(1) * t184 + qJDD(1) * t156 + t147;
t149 = pkin(3) * t247 + qJDD(4) - t150;
t159 = t232 + t182;
t229 = qJDD(1) * t159 + t149 + t270;
t164 = -t182 - t281;
t228 = -qJDD(1) * t164 - t150 + t270;
t227 = -t178 * pkin(2) - qJ(3) * t177 + t244;
t154 = qJD(1) * t232 + qJD(4) + t171;
t163 = -qJ(3) * t269 - t171;
t226 = t163 * MDP(11) + (-qJD(4) - t154) * MDP(15);
t225 = (-pkin(2) - qJ(4)) * t217 + t245;
t148 = qJDD(2) + (qJD(1) * t235 + qJDD(1) * t225) * t216;
t160 = t225 * t216;
t176 = t235 * t216;
t224 = (-qJD(1) * t176 - qJDD(1) * t160 - t148) * t216;
t223 = -g(1) * t179 - g(2) * t177 + g(3) * t278 + t246;
t214 = t218 ^ 2;
t213 = t217 ^ 2;
t187 = t211 * t251;
t186 = t273 * t279;
t181 = -t201 + t258;
t170 = qJD(1) * t258 - t192;
t162 = qJD(1) * t166 + qJD(2);
t161 = qJD(1) * t240 + t265;
t157 = t217 * t257 - t272;
t155 = qJD(1) * t160 + qJD(2);
t151 = qJD(1) * t222 + t265;
t1 = [qJDD(1) * MDP(1) + (-g(2) * t220 + t282) * MDP(2) + t237 * MDP(3) + (t231 * t217 + (qJDD(1) * t181 + t157 - t189) * t218 + t238) * MDP(4) + (-t231 * t215 + (-qJDD(1) * t182 - t158 - t191) * t218 - t239) * MDP(5) + (t213 * t251 + t187 + (-t157 * t215 + t158 * t217 + (-t181 * t215 + t182 * t217) * qJDD(1) - t237) * t216) * MDP(6) + (t158 * t182 + t157 * t181 - g(1) * t244 - g(2) * t271 + (-t200 * pkin(1) + (-t170 * t215 + t171 * t217) * qJD(2)) * t216) * MDP(7) + (t187 + (t215 * t243 + t217 * t228 - t237) * t216) * MDP(8) + (-t249 * t285 + t243 * t218 + (qJD(2) * t252 + t217 * t242) * t216 - t238) * MDP(9) + (t218 * t228 - t242 * t279 + t263 * t280 + t239) * MDP(10) + (t153 * t166 + t150 * t164 - t163 * t185 + t152 * t165 - g(1) * t227 - g(2) * t233 + (qJD(2) * t161 - qJD(3) * t162) * t279) * MDP(11) + (t215 * t230 + t217 * t229 - t237) * t216 * MDP(12) + (t215 * t224 + t218 * t229 + t239) * MDP(13) + (t217 * t224 - t218 * t230 + t238) * MDP(14) + (t148 * t160 + t155 * t176 + t147 * t156 + t151 * t184 + t149 * t159 + t154 * t185 - g(1) * (pkin(3) * t276 - qJ(4) * t178 + t227) - g(2) * (pkin(3) * t277 + qJ(4) * t180 + t233)) * MDP(15); (MDP(6) + t267) * (-t211 - t213) * t221 * t212 + (-MDP(4) + t266) * (-t186 + t247) + ((MDP(5) - t259) * (t255 + t262) - t284 * t282 + (-pkin(1) * MDP(7) + t234 * MDP(11) + (-qJ(4) * t217 + t234) * MDP(15)) * qJDD(1) + ((-t171 * MDP(7) + t226) * t217 + (t170 * MDP(7) + (-qJD(3) - t161) * MDP(11) + (-qJD(3) - t151) * MDP(15)) * t215) * qJD(1)) * t216 + t284 * (g(2) * t276 - g(3) * t218 + qJDD(2)); (t162 * t253 + t223) * MDP(11) + (t155 * t253 + t196 + t223) * MDP(15) + t259 * (-t214 - t280) * t221 + t266 * (t221 * t256 + t260) + t267 * (-t255 + t262) * t216 + (t226 * qJD(1) + (MDP(11) * t254 + MDP(15) * t241) * qJDD(1)) * t218; (t186 + t247) * MDP(12) + MDP(13) * t260 + (-g(3) * t279 - g(1) * t180 - g(2) * t178 + (t151 * t218 + t155 * t278) * qJD(1) + t149) * MDP(15) + (-MDP(13) * t256 + (-t212 * t213 - t214) * MDP(14)) * t221;];
tau  = t1;
