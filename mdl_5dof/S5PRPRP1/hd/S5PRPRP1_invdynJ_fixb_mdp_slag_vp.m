% Calculate vector of inverse dynamics joint torques for
% S5PRPRP1
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
% MDP [19x1]
%   Minimal dynamic parameter vector (fixed base model)
%   see S5PRPRP1_convert_par2_MPV_fixb.m
% 
% Output:
% tau [5x1]
%   joint torques of inverse dynamics (contains inertial, gravitational Coriolis and centrifugal forces)

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-05 15:29
% Revision: 77da58f92bca3eff71542919beafa37024070d86 (2019-12-05)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tau = S5PRPRP1_invdynJ_fixb_mdp_slag_vp(qJ, qJD, qJDD, g, pkin, MDP)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(5,1),zeros(3,1),zeros(8,1),zeros(19,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5PRPRP1_invdynJ_fixb_mdp_slag_vp: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5PRPRP1_invdynJ_fixb_mdp_slag_vp: qJD has to be [5x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [5 1]), ...
  'S5PRPRP1_invdynJ_fixb_mdp_slag_vp: qJDD has to be [5x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S5PRPRP1_invdynJ_fixb_mdp_slag_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S5PRPRP1_invdynJ_fixb_mdp_slag_vp: pkin has to be [8x1] (double)');
assert(isreal(MDP) && all(size(MDP) == [19 1]), ...
  'S5PRPRP1_invdynJ_fixb_mdp_slag_vp: MDP has to be [19x1] (double)'); 

%% Symbolic Calculation
% From invdyn_joint_fixb_mdp_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-05 15:28:45
% EndTime: 2019-12-05 15:28:48
% DurationCPUTime: 1.18s
% Computational Cost: add. (973->206), mult. (1965->244), div. (0->0), fcn. (1389->8), ass. (0->98)
t219 = pkin(7) + qJ(2);
t213 = sin(t219);
t215 = cos(t219);
t284 = g(1) * t213 - g(2) * t215;
t272 = qJDD(2) * pkin(2);
t235 = -qJDD(3) + t272 + t284;
t244 = g(1) * t215 + g(2) * t213;
t259 = qJD(2) * qJD(3);
t283 = qJ(3) * qJDD(2) + t259;
t220 = sin(pkin(8));
t221 = cos(pkin(8));
t223 = sin(qJ(4));
t277 = cos(qJ(4));
t185 = t277 * t220 + t223 * t221;
t180 = t185 * qJD(2);
t282 = MDP(14) + MDP(16);
t281 = -MDP(15) + MDP(18);
t253 = t277 * t221;
t232 = -t223 * t220 + t253;
t273 = pkin(6) + qJ(3);
t193 = t273 * t220;
t194 = t273 * t221;
t233 = -t277 * t193 - t223 * t194;
t151 = t232 * qJD(3) + t233 * qJD(4);
t167 = -t223 * t193 + t277 * t194;
t218 = pkin(8) + qJ(4);
t212 = sin(t218);
t280 = -qJD(4) * t151 - qJDD(4) * t167 - t212 * t284;
t214 = cos(t218);
t209 = t221 * qJDD(1);
t168 = t209 + (-t273 * qJDD(2) - t259) * t220;
t258 = qJDD(2) * t221;
t174 = qJ(3) * t258 + t220 * qJDD(1) + t221 * t259;
t169 = pkin(6) * t258 + t174;
t211 = t221 * qJD(1);
t175 = -qJD(2) * t193 + t211;
t250 = qJD(4) * t277;
t255 = t223 * t168 + t277 * t169 + t175 * t250;
t279 = -g(3) * t212 - t244 * t214 + t255;
t278 = t180 ^ 2;
t271 = qJDD(4) * pkin(4);
t245 = qJD(2) * t253;
t265 = qJD(2) * t220;
t252 = t223 * t265;
t178 = -t245 + t252;
t270 = t178 * t180;
t269 = t220 * MDP(6);
t264 = qJD(2) * t221;
t187 = qJ(3) * t264 + t220 * qJD(1);
t176 = pkin(6) * t264 + t187;
t268 = t223 * t176;
t183 = t185 * qJD(4);
t248 = qJDD(2) * t277;
t257 = qJDD(2) * t223;
t240 = t220 * t257 - t221 * t248;
t158 = qJD(2) * t183 + t240;
t262 = qJD(4) * t223;
t251 = t220 * t262;
t182 = -t221 * t250 + t251;
t267 = -t185 * t158 + t182 * t178;
t266 = t220 ^ 2 + t221 ^ 2;
t154 = t223 * t175 + t277 * t176;
t263 = qJD(4) * t154;
t153 = t277 * t175 - t268;
t261 = qJD(5) - t153;
t256 = qJDD(4) * qJ(5);
t254 = qJD(4) * t245 + t220 * t248 + t221 * t257;
t206 = pkin(3) * t221 + pkin(2);
t247 = t153 + t268;
t246 = -t277 * t168 + t223 * t169 + t175 * t262 + t176 * t250;
t242 = pkin(4) * t214 + qJ(5) * t212;
t157 = qJD(2) * t251 - t254;
t239 = t157 * t232 + t180 * t183;
t238 = (-qJ(3) * t265 + t211) * t220 - t187 * t221;
t160 = -qJD(4) * t183 + qJDD(4) * t232;
t234 = t206 + t242;
t192 = -t206 * qJD(2) + qJD(3);
t188 = -t206 * qJDD(2) + qJDD(3);
t230 = -g(3) * t214 + t244 * t212 - t246;
t173 = -t220 * t283 + t209;
t229 = -t173 * t220 + t174 * t221 - t244;
t152 = t185 * qJD(3) + t167 * qJD(4);
t228 = -qJD(4) * t152 + qJDD(4) * t233 + t284 * t214;
t227 = pkin(4) * t158 + qJ(5) * t157 + t188;
t149 = pkin(4) * t178 - qJ(5) * t180 + t192;
t226 = t149 * t180 + qJDD(5) - t230;
t177 = t178 ^ 2;
t159 = -qJD(4) * t182 + qJDD(4) * t185;
t156 = -pkin(4) * t232 - qJ(5) * t185 - t206;
t155 = pkin(4) * t180 + qJ(5) * t178;
t150 = qJD(4) * qJ(5) + t154;
t148 = -qJD(4) * pkin(4) + t261;
t146 = pkin(4) * t183 + qJ(5) * t182 - qJD(5) * t185;
t145 = (t178 - t252) * qJD(4) + t254;
t143 = -qJD(5) * t180 + t227;
t142 = qJDD(5) + t246 - t271;
t141 = t256 + (qJD(5) - t268) * qJD(4) + t255;
t1 = [(qJDD(1) - g(3)) * MDP(1) + (t173 * t221 + t174 * t220 - g(3)) * MDP(8) + (t239 + t267) * MDP(17) + (t141 * t185 - t142 * t232 + t148 * t183 - t150 * t182 - g(3)) * MDP(19) + t281 * t159 + t282 * t160; qJDD(2) * MDP(2) + t284 * MDP(3) + t244 * MDP(4) + (t266 * t283 + t229) * MDP(7) + (t235 * pkin(2) + t229 * qJ(3) - t238 * qJD(3)) * MDP(8) + (-t157 * t185 - t180 * t182) * MDP(9) + (-t239 + t267) * MDP(10) + t159 * MDP(11) + t160 * MDP(12) + (-t158 * t206 + t183 * t192 - t188 * t232 + t228) * MDP(14) + (t157 * t206 - t182 * t192 + t185 * t188 + t280) * MDP(15) + (-t143 * t232 + t146 * t178 + t149 * t183 + t156 * t158 + t228) * MDP(16) + (t141 * t232 + t142 * t185 - t148 * t182 - t150 * t183 - t151 * t178 + t152 * t180 + t157 * t233 - t158 * t167 - t244) * MDP(17) + (-t143 * t185 - t146 * t180 + t149 * t182 + t156 * t157 - t280) * MDP(18) + (t141 * t167 - t142 * t233 + t143 * t156 + t149 * t146 + t148 * t152 + t150 * t151 + (-g(1) * t273 - g(2) * t234) * t215 + (g(1) * t234 - g(2) * t273) * t213) * MDP(19) + (t221 * MDP(5) - t269) * (t235 + t272); -MDP(5) * t258 + qJDD(2) * t269 - t266 * MDP(7) * qJD(2) ^ 2 + (t238 * qJD(2) - t235) * MDP(8) + (-t177 - t278) * MDP(17) + (t150 * t178 + (-qJD(5) - t148) * t180 + t227 - t284) * MDP(19) + t282 * (0.2e1 * t180 * qJD(4) + t240) + t281 * ((t178 + t252) * qJD(4) - t254); MDP(9) * t270 + (-t177 + t278) * MDP(10) + t145 * MDP(11) - t240 * MDP(12) + qJDD(4) * MDP(13) + (-t180 * t192 + t230 + t263) * MDP(14) + (t247 * qJD(4) + t178 * t192 - t279) * MDP(15) + (-t155 * t178 - t226 + t263 + 0.2e1 * t271) * MDP(16) + (pkin(4) * t157 - qJ(5) * t158 + (t150 - t154) * t180 + (t148 - t261) * t178) * MDP(17) + (0.2e1 * t256 - t149 * t178 + t155 * t180 + (0.2e1 * qJD(5) - t247) * qJD(4) + t279) * MDP(18) + (-t142 * pkin(4) - g(3) * t242 + t141 * qJ(5) - t148 * t154 - t149 * t155 + t261 * t150 + t244 * (pkin(4) * t212 - qJ(5) * t214)) * MDP(19); (-qJDD(4) + t270) * MDP(16) + t145 * MDP(17) + (-qJD(4) ^ 2 - t278) * MDP(18) + (-qJD(4) * t150 + t226 - t271) * MDP(19);];
tau = t1;
