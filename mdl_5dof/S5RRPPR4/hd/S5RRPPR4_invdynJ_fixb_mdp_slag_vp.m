% Calculate vector of inverse dynamics joint torques for
% S5RRPPR4
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
%   pkin=[a2,a3,a4,a5,d1,d2,d5,theta4]';
% MDP [19x1]
%   Minimal dynamic parameter vector (fixed base model)
%   see S5RRPPR4_convert_par2_MPV_fixb.m
% 
% Output:
% tau [5x1]
%   joint torques of inverse dynamics (contains inertial, gravitational Coriolis and centrifugal forces)

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-31 19:28
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tau = S5RRPPR4_invdynJ_fixb_mdp_slag_vp(qJ, qJD, qJDD, g, pkin, MDP)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(5,1),zeros(3,1),zeros(8,1),zeros(19,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RRPPR4_invdynJ_fixb_mdp_slag_vp: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RRPPR4_invdynJ_fixb_mdp_slag_vp: qJD has to be [5x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [5 1]), ...
  'S5RRPPR4_invdynJ_fixb_mdp_slag_vp: qJDD has to be [5x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S5RRPPR4_invdynJ_fixb_mdp_slag_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S5RRPPR4_invdynJ_fixb_mdp_slag_vp: pkin has to be [8x1] (double)');
assert(isreal(MDP) && all(size(MDP) == [19 1]), ...
  'S5RRPPR4_invdynJ_fixb_mdp_slag_vp: MDP has to be [19x1] (double)'); 

%% Symbolic Calculation
% From invdyn_joint_fixb_mdp_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 19:27:52
% EndTime: 2019-12-31 19:27:55
% DurationCPUTime: 1.11s
% Computational Cost: add. (945->193), mult. (1168->237), div. (0->0), fcn. (618->10), ass. (0->100)
t219 = qJ(1) + qJ(2);
t210 = sin(t219);
t211 = cos(t219);
t220 = sin(pkin(8));
t221 = cos(pkin(8));
t170 = -t210 * t221 + t211 * t220;
t171 = t210 * t220 + t211 * t221;
t280 = g(1) * t170 + g(2) * t171;
t226 = cos(qJ(2));
t272 = pkin(1) * qJD(1);
t254 = t226 * t272;
t223 = sin(qJ(2));
t256 = t223 * t272;
t168 = t220 * t254 - t221 * t256;
t215 = qJD(1) + qJD(2);
t279 = (t220 * qJD(3) - t168) * t215;
t278 = -g(1) * t211 - g(2) * t210;
t214 = qJDD(1) + qJDD(2);
t228 = -pkin(2) - pkin(3);
t253 = qJD(2) * t272;
t271 = pkin(1) * qJDD(1);
t262 = -t223 * t253 + t226 * t271;
t250 = -qJDD(3) + t262;
t154 = t228 * t214 - t250;
t206 = t214 * qJ(3);
t208 = t215 * qJD(3);
t265 = t223 * t271 + t226 * t253;
t155 = t206 + t208 + t265;
t145 = t220 * t154 + t221 * t155;
t235 = -g(1) * t171 + g(2) * t170 + t145;
t243 = qJD(3) - t254;
t167 = t228 * t215 + t243;
t178 = qJ(3) * t215 + t256;
t148 = t167 * t221 - t178 * t220;
t146 = pkin(4) * t215 - t148;
t203 = -pkin(1) * t226 - pkin(2);
t192 = -pkin(3) + t203;
t195 = pkin(1) * t223 + qJ(3);
t161 = t192 * t221 - t195 * t220;
t156 = pkin(4) - t161;
t162 = t220 * t192 + t221 * t195;
t157 = -pkin(7) + t162;
t259 = qJD(2) * t226;
t189 = pkin(1) * t259 + qJD(3);
t260 = qJD(2) * t223;
t255 = pkin(1) * t260;
t165 = t189 * t221 + t220 * t255;
t277 = -qJDD(5) * t157 + (-t156 * t215 - t146 - t165) * qJD(5);
t224 = sin(qJ(1));
t276 = pkin(1) * t224;
t275 = pkin(2) * t214;
t149 = t220 * t167 + t221 * t178;
t270 = t149 * t221;
t164 = t189 * t220 - t221 * t255;
t269 = t164 * t215;
t222 = sin(qJ(5));
t225 = cos(qJ(5));
t268 = t222 * t225;
t267 = qJDD(4) + g(3);
t266 = t280 * t222;
t264 = t211 * pkin(2) + t210 * qJ(3);
t263 = g(1) * t210 - g(2) * t211;
t180 = t221 * qJ(3) + t220 * t228;
t217 = t222 ^ 2;
t261 = -t225 ^ 2 + t217;
t258 = qJD(5) * t215;
t227 = cos(qJ(1));
t252 = t227 * pkin(1) + t264;
t251 = t215 * t260;
t194 = t211 * qJ(3);
t249 = -pkin(2) * t210 + t194;
t144 = t154 * t221 - t220 * t155;
t169 = (t220 * t223 + t221 * t226) * t272;
t247 = t221 * qJD(3) - t169;
t246 = -t262 - t263;
t245 = t265 + t278;
t244 = t228 * t210 + t194;
t229 = qJD(5) ^ 2;
t181 = -qJDD(5) * t222 - t225 * t229;
t182 = qJDD(5) * t225 - t222 * t229;
t240 = 0.2e1 * (t214 * t268 - t261 * t258) * MDP(14) + (t214 * t217 + 0.2e1 * t258 * t268) * MDP(13) + t181 * MDP(15) + t214 * MDP(4) - t182 * MDP(16);
t179 = -qJ(3) * t220 + t221 * t228;
t174 = pkin(4) - t179;
t175 = -pkin(7) + t180;
t239 = -t174 * t214 + t175 * t229;
t238 = qJDD(3) + t246;
t236 = -MDP(18) * t225 + MDP(19) * t222 - MDP(10);
t234 = -t156 * t214 + t157 * t229 - t269;
t233 = pkin(7) * t214 + t146 * t215 - t235;
t232 = t215 * t254 - t245;
t231 = -t280 + t279;
t230 = -qJDD(5) * t175 + (-t174 * t215 - t146 - t247) * qJD(5);
t213 = t215 ^ 2;
t197 = t211 * pkin(3);
t185 = t215 * t256;
t173 = -pkin(2) * t215 + t243;
t166 = -t250 - t275;
t142 = pkin(4) * t214 - t144;
t141 = t142 * t225;
t1 = [qJDD(1) * MDP(1) + (t145 * t162 + t149 * t165 + t144 * t161 - t148 * t164 - g(1) * (t244 - t276) - g(2) * (t197 + t252)) * MDP(12) + t141 * MDP(18) + t266 * MDP(19) + t240 + ((t214 * t226 - t251) * pkin(1) - t246) * MDP(5) + ((-t214 * t223 - t215 * t259) * pkin(1) - t245) * MDP(6) + (g(1) * t224 - g(2) * t227) * MDP(2) + (g(1) * t227 + g(2) * t224) * MDP(3) + (-t161 * t214 - t144 + t269 - t280) * MDP(10) + (t155 * t195 + t178 * t189 + t166 * t203 + t173 * t255 - g(1) * (t249 - t276) - g(2) * t252) * MDP(9) + (t189 * t215 + t195 * t214 + t155 + t278) * MDP(8) + (t162 * t214 + t165 * t215 + t235) * MDP(11) + ((-t234 - t280) * MDP(18) + t277 * MDP(19)) * t225 + ((-t142 + t234) * MDP(19) + t277 * MDP(18)) * t222 + (-pkin(1) * t251 + (pkin(2) - t203) * t214 - t238) * MDP(7); (t185 - t246) * MDP(5) + t232 * MDP(6) + (t185 - t238 + 0.2e1 * t275) * MDP(7) + (0.2e1 * t206 + 0.2e1 * t208 - t232) * MDP(8) + (t155 * qJ(3) + t178 * qJD(3) - t166 * pkin(2) - g(1) * t249 - g(2) * t264 + (-t173 * t223 - t178 * t226) * t272) * MDP(9) + (-t179 * t214 - t144 + t231) * MDP(10) + (t180 * t214 + t247 * t215 + t235) * MDP(11) + (t145 * t180 + t144 * t179 - t149 * t169 + t148 * t168 - g(1) * t244 - g(2) * (t197 + t264) + (-t148 * t220 + t270) * qJD(3)) * MDP(12) + (t141 + t230 * t222 + (t231 - t239) * t225) * MDP(18) + (t230 * t225 + (-t142 + t239 - t279) * t222 + t266) * MDP(19) + t240; t238 * MDP(9) + (t144 * t221 - t263) * MDP(12) + (t145 * MDP(12) + t181 * MDP(18) - t182 * MDP(19)) * t220 + (t220 * MDP(11) - pkin(2) * MDP(9) + t236 * t221 - MDP(7)) * t214 + (-t178 * MDP(9) - MDP(12) * t270 + (-t221 * MDP(11) - MDP(8)) * t215 + (t148 * MDP(12) + t236 * t215) * t220 + 0.2e1 * (MDP(18) * t222 + MDP(19) * t225) * t221 * qJD(5)) * t215; t267 * MDP(12) + t182 * MDP(18) + t181 * MDP(19); qJDD(5) * MDP(17) + t261 * MDP(14) * t213 + (-t214 * MDP(16) + t267 * MDP(18) + t233 * MDP(19)) * t225 + (-t213 * t225 * MDP(13) - t214 * MDP(15) + t233 * MDP(18) - t267 * MDP(19)) * t222;];
tau = t1;
