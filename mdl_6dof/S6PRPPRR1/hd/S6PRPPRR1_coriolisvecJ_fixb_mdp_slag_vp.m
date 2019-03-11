% Calculate minimal parameter regressor of Coriolis joint torque vector for
% S6PRPPRR1
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% qJD [6x1]
%   Generalized joint velocities
% pkin [12x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,alpha2,d2,d5,d6,theta1,theta3,theta4]';
% MDP [23x1]
%   Minimal dynamic parameter vector (fixed base model)
%   see S6PRPPRR1_convert_par2_MPV_fixb.m
% 
% Output:
% tauc [6x1]
%   joint torques required to compensate Coriolis and centrifugal load

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-08 19:16
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauc = S6PRPPRR1_coriolisvecJ_fixb_mdp_slag_vp(qJ, qJD, pkin, MDP)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(12,1),zeros(23,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6PRPPRR1_coriolisvecJ_fixb_mdp_slag_vp: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6PRPPRR1_coriolisvecJ_fixb_mdp_slag_vp: qJD has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [12 1]), ...
  'S6PRPPRR1_coriolisvecJ_fixb_mdp_slag_vp: pkin has to be [12x1] (double)');
assert(isreal(MDP) && all(size(MDP) == [23 1]), ...
  'S6PRPPRR1_coriolisvecJ_fixb_mdp_slag_vp: MDP has to be [23x1] (double)'); 

%% Symbolic Calculation
% From coriolisvec_joint_fixb_mdp_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-08 19:16:17
% EndTime: 2019-03-08 19:16:22
% DurationCPUTime: 2.20s
% Computational Cost: add. (1578->249), mult. (4167->355), div. (0->0), fcn. (3419->12), ass. (0->126)
t239 = sin(pkin(12));
t242 = cos(pkin(12));
t310 = t242 * MDP(6) - t239 * MDP(7);
t249 = cos(qJ(5));
t282 = qJD(2) * t242;
t231 = t249 * t282;
t246 = sin(qJ(5));
t292 = t239 * t246;
t272 = qJD(2) * t292;
t213 = -t231 + t272;
t212 = qJD(6) + t213;
t240 = sin(pkin(11));
t247 = sin(qJ(2));
t241 = sin(pkin(6));
t283 = qJD(1) * t241;
t274 = t247 * t283;
t228 = t240 * t274;
t243 = cos(pkin(11));
t250 = cos(qJ(2));
t273 = t250 * t283;
t266 = t243 * t273;
t207 = -t228 + t266;
t307 = t207 - qJD(4);
t284 = t239 ^ 2 + t242 ^ 2;
t309 = MDP(8) * t284;
t223 = qJD(2) * t266;
t194 = t223 + (qJD(4) - t228) * qJD(2);
t221 = -t249 * t242 + t292;
t308 = t194 * t221;
t208 = (t240 * t247 - t243 * t250) * t241;
t224 = qJD(2) * pkin(2) + t273;
t229 = t243 * t274;
t196 = t240 * t224 + t229;
t193 = qJD(2) * qJ(4) + t196;
t244 = cos(pkin(6));
t232 = qJD(1) * t244 + qJD(3);
t226 = t242 * t232;
t175 = t226 + (-pkin(8) * qJD(2) - t193) * t239;
t180 = t242 * t193 + t239 * t232;
t176 = pkin(8) * t282 + t180;
t164 = t175 * t246 + t176 * t249;
t222 = t239 * t249 + t242 * t246;
t158 = qJD(5) * t164 + t194 * t222;
t215 = t222 * qJD(2);
t306 = t212 * (pkin(5) * t215 + t212 * pkin(9)) + t158;
t163 = t175 * t249 - t176 * t246;
t157 = qJD(5) * t163 - t308;
t161 = -qJD(5) * pkin(5) - t163;
t195 = t224 * t243 - t228;
t265 = qJD(4) - t195;
t275 = -pkin(4) * t242 - pkin(3);
t186 = qJD(2) * t275 + t265;
t165 = pkin(5) * t213 - pkin(9) * t215 + t186;
t304 = pkin(2) * t243;
t227 = t275 - t304;
t182 = pkin(5) * t221 - pkin(9) * t222 + t227;
t233 = pkin(2) * t240 + qJ(4);
t303 = pkin(8) + t233;
t218 = t303 * t239;
t219 = t303 * t242;
t184 = -t218 * t246 + t219 * t249;
t217 = t222 * qJD(5);
t211 = qJD(2) * t217;
t216 = t221 * qJD(5);
t258 = -t218 * t249 - t219 * t246;
t288 = -qJD(5) * t258 - t307 * t221;
t305 = -(qJD(6) * t165 + t157) * t221 + t158 * t222 - t161 * t216 + (-qJD(6) * t182 + t288) * t212 - t184 * t211;
t245 = sin(qJ(6));
t280 = qJD(6) * t245;
t230 = qJD(5) * t231;
t210 = -qJD(5) * t272 + t230;
t248 = cos(qJ(6));
t279 = t248 * qJD(5);
t285 = qJD(6) * t279 + t248 * t210;
t177 = -t215 * t280 + t285;
t302 = t177 * t245;
t301 = t182 * t211;
t209 = (t240 * t250 + t243 * t247) * t241;
t205 = qJD(2) * t209;
t197 = qJD(1) * t205;
t300 = t197 * t208;
t294 = t215 * t245;
t199 = -t279 + t294;
t299 = t199 * t212;
t298 = t199 * t215;
t201 = qJD(5) * t245 + t215 * t248;
t297 = t201 * t212;
t296 = t201 * t215;
t295 = t210 * t245;
t290 = t245 * t211;
t203 = t248 * t211;
t289 = t177 * t221 + t201 * t217;
t287 = qJD(5) * t184 - t307 * t222;
t204 = t240 * t273 + t229;
t286 = pkin(5) * t217 + pkin(9) * t216 - t204;
t281 = qJD(6) * t222;
t278 = MDP(13) * qJD(5);
t277 = t222 * t290;
t276 = t222 * t203;
t270 = t284 * t194;
t269 = t248 * t212;
t162 = qJD(5) * pkin(9) + t164;
t156 = t162 * t248 + t165 * t245;
t264 = t162 * t245 - t165 * t248;
t191 = -t209 * t239 + t242 * t244;
t192 = t209 * t242 + t239 * t244;
t168 = t191 * t246 + t192 * t249;
t263 = t168 * t248 + t208 * t245;
t262 = -t168 * t245 + t208 * t248;
t178 = t201 * qJD(6) + t295;
t261 = -t178 * t221 - t199 * t217;
t260 = (-t193 * t239 + t226) * t239 - t180 * t242;
t259 = t191 * t249 - t192 * t246;
t256 = t203 + (-t213 * t245 - t280) * t212;
t255 = t216 * t245 - t248 * t281;
t254 = t216 * t248 + t222 * t280;
t253 = -pkin(9) * t211 + (t161 + t163) * t212;
t251 = qJD(2) ^ 2;
t206 = qJD(2) * t208;
t198 = -qJD(2) * t228 + t223;
t190 = -qJD(2) * pkin(3) + t265;
t169 = pkin(5) * t211 - pkin(9) * t210 + t197;
t166 = t248 * t169;
t160 = qJD(5) * t168 - t206 * t222;
t159 = qJD(5) * t259 + t206 * t221;
t1 = [(-t195 * t205 - t196 * t206 + t198 * t209 + t300) * MDP(5) + (t190 * t205 + t300 + t260 * t206 + (-t191 * t239 + t192 * t242) * t194) * MDP(9) + (-qJD(5) * t160 + t205 * t213 + t208 * t211) * MDP(15) + (-qJD(5) * t159 + t205 * t215 + t208 * t210) * MDP(16) + ((-qJD(6) * t263 - t159 * t245 + t205 * t248) * t212 + t262 * t211 + t160 * t199 - t259 * t178) * MDP(22) + (-(qJD(6) * t262 + t159 * t248 + t205 * t245) * t212 - t263 * t211 + t160 * t201 - t259 * t177) * MDP(23) + (-MDP(3) * t247 - MDP(4) * t250) * t251 * t241 + (-t205 * t310 - t206 * t309) * qJD(2); (t195 * t204 - t196 * t207 + (-t197 * t243 + t198 * t240) * pkin(2)) * MDP(5) + (-t307 * qJD(2) * t284 + t270) * MDP(8) + (t197 * (-pkin(3) - t304) - t190 * t204 + t233 * t270 + t307 * t260) * MDP(9) + (t210 * t222 - t215 * t216) * MDP(10) + (-t210 * t221 - t211 * t222 + t213 * t216 - t215 * t217) * MDP(11) - t216 * qJD(5) * MDP(12) - t217 * t278 + (-qJD(5) * t287 + t186 * t217 + t197 * t221 - t204 * t213 + t211 * t227) * MDP(15) + (qJD(5) * t288 - t186 * t216 + t197 * t222 - t204 * t215 + t210 * t227) * MDP(16) + (t177 * t222 * t248 - t201 * t254) * MDP(17) + (-(-t199 * t248 - t201 * t245) * t216 + (-t302 - t178 * t248 + (t199 * t245 - t201 * t248) * qJD(6)) * t222) * MDP(18) + (-t212 * t254 + t276 + t289) * MDP(19) + (t212 * t255 + t261 - t277) * MDP(20) + (t211 * t221 + t212 * t217) * MDP(21) + (-t264 * t217 + t166 * t221 - t258 * t178 + t287 * t199 + (t301 + t286 * t212 + (t161 * t222 - t162 * t221 - t184 * t212) * qJD(6)) * t248 + t305 * t245) * MDP(22) + (-t156 * t217 - t258 * t177 + t287 * t201 + (-t301 - (-qJD(6) * t162 + t169) * t221 - t161 * t281 + (qJD(6) * t184 - t286) * t212) * t245 + t305 * t248) * MDP(23) + t310 * (qJD(2) * t204 - t197); (-t261 - t277) * MDP(22) + (-t276 + t289) * MDP(23) + (MDP(22) * t255 + MDP(23) * t254) * t212 + (-MDP(15) * t217 + MDP(16) * t216) * qJD(5); t215 * qJD(5) * MDP(15) + (-qJD(5) * t213 + t230) * MDP(16) + (t256 - t298) * MDP(22) + (-t212 ^ 2 * t248 - t290 - t296) * MDP(23) - t251 * t309 + ((qJD(1) * t209 + t260) * MDP(9) + (MDP(15) * t222 - MDP(16) * t292) * qJD(5)) * qJD(2); -t213 ^ 2 * MDP(11) + (t230 + (t213 - t272) * qJD(5)) * MDP(12) + (t186 * t213 + t308) * MDP(16) + (t201 * t269 + t302) * MDP(17) + ((t177 - t299) * t248 + (-t178 - t297) * t245) * MDP(18) + (t212 * t269 + t290 - t296) * MDP(19) + (t256 + t298) * MDP(20) + (-pkin(5) * t178 - t164 * t199 + t253 * t245 - t306 * t248) * MDP(22) + (-pkin(5) * t177 - t164 * t201 + t306 * t245 + t253 * t248) * MDP(23) + (-MDP(15) * t194 - qJD(2) * t278) * t222 + (MDP(10) * t213 + t215 * MDP(11) - t186 * MDP(15) - t212 * MDP(21) + MDP(22) * t264 + t156 * MDP(23) + t278) * t215; t201 * t199 * MDP(17) + (-t199 ^ 2 + t201 ^ 2) * MDP(18) + (t285 + t299) * MDP(19) + (-t295 + t297) * MDP(20) + t211 * MDP(21) + (t156 * t212 - t157 * t245 - t161 * t201 + t166) * MDP(22) + (-t157 * t248 + t161 * t199 - t169 * t245 - t212 * t264) * MDP(23) + (-MDP(19) * t294 - MDP(20) * t201 - MDP(22) * t156 + MDP(23) * t264) * qJD(6);];
tauc  = t1;
