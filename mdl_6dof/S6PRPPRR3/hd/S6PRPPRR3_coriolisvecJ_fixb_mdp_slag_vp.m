% Calculate minimal parameter regressor of Coriolis joint torque vector for
% S6PRPPRR3
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% qJD [6x1]
%   Generalized joint velocities
% pkin [11x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,alpha2,d2,d5,d6,theta1,theta4]';
% MDP [24x1]
%   Minimal dynamic parameter vector (fixed base model)
%   see S6PRPPRR3_convert_par2_MPV_fixb.m
% 
% Output:
% tauc [6x1]
%   joint torques required to compensate Coriolis and centrifugal load

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-08 19:23
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauc = S6PRPPRR3_coriolisvecJ_fixb_mdp_slag_vp(qJ, qJD, pkin, MDP)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(11,1),zeros(24,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6PRPPRR3_coriolisvecJ_fixb_mdp_slag_vp: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6PRPPRR3_coriolisvecJ_fixb_mdp_slag_vp: qJD has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [11 1]), ...
  'S6PRPPRR3_coriolisvecJ_fixb_mdp_slag_vp: pkin has to be [11x1] (double)');
assert(isreal(MDP) && all(size(MDP) == [24 1]), ...
  'S6PRPPRR3_coriolisvecJ_fixb_mdp_slag_vp: MDP has to be [24x1] (double)'); 

%% Symbolic Calculation
% From coriolisvec_joint_fixb_mdp_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-08 19:23:06
% EndTime: 2019-03-08 19:23:11
% DurationCPUTime: 2.17s
% Computational Cost: add. (1031->261), mult. (2329->393), div. (0->0), fcn. (1631->10), ass. (0->128)
t221 = sin(qJ(5));
t224 = cos(qJ(5));
t279 = MDP(17) * t224;
t308 = MDP(16) * t221 + t279;
t225 = cos(qJ(2));
t216 = sin(pkin(6));
t277 = qJD(1) * t216;
t257 = t225 * t277;
t191 = (qJD(3) + t257) * qJD(2);
t215 = sin(pkin(11));
t217 = cos(pkin(11));
t222 = sin(qJ(2));
t276 = qJD(2) * t216;
t256 = t222 * t276;
t242 = qJD(1) * t256;
t172 = t191 * t215 - t217 * t242;
t258 = t222 * t277;
t178 = t215 * t257 - t217 * t258;
t246 = qJD(3) * t215 - t178;
t307 = t246 * qJD(2) + t172;
t227 = qJD(5) ^ 2;
t228 = qJD(2) ^ 2;
t306 = t215 * (t227 + t228);
t226 = -pkin(2) - pkin(3);
t283 = t217 * qJ(3) + t215 * t226;
t196 = -pkin(8) + t283;
t305 = -t196 * t227 + t307;
t304 = t224 * MDP(16);
t213 = t221 ^ 2;
t303 = (-t224 ^ 2 + t213) * MDP(12);
t302 = qJD(2) * pkin(2);
t301 = qJD(5) * pkin(9);
t239 = qJD(3) - t257;
t188 = t226 * qJD(2) + t239;
t200 = qJD(2) * qJ(3) + t258;
t171 = t215 * t188 + t217 * t200;
t169 = -qJD(2) * pkin(8) + t171;
t218 = cos(pkin(6));
t205 = -qJD(1) * t218 + qJD(4);
t164 = t169 * t224 + t205 * t221;
t173 = t217 * t191 + t215 * t242;
t158 = qJD(5) * t164 + t173 * t221;
t220 = sin(qJ(6));
t300 = t158 * t220;
t223 = cos(qJ(6));
t299 = t158 * t223;
t263 = qJD(2) * qJD(5);
t251 = t224 * t263;
t271 = qJD(6) * t220;
t254 = t221 * t271;
t262 = qJD(5) * qJD(6);
t284 = qJD(2) * t254 + t223 * t262;
t176 = -t223 * t251 + t284;
t298 = t176 * t220;
t297 = t176 * t224;
t207 = t220 * t262;
t270 = qJD(6) * t223;
t253 = t221 * t270;
t266 = t220 * qJD(5);
t177 = -t207 + (t224 * t266 + t253) * qJD(2);
t296 = t177 * t224;
t264 = t223 * qJD(5);
t267 = t220 * qJD(2);
t193 = t221 * t267 + t264;
t275 = qJD(2) * t224;
t206 = qJD(6) + t275;
t295 = t193 * t206;
t294 = t193 * t221;
t288 = t221 * t223;
t194 = qJD(2) * t288 - t266;
t293 = t194 * t206;
t292 = t194 * t221;
t291 = t215 * t225;
t290 = t220 * t206;
t289 = t220 * t224;
t287 = t223 * t206;
t286 = t223 * t224;
t240 = -pkin(5) * t221 + pkin(9) * t224;
t231 = t240 * qJD(5);
t285 = -t231 - t246;
t278 = MDP(7) * qJD(2);
t273 = qJD(5) * t196;
t272 = qJD(5) * t224;
t269 = qJD(6) * t224;
t237 = t169 * t221 - t205 * t224;
t161 = -qJD(5) * pkin(5) + t237;
t268 = t161 * qJD(6);
t265 = t221 * MDP(22);
t261 = 0.2e1 * t263;
t260 = t206 * t289;
t259 = t206 * t286;
t255 = t206 * t271;
t252 = t206 * t270;
t162 = t164 + t301;
t250 = t196 * t206 + t162;
t170 = t188 * t217 - t215 * t200;
t248 = -t215 * qJ(3) + t217 * t226;
t184 = (t215 * t222 + t217 * t225) * t216;
t181 = qJD(1) * t184;
t245 = qJD(3) * t217 - t181;
t244 = 0.2e1 * t251;
t243 = t221 * t252;
t195 = pkin(4) - t248;
t241 = pkin(5) * t224 + pkin(9) * t221;
t168 = qJD(2) * pkin(4) - t170;
t165 = t241 * qJD(2) + t168;
t156 = t162 * t223 + t165 * t220;
t238 = t162 * t220 - t165 * t223;
t236 = t170 * t215 - t171 * t217;
t185 = (t217 * t222 - t291) * t216;
t175 = t185 * t224 - t218 * t221;
t235 = t175 * t223 + t184 * t220;
t234 = -t175 * t220 + t184 * t223;
t174 = t185 * t221 + t218 * t224;
t232 = t223 * t213 * t263 + t206 * t254;
t230 = qJD(5) * (-qJD(2) * t195 - t168 - t245);
t157 = -t237 * qJD(5) + t173 * t224;
t229 = -t161 * qJD(5) - qJD(6) * t165 - t245 * t206 - t157;
t199 = t240 * qJD(2);
t192 = t239 - t302;
t182 = t195 + t241;
t180 = qJD(2) * t184;
t179 = -t217 * t256 + t276 * t291;
t167 = qJD(2) * t231 + t172;
t166 = t223 * t167;
t160 = -t174 * qJD(5) + t180 * t224;
t159 = t175 * qJD(5) + t180 * t221;
t1 = [(-t170 * t179 + t171 * t180 + t172 * t184 + t173 * t185) * MDP(10) + ((-t235 * qJD(6) - t160 * t220 + t179 * t223) * t206 - t159 * t193 - t174 * t177) * MDP(23) + (-(t234 * qJD(6) + t160 * t223 + t179 * t220) * t206 - t159 * t194 + t174 * t176) * MDP(24) + (-MDP(16) * t159 - MDP(17) * t160) * qJD(5) + (t180 * MDP(9) + (-t221 * MDP(17) + MDP(8) + t304) * t179 + (-t184 * t279 + (-t184 * MDP(16) - t234 * MDP(23) + t235 * MDP(24)) * t221) * qJD(5)) * qJD(2) + ((t200 * t278 + (-MDP(4) + MDP(6)) * t228) * t225 + ((-MDP(3) - MDP(5)) * t228 + (t191 + (t192 - t257) * qJD(2)) * MDP(7)) * t222) * t216; 0.2e1 * qJD(2) * qJD(3) * MDP(6) + (qJ(3) * t191 + qJD(3) * t200 + (-t200 * t225 + (-t192 - t302) * t222) * t277) * MDP(7) + t307 * MDP(8) + (t245 * qJD(2) + t173) * MDP(9) + (-t236 * qJD(3) + t170 * t178 - t171 * t181 - t172 * t248 + t173 * t283) * MDP(10) - t261 * t303 - t227 * t224 * MDP(13) + t305 * t304 + t230 * t279 + (-t176 * t288 + (t224 * t264 - t254) * t194) * MDP(18) + (-t193 * t223 - t194 * t220) * t272 * MDP(19) + (t297 + (-t259 + t292) * qJD(5) + t232) * MDP(20) + (t243 + t296 + (-t294 + (-qJD(2) * t213 + t206 * t224) * t220) * qJD(5)) * MDP(21) + (-t206 - t275) * qJD(5) * t265 + ((-t182 * t271 - t285 * t223) * t206 + (-t193 * t273 + t229 * t220 - t250 * t270 + t166) * t224) * MDP(23) + ((-t182 * t270 + t285 * t220) * t206 + (-t194 * t273 + (t250 * qJD(6) - t167) * t220 + t229 * t223) * t224) * MDP(24) + (MDP(11) * t244 + t227 * MDP(14) + t230 * MDP(16) - t305 * MDP(17) + (t298 - t177 * t223 + (t193 * t220 - t194 * t223) * qJD(6)) * MDP(19) + (-t223 * t268 - t300 - t196 * t177 - t245 * t193 + (t196 * t290 - (t182 * t223 - t196 * t289) * qJD(2) + t238) * qJD(5)) * MDP(23) + (t220 * t268 - t299 + t196 * t176 - t245 * t194 + (t196 * t287 + (t182 * t220 + t196 * t286) * qJD(2) + t156) * qJD(5)) * MDP(24)) * t221; (-t200 + t258) * t278 + (t236 * qJD(2) - t172 * t217 + t173 * t215) * MDP(10) + (t217 * t221 * t261 - t224 * t306) * MDP(16) + (t217 * t244 + t221 * t306) * MDP(17) + (t217 * t255 + ((t221 * t266 - t223 * t269) * t206 - t193 * t272 - t221 * t177) * t215 + (-(t215 * t223 - t217 * t289) * t206 + (-(-t215 * t289 - t217 * t223) * qJD(5) + t217 * t193) * t221) * qJD(2)) * MDP(23) + (t217 * t252 + (-(-t220 * t269 - t221 * t264) * t206 - t194 * t272 + t221 * t176) * t215 + ((t215 * t220 + t217 * t286) * t206 + ((t215 * t286 - t217 * t220) * qJD(5) + t217 * t194) * t221) * qJD(2)) * MDP(24) + (-t215 * MDP(8) - t217 * MDP(9) - MDP(6)) * t228; (-t243 + t296) * MDP(23) + (t232 - t297) * MDP(24) - t308 * t227 + ((t213 * t267 - t260 - t294) * MDP(23) + (-t259 - t292) * MDP(24)) * qJD(5); (-t194 * t287 + t298) * MDP(18) + ((t176 + t295) * t223 + (t177 + t293) * t220) * MDP(19) + (t252 + (t259 + (-t194 - t266) * t221) * qJD(2)) * MDP(20) + (-t255 + (-t260 + (t193 - t264) * t221) * qJD(2)) * MDP(21) + t206 * qJD(2) * t265 + (pkin(5) * t177 - t299 - (t199 * t223 + t220 * t237) * t206 + t164 * t193 + (-pkin(9) * t287 + t161 * t220) * qJD(6) + (-t238 * t221 + (t161 * t224 + t221 * t301) * t220) * qJD(2)) * MDP(23) + (-pkin(5) * t176 + t300 + (t199 * t220 - t223 * t237) * t206 + t164 * t194 + (pkin(9) * t290 + t161 * t223) * qJD(6) + (t161 * t286 + (pkin(9) * t264 - t156) * t221) * qJD(2)) * MDP(24) + (-MDP(11) * t221 * t224 + t303) * t228 + t308 * (qJD(2) * t168 - t173); t194 * t193 * MDP(18) + (-t193 ^ 2 + t194 ^ 2) * MDP(19) + (t284 - t295) * MDP(20) + (-t207 - t293) * MDP(21) + (t156 * t206 - t157 * t220 + t161 * t194 + t166) * MDP(23) + (-t157 * t223 - t161 * t193 - t167 * t220 - t206 * t238) * MDP(24) + (-t156 * MDP(23) + t238 * MDP(24)) * qJD(6) + (MDP(21) * t253 + (-t265 + (-t223 * MDP(20) + t220 * MDP(21)) * t224) * qJD(5)) * qJD(2);];
tauc  = t1;
