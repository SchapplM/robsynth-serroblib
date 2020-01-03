% Calculate minimal parameter regressor of Coriolis joint torque vector for
% S5RPPRR4
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% qJD [5x1]
%   Generalized joint velocities
% pkin [9x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d1,d4,d5,theta2,theta3]';
% MDP [25x1]
%   Minimal dynamic parameter vector (fixed base model)
%   see S5RPPRR4_convert_par2_MPV_fixb.m
% 
% Output:
% tauc [5x1]
%   joint torques required to compensate Coriolis and centrifugal load

% Quelle: HybrDyn-Toolbox
% Datum: 2020-01-03 11:32
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauc = S5RPPRR4_coriolisvecJ_fixb_mdp_slag_vp(qJ, qJD, pkin, MDP)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(9,1),zeros(25,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RPPRR4_coriolisvecJ_fixb_mdp_slag_vp: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RPPRR4_coriolisvecJ_fixb_mdp_slag_vp: qJD has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S5RPPRR4_coriolisvecJ_fixb_mdp_slag_vp: pkin has to be [9x1] (double)');
assert(isreal(MDP) && all(size(MDP) == [25 1]), ...
  'S5RPPRR4_coriolisvecJ_fixb_mdp_slag_vp: MDP has to be [25x1] (double)'); 

%% Symbolic Calculation
% From coriolisvec_joint_fixb_mdp_matlab.m
% OptimizationMode: 2
% StartTime: 2020-01-03 11:31:30
% EndTime: 2020-01-03 11:31:39
% DurationCPUTime: 2.24s
% Computational Cost: add. (1372->219), mult. (3986->330), div. (0->0), fcn. (3046->8), ass. (0->117)
t254 = sin(pkin(9));
t256 = cos(pkin(9));
t259 = sin(qJ(4));
t261 = cos(qJ(4));
t235 = t254 * t261 + t256 * t259;
t231 = t235 * qJD(4);
t255 = sin(pkin(8));
t215 = t255 * t231;
t205 = qJD(1) * t215;
t300 = qJD(1) * t255;
t287 = t254 * t300;
t277 = t259 * t287;
t237 = qJD(4) * t277;
t306 = t256 * t261;
t289 = t255 * t306;
t278 = qJD(4) * t289;
t206 = qJD(1) * t278 - t237;
t213 = qJD(1) * t289 - t277;
t258 = sin(qJ(5));
t260 = cos(qJ(5));
t293 = qJD(5) * t258;
t266 = qJD(1) * t235;
t210 = t255 * t266;
t305 = t260 * t210;
t172 = -qJD(5) * t305 - t260 * t205 - t258 * t206 - t213 * t293;
t271 = -t210 * t258 + t260 * t213;
t173 = t271 * qJD(5) - t205 * t258 + t260 * t206;
t179 = t213 * t258 + t305;
t257 = cos(pkin(8));
t299 = qJD(1) * t257;
t320 = qJD(4) - t299;
t243 = -qJD(5) - t320;
t310 = t179 * t243;
t311 = t271 * t243;
t321 = t179 * t271 * MDP(19) + (-t173 - t311) * MDP(22) + (-t179 ^ 2 + t271 ^ 2) * MDP(20) + (t172 - t310) * MDP(21);
t238 = -pkin(2) * t257 - qJ(3) * t255 - pkin(1);
t226 = t238 * qJD(1) + qJD(2);
t218 = t256 * t226;
t265 = -pkin(6) * t255 * t256 + (-qJ(2) * t254 - pkin(3)) * t257;
t185 = t265 * qJD(1) + t218;
t288 = qJ(2) * t299;
t195 = t254 * t226 + t256 * t288;
t190 = -pkin(6) * t287 + t195;
t273 = -t185 * t259 - t190 * t261;
t296 = qJD(3) * t255;
t297 = qJD(2) * t257;
t228 = -t254 * t297 - t256 * t296;
t222 = t228 * qJD(1);
t229 = -t254 * t296 + t256 * t297;
t223 = t229 * qJD(1);
t280 = t261 * t222 - t259 * t223;
t264 = t273 * qJD(4) + t280;
t165 = pkin(7) * t205 + t264;
t171 = -pkin(7) * t210 - t273;
t169 = t171 * t293;
t244 = qJ(2) * t300 + qJD(3);
t227 = pkin(3) * t287 + t244;
t188 = pkin(4) * t210 + t227;
t319 = t188 * t179 + t169 + (t171 * t243 - t165) * t258;
t316 = qJD(5) + t243;
t294 = qJD(4) * t261;
t295 = qJD(4) * t259;
t268 = -t185 * t294 + t190 * t295 - t259 * t222 - t261 * t223;
t164 = -pkin(7) * t206 - t268;
t283 = -t258 * t164 + t260 * t165;
t315 = -t188 * t271 + t283;
t314 = pkin(4) * t213;
t313 = qJ(2) * t257;
t312 = t171 * t260;
t309 = t210 * t320;
t308 = t213 * t320;
t307 = t254 * t255;
t304 = -t257 * t266 + t231;
t234 = -t254 * t259 + t306;
t303 = t320 * t234;
t302 = t254 * t238 + t256 * t313;
t236 = pkin(3) * t307 + t255 * qJ(2);
t252 = t255 ^ 2;
t253 = t257 ^ 2;
t301 = t252 + t253;
t298 = qJD(2) * t255;
t291 = qJD(1) * qJD(2);
t290 = 0.2e1 * qJD(2) * t252;
t286 = qJ(2) * t291;
t282 = t261 * t185 - t190 * t259;
t170 = -pkin(7) * t213 + t282;
t168 = pkin(4) * t320 + t170;
t285 = pkin(4) * t243 - t168;
t279 = t261 * t228 - t229 * t259;
t276 = qJD(5) * t235 + t304;
t275 = qJD(5) * t234 + t303;
t274 = -t168 * t258 - t312;
t233 = t256 * t238;
t191 = t233 + t265;
t196 = -pkin(6) * t307 + t302;
t272 = -t191 * t259 - t196 * t261;
t270 = -t222 * t256 - t223 * t254;
t224 = t235 * t255;
t225 = t234 * t255;
t186 = t224 * t260 + t225 * t258;
t187 = -t224 * t258 + t225 * t260;
t267 = t191 * t294 - t196 * t295 + t259 * t228 + t261 * t229;
t262 = qJD(1) ^ 2;
t249 = t255 * t291;
t246 = t252 * t286;
t216 = -t295 * t307 + t278;
t197 = pkin(4) * t216 + t298;
t194 = -t254 * t288 + t218;
t193 = pkin(4) * t206 + t249;
t192 = pkin(4) * t224 + t236;
t177 = -pkin(7) * t224 - t272;
t176 = -pkin(4) * t257 - pkin(7) * t225 + t191 * t261 - t196 * t259;
t175 = t187 * qJD(5) - t215 * t258 + t260 * t216;
t174 = -t186 * qJD(5) - t215 * t260 - t216 * t258;
t167 = pkin(7) * t215 + t272 * qJD(4) + t279;
t166 = -pkin(7) * t216 + t267;
t1 = [0.2e1 * t301 * MDP(6) * t291 + 0.2e1 * (t253 * t286 + t246) * MDP(7) + (-t222 * t257 + (-t228 * t257 + t254 * t290) * qJD(1)) * MDP(8) + (t223 * t257 + (t229 * t257 + t256 * t290) * qJD(1)) * MDP(9) + ((-t228 * t256 - t229 * t254) * qJD(1) + t270) * t255 * MDP(10) + (t223 * t302 + t195 * t229 + t222 * (-t254 * t313 + t233) + t194 * t228 + t246 + t244 * t298) * MDP(11) + (-t205 * t225 - t213 * t215) * MDP(12) + (t205 * t224 - t206 * t225 + t210 * t215 - t213 * t216) * MDP(13) + (t205 * t257 - t215 * t320) * MDP(14) + (t206 * t257 - t216 * t320) * MDP(15) + (t279 * t320 - t280 * t257 + t236 * t206 + t227 * t216 + (-t257 * t273 + t272 * t320) * qJD(4) + (qJD(1) * t224 + t210) * t298) * MDP(17) + (-t267 * t320 - t268 * t257 - t236 * t205 - t227 * t215 + (qJD(1) * t225 + t213) * t298) * MDP(18) + (t172 * t187 + t174 * t271) * MDP(19) + (-t172 * t186 - t173 * t187 - t174 * t179 - t175 * t271) * MDP(20) + (-t172 * t257 - t174 * t243) * MDP(21) + (t173 * t257 + t175 * t243) * MDP(22) + (-(-t166 * t258 + t167 * t260) * t243 - t283 * t257 + t197 * t179 + t192 * t173 + t193 * t186 + t188 * t175 + (-(-t176 * t258 - t177 * t260) * t243 - t274 * t257) * qJD(5)) * MDP(24) + (-t169 * t257 + t192 * t172 + t188 * t174 + t197 * t271 + t193 * t187 + ((-qJD(5) * t177 + t167) * t243 + t165 * t257) * t258 + ((qJD(5) * t176 + t166) * t243 + (qJD(5) * t168 + t164) * t257) * t260) * MDP(25); ((-t244 * t255 + (t194 * t254 - t195 * t256) * t257) * qJD(1) - t270) * MDP(11) + (-t210 * t300 - t304 * t320) * MDP(17) + (-t213 * t300 - t303 * t320) * MDP(18) + (-t179 * t300 + (t258 * t275 + t260 * t276) * t243) * MDP(24) + (-t271 * t300 + (-t258 * t276 + t260 * t275) * t243) * MDP(25) - (MDP(7) * qJ(2) + MDP(8) * t254 + MDP(9) * t256 + MDP(6)) * t301 * t262; t249 * MDP(11) + (-t237 + t308) * MDP(17) - MDP(18) * t309 + (t173 - t311) * MDP(24) + (t172 + t310) * MDP(25) + (-t254 ^ 2 - t256 ^ 2) * t262 * MDP(10) * t252 + ((-MDP(8) * t256 + MDP(9) * t254) * t262 * t257 + ((t194 * t256 + t195 * t254) * MDP(11) + (MDP(17) * t306 - MDP(18) * t235) * qJD(4)) * qJD(1)) * t255; t213 * t210 * MDP(12) + (-t210 ^ 2 + t213 ^ 2) * MDP(13) + (-t205 + t309) * MDP(14) + (-t206 + t308) * MDP(15) + (-t227 * t213 - t273 * t320 + t264) * MDP(17) + (t227 * t210 + t282 * t320 + t268) * MDP(18) + ((-t170 * t258 - t312) * t243 - t179 * t314 + (t258 * t285 - t312) * qJD(5) + t315) * MDP(24) + (-t271 * t314 + (qJD(5) * t285 - t170 * t243 - t164) * t260 + t319) * MDP(25) + t321; (t316 * t274 + t315) * MDP(24) + ((-t316 * t168 - t164) * t260 + t319) * MDP(25) + t321;];
tauc = t1;
