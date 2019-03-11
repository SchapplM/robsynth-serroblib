% Calculate minimal parameter regressor of Coriolis joint torque vector for
% S6RPPRPR3
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% qJD [6x1]
%   Generalized joint velocities
% pkin [10x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d4,d6,theta2,theta5]';
% MDP [23x1]
%   Minimal dynamic parameter vector (fixed base model)
%   see S6RPPRPR3_convert_par2_MPV_fixb.m
% 
% Output:
% tauc [6x1]
%   joint torques required to compensate Coriolis and centrifugal load

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-09 01:45
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauc = S6RPPRPR3_coriolisvecJ_fixb_mdp_slag_vp(qJ, qJD, pkin, MDP)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(10,1),zeros(23,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPPRPR3_coriolisvecJ_fixb_mdp_slag_vp: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RPPRPR3_coriolisvecJ_fixb_mdp_slag_vp: qJD has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S6RPPRPR3_coriolisvecJ_fixb_mdp_slag_vp: pkin has to be [10x1] (double)');
assert(isreal(MDP) && all(size(MDP) == [23 1]), ...
  'S6RPPRPR3_coriolisvecJ_fixb_mdp_slag_vp: MDP has to be [23x1] (double)'); 

%% Symbolic Calculation
% From coriolisvec_joint_fixb_mdp_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 01:45:01
% EndTime: 2019-03-09 01:45:04
% DurationCPUTime: 1.92s
% Computational Cost: add. (1629->236), mult. (3393->330), div. (0->0), fcn. (2233->8), ass. (0->112)
t236 = cos(qJ(6));
t270 = t236 * qJD(4);
t231 = sin(pkin(10));
t237 = cos(qJ(4));
t299 = cos(pkin(10));
t258 = t299 * t237;
t235 = sin(qJ(4));
t277 = qJD(1) * t235;
t203 = qJD(1) * t258 - t231 * t277;
t234 = sin(qJ(6));
t287 = t203 * t234;
t189 = -t270 + t287;
t245 = t231 * t237 + t235 * t299;
t201 = t245 * qJD(1);
t308 = qJD(6) + t201;
t312 = t189 * t308;
t191 = qJD(4) * t234 + t203 * t236;
t311 = t191 * t308;
t216 = sin(pkin(9)) * pkin(1) + qJ(3);
t210 = qJD(1) * t216;
t310 = t210 * MDP(7);
t300 = 2 * qJD(3);
t256 = t308 * t236;
t196 = t203 * qJD(4);
t289 = t196 * t234;
t309 = -t256 * t308 - t289;
t307 = MDP(8) * t235;
t306 = MDP(9) * (t235 ^ 2 - t237 ^ 2);
t265 = qJD(1) * qJD(5);
t214 = -cos(pkin(9)) * pkin(1) - pkin(2) - pkin(7);
t207 = qJD(1) * t214 + qJD(3);
t255 = qJ(5) * qJD(1) - t207;
t302 = -t237 * qJD(2) + t235 * t255;
t304 = t302 * qJD(4) - t237 * t265;
t275 = qJD(2) * t235;
t180 = -t235 * t265 + (-t237 * t255 - t275) * qJD(4);
t157 = t180 * t231 - t299 * t304;
t220 = pkin(4) * t231 + pkin(8);
t276 = qJD(1) * t237;
t303 = t308 * (pkin(4) * t276 + pkin(5) * t203 + pkin(8) * t201 + qJD(6) * t220) + t157;
t158 = t299 * t180 + t304 * t231;
t186 = -qJ(5) * t276 + t237 * t207 - t275;
t185 = qJD(4) * pkin(4) + t186;
t286 = t231 * t302;
t161 = t185 * t299 + t286;
t159 = -qJD(4) * pkin(5) - t161;
t283 = qJ(5) - t214;
t257 = t283 * t237;
t192 = -qJD(4) * t257 - qJD(5) * t235;
t274 = qJD(4) * t235;
t243 = -qJD(5) * t237 + t274 * t283;
t168 = t192 * t299 + t231 * t243;
t200 = pkin(4) * t277 + qJD(5) + t210;
t169 = pkin(5) * t201 - pkin(8) * t203 + t200;
t209 = -t231 * t235 + t258;
t262 = t235 * pkin(4) + t216;
t177 = pkin(5) * t245 - pkin(8) * t209 + t262;
t205 = t245 * qJD(4);
t206 = t283 * t235;
t179 = -t206 * t299 - t231 * t257;
t293 = t179 * t196;
t298 = t157 * t209;
t301 = -t159 * t205 - (qJD(6) * t177 + t168) * t308 - (qJD(6) * t169 + t158) * t245 - t293 + t298;
t297 = t159 * t209;
t272 = qJD(6) * t234;
t197 = qJD(1) * t205;
t280 = qJD(6) * t270 - t236 * t197;
t171 = -t203 * t272 + t280;
t296 = t171 * t209;
t295 = t171 * t234;
t294 = t177 * t196;
t292 = t189 * t203;
t291 = t191 * t203;
t290 = t196 * t245;
t288 = t197 * t234;
t238 = qJD(4) ^ 2;
t285 = t235 * t238;
t193 = t236 * t196;
t284 = t237 * t238;
t202 = -qJD(4) * t258 + t231 * t274;
t282 = t171 * t245 - t191 * t202;
t183 = t299 * t302;
t162 = t231 * t185 - t183;
t228 = qJD(3) * qJD(1);
t266 = qJD(1) * qJD(4);
t260 = t237 * t266;
t281 = pkin(4) * t260 + t228;
t273 = qJD(4) * t237;
t271 = qJD(6) * t236;
t268 = pkin(4) * t273 + qJD(3);
t264 = t209 * t289;
t263 = t209 * t193;
t160 = qJD(4) * pkin(8) + t162;
t156 = t160 * t236 + t169 * t234;
t250 = t160 * t234 - t169 * t236;
t172 = qJD(6) * t191 - t288;
t249 = -t172 * t245 + t189 * t202;
t248 = t193 + (-t201 * t234 - t272) * t308;
t247 = t205 * t234 - t209 * t271;
t246 = t205 * t236 + t209 * t272;
t164 = t186 * t299 + t286;
t242 = -t220 * t196 + (t159 + t164) * t308;
t240 = t158 * t245 - t161 * t205 - t162 * t202 - t298;
t239 = qJD(1) ^ 2;
t221 = -pkin(4) * t299 - pkin(5);
t178 = -t206 * t231 + t257 * t299;
t173 = -pkin(5) * t202 + pkin(8) * t205 + t268;
t170 = pkin(5) * t196 + pkin(8) * t197 + t281;
t167 = t192 * t231 - t243 * t299;
t166 = t236 * t170;
t163 = t186 * t231 - t183;
t1 = [0.2e1 * MDP(6) * t228 + t300 * t310 - 0.2e1 * t260 * t307 + 0.2e1 * t266 * t306 - MDP(10) * t285 - MDP(11) * t284 + (t210 * t273 - t214 * t285 + (t216 * t273 + t235 * t300) * qJD(1)) * MDP(13) + (-t210 * t274 - t214 * t284 + (-t216 * t274 + t237 * t300) * qJD(1)) * MDP(14) + (t167 * t203 - t168 * t201 - t178 * t197 - t240 - t293) * MDP(15) + (t157 * t178 + t158 * t179 - t161 * t167 + t162 * t168 + t200 * t268 + t262 * t281) * MDP(16) + (-t191 * t246 + t236 * t296) * MDP(17) + ((t189 * t236 + t191 * t234) * t205 + (-t295 - t172 * t236 + (t189 * t234 - t191 * t236) * qJD(6)) * t209) * MDP(18) + (-t246 * t308 + t263 + t282) * MDP(19) + (t247 * t308 + t249 - t264) * MDP(20) + (-t202 * t308 + t290) * MDP(21) + (t250 * t202 + t166 * t245 + t167 * t189 + t178 * t172 + (t173 * t308 + t294 + (-t160 * t245 - t179 * t308 + t297) * qJD(6)) * t236 + t301 * t234) * MDP(22) + (t156 * t202 + t167 * t191 + t178 * t171 + (-(-qJD(6) * t179 + t173) * t308 - t294 - (-qJD(6) * t160 + t170) * t245 - qJD(6) * t297) * t234 + t301 * t236) * MDP(23); (-t196 * t209 - t197 * t245 + t201 * t205 - t202 * t203) * MDP(15) + (t157 * t245 + t158 * t209 + t161 * t202 - t162 * t205) * MDP(16) + (-t249 - t264) * MDP(22) + (-t263 + t282) * MDP(23) + (-MDP(13) * t237 + MDP(14) * t235) * t238 + (MDP(22) * t247 + MDP(23) * t246) * t308; -t239 * MDP(6) - qJD(1) * t310 + (t197 * t209 + t201 * t202 + t203 * t205 - t290) * MDP(15) + (-qJD(1) * t200 + t240) * MDP(16) + (-t172 * t209 + t189 * t205 - t245 * t289) * MDP(22) + (t191 * t205 - t193 * t245 - t296) * MDP(23) + ((-qJD(1) * t236 + t202 * t234 - t245 * t271) * MDP(22) + (qJD(1) * t234 + t202 * t236 + t245 * t272) * MDP(23)) * t308 + (MDP(13) * t235 + MDP(14) * t237) * (-t238 - t239); ((t162 - t163) * t203 - (t161 - t164) * t201 + (-t196 * t231 + t197 * t299) * pkin(4)) * MDP(15) + (t161 * t163 - t162 * t164 + (-t157 * t299 + t158 * t231 - t200 * t276) * pkin(4)) * MDP(16) + (t191 * t256 + t295) * MDP(17) + ((t171 - t312) * t236 + (-t172 - t311) * t234) * MDP(18) + (-t291 - t309) * MDP(19) + (t248 + t292) * MDP(20) - t308 * t203 * MDP(21) + (-t163 * t189 + t221 * t172 + t203 * t250 + t242 * t234 - t303 * t236) * MDP(22) + (t156 * t203 - t163 * t191 + t221 * t171 + t303 * t234 + t242 * t236) * MDP(23) + (t237 * t307 - t306) * t239 + (-MDP(13) * t276 + MDP(14) * t277) * t210; (-t201 ^ 2 - t203 ^ 2) * MDP(15) + (t161 * t203 + t162 * t201 + t281) * MDP(16) + (t248 - t292) * MDP(22) + (-t291 + t309) * MDP(23); t191 * t189 * MDP(17) + (-t189 ^ 2 + t191 ^ 2) * MDP(18) + (t280 + t312) * MDP(19) + (t288 + t311) * MDP(20) + t196 * MDP(21) + (t156 * t308 - t158 * t234 - t159 * t191 + t166) * MDP(22) + (-t158 * t236 + t159 * t189 - t170 * t234 - t250 * t308) * MDP(23) + (-MDP(19) * t287 - MDP(20) * t191 - MDP(22) * t156 + MDP(23) * t250) * qJD(6);];
tauc  = t1;
