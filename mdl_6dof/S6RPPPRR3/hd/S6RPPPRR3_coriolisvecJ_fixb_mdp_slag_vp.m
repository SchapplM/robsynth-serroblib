% Calculate minimal parameter regressor of Coriolis joint torque vector for
% S6RPPPRR3
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% qJD [6x1]
%   Generalized joint velocities
% pkin [10x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d5,d6,theta3,theta4]';
% MDP [27x1]
%   Minimal dynamic parameter vector (fixed base model)
%   see S6RPPPRR3_convert_par2_MPV_fixb.m
% 
% Output:
% tauc [6x1]
%   joint torques required to compensate Coriolis and centrifugal load

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-09 01:34
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauc = S6RPPPRR3_coriolisvecJ_fixb_mdp_slag_vp(qJ, qJD, pkin, MDP)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(10,1),zeros(27,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPPPRR3_coriolisvecJ_fixb_mdp_slag_vp: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RPPPRR3_coriolisvecJ_fixb_mdp_slag_vp: qJD has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S6RPPPRR3_coriolisvecJ_fixb_mdp_slag_vp: pkin has to be [10x1] (double)');
assert(isreal(MDP) && all(size(MDP) == [27 1]), ...
  'S6RPPPRR3_coriolisvecJ_fixb_mdp_slag_vp: MDP has to be [27x1] (double)'); 

%% Symbolic Calculation
% From coriolisvec_joint_fixb_mdp_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 01:34:02
% EndTime: 2019-03-09 01:34:06
% DurationCPUTime: 1.97s
% Computational Cost: add. (1508->234), mult. (3216->331), div. (0->0), fcn. (2244->8), ass. (0->115)
t234 = cos(qJ(6));
t265 = t234 * qJD(5);
t227 = sin(pkin(10));
t229 = cos(pkin(10));
t233 = sin(qJ(5));
t235 = cos(qJ(5));
t203 = t227 * t235 + t229 * t233;
t199 = qJD(1) * t203;
t232 = sin(qJ(6));
t285 = t199 * t232;
t177 = -t265 - t285;
t270 = qJD(1) * t235;
t209 = t229 * t270;
t271 = qJD(1) * t233;
t198 = t227 * t271 - t209;
t264 = qJD(6) - t198;
t305 = t177 * t264;
t179 = qJD(5) * t232 - t199 * t234;
t304 = t179 * t264;
t303 = t229 * MDP(10) - t227 * MDP(11) + MDP(7);
t255 = t264 * t234;
t201 = t203 * qJD(5);
t193 = qJD(1) * t201;
t281 = t232 * t193;
t302 = t264 * t255 - t281;
t280 = t235 * t229;
t284 = t227 * t233;
t202 = -t280 + t284;
t230 = cos(pkin(9));
t211 = qJD(2) * t230 - qJD(4);
t206 = t211 * qJD(1);
t300 = t202 * t206;
t299 = t203 * t206;
t275 = t227 ^ 2 + t229 ^ 2;
t257 = t206 * t275;
t237 = qJD(1) ^ 2;
t298 = t275 * MDP(12) * t237;
t296 = qJ(2) * MDP(6) + t230 * MDP(8) + MDP(5);
t236 = -pkin(1) - pkin(2);
t208 = qJD(1) * t236 + qJD(2);
t228 = sin(pkin(9));
t272 = qJD(1) * t230;
t195 = qJ(2) * t272 + t228 * t208;
t183 = -qJD(1) * qJ(4) + t195;
t220 = t229 * qJD(3);
t169 = t220 + (pkin(7) * qJD(1) - t183) * t227;
t174 = t227 * qJD(3) + t229 * t183;
t273 = qJD(1) * t229;
t170 = -pkin(7) * t273 + t174;
t157 = t169 * t233 + t170 * t235;
t153 = qJD(5) * t157 + t299;
t295 = t264 * (-pkin(5) * t199 + pkin(8) * t264) + t153;
t156 = t169 * t235 - t170 * t233;
t152 = qJD(5) * t156 - t300;
t154 = -qJD(5) * pkin(5) - t156;
t276 = t230 * qJ(2) + t228 * t236;
t204 = -qJ(4) + t276;
t292 = pkin(7) - t204;
t188 = t292 * t227;
t189 = t292 * t229;
t245 = t188 * t235 + t189 * t233;
t158 = qJD(5) * t245 - t202 * t211;
t266 = t228 * qJD(1);
t194 = -qJ(2) * t266 + t208 * t230;
t180 = qJD(1) * pkin(3) + qJD(4) - t194;
t175 = pkin(4) * t273 + t180;
t160 = -pkin(5) * t198 + pkin(8) * t199 + t175;
t163 = t188 * t233 - t189 * t235;
t256 = -t228 * qJ(2) + t230 * t236;
t251 = pkin(3) - t256;
t197 = t229 * pkin(4) + t251;
t167 = -pkin(5) * t202 + pkin(8) * t203 + t197;
t259 = qJD(5) * t284;
t200 = -qJD(5) * t280 + t259;
t293 = -t153 * t203 + t154 * t200 + t163 * t193 - (qJD(6) * t167 + t158) * t264 + (qJD(6) * t160 + t152) * t202;
t267 = qJD(6) * t232;
t207 = qJD(1) * t259;
t192 = -qJD(5) * t209 + t207;
t279 = qJD(6) * t265 + t234 * t192;
t164 = t199 * t267 + t279;
t290 = t164 * t232;
t289 = t167 * t193;
t288 = t177 * t199;
t287 = t179 * t199;
t286 = t192 * t232;
t182 = t234 * t193;
t278 = t228 * t201 - t202 * t272;
t191 = t202 * t228;
t277 = -qJD(5) * t191 - t230 * t199;
t269 = qJD(2) * t228;
t268 = qJD(6) * t203;
t263 = qJD(1) * qJD(2);
t261 = t203 * t281;
t260 = t203 * t182;
t212 = t228 * t263;
t250 = -qJD(6) * t230 - t278;
t155 = qJD(5) * pkin(8) + t157;
t151 = t155 * t234 + t160 * t232;
t249 = t155 * t232 - t160 * t234;
t248 = -t164 * t202 - t179 * t201;
t165 = t179 * qJD(6) + t286;
t247 = t165 * t202 + t177 * t201;
t246 = -(-t183 * t227 + t220) * t227 + t174 * t229;
t244 = t194 * t228 - t195 * t230;
t243 = -qJD(6) * t191 + t266;
t242 = -t182 + (t198 * t232 - t267) * t264;
t241 = -t232 * t200 + t234 * t268;
t240 = t234 * t200 + t203 * t267;
t238 = pkin(8) * t193 + (t154 + t156) * t264;
t190 = t203 * t228;
t168 = -pkin(5) * t201 - pkin(8) * t200 + t269;
t166 = -pkin(5) * t193 - pkin(8) * t192 + t212;
t161 = t234 * t166;
t159 = qJD(5) * t163 + t203 * t211;
t1 = [((-t228 * t256 + t230 * t276) * qJD(1) - t244) * qJD(2) * MDP(9) + (t246 * t211 + t204 * t257 + (qJD(1) * t251 + t180) * t269) * MDP(13) + (-t192 * t203 - t199 * t200) * MDP(14) + (t192 * t202 - t193 * t203 + t198 * t200 - t199 * t201) * MDP(15) + (-t175 * t201 - t193 * t197 + (-qJD(1) * t202 - t198) * t269) * MDP(19) + (t175 * t200 + t192 * t197 - 0.2e1 * t199 * t269) * MDP(20) + (-t164 * t203 * t234 + t179 * t240) * MDP(21) + ((-t177 * t234 - t179 * t232) * t200 + (t290 + t165 * t234 + (-t177 * t232 + t179 * t234) * qJD(6)) * t203) * MDP(22) + (t240 * t264 + t248 + t260) * MDP(23) + (t241 * t264 + t247 - t261) * MDP(24) + (t193 * t202 - t201 * t264) * MDP(25) + (t249 * t201 + t159 * t177 - t161 * t202 - t245 * t165 + (-t289 + t168 * t264 + (-t154 * t203 + t155 * t202 - t163 * t264) * qJD(6)) * t234 + t293 * t232) * MDP(26) + (t151 * t201 + t159 * t179 - t245 * t164 + (-(-qJD(6) * t163 + t168) * t264 + t289 + (-qJD(6) * t155 + t166) * t202 + t154 * t268) * t232 + t293 * t234) * MDP(27) - 0.2e1 * MDP(12) * t257 + (MDP(16) * t200 + MDP(17) * t201 - MDP(19) * t159 - MDP(20) * t158) * qJD(5) + 0.2e1 * t296 * t263 + 0.2e1 * t303 * t212; t244 * qJD(1) * MDP(9) + t230 * t298 + (t228 * t257 + (-t180 * t228 + (-t246 - t269) * t230) * qJD(1)) * MDP(13) + (-qJD(5) * t277 + t193 * t230 + t198 * t266) * MDP(19) + (qJD(5) * t278 - t192 * t230 + t199 * t266) * MDP(20) + (-(t191 * t232 - t230 * t234) * t193 + t190 * t165 - (t232 * t250 + t234 * t243) * t264 + t277 * t177) * MDP(26) + ((-t191 * t234 - t230 * t232) * t193 + t190 * t164 - (-t232 * t243 + t234 * t250) * t264 + t277 * t179) * MDP(27) + (-t228 * t303 - t296) * t237; (t247 + t261) * MDP(26) + (-t248 + t260) * MDP(27) - (MDP(26) * t241 - MDP(27) * t240) * t264 + (-t201 * MDP(19) + t200 * MDP(20)) * qJD(5); (qJD(1) * t246 + t212) * MDP(13) + t207 * MDP(20) + (t242 + t288) * MDP(26) + (t287 - t302) * MDP(27) - t298 + ((-t227 * t270 - t229 * t271 - t199) * MDP(19) + (t198 - t209) * MDP(20)) * qJD(5); -t198 ^ 2 * MDP(15) + (t207 + (-t198 - t209) * qJD(5)) * MDP(16) - MDP(19) * t299 + (-t175 * t198 + t300) * MDP(20) + (t179 * t255 + t290) * MDP(21) + ((t164 - t305) * t234 + (-t165 - t304) * t232) * MDP(22) + (t287 + t302) * MDP(23) + (t242 - t288) * MDP(24) + (-pkin(5) * t165 - t157 * t177 + t238 * t232 - t295 * t234) * MDP(26) + (-pkin(5) * t164 - t157 * t179 + t295 * t232 + t238 * t234) * MDP(27) + (MDP(14) * t198 + t199 * MDP(15) + t175 * MDP(19) + MDP(25) * t264 - MDP(26) * t249 - t151 * MDP(27)) * t199; t179 * t177 * MDP(21) + (-t177 ^ 2 + t179 ^ 2) * MDP(22) + (t279 + t305) * MDP(23) + (-t286 + t304) * MDP(24) - t193 * MDP(25) + (t151 * t264 - t152 * t232 - t154 * t179 + t161) * MDP(26) + (-t152 * t234 + t154 * t177 - t166 * t232 - t249 * t264) * MDP(27) + (MDP(23) * t285 - MDP(24) * t179 - MDP(26) * t151 + MDP(27) * t249) * qJD(6);];
tauc  = t1;
