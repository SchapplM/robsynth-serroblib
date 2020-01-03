% Calculate minimal parameter regressor of Coriolis joint torque vector for
% S5RPRRR6
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% qJD [5x1]
%   Generalized joint velocities
% pkin [9x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d1,d3,d4,d5,theta2]';
% MDP [25x1]
%   Minimal dynamic parameter vector (fixed base model)
%   see S5RPRRR6_convert_par2_MPV_fixb.m
% 
% Output:
% tauc [5x1]
%   joint torques required to compensate Coriolis and centrifugal load

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-31 19:02
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauc = S5RPRRR6_coriolisvecJ_fixb_mdp_slag_vp(qJ, qJD, pkin, MDP)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(9,1),zeros(25,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RPRRR6_coriolisvecJ_fixb_mdp_slag_vp: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RPRRR6_coriolisvecJ_fixb_mdp_slag_vp: qJD has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S5RPRRR6_coriolisvecJ_fixb_mdp_slag_vp: pkin has to be [9x1] (double)');
assert(isreal(MDP) && all(size(MDP) == [25 1]), ...
  'S5RPRRR6_coriolisvecJ_fixb_mdp_slag_vp: MDP has to be [25x1] (double)'); 

%% Symbolic Calculation
% From coriolisvec_joint_fixb_mdp_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 19:01:44
% EndTime: 2019-12-31 19:01:49
% DurationCPUTime: 2.17s
% Computational Cost: add. (1695->224), mult. (3964->321), div. (0->0), fcn. (2734->8), ass. (0->112)
t238 = sin(qJ(4));
t241 = cos(qJ(3));
t307 = cos(qJ(4));
t271 = qJD(1) * t307;
t239 = sin(qJ(3));
t285 = qJD(1) * t239;
t313 = -t238 * t285 + t241 * t271;
t226 = sin(pkin(9)) * pkin(1) + pkin(6);
t306 = pkin(7) + t226;
t232 = qJD(3) + qJD(4);
t312 = MDP(5) * t241;
t311 = (t239 ^ 2 - t241 ^ 2) * MDP(6);
t267 = t306 * qJD(1);
t203 = t241 * qJD(2) - t267 * t239;
t310 = -t239 * MDP(10) - t241 * MDP(11);
t292 = t238 * t241;
t215 = t307 * t239 + t292;
t309 = qJD(1) * t215;
t204 = qJD(2) * t239 + t267 * t241;
t197 = t203 * qJD(3);
t198 = t204 * qJD(3);
t305 = qJD(3) * pkin(3);
t199 = t203 + t305;
t270 = qJD(4) * t307;
t283 = qJD(4) * t238;
t160 = t307 * t197 - t238 * t198 + t199 * t270 - t204 * t283;
t273 = t307 * t204;
t177 = t238 * t199 + t273;
t266 = t238 * t197 + t307 * t198;
t161 = t177 * qJD(4) + t266;
t268 = qJD(3) * t306;
t205 = t239 * t268;
t206 = t241 * t268;
t212 = t306 * t239;
t213 = t306 * t241;
t250 = -t307 * t212 - t238 * t213;
t167 = t250 * qJD(4) - t307 * t205 - t238 * t206;
t293 = t238 * t204;
t176 = t307 * t199 - t293;
t174 = -t232 * pkin(4) - t176;
t210 = -qJD(1) * t292 - t239 * t271;
t227 = -cos(pkin(9)) * pkin(1) - pkin(2);
t216 = -pkin(3) * t241 + t227;
t211 = t216 * qJD(1);
t184 = -pkin(4) * t313 + pkin(8) * t210 + t211;
t249 = -t238 * t239 + t307 * t241;
t187 = -pkin(4) * t249 - pkin(8) * t215 + t216;
t189 = -t238 * t212 + t307 * t213;
t195 = t232 * t215;
t191 = t195 * qJD(1);
t194 = t232 * t249;
t207 = qJD(5) - t313;
t308 = t161 * t215 + t174 * t194 - t189 * t191 - (qJD(5) * t187 + t167) * t207 + (qJD(5) * t184 + t160) * t249;
t237 = sin(qJ(5));
t282 = qJD(5) * t237;
t190 = t313 * t232;
t240 = cos(qJ(5));
t281 = qJD(5) * t240;
t288 = t240 * t190 + t232 * t281;
t170 = t210 * t282 + t288;
t304 = t170 * t237;
t303 = t174 * t313;
t302 = t174 * t215;
t301 = t187 * t191;
t300 = t190 * t237;
t299 = t191 * t237;
t298 = t191 * t240;
t295 = t210 * t237;
t200 = -t240 * t232 - t295;
t297 = t200 * t207;
t255 = t210 * t240 - t232 * t237;
t296 = t255 * t207;
t294 = t215 * t240;
t242 = qJD(3) ^ 2;
t291 = t239 * t242;
t290 = t241 * t242;
t289 = -t170 * t249 - t195 * t255;
t218 = qJD(1) * t227;
t278 = qJD(1) * qJD(3);
t277 = pkin(3) * t285;
t276 = t239 * t305;
t275 = t215 * t299;
t274 = t191 * t294;
t269 = t239 * t278;
t265 = t207 * t240;
t192 = -pkin(4) * t210 - pkin(8) * t313;
t229 = pkin(3) * t238 + pkin(8);
t262 = qJD(5) * t229 + t192 + t277;
t181 = t238 * t203 + t273;
t261 = pkin(3) * t283 - t181;
t175 = t232 * pkin(8) + t177;
t163 = t175 * t240 + t184 * t237;
t260 = t161 * t237 - t163 * t210 + t174 * t281;
t171 = -t255 * qJD(5) + t300;
t258 = t171 * t249 - t195 * t200;
t257 = -t191 * t229 - t303;
t256 = t175 * t237 - t184 * t240;
t254 = 0.2e1 * qJD(3) * t218;
t182 = t307 * t203 - t293;
t253 = -pkin(3) * t270 + t182;
t252 = -t161 * t240 + t174 * t282 - t210 * t256;
t251 = t211 * t210 - t266;
t248 = -t194 * t237 - t215 * t281;
t247 = -t194 * t240 + t215 * t282;
t245 = ((t170 - t297) * t240 + (-t171 + t296) * t237) * MDP(20) + (-t255 * t265 + t304) * MDP(19) + (-t207 ^ 2 * t237 - t200 * t210 + t298) * MDP(22) + (t207 * t265 - t210 * t255 + t299) * MDP(21) + t190 * MDP(14) + (t210 ^ 2 - t313 ^ 2) * MDP(13) + (MDP(12) * t313 + t207 * MDP(23)) * t210 + (-t313 * MDP(14) + (-t210 - t309) * MDP(15)) * t232;
t244 = -t211 * t313 - t160;
t230 = -t307 * pkin(3) - pkin(4);
t169 = pkin(4) * t195 - pkin(8) * t194 + t276;
t168 = t189 * qJD(4) - t238 * t205 + t307 * t206;
t165 = pkin(3) * t269 + pkin(4) * t191 - pkin(8) * t190;
t164 = t240 * t165;
t1 = [0.2e1 * t269 * t312 - 0.2e1 * t278 * t311 + MDP(7) * t290 - MDP(8) * t291 + (-t226 * t290 + t239 * t254) * MDP(10) + (t226 * t291 + t241 * t254) * MDP(11) + (t190 * t215 - t210 * t194) * MDP(12) + (t190 * t249 - t191 * t215 + t194 * t313 + t195 * t210) * MDP(13) + (t191 * t216 + t195 * t211 + (-qJD(1) * t249 - t313) * t276) * MDP(17) + (t190 * t216 + t194 * t211 + (-t210 + t309) * t276) * MDP(18) + (t170 * t294 + t247 * t255) * MDP(19) + ((-t200 * t240 + t237 * t255) * t194 + (-t304 - t171 * t240 + (t200 * t237 + t240 * t255) * qJD(5)) * t215) * MDP(20) + (-t247 * t207 + t274 + t289) * MDP(21) + (t248 * t207 + t258 - t275) * MDP(22) + (-t191 * t249 + t195 * t207) * MDP(23) + (-t256 * t195 - t164 * t249 + t168 * t200 - t250 * t171 + (t169 * t207 + t301 + (t175 * t249 - t189 * t207 + t302) * qJD(5)) * t240 + t308 * t237) * MDP(24) + (-t163 * t195 - t168 * t255 - t250 * t170 + (-(-qJD(5) * t189 + t169) * t207 - t301 + (-qJD(5) * t175 + t165) * t249 - qJD(5) * t302) * t237 + t308 * t240) * MDP(25) + (t194 * MDP(14) - t195 * MDP(15) - t168 * MDP(17) - t167 * MDP(18)) * t232; (-t258 - t275) * MDP(24) + (-t274 + t289) * MDP(25) + t310 * t242 + (-t195 * MDP(17) - t194 * MDP(18)) * t232 + (t248 * MDP(24) + t247 * MDP(25)) * t207; t245 + (t230 * t171 + t257 * t237 + t261 * t200 + (t253 * t237 - t262 * t240) * t207 + t252) * MDP(24) + (t230 * t170 + t257 * t240 - t261 * t255 + (t262 * t237 + t253 * t240) * t207 + t260) * MDP(25) + (t313 * t277 + t181 * t232 + (-t273 + (-pkin(3) * t232 - t199) * t238) * qJD(4) + t251) * MDP(17) + (t182 * t232 + (t210 * t285 - t232 * t270) * pkin(3) + t244) * MDP(18) + (t310 * t218 + (-t239 * t312 + t311) * qJD(1)) * qJD(1); (t251 + (-qJD(4) + t232) * t177) * MDP(17) + (t176 * t232 + t244) * MDP(18) + (-pkin(4) * t171 - (-t176 * t237 + t192 * t240) * t207 - t177 * t200 - t237 * t303 + (-t207 * t281 - t299) * pkin(8) + t252) * MDP(24) + (-pkin(4) * t170 + (t176 * t240 + t192 * t237) * t207 + t177 * t255 - t240 * t303 + (t207 * t282 - t298) * pkin(8) + t260) * MDP(25) + t245; -t255 * t200 * MDP(19) + (-t200 ^ 2 + t255 ^ 2) * MDP(20) + (t288 + t297) * MDP(21) + (-t296 - t300) * MDP(22) + t191 * MDP(23) + (-t160 * t237 + t163 * t207 + t174 * t255 + t164) * MDP(24) + (-t160 * t240 - t165 * t237 + t174 * t200 - t207 * t256) * MDP(25) + (MDP(21) * t295 + t255 * MDP(22) - t163 * MDP(24) + t256 * MDP(25)) * qJD(5);];
tauc = t1;
