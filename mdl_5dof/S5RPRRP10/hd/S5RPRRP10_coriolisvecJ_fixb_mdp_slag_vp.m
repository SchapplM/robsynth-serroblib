% Calculate minimal parameter regressor of Coriolis joint torque vector for
% S5RPRRP10
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% qJD [5x1]
%   Generalized joint velocities
% pkin [8x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d1,d3,d4,theta2]';
% MDP [23x1]
%   Minimal dynamic parameter vector (fixed base model)
%   see S5RPRRP10_convert_par2_MPV_fixb.m
% 
% Output:
% tauc [5x1]
%   joint torques required to compensate Coriolis and centrifugal load

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-31 18:52
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauc = S5RPRRP10_coriolisvecJ_fixb_mdp_slag_vp(qJ, qJD, pkin, MDP)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(8,1),zeros(23,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RPRRP10_coriolisvecJ_fixb_mdp_slag_vp: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RPRRP10_coriolisvecJ_fixb_mdp_slag_vp: qJD has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S5RPRRP10_coriolisvecJ_fixb_mdp_slag_vp: pkin has to be [8x1] (double)');
assert(isreal(MDP) && all(size(MDP) == [23 1]), ...
  'S5RPRRP10_coriolisvecJ_fixb_mdp_slag_vp: MDP has to be [23x1] (double)'); 

%% Symbolic Calculation
% From coriolisvec_joint_fixb_mdp_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 18:52:08
% EndTime: 2019-12-31 18:52:14
% DurationCPUTime: 2.34s
% Computational Cost: add. (1923->253), mult. (5020->346), div. (0->0), fcn. (3624->6), ass. (0->105)
t240 = sin(pkin(8));
t291 = pkin(6) + qJ(2);
t226 = t291 * t240;
t223 = qJD(1) * t226;
t241 = cos(pkin(8));
t227 = t291 * t241;
t224 = qJD(1) * t227;
t243 = sin(qJ(3));
t245 = cos(qJ(3));
t195 = -t243 * t223 + t245 * t224;
t301 = qJD(3) * t195;
t300 = (t240 ^ 2 + t241 ^ 2) * (qJ(2) * MDP(7) + MDP(6));
t222 = t240 * t245 + t241 * t243;
t217 = t222 * qJD(1);
t244 = cos(qJ(4));
t270 = qJD(1) * t245;
t233 = t241 * t270;
t232 = qJD(3) * t233;
t271 = qJD(1) * t243;
t263 = t240 * t271;
t254 = t263 * qJD(3) - t232;
t267 = t244 * qJD(3);
t242 = sin(qJ(4));
t269 = qJD(4) * t242;
t170 = -qJD(4) * t267 + t217 * t269 + t244 * t254;
t203 = qJD(3) * t242 + t217 * t244;
t219 = t222 * qJD(3);
t207 = qJD(1) * t219;
t216 = t233 - t263;
t235 = -pkin(2) * t241 - pkin(1);
t225 = qJD(1) * t235 + qJD(2);
t174 = -pkin(3) * t216 - pkin(7) * t217 + t225;
t190 = qJD(3) * pkin(7) + t195;
t162 = t174 * t242 + t190 * t244;
t221 = t240 * t243 - t245 * t241;
t249 = t221 * qJD(2);
t293 = -t223 * t245 - t243 * t224;
t166 = -qJD(1) * t249 + qJD(3) * t293;
t179 = t207 * pkin(3) + pkin(7) * t254;
t178 = t244 * t179;
t247 = -qJD(4) * t162 - t166 * t242 + t178;
t150 = pkin(4) * t207 + qJ(5) * t170 - qJD(5) * t203 + t247;
t201 = t217 * t242 - t267;
t156 = -qJ(5) * t201 + t162;
t210 = qJD(4) - t216;
t299 = t156 * t210 + t150;
t171 = qJD(4) * t203 - t242 * t254;
t268 = qJD(4) * t244;
t250 = t244 * t166 + t174 * t268 + t242 * t179 - t190 * t269;
t151 = -qJ(5) * t171 - qJD(5) * t201 + t250;
t161 = t244 * t174 - t190 * t242;
t155 = -qJ(5) * t203 + t161;
t154 = pkin(4) * t210 + t155;
t298 = -t154 * t210 + t151;
t297 = t240 * t270 + t241 * t271;
t198 = t226 * t245 + t243 * t227;
t292 = t203 ^ 2;
t290 = -qJ(5) - pkin(7);
t288 = t170 * t242;
t287 = t201 * t216;
t286 = t203 * t210;
t285 = t210 * t242;
t284 = t216 * t242;
t283 = t222 * t242;
t282 = t222 * t244;
t199 = -t226 * t243 + t227 * t245;
t196 = t244 * t199;
t206 = t244 * t207;
t279 = t154 - t155;
t278 = -t242 * t171 - t201 * t268;
t277 = t210 * t284 + t206;
t191 = pkin(3) * t217 - pkin(7) * t216;
t276 = t242 * t191 + t244 * t293;
t193 = pkin(3) * t221 - pkin(7) * t222 + t235;
t275 = t242 * t193 + t196;
t259 = qJD(4) * t290;
t274 = qJ(5) * t284 + qJD(5) * t244 + t242 * t259 - t276;
t183 = t244 * t191;
t273 = -pkin(4) * t217 - t183 + (qJ(5) * t216 + t259) * t244 + (-qJD(5) + t293) * t242;
t265 = qJD(1) * qJD(2);
t175 = -t198 * qJD(3) - t249;
t218 = t221 * qJD(3);
t192 = pkin(3) * t219 + pkin(7) * t218;
t264 = t244 * t175 + t242 * t192 + t193 * t268;
t260 = t222 * t268;
t258 = t210 * t244;
t167 = qJD(2) * t297 + t301;
t256 = qJ(5) * t218 - qJD(5) * t222;
t189 = -qJD(3) * pkin(3) - t293;
t253 = -t218 * t242 + t260;
t252 = -t218 * t244 - t222 * t269;
t158 = pkin(4) * t171 + t167;
t251 = -pkin(7) * t207 + t189 * t210;
t176 = qJD(2) * t222 + qJD(3) * t199;
t229 = t290 * t244;
t228 = t290 * t242;
t200 = t201 ^ 2;
t186 = t244 * t193;
t184 = t244 * t192;
t164 = pkin(4) * t201 + qJD(5) + t189;
t163 = -qJ(5) * t283 + t275;
t159 = pkin(4) * t221 - qJ(5) * t282 - t199 * t242 + t186;
t153 = -qJ(5) * t260 + (-qJD(4) * t199 + t256) * t242 + t264;
t152 = pkin(4) * t219 - t175 * t242 + t184 + t256 * t244 + (-t196 + (qJ(5) * t222 - t193) * t242) * qJD(4);
t1 = [(-t217 * t218 - t222 * t254) * MDP(8) + (-t222 * t207 - t218 * t216 - t217 * t219 + t221 * t254) * MDP(9) + (t207 * t235 + t219 * t225) * MDP(13) + (-t225 * t218 - t235 * t254) * MDP(14) + (-t170 * t282 + t203 * t252) * MDP(15) + (-(-t201 * t244 - t203 * t242) * t218 + (t288 - t171 * t244 + (t201 * t242 - t203 * t244) * qJD(4)) * t222) * MDP(16) + (-t170 * t221 + t203 * t219 + t206 * t222 + t210 * t252) * MDP(17) + (-t171 * t221 - t201 * t219 - t207 * t283 - t210 * t253) * MDP(18) + (t207 * t221 + t210 * t219) * MDP(19) + ((-t199 * t268 + t184) * t210 + t186 * t207 + (-t190 * t268 + t178) * t221 + t161 * t219 + t176 * t201 + t198 * t171 + t189 * t260 + ((-qJD(4) * t193 - t175) * t210 - t199 * t207 + (-qJD(4) * t174 - t166) * t221 + t167 * t222 - t189 * t218) * t242) * MDP(20) + (-(-t199 * t269 + t264) * t210 - t275 * t207 - t250 * t221 - t162 * t219 + t176 * t203 - t198 * t170 + t167 * t282 + t252 * t189) * MDP(21) + (-t152 * t203 - t153 * t201 + t159 * t170 - t163 * t171 - (-t154 * t244 - t156 * t242) * t218 + (-t150 * t244 - t151 * t242 + (t154 * t242 - t156 * t244) * qJD(4)) * t222) * MDP(22) + (t151 * t163 + t156 * t153 + t150 * t159 + t154 * t152 + t158 * (pkin(4) * t283 + t198) + t164 * (pkin(4) * t253 + t176)) * MDP(23) + 0.2e1 * t265 * t300 + (-t218 * MDP(10) - t219 * MDP(11) - t176 * MDP(13) - MDP(14) * t175) * qJD(3); t232 * MDP(14) + t277 * MDP(20) + t278 * MDP(22) + (-MDP(20) * t201 - MDP(21) * t203 - MDP(23) * t164) * t217 + (-qJD(4) * t210 * MDP(20) - t207 * MDP(21) + MDP(22) * t286 + MDP(23) * t298) * t242 + ((t170 + t287) * MDP(22) + t299 * MDP(23) - t210 ^ 2 * MDP(21)) * t244 + ((t217 + t297) * MDP(13) + (t216 - t263) * MDP(14)) * qJD(3) - qJD(1) ^ 2 * t300; -t216 ^ 2 * MDP(9) + (t232 + (-t216 - t263) * qJD(3)) * MDP(10) + (-t167 + t301) * MDP(13) + (-t216 * t225 + t221 * t265) * MDP(14) + (t203 * t258 - t288) * MDP(15) + ((-t170 + t287) * t244 - t203 * t285 + t278) * MDP(16) + (t207 * t242 + t210 * t258) * MDP(17) + (-t210 * t269 + t277) * MDP(18) + (-pkin(3) * t171 - t167 * t244 - t195 * t201 + (-pkin(7) * t268 - t183) * t210 + (t210 * t293 + t251) * t242) * MDP(20) + (pkin(3) * t170 + t167 * t242 - t195 * t203 + (pkin(7) * t269 + t276) * t210 + t251 * t244) * MDP(21) + (t170 * t228 + t171 * t229 - t274 * t201 - t273 * t203 - t242 * t299 + t244 * t298) * MDP(22) + (-t151 * t229 + t150 * t228 + t158 * (-pkin(4) * t244 - pkin(3)) + (pkin(4) * t285 - t195) * t164 + t274 * t156 + t273 * t154) * MDP(23) + (-t225 * MDP(13) - t203 * MDP(17) + t201 * MDP(18) - t210 * MDP(19) - t161 * MDP(20) + t162 * MDP(21) - t216 * MDP(8) + MDP(9) * t217) * t217; t203 * t201 * MDP(15) + (-t200 + t292) * MDP(16) + (t201 * t210 - t170) * MDP(17) + (-t171 + t286) * MDP(18) + t207 * MDP(19) + (t162 * t210 - t189 * t203 + t247) * MDP(20) + (t161 * t210 + t189 * t201 - t250) * MDP(21) + (pkin(4) * t170 - t201 * t279) * MDP(22) + (t279 * t156 + (-t164 * t203 + t150) * pkin(4)) * MDP(23); (-t200 - t292) * MDP(22) + (t154 * t203 + t156 * t201 + t158) * MDP(23);];
tauc = t1;
