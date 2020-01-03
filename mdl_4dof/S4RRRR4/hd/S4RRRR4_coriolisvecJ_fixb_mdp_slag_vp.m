% Calculate minimal parameter regressor of Coriolis joint torque vector for
% S4RRRR4
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [4x1]
%   Generalized joint coordinates (joint angles)
% qJD [4x1]
%   Generalized joint velocities
% pkin [7x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,d1,d2,d3,d4]';
% MDP [24x1]
%   Minimal dynamic parameter vector (fixed base model)
%   see S4RRRR4_convert_par2_MPV_fixb.m
% 
% Output:
% tauc [4x1]
%   joint torques required to compensate Coriolis and centrifugal load

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-31 17:26
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauc = S4RRRR4_coriolisvecJ_fixb_mdp_slag_vp(qJ, qJD, pkin, MDP)
%% Coder Information
%#codegen
%$cgargs {zeros(4,1),zeros(4,1),zeros(7,1),zeros(24,1)}
assert(isreal(qJ) && all(size(qJ) == [4 1]), ...
  'S4RRRR4_coriolisvecJ_fixb_mdp_slag_vp: qJ has to be [4x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [4 1]), ...
  'S4RRRR4_coriolisvecJ_fixb_mdp_slag_vp: qJD has to be [4x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [7 1]), ...
  'S4RRRR4_coriolisvecJ_fixb_mdp_slag_vp: pkin has to be [7x1] (double)');
assert(isreal(MDP) && all(size(MDP) == [24 1]), ...
  'S4RRRR4_coriolisvecJ_fixb_mdp_slag_vp: MDP has to be [24x1] (double)'); 

%% Symbolic Calculation
% From coriolisvec_joint_fixb_mdp_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 17:26:18
% EndTime: 2019-12-31 17:26:22
% DurationCPUTime: 1.78s
% Computational Cost: add. (1276->212), mult. (3236->310), div. (0->0), fcn. (2234->6), ass. (0->104)
t212 = sin(qJ(3));
t215 = cos(qJ(2));
t269 = cos(qJ(3));
t240 = qJD(1) * t269;
t213 = sin(qJ(2));
t250 = qJD(1) * t213;
t275 = -t212 * t250 + t215 * t240;
t208 = qJD(2) + qJD(3);
t246 = qJD(1) * qJD(2);
t274 = -0.2e1 * t246;
t273 = MDP(5) * (t213 ^ 2 - t215 ^ 2);
t257 = t212 * t215;
t189 = t269 * t213 + t257;
t272 = qJD(1) * t189;
t270 = pkin(5) + pkin(6);
t195 = t270 * t213;
t190 = qJD(1) * t195;
t268 = qJD(2) * pkin(2);
t185 = -t190 + t268;
t243 = qJD(2) * t270;
t231 = qJD(1) * t243;
t186 = t213 * t231;
t187 = t215 * t231;
t196 = t270 * t215;
t192 = qJD(1) * t196;
t239 = qJD(3) * t269;
t249 = qJD(3) * t212;
t147 = t185 * t239 - t269 * t186 - t212 * t187 - t192 * t249;
t242 = t269 * t192;
t169 = t212 * t185 + t242;
t236 = -t212 * t186 + t269 * t187;
t148 = t169 * qJD(3) + t236;
t191 = t213 * t243;
t193 = t215 * t243;
t223 = -t269 * t195 - t212 * t196;
t152 = t223 * qJD(3) - t269 * t191 - t212 * t193;
t184 = -qJD(1) * t257 - t213 * t240;
t207 = -pkin(2) * t215 - pkin(1);
t194 = t207 * qJD(1);
t159 = -pkin(3) * t275 + pkin(7) * t184 + t194;
t258 = t212 * t192;
t168 = t269 * t185 - t258;
t161 = -t208 * pkin(3) - t168;
t173 = t208 * t189;
t165 = t173 * qJD(1);
t222 = -t212 * t213 + t269 * t215;
t167 = -pkin(3) * t222 - pkin(7) * t189 + t207;
t172 = t208 * t222;
t178 = -t212 * t195 + t269 * t196;
t180 = qJD(4) - t275;
t271 = t148 * t189 + t161 * t172 - t178 * t165 - (qJD(4) * t167 + t152) * t180 + (qJD(4) * t159 + t147) * t222;
t211 = sin(qJ(4));
t248 = qJD(4) * t211;
t164 = t275 * t208;
t214 = cos(qJ(4));
t247 = qJD(4) * t214;
t252 = t214 * t164 + t208 * t247;
t150 = t184 * t248 + t252;
t267 = t150 * t211;
t266 = t161 * t275;
t265 = t161 * t189;
t264 = t164 * t211;
t263 = t167 * t165;
t260 = t184 * t211;
t174 = -t214 * t208 - t260;
t262 = t174 * t180;
t227 = t184 * t214 - t208 * t211;
t261 = t227 * t180;
t259 = t211 * t165;
t216 = qJD(2) ^ 2;
t256 = t213 * t216;
t255 = t214 * t165;
t254 = t215 * t216;
t217 = qJD(1) ^ 2;
t253 = t215 * t217;
t245 = pkin(2) * t250;
t244 = t213 * t268;
t238 = t213 * t246;
t237 = pkin(1) * t274;
t235 = t214 * t180;
t166 = -pkin(3) * t184 - pkin(7) * t275;
t205 = pkin(2) * t212 + pkin(7);
t232 = qJD(4) * t205 + t166 + t245;
t170 = -t212 * t190 + t242;
t230 = pkin(2) * t249 - t170;
t162 = t208 * pkin(7) + t169;
t143 = t159 * t211 + t162 * t214;
t229 = -t143 * t184 + t148 * t211 + t161 * t247;
t142 = t159 * t214 - t162 * t211;
t228 = -t165 * t205 - t266;
t171 = -t269 * t190 - t258;
t226 = -pkin(2) * t239 + t171;
t225 = t142 * t184 - t148 * t214 + t161 * t248;
t224 = t194 * t184 - t236;
t221 = t172 * t214 - t189 * t248;
t151 = -t227 * qJD(4) + t264;
t219 = ((t150 - t262) * t214 + (-t151 + t261) * t211) * MDP(19) + (-t227 * t235 + t267) * MDP(18) + (-t180 ^ 2 * t211 - t174 * t184 + t255) * MDP(21) + (t180 * t235 - t184 * t227 + t259) * MDP(20) + t164 * MDP(13) + (t184 ^ 2 - t275 ^ 2) * MDP(12) + (MDP(11) * t275 + t180 * MDP(22)) * t184 + (-t275 * MDP(13) + (-t184 - t272) * MDP(14)) * t208;
t218 = -t194 * t275 - t147;
t206 = -t269 * pkin(2) - pkin(3);
t153 = t178 * qJD(3) - t212 * t191 + t269 * t193;
t149 = pkin(3) * t173 - pkin(7) * t172 + t244;
t145 = pkin(2) * t238 + pkin(3) * t165 - pkin(7) * t164;
t144 = t214 * t145;
t1 = [0.2e1 * t215 * MDP(4) * t238 + t273 * t274 + MDP(6) * t254 - MDP(7) * t256 + (-pkin(5) * t254 + t213 * t237) * MDP(9) + (pkin(5) * t256 + t215 * t237) * MDP(10) + (t164 * t189 - t172 * t184) * MDP(11) + (t164 * t222 - t165 * t189 + t172 * t275 + t173 * t184) * MDP(12) + (t165 * t207 + t173 * t194 + (-qJD(1) * t222 - t275) * t244) * MDP(16) + (t164 * t207 + t172 * t194 + (-t184 + t272) * t244) * MDP(17) + (t150 * t189 * t214 - t221 * t227) * MDP(18) + ((-t174 * t214 + t211 * t227) * t172 + (-t267 - t151 * t214 + (t174 * t211 + t214 * t227) * qJD(4)) * t189) * MDP(19) + (-t150 * t222 - t173 * t227 + t221 * t180 + t189 * t255) * MDP(20) + (-t189 * t259 + t151 * t222 - t173 * t174 + (-t172 * t211 - t189 * t247) * t180) * MDP(21) + (-t165 * t222 + t173 * t180) * MDP(22) + (t142 * t173 - t144 * t222 - t223 * t151 + t153 * t174 + (t149 * t180 + t263 + (t162 * t222 - t178 * t180 + t265) * qJD(4)) * t214 + t271 * t211) * MDP(23) + (-t143 * t173 - t223 * t150 - t153 * t227 + (-(-qJD(4) * t178 + t149) * t180 - t263 + (-qJD(4) * t162 + t145) * t222 - qJD(4) * t265) * t211 + t271 * t214) * MDP(24) + (t172 * MDP(13) - t173 * MDP(14) - t153 * MDP(16) - t152 * MDP(17)) * t208; (t206 * t151 + t228 * t211 + t230 * t174 + (t226 * t211 - t232 * t214) * t180 + t225) * MDP(23) + (t206 * t150 + t228 * t214 - t230 * t227 + (t232 * t211 + t226 * t214) * t180 + t229) * MDP(24) + (t171 * t208 + (t184 * t250 - t208 * t239) * pkin(2) + t218) * MDP(17) + (t275 * t245 + t170 * t208 + (-t242 + (-pkin(2) * t208 - t185) * t212) * qJD(3) + t224) * MDP(16) - t213 * MDP(4) * t253 + t219 + t217 * t273 + (t217 * t213 * MDP(9) + MDP(10) * t253) * pkin(1); (t224 + (-qJD(3) + t208) * t169) * MDP(16) + (t168 * t208 + t218) * MDP(17) + (-pkin(3) * t151 - (t166 * t214 - t168 * t211) * t180 - t169 * t174 - t211 * t266 + (-t180 * t247 - t259) * pkin(7) + t225) * MDP(23) + (-pkin(3) * t150 + (t166 * t211 + t168 * t214) * t180 + t169 * t227 - t214 * t266 + (t180 * t248 - t255) * pkin(7) + t229) * MDP(24) + t219; -t227 * t174 * MDP(18) + (-t174 ^ 2 + t227 ^ 2) * MDP(19) + (t252 + t262) * MDP(20) + (-t261 - t264) * MDP(21) + t165 * MDP(22) + (t143 * t180 - t147 * t211 + t161 * t227 + t144) * MDP(23) + (t142 * t180 - t145 * t211 - t147 * t214 + t161 * t174) * MDP(24) + (MDP(20) * t260 + t227 * MDP(21) - t143 * MDP(23) - t142 * MDP(24)) * qJD(4);];
tauc = t1;
