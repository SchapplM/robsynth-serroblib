% Calculate Coriolis joint torque vector for
% S4RRRP6
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [4x1]
%   Generalized joint coordinates (joint angles)
% qJD [4x1]
%   Generalized joint velocities
% pkin [6x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,d1,d2,d3]';
% MDP [21x1]
%   Minimal dynamic parameter vector (fixed base model)
%   see S4RRRP6_convert_par2_MPV_fixb.m
% 
% Output:
% tauc [4x1]
%   joint torques required to compensate Coriolis and centrifugal load

% Quelle: HybrDyn-Toolbox
% Datum: 2021-01-15 14:39
% Revision: 24b2e7d74a0c1a3b64fa2f8f5ad758691ad61af3 (2021-01-15)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauc = S4RRRP6_coriolisvecJ_fixb_mdp_slag_vp(qJ, qJD, pkin, MDP)
%% Coder Information
%#codegen
%$cgargs {zeros(4,1),zeros(4,1),zeros(6,1),zeros(21,1)}
assert(isreal(qJ) && all(size(qJ) == [4 1]), ...
  'S4RRRP6_coriolisvecJ_fixb_mdp_slag_vp: qJ has to be [4x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [4 1]), ...
  'S4RRRP6_coriolisvecJ_fixb_mdp_slag_vp: qJD has to be [4x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [6 1]), ...
  'S4RRRP6_coriolisvecJ_fixb_mdp_slag_vp: pkin has to be [6x1] (double)');
assert(isreal(MDP) && all(size(MDP) == [21 1]), ...
  'S4RRRP6_coriolisvecJ_fixb_mdp_slag_vp: MDP has to be [21x1] (double)'); 

%% Symbolic Calculation
% From coriolisvec_joint_fixb_mdp_matlab.m
% OptimizationMode: 2
% StartTime: 2021-01-15 14:39:20
% EndTime: 2021-01-15 14:39:25
% DurationCPUTime: 1.68s
% Computational Cost: add. (1004->253), mult. (2551->362), div. (0->0), fcn. (1471->4), ass. (0->107)
t225 = (qJD(1) * qJD(2));
t271 = -2 * t225;
t201 = sin(qJ(2));
t270 = t201 * MDP(4);
t198 = t201 ^ 2;
t203 = cos(qJ(2));
t269 = (-t203 ^ 2 + t198) * MDP(5);
t180 = -pkin(2) * t203 - pkin(6) * t201 - pkin(1);
t165 = t180 * qJD(1);
t233 = qJD(1) * t203;
t196 = pkin(5) * t233;
t186 = qJD(2) * pkin(6) + t196;
t200 = sin(qJ(3));
t202 = cos(qJ(3));
t148 = t202 * t165 - t186 * t200;
t232 = qJD(2) * t200;
t234 = qJD(1) * t201;
t175 = t202 * t234 + t232;
t142 = -qJ(4) * t175 + t148;
t191 = -qJD(3) + t233;
t141 = -pkin(3) * t191 + t142;
t268 = t141 - t142;
t229 = qJD(3) * t200;
t222 = t201 * t229;
t219 = t203 * t225;
t224 = qJD(2) * qJD(3);
t238 = (t219 + t224) * t202;
t152 = qJD(1) * t222 - t238;
t223 = t200 * t234;
t226 = t202 * qJD(2);
t173 = t223 - t226;
t267 = t173 * t191 - t152;
t228 = qJD(3) * t202;
t221 = t201 * t228;
t230 = qJD(2) * t203;
t206 = t200 * t230 + t221;
t153 = t206 * qJD(1) + t200 * t224;
t266 = -t175 * t191 + t153;
t265 = pkin(3) * t173;
t264 = pkin(3) * t200;
t263 = pkin(5) * t200;
t262 = qJ(4) + pkin(6);
t261 = qJ(4) * t201;
t250 = t200 * t165;
t149 = t186 * t202 + t250;
t143 = -qJ(4) * t173 + t149;
t260 = t143 * t191;
t145 = pkin(3) * t153 + pkin(5) * t219;
t259 = t145 * t200;
t258 = t145 * t202;
t257 = t152 * t200;
t254 = t175 * t200;
t185 = -qJD(2) * pkin(2) + pkin(5) * t234;
t253 = t185 * t200;
t252 = t185 * t202;
t251 = t191 * t202;
t249 = t200 * t203;
t248 = t201 * t202;
t204 = qJD(2) ^ 2;
t247 = t201 * t204;
t246 = t202 * t203;
t245 = t203 * t204;
t205 = qJD(1) ^ 2;
t244 = t203 * t205;
t208 = pkin(3) * t201 - qJ(4) * t246;
t218 = qJD(3) * t262;
t211 = pkin(2) * t201 - pkin(6) * t203;
t177 = t211 * qJD(1);
t237 = pkin(5) * t223 + t202 * t177;
t243 = t208 * qJD(1) + qJD(4) * t200 + t202 * t218 + t237;
t161 = t200 * t177;
t227 = qJD(4) * t202;
t242 = t161 + (-pkin(5) * t248 - qJ(4) * t249) * qJD(1) + t200 * t218 - t227;
t178 = t211 * qJD(2);
t166 = qJD(1) * t178;
t220 = t201 * t225;
t213 = pkin(5) * t220;
t241 = t202 * t166 + t200 * t213;
t240 = t200 * t178 + t180 * t228;
t231 = qJD(2) * t201;
t239 = t202 * t178 + t231 * t263;
t192 = pkin(5) * t246;
t236 = t200 * t180 + t192;
t217 = -qJD(4) - t265;
t216 = pkin(1) * t271;
t215 = t173 + t226;
t214 = -t175 + t232;
t212 = t165 * t228 + t200 * t166 - t186 * t229 - t202 * t213;
t210 = qJ(4) * t152 + t241;
t209 = qJD(1) * t198 - t191 * t203;
t207 = -qJ(4) * t153 + t212;
t194 = -pkin(3) * t202 - pkin(2);
t184 = t262 * t202;
t183 = t262 * t200;
t179 = (pkin(5) + t264) * t201;
t172 = t202 * t180;
t170 = t173 ^ 2;
t168 = t233 * t264 + t196;
t154 = t206 * pkin(3) + pkin(5) * t230;
t151 = t185 - t217;
t150 = -t200 * t261 + t236;
t147 = -qJ(4) * t248 + t172 + (-pkin(3) - t263) * t203;
t140 = (-pkin(5) * qJD(2) - qJ(4) * qJD(3)) * t248 + (-qJD(4) * t201 + (-pkin(5) * qJD(3) - qJ(4) * qJD(2)) * t203) * t200 + t240;
t139 = -t201 * t227 + t208 * qJD(2) + (-t192 + (-t180 + t261) * t200) * qJD(3) + t239;
t138 = -qJD(4) * t173 + t207;
t137 = pkin(3) * t220 - t149 * qJD(3) - qJD(4) * t175 + t210;
t1 = [0.2e1 * t219 * t270 + t269 * t271 + MDP(6) * t245 - MDP(7) * t247 + (-pkin(5) * t245 + t201 * t216) * MDP(9) + (pkin(5) * t247 + t203 * t216) * MDP(10) + (-t152 * t248 + (t203 * t226 - t222) * t175) * MDP(11) + ((-t173 * t202 - t254) * t230 + (t257 - t153 * t202 + (t173 * t200 - t175 * t202) * qJD(3)) * t201) * MDP(12) + (t191 * t222 + t152 * t203 + (t175 * t201 + t209 * t202) * qJD(2)) * MDP(13) + (t191 * t221 + t153 * t203 + (-t173 * t201 - t209 * t200) * qJD(2)) * MDP(14) + (-t191 - t233) * MDP(15) * t231 + (-(-t180 * t229 + t239) * t191 + (t185 * t228 + pkin(5) * t153 + (qJD(1) * t172 + t148) * qJD(2)) * t201 + ((pkin(5) * t173 + t253) * qJD(2) + (t250 + (pkin(5) * t191 + t186) * t202) * qJD(3) - t241) * t203) * MDP(16) + ((-pkin(5) * t203 * t229 + t240) * t191 + t212 * t203 + (-pkin(5) * t152 - t185 * t229) * t201 + ((pkin(5) * t175 + t252) * t203 + (-t236 * qJD(1) - t149 + (-t191 + t233) * t202 * pkin(5)) * t201) * qJD(2)) * MDP(17) + (-t139 * t191 + t153 * t179 + t154 * t173 + (t151 * t232 - t137) * t203 + (t151 * t228 + t259 + (qJD(1) * t147 + t141) * qJD(2)) * t201) * MDP(18) + (t140 * t191 - t152 * t179 + t154 * t175 + (t151 * t226 + t138) * t203 + (-t151 * t229 + t258 + (-qJD(1) * t150 - t143) * qJD(2)) * t201) * MDP(19) + (-t139 * t175 - t140 * t173 + t147 * t152 - t150 * t153 + (-t141 * t202 - t143 * t200) * t230 + (-t137 * t202 - t138 * t200 + (t141 * t200 - t143 * t202) * qJD(3)) * t201) * MDP(20) + (t137 * t147 + t138 * t150 + t139 * t141 + t140 * t143 + t145 * t179 + t151 * t154) * MDP(21); -t244 * t270 + t205 * t269 + (-t175 * t251 - t257) * MDP(11) + (-t200 * t266 + t202 * t267) * MDP(12) + (-t191 * t228 + (t191 * t246 + t214 * t201) * qJD(1)) * MDP(13) + (t191 * t229 + (-t191 * t249 + t215 * t201) * qJD(1)) * MDP(14) + t191 * MDP(15) * t234 + (-pkin(2) * t153 + t237 * t191 + (pkin(6) * t251 + t253) * qJD(3) + ((-pkin(6) * t232 - t148) * t201 + (-t215 * pkin(5) - t253) * t203) * qJD(1)) * MDP(16) + (pkin(2) * t152 - t161 * t191 + (-pkin(6) * t191 * t200 + t252) * qJD(3) + (-t185 * t246 + (-pkin(6) * t226 + t149) * t201 + (t191 * t248 + t214 * t203) * pkin(5)) * qJD(1)) * MDP(17) + (-t258 + t153 * t194 - t168 * t173 + t243 * t191 + (t151 + t265) * t229 + (-t151 * t249 + (-qJD(2) * t183 - t141) * t201) * qJD(1)) * MDP(18) + (t259 - t152 * t194 - t168 * t175 - t242 * t191 + (pkin(3) * t254 + t151 * t202) * qJD(3) + (-t151 * t246 + (-qJD(2) * t184 + t143) * t201) * qJD(1)) * MDP(19) + (-t152 * t183 - t153 * t184 + t243 * t175 + t242 * t173 + (t191 * t141 + t138) * t202 + (-t137 + t260) * t200) * MDP(20) + (-t137 * t183 + t138 * t184 + t145 * t194 + (pkin(3) * t229 - t168) * t151 - t242 * t143 - t243 * t141) * MDP(21) + (t205 * t201 * MDP(9) + MDP(10) * t244) * pkin(1); -t170 * MDP(12) + t238 * MDP(13) + (-t149 * t191 + t241) * MDP(16) + (-t148 * t191 - t212) * MDP(17) + (t210 - t260) * MDP(18) + (-t142 * t191 - t207) * MDP(19) + t152 * pkin(3) * MDP(20) + (pkin(3) * t137 + t143 * t268) * MDP(21) + (-t191 * MDP(13) + t185 * MDP(17) + (qJD(4) + t151) * MDP(19) - t268 * MDP(20)) * t173 + (-MDP(14) * t249 + (0.2e1 * MDP(18) * pkin(3) + MDP(15)) * t201) * t225 + (-MDP(13) * t223 - t175 * MDP(14) + (-MDP(16) - MDP(18)) * t149) * qJD(3) + (t173 * MDP(11) - t191 * MDP(14) - t185 * MDP(16) + (-t151 + t217) * MDP(18) - pkin(3) * t151 * MDP(21) + (-MDP(19) * pkin(3) + MDP(12)) * t175) * t175; t266 * MDP(18) + t267 * MDP(19) + (-t175 ^ 2 - t170) * MDP(20) + (t141 * t175 + t143 * t173 + t145) * MDP(21);];
tauc = t1;
