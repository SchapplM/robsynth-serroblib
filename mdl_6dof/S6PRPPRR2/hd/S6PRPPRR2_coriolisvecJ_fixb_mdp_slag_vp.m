% Calculate minimal parameter regressor of Coriolis joint torque vector for
% S6PRPPRR2
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% qJD [6x1]
%   Generalized joint velocities
% pkin [11x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,alpha2,d2,d5,d6,theta1,theta3]';
% MDP [22x1]
%   Minimal dynamic parameter vector (fixed base model)
%   see S6PRPPRR2_convert_par2_MPV_fixb.m
% 
% Output:
% tauc [6x1]
%   joint torques required to compensate Coriolis and centrifugal load

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-08 19:20
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauc = S6PRPPRR2_coriolisvecJ_fixb_mdp_slag_vp(qJ, qJD, pkin, MDP)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(11,1),zeros(22,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6PRPPRR2_coriolisvecJ_fixb_mdp_slag_vp: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6PRPPRR2_coriolisvecJ_fixb_mdp_slag_vp: qJD has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [11 1]), ...
  'S6PRPPRR2_coriolisvecJ_fixb_mdp_slag_vp: pkin has to be [11x1] (double)');
assert(isreal(MDP) && all(size(MDP) == [22 1]), ...
  'S6PRPPRR2_coriolisvecJ_fixb_mdp_slag_vp: MDP has to be [22x1] (double)'); 

%% Symbolic Calculation
% From coriolisvec_joint_fixb_mdp_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-08 19:20:04
% EndTime: 2019-03-08 19:20:08
% DurationCPUTime: 2.00s
% Computational Cost: add. (868->229), mult. (2144->329), div. (0->0), fcn. (1566->10), ass. (0->117)
t190 = sin(pkin(11));
t191 = sin(pkin(6));
t192 = cos(pkin(11));
t196 = sin(qJ(2));
t199 = cos(qJ(2));
t166 = (-t190 * t196 + t192 * t199) * t191;
t164 = qJD(2) * t166;
t246 = qJD(1) * t191;
t230 = t199 * t246;
t178 = qJD(2) * pkin(2) + t230;
t231 = t196 * t246;
t157 = t178 * t190 + t192 * t231;
t151 = qJD(2) * qJ(4) + t157;
t167 = (t190 * t199 + t192 * t196) * t191;
t162 = qJD(1) * t167;
t271 = -t151 + t162;
t198 = cos(qJ(5));
t195 = sin(qJ(5));
t247 = MDP(15) * t195;
t270 = -MDP(14) * t198 + t247;
t197 = cos(qJ(6));
t242 = qJD(2) * t198;
t227 = t197 * t242;
t194 = sin(qJ(6));
t240 = qJD(5) * t194;
t174 = t227 + t240;
t269 = t174 * MDP(16);
t189 = t198 ^ 2;
t268 = (t195 ^ 2 - t189) * MDP(10);
t267 = pkin(2) * t190;
t179 = t190 * t231;
t156 = t178 * t192 - t179;
t214 = qJD(4) - t156;
t149 = (-pkin(3) - pkin(8)) * qJD(2) + t214;
t193 = cos(pkin(6));
t182 = qJD(1) * t193 + qJD(3);
t145 = t149 * t195 + t182 * t198;
t163 = qJD(2) * t167;
t158 = qJD(1) * t163;
t139 = qJD(5) * t145 - t158 * t198;
t266 = t139 * t194;
t265 = t139 * t197;
t264 = t158 * t166;
t235 = t197 * qJD(5);
t186 = qJD(6) * t235;
t237 = qJD(6) * t194;
t225 = t198 * t237;
t204 = -t195 * t235 - t225;
t160 = t204 * qJD(2) + t186;
t263 = t160 * t194;
t233 = qJD(2) * qJD(5);
t223 = t195 * t233;
t181 = t194 * t223;
t161 = t174 * qJD(6) - t181;
t262 = t161 * t195;
t228 = t194 * t242;
t172 = t228 - t235;
t243 = qJD(2) * t195;
t184 = qJD(6) + t243;
t261 = t172 * t184;
t260 = t172 * t198;
t259 = t174 * t184;
t258 = t184 * t194;
t257 = t184 * t197;
t254 = t194 * t195;
t253 = t195 * t197;
t238 = qJD(5) * t198;
t252 = t160 * t195 + t174 * t238;
t165 = t192 * t230 - t179;
t216 = pkin(5) * t198 + pkin(9) * t195;
t171 = t216 * qJD(5) + qJD(4);
t251 = t165 - t171;
t159 = t164 * qJD(1);
t200 = qJD(5) ^ 2;
t201 = qJD(2) ^ 2;
t249 = -t200 - t201;
t245 = qJD(2) * t189;
t232 = -pkin(2) * t192 - pkin(3);
t183 = -pkin(8) + t232;
t241 = qJD(5) * t183;
t236 = qJD(6) * t197;
t234 = qJD(4) - t165;
t229 = t194 * t245;
t226 = t184 * t240;
t224 = t184 * t236;
t222 = MDP(20) * t238;
t143 = qJD(5) * pkin(9) + t145;
t221 = t183 * t184 + t143;
t218 = t184 * t225;
t185 = qJ(4) + t267;
t215 = qJD(2) * t185 - t271;
t208 = pkin(5) * t195 - pkin(9) * t198 + qJ(4);
t148 = t208 * qJD(2) + t157;
t137 = t143 * t197 + t148 * t194;
t213 = t143 * t194 - t148 * t197;
t144 = t149 * t198 - t182 * t195;
t153 = -t166 * t195 + t193 * t198;
t212 = t153 * t197 + t167 * t194;
t211 = -t153 * t194 + t167 * t197;
t210 = -t166 * t198 - t193 * t195;
t209 = -t184 * t195 + t245;
t207 = t195 * t226 - t198 * t224;
t142 = -qJD(5) * pkin(5) - t144;
t206 = -pkin(9) * t238 + t142 * t195;
t205 = (-MDP(21) * t197 + MDP(22) * t194) * t184;
t154 = qJD(4) * qJD(2) + t159;
t203 = t234 * qJD(2) - t183 * t200 + t154;
t138 = qJD(5) * t144 + t158 * t195;
t202 = -qJD(5) * t142 - qJD(6) * t148 + t162 * t184 - t138;
t177 = t216 * qJD(2);
t169 = t208 + t267;
t150 = -qJD(2) * pkin(3) + t214;
t147 = t171 * qJD(2) + t159;
t146 = t197 * t147;
t141 = t153 * qJD(5) - t163 * t198;
t140 = t210 * qJD(5) + t163 * t195;
t1 = [(-t156 * t163 + t157 * t164 + t159 * t167 - t264) * MDP(5) + (t150 * t163 + t151 * t164 + t154 * t167 - t264) * MDP(8) + ((-t212 * qJD(6) - t140 * t194 + t164 * t197) * t184 + t141 * t172 - t210 * t161) * MDP(21) + (-(t211 * qJD(6) + t140 * t197 + t164 * t194) * t184 + t141 * t174 - t210 * t160) * MDP(22) + (-MDP(3) * t196 - MDP(4) * t199) * t201 * t191 + (-MDP(14) * t141 - MDP(15) * t140) * qJD(5) + (t163 * MDP(6) + (MDP(14) * t195 + MDP(15) * t198 + MDP(7)) * t164 + (-t167 * t247 + (t167 * MDP(14) + t211 * MDP(21) - t212 * MDP(22)) * t198) * qJD(5)) * qJD(2); (t156 * t162 - t157 * t165 + (-t158 * t192 + t159 * t190) * pkin(2)) * MDP(5) + ((0.2e1 * qJD(4) - t165) * qJD(2) + t159) * MDP(7) + (-t150 * t162 + t234 * t151 + t154 * t185 + t158 * t232) * MDP(8) + 0.2e1 * t233 * t268 + t215 * t238 * MDP(14) + t204 * t269 + (t209 * t235 - t218 + t252) * MDP(18) + (-t262 + (-t229 - t260) * qJD(5) + t207) * MDP(19) + t243 * t222 + (-0.2e1 * MDP(9) * t223 - t200 * MDP(12) + t203 * MDP(15) + t160 * t197 * MDP(16) + (-t263 - t161 * t197 + (t172 * t194 - t174 * t197) * qJD(6)) * MDP(17) + (t142 * t236 + t266 - t183 * t161 + t162 * t172 + (-t183 * t258 + (t169 * t197 - t183 * t254) * qJD(2) - t213) * qJD(5)) * MDP(21) + (-t142 * t237 + t265 - t183 * t160 + t162 * t174 + (-t183 * t257 - (t169 * t194 + t183 * t253) * qJD(2) - t137) * qJD(5)) * MDP(22)) * t198 + (t222 + (-t169 * t237 - t251 * t197) * MDP(21) + (-t169 * t236 + t251 * t194) * MDP(22)) * t184 + (-t200 * MDP(11) + t203 * MDP(14) + (t172 * t241 + t202 * t194 - t221 * t236 + t146) * MDP(21) + (t174 * t241 + (t221 * qJD(6) - t147) * t194 + t202 * t197) * MDP(22) + (-t215 * MDP(15) + (t172 * t197 + t174 * t194) * MDP(17)) * qJD(5)) * t195; (t207 + t262) * MDP(21) + (t218 + t252) * MDP(22) + t270 * t200 + ((-t229 + t260) * MDP(21) - t209 * MDP(22) * t197) * qJD(5); -t201 * MDP(7) + (MDP(8) * t271 + t205) * qJD(2) + (t249 * MDP(15) + (-t161 - t226) * MDP(21) + (-t184 * t235 - t160) * MDP(22)) * t198 + (t249 * MDP(14) + qJD(6) * t205 + ((t172 - t228) * MDP(21) + (t174 - t227) * MDP(22)) * qJD(5)) * t195; (t174 * t257 + t263) * MDP(16) + ((t160 - t261) * t197 + (-t161 - t259) * t194) * MDP(17) + (t224 + (t184 * t253 + (-t174 + t240) * t198) * qJD(2)) * MDP(18) + (-t184 * t237 + (-t184 * t254 + (t172 + t235) * t198) * qJD(2)) * MDP(19) - t184 * MDP(20) * t242 + (-pkin(5) * t161 - t265 - (-t144 * t194 + t177 * t197) * t184 - t145 * t172 + (-pkin(9) * t257 + t142 * t194) * qJD(6) + (t206 * t194 + t198 * t213) * qJD(2)) * MDP(21) + (-pkin(5) * t160 + t266 + (t144 * t197 + t177 * t194) * t184 - t145 * t174 + (pkin(9) * t258 + t142 * t197) * qJD(6) + (t137 * t198 + t206 * t197) * qJD(2)) * MDP(22) + t270 * (qJD(2) * t151 - t158) + (t198 * t195 * MDP(9) - t268) * t201; t172 * t269 + (-t172 ^ 2 + t174 ^ 2) * MDP(17) + (-t197 * t223 + t186 + t261) * MDP(18) + (t181 + t259) * MDP(19) + qJD(2) * t222 + (t137 * t184 - t194 * t138 - t142 * t174 + t146) * MDP(21) + (-t197 * t138 + t142 * t172 - t194 * t147 - t184 * t213) * MDP(22) + (-MDP(18) * t228 - t174 * MDP(19) - t137 * MDP(21) + t213 * MDP(22)) * qJD(6);];
tauc  = t1;
