% Calculate minimal parameter regressor of Coriolis joint torque vector for
% S5PPRRR4
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% qJD [5x1]
%   Generalized joint velocities
% pkin [11x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,alpha2,alpha3,d3,d4,d5,theta1,theta2]';
% MDP [19x1]
%   Minimal dynamic parameter vector (fixed base model)
%   see S5PPRRR4_convert_par2_MPV_fixb.m
% 
% Output:
% tauc [5x1]
%   joint torques required to compensate Coriolis and centrifugal load

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-05 15:20
% Revision: 77da58f92bca3eff71542919beafa37024070d86 (2019-12-05)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauc = S5PPRRR4_coriolisvecJ_fixb_mdp_slag_vp(qJ, qJD, pkin, MDP)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(11,1),zeros(19,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5PPRRR4_coriolisvecJ_fixb_mdp_slag_vp: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5PPRRR4_coriolisvecJ_fixb_mdp_slag_vp: qJD has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [11 1]), ...
  'S5PPRRR4_coriolisvecJ_fixb_mdp_slag_vp: pkin has to be [11x1] (double)');
assert(isreal(MDP) && all(size(MDP) == [19 1]), ...
  'S5PPRRR4_coriolisvecJ_fixb_mdp_slag_vp: MDP has to be [19x1] (double)'); 

%% Symbolic Calculation
% From coriolisvec_joint_fixb_mdp_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-05 15:19:48
% EndTime: 2019-12-05 15:19:54
% DurationCPUTime: 1.90s
% Computational Cost: add. (869->212), mult. (2423->347), div. (0->0), fcn. (2038->12), ass. (0->112)
t191 = cos(pkin(5));
t187 = sin(pkin(6));
t197 = cos(qJ(3));
t255 = t187 * t197;
t186 = sin(pkin(11));
t188 = sin(pkin(5));
t194 = sin(qJ(3));
t189 = cos(pkin(11));
t190 = cos(pkin(6));
t254 = t189 * t190;
t265 = (-t186 * t194 + t197 * t254) * t188;
t268 = t191 * t255 + t265;
t196 = cos(qJ(4));
t246 = MDP(12) * t196;
t180 = qJD(1) * t191 + qJD(2);
t267 = qJD(1) * t265 + t180 * t255;
t193 = sin(qJ(4));
t266 = MDP(6) * t193;
t184 = t193 ^ 2;
t264 = (-t196 ^ 2 + t184) * MDP(7);
t245 = qJD(1) * t188;
t225 = t189 * t245;
t216 = t190 * t225;
t226 = t186 * t245;
t256 = t187 * t194;
t229 = t180 * t256;
t242 = qJD(3) * t197;
t152 = t226 * t242 + (t194 * t216 + t229) * qJD(3);
t202 = (t186 * t197 + t194 * t254) * t188;
t156 = qJD(1) * t202 + t229;
t263 = qJD(3) * t156 - t152;
t155 = -t194 * t226 + t197 * (t180 * t187 + t216);
t262 = qJD(3) * pkin(3);
t154 = qJD(3) * pkin(8) + t156;
t165 = t180 * t190 - t187 * t225;
t145 = t154 * t196 + t165 * t193;
t151 = t267 * qJD(3);
t139 = qJD(4) * t145 + t151 * t193;
t192 = sin(qJ(5));
t261 = t139 * t192;
t195 = cos(qJ(5));
t260 = t139 * t195;
t238 = qJD(5) * t192;
t221 = t193 * t238;
t231 = qJD(3) * qJD(4);
t219 = t196 * t231;
t230 = qJD(4) * qJD(5);
t248 = (t219 + t230) * t195;
t163 = -qJD(3) * t221 + t248;
t259 = t163 * t192;
t234 = t192 * qJD(3);
t223 = t193 * t234;
t232 = t195 * qJD(4);
t171 = t223 - t232;
t243 = qJD(3) * t196;
t181 = -qJD(5) + t243;
t258 = t171 * t181;
t233 = t195 * qJD(3);
t241 = qJD(4) * t192;
t173 = t193 * t233 + t241;
t257 = t173 * t181;
t253 = t192 * t181;
t252 = t192 * t196;
t251 = t195 * t181;
t250 = t195 * t196;
t215 = pkin(4) * t193 - pkin(9) * t196;
t176 = t215 * qJD(4);
t249 = t156 - t176;
t240 = qJD(4) * t193;
t239 = qJD(4) * t196;
t237 = qJD(5) * t195;
t236 = qJD(5) * t197;
t211 = t154 * t193 - t165 * t196;
t142 = -qJD(4) * pkin(4) + t211;
t235 = t142 * qJD(5);
t224 = t187 * t242;
t222 = t181 * t237;
t220 = t193 * t237;
t218 = MDP(17) * t240;
t143 = qJD(4) * pkin(9) + t145;
t178 = -pkin(4) * t196 - pkin(9) * t193 - pkin(3);
t150 = qJD(3) * t178 - t155;
t137 = t143 * t195 + t150 * t192;
t214 = t143 * t192 - t150 * t195;
t160 = t191 * t256 + t202;
t166 = -t187 * t188 * t189 + t190 * t191;
t149 = t160 * t196 + t166 * t193;
t213 = t149 * t195 - t192 * t268;
t212 = -t149 * t192 - t195 * t268;
t148 = t160 * t193 - t166 * t196;
t210 = qJD(3) * t184 - t181 * t196;
t198 = qJD(4) ^ 2;
t209 = pkin(8) * t198 - t263;
t153 = -t155 - t262;
t208 = qJD(4) * (t153 + t155 - t262);
t168 = t190 * t193 + t196 * t256;
t167 = -t190 * t196 + t193 * t256;
t206 = -MDP(11) * t196 + MDP(12) * t193 - MDP(4);
t138 = -qJD(4) * t211 + t151 * t196;
t201 = qJD(4) * t142 + qJD(5) * t150 - t155 * t181 + t138;
t199 = qJD(3) ^ 2;
t175 = t215 * qJD(3);
t164 = t192 * t230 + (t192 * t239 + t220) * qJD(3);
t162 = qJD(4) * t168 + t193 * t224;
t161 = -qJD(4) * t167 + t196 * t224;
t158 = t160 * qJD(3);
t157 = t268 * qJD(3);
t147 = qJD(3) * t176 + t152;
t146 = t195 * t147;
t141 = -qJD(4) * t148 + t157 * t196;
t140 = qJD(4) * t149 + t157 * t193;
t1 = [(-(-qJD(5) * t213 - t141 * t192 + t158 * t195) * t181 + t140 * t171 + t148 * t164) * MDP(18) + ((qJD(5) * t212 + t141 * t195 + t158 * t192) * t181 + t140 * t173 + t148 * t163) * MDP(19) + (-MDP(11) * t140 - MDP(12) * t141) * qJD(4) + (-t157 * MDP(5) + t206 * t158 + (-t268 * t246 + (-MDP(11) * t268 + MDP(18) * t212 - MDP(19) * t213) * t193) * qJD(4)) * qJD(3); (-(-t161 * t192 - t168 * t237) * t181 + t162 * t171 + t167 * t164) * MDP(18) + ((t161 * t195 - t168 * t238) * t181 + t162 * t173 + t167 * t163) * MDP(19) + ((-(t192 * t236 + t194 * t233) * MDP(18) + (t194 * t234 - t195 * t236) * MDP(19)) * t181 + (-t197 * MDP(5) + t194 * t206) * t199) * t187 + (-t162 * MDP(11) - t161 * MDP(12) + (-t246 * t255 + (-MDP(11) * t255 + (-t168 * t192 - t195 * t255) * MDP(18) - (t168 * t195 - t192 * t255) * MDP(19)) * t193) * qJD(3)) * qJD(4); t263 * MDP(4) + (t155 - t267) * qJD(3) * MDP(5) + 0.2e1 * t219 * t266 - 0.2e1 * t231 * t264 + (t193 * t208 - t196 * t209) * MDP(11) + (t193 * t209 + t196 * t208) * MDP(12) + (t163 * t193 * t195 + (t196 * t232 - t221) * t173) * MDP(13) + ((-t171 * t195 - t173 * t192) * t239 + (-t259 - t164 * t195 + (t171 * t192 - t173 * t195) * qJD(5)) * t193) * MDP(14) + (t181 * t221 - t163 * t196 + (t173 * t193 + t210 * t195) * qJD(4)) * MDP(15) + (t181 * t220 + t164 * t196 + (-t171 * t193 - t192 * t210) * qJD(4)) * MDP(16) + (-t181 - t243) * t218 + ((t178 * t238 + t249 * t195) * t181 + (t143 * t237 - t146 + (qJD(4) * t171 + t222) * pkin(8) + t201 * t192) * t196 + (t195 * t235 + pkin(8) * t164 + t261 - t155 * t171 + (-pkin(8) * t253 + (-pkin(8) * t252 + t178 * t195) * qJD(3) - t214) * qJD(4)) * t193) * MDP(18) + ((t178 * t237 - t192 * t249) * t181 + (qJD(4) * pkin(8) * t173 + (t147 + (-pkin(8) * t181 - t143) * qJD(5)) * t192 + t201 * t195) * t196 + (-t192 * t235 + pkin(8) * t163 + t260 - t155 * t173 + (-pkin(8) * t251 - (pkin(8) * t250 + t178 * t192) * qJD(3) - t137) * qJD(4)) * t193) * MDP(19) + (t196 * MDP(8) - t193 * MDP(9)) * t198; (-t173 * t251 + t259) * MDP(13) + ((t163 + t258) * t195 + (-t164 + t257) * t192) * MDP(14) + (-t222 + (t181 * t250 + (-t173 + t241) * t193) * qJD(3)) * MDP(15) + (t181 * t238 + (-t181 * t252 + (t171 + t232) * t193) * qJD(3)) * MDP(16) + t181 * t193 * qJD(3) * MDP(17) + (-pkin(4) * t164 - t260 + (t175 * t195 + t192 * t211) * t181 - t145 * t171 + (pkin(9) * t251 + t142 * t192) * qJD(5) + (t214 * t193 + (-pkin(9) * t240 - t142 * t196) * t192) * qJD(3)) * MDP(18) + (-pkin(4) * t163 + t261 - (t175 * t192 - t195 * t211) * t181 - t145 * t173 + (-pkin(9) * t253 + t142 * t195) * qJD(5) + (-t142 * t250 + (-pkin(9) * t232 + t137) * t193) * qJD(3)) * MDP(19) + (-t196 * t266 + t264) * t199 + (t193 * MDP(11) + t246) * (-qJD(3) * t153 - t151); t173 * t171 * MDP(13) + (-t171 ^ 2 + t173 ^ 2) * MDP(14) + (t248 - t258) * MDP(15) + (-t192 * t219 - t257) * MDP(16) + qJD(3) * t218 + (-t137 * t181 - t138 * t192 - t142 * t173 + t146) * MDP(18) + (-t195 * t138 + t142 * t171 - t192 * t147 + t181 * t214) * MDP(19) + (-MDP(15) * t223 - MDP(16) * t173 - MDP(18) * t137 + MDP(19) * t214) * qJD(5);];
tauc = t1;
