% Calculate Coriolis joint torque vector for
% S5RPRRP1
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% qJD [5x1]
%   Generalized joint velocities
% pkin [7x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d1,d3,d4]';
% MDP [24x1]
%   Minimal dynamic parameter vector (fixed base model)
%   see S5RPRRP1_convert_par2_MPV_fixb.m
% 
% Output:
% tauc [5x1]
%   joint torques required to compensate Coriolis and centrifugal load

% Quelle: HybrDyn-Toolbox
% Datum: 2021-01-15 12:27
% Revision: d12c3222fdeb2c5f3b3c8fa5751e113be2fc3aae (2021-01-15)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauc = S5RPRRP1_coriolisvecJ_fixb_mdp_slag_vp(qJ, qJD, pkin, MDP)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(7,1),zeros(24,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RPRRP1_coriolisvecJ_fixb_mdp_slag_vp: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RPRRP1_coriolisvecJ_fixb_mdp_slag_vp: qJD has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [7 1]), ...
  'S5RPRRP1_coriolisvecJ_fixb_mdp_slag_vp: pkin has to be [7x1] (double)');
assert(isreal(MDP) && all(size(MDP) == [24 1]), ...
  'S5RPRRP1_coriolisvecJ_fixb_mdp_slag_vp: MDP has to be [24x1] (double)'); 

%% Symbolic Calculation
% From coriolisvec_joint_fixb_mdp_matlab.m
% OptimizationMode: 2
% StartTime: 2021-01-15 12:27:24
% EndTime: 2021-01-15 12:27:28
% DurationCPUTime: 1.16s
% Computational Cost: add. (1329->195), mult. (2793->248), div. (0->0), fcn. (1695->4), ass. (0->99)
t218 = cos(qJ(4));
t211 = qJD(3) + qJD(4);
t219 = cos(qJ(3));
t238 = t211 * t219;
t236 = t218 * t238;
t216 = sin(qJ(4));
t217 = sin(qJ(3));
t259 = qJD(1) * t217;
t247 = t216 * t259;
t264 = t211 * t247;
t160 = qJD(1) * t236 - t264;
t220 = -pkin(1) - pkin(6);
t199 = t220 * qJD(1) + qJD(2);
t258 = qJD(1) * t219;
t176 = -pkin(7) * t258 + t219 * t199;
t171 = qJD(3) * pkin(3) + t176;
t245 = pkin(7) * qJD(1) - t199;
t257 = qJD(3) * t217;
t172 = t245 * t257;
t256 = qJD(3) * t219;
t173 = t245 * t256;
t175 = -pkin(7) * t259 + t199 * t217;
t255 = qJD(4) * t216;
t230 = -(qJD(4) * t171 - t173) * t218 - t216 * t172 + t175 * t255;
t282 = -qJ(5) * t160 - t230;
t237 = MDP(22) * t211;
t281 = (t219 * MDP(12) - t217 * MDP(13)) * qJ(2) - t219 * t217 * MDP(7) + (t217 ^ 2 - t219 ^ 2) * MDP(8);
t248 = t218 * t258;
t184 = -t247 + t248;
t177 = t184 * qJ(5);
t168 = t216 * t175;
t243 = t218 * t171 - t168;
t146 = -t177 + t243;
t279 = MDP(12) * t217 + MDP(13) * t219;
t267 = t216 * t219;
t189 = t217 * t218 + t267;
t182 = t189 * qJD(1);
t277 = t184 ^ 2;
t276 = pkin(3) * t211;
t275 = pkin(7) - t220;
t161 = t211 * t189;
t159 = t161 * qJD(1);
t273 = qJ(5) * t159;
t271 = qJ(5) * t182;
t269 = t184 * t211;
t195 = pkin(3) * t259 + qJD(1) * qJ(2);
t268 = t195 * t184;
t169 = t218 * t175;
t145 = pkin(4) * t211 + t146;
t266 = t145 - t146;
t265 = t218 * t176 - t168;
t250 = pkin(3) * t258;
t192 = qJD(1) * qJD(2) + qJD(3) * t250;
t205 = t217 * pkin(3) + qJ(2);
t254 = qJD(4) * t218;
t246 = -pkin(4) * t182 - qJD(5);
t163 = -t246 + t195;
t251 = qJD(5) + t163;
t200 = pkin(3) * t256 + qJD(2);
t194 = t275 * t219;
t244 = -qJ(2) * MDP(6) - MDP(5);
t242 = t218 * t172 + t216 * t173;
t240 = -t176 * t216 - t169;
t155 = pkin(4) * t160 + t192;
t190 = -t216 * t217 + t218 * t219;
t235 = -t159 * t190 - t161 * t184;
t234 = -t171 * t216 - t169;
t193 = t275 * t217;
t233 = t193 * t218 + t194 * t216;
t181 = t182 ^ 2;
t232 = t182 * t184 * MDP(14) + (-t211 * t248 + t264 + t269) * MDP(17) + (-t181 + t277) * MDP(15);
t187 = t275 * t257;
t188 = qJD(3) * t194;
t231 = -t216 * t187 + t218 * t188 - t193 * t255 + t194 * t254;
t141 = -qJD(5) * t182 + t282;
t228 = t234 * qJD(4) + t242;
t225 = t228 + t273;
t142 = -qJD(5) * t184 + t225;
t147 = -t234 - t271;
t162 = -t216 * t257 - t217 * t255 + t236;
t229 = t141 * t189 + t142 * t190 - t145 * t161 + t147 * t162;
t227 = t233 * qJD(4) + t218 * t187 + t188 * t216;
t226 = t195 * t182 + t230;
t224 = (-t169 + (-t171 - t276) * t216) * qJD(4) + t242;
t223 = t251 * t182 - t282;
t222 = qJD(1) ^ 2;
t221 = qJD(3) ^ 2;
t206 = pkin(3) * t218 + pkin(4);
t196 = t254 * t276;
t174 = pkin(4) * t189 + t205;
t166 = pkin(4) * t184 + t250;
t158 = pkin(4) * t162 + t200;
t157 = -qJ(5) * t189 - t233;
t156 = -qJ(5) * t190 + t193 * t216 - t194 * t218;
t149 = -t177 + t265;
t148 = t240 + t271;
t144 = qJ(5) * t161 - qJD(5) * t190 + t227;
t143 = -qJ(5) * t162 - qJD(5) * t189 - t231;
t1 = [t235 * MDP(14) + (t159 * t189 - t160 * t190 + t161 * t182 - t162 * t184) * MDP(15) + (t205 * t160 + t195 * t162 + t200 * t182 + t192 * t189) * MDP(19) + (-t205 * t159 - t195 * t161 + t200 * t184 + t192 * t190) * MDP(20) + (t155 * t189 + t158 * t182 + t160 * t174 + t162 * t163) * MDP(21) + (t155 * t190 + t158 * t184 - t159 * t174 - t161 * t163) * MDP(22) + (-t143 * t182 - t144 * t184 + t156 * t159 - t157 * t160 - t229) * MDP(23) + (t141 * t157 + t142 * t156 + t143 * t147 + t144 * t145 + t155 * t174 + t158 * t163) * MDP(24) + (-t161 * MDP(16) - t162 * MDP(17) + t227 * MDP(19) + t231 * MDP(20) + t144 * MDP(21) - t143 * MDP(22)) * t211 + ((-MDP(13) * t220 - MDP(10)) * t219 + (-MDP(12) * t220 - MDP(9)) * t217) * t221 + (0.2e1 * (-t244 + t279) * qJD(2) + 0.2e1 * t281 * qJD(3)) * qJD(1); (-t160 * t189 - t162 * t182 - t235) * MDP(23) + (-qJD(1) * t163 + t229) * MDP(24) + t244 * t222 + (MDP(19) + MDP(21)) * (-qJD(1) * t182 - t161 * t211) + (MDP(20) + MDP(22)) * (-qJD(1) * t184 - t162 * t211) + t279 * (-t221 - t222); (-t182 * t250 - t240 * t211 + t224 - t268) * MDP(19) + (-t184 * t250 + t265 * t211 - t196 + t226) * MDP(20) + (-t148 * t211 - t166 * t182 - t251 * t184 + t224 + t273) * MDP(21) + (t149 * t211 - t166 * t184 - t196 + t223) * MDP(22) + (t159 * t206 + (t147 + t148) * t184 + (-t145 + t149) * t182 + (-t160 * t216 + (-t182 * t218 + t184 * t216) * qJD(4)) * pkin(3)) * MDP(23) + (t142 * t206 - t145 * t148 - t147 * t149 - t163 * t166 + (t141 * t216 + (-t145 * t216 + t147 * t218) * qJD(4)) * pkin(3)) * MDP(24) + t232 - t281 * t222; (-t234 * t211 + t228 - t268) * MDP(19) + (t243 * t211 + t226) * MDP(20) + (t147 * t211 + (-t163 + t246) * t184 + t225) * MDP(21) + (-pkin(4) * t277 + t146 * t211 + t223) * MDP(22) + (pkin(4) * t159 - t266 * t182) * MDP(23) + (t266 * t147 + (-t163 * t184 + t142) * pkin(4)) * MDP(24) + t232; (-t264 + t269) * MDP(21) - t182 * t237 + (-t181 - t277) * MDP(23) + (t145 * t184 + t147 * t182 + t155) * MDP(24) + (-t237 * t267 + (MDP(21) * t238 - t217 * t237) * t218) * qJD(1);];
tauc = t1;
