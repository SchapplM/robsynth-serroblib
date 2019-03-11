% Calculate minimal parameter regressor of Coriolis joint torque vector for
% S6RPPPRR2
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% qJD [6x1]
%   Generalized joint velocities
% pkin [10x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d5,d6,theta2,theta4]';
% MDP [25x1]
%   Minimal dynamic parameter vector (fixed base model)
%   see S6RPPPRR2_convert_par2_MPV_fixb.m
% 
% Output:
% tauc [6x1]
%   joint torques required to compensate Coriolis and centrifugal load

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-09 01:32
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauc = S6RPPPRR2_coriolisvecJ_fixb_mdp_slag_vp(qJ, qJD, pkin, MDP)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(10,1),zeros(25,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPPPRR2_coriolisvecJ_fixb_mdp_slag_vp: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RPPPRR2_coriolisvecJ_fixb_mdp_slag_vp: qJD has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S6RPPPRR2_coriolisvecJ_fixb_mdp_slag_vp: pkin has to be [10x1] (double)');
assert(isreal(MDP) && all(size(MDP) == [25 1]), ...
  'S6RPPPRR2_coriolisvecJ_fixb_mdp_slag_vp: MDP has to be [25x1] (double)'); 

%% Symbolic Calculation
% From coriolisvec_joint_fixb_mdp_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 01:32:05
% EndTime: 2019-03-09 01:32:09
% DurationCPUTime: 1.80s
% Computational Cost: add. (1243->211), mult. (2717->295), div. (0->0), fcn. (1930->8), ass. (0->100)
t271 = 2 * qJD(3);
t190 = sin(pkin(9)) * pkin(1) + qJ(3);
t267 = qJD(1) * t190;
t270 = t267 * MDP(7);
t209 = cos(qJ(6));
t235 = t209 * qJD(5);
t203 = sin(pkin(10));
t208 = sin(qJ(5));
t241 = qJD(1) * t208;
t227 = t203 * t241;
t205 = cos(pkin(10));
t210 = cos(qJ(5));
t240 = qJD(1) * t210;
t229 = t205 * t240;
t175 = -t227 + t229;
t207 = sin(qJ(6));
t250 = t175 * t207;
t158 = -t235 + t250;
t180 = t203 * t210 + t205 * t208;
t215 = qJD(1) * t180;
t265 = qJD(6) + t215;
t269 = t158 * t265;
t160 = qJD(5) * t207 + t175 * t209;
t268 = t160 * t265;
t225 = t209 * t265;
t184 = qJD(5) * t227;
t238 = qJD(5) * t210;
t228 = t205 * t238;
t166 = qJD(1) * t228 - t184;
t247 = t207 * t166;
t266 = -t225 * t265 - t247;
t244 = t203 ^ 2 + t205 ^ 2;
t263 = MDP(10) * t244;
t189 = -cos(pkin(9)) * pkin(1) - pkin(2) - qJ(4);
t262 = qJD(1) * t189;
t261 = t203 * MDP(8) + t205 * MDP(9) + MDP(6);
t178 = qJD(3) + t262;
t161 = -qJD(2) * t203 + t205 * t178;
t156 = -pkin(7) * qJD(1) * t205 + t161;
t162 = t205 * qJD(2) + t203 * t178;
t242 = qJD(1) * t203;
t157 = -pkin(7) * t242 + t162;
t140 = t156 * t208 + t157 * t210;
t179 = t203 * t208 - t210 * t205;
t214 = t179 * qJD(4);
t136 = -qJD(1) * t214 + qJD(5) * t140;
t260 = t265 * (pkin(5) * t175 + t265 * pkin(8)) + t136;
t139 = t156 * t210 - t157 * t208;
t216 = t180 * qJD(4);
t135 = -qJD(1) * t216 + qJD(5) * t139;
t137 = -qJD(5) * pkin(5) - t139;
t258 = -pkin(7) + t189;
t171 = t258 * t203;
t172 = t258 * t205;
t151 = t171 * t208 - t172 * t210;
t142 = -qJD(5) * t151 - t216;
t181 = qJD(4) + t267;
t170 = pkin(4) * t242 + t181;
t144 = pkin(5) * t215 - pkin(8) * t175 + t170;
t182 = t203 * pkin(4) + t190;
t150 = pkin(5) * t180 + pkin(8) * t179 + t182;
t152 = t171 * t210 + t172 * t208;
t239 = qJD(5) * t208;
t176 = -t203 * t238 - t205 * t239;
t259 = -(qJD(6) * t150 + t142) * t265 - t152 * t166 - (qJD(6) * t144 + t135) * t180 - t136 * t179 + t137 * t176;
t257 = t137 * t179;
t237 = qJD(6) * t207;
t165 = qJD(5) * t215;
t245 = qJD(6) * t235 - t209 * t165;
t145 = -t175 * t237 + t245;
t256 = t145 * t179;
t255 = t145 * t207;
t254 = t150 * t166;
t253 = t158 * t175;
t252 = t160 * t175;
t251 = t165 * t207;
t164 = t209 * t166;
t177 = -t203 * t239 + t228;
t246 = t145 * t180 + t160 * t177;
t236 = qJD(6) * t209;
t233 = qJD(1) * qJD(4);
t201 = qJD(3) * qJD(1);
t231 = t179 * t247;
t230 = t179 * t164;
t138 = qJD(5) * pkin(8) + t140;
t134 = t138 * t209 + t144 * t207;
t222 = t138 * t207 - t144 * t209;
t146 = t160 * qJD(6) - t251;
t221 = -t146 * t180 - t158 * t177;
t220 = t161 * t205 + t162 * t203;
t219 = t164 + (-t207 * t215 - t237) * t265;
t218 = -t176 * t207 + t179 * t236;
t217 = t176 * t209 + t179 * t237;
t212 = -pkin(8) * t166 + (t137 + t139) * t265;
t211 = qJD(1) ^ 2;
t153 = pkin(5) * t177 - pkin(8) * t176 + qJD(3);
t148 = pkin(5) * t166 + pkin(8) * t165 + t201;
t147 = t209 * t148;
t143 = qJD(5) * t152 - t214;
t1 = [t270 * t271 + 0.2e1 * t233 * t263 + ((t181 + t267) * qJD(3) + (-t244 * t262 - t220) * qJD(4)) * MDP(11) + (t165 * t179 + t175 * t176) * MDP(12) + (t165 * t180 + t166 * t179 - t175 * t177 - t176 * t215) * MDP(13) + (t166 * t182 + t170 * t177 + t215 * t271) * MDP(17) + (-t165 * t182 + t170 * t176 + (-qJD(1) * t179 + t175) * qJD(3)) * MDP(18) + (t217 * t160 - t209 * t256) * MDP(19) + ((-t158 * t209 - t160 * t207) * t176 + (t255 + t146 * t209 + (-t158 * t207 + t160 * t209) * qJD(6)) * t179) * MDP(20) + (t217 * t265 - t230 + t246) * MDP(21) + (t218 * t265 + t221 + t231) * MDP(22) + (t166 * t180 + t177 * t265) * MDP(23) + (-t222 * t177 + t143 * t158 + t151 * t146 + t147 * t180 + (t254 + t153 * t265 + (-t138 * t180 - t152 * t265 - t257) * qJD(6)) * t209 + t259 * t207) * MDP(24) + (-t134 * t177 + t143 * t160 + t151 * t145 + (-(-qJD(6) * t152 + t153) * t265 - t254 - (-qJD(6) * t138 + t148) * t180 + qJD(6) * t257) * t207 + t259 * t209) * MDP(25) + (MDP(14) * t176 - MDP(15) * t177 - MDP(17) * t143 - MDP(18) * t142) * qJD(5) + 0.2e1 * t261 * t201; (-t221 + t231) * MDP(24) + (t230 + t246) * MDP(25) + (MDP(24) * t218 - MDP(25) * t217) * t265 + (-t177 * MDP(17) - t176 * MDP(18)) * qJD(5); (t146 * t179 - t158 * t176 - t180 * t247) * MDP(24) + (-t160 * t176 - t180 * t164 + t256) * MDP(25) + ((-t177 * t207 - t180 * t236) * MDP(24) + (-t177 * t209 + t180 * t237) * MDP(25)) * t265 + (MDP(17) * t176 - MDP(18) * t177) * qJD(5) - t261 * t211 + (-t270 + (-t244 * qJD(4) - t181) * MDP(11) - t215 * MDP(17) - t175 * MDP(18) + (-MDP(24) * t209 + MDP(25) * t207) * t265) * qJD(1); (qJD(1) * t220 + t201) * MDP(11) - t184 * MDP(17) + (t219 - t253) * MDP(24) + (-t252 + t266) * MDP(25) - t211 * t263 + ((t175 + t229) * MDP(17) + (-t203 * t240 - t205 * t241 - t215) * MDP(18)) * qJD(5); t175 ^ 2 * MDP(13) + (t184 + (t175 - t229) * qJD(5)) * MDP(15) + (-t170 * t175 + t179 * t233) * MDP(17) + (t160 * t225 + t255) * MDP(19) + ((t145 - t269) * t209 + (-t146 - t268) * t207) * MDP(20) + (-t252 - t266) * MDP(21) + (t219 + t253) * MDP(22) - t265 * t175 * MDP(23) + (-pkin(5) * t146 - t140 * t158 + t175 * t222 + t212 * t207 - t260 * t209) * MDP(24) + (-pkin(5) * t145 + t134 * t175 - t140 * t160 + t260 * t207 + t212 * t209) * MDP(25) + (t175 * MDP(12) + (qJD(4) + t170) * MDP(18) - MDP(13) * t215) * t215; t160 * t158 * MDP(19) + (-t158 ^ 2 + t160 ^ 2) * MDP(20) + (t245 + t269) * MDP(21) + (t251 + t268) * MDP(22) + t166 * MDP(23) + (t134 * t265 - t135 * t207 - t137 * t160 + t147) * MDP(24) + (-t135 * t209 + t137 * t158 - t148 * t207 - t222 * t265) * MDP(25) + (-MDP(21) * t250 - MDP(22) * t160 - MDP(24) * t134 + MDP(25) * t222) * qJD(6);];
tauc  = t1;
