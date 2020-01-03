% Calculate minimal parameter regressor of Coriolis joint torque vector for
% S5RPPPR2
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% qJD [5x1]
%   Generalized joint velocities
% pkin [9x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d1,d5,theta2,theta3,theta4]';
% MDP [22x1]
%   Minimal dynamic parameter vector (fixed base model)
%   see S5RPPPR2_convert_par2_MPV_fixb.m
% 
% Output:
% tauc [5x1]
%   joint torques required to compensate Coriolis and centrifugal load

% Quelle: HybrDyn-Toolbox
% Datum: 2020-01-03 11:23
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauc = S5RPPPR2_coriolisvecJ_fixb_mdp_slag_vp(qJ, qJD, pkin, MDP)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(9,1),zeros(22,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RPPPR2_coriolisvecJ_fixb_mdp_slag_vp: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RPPPR2_coriolisvecJ_fixb_mdp_slag_vp: qJD has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S5RPPPR2_coriolisvecJ_fixb_mdp_slag_vp: pkin has to be [9x1] (double)');
assert(isreal(MDP) && all(size(MDP) == [22 1]), ...
  'S5RPPPR2_coriolisvecJ_fixb_mdp_slag_vp: MDP has to be [22x1] (double)'); 

%% Symbolic Calculation
% From coriolisvec_joint_fixb_mdp_matlab.m
% OptimizationMode: 2
% StartTime: 2020-01-03 11:23:02
% EndTime: 2020-01-03 11:23:08
% DurationCPUTime: 1.66s
% Computational Cost: add. (816->199), mult. (2353->330), div. (0->0), fcn. (1800->8), ass. (0->107)
t213 = cos(pkin(8));
t212 = cos(pkin(9));
t211 = sin(pkin(7));
t255 = qJD(1) * t211;
t237 = t212 * t255;
t209 = sin(pkin(9));
t214 = cos(pkin(7));
t254 = qJD(1) * t214;
t239 = t209 * t254;
t178 = t213 * t237 - t239;
t215 = sin(qJ(5));
t216 = cos(qJ(5));
t210 = sin(pkin(8));
t256 = qJD(1) * t210;
t242 = t216 * t256;
t273 = -t215 * t178 + t211 * t242;
t261 = t214 * t212;
t262 = t211 * t213;
t185 = t209 * t262 + t261;
t175 = t185 * qJD(1);
t171 = qJD(5) + t175;
t272 = qJD(5) - t171;
t270 = qJ(2) * t214;
t252 = qJD(3) * t211;
t236 = t213 * t252;
t253 = qJD(2) * t214;
t187 = t210 * t253 + t236;
t181 = t187 * qJD(1);
t269 = t181 * t213;
t196 = qJ(2) * t255 + qJD(3);
t268 = t196 * t211;
t207 = t211 ^ 2;
t217 = qJD(1) ^ 2;
t267 = t207 * t217;
t266 = t209 * t210;
t265 = t210 * t211;
t264 = t210 * t215;
t263 = t210 * t216;
t191 = -t214 * pkin(2) - t211 * qJ(3) - pkin(1);
t184 = t191 * qJD(1) + qJD(2);
t245 = qJ(2) * t254;
t161 = t210 * t184 + t213 * t245;
t151 = -qJ(4) * t254 + t161;
t229 = pkin(3) * t210 - qJ(4) * t213;
t166 = t229 * t255 + t196;
t141 = t212 * t151 + t209 * t166;
t258 = t210 * t191 + t213 * t270;
t165 = -t214 * qJ(4) + t258;
t172 = (qJ(2) + t229) * t211;
t259 = t212 * t165 + t209 * t172;
t208 = t214 ^ 2;
t257 = t207 + t208;
t251 = qJD(5) * t171;
t249 = qJD(1) * qJD(2);
t248 = 0.2e1 * qJD(2) * t207;
t206 = t210 ^ 2;
t247 = t206 * t267;
t244 = t210 * t255;
t243 = t215 * t256;
t241 = t213 * t255;
t240 = t209 * t251;
t238 = t210 * t252;
t199 = t213 * t253;
t235 = qJ(2) * t249;
t193 = qJD(1) * t199;
t220 = -t214 * qJD(4) - t238;
t167 = t220 * qJD(1) + t193;
t223 = (-qJD(4) * t213 + qJD(2)) * t211;
t219 = qJD(1) * t223;
t146 = t212 * t167 + t209 * t219;
t233 = -t215 * t146 + t216 * t181;
t160 = t213 * t184 - t210 * t245;
t232 = t213 * t191 - t210 * t270;
t230 = t214 * pkin(3) - t232;
t152 = t273 * qJD(5);
t150 = pkin(3) * t254 + qJD(4) - t160;
t137 = t175 * pkin(4) - t178 * pkin(6) + t150;
t139 = pkin(6) * t244 + t141;
t228 = t216 * t137 - t215 * t139;
t227 = -t215 * t137 - t216 * t139;
t226 = t216 * t146 + t215 * t181;
t140 = -t209 * t151 + t212 * t166;
t225 = -t209 * t165 + t212 * t172;
t182 = -qJD(1) * t238 + t193;
t224 = t182 * t210 - t269;
t186 = -t214 * t209 + t212 * t262;
t222 = -t186 * t215 + t211 * t263;
t163 = t186 * t216 + t211 * t264;
t221 = t212 * t264 + t213 * t216;
t156 = t216 * t178 + t211 * t243;
t218 = (-t212 * t263 + t213 * t215) * t171;
t202 = t207 * t235;
t188 = t199 - t238;
t180 = (t209 * t211 + t213 * t261) * qJD(1);
t177 = t213 * t239 - t237;
t170 = t199 + t220;
t158 = t163 * qJD(5);
t157 = t222 * qJD(5);
t153 = t156 * qJD(5);
t148 = t212 * t170 + t209 * t223;
t147 = t209 * t170 - t212 * t223;
t145 = t209 * t167 - t212 * t219;
t144 = pkin(6) * t265 + t259;
t143 = -pkin(4) * t265 - t225;
t142 = t185 * pkin(4) - t186 * pkin(6) + t230;
t138 = -pkin(4) * t244 - t140;
t1 = [0.2e1 * t257 * MDP(6) * t249 + 0.2e1 * (t208 * t235 + t202) * MDP(7) + (t181 * t214 + (t187 * t214 + t210 * t248) * qJD(1)) * MDP(8) + (t182 * t214 + (t188 * t214 + t213 * t248) * qJD(1)) * MDP(9) + ((t187 * t213 - t188 * t210) * qJD(1) - t224) * t211 * MDP(10) + (qJD(2) * t268 - t160 * t187 + t161 * t188 - t181 * t232 + t182 * t258 + t202) * MDP(11) + (t187 * t175 + t181 * t185 + (-qJD(1) * t147 - t145) * t265) * MDP(12) + (t187 * t178 + t181 * t186 + (-qJD(1) * t148 - t146) * t265) * MDP(13) + (t145 * t186 - t146 * t185 + t147 * t178 - t148 * t175) * MDP(14) + (-t140 * t147 + t141 * t148 - t145 * t225 + t146 * t259 + t150 * t187 + t181 * t230) * MDP(15) + (t152 * t163 + t156 * t157) * MDP(16) + (t152 * t222 - t163 * t153 - t156 * t158 + t157 * t273) * MDP(17) + (t152 * t185 + t157 * t171) * MDP(18) + (-t153 * t185 - t158 * t171) * MDP(19) + ((-t215 * t148 + t216 * t187) * t171 + t233 * t185 - t147 * t273 + t143 * t153 - t145 * t222 + t138 * t158 + ((-t142 * t215 - t144 * t216) * t171 + t227 * t185) * qJD(5)) * MDP(21) + (-(t216 * t148 + t215 * t187) * t171 - t226 * t185 + t147 * t156 + t143 * t152 + t145 * t163 + t138 * t157 + (-(t142 * t216 - t144 * t215) * t171 - t228 * t185) * qJD(5)) * MDP(22); ((-t268 + (t160 * t210 - t161 * t213) * t214) * qJD(1) + t224) * MDP(11) + (t180 * t175 - t177 * t178) * MDP(14) + (t140 * t177 - t141 * t180 - t269 + (t145 * t209 + t146 * t212 - t150 * t254) * t210) * MDP(15) + (t153 * t266 - (-t215 * t180 + t214 * t242) * t171 + t177 * t273 + qJD(5) * t218) * MDP(21) + (t152 * t266 + (t216 * t180 + t214 * t243) * t171 - t177 * t156 + t221 * t251) * MDP(22) + ((-t175 * t214 + t177 * t211) * MDP(12) + (-t178 * t214 + t180 * t211) * MDP(13)) * t256 - (qJ(2) * MDP(7) + t210 * MDP(8) + t213 * MDP(9) + MDP(6)) * t257 * t217; (-t213 ^ 2 - t206) * MDP(10) * t267 + (t160 * t213 + t161 * t210 + qJD(2)) * MDP(11) * t255 + (-t175 * t241 - t209 * t247) * MDP(12) + (-t178 * t241 - t212 * t247) * MDP(13) + (-t175 * t212 + t178 * t209) * MDP(14) * t244 + (-t145 * t212 + t146 * t209 + (-t150 * t213 + (-t140 * t209 + t141 * t212) * t210) * t255) * MDP(15) + (-t216 * t240 - t212 * t153 + (-t221 * t171 - t266 * t273) * t255) * MDP(21) + (t215 * t240 - t212 * t152 + (t156 * t266 + t218) * t255) * MDP(22) + (-t213 * MDP(8) + t210 * MDP(9)) * t211 * t214 * t217; (-MDP(14) * t175 + t141 * MDP(15)) * t175 + (-MDP(14) * t178 + t140 * MDP(15) + MDP(21) * t273 - t156 * MDP(22)) * t178 + (MDP(15) * t236 + (MDP(15) * t253 + (MDP(12) * t178 - MDP(13) * t175) * t211) * t210) * qJD(1) - (MDP(21) * t215 + MDP(22) * t216) * t171 ^ 2; -t273 ^ 2 * MDP(17) + (-t171 * t273 + t152) * MDP(18) + (t272 * t227 + t233) * MDP(21) + (-t138 * t273 - t228 * t272 - t226) * MDP(22) + (-MDP(16) * t273 + MDP(17) * t156 - MDP(19) * t272 - t138 * MDP(21)) * t156;];
tauc = t1;
