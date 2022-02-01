% Calculate Coriolis joint torque vector for
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
% Datum: 2022-01-23 09:00
% Revision: fd3771346c4aea32fdeb66112c511235427c26a7 (2022-01-20)
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
% StartTime: 2022-01-23 08:59:49
% EndTime: 2022-01-23 08:59:52
% DurationCPUTime: 1.71s
% Computational Cost: add. (816->199), mult. (2353->331), div. (0->0), fcn. (1800->8), ass. (0->106)
t213 = cos(pkin(8));
t212 = cos(pkin(9));
t211 = sin(pkin(7));
t255 = qJD(1) * t211;
t236 = t212 * t255;
t209 = sin(pkin(9));
t214 = cos(pkin(7));
t254 = qJD(1) * t214;
t237 = t209 * t254;
t178 = t213 * t236 - t237;
t215 = sin(qJ(5));
t216 = cos(qJ(5));
t210 = sin(pkin(8));
t256 = qJD(1) * t210;
t242 = t216 * t256;
t272 = -t178 * t215 + t211 * t242;
t261 = t211 * t213;
t185 = t209 * t261 + t212 * t214;
t175 = t185 * qJD(1);
t171 = qJD(5) + t175;
t271 = qJD(5) - t171;
t252 = qJD(3) * t211;
t239 = t213 * t252;
t253 = qJD(2) * t214;
t187 = t210 * t253 + t239;
t181 = t187 * qJD(1);
t268 = t181 * t213;
t196 = qJ(2) * t255 + qJD(3);
t267 = t196 * t211;
t207 = t211 ^ 2;
t217 = qJD(1) ^ 2;
t266 = t207 * t217;
t265 = t209 * t210;
t264 = t210 * t211;
t263 = t210 * t215;
t262 = t210 * t216;
t260 = t213 * t214;
t191 = -pkin(2) * t214 - qJ(3) * t211 - pkin(1);
t184 = t191 * qJD(1) + qJD(2);
t245 = qJ(2) * t254;
t161 = t210 * t184 + t213 * t245;
t151 = -qJ(4) * t254 + t161;
t229 = pkin(3) * t210 - qJ(4) * t213;
t166 = t229 * t255 + t196;
t141 = t212 * t151 + t209 * t166;
t258 = qJ(2) * t260 + t210 * t191;
t165 = -qJ(4) * t214 + t258;
t172 = (qJ(2) + t229) * t211;
t259 = t212 * t165 + t209 * t172;
t208 = t214 ^ 2;
t257 = t207 + t208;
t251 = qJD(5) * t171;
t249 = qJD(1) * qJD(2);
t248 = 0.2e1 * qJD(2) * t207;
t206 = t210 ^ 2;
t247 = t206 * t266;
t244 = t210 * t255;
t243 = t215 * t256;
t241 = t213 * t255;
t240 = t210 * t252;
t238 = t209 * t251;
t199 = t213 * t253;
t235 = qJ(2) * t249;
t193 = qJD(1) * t199;
t220 = -qJD(4) * t214 - t240;
t167 = t220 * qJD(1) + t193;
t223 = (-qJD(4) * t213 + qJD(2)) * t211;
t219 = qJD(1) * t223;
t146 = t212 * t167 + t209 * t219;
t233 = -t215 * t146 + t216 * t181;
t160 = t184 * t213 - t210 * t245;
t232 = -t210 * t214 * qJ(2) + t191 * t213;
t230 = t214 * pkin(3) - t232;
t152 = t272 * qJD(5);
t150 = pkin(3) * t254 + qJD(4) - t160;
t137 = pkin(4) * t175 - pkin(6) * t178 + t150;
t139 = pkin(6) * t244 + t141;
t228 = t137 * t216 - t139 * t215;
t227 = -t137 * t215 - t139 * t216;
t226 = t216 * t146 + t215 * t181;
t140 = -t151 * t209 + t166 * t212;
t225 = -t165 * t209 + t172 * t212;
t182 = -qJD(1) * t240 + t193;
t224 = t182 * t210 - t268;
t186 = -t209 * t214 + t212 * t261;
t222 = -t186 * t215 + t211 * t262;
t163 = t186 * t216 + t211 * t263;
t221 = t212 * t263 + t213 * t216;
t156 = t178 * t216 + t211 * t243;
t218 = (-t212 * t262 + t213 * t215) * t171;
t202 = t207 * t235;
t188 = t199 - t240;
t180 = (t209 * t211 + t212 * t260) * qJD(1);
t177 = t213 * t237 - t236;
t170 = t199 + t220;
t158 = t163 * qJD(5);
t157 = t222 * qJD(5);
t153 = t156 * qJD(5);
t148 = t212 * t170 + t209 * t223;
t147 = t170 * t209 - t212 * t223;
t145 = t167 * t209 - t212 * t219;
t144 = pkin(6) * t264 + t259;
t143 = -pkin(4) * t264 - t225;
t142 = pkin(4) * t185 - pkin(6) * t186 + t230;
t138 = -pkin(4) * t244 - t140;
t1 = [0.2e1 * t257 * MDP(6) * t249 + 0.2e1 * (t208 * t235 + t202) * MDP(7) + (t181 * t214 + (t187 * t214 + t210 * t248) * qJD(1)) * MDP(8) + (t182 * t214 + (t188 * t214 + t213 * t248) * qJD(1)) * MDP(9) + ((t187 * t213 - t188 * t210) * qJD(1) - t224) * t211 * MDP(10) + (qJD(2) * t267 - t160 * t187 + t161 * t188 - t181 * t232 + t182 * t258 + t202) * MDP(11) + (t175 * t187 + t181 * t185 + (-qJD(1) * t147 - t145) * t264) * MDP(12) + (t178 * t187 + t181 * t186 + (-qJD(1) * t148 - t146) * t264) * MDP(13) + (t145 * t186 - t146 * t185 + t147 * t178 - t148 * t175) * MDP(14) + (-t140 * t147 + t141 * t148 - t145 * t225 + t146 * t259 + t150 * t187 + t181 * t230) * MDP(15) + (t152 * t163 + t156 * t157) * MDP(16) + (t152 * t222 - t153 * t163 - t156 * t158 + t157 * t272) * MDP(17) + (t152 * t185 + t157 * t171) * MDP(18) + (-t153 * t185 - t158 * t171) * MDP(19) + ((-t148 * t215 + t187 * t216) * t171 + t233 * t185 - t147 * t272 + t143 * t153 - t145 * t222 + t138 * t158 + ((-t142 * t215 - t144 * t216) * t171 + t227 * t185) * qJD(5)) * MDP(21) + (-(t148 * t216 + t187 * t215) * t171 - t226 * t185 + t147 * t156 + t143 * t152 + t145 * t163 + t138 * t157 + (-(t142 * t216 - t144 * t215) * t171 - t228 * t185) * qJD(5)) * MDP(22); ((-t267 + (t160 * t210 - t161 * t213) * t214) * qJD(1) + t224) * MDP(11) + (t175 * t180 - t177 * t178) * MDP(14) + (t140 * t177 - t141 * t180 - t268 + (t145 * t209 + t146 * t212 - t150 * t254) * t210) * MDP(15) + (t153 * t265 - (-t180 * t215 + t214 * t242) * t171 + t177 * t272 + qJD(5) * t218) * MDP(21) + (t152 * t265 + (t180 * t216 + t214 * t243) * t171 - t177 * t156 + t221 * t251) * MDP(22) + ((-t175 * t214 + t177 * t211) * MDP(12) + (-t178 * t214 + t180 * t211) * MDP(13)) * t256 - (qJ(2) * MDP(7) + t210 * MDP(8) + t213 * MDP(9) + MDP(6)) * t257 * t217; (-t213 ^ 2 - t206) * MDP(10) * t266 + (t160 * t213 + t161 * t210 + qJD(2)) * MDP(11) * t255 + (-t175 * t241 - t209 * t247) * MDP(12) + (-t178 * t241 - t212 * t247) * MDP(13) + (-t175 * t212 + t178 * t209) * MDP(14) * t244 + (-t145 * t212 + t146 * t209 + (-t150 * t213 + (-t140 * t209 + t141 * t212) * t210) * t255) * MDP(15) + (-t216 * t238 - t212 * t153 + (-t221 * t171 - t265 * t272) * t255) * MDP(21) + (t215 * t238 - t212 * t152 + (t156 * t265 + t218) * t255) * MDP(22) + (-t213 * MDP(8) + t210 * MDP(9)) * t211 * t214 * t217; (-t175 * MDP(14) + t141 * MDP(15)) * t175 + (-t178 * MDP(14) + t140 * MDP(15) + MDP(21) * t272 - t156 * MDP(22)) * t178 + (MDP(15) * t239 + (MDP(15) * t253 + (MDP(12) * t178 - MDP(13) * t175) * t211) * t210) * qJD(1) - (MDP(21) * t215 + MDP(22) * t216) * t171 ^ 2; -t272 ^ 2 * MDP(17) + (-t171 * t272 + t152) * MDP(18) + (t271 * t227 + t233) * MDP(21) + (-t138 * t272 - t228 * t271 - t226) * MDP(22) + (-MDP(16) * t272 + MDP(17) * t156 - MDP(19) * t271 - t138 * MDP(21)) * t156;];
tauc = t1;
