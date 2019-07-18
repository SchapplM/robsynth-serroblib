% Calculate minimal parameter regressor of Coriolis joint torque vector for
% S5PRRRR1
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% qJD [5x1]
%   Generalized joint velocities
% pkin [2x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a4]';
% MDP [25x1]
%   Minimal dynamic parameter vector (fixed base model)
%   see S5PRRRR1_convert_par2_MPV_fixb.m
% 
% Output:
% tauc [5x1]
%   joint torques required to compensate Coriolis and centrifugal load

% Quelle: HybrDyn-Toolbox
% Datum: 2019-07-18 13:29
% Revision: 08c8d617a845f5dd194efdf9aca2774760f7818f (2019-07-16)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauc = S5PRRRR1_coriolisvecJ_fixb_mdp_slag_vp(qJ, qJD, pkin, MDP)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(2,1),zeros(25,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5PRRRR1_coriolisvecJ_fixb_mdp_slag_vp: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5PRRRR1_coriolisvecJ_fixb_mdp_slag_vp: qJD has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [2 1]), ...
  'S5PRRRR1_coriolisvecJ_fixb_mdp_slag_vp: pkin has to be [2x1] (double)');
assert(isreal(MDP) && all(size(MDP) == [25 1]), ...
  'S5PRRRR1_coriolisvecJ_fixb_mdp_slag_vp: MDP has to be [25x1] (double)'); 

%% Symbolic Calculation
% From coriolisvec_joint_fixb_mdp_matlab.m
% OptimizationMode: 2
% StartTime: 2019-07-18 13:28:54
% EndTime: 2019-07-18 13:28:59
% DurationCPUTime: 2.13s
% Computational Cost: add. (995->206), mult. (2786->318), div. (0->0), fcn. (2200->8), ass. (0->97)
t197 = qJD(3) + qJD(4);
t205 = cos(qJ(3));
t271 = cos(qJ(4));
t235 = qJD(2) * t271;
t194 = t205 * t235;
t201 = sin(qJ(4));
t202 = sin(qJ(3));
t249 = qJD(2) * t202;
t277 = -t201 * t249 + t194;
t155 = t277 * t197;
t258 = t201 * t205;
t182 = t202 * t271 + t258;
t275 = t197 * t182;
t156 = t275 * qJD(2);
t276 = -0.2e1 * qJD(2) * qJD(3);
t175 = qJD(5) - t277;
t245 = qJD(3) * t202;
t203 = sin(qJ(2));
t252 = qJD(1) * t203;
t266 = t156 * t205;
t274 = (t175 * t245 - t266) * pkin(2) - t175 * t252;
t273 = -t202 * t205 * MDP(5) + (t202 ^ 2 - t205 ^ 2) * MDP(6);
t206 = cos(qJ(2));
t246 = qJD(2) * t206;
t171 = (-qJD(3) * t203 * t205 - t202 * t246) * qJD(1);
t239 = t202 * t252;
t186 = qJD(3) * pkin(2) - t239;
t228 = qJD(3) * t239;
t251 = qJD(1) * t205;
t237 = t203 * t251;
t229 = t201 * t237;
t234 = qJD(4) * t271;
t250 = qJD(1) * t206;
t143 = -qJD(4) * t229 + t201 * t171 + t186 * t234 + t194 * t250 - t271 * t228;
t254 = -t271 * t171 - t201 * t228;
t260 = t186 * t201;
t144 = (t201 * t246 + t203 * t234) * t251 + qJD(4) * t260 + t254;
t218 = -t201 * t202 + t205 * t271;
t157 = t197 * t218;
t164 = -t186 * t271 + t229;
t214 = t206 * t218;
t170 = qJD(1) * t214;
t247 = qJD(2) * t205;
t187 = -pkin(2) * t247 - t250;
t272 = (qJD(5) * t187 + t143) * t218 + t144 * t182 + t164 * t157 + (pkin(2) * qJD(5) * t205 + t170) * t175;
t178 = -t201 * t247 - t202 * t235;
t200 = sin(qJ(5));
t244 = qJD(5) * t200;
t204 = cos(qJ(5));
t243 = qJD(5) * t204;
t255 = t155 * t204 + t197 * t243;
t145 = t178 * t244 + t255;
t270 = t145 * t200;
t269 = t155 * t200;
t268 = t156 * t200;
t267 = t156 * t204;
t262 = t178 * t200;
t161 = -t197 * t204 - t262;
t265 = t161 * t175;
t222 = t178 * t204 - t197 * t200;
t264 = t222 * t175;
t263 = t164 * t182;
t261 = t182 * t204;
t259 = t187 * t178;
t207 = qJD(3) ^ 2;
t257 = t202 * t207;
t256 = t205 * t207;
t248 = qJD(2) * t203;
t233 = t164 * (-t175 - t277);
t232 = t175 * t204;
t230 = t206 * t276;
t188 = t271 * t237;
t165 = t188 + t260;
t153 = t165 * t204 + t187 * t200;
t227 = t144 * t200 - t153 * t178 + t164 * t243;
t169 = t182 * t252;
t226 = -t164 * t277 - t169 * t175;
t225 = t165 * t200 - t187 * t204;
t173 = t218 * t203;
t224 = t173 * t204 - t200 * t206;
t223 = -t173 * t200 - t204 * t206;
t221 = qJD(5) * t201 + t249;
t219 = -t144 * t204 + t164 * t244 - t178 * t225;
t216 = t157 * t204 - t182 * t244;
t215 = t206 * t182;
t146 = -qJD(5) * t222 + t269;
t211 = ((t145 - t265) * t204 + (-t146 + t264) * t200) * MDP(20) + (-t222 * t232 + t270) * MDP(19) + (-t175 ^ 2 * t200 - t161 * t178 + t267) * MDP(22) + (t175 * t232 - t178 * t222 + t268) * MDP(21) + (-t178 * t197 - t156) * MDP(15) + (t178 ^ 2 - t277 ^ 2) * MDP(13) + (MDP(12) * t277 + MDP(23) * t175) * t178;
t209 = -t187 * t277 - t143;
t208 = qJD(2) ^ 2;
t180 = (pkin(2) * t245 + t252) * qJD(2);
t174 = t204 * t180;
t172 = t182 * t203;
t168 = qJD(1) * t215;
t167 = t201 * t239 - t188;
t148 = qJD(2) * t215 + t157 * t203;
t147 = qJD(2) * t214 - t203 * t275;
t1 = [(t202 * t230 - t203 * t256) * MDP(10) + (t203 * t257 + t205 * t230) * MDP(11) + (-t148 * t197 - t206 * t156 - t248 * t277) * MDP(17) + (-t147 * t197 - t155 * t206 - t178 * t248) * MDP(18) + ((-qJD(5) * t224 - t147 * t200 + t204 * t248) * t175 + t223 * t156 + t148 * t161 + t172 * t146) * MDP(24) + (-(qJD(5) * t223 + t147 * t204 + t200 * t248) * t175 - t224 * t156 - t148 * t222 + t172 * t145) * MDP(25) + (-t206 * MDP(4) + (-MDP(10) * t205 + MDP(11) * t202 - MDP(3)) * t203) * t208; MDP(7) * t256 - MDP(8) * t257 + (t155 * t182 - t157 * t178) * MDP(12) + (t155 * t218 - t156 * t182 + t157 * t277 + t178 * t275) * MDP(13) + (t277 * t252 + t275 * t187 - t180 * t218 + (-t245 * t277 - t266) * pkin(2)) * MDP(17) + (t178 * t252 + t157 * t187 + t180 * t182 + (-t155 * t205 - t178 * t245) * pkin(2)) * MDP(18) + (t145 * t261 - t216 * t222) * MDP(19) + ((-t161 * t204 + t200 * t222) * t157 + (-t270 - t146 * t204 + (t161 * t200 + t204 * t222) * qJD(5)) * t182) * MDP(20) + (-t145 * t218 + t156 * t261 + t175 * t216 - t222 * t275) * MDP(21) + (-t182 * t268 + t146 * t218 - t275 * t161 + (-t157 * t200 - t182 * t243) * t175) * MDP(22) + (-t156 * t218 + t175 * t275) * MDP(23) + (-t225 * t275 - t168 * t161 - t174 * t218 + ((t165 * t218 + t263) * qJD(5) + t274) * t204 + t272 * t200) * MDP(24) + (-t153 * t275 + t168 * t222 + t272 * t204 + ((-qJD(5) * t165 + t180) * t218 - qJD(5) * t263 - t274) * t200) * MDP(25) + t273 * t276 + (MDP(14) * t157 - MDP(15) * t275 + MDP(17) * t168 + MDP(18) * t170) * t197; (-t167 * t222 + t226 * t204 + (-t271 * t145 + (-qJD(4) * t222 - t267) * t201 + (t200 * t221 - t204 * t234) * t175) * pkin(2) + t227) * MDP(25) + (-t167 * t197 + t259 + (pkin(2) * t202 * t277 - t250 * t258) * qJD(2) + (-t188 + (-pkin(2) * t197 - t186) * t201) * qJD(4) - t254) * MDP(17) + (-t169 * t197 + (t178 * t249 - t197 * t234) * pkin(2) + t209) * MDP(18) + t211 + (t167 * t161 + t226 * t200 + (-t271 * t146 + (qJD(4) * t161 - t268) * t201 + (-t200 * t234 - t204 * t221) * t175) * pkin(2) + t219) * MDP(24) + t273 * t208; (t165 * t197 - t144 + t259) * MDP(17) + (-t164 * t197 + t209) * MDP(18) + (-t161 * t165 + t200 * t233 + t219) * MDP(24) + (t165 * t222 + t204 * t233 + t227) * MDP(25) + t211; -t222 * t161 * MDP(19) + (-t161 ^ 2 + t222 ^ 2) * MDP(20) + (t255 + t265) * MDP(21) + (-t264 - t269) * MDP(22) + t156 * MDP(23) + (-t143 * t200 + t153 * t175 + t164 * t222 + t174) * MDP(24) + (-t143 * t204 + t161 * t164 - t175 * t225 - t180 * t200) * MDP(25) + (MDP(21) * t262 + MDP(22) * t222 - MDP(24) * t153 + MDP(25) * t225) * qJD(5);];
tauc  = t1;
