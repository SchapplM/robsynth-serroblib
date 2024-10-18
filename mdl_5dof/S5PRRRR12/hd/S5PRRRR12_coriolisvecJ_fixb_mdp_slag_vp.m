% Calculate Coriolis joint torque vector for
% S5PRRRR12
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% qJD [5x1]
%   Generalized joint velocities
% pkin [11x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,alpha2,alpha5,d2,d3,d4,d5,theta1]';
% MDP [17x1]
%   Minimal dynamic parameter vector (fixed base model)
%   see S5PRRRR12_convert_par2_MPV_fixb.m
% 
% Output:
% tauc [5x1]
%   joint torques required to compensate Coriolis and centrifugal load

% Quelle: HybrDyn-Toolbox
% Datum: 2024-09-28 18:09
% Revision: 4128c0dfad36ceb08eaa7a9fbc8f797c888f59be (2024-07-29)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauc = S5PRRRR12_coriolisvecJ_fixb_mdp_slag_vp(qJ, qJD, pkin, MDP)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(11,1),zeros(17,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5PRRRR12_coriolisvecJ_fixb_mdp_slag_vp: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5PRRRR12_coriolisvecJ_fixb_mdp_slag_vp: qJD has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [11 1]), ...
  'S5PRRRR12_coriolisvecJ_fixb_mdp_slag_vp: pkin has to be [11x1] (double)');
assert(isreal(MDP) && all(size(MDP) == [17 1]), ...
  'S5PRRRR12_coriolisvecJ_fixb_mdp_slag_vp: MDP has to be [17x1] (double)'); 

%% Symbolic Calculation
% From coriolisvec_joint_fixb_mdp_matlab.m
% OptimizationMode: 2
% StartTime: 2024-09-28 18:08:21
% EndTime: 2024-09-28 18:08:23
% DurationCPUTime: 0.84s
% Computational Cost: add. (1445->168), mult. (3647->270), div. (0->0), fcn. (2816->12), ass. (0->108)
t210 = cos(qJ(2));
t200 = sin(pkin(5));
t254 = qJD(1) * t200;
t241 = t210 * t254;
t187 = qJD(2) * pkin(2) + t241;
t209 = cos(qJ(3));
t205 = sin(qJ(3));
t206 = sin(qJ(2));
t242 = t206 * t254;
t229 = t205 * t242;
t175 = t209 * t187 - t229;
t196 = qJD(2) + qJD(3);
t172 = pkin(3) * t196 + t175;
t204 = sin(qJ(4));
t176 = t187 * t205 + t209 * t242;
t208 = cos(qJ(4));
t271 = t176 * t208;
t154 = -t172 * t204 - t271;
t193 = qJD(4) + t196;
t199 = sin(pkin(6));
t269 = t193 * t199;
t151 = pkin(10) * t269 - t154;
t201 = cos(pkin(6));
t268 = t193 * t201;
t188 = qJD(5) + t268;
t247 = qJD(5) - t188;
t280 = t247 * t151;
t228 = qJD(2) * t241;
t252 = qJD(3) * t187;
t256 = t196 * t229;
t159 = (t228 + t252) * t209 - t256;
t259 = t206 * t209;
t222 = t205 * t210 + t259;
t219 = t222 * qJD(2);
t212 = (-qJD(3) * t259 - t219) * t254;
t160 = -t205 * t252 + t212;
t251 = qJD(4) * t204;
t234 = t160 * t204 - t176 * t251;
t144 = (qJD(4) * t172 + t159) * t208 + t234;
t272 = t176 * t204;
t153 = t208 * t172 - t272;
t275 = pkin(4) * t193;
t152 = t153 + t275;
t274 = t152 * t201;
t279 = -t247 * t274 - t144;
t180 = (-t205 * t206 + t209 * t210) * t200;
t203 = sin(qJ(5));
t207 = cos(qJ(5));
t277 = MDP(12) * (t203 ^ 2 - t207 ^ 2);
t181 = t222 * t200;
t163 = t180 * t208 - t181 * t204;
t202 = cos(pkin(5));
t276 = -(t163 * t201 + t199 * t202) * t188 + (-t163 * t199 + t201 * t202) * t269;
t194 = t199 * pkin(10);
t164 = t180 * t204 + t181 * t208;
t273 = t164 * t188;
t195 = t199 ^ 2;
t270 = t193 * t195;
t267 = t195 * t203;
t266 = t195 * t207;
t264 = t201 * t203;
t263 = t201 * t207;
t235 = -t159 * t204 + t208 * t160;
t145 = t154 * qJD(4) + t235;
t262 = t203 * t145;
t261 = t204 * t205;
t260 = t205 * t208;
t177 = t222 * t254;
t178 = qJD(1) * t180;
t191 = pkin(2) * t209 + pkin(3);
t250 = qJD(4) * t208;
t258 = -t177 * t208 - t178 * t204 + t191 * t251 - (-t205 * t250 + (-t204 * t209 - t260) * qJD(3)) * pkin(2);
t257 = -t177 * t204 + t178 * t208 - t191 * t250 - (-t205 * t251 + (t208 * t209 - t261) * qJD(3)) * pkin(2);
t253 = qJD(1) * t202;
t249 = qJD(5) * t203;
t248 = qJD(5) * t207;
t245 = t193 * t266;
t165 = t196 * t180;
t166 = (-t222 * qJD(3) - t219) * t200;
t147 = -t164 * qJD(4) - t165 * t204 + t166 * t208;
t244 = t147 * t270;
t143 = t145 * t263;
t150 = -t152 * t199 + t201 * t253;
t220 = t199 * t253 + t274;
t240 = t199 * t249;
t243 = (-t203 * t144 + t143 + (-t151 * t207 - t220 * t203) * qJD(5)) * t201 + t145 * t266 + t150 * t240;
t239 = t199 * t248;
t238 = -pkin(3) * t193 - t172;
t237 = -(t201 * t262 + t207 * t144 + (-t151 * t203 + t220 * t207) * qJD(5)) * t201 + t150 * t239;
t236 = t154 * t193 - t145;
t230 = (MDP(13) * t239 - MDP(14) * t240) * (t188 + t268) + 0.2e1 * (MDP(11) * t203 * t207 - t277) * qJD(5) * t270;
t155 = -t175 * t204 - t271;
t227 = pkin(3) * t251 + t155;
t156 = t175 * t208 - t272;
t226 = -pkin(3) * t250 + t156;
t225 = t247 * t253;
t224 = (-pkin(2) * t196 - t187) * qJD(3);
t182 = -pkin(2) * t261 + t191 * t208 + pkin(4);
t217 = -t182 * t249 - t258 * t207;
t216 = -t182 * t248 + t258 * t203;
t190 = pkin(3) * t208 + pkin(4);
t214 = -t190 * t249 - t227 * t207;
t213 = -t190 * t248 + t227 * t203;
t192 = t193 ^ 2;
t189 = pkin(3) * t204 + t194;
t179 = pkin(2) * t260 + t191 * t204 + t194;
t146 = t163 * qJD(4) + t165 * t208 + t166 * t204;
t1 = [((-t146 * t203 + t147 * t263) * t188 + t207 * t244 + (t276 * t203 - t207 * t273) * qJD(5)) * MDP(16) + (-(t146 * t207 + t147 * t264) * t188 - t203 * t244 + (t203 * t273 + t276 * t207) * qJD(5)) * MDP(17) + (-t206 * MDP(3) - t210 * MDP(4)) * t200 * qJD(2) ^ 2 + (t166 * MDP(6) - t165 * MDP(7)) * t196 + (-t146 * MDP(10) + t147 * MDP(9)) * t193; (t177 * t196 + t205 * t224 + t212) * MDP(6) + (t178 * t196 + (t224 - t228) * t209 + t256) * MDP(7) + (-t258 * t193 + t145) * MDP(9) + (t257 * t193 - t144) * MDP(10) + ((-t179 * t248 + t217 * t201 + t257 * t203) * t188 + t217 * t270 + t243) * MDP(16) + ((t179 * t249 + t216 * t201 + t257 * t207) * t188 + (t216 * t193 - t262) * t195 + t237) * MDP(17) + t230; (t176 * t196 + t160) * MDP(6) + (t175 * t196 - t159) * MDP(7) + (-t155 * t193 + (t238 * t204 - t271) * qJD(4) + t235) * MDP(9) + (t156 * t193 + (t238 * qJD(4) - t159) * t208 - t234) * MDP(10) + ((-t189 * t248 + t214 * t201 + t226 * t203) * t188 + t214 * t270 + t243) * MDP(16) + ((t189 * t249 + t213 * t201 + t226 * t207) * t188 + (t213 * t193 - t262) * t195 + t237) * MDP(17) + t230; -t236 * MDP(9) + (t153 * t193 - t144) * MDP(10) + (-(-t153 * t203 + t154 * t263) * t188 - t154 * t245 + ((-pkin(4) * t264 - t207 * t194) * t188 - t267 * t275) * qJD(5) + t243) * MDP(16) + ((t153 * t207 + t154 * t264) * t188 + t236 * t267 + (-(pkin(4) * t263 - t203 * t194) * t188 - pkin(4) * t245) * qJD(5) + t237) * MDP(17) + t230; t143 * MDP(16) + t195 * t192 * t277 + (t279 * MDP(17) - MDP(16) * t280 + (-MDP(17) * t225 + (t247 * MDP(13) - t150 * MDP(17)) * t193) * t199) * t207 + (-t192 * MDP(11) * t266 + t279 * MDP(16) + (-t201 * t145 + t280) * MDP(17) + (-MDP(16) * t225 + (-t247 * MDP(14) - t150 * MDP(16)) * t193) * t199) * t203;];
tauc = t1;
