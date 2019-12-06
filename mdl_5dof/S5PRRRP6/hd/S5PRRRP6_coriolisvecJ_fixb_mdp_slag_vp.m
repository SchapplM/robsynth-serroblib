% Calculate minimal parameter regressor of Coriolis joint torque vector for
% S5PRRRP6
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% qJD [5x1]
%   Generalized joint velocities
% pkin [8x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d2,d3,d4,theta1]';
% MDP [22x1]
%   Minimal dynamic parameter vector (fixed base model)
%   see S5PRRRP6_convert_par2_MPV_fixb.m
% 
% Output:
% tauc [5x1]
%   joint torques required to compensate Coriolis and centrifugal load

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-05 16:53
% Revision: 77da58f92bca3eff71542919beafa37024070d86 (2019-12-05)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauc = S5PRRRP6_coriolisvecJ_fixb_mdp_slag_vp(qJ, qJD, pkin, MDP)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(8,1),zeros(22,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5PRRRP6_coriolisvecJ_fixb_mdp_slag_vp: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5PRRRP6_coriolisvecJ_fixb_mdp_slag_vp: qJD has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S5PRRRP6_coriolisvecJ_fixb_mdp_slag_vp: pkin has to be [8x1] (double)');
assert(isreal(MDP) && all(size(MDP) == [22 1]), ...
  'S5PRRRP6_coriolisvecJ_fixb_mdp_slag_vp: MDP has to be [22x1] (double)'); 

%% Symbolic Calculation
% From coriolisvec_joint_fixb_mdp_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-05 16:52:21
% EndTime: 2019-12-05 16:52:26
% DurationCPUTime: 1.48s
% Computational Cost: add. (1262->197), mult. (3053->277), div. (0->0), fcn. (2022->6), ass. (0->95)
t219 = sin(qJ(4));
t220 = sin(qJ(3));
t222 = cos(qJ(4));
t223 = cos(qJ(3));
t196 = t219 * t223 + t220 * t222;
t216 = qJD(3) + qJD(4);
t282 = t216 * t196;
t164 = t282 * qJD(2);
t281 = MDP(17) + MDP(19);
t255 = qJD(2) * t223;
t243 = t222 * t255;
t257 = qJD(2) * t220;
t244 = t219 * t257;
t189 = -t243 + t244;
t191 = t196 * qJD(2);
t280 = (MDP(17) * t189 + MDP(18) * t191) * pkin(3);
t269 = t219 * t220;
t195 = -t222 * t223 + t269;
t221 = sin(qJ(2));
t188 = t195 * t221;
t279 = (t220 ^ 2 - t223 ^ 2) * MDP(6);
t259 = qJD(1) * t221;
t203 = qJD(2) * pkin(6) + t259;
t242 = pkin(7) * qJD(2) + t203;
t185 = t242 * t220;
t251 = qJD(4) * t222;
t186 = t242 * t223;
t271 = t186 * t219;
t278 = -pkin(3) * t251 - t185 * t222 - t271;
t277 = -MDP(10) * t223 + MDP(11) * t220;
t249 = MDP(18) - MDP(21);
t275 = t191 ^ 2;
t274 = -pkin(7) - pkin(6);
t273 = qJD(2) * pkin(2);
t214 = -pkin(3) * t223 - pkin(2);
t224 = cos(qJ(2));
t258 = qJD(1) * t224;
t193 = t214 * qJD(2) - t258;
t157 = pkin(4) * t189 - qJ(5) * t191 + t193;
t272 = t157 * t191;
t270 = t186 * t222;
t245 = qJD(3) * t274;
t197 = t220 * t245;
t198 = t223 * t245;
t229 = t195 * t224;
t199 = t274 * t220;
t200 = t274 * t223;
t235 = t199 * t222 + t200 * t219;
t268 = qJD(1) * t229 + t235 * qJD(4) + t197 * t222 + t198 * t219;
t178 = t199 * t219 - t200 * t222;
t230 = t196 * t224;
t267 = -qJD(1) * t230 + t178 * qJD(4) + t197 * t219 - t198 * t222;
t266 = qJD(5) - t278;
t265 = t216 * t243;
t254 = qJD(3) * t220;
t248 = pkin(3) * t254;
t256 = qJD(2) * t221;
t194 = qJD(1) * t256 + qJD(2) * t248;
t260 = MDP(11) * t223;
t253 = qJD(3) * t223;
t252 = qJD(4) * t219;
t181 = qJD(3) * pkin(3) - t185;
t158 = t181 * t222 - t271;
t250 = qJD(5) - t158;
t153 = t265 + (t189 - t244) * t216;
t246 = t153 * MDP(14) + (t191 * t216 - t164) * MDP(15) + (-t189 ^ 2 + t275) * MDP(13);
t236 = t216 * t269;
t171 = -t222 * t253 - t223 * t251 + t236;
t240 = -pkin(4) * t282 - qJ(5) * t171 + qJD(5) * t196 - t248 + t259;
t173 = -t203 * t254 + (-pkin(7) * t254 + t223 * t258) * qJD(2);
t174 = -t203 * t253 + (-pkin(7) * t253 - t220 * t258) * qJD(2);
t238 = t222 * t173 + t219 * t174 + t181 * t251 - t186 * t252;
t144 = t219 * t173 - t222 * t174 + t181 * t252 + t186 * t251;
t237 = pkin(3) * t252 + t185 * t219 - t270;
t165 = pkin(4) * t191 + qJ(5) * t189;
t215 = t216 * qJD(5);
t143 = t215 + t238;
t159 = t181 * t219 + t270;
t233 = -0.2e1 * t273;
t232 = t158 * t216 - t238;
t231 = t159 * t216 - t144;
t228 = t191 * MDP(12) + t193 * MDP(18) - t157 * MDP(21);
t226 = qJD(2) ^ 2;
t225 = qJD(3) ^ 2;
t213 = -pkin(3) * t222 - pkin(4);
t210 = pkin(3) * t219 + qJ(5);
t187 = t196 * t221;
t166 = pkin(4) * t195 - qJ(5) * t196 + t214;
t163 = qJD(2) * t236 - t265;
t162 = pkin(3) * t257 + t165;
t155 = qJ(5) * t216 + t159;
t154 = -pkin(4) * t216 + t250;
t148 = qJD(2) * t230 - t188 * t216;
t145 = pkin(4) * t164 + qJ(5) * t163 - qJD(5) * t191 + t194;
t1 = [(t148 * t191 - t163 * t187 + t164 * t188) * MDP(20) + (-t143 * t188 + t144 * t187 + t148 * t154) * MDP(22) + (-t189 * MDP(20) + t155 * MDP(22) - t249 * t216) * (-qJD(2) * t229 - t221 * t282) + (-t145 * MDP(22) - t226 * MDP(4) - t281 * t164 + t249 * t163 + 0.2e1 * (-MDP(10) * t220 - t260) * qJD(3) * qJD(2)) * t224 + (-t226 * MDP(3) + (t157 * MDP(22) + t249 * t191) * qJD(2) + t277 * (t225 + t226)) * t221 + t281 * (-t148 * t216 + t189 * t256); -t163 * t196 * MDP(12) + (t163 * t195 - t164 * t196) * MDP(13) + (t164 * t214 + t193 * t282 + t194 * t195) * MDP(17) + (-t163 * t214 - t171 * t193 + t194 * t196) * MDP(18) + (t145 * t195 + t157 * t282 + t164 * t166) * MDP(19) + (-t143 * t195 + t144 * t196 - t154 * t171 - t155 * t282 + t163 * t235 - t164 * t178) * MDP(20) + (-t145 * t196 + t157 * t171 + t163 * t166) * MDP(21) + (t143 * t178 - t144 * t235 + t145 * t166 + t267 * t154 + t268 * t155 - t240 * t157) * MDP(22) + (-t171 * MDP(12) - MDP(13) * t282 - MDP(18) * t259 + t267 * MDP(20) + t240 * MDP(21)) * t191 + (t171 * MDP(13) - MDP(17) * t259 - t240 * MDP(19) - t268 * MDP(20)) * t189 + (t223 * MDP(7) - t220 * MDP(8) + pkin(6) * t277) * t225 + (-t171 * MDP(14) - MDP(15) * t282 - t249 * t268 - t267 * t281) * t216 + (-0.2e1 * qJD(2) * t279 + t233 * t260 + (t233 * MDP(10) + 0.2e1 * MDP(5) * t255 + t280) * t220) * qJD(3); -t238 * MDP(18) + (-t163 * t213 - t164 * t210) * MDP(20) + t143 * MDP(21) + (t143 * t210 + t237 * t154 + t266 * t155 - t157 * t162) * MDP(22) + (-t193 * MDP(17) - t157 * MDP(19) + (t155 + t237) * MDP(20) + t162 * MDP(21)) * t191 + (-t220 * t223 * MDP(5) + t279) * t226 + (-t162 * MDP(19) + (t154 - t266) * MDP(20) + t228) * t189 + (MDP(18) * t278 + t266 * MDP(21) - t237 * t281) * t216 + (t273 * t260 + (MDP(10) * t273 - t280) * t220) * qJD(2) + t246 + (MDP(22) * t213 - t281) * t144; (-t191 * t193 + t231) * MDP(17) + t232 * MDP(18) + (t231 - t272) * MDP(19) + (pkin(4) * t163 - qJ(5) * t164 - (-t155 + t159) * t191) * MDP(20) + (t165 * t191 + 0.2e1 * t215 - t232) * MDP(21) + (-pkin(4) * t144 + qJ(5) * t143 - t154 * t159 + t250 * t155 - t157 * t165) * MDP(22) + (-t165 * MDP(19) + (t154 - t250) * MDP(20) + t228) * t189 + t246; t191 * t189 * MDP(19) + t153 * MDP(20) + (-t216 ^ 2 - t275) * MDP(21) + (-t155 * t216 + t144 + t272) * MDP(22);];
tauc = t1;
