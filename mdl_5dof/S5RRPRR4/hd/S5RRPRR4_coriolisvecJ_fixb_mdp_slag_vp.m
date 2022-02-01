% Calculate Coriolis joint torque vector for
% S5RRPRR4
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% qJD [5x1]
%   Generalized joint velocities
% pkin [9x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d1,d2,d4,d5,theta3]';
% MDP [21x1]
%   Minimal dynamic parameter vector (fixed base model)
%   see S5RRPRR4_convert_par2_MPV_fixb.m
% 
% Output:
% tauc [5x1]
%   joint torques required to compensate Coriolis and centrifugal load

% Quelle: HybrDyn-Toolbox
% Datum: 2022-01-20 10:49
% Revision: fd3771346c4aea32fdeb66112c511235427c26a7 (2022-01-20)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauc = S5RRPRR4_coriolisvecJ_fixb_mdp_slag_vp(qJ, qJD, pkin, MDP)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(9,1),zeros(21,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RRPRR4_coriolisvecJ_fixb_mdp_slag_vp: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RRPRR4_coriolisvecJ_fixb_mdp_slag_vp: qJD has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S5RRPRR4_coriolisvecJ_fixb_mdp_slag_vp: pkin has to be [9x1] (double)');
assert(isreal(MDP) && all(size(MDP) == [21 1]), ...
  'S5RRPRR4_coriolisvecJ_fixb_mdp_slag_vp: MDP has to be [21x1] (double)'); 

%% Symbolic Calculation
% From coriolisvec_joint_fixb_mdp_matlab.m
% OptimizationMode: 2
% StartTime: 2022-01-20 10:48:35
% EndTime: 2022-01-20 10:48:38
% DurationCPUTime: 0.88s
% Computational Cost: add. (927->153), mult. (1769->224), div. (0->0), fcn. (1106->8), ass. (0->102)
t222 = sin(qJ(5));
t223 = sin(qJ(4));
t225 = cos(qJ(5));
t226 = cos(qJ(4));
t196 = t222 * t226 + t223 * t225;
t217 = qJD(1) + qJD(2);
t182 = t196 * t217;
t216 = qJD(4) + qJD(5);
t283 = qJD(5) - t216;
t287 = MDP(8) * t223;
t286 = (t223 ^ 2 - t226 ^ 2) * MDP(9);
t227 = cos(qJ(2));
t278 = pkin(1) * qJD(1);
t252 = t227 * t278;
t197 = pkin(2) * t217 + t252;
t221 = cos(pkin(9));
t224 = sin(qJ(2));
t253 = t224 * t278;
t206 = t221 * t253;
t220 = sin(pkin(9));
t173 = t220 * t197 + t206;
t245 = t173 + (pkin(7) + pkin(8)) * t217;
t154 = t226 * qJD(3) - t245 * t223;
t285 = MDP(13) * t223 + MDP(14) * t226;
t284 = MDP(5) * t224 + MDP(6) * t227;
t155 = t223 * qJD(3) + t245 * t226;
t282 = pkin(2) * t221;
t281 = pkin(4) * t226;
t212 = pkin(1) * t227 + pkin(2);
t272 = t221 * t224;
t263 = pkin(1) * t272 + t220 * t212;
t184 = pkin(7) + t263;
t280 = -pkin(8) - t184;
t209 = pkin(2) * t220 + pkin(7);
t279 = -pkin(8) - t209;
t277 = pkin(1) * qJD(2);
t205 = t220 * t253;
t187 = t221 * t252 - t205;
t274 = t187 * t216;
t273 = t220 * t224;
t271 = t222 * t223;
t268 = t225 * t155;
t267 = t225 * t226;
t172 = t221 * t197 - t205;
t248 = -pkin(3) - t281;
t160 = t248 * t217 - t172;
t195 = -t267 + t271;
t165 = t216 * t195;
t186 = (t220 * t227 + t272) * t277;
t178 = qJD(1) * t186;
t259 = qJD(4) * t223;
t251 = pkin(4) * t259;
t167 = t217 * t251 + t178;
t266 = -t160 * t165 + t167 * t196;
t166 = t216 * t196;
t265 = t160 * t166 + t167 * t195;
t168 = -pkin(3) * t217 - t172;
t258 = qJD(4) * t226;
t264 = t168 * t258 + t178 * t223;
t228 = qJD(4) ^ 2;
t256 = t228 * MDP(11);
t254 = pkin(4) * t217 * t223;
t250 = t217 * t271;
t249 = t217 * t267;
t247 = t217 * t258;
t153 = qJD(4) * pkin(4) + t154;
t246 = -pkin(4) * t216 - t153;
t244 = qJD(4) * t280;
t243 = qJD(4) * t279;
t241 = -pkin(1) * t273 + t212 * t221;
t183 = -pkin(3) - t241;
t185 = t220 * t252 + t206;
t240 = -t185 + t251;
t236 = t184 * t228 + t186 * t217;
t235 = -t185 * t217 + t209 * t228;
t188 = (t221 * t227 - t273) * t277;
t234 = qJD(4) * (t183 * t217 - t188);
t233 = qJD(4) * ((-pkin(3) - t282) * t217 + t187);
t179 = qJD(1) * t188;
t143 = qJD(4) * t154 + t226 * t179;
t144 = -qJD(4) * t155 - t223 * t179;
t232 = -t222 * t143 + t225 * t144 - t160 * t182;
t156 = qJD(5) * t249 - t216 * t250 + t225 * t247;
t180 = -t249 + t250;
t231 = t182 * t180 * MDP(15) + (-t180 ^ 2 + t182 ^ 2) * MDP(16) + (t180 * t216 + t156) * MDP(17);
t157 = t166 * t217;
t230 = (-t156 * t195 - t157 * t196 + t165 * t180 - t166 * t182) * MDP(16) + (t156 * t196 - t165 * t182) * MDP(15) - 0.2e1 * t217 * qJD(4) * t286 + 0.2e1 * t247 * t287 + t228 * t226 * MDP(10) + (-t165 * MDP(17) - t166 * MDP(18)) * t216;
t229 = t160 * t180 + (t283 * t155 - t144) * t222;
t214 = t226 * pkin(8);
t200 = t248 - t282;
t193 = t209 * t226 + t214;
t192 = t279 * t223;
t190 = t226 * t243;
t189 = t223 * t243;
t176 = t183 - t281;
t175 = t186 + t251;
t171 = t184 * t226 + t214;
t170 = t280 * t223;
t163 = t168 * t259;
t152 = -t223 * t188 + t226 * t244;
t151 = t226 * t188 + t223 * t244;
t1 = [(-t172 * t186 + t173 * t188 - t178 * t241 + t179 * t263) * MDP(7) + t163 * MDP(13) + t264 * MDP(14) + (t175 * t180 + t176 * t157 + (-t222 * t151 + t225 * t152 + (-t170 * t222 - t171 * t225) * qJD(5)) * t216 + t265) * MDP(20) + (t175 * t182 + t176 * t156 - (t225 * t151 + t222 * t152 + (t170 * t225 - t171 * t222) * qJD(5)) * t216 + t266) * MDP(21) + ((-t178 - t236) * MDP(13) + MDP(14) * t234) * t226 + (MDP(13) * t234 + t236 * MDP(14) - t256) * t223 + t230 + t284 * t277 * (-qJD(1) - t217); (t172 * t185 - t173 * t187 + (-t178 * t221 + t179 * t220) * pkin(2)) * MDP(7) - t223 * t256 + (t163 + t223 * t233 + (-t178 - t235) * t226) * MDP(13) + (t235 * t223 + t226 * t233 + t264) * MDP(14) + (t200 * t157 + (-t222 * t189 + t225 * t190 + (-t192 * t222 - t193 * t225) * qJD(5)) * t216 + t196 * t274 + t240 * t180 + t265) * MDP(20) + (t200 * t156 - (t225 * t189 + t222 * t190 + (t192 * t225 - t193 * t222) * qJD(5)) * t216 - t195 * t274 + t240 * t182 + t266) * MDP(21) + t230 + t284 * (-qJD(2) + t217) * t278; -t285 * t228 + (-MDP(20) * t166 + MDP(21) * t165) * t216; (-t180 * t254 - (-t154 * t222 - t268) * t216 + (t246 * t222 - t268) * qJD(5) + t232) * MDP(20) + (-t182 * t254 + (t246 * qJD(5) + t154 * t216 - t143) * t225 + t229) * MDP(21) + t231 + t285 * (-t168 * t217 - t179) + (-t226 * t287 + t286) * t217 ^ 2; (t232 + t283 * (-t153 * t222 - t268)) * MDP(20) + ((-t153 * t283 - t143) * t225 + t229) * MDP(21) + t231;];
tauc = t1;
