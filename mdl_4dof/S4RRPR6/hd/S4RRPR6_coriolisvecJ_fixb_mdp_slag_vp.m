% Calculate Coriolis joint torque vector for
% S4RRPR6
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [4x1]
%   Generalized joint coordinates (joint angles)
% qJD [4x1]
%   Generalized joint velocities
% pkin [7x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,d1,d2,d4,theta3]';
% MDP [21x1]
%   Minimal dynamic parameter vector (fixed base model)
%   see S4RRPR6_convert_par2_MPV_fixb.m
% 
% Output:
% tauc [4x1]
%   joint torques required to compensate Coriolis and centrifugal load

% Quelle: HybrDyn-Toolbox
% Datum: 2021-01-15 10:46
% Revision: d12c3222fdeb2c5f3b3c8fa5751e113be2fc3aae (2021-01-15)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauc = S4RRPR6_coriolisvecJ_fixb_mdp_slag_vp(qJ, qJD, pkin, MDP)
%% Coder Information
%#codegen
%$cgargs {zeros(4,1),zeros(4,1),zeros(7,1),zeros(21,1)}
assert(isreal(qJ) && all(size(qJ) == [4 1]), ...
  'S4RRPR6_coriolisvecJ_fixb_mdp_slag_vp: qJ has to be [4x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [4 1]), ...
  'S4RRPR6_coriolisvecJ_fixb_mdp_slag_vp: qJD has to be [4x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [7 1]), ...
  'S4RRPR6_coriolisvecJ_fixb_mdp_slag_vp: pkin has to be [7x1] (double)');
assert(isreal(MDP) && all(size(MDP) == [21 1]), ...
  'S4RRPR6_coriolisvecJ_fixb_mdp_slag_vp: MDP has to be [21x1] (double)'); 

%% Symbolic Calculation
% From coriolisvec_joint_fixb_mdp_matlab.m
% OptimizationMode: 2
% StartTime: 2021-01-15 10:46:27
% EndTime: 2021-01-15 10:46:32
% DurationCPUTime: 1.37s
% Computational Cost: add. (814->175), mult. (2223->251), div. (0->0), fcn. (1551->6), ass. (0->89)
t227 = sin(pkin(7));
t228 = cos(pkin(7));
t232 = cos(qJ(2));
t249 = qJD(1) * t232;
t242 = t228 * t249;
t230 = sin(qJ(2));
t250 = qJD(1) * t230;
t197 = t227 * t250 - t242;
t231 = cos(qJ(4));
t188 = t231 * t197;
t209 = t227 * t232 + t228 * t230;
t199 = t209 * qJD(2);
t191 = qJD(1) * t199;
t248 = qJD(2) * t230;
t241 = qJD(1) * t248;
t218 = t227 * t241;
t192 = qJD(2) * t242 - t218;
t229 = sin(qJ(4));
t247 = qJD(4) * t229;
t252 = qJD(1) * t209;
t138 = -qJD(4) * t188 - t229 * t191 + t231 * t192 - t247 * t252;
t161 = -t229 * t252 - t188;
t236 = t197 * t229 - t231 * t252;
t235 = t236 * qJD(4) - t231 * t191 - t192 * t229;
t224 = qJD(2) + qJD(4);
t257 = t161 * t224;
t258 = t236 * t224;
t270 = t161 * t236 * MDP(15) + (-t161 ^ 2 + t236 ^ 2) * MDP(16) + (t138 - t257) * MDP(17) + (t235 - t258) * MDP(18);
t259 = -qJ(3) - pkin(5);
t216 = t259 * t230;
t213 = qJD(1) * t216;
t207 = qJD(2) * pkin(2) + t213;
t217 = t259 * t232;
t214 = qJD(1) * t217;
t256 = t228 * t214;
t164 = t227 * t207 - t256;
t261 = pkin(6) * t197;
t149 = t164 - t261;
t223 = -pkin(2) * t232 - pkin(1);
t251 = qJD(1) * t223;
t215 = qJD(3) + t251;
t169 = pkin(3) * t197 + t215;
t269 = t149 * t247 - t169 * t161;
t240 = qJD(2) * t259;
t195 = qJD(3) * t232 + t230 * t240;
t181 = t195 * qJD(1);
t196 = -qJD(3) * t230 + t232 * t240;
t182 = t196 * qJD(1);
t150 = -t181 * t227 + t228 * t182;
t143 = -pkin(6) * t192 + t150;
t151 = t228 * t181 + t227 * t182;
t144 = -pkin(6) * t191 + t151;
t268 = t231 * t143 - t229 * t144 + t169 * t236;
t266 = t232 * MDP(4) - pkin(1) * MDP(9);
t265 = qJD(4) - t224;
t264 = pkin(1) * t232 * MDP(10) + (t230 ^ 2 - t232 ^ 2) * MDP(5);
t262 = pkin(2) * t227;
t260 = pkin(6) * t252;
t203 = t227 * t214;
t155 = t228 * t195 + t227 * t196;
t168 = t228 * t213 + t203;
t172 = t227 * t216 - t228 * t217;
t245 = 0.2e1 * qJD(1);
t243 = pkin(2) * t250;
t154 = -t195 * t227 + t228 * t196;
t163 = t228 * t207 + t203;
t167 = -t213 * t227 + t256;
t171 = t228 * t216 + t217 * t227;
t148 = qJD(2) * pkin(3) + t163 - t260;
t237 = -t229 * t148 - t231 * t149;
t208 = t227 * t230 - t228 * t232;
t165 = t208 * t231 + t209 * t229;
t166 = -t208 * t229 + t209 * t231;
t222 = pkin(2) * t228 + pkin(3);
t220 = pkin(2) * t241;
t202 = t208 * qJD(2);
t184 = pkin(3) * t208 + t223;
t177 = pkin(2) * t248 + pkin(3) * t199;
t176 = pkin(3) * t252 + t243;
t170 = pkin(3) * t191 + t220;
t157 = -pkin(6) * t208 + t172;
t156 = -pkin(6) * t209 + t171;
t153 = t168 - t260;
t152 = t167 + t261;
t146 = -pkin(6) * t199 + t155;
t145 = pkin(6) * t202 + t154;
t141 = qJD(4) * t166 + t231 * t199 - t202 * t229;
t140 = -qJD(4) * t165 - t199 * t229 - t202 * t231;
t1 = [(t191 * t223 + t199 * t215) * MDP(11) + (t192 * t223 - t202 * t215) * MDP(12) + (-t150 * t209 - t151 * t208 - t154 * t252 - t155 * t197 + t163 * t202 - t164 * t199 - t171 * t192 - t172 * t191) * MDP(13) + (t150 * t171 + t151 * t172 + t154 * t163 + t155 * t164) * MDP(14) + (t138 * t166 - t140 * t236) * MDP(15) + (-t138 * t165 + t140 * t161 + t141 * t236 + t166 * t235) * MDP(16) + (t169 * t141 - t161 * t177 + t170 * t165 - t184 * t235) * MDP(20) + (t184 * t138 + t169 * t140 + t170 * t166 - t177 * t236) * MDP(21) + (t140 * MDP(17) - t141 * MDP(18) + (t145 * t231 - t146 * t229) * MDP(20) + (-t145 * t229 - t146 * t231) * MDP(21) + ((-t156 * t229 - t157 * t231) * MDP(20) + (-t156 * t231 + t157 * t229) * MDP(21)) * qJD(4)) * t224 + (t154 * MDP(11) - t155 * MDP(12) - t264 * t245 + (t266 * t245 + ((qJD(1) * t208 + t197) * MDP(11) + 0.2e1 * t252 * MDP(12) + (t215 + t251) * MDP(14)) * pkin(2)) * t230 + (t232 * MDP(6) - t230 * MDP(7) + (MDP(10) * t230 - MDP(9) * t232) * pkin(5)) * qJD(2)) * qJD(2); (-qJD(2) * t167 - t197 * t243 - t215 * t252 + t150) * MDP(11) + (qJD(2) * t168 + t197 * t215 - t243 * t252 - t151) * MDP(12) + ((t164 + t167) * t252 + (-t163 + t168) * t197 + (-t191 * t227 - t192 * t228) * pkin(2)) * MDP(13) + (-t163 * t167 - t164 * t168 + (t150 * t228 + t151 * t227 - t215 * t250) * pkin(2)) * MDP(14) + (t176 * t161 - (t152 * t231 - t153 * t229) * t224 + ((-t222 * t229 - t231 * t262) * t224 + t237) * qJD(4) + t268) * MDP(20) + (-t231 * t144 - t229 * t143 + t176 * t236 + (t152 * t229 + t153 * t231) * t224 + (-(t222 * t231 - t229 * t262) * t224 - t231 * t148) * qJD(4) + t269) * MDP(21) + (-t266 * t230 + t264) * qJD(1) ^ 2 + t270; -t218 * MDP(12) + (-t197 ^ 2 - t252 ^ 2) * MDP(13) + (t163 * t252 + t164 * t197 + t220) * MDP(14) + (-t235 - t258) * MDP(20) + (t138 + t257) * MDP(21) + ((t227 * t249 + t228 * t250 + t252) * MDP(11) + (-t197 + t242) * MDP(12)) * qJD(2); (t265 * t237 + t268) * MDP(20) + ((-t149 * t224 - t143) * t229 + (-t265 * t148 - t144) * t231 + t269) * MDP(21) + t270;];
tauc = t1;
