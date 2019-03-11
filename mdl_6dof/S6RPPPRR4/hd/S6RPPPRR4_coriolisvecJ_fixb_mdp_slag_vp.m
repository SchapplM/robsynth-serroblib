% Calculate minimal parameter regressor of Coriolis joint torque vector for
% S6RPPPRR4
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% qJD [6x1]
%   Generalized joint velocities
% pkin [9x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d5,d6,theta3]';
% MDP [26x1]
%   Minimal dynamic parameter vector (fixed base model)
%   see S6RPPPRR4_convert_par2_MPV_fixb.m
% 
% Output:
% tauc [6x1]
%   joint torques required to compensate Coriolis and centrifugal load

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-09 01:36
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauc = S6RPPPRR4_coriolisvecJ_fixb_mdp_slag_vp(qJ, qJD, pkin, MDP)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(9,1),zeros(26,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPPPRR4_coriolisvecJ_fixb_mdp_slag_vp: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RPPPRR4_coriolisvecJ_fixb_mdp_slag_vp: qJD has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S6RPPPRR4_coriolisvecJ_fixb_mdp_slag_vp: pkin has to be [9x1] (double)');
assert(isreal(MDP) && all(size(MDP) == [26 1]), ...
  'S6RPPPRR4_coriolisvecJ_fixb_mdp_slag_vp: MDP has to be [26x1] (double)'); 

%% Symbolic Calculation
% From coriolisvec_joint_fixb_mdp_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 01:35:55
% EndTime: 2019-03-09 01:35:58
% DurationCPUTime: 1.52s
% Computational Cost: add. (843->233), mult. (1587->334), div. (0->0), fcn. (876->6), ass. (0->113)
t181 = sin(qJ(5));
t235 = qJD(1) * t181;
t165 = -qJD(6) + t235;
t180 = sin(qJ(6));
t182 = cos(qJ(6));
t263 = t165 * (MDP(25) * t182 - MDP(26) * t180);
t262 = qJ(2) * MDP(6) + MDP(5);
t261 = MDP(13) * t181;
t183 = cos(qJ(5));
t176 = t183 ^ 2;
t260 = MDP(14) * (t181 ^ 2 - t176);
t259 = MDP(7) - MDP(10);
t184 = -pkin(1) - pkin(2);
t163 = qJD(1) * t184 + qJD(2);
t177 = sin(pkin(9));
t178 = cos(pkin(9));
t237 = qJD(1) * qJ(2);
t150 = t163 * t178 - t177 * t237;
t147 = qJD(1) * pkin(3) + qJD(4) - t150;
t143 = qJD(1) * pkin(7) + t147;
t140 = qJD(3) * t183 + t143 * t181;
t222 = qJD(1) * qJD(2);
t211 = t177 * t222;
t136 = t140 * qJD(5) - t183 * t211;
t257 = t136 * t180;
t256 = t136 * t182;
t226 = qJD(6) * t183;
t208 = qJD(1) * t226;
t221 = qJD(1) * qJD(5);
t209 = t182 * t221;
t220 = qJD(5) * qJD(6);
t145 = t180 * t208 + t181 * t209 + t182 * t220;
t255 = t145 * t180;
t254 = t145 * t181;
t210 = t181 * t221;
t241 = (-t210 - t220) * t180;
t146 = t182 * t208 + t241;
t253 = t146 * t181;
t234 = qJD(1) * t183;
t216 = t180 * t234;
t230 = qJD(5) * t182;
t155 = t216 + t230;
t252 = t155 * t165;
t251 = t155 * t183;
t215 = t182 * t234;
t224 = t180 * qJD(5);
t156 = t215 - t224;
t250 = t156 * t165;
t249 = t156 * t183;
t248 = t165 * t180;
t247 = t165 * t182;
t246 = t177 * t163;
t171 = t178 * qJ(2);
t244 = t180 * t181;
t243 = t181 * t182;
t242 = t183 * t145;
t240 = t177 * t184 + t171;
t185 = qJD(5) ^ 2;
t186 = qJD(1) ^ 2;
t238 = t185 + t186;
t236 = qJD(1) * t176;
t233 = qJD(2) * t177;
t203 = -t177 * qJ(2) + t178 * t184;
t199 = pkin(3) - t203;
t153 = pkin(7) + t199;
t232 = qJD(5) * t153;
t231 = qJD(5) * t181;
t229 = qJD(5) * t183;
t228 = qJD(6) * t180;
t227 = qJD(6) * t182;
t225 = t178 * qJD(2);
t223 = t186 * MDP(11);
t218 = t165 * t244;
t217 = t165 * t243;
t154 = -qJ(4) + t240;
t214 = t165 * t228;
t213 = t180 * t226;
t212 = t165 * t227;
t207 = t178 * t222;
t206 = MDP(24) * t229;
t205 = t178 * t238;
t138 = qJD(5) * pkin(8) + t140;
t204 = t153 * t165 + t138;
t151 = t178 * t237 + t246;
t174 = qJD(1) * qJ(4);
t148 = t151 - t174;
t202 = t148 + t233;
t201 = -0.2e1 * t183 * t221;
t200 = t183 * t212;
t198 = -pkin(5) * t183 - pkin(8) * t181;
t197 = -pkin(5) * t181 + pkin(8) * t183;
t141 = t246 - t174 + (t197 + t171) * qJD(1);
t134 = t138 * t182 + t141 * t180;
t196 = t138 * t180 - t141 * t182;
t195 = t150 * t177 - t151 * t178;
t194 = qJD(3) * t181 - t143 * t183;
t192 = -t165 * t213 + t176 * t209;
t137 = -qJD(5) * pkin(5) + t194;
t191 = pkin(8) * t229 - t137 * t181;
t173 = qJD(1) * qJD(4);
t158 = -t173 + t207;
t164 = -qJD(4) + t225;
t190 = -qJD(1) * t164 - t153 * t185 - t158;
t189 = qJD(1) * t154 + t148 - t233;
t188 = qJD(5) * t198 + t225;
t135 = -qJD(5) * t194 + t181 * t211;
t187 = qJD(5) * t137 + qJD(6) * t141 + t165 * t233 + t135;
t157 = t198 * qJD(1);
t152 = -qJD(4) + t188;
t149 = t197 + t154;
t144 = qJD(1) * t188 - t173;
t142 = t182 * t144;
t1 = [0.2e1 * MDP(8) * t207 + ((-t177 * t203 + t178 * t240) * qJD(1) - t195) * qJD(2) * MDP(9) + (t173 + (-t164 - t225) * qJD(1)) * MDP(11) + (t148 * t164 + t158 * t154 + (qJD(1) * t199 + t147) * t233) * MDP(12) + t201 * t261 + 0.2e1 * t221 * t260 + (t181 * t190 - t189 * t229) * MDP(18) + (t183 * t190 + t189 * t231) * MDP(19) + (-t182 * t242 + (-t181 * t230 - t213) * t156) * MDP(20) + ((t155 * t182 + t156 * t180) * t231 + (t255 - t146 * t182 + (t155 * t180 - t156 * t182) * qJD(6)) * t183) * MDP(21) + (-t254 + (-t217 + t249) * qJD(5) + t192) * MDP(22) + (-t200 - t253 + (-t251 + (t165 * t181 - t236) * t180) * qJD(5)) * MDP(23) + (t165 + t235) * t206 + (-(-t149 * t228 + t182 * t152) * t165 + (-t155 * t232 + t180 * t187 + t204 * t227 - t142) * t181 + (t155 * t233 - t137 * t227 - t257 + t153 * t146 + (t153 * t248 - (t149 * t182 - t153 * t244) * qJD(1) + t196) * qJD(5)) * t183) * MDP(25) + ((t149 * t227 + t180 * t152) * t165 + (-t156 * t232 + (-qJD(6) * t204 + t144) * t180 + t187 * t182) * t181 + (t156 * t233 + t137 * t228 - t256 - t153 * t145 + (t153 * t247 + (t149 * t180 + t153 * t243) * qJD(1) + t134) * qJD(5)) * t183) * MDP(26) + 0.2e1 * t262 * t222 + 0.2e1 * t259 * t211 + (MDP(15) * t181 + MDP(16) * t183) * t185; t195 * qJD(1) * MDP(9) + t178 * t223 + (t158 * t177 + (-t148 * t178 + (-t147 - t225) * t177) * qJD(1)) * MDP(12) + (t177 * t201 + t181 * t205) * MDP(18) + (0.2e1 * t177 * t210 + t183 * t205) * MDP(19) + (t177 * t214 + (-(t181 * t227 + t183 * t224) * t165 + t155 * t231 - t183 * t146) * t178 + ((-t177 * t244 + t178 * t182) * t165 + (-(t177 * t182 + t178 * t244) * qJD(5) - t177 * t155) * t183) * qJD(1)) * MDP(25) + (t177 * t212 + ((t181 * t228 - t182 * t229) * t165 + t156 * t231 + t242) * t178 + (-(t177 * t243 + t178 * t180) * t165 + ((t177 * t180 - t178 * t243) * qJD(5) - t177 * t156) * t183) * qJD(1)) * MDP(26) + (-t178 * MDP(8) - t259 * t177 - t262) * t186; (t200 - t253) * MDP(25) + (t192 + t254) * MDP(26) + (-MDP(18) * t183 + MDP(19) * t181) * t185 + ((t180 * t236 - t218 - t251) * MDP(25) + (-t217 - t249) * MDP(26)) * qJD(5); -t223 + (MDP(12) * t202 - t263) * qJD(1) + (-t238 * MDP(19) + (t165 * t224 + t146) * MDP(25) + (t165 * t230 - t145) * MDP(26)) * t183 + (-t238 * MDP(18) + qJD(6) * t263 + ((-t155 + t216) * MDP(25) + (-t156 + t215) * MDP(26)) * qJD(5)) * t181; -t202 * t235 * MDP(19) + (t156 * t247 + t255) * MDP(20) + ((t145 - t252) * t182 + (t146 - t250) * t180) * MDP(21) + (-t212 + (t217 + (-t156 - t224) * t183) * qJD(1)) * MDP(22) + (t214 + (-t218 + (t155 - t230) * t183) * qJD(1)) * MDP(23) + (pkin(5) * t146 - t256 + (t157 * t182 + t180 * t194) * t165 + t140 * t155 + (pkin(8) * t247 + t137 * t180) * qJD(6) + (t180 * t191 - t183 * t196) * qJD(1)) * MDP(25) + (-pkin(5) * t145 + t257 - (t157 * t180 - t182 * t194) * t165 + t140 * t156 + (-pkin(8) * t248 + t137 * t182) * qJD(6) + (-t134 * t183 + t182 * t191) * qJD(1)) * MDP(26) + (MDP(18) * t202 - t165 * MDP(24)) * t234 + (t183 * t261 - t260) * t186; t156 * t155 * MDP(20) + (-t155 ^ 2 + t156 ^ 2) * MDP(21) + (t145 + t252) * MDP(22) + (t241 + t250) * MDP(23) - qJD(1) * t206 + (-t134 * t165 - t180 * t135 + t137 * t156 + t142) * MDP(25) + (-t182 * t135 - t137 * t155 - t180 * t144 + t165 * t196) * MDP(26) + (MDP(23) * t215 - MDP(25) * t134 + MDP(26) * t196) * qJD(6);];
tauc  = t1;
