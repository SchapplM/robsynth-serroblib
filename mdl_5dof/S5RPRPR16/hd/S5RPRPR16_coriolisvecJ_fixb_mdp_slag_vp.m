% Calculate minimal parameter regressor of Coriolis joint torque vector for
% S5RPRPR16
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% qJD [5x1]
%   Generalized joint velocities
% pkin [7x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d1,d3,d5]';
% MDP [24x1]
%   Minimal dynamic parameter vector (fixed base model)
%   see S5RPRPR16_convert_par2_MPV_fixb.m
% 
% Output:
% tauc [5x1]
%   joint torques required to compensate Coriolis and centrifugal load

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-31 18:39
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauc = S5RPRPR16_coriolisvecJ_fixb_mdp_slag_vp(qJ, qJD, pkin, MDP)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(7,1),zeros(24,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RPRPR16_coriolisvecJ_fixb_mdp_slag_vp: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RPRPR16_coriolisvecJ_fixb_mdp_slag_vp: qJD has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [7 1]), ...
  'S5RPRPR16_coriolisvecJ_fixb_mdp_slag_vp: pkin has to be [7x1] (double)');
assert(isreal(MDP) && all(size(MDP) == [24 1]), ...
  'S5RPRPR16_coriolisvecJ_fixb_mdp_slag_vp: MDP has to be [24x1] (double)'); 

%% Symbolic Calculation
% From coriolisvec_joint_fixb_mdp_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 18:39:35
% EndTime: 2019-12-31 18:39:40
% DurationCPUTime: 1.64s
% Computational Cost: add. (658->216), mult. (1348->310), div. (0->0), fcn. (672->4), ass. (0->109)
t193 = -pkin(1) - pkin(6);
t173 = qJD(1) * t193 + qJD(2);
t191 = cos(qJ(3));
t150 = (pkin(4) * qJD(1) - t173) * t191;
t223 = qJD(4) + t150;
t246 = t173 * t191;
t151 = (-qJD(4) - t246) * qJD(3);
t189 = sin(qJ(3));
t257 = qJD(3) * pkin(3);
t208 = -qJD(4) + t257;
t152 = -t208 - t246;
t206 = -t152 + t246;
t163 = t189 * t173;
t234 = qJD(3) * qJ(4);
t154 = -t163 - t234;
t250 = t154 * t191;
t263 = (t189 * t206 + t250) * qJD(3) + t151 * t189;
t186 = t189 ^ 2;
t187 = t191 ^ 2;
t262 = MDP(8) * (t186 - t187);
t261 = qJ(2) * MDP(6) + MDP(5);
t260 = t189 * pkin(3) + qJ(2);
t192 = -pkin(3) - pkin(7);
t229 = qJD(3) * t192;
t258 = pkin(4) - t193;
t255 = qJ(4) * t191;
t188 = sin(qJ(5));
t227 = qJD(5) * t188;
t222 = qJD(1) * qJD(3);
t212 = t191 * t222;
t190 = cos(qJ(5));
t236 = qJD(1) * t189;
t215 = t190 * t236;
t240 = qJD(5) * t215 + t188 * t212;
t138 = -qJD(3) * t227 + t240;
t254 = t138 * t190;
t141 = (qJD(4) - t150) * qJD(3);
t253 = t141 * t188;
t252 = t141 * t190;
t233 = qJD(3) * t188;
t158 = -t215 + t233;
t235 = qJD(1) * t191;
t177 = qJD(5) + t235;
t249 = t158 * t177;
t248 = t158 * t189;
t231 = qJD(3) * t190;
t160 = t188 * t236 + t231;
t247 = t160 * t177;
t245 = t177 * t191;
t244 = t177 * t192;
t243 = t190 * t191;
t195 = qJD(1) ^ 2;
t242 = t191 * t195;
t194 = qJD(3) ^ 2;
t241 = t193 * t194;
t162 = pkin(3) * t235 + qJ(4) * t236;
t239 = pkin(3) * t236 + qJD(1) * qJ(2);
t237 = t194 + t195;
t232 = qJD(3) * t189;
t230 = qJD(3) * t191;
t228 = qJD(4) * t191;
t226 = qJD(5) * t190;
t149 = -pkin(4) * t236 + t163;
t144 = t149 + t234;
t225 = t144 * qJD(5);
t224 = t190 * MDP(24);
t184 = qJD(1) * qJD(2);
t221 = 0.2e1 * qJD(1);
t219 = t189 * t242;
t213 = t189 * t222;
t218 = pkin(3) * t212 + qJ(4) * t213 + t184;
t217 = t177 * t227;
t216 = t177 * t226;
t214 = pkin(3) * t230 + qJ(4) * t232 + qJD(2);
t211 = MDP(22) * t232;
t210 = qJD(3) * t258;
t201 = (qJD(3) * pkin(7) - qJD(4)) * t191;
t135 = qJD(1) * t201 + t218;
t161 = t173 * t232;
t145 = -pkin(4) * t213 + t161;
t207 = -t135 * t188 + t190 * t145;
t153 = -qJ(4) * t235 + t239;
t164 = -t255 + t260;
t205 = qJD(1) * t164 + t153;
t204 = t177 + t235;
t202 = pkin(7) * t189 - t255;
t136 = t229 + t223;
t143 = qJD(1) * t202 + t239;
t133 = t136 * t190 - t143 * t188;
t134 = t136 * t188 + t143 * t190;
t155 = t202 + t260;
t166 = t258 * t191;
t200 = t155 * t190 + t166 * t188;
t199 = -qJD(1) * t186 + t245;
t198 = t177 * t188;
t197 = t154 * MDP(17) - t160 * MDP(24);
t137 = -qJD(1) * t228 + t218;
t147 = t214 - t228;
t196 = -qJD(1) * t147 - t137 + t241;
t171 = t190 * t212;
t169 = t188 * t213;
t165 = t258 * t189;
t157 = t191 * t210;
t156 = t189 * t210;
t148 = pkin(7) * t235 + t162;
t146 = t153 * t235;
t140 = t201 + t214;
t139 = qJD(5) * t160 - t171;
t1 = [-0.2e1 * t189 * MDP(7) * t212 + 0.2e1 * t222 * t262 + (-t189 * t241 + (qJ(2) * t230 + qJD(2) * t189) * t221) * MDP(12) + (-t191 * t241 + (-qJ(2) * t232 + qJD(2) * t191) * t221) * MDP(13) + t263 * MDP(14) + (t189 * t196 - t205 * t230) * MDP(15) + (t191 * t196 + t205 * t232) * MDP(16) + (t137 * t164 + t147 * t153 - t193 * t263) * MDP(17) + (t138 * t188 * t189 + (t188 * t230 + t189 * t226) * t160) * MDP(18) + ((-t158 * t188 + t160 * t190) * t230 + (t254 - t139 * t188 + (-t158 * t190 - t160 * t188) * qJD(5)) * t189) * MDP(19) + (t189 * t216 + t138 * t191 + (-t160 * t189 + t188 * t199) * qJD(3)) * MDP(20) + (-t189 * t217 - t139 * t191 + (t190 * t199 + t248) * qJD(3)) * MDP(21) - t204 * t211 + ((-t140 * t188 - t190 * t156) * t177 - t157 * t158 - t165 * t139 + (-t144 * t231 + t207) * t191 + (-t134 * t191 - t177 * t200) * qJD(5) + (t188 * t225 - t252 + (-(-t155 * t188 + t190 * t166) * qJD(1) - t133) * qJD(3)) * t189) * MDP(23) + (-t165 * t138 - t157 * t160 + (-(qJD(5) * t166 + t140) * t177 - (qJD(5) * t136 + t135) * t191) * t190 + (-(-qJD(5) * t155 - t156) * t177 + (qJD(3) * t144 + qJD(5) * t143 - t145) * t191) * t188 + (t190 * t225 + t253 + (qJD(1) * t200 + t134) * qJD(3)) * t189) * MDP(24) + 0.2e1 * t261 * t184 + (-MDP(10) * t191 - t189 * MDP(9)) * t194; -t153 * qJD(1) * MDP(17) - t261 * t195 + (t188 * MDP(23) + t224) * t177 * (qJD(5) * t191 + qJD(1)) + ((-MDP(13) + MDP(16)) * t237 + (t158 * MDP(23) - t197) * qJD(3)) * t191 + (-t151 * MDP(17) + t139 * MDP(23) + t138 * MDP(24) + (-MDP(12) + MDP(15)) * t237 + (-MDP(17) * t206 + (MDP(23) * t190 - MDP(24) * t188) * t204) * qJD(3)) * t189; MDP(7) * t219 - t195 * t262 + ((-t154 - t234) * t191 + (t152 + t208) * t189) * qJD(1) * MDP(14) + (t162 * t236 + t146) * MDP(15) + (0.2e1 * qJD(3) * qJD(4) + (-t153 * t189 + t162 * t191) * qJD(1)) * MDP(16) + (-qJ(4) * t151 - qJD(4) * t154 - t153 * t162 + (t250 + (-t152 - t257) * t189) * t173) * MDP(17) + (-t160 * t198 + t254) * MDP(18) + ((-t139 - t247) * t190 + (-t138 + t249) * t188) * MDP(19) + (-t217 + (-t188 * t245 + (t160 - t231) * t189) * qJD(1)) * MDP(20) + (-t216 + t169 + (-t177 * t243 - t248) * qJD(1)) * MDP(21) + t177 * MDP(22) * t236 + (qJ(4) * t139 + t253 - (-t148 * t188 + t149 * t190) * t177 + t223 * t158 + (t144 * t190 - t188 * t244) * qJD(5) + (t144 * t243 + (-t190 * t229 + t133) * t189) * qJD(1)) * MDP(23) + (qJ(4) * t138 + t252 + (t148 * t190 + t149 * t188) * t177 + t223 * t160 + (-t144 * t188 - t190 * t244) * qJD(5) + (-t134 * t189 + (-t144 * t191 + t189 * t229) * t188) * qJD(1)) * MDP(24) + (MDP(13) * t189 * t195 - MDP(12) * t242) * qJ(2); -MDP(15) * t219 + (-t187 * t195 - t194) * MDP(16) + (t146 + t161) * MDP(17) + t169 * MDP(24) + ((-t158 - t215) * MDP(23) + t197) * qJD(3) + (-MDP(23) * t198 - t177 * t224) * t177; t160 * t158 * MDP(18) + (-t158 ^ 2 + t160 ^ 2) * MDP(19) + (t240 + t249) * MDP(20) + (t171 + t247) * MDP(21) - qJD(1) * t211 + (t134 * t177 - t144 * t160 + t207) * MDP(23) + (t133 * t177 - t135 * t190 + t144 * t158 - t145 * t188) * MDP(24) + (-MDP(20) * t233 - MDP(21) * t160 - MDP(23) * t134 - MDP(24) * t133) * qJD(5);];
tauc = t1;
