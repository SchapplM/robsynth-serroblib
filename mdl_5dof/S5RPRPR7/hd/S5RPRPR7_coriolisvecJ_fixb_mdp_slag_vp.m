% Calculate minimal parameter regressor of Coriolis joint torque vector for
% S5RPRPR7
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% qJD [5x1]
%   Generalized joint velocities
% pkin [9x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d1,d3,d5,theta2,theta4]';
% MDP [20x1]
%   Minimal dynamic parameter vector (fixed base model)
%   see S5RPRPR7_convert_par2_MPV_fixb.m
% 
% Output:
% tauc [5x1]
%   joint torques required to compensate Coriolis and centrifugal load

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-31 18:20
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauc = S5RPRPR7_coriolisvecJ_fixb_mdp_slag_vp(qJ, qJD, pkin, MDP)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(9,1),zeros(20,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RPRPR7_coriolisvecJ_fixb_mdp_slag_vp: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RPRPR7_coriolisvecJ_fixb_mdp_slag_vp: qJD has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S5RPRPR7_coriolisvecJ_fixb_mdp_slag_vp: pkin has to be [9x1] (double)');
assert(isreal(MDP) && all(size(MDP) == [20 1]), ...
  'S5RPRPR7_coriolisvecJ_fixb_mdp_slag_vp: MDP has to be [20x1] (double)'); 

%% Symbolic Calculation
% From coriolisvec_joint_fixb_mdp_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 18:19:54
% EndTime: 2019-12-31 18:19:57
% DurationCPUTime: 1.31s
% Computational Cost: add. (1215->200), mult. (2953->290), div. (0->0), fcn. (2016->8), ass. (0->102)
t191 = sin(pkin(8)) * pkin(1) + pkin(6);
t247 = qJ(4) + t191;
t204 = sin(qJ(3));
t265 = MDP(5) * t204;
t206 = cos(qJ(3));
t264 = (t204 ^ 2 - t206 ^ 2) * MDP(6);
t226 = t247 * qJD(1);
t166 = t206 * qJD(2) - t226 * t204;
t167 = t204 * qJD(2) + t226 * t206;
t236 = qJD(1) * qJD(4);
t263 = -t167 * qJD(3) - t204 * t236;
t157 = t166 * qJD(3) + t206 * t236;
t199 = sin(pkin(9));
t201 = cos(pkin(9));
t133 = t157 * t199 - t201 * t263;
t242 = qJD(1) * t206;
t243 = qJD(1) * t204;
t176 = -t199 * t243 + t201 * t242;
t174 = qJD(5) - t176;
t183 = t199 * t206 + t201 * t204;
t178 = t183 * qJD(1);
t190 = pkin(3) * t199 + pkin(7);
t262 = (pkin(3) * t243 + pkin(4) * t178 - pkin(7) * t176 + qJD(5) * t190) * t174 + t133;
t134 = t201 * t157 + t263 * t199;
t260 = qJD(3) * pkin(3);
t161 = t166 + t260;
t253 = t167 * t199;
t138 = t161 * t201 - t253;
t136 = -qJD(3) * pkin(4) - t138;
t228 = qJD(3) * t247;
t168 = qJD(4) * t206 - t204 * t228;
t211 = -qJD(4) * t204 - t206 * t228;
t143 = t201 * t168 + t199 * t211;
t193 = -cos(pkin(8)) * pkin(1) - pkin(2);
t218 = -pkin(3) * t206 + t193;
t213 = t218 * qJD(1);
t175 = qJD(4) + t213;
t144 = -pkin(4) * t176 - pkin(7) * t178 + t175;
t182 = t199 * t204 - t201 * t206;
t152 = pkin(4) * t182 - pkin(7) * t183 + t218;
t180 = t182 * qJD(3);
t181 = t247 * t206;
t229 = t247 * t204;
t155 = t201 * t181 - t199 * t229;
t177 = t183 * qJD(3);
t171 = qJD(1) * t177;
t222 = t133 * t183 - t155 * t171;
t261 = -t136 * t180 - (qJD(5) * t152 + t143) * t174 - (qJD(5) * t144 + t134) * t182 + t222;
t203 = sin(qJ(5));
t240 = qJD(5) * t203;
t237 = qJD(1) * qJD(3);
t231 = t206 * t237;
t232 = t204 * t237;
t172 = -t199 * t232 + t201 * t231;
t205 = cos(qJ(5));
t238 = t205 * qJD(3);
t245 = qJD(5) * t238 + t205 * t172;
t146 = -t178 * t240 + t245;
t259 = t146 * t203;
t258 = t152 * t171;
t250 = t178 * t203;
t163 = -t238 + t250;
t257 = t163 * t174;
t256 = t163 * t178;
t165 = qJD(3) * t203 + t178 * t205;
t255 = t165 * t174;
t254 = t165 * t178;
t252 = t171 * t203;
t251 = t172 * t203;
t159 = t201 * t167;
t207 = qJD(3) ^ 2;
t249 = t204 * t207;
t169 = t205 * t171;
t248 = t206 * t207;
t246 = t146 * t182 + t165 * t177;
t139 = t199 * t161 + t159;
t185 = qJD(1) * t193;
t241 = qJD(5) * t183;
t235 = t204 * t260;
t234 = t183 * t252;
t233 = t183 * t169;
t227 = t174 * t205;
t137 = qJD(3) * pkin(7) + t139;
t132 = t137 * t205 + t144 * t203;
t221 = t137 * t203 - t144 * t205;
t147 = qJD(5) * t165 + t251;
t220 = -t147 * t182 - t163 * t177;
t217 = 0.2e1 * qJD(3) * t185;
t216 = t169 + (t176 * t203 - t240) * t174;
t215 = t180 * t203 - t205 * t241;
t214 = t180 * t205 + t183 * t240;
t141 = t166 * t201 - t253;
t210 = -t190 * t171 + (t136 + t141) * t174;
t192 = -pkin(3) * t201 - pkin(4);
t188 = pkin(3) * t232;
t154 = t181 * t199 + t201 * t229;
t150 = pkin(4) * t177 + pkin(7) * t180 + t235;
t148 = pkin(4) * t171 - pkin(7) * t172 + t188;
t145 = t205 * t148;
t142 = t168 * t199 - t201 * t211;
t140 = t166 * t199 + t159;
t1 = [0.2e1 * t231 * t265 - 0.2e1 * t237 * t264 + MDP(7) * t248 - MDP(8) * t249 + (-t191 * t248 + t204 * t217) * MDP(10) + (t191 * t249 + t206 * t217) * MDP(11) + (-t134 * t182 + t138 * t180 - t139 * t177 + t142 * t178 + t143 * t176 + t154 * t172 + t222) * MDP(12) + (t133 * t154 + t134 * t155 - t138 * t142 + t139 * t143 + (t175 + t213) * t235) * MDP(13) + (t146 * t183 * t205 - t214 * t165) * MDP(14) + (-(-t163 * t205 - t165 * t203) * t180 + (-t259 - t147 * t205 + (t163 * t203 - t165 * t205) * qJD(5)) * t183) * MDP(15) + (-t214 * t174 + t233 + t246) * MDP(16) + (t215 * t174 + t220 - t234) * MDP(17) + (t171 * t182 + t174 * t177) * MDP(18) + (-t221 * t177 + t142 * t163 + t145 * t182 + t154 * t147 + (t150 * t174 + t258 + (t136 * t183 - t137 * t182 - t155 * t174) * qJD(5)) * t205 + t261 * t203) * MDP(19) + (-t132 * t177 + t142 * t165 + t154 * t146 + (-(-qJD(5) * t155 + t150) * t174 - t258 - (-qJD(5) * t137 + t148) * t182 - t136 * t241) * t203 + t261 * t205) * MDP(20); (-t171 * t183 + t172 * t182 - t176 * t180 + t177 * t178) * MDP(12) + (t133 * t182 + t134 * t183 - t138 * t177 - t139 * t180) * MDP(13) + (-t220 - t234) * MDP(19) + (-t233 + t246) * MDP(20) + (-MDP(10) * t204 - MDP(11) * t206) * t207 + (t215 * MDP(19) + t214 * MDP(20)) * t174; ((t139 - t140) * t178 + (t138 - t141) * t176 + (-t171 * t199 - t172 * t201) * pkin(3)) * MDP(12) + (t138 * t140 - t139 * t141 + (-t133 * t201 + t134 * t199 - t175 * t243) * pkin(3)) * MDP(13) + (t165 * t227 + t259) * MDP(14) + ((t146 - t257) * t205 + (-t147 - t255) * t203) * MDP(15) + (t174 * t227 + t252 - t254) * MDP(16) + (t216 + t256) * MDP(17) - t174 * t178 * MDP(18) + (-t140 * t163 + t192 * t147 + t178 * t221 + t210 * t203 - t262 * t205) * MDP(19) + (t132 * t178 - t140 * t165 + t192 * t146 + t262 * t203 + t210 * t205) * MDP(20) + (-t206 * t265 + t264) * qJD(1) ^ 2 + (-MDP(10) * t243 - MDP(11) * t242) * t185; (-t176 ^ 2 - t178 ^ 2) * MDP(12) + (t138 * t178 - t139 * t176 + t188) * MDP(13) + (t216 - t256) * MDP(19) + (-t174 ^ 2 * t205 - t252 - t254) * MDP(20); t165 * t163 * MDP(14) + (-t163 ^ 2 + t165 ^ 2) * MDP(15) + (t245 + t257) * MDP(16) + (-t251 + t255) * MDP(17) + t171 * MDP(18) + (t132 * t174 - t134 * t203 - t136 * t165 + t145) * MDP(19) + (-t134 * t205 + t136 * t163 - t148 * t203 - t174 * t221) * MDP(20) + (-MDP(16) * t250 - t165 * MDP(17) - t132 * MDP(19) + t221 * MDP(20)) * qJD(5);];
tauc = t1;
