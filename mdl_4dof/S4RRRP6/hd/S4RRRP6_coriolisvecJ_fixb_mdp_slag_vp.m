% Calculate minimal parameter regressor of Coriolis joint torque vector for
% S4RRRP6
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [4x1]
%   Generalized joint coordinates (joint angles)
% qJD [4x1]
%   Generalized joint velocities
% pkin [6x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,d1,d2,d3]';
% MDP [19x1]
%   Minimal dynamic parameter vector (fixed base model)
%   see S4RRRP6_convert_par2_MPV_fixb.m
% 
% Output:
% tauc [4x1]
%   joint torques required to compensate Coriolis and centrifugal load

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-31 17:19
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauc = S4RRRP6_coriolisvecJ_fixb_mdp_slag_vp(qJ, qJD, pkin, MDP)
%% Coder Information
%#codegen
%$cgargs {zeros(4,1),zeros(4,1),zeros(6,1),zeros(19,1)}
assert(isreal(qJ) && all(size(qJ) == [4 1]), ...
  'S4RRRP6_coriolisvecJ_fixb_mdp_slag_vp: qJ has to be [4x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [4 1]), ...
  'S4RRRP6_coriolisvecJ_fixb_mdp_slag_vp: qJD has to be [4x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [6 1]), ...
  'S4RRRP6_coriolisvecJ_fixb_mdp_slag_vp: pkin has to be [6x1] (double)');
assert(isreal(MDP) && all(size(MDP) == [19 1]), ...
  'S4RRRP6_coriolisvecJ_fixb_mdp_slag_vp: MDP has to be [19x1] (double)'); 

%% Symbolic Calculation
% From coriolisvec_joint_fixb_mdp_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 17:19:13
% EndTime: 2019-12-31 17:19:16
% DurationCPUTime: 1.22s
% Computational Cost: add. (725->203), mult. (1871->309), div. (0->0), fcn. (1072->4), ass. (0->99)
t204 = (qJD(1) * qJD(2));
t244 = -2 * t204;
t180 = sin(qJ(2));
t243 = t180 * MDP(4);
t177 = t180 ^ 2;
t182 = cos(qJ(2));
t242 = (-t182 ^ 2 + t177) * MDP(5);
t181 = cos(qJ(3));
t179 = sin(qJ(3));
t212 = qJD(2) * t179;
t214 = qJD(1) * t180;
t157 = t181 * t214 + t212;
t241 = t157 ^ 2;
t240 = pkin(5) * t179;
t239 = -qJ(4) - pkin(6);
t238 = qJ(4) * t180;
t209 = qJD(3) * t179;
t199 = t180 * t209;
t196 = t182 * t204;
t203 = qJD(2) * qJD(3);
t218 = (t196 + t203) * t181;
t137 = qJD(1) * t199 - t218;
t237 = t137 * t179;
t213 = qJD(1) * t182;
t171 = -qJD(3) + t213;
t236 = t157 * t171;
t165 = -qJD(2) * pkin(2) + pkin(5) * t214;
t235 = t165 * t179;
t234 = t165 * t181;
t233 = t171 * t179;
t232 = t171 * t181;
t161 = -pkin(2) * t182 - pkin(6) * t180 - pkin(1);
t149 = t161 * qJD(1);
t231 = t179 * t149;
t230 = t179 * t182;
t229 = t180 * t181;
t183 = qJD(2) ^ 2;
t228 = t180 * t183;
t227 = t181 * t182;
t226 = t182 * t183;
t184 = qJD(1) ^ 2;
t225 = t182 * t184;
t175 = pkin(5) * t213;
t166 = qJD(2) * pkin(6) + t175;
t133 = t181 * t149 - t166 * t179;
t127 = -qJ(4) * t157 + t133;
t126 = -pkin(3) * t171 + t127;
t224 = t126 - t127;
t189 = pkin(2) * t180 - pkin(6) * t182;
t160 = t189 * qJD(2);
t150 = qJD(1) * t160;
t197 = t180 * t204;
t190 = pkin(5) * t197;
t223 = -t181 * t150 - t179 * t190;
t159 = t189 * qJD(1);
t145 = t179 * t159;
t194 = qJD(3) * t239;
t207 = qJD(4) * t181;
t222 = t179 * t194 + t207 - t145 - (-pkin(5) * t229 - qJ(4) * t230) * qJD(1);
t186 = pkin(3) * t180 - qJ(4) * t227;
t201 = t179 * t214;
t217 = pkin(5) * t201 + t181 * t159;
t221 = -t186 * qJD(1) - qJD(4) * t179 + t181 * t194 - t217;
t208 = qJD(3) * t181;
t220 = t179 * t160 + t161 * t208;
t211 = qJD(2) * t180;
t219 = t181 * t160 + t211 * t240;
t172 = pkin(5) * t227;
t216 = t179 * t161 + t172;
t210 = qJD(2) * t182;
t206 = t165 * qJD(3);
t205 = t181 * qJD(2);
t202 = pkin(3) * t179 + pkin(5);
t200 = t179 * t210;
t198 = t180 * t208;
t195 = t179 * t203;
t138 = t195 + (t198 + t200) * qJD(1);
t130 = pkin(3) * t138 + pkin(5) * t196;
t193 = pkin(1) * t244;
t155 = t201 - t205;
t192 = t155 + t205;
t191 = -t157 + t212;
t134 = t166 * t181 + t231;
t188 = qJD(1) * t177 - t171 * t182;
t187 = -t149 * t208 - t179 * t150 + t166 * t209;
t185 = -t134 * qJD(3) - t223;
t164 = t239 * t181;
t163 = t239 * t179;
t154 = t181 * t161;
t152 = t155 ^ 2;
t136 = pkin(3) * t155 + qJD(4) + t165;
t135 = -t179 * t238 + t216;
t132 = -qJ(4) * t229 + t154 + (-pkin(3) - t240) * t182;
t128 = -qJ(4) * t155 + t134;
t125 = (-pkin(5) * qJD(2) - qJ(4) * qJD(3)) * t229 + (-qJD(4) * t180 + (-pkin(5) * qJD(3) - qJ(4) * qJD(2)) * t182) * t179 + t220;
t124 = -t180 * t207 + t186 * qJD(2) + (-t172 + (-t161 + t238) * t179) * qJD(3) + t219;
t123 = -qJ(4) * t138 - qJD(4) * t155 - t181 * t190 - t187;
t122 = pkin(3) * t197 + qJ(4) * t137 - qJD(4) * t157 + t185;
t1 = [0.2e1 * t196 * t243 + t242 * t244 + MDP(6) * t226 - MDP(7) * t228 + (-pkin(5) * t226 + t180 * t193) * MDP(9) + (pkin(5) * t228 + t182 * t193) * MDP(10) + (-t137 * t229 + (t182 * t205 - t199) * t157) * MDP(11) + ((-t155 * t181 - t157 * t179) * t210 + (t237 - t138 * t181 + (t155 * t179 - t157 * t181) * qJD(3)) * t180) * MDP(12) + (t171 * t199 + t137 * t182 + (t157 * t180 + t188 * t181) * qJD(2)) * MDP(13) + (t171 * t198 + t138 * t182 + (-t155 * t180 - t188 * t179) * qJD(2)) * MDP(14) + (-t171 - t213) * MDP(15) * t211 + (-(-t161 * t209 + t219) * t171 + (t181 * t206 + pkin(5) * t138 + (t154 * qJD(1) + t133) * qJD(2)) * t180 + ((pkin(5) * t155 + t235) * qJD(2) + (t231 + (pkin(5) * t171 + t166) * t181) * qJD(3) + t223) * t182) * MDP(16) + ((-pkin(5) * t182 * t209 + t220) * t171 - t187 * t182 + (-pkin(5) * t137 - t179 * t206) * t180 + ((pkin(5) * t157 + t234) * t182 + (-pkin(5) * t232 - t216 * qJD(1) - t134) * t180) * qJD(2)) * MDP(17) + (-t124 * t157 - t125 * t155 + t132 * t137 - t135 * t138 + (-t126 * t181 - t128 * t179) * t210 + (-t122 * t181 - t123 * t179 + (t126 * t179 - t128 * t181) * qJD(3)) * t180) * MDP(18) + (t122 * t132 + t123 * t135 + t126 * t124 + t128 * t125 + t136 * t202 * t210 + (t136 * pkin(3) * t208 + t130 * t202) * t180) * MDP(19); -t225 * t243 + t184 * t242 + (-t157 * t232 - t237) * MDP(11) + ((t171 * t155 - t137) * t181 + (-t138 + t236) * t179) * MDP(12) + (-t171 * t208 + (t171 * t227 + t191 * t180) * qJD(1)) * MDP(13) + (t171 * t209 + (-t171 * t230 + t192 * t180) * qJD(1)) * MDP(14) + t171 * MDP(15) * t214 + (-pkin(2) * t138 + t217 * t171 + (pkin(6) * t232 + t235) * qJD(3) + ((-pkin(6) * t212 - t133) * t180 + (-t192 * pkin(5) - t235) * t182) * qJD(1)) * MDP(16) + (pkin(2) * t137 - t145 * t171 + (-pkin(6) * t233 + t234) * qJD(3) + (-t165 * t227 + (-pkin(6) * t205 + t134) * t180 + (t171 * t229 + t191 * t182) * pkin(5)) * qJD(1)) * MDP(17) + (t137 * t163 + t138 * t164 - t221 * t157 - t222 * t155 + (t171 * t126 + t123) * t181 + (t171 * t128 - t122) * t179) * MDP(18) + (-t123 * t164 + t122 * t163 + t130 * (-pkin(3) * t181 - pkin(2)) + (-pkin(3) * t233 - t175) * t136 + t222 * t128 + t221 * t126) * MDP(19) + (t184 * t180 * MDP(9) + MDP(10) * t225) * pkin(1); (-t152 + t241) * MDP(12) + t218 * MDP(13) + (-t195 - t236) * MDP(14) + (-t134 * t171 - t157 * t165 + t185) * MDP(16) + (-t133 * t171 + t187) * MDP(17) + t224 * MDP(19) * t128 + (t137 * MDP(18) + (-t136 * t157 + t122) * MDP(19)) * pkin(3) + (t157 * MDP(11) - t171 * MDP(13) + t165 * MDP(17) - t224 * MDP(18)) * t155 + (-MDP(14) * t200 + ((-t179 * MDP(13) - t181 * MDP(14)) * qJD(3) + (pkin(5) * t181 * MDP(17) + MDP(15)) * qJD(2)) * t180) * qJD(1); (-t152 - t241) * MDP(18) + (t126 * t157 + t128 * t155 + t130) * MDP(19);];
tauc = t1;
