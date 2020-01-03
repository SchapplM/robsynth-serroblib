% Calculate minimal parameter regressor of Coriolis joint torque vector for
% S4RRRP4
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
%   see S4RRRP4_convert_par2_MPV_fixb.m
% 
% Output:
% tauc [4x1]
%   joint torques required to compensate Coriolis and centrifugal load

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-31 17:15
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauc = S4RRRP4_coriolisvecJ_fixb_mdp_slag_vp(qJ, qJD, pkin, MDP)
%% Coder Information
%#codegen
%$cgargs {zeros(4,1),zeros(4,1),zeros(6,1),zeros(19,1)}
assert(isreal(qJ) && all(size(qJ) == [4 1]), ...
  'S4RRRP4_coriolisvecJ_fixb_mdp_slag_vp: qJ has to be [4x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [4 1]), ...
  'S4RRRP4_coriolisvecJ_fixb_mdp_slag_vp: qJD has to be [4x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [6 1]), ...
  'S4RRRP4_coriolisvecJ_fixb_mdp_slag_vp: pkin has to be [6x1] (double)');
assert(isreal(MDP) && all(size(MDP) == [19 1]), ...
  'S4RRRP4_coriolisvecJ_fixb_mdp_slag_vp: MDP has to be [19x1] (double)'); 

%% Symbolic Calculation
% From coriolisvec_joint_fixb_mdp_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 17:15:43
% EndTime: 2019-12-31 17:15:45
% DurationCPUTime: 0.81s
% Computational Cost: add. (658->141), mult. (1759->200), div. (0->0), fcn. (1112->4), ass. (0->88)
t175 = cos(qJ(2));
t220 = pkin(5) + pkin(6);
t159 = t220 * t175;
t155 = qJD(1) * t159;
t172 = sin(qJ(3));
t141 = t172 * t155;
t173 = sin(qJ(2));
t158 = t220 * t173;
t153 = qJD(1) * t158;
t219 = qJD(2) * pkin(2);
t147 = -t153 + t219;
t174 = cos(qJ(3));
t191 = t174 * t147 - t141;
t152 = t172 * t175 + t174 * t173;
t205 = qJD(1) * t152;
t216 = t205 * qJ(4);
t226 = t216 - t191;
t201 = qJD(1) * qJD(2);
t225 = -0.2e1 * t201;
t224 = t173 * MDP(4);
t198 = qJD(2) * t220;
t186 = qJD(1) * t198;
t148 = t173 * t186;
t223 = (qJD(3) * t147 - t148) * t174;
t222 = (t173 ^ 2 - t175 ^ 2) * MDP(5);
t169 = qJD(2) + qJD(3);
t221 = t205 ^ 2;
t212 = t174 * t175;
t196 = qJD(1) * t212;
t204 = qJD(1) * t173;
t197 = t172 * t204;
t137 = -t196 + t197;
t168 = -t175 * pkin(2) - pkin(1);
t157 = t168 * qJD(1);
t131 = t137 * pkin(3) + qJD(4) + t157;
t218 = t131 * t205;
t217 = t137 * qJ(4);
t215 = t157 * t205;
t214 = t172 * t173;
t176 = qJD(2) ^ 2;
t213 = t173 * t176;
t145 = t174 * t155;
t211 = t175 * t176;
t177 = qJD(1) ^ 2;
t210 = t175 * t177;
t119 = t169 * pkin(3) - t226;
t209 = t119 + t226;
t208 = -t174 * t153 - t141;
t195 = t175 * t201;
t207 = -qJD(3) * t196 - t174 * t195;
t203 = qJD(3) * t172;
t202 = qJD(3) * t174;
t200 = pkin(2) * t204;
t199 = t173 * t219;
t194 = -pkin(2) * t169 - t147;
t130 = t169 * t152;
t128 = t130 * qJD(1);
t193 = t128 * pkin(3) + qJD(1) * t199;
t192 = pkin(1) * t225;
t149 = t175 * t186;
t190 = t172 * t148 - t174 * t149;
t189 = -t172 * t149 - t155 * t203;
t188 = t172 * t153 - t145;
t185 = t169 * t214;
t184 = -t172 * t147 - t145;
t183 = t172 * t158 - t174 * t159;
t182 = t157 * t137 - t189;
t136 = t137 ^ 2;
t181 = t205 * t137 * MDP(11) + (-t136 + t221) * MDP(12) + (-t207 + (t137 - t197) * t169) * MDP(13);
t154 = t173 * t198;
t156 = t175 * t198;
t180 = -t174 * t154 - t172 * t156 - t158 * t202 - t159 * t203;
t179 = t184 * qJD(3) + t190;
t178 = t183 * qJD(3) + t172 * t154 - t174 * t156;
t167 = t174 * pkin(2) + pkin(3);
t151 = -t212 + t214;
t129 = -qJD(2) * t212 - t175 * t202 + t185;
t127 = qJD(1) * t185 + t207;
t126 = -t151 * qJ(4) - t183;
t125 = -t152 * qJ(4) - t174 * t158 - t172 * t159;
t124 = -t216 + t208;
t123 = t188 + t217;
t121 = -t184 - t217;
t116 = t129 * qJ(4) - t152 * qJD(4) + t178;
t115 = -t130 * qJ(4) - t151 * qJD(4) + t180;
t114 = t127 * qJ(4) - qJD(4) * t205 + t179;
t113 = -t128 * qJ(4) - t137 * qJD(4) + t189 + t223;
t1 = [0.2e1 * t195 * t224 + t222 * t225 + MDP(6) * t211 - MDP(7) * t213 + (-pkin(5) * t211 + t173 * t192) * MDP(9) + (pkin(5) * t213 + t175 * t192) * MDP(10) + (-t127 * t152 - t129 * t205) * MDP(11) + (t127 * t151 - t152 * t128 + t129 * t137 - t130 * t205) * MDP(12) + (t168 * t128 + t157 * t130 + (qJD(1) * t151 + t137) * t199) * MDP(16) + (-t168 * t127 - t157 * t129 + 0.2e1 * t205 * t199) * MDP(17) + (-t113 * t151 - t114 * t152 - t115 * t137 - t116 * t205 + t119 * t129 - t121 * t130 + t125 * t127 - t126 * t128) * MDP(18) + (t113 * t126 + t121 * t115 + t114 * t125 + t119 * t116 + t193 * (t151 * pkin(3) + t168) + t131 * (t130 * pkin(3) + t199)) * MDP(19) + (-t129 * MDP(13) - t130 * MDP(14) + t178 * MDP(16) - t180 * MDP(17)) * t169; -t210 * t224 + t177 * t222 + (-t137 * t200 - t215 - t188 * t169 + (t194 * t172 - t145) * qJD(3) + t190) * MDP(16) + (-t205 * t200 + t208 * t169 + (t194 * qJD(3) + t148) * t174 + t182) * MDP(17) + (t167 * t127 + (t121 + t123) * t205 + (-t119 + t124) * t137 + (-t128 * t172 + (-t137 * t174 + t172 * t205) * qJD(3)) * pkin(2)) * MDP(18) + (-pkin(3) * t218 + t114 * t167 - t119 * t123 - t121 * t124 + (-t131 * t204 + t113 * t172 + (-t119 * t172 + t121 * t174) * qJD(3)) * pkin(2)) * MDP(19) + t181 + (t177 * t173 * MDP(9) + MDP(10) * t210) * pkin(1); (-t184 * t169 + t179 - t215) * MDP(16) + (t191 * t169 + t182 - t223) * MDP(17) + (pkin(3) * t127 - t209 * t137) * MDP(18) + (t209 * t121 + (t114 - t218) * pkin(3)) * MDP(19) + t181; (-t136 - t221) * MDP(18) + (t119 * t205 + t121 * t137 + t193) * MDP(19);];
tauc = t1;
