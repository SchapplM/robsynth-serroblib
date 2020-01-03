% Calculate minimal parameter regressor of Coriolis joint torque vector for
% S5RPRRP5
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% qJD [5x1]
%   Generalized joint velocities
% pkin [8x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d1,d3,d4,theta2]';
% MDP [18x1]
%   Minimal dynamic parameter vector (fixed base model)
%   see S5RPRRP5_convert_par2_MPV_fixb.m
% 
% Output:
% tauc [5x1]
%   joint torques required to compensate Coriolis and centrifugal load

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-31 18:41
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauc = S5RPRRP5_coriolisvecJ_fixb_mdp_slag_vp(qJ, qJD, pkin, MDP)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(8,1),zeros(18,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RPRRP5_coriolisvecJ_fixb_mdp_slag_vp: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RPRRP5_coriolisvecJ_fixb_mdp_slag_vp: qJD has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S5RPRRP5_coriolisvecJ_fixb_mdp_slag_vp: pkin has to be [8x1] (double)');
assert(isreal(MDP) && all(size(MDP) == [18 1]), ...
  'S5RPRRP5_coriolisvecJ_fixb_mdp_slag_vp: MDP has to be [18x1] (double)'); 

%% Symbolic Calculation
% From coriolisvec_joint_fixb_mdp_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 18:40:58
% EndTime: 2019-12-31 18:41:01
% DurationCPUTime: 0.62s
% Computational Cost: add. (810->140), mult. (1649->181), div. (0->0), fcn. (842->6), ass. (0->73)
t163 = sin(qJ(3));
t165 = cos(qJ(3));
t206 = pkin(1) * sin(pkin(8));
t180 = qJD(1) * t206;
t173 = qJD(3) * t180;
t151 = cos(pkin(8)) * pkin(1) + pkin(2);
t148 = t151 * qJD(1);
t190 = qJD(3) * t148;
t132 = -t163 * t173 + t165 * t190;
t213 = qJD(2) * qJD(4) + t132;
t162 = sin(qJ(4));
t158 = t162 ^ 2;
t164 = cos(qJ(4));
t159 = t164 ^ 2;
t194 = t162 * t164;
t212 = MDP(8) * t194 - (t158 - t159) * MDP(9);
t135 = t148 * t163 + t165 * t180;
t157 = qJD(1) + qJD(3);
t129 = pkin(7) * t157 + t135;
t202 = t129 * t162;
t119 = qJD(2) * t164 - t202;
t211 = qJD(5) - t119;
t116 = -qJD(4) * pkin(4) + t211;
t120 = qJD(2) * t162 + t129 * t164;
t117 = qJD(4) * qJ(5) + t120;
t209 = -t151 * t165 + t163 * t206;
t208 = (t158 + t159) * t157;
t182 = MDP(13) + MDP(15);
t207 = -MDP(14) + MDP(17);
t205 = pkin(3) * t157;
t166 = qJD(4) ^ 2;
t204 = pkin(7) * t166;
t203 = MDP(18) * pkin(7);
t201 = t135 * t157;
t170 = t151 * t163 + t165 * t206;
t138 = t170 * qJD(3);
t200 = t138 * t157;
t140 = pkin(7) + t170;
t199 = t140 * t166;
t145 = -pkin(4) * t164 - qJ(5) * t162 - pkin(3);
t198 = t145 * t157;
t134 = t148 * t165 - t163 * t180;
t188 = qJD(4) * t162;
t193 = t134 * t188 + t164 * t201;
t192 = t213 * t164;
t133 = t163 * t190 + t165 * t173;
t187 = qJD(4) * t164;
t186 = qJD(5) * t162;
t185 = t140 * MDP(18);
t184 = t166 * MDP(11);
t109 = (qJD(5) - t202) * qJD(4) + t192;
t111 = t129 * t187 + t213 * t162;
t178 = t109 * t164 + t111 * t162 + t116 * t187;
t176 = t119 + t202;
t130 = t145 + t209;
t137 = t209 * qJD(3);
t175 = t130 * t157 + t137;
t174 = (-pkin(3) + t209) * t157 + t137;
t172 = pkin(4) * t162 - qJ(5) * t164;
t171 = t199 + t200;
t112 = (t172 * qJD(4) - t186) * t157 + t133;
t141 = pkin(4) * t188 - qJ(5) * t187 - t186;
t169 = -t141 * t157 - t112 - t204;
t118 = t138 + t141;
t168 = -t118 * t157 - t112 - t199;
t128 = -t134 - t205;
t167 = t166 * t164 * MDP(10) + (t128 * t187 + t133 * t162) * MDP(14) + 0.2e1 * t212 * qJD(4) * t157;
t156 = t157 ^ 2;
t142 = t172 * t157;
t121 = t128 * t188;
t115 = -t134 + t198;
t113 = t115 * t188;
t1 = [(-t133 - t200) * MDP(6) + (t137 * t157 - t132) * MDP(7) + t121 * MDP(13) + t113 * MDP(15) + (-t208 * t137 + t178) * MDP(16) + (t112 * t130 + t115 * t118) * MDP(18) + (-t184 + t171 * MDP(14) + t168 * MDP(17) + (t111 * t140 - t116 * t137) * MDP(18) + (t174 * MDP(13) + t175 * MDP(15) + (-MDP(16) - t185) * t117) * qJD(4)) * t162 + ((-t133 - t171) * MDP(13) + t168 * MDP(15) + (t109 * t140 - t117 * t137) * MDP(18) + (t174 * MDP(14) + (-t115 - t175) * MDP(17) + t116 * t185) * qJD(4)) * t164 + t167; (t109 * t162 - t111 * t164 + (t116 * t162 + t117 * t164) * qJD(4)) * MDP(18) + (-t182 * t162 + t207 * t164) * t166; (-t133 + t201) * MDP(6) + (t134 * t157 - t132) * MDP(7) + (t121 + t193) * MDP(13) + (t113 + t193) * MDP(15) + (-t208 * t134 + t178) * MDP(16) + (t112 * t145 + (-t135 + t141) * t115) * MDP(18) + (-t184 + (-t201 + t204) * MDP(14) + (t169 + t201) * MDP(17) + (pkin(7) * t111 - t116 * t134) * MDP(18) + ((-MDP(13) * pkin(3) + MDP(15) * t145) * t157 + (-MDP(16) - t203) * t117) * qJD(4)) * t162 + ((-t133 - t204) * MDP(13) + t169 * MDP(15) + (pkin(7) * t109 - t117 * t134) * MDP(18) + ((t134 - t205) * MDP(14) + (-t115 - t134 - t198) * MDP(17) + t116 * t203) * qJD(4)) * t164 + t167; (qJ(5) * t109 - t115 * t142 - t116 * t120 + t211 * t117) * MDP(18) - t212 * t156 + (t176 * MDP(14) + (0.2e1 * qJD(5) - t176) * MDP(17) + t182 * t120) * qJD(4) + ((-t128 * MDP(13) - t115 * MDP(15) + t142 * MDP(17)) * t162 + (-t128 * MDP(14) + t142 * MDP(15) + t115 * MDP(17)) * t164) * t157 + t207 * t192 + (-pkin(4) * MDP(18) - t182) * t111; -t156 * MDP(15) * t194 + (-t156 * t158 - t166) * MDP(17) + (t115 * t157 * t162 - t117 * qJD(4) + t111) * MDP(18);];
tauc = t1;
