% Calculate minimal parameter regressor of Coriolis joint torque vector for
% S5RRPPR1
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% qJD [5x1]
%   Generalized joint velocities
% pkin [9x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d1,d2,d5,theta3,theta4]';
% MDP [18x1]
%   Minimal dynamic parameter vector (fixed base model)
%   see S5RRPPR1_convert_par2_MPV_fixb.m
% 
% Output:
% tauc [5x1]
%   joint torques required to compensate Coriolis and centrifugal load

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-05 18:19
% Revision: 77da58f92bca3eff71542919beafa37024070d86 (2019-12-05)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauc = S5RRPPR1_coriolisvecJ_fixb_mdp_slag_vp(qJ, qJD, pkin, MDP)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(9,1),zeros(18,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RRPPR1_coriolisvecJ_fixb_mdp_slag_vp: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RRPPR1_coriolisvecJ_fixb_mdp_slag_vp: qJD has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S5RRPPR1_coriolisvecJ_fixb_mdp_slag_vp: pkin has to be [9x1] (double)');
assert(isreal(MDP) && all(size(MDP) == [18 1]), ...
  'S5RRPPR1_coriolisvecJ_fixb_mdp_slag_vp: MDP has to be [18x1] (double)'); 

%% Symbolic Calculation
% From coriolisvec_joint_fixb_mdp_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-05 18:18:29
% EndTime: 2019-12-05 18:18:33
% DurationCPUTime: 0.67s
% Computational Cost: add. (601->117), mult. (1272->168), div. (0->0), fcn. (793->8), ass. (0->74)
t167 = sin(pkin(9));
t169 = cos(pkin(9));
t190 = t167 ^ 2 + t169 ^ 2;
t168 = sin(pkin(8));
t172 = sin(qJ(2));
t203 = pkin(1) * qJD(1);
t188 = t172 * t203;
t153 = t168 * t188;
t170 = cos(pkin(8));
t174 = cos(qJ(2));
t202 = pkin(1) * qJD(2);
t155 = t170 * t174 * t202;
t129 = qJD(1) * t155 - qJD(2) * t153;
t166 = qJD(1) + qJD(2);
t121 = qJD(4) * t166 + t129;
t211 = t121 * t190;
t210 = t172 * MDP(5) + t174 * MDP(6);
t171 = sin(qJ(5));
t173 = cos(qJ(5));
t146 = t167 * t173 + t169 * t171;
t132 = t146 * t166;
t209 = t166 * t190;
t187 = t174 * t203;
t139 = t170 * t187 - t153;
t206 = t139 - qJD(4);
t205 = pkin(2) * t170;
t204 = pkin(4) * t169;
t201 = t166 * t167;
t200 = t167 * t171;
t198 = t169 * MDP(8);
t196 = t169 * t173;
t195 = t170 * t172;
t147 = t166 * pkin(2) + t187;
t122 = t147 * t170 - t153;
t177 = qJD(4) - t122;
t184 = -pkin(3) - t204;
t111 = t184 * t166 + t177;
t138 = (t168 * t174 + t195) * t202;
t128 = qJD(1) * t138;
t145 = -t196 + t200;
t140 = t145 * qJD(5);
t194 = -t111 * t140 + t128 * t146;
t141 = t146 * qJD(5);
t193 = t111 * t141 + t128 * t145;
t154 = t170 * t188;
t123 = t168 * t147 + t154;
t159 = pkin(1) * t174 + pkin(2);
t191 = pkin(1) * t195 + t168 * t159;
t156 = t168 * t172 * pkin(1);
t186 = t166 * t200;
t185 = t166 * t196;
t182 = t159 * t170 - t156;
t148 = qJD(5) * t185;
t125 = -qJD(5) * t186 + t148;
t126 = t166 * t141;
t130 = -t185 + t186;
t181 = (-t125 * t145 - t126 * t146 + t130 * t140 - t132 * t141) * MDP(13) + (t125 * t146 - t132 * t140) * MDP(12) + (-t140 * MDP(14) - t141 * MDP(15)) * qJD(5);
t180 = -pkin(3) - t182;
t176 = t190 * (qJ(4) * t166 + t123);
t175 = -qJD(2) * t156 + t155;
t162 = t169 * pkin(7);
t158 = pkin(2) * t168 + qJ(4);
t149 = t184 - t205;
t143 = t158 * t169 + t162;
t142 = (-pkin(7) - t158) * t167;
t137 = t168 * t187 + t154;
t136 = qJ(4) + t191;
t135 = qJD(4) + t175;
t127 = t180 - t204;
t124 = t128 * t167;
t120 = t136 * t169 + t162;
t119 = (-pkin(7) - t136) * t167;
t115 = -pkin(3) * t166 + t177;
t1 = [(-t122 * t138 + t123 * t175 - t128 * t182 + t129 * t191) * MDP(7) + (-t138 * t166 - t128) * t198 + (t138 * t201 + t124) * MDP(9) + (t135 * t209 + t211) * MDP(10) + (t115 * t138 + t128 * t180 + t176 * t135 + t136 * t211) * MDP(11) + (t127 * t126 + t138 * t130 + ((-t119 * t171 - t120 * t173) * qJD(5) - t146 * t135) * qJD(5) + t193) * MDP(17) + (t127 * t125 + t138 * t132 + ((-t119 * t173 + t120 * t171) * qJD(5) + t145 * t135) * qJD(5) + t194) * MDP(18) + t181 + t210 * (-qJD(1) - t166) * t202; (t122 * t137 - t123 * t139 + (-t128 * t170 + t129 * t168) * pkin(2)) * MDP(7) + (t137 * t166 - t128) * t198 + (-t137 * t201 + t124) * MDP(9) + (-t206 * t209 + t211) * MDP(10) + (t128 * (-pkin(3) - t205) - t115 * t137 + t158 * t211 - t206 * t176) * MDP(11) + (t149 * t126 - t137 * t130 + ((-t142 * t171 - t143 * t173) * qJD(5) + t206 * t146) * qJD(5) + t193) * MDP(17) + (t149 * t125 - t137 * t132 + ((-t142 * t173 + t143 * t171) * qJD(5) - t206 * t145) * qJD(5) + t194) * MDP(18) + t181 + t210 * (-qJD(2) + t166) * t203; (-MDP(17) * t141 + MDP(18) * t140) * qJD(5); (-t176 * t166 + t128) * MDP(11) + t148 * MDP(18) - t190 * MDP(10) * t166 ^ 2 + (0.2e1 * t132 * MDP(17) + (-t130 - t186) * MDP(18)) * qJD(5); t132 * t130 * MDP(12) + (-t130 ^ 2 + t132 ^ 2) * MDP(13) + (t148 + (t130 - t186) * qJD(5)) * MDP(14) + (-t111 * t132 - t146 * t121) * MDP(17) + (t111 * t130 + t145 * t121) * MDP(18);];
tauc = t1;
