% Calculate Coriolis joint torque vector for
% S5PRPRP2
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% qJD [5x1]
%   Generalized joint velocities
% pkin [8x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d2,d4,theta1,theta3]';
% MDP [18x1]
%   Minimal dynamic parameter vector (fixed base model)
%   see S5PRPRP2_convert_par2_MPV_fixb.m
% 
% Output:
% tauc [5x1]
%   joint torques required to compensate Coriolis and centrifugal load

% Quelle: HybrDyn-Toolbox
% Datum: 2021-01-15 15:05
% Revision: 24b2e7d74a0c1a3b64fa2f8f5ad758691ad61af3 (2021-01-15)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauc = S5PRPRP2_coriolisvecJ_fixb_mdp_slag_vp(qJ, qJD, pkin, MDP)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(8,1),zeros(18,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5PRPRP2_coriolisvecJ_fixb_mdp_slag_vp: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5PRPRP2_coriolisvecJ_fixb_mdp_slag_vp: qJD has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S5PRPRP2_coriolisvecJ_fixb_mdp_slag_vp: pkin has to be [8x1] (double)');
assert(isreal(MDP) && all(size(MDP) == [18 1]), ...
  'S5PRPRP2_coriolisvecJ_fixb_mdp_slag_vp: MDP has to be [18x1] (double)'); 

%% Symbolic Calculation
% From coriolisvec_joint_fixb_mdp_matlab.m
% OptimizationMode: 2
% StartTime: 2021-01-15 15:04:54
% EndTime: 2021-01-15 15:04:57
% DurationCPUTime: 0.82s
% Computational Cost: add. (533->140), mult. (1326->214), div. (0->0), fcn. (786->4), ass. (0->78)
t137 = sin(qJ(4));
t133 = t137 ^ 2;
t138 = cos(qJ(4));
t134 = t138 ^ 2;
t191 = MDP(9) * (t133 - t134);
t135 = sin(pkin(8));
t136 = cos(pkin(8));
t122 = -pkin(3) * t136 - pkin(6) * t135 - pkin(2);
t153 = qJ(5) * t135 - t122;
t182 = t136 * t138;
t190 = -qJ(3) * t182 + t153 * t137;
t113 = t122 * qJD(2) + qJD(3);
t176 = qJD(2) * t136;
t120 = qJ(3) * t176 + qJD(1) * t135;
t144 = -t113 * t137 - t120 * t138;
t177 = qJD(2) * t135;
t160 = t137 * t177;
t102 = -qJ(5) * t160 - t144;
t172 = qJD(4) * t138;
t161 = qJ(5) * t172;
t169 = qJD(2) * qJD(3);
t155 = t136 * t169;
t173 = qJD(4) * t137;
t163 = -t113 * t172 + t120 * t173 - t138 * t155;
t171 = qJD(5) * t137;
t96 = (-t161 - t171) * t177 - t163;
t168 = qJD(2) * qJD(4);
t154 = t135 * t168;
t148 = t137 * t154;
t123 = qJ(5) * t148;
t175 = qJD(3) * t137;
t156 = t136 * t175;
t184 = t135 * t138;
t141 = -qJD(5) * t184 - t156;
t142 = t144 * qJD(4);
t97 = t141 * qJD(2) + t123 + t142;
t152 = t138 * t113 - t120 * t137;
t159 = t138 * t177;
t101 = -qJ(5) * t159 + t152;
t128 = -qJD(4) + t176;
t98 = -pkin(4) * t128 + t101;
t189 = t96 * t137 + t97 * t138 + (t102 * t138 - t137 * t98) * qJD(4);
t132 = t136 ^ 2;
t131 = t135 ^ 2;
t188 = 0.2e1 * t131;
t187 = -t101 + t98;
t186 = t128 * t136;
t139 = qJD(2) ^ 2;
t185 = t131 * t139;
t183 = t136 * t137;
t174 = qJD(3) * t138;
t181 = t122 * t172 + t136 * t174;
t147 = t138 * t154;
t112 = pkin(4) * t147 + t135 * t169;
t180 = t131 + t132;
t178 = qJD(2) * t131;
t118 = (pkin(4) * t137 + qJ(3)) * t135;
t130 = t136 * qJD(1);
t109 = qJD(2) * t118 + qJD(5) - t130;
t170 = qJD(5) + t109;
t167 = MDP(13) + MDP(15);
t166 = MDP(14) + MDP(16);
t165 = pkin(4) * t185;
t162 = qJ(3) * t173;
t158 = t135 * t173;
t157 = t135 * t172;
t151 = t180 * qJD(2);
t150 = t138 * t131 * t137 * MDP(8);
t149 = t128 * t158;
t119 = qJ(3) * t177 - t130;
t143 = t119 * t135 + t120 * t136;
t121 = t136 * t148;
t117 = (pkin(4) * t172 + qJD(3)) * t135;
t116 = t128 * t159;
t105 = -t153 * t138 + (-qJ(3) * t137 - pkin(4)) * t136;
t100 = qJD(4) * t190 + t141;
t99 = -t135 * t171 + (-qJ(3) * t183 - qJ(5) * t184) * qJD(4) + t181;
t1 = [t167 * (t128 - t176) * t157 + t166 * (t121 - t149) + (-t112 * t136 + (-t137 * t97 + t138 * t96 + (-t102 * t137 - t138 * t98) * qJD(4)) * t135) * MDP(18); (t121 + t149) * MDP(10) + (t128 + t176) * MDP(11) * t157 + ((t186 + (t188 + t132) * qJD(2)) * t175 + ((t113 * t136 + t122 * t128) * t137 + ((t178 + t186) * qJ(3) + t143) * t138) * qJD(4)) * MDP(13) + ((-t136 * t162 + t181) * t128 - t163 * t136 - t119 * t158 + (-t162 + 0.2e1 * t174) * t178) * MDP(14) + (-t100 * t128 - t97 * t136 + (t109 * t172 + t112 * t137 + (t117 * t137 + t118 * t172) * qJD(2)) * t135) * MDP(15) + (t99 * t128 + t96 * t136 + (-t109 * t173 + t112 * t138 + (t117 * t138 - t118 * t173) * qJD(2)) * t135) * MDP(16) + ((-t100 * t138 - t137 * t99 + (t105 * t137 + t138 * t190) * qJD(4)) * qJD(2) - t189) * t135 * MDP(17) + (t100 * t98 + t102 * t99 + t105 * t97 + t109 * t117 + t112 * t118 - t190 * t96) * MDP(18) + (t188 * t191 - 0.2e1 * t150) * t168 + (0.2e1 * MDP(6) * t151 + (qJ(3) * t151 + t143) * MDP(7)) * qJD(3); t189 * MDP(18) - t180 * MDP(6) * t139 + (-t143 * MDP(7) + (-t102 * t182 - t109 * t135 + t98 * t183) * MDP(18)) * qJD(2) + (t167 * t137 + t166 * t138) * (-t128 ^ 2 - t185); t139 * t150 - t185 * t191 + (-t116 - t147) * MDP(11) + (t144 * t128 + t142 + (-t119 * t184 - t156) * qJD(2)) * MDP(13) + (-t152 * t128 + t163) * MDP(14) + (-t102 * t128 + t123 + (-qJD(4) * t113 - t155) * t137 + (-qJD(4) * t120 - t137 * t165 - t170 * t177) * t138) * MDP(15) + (-t134 * t165 - t101 * t128 + (t170 * t137 + t161) * t177 + t163) * MDP(16) + (t187 * t102 + (-t109 * t159 + t97) * pkin(4)) * MDP(18) + ((-qJD(4) - t128) * MDP(10) + t119 * MDP(14) + (pkin(4) * qJD(4) - t187) * MDP(17)) * t160; -t116 * MDP(15) + t112 * MDP(18) + (-t133 - t134) * MDP(17) * t185 + ((MDP(15) * qJD(4) + MDP(18) * t98) * t138 + ((-qJD(4) + t128) * MDP(16) + t102 * MDP(18)) * t137) * t177;];
tauc = t1;
