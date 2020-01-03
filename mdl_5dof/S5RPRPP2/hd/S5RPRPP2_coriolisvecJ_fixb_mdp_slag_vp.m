% Calculate minimal parameter regressor of Coriolis joint torque vector for
% S5RPRPP2
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% qJD [5x1]
%   Generalized joint velocities
% pkin [7x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d1,d3,theta2]';
% MDP [19x1]
%   Minimal dynamic parameter vector (fixed base model)
%   see S5RPRPP2_convert_par2_MPV_fixb.m
% 
% Output:
% tauc [5x1]
%   joint torques required to compensate Coriolis and centrifugal load

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-31 18:11
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauc = S5RPRPP2_coriolisvecJ_fixb_mdp_slag_vp(qJ, qJD, pkin, MDP)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(7,1),zeros(19,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RPRPP2_coriolisvecJ_fixb_mdp_slag_vp: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RPRPP2_coriolisvecJ_fixb_mdp_slag_vp: qJD has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [7 1]), ...
  'S5RPRPP2_coriolisvecJ_fixb_mdp_slag_vp: pkin has to be [7x1] (double)');
assert(isreal(MDP) && all(size(MDP) == [19 1]), ...
  'S5RPRPP2_coriolisvecJ_fixb_mdp_slag_vp: MDP has to be [19x1] (double)'); 

%% Symbolic Calculation
% From coriolisvec_joint_fixb_mdp_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 18:11:13
% EndTime: 2019-12-31 18:11:15
% DurationCPUTime: 0.66s
% Computational Cost: add. (518->141), mult. (1139->188), div. (0->0), fcn. (517->4), ass. (0->77)
t166 = -MDP(12) - MDP(16);
t188 = -MDP(10) + t166;
t132 = sin(pkin(7)) * pkin(1) + pkin(6);
t127 = t132 * qJD(1);
t146 = sin(qJ(3));
t124 = t146 * t127;
t147 = cos(qJ(3));
t116 = t147 * qJD(2) - t124;
t187 = qJD(4) - t116;
t105 = qJ(5) * qJD(1) * t146 + t116;
t171 = qJD(4) - t105;
t184 = pkin(3) + pkin(4);
t100 = -t184 * qJD(3) + t171;
t142 = t146 ^ 2;
t143 = t147 ^ 2;
t186 = (t142 - t143) * MDP(6);
t117 = t146 * qJD(2) + t147 * t127;
t177 = qJD(1) * t147;
t106 = -qJ(5) * t177 + t117;
t141 = qJD(3) * qJ(4);
t103 = t106 + t141;
t111 = t141 + t117;
t107 = -qJD(3) * pkin(3) + t187;
t185 = MDP(11) - MDP(14);
t183 = qJ(4) * t147;
t182 = qJ(5) - t132;
t168 = qJD(2) * qJD(3);
t175 = qJD(3) * t147;
t114 = t127 * t175 + t146 * t168;
t169 = qJD(1) * qJD(3);
t163 = t147 * t169;
t174 = qJD(4) * t146;
t181 = qJ(4) * t163 + qJD(1) * t174;
t136 = t147 * t168;
t140 = qJD(3) * qJD(4);
t180 = t136 + 0.2e1 * t140;
t133 = -cos(pkin(7)) * pkin(1) - pkin(2);
t153 = qJ(4) * t146 - t133;
t115 = t184 * t147 + t153;
t178 = qJD(1) * t115;
t120 = -pkin(3) * t147 - t153;
t112 = qJD(1) * t120;
t128 = qJD(1) * t133;
t176 = qJD(3) * t146;
t173 = qJD(5) * t146;
t172 = qJD(5) * t147;
t101 = qJD(5) + t178;
t170 = qJD(5) + t101;
t167 = MDP(10) + MDP(12);
t165 = MDP(14) + MDP(17);
t164 = t146 * t169;
t122 = t182 * t147;
t160 = t101 + t178;
t102 = -t184 * t164 + t181;
t152 = -t184 * t146 + t183;
t108 = t152 * qJD(3) + t174;
t159 = qJD(1) * t108 + t102;
t158 = 0.2e1 * t112;
t113 = pkin(3) * t164 - t181;
t154 = pkin(3) * t146 - t183;
t119 = t154 * qJD(3) - t174;
t157 = -qJD(1) * t119 - t113;
t156 = 0.2e1 * t128;
t155 = t132 * MDP(15) + MDP(13);
t104 = -t127 * t176 + t136 + t140;
t151 = (t100 * t146 + t103 * t147) * MDP(19);
t150 = qJD(1) ^ 2;
t149 = qJD(3) ^ 2;
t130 = qJ(5) * t164;
t126 = t154 * qJD(1);
t121 = t182 * t146;
t118 = t152 * qJD(1);
t110 = -qJD(3) * t122 - t173;
t109 = t182 * t176 - t172;
t99 = (-qJ(5) * t175 - t173) * qJD(1) + t114;
t98 = -qJD(1) * t172 + t104 + t130;
t1 = [(t112 * t119 + t113 * t120) * MDP(15) + (t100 * t110 + t101 * t108 + t102 * t115 + t103 * t109 - t121 * t99 - t122 * t98) * MDP(19) + (t157 * MDP(12) + t159 * MDP(16) + (-qJD(1) * t109 - t98) * MDP(18) + (-t167 * t132 + MDP(7)) * t149 + t155 * t104) * t147 + (t157 * MDP(14) + t159 * MDP(17) + (-qJD(1) * t110 - t99) * MDP(18) + (t185 * t132 - MDP(8)) * t149 + t155 * t114) * t146 + (-t110 * MDP(16) + t109 * MDP(17) - 0.2e1 * qJD(1) * t186 + (t156 * MDP(11) - t158 * MDP(14) + t160 * MDP(17) + (qJD(1) * t121 - t100) * MDP(18) + t155 * t107) * t147 + (0.2e1 * MDP(5) * t177 + t156 * MDP(10) + t158 * MDP(12) - t160 * MDP(16) + (-qJD(1) * t122 + t103) * MDP(18) - t155 * t111) * t146) * qJD(3); (t104 * t146 - t114 * t147) * MDP(15) + (t146 * t98 - t147 * t99) * MDP(19) + ((t107 * t146 + t111 * t147) * MDP(15) + t151) * qJD(3) + ((-MDP(11) + t165) * t147 + t188 * t146) * t149; -t136 * MDP(11) + t180 * MDP(14) + (qJ(4) * t104 - t107 * t117 + t187 * t111 - t112 * t126) * MDP(15) + (t130 + t180) * MDP(17) + (qJ(4) * t98 - t100 * t106 - t101 * t118 + t171 * t103 - t184 * t99) * MDP(19) + (-t146 * t147 * MDP(5) + t186) * t150 + (t106 * MDP(16) + (-t105 - t124) * MDP(17) + t167 * t117 + t185 * (t116 + t124)) * qJD(3) + ((-t128 * MDP(10) - t112 * MDP(12) + t126 * MDP(14) + t170 * MDP(16) - t118 * MDP(17)) * t146 + (-t128 * MDP(11) + t126 * MDP(12) + t112 * MDP(14) + (qJ(5) * qJD(3) - t118) * MDP(16) - t170 * MDP(17)) * t147) * qJD(1) + (-pkin(3) * MDP(15) + t188) * t114; (-qJD(3) * t111 + t114) * MDP(15) + (-qJ(5) * t163 - qJD(3) * t103 + t114) * MDP(19) + t165 * (-t142 * t150 - t149) + (t166 * t150 * t147 + (t112 * MDP(15) - t170 * MDP(19)) * qJD(1)) * t146; t181 * MDP(19) + (-t142 - t143) * MDP(18) * t150 + (t151 + (0.2e1 * t147 * MDP(17) + (-t184 * MDP(19) - 0.2e1 * MDP(16)) * t146) * qJD(3)) * qJD(1);];
tauc = t1;
