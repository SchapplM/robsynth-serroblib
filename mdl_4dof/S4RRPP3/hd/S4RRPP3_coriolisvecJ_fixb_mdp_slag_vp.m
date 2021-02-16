% Calculate Coriolis joint torque vector for
% S4RRPP3
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [4x1]
%   Generalized joint coordinates (joint angles)
% qJD [4x1]
%   Generalized joint velocities
% pkin [6x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,d1,d2,theta3]';
% MDP [18x1]
%   Minimal dynamic parameter vector (fixed base model)
%   see S4RRPP3_convert_par2_MPV_fixb.m
% 
% Output:
% tauc [4x1]
%   joint torques required to compensate Coriolis and centrifugal load

% Quelle: HybrDyn-Toolbox
% Datum: 2021-01-15 10:36
% Revision: d12c3222fdeb2c5f3b3c8fa5751e113be2fc3aae (2021-01-15)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauc = S4RRPP3_coriolisvecJ_fixb_mdp_slag_vp(qJ, qJD, pkin, MDP)
%% Coder Information
%#codegen
%$cgargs {zeros(4,1),zeros(4,1),zeros(6,1),zeros(18,1)}
assert(isreal(qJ) && all(size(qJ) == [4 1]), ...
  'S4RRPP3_coriolisvecJ_fixb_mdp_slag_vp: qJ has to be [4x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [4 1]), ...
  'S4RRPP3_coriolisvecJ_fixb_mdp_slag_vp: qJD has to be [4x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [6 1]), ...
  'S4RRPP3_coriolisvecJ_fixb_mdp_slag_vp: pkin has to be [6x1] (double)');
assert(isreal(MDP) && all(size(MDP) == [18 1]), ...
  'S4RRPP3_coriolisvecJ_fixb_mdp_slag_vp: MDP has to be [18x1] (double)'); 

%% Symbolic Calculation
% From coriolisvec_joint_fixb_mdp_matlab.m
% OptimizationMode: 2
% StartTime: 2021-01-15 10:36:07
% EndTime: 2021-01-15 10:36:10
% DurationCPUTime: 0.75s
% Computational Cost: add. (643->151), mult. (1714->202), div. (0->0), fcn. (1035->4), ass. (0->71)
t157 = sin(pkin(6));
t158 = cos(pkin(6));
t159 = sin(qJ(2));
t160 = cos(qJ(2));
t141 = t157 * t160 + t158 * t159;
t181 = qJD(1) * t141;
t190 = t181 * MDP(12);
t189 = (t159 ^ 2 - t160 ^ 2) * MDP(5);
t174 = MDP(11) + MDP(15);
t130 = t181 ^ 2;
t187 = -qJ(3) - pkin(5);
t168 = qJD(2) * t187;
t128 = qJD(3) * t160 + t159 * t168;
t122 = t128 * qJD(1);
t164 = -qJD(3) * t159 + t160 * t168;
t123 = t164 * qJD(1);
t102 = t157 * t122 - t158 * t123;
t146 = t187 * t160;
t170 = t187 * t159;
t117 = -t146 * t157 - t158 * t170;
t186 = t102 * t117;
t144 = qJD(1) * t146;
t185 = t144 * t157;
t137 = t158 * t144;
t184 = t158 * t160;
t183 = t160 * MDP(4);
t103 = t158 * t122 + t157 * t123;
t143 = qJD(1) * t170;
t139 = qJD(2) * pkin(2) + t143;
t114 = t157 * t139 - t137;
t154 = -pkin(2) * t160 - pkin(1);
t180 = qJD(1) * t154;
t179 = qJD(1) * t159;
t178 = qJD(2) * t159;
t177 = t160 * MDP(10);
t116 = t143 * t158 + t185;
t176 = qJD(4) - t116;
t175 = qJD(2) * qJD(4);
t173 = MDP(12) - MDP(17);
t172 = 0.2e1 * qJD(1);
t171 = qJD(1) * t184;
t169 = qJD(1) * t178;
t167 = qJD(2) * t116 - t103;
t113 = t139 * t158 + t185;
t133 = t141 * qJD(2);
t125 = qJD(1) * t133;
t147 = t157 * t169;
t126 = qJD(2) * t171 - t147;
t150 = pkin(2) * t169;
t166 = pkin(3) * t125 - qJ(4) * t126 + t150;
t145 = qJD(3) + t180;
t107 = t128 * t157 - t158 * t164;
t108 = t158 * t128 + t157 * t164;
t118 = -t158 * t146 + t157 * t170;
t131 = t157 * t179 - t171;
t163 = t102 * t141 + t107 * t181 - t108 * t131 + t117 * t126 - t118 * t125;
t161 = qJD(2) ^ 2;
t153 = -pkin(2) * t158 - pkin(3);
t151 = pkin(2) * t157 + qJ(4);
t140 = t157 * t159 - t184;
t136 = qJD(2) * t184 - t157 * t178;
t115 = t143 * t157 - t137;
t112 = pkin(3) * t140 - qJ(4) * t141 + t154;
t110 = qJD(2) * qJ(4) + t114;
t109 = -qJD(2) * pkin(3) + qJD(4) - t113;
t106 = pkin(2) * t179 + pkin(3) * t181 + qJ(4) * t131;
t105 = pkin(3) * t131 - qJ(4) * t181 + t145;
t101 = t175 + t103;
t100 = pkin(2) * t178 + pkin(3) * t133 - qJ(4) * t136 - qJD(4) * t141;
t98 = -qJD(4) * t181 + t166;
t1 = [(t125 * t154 + t133 * t145) * MDP(11) + (t126 * t154 + t136 * t145) * MDP(12) + (-t103 * t140 - t113 * t136 - t114 * t133 + t163) * MDP(13) + (t103 * t118 - t107 * t113 + t108 * t114 + t186) * MDP(14) + (t100 * t131 + t105 * t133 + t112 * t125 + t140 * t98) * MDP(15) + (-t101 * t140 + t109 * t136 - t110 * t133 + t163) * MDP(16) + (-t100 * t181 - t105 * t136 - t112 * t126 - t141 * t98) * MDP(17) + (t100 * t105 + t101 * t118 + t107 * t109 + t108 * t110 + t112 * t98 + t186) * MDP(18) + (t160 * MDP(6) - t159 * MDP(7) + (MDP(10) * t159 - MDP(9) * t160) * pkin(5)) * t161 + (-t173 * t108 - t174 * t107 + (-pkin(1) * t177 - t189) * t172 + ((-pkin(1) * MDP(9) + t183) * t172 + ((qJD(1) * t140 + t131) * MDP(11) + 0.2e1 * t190 + (t145 + t180) * MDP(14)) * pkin(2)) * t159) * qJD(2); t167 * MDP(12) + (t113 * t115 - t114 * t116) * MDP(14) + (-t125 * t151 + t153 * t126) * MDP(16) + (-t167 + 0.2e1 * t175) * MDP(17) + (t101 * t151 + t102 * t153 - t105 * t106 - t109 * t115 + t110 * t176) * MDP(18) + (-t145 * MDP(11) + (t114 - t115) * MDP(13) - t105 * MDP(15) + (t110 - t115) * MDP(16) + t106 * MDP(17)) * t181 + (t145 * MDP(12) + (-t113 + t116) * MDP(13) - t106 * MDP(15) + (t109 - t176) * MDP(16) - t105 * MDP(17)) * t131 + ((-t125 * t157 - t126 * t158) * MDP(13) + (-t102 * t158 + t103 * t157) * MDP(14) + (-MDP(11) * t131 - MDP(14) * t145 - t190) * t179) * pkin(2) + (-t159 * t183 + t189 + (t159 * MDP(9) + t177) * pkin(1)) * qJD(1) ^ 2 + t174 * (qJD(2) * t115 - t102); (t113 * t181 + t114 * t131 + t150) * MDP(14) + (t110 * t131 + (-qJD(4) - t109) * t181 + t166) * MDP(18) + t173 * (-t147 + (-t131 + t171) * qJD(2)) + (MDP(13) + MDP(16)) * (-t131 ^ 2 - t130) + 0.2e1 * t174 * t181 * qJD(2); t181 * t131 * MDP(15) + (-t147 + (t131 + t171) * qJD(2)) * MDP(16) + (-t130 - t161) * MDP(17) + (-qJD(2) * t110 + t105 * t181 + t102) * MDP(18);];
tauc = t1;
