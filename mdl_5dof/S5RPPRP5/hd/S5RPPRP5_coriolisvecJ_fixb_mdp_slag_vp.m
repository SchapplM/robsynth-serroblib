% Calculate minimal parameter regressor of Coriolis joint torque vector for
% S5RPPRP5
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% qJD [5x1]
%   Generalized joint velocities
% pkin [7x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d1,d4,theta2]';
% MDP [22x1]
%   Minimal dynamic parameter vector (fixed base model)
%   see S5RPPRP5_convert_par2_MPV_fixb.m
% 
% Output:
% tauc [5x1]
%   joint torques required to compensate Coriolis and centrifugal load

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-31 17:54
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauc = S5RPPRP5_coriolisvecJ_fixb_mdp_slag_vp(qJ, qJD, pkin, MDP)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(7,1),zeros(22,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RPPRP5_coriolisvecJ_fixb_mdp_slag_vp: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RPPRP5_coriolisvecJ_fixb_mdp_slag_vp: qJD has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [7 1]), ...
  'S5RPPRP5_coriolisvecJ_fixb_mdp_slag_vp: pkin has to be [7x1] (double)');
assert(isreal(MDP) && all(size(MDP) == [22 1]), ...
  'S5RPPRP5_coriolisvecJ_fixb_mdp_slag_vp: MDP has to be [22x1] (double)'); 

%% Symbolic Calculation
% From coriolisvec_joint_fixb_mdp_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 17:53:44
% EndTime: 2019-12-31 17:53:47
% DurationCPUTime: 0.77s
% Computational Cost: add. (634->159), mult. (1594->216), div. (0->0), fcn. (1029->4), ass. (0->76)
t182 = MDP(17) + MDP(19);
t181 = MDP(18) - MDP(21);
t158 = sin(pkin(7));
t156 = t158 ^ 2;
t159 = cos(pkin(7));
t157 = t159 ^ 2;
t196 = t156 + t157;
t202 = (MDP(6) + MDP(9)) * t196;
t161 = sin(qJ(4));
t193 = qJD(1) * t159;
t176 = t161 * t193;
t162 = cos(qJ(4));
t194 = qJD(1) * t158;
t178 = t162 * t194;
t129 = -t176 + t178;
t201 = t129 ^ 2;
t199 = -pkin(6) + qJ(2);
t139 = t199 * t159;
t136 = qJD(1) * t139;
t198 = t136 * t161;
t197 = t159 * MDP(8);
t134 = t158 * t161 + t159 * t162;
t195 = qJD(1) * t134;
t192 = qJD(3) * t158;
t138 = t199 * t158;
t166 = t138 * t162 - t139 * t161;
t105 = t134 * qJD(2) + t166 * qJD(4);
t191 = qJD(4) * t105;
t116 = t138 * t161 + t139 * t162;
t135 = t158 * t162 - t159 * t161;
t106 = -t135 * qJD(2) + t116 * qJD(4);
t190 = qJD(4) * t106;
t189 = qJD(4) * t129;
t188 = qJD(4) * t161;
t187 = qJD(4) * t162;
t186 = t156 * MDP(10);
t145 = qJ(2) * t194 + qJD(3);
t132 = -pkin(6) * t194 + t145;
t113 = t132 * t162 - t198;
t185 = qJD(5) - t113;
t184 = qJD(1) * qJD(2);
t183 = qJD(1) * qJD(3);
t180 = -t159 * pkin(2) - t158 * qJ(3) - pkin(1);
t173 = t162 * t184;
t174 = t161 * t184;
t179 = t132 * t187 + t158 * t174 + t159 * t173;
t177 = t158 * t187;
t175 = qJ(2) * t184;
t172 = t158 * t183;
t171 = t113 + t198;
t131 = t159 * pkin(3) - t180;
t170 = 0.2e1 * t195;
t104 = t132 * t188 + t136 * t187 - t158 * t173 + t159 * t174;
t114 = t132 * t161 + t136 * t162;
t109 = qJD(4) * qJ(5) + t114;
t169 = qJD(4) * t109 - t104;
t124 = -qJD(1) * pkin(1) - pkin(2) * t193 - qJ(3) * t194 + qJD(2);
t167 = t157 * t175;
t118 = pkin(3) * t193 - t124;
t125 = t134 * qJD(4);
t121 = qJD(1) * t125;
t144 = qJD(4) * t176;
t122 = qJD(1) * t177 - t144;
t165 = pkin(4) * t122 + qJ(5) * t121 - qJD(5) * t129;
t164 = qJD(1) ^ 2;
t163 = qJD(4) ^ 2;
t146 = t156 * t175;
t126 = -t159 * t188 + t177;
t112 = pkin(4) * t129 + qJ(5) * t195;
t108 = -qJD(4) * pkin(4) + t185;
t107 = pkin(4) * t134 - qJ(5) * t135 + t131;
t103 = pkin(4) * t126 + qJ(5) * t125 - qJD(5) * t135 + t192;
t102 = pkin(4) * t195 - qJ(5) * t129 + t118;
t101 = (qJD(5) - t198) * qJD(4) + t179;
t100 = t165 + t172;
t1 = [0.2e1 * (t146 + t167) * MDP(7) + 0.2e1 * t172 * t197 + 0.2e1 * t183 * t186 + (0.2e1 * t167 + t146 + (t145 * qJD(2) + (-qJD(1) * t180 - t124) * qJD(3)) * t158) * MDP(11) + (-t121 * t135 - t125 * t129) * MDP(12) + (t121 * t134 - t122 * t135 + t125 * t195 - t126 * t129) * MDP(13) + (t118 * t126 + t122 * t131 + t170 * t192 - t190) * MDP(17) + (-t191 - t118 * t125 - t121 * t131 + (qJD(1) * t135 + t129) * t192) * MDP(18) + (t100 * t134 + t102 * t126 + t103 * t195 + t107 * t122 - t190) * MDP(19) + (-t101 * t134 + t104 * t135 - t105 * t195 + t106 * t129 - t108 * t125 - t109 * t126 - t116 * t122 + t121 * t166) * MDP(20) + (-t100 * t135 + t102 * t125 - t103 * t129 + t107 * t121 + t191) * MDP(21) + (t100 * t107 + t101 * t116 + t102 * t103 - t104 * t166 + t105 * t109 + t106 * t108) * MDP(22) + 0.2e1 * t184 * t202 + (-t125 * MDP(14) - t126 * MDP(15)) * qJD(4); (t195 ^ 2 + t201) * MDP(20) + (t108 * t129 - t109 * t195 - t165) * MDP(22) + t181 * t170 * qJD(4) + ((-qJD(3) - t145) * MDP(11) - qJD(3) * MDP(22) - t182 * t187) * t194 + t182 * (t144 - t189) + (-t202 + (-t157 * MDP(11) - t196 * MDP(7)) * qJ(2)) * t164; -t164 * t186 + t182 * (-t161 * t163 - t194 * t195) + ((-t122 + t189) * MDP(20) + (qJD(4) * t108 + t101) * MDP(22)) * t161 + ((-qJD(4) * t195 + t121) * MDP(20) + t169 * MDP(22) - t181 * t163) * t162 + (-t164 * t197 + ((qJD(2) + t124) * MDP(11) - t102 * MDP(22) - t181 * t129) * qJD(1)) * t158; t144 * MDP(15) + (pkin(4) * t121 - qJ(5) * t122) * MDP(20) + (qJ(5) * t101 - t102 * t112 - t108 * t114 + t185 * t109) * MDP(22) + (-t118 * MDP(17) - t102 * MDP(19) + (t109 - t114) * MDP(20) + t112 * MDP(21) + MDP(13) * t129) * t129 + (t129 * MDP(12) + t118 * MDP(18) - t112 * MDP(19) + (t108 - t185) * MDP(20) - t102 * MDP(21) - MDP(13) * t195) * t195 + ((t129 - t178) * MDP(15) + t171 * MDP(18) + (0.2e1 * qJD(5) - t171) * MDP(21) + t182 * t114) * qJD(4) - t181 * t179 + (-pkin(4) * MDP(22) - t182) * t104; t129 * t195 * MDP(19) + (-t163 - t201) * MDP(21) + (t102 * t129 - t169) * MDP(22);];
tauc = t1;
