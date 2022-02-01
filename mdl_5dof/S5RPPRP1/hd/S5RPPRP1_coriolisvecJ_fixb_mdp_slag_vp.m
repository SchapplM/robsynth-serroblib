% Calculate Coriolis joint torque vector for
% S5RPPRP1
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% qJD [5x1]
%   Generalized joint velocities
% pkin [8x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d1,d4,theta2,theta3]';
% MDP [18x1]
%   Minimal dynamic parameter vector (fixed base model)
%   see S5RPPRP1_convert_par2_MPV_fixb.m
% 
% Output:
% tauc [5x1]
%   joint torques required to compensate Coriolis and centrifugal load

% Quelle: HybrDyn-Toolbox
% Datum: 2022-01-23 09:13
% Revision: fd3771346c4aea32fdeb66112c511235427c26a7 (2022-01-20)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauc = S5RPPRP1_coriolisvecJ_fixb_mdp_slag_vp(qJ, qJD, pkin, MDP)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(8,1),zeros(18,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RPPRP1_coriolisvecJ_fixb_mdp_slag_vp: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RPPRP1_coriolisvecJ_fixb_mdp_slag_vp: qJD has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S5RPPRP1_coriolisvecJ_fixb_mdp_slag_vp: pkin has to be [8x1] (double)');
assert(isreal(MDP) && all(size(MDP) == [18 1]), ...
  'S5RPPRP1_coriolisvecJ_fixb_mdp_slag_vp: MDP has to be [18x1] (double)'); 

%% Symbolic Calculation
% From coriolisvec_joint_fixb_mdp_matlab.m
% OptimizationMode: 2
% StartTime: 2022-01-23 09:12:47
% EndTime: 2022-01-23 09:12:48
% DurationCPUTime: 0.85s
% Computational Cost: add. (661->143), mult. (1551->219), div. (0->0), fcn. (914->6), ass. (0->81)
t145 = sin(qJ(4));
t139 = t145 ^ 2;
t146 = cos(qJ(4));
t140 = t146 ^ 2;
t200 = MDP(9) * (t139 - t140);
t134 = sin(pkin(7)) * pkin(1) + qJ(3);
t141 = sin(pkin(8));
t143 = cos(pkin(8));
t123 = -cos(pkin(7)) * pkin(1) - pkin(3) * t143 - pkin(6) * t141 - pkin(2);
t161 = qJ(5) * t141 - t123;
t192 = t143 * t146;
t199 = -t134 * t192 + t161 * t145;
t180 = qJD(4) * t146;
t170 = qJ(5) * t180;
t117 = t123 * qJD(1) + qJD(3);
t128 = t134 * qJD(1);
t120 = qJD(2) * t141 + t128 * t143;
t177 = qJD(1) * qJD(3);
t163 = t143 * t177;
t181 = qJD(4) * t145;
t171 = -t117 * t180 + t120 * t181 - t146 * t163;
t179 = qJD(5) * t145;
t186 = qJD(1) * t141;
t100 = (-t170 - t179) * t186 - t171;
t176 = qJD(1) * qJD(4);
t162 = t141 * t176;
t156 = t145 * t162;
t127 = qJ(5) * t156;
t183 = qJD(3) * t145;
t167 = t143 * t183;
t194 = t141 * t146;
t149 = -qJD(5) * t194 - t167;
t152 = -t117 * t145 - t120 * t146;
t150 = t152 * qJD(4);
t101 = t149 * qJD(1) + t127 + t150;
t160 = t146 * t117 - t120 * t145;
t168 = t146 * t186;
t105 = -qJ(5) * t168 + t160;
t185 = qJD(1) * t143;
t133 = -qJD(4) + t185;
t102 = -pkin(4) * t133 + t105;
t184 = qJD(1) * t145;
t169 = t141 * t184;
t106 = -qJ(5) * t169 - t152;
t198 = -(t102 * t145 - t106 * t146) * qJD(4) + t100 * t145 + t101 * t146;
t138 = t143 ^ 2;
t137 = t141 ^ 2;
t197 = 0.2e1 * t137;
t196 = t133 * t143;
t147 = qJD(1) ^ 2;
t195 = t137 * t147;
t193 = t143 * t145;
t191 = t102 - t105;
t182 = qJD(3) * t146;
t190 = t123 * t180 + t143 * t182;
t155 = t146 * t162;
t122 = pkin(4) * t155 + t141 * t177;
t189 = t137 + t138;
t187 = qJD(1) * t137;
t136 = t143 * qJD(2);
t111 = qJD(5) - t136 + (pkin(4) * t184 + t128) * t141;
t178 = qJD(5) + t111;
t175 = MDP(13) + MDP(15);
t174 = MDP(14) + MDP(16);
t173 = pkin(4) * t195;
t166 = t134 * t181;
t165 = t141 * t181;
t164 = t141 * t180;
t159 = qJD(1) * t189;
t158 = t146 * t137 * t145 * MDP(8);
t157 = t133 * t165;
t119 = t128 * t141 - t136;
t151 = t119 * t141 + t120 * t143;
t126 = t143 * t156;
t125 = (pkin(4) * t180 + qJD(3)) * t141;
t124 = t133 * t168;
t121 = (pkin(4) * t145 + t134) * t141;
t107 = -t161 * t146 + (-t134 * t145 - pkin(4)) * t143;
t104 = t199 * qJD(4) + t149;
t103 = -t141 * t179 + (-qJ(5) * t194 - t134 * t193) * qJD(4) + t190;
t1 = [(t126 + t157) * MDP(10) + (t133 + t185) * MDP(11) * t164 + ((t196 + (t197 + t138) * qJD(1)) * t183 + ((t117 * t143 + t123 * t133) * t145 + ((t187 + t196) * t134 + t151) * t146) * qJD(4)) * MDP(13) + ((-t143 * t166 + t190) * t133 - t171 * t143 - t119 * t165 + (-t166 + 0.2e1 * t182) * t187) * MDP(14) + (-t101 * t143 - t104 * t133 + (t111 * t180 + t122 * t145 + (t121 * t180 + t125 * t145) * qJD(1)) * t141) * MDP(15) + (t100 * t143 + t103 * t133 + (-t111 * t181 + t122 * t146 + (-t121 * t181 + t125 * t146) * qJD(1)) * t141) * MDP(16) + ((-t103 * t145 - t104 * t146 + (t107 * t145 + t146 * t199) * qJD(4)) * qJD(1) - t198) * t141 * MDP(17) + (-t100 * t199 + t101 * t107 + t102 * t104 + t103 * t106 + t111 * t125 + t121 * t122) * MDP(18) + (t197 * t200 - 0.2e1 * t158) * t176 + (0.2e1 * MDP(6) * t159 + (t134 * t159 + t151) * MDP(7)) * qJD(3); t175 * (t133 - t185) * t164 + t174 * (t126 - t157) + (-t122 * t143 + (t100 * t146 - t101 * t145 + (-t102 * t146 - t106 * t145) * qJD(4)) * t141) * MDP(18); t198 * MDP(18) - t189 * MDP(6) * t147 + (-t151 * MDP(7) + (t102 * t193 - t106 * t192 - t111 * t141) * MDP(18)) * qJD(1) + (t175 * t145 + t174 * t146) * (-t133 ^ 2 - t195); t147 * t158 - t195 * t200 + (-t124 - t155) * MDP(11) + (t152 * t133 + t150 + (-t119 * t194 - t167) * qJD(1)) * MDP(13) + (-t160 * t133 + t171) * MDP(14) + (-t106 * t133 + t127 + (-qJD(4) * t117 - t163) * t145 + (-qJD(4) * t120 - t145 * t173 - t178 * t186) * t146) * MDP(15) + (-t140 * t173 - t105 * t133 + (t178 * t145 + t170) * t186 + t171) * MDP(16) + (t191 * t106 + (-t111 * t168 + t101) * pkin(4)) * MDP(18) + ((-qJD(4) - t133) * MDP(10) + t119 * MDP(14) + (pkin(4) * qJD(4) - t191) * MDP(17)) * t169; -t124 * MDP(15) + t122 * MDP(18) + (-t139 - t140) * MDP(17) * t195 + ((MDP(15) * qJD(4) + MDP(18) * t102) * t146 + ((-qJD(4) + t133) * MDP(16) + t106 * MDP(18)) * t145) * t186;];
tauc = t1;
