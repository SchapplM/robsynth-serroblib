% Calculate Coriolis joint torque vector for
% S5RPPRP4
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% qJD [5x1]
%   Generalized joint velocities
% pkin [7x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d1,d4,theta3]';
% MDP [18x1]
%   Minimal dynamic parameter vector (fixed base model)
%   see S5RPPRP4_convert_par2_MPV_fixb.m
% 
% Output:
% tauc [5x1]
%   joint torques required to compensate Coriolis and centrifugal load

% Quelle: HybrDyn-Toolbox
% Datum: 2021-01-15 17:13
% Revision: 24b2e7d74a0c1a3b64fa2f8f5ad758691ad61af3 (2021-01-15)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauc = S5RPPRP4_coriolisvecJ_fixb_mdp_slag_vp(qJ, qJD, pkin, MDP)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(7,1),zeros(18,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RPPRP4_coriolisvecJ_fixb_mdp_slag_vp: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RPPRP4_coriolisvecJ_fixb_mdp_slag_vp: qJD has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [7 1]), ...
  'S5RPPRP4_coriolisvecJ_fixb_mdp_slag_vp: pkin has to be [7x1] (double)');
assert(isreal(MDP) && all(size(MDP) == [18 1]), ...
  'S5RPPRP4_coriolisvecJ_fixb_mdp_slag_vp: MDP has to be [18x1] (double)'); 

%% Symbolic Calculation
% From coriolisvec_joint_fixb_mdp_matlab.m
% OptimizationMode: 2
% StartTime: 2021-01-15 17:13:18
% EndTime: 2021-01-15 17:13:20
% DurationCPUTime: 0.52s
% Computational Cost: add. (515->121), mult. (983->177), div. (0->0), fcn. (465->4), ass. (0->74)
t128 = sin(pkin(7));
t134 = qJD(4) ^ 2;
t135 = qJD(1) ^ 2;
t183 = t128 * (t134 + t135);
t131 = sin(qJ(4));
t126 = t131 ^ 2;
t132 = cos(qJ(4));
t127 = t132 ^ 2;
t182 = (t126 - t127) * MDP(9);
t181 = qJ(2) * MDP(6) + MDP(5);
t180 = MDP(17) * (t126 + t127);
t129 = cos(pkin(7));
t162 = qJD(2) * t129;
t145 = -qJD(5) + t162;
t139 = t145 * t132;
t161 = qJD(4) * t131;
t154 = qJ(5) * t161;
t158 = qJD(4) * qJD(3);
t133 = -pkin(1) - pkin(2);
t114 = t133 * qJD(1) + qJD(2);
t165 = qJD(1) * qJ(2);
t107 = t128 * t114 + t129 * t165;
t103 = -qJD(1) * pkin(6) + t107;
t99 = t103 * t161;
t91 = t132 * t158 - t99 + (t139 + t154) * qJD(1);
t174 = t103 * t132;
t140 = -qJD(3) * t131 - t174;
t160 = qJD(1) * qJD(2);
t159 = qJD(1) * qJD(4);
t152 = t132 * t159;
t164 = qJD(1) * t131;
t170 = qJ(5) * t152 + qJD(5) * t164;
t173 = t129 * t131;
t92 = t140 * qJD(4) - t160 * t173 + t170;
t177 = qJD(4) * pkin(4);
t125 = t132 * qJD(3);
t150 = -t103 * t131 + t125;
t96 = qJ(5) * t164 + t150;
t93 = t96 + t177;
t163 = qJD(1) * t132;
t97 = -qJ(5) * t163 - t140;
t179 = -(t131 * t97 + t132 * t93) * qJD(4) - t92 * t131 + t91 * t132;
t178 = t93 - t96;
t176 = t132 * t97;
t172 = t132 * t135;
t169 = t129 * qJ(2) + t128 * t133;
t111 = -pkin(6) + t169;
t171 = qJ(5) - t111;
t157 = MDP(13) + MDP(15);
t156 = MDP(14) + MDP(16);
t153 = t131 * t159;
t116 = t128 * t160;
t149 = -t128 * qJ(2) + t129 * t133;
t110 = pkin(3) - t149;
t108 = pkin(4) * t132 + t110;
t106 = t114 * t129 - t128 * t165;
t102 = qJD(1) * pkin(3) - t106;
t98 = pkin(4) * t163 + qJD(5) + t102;
t151 = -qJD(1) * t108 - t98;
t148 = qJD(4) * t171;
t109 = -pkin(4) * t153 + t116;
t113 = -pkin(4) * t161 + qJD(2) * t128;
t147 = qJD(1) * t113 + t109;
t146 = 0.2e1 * t152;
t142 = -t131 * t93 + t176;
t141 = t106 * t128 - t107 * t129;
t138 = (t102 - t162) * qJD(1);
t137 = -t111 * t134 + 0.2e1 * t116;
t136 = qJD(4) * (-qJD(1) * t110 - t102 - t162);
t105 = t171 * t132;
t104 = t171 * t131;
t95 = -t145 * t131 + t132 * t148;
t94 = t131 * t148 + t139;
t1 = [((-t128 * t149 + t129 * t169) * qJD(1) - t141) * qJD(2) * MDP(7) + t131 * MDP(8) * t146 - 0.2e1 * t159 * t182 + (t131 * t136 + t137 * t132) * MDP(13) + (-t137 * t131 + t132 * t136) * MDP(14) + (t147 * t132 + (t151 * t131 + t95) * qJD(4)) * MDP(15) + (-t147 * t131 + (t151 * t132 - t94) * qJD(4)) * MDP(16) + ((t131 * t95 - t132 * t94 + (t104 * t132 - t105 * t131) * qJD(4)) * qJD(1) - t179) * MDP(17) + (t104 * t92 - t105 * t91 + t108 * t109 + t113 * t98 + t93 * t95 + t94 * t97) * MDP(18) + 0.2e1 * t181 * t160 + (-t132 * MDP(10) + t131 * MDP(11)) * t134; t156 * (t129 * t146 + t131 * t183) + t157 * (0.2e1 * t129 * t153 - t132 * t183) + (t129 * t180 - t181) * t135 + t141 * MDP(7) * qJD(1) + (-t109 * t129 + t179 * t128 + (-t128 * t98 - t129 * t176 + t93 * t173) * qJD(1)) * MDP(18); (t142 * qJD(4) + t91 * t131 + t92 * t132) * MDP(18) + (-t157 * t131 - t156 * t132) * t134; t135 * t182 + (t99 + t150 * qJD(4) + (t138 - t158) * t132) * MDP(14) + ((t97 - t174) * qJD(4) + t170) * MDP(15) + (-pkin(4) * t126 * t135 + t99 + (t96 - t125) * qJD(4) + (-t154 + (-t145 + t98) * t132) * qJD(1)) * MDP(16) + (t177 - t178) * MDP(17) * t163 + (t178 * t97 + (t98 * t164 + t92) * pkin(4)) * MDP(18) + (-MDP(8) * t172 + MDP(13) * t138 + (pkin(4) * t172 - t158 + (t98 - t162) * qJD(1)) * MDP(15)) * t131; t116 * MDP(18) - t135 * t180 + (t142 * MDP(18) + (-0.2e1 * t132 * MDP(16) + (-MDP(18) * pkin(4) - 0.2e1 * MDP(15)) * t131) * qJD(4)) * qJD(1);];
tauc = t1;
