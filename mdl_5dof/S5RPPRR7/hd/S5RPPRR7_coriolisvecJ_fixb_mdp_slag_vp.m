% Calculate minimal parameter regressor of Coriolis joint torque vector for
% S5RPPRR7
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% qJD [5x1]
%   Generalized joint velocities
% pkin [8x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d1,d4,d5,theta2]';
% MDP [21x1]
%   Minimal dynamic parameter vector (fixed base model)
%   see S5RPPRR7_convert_par2_MPV_fixb.m
% 
% Output:
% tauc [5x1]
%   joint torques required to compensate Coriolis and centrifugal load

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-31 18:00
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauc = S5RPPRR7_coriolisvecJ_fixb_mdp_slag_vp(qJ, qJD, pkin, MDP)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(8,1),zeros(21,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RPPRR7_coriolisvecJ_fixb_mdp_slag_vp: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RPPRR7_coriolisvecJ_fixb_mdp_slag_vp: qJD has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S5RPPRR7_coriolisvecJ_fixb_mdp_slag_vp: pkin has to be [8x1] (double)');
assert(isreal(MDP) && all(size(MDP) == [21 1]), ...
  'S5RPPRR7_coriolisvecJ_fixb_mdp_slag_vp: MDP has to be [21x1] (double)'); 

%% Symbolic Calculation
% From coriolisvec_joint_fixb_mdp_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 17:59:55
% EndTime: 2019-12-31 17:59:58
% DurationCPUTime: 0.86s
% Computational Cost: add. (519->163), mult. (1105->245), div. (0->0), fcn. (619->6), ass. (0->84)
t130 = cos(qJ(4));
t186 = MDP(8) * t130;
t124 = t130 ^ 2;
t128 = sin(qJ(4));
t185 = (t128 ^ 2 - t124) * MDP(9);
t184 = 2 * qJD(3);
t127 = sin(qJ(5));
t118 = -cos(pkin(8)) * pkin(1) - pkin(2) - pkin(6);
t110 = t118 * qJD(1) + qJD(3);
t137 = qJD(2) * t128 - t110 * t130;
t96 = -qJD(4) * pkin(4) + t137;
t183 = t127 * t96;
t129 = cos(qJ(5));
t182 = t129 * t96;
t103 = qJD(2) * t130 + t110 * t128;
t99 = qJD(4) * t103;
t181 = t99 * t127;
t180 = t99 * t129;
t155 = t129 * qJD(4);
t121 = qJD(5) * t155;
t158 = qJD(5) * t127;
t147 = t130 * t158;
t133 = -t128 * t155 - t147;
t100 = t133 * qJD(1) + t121;
t166 = t129 * t130;
t146 = qJD(1) * t166;
t161 = qJD(4) * t127;
t114 = t146 + t161;
t159 = qJD(4) * t130;
t179 = t100 * t128 + t114 * t159;
t178 = t100 * t127;
t152 = qJD(1) * qJD(4);
t145 = t128 * t152;
t117 = t127 * t145;
t101 = qJD(5) * t114 - t117;
t177 = t101 * t128;
t120 = sin(pkin(8)) * pkin(1) + qJ(3);
t108 = pkin(4) * t128 - pkin(7) * t130 + t120;
t105 = t108 * qJD(1);
t176 = t105 * t127;
t175 = t105 * t129;
t156 = t127 * qJD(1);
t150 = t130 * t156;
t112 = t150 - t155;
t162 = qJD(1) * t128;
t119 = qJD(5) + t162;
t174 = t112 * t119;
t173 = t112 * t130;
t172 = t114 * t119;
t116 = qJD(1) * t120;
t171 = t116 * MDP(7);
t170 = t119 * t127;
t169 = t127 * t128;
t131 = qJD(4) ^ 2;
t168 = t128 * t131;
t167 = t129 * t119;
t165 = t130 * t131;
t132 = qJD(1) ^ 2;
t163 = -t131 - t132;
t160 = qJD(4) * t128;
t157 = qJD(5) * t129;
t154 = t130 * MDP(13);
t153 = t130 * MDP(19);
t151 = t124 * t156;
t149 = t119 * t161;
t148 = t119 * t157;
t144 = qJD(4) * t153;
t140 = pkin(4) * t130 + pkin(7) * t128;
t111 = t140 * qJD(4) + qJD(3);
t106 = t111 * qJD(1);
t98 = t137 * qJD(4);
t143 = t129 * t106 + t127 * t98;
t97 = qJD(4) * pkin(7) + t103;
t142 = t118 * t119 + t97;
t141 = t119 * t147;
t93 = -t127 * t97 + t175;
t94 = t129 * t97 + t176;
t139 = -t106 * t127 + t129 * t98;
t138 = qJD(1) * t124 - t119 * t128;
t136 = t128 * t149 - t130 * t148;
t135 = -pkin(7) * t159 + t128 * t96;
t134 = (-t129 * MDP(20) + t127 * MDP(21)) * t119;
t115 = t140 * qJD(1);
t1 = [qJD(1) * MDP(6) * t184 + t171 * t184 - 0.2e1 * t145 * t186 + 0.2e1 * t152 * t185 - MDP(10) * t168 - MDP(11) * t165 + (t116 * t159 - t118 * t168 + (t120 * t159 + t128 * t184) * qJD(1)) * MDP(13) + (-t116 * t160 - t118 * t165 + (-t120 * t160 + t130 * t184) * qJD(1)) * MDP(14) + (t100 * t166 + t133 * t114) * MDP(15) + ((t112 * t129 + t114 * t127) * t160 + (-t178 - t101 * t129 + (t112 * t127 - t114 * t129) * qJD(5)) * t130) * MDP(16) + (t138 * t155 - t141 + t179) * MDP(17) + (-t177 + (-t151 - t173) * qJD(4) + t136) * MDP(18) + (t119 + t162) * t144 + ((-t108 * t158 + t111 * t129) * t119 + ((t112 * t118 - t183) * qJD(4) + (-t142 * t129 - t176) * qJD(5) + t143) * t128 + (t96 * t157 - t118 * t101 + t181 + (-t118 * t170 + (t108 * t129 - t118 * t169) * qJD(1) + t93) * qJD(4)) * t130) * MDP(20) + (-(t108 * t157 + t111 * t127) * t119 + ((t114 * t118 - t182) * qJD(4) + (t142 * t127 - t175) * qJD(5) + t139) * t128 + (-t96 * t158 - t118 * t100 + t180 + (-t118 * t167 - (t118 * t128 * t129 + t108 * t127) * qJD(1) - t94) * qJD(4)) * t130) * MDP(21); (t136 + t177) * MDP(20) + (t141 + t179) * MDP(21) + (t128 * MDP(14) - t154) * t131 + ((-t151 + t173) * MDP(20) - t138 * MDP(21) * t129) * qJD(4); -t132 * MDP(6) + (t134 - t171) * qJD(1) + (t163 * MDP(14) + (-t101 - t149) * MDP(20) + (-t119 * t155 - t100) * MDP(21)) * t130 + (t163 * MDP(13) + qJD(5) * t134 + ((t112 - t150) * MDP(20) + (t114 - t146) * MDP(21)) * qJD(4)) * t128; t116 * t162 * MDP(14) + (t114 * t167 + t178) * MDP(15) + ((t100 - t174) * t129 + (-t101 - t172) * t127) * MDP(16) + t148 * MDP(17) - t119 * t158 * MDP(18) + (-pkin(4) * t101 - t180 - (t115 * t129 + t127 * t137) * t119 - t103 * t112 + (-pkin(7) * t167 + t183) * qJD(5)) * MDP(20) + (-pkin(4) * t100 + t181 + (t115 * t127 - t129 * t137) * t119 - t103 * t114 + (pkin(7) * t170 + t182) * qJD(5)) * MDP(21) + (t128 * t186 - t185) * t132 + (-t116 * t154 + (t128 * t167 + (-t114 + t161) * t130) * MDP(17) + (-t119 * t169 + (t112 + t155) * t130) * MDP(18) - t119 * t153 + (t135 * t127 - t93 * t130) * MDP(20) + (t135 * t129 + t94 * t130) * MDP(21)) * qJD(1); t114 * t112 * MDP(15) + (-t112 ^ 2 + t114 ^ 2) * MDP(16) + (-t129 * t145 + t121 + t174) * MDP(17) + (t117 + t172) * MDP(18) + qJD(1) * t144 + (-t114 * t96 + t119 * t94 + t143) * MDP(20) + (t112 * t96 + t119 * t93 + t139) * MDP(21) + (-MDP(17) * t150 - t114 * MDP(18) - t94 * MDP(20) - t93 * MDP(21)) * qJD(5);];
tauc = t1;
