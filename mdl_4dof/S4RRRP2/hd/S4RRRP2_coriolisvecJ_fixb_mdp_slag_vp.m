% Calculate Coriolis joint torque vector for
% S4RRRP2
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
% MDP [17x1]
%   Minimal dynamic parameter vector (fixed base model)
%   see S4RRRP2_convert_par2_MPV_fixb.m
% 
% Output:
% tauc [4x1]
%   joint torques required to compensate Coriolis and centrifugal load

% Quelle: HybrDyn-Toolbox
% Datum: 2021-01-15 11:04
% Revision: d12c3222fdeb2c5f3b3c8fa5751e113be2fc3aae (2021-01-15)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauc = S4RRRP2_coriolisvecJ_fixb_mdp_slag_vp(qJ, qJD, pkin, MDP)
%% Coder Information
%#codegen
%$cgargs {zeros(4,1),zeros(4,1),zeros(6,1),zeros(17,1)}
assert(isreal(qJ) && all(size(qJ) == [4 1]), ...
  'S4RRRP2_coriolisvecJ_fixb_mdp_slag_vp: qJ has to be [4x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [4 1]), ...
  'S4RRRP2_coriolisvecJ_fixb_mdp_slag_vp: qJD has to be [4x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [6 1]), ...
  'S4RRRP2_coriolisvecJ_fixb_mdp_slag_vp: pkin has to be [6x1] (double)');
assert(isreal(MDP) && all(size(MDP) == [17 1]), ...
  'S4RRRP2_coriolisvecJ_fixb_mdp_slag_vp: MDP has to be [17x1] (double)'); 

%% Symbolic Calculation
% From coriolisvec_joint_fixb_mdp_matlab.m
% OptimizationMode: 2
% StartTime: 2021-01-15 11:04:30
% EndTime: 2021-01-15 11:04:32
% DurationCPUTime: 0.56s
% Computational Cost: add. (498->140), mult. (896->200), div. (0->0), fcn. (399->4), ass. (0->81)
t135 = sin(qJ(3));
t133 = t135 ^ 2;
t137 = cos(qJ(3));
t134 = t137 ^ 2;
t184 = (t133 - t134) * MDP(8);
t138 = cos(qJ(2));
t183 = pkin(1) * t138;
t182 = pkin(3) * t137;
t181 = -qJ(4) - pkin(6);
t180 = MDP(17) * pkin(3);
t179 = pkin(1) * qJD(1);
t178 = pkin(1) * qJD(2);
t177 = qJD(3) * pkin(3);
t132 = qJD(1) + qJD(2);
t136 = sin(qJ(2));
t158 = t136 * t179;
t115 = pkin(6) * t132 + t158;
t149 = qJ(4) * t132 + t115;
t100 = t149 * t135;
t97 = -t100 + t177;
t176 = t100 + t97;
t128 = -pkin(2) - t182;
t157 = t138 * t179;
t103 = t128 * t132 + qJD(4) - t157;
t153 = qJD(1) * t178;
t124 = t136 * t153;
t165 = qJD(3) * t135;
t152 = t132 * t165;
t107 = pkin(3) * t152 + t124;
t164 = qJD(3) * t137;
t175 = t103 * t164 + t107 * t135;
t131 = t132 ^ 2;
t174 = t131 * t135;
t173 = t132 * t137;
t139 = qJD(3) ^ 2;
t172 = t137 * t139;
t171 = t138 * MDP(6);
t126 = pkin(1) * t136 + pkin(6);
t170 = -qJ(4) - t126;
t116 = -pkin(2) * t132 - t157;
t169 = t116 * t164 + t135 * t124;
t146 = qJD(3) * t157;
t168 = t135 * t146 + t158 * t173;
t166 = MDP(15) * t132;
t163 = t137 * MDP(12);
t162 = -qJD(4) - t103;
t161 = qJD(3) * MDP(14);
t160 = qJD(3) * MDP(15);
t159 = -MDP(13) - MDP(15);
t156 = t138 * t178;
t155 = t103 * t180;
t154 = pkin(3) * t133 * MDP(15);
t151 = qJD(3) * t181;
t150 = t103 * t165 - t107 * t137;
t148 = qJD(3) * t170;
t147 = (-t133 - t134) * MDP(16);
t145 = t138 * t153;
t101 = t149 * t137;
t144 = -t101 * t137 + t135 * t97;
t143 = qJD(3) * t149;
t142 = 0.2e1 * t137 * MDP(7) * t152 - 0.2e1 * t132 * qJD(3) * t184 - t124 * MDP(5) + MDP(9) * t172;
t141 = (-pkin(2) - t183) * t132 - t156;
t140 = qJD(4) * t132 + t145;
t130 = t137 * qJ(4);
t129 = t137 * qJD(4);
t122 = pkin(6) * t137 + t130;
t121 = t181 * t135;
t119 = t137 * t146;
t117 = t128 - t183;
t113 = pkin(3) * t165 + t136 * t178;
t111 = t126 * t137 + t130;
t110 = t170 * t135;
t108 = t116 * t165;
t106 = -qJD(4) * t135 + t137 * t151;
t105 = t135 * t151 + t129;
t96 = (-qJD(4) - t156) * t135 + t137 * t148;
t95 = t135 * t148 + t137 * t156 + t129;
t94 = -t140 * t135 - t137 * t143;
t93 = -t135 * t143 + t140 * t137;
t92 = t93 * t137;
t1 = [(-t126 * t172 + t108) * MDP(12) + t169 * MDP(13) + (-t113 * t173 + t150) * MDP(14) + t175 * MDP(15) + (t95 * t173 + t92) * MDP(16) + (t101 * t95 + t103 * t113 + t107 * t117 + t110 * t94 + t111 * t93 + t96 * t97) * MDP(17) + (t113 * t166 + (-t132 * t96 - t94) * MDP(16) + (t126 * MDP(13) - MDP(10)) * t139) * t135 + ((-qJD(1) - t132) * t171 + (-qJD(1) * t163 + (t135 * MDP(13) - MDP(5) - t163) * t132) * t136) * t178 + (t96 * MDP(14) - t95 * MDP(15) + (t141 * MDP(13) + t117 * t166 + (-t110 * t132 - t97) * MDP(16)) * t137 + (t141 * MDP(12) + t117 * t132 * MDP(14) + (-t111 * t132 - t101) * MDP(16)) * t135) * qJD(3) + t142; (-pkin(6) * t172 + t108 + t168) * MDP(12) + (t119 + t169) * MDP(13) + (qJD(3) * t106 + t150 + t168) * MDP(14) + (-qJD(3) * t105 + t119 + t175) * MDP(15) + (-t97 * t164 + t92) * MDP(16) + (t101 * t105 + t106 * t97 + t107 * t128 + t121 * t94 + t122 * t93) * MDP(17) + ((-qJD(3) * t101 - t94) * MDP(16) + qJD(3) * t155 + (pkin(6) * MDP(13) - MDP(10)) * t139) * t135 + ((-t103 * t136 + t144 * t138) * MDP(17) + (-t136 * t163 - t171) * qJD(2)) * t179 + ((t105 * t137 - t106 * t135) * MDP(16) + ((MDP(6) + t147) * t138 + (t159 * t135 + MDP(5)) * t136) * t179 + (t154 + (-pkin(2) * MDP(13) + t128 * MDP(15) - t121 * MDP(16)) * t137 + (-pkin(2) * MDP(12) + (t128 - t182) * MDP(14) - t122 * MDP(16)) * t135) * qJD(3)) * t132 + t142; -t100 * t160 + t94 * t180 + (-t154 + t184) * t131 + (t176 * MDP(17) + t161) * t101 + ((-MDP(12) - MDP(14)) * t145 + (-t116 * MDP(12) + t162 * MDP(14) - t155) * t132 + t149 * t160) * t135 + (-MDP(7) * t174 + (pkin(3) * t174 - qJD(3) * t115) * MDP(14) + t159 * t145 + (-t116 * MDP(13) - qJ(4) * t161 + t162 * MDP(15) + (t176 - t177) * MDP(16)) * t132) * t137; t107 * MDP(17) + t131 * t147 + (t144 * MDP(17) + 0.2e1 * (t135 * MDP(14) + t137 * MDP(15)) * qJD(3)) * t132;];
tauc = t1;
