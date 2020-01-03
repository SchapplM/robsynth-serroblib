% Calculate minimal parameter regressor of Coriolis joint torque vector for
% S5PRRPP4
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% qJD [5x1]
%   Generalized joint velocities
% pkin [7x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d2,d3,theta1]';
% MDP [19x1]
%   Minimal dynamic parameter vector (fixed base model)
%   see S5PRRPP4_convert_par2_MPV_fixb.m
% 
% Output:
% tauc [5x1]
%   joint torques required to compensate Coriolis and centrifugal load

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-31 17:41
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauc = S5PRRPP4_coriolisvecJ_fixb_mdp_slag_vp(qJ, qJD, pkin, MDP)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(7,1),zeros(19,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5PRRPP4_coriolisvecJ_fixb_mdp_slag_vp: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5PRRPP4_coriolisvecJ_fixb_mdp_slag_vp: qJD has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [7 1]), ...
  'S5PRRPP4_coriolisvecJ_fixb_mdp_slag_vp: pkin has to be [7x1] (double)');
assert(isreal(MDP) && all(size(MDP) == [19 1]), ...
  'S5PRRPP4_coriolisvecJ_fixb_mdp_slag_vp: MDP has to be [19x1] (double)'); 

%% Symbolic Calculation
% From coriolisvec_joint_fixb_mdp_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 17:41:14
% EndTime: 2019-12-31 17:41:16
% DurationCPUTime: 0.69s
% Computational Cost: add. (408->148), mult. (945->196), div. (0->0), fcn. (407->2), ass. (0->71)
t135 = cos(qJ(3));
t134 = sin(qJ(3));
t161 = qJD(2) * t134;
t113 = -pkin(6) * t161 + t135 * qJD(1);
t175 = qJD(4) - t113;
t102 = qJ(5) * t161 + t113;
t155 = qJD(4) - t102;
t171 = pkin(3) + pkin(4);
t95 = -t171 * qJD(3) + t155;
t132 = t134 ^ 2;
t133 = t135 ^ 2;
t174 = (t132 - t133) * MDP(6);
t160 = qJD(2) * t135;
t114 = pkin(6) * t160 + t134 * qJD(1);
t104 = -qJ(5) * t160 + t114;
t131 = qJD(3) * qJ(4);
t99 = t104 + t131;
t110 = t131 + t114;
t107 = -qJD(3) * pkin(3) + t175;
t173 = MDP(10) + MDP(12);
t172 = -2 * pkin(2);
t170 = pkin(6) - qJ(5);
t169 = pkin(6) * qJD(3);
t168 = qJ(4) * t135;
t167 = t135 * MDP(5);
t148 = qJ(4) * t134 + pkin(2);
t109 = t171 * t135 + t148;
t97 = qJD(2) * t109 + qJD(5);
t166 = qJD(5) + t97;
t153 = qJD(2) * qJD(3);
t149 = t135 * t153;
t158 = qJD(4) * t134;
t165 = qJ(4) * t149 + qJD(2) * t158;
t154 = qJD(1) * qJD(3);
t111 = pkin(6) * t149 + t134 * t154;
t124 = t135 * t154;
t130 = qJD(3) * qJD(4);
t164 = t124 + t130;
t162 = qJ(5) * qJD(3);
t159 = qJD(3) * t134;
t157 = qJD(5) * t134;
t156 = qJD(5) * t135;
t152 = -MDP(12) - MDP(16);
t151 = MDP(14) + MDP(17);
t117 = t170 * t135;
t150 = t134 * t153;
t93 = -t171 * t150 + t165;
t140 = -t171 * t134 + t168;
t96 = t140 * qJD(3) + t158;
t147 = qJD(2) * t96 + t93;
t142 = pkin(3) * t134 - t168;
t106 = t142 * qJD(3) - t158;
t98 = pkin(3) * t150 - t165;
t145 = -qJD(2) * t106 - t98;
t144 = qJD(3) * t113 - t124;
t115 = -pkin(3) * t135 - t148;
t139 = (t134 * t95 + t135 * t99) * MDP(19);
t138 = qJD(2) ^ 2;
t137 = qJD(3) ^ 2;
t129 = 0.2e1 * t130;
t119 = qJ(5) * t150;
t116 = t170 * t134;
t112 = t142 * qJD(2);
t108 = qJD(2) * t115;
t105 = qJD(3) * t117 - t157;
t103 = -t170 * t159 - t156;
t101 = -pkin(6) * t150 + t164;
t100 = t140 * qJD(2);
t94 = (-t135 * t162 - t157) * qJD(2) + t111;
t92 = t119 + (-pkin(6) * t159 - t156) * qJD(2) + t164;
t1 = [(t101 * t134 - t111 * t135) * MDP(15) + (t134 * t92 - t135 * t94) * MDP(19) + ((t107 * t134 + t110 * t135) * MDP(15) + t139) * qJD(3) + ((-MDP(11) + t151) * t135 + (-MDP(10) + t152) * t134) * t137; (t106 * t108 + t115 * t98) * MDP(15) + (t103 * t99 + t105 * t95 + t109 * t93 + t116 * t94 + t117 * t92 + t96 * t97) * MDP(19) + (t137 * MDP(7) + t145 * MDP(12) + t101 * MDP(13) + t147 * MDP(16) + (-qJD(2) * t103 - t92) * MDP(18) + (t101 * MDP(15) - t173 * t137) * pkin(6)) * t135 + (-t137 * MDP(8) + t111 * MDP(13) + t145 * MDP(14) + t147 * MDP(17) + (-qJD(2) * t105 - t94) * MDP(18) + (t111 * MDP(15) + (MDP(11) - MDP(14)) * t137) * pkin(6)) * t134 + ((-t134 * t97 - t105) * MDP(16) + (t135 * t97 + t103) * MDP(17) + (t134 * t99 - t135 * t95) * MDP(18) + (t134 * MDP(12) - t135 * MDP(14)) * t108 + (-0.2e1 * t174 + (MDP(11) * t172 - t115 * MDP(14) + t109 * MDP(17) - t116 * MDP(18)) * t135 + (MDP(10) * t172 + t115 * MDP(12) - t109 * MDP(16) + t117 * MDP(18) + 0.2e1 * t167) * t134) * qJD(2) + (MDP(15) * pkin(6) + MDP(13)) * (t107 * t135 - t110 * t134)) * qJD(3); t144 * MDP(11) + (t129 - t144) * MDP(14) + (-pkin(3) * t111 + qJ(4) * t101 - t107 * t114 - t108 * t112 + t175 * t110) * MDP(15) + (qJD(3) * t104 - t111) * MDP(16) + (-qJD(3) * t102 + t119 + t124 + t129) * MDP(17) + (qJ(4) * t92 - t100 * t97 - t104 * t95 + t155 * t99 - t171 * t94) * MDP(19) + (-t134 * t167 + t174 + (t134 * MDP(10) + t135 * MDP(11)) * pkin(2)) * t138 + ((t112 * MDP(12) + t108 * MDP(14) + (-t100 + t162) * MDP(16) - t166 * MDP(17)) * t135 + (MDP(11) * t169 - t108 * MDP(12) + (t112 - t169) * MDP(14) + t166 * MDP(16) + (-t100 - t169) * MDP(17)) * t134) * qJD(2) + t173 * (qJD(3) * t114 - t111); (-qJD(3) * t110 + t111) * MDP(15) + (-qJ(5) * t149 - qJD(3) * t99 + t111) * MDP(19) + t151 * (-t132 * t138 - t137) + (t152 * t138 * t135 + (t108 * MDP(15) - t166 * MDP(19)) * qJD(2)) * t134; t165 * MDP(19) + (-t132 - t133) * MDP(18) * t138 + (t139 + (0.2e1 * t135 * MDP(17) + (-t171 * MDP(19) - 0.2e1 * MDP(16)) * t134) * qJD(3)) * qJD(2);];
tauc = t1;
