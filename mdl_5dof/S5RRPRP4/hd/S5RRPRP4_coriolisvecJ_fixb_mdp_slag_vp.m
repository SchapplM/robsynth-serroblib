% Calculate minimal parameter regressor of Coriolis joint torque vector for
% S5RRPRP4
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% qJD [5x1]
%   Generalized joint velocities
% pkin [7x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d1,d2,d4]';
% MDP [20x1]
%   Minimal dynamic parameter vector (fixed base model)
%   see S5RRPRP4_convert_par2_MPV_fixb.m
% 
% Output:
% tauc [5x1]
%   joint torques required to compensate Coriolis and centrifugal load

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-31 19:53
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauc = S5RRPRP4_coriolisvecJ_fixb_mdp_slag_vp(qJ, qJD, pkin, MDP)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(7,1),zeros(20,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RRPRP4_coriolisvecJ_fixb_mdp_slag_vp: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RRPRP4_coriolisvecJ_fixb_mdp_slag_vp: qJD has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [7 1]), ...
  'S5RRPRP4_coriolisvecJ_fixb_mdp_slag_vp: pkin has to be [7x1] (double)');
assert(isreal(MDP) && all(size(MDP) == [20 1]), ...
  'S5RRPRP4_coriolisvecJ_fixb_mdp_slag_vp: MDP has to be [20x1] (double)'); 

%% Symbolic Calculation
% From coriolisvec_joint_fixb_mdp_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 19:52:57
% EndTime: 2019-12-31 19:53:00
% DurationCPUTime: 0.72s
% Computational Cost: add. (683->158), mult. (1021->196), div. (0->0), fcn. (401->4), ass. (0->77)
t139 = sin(qJ(4));
t137 = t139 ^ 2;
t141 = cos(qJ(4));
t138 = t141 ^ 2;
t179 = t139 * t141;
t195 = MDP(10) * t179 - (t137 - t138) * MDP(11);
t136 = qJD(1) + qJD(2);
t143 = -pkin(2) - pkin(7);
t142 = cos(qJ(2));
t189 = pkin(1) * qJD(1);
t163 = t142 * t189;
t151 = qJD(3) - t163;
t110 = t143 * t136 + t151;
t181 = t110 * t141;
t156 = qJD(5) - t181;
t187 = qJD(4) * pkin(4);
t104 = t156 - t187;
t170 = qJD(4) * qJ(5);
t182 = t110 * t139;
t105 = t170 + t182;
t140 = sin(qJ(2));
t188 = pkin(1) * qJD(2);
t160 = qJD(1) * t188;
t153 = t140 * t160;
t125 = t139 * t153;
t98 = t125 + (qJD(5) + t181) * qJD(4);
t126 = t141 * t153;
t168 = qJD(4) * t139;
t100 = t110 * t168 - t126;
t99 = t100 * t141;
t194 = (t104 * t139 + t105 * t141) * qJD(4) + t139 * t98 - t99;
t124 = t139 * pkin(4) - qJ(5) * t141 + qJ(3);
t193 = t124 * t136;
t192 = -MDP(5) + MDP(7);
t135 = t136 ^ 2;
t190 = pkin(1) * t140;
t164 = t140 * t189;
t103 = t164 + t193;
t167 = qJD(4) * t141;
t129 = t142 * t160;
t149 = pkin(4) * t141 + qJ(5) * t139;
t154 = -qJD(5) * t141 + qJD(3);
t97 = t129 + (t149 * qJD(4) + t154) * t136;
t185 = t103 * t167 + t97 * t139;
t184 = qJ(3) * t136;
t183 = t103 * t136;
t159 = -pkin(1) * t142 - pkin(2);
t130 = -pkin(7) + t159;
t144 = qJD(4) ^ 2;
t180 = t130 * t144;
t178 = t139 * t144;
t177 = -0.2e1 * t195 * qJD(4) * t136;
t171 = qJD(3) * t136;
t118 = t129 + t171;
t123 = t164 + t184;
t176 = t118 * t139 + t123 * t167;
t173 = t137 + t138;
t172 = qJD(2) * t140;
t166 = MDP(15) + MDP(17);
t165 = -MDP(16) + MDP(19);
t162 = pkin(1) * t172;
t161 = t142 * t188;
t155 = t143 * MDP(20) - MDP(18);
t109 = pkin(4) * t167 + qJ(5) * t168 + t154;
t152 = -t109 + t163;
t150 = -t130 * t178 + t162 * t167;
t147 = -t164 + t184;
t146 = -t164 + t193;
t131 = qJ(3) + t190;
t128 = qJD(3) + t161;
t120 = -pkin(2) * t136 + t151;
t117 = t124 + t190;
t115 = t149 * t136;
t114 = t118 * t141;
t106 = t109 + t161;
t101 = t103 * t168;
t1 = [(-t136 * t161 - t129) * MDP(6) + (t129 + (qJD(3) + t128) * t136) * MDP(8) + (t118 * t131 + t123 * t128 + (t159 * qJD(1) + t120) * t162) * MDP(9) - MDP(12) * t178 - t144 * t141 * MDP(13) + ((t128 * t139 + t131 * t167) * t136 + t150 + t176) * MDP(15) + (t114 + (t128 * t136 - t180) * t141 + (-t131 * t136 - t123 - t162) * t168) * MDP(16) + ((t106 * t139 + t117 * t167) * t136 + t150 + t185) * MDP(17) + (-t173 * t136 * t162 - t194) * MDP(18) + (t101 + (t117 * t136 + t162) * t168 + (-t106 * t136 + t180 - t97) * t141) * MDP(19) + (t103 * t106 + t117 * t97 + (-t104 * t141 + t105 * t139) * t162 + t194 * t130) * MDP(20) + t177 + t192 * (qJD(1) + t136) * t162; -t129 * MDP(6) + (t129 + 0.2e1 * t171) * MDP(8) + (qJ(3) * t118 + qJD(3) * t123) * MDP(9) + t176 * MDP(15) + t114 * MDP(16) + t185 * MDP(17) + t99 * MDP(18) + t101 * MDP(19) + (t103 * t109 + t124 * t97) * MDP(20) + t192 * (qJD(2) - t136) * t164 + ((-pkin(2) * t172 - t120 * t140 - t123 * t142) * MDP(9) - t103 * t142 * MDP(20) + ((MDP(6) - MDP(8)) * t142 + t173 * MDP(18) * t140) * t136) * t189 + (-t97 * MDP(19) + (-t100 * t143 + t104 * t164) * MDP(20) + (t165 * t143 - MDP(13)) * t144 + (t151 * MDP(16) + t152 * MDP(19)) * t136 + (t147 * MDP(15) + t146 * MDP(17) + t155 * t105) * qJD(4)) * t141 + (-t98 * MDP(18) + (-t105 * t164 + t143 * t98) * MDP(20) + (-t166 * t143 - MDP(12)) * t144 + (t151 * MDP(15) - t152 * MDP(17)) * t136 + ((-t123 - t147) * MDP(16) + t146 * MDP(19) + t155 * t104) * qJD(4)) * t139 + t177; -t135 * MDP(8) + (-t123 * t136 + t153) * MDP(9) + (-t183 + t194) * MDP(20) + (-t166 * t139 + t165 * t141) * (t135 + t144); -t125 * MDP(16) + (0.2e1 * qJD(4) * qJD(5) + t125) * MDP(19) + (-pkin(4) * t100 + qJ(5) * t98 - t103 * t115 - t104 * t182 + t156 * t105) * MDP(20) + t166 * t126 + t195 * t135 + ((-t123 * MDP(15) - t103 * MDP(17) + (t105 - t170) * MDP(18) + t115 * MDP(19)) * t141 + (t123 * MDP(16) - t115 * MDP(17) + (-qJD(5) + t104 + t187) * MDP(18) - t103 * MDP(19)) * t139) * t136; t135 * MDP(17) * t179 + (-t135 * t138 - t144) * MDP(19) + (-t105 * qJD(4) + t141 * t183 + t100) * MDP(20);];
tauc = t1;
