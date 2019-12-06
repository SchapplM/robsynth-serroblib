% Calculate minimal parameter regressor of Coriolis joint torque vector for
% S5RRPRP1
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% qJD [5x1]
%   Generalized joint velocities
% pkin [8x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d1,d2,d4,theta3]';
% MDP [16x1]
%   Minimal dynamic parameter vector (fixed base model)
%   see S5RRPRP1_convert_par2_MPV_fixb.m
% 
% Output:
% tauc [5x1]
%   joint torques required to compensate Coriolis and centrifugal load

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-05 18:22
% Revision: 77da58f92bca3eff71542919beafa37024070d86 (2019-12-05)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauc = S5RRPRP1_coriolisvecJ_fixb_mdp_slag_vp(qJ, qJD, pkin, MDP)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(8,1),zeros(16,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RRPRP1_coriolisvecJ_fixb_mdp_slag_vp: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RRPRP1_coriolisvecJ_fixb_mdp_slag_vp: qJD has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S5RRPRP1_coriolisvecJ_fixb_mdp_slag_vp: pkin has to be [8x1] (double)');
assert(isreal(MDP) && all(size(MDP) == [16 1]), ...
  'S5RRPRP1_coriolisvecJ_fixb_mdp_slag_vp: MDP has to be [16x1] (double)'); 

%% Symbolic Calculation
% From coriolisvec_joint_fixb_mdp_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-05 18:22:25
% EndTime: 2019-12-05 18:22:28
% DurationCPUTime: 0.66s
% Computational Cost: add. (695->131), mult. (1363->191), div. (0->0), fcn. (735->6), ass. (0->81)
t152 = sin(qJ(4));
t148 = t152 ^ 2;
t154 = cos(qJ(4));
t149 = t154 ^ 2;
t197 = (t148 - t149) * MDP(9);
t147 = qJD(1) + qJD(2);
t155 = cos(qJ(2));
t192 = pkin(1) * qJD(1);
t172 = t155 * t192;
t132 = pkin(2) * t147 + t172;
t151 = cos(pkin(8));
t153 = sin(qJ(2));
t173 = t153 * t192;
t136 = t151 * t173;
t150 = sin(pkin(8));
t115 = t150 * t132 + t136;
t168 = t115 + (pkin(7) + qJ(5)) * t147;
t104 = t154 * qJD(3) - t168 * t152;
t196 = MDP(5) * t153 + MDP(6) * t155;
t105 = qJD(3) * t152 + t168 * t154;
t195 = qJD(4) * t105;
t194 = pkin(2) * t151;
t193 = pkin(4) * t154;
t191 = pkin(1) * qJD(2);
t190 = qJD(4) * pkin(4);
t187 = t150 * t153;
t186 = t151 * t153;
t156 = qJD(4) ^ 2;
t185 = t154 * t156;
t142 = pkin(1) * t155 + pkin(2);
t180 = pkin(1) * t186 + t150 * t142;
t123 = pkin(7) + t180;
t184 = -qJ(5) - t123;
t139 = pkin(2) * t150 + pkin(7);
t183 = -qJ(5) - t139;
t103 = t104 + t190;
t182 = t103 - t104;
t135 = t150 * t173;
t114 = t132 * t151 - t135;
t110 = -pkin(3) * t147 - t114;
t158 = t150 * t155 + t186;
t125 = t158 * t191;
t120 = qJD(1) * t125;
t174 = qJD(4) * t154;
t181 = t110 * t174 + t120 * t152;
t179 = -t148 - t149;
t170 = -pkin(3) - t193;
t106 = t170 * t147 + qJD(5) - t114;
t177 = MDP(16) * t106;
t175 = qJD(4) * t152;
t169 = t147 * t175;
t171 = 0.2e1 * t154 * MDP(8) * t169 - 0.2e1 * qJD(4) * t147 * t197 + MDP(10) * t185;
t127 = (t151 * t155 - t187) * t191;
t121 = qJD(1) * t127;
t167 = -t110 * t147 - t121;
t166 = -pkin(1) * t187 + t142 * t151;
t165 = qJD(4) * t184;
t164 = qJD(4) * t183;
t163 = qJD(5) * t147 + t121;
t122 = -pkin(3) - t166;
t124 = t150 * t172 + t136;
t159 = -t124 * t147 + t139 * t156;
t126 = t151 * t172 - t135;
t157 = qJD(4) * ((-pkin(3) - t194) * t147 + t126);
t109 = pkin(4) * t169 + t120;
t146 = t147 ^ 2;
t145 = t154 * qJ(5);
t143 = t154 * qJD(5);
t130 = t139 * t154 + t145;
t129 = t183 * t152;
t118 = -qJD(5) * t152 + t154 * t164;
t117 = t152 * t164 + t143;
t113 = t123 * t154 + t145;
t112 = t184 * t152;
t107 = t110 * t175;
t102 = (-qJD(5) - t127) * t152 + t154 * t165;
t101 = t127 * t154 + t152 * t165 + t143;
t100 = -t163 * t152 - t195;
t99 = qJD(4) * t104 + t163 * t154;
t98 = t99 * t154;
t1 = [(-t114 * t125 + t115 * t127 - t120 * t166 + t121 * t180) * MDP(7) + (-t120 * t154 - t123 * t185 + t107) * MDP(13) + (-t127 * t174 + t181) * MDP(14) + (-t103 * t174 + t98) * MDP(15) + (t99 * t113 + t105 * t101 + t100 * t112 + t103 * t102 + t109 * (t122 - t193)) * MDP(16) + (-t100 * MDP(15) + (MDP(14) * t123 - MDP(11)) * t156 + (-MDP(13) * t127 - MDP(15) * t105 + pkin(4) * t177) * qJD(4)) * t152 + (-qJD(1) * t196 + t158 * t177) * t191 + ((t122 * t175 - t125 * t154) * MDP(13) + (t122 * t174 + t125 * t152) * MDP(14) + (t101 * t154 - t102 * t152 - t112 * t174 - t113 * t175) * MDP(15) - t196 * t191) * t147 + t171; (t114 * t124 - t115 * t126 + (-t120 * t151 + t121 * t150) * pkin(2)) * MDP(7) - t156 * t152 * MDP(11) + (t107 + t152 * t157 + (-t120 - t159) * t154) * MDP(13) + (t159 * t152 + t154 * t157 + t181) * MDP(14) + (-t100 * t152 + t98 + (-t103 * t154 - t105 * t152) * qJD(4) + (t117 * t154 - t118 * t152 + t179 * t126 + (-t129 * t154 - t130 * t152) * qJD(4)) * t147) * MDP(15) + (t99 * t130 + t100 * t129 + t109 * (t170 - t194) + (pkin(4) * t175 - t124) * t106 + (-t126 * t154 + t117) * t105 + (t126 * t152 + t118) * t103) * MDP(16) + t171 + t196 * (-qJD(2) + t147) * t192; (-t156 * MDP(14) + (t100 + t195) * MDP(16)) * t154 + (-t156 * MDP(13) + (-qJD(4) * t103 + t99) * MDP(16)) * t152; t146 * t197 + t167 * t152 * MDP(13) + (t182 * t105 + (-t106 * t147 * t152 + t100) * pkin(4)) * MDP(16) + (-t152 * t146 * MDP(8) + t167 * MDP(14) + (t182 - t190) * t147 * MDP(15)) * t154; ((t103 * t152 - t105 * t154) * t147 + t109) * MDP(16) + t179 * MDP(15) * t146;];
tauc = t1;
