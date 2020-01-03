% Calculate minimal parameter regressor of Coriolis joint torque vector for
% S5RPRPP5
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% qJD [5x1]
%   Generalized joint velocities
% pkin [6x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d1,d3]';
% MDP [21x1]
%   Minimal dynamic parameter vector (fixed base model)
%   see S5RPRPP5_convert_par2_MPV_fixb.m
% 
% Output:
% tauc [5x1]
%   joint torques required to compensate Coriolis and centrifugal load

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-31 18:16
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauc = S5RPRPP5_coriolisvecJ_fixb_mdp_slag_vp(qJ, qJD, pkin, MDP)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(6,1),zeros(21,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RPRPP5_coriolisvecJ_fixb_mdp_slag_vp: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RPRPP5_coriolisvecJ_fixb_mdp_slag_vp: qJD has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [6 1]), ...
  'S5RPRPP5_coriolisvecJ_fixb_mdp_slag_vp: pkin has to be [6x1] (double)');
assert(isreal(MDP) && all(size(MDP) == [21 1]), ...
  'S5RPRPP5_coriolisvecJ_fixb_mdp_slag_vp: MDP has to be [21x1] (double)'); 

%% Symbolic Calculation
% From coriolisvec_joint_fixb_mdp_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 18:16:44
% EndTime: 2019-12-31 18:16:46
% DurationCPUTime: 0.75s
% Computational Cost: add. (493->151), mult. (930->190), div. (0->0), fcn. (344->2), ass. (0->73)
t144 = -pkin(3) - pkin(4);
t145 = -pkin(1) - pkin(6);
t130 = t145 * qJD(1) + qJD(2);
t143 = cos(qJ(3));
t172 = qJ(5) * qJD(1);
t114 = (t130 + t172) * t143;
t167 = qJD(4) - t114;
t105 = qJD(3) * t144 + t167;
t142 = sin(qJ(3));
t140 = t142 ^ 2;
t141 = t143 ^ 2;
t185 = (t140 - t141) * MDP(8);
t183 = -MDP(17) * t145 + MDP(15);
t182 = qJD(3) * pkin(3);
t181 = qJ(4) * t142;
t124 = t142 * t130;
t139 = qJD(3) * qJ(4);
t119 = t124 + t139;
t180 = t119 * t143;
t179 = t130 * t143;
t178 = t143 * MDP(7);
t177 = qJ(5) + t145;
t171 = qJD(3) * t142;
t121 = t130 * t171;
t160 = qJD(3) * t172;
t176 = t142 * t160 + t121;
t170 = qJD(3) * t143;
t169 = qJD(4) * t143;
t157 = qJ(4) * t143 - qJ(2);
t127 = pkin(3) * t142 - t157;
t118 = qJD(1) * t127;
t168 = t118 * MDP(14);
t120 = t142 * t144 + t157;
t108 = qJD(1) * t120 + qJD(5);
t166 = qJD(5) + t108;
t165 = qJD(1) * qJD(5);
t164 = MDP(14) + MDP(18);
t163 = MDP(16) + MDP(19);
t138 = qJD(3) * qJD(4);
t122 = t130 * t170;
t161 = t142 * t165 + t143 * t160 + t122;
t101 = t138 + t161;
t136 = t142 * t172;
t109 = t136 + t119;
t162 = t101 * t142 + t105 * t171 + t109 * t170;
t159 = -qJ(2) * MDP(6) - MDP(5);
t158 = -qJD(4) + t182;
t117 = -t158 - t179;
t156 = -t117 + t179;
t155 = qJD(3) * t177;
t152 = pkin(3) * t143 + t181;
t151 = t118 * MDP(16) - t108 * MDP(19);
t150 = t144 * t143 - t181;
t149 = t152 * qJD(3) + qJD(2);
t148 = t150 * qJD(3) - qJD(2);
t147 = qJD(1) ^ 2;
t146 = qJD(3) ^ 2;
t137 = 0.2e1 * t138;
t135 = qJD(1) * t169;
t126 = t177 * t143;
t125 = t177 * t142;
t123 = t152 * qJD(1);
t116 = t122 + t138;
t115 = t150 * qJD(1);
t113 = t124 + t136;
t112 = t149 - t169;
t111 = qJD(5) * t142 + t143 * t155;
t110 = -qJD(5) * t143 + t142 * t155;
t107 = t148 + t169;
t106 = t149 * qJD(1) - t135;
t103 = -t143 * t165 + t176;
t100 = t148 * qJD(1) + t135;
t1 = [(t106 * t127 + t112 * t118) * MDP(17) + t162 * MDP(20) + (t100 * t120 + t101 * t125 - t103 * t126 + t105 * t110 + t107 * t108 + t109 * t111) * MDP(21) + (-t110 * MDP(18) + t111 * MDP(19)) * qJD(3) + (-t106 * MDP(16) + t100 * MDP(19) - t103 * MDP(20) + (-MDP(10) + (-MDP(13) + MDP(16)) * t145) * t146 + (-t108 * MDP(18) - t119 * t183 + t168) * qJD(3)) * t143 + (t106 * MDP(14) - t100 * MDP(18) + (-MDP(9) + (-MDP(12) - MDP(14)) * t145) * t146 - t183 * t116 + (t183 * t156 + t151) * qJD(3)) * t142 + ((-t110 * t143 + t111 * t142) * MDP(20) + (t142 * MDP(14) - t143 * MDP(16)) * t112 + (-t142 * MDP(18) + t143 * MDP(19)) * t107 + 0.2e1 * (t142 * MDP(12) + t143 * MDP(13) - t159) * qJD(2) + (0.2e1 * t185 + (0.2e1 * qJ(2) * MDP(12) + t127 * MDP(14) - t120 * MDP(18) + t125 * MDP(20)) * t143 + (-0.2e1 * qJ(2) * MDP(13) + t127 * MDP(16) - t120 * MDP(19) - t126 * MDP(20) - 0.2e1 * t178) * t142) * qJD(3)) * qJD(1); (qJD(1) * t108 - t103 * t143 + t162) * MDP(21) + t159 * t147 + ((-MDP(13) + t163) * t143 - (MDP(12) + t164) * t142) * (t146 + t147) + (-qJD(1) * t118 + t116 * t142 + (-t156 * t142 + t180) * qJD(3)) * MDP(17); t137 * MDP(16) + (qJD(3) * t113 - t176) * MDP(18) + (-qJD(3) * t114 + t137 + t161) * MDP(19) + (qJ(4) * t101 + t103 * t144 - t105 * t113 - t108 * t115 + t167 * t109) * MDP(21) + (t142 * t178 - t185 + (-t143 * MDP(12) + t142 * MDP(13)) * qJ(2)) * t147 + ((-t168 + (t119 - t139) * MDP(15) + t123 * MDP(16) + t166 * MDP(18) - t115 * MDP(19) + (-t109 + t113 + t139) * MDP(20)) * t143 + (-t123 * MDP(14) + (t117 + t158) * MDP(15) + t115 * MDP(18) - t151) * t142) * qJD(1) + (qJ(4) * t116 + qJD(4) * t119 - t118 * t123 + (-t180 + (-t117 - t182) * t142) * t130) * MDP(17); (-qJD(3) * t119 + t121) * MDP(17) + (-qJD(3) * t109 + t176) * MDP(21) + t163 * (-t141 * t147 - t146) + (t164 * t147 * t142 + (t118 * MDP(17) - t166 * MDP(21)) * qJD(1)) * t143; t135 * MDP(21) + (-t140 - t141) * MDP(20) * t147 + ((t105 * t143 - t109 * t142 - qJD(2)) * MDP(21) + ((-MDP(21) * qJ(4) - 0.2e1 * MDP(19)) * t142 + (t144 * MDP(21) - 0.2e1 * MDP(18)) * t143) * qJD(3)) * qJD(1);];
tauc = t1;
