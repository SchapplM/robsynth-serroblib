% Calculate minimal parameter regressor of Coriolis joint torque vector for
% S5PPRRP3
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% qJD [5x1]
%   Generalized joint velocities
% pkin [8x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d3,d4,theta1,theta2]';
% MDP [16x1]
%   Minimal dynamic parameter vector (fixed base model)
%   see S5PPRRP3_convert_par2_MPV_fixb.m
% 
% Output:
% tauc [5x1]
%   joint torques required to compensate Coriolis and centrifugal load

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-05 15:11
% Revision: 77da58f92bca3eff71542919beafa37024070d86 (2019-12-05)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauc = S5PPRRP3_coriolisvecJ_fixb_mdp_slag_vp(qJ, qJD, pkin, MDP)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(8,1),zeros(16,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5PPRRP3_coriolisvecJ_fixb_mdp_slag_vp: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5PPRRP3_coriolisvecJ_fixb_mdp_slag_vp: qJD has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S5PPRRP3_coriolisvecJ_fixb_mdp_slag_vp: pkin has to be [8x1] (double)');
assert(isreal(MDP) && all(size(MDP) == [16 1]), ...
  'S5PPRRP3_coriolisvecJ_fixb_mdp_slag_vp: MDP has to be [16x1] (double)'); 

%% Symbolic Calculation
% From coriolisvec_joint_fixb_mdp_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-05 15:11:18
% EndTime: 2019-12-05 15:11:20
% DurationCPUTime: 0.69s
% Computational Cost: add. (436->119), mult. (1098->172), div. (0->0), fcn. (705->6), ass. (0->67)
t139 = sin(qJ(3));
t141 = cos(qJ(3));
t136 = sin(pkin(8));
t175 = qJD(1) * t136;
t159 = t141 * t175;
t125 = qJD(2) * t139 + t159;
t119 = qJD(3) * pkin(6) + t125;
t138 = sin(qJ(4));
t140 = cos(qJ(4));
t137 = cos(pkin(8));
t174 = qJD(1) * t137;
t146 = t119 * t138 + t140 * t174;
t188 = qJD(5) + t146;
t103 = -qJD(4) * pkin(4) + t188;
t187 = -0.2e1 * qJD(4);
t142 = qJD(4) ^ 2;
t143 = qJD(3) ^ 2;
t186 = (t142 + t143) * t139;
t134 = t138 ^ 2;
t135 = t140 ^ 2;
t185 = (t134 - t135) * MDP(7);
t158 = t138 * t174;
t108 = t119 * t140 - t158;
t104 = qJD(4) * qJ(5) + t108;
t145 = -qJD(2) * t141 + t139 * t175;
t163 = MDP(11) + MDP(13);
t184 = -MDP(5) + (t134 + t135) * MDP(14);
t183 = pkin(6) * t142;
t182 = qJD(3) * pkin(3);
t181 = t136 * t141;
t180 = t138 * t140;
t168 = qJD(3) * t139;
t120 = qJD(2) * t168 + qJD(3) * t159;
t128 = -pkin(4) * t140 - qJ(5) * t138 - pkin(3);
t169 = qJD(3) * t128;
t109 = t145 + t169;
t172 = qJD(3) * t109;
t148 = pkin(4) * t138 - qJ(5) * t140;
t117 = qJD(4) * t148 - qJD(5) * t138;
t171 = qJD(3) * t117;
t170 = qJD(3) * t125;
t167 = qJD(3) * t140;
t165 = qJD(4) * t138;
t164 = qJD(4) * t140;
t162 = MDP(12) - MDP(15);
t161 = t138 * t181;
t160 = t140 * t181;
t121 = t145 * qJD(3);
t101 = -qJD(4) * t158 + t119 * t164 - t138 * t121;
t157 = t136 * t168;
t156 = -t120 - t183;
t118 = t145 - t182;
t155 = t118 - t182;
t154 = pkin(6) * MDP(16) + MDP(14);
t152 = t109 + t169;
t150 = qJD(3) * t141 * t187;
t149 = t140 * t157;
t102 = t171 + t120;
t147 = -t102 - t171 - t183;
t122 = t137 * t140 + t161;
t126 = t148 * qJD(3);
t123 = -t137 * t138 + t160;
t114 = t140 * t121;
t106 = -qJD(4) * t122 - t149;
t105 = -qJD(4) * t160 + (qJD(4) * t137 + t157) * t138;
t100 = -t114 + (qJD(5) - t146) * qJD(4);
t1 = [(t100 * t123 + t101 * t122 - t103 * t105 + t104 * t106) * MDP(16) - t162 * (-t143 * t161 + (t106 - t149) * qJD(4)) + (t102 * t139 * MDP(16) + (MDP(5) * t139 + (-t140 * t163 - MDP(4)) * t141) * t143) * t136 + ((-t105 * t138 + t106 * t140 + t122 * t164 - t123 * t165) * MDP(14) + t109 * MDP(16) * t181) * qJD(3) + t163 * (qJD(4) * t105 + t157 * t165); t163 * (t138 * t150 - t140 * t186) + t162 * (t138 * t186 + t140 * t150) + ((-t102 + (t103 * t138 + t104 * t140) * qJD(3)) * MDP(16) + t184 * t143) * t141 + (-t143 * MDP(4) + (t100 * t140 + t101 * t138 + t103 * t164 - t104 * t165 + t172) * MDP(16)) * t139; -t120 * MDP(4) + (t102 * t128 + (t117 - t125) * t109) * MDP(16) + (t125 * MDP(4) + t185 * t187 + (MDP(5) + t184) * t145) * qJD(3) + (t142 * MDP(8) + t156 * MDP(11) + t147 * MDP(13) + t100 * MDP(14) + (pkin(6) * t100 + t104 * t145) * MDP(16) + ((-t145 + t155) * MDP(12) + (t145 - t152) * MDP(15) + t154 * t103) * qJD(4)) * t140 + (-t142 * MDP(9) + (-t156 - t170) * MDP(12) + t101 * MDP(14) + (t147 + t170) * MDP(15) + (pkin(6) * t101 + t103 * t145) * MDP(16) + (MDP(11) * t155 + MDP(13) * t152 + 0.2e1 * MDP(6) * t167 - t104 * t154) * qJD(4)) * t138 + t163 * (t125 * t167 - t145 * t165); (qJ(5) * t100 - t103 * t108 + t188 * t104 - t109 * t126) * MDP(16) + t162 * t114 + (-MDP(6) * t180 + t185) * t143 + (0.2e1 * qJD(5) * MDP(15) + t163 * t108) * qJD(4) + ((-t118 * MDP(11) - t109 * MDP(13) + t126 * MDP(15)) * t138 + (-t118 * MDP(12) + t126 * MDP(13) + t109 * MDP(15)) * t140) * qJD(3) + (-MDP(16) * pkin(4) - t163) * t101; -t143 * MDP(13) * t180 + (-t134 * t143 - t142) * MDP(15) + (-qJD(4) * t104 + t138 * t172 + t101) * MDP(16);];
tauc = t1;
