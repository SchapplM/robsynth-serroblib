% Calculate minimal parameter regressor of Coriolis joint torque vector for
% S5RPPRP6
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
% MDP [21x1]
%   Minimal dynamic parameter vector (fixed base model)
%   see S5RPPRP6_convert_par2_MPV_fixb.m
% 
% Output:
% tauc [5x1]
%   joint torques required to compensate Coriolis and centrifugal load

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-31 17:55
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauc = S5RPPRP6_coriolisvecJ_fixb_mdp_slag_vp(qJ, qJD, pkin, MDP)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(7,1),zeros(21,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RPPRP6_coriolisvecJ_fixb_mdp_slag_vp: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RPPRP6_coriolisvecJ_fixb_mdp_slag_vp: qJD has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [7 1]), ...
  'S5RPPRP6_coriolisvecJ_fixb_mdp_slag_vp: pkin has to be [7x1] (double)');
assert(isreal(MDP) && all(size(MDP) == [21 1]), ...
  'S5RPPRP6_coriolisvecJ_fixb_mdp_slag_vp: MDP has to be [21x1] (double)'); 

%% Symbolic Calculation
% From coriolisvec_joint_fixb_mdp_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 17:55:18
% EndTime: 2019-12-31 17:55:20
% DurationCPUTime: 0.78s
% Computational Cost: add. (739->151), mult. (1583->195), div. (0->0), fcn. (987->4), ass. (0->67)
t148 = sin(pkin(7));
t149 = cos(pkin(7));
t174 = t148 ^ 2 + t149 ^ 2;
t184 = t174 * qJD(3);
t151 = sin(qJ(4));
t173 = qJD(1) * t148;
t163 = t151 * t173;
t152 = cos(qJ(4));
t175 = t152 * t149;
t165 = qJD(1) * t175;
t122 = -t163 + t165;
t127 = t148 * t152 + t149 * t151;
t183 = qJD(1) * t127;
t156 = t127 * qJD(3);
t182 = qJD(1) * t156;
t167 = MDP(16) + MDP(18);
t181 = MDP(17) - MDP(20);
t180 = t122 ^ 2;
t150 = -pkin(1) - qJ(3);
t179 = -pkin(6) + t150;
t134 = t150 * qJD(1) + qJD(2);
t162 = -pkin(6) * qJD(1) + t134;
t117 = t162 * t148;
t178 = t117 * t151;
t177 = t148 * MDP(7);
t176 = t149 * MDP(8);
t139 = t148 * pkin(3) + qJ(2);
t170 = qJD(4) * t151;
t169 = qJD(4) * t152;
t142 = qJD(1) * qJ(2) + qJD(3);
t118 = t162 * t149;
t104 = t118 * t152 - t178;
t168 = qJD(5) - t104;
t164 = t149 * t169;
t130 = pkin(3) * t173 + t142;
t161 = t104 + t178;
t96 = t122 * qJD(3) + t117 * t169 + t118 * t170;
t115 = qJD(4) * t183;
t124 = -t148 * t169 - t149 * t170;
t126 = t148 * t151 - t175;
t159 = t115 * t126 + t122 * t124;
t105 = t117 * t152 + t118 * t151;
t128 = t179 * t148;
t129 = t179 * t149;
t111 = t128 * t152 + t129 * t151;
t110 = t128 * t151 - t129 * t152;
t133 = qJD(4) * t163;
t116 = qJD(1) * t164 - t133;
t146 = qJD(1) * qJD(2);
t158 = pkin(4) * t116 + qJ(5) * t115 + t146;
t103 = pkin(4) * t183 - qJ(5) * t122 + t130;
t157 = -t103 * t122 - t96;
t102 = qJD(4) * qJ(5) + t105;
t125 = -t148 * t170 + t164;
t114 = t118 * t169;
t95 = t114 + (qJD(5) - t178) * qJD(4) - t182;
t99 = -qJD(4) * pkin(4) + t168;
t154 = t102 * t125 - t124 * t99 + t126 * t96 + t127 * t95;
t153 = qJD(1) ^ 2;
t119 = t183 ^ 2;
t109 = pkin(4) * t122 + qJ(5) * t183;
t108 = pkin(4) * t127 + qJ(5) * t126 + t139;
t101 = -t126 * qJD(3) + t111 * qJD(4);
t100 = -t110 * qJD(4) - t156;
t98 = pkin(4) * t125 - qJ(5) * t124 + qJD(5) * t126 + qJD(2);
t97 = -qJD(5) * t122 + t158;
t1 = [(qJD(2) * t142 - t184 * t134) * MDP(10) + t159 * MDP(11) + (t115 * t127 + t116 * t126 - t122 * t125 - t124 * t183) * MDP(12) + (qJD(2) * t183 + t116 * t139 + t125 * t130) * MDP(16) + (qJD(2) * t122 - t115 * t139 + t124 * t130) * MDP(17) + (t103 * t125 + t108 * t116 + t127 * t97 + t183 * t98) * MDP(18) + (-t100 * t183 + t101 * t122 - t110 * t115 - t111 * t116 - t154) * MDP(19) + (-t103 * t124 + t108 * t115 - t122 * t98 + t126 * t97) * MDP(20) + (t100 * t102 + t101 * t99 + t103 * t98 + t108 * t97 + t110 * t96 + t111 * t95) * MDP(21) + (t124 * MDP(13) - t125 * MDP(14) - t100 * t181 - t167 * t101) * qJD(4) + ((t127 * MDP(16) - t126 * MDP(17) + 0.2e1 * t177 + 0.2e1 * t176 + (2 * MDP(5)) + (MDP(10) + (2 * MDP(6))) * qJ(2)) * qJD(2) + (-MDP(10) * t150 + (2 * MDP(9))) * t184) * qJD(1); (-t116 * t127 - t125 * t183 - t159) * MDP(19) + t154 * MDP(21) + (t167 * t124 - t125 * t181) * qJD(4) + ((-t142 - t184) * MDP(10) - t103 * MDP(21) - t181 * t122 - t167 * t183) * qJD(1) + (-qJ(2) * MDP(6) - MDP(5) - t176 - t177) * t153; (t174 * t134 * qJD(1) + t146) * MDP(10) + (-t119 - t180) * MDP(19) + (t102 * t183 + (-qJD(5) - t99) * t122 + t158) * MDP(21) - t174 * MDP(9) * t153 - t167 * t133 + (-0.2e1 * t181 * t183 + t167 * (t122 + t165)) * qJD(4); (-t119 + t180) * MDP(12) + t133 * MDP(14) + (-t122 * t130 - t96) * MDP(16) - t114 * MDP(17) + t157 * MDP(18) + (pkin(4) * t115 - qJ(5) * t116 + (t102 - t105) * t122) * MDP(19) + (t109 * t122 + t114) * MDP(20) + (-pkin(4) * t96 + qJ(5) * t95 + t168 * t102 - t103 * t109 - t105 * t99) * MDP(21) + (t122 * MDP(11) + t130 * MDP(17) - t109 * MDP(18) + (t99 - t168) * MDP(19) - t103 * MDP(20)) * t183 + ((t122 - t165) * MDP(14) + t161 * MDP(17) + (0.2e1 * qJD(5) - t161) * MDP(20) + t167 * t105) * qJD(4) + t181 * t182; t122 * t183 * MDP(18) + (-qJD(4) ^ 2 - t180) * MDP(20) + (-qJD(4) * t102 - t157) * MDP(21);];
tauc = t1;
