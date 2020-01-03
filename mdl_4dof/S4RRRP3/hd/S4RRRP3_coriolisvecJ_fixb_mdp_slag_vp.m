% Calculate minimal parameter regressor of Coriolis joint torque vector for
% S4RRRP3
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
%   see S4RRRP3_convert_par2_MPV_fixb.m
% 
% Output:
% tauc [4x1]
%   joint torques required to compensate Coriolis and centrifugal load

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-31 17:14
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauc = S4RRRP3_coriolisvecJ_fixb_mdp_slag_vp(qJ, qJD, pkin, MDP)
%% Coder Information
%#codegen
%$cgargs {zeros(4,1),zeros(4,1),zeros(6,1),zeros(17,1)}
assert(isreal(qJ) && all(size(qJ) == [4 1]), ...
  'S4RRRP3_coriolisvecJ_fixb_mdp_slag_vp: qJ has to be [4x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [4 1]), ...
  'S4RRRP3_coriolisvecJ_fixb_mdp_slag_vp: qJD has to be [4x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [6 1]), ...
  'S4RRRP3_coriolisvecJ_fixb_mdp_slag_vp: pkin has to be [6x1] (double)');
assert(isreal(MDP) && all(size(MDP) == [17 1]), ...
  'S4RRRP3_coriolisvecJ_fixb_mdp_slag_vp: MDP has to be [17x1] (double)'); 

%% Symbolic Calculation
% From coriolisvec_joint_fixb_mdp_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 17:14:16
% EndTime: 2019-12-31 17:14:18
% DurationCPUTime: 0.56s
% Computational Cost: add. (440->117), mult. (817->165), div. (0->0), fcn. (347->4), ass. (0->62)
t120 = sin(qJ(3));
t118 = t120 ^ 2;
t122 = cos(qJ(3));
t119 = t122 ^ 2;
t154 = t120 * t122;
t168 = MDP(7) * t154 - (t118 - t119) * MDP(8);
t117 = qJD(1) + qJD(2);
t132 = -qJD(3) * pkin(3) + qJD(4);
t121 = sin(qJ(2));
t159 = pkin(1) * qJD(1);
t138 = t121 * t159;
t105 = pkin(6) * t117 + t138;
t156 = t105 * t120;
t95 = t132 + t156;
t147 = qJD(3) * qJ(4);
t96 = t105 * t122 + t147;
t167 = (-MDP(6) + (t118 + t119) * MDP(15)) * t117 + (t120 * t95 + t122 * t96) * MDP(17);
t123 = cos(qJ(2));
t158 = pkin(1) * qJD(2);
t136 = qJD(1) * t158;
t129 = t123 * t136;
t109 = t120 * t129;
t145 = qJD(3) * t122;
t91 = t105 * t145 + t109;
t166 = -qJD(3) * t96 + t91;
t164 = MDP(12) * t122 + MDP(5);
t161 = 0.2e1 * qJD(3);
t160 = pkin(1) * t123;
t107 = -pkin(3) * t122 - qJ(4) * t120 - pkin(2);
t155 = t107 * t117;
t137 = t123 * t159;
t146 = qJD(3) * t120;
t153 = t122 * t117 * t138 + t137 * t146;
t149 = MDP(13) * t117;
t113 = pkin(1) * t121 + pkin(6);
t148 = MDP(17) * t113;
t144 = t120 * qJD(4);
t124 = qJD(3) ^ 2;
t143 = t124 * MDP(10);
t141 = -MDP(12) - MDP(14);
t140 = MDP(13) - MDP(16);
t127 = pkin(3) * t120 - qJ(4) * t122;
t128 = t121 * t136;
t85 = t128 + (t127 * qJD(3) - t144) * t117;
t97 = pkin(3) * t146 - qJ(4) * t145 - t144;
t93 = t121 * t158 + t97;
t134 = -t117 * t93 - t85;
t133 = -t117 * t97 - t85;
t131 = t141 * t124;
t130 = t140 * t124;
t106 = -pkin(2) * t117 - t137;
t110 = t122 * t129;
t88 = t110 + (qJD(4) - t156) * qJD(3);
t125 = t124 * t122 * MDP(9) + (t106 * t145 + t120 * t128) * MDP(13) + (t91 * t120 + t88 * t122 + t95 * t145) * MDP(15) + t168 * t117 * t161;
t116 = t117 ^ 2;
t114 = -pkin(2) - t160;
t102 = t107 - t160;
t100 = t106 * t146;
t98 = t127 * t117;
t90 = -t137 + t155;
t87 = t90 * t146;
t1 = [t100 * MDP(12) + t87 * MDP(14) + (t102 * t85 + t90 * t93) * MDP(17) + (t134 * MDP(14) + (t88 * MDP(17) + t131) * t113 + (t114 * t149 + (-t102 * t117 - t90) * MDP(16) + t95 * t148) * qJD(3)) * t122 + (-t143 + t134 * MDP(16) + (t91 * MDP(17) + t130) * t113 + ((-MDP(15) - t148) * t96 + (MDP(12) * t114 + MDP(14) * t102) * t117) * qJD(3)) * t120 + ((t120 * t149 + t164 * (-qJD(1) - t117)) * t121 + (-qJD(1) * MDP(6) + (t141 * t120 - t140 * t122) * qJD(3) + t167) * t123) * t158 + t125; (t100 + t153) * MDP(12) + (t87 + t153) * MDP(14) + (t107 * t85 + t90 * t97) * MDP(17) + (t133 * MDP(14) + (-pkin(2) * t149 + (-t90 - t155) * MDP(16)) * qJD(3) + ((qJD(3) * t95 + t88) * MDP(17) + t131) * pkin(6)) * t122 + (-t143 + t133 * MDP(16) + (-t96 * MDP(15) + (-MDP(12) * pkin(2) + MDP(14) * t107) * t117) * qJD(3) + (MDP(17) * t166 + t130) * pkin(6)) * t120 + ((-t90 * MDP(17) - t164 * qJD(2) + (-t140 * t120 + MDP(5)) * t117) * t121 + (-qJD(2) * MDP(6) + t140 * t145 - t167) * t123) * t159 + t125; -t110 * MDP(13) + (qJD(4) * t161 + t110) * MDP(16) + (-t91 * pkin(3) + t88 * qJ(4) + t96 * qJD(4) - t90 * t98 + (t120 * t96 - t122 * t95) * t105) * MDP(17) + t141 * t109 - t168 * t116 + ((-t106 * MDP(12) - t90 * MDP(14) + (t96 - t147) * MDP(15) + t98 * MDP(16)) * t120 + (-t106 * MDP(13) + t98 * MDP(14) + (t132 - t95) * MDP(15) + t90 * MDP(16)) * t122) * t117; -t116 * MDP(14) * t154 + (-t116 * t118 - t124) * MDP(16) + (t117 * t120 * t90 + t166) * MDP(17);];
tauc = t1;
