% Calculate minimal parameter regressor of Coriolis joint torque vector for
% S4RRRR1
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [4x1]
%   Generalized joint coordinates (joint angles)
% qJD [4x1]
%   Generalized joint velocities
% pkin [7x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,d1,d2,d3,d4]';
% MDP [16x1]
%   Minimal dynamic parameter vector (fixed base model)
%   see S4RRRR1_convert_par2_MPV_fixb.m
% 
% Output:
% tauc [4x1]
%   joint torques required to compensate Coriolis and centrifugal load

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-31 17:22
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauc = S4RRRR1_coriolisvecJ_fixb_mdp_slag_vp(qJ, qJD, pkin, MDP)
%% Coder Information
%#codegen
%$cgargs {zeros(4,1),zeros(4,1),zeros(7,1),zeros(16,1)}
assert(isreal(qJ) && all(size(qJ) == [4 1]), ...
  'S4RRRR1_coriolisvecJ_fixb_mdp_slag_vp: qJ has to be [4x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [4 1]), ...
  'S4RRRR1_coriolisvecJ_fixb_mdp_slag_vp: qJD has to be [4x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [7 1]), ...
  'S4RRRR1_coriolisvecJ_fixb_mdp_slag_vp: pkin has to be [7x1] (double)');
assert(isreal(MDP) && all(size(MDP) == [16 1]), ...
  'S4RRRR1_coriolisvecJ_fixb_mdp_slag_vp: MDP has to be [16x1] (double)'); 

%% Symbolic Calculation
% From coriolisvec_joint_fixb_mdp_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 17:22:15
% EndTime: 2019-12-31 17:22:16
% DurationCPUTime: 0.41s
% Computational Cost: add. (338->90), mult. (697->135), div. (0->0), fcn. (340->6), ass. (0->56)
t105 = sin(qJ(4));
t108 = cos(qJ(4));
t157 = t105 * t108 * MDP(10) - (t105 ^ 2 - t108 ^ 2) * MDP(11);
t134 = -qJD(2) - qJD(3);
t106 = sin(qJ(3));
t110 = cos(qJ(2));
t107 = sin(qJ(2));
t109 = cos(qJ(3));
t146 = t107 * t109;
t117 = t106 * t110 + t146;
t141 = qJD(3) * t106;
t150 = pkin(1) * qJD(1);
t156 = -pkin(2) * t141 + t117 * t150;
t155 = t105 * MDP(15) + t108 * MDP(16);
t101 = qJD(1) - t134;
t153 = t101 * pkin(3);
t138 = qJD(4) * t108;
t140 = qJD(3) * t109;
t127 = t107 * t140;
t112 = (t117 * qJD(2) + t127) * pkin(1);
t102 = qJD(1) + qJD(2);
t129 = t110 * t150;
t93 = t102 * pkin(2) + t129;
t128 = t93 * t141;
t79 = qJD(1) * t112 + t128;
t131 = t107 * t150;
t124 = t106 * t131;
t86 = t109 * t93 - t124;
t84 = -t86 - t153;
t152 = t79 * t105 + t84 * t138;
t151 = t134 * t124;
t99 = t110 * pkin(1) + pkin(2);
t149 = (t99 * t141 + t112) * t101;
t148 = (t106 * t93 + t109 * t131) * t101;
t147 = t106 * t107;
t111 = qJD(4) ^ 2;
t145 = t108 * t111;
t144 = t109 * MDP(8);
t143 = t109 * t110;
t139 = qJD(4) * t105;
t135 = t111 * MDP(13);
t133 = 0.2e1 * t157 * qJD(4) * t101 + MDP(12) * t145;
t123 = MDP(8) * t127;
t121 = pkin(7) * t111 - t148;
t120 = qJD(4) * (t86 - t153);
t119 = t111 * (pkin(1) * t146 + t106 * t99 + pkin(7)) + t149;
t116 = t143 - t147;
t80 = t99 * t140 + (t116 * qJD(2) - t107 * t141) * pkin(1);
t118 = qJD(4) * (t101 * (pkin(1) * t147 - t109 * t99 - pkin(3)) - t80);
t115 = -t93 * t140 - t151;
t114 = -t106 * MDP(8) - t109 * MDP(9) - MDP(6);
t82 = t84 * t139;
t113 = t82 * MDP(15) + t152 * MDP(16) + t133;
t98 = -t109 * pkin(2) - pkin(3);
t97 = t106 * pkin(2) + pkin(7);
t1 = [(-t128 - t149) * MDP(8) + (-t80 * t101 + t115) * MDP(9) + ((-t119 - t79) * MDP(15) + MDP(16) * t118) * t108 + (MDP(15) * t118 + t119 * MDP(16) - t135) * t105 + (-qJD(1) * t123 + ((-t107 * MDP(5) - t110 * MDP(6)) * t102 + ((-MDP(5) - t144) * t107 + t114 * t110) * qJD(1)) * qJD(2)) * pkin(1) + t113; -MDP(8) * t128 + t115 * MDP(9) - t105 * t135 + (-t79 * t108 - t97 * t145 + t82) * MDP(15) + (t111 * t105 * t97 + t152) * MDP(16) + (t156 * MDP(8) + (t156 * t108 + t98 * t139) * MDP(15) + (-t156 * t105 + t98 * t138) * MDP(16)) * t101 + ((t102 * MDP(6) + t114 * qJD(2)) * t110 + ((-qJD(2) + t102) * MDP(5) + t134 * t144) * t107) * t150 + t133 + (MDP(9) * t101 + t155 * qJD(4)) * (-pkin(2) * t140 + t116 * t150); (-t128 + t148) * MDP(8) + (t86 * t101 + t115) * MDP(9) + (-t123 + (-t117 * MDP(8) - MDP(9) * t143) * qJD(2)) * t150 + ((-t121 - t79) * MDP(15) + MDP(16) * t120) * t108 + (MDP(15) * t120 + t121 * MDP(16) - t135) * t105 + t113; -t157 * t101 ^ 2 + t155 * (-t101 * t84 - (qJD(2) * t129 + qJD(3) * t93) * t109 - t151);];
tauc = t1;
