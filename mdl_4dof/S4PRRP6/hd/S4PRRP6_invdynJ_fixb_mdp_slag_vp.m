% Calculate vector of inverse dynamics joint torques for
% S4PRRP6
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [4x1]
%   Generalized joint coordinates (joint angles)
% qJD [4x1]
%   Generalized joint velocities
% qJDD [4x1]
%   Generalized joint accelerations
% g [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [6x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,d2,d3,theta1]';
% MDP [15x1]
%   Minimal dynamic parameter vector (fixed base model)
%   see S4PRRP6_convert_par2_MPV_fixb.m
% 
% Output:
% tau [4x1]
%   joint torques of inverse dynamics (contains inertial, gravitational Coriolis and centrifugal forces)

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-31 16:30
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tau = S4PRRP6_invdynJ_fixb_mdp_slag_vp(qJ, qJD, qJDD, g, pkin, MDP)
%% Coder Information
%#codegen
%$cgargs {zeros(4,1),zeros(4,1),zeros(4,1),zeros(3,1),zeros(6,1),zeros(15,1)}
assert(isreal(qJ) && all(size(qJ) == [4 1]), ...
  'S4PRRP6_invdynJ_fixb_mdp_slag_vp: qJ has to be [4x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [4 1]), ...
  'S4PRRP6_invdynJ_fixb_mdp_slag_vp: qJD has to be [4x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [4 1]), ...
  'S4PRRP6_invdynJ_fixb_mdp_slag_vp: qJDD has to be [4x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S4PRRP6_invdynJ_fixb_mdp_slag_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [6 1]), ...
  'S4PRRP6_invdynJ_fixb_mdp_slag_vp: pkin has to be [6x1] (double)');
assert(isreal(MDP) && all(size(MDP) == [15 1]), ...
  'S4PRRP6_invdynJ_fixb_mdp_slag_vp: MDP has to be [15x1] (double)'); 

%% Symbolic Calculation
% From invdyn_joint_fixb_mdp_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 16:30:45
% EndTime: 2019-12-31 16:30:46
% DurationCPUTime: 0.77s
% Computational Cost: add. (392->158), mult. (843->208), div. (0->0), fcn. (502->6), ass. (0->79)
t125 = sin(pkin(6));
t126 = cos(pkin(6));
t144 = g(1) * t126 + g(2) * t125;
t127 = sin(qJ(3));
t129 = cos(qJ(3));
t140 = t129 * pkin(3) + t127 * qJ(4) + pkin(2);
t187 = t140 * qJD(2);
t143 = pkin(3) * t127 - qJ(4) * t129;
t101 = t143 * qJD(3) - t127 * qJD(4);
t186 = qJD(2) * t101 - t140 * qJDD(2);
t128 = sin(qJ(2));
t158 = qJD(1) * qJD(2);
t185 = (t144 + t158) * t128;
t117 = qJD(2) * pkin(5) + t128 * qJD(1);
t161 = qJD(3) * t129;
t175 = qJDD(3) * pkin(3);
t184 = t117 * t161 - t175;
t181 = g(3) * t128;
t180 = qJD(2) * pkin(2);
t179 = qJD(3) * pkin(3);
t178 = pkin(5) * qJDD(3);
t176 = qJDD(2) * pkin(5);
t162 = qJD(3) * qJ(4);
t170 = t129 * t117;
t107 = t162 + t170;
t174 = t107 * t129;
t173 = t127 * t117;
t130 = cos(qJ(2));
t172 = t127 * t130;
t171 = t128 * t129;
t169 = t129 * t130;
t168 = qJDD(1) - g(3);
t123 = t127 ^ 2;
t124 = t129 ^ 2;
t167 = t123 - t124;
t166 = t123 + t124;
t131 = qJD(3) ^ 2;
t132 = qJD(2) ^ 2;
t165 = t131 + t132;
t163 = qJD(2) * t127;
t160 = t107 * qJD(3);
t159 = t130 * qJD(1);
t157 = qJD(2) * qJD(3);
t156 = qJDD(2) * t127;
t155 = qJDD(2) * t129;
t154 = qJDD(3) * qJ(4);
t153 = qJDD(3) * t127;
t152 = t127 * t132 * t129;
t151 = t127 * t157;
t122 = t128 * t158;
t118 = -t159 - t180;
t150 = t118 - t180;
t98 = -t159 - t187;
t149 = t98 - t187;
t148 = qJD(4) + t173;
t147 = -t130 * qJDD(1) + t122;
t146 = t127 * qJD(3) * t159 + t129 * t122 + t144 * t171;
t145 = pkin(5) * t131 + g(3) * t130;
t109 = t128 * qJDD(1) + t130 * t158 + t176;
t100 = t129 * t109;
t96 = t154 + t100 + (qJD(4) - t173) * qJD(3);
t99 = t127 * t109;
t97 = qJDD(4) + t99 + t184;
t142 = t97 * t127 + t96 * t129;
t141 = t166 * MDP(13) - MDP(4);
t104 = t125 * t169 - t126 * t127;
t106 = t125 * t127 + t126 * t169;
t139 = g(1) * t106 + g(2) * t104 - t100;
t103 = t125 * t172 + t126 * t129;
t105 = -t125 * t129 + t126 * t172;
t138 = g(1) * t105 + g(2) * t103 + t127 * t181 - t99;
t136 = -qJDD(4) + t138;
t135 = 0.2e1 * qJDD(2) * pkin(2) - t145 - t147;
t93 = t147 + t186;
t134 = -t145 - t186 - t93;
t102 = t148 - t179;
t133 = (t102 * t129 - t107 * t127) * qJD(3) + t142;
t111 = t143 * qJD(2);
t1 = [t168 * MDP(1) - g(3) * MDP(15) + (MDP(10) + MDP(12)) * ((-0.2e1 * t151 + t155) * t130 + (-t165 * t129 - t153) * t128) + (MDP(11) - MDP(14)) * ((-qJDD(3) * t128 - 0.2e1 * t130 * t157) * t129 + (-t130 * qJDD(2) + t165 * t128) * t127) + (qJDD(2) * MDP(3) + (qJD(2) * t174 + t102 * t163 - t93) * MDP(15) + t141 * t132) * t130 + (-t132 * MDP(3) + (qJD(2) * t98 + t102 * t161 - t127 * t160 + t142) * MDP(15) + t141 * qJDD(2)) * t128; qJDD(2) * MDP(2) + (t144 * t128 + t168 * t130) * MDP(3) + (-t168 * t128 + t144 * t130) * MDP(4) + (t123 * qJDD(2) + 0.2e1 * t129 * t151) * MDP(5) + 0.2e1 * (t127 * t155 - t167 * t157) * MDP(6) + (t131 * t129 + t153) * MDP(7) + (qJDD(3) * t129 - t131 * t127) * MDP(8) + ((t150 * qJD(3) - t178) * t127 + t135 * t129 + t146) * MDP(10) + ((-t178 + (t150 + t159) * qJD(3)) * t129 + (-t135 - t185) * t127) * MDP(11) + ((t149 * qJD(3) - t178) * t127 + t134 * t129 + t146) * MDP(12) + (-t181 + t166 * t176 + (-t166 * t158 - t144) * t130 + t133) * MDP(13) + ((t178 + (-t149 - t159) * qJD(3)) * t129 + (t134 + t185) * t127) * MDP(14) + (t98 * t101 - t93 * t140 + t133 * pkin(5) + (-g(3) * t140 + (-t102 * t127 - t174) * qJD(1) - t144 * pkin(5)) * t130 + (-g(3) * pkin(5) - t98 * qJD(1) + t144 * t140) * t128) * MDP(15); -MDP(5) * t152 + t167 * MDP(6) * t132 + MDP(7) * t156 + MDP(8) * t155 + qJDD(3) * MDP(9) + (-t118 * t163 + t138) * MDP(10) + ((-qJD(2) * t118 + t181) * t129 + t139) * MDP(11) + (0.2e1 * t175 + (t111 * t129 - t127 * t98) * qJD(2) + t136) * MDP(12) + (-t143 * qJDD(2) + ((t107 - t162) * t127 + (qJD(4) - t102 - t179) * t129) * qJD(2)) * MDP(13) + (-g(3) * t171 + 0.2e1 * t154 + 0.2e1 * qJD(3) * qJD(4) + (t111 * t127 + t129 * t98) * qJD(2) - t139) * MDP(14) + (t96 * qJ(4) - t97 * pkin(3) - t98 * t111 - t102 * t170 - g(1) * (-t105 * pkin(3) + t106 * qJ(4)) - g(2) * (-t103 * pkin(3) + t104 * qJ(4)) + t143 * t181 + t148 * t107) * MDP(15); (-qJDD(3) - t152) * MDP(12) + MDP(13) * t156 + (-t123 * t132 - t131) * MDP(14) + (t98 * t163 - t136 - t160 + t184) * MDP(15);];
tau = t1;
