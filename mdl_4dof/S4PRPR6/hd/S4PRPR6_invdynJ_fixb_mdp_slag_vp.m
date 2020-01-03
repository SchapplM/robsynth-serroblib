% Calculate vector of inverse dynamics joint torques for
% S4PRPR6
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
% pkin [7x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,d2,d4,theta1,theta3]';
% MDP [15x1]
%   Minimal dynamic parameter vector (fixed base model)
%   see S4PRPR6_convert_par2_MPV_fixb.m
% 
% Output:
% tau [4x1]
%   joint torques of inverse dynamics (contains inertial, gravitational Coriolis and centrifugal forces)

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-31 16:24
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tau = S4PRPR6_invdynJ_fixb_mdp_slag_vp(qJ, qJD, qJDD, g, pkin, MDP)
%% Coder Information
%#codegen
%$cgargs {zeros(4,1),zeros(4,1),zeros(4,1),zeros(3,1),zeros(7,1),zeros(15,1)}
assert(isreal(qJ) && all(size(qJ) == [4 1]), ...
  'S4PRPR6_invdynJ_fixb_mdp_slag_vp: qJ has to be [4x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [4 1]), ...
  'S4PRPR6_invdynJ_fixb_mdp_slag_vp: qJD has to be [4x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [4 1]), ...
  'S4PRPR6_invdynJ_fixb_mdp_slag_vp: qJDD has to be [4x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S4PRPR6_invdynJ_fixb_mdp_slag_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [7 1]), ...
  'S4PRPR6_invdynJ_fixb_mdp_slag_vp: pkin has to be [7x1] (double)');
assert(isreal(MDP) && all(size(MDP) == [15 1]), ...
  'S4PRPR6_invdynJ_fixb_mdp_slag_vp: MDP has to be [15x1] (double)'); 

%% Symbolic Calculation
% From invdyn_joint_fixb_mdp_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 16:24:44
% EndTime: 2019-12-31 16:24:45
% DurationCPUTime: 0.81s
% Computational Cost: add. (353->137), mult. (781->190), div. (0->0), fcn. (553->10), ass. (0->73)
t131 = sin(qJ(2));
t127 = sin(pkin(6));
t129 = cos(pkin(6));
t148 = g(1) * t129 + g(2) * t127;
t141 = t148 * t131;
t133 = cos(qJ(2));
t175 = g(3) * t133;
t137 = t141 - t175;
t126 = sin(pkin(7));
t128 = cos(pkin(7));
t170 = t128 * MDP(5);
t184 = -t126 * MDP(6) + t170;
t130 = sin(qJ(4));
t132 = cos(qJ(4));
t108 = t126 * t132 + t128 * t130;
t181 = t108 * qJD(4);
t183 = qJD(2) * t181;
t182 = t148 * t133;
t167 = t126 ^ 2 + t128 ^ 2;
t165 = qJD(1) * t133;
t149 = qJD(3) - t165;
t179 = qJD(4) ^ 2;
t176 = g(3) * t131;
t174 = pkin(5) + qJ(3);
t173 = qJDD(2) * pkin(2);
t171 = t127 * t133;
t169 = t129 * t133;
t168 = qJDD(1) - g(3);
t166 = qJD(1) * t131;
t163 = qJD(2) * t130;
t162 = qJD(2) * t132;
t161 = qJD(1) * qJD(2);
t160 = qJDD(2) * t130;
t159 = qJDD(2) * t132;
t107 = t126 * t130 - t132 * t128;
t158 = qJDD(4) * t107;
t157 = qJDD(4) * t108;
t154 = t128 * t162;
t156 = qJD(4) * t154 + t126 * t159 + t128 * t160;
t119 = -pkin(3) * t128 - pkin(2);
t155 = t126 * t163;
t99 = qJDD(2) * qJ(3) + qJDD(1) * t131 + (qJD(3) + t165) * qJD(2);
t153 = t167 * t99;
t152 = pkin(5) * qJDD(2) + t99;
t151 = t167 * t133;
t150 = t167 * qJDD(2);
t116 = t128 * t159;
t147 = -t126 * t160 + t116;
t109 = t174 * t126;
t110 = t174 * t128;
t146 = -t109 * t132 - t110 * t130;
t145 = -t109 * t130 + t110 * t132;
t144 = -qJDD(1) * t133 + t131 * t161 + qJDD(3);
t104 = t107 * qJD(4);
t139 = -t126 * t162 - t128 * t163;
t136 = t149 * t167;
t135 = t153 - t176 - t182;
t134 = qJD(2) ^ 2;
t125 = pkin(7) + qJ(4);
t122 = cos(t125);
t121 = sin(t125);
t113 = qJD(2) * qJ(3) + t166;
t111 = -qJD(2) * pkin(2) + t149;
t106 = t119 * qJD(2) + t149;
t103 = t108 * qJD(2);
t101 = -t154 + t155;
t100 = t144 - t173;
t96 = t119 * qJDD(2) + t144;
t95 = t152 * t128;
t94 = t152 * t126;
t93 = -t147 + t183;
t92 = -qJD(4) * t155 + t156;
t1 = [t168 * MDP(1) + (-qJDD(2) * t131 - t133 * t134) * MDP(4) + (t131 * t150 + t134 * t151) * MDP(7) + (-t100 * t133 - g(3) + t131 * t153 + (t111 * t131 + t113 * t151) * qJD(2)) * MDP(8) + ((-t93 - t183) * t133 + (qJD(2) * t101 + t107 * t179 - t157) * t131) * MDP(14) + ((qJD(2) * t104 - t92) * t133 + (qJD(2) * t103 + t108 * t179 + t158) * t131) * MDP(15) + (MDP(3) + t184) * (qJDD(2) * t133 - t131 * t134); qJDD(2) * MDP(2) + (t168 * t133 + t141) * MDP(3) + (-t168 * t131 + t182) * MDP(4) + (qJ(3) * t150 + t136 * qJD(2) + t135) * MDP(7) + (-t111 * t166 + (-t100 + t137) * pkin(2) + t135 * qJ(3) + t136 * t113) * MDP(8) + (-t103 * t104 + t108 * t92) * MDP(9) + (t101 * t104 - t103 * t181 - t107 * t92 - t108 * t93) * MDP(10) + (-qJD(4) * t104 + t157) * MDP(11) + (-qJD(4) * t181 - t158) * MDP(12) + (t146 * qJDD(4) + t119 * t93 + t96 * t107 + t106 * t181 + (-t108 * qJD(3) - t145 * qJD(4)) * qJD(4) + t137 * t122 + (-t131 * t101 + t133 * t181) * qJD(1)) * MDP(14) + (-t145 * qJDD(4) + t119 * t92 + t96 * t108 - t106 * t104 + (t107 * qJD(3) - t146 * qJD(4)) * qJD(4) - t137 * t121 + (-t131 * t103 - t133 * t104) * qJD(1)) * MDP(15) + t184 * ((t148 + t161) * t131 - t100 + t173 - t175); (-t167 * qJD(2) * t113 - t137 + t144) * MDP(8) - t116 * MDP(14) + t156 * MDP(15) - t167 * MDP(7) * t134 + (-t170 - pkin(2) * MDP(8) + (t130 * MDP(14) + MDP(6)) * t126) * qJDD(2) + ((t103 - t139) * MDP(14) + (-t101 - t155) * MDP(15)) * qJD(4); t103 * t101 * MDP(9) + (-t101 ^ 2 + t103 ^ 2) * MDP(10) + t156 * MDP(11) + t147 * MDP(12) + qJDD(4) * MDP(13) + (-t130 * t95 - t132 * t94 - t106 * t103 - g(1) * (-t121 * t169 + t122 * t127) - g(2) * (-t121 * t171 - t122 * t129) + t121 * t176) * MDP(14) + (-t132 * t95 + t130 * t94 + t106 * t101 - g(1) * (-t121 * t127 - t122 * t169) - g(2) * (t121 * t129 - t122 * t171) + t122 * t176) * MDP(15) + ((t101 - t155) * MDP(11) + (t103 + t139) * MDP(12)) * qJD(4);];
tau = t1;
