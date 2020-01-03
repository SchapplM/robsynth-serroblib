% Calculate vector of inverse dynamics joint torques for
% S4RPRR2
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
%   pkin=[a2,a3,a4,d1,d3,d4,theta2]';
% MDP [14x1]
%   Minimal dynamic parameter vector (fixed base model)
%   see S4RPRR2_convert_par2_MPV_fixb.m
% 
% Output:
% tau [4x1]
%   joint torques of inverse dynamics (contains inertial, gravitational Coriolis and centrifugal forces)

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-31 16:48
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tau = S4RPRR2_invdynJ_fixb_mdp_slag_vp(qJ, qJD, qJDD, g, pkin, MDP)
%% Coder Information
%#codegen
%$cgargs {zeros(4,1),zeros(4,1),zeros(4,1),zeros(3,1),zeros(7,1),zeros(14,1)}
assert(isreal(qJ) && all(size(qJ) == [4 1]), ...
  'S4RPRR2_invdynJ_fixb_mdp_slag_vp: qJ has to be [4x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [4 1]), ...
  'S4RPRR2_invdynJ_fixb_mdp_slag_vp: qJD has to be [4x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [4 1]), ...
  'S4RPRR2_invdynJ_fixb_mdp_slag_vp: qJDD has to be [4x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S4RPRR2_invdynJ_fixb_mdp_slag_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [7 1]), ...
  'S4RPRR2_invdynJ_fixb_mdp_slag_vp: pkin has to be [7x1] (double)');
assert(isreal(MDP) && all(size(MDP) == [14 1]), ...
  'S4RPRR2_invdynJ_fixb_mdp_slag_vp: MDP has to be [14x1] (double)'); 

%% Symbolic Calculation
% From invdyn_joint_fixb_mdp_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 16:48:15
% EndTime: 2019-12-31 16:48:15
% DurationCPUTime: 0.35s
% Computational Cost: add. (346->91), mult. (605->124), div. (0->0), fcn. (328->10), ass. (0->57)
t124 = cos(pkin(7));
t115 = t124 * pkin(1) + pkin(2);
t105 = t115 * qJDD(1);
t117 = qJ(1) + pkin(7) + qJ(3);
t114 = cos(t117);
t126 = sin(qJ(3));
t129 = cos(qJ(3));
t107 = t115 * qJD(1);
t154 = qJD(3) * t107;
t167 = -g(2) * t114 + t129 * t105 - t126 * t154;
t119 = qJDD(1) + qJDD(3);
t162 = t119 * pkin(3);
t148 = qJD(1) * qJD(3) * t129;
t123 = sin(pkin(7));
t164 = pkin(1) * t123;
t165 = (qJDD(1) * t126 + t148) * t164;
t166 = -t162 + t165 - t167;
t156 = t126 * t115 + t129 * t164;
t113 = sin(t117);
t150 = qJD(1) * t164;
t146 = t126 * t150;
t133 = g(1) * t114 + g(2) * t113 - (qJDD(1) * t164 + t154) * t129 + qJD(3) * t146 - t126 * t105;
t110 = g(1) * t113;
t120 = qJD(1) + qJD(3);
t161 = (t126 * t107 + t129 * t150) * t120;
t160 = t156 * qJD(3) * t120;
t125 = sin(qJ(4));
t128 = cos(qJ(4));
t91 = t129 * t107 - t146;
t89 = -t120 * pkin(3) - t91;
t159 = t89 * qJD(4) * t125 + t128 * t110;
t158 = t129 * t115;
t157 = qJDD(2) - g(3);
t121 = t125 ^ 2;
t155 = -t128 ^ 2 + t121;
t153 = qJD(4) * t120;
t152 = t128 * qJD(4);
t151 = t166 * t125 + t89 * t152;
t131 = qJD(4) ^ 2;
t101 = qJDD(4) * t125 + t131 * t128;
t102 = qJDD(4) * t128 - t131 * t125;
t145 = t102 * MDP(11) + t119 * MDP(5) + 0.2e1 * (t125 * t119 * t128 - t155 * t153) * MDP(9) + (0.2e1 * t125 * t120 * t152 + t121 * t119) * MDP(8) + t101 * MDP(10);
t127 = sin(qJ(1));
t130 = cos(qJ(1));
t144 = g(1) * t127 - g(2) * t130;
t142 = -t126 * t164 + t158;
t139 = pkin(6) * t131 - t161 - t162;
t96 = -pkin(3) - t142;
t97 = pkin(6) + t156;
t138 = t119 * t96 + t131 * t97 + t160;
t137 = -t119 * pkin(6) - t89 * t120 + t133;
t136 = t110 + t167;
t94 = t142 * qJD(3);
t135 = -qJDD(4) * t97 + (t120 * t96 - t94) * qJD(4);
t134 = -pkin(3) * t153 - pkin(6) * qJDD(4) + qJD(4) * t91;
t118 = t120 ^ 2;
t1 = [qJDD(1) * MDP(1) + t144 * MDP(2) + (g(1) * t130 + g(2) * t127) * MDP(3) + (t144 + (t123 ^ 2 + t124 ^ 2) * qJDD(1) * pkin(1)) * pkin(1) * MDP(4) + (t119 * t158 - t160 + (-t148 + (-qJDD(1) - t119) * t126) * t164 + t136) * MDP(6) + (-t156 * t119 - t94 * t120 + t133) * MDP(7) + (t135 * t125 + (-t138 - t166) * t128 + t159) * MDP(13) + (t135 * t128 + (t138 - t110) * t125 + t151) * MDP(14) + t145; t102 * MDP(13) - t101 * MDP(14) + t157 * MDP(4); (t136 + t161 - t165) * MDP(6) + (t91 * t120 + t133) * MDP(7) + t159 * MDP(13) + t151 * MDP(14) + (t134 * MDP(13) + (t139 - t110) * MDP(14)) * t125 + ((-t139 - t166) * MDP(13) + t134 * MDP(14)) * t128 + t145; qJDD(4) * MDP(12) + t155 * MDP(9) * t118 + (t119 * MDP(11) + t157 * MDP(13) + t137 * MDP(14)) * t128 + (-t118 * t128 * MDP(8) + t119 * MDP(10) + t137 * MDP(13) - t157 * MDP(14)) * t125;];
tau = t1;
