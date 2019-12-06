% Calculate vector of inverse dynamics joint torques for
% S5PPPRR1
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% qJD [5x1]
%   Generalized joint velocities
% qJDD [5x1]
%   Generalized joint accelerations
% g [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [9x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d4,d5,theta1,theta2,theta3]';
% MDP [13x1]
%   Minimal dynamic parameter vector (fixed base model)
%   see S5PPPRR1_convert_par2_MPV_fixb.m
% 
% Output:
% tau [5x1]
%   joint torques of inverse dynamics (contains inertial, gravitational Coriolis and centrifugal forces)

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-05 14:58
% Revision: 77da58f92bca3eff71542919beafa37024070d86 (2019-12-05)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tau = S5PPPRR1_invdynJ_fixb_mdp_slag_vp(qJ, qJD, qJDD, g, pkin, MDP)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(5,1),zeros(3,1),zeros(9,1),zeros(13,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5PPPRR1_invdynJ_fixb_mdp_slag_vp: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5PPPRR1_invdynJ_fixb_mdp_slag_vp: qJD has to be [5x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [5 1]), ...
  'S5PPPRR1_invdynJ_fixb_mdp_slag_vp: qJDD has to be [5x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S5PPPRR1_invdynJ_fixb_mdp_slag_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S5PPPRR1_invdynJ_fixb_mdp_slag_vp: pkin has to be [9x1] (double)');
assert(isreal(MDP) && all(size(MDP) == [13 1]), ...
  'S5PPPRR1_invdynJ_fixb_mdp_slag_vp: MDP has to be [13x1] (double)'); 

%% Symbolic Calculation
% From invdyn_joint_fixb_mdp_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-05 14:58:10
% EndTime: 2019-12-05 14:58:11
% DurationCPUTime: 0.45s
% Computational Cost: add. (311->92), mult. (672->149), div. (0->0), fcn. (589->12), ass. (0->56)
t117 = sin(pkin(8));
t118 = sin(pkin(7));
t121 = cos(pkin(7));
t166 = (g(1) * t121 + g(2) * t118) * t117;
t127 = qJD(4) ^ 2;
t165 = 2 * qJDD(4);
t122 = sin(qJ(5));
t114 = t122 ^ 2;
t124 = cos(qJ(5));
t164 = (-t124 ^ 2 + t114) * MDP(8);
t116 = sin(pkin(9));
t119 = cos(pkin(9));
t145 = qJDD(1) * t117;
t100 = qJDD(2) * t119 - t116 * t145;
t101 = qJDD(2) * t116 + t119 * t145;
t113 = pkin(9) + qJ(4);
t110 = sin(t113);
t111 = cos(t113);
t123 = sin(qJ(4));
t125 = cos(qJ(4));
t120 = cos(pkin(8));
t150 = t120 * t121;
t151 = t118 * t120;
t157 = g(3) * t117;
t163 = g(1) * (t110 * t118 + t111 * t150) + g(2) * (-t110 * t121 + t111 * t151) + t111 * t157 - t100 * t123 - t101 * t125;
t162 = -g(1) * (-t110 * t150 + t111 * t118) - g(2) * (-t110 * t151 - t111 * t121) + t110 * t157 + t100 * t125 - t101 * t123;
t104 = t116 * t123 - t119 * t125;
t98 = qJD(4) * t104;
t160 = qJD(4) * t117;
t156 = qJD(4) * pkin(4);
t155 = pkin(6) * qJDD(5);
t97 = t104 * t117;
t153 = qJD(5) * t97;
t149 = t124 * MDP(7);
t148 = qJDD(1) - g(3);
t144 = qJDD(5) * t122;
t143 = qJDD(5) * t124;
t142 = -g(1) * t118 + g(2) * t121;
t141 = -0.2e1 * t156;
t105 = t116 * t125 + t119 * t123;
t96 = t105 * t117;
t138 = -qJDD(4) * t96 + t98 * t160;
t137 = -qJDD(4) * t104 - t127 * t105;
t90 = t105 * t160;
t136 = qJD(4) * t96 + qJD(5) * t120 + t90;
t126 = qJD(5) ^ 2;
t133 = t105 * t126 - t137;
t132 = qJDD(5) * t120 - t138;
t131 = 0.2e1 * qJD(5) * t98 - qJDD(5) * t105;
t108 = -qJDD(1) * t120 + qJDD(3);
t130 = -g(3) * t120 - t108 + t166;
t129 = -(qJDD(4) * pkin(6)) + qJD(4) * t156 + t163;
t128 = pkin(4) * t165 - pkin(6) * t126 + t162;
t107 = -t122 * t126 + t143;
t106 = t124 * t126 + t144;
t1 = [t148 * MDP(1) + (-g(3) + (t117 ^ 2 + t120 ^ 2) * qJDD(1)) * MDP(2) + (-t108 * t120 - g(3) + (-t100 * t116 + t101 * t119) * t117) * MDP(3) + t138 * MDP(5) + (qJD(4) * t90 + qJDD(4) * t97) * MDP(6) + (t97 * t144 - t132 * t124 + (t136 * t122 + t124 * t153) * qJD(5)) * MDP(12) + (t97 * t143 + t132 * t122 + (-t122 * t153 + t136 * t124) * qJD(5)) * MDP(13); (qJDD(2) + t142) * MDP(2) + (t100 * t119 + t101 * t116 + t142) * MDP(3) + t137 * MDP(5) + (qJD(4) * t98 - qJDD(4) * t105) * MDP(6) + (-t133 * MDP(12) + t131 * MDP(13)) * t124 + (t131 * MDP(12) + t133 * MDP(13)) * t122; (-t148 * t120 + qJDD(3) - t166) * MDP(3) + t107 * MDP(12) - t106 * MDP(13); t162 * MDP(5) + t163 * MDP(6) + t106 * MDP(9) + t107 * MDP(10) + (t114 * MDP(7) + MDP(4)) * qJDD(4) - 0.2e1 * qJD(5) * qJD(4) * t164 + (t128 * MDP(12) + (t141 * qJD(5) - t155) * MDP(13)) * t124 + (t124 * MDP(8) * t165 - MDP(12) * t155 - t128 * MDP(13) + (t141 * MDP(12) + 0.2e1 * qJD(4) * t149) * qJD(5)) * t122; qJDD(5) * MDP(11) + t127 * t164 + (qJDD(4) * MDP(10) - t130 * MDP(12) + t129 * MDP(13)) * t124 + (t129 * MDP(12) + t130 * MDP(13) + qJDD(4) * MDP(9) - t127 * t149) * t122;];
tau = t1;
