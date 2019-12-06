% Calculate vector of inverse dynamics joint torques for
% S5PPPRR2
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
%   see S5PPPRR2_convert_par2_MPV_fixb.m
% 
% Output:
% tau [5x1]
%   joint torques of inverse dynamics (contains inertial, gravitational Coriolis and centrifugal forces)

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-05 14:59
% Revision: 77da58f92bca3eff71542919beafa37024070d86 (2019-12-05)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tau = S5PPPRR2_invdynJ_fixb_mdp_slag_vp(qJ, qJD, qJDD, g, pkin, MDP)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(5,1),zeros(3,1),zeros(9,1),zeros(13,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5PPPRR2_invdynJ_fixb_mdp_slag_vp: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5PPPRR2_invdynJ_fixb_mdp_slag_vp: qJD has to be [5x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [5 1]), ...
  'S5PPPRR2_invdynJ_fixb_mdp_slag_vp: qJDD has to be [5x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S5PPPRR2_invdynJ_fixb_mdp_slag_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S5PPPRR2_invdynJ_fixb_mdp_slag_vp: pkin has to be [9x1] (double)');
assert(isreal(MDP) && all(size(MDP) == [13 1]), ...
  'S5PPPRR2_invdynJ_fixb_mdp_slag_vp: MDP has to be [13x1] (double)'); 

%% Symbolic Calculation
% From invdyn_joint_fixb_mdp_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-05 14:59:40
% EndTime: 2019-12-05 14:59:42
% DurationCPUTime: 0.62s
% Computational Cost: add. (301->115), mult. (690->190), div. (0->0), fcn. (613->10), ass. (0->53)
t140 = -2 * qJD(4) * qJD(5);
t125 = qJD(4) ^ 2;
t120 = sin(qJ(5));
t112 = t120 ^ 2;
t122 = cos(qJ(5));
t163 = (-t122 ^ 2 + t112) * MDP(8);
t114 = sin(pkin(9));
t116 = sin(pkin(7));
t118 = cos(pkin(8));
t117 = cos(pkin(9));
t119 = cos(pkin(7));
t151 = t119 * t117;
t101 = t114 * t116 + t118 * t151;
t123 = cos(qJ(4));
t115 = sin(pkin(8));
t121 = sin(qJ(4));
t155 = t115 * t121;
t102 = t117 * t155 + t118 * t123;
t144 = qJDD(1) * t115;
t105 = qJDD(2) * t114 + t117 * t144;
t108 = -qJDD(1) * t118 + qJDD(3);
t154 = t115 * t123;
t152 = t119 * t114;
t153 = t116 * t118;
t99 = t117 * t153 - t152;
t162 = -g(1) * (-t101 * t121 + t119 * t154) - g(2) * (t116 * t154 - t121 * t99) + g(3) * t102 - t121 * t105 + t123 * t108;
t161 = -t123 * t105 - t121 * t108 + g(1) * (t101 * t123 + t119 * t155) + g(2) * (t116 * t155 + t123 * t99);
t160 = 2 * qJD(4);
t159 = qJD(4) * pkin(4);
t157 = t114 * t115;
t156 = t114 * t123;
t150 = t122 * t125;
t148 = qJDD(1) - g(3);
t103 = t117 * t154 - t118 * t121;
t146 = qJD(5) * t103;
t145 = t102 * qJD(4);
t143 = t117 * qJDD(5);
t142 = t122 * qJDD(4);
t141 = t123 * qJDD(5);
t139 = -g(1) * t116 + g(2) * t119;
t138 = t114 * t121 * t160;
t137 = -2 * t159;
t132 = qJDD(4) * t102 + t103 * t125;
t131 = qJDD(4) * t121 + t123 * t125;
t104 = -qJDD(2) * t117 + t114 * t144;
t130 = g(1) * (-t116 * t117 + t118 * t152) + g(2) * (t114 * t153 + t151) - t104;
t128 = -qJD(5) * t157 + 0.2e1 * t145;
t127 = -qJDD(4) * pkin(6) + (qJD(4) * t159) + t161;
t124 = qJD(5) ^ 2;
t126 = 0.2e1 * qJDD(4) * pkin(4) - pkin(6) * t124 + t162;
t93 = t103 * t122 + t120 * t157;
t92 = -t103 * t120 + t122 * t157;
t1 = [t148 * MDP(1) + (-g(3) + (t115 ^ 2 + t118 ^ 2) * qJDD(1)) * MDP(2) + (-t108 * t118 - g(3) + (t104 * t114 + t105 * t117) * t115) * MDP(3) - t132 * MDP(5) + (qJD(4) * t145 - qJDD(4) * t103) * MDP(6) + (t92 * qJDD(5) - t132 * t122 + (t120 * t128 - t122 * t146) * qJD(5)) * MDP(12) + (-t93 * qJDD(5) + t132 * t120 + (t120 * t146 + t122 * t128) * qJD(5)) * MDP(13); (qJDD(2) + t139) * MDP(2) + (-t104 * t117 + t139) * MDP(3) + (-t122 * t143 + (t120 * t138 + (t117 * t120 - t122 * t156) * qJD(5)) * qJD(5)) * MDP(12) + (t120 * t143 + (t122 * t138 + (t117 * t122 + t120 * t156) * qJD(5)) * qJD(5)) * MDP(13) + (t105 * MDP(3) - t131 * MDP(5) + (-t123 * qJDD(4) + t121 * t125) * MDP(6) + (-t120 * t141 - t122 * t131) * MDP(12) + (t120 * t131 - t122 * t141) * MDP(13)) * t114; (qJDD(3) - t148 * t118 + (-g(1) * t119 - g(2) * t116) * t115) * MDP(3) + (qJDD(4) * MDP(5) - t125 * MDP(6) + (t120 * t140 + t142) * MDP(12) + (-qJDD(4) * t120 + t122 * t140) * MDP(13)) * t123 + (-t125 * MDP(5) - qJDD(4) * MDP(6) + (-qJDD(5) * t120 - t122 * t124 - t150) * MDP(12) + (-qJDD(5) * t122 + (t124 + t125) * t120) * MDP(13)) * t121; t162 * MDP(5) + (g(3) * t103 + t161) * MDP(6) + (MDP(7) * t112 + MDP(4)) * qJDD(4) + t163 * t140 + ((t124 * MDP(9)) + t126 * MDP(12) + (-MDP(13) * pkin(6) + MDP(10)) * qJDD(5) + t137 * MDP(13) * qJD(5)) * t122 + (0.2e1 * MDP(8) * t142 - t124 * MDP(10) - t126 * MDP(13) + (-MDP(12) * pkin(6) + MDP(9)) * qJDD(5) + (MDP(7) * t122 * t160 + MDP(12) * t137) * qJD(5)) * t120; qJDD(5) * MDP(11) + t125 * t163 + (-MDP(12) * t92 + MDP(13) * t93) * g(3) + (qJDD(4) * MDP(10) - MDP(12) * t130 + MDP(13) * t127) * t122 + (MDP(12) * t127 + MDP(13) * t130 - MDP(7) * t150 + qJDD(4) * MDP(9)) * t120;];
tau = t1;
