% Calculate vector of inverse dynamics joint torques for
% S4RPPR3
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
%   pkin=[a2,a3,a4,d1,d4,theta2,theta3]';
% MDP [15x1]
%   Minimal dynamic parameter vector (fixed base model)
%   see S4RPPR3_convert_par2_MPV_fixb.m
% 
% Output:
% tau [4x1]
%   joint torques of inverse dynamics (contains inertial, gravitational Coriolis and centrifugal forces)

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-31 16:38
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tau = S4RPPR3_invdynJ_fixb_mdp_slag_vp(qJ, qJD, qJDD, g, pkin, MDP)
%% Coder Information
%#codegen
%$cgargs {zeros(4,1),zeros(4,1),zeros(4,1),zeros(3,1),zeros(7,1),zeros(15,1)}
assert(isreal(qJ) && all(size(qJ) == [4 1]), ...
  'S4RPPR3_invdynJ_fixb_mdp_slag_vp: qJ has to be [4x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [4 1]), ...
  'S4RPPR3_invdynJ_fixb_mdp_slag_vp: qJD has to be [4x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [4 1]), ...
  'S4RPPR3_invdynJ_fixb_mdp_slag_vp: qJDD has to be [4x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S4RPPR3_invdynJ_fixb_mdp_slag_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [7 1]), ...
  'S4RPPR3_invdynJ_fixb_mdp_slag_vp: pkin has to be [7x1] (double)');
assert(isreal(MDP) && all(size(MDP) == [15 1]), ...
  'S4RPPR3_invdynJ_fixb_mdp_slag_vp: MDP has to be [15x1] (double)'); 

%% Symbolic Calculation
% From invdyn_joint_fixb_mdp_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 16:37:58
% EndTime: 2019-12-31 16:38:00
% DurationCPUTime: 0.70s
% Computational Cost: add. (322->114), mult. (620->157), div. (0->0), fcn. (423->12), ass. (0->61)
t125 = sin(pkin(7));
t127 = cos(pkin(7));
t151 = t125 ^ 2 + t127 ^ 2;
t161 = qJD(1) * t151;
t126 = sin(pkin(6));
t110 = pkin(1) * t126 + qJ(3);
t160 = t110 * t161;
t124 = qJ(1) + pkin(6);
t118 = sin(t124);
t120 = cos(t124);
t159 = g(1) * t120 + g(2) * t118;
t99 = qJD(1) * qJD(3) + qJDD(1) * t110;
t156 = pkin(5) + t110;
t87 = t125 * qJDD(2) + t127 * t99;
t153 = pkin(5) * qJDD(1);
t152 = t127 * MDP(5);
t129 = sin(qJ(4));
t150 = qJD(1) * t129;
t131 = cos(qJ(4));
t149 = qJD(1) * t131;
t128 = cos(pkin(6));
t112 = -pkin(1) * t128 - pkin(2);
t147 = qJDD(1) * t112;
t146 = qJDD(1) * t129;
t145 = qJDD(1) * t131;
t144 = t125 * t145 + (qJD(4) * t149 + t146) * t127;
t143 = t125 * t150;
t142 = g(1) * t118 - g(2) * t120;
t130 = sin(qJ(1));
t132 = cos(qJ(1));
t141 = g(1) * t130 - g(2) * t132;
t109 = t127 * t145;
t140 = -t125 * t146 + t109;
t114 = t127 * qJDD(2);
t86 = -t125 * t99 + t114;
t139 = -t125 * t86 + t127 * t87;
t97 = t156 * t125;
t98 = t156 * t127;
t138 = t129 * t98 + t131 * t97;
t137 = -t129 * t97 + t131 * t98;
t102 = t125 * t131 + t127 * t129;
t101 = t125 * t129 - t131 * t127;
t104 = -pkin(3) * t127 + t112;
t96 = t102 * qJD(4);
t136 = -t125 * t149 - t127 * t150;
t123 = pkin(7) + qJ(4);
t119 = cos(t123);
t117 = sin(t123);
t103 = qJDD(3) + t147;
t95 = t101 * qJD(4);
t94 = t102 * qJD(1);
t93 = t101 * qJD(1);
t92 = t104 * qJD(1) + qJD(3);
t90 = t104 * qJDD(1) + qJDD(3);
t83 = t127 * t153 + t87;
t82 = t114 + (-t99 - t153) * t125;
t81 = -qJD(4) * t96 - qJDD(4) * t101;
t80 = -qJD(4) * t95 + qJDD(4) * t102;
t79 = qJD(1) * t96 - t140;
t78 = -qJD(4) * t143 + t144;
t1 = [qJDD(1) * MDP(1) + t141 * MDP(2) + (g(1) * t132 + g(2) * t130) * MDP(3) + (t141 + (t126 ^ 2 + t128 ^ 2) * qJDD(1) * pkin(1)) * pkin(1) * MDP(4) + (t99 * t151 + t139 - t159) * MDP(7) + (t103 * t112 - g(1) * (-pkin(1) * t130 - pkin(2) * t118 + qJ(3) * t120) - g(2) * (pkin(1) * t132 + pkin(2) * t120 + qJ(3) * t118) + t139 * t110 + t160 * qJD(3)) * MDP(8) + (t102 * t78 - t94 * t95) * MDP(9) + (-t101 * t78 - t102 * t79 + t93 * t95 - t94 * t96) * MDP(10) + t80 * MDP(11) + t81 * MDP(12) + (t104 * t79 + t90 * t101 + t92 * t96 - t138 * qJDD(4) + t142 * t119 + (-t102 * qJD(3) - t137 * qJD(4)) * qJD(4)) * MDP(14) + (t104 * t78 + t90 * t102 - t92 * t95 - t137 * qJDD(4) - t142 * t117 + (t101 * qJD(3) + t138 * qJD(4)) * qJD(4)) * MDP(15) + (-t125 * MDP(6) + t152) * (-t103 + t142 - t147); (qJDD(2) - g(3)) * MDP(4) + (t125 * t87 + t127 * t86 - g(3)) * MDP(8) + t81 * MDP(14) - t80 * MDP(15); (qJDD(3) - t142) * MDP(8) - t109 * MDP(14) + t144 * MDP(15) + (-t152 + t112 * MDP(8) + (t129 * MDP(14) + MDP(6)) * t125) * qJDD(1) + ((-t136 + t94) * MDP(14) + (-t93 - t143) * MDP(15)) * qJD(4) + (-MDP(7) * t161 - t160 * MDP(8)) * qJD(1); t94 * t93 * MDP(9) + (-t93 ^ 2 + t94 ^ 2) * MDP(10) + t144 * MDP(11) + t140 * MDP(12) + qJDD(4) * MDP(13) + (-g(3) * t119 + t159 * t117 - t129 * t83 + t131 * t82 - t92 * t94) * MDP(14) + (g(3) * t117 + t159 * t119 - t129 * t82 - t131 * t83 + t92 * t93) * MDP(15) + ((t93 - t143) * MDP(11) + (t136 + t94) * MDP(12)) * qJD(4);];
tau = t1;
