% Calculate vector of inverse dynamics joint torques for
% S4RRRP1
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
%   pkin=[a2,a3,a4,d1,d2,d3]';
% MDP [10x1]
%   Minimal dynamic parameter vector (fixed base model)
%   see S4RRRP1_convert_par2_MPV_fixb.m
% 
% Output:
% tau [4x1]
%   joint torques of inverse dynamics (contains inertial, gravitational Coriolis and centrifugal forces)

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-08 18:36
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tau = S4RRRP1_invdynJ_fixb_mdp_slag_vp(qJ, qJD, qJDD, g, pkin, MDP)
%% Coder Information
%#codegen
%$cgargs {zeros(4,1),zeros(4,1),zeros(4,1),zeros(3,1),zeros(6,1),zeros(10,1)}
assert(isreal(qJ) && all(size(qJ) == [4 1]), ...
  'S4RRRP1_invdynJ_fixb_mdp_slag_vp: qJ has to be [4x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [4 1]), ...
  'S4RRRP1_invdynJ_fixb_mdp_slag_vp: qJD has to be [4x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [4 1]), ...
  'S4RRRP1_invdynJ_fixb_mdp_slag_vp: qJDD has to be [4x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S4RRRP1_invdynJ_fixb_mdp_slag_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [6 1]), ...
  'S4RRRP1_invdynJ_fixb_mdp_slag_vp: pkin has to be [6x1] (double)');
assert(isreal(MDP) && all(size(MDP) == [10 1]), ...
  'S4RRRP1_invdynJ_fixb_mdp_slag_vp: MDP has to be [10x1] (double)'); 

%% Symbolic Calculation
% From invdyn_joint_fixb_mdp_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-08 18:36:05
% EndTime: 2019-03-08 18:36:06
% DurationCPUTime: 0.43s
% Computational Cost: add. (357->114), mult. (569->148), div. (0->0), fcn. (291->10), ass. (0->66)
t119 = cos(qJ(3));
t120 = cos(qJ(2));
t142 = qJD(2) * t120;
t135 = qJD(1) * t142;
t117 = sin(qJ(2));
t138 = qJDD(1) * t117;
t124 = t135 + t138;
t114 = qJD(1) + qJD(2);
t150 = pkin(1) * qJD(1);
t91 = t114 * pkin(2) + t120 * t150;
t155 = (t124 * pkin(1) + qJD(3) * t91) * t119;
t115 = qJ(1) + qJ(2);
t112 = qJ(3) + t115;
t104 = cos(t112);
t103 = sin(t112);
t96 = g(1) * t103;
t154 = -g(2) * t104 + t96;
t110 = sin(t115);
t100 = g(1) * t110;
t151 = t120 * pkin(1);
t111 = cos(t115);
t149 = g(1) * t111 + g(2) * t110;
t113 = qJDD(1) + qJDD(2);
t108 = qJDD(3) + t113;
t94 = t108 * MDP(7);
t148 = t113 * MDP(4) + t94;
t147 = pkin(2) * t111 + pkin(3) * t104;
t116 = sin(qJ(3));
t146 = MDP(8) * t116;
t144 = t116 * t117;
t143 = t117 * t119;
t141 = qJD(3) * t116;
t140 = qJD(3) * t117;
t139 = qJD(3) * t119;
t137 = t117 * t150;
t136 = t117 * t139;
t134 = qJD(1) * t140;
t131 = t116 * t137;
t106 = qJDD(1) * t151;
t86 = t113 * pkin(2) - qJD(2) * t137 + t106;
t133 = -qJD(3) * t131 + t116 * t86;
t132 = qJD(2) * (-qJD(1) - t114);
t130 = -g(2) * t111 + t100 + t106;
t129 = t119 * t86 - t91 * t141;
t128 = g(1) * t104 + g(2) * t103 - t133;
t127 = -t116 * t120 - t143;
t126 = t119 * t120 - t144;
t83 = t119 * t91 - t131;
t125 = -MDP(9) * t119 - MDP(6) - t146;
t123 = t129 + t154;
t122 = -t91 * t139 + t128;
t121 = cos(qJ(1));
t118 = sin(qJ(1));
t109 = qJD(3) + t114;
t105 = pkin(2) + t151;
t93 = t119 * t105;
t89 = pkin(1) * t143 + t116 * t105;
t88 = t126 * t150;
t87 = t127 * t150;
t84 = t116 * t91 + t119 * t137;
t82 = t109 * pkin(3) + t83;
t81 = -t105 * t141 + (t127 * qJD(2) - t136) * pkin(1);
t80 = t105 * t139 + (t126 * qJD(2) - t116 * t140) * pkin(1);
t79 = t133 + t155;
t78 = t108 * pkin(3) + (-t116 * t138 + (-t116 * t142 - t136) * qJD(1)) * pkin(1) + t129;
t1 = [qJDD(1) * MDP(1) + (g(1) * t118 - g(2) * t121) * MDP(2) + (g(1) * t121 + g(2) * t118) * MDP(3) + ((t113 * t120 + t117 * t132) * pkin(1) + t130) * MDP(5) + (((-qJDD(1) - t113) * t117 + t120 * t132) * pkin(1) + t149) * MDP(6) + (t93 * t108 + t81 * t109 + (-t119 * t134 + (-t135 + (-qJDD(1) - t108) * t117) * t116) * pkin(1) + t123) * MDP(8) + (-t89 * t108 - t80 * t109 + t128 - t155) * MDP(9) + (t79 * t89 + t84 * t80 + t78 * (-pkin(1) * t144 + pkin(3) + t93) + t82 * t81 - g(1) * (-t118 * pkin(1) - pkin(2) * t110 - pkin(3) * t103) - g(2) * (t121 * pkin(1) + t147)) * MDP(10) + t148; t130 * MDP(5) + t149 * MDP(6) + (-t87 * t109 + t123) * MDP(8) + (t88 * t109 + t122) * MDP(9) + (-g(2) * t147 - t82 * t87 - t84 * t88 + (t78 + t96) * pkin(3)) * MDP(10) + ((t108 * t119 - t109 * t141) * MDP(8) + (-t108 * t116 - t109 * t139) * MDP(9) + (t79 * t116 + t78 * t119 + t84 * t139 - t82 * t141 + t100) * MDP(10)) * pkin(2) + (t125 * t138 + (((-qJD(2) + t114) * MDP(5) - MDP(8) * t139) * t117 + (t114 * MDP(6) + t125 * qJD(2)) * t120) * qJD(1)) * pkin(1) + t148; t94 + t123 * MDP(8) + (t83 * t109 + t122) * MDP(9) + (t109 * MDP(8) + (t82 - t83) * MDP(10)) * t84 + (t78 + t154) * MDP(10) * pkin(3) + (-t124 * t146 + (-MDP(8) * t134 - t124 * MDP(9)) * t119) * pkin(1); (qJDD(4) - g(3)) * MDP(10);];
tau  = t1;
