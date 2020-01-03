% Calculate vector of inverse dynamics joint torques for
% S5PPRRP4
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
% pkin [7x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d3,d4,theta1]';
% MDP [14x1]
%   Minimal dynamic parameter vector (fixed base model)
%   see S5PPRRP4_convert_par2_MPV_fixb.m
% 
% Output:
% tau [5x1]
%   joint torques of inverse dynamics (contains inertial, gravitational Coriolis and centrifugal forces)

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-31 17:34
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tau = S5PPRRP4_invdynJ_fixb_mdp_slag_vp(qJ, qJD, qJDD, g, pkin, MDP)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(5,1),zeros(3,1),zeros(7,1),zeros(14,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5PPRRP4_invdynJ_fixb_mdp_slag_vp: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5PPRRP4_invdynJ_fixb_mdp_slag_vp: qJD has to be [5x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [5 1]), ...
  'S5PPRRP4_invdynJ_fixb_mdp_slag_vp: qJDD has to be [5x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S5PPRRP4_invdynJ_fixb_mdp_slag_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [7 1]), ...
  'S5PPRRP4_invdynJ_fixb_mdp_slag_vp: pkin has to be [7x1] (double)');
assert(isreal(MDP) && all(size(MDP) == [14 1]), ...
  'S5PPRRP4_invdynJ_fixb_mdp_slag_vp: MDP has to be [14x1] (double)'); 

%% Symbolic Calculation
% From invdyn_joint_fixb_mdp_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 17:34:44
% EndTime: 2019-12-31 17:34:45
% DurationCPUTime: 0.61s
% Computational Cost: add. (393->137), mult. (797->179), div. (0->0), fcn. (484->6), ass. (0->69)
t155 = qJ(5) + pkin(6);
t109 = sin(qJ(4));
t105 = t109 ^ 2;
t111 = cos(qJ(4));
t106 = t111 ^ 2;
t156 = (t105 - t106) * MDP(7);
t126 = (t105 + t106) * MDP(13);
t152 = qJD(4) * pkin(4);
t110 = sin(qJ(3));
t142 = qJD(2) * t110;
t127 = t155 * qJD(3) + t142;
t81 = -qJD(1) * t111 - t127 * t109;
t80 = t81 + t152;
t154 = t80 - t81;
t153 = qJD(3) * pkin(3);
t151 = t109 * t80;
t150 = cos(pkin(7));
t149 = sin(pkin(7));
t102 = pkin(4) * t111 + pkin(3);
t112 = cos(qJ(3));
t141 = qJD(2) * t112;
t85 = -t102 * qJD(3) + qJD(5) - t141;
t148 = qJD(3) * t85;
t146 = qJDD(1) - g(3);
t113 = qJD(4) ^ 2;
t114 = qJD(3) ^ 2;
t143 = t113 + t114;
t140 = qJD(3) * t111;
t139 = qJD(1) * qJD(4);
t138 = qJD(2) * qJD(3);
t137 = qJD(3) * qJD(4);
t136 = qJDD(2) * t110;
t135 = qJDD(3) * t109;
t134 = qJDD(3) * t111;
t133 = t112 * qJDD(2);
t132 = qJDD(4) * MDP(11);
t131 = qJDD(4) * MDP(12);
t130 = MDP(14) * t141;
t129 = t109 * t137;
t100 = t110 * t138;
t128 = qJD(4) * t155;
t125 = t100 - t133;
t88 = -t149 * t110 - t150 * t112;
t89 = t150 * t110 - t149 * t112;
t124 = g(1) * t89 - g(2) * t88;
t123 = g(1) * t88 + g(2) * t89;
t122 = -g(1) * t149 + g(2) * t150;
t103 = t109 * qJD(1);
t82 = t127 * t111 - t103;
t121 = -t111 * t82 + t151;
t120 = -MDP(5) + t126;
t95 = -t141 - t153;
t119 = t141 + t95 - t153;
t87 = qJDD(3) * pkin(6) + t112 * t138 + t136;
t118 = -qJ(5) * qJDD(3) - qJD(3) * qJD(5) - t87;
t117 = -t127 * qJD(4) - qJDD(1);
t116 = -qJD(3) * t95 - t123 - t87;
t79 = pkin(4) * t129 - t102 * qJDD(3) + qJDD(5) + t125;
t115 = 0.2e1 * qJDD(3) * pkin(3) - pkin(6) * t113 + t100 + t124 - t125;
t101 = t109 * t139;
t93 = t155 * t111;
t92 = t155 * t109;
t91 = qJDD(4) * t111 - t109 * t113;
t90 = qJDD(4) * t109 + t111 * t113;
t84 = -qJD(5) * t109 - t111 * t128;
t83 = qJD(5) * t111 - t109 * t128;
t78 = t117 * t109 + (-t118 - t139) * t111;
t77 = qJDD(4) * pkin(4) + t118 * t109 + t117 * t111 + t101;
t1 = [-t91 * MDP(11) + t90 * MDP(12) + (t121 * qJD(4) - t109 * t78 - t111 * t77 - g(3)) * MDP(14) + (MDP(1) + MDP(2)) * t146; (qJDD(2) + t122) * MDP(2) + t122 * MDP(14) + (qJDD(3) * MDP(4) + (-0.2e1 * t129 + t134) * MDP(11) + (-0.2e1 * t111 * t137 - t135) * MDP(12) + (-qJD(3) * t151 + t82 * t140 - t79) * MDP(14) + t120 * t114) * t112 + (MDP(14) * t148 - t114 * MDP(4) + (-t143 * MDP(11) - t131 + (-qJD(4) * t80 + t78) * MDP(14)) * t111 + (-t132 + t143 * MDP(12) + (-qJD(4) * t82 - t77) * MDP(14)) * t109 + t120 * qJDD(3)) * t110; (t124 + t133) * MDP(4) + (-t123 - t136) * MDP(5) + t90 * MDP(8) + t91 * MDP(9) + t123 * MDP(13) + (t78 * t93 + t82 * t83 - t77 * t92 + t80 * t84 - t79 * t102 - t85 * t142 - g(1) * (-t102 * t89 - t155 * t88) - g(2) * (t102 * t88 - t155 * t89)) * MDP(14) + (t105 * MDP(6) + MDP(3)) * qJDD(3) + (-0.2e1 * qJD(4) * t156 - t126 * t141) * qJD(3) + (t115 * MDP(11) - pkin(6) * t131 + (qJD(3) * t83 + qJDD(3) * t93 + t78) * MDP(13) - t82 * t130 + (t119 * MDP(12) + (qJD(3) * t92 - t80) * MDP(13)) * qJD(4)) * t111 + (0.2e1 * MDP(7) * t134 - pkin(6) * t132 - t115 * MDP(12) + (-qJD(3) * t84 + qJDD(3) * t92 - t77) * MDP(13) + t80 * t130 + (0.2e1 * MDP(6) * t140 + t119 * MDP(11) + (-qJD(3) * t93 - t82) * MDP(13) + pkin(4) * t85 * MDP(14)) * qJD(4)) * t109; MDP(8) * t135 + MDP(9) * t134 + qJDD(4) * MDP(10) + (-t103 * qJD(4) + t116 * t109 - t146 * t111 + t101) * MDP(11) + (t146 * t109 + t116 * t111) * MDP(12) + (-pkin(4) * t135 + (-t152 + t154) * t140) * MDP(13) + (t154 * t82 + (g(3) * t111 + t77 + (-t123 - t148) * t109) * pkin(4)) * MDP(14) + (-t109 * t111 * MDP(6) + t156) * t114; (t121 * qJD(3) - t124 + t79) * MDP(14) - t114 * t126;];
tau = t1;
