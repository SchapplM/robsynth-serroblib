% Calculate vector of inverse dynamics joint torques for
% S4PRRP3
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
%   see S4PRRP3_convert_par2_MPV_fixb.m
% 
% Output:
% tau [4x1]
%   joint torques of inverse dynamics (contains inertial, gravitational Coriolis and centrifugal forces)

% Quelle: HybrDyn-Toolbox
% Datum: 2021-01-14 22:27
% Revision: beb2ba9bd8c5bd556f42a244985f3dab86917626 (2021-01-14)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tau = S4PRRP3_invdynJ_fixb_mdp_slag_vp(qJ, qJD, qJDD, g, pkin, MDP)
%% Coder Information
%#codegen
%$cgargs {zeros(4,1),zeros(4,1),zeros(4,1),zeros(3,1),zeros(6,1),zeros(15,1)}
assert(isreal(qJ) && all(size(qJ) == [4 1]), ...
  'S4PRRP3_invdynJ_fixb_mdp_slag_vp: qJ has to be [4x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [4 1]), ...
  'S4PRRP3_invdynJ_fixb_mdp_slag_vp: qJD has to be [4x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [4 1]), ...
  'S4PRRP3_invdynJ_fixb_mdp_slag_vp: qJDD has to be [4x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S4PRRP3_invdynJ_fixb_mdp_slag_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [6 1]), ...
  'S4PRRP3_invdynJ_fixb_mdp_slag_vp: pkin has to be [6x1] (double)');
assert(isreal(MDP) && all(size(MDP) == [15 1]), ...
  'S4PRRP3_invdynJ_fixb_mdp_slag_vp: MDP has to be [15x1] (double)'); 

%% Symbolic Calculation
% From invdyn_joint_fixb_mdp_matlab.m
% OptimizationMode: 2
% StartTime: 2021-01-14 22:27:07
% EndTime: 2021-01-14 22:27:09
% DurationCPUTime: 0.63s
% Computational Cost: add. (325->139), mult. (637->176), div. (0->0), fcn. (318->4), ass. (0->70)
t154 = qJDD(1) - g(3);
t108 = sin(qJ(3));
t105 = t108 ^ 2;
t109 = cos(qJ(3));
t106 = t109 ^ 2;
t153 = (t105 - t106) * MDP(6);
t130 = qJD(1) * qJD(3);
t142 = qJ(4) + pkin(5);
t113 = -t142 * qJDD(2) - t130;
t152 = -qJD(2) * qJD(4) + t113;
t151 = -2 * pkin(2);
t150 = 0.2e1 * t109;
t104 = pkin(6) + qJ(2);
t99 = sin(t104);
t149 = g(1) * t99;
t148 = g(2) * t99;
t147 = pkin(3) * t105;
t146 = pkin(3) * t109;
t100 = cos(t104);
t145 = g(1) * t100;
t144 = g(2) * t100;
t143 = g(3) * t109;
t140 = qJD(3) * pkin(3);
t102 = t109 * qJD(1);
t124 = qJD(2) * t142;
t81 = -t108 * t124 + t102;
t79 = t81 + t140;
t141 = t79 - t81;
t139 = t109 * t99;
t138 = pkin(5) * qJDD(3);
t98 = pkin(2) + t146;
t137 = qJD(2) * t98;
t136 = qJDD(3) * pkin(3);
t135 = t100 * t108;
t111 = qJD(2) ^ 2;
t134 = t109 * t111;
t132 = qJD(2) * t108;
t131 = qJDD(2) * t98;
t129 = qJD(2) * qJD(3);
t127 = qJDD(2) * t108;
t126 = qJDD(2) * t109;
t89 = t142 * t109;
t125 = t108 * t129;
t123 = qJD(3) * t142;
t80 = pkin(3) * t125 + qJDD(4) - t131;
t122 = t80 - t131;
t121 = t109 * t123;
t120 = -t144 + t149;
t119 = t145 + t148;
t110 = qJD(3) ^ 2;
t118 = pkin(5) * t110 + qJDD(2) * t151;
t117 = pkin(2) * t111 - pkin(5) * qJDD(2);
t83 = t108 * qJD(1) + qJD(2) * t89;
t116 = t108 * t79 - t109 * t83;
t101 = t109 * qJDD(1);
t115 = g(1) * t135 + t108 * t148 + t101 - t143;
t96 = pkin(5) * t125;
t114 = g(2) * t139 - t154 * t108 + t109 * t145 + t96;
t85 = qJD(4) - t137;
t112 = (-qJD(4) - t85) * qJD(2) + t113;
t93 = g(1) * t139;
t92 = g(2) * t135;
t88 = t142 * t108;
t87 = qJDD(3) * t109 - t108 * t110;
t86 = qJDD(3) * t108 + t109 * t110;
t84 = -t108 * qJD(4) - t121;
t82 = t109 * qJD(4) - t108 * t123;
t78 = -t96 + (-qJ(4) * t129 + qJDD(1)) * t108 - t152 * t109;
t77 = -qJD(2) * t121 + t152 * t108 + t101 + t136;
t1 = [t154 * MDP(1) + (-t116 * qJD(3) + t78 * t108 + t77 * t109 - g(3)) * MDP(15) + (MDP(10) + MDP(12)) * t87 + (-MDP(11) - MDP(13)) * t86; t120 * MDP(3) + t86 * MDP(7) + t87 * MDP(8) + t93 * MDP(10) + t92 * MDP(11) + (-qJDD(3) * t88 + t93) * MDP(12) + (-qJDD(3) * t89 + t92) * MDP(13) + (t78 * t89 + t83 * t82 - t77 * t88 + t79 * t84 - t80 * t98 - g(1) * (t100 * t142 - t98 * t99) - g(2) * (t100 * t98 + t142 * t99)) * MDP(15) + (t105 * MDP(5) + MDP(2)) * qJDD(2) + (t84 * MDP(12) - t82 * MDP(13) + (MDP(13) * t147 - 0.2e1 * t153) * qJD(2)) * qJD(3) + ((-t118 - t144) * MDP(10) - MDP(11) * t138 + (-t122 - t144) * MDP(12) + (qJD(2) * t82 + qJDD(2) * t89 + t78) * MDP(14) + (qJD(2) * MDP(11) * t151 + (t85 - t137) * MDP(13) + (qJD(2) * t88 - t79) * MDP(14)) * qJD(3)) * t109 + (0.2e1 * MDP(6) * t126 - MDP(10) * t138 + (t118 - t149) * MDP(11) + (t122 - t149) * MDP(13) + (-qJD(2) * t84 + qJDD(2) * t88 - t77) * MDP(14) + (-t83 * MDP(14) + (MDP(15) * pkin(3) + MDP(12)) * t85 + (MDP(5) * t150 + MDP(10) * t151 + (-t98 - t146) * MDP(12) - t89 * MDP(14)) * qJD(2)) * qJD(3)) * t108 + (MDP(4) - MDP(14)) * t119; -t108 * MDP(5) * t134 + t111 * t153 + MDP(7) * t127 + MDP(8) * t126 + qJDD(3) * MDP(9) + (t117 * t108 + t115) * MDP(10) + ((-pkin(5) * t132 + t102) * qJD(3) + (t117 - t130) * t109 + t114) * MDP(11) + (0.2e1 * t136 + (-t109 * t124 + t83) * qJD(3) + (pkin(3) * t134 + t112) * t108 + t115) * MDP(12) + (-t111 * t147 + (qJ(4) * t132 + t81) * qJD(3) + t112 * t109 + t114) * MDP(13) + (-pkin(3) * t127 + (-t140 + t141) * t109 * qJD(2)) * MDP(14) + (t141 * t83 + (-t143 + t77 + (-qJD(2) * t85 + t119) * t108) * pkin(3)) * MDP(15); (0.2e1 * t125 - t126) * MDP(12) + (t129 * t150 + t127) * MDP(13) + (t116 * qJD(2) - t120 + t80) * MDP(15) + (-t105 - t106) * MDP(14) * t111;];
tau = t1;
