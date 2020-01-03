% Calculate vector of inverse dynamics joint torques for
% S4RRPR5
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
%   pkin=[a2,a3,a4,d1,d2,d4]';
% MDP [16x1]
%   Minimal dynamic parameter vector (fixed base model)
%   see S4RRPR5_convert_par2_MPV_fixb.m
% 
% Output:
% tau [4x1]
%   joint torques of inverse dynamics (contains inertial, gravitational Coriolis and centrifugal forces)

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-31 17:03
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tau = S4RRPR5_invdynJ_fixb_mdp_slag_vp(qJ, qJD, qJDD, g, pkin, MDP)
%% Coder Information
%#codegen
%$cgargs {zeros(4,1),zeros(4,1),zeros(4,1),zeros(3,1),zeros(6,1),zeros(16,1)}
assert(isreal(qJ) && all(size(qJ) == [4 1]), ...
  'S4RRPR5_invdynJ_fixb_mdp_slag_vp: qJ has to be [4x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [4 1]), ...
  'S4RRPR5_invdynJ_fixb_mdp_slag_vp: qJD has to be [4x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [4 1]), ...
  'S4RRPR5_invdynJ_fixb_mdp_slag_vp: qJDD has to be [4x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S4RRPR5_invdynJ_fixb_mdp_slag_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [6 1]), ...
  'S4RRPR5_invdynJ_fixb_mdp_slag_vp: pkin has to be [6x1] (double)');
assert(isreal(MDP) && all(size(MDP) == [16 1]), ...
  'S4RRPR5_invdynJ_fixb_mdp_slag_vp: MDP has to be [16x1] (double)'); 

%% Symbolic Calculation
% From invdyn_joint_fixb_mdp_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 17:03:35
% EndTime: 2019-12-31 17:03:36
% DurationCPUTime: 0.56s
% Computational Cost: add. (361->127), mult. (482->153), div. (0->0), fcn. (214->8), ass. (0->63)
t116 = qJD(1) + qJD(2);
t119 = qJ(1) + qJ(2);
t112 = sin(t119);
t113 = cos(t119);
t154 = g(1) * t113 + g(2) * t112;
t124 = cos(qJ(2));
t150 = qJD(2) * t124;
t97 = pkin(1) * t150 + qJD(3);
t169 = t116 * t97 - t154;
t121 = sin(qJ(2));
t144 = qJDD(1) * t121;
t151 = qJD(1) * t124;
t168 = -((qJD(2) - t116) * t151 + t144) * pkin(1) + t154;
t167 = -g(1) * t112 + g(2) * t113;
t160 = pkin(1) * qJD(1);
t141 = t121 * t160;
t155 = qJ(3) * t116;
t96 = t141 + t155;
t166 = -t116 * t96 + t167;
t114 = t116 ^ 2;
t126 = -pkin(2) - pkin(6);
t164 = pkin(1) * t121;
t163 = pkin(1) * t124;
t115 = qJDD(1) + qJDD(2);
t162 = pkin(2) * t115;
t120 = sin(qJ(4));
t123 = cos(qJ(4));
t147 = qJD(4) * t123;
t149 = qJD(3) * t116;
t156 = qJ(3) * t115;
t89 = t156 + t149 + (qJD(1) * t150 + t144) * pkin(1);
t161 = t89 * t120 + t96 * t147;
t140 = qJD(2) * t164;
t157 = qJD(1) * t140 - qJDD(1) * t163;
t127 = qJD(4) ^ 2;
t153 = -t114 - t127;
t118 = t123 ^ 2;
t152 = t120 ^ 2 - t118;
t148 = qJD(4) * t116;
t108 = -pkin(2) - t163;
t99 = -pkin(6) + t108;
t146 = qJDD(4) * t99;
t143 = qJDD(4) * t120;
t142 = qJDD(4) * t126;
t139 = qJDD(3) + t157;
t138 = t157 + t167;
t111 = qJDD(4) * t123;
t137 = t115 * MDP(4) + 0.2e1 * (-t115 * t120 * t123 + t152 * t148) * MDP(11) + (-0.2e1 * t116 * t120 * t147 + t115 * t118) * MDP(10) + (-t123 * t127 - t143) * MDP(13) + (-t120 * t127 + t111) * MDP(12);
t136 = -pkin(1) * t151 + qJD(3);
t122 = sin(qJ(1));
t125 = cos(qJ(1));
t134 = g(1) * t122 - g(2) * t125;
t91 = t139 - t162;
t133 = -t141 + t155;
t132 = -t126 * t115 - t139 - t166;
t131 = t116 * t141 - t138;
t130 = -g(1) * (-pkin(2) * t112 + t113 * qJ(3)) - g(2) * (t113 * pkin(2) + t112 * qJ(3));
t102 = qJ(3) + t164;
t129 = t102 * t115 - t127 * t99 + t169;
t128 = t136 * t116 - t126 * t127 - t154 + t156;
t93 = -pkin(2) * t116 + t136;
t87 = t89 * t123;
t1 = [t161 * MDP(15) + t87 * MDP(16) + (qJDD(3) + t138) * MDP(7) + (g(1) * t125 + g(2) * t122) * MDP(3) + (t149 + t169) * MDP(8) + ((-pkin(2) + t108) * MDP(7) + (qJ(3) + t102) * MDP(8)) * t115 - t138 * MDP(5) + t154 * MDP(6) + ((qJD(4) * t140 + t102 * t148 + t146) * MDP(15) + t129 * MDP(16)) * t123 + (t129 * MDP(15) + (-t146 + (-t102 * t116 - t140 - t96) * qJD(4)) * MDP(16)) * t120 + (t89 * t102 + t91 * t108 + t96 * t97 + t130) * MDP(9) + t137 + qJDD(1) * MDP(1) + t134 * MDP(2) + (t115 * t124 * MDP(5) + t134 * MDP(9) + ((-qJDD(1) - t115) * MDP(6) + qJDD(1) * MDP(8)) * t121 + (((-qJD(1) - t116) * MDP(6) + qJD(1) * MDP(8)) * t124 + (t93 * MDP(9) + (-MDP(5) + MDP(7)) * t116) * t121) * qJD(2)) * pkin(1); t131 * MDP(5) + t168 * MDP(6) + (qJDD(3) - t131 - 0.2e1 * t162) * MDP(7) + (0.2e1 * t149 + 0.2e1 * t156 - t168) * MDP(8) + (t89 * qJ(3) + t96 * qJD(3) - t91 * pkin(2) + (-t121 * t93 - t124 * t96) * t160 + t130) * MDP(9) + ((t133 * qJD(4) + t142) * t123 + t128 * t120 + t161) * MDP(15) + (t87 + (-t142 + (-t133 - t96) * qJD(4)) * t120 + t128 * t123) * MDP(16) + t137; t115 * MDP(7) - t114 * MDP(8) + (t91 + t166) * MDP(9) + (t153 * t120 + t111) * MDP(15) + (t153 * t123 - t143) * MDP(16); qJDD(4) * MDP(14) - t152 * MDP(11) * t114 + (t115 * MDP(12) - t132 * MDP(15) + g(3) * MDP(16)) * t123 + (t123 * t114 * MDP(10) - t115 * MDP(13) + g(3) * MDP(15) + t132 * MDP(16)) * t120;];
tau = t1;
