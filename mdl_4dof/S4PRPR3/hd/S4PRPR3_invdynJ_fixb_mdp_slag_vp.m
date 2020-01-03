% Calculate vector of inverse dynamics joint torques for
% S4PRPR3
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
%   see S4PRPR3_convert_par2_MPV_fixb.m
% 
% Output:
% tau [4x1]
%   joint torques of inverse dynamics (contains inertial, gravitational Coriolis and centrifugal forces)

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-31 16:21
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tau = S4PRPR3_invdynJ_fixb_mdp_slag_vp(qJ, qJD, qJDD, g, pkin, MDP)
%% Coder Information
%#codegen
%$cgargs {zeros(4,1),zeros(4,1),zeros(4,1),zeros(3,1),zeros(7,1),zeros(15,1)}
assert(isreal(qJ) && all(size(qJ) == [4 1]), ...
  'S4PRPR3_invdynJ_fixb_mdp_slag_vp: qJ has to be [4x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [4 1]), ...
  'S4PRPR3_invdynJ_fixb_mdp_slag_vp: qJD has to be [4x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [4 1]), ...
  'S4PRPR3_invdynJ_fixb_mdp_slag_vp: qJDD has to be [4x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S4PRPR3_invdynJ_fixb_mdp_slag_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [7 1]), ...
  'S4PRPR3_invdynJ_fixb_mdp_slag_vp: pkin has to be [7x1] (double)');
assert(isreal(MDP) && all(size(MDP) == [15 1]), ...
  'S4PRPR3_invdynJ_fixb_mdp_slag_vp: MDP has to be [15x1] (double)'); 

%% Symbolic Calculation
% From invdyn_joint_fixb_mdp_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 16:20:57
% EndTime: 2019-12-31 16:20:58
% DurationCPUTime: 0.46s
% Computational Cost: add. (271->105), mult. (534->147), div. (0->0), fcn. (378->8), ass. (0->58)
t135 = qJD(2) * qJD(3);
t148 = qJ(3) * qJDD(2) + t135;
t115 = pkin(6) + qJ(2);
t109 = sin(t115);
t111 = cos(t115);
t129 = g(1) * t111 + g(2) * t109;
t128 = g(1) * t109 - g(2) * t111;
t147 = t128 - qJDD(3);
t144 = pkin(5) + qJ(3);
t143 = qJDD(2) * pkin(2);
t117 = cos(pkin(7));
t142 = t117 * MDP(5);
t116 = sin(pkin(7));
t141 = t116 ^ 2 + t117 ^ 2;
t140 = qJD(2) * t116;
t139 = qJD(2) * t117;
t118 = sin(qJ(4));
t138 = qJD(2) * t118;
t119 = cos(qJ(4));
t137 = qJD(2) * t119;
t134 = qJDD(2) * t117;
t133 = qJDD(2) * t118;
t132 = qJDD(2) * t119;
t131 = t116 * t132 + (qJD(4) * t137 + t133) * t117;
t79 = qJ(3) * t134 + t116 * qJDD(1) + t117 * t135;
t102 = -pkin(3) * t117 - pkin(2);
t130 = t116 * t138;
t97 = t117 * t132;
t127 = -t116 * t133 + t97;
t92 = t144 * t116;
t93 = t144 * t117;
t126 = t118 * t93 + t119 * t92;
t125 = -t118 * t92 + t119 * t93;
t87 = t116 * t119 + t117 * t118;
t86 = t116 * t118 - t119 * t117;
t124 = t143 + t147;
t85 = t87 * qJD(4);
t123 = -t116 * t137 - t117 * t138;
t105 = t117 * qJDD(1);
t78 = -t148 * t116 + t105;
t121 = -t116 * t78 + t117 * t79 - t129;
t114 = pkin(7) + qJ(4);
t110 = cos(t114);
t108 = sin(t114);
t91 = t102 * qJD(2) + qJD(3);
t90 = t102 * qJDD(2) + qJDD(3);
t89 = qJ(3) * t139 + t116 * qJD(1);
t88 = -qJ(3) * t140 + t117 * qJD(1);
t84 = t86 * qJD(4);
t83 = t87 * qJD(2);
t82 = t86 * qJD(2);
t77 = pkin(5) * t134 + t79;
t76 = t105 + (-t144 * qJDD(2) - t135) * t116;
t75 = -qJD(4) * t85 - qJDD(4) * t86;
t74 = -qJD(4) * t84 + qJDD(4) * t87;
t73 = qJD(2) * t85 - t127;
t72 = -qJD(4) * t130 + t131;
t1 = [(qJDD(1) - g(3)) * MDP(1) + (t116 * t79 + t117 * t78 - g(3)) * MDP(8) + t75 * MDP(14) - t74 * MDP(15); qJDD(2) * MDP(2) + t128 * MDP(3) + t129 * MDP(4) + (t148 * t141 + t121) * MDP(7) + ((-t116 * t88 + t117 * t89) * qJD(3) + t124 * pkin(2) + t121 * qJ(3)) * MDP(8) + (t72 * t87 - t83 * t84) * MDP(9) + (-t72 * t86 - t73 * t87 + t82 * t84 - t83 * t85) * MDP(10) + t74 * MDP(11) + t75 * MDP(12) + (t102 * t73 + t90 * t86 + t91 * t85 - t126 * qJDD(4) + t128 * t110 + (-t87 * qJD(3) - t125 * qJD(4)) * qJD(4)) * MDP(14) + (t102 * t72 + t90 * t87 - t91 * t84 - t125 * qJDD(4) - t128 * t108 + (t86 * qJD(3) + t126 * qJD(4)) * qJD(4)) * MDP(15) + (-MDP(6) * t116 + t142) * (t124 + t143); (-t89 * t139 + t88 * t140 - t147) * MDP(8) - t97 * MDP(14) + t131 * MDP(15) - t141 * MDP(7) * qJD(2) ^ 2 + (-t142 - pkin(2) * MDP(8) + (t118 * MDP(14) + MDP(6)) * t116) * qJDD(2) + ((-t123 + t83) * MDP(14) + (-t82 - t130) * MDP(15)) * qJD(4); t83 * t82 * MDP(9) + (-t82 ^ 2 + t83 ^ 2) * MDP(10) + t131 * MDP(11) + t127 * MDP(12) + qJDD(4) * MDP(13) + (-g(3) * t110 + t129 * t108 - t118 * t77 + t119 * t76 - t91 * t83) * MDP(14) + (g(3) * t108 + t129 * t110 - t118 * t76 - t119 * t77 + t91 * t82) * MDP(15) + ((t82 - t130) * MDP(11) + (t123 + t83) * MDP(12)) * qJD(4);];
tau = t1;
