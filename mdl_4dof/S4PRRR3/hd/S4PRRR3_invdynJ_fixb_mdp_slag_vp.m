% Calculate vector of inverse dynamics joint torques for
% S4PRRR3
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
%   pkin=[a2,a3,a4,d2,d3,d4,theta1]';
% MDP [14x1]
%   Minimal dynamic parameter vector (fixed base model)
%   see S4PRRR3_convert_par2_MPV_fixb.m
% 
% Output:
% tau [4x1]
%   joint torques of inverse dynamics (contains inertial, gravitational Coriolis and centrifugal forces)

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-31 16:31
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tau = S4PRRR3_invdynJ_fixb_mdp_slag_vp(qJ, qJD, qJDD, g, pkin, MDP)
%% Coder Information
%#codegen
%$cgargs {zeros(4,1),zeros(4,1),zeros(4,1),zeros(3,1),zeros(7,1),zeros(14,1)}
assert(isreal(qJ) && all(size(qJ) == [4 1]), ...
  'S4PRRR3_invdynJ_fixb_mdp_slag_vp: qJ has to be [4x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [4 1]), ...
  'S4PRRR3_invdynJ_fixb_mdp_slag_vp: qJD has to be [4x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [4 1]), ...
  'S4PRRR3_invdynJ_fixb_mdp_slag_vp: qJDD has to be [4x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S4PRRR3_invdynJ_fixb_mdp_slag_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [7 1]), ...
  'S4PRRR3_invdynJ_fixb_mdp_slag_vp: pkin has to be [7x1] (double)');
assert(isreal(MDP) && all(size(MDP) == [14 1]), ...
  'S4PRRR3_invdynJ_fixb_mdp_slag_vp: MDP has to be [14x1] (double)'); 

%% Symbolic Calculation
% From invdyn_joint_fixb_mdp_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 16:31:40
% EndTime: 2019-12-31 16:31:41
% DurationCPUTime: 0.33s
% Computational Cost: add. (242->80), mult. (348->116), div. (0->0), fcn. (168->8), ass. (0->42)
t100 = cos(qJ(3));
t118 = pkin(2) * t100;
t93 = pkin(7) + qJ(2);
t90 = qJ(3) + t93;
t83 = cos(t90);
t124 = -g(2) * t83 + qJDD(2) * t118;
t98 = sin(qJ(3));
t121 = pkin(2) * t98;
t110 = qJD(2) * t121;
t92 = qJDD(2) + qJDD(3);
t120 = pkin(3) * t92;
t123 = qJD(3) * t110 - t120 - t124;
t82 = sin(t90);
t122 = g(1) * t83 + g(2) * t82;
t80 = g(1) * t82;
t97 = sin(qJ(4));
t95 = t97 ^ 2;
t99 = cos(qJ(4));
t117 = -t99 ^ 2 + t95;
t94 = qJD(2) + qJD(3);
t116 = qJD(4) * t94;
t115 = t99 * qJD(4);
t113 = qJDD(1) - g(3);
t112 = qJD(2) * t100;
t111 = qJDD(2) * t98;
t109 = pkin(2) * t112;
t101 = qJD(4) ^ 2;
t84 = pkin(6) + t121;
t85 = -pkin(3) - t118;
t107 = t101 * t84 + t85 * t92;
t106 = -qJDD(4) * t84 + t85 * t116;
t74 = -pkin(3) * t94 - t109;
t105 = -t74 * t94 - pkin(6) * t92 - (qJD(3) * t112 + t111) * pkin(2) + t122;
t104 = pkin(6) * t101 - t94 * t110 - t120;
t103 = -pkin(3) * t116 - pkin(6) * qJDD(4) + qJD(4) * t109;
t75 = qJDD(4) * t97 + t101 * t99;
t76 = qJDD(4) * t99 - t101 * t97;
t102 = t92 * MDP(5) + (t80 + t124) * MDP(6) + t122 * MDP(7) + (0.2e1 * t94 * t97 * t115 + t92 * t95) * MDP(8) + 0.2e1 * (t92 * t97 * t99 - t117 * t116) * MDP(9) + t75 * MDP(10) + t76 * MDP(11) + (t74 * qJD(4) * t97 + t99 * t80) * MDP(13) + (t74 * t115 + t123 * t97) * MDP(14);
t91 = t94 ^ 2;
t89 = cos(t93);
t88 = sin(t93);
t1 = [t113 * MDP(1) + t76 * MDP(13) - t75 * MDP(14); (g(1) * t88 - g(2) * t89) * MDP(3) + (g(1) * t89 + g(2) * t88) * MDP(4) + (t100 * t92 * MDP(6) + (-qJDD(2) - t92) * MDP(7) * t98 + ((-qJD(2) * MDP(6) + (-t99 * MDP(13) + t97 * MDP(14) - MDP(6)) * t94) * t98 + ((-qJD(2) - t94) * MDP(7) + (-t97 * MDP(13) - t99 * MDP(14)) * qJD(4)) * t100) * qJD(3)) * pkin(2) + t102 + ((-t107 - t123) * MDP(13) + t106 * MDP(14)) * t99 + (t106 * MDP(13) + (t107 - t80) * MDP(14)) * t97 + qJDD(2) * MDP(2); (-MDP(7) * t111 + (t98 * MDP(6) + t100 * MDP(7)) * qJD(2) * (-qJD(3) + t94)) * pkin(2) + (t103 * MDP(13) + (t104 - t80) * MDP(14)) * t97 + ((-t104 - t123) * MDP(13) + t103 * MDP(14)) * t99 + t102; qJDD(4) * MDP(12) + t117 * MDP(9) * t91 + (t92 * MDP(11) + t113 * MDP(13) + t105 * MDP(14)) * t99 + (-t91 * t99 * MDP(8) + t92 * MDP(10) + t105 * MDP(13) - t113 * MDP(14)) * t97;];
tau = t1;
