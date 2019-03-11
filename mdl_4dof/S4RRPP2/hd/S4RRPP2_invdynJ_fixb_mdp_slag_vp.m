% Calculate vector of inverse dynamics joint torques for
% S4RRPP2
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
% pkin [5x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,d1,d2]';
% MDP [12x1]
%   Minimal dynamic parameter vector (fixed base model)
%   see S4RRPP2_convert_par2_MPV_fixb.m
% 
% Output:
% tau [4x1]
%   joint torques of inverse dynamics (contains inertial, gravitational Coriolis and centrifugal forces)

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-08 18:34
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tau = S4RRPP2_invdynJ_fixb_mdp_slag_vp(qJ, qJD, qJDD, g, pkin, MDP)
%% Coder Information
%#codegen
%$cgargs {zeros(4,1),zeros(4,1),zeros(4,1),zeros(3,1),zeros(5,1),zeros(12,1)}
assert(isreal(qJ) && all(size(qJ) == [4 1]), ...
  'S4RRPP2_invdynJ_fixb_mdp_slag_vp: qJ has to be [4x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [4 1]), ...
  'S4RRPP2_invdynJ_fixb_mdp_slag_vp: qJD has to be [4x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [4 1]), ...
  'S4RRPP2_invdynJ_fixb_mdp_slag_vp: qJDD has to be [4x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S4RRPP2_invdynJ_fixb_mdp_slag_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [5 1]), ...
  'S4RRPP2_invdynJ_fixb_mdp_slag_vp: pkin has to be [5x1] (double)');
assert(isreal(MDP) && all(size(MDP) == [12 1]), ...
  'S4RRPP2_invdynJ_fixb_mdp_slag_vp: MDP has to be [12x1] (double)'); 

%% Symbolic Calculation
% From invdyn_joint_fixb_mdp_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-08 18:34:06
% EndTime: 2019-03-08 18:34:07
% DurationCPUTime: 0.29s
% Computational Cost: add. (336->99), mult. (367->103), div. (0->0), fcn. (142->6), ass. (0->57)
t111 = qJ(1) + qJ(2);
t106 = sin(t111);
t107 = cos(t111);
t145 = -g(1) * t107 - g(2) * t106;
t144 = MDP(8) + MDP(11);
t116 = -pkin(2) - pkin(3);
t113 = sin(qJ(1));
t143 = g(1) * t113;
t109 = qJDD(1) + qJDD(2);
t101 = t109 * qJ(3);
t110 = qJD(1) + qJD(2);
t103 = t110 * qJD(3);
t112 = sin(qJ(2));
t114 = cos(qJ(2));
t137 = pkin(1) * qJD(2);
t128 = qJD(1) * t137;
t134 = pkin(1) * qJDD(1);
t140 = t112 * t134 + t114 * t128;
t77 = t101 + t103 + t140;
t138 = pkin(1) * qJD(1);
t130 = t112 * t138;
t83 = qJ(3) * t110 + t130;
t85 = t114 * t137 + qJD(3);
t91 = pkin(1) * t112 + qJ(3);
t142 = t77 * t91 + t83 * t85;
t141 = t77 * qJ(3) + t83 * qJD(3);
t139 = t107 * pkin(2) + t106 * qJ(3);
t136 = t114 * t83;
t135 = t112 * t128 - t114 * t134;
t133 = -MDP(7) - MDP(10);
t132 = MDP(9) + MDP(12);
t115 = cos(qJ(1));
t131 = t115 * pkin(1) + t139;
t129 = t114 * t138;
t98 = -pkin(1) * t114 - pkin(2);
t127 = t140 + t145;
t126 = g(1) * t106 - g(2) * t107 - t135;
t84 = t110 * t130;
t105 = t109 * pkin(2);
t78 = -t105 + qJDD(3) + t135;
t90 = t107 * qJ(3);
t125 = g(1) * (-pkin(2) * t106 + t90);
t124 = qJD(3) - t129;
t123 = -qJDD(3) + t126;
t122 = t77 + t145;
t121 = g(1) * (t116 * t106 + t90);
t120 = -t105 - t123;
t104 = t109 * pkin(3);
t119 = -t104 + t120;
t118 = t110 * t129 - t127;
t102 = t109 * MDP(4);
t92 = t107 * pkin(3);
t88 = -pkin(3) + t98;
t82 = -pkin(2) * t110 + t124;
t80 = t116 * t110 + t124;
t76 = -t104 + t78;
t1 = [qJDD(1) * MDP(1) + (-g(2) * t115 + t143) * MDP(2) + (g(1) * t115 + g(2) * t113) * MDP(3) + t102 + t126 * MDP(5) - t127 * MDP(6) + (-t109 * t98 - t120) * MDP(7) + (-g(2) * t131 + t78 * t98 - t125 + t142) * MDP(9) + (-t109 * t88 - t119) * MDP(10) + (t76 * t88 - t121 - g(2) * (t92 + t131) + t142) * MDP(12) + ((t114 * MDP(5) - t112 * MDP(6)) * t109 + t132 * t143 + (-t110 * t114 * MDP(6) + (t80 * MDP(12) + t82 * MDP(9) + (-MDP(5) + t133) * t110) * t112) * qJD(2)) * pkin(1) + t144 * (t109 * t91 + t110 * t85 + t122); t102 + (t84 + t126) * MDP(5) + t118 * MDP(6) + (0.2e1 * t105 + t84 + t123) * MDP(7) + (0.2e1 * t101 + 0.2e1 * t103 - t118) * MDP(8) + (-t78 * pkin(2) - t125 - g(2) * t139 + (-t112 * t82 - t136) * t138 + t141) * MDP(9) + (-t109 * t116 - t119 + t84) * MDP(10) + (t124 * t110 + t101 + t122) * MDP(11) + (t76 * t116 - t121 - g(2) * (t92 + t139) + (-t112 * t80 - t136) * t138 + t141) * MDP(12); t120 * MDP(9) + t119 * MDP(12) + t133 * t109 + (-t144 * t110 - t132 * t83) * t110; (qJDD(4) + g(3)) * MDP(12);];
tau  = t1;
