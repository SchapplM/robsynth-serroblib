% Calculate vector of inverse dynamics joint torques for
% S4RPRP1
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
%   pkin=[a2,a3,a4,d1,d3,theta2]';
% MDP [10x1]
%   Minimal dynamic parameter vector (fixed base model)
%   see S4RPRP1_convert_par2_MPV_fixb.m
% 
% Output:
% tau [4x1]
%   joint torques of inverse dynamics (contains inertial, gravitational Coriolis and centrifugal forces)

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-08 18:30
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tau = S4RPRP1_invdynJ_fixb_mdp_slag_vp(qJ, qJD, qJDD, g, pkin, MDP)
%% Coder Information
%#codegen
%$cgargs {zeros(4,1),zeros(4,1),zeros(4,1),zeros(3,1),zeros(6,1),zeros(10,1)}
assert(isreal(qJ) && all(size(qJ) == [4 1]), ...
  'S4RPRP1_invdynJ_fixb_mdp_slag_vp: qJ has to be [4x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [4 1]), ...
  'S4RPRP1_invdynJ_fixb_mdp_slag_vp: qJD has to be [4x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [4 1]), ...
  'S4RPRP1_invdynJ_fixb_mdp_slag_vp: qJDD has to be [4x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S4RPRP1_invdynJ_fixb_mdp_slag_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [6 1]), ...
  'S4RPRP1_invdynJ_fixb_mdp_slag_vp: pkin has to be [6x1] (double)');
assert(isreal(MDP) && all(size(MDP) == [10 1]), ...
  'S4RPRP1_invdynJ_fixb_mdp_slag_vp: MDP has to be [10x1] (double)'); 

%% Symbolic Calculation
% From invdyn_joint_fixb_mdp_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-08 18:29:44
% EndTime: 2019-03-08 18:29:44
% DurationCPUTime: 0.24s
% Computational Cost: add. (309->84), mult. (468->99), div. (0->0), fcn. (238->10), ass. (0->53)
t114 = cos(pkin(6));
t104 = pkin(1) * t114 + pkin(2);
t113 = sin(pkin(6));
t143 = pkin(1) * t113;
t132 = qJD(3) * t143;
t146 = -qJD(1) * t132 + t104 * qJDD(1);
t111 = qJ(1) + pkin(6);
t108 = qJ(3) + t111;
t102 = sin(t108);
t103 = cos(t108);
t145 = -g(1) * t103 - g(2) * t102;
t144 = g(1) * t102 - g(2) * t103;
t109 = qJDD(1) + qJDD(3);
t142 = pkin(3) * t109;
t116 = sin(qJ(1));
t140 = t116 * pkin(1);
t118 = cos(qJ(1));
t139 = t118 * pkin(1);
t138 = t103 * pkin(3) + t102 * qJ(4);
t115 = sin(qJ(3));
t91 = t104 * qJD(1);
t137 = t115 * t91;
t117 = cos(qJ(3));
t136 = t115 * t104 + t117 * t143;
t133 = qJD(1) * t143;
t77 = -t115 * t133 + t117 * t91;
t135 = qJD(4) - t77;
t134 = qJD(3) * t117;
t131 = qJDD(1) * t143;
t130 = -t102 * pkin(3) + t103 * qJ(4);
t129 = t146 * t115 + t117 * t131 + t91 * t134;
t128 = -qJD(3) * t137 - t115 * t131 + t146 * t117;
t126 = t129 + t145;
t125 = t104 * t134 - t115 * t132;
t105 = t109 * qJ(4);
t110 = qJD(1) + qJD(3);
t107 = t110 * qJD(4);
t73 = t105 + t107 + t129;
t124 = t104 * t117 - t115 * t143;
t123 = t128 + t144;
t78 = t117 * t133 + t137;
t74 = qJDD(4) - t128 - t142;
t122 = t110 * t77 - t126;
t121 = t110 * t78 + t123;
t80 = t136 * qJD(3);
t120 = -t80 * t110 + t123;
t106 = t109 * MDP(5);
t82 = -pkin(3) - t124;
t81 = qJ(4) + t136;
t79 = qJD(4) + t125;
t76 = qJ(4) * t110 + t78;
t75 = -pkin(3) * t110 + t135;
t1 = [(g(1) * t116 - g(2) * t118) * MDP(2) + (g(1) * t118 + g(2) * t116) * MDP(3) + (g(1) * t140 - g(2) * t139) * MDP(4) + t106 + (t124 * t109 + t120) * MDP(6) + (-t136 * t109 - t125 * t110 - t126) * MDP(7) + (-qJDD(4) + (pkin(3) - t82) * t109 + t120) * MDP(8) + (t109 * t81 + t110 * t79 + t145 + t73) * MDP(9) + (t73 * t81 + t76 * t79 + t74 * t82 + t75 * t80 - g(1) * (-pkin(2) * sin(t111) - t140 + t130) - g(2) * (pkin(2) * cos(t111) + t139 + t138)) * MDP(10) + (MDP(1) + (t113 ^ 2 + t114 ^ 2) * MDP(4) * pkin(1) ^ 2) * qJDD(1); (MDP(4) + MDP(10)) * (qJDD(2) - g(3)); t106 + t121 * MDP(6) + t122 * MDP(7) + (-qJDD(4) + t121 + 0.2e1 * t142) * MDP(8) + (0.2e1 * t105 + 0.2e1 * t107 - t122) * MDP(9) + (-t74 * pkin(3) - g(1) * t130 - g(2) * t138 + t73 * qJ(4) + t135 * t76 - t75 * t78) * MDP(10); -t109 * MDP(8) - t110 ^ 2 * MDP(9) + (-t110 * t76 - t144 + t74) * MDP(10);];
tau  = t1;
