% Calculate vector of inverse dynamics joint torques for
% S5PPRPR5
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
%   pkin=[a2,a3,a4,a5,d3,d5,theta1]';
% MDP [15x1]
%   Minimal dynamic parameter vector (fixed base model)
%   see S5PPRPR5_convert_par2_MPV_fixb.m
% 
% Output:
% tau [5x1]
%   joint torques of inverse dynamics (contains inertial, gravitational Coriolis and centrifugal forces)

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-31 17:33
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tau = S5PPRPR5_invdynJ_fixb_mdp_slag_vp(qJ, qJD, qJDD, g, pkin, MDP)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(5,1),zeros(3,1),zeros(7,1),zeros(15,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5PPRPR5_invdynJ_fixb_mdp_slag_vp: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5PPRPR5_invdynJ_fixb_mdp_slag_vp: qJD has to be [5x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [5 1]), ...
  'S5PPRPR5_invdynJ_fixb_mdp_slag_vp: qJDD has to be [5x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S5PPRPR5_invdynJ_fixb_mdp_slag_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [7 1]), ...
  'S5PPRPR5_invdynJ_fixb_mdp_slag_vp: pkin has to be [7x1] (double)');
assert(isreal(MDP) && all(size(MDP) == [15 1]), ...
  'S5PPRPR5_invdynJ_fixb_mdp_slag_vp: MDP has to be [15x1] (double)'); 

%% Symbolic Calculation
% From invdyn_joint_fixb_mdp_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 17:33:32
% EndTime: 2019-12-31 17:33:33
% DurationCPUTime: 0.29s
% Computational Cost: add. (206->96), mult. (372->118), div. (0->0), fcn. (235->6), ass. (0->49)
t111 = qJD(3) * qJ(4);
t84 = sin(qJ(3));
t114 = t84 * qJD(2);
t73 = t111 + t114;
t87 = -pkin(3) - pkin(6);
t122 = (t111 + t73 - t114) * qJD(5) + qJDD(5) * t87;
t83 = sin(qJ(5));
t85 = cos(qJ(5));
t121 = t83 * t85;
t81 = t85 ^ 2;
t120 = t83 ^ 2 - t81;
t88 = qJD(5) ^ 2;
t89 = qJD(3) ^ 2;
t119 = t88 + t89;
t118 = cos(pkin(7));
t117 = sin(pkin(7));
t116 = qJDD(3) * pkin(3);
t115 = t73 * qJD(3);
t86 = cos(qJ(3));
t113 = t86 * qJD(2);
t112 = qJDD(1) - g(3);
t110 = qJDD(3) * t84;
t109 = qJDD(5) * t83;
t107 = t84 * qJDD(2);
t106 = t86 * qJDD(2);
t105 = qJD(3) * qJD(5);
t104 = qJDD(3) * qJ(4);
t103 = qJD(4) - t113;
t65 = -t117 * t84 - t118 * t86;
t66 = -t117 * t86 + t118 * t84;
t102 = g(1) * t66 - g(2) * t65;
t101 = g(1) * t65 + g(2) * t66;
t100 = -g(1) * t117 + g(2) * t118;
t99 = (-qJD(3) * pkin(3) + t103) * t84 + t73 * t86;
t77 = qJD(3) * t114;
t98 = qJDD(4) + t77 - t106;
t96 = t119 * t86 + t110;
t95 = -qJDD(5) * t86 + 0.2e1 * t84 * t105;
t94 = -t101 - t107;
t93 = t102 + t106;
t92 = -t87 * qJDD(3) + t102 + t115 - t98;
t91 = qJDD(4) - t93;
t63 = t104 + t107 + (qJD(4) + t113) * qJD(3);
t90 = t103 * qJD(3) - t87 * t88 + t101 + t104 + t63;
t78 = qJDD(5) * t85;
t69 = -t88 * t83 + t78;
t68 = t88 * t85 + t109;
t64 = t98 - t116;
t1 = [t68 * MDP(14) + t69 * MDP(15) + (MDP(1) + MDP(2) + MDP(8)) * t112; (qJDD(2) + t100) * MDP(2) + (t99 * qJD(3) + t63 * t84 - t64 * t86 + t100) * MDP(8) + (t96 * t83 + t95 * t85) * MDP(14) + (-t95 * t83 + t96 * t85) * MDP(15) + (-MDP(4) + MDP(6)) * (-t86 * qJDD(3) + t89 * t84) + (-MDP(5) + MDP(7)) * (t89 * t86 + t110); qJDD(3) * MDP(3) + t93 * MDP(4) + t94 * MDP(5) + (t91 - 0.2e1 * t116) * MDP(6) + (0.2e1 * qJD(3) * qJD(4) + 0.2e1 * t104 - t94) * MDP(7) + (t63 * qJ(4) + t73 * qJD(4) - t64 * pkin(3) - g(1) * (-t66 * pkin(3) - t65 * qJ(4)) - g(2) * (t65 * pkin(3) - t66 * qJ(4)) - t99 * qJD(2)) * MDP(8) + (t81 * qJDD(3) - 0.2e1 * t105 * t121) * MDP(9) + 0.2e1 * (-qJDD(3) * t121 + t120 * t105) * MDP(10) + t69 * MDP(11) - t68 * MDP(12) + (t122 * t85 + t90 * t83) * MDP(14) + (-t122 * t83 + t90 * t85) * MDP(15); qJDD(3) * MDP(6) - t89 * MDP(7) + (-t115 + t77 + t91 - t116) * MDP(8) + (-t119 * t83 + t78) * MDP(14) + (-t119 * t85 - t109) * MDP(15); qJDD(5) * MDP(13) - t120 * MDP(10) * t89 + (qJDD(3) * MDP(11) - t92 * MDP(14) + t112 * MDP(15)) * t85 + (t85 * t89 * MDP(9) - qJDD(3) * MDP(12) + t112 * MDP(14) + t92 * MDP(15)) * t83;];
tau = t1;
