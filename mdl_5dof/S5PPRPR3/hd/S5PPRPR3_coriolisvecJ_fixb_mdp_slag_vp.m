% Calculate minimal parameter regressor of Coriolis joint torque vector for
% S5PPRPR3
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% qJD [5x1]
%   Generalized joint velocities
% pkin [9x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d3,d5,theta1,theta2,theta4]';
% MDP [13x1]
%   Minimal dynamic parameter vector (fixed base model)
%   see S5PPRPR3_convert_par2_MPV_fixb.m
% 
% Output:
% tauc [5x1]
%   joint torques required to compensate Coriolis and centrifugal load

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-05 15:05
% Revision: 77da58f92bca3eff71542919beafa37024070d86 (2019-12-05)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauc = S5PPRPR3_coriolisvecJ_fixb_mdp_slag_vp(qJ, qJD, pkin, MDP)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(9,1),zeros(13,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5PPRPR3_coriolisvecJ_fixb_mdp_slag_vp: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5PPRPR3_coriolisvecJ_fixb_mdp_slag_vp: qJD has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S5PPRPR3_coriolisvecJ_fixb_mdp_slag_vp: pkin has to be [9x1] (double)');
assert(isreal(MDP) && all(size(MDP) == [13 1]), ...
  'S5PPRPR3_coriolisvecJ_fixb_mdp_slag_vp: MDP has to be [13x1] (double)'); 

%% Symbolic Calculation
% From coriolisvec_joint_fixb_mdp_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-05 15:05:23
% EndTime: 2019-12-05 15:05:25
% DurationCPUTime: 0.31s
% Computational Cost: add. (215->63), mult. (590->119), div. (0->0), fcn. (444->8), ass. (0->46)
t100 = sin(qJ(3));
t102 = cos(qJ(3));
t96 = sin(pkin(8));
t117 = qJD(1) * t96;
t127 = t102 * qJD(2) - t100 * t117;
t101 = cos(qJ(5));
t99 = sin(qJ(5));
t126 = (-t101 ^ 2 + t99 ^ 2) * MDP(8);
t95 = sin(pkin(9));
t97 = cos(pkin(9));
t84 = t95 * t100 - t97 * t102;
t87 = t100 * qJD(2) + t102 * t117;
t125 = t95 * t87;
t124 = t97 * t87;
t103 = qJD(5) ^ 2;
t85 = t97 * t100 + t95 * t102;
t122 = t103 * t85;
t118 = t99 * MDP(7);
t116 = qJD(3) * t96;
t73 = t84 * t116;
t115 = t73 * qJD(3);
t114 = t99 * MDP(12);
t113 = t99 * qJD(5);
t112 = t101 * qJD(5);
t105 = t87 * qJD(3);
t83 = t127 * qJD(3);
t67 = t97 * t105 + t95 * t83;
t110 = t103 * (t95 * pkin(3) + pkin(6)) + t67;
t80 = qJD(3) * pkin(3) + t127;
t69 = t97 * t80 - t125;
t65 = -qJD(3) * pkin(4) - t69;
t72 = t127 * t97 - t125;
t109 = qJD(5) * (t65 + t72);
t68 = -t95 * t105 + t97 * t83;
t108 = -t65 * qJD(3) - t68;
t74 = t85 * t116;
t77 = t85 * t96;
t107 = t77 * qJD(3) + qJD(5) * cos(pkin(8)) + t74;
t104 = qJD(3) ^ 2;
t91 = -t97 * pkin(3) - pkin(4);
t82 = t84 * qJD(3);
t81 = t85 * qJD(3);
t78 = t84 * t96;
t71 = t127 * t95 + t124;
t70 = t95 * t80 + t124;
t1 = [(t67 * t77 - t68 * t78 + t69 * t73 - t70 * t74) * MDP(6) + (t101 * t115 + (t107 * t99 + t78 * t112) * qJD(5)) * MDP(12) + (-t99 * t115 + (t107 * t101 - t78 * t113) * qJD(5)) * MDP(13) + (-MDP(4) * t102 + MDP(5) * t100) * t96 * t104; (t67 * t84 + t68 * t85 - t69 * t81 - t70 * t82) * MDP(6) + (-t101 * t122 + t82 * t113) * MDP(12) + (t82 * t112 + t99 * t122) * MDP(13) + (-t100 * MDP(4) - t102 * MDP(5)) * t104 + ((-t101 * t81 + t84 * t113) * MDP(12) + (t84 * t112 + t81 * t99) * MDP(13)) * qJD(3); (t69 * t71 - t70 * t72 + (-t67 * t97 + t68 * t95) * pkin(3)) * MDP(6) + (-t103 * MDP(10) + MDP(12) * t109 + t110 * MDP(13)) * t99 + (-t110 * MDP(12) + MDP(13) * t109 + t103 * MDP(9)) * t101 + ((t101 * MDP(12) - t99 * MDP(13)) * t71 + (-0.2e1 * t126 + t91 * t114 + (t91 * MDP(13) + 0.2e1 * t118) * t101) * qJD(5)) * qJD(3); (-t101 * MDP(13) - t114) * t103; t108 * t114 + t104 * t126 + (t108 * MDP(13) - t104 * t118) * t101;];
tauc = t1;
