% Calculate Coriolis joint torque vector for
% S5RPPRP3
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% qJD [5x1]
%   Generalized joint velocities
% pkin [7x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d1,d4,theta2]';
% MDP [18x1]
%   Minimal dynamic parameter vector (fixed base model)
%   see S5RPPRP3_convert_par2_MPV_fixb.m
% 
% Output:
% tauc [5x1]
%   joint torques required to compensate Coriolis and centrifugal load

% Quelle: HybrDyn-Toolbox
% Datum: 2021-01-15 17:04
% Revision: 24b2e7d74a0c1a3b64fa2f8f5ad758691ad61af3 (2021-01-15)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauc = S5RPPRP3_coriolisvecJ_fixb_mdp_slag_vp(qJ, qJD, pkin, MDP)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(7,1),zeros(18,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RPPRP3_coriolisvecJ_fixb_mdp_slag_vp: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RPPRP3_coriolisvecJ_fixb_mdp_slag_vp: qJD has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [7 1]), ...
  'S5RPPRP3_coriolisvecJ_fixb_mdp_slag_vp: pkin has to be [7x1] (double)');
assert(isreal(MDP) && all(size(MDP) == [18 1]), ...
  'S5RPPRP3_coriolisvecJ_fixb_mdp_slag_vp: MDP has to be [18x1] (double)'); 

%% Symbolic Calculation
% From coriolisvec_joint_fixb_mdp_matlab.m
% OptimizationMode: 2
% StartTime: 2021-01-15 17:04:43
% EndTime: 2021-01-15 17:04:45
% DurationCPUTime: 0.36s
% Computational Cost: add. (363->100), mult. (692->138), div. (0->0), fcn. (294->4), ass. (0->57)
t85 = -cos(pkin(7)) * pkin(1) - pkin(2) - pkin(6);
t120 = qJ(5) - t85;
t95 = sin(qJ(4));
t75 = t120 * t95;
t91 = t95 ^ 2;
t96 = cos(qJ(4));
t92 = t96 ^ 2;
t129 = MDP(9) * (t91 - t92);
t128 = 2 * qJD(1);
t77 = qJD(1) * t85 + qJD(3);
t127 = t95 * t77;
t73 = t96 * t77;
t98 = qJD(1) ^ 2;
t126 = t98 * t95;
t122 = qJD(4) * pkin(4);
t102 = -t95 * qJD(2) + t73;
t112 = qJ(5) * qJD(1);
t69 = -t112 * t96 + t102;
t68 = t69 + t122;
t125 = t68 - t69;
t115 = t96 * qJD(4);
t108 = pkin(4) * t115;
t90 = qJD(3) * qJD(1);
t78 = qJD(1) * t108 + t90;
t121 = t96 * MDP(8);
t87 = sin(pkin(7)) * pkin(1) + qJ(3);
t79 = t95 * pkin(4) + t87;
t119 = qJD(1) * t79;
t82 = qJD(1) * t87;
t74 = qJD(5) + t119;
t118 = t74 * qJD(1);
t117 = t82 * qJD(1);
t116 = t95 * qJD(5);
t114 = t96 * qJD(5);
t113 = qJD(5) + t74;
t111 = qJD(2) * qJD(4);
t110 = MDP(13) + MDP(15);
t109 = MDP(14) + MDP(16);
t107 = qJ(5) * t115;
t106 = t95 * t112;
t76 = t120 * t96;
t105 = t74 + t119;
t83 = qJD(3) + t108;
t104 = qJD(1) * t83 + t78;
t103 = 0.2e1 * t82;
t100 = -t96 * qJD(2) - t127;
t70 = -t100 - t106;
t101 = t68 * t96 + t70 * t95;
t97 = qJD(4) ^ 2;
t99 = -t85 * t97 + 0.2e1 * t90;
t88 = t95 * t111;
t84 = qJD(4) * t106;
t72 = -qJD(4) * t76 - t116;
t71 = qJD(4) * t75 - t114;
t67 = -qJD(1) * t114 + qJD(4) * t100 + t84;
t66 = t77 * t115 - t88 + (-t107 - t116) * qJD(1);
t1 = [(-t66 * t75 - t67 * t76 + t68 * t71 + t70 * t72 + t74 * t83 + t78 * t79) * MDP(18) + ((MDP(6) * t128) + MDP(7) * t103) * qJD(3) + (-t97 * MDP(11) + t99 * MDP(14) + t104 * MDP(16) + (-qJD(1) * t71 - t67) * MDP(17)) * t96 + (-t97 * MDP(10) + t99 * MDP(13) + t104 * MDP(15) + (-qJD(1) * t72 - t66) * MDP(17)) * t95 + (t71 * MDP(15) - t72 * MDP(16) + t128 * t129 + (t103 * MDP(13) + t105 * MDP(15) + (qJD(1) * t75 - t70) * MDP(17)) * t96 + (-0.2e1 * qJD(1) * t121 - t103 * MDP(14) - t105 * MDP(16) + (-qJD(1) * t76 + t68) * MDP(17)) * t95) * qJD(4); (-qJD(4) * t101 + t66 * t96 - t67 * t95) * MDP(18) + (t109 * t95 - t110 * t96) * t97; -(t98 * MDP(6)) - MDP(7) * t117 + (-t118 + t66 * t95 + t67 * t96 + (-t68 * t95 + t70 * t96) * qJD(4)) * MDP(18) + (t109 * t96 + t110 * t95) * (-t97 - t98); t121 * t126 - t98 * t129 - t96 * MDP(13) * t117 + (t95 * t117 + t88 + (t102 - t73) * qJD(4)) * MDP(14) + (t84 + (t70 - t127) * qJD(4) + (-pkin(4) * t126 - qJD(1) * t113 - t111) * t96) * MDP(15) + (-t92 * t98 * pkin(4) + t88 + (t69 - t73) * qJD(4) + (t113 * t95 + t107) * qJD(1)) * MDP(16) + (t122 - t125) * t95 * qJD(1) * MDP(17) + (t125 * t70 + (-t118 * t96 + t67) * pkin(4)) * MDP(18); t78 * MDP(18) + (-t91 - t92) * MDP(17) * t98 + (t101 * MDP(18) + 0.2e1 * (t96 * MDP(15) - t95 * MDP(16)) * qJD(4)) * qJD(1);];
tauc = t1;
