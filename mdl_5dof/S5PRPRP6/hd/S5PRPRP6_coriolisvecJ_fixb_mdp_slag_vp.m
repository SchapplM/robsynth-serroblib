% Calculate minimal parameter regressor of Coriolis joint torque vector for
% S5PRPRP6
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% qJD [5x1]
%   Generalized joint velocities
% pkin [7x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d2,d4,theta1]';
% MDP [18x1]
%   Minimal dynamic parameter vector (fixed base model)
%   see S5PRPRP6_convert_par2_MPV_fixb.m
% 
% Output:
% tauc [5x1]
%   joint torques required to compensate Coriolis and centrifugal load

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-05 15:41
% Revision: 77da58f92bca3eff71542919beafa37024070d86 (2019-12-05)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauc = S5PRPRP6_coriolisvecJ_fixb_mdp_slag_vp(qJ, qJD, pkin, MDP)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(7,1),zeros(18,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5PRPRP6_coriolisvecJ_fixb_mdp_slag_vp: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5PRPRP6_coriolisvecJ_fixb_mdp_slag_vp: qJD has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [7 1]), ...
  'S5PRPRP6_coriolisvecJ_fixb_mdp_slag_vp: pkin has to be [7x1] (double)');
assert(isreal(MDP) && all(size(MDP) == [18 1]), ...
  'S5PRPRP6_coriolisvecJ_fixb_mdp_slag_vp: MDP has to be [18x1] (double)'); 

%% Symbolic Calculation
% From coriolisvec_joint_fixb_mdp_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-05 15:41:21
% EndTime: 2019-12-05 15:41:23
% DurationCPUTime: 0.45s
% Computational Cost: add. (349->112), mult. (703->147), div. (0->0), fcn. (318->4), ass. (0->57)
t91 = sin(qJ(4));
t118 = qJD(4) * t91;
t94 = cos(qJ(2));
t116 = t94 * qJD(1);
t102 = qJD(3) - t116;
t95 = -pkin(2) - pkin(6);
t79 = t95 * qJD(2) + t102;
t92 = sin(qJ(2));
t117 = t92 * qJD(1);
t110 = qJD(2) * t117;
t93 = cos(qJ(4));
t88 = t93 * t110;
t71 = t79 * t118 - t88;
t114 = qJD(4) * qJ(5);
t128 = t79 * t91;
t76 = t114 + t128;
t132 = qJD(4) * t76 - t71;
t101 = pkin(4) * t93 + qJ(5) * t91;
t77 = t101 * qJD(4) - t93 * qJD(5) + qJD(3);
t106 = -t77 + t116;
t89 = t91 ^ 2;
t90 = t93 ^ 2;
t131 = (t89 - t90) * MDP(9);
t130 = (t89 + t90) * MDP(16);
t129 = 0.2e1 * qJD(4);
t127 = t79 * t93;
t126 = t91 * t93;
t96 = qJD(4) ^ 2;
t125 = t95 * t96;
t97 = qJD(2) ^ 2;
t122 = t96 + t97;
t81 = pkin(4) * t91 - qJ(5) * t93 + qJ(3);
t121 = qJD(2) * t81;
t120 = qJD(2) * t93;
t115 = qJD(2) * qJ(3);
t113 = MDP(13) + MDP(15);
t112 = MDP(14) - MDP(17);
t111 = t122 * t94;
t109 = qJD(2) * qJD(4) * t92;
t69 = (t77 + t116) * qJD(2);
t108 = -t69 + t125;
t82 = (qJD(3) + t116) * qJD(2);
t107 = t82 - t125;
t86 = t115 + t117;
t105 = -t86 + t117;
t104 = qJD(4) * pkin(4) - qJD(5);
t103 = MDP(18) * t95 - MDP(16);
t75 = t117 + t121;
t100 = -t117 + t75 + t121;
t99 = -t105 + t115;
t87 = t91 * t110;
t70 = t87 + (qJD(5) + t127) * qJD(4);
t74 = -t104 - t127;
t98 = t74 * t118 + t132 * t93 + t70 * t91;
t85 = -qJD(2) * pkin(2) + t102;
t80 = t101 * qJD(2);
t1 = [t112 * (-0.2e1 * t91 * t109 + t93 * t111) + t113 * (0.2e1 * t93 * t109 + t91 * t111) + (t86 * qJD(2) * MDP(7) + (qJD(2) * t75 - t98) * MDP(18) + (-MDP(4) + MDP(6)) * t97) * t94 + (t69 * MDP(18) + t82 * MDP(7) + ((t85 - t116) * MDP(7) + (-t74 * t93 + t76 * t91) * MDP(18)) * qJD(2) + (-MDP(3) + MDP(5) - t130) * t97) * t92; (qJ(3) * t82 + t102 * t86 - t85 * t117) * MDP(7) + (-t106 * t75 + t69 * t81) * MDP(18) + (0.2e1 * qJD(3) * MDP(6) + t129 * t131 + (-pkin(2) * MDP(7) + t130) * t117) * qJD(2) + (-t96 * MDP(11) + t107 * MDP(14) + t71 * MDP(16) + t108 * MDP(17) + (t74 * t117 - t71 * t95) * MDP(18) + (t102 * MDP(14) + t106 * MDP(17)) * qJD(2) + (t99 * MDP(13) + t100 * MDP(15) + t103 * t76) * qJD(4)) * t93 + (-t96 * MDP(10) + t107 * MDP(13) - t108 * MDP(15) - t70 * MDP(16) + (-t76 * t117 + t70 * t95) * MDP(18) + (t102 * MDP(13) - t106 * MDP(15)) * qJD(2) + (-t99 * MDP(14) + t100 * MDP(17) - 0.2e1 * MDP(8) * t120 + t103 * t74) * qJD(4)) * t91; -t97 * MDP(6) + t98 * MDP(18) + (-t75 * MDP(18) + t105 * MDP(7)) * qJD(2) + (-t112 * t93 - t113 * t91) * t122; -t87 * MDP(14) + (qJD(5) * t129 + t87) * MDP(17) + (-t74 * t128 - pkin(4) * t71 + qJ(5) * t70 - t75 * t80 + (qJD(5) - t127) * t76) * MDP(18) + t113 * t88 + (MDP(8) * t126 - t131) * t97 + ((-t86 * MDP(13) - t75 * MDP(15) + (t76 - t114) * MDP(16) + t80 * MDP(17)) * t93 + (t86 * MDP(14) - t80 * MDP(15) + (t104 + t74) * MDP(16) - t75 * MDP(17)) * t91) * qJD(2); t97 * MDP(15) * t126 + (-t90 * t97 - t96) * MDP(17) + (t75 * t120 - t132) * MDP(18);];
tauc = t1;
