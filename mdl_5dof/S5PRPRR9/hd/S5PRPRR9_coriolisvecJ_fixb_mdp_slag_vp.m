% Calculate minimal parameter regressor of Coriolis joint torque vector for
% S5PRPRR9
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% qJD [5x1]
%   Generalized joint velocities
% pkin [8x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d2,d4,d5,theta1]';
% MDP [17x1]
%   Minimal dynamic parameter vector (fixed base model)
%   see S5PRPRR9_convert_par2_MPV_fixb.m
% 
% Output:
% tauc [5x1]
%   joint torques required to compensate Coriolis and centrifugal load

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-31 17:40
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauc = S5PRPRR9_coriolisvecJ_fixb_mdp_slag_vp(qJ, qJD, pkin, MDP)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(8,1),zeros(17,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5PRPRR9_coriolisvecJ_fixb_mdp_slag_vp: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5PRPRR9_coriolisvecJ_fixb_mdp_slag_vp: qJD has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S5PRPRR9_coriolisvecJ_fixb_mdp_slag_vp: pkin has to be [8x1] (double)');
assert(isreal(MDP) && all(size(MDP) == [17 1]), ...
  'S5PRPRR9_coriolisvecJ_fixb_mdp_slag_vp: MDP has to be [17x1] (double)'); 

%% Symbolic Calculation
% From coriolisvec_joint_fixb_mdp_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 17:39:46
% EndTime: 2019-12-31 17:39:47
% DurationCPUTime: 0.26s
% Computational Cost: add. (245->62), mult. (421->98), div. (0->0), fcn. (170->4), ass. (0->43)
t101 = qJ(3) * qJD(2);
t81 = sin(qJ(4));
t107 = qJD(4) * t81;
t84 = -pkin(2) - pkin(3);
t72 = t84 * qJD(2) + qJD(3);
t83 = cos(qJ(4));
t95 = qJD(4) * t101;
t99 = qJD(2) * qJD(3);
t60 = t72 * t107 + t81 * t99 + t83 * t95;
t98 = qJD(2) - qJD(4);
t118 = -t60 - (t83 * t101 + t81 * t72) * t98;
t91 = t83 * qJ(3) + t81 * t84;
t117 = t60 + (t81 * qJD(3) + t91 * qJD(4)) * t98;
t80 = sin(qJ(5));
t82 = cos(qJ(5));
t116 = (t80 ^ 2 - t82 ^ 2) * MDP(12);
t115 = 0.2e1 * qJD(3);
t114 = t98 * pkin(4);
t111 = t83 * t72;
t110 = qJD(4) * t111 + t83 * t99;
t108 = MDP(17) * t80;
t106 = t98 * qJD(5) * t116;
t105 = t82 * MDP(11);
t104 = t82 * MDP(17);
t85 = qJD(5) ^ 2;
t103 = t85 * MDP(13);
t102 = t85 * MDP(14);
t100 = MDP(17) * qJD(5);
t97 = t98 * t105;
t59 = -t81 * t95 + t110;
t65 = -t81 * t101 + t111;
t61 = -t65 + t114;
t96 = t61 * t98 - t59;
t94 = t61 + t65 + t114;
t93 = t98 * MDP(16);
t90 = -t81 * qJ(3) + t83 * t84;
t63 = t83 * qJD(3) + t90 * qJD(4);
t92 = -(pkin(4) - t90) * t98 - t61 - t63;
t89 = t80 * MDP(16) + t104;
t88 = pkin(7) * t85 - t118;
t87 = (-pkin(7) + t91) * t85 - t117;
t76 = t98 ^ 2;
t1 = [t89 * t85; t117 * MDP(9) + (t63 * t98 + t110) * MDP(10) - 0.2e1 * t106 + (MDP(6) * t115 + (-MDP(10) * t107 + MDP(7) * t115) * qJ(3)) * qJD(2) + (-t87 * MDP(16) + t92 * t100 - t103) * t82 + (t102 + t87 * MDP(17) + (t92 * MDP(16) + 0.2e1 * t97) * qJD(5)) * t80; (-qJ(3) * MDP(7) - MDP(6)) * qJD(2) ^ 2 + (-MDP(16) * t82 + t108) * t85 * t81 + (t98 * t104 + t80 * t93) * t83 * qJD(5) - ((t98 * MDP(10) - t89 * qJD(5)) * t83 + (t82 * t93 + (MDP(9) - t108) * t98) * t81) * t98; t118 * MDP(9) + (-t65 * t98 - t59) * MDP(10) + 0.2e1 * t106 + (-t88 * MDP(16) + t94 * t100 + t103) * t82 + (-t102 + t88 * MDP(17) + (t94 * MDP(16) - 0.2e1 * t97) * qJD(5)) * t80; t96 * t104 + t76 * t116 + (t96 * MDP(16) - t76 * t105) * t80;];
tauc = t1;
