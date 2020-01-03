% Calculate minimal parameter regressor of Coriolis joint torque vector for
% S4RPPR6
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [4x1]
%   Generalized joint coordinates (joint angles)
% qJD [4x1]
%   Generalized joint velocities
% pkin [6x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,d1,d4,theta2]';
% MDP [18x1]
%   Minimal dynamic parameter vector (fixed base model)
%   see S4RPPR6_convert_par2_MPV_fixb.m
% 
% Output:
% tauc [4x1]
%   joint torques required to compensate Coriolis and centrifugal load

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-31 16:40
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauc = S4RPPR6_coriolisvecJ_fixb_mdp_slag_vp(qJ, qJD, pkin, MDP)
%% Coder Information
%#codegen
%$cgargs {zeros(4,1),zeros(4,1),zeros(6,1),zeros(18,1)}
assert(isreal(qJ) && all(size(qJ) == [4 1]), ...
  'S4RPPR6_coriolisvecJ_fixb_mdp_slag_vp: qJ has to be [4x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [4 1]), ...
  'S4RPPR6_coriolisvecJ_fixb_mdp_slag_vp: qJD has to be [4x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [6 1]), ...
  'S4RPPR6_coriolisvecJ_fixb_mdp_slag_vp: pkin has to be [6x1] (double)');
assert(isreal(MDP) && all(size(MDP) == [18 1]), ...
  'S4RPPR6_coriolisvecJ_fixb_mdp_slag_vp: MDP has to be [18x1] (double)'); 

%% Symbolic Calculation
% From coriolisvec_joint_fixb_mdp_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 16:40:45
% EndTime: 2019-12-31 16:40:46
% DurationCPUTime: 0.34s
% Computational Cost: add. (162->71), mult. (454->119), div. (0->0), fcn. (270->4), ass. (0->43)
t78 = sin(pkin(6));
t76 = t78 ^ 2;
t79 = cos(pkin(6));
t77 = t79 ^ 2;
t100 = t76 + t77;
t105 = (MDP(6) + MDP(9)) * t100;
t81 = cos(qJ(4));
t103 = t78 * t81;
t80 = sin(qJ(4));
t102 = t79 * t80;
t101 = -pkin(5) + qJ(2);
t99 = t79 * MDP(8);
t64 = t78 * t80 + t79 * t81;
t98 = qJD(1) * t64;
t97 = qJD(1) * t79;
t96 = qJD(3) * t78;
t95 = t76 * MDP(10);
t94 = t78 * qJD(1);
t93 = qJD(1) * qJD(2);
t92 = t81 * t94;
t91 = qJD(4) * t103;
t90 = t80 * t97;
t88 = qJ(3) * t78 + pkin(1);
t87 = qJ(2) * t93;
t65 = -t102 + t103;
t86 = -pkin(2) * t79 - t88;
t85 = t65 * qJD(2);
t84 = t64 * qJD(2);
t56 = t64 * qJD(4);
t61 = (pkin(2) + pkin(3)) * t79 + t88;
t83 = qJD(1) ^ 2;
t72 = t76 * t87;
t71 = qJ(2) * t94 + qJD(3);
t70 = qJD(4) * t90;
t69 = t101 * t79;
t68 = t101 * t78;
t60 = -t90 + t92;
t57 = -qJD(4) * t102 + t91;
t55 = t86 * qJD(1) + qJD(2);
t54 = qJD(1) * t91 - t70;
t53 = qJD(1) * t56;
t52 = t61 * qJD(1) - qJD(2);
t1 = [0.2e1 * (t77 * t87 + t72) * MDP(7) + (t72 + (t71 * qJD(2) - t55 * qJD(3)) * t78 + (0.2e1 * t77 * qJ(2) * qJD(2) - t86 * t96) * qJD(1)) * MDP(11) + (-t53 * t65 - t56 * t60) * MDP(12) + (t53 * t64 - t54 * t65 + t56 * t98 - t57 * t60) * MDP(13) - t56 * qJD(4) * MDP(14) - t57 * qJD(4) * MDP(15) + (t52 * t57 + t61 * t54 + ((-t68 * t80 - t69 * t81) * qJD(4) + t85) * qJD(4) + 0.2e1 * t98 * t96) * MDP(17) + (-t52 * t56 - t61 * t53 + ((-t68 * t81 + t69 * t80) * qJD(4) - t84) * qJD(4) + (qJD(1) * t65 + t60) * t96) * MDP(18) + 0.2e1 * (t99 * t78 + t95) * qJD(1) * qJD(3) + 0.2e1 * t93 * t105; t70 * MDP(17) + (-qJD(3) - t71) * MDP(11) * t94 + ((-t60 - t92) * MDP(17) + (t80 * t94 + t81 * t97 + t98) * MDP(18)) * qJD(4) + (-t105 + (-t77 * MDP(11) - t100 * MDP(7)) * qJ(2)) * t83; -t83 * t95 + (-t80 * MDP(17) - t81 * MDP(18)) * qJD(4) ^ 2 + (-t83 * t99 + ((qJD(2) + t55) * MDP(11) - t98 * MDP(17) - t60 * MDP(18)) * qJD(1)) * t78; t60 * t98 * MDP(12) + (t60 ^ 2 - t98 ^ 2) * MDP(13) + (t70 + (t60 - t92) * qJD(4)) * MDP(15) + (qJD(1) * t85 - t52 * t60) * MDP(17) + (-qJD(1) * t84 + t52 * t98) * MDP(18);];
tauc = t1;
