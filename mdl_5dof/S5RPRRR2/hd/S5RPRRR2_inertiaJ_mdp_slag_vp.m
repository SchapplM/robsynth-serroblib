% Calculate joint inertia matrix for
% S5RPRRR2
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% pkin [9x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d1,d3,d4,d5,theta2]';
% MDP [28x1]
%   Minimal dynamic parameter vector (fixed base model)
%   see S5RPRRR2_convert_par2_MPV_fixb.m
% 
% Output:
% Mq [5x5]
%   inertia matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-05 18:13
% Revision: 77da58f92bca3eff71542919beafa37024070d86 (2019-12-05)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Mq = S5RPRRR2_inertiaJ_mdp_slag_vp(qJ, pkin, MDP)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(9,1),zeros(28,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RPRRR2_inertiaJ_mdp_slag_vp: qJ has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S5RPRRR2_inertiaJ_mdp_slag_vp: pkin has to be [9x1] (double)');
assert(isreal(MDP) && all(size(MDP) == [28 1]), ...
  'S5RPRRR2_inertiaJ_mdp_slag_vp: MDP has to be [28x1] (double)'); 

%% Symbolic Calculation
% From inertia_joint_joint_fixb_mdp_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-05 18:12:15
% EndTime: 2019-12-05 18:12:16
% DurationCPUTime: 0.17s
% Computational Cost: add. (471->79), mult. (869->110), div. (0->0), fcn. (1012->8), ass. (0->48)
t86 = sin(pkin(9));
t87 = cos(pkin(9));
t90 = sin(qJ(3));
t93 = cos(qJ(3));
t76 = t90 * t86 - t93 * t87;
t77 = t93 * t86 + t90 * t87;
t89 = sin(qJ(4));
t92 = cos(qJ(4));
t69 = t92 * t76 + t89 * t77;
t82 = -t87 * pkin(2) - pkin(1);
t71 = t76 * pkin(3) + t82;
t118 = 0.2e1 * t69 * pkin(4) + 0.2e1 * t71;
t117 = 0.2e1 * t71;
t116 = 0.2e1 * t82;
t115 = pkin(3) * t89;
t114 = pkin(1) * MDP(7);
t113 = pkin(6) + qJ(2);
t111 = t86 * MDP(5);
t110 = t87 * MDP(4);
t70 = -t89 * t76 + t92 * t77;
t88 = sin(qJ(5));
t91 = cos(qJ(5));
t60 = t91 * t69 + t88 * t70;
t109 = t60 * MDP(27);
t108 = t69 * MDP(20);
t83 = t92 * pkin(3) + pkin(4);
t80 = t91 * t83;
t107 = (-t88 * t115 + t80) * MDP(27);
t106 = (-t91 * t115 - t88 * t83) * MDP(28);
t105 = t76 * MDP(13);
t104 = t88 * MDP(28);
t103 = t92 * MDP(20);
t102 = MDP(19) + MDP(26);
t78 = t113 * t86;
t79 = t113 * t87;
t100 = -t93 * t78 - t90 * t79;
t64 = -t77 * pkin(7) + t100;
t97 = t90 * t78 - t93 * t79;
t65 = -t76 * pkin(7) - t97;
t101 = t92 * t64 - t89 * t65;
t54 = -t70 * pkin(8) + t101;
t98 = -t89 * t64 - t92 * t65;
t55 = -t69 * pkin(8) - t98;
t61 = -t88 * t69 + t91 * t70;
t99 = t61 * MDP(24) - t60 * MDP(25) + (t91 * t54 - t88 * t55) * MDP(27) + (-t88 * t54 - t91 * t55) * MDP(28);
t96 = (t91 * MDP(27) - t104) * pkin(4);
t95 = t70 * MDP(17) - t69 * MDP(18) + t101 * MDP(20) + t98 * MDP(21) + t99;
t1 = [t105 * t116 + t108 * t117 + t109 * t118 + MDP(1) + (0.2e1 * t110 - 0.2e1 * t111 + t114) * pkin(1) + (MDP(14) * t116 + MDP(8) * t77 - 0.2e1 * t76 * MDP(9)) * t77 + (MDP(15) * t70 - 0.2e1 * t69 * MDP(16) + MDP(21) * t117) * t70 + (MDP(22) * t61 - 0.2e1 * t60 * MDP(23) + MDP(28) * t118) * t61 + (MDP(7) * qJ(2) + (2 * MDP(6))) * (t86 ^ 2 + t87 ^ 2) * qJ(2); t77 * MDP(14) + t70 * MDP(21) + t61 * MDP(28) + t105 + t108 + t109 - t110 + t111 - t114; MDP(7); t77 * MDP(10) - t76 * MDP(11) + t100 * MDP(13) + t97 * MDP(14) + t95; 0; MDP(12) + 0.2e1 * (-t89 * MDP(21) + t103) * pkin(3) + 0.2e1 * t107 + 0.2e1 * t106 + t102; t95; 0; (t91 * pkin(4) + t80) * MDP(27) + (-pkin(4) - t83) * t104 + (t103 + (-MDP(27) * t88 - MDP(28) * t91 - MDP(21)) * t89) * pkin(3) + t102; t102 + 0.2e1 * t96; t99; 0; MDP(26) + t106 + t107; MDP(26) + t96; MDP(26);];
%% Postprocessing: Reshape Output
% From vec2symmat_5_matlab.m
res = [t1(1), t1(2), t1(4), t1(7), t1(11); t1(2), t1(3), t1(5), t1(8), t1(12); t1(4), t1(5), t1(6), t1(9), t1(13); t1(7), t1(8), t1(9), t1(10), t1(14); t1(11), t1(12), t1(13), t1(14), t1(15);];
Mq = res;
