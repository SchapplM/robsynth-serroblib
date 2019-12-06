% Calculate joint inertia matrix for
% S5RRPPR1
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% pkin [9x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d1,d2,d5,theta3,theta4]';
% MDP [18x1]
%   Minimal dynamic parameter vector (fixed base model)
%   see S5RRPPR1_convert_par2_MPV_fixb.m
% 
% Output:
% Mq [5x5]
%   inertia matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-05 18:19
% Revision: 77da58f92bca3eff71542919beafa37024070d86 (2019-12-05)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Mq = S5RRPPR1_inertiaJ_mdp_slag_vp(qJ, pkin, MDP)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(9,1),zeros(18,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RRPPR1_inertiaJ_mdp_slag_vp: qJ has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S5RRPPR1_inertiaJ_mdp_slag_vp: pkin has to be [9x1] (double)');
assert(isreal(MDP) && all(size(MDP) == [18 1]), ...
  'S5RRPPR1_inertiaJ_mdp_slag_vp: MDP has to be [18x1] (double)'); 

%% Symbolic Calculation
% From inertia_joint_joint_fixb_mdp_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-05 18:18:28
% EndTime: 2019-12-05 18:18:29
% DurationCPUTime: 0.23s
% Computational Cost: add. (214->62), mult. (363->82), div. (0->0), fcn. (307->8), ass. (0->42)
t87 = sin(pkin(9));
t89 = cos(pkin(9));
t106 = t87 ^ 2 + t89 ^ 2;
t88 = sin(pkin(8));
t79 = t88 * pkin(2) + qJ(4);
t107 = t106 * t79;
t91 = sin(qJ(5));
t93 = cos(qJ(5));
t70 = t91 * t87 - t93 * t89;
t71 = t93 * t87 + t91 * t89;
t98 = t70 * MDP(17) + t71 * MDP(18);
t115 = -t89 * MDP(8) + t87 * MDP(9);
t114 = 2 * MDP(10);
t92 = sin(qJ(2));
t113 = pkin(1) * t92;
t112 = t89 * pkin(4);
t94 = cos(qJ(2));
t82 = t94 * pkin(1) + pkin(2);
t90 = cos(pkin(8));
t62 = t90 * t113 + t88 * t82;
t59 = qJ(4) + t62;
t110 = t106 * t59;
t108 = t71 * MDP(14) - t70 * MDP(15);
t61 = -t88 * t113 + t90 * t82;
t60 = -pkin(3) - t61;
t104 = t60 * MDP(11);
t81 = -t90 * pkin(2) - pkin(3);
t103 = t81 * MDP(11);
t102 = MDP(4) + (MDP(12) * t71 - 0.2e1 * MDP(13) * t70) * t71;
t101 = t106 * MDP(11);
t100 = 0.2e1 * t115;
t99 = t115 + t98;
t97 = (t94 * MDP(5) - t92 * MDP(6)) * pkin(1);
t96 = 0.2e1 * t98;
t84 = t89 * pkin(7);
t72 = t81 - t112;
t64 = t89 * t79 + t84;
t63 = (-pkin(7) - t79) * t87;
t55 = t60 - t112;
t54 = t89 * t59 + t84;
t53 = (-pkin(7) - t59) * t87;
t1 = [MDP(1) + (t61 ^ 2 + t62 ^ 2) * MDP(7) + t110 * t114 + t59 ^ 2 * t101 + (t100 + t104) * t60 + t55 * t96 + 0.2e1 * t97 + t102; (t107 + t110) * MDP(10) + (t59 * t107 + t60 * t81) * MDP(11) + (t61 * t90 + t62 * t88) * MDP(7) * pkin(2) + t97 + t102 + t98 * (t55 + t72) + t115 * (t60 + t81); t107 * t114 + (t88 ^ 2 + t90 ^ 2) * MDP(7) * pkin(2) ^ 2 + t79 ^ 2 * t101 + (t100 + t103) * t81 + t72 * t96 + t102; 0; 0; MDP(7) + t101; t99 + t104; t99 + t103; 0; MDP(11); (t93 * t53 - t91 * t54) * MDP(17) + (-t91 * t53 - t93 * t54) * MDP(18) + t108; (t93 * t63 - t91 * t64) * MDP(17) + (-t91 * t63 - t93 * t64) * MDP(18) + t108; -t98; 0; MDP(16);];
%% Postprocessing: Reshape Output
% From vec2symmat_5_matlab.m
res = [t1(1), t1(2), t1(4), t1(7), t1(11); t1(2), t1(3), t1(5), t1(8), t1(12); t1(4), t1(5), t1(6), t1(9), t1(13); t1(7), t1(8), t1(9), t1(10), t1(14); t1(11), t1(12), t1(13), t1(14), t1(15);];
Mq = res;
