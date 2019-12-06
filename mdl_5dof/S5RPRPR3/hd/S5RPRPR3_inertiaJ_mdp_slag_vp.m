% Calculate joint inertia matrix for
% S5RPRPR3
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% pkin [9x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d1,d3,d5,theta2,theta4]';
% MDP [18x1]
%   Minimal dynamic parameter vector (fixed base model)
%   see S5RPRPR3_convert_par2_MPV_fixb.m
% 
% Output:
% Mq [5x5]
%   inertia matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-05 17:52
% Revision: 77da58f92bca3eff71542919beafa37024070d86 (2019-12-05)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Mq = S5RPRPR3_inertiaJ_mdp_slag_vp(qJ, pkin, MDP)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(9,1),zeros(18,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RPRPR3_inertiaJ_mdp_slag_vp: qJ has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S5RPRPR3_inertiaJ_mdp_slag_vp: pkin has to be [9x1] (double)');
assert(isreal(MDP) && all(size(MDP) == [18 1]), ...
  'S5RPRPR3_inertiaJ_mdp_slag_vp: MDP has to be [18x1] (double)'); 

%% Symbolic Calculation
% From inertia_joint_joint_fixb_mdp_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-05 17:51:34
% EndTime: 2019-12-05 17:51:34
% DurationCPUTime: 0.17s
% Computational Cost: add. (220->70), mult. (420->98), div. (0->0), fcn. (314->8), ass. (0->53)
t92 = sin(qJ(5));
t94 = cos(qJ(5));
t124 = t94 * MDP(17) - t92 * MDP(18);
t90 = cos(pkin(9));
t123 = 0.2e1 * t90;
t122 = 2 * MDP(10);
t121 = 0.2e1 * MDP(17);
t120 = 0.2e1 * MDP(18);
t89 = sin(pkin(8));
t119 = pkin(1) * t89;
t91 = cos(pkin(8));
t81 = t91 * pkin(1) + pkin(2);
t93 = sin(qJ(3));
t95 = cos(qJ(3));
t69 = -t93 * t119 + t95 * t81;
t68 = -pkin(3) - t69;
t118 = pkin(3) - t68;
t70 = -t95 * t119 - t93 * t81;
t67 = qJ(4) - t70;
t117 = t67 * t90;
t88 = sin(pkin(9));
t86 = t88 ^ 2;
t64 = t86 * t67;
t116 = t86 * t94;
t74 = -t90 * pkin(4) - t88 * pkin(7) - pkin(3);
t58 = -t69 + t74;
t57 = t94 * t117 + t92 * t58;
t115 = t67 * t116 + t57 * t90;
t109 = qJ(4) * t90;
t63 = t94 * t109 + t92 * t74;
t83 = t86 * qJ(4);
t114 = t63 * t90 + t94 * t83;
t87 = t90 ^ 2;
t113 = t87 * t67 + t64;
t112 = t87 * qJ(4) + t83;
t111 = t86 + t87;
t110 = pkin(3) * MDP(11);
t108 = t69 * MDP(6);
t107 = t70 * MDP(7);
t84 = t88 * MDP(9);
t106 = t68 * MDP(11);
t103 = t92 * t88 * MDP(15);
t77 = t94 * t88 * MDP(14);
t102 = t111 * MDP(11);
t101 = MDP(8) * t123 - 0.2e1 * t84;
t100 = t94 ^ 2 * t86 * MDP(12) - 0.2e1 * t92 * MDP(13) * t116 + t87 * MDP(16) + t103 * t123 - 0.2e1 * t90 * t77 + MDP(5);
t99 = -t90 * MDP(16) - t103 + t77;
t98 = t84 + (-MDP(8) - t124) * t90;
t78 = t92 * t83;
t62 = -t92 * t109 + t94 * t74;
t60 = t92 * t64;
t56 = -t92 * t117 + t94 * t58;
t1 = [MDP(1) + (t89 ^ 2 + t91 ^ 2) * MDP(4) * pkin(1) ^ 2 + t67 ^ 2 * t102 + (-t101 + t106) * t68 + 0.2e1 * t108 + 0.2e1 * t107 + t113 * t122 + (-t56 * t90 + t60) * t121 + t115 * t120 + t100; 0; MDP(4) + t102; t108 + t107 + (t112 + t113) * MDP(10) + (t111 * t67 * qJ(4) - t68 * pkin(3)) * MDP(11) + (t60 + t78) * MDP(17) + (t114 + t115) * MDP(18) - t118 * t84 + (t118 * MDP(8) + (-t56 - t62) * MDP(17)) * t90 + t100; 0; qJ(4) ^ 2 * t102 + (t101 + t110) * pkin(3) + t112 * t122 + (-t62 * t90 + t78) * t121 + t114 * t120 + t100; t98 + t106; 0; t98 - t110; MDP(11); t56 * MDP(17) - t57 * MDP(18) + t99; (-MDP(17) * t92 - MDP(18) * t94) * t88; t62 * MDP(17) - t63 * MDP(18) + t99; t124; MDP(16);];
%% Postprocessing: Reshape Output
% From vec2symmat_5_matlab.m
res = [t1(1), t1(2), t1(4), t1(7), t1(11); t1(2), t1(3), t1(5), t1(8), t1(12); t1(4), t1(5), t1(6), t1(9), t1(13); t1(7), t1(8), t1(9), t1(10), t1(14); t1(11), t1(12), t1(13), t1(14), t1(15);];
Mq = res;
