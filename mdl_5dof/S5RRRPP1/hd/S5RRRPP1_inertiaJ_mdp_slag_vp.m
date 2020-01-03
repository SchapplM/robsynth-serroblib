% Calculate joint inertia matrix for
% S5RRRPP1
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% pkin [8x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d1,d2,d3,theta4]';
% MDP [19x1]
%   Minimal dynamic parameter vector (fixed base model)
%   see S5RRRPP1_convert_par2_MPV_fixb.m
% 
% Output:
% Mq [5x5]
%   inertia matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-31 20:50
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Mq = S5RRRPP1_inertiaJ_mdp_slag_vp(qJ, pkin, MDP)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(8,1),zeros(19,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RRRPP1_inertiaJ_mdp_slag_vp: qJ has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S5RRRPP1_inertiaJ_mdp_slag_vp: pkin has to be [8x1] (double)');
assert(isreal(MDP) && all(size(MDP) == [19 1]), ...
  'S5RRRPP1_inertiaJ_mdp_slag_vp: MDP has to be [19x1] (double)'); 

%% Symbolic Calculation
% From inertia_joint_joint_fixb_mdp_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 20:49:58
% EndTime: 2019-12-31 20:49:59
% DurationCPUTime: 0.31s
% Computational Cost: add. (430->93), mult. (723->121), div. (0->0), fcn. (674->6), ass. (0->45)
t149 = 2 * MDP(14) + 2 * MDP(17);
t109 = sin(pkin(8));
t110 = cos(pkin(8));
t111 = sin(qJ(3));
t113 = cos(qJ(3));
t93 = t109 * t111 - t110 * t113;
t94 = t109 * t113 + t110 * t111;
t134 = t93 * MDP(16) - t94 * MDP(18);
t112 = sin(qJ(2));
t124 = t112 * pkin(1) + pkin(7);
t117 = (-qJ(4) - t124) * t111;
t108 = t113 * qJ(4);
t120 = t113 * t124;
t92 = t120 + t108;
t80 = t109 * t92 - t110 * t117;
t82 = t109 * t117 + t110 * t92;
t147 = t80 * t94 - t82 * t93;
t122 = (-qJ(4) - pkin(7)) * t111;
t97 = t113 * pkin(7) + t108;
t85 = t109 * t97 - t110 * t122;
t87 = t109 * t122 + t110 * t97;
t146 = t85 * t94 - t87 * t93;
t118 = t113 * MDP(12) - t111 * MDP(13);
t114 = cos(qJ(2));
t143 = t114 * pkin(1);
t133 = MDP(15) * pkin(3);
t104 = -t113 * pkin(3) - pkin(2);
t83 = t93 * pkin(4) - t94 * qJ(5) + t104;
t77 = t83 - t143;
t132 = t77 * MDP(19);
t131 = t83 * MDP(19);
t130 = t111 * MDP(12);
t127 = t80 ^ 2 + t82 ^ 2;
t126 = t85 ^ 2 + t87 ^ 2;
t125 = MDP(4) + (MDP(7) * t111 + 0.2e1 * MDP(8) * t113) * t111;
t123 = t80 * t85 + t82 * t87;
t101 = t110 * pkin(3) + pkin(4);
t99 = t109 * pkin(3) + qJ(5);
t121 = t113 * MDP(10) + t111 * MDP(9) + (-t101 * t94 - t99 * t93) * MDP(17) + (-t109 * t93 - t110 * t94) * pkin(3) * MDP(14);
t119 = 0.2e1 * t134;
t116 = (t114 * MDP(5) - t112 * MDP(6)) * pkin(1);
t103 = -pkin(2) - t143;
t96 = t104 - t143;
t90 = t94 * MDP(17);
t1 = [MDP(1) + (t96 ^ 2 + t127) * MDP(15) + t127 * MDP(19) + (t119 + t132) * t77 + t125 - 0.2e1 * t118 * t103 + 0.2e1 * t116 + t149 * t147; (t96 * t104 + t123) * MDP(15) + (t77 * t83 + t123) * MDP(19) + t116 + t125 + t118 * (pkin(2) - t103) + t134 * (t77 + t83) + (MDP(14) + MDP(17)) * (t146 + t147); (t104 ^ 2 + t126) * MDP(15) + t126 * MDP(19) + (t119 + t131) * t83 + 0.2e1 * t118 * pkin(2) + t125 + t149 * t146; -t124 * t130 - MDP(13) * t120 - t80 * MDP(16) + t82 * MDP(18) + (-t80 * t101 + t82 * t99) * MDP(19) + (t109 * t82 - t110 * t80) * t133 + t121; -t85 * MDP(16) + t87 * MDP(18) + (-t85 * t101 + t87 * t99) * MDP(19) + (-t113 * MDP(13) - t130) * pkin(7) + (t109 * t87 - t110 * t85) * t133 + t121; MDP(11) + (t101 ^ 2 + t99 ^ 2) * MDP(19) + (t109 ^ 2 + t110 ^ 2) * MDP(15) * pkin(3) ^ 2 + 0.2e1 * t101 * MDP(16) + 0.2e1 * t99 * MDP(18); t96 * MDP(15) + t132 + t134; t104 * MDP(15) + t131 + t134; 0; MDP(15) + MDP(19); t80 * MDP(19) + t90; t85 * MDP(19) + t90; -t101 * MDP(19) - MDP(16); 0; MDP(19);];
%% Postprocessing: Reshape Output
% From vec2symmat_5_matlab.m
res = [t1(1), t1(2), t1(4), t1(7), t1(11); t1(2), t1(3), t1(5), t1(8), t1(12); t1(4), t1(5), t1(6), t1(9), t1(13); t1(7), t1(8), t1(9), t1(10), t1(14); t1(11), t1(12), t1(13), t1(14), t1(15);];
Mq = res;
