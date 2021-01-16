% Calculate joint inertia matrix for
% S5RRPPR6
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% pkin [9x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d1,d2,d5,theta3,theta4]';
% MDP [25x1]
%   Minimal dynamic parameter vector (fixed base model)
%   see S5RRPPR6_convert_par2_MPV_fixb.m
% 
% Output:
% Mq [5x5]
%   inertia matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2021-01-15 19:48
% Revision: 24b2e7d74a0c1a3b64fa2f8f5ad758691ad61af3 (2021-01-15)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Mq = S5RRPPR6_inertiaJ_mdp_slag_vp(qJ, pkin, MDP)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(9,1),zeros(25,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RRPPR6_inertiaJ_mdp_slag_vp: qJ has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S5RRPPR6_inertiaJ_mdp_slag_vp: pkin has to be [9x1] (double)');
assert(isreal(MDP) && all(size(MDP) == [25 1]), ...
  'S5RRPPR6_inertiaJ_mdp_slag_vp: MDP has to be [25x1] (double)'); 

%% Symbolic Calculation
% From inertia_joint_joint_fixb_mdp_matlab.m
% OptimizationMode: 2
% StartTime: 2021-01-15 19:47:34
% EndTime: 2021-01-15 19:47:37
% DurationCPUTime: 0.39s
% Computational Cost: add. (527->110), mult. (999->168), div. (0->0), fcn. (1062->8), ass. (0->59)
t100 = cos(pkin(9));
t98 = sin(pkin(9));
t127 = t100 ^ 2 + t98 ^ 2;
t136 = t127 * MDP(17);
t101 = cos(pkin(8));
t103 = sin(qJ(2));
t128 = -qJ(3) - pkin(6);
t117 = t128 * t103;
t105 = cos(qJ(2));
t89 = t128 * t105;
t99 = sin(pkin(8));
t78 = -t101 * t117 - t89 * t99;
t135 = t78 ^ 2;
t94 = -pkin(2) * t101 - pkin(3);
t88 = -pkin(4) * t100 + t94;
t134 = 0.2e1 * t88;
t95 = -pkin(2) * t105 - pkin(1);
t133 = 0.2e1 * t95;
t132 = 0.2e1 * t105;
t131 = -2 * MDP(20);
t91 = pkin(2) * t99 + qJ(4);
t130 = pkin(7) + t91;
t85 = t101 * t103 + t105 * t99;
t129 = t85 * t98;
t83 = -t101 * t105 + t103 * t99;
t76 = pkin(3) * t83 - qJ(4) * t85 + t95;
t80 = -t101 * t89 + t99 * t117;
t68 = t100 * t80 + t98 * t76;
t126 = t100 * t85;
t125 = MDP(16) * t98;
t124 = MDP(18) * t94;
t102 = sin(qJ(5));
t104 = cos(qJ(5));
t86 = t100 * t102 + t104 * t98;
t123 = MDP(19) * t86;
t122 = MDP(23) * t83;
t113 = t104 * t100 - t102 * t98;
t121 = MDP(24) * t113;
t70 = t86 * t85;
t120 = t70 * MDP(22);
t71 = t113 * t85;
t119 = t71 * MDP(21);
t118 = MDP(15) * t100;
t67 = t100 * t76 - t80 * t98;
t116 = t127 * MDP(18);
t115 = t100 * t68 - t67 * t98;
t114 = t100 * t67 + t68 * t98;
t65 = pkin(4) * t83 - pkin(7) * t126 + t67;
t66 = -pkin(7) * t129 + t68;
t112 = (-t102 * t66 + t104 * t65) * MDP(24) - (t102 * t65 + t104 * t66) * MDP(25);
t111 = t70 * MDP(24) + t71 * MDP(25);
t110 = -MDP(25) * t86 + t121;
t109 = t98 * MDP(15) + t100 * MDP(16);
t81 = t130 * t98;
t82 = t130 * t100;
t108 = t86 * MDP(21) + t113 * MDP(22) + (-t102 * t82 - t104 * t81) * MDP(24) - (-t102 * t81 + t104 * t82) * MDP(25);
t107 = t110 + t118 - t125;
t69 = pkin(4) * t129 + t78;
t1 = [MDP(1) + pkin(1) * MDP(9) * t132 + t85 * MDP(12) * t133 + (t80 ^ 2 + t95 ^ 2 + t135) * MDP(14) + (t67 ^ 2 + t68 ^ 2 + t135) * MDP(18) + (MDP(19) * t71 + t70 * t131) * t71 + (-0.2e1 * pkin(1) * MDP(10) + MDP(4) * t103 + MDP(5) * t132) * t103 + (MDP(11) * t133 + 0.2e1 * t119 - 0.2e1 * t120 + t122) * t83 + 0.2e1 * t111 * t69 + 0.2e1 * (-t80 * MDP(13) + t67 * MDP(15) - t68 * MDP(16) + t112) * t83 + 0.2e1 * (-t114 * MDP(17) + (MDP(13) + t109) * t78) * t85; t103 * MDP(6) + t105 * MDP(7) - t78 * MDP(11) - t80 * MDP(12) + (-t100 * t78 + t94 * t129) * MDP(15) + (t94 * t126 + t78 * t98) * MDP(16) + t115 * MDP(17) + (t115 * t91 + t78 * t94) * MDP(18) + t71 * t123 + (t113 * t71 - t70 * t86) * MDP(20) + (-t113 * t69 + t70 * t88) * MDP(24) + (t69 * t86 + t71 * t88) * MDP(25) + (-t109 * t91 + t108) * t83 + (-MDP(10) * t105 - MDP(9) * t103) * pkin(6) + ((-t101 * t85 - t83 * t99) * MDP(13) + (-t101 * t78 + t80 * t99) * MDP(14)) * pkin(2); -t121 * t134 + MDP(8) + (-0.2e1 * t118 + t124 + 0.2e1 * t125) * t94 + (MDP(25) * t134 - t113 * t131 + t123) * t86 + (t116 * t91 + 0.2e1 * t136) * t91 + (0.2e1 * MDP(11) * t101 - 0.2e1 * MDP(12) * t99 + (t101 ^ 2 + t99 ^ 2) * MDP(14) * pkin(2)) * pkin(2); t95 * MDP(14) + t114 * MDP(18) + (MDP(12) - t136) * t85 + (MDP(11) + t107) * t83; 0; MDP(14) + t116; MDP(18) * t78 + t109 * t85 + t111; -t107 + t124; 0; MDP(18); t112 + t119 - t120 + t122; t108; t110; 0; MDP(23);];
%% Postprocessing: Reshape Output
% From vec2symmat_5_matlab.m
res = [t1(1), t1(2), t1(4), t1(7), t1(11); t1(2), t1(3), t1(5), t1(8), t1(12); t1(4), t1(5), t1(6), t1(9), t1(13); t1(7), t1(8), t1(9), t1(10), t1(14); t1(11), t1(12), t1(13), t1(14), t1(15);];
Mq = res;
