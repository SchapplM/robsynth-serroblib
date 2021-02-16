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
% MDP [21x1]
%   Minimal dynamic parameter vector (fixed base model)
%   see S5RRRPP1_convert_par2_MPV_fixb.m
% 
% Output:
% Mq [5x5]
%   inertia matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2021-01-15 22:15
% Revision: 24b2e7d74a0c1a3b64fa2f8f5ad758691ad61af3 (2021-01-15)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Mq = S5RRRPP1_inertiaJ_mdp_slag_vp(qJ, pkin, MDP)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(8,1),zeros(21,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RRRPP1_inertiaJ_mdp_slag_vp: qJ has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S5RRRPP1_inertiaJ_mdp_slag_vp: pkin has to be [8x1] (double)');
assert(isreal(MDP) && all(size(MDP) == [21 1]), ...
  'S5RRRPP1_inertiaJ_mdp_slag_vp: MDP has to be [21x1] (double)'); 

%% Symbolic Calculation
% From inertia_joint_joint_fixb_mdp_matlab.m
% OptimizationMode: 2
% StartTime: 2021-01-15 22:14:50
% EndTime: 2021-01-15 22:14:53
% DurationCPUTime: 0.37s
% Computational Cost: add. (488->105), mult. (819->137), div. (0->0), fcn. (764->6), ass. (0->50)
t160 = 2 * MDP(16) + 2 * MDP(19);
t158 = MDP(15) - MDP(20);
t146 = MDP(17) * pkin(3);
t116 = sin(pkin(8));
t117 = cos(pkin(8));
t118 = sin(qJ(3));
t120 = cos(qJ(3));
t100 = t116 * t118 - t117 * t120;
t101 = t116 * t120 + t117 * t118;
t115 = t120 * qJ(4);
t104 = pkin(7) * t120 + t115;
t132 = (-qJ(4) - pkin(7)) * t118;
t90 = t104 * t116 - t117 * t132;
t92 = t117 * t104 + t116 * t132;
t157 = -t92 * t100 + t90 * t101;
t119 = sin(qJ(2));
t134 = pkin(1) * t119 + pkin(7);
t124 = (-qJ(4) - t134) * t118;
t128 = t120 * t134;
t99 = t128 + t115;
t85 = t116 * t99 - t117 * t124;
t87 = t116 * t124 + t117 * t99;
t156 = -t87 * t100 + t85 * t101;
t125 = MDP(12) * t120 - MDP(13) * t118;
t121 = cos(qJ(2));
t153 = pkin(1) * t121;
t111 = -pkin(3) * t120 - pkin(2);
t88 = pkin(4) * t100 - qJ(5) * t101 + t111;
t82 = t88 - t153;
t147 = t82 + t88;
t103 = t111 - t153;
t143 = t103 + t111;
t142 = MDP(12) * t118;
t139 = t85 ^ 2 + t87 ^ 2;
t138 = t90 ^ 2 + t92 ^ 2;
t137 = 0.2e1 * t100;
t136 = 0.2e1 * t101;
t135 = MDP(4) + (MDP(7) * t118 + 0.2e1 * MDP(8) * t120) * t118;
t133 = t85 * t90 + t87 * t92;
t131 = t158 * t101 + (MDP(14) + MDP(18)) * t100;
t106 = pkin(3) * t116 + qJ(5);
t108 = pkin(3) * t117 + pkin(4);
t130 = t120 * MDP(10) + t118 * MDP(9) + (-t100 * t106 - t101 * t108) * MDP(19) + (-t100 * t116 - t101 * t117) * pkin(3) * MDP(16);
t129 = -MDP(21) * t108 - MDP(18);
t127 = -MDP(14) + t129;
t126 = MDP(21) * t106 - t158;
t123 = (t121 * MDP(5) - t119 * MDP(6)) * pkin(1);
t110 = -pkin(2) - t153;
t96 = t101 * MDP(19);
t1 = [MDP(1) + (t103 ^ 2 + t139) * MDP(17) + (t82 ^ 2 + t139) * MDP(21) + (t103 * MDP(15) - t82 * MDP(20)) * t136 + (t103 * MDP(14) + t82 * MDP(18)) * t137 + t135 - 0.2e1 * t125 * t110 + 0.2e1 * t123 + t160 * t156; (t103 * t111 + t133) * MDP(17) + (t82 * t88 + t133) * MDP(21) + t123 + (t143 * MDP(15) - t147 * MDP(20)) * t101 + (t143 * MDP(14) + t147 * MDP(18)) * t100 + t135 + t125 * (pkin(2) - t110) + (MDP(16) + MDP(19)) * (t157 + t156); (t111 ^ 2 + t138) * MDP(17) + (t88 ^ 2 + t138) * MDP(21) + (MDP(15) * t111 - MDP(20) * t88) * t136 + (MDP(14) * t111 + MDP(18) * t88) * t137 + 0.2e1 * t125 * pkin(2) + t135 + t160 * t157; -t134 * t142 - MDP(13) * t128 + t126 * t87 + t127 * t85 + (t116 * t87 - t117 * t85) * t146 + t130; t126 * t92 + t127 * t90 + (-MDP(13) * t120 - t142) * pkin(7) + (t116 * t92 - t117 * t90) * t146 + t130; MDP(11) + (t106 ^ 2 + t108 ^ 2) * MDP(21) + 0.2e1 * MDP(18) * t108 + 0.2e1 * MDP(20) * t106 + (0.2e1 * MDP(14) * t117 - 0.2e1 * MDP(15) * t116 + (t116 ^ 2 + t117 ^ 2) * t146) * pkin(3); t103 * MDP(17) + t82 * MDP(21) + t131; t111 * MDP(17) + t88 * MDP(21) + t131; 0; MDP(17) + MDP(21); t85 * MDP(21) + t96; t90 * MDP(21) + t96; t129; 0; MDP(21);];
%% Postprocessing: Reshape Output
% From vec2symmat_5_matlab.m
res = [t1(1), t1(2), t1(4), t1(7), t1(11); t1(2), t1(3), t1(5), t1(8), t1(12); t1(4), t1(5), t1(6), t1(9), t1(13); t1(7), t1(8), t1(9), t1(10), t1(14); t1(11), t1(12), t1(13), t1(14), t1(15);];
Mq = res;
