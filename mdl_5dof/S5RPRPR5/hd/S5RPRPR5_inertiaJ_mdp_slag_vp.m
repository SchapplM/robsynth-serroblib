% Calculate joint inertia matrix for
% S5RPRPR5
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% pkin [9x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d1,d3,d5,theta2,theta4]';
% MDP [24x1]
%   Minimal dynamic parameter vector (fixed base model)
%   see S5RPRPR5_convert_par2_MPV_fixb.m
% 
% Output:
% Mq [5x5]
%   inertia matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2021-01-15 11:56
% Revision: d12c3222fdeb2c5f3b3c8fa5751e113be2fc3aae (2021-01-15)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Mq = S5RPRPR5_inertiaJ_mdp_slag_vp(qJ, pkin, MDP)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(9,1),zeros(24,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RPRPR5_inertiaJ_mdp_slag_vp: qJ has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S5RPRPR5_inertiaJ_mdp_slag_vp: pkin has to be [9x1] (double)');
assert(isreal(MDP) && all(size(MDP) == [24 1]), ...
  'S5RPRPR5_inertiaJ_mdp_slag_vp: MDP has to be [24x1] (double)'); 

%% Symbolic Calculation
% From inertia_joint_joint_fixb_mdp_matlab.m
% OptimizationMode: 2
% StartTime: 2021-01-15 11:55:15
% EndTime: 2021-01-15 11:55:17
% DurationCPUTime: 0.43s
% Computational Cost: add. (479->93), mult. (978->143), div. (0->0), fcn. (986->8), ass. (0->46)
t103 = sin(qJ(5));
t105 = cos(qJ(5));
t102 = cos(pkin(8));
t101 = cos(pkin(9));
t106 = cos(qJ(3));
t100 = sin(pkin(8));
t126 = qJ(4) * t100;
t104 = sin(qJ(3));
t127 = qJ(2) * t104;
t92 = -pkin(2) * t102 - t100 * pkin(6) - pkin(1);
t88 = t106 * t92;
t75 = -t106 * t126 + t88 + (-pkin(3) - t127) * t102;
t117 = t106 * t102 * qJ(2);
t80 = t117 + (t92 - t126) * t104;
t99 = sin(pkin(9));
t65 = t101 * t75 - t99 * t80;
t124 = t101 * t106;
t125 = t100 * t104;
t85 = t100 * t124 - t99 * t125;
t63 = -t102 * pkin(4) - t85 * pkin(7) + t65;
t66 = t101 * t80 + t99 * t75;
t90 = t101 * t104 + t99 * t106;
t84 = t90 * t100;
t64 = -t84 * pkin(7) + t66;
t69 = t103 * t85 + t105 * t84;
t70 = -t103 * t84 + t105 * t85;
t139 = -(t103 * t63 + t105 * t64) * MDP(24) + (-t103 * t64 + t105 * t63) * MDP(23) + t70 * MDP(20) - t69 * MDP(21);
t138 = -(t104 * MDP(10) - t106 * MDP(9)) * t100 - t66 * MDP(15) + (-t102 * t127 + t88) * MDP(12) - (t104 * t92 + t117) * MDP(13) + t65 * MDP(14) + t139;
t134 = MDP(17) * pkin(3);
t89 = -t99 * t104 + t124;
t128 = (-t103 * t90 + t105 * t89) * MDP(23) - (t103 * t89 + t105 * t90) * MDP(24);
t132 = -t106 * MDP(12) + t104 * MDP(13) - t89 * MDP(14) + t90 * MDP(15) - t128;
t130 = pkin(3) * t99;
t91 = pkin(3) * t125 + t100 * qJ(2);
t95 = t101 * pkin(3) + pkin(4);
t120 = (-t103 * t130 + t105 * t95) * MDP(23);
t119 = (t103 * t95 + t105 * t130) * MDP(24);
t118 = MDP(11) + MDP(22);
t115 = t84 * MDP(14) + t85 * MDP(15);
t113 = t69 * MDP(23) + t70 * MDP(24);
t111 = t101 * MDP(14) - t99 * MDP(15);
t110 = MDP(22) - t119 + t120;
t107 = qJ(2) ^ 2;
t98 = t102 ^ 2;
t97 = t100 ^ 2;
t1 = [MDP(1) + (pkin(1) ^ 2 + t97 * t107) * MDP(6) + (t65 ^ 2 + t66 ^ 2 + t91 ^ 2) * MDP(17) + (t107 * MDP(6) + t118) * t98 + (MDP(18) * t70 - 0.2e1 * t69 * MDP(19)) * t70 + (MDP(7) * t106 - 0.2e1 * MDP(8) * t104) * t97 * t106 + 0.2e1 * (-t65 * t85 - t66 * t84) * MDP(16) + 0.2e1 * t115 * t91 + 0.2e1 * t113 * (t84 * pkin(4) + t91) + 0.2e1 * (t98 * MDP(5) + (t104 * MDP(12) + t106 * MDP(13) + MDP(5)) * t97) * qJ(2) + 0.2e1 * (pkin(1) * MDP(4) - t138) * t102; -pkin(1) * MDP(6) + (-t90 * t84 - t89 * t85) * MDP(16) + (t65 * t89 + t66 * t90) * MDP(17) + (-MDP(4) + t132) * t102; MDP(6) + (t89 ^ 2 + t90 ^ 2) * MDP(17); (-MDP(11) - t110) * t102 + ((-t101 * t85 - t84 * t99) * MDP(16) + (t101 * t65 + t66 * t99) * MDP(17) - t111 * t102) * pkin(3) + t138; (t101 * t89 + t90 * t99) * t134 - t132; 0.2e1 * t120 - 0.2e1 * t119 + t118 + (0.2e1 * t111 + (t101 ^ 2 + t99 ^ 2) * t134) * pkin(3); t91 * MDP(17) + t113 + t115; 0; 0; MDP(17); -t102 * MDP(22) + t139; t128; t110; 0; MDP(22);];
%% Postprocessing: Reshape Output
% From vec2symmat_5_matlab.m
res = [t1(1), t1(2), t1(4), t1(7), t1(11); t1(2), t1(3), t1(5), t1(8), t1(12); t1(4), t1(5), t1(6), t1(9), t1(13); t1(7), t1(8), t1(9), t1(10), t1(14); t1(11), t1(12), t1(13), t1(14), t1(15);];
Mq = res;
