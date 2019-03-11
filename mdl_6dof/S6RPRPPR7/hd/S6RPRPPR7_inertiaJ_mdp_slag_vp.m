% Calculate joint inertia matrix for
% S6RPRPPR7
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% pkin [9x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d3,d6,theta4]';
% MDP [26x1]
%   Minimal dynamic parameter vector (fixed base model)
%   see S6RPRPPR7_convert_par2_MPV_fixb.m
% 
% Output:
% Mq [6x6]
%   inertia matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-09 02:57
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Mq = S6RPRPPR7_inertiaJ_mdp_slag_vp(qJ, pkin, MDP)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(9,1),zeros(26,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPRPPR7_inertiaJ_mdp_slag_vp: qJ has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S6RPRPPR7_inertiaJ_mdp_slag_vp: pkin has to be [9x1] (double)');
assert(isreal(MDP) && all(size(MDP) == [26 1]), ...
  'S6RPRPPR7_inertiaJ_mdp_slag_vp: MDP has to be [26x1] (double)'); 

%% Symbolic Calculation
% From inertia_joint_joint_fixb_mdp_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 02:57:27
% EndTime: 2019-03-09 02:57:28
% DurationCPUTime: 0.33s
% Computational Cost: add. (411->103), mult. (657->139), div. (0->0), fcn. (664->6), ass. (0->52)
t98 = sin(qJ(6));
t99 = cos(qJ(6));
t108 = MDP(25) * t99 - t98 * MDP(26);
t129 = MDP(15) + MDP(19);
t100 = cos(qJ(3));
t127 = sin(qJ(3));
t96 = sin(pkin(9));
t97 = cos(pkin(9));
t80 = t100 * t97 - t96 * t127;
t128 = t80 ^ 2;
t81 = -t96 * t100 - t97 * t127;
t78 = t81 ^ 2;
t132 = t78 + t128;
t101 = -pkin(1) - pkin(7);
t131 = -qJ(4) + t101;
t130 = MDP(14) + MDP(16);
t126 = (pkin(1) * MDP(6));
t112 = t131 * t100;
t84 = t131 * t127;
t74 = t96 * t112 + t97 * t84;
t70 = t81 * pkin(5) + t74;
t125 = t70 * t81;
t124 = t98 * t99;
t123 = MDP(15) * pkin(3);
t92 = t127 * pkin(3) + qJ(2);
t122 = MDP(17) * t81;
t86 = pkin(3) * t96 + qJ(5);
t121 = MDP(19) * t86;
t91 = -pkin(3) * t97 - pkin(4);
t120 = MDP(19) * t91;
t119 = MDP(25) * t98;
t117 = MDP(26) * t99;
t115 = MDP(21) * t124;
t72 = -t97 * t112 + t84 * t96;
t114 = t72 ^ 2 + t74 ^ 2;
t113 = t127 * MDP(13);
t111 = MDP(17) + t120;
t110 = -qJ(5) * t80 + t92;
t109 = t72 * t80 + t74 * t81;
t107 = t117 + t119;
t106 = MDP(16) + t108;
t105 = -MDP(18) - t107;
t104 = (-MDP(22) * t98 - MDP(23) * t99) * t81;
t103 = t99 * MDP(22) - t98 * MDP(23) + t108 * (-pkin(8) + t91);
t95 = t99 ^ 2;
t94 = t98 ^ 2;
t71 = -pkin(4) * t81 + t110;
t69 = pkin(5) * t80 + t72;
t68 = (-pkin(4) - pkin(8)) * t81 + t110;
t66 = t68 * t99 + t69 * t98;
t65 = -t68 * t98 + t69 * t99;
t1 = [MDP(1) + (t92 ^ 2 + t114) * MDP(15) + 0.2e1 * t71 * t122 + (t71 ^ 2 + t114) * MDP(19) + t128 * MDP(24) + (MDP(20) * t94 + 0.2e1 * t115) * t78 + (MDP(7) * t100 - 0.2e1 * t127 * MDP(8)) * t100 + ((-2 * MDP(4) + t126) * pkin(1)) + (0.2e1 * t127 * MDP(12) + 0.2e1 * t100 * MDP(13) + MDP(6) * qJ(2) + (2 * MDP(5))) * qJ(2) + 0.2e1 * (t99 * t125 + t65 * t80) * MDP(25) + 0.2e1 * (-t98 * t125 - t66 * t80) * MDP(26) + 0.2e1 * (-t71 * MDP(18) + t104) * t80 + 0.2e1 * t130 * t109; -t129 * t109 + MDP(4) - t126 + (-t130 - t108) * t132; t129 * t132 + MDP(6); -t127 * MDP(10) - t101 * t113 + t72 * MDP(17) + t74 * MDP(18) + (t72 * t91 + t74 * t86) * MDP(19) + t107 * t70 + (MDP(12) * t101 + MDP(9)) * t100 + (t91 * MDP(16) + t103) * t80 + (-MDP(20) * t124 + (t94 - t95) * MDP(21) + t106 * t86) * t81 + ((-t80 * t97 + t81 * t96) * MDP(14) + (-t72 * t97 + t74 * t96) * MDP(15)) * pkin(3); t100 * MDP(12) - t113 + (t97 * t123 - t111) * t80 + (-t96 * t123 + t105 - t121) * t81; -0.2e1 * t115 + t95 * MDP(20) + MDP(11) + (0.2e1 * MDP(17) + t120) * t91 + (t96 ^ 2 + t97 ^ 2) * MDP(15) * pkin(3) ^ 2 + (0.2e1 * MDP(18) + 0.2e1 * t117 + 0.2e1 * t119 + t121) * t86; MDP(15) * t92 + MDP(19) * t71 + t105 * t80 + t122; 0; 0; t129; MDP(19) * t72 + t106 * t80; -t80 * MDP(19); t111; 0; MDP(19); MDP(24) * t80 + MDP(25) * t65 - MDP(26) * t66 + t104; -t108 * t80; t103; -t107; t108; MDP(24);];
%% Postprocessing: Reshape Output
% From vec2symmat_6_matlab.m
res = [t1(1) t1(2) t1(4) t1(7) t1(11) t1(16); t1(2) t1(3) t1(5) t1(8) t1(12) t1(17); t1(4) t1(5) t1(6) t1(9) t1(13) t1(18); t1(7) t1(8) t1(9) t1(10) t1(14) t1(19); t1(11) t1(12) t1(13) t1(14) t1(15) t1(20); t1(16) t1(17) t1(18) t1(19) t1(20) t1(21);];
Mq  = res;
