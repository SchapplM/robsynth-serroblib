% Calculate joint inertia matrix for
% S6PRRPPR3
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% pkin [10x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,alpha2,d2,d3,d6,theta1]';
% MDP [26x1]
%   Minimal dynamic parameter vector (fixed base model)
%   see S6PRRPPR3_convert_par2_MPV_fixb.m
% 
% Output:
% Mq [6x6]
%   inertia matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-08 21:11
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Mq = S6PRRPPR3_inertiaJ_mdp_slag_vp(qJ, pkin, MDP)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(10,1),zeros(26,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6PRRPPR3_inertiaJ_mdp_slag_vp: qJ has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S6PRRPPR3_inertiaJ_mdp_slag_vp: pkin has to be [10x1] (double)');
assert(isreal(MDP) && all(size(MDP) == [26 1]), ...
  'S6PRRPPR3_inertiaJ_mdp_slag_vp: MDP has to be [26x1] (double)'); 

%% Symbolic Calculation
% From inertia_joint_joint_fixb_mdp_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-08 21:10:58
% EndTime: 2019-03-08 21:10:59
% DurationCPUTime: 0.33s
% Computational Cost: add. (247->112), mult. (468->150), div. (0->0), fcn. (408->8), ass. (0->50)
t112 = MDP(13) - MDP(18);
t123 = MDP(15) * pkin(8) + t112;
t122 = pkin(8) - qJ(5);
t109 = MDP(15) + MDP(19);
t94 = -pkin(3) - pkin(4);
t86 = sin(pkin(6));
t90 = sin(qJ(2));
t121 = t86 * t90;
t93 = cos(qJ(2));
t120 = t86 * t93;
t88 = sin(qJ(6));
t91 = cos(qJ(6));
t119 = t88 * t91;
t92 = cos(qJ(3));
t118 = t88 * t92;
t117 = t91 * t92;
t115 = cos(pkin(6));
t89 = sin(qJ(3));
t71 = -t92 * pkin(3) - t89 * qJ(4) - pkin(2);
t70 = t92 * pkin(4) - t71;
t114 = t70 * MDP(19);
t113 = pkin(8) ^ 2 * MDP(15);
t111 = MDP(14) - MDP(11);
t110 = MDP(14) + MDP(16);
t108 = 0.2e1 * t92;
t107 = MDP(21) * t119;
t105 = -pkin(3) * MDP(15) - MDP(12);
t104 = t94 * MDP(19) + MDP(17);
t102 = -MDP(22) * t91 + MDP(23) * t88;
t101 = t91 * MDP(25) - t88 * MDP(26);
t100 = -t88 * MDP(25) - t91 * MDP(26);
t99 = MDP(16) + t101;
t98 = t104 + t105;
t97 = -t88 * MDP(22) - t91 * MDP(23) + t100 * (-pkin(9) + t94);
t95 = qJ(4) ^ 2;
t87 = qJ(4) + pkin(5);
t85 = t92 ^ 2;
t84 = t91 ^ 2;
t83 = t89 ^ 2;
t82 = t88 ^ 2;
t73 = t122 * t92;
t72 = t122 * t89;
t68 = t115 * t89 + t92 * t121;
t67 = -t115 * t92 + t89 * t121;
t65 = t89 * pkin(5) + t92 * pkin(9) + t70;
t64 = t88 * t120 + t67 * t91;
t63 = t91 * t120 - t67 * t88;
t62 = t88 * t65 + t91 * t72;
t61 = t91 * t65 - t88 * t72;
t1 = [MDP(1) + t109 * (t86 ^ 2 * t93 ^ 2 + t67 ^ 2 + t68 ^ 2); (t67 * t72 + t68 * t73) * MDP(19) + (-t68 * t118 + t63 * t89) * MDP(25) + (-t68 * t117 - t64 * t89) * MDP(26) + (-t90 * MDP(4) + (-MDP(15) * t71 + t114 + MDP(3) + (MDP(10) + MDP(12) - MDP(17)) * t92 + (-MDP(11) + t110) * t89) * t93) * t86 + t123 * (t67 * t89 + t68 * t92); MDP(2) + t71 ^ 2 * MDP(15) + (t70 ^ 2 + t72 ^ 2 + t73 ^ 2) * MDP(19) + (t84 * MDP(20) - 0.2e1 * t107 + t113) * t85 + (MDP(24) + MDP(5) + t113) * t83 + (pkin(2) * MDP(10) - t71 * MDP(12) - t70 * MDP(17)) * t108 + (-0.2e1 * pkin(2) * MDP(11) - 0.2e1 * t71 * MDP(14) + 0.2e1 * t70 * MDP(16) + (MDP(6) + t102) * t108) * t89 + 0.2e1 * (-t72 * t89 - t73 * t92) * MDP(18) + 0.2e1 * (-t73 * t118 + t61 * t89) * MDP(25) + 0.2e1 * (-t73 * t117 - t62 * t89) * MDP(26) + 0.2e1 * (t83 + t85) * MDP(13) * pkin(8); (-MDP(10) + t98) * t67 + (t109 * qJ(4) + t111 + t99) * t68; t104 * t72 + (qJ(4) * MDP(19) + t99) * t73 + (-pkin(3) * MDP(13) - t94 * MDP(18) + MDP(7) + t97) * t89 + (MDP(8) + MDP(20) * t119 + (-t82 + t84) * MDP(21) + t100 * t87 + t112 * qJ(4)) * t92 + ((qJ(4) * MDP(15) + t111) * t92 + (-MDP(10) + t105) * t89) * pkin(8); MDP(9) + 0.2e1 * pkin(3) * MDP(12) + (pkin(3) ^ 2 + t95) * MDP(15) + 0.2e1 * t94 * MDP(17) + (t94 ^ 2 + t95) * MDP(19) + t82 * MDP(20) + 0.2e1 * t107 + 0.2e1 * t110 * qJ(4) + 0.2e1 * t101 * t87; t109 * t67; t72 * MDP(19) + (t100 + t123) * t89; t98; t109; MDP(19) * t120; -t92 * MDP(17) + t99 * t89 + t114; 0; 0; MDP(19); t63 * MDP(25) - t64 * MDP(26); t89 * MDP(24) + t61 * MDP(25) - t62 * MDP(26) + t102 * t92; t97; t100; t101; MDP(24);];
%% Postprocessing: Reshape Output
% From vec2symmat_6_matlab.m
res = [t1(1) t1(2) t1(4) t1(7) t1(11) t1(16); t1(2) t1(3) t1(5) t1(8) t1(12) t1(17); t1(4) t1(5) t1(6) t1(9) t1(13) t1(18); t1(7) t1(8) t1(9) t1(10) t1(14) t1(19); t1(11) t1(12) t1(13) t1(14) t1(15) t1(20); t1(16) t1(17) t1(18) t1(19) t1(20) t1(21);];
Mq  = res;
