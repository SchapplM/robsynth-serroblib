% Calculate joint inertia matrix for
% S6RPRPRR1
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% pkin [11x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d3,d5,d6,theta2,theta4]';
% MDP [27x1]
%   Minimal dynamic parameter vector (fixed base model)
%   see S6RPRPRR1_convert_par2_MPV_fixb.m
% 
% Output:
% Mq [6x6]
%   inertia matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-09 03:35
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Mq = S6RPRPRR1_inertiaJ_mdp_slag_vp(qJ, pkin, MDP)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(11,1),zeros(27,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPRPRR1_inertiaJ_mdp_slag_vp: qJ has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [11 1]), ...
  'S6RPRPRR1_inertiaJ_mdp_slag_vp: pkin has to be [11x1] (double)');
assert(isreal(MDP) && all(size(MDP) == [27 1]), ...
  'S6RPRPRR1_inertiaJ_mdp_slag_vp: MDP has to be [27x1] (double)'); 

%% Symbolic Calculation
% From inertia_joint_joint_fixb_mdp_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 03:35:20
% EndTime: 2019-03-09 03:35:21
% DurationCPUTime: 0.40s
% Computational Cost: add. (661->104), mult. (1194->157), div. (0->0), fcn. (1353->10), ass. (0->58)
t114 = sin(qJ(6));
t117 = cos(qJ(6));
t124 = t117 * MDP(26) - t114 * MDP(27);
t115 = sin(qJ(5));
t141 = cos(qJ(5));
t110 = sin(pkin(11));
t112 = cos(pkin(11));
t116 = sin(qJ(3));
t118 = cos(qJ(3));
t98 = -t110 * t116 + t112 * t118;
t99 = t110 * t118 + t112 * t116;
t86 = t115 * t98 + t141 * t99;
t83 = t86 * MDP(20);
t85 = t115 * t99 - t141 * t98;
t122 = -(MDP(19) + t124) * t85 - t83;
t133 = t114 * MDP(23) + t117 * MDP(24);
t113 = cos(pkin(10));
t104 = -t113 * pkin(1) - pkin(2);
t100 = -t118 * pkin(3) + t104;
t87 = -t98 * pkin(4) + t100;
t143 = 0.2e1 * t87;
t140 = pkin(3) * t110;
t111 = sin(pkin(10));
t102 = t111 * pkin(1) + pkin(7);
t134 = qJ(4) + t102;
t96 = t134 * t116;
t97 = t134 * t118;
t79 = -t110 * t96 + t112 * t97;
t78 = -t110 * t97 - t112 * t96;
t72 = -t99 * pkin(8) + t78;
t73 = t98 * pkin(8) + t79;
t68 = t115 * t73 - t141 * t72;
t65 = t68 * t114;
t139 = t68 * t117;
t138 = t114 * t117;
t137 = t85 * MDP(25);
t103 = t112 * pkin(3) + pkin(4);
t93 = t141 * t103 - t115 * t140;
t136 = t93 * MDP(19);
t94 = -t115 * t103 - t141 * t140;
t135 = t94 * MDP(20);
t130 = t118 * MDP(10);
t129 = MDP(22) * t138;
t108 = t114 ^ 2;
t128 = t108 * MDP(21) + MDP(18) + 0.2e1 * t129;
t127 = -pkin(5) * t86 - pkin(9) * t85;
t91 = -pkin(5) - t93;
t92 = pkin(9) - t94;
t126 = -t85 * t92 + t86 * t91;
t125 = MDP(23) * t117 - MDP(24) * t114;
t123 = -MDP(26) * t114 - MDP(27) * t117;
t109 = t117 ^ 2;
t69 = t115 * t72 + t141 * t73;
t121 = -t68 * MDP(19) - t69 * MDP(20) + ((-t108 + t109) * MDP(22) + MDP(21) * t138 + MDP(16)) * t86 + (-MDP(17) + t133) * t85;
t70 = t85 * pkin(5) - t86 * pkin(9) + t87;
t64 = t114 * t70 + t117 * t69;
t63 = -t114 * t69 + t117 * t70;
t1 = [MDP(1) - 0.2e1 * t104 * t130 + (t100 ^ 2 + t78 ^ 2 + t79 ^ 2) * MDP(13) + t83 * t143 + (t111 ^ 2 + t113 ^ 2) * MDP(4) * pkin(1) ^ 2 + (t109 * MDP(21) + MDP(14) - 0.2e1 * t129) * t86 ^ 2 + (0.2e1 * t104 * MDP(11) + MDP(5) * t116 + 0.2e1 * t118 * MDP(6)) * t116 + (MDP(19) * t143 + t137 + 0.2e1 * (-MDP(15) + t125) * t86) * t85 + 0.2e1 * (-t78 * t99 + t79 * t98) * MDP(12) + 0.2e1 * (t63 * t85 + t86 * t65) * MDP(26) + 0.2e1 * (t86 * t139 - t64 * t85) * MDP(27); (t78 * t98 + t79 * t99) * MDP(13); MDP(4) + (t98 ^ 2 + t99 ^ 2) * MDP(13); t116 * MDP(7) + t118 * MDP(8) + (t126 * t114 - t139) * MDP(26) + (t126 * t117 + t65) * MDP(27) + (-t116 * MDP(10) - t118 * MDP(11)) * t102 + ((t110 * t98 - t112 * t99) * MDP(12) + (t110 * t79 + t112 * t78) * MDP(13)) * pkin(3) + t121; t130 - t116 * MDP(11) + (t110 * t99 + t112 * t98) * MDP(13) * pkin(3) + t122; MDP(9) - 0.2e1 * t124 * t91 + (t110 ^ 2 + t112 ^ 2) * MDP(13) * pkin(3) ^ 2 + 0.2e1 * t136 + 0.2e1 * t135 + t128; t100 * MDP(13) - t122; 0; 0; MDP(13); (t127 * t114 - t139) * MDP(26) + (t127 * t117 + t65) * MDP(27) + t121; t122; t128 + t135 + t136 + t124 * (pkin(5) - t91); 0; 0.2e1 * pkin(5) * t124 + t128; t63 * MDP(26) - t64 * MDP(27) + t125 * t86 + t137; t123 * t86; t123 * t92 + t133; t124; t123 * pkin(9) + t133; MDP(25);];
%% Postprocessing: Reshape Output
% From vec2symmat_6_matlab.m
res = [t1(1) t1(2) t1(4) t1(7) t1(11) t1(16); t1(2) t1(3) t1(5) t1(8) t1(12) t1(17); t1(4) t1(5) t1(6) t1(9) t1(13) t1(18); t1(7) t1(8) t1(9) t1(10) t1(14) t1(19); t1(11) t1(12) t1(13) t1(14) t1(15) t1(20); t1(16) t1(17) t1(18) t1(19) t1(20) t1(21);];
Mq  = res;
